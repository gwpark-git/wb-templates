#!/usr/bin/env python3
"""Optional fresh non-cake comparison against an explicitly supplied public mBLA checkout.

Requires NumPy/SciPy. Does not implement cake or establish a historical tolerance.
The caller must review the supplied Python source before importing it.
"""
import argparse
import contextlib
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--python-source', required=True, type=Path, help='reviewed checkout scripts directory')
    p.add_argument('--exe', required=True, type=Path)
    p.add_argument('--radial-grid', choices=['native','cpp-log'], default='native')
    p.add_argument('--inlet-velocity-order', type=int, choices=(1, 2), default=2,
                   help='C++ inlet axial velocity reconstruction order with prescribed inlet pressure (default: 2)')
    p.add_argument('--nr', type=int, default=120)
    p.add_argument('--work', required=True, type=Path)
    a = p.parse_args()
    import numpy as np
    import scipy
    required = ['sol_solvent.py', 'sol_CT.py', 'sol_GT.py', 'transport_properties_SPHS.py',
                'osmotic_pressure_CS.py', 'membrane_geometry_functions.py', 'aux_functions.py']
    identity = {name: hashlib.sha256((a.python_source/name).read_bytes()).hexdigest() for name in required}
    sys.path.insert(0, str(a.python_source.resolve()))
    import sol_solvent as PS
    import sol_CT as CT
    import sol_GT as GT
    import transport_properties_SPHS as TP
    import osmotic_pressure_CS as OS
    a.work = a.work.resolve()
    a.work.mkdir(parents=True, exist_ok=True)
    config = (Path(__file__).resolve().parents[1]/'examples/noncake.inp').read_text()
    if not 3 <= a.nr <= 1000: raise ValueError('comparison Nr must lie in [3,1000]')
    config = config.replace('ref_DTP=5000', 'ref_DTP=1000').replace('Nr=120', 'Nr='+str(a.nr))
    config += f'inlet_velocity_order={a.inlet_velocity_order}\n'
    (a.work/'case.inp').write_text(config)
    with (a.work/'cpp.log').open('w') as log:
        subprocess.run([str(a.exe.resolve()), 'case.inp'], cwd=a.work,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=35)
    cpp = np.loadtxt(a.work/'result.dat')
    # Explicit common physical inputs and pressure convention. The GSL constant
    # is retained for this comparison; do not silently substitute exact SI k_B.
    nz, nr = 25, a.nr
    R, L, lp, eta, temp, particle, feed = 0.0005, 0.5, 6.7e-10, 0.001, 293.15, 1e-8, 0.001
    dlp, tmp, pout = 130., 1000., 101325.
    pin, pper = pout+dlp, pout+dlp/2-tmp
    k=4*math.sqrt(L*L*lp*eta/R**3)
    ustar=R*R*dlp/(4*eta*L)
    pre={'k':k,'R':R,'L':L,'Lp':lp,'eta0':eta,'membrane_geometry':'HF','lam1':4.,'lam2':4.,
         'define_permeability':'Lp','h':R/2,'kappa_Darcy':lp*eta*R*math.log(1.5),
         'BC_inlet':'pressure','DLP':dlp,'phi_freeze':0.494,'Nz':nz}
    with contextlib.redirect_stdout(io.StringIO()):
        ps=PS.get_cond(pre,pin,pout,pper,ustar)
        ct=CT.get_cond(ps,feed,particle,particle,4*math.pi*particle**3/3,
                       1.3806504e-23*temp,L/(nz-1),0.)
        cond=GT.get_cond(ct,nr,0.1)
    z=np.linspace(0,1,nz)
    y=GT.gen_y_div_R_arr(cond) if a.radial_grid == 'native' else np.logspace(-6,0,nr)-1e-6
    relative=np.ones(nz)
    gp,gm=np.zeros(nz),np.zeros(nz)
    fields=[[np.zeros(nr) for _ in range(nz)] for _ in range(3)]

    def pressure_state(profile):
        pi=np.asarray(OS.Pi_CS(profile*feed,cond))/dlp
        CT.gen_gpm_arr(1.,z,pi,k,gp);CT.gen_gpm_arr(-1.,z,pi,k,gm)
        cond['Gk']=CT.get_Gk_boost(k,1/(nz-1),gp[-1],gm[-1],cond['denom_Gk_BC_specific'])
        cond['Bp']=CT.get_Bpm_conv(1.,cond);cond['Bm']=CT.get_Bpm_conv(-1.,cond)
        return pi

    start=time.monotonic()
    diagnostic=''
    for iteration in range(400):
        pi=pressure_state(relative)
        candidate=relative.copy()
        capture=io.StringIO()
        with contextlib.redirect_stdout(capture):
            GT.gen_new_phiw_div_phib_arr(1,candidate,cond,TP.Dc_short_div_D0_SPHS,TP.eta_div_eta0_SPHS,
                z,relative,pi,0.1,gp,gm,y,fields[0],fields[1],fields[2])
        diagnostic=capture.getvalue()
        residual=float(np.max(np.abs(candidate-relative)))
        relative=candidate
        if not np.isfinite(relative).all(): raise RuntimeError('nonfinite Python iterate')
        if residual < 1e-6 and not diagnostic: break
        if time.monotonic()-start > 50: raise TimeoutError('bounded Python comparison exceeded 50 seconds')
    if residual >= 1e-6 or diagnostic: raise RuntimeError('Python comparison did not converge cake-free')
    pi=pressure_state(relative)
    python_p=np.array([GT.get_P_conv(0.,zi,cond,gp[i],gm[i])*dlp for i,zi in enumerate(z)])
    python_v=np.array([GT.get_v_conv(1.,zi,pi[i],cond,gp[i],gm[i])*cond['vw0'] for i,zi in enumerate(z)])
    python_u=np.array([GT.get_uZ_out(zi,k,cond['Bp'],cond['Bm'],gp[i],gm[i])*ustar for i,zi in enumerate(z)])
    assert np.max(cpp[:,7]) == 0 and np.max(relative*feed) < 0.494
    observed={
        'wall_concentration_max_relative_difference':float(np.max(np.abs(cpp[:,1]/(relative*feed)-1))),
        'wall_permeate_max_relative_difference':float(np.max(np.abs(cpp[:,3]/python_v-1))),
        'axial_center_factor_max_relative_difference':float(np.max(np.abs(cpp[:,4]/python_u-1))),
        'pressure_max_absolute_difference_Pa':float(np.max(np.abs(cpp[:,2]-python_p))),
    }
    # Prospectively chosen 1% current comparison allowance, separate from the
    # researcher's established historical verification with unrecovered tolerance.
    accepted=all(observed[name] < 0.01 for name in observed if 'relative_difference' in name)
    result={'status':'passed' if accepted else 'not_accepted','scope':'new non-cake comparison only',
        'python_module_sha256':identity,'numpy':np.__version__,'scipy':scipy.__version__,
        'cpp_inlet_velocity_order':a.inlet_velocity_order,
        'radial_grid':a.radial_grid,'Nr':nr,'python_iterations':iteration+1,'python_seconds':time.monotonic()-start,
        'python_final_max_relaxed_increment':residual,'observed':observed,
        'new_acceptance':'<1% maximum relative difference in wall concentration, permeate and axial center factor',
        'differences':'Python radial grid is explicitly named; C++ uses logarithmic grid. Python integral axial pressure differs from C++ finite-difference pressure.',
        'history':'Researcher-established prior non-cake agreement remains valid reported evidence; historical conditions/tolerance not recovered here.',
        'cake':'Python source has no cake extension; no cake parity claim'}
    (a.work/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savetxt(a.work/'python-noncake.dat',np.column_stack([z*L,relative*feed,python_p,python_v,python_u]))
    print(json.dumps(result,indent=2))
    return 0 if accepted else 1


if __name__ == '__main__':
    sys.exit(main())
