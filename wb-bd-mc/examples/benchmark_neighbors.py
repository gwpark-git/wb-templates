"""Bounded CPU search comparison; accepted pairs are never written to disk."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

import numpy as np


def worker(n,backend,repeats):
    from bdmc import _core
    L=(n/.3)**(1/3)
    positions=np.random.default_rng(1900+n).uniform(0,L,size=(n,3))
    times=[]
    for _ in range(repeats):
        start=time.perf_counter()
        result=_core.neighbor_pairs(positions,[L]*3,1.,backend=backend)
        times.append(time.perf_counter()-start)
    digest=hashlib.sha256()
    for key in ('pair_indices','distances','vectors'):
        digest.update(result[key].tobytes())
    return {'n':n,'backend':backend,'density':.3,'cutoff':1.,'box_length':L,
            'repeats':repeats,'seconds':times,'median_seconds':statistics.median(times),
            'accepted_pairs':len(result['distances']),'output_sha256':digest.hexdigest(),
            'stats':dict(result['stats']),
            'process_peak_rss_native_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'rss_note':'whole child process including Python/NumPy; Linux native units KiB, macOS bytes',
            'core_binary_sha256':hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--sizes',nargs='+',type=int,default=[256,1024,4096])
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--worker',choices=['cell_list','brute_force'])
    parser.add_argument('--n',type=int)
    args=parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.n,args.worker,args.repeats)))
        raise SystemExit(0)
    if args.output is None or args.output.exists():
        parser.error('a fresh --output JSON path is required')
    if args.repeats<1 or args.repeats>5 or any(n<2 or n>8192 for n in args.sizes) or len(args.sizes)>5:
        parser.error('bounded benchmark: 1..5 repeats, <=5 sizes, 2..8192 particles')
    records=[]
    for n in args.sizes:
        pair=[]
        for backend in ('cell_list','brute_force'):
            process=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',backend,
                '--n',str(n),'--repeats',str(args.repeats)],capture_output=True,text=True,check=True,timeout=60)
            pair.append(json.loads(process.stdout))
        if pair[0]['output_sha256']!=pair[1]['output_sha256']:
            raise AssertionError(f'pair geometry mismatch at N={n}')
        records.append({'n':n,'same_pair_geometry':True,'cell_list':pair[0],'brute_force':pair[1],
                        'measured_speed_ratio_brute_over_cells':pair[1]['median_seconds']/pair[0]['median_seconds'],
                        'legacy_3d_pair_payload_lower_bound_bytes':32*n*n})
    report={'schema':'bdmc.neighbor-benchmark.v1','platform':platform.platform(),'records':records,
            'limits':['CPU toy search microbenchmark; not whole-simulation or GPU/MPI speedup',
                      'uniform independent positions at fixed density; network/dense regimes may differ',
                      'both algorithms sort/copy accepted output; timing includes binding overhead',
                      'legacy lower bound is static estimate, not measured legacy allocation']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    for row in records:
        print(row['n'],'pair match',row['same_pair_geometry'],'speed ratio',row['measured_speed_ratio_brute_over_cells'])
