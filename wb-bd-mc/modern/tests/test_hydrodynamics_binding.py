"""Direct NumPy HI contract and dimensional analytical checks."""
import unittest
import numpy as np

from bdmc.hydrodynamics import (
    HydrodynamicSpec, HydrodynamicSimulation, StokesianSpec,
    mobility, solve_stokesian,
)


class HydrodynamicsBindingTests(unittest.TestCase):
    def test_literal_mobility_and_owned_arrays(self):
        q=np.array([[0.,0,0],[4,0,0]])
        s=HydrodynamicSpec(n_particles=2,viscosity=1/(6*np.pi),seed=810729)
        m=mobility(q,s)
        self.assertEqual(m.shape,(6,6))
        self.assertAlmostEqual(m[0,3],.359375,places=14)
        sim=HydrodynamicSimulation(s,q)
        q[:]=999
        first=sim.observe()
        self.assertEqual(first['positions'][0,0],0)
        first['positions'][:]=888
        first['unwrapped_positions'][:]=777
        self.assertEqual(sim.observe()['positions'][0,0],0)
        for key in ('positions','unwrapped_positions','forces','pair_virial'):
            self.assertEqual(first[key].dtype,np.float64)
            self.assertTrue(first[key].flags.c_contiguous)
            self.assertTrue(first[key].flags.owndata)
        self.assertIn('no periodic',sim.metadata['boundary'])

    def test_batch_replay_and_spec_ownership(self):
        q=np.array([[0.,0,0],[4,0,0]])
        s=HydrodynamicSpec(n_particles=2,seed=810729,dt=.002)
        a,b=HydrodynamicSimulation(s,q),HydrodynamicSimulation(s,q)
        s.dt=999
        copied=a.spec
        copied.dt=888
        out=a.sample(7,stride=3)
        np.testing.assert_array_equal(out['steps'],[0,3,6,7])
        np.testing.assert_allclose(out['times'],np.array([0,3,6,7])*.002,rtol=0,atol=0)
        self.assertEqual(out['positions'].shape,(4,2,3))
        b.advance(7)
        np.testing.assert_array_equal(a.observe()['positions'],b.observe()['positions'])
        out['positions'][:]=111
        out['unwrapped_positions'][:]=222
        np.testing.assert_array_equal(a.observe()['positions'],b.observe()['positions'])
        empty=a.sample(0,include_initial=False)
        self.assertEqual(empty['positions'].shape,(0,2,3))

    def test_free_space_force_cutoff_and_shear(self):
        q=np.array([[0.,0,0],[4,0,0]])
        s=HydrodynamicSpec(n_particles=2,thermal_energy=0,epsilon=6,sigma=8,seed=810729)
        frame=HydrodynamicSimulation(s,q).observe()
        np.testing.assert_allclose(frame['forces'],[[-.5625,0,0],[.5625,0,0]],rtol=0,atol=1e-15)
        self.assertEqual(frame['potential_energy'],1.25)
        # Large centre separation stays outside sigma; no periodic images.
        far=HydrodynamicSimulation(s,[[0,0,0],[100,0,0]]).observe()
        np.testing.assert_array_equal(far['forces'],np.zeros((2,3)))
        one=HydrodynamicSpec(thermal_energy=0,dt=.01,shear_rate=.3,
            external_forces=[[2,0,0]],background_velocity=[[.2,0,0]])
        sim=HydrodynamicSimulation(one,[[0,2,0]])
        sim.advance(1)
        self.assertAlmostEqual(sim.observe()['positions'][0,0],.01*(2/(6*np.pi)+.2+.6))

    def test_athermal_stokesian_response(self):
        radius,eta=.7,1.3
        s=StokesianSpec(radius=radius,viscosity=eta)
        f=np.array([[1.,2,3]])
        t=np.array([[2.,-1,.5]])
        E=np.diag([1.,-.3,-.7]);E[0,1]=E[1,0]=.2
        r=solve_stokesian([[0,0,0]],s,forces=f,torques=t,strain=E)
        np.testing.assert_allclose(r['translation'],f/(6*np.pi*eta*radius),rtol=2e-14,atol=2e-14)
        np.testing.assert_allclose(r['rotation'],t/(8*np.pi*eta*radius**3),rtol=2e-14,atol=2e-14)
        np.testing.assert_allclose(r['stresslet'][0],(20*np.pi/3)*eta*radius**3*E,rtol=2e-14,atol=2e-14)
        self.assertLess(r['relative_residual'],1e-9)
        self.assertIn('athermal',r['metadata']['thermal'])
        self.assertEqual(r['stresslet'].shape,(1,3,3))
        self.assertTrue(r['stresslet'].flags.owndata)

    def test_rejections_and_state_transaction(self):
        with self.assertRaises(ValueError): HydrodynamicSpec(n_particles=129)
        with self.assertRaises(ValueError): HydrodynamicSpec(seed=2**32)
        with self.assertRaises(ValueError): HydrodynamicSpec(mobility='periodic_rpy')
        with self.assertRaises(ValueError): HydrodynamicSpec(external_forces=[1,2,3])
        with self.assertRaises(ValueError): HydrodynamicSimulation(HydrodynamicSpec(),[0,0,0])
        with self.assertRaises(ValueError): mobility([[0,0,np.nan]],HydrodynamicSpec())
        with self.assertRaises(ValueError): solve_stokesian([[0,0,0]],strain=np.eye(3))
        with self.assertRaises(ValueError): solve_stokesian([[0,0,0]],strain=np.zeros((1,9)))
        with self.assertRaises(ValueError): solve_stokesian([[0,0,0],[2.00001,0,0]])
        s=HydrodynamicSpec(n_particles=2,viscosity=1/(6*np.pi),dt=.1,
                           thermal_energy=1,seed=3519,mobility='oseen')
        sim=HydrodynamicSimulation(s,[[0,0,0],[2.000001,0,0]])
        initial=sim.observe()
        for _ in range(8):
            with self.assertRaises(ValueError): sim.advance(1)
            np.testing.assert_array_equal(sim.observe()['positions'],initial['positions'])
            self.assertEqual(sim.observe()['step'],0)
        with self.assertRaises(ValueError): sim.sample(2**64-1,stride=1)


if __name__=='__main__':
    unittest.main()
