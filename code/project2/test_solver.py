"""Independent operator, boundary, constitutive, and nullspace regression tests."""
import unittest
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
import sympy as sy
from run_all import assemble, metrics, forcing, layered, graph_neumann, neumann_solve


class SolverTests(unittest.TestCase):
    def test_constant_operator_against_kronecker(self):
        N=6; m=N-1
        A,b,_=assemble(N,lambda x,y:1,lambda x,y:0,lambda x,y:0)
        T=sp.diags([-np.ones(m-1),2*np.ones(m),-np.ones(m-1)],[-1,0,1])
        oracle=N*N*(sp.kron(T,sp.eye(m))+sp.kron(sp.eye(m),T))
        self.assertLess(sp.linalg.norm(A-oracle),1e-13)
        np.testing.assert_array_equal(b,0)

    def test_zero_solution_and_integer_inputs(self):
        A,b,_=assemble(4,lambda x,y:1,lambda x,y:0,lambda x,y:0)
        self.assertEqual(metrics(A,b,splu(A).solve(b)),(0.,0.))

    def test_affine_nonzero_boundary(self):
        A,b,d=assemble(8,lambda x,y:1+x,lambda x,y:-1,lambda x,y:x+2*y)
        self.assertEqual((A[0,0],A[0,1],A[0,7],b[0]),(288,-72,-76,25))
        ref=(d[0]+2*d[1])[1:-1,1:-1].ravel()
        np.testing.assert_allclose(splu(A).solve(b),ref,atol=2e-14,rtol=0)

    def test_callable_forcing_against_symbolic_derivative(self):
        x,y=sy.symbols('x y'); u=sy.sin(sy.pi*x)*sy.sin(sy.pi*y)
        oracle=sy.lambdify((x,y),-sy.diff((1+x)*sy.diff(u,x),x)-sy.diff((1+x)*sy.diff(u,y),y),'numpy')
        xx=np.array([.13,.37,.81]);yy=np.array([.21,.55,.72])
        np.testing.assert_allclose(forcing(xx,yy),oracle(xx,yy),rtol=1e-14)

    def test_harmonic_interface_exact_and_arithmetic_falsified(self):
        errors={}
        for rule in ['harmonic','arithmetic']:
            A,b,d=assemble(15,lambda x,y:np.where(x<.5,1.,10.),lambda x,y:0,layered,rule)
            errors[rule]=np.max(abs(splu(A).solve(b)-layered(d[0][1:-1,1:-1],d[1][1:-1,1:-1]).ravel()))
            self.assertAlmostEqual(d[3][7,7],20/11 if rule=='harmonic' else 11/2)
        self.assertLess(errors['harmonic'],1e-12)
        self.assertGreater(errors['arithmetic'],1e-3)

    def test_neumann_nullspace_and_mean_gauge(self):
        A,b,ref=graph_neumann(8); u,lag,defect=neumann_solve(A,b)
        np.testing.assert_allclose(A@np.ones(len(b)),0,atol=1e-12)
        self.assertLess(abs(np.mean(u)),1e-13)
        self.assertLess(np.linalg.norm(A@u-b),1e-10)

    def test_incompatible_neumann_is_rejected(self):
        A,b,_=graph_neumann(8)
        with self.assertRaises(ValueError): neumann_solve(A,b+1)

    def test_smallest_grid(self):
        A,b,_=assemble(2,lambda x,y:1,lambda x,y:0,lambda x,y:1)
        np.testing.assert_allclose(splu(A).solve(b),[1.])

    def test_invalid_grid(self):
        for N in [0,1,3.5]:
            with self.assertRaises(ValueError): assemble(N,lambda x,y:1,lambda x,y:0,lambda x,y:0)

    def test_invalid_coefficient_and_data(self):
        for a in [0,-1,np.nan,np.inf]:
            with self.assertRaises(ValueError): assemble(4,lambda x,y:a,lambda x,y:0,lambda x,y:0)
        with self.assertRaises(ValueError): assemble(4,lambda x,y:1,lambda x,y:np.nan,lambda x,y:0)


if __name__=='__main__': unittest.main(verbosity=2)
