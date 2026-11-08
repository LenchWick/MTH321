"""Independent analytic special-case checks; no full Monte Carlo run required."""
import unittest
import numpy as np
from kernels import gbm, sv


class AnalyticChecks(unittest.TestCase):
    def test_exact_transition_composition(self):
        dw = np.array([[.2, -.1], [-.3, .4], [.1, -.2]])
        h = .1
        oracle = 100*np.exp((.2-.09/2)*3*h+.3*dw.sum(axis=0))
        np.testing.assert_allclose(gbm(dw, h, 'exact'), oracle, rtol=1e-13)

    def test_zero_diffusion_deterministic_recurrence(self):
        dw = np.array([[.2], [-.3], [.1]])
        for method in ['EM', 'Milstein']:
            np.testing.assert_allclose(gbm(dw, .1, method, sigma=0),
                                       [100*(1+.2*.1)**3], rtol=1e-13)

    def test_one_step_milstein_independent_formula(self):
        h, w = .25, -.2
        oracle = 100*(1+.2*h+.3*w+.045*(w*w-h))
        np.testing.assert_allclose(gbm(np.array([[w]]), h, 'Milstein'), [oracle])

    def test_frozen_factor_exact_gbm(self):
        dw = np.array([[.2, -.1], [-.3, .4], [.1, -.2]])
        for balanced in [False, True]:
            price, factor = sv(dw, dw, .1, mu=.2, theta=0, xi=0,
                               balanced=balanced)
            oracle = 100*np.exp(.155*.3+.3*dw.sum(axis=0))
            np.testing.assert_allclose(price, oracle, rtol=1e-13)
            np.testing.assert_equal(factor, np.zeros(2))

    def test_no_shocks_log_scheme(self):
        dw = np.zeros((10, 2))
        price, _ = sv(dw, dw, .1, mu=.05, theta=0, xi=0)
        np.testing.assert_allclose(price, 100*np.exp(.05-.3**2/2), rtol=1e-13)

    def test_em_negative_tail_and_milstein_stress(self):
        dw = np.array([[-1.5]])
        self.assertLess(gbm(dw, 1, 'EM', mu=.05, sigma=1.5)[0], 0)
        self.assertLess(gbm(np.array([[-1/1.5]]), 1, 'Milstein',
                            mu=.05, sigma=1.5)[0], 0)

    def test_milstein_minimum_condition(self):
        # Evaluate exactly at the quadratic minimizer; no random sampling.
        for h in [.125, .25, .5, 1.]:
            oracle = 100*(.5+(.05-.5*1.5**2)*h)
            actual = gbm(np.array([[-1/1.5]]), h, 'Milstein', mu=.05, sigma=1.5)
            np.testing.assert_allclose(actual, [oracle], atol=1e-12)

    def test_balanced_stationary_variance_identity(self):
        kappa, xi, theta = 2., .6, -.2
        v = xi**2*(1+theta**2)/(2*kappa-xi**2)
        for h in [.125, .25, .5, 1.]:
            a = np.exp(-kappa*h)
            r = xi**2*(1-a*a)/(2*kappa)
            self.assertAlmostEqual((a*a+r)*v+r*(1+theta**2), v, places=13)


if __name__ == '__main__':
    unittest.main(verbosity=2)
