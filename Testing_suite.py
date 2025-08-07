import unittest
import numpy as np
import astropy.units as u
import lcEnhance.LCE as lce

class TestTransitSimulation(unittest.TestCase):
    def setUp(self):
        self.planet = lce.Exoplanet(1, 0.05)
        self.star = lce.Star(1)
        self.lc = lce.LightCurveExoplanet(self.planet, self.star, ticksinper=100, noise=0, numper=1)

    def test_depth_calculation(self):
        expected_depth = float(((self.planet.radius / self.star.radius).to(''))**2)
        print(self.lc.depth)
        self.assertAlmostEqual(self.lc.depth, expected_depth, places=5)

    def test_duration_calculation(self):
        expected_duration = float((self.star.radius / (self.planet.a * 2 * np.pi)).to('') * 100)
        self.assertAlmostEqual(self.lc.duration, expected_duration, places=5)
        self.lc.plot_transit()

    def test_flux_array_shape(self):
        self.assertEqual(len(self.lc.fluxOG), 100)

    def test_noise_application(self):
        noisy_lc = lce.LightCurveExoplanet(self.planet, self.star, noise=0.005)
        std = np.std(noisy_lc.fluxOG - 1)
        self.assertTrue(0.003 < std < 0.007)

    def test_plot_transit_bounds(self):
        self.lc.plot_transit()
        self.assertTrue(0 <= self.lc.lb < self.lc.ub <= self.lc.length)

    def test_phase_folded_plot_does_not_crash(self):
        try:
            self.lc.plot_transit(phase_flag=True)
        except Exception as e:
            self.fail(f"Phase-folded plot raised an exception: {e}")

if __name__ == '__main__':
    unittest.main()
