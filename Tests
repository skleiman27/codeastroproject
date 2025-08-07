import unittest
import numpy as np
import astropy.units as u

class Exoplanet(object):
    def __init__(self, rad, a, rad_unit="R_J", a_unit="AU"):
        if rad_unit == "R_J":
            self.radius = rad * u.R_jup
        elif rad_unit == "R_earth":
            self.radius = rad * u.R_earth
        elif rad_unit == "m":
            self.radius = rad * u.m
        else:
            raise ValueError("rad_unit must be 'R_J', 'R_earth', or 'm'.")
        if a_unit == "AU":
            self.a = a * u.au
        elif a_unit == "m":
            self.a = a * u.m
        else:
            raise ValueError("a_unit must be either 'AU' or 'm'.")

class Star(object):
    def __init__(self, rad):
        self.radius = rad * u.R_sun

class LightCurveExoplanet(object):
    def __init__(self, planet, star, ticksinper=100, noise=0.001, numper=1):
        self.ticksinper = ticksinper
        self.length = self.ticksinper * numper
        self.depth = float(((planet.radius / star.radius).to(''))**2)
        self.duration = float((star.radius / (planet.a * 2 * np.pi)).to('') * self.ticksinper)
        self.noise = noise
        self.location = np.random.randint(0, ticksinper)
        self.flux = np.ones(self.length) + np.random.normal(loc=0, scale=self.noise, size=self.length)
        self.per = 0
        self.numper = numper

    def new_transit(self):
        self.depth = np.random.normal(loc=self.depth, scale=0.001, size=1)
        self.location = self.location + self.ticksinper

    def plot_transit(self, phase_flag=False):
        lowerbound = int(self.location - self.duration / 2)
        upperbound = int(self.location + self.duration / 2)
        if lowerbound < 0:
            self.ub = self.length + lowerbound
            self.lb = upperbound
            self.flux[0:self.lb] -= self.depth
            self.flux[self.ub:self.length] -= self.depth
            self.overflow_flag = True
        elif upperbound > self.length:
            self.lb = int(upperbound - self.length)
            self.ub = int(lowerbound)
            self.flux[0:self.lb] -= self.depth
            self.flux[self.ub:self.length] -= self.depth
            self.overflow_flag = True
        else:
            self.ub = int(upperbound)
            self.lb = int(lowerbound)
            self.flux[self.lb:self.ub] -= self.depth
            self.overflow_flag = False
        if self.lb == 0:
            self.lb = 1
        if self.ub == self.length:
            self.ub = self.ub - 1

class TestTransitSimulation(unittest.TestCase):
    def setUp(self):
        self.planet = Exoplanet(rad=1, a=0.05)
        self.star = Star(rad=1)
        self.lc = LightCurveExoplanet(self.planet, self.star, ticksinper=100, noise=0, numper=1)

    def test_depth_calculation(self):
        expected_depth = float(((self.planet.radius / self.star.radius).to(''))**2)
        self.assertAlmostEqual(self.lc.depth, expected_depth, places=5)

    def test_duration_calculation(self):
        expected_duration = float((self.star.radius / (self.planet.a * 2 * np.pi)).to('') * 100)
        self.assertAlmostEqual(self.lc.duration, expected_duration, places=5)

    def test_flux_array_shape(self):
        self.assertEqual(len(self.lc.flux), 100)

    def test_noise_application(self):
        noisy_lc = LightCurveExoplanet(self.planet, self.star, noise=0.005)
        std = np.std(noisy_lc.flux - 1)
        self.assertTrue(0.003 < std < 0.007)

    def test_plot_transit_bounds(self):
        self.lc.plot_transit()
        self.assertTrue(0 <= self.lc.lb < self.lc.ub <= self.lc.length)

    def test_transit_overflow_flag(self):
        self.lc.location = 2
        self.lc.plot_transit()
        self.assertTrue(self.lc.overflow_flag)

    def test_transit_no_overflow(self):
        self.lc.location = 50
        self.lc.plot_transit()
        self.assertFalse(self.lc.overflow_flag)

    def test_new_transit_period_increments(self):
        loc_before = self.lc.location
        self.lc.new_transit()
        self.assertEqual(self.lc.location, loc_before + 100)

    def test_phase_folded_plot_does_not_crash(self):
        try:
            self.lc.plot_transit(phase_flag=True)
        except Exception as e:
            self.fail(f"Phase-folded plot raised an exception: {e}")

if __name__ == '__main__':
    unittest.main()
