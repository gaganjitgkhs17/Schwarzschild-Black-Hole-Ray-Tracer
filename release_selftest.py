"""Portable standalone smoke tests. Run: python release_selftest.py"""
import math
import unittest
import numpy as np
import bh_engine_v010s as e

class ReleaseSmokeTests(unittest.TestCase):
    def test_constants(self):
        self.assertGreater(e.M, 0)
        self.assertGreater(e.R_OBSERVER, 2 * e.M)
        self.assertGreaterEqual(e.DISK_INNER, 6 * e.M)
        self.assertGreater(e.DISK_OUTER, e.DISK_INNER)

    def test_camera_impact_parameter(self):
        x = np.array([-4., 0., 5.]); y = np.array([3., 0., -2.])
        rho = np.hypot(x, y)
        expected = (e.R_OBSERVER / np.sqrt(1-2*e.M/e.R_OBSERVER)
                    *rho/np.hypot(e.FOCAL_LENGTH,rho))
        np.testing.assert_allclose(e.calculate_impact_parameter(x,y), expected,
                                   rtol=0, atol=1e-12)

    def test_critical_impact_parameter(self):
        bc = 3 * math.sqrt(3) * e.M
        angle = math.asin(bc * math.sqrt(1-2*e.M/e.R_OBSERVER)/e.R_OBSERVER)
        x = np.array([e.FOCAL_LENGTH * math.tan(angle)])
        np.testing.assert_allclose(e.calculate_impact_parameter(x,np.zeros(1)),
                                   bc, atol=1e-12, rtol=0)

    def test_ray_radiance(self):
        x = np.array([-8., -4., 0., 8.]); y = np.array([-3., 2., 0., 3.])
        rgb, hits = e.trace_radiance(x, y)
        self.assertEqual(rgb.shape, (4, 3))
        self.assertEqual(hits.shape, (4,))
        self.assertTrue(np.isfinite(rgb).all())
        self.assertTrue((rgb >= 0).all())

    def test_tile_independence(self):
        a, ah = e.render_linear(12, 2, 1, .008, verbose=False)
        b, bh = e.render_linear(12, 2, 5, .008, verbose=False)
        np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(ah,bh)

    def test_supersampling_alignment(self):
        a, ah = e.render_linear(12, 2, 3, .008, verbose=False)
        b, bh = e.render_linear(24, 1, 5, .008, verbose=False)
        np.testing.assert_allclose(a,b.reshape(12,2,12,2,3).mean((1,3)),atol=1e-12,rtol=0)
        np.testing.assert_array_equal(ah,bh.reshape(12,2,12,2).mean((1,3)))

if __name__ == '__main__':
    unittest.main(verbosity=2)
