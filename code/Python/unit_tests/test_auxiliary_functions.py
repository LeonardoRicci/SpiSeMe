import unittest
import numpy as np
import sys
sys.path.insert(0, "../spiSeMe")
import spiSeMe_distribution_test
import spiSeMe_event_autocorrelation
import spiSeMe_event_cross_correlation

class TestStringMethods(unittest.TestCase):
	def setUp(self):
		self.data_A = np.loadtxt('test_data_A.dat')
		self.data_B = np.loadtxt('test_data_B.dat')
		self.data_Z = [1, 2, 3, 4, 5, 6]

	def test_aux_distrtest_arg(self):
		print("Distribution Test: Intentionally providing bad arguments...")
		p, d = spiSeMe_distribution_test.spiSeMe_distribution_test(self.data_Z, self.data_A)
		self.assertFalse(p)
		self.assertFalse(d)
		p, d = spiSeMe_distribution_test.spiSeMe_distribution_test(self.data_A, self.data_Z)
		self.assertFalse(p)
		self.assertFalse(d)
	def test_aux_distrdtest_run(self):
		print("Distribution Test: Checking output formats...")
		self.p, self.d = spiSeMe_distribution_test.spiSeMe_distribution_test(self.data_A, self.data_A)
		self.assertEqual(self.p, 1.0)
		self.assertEqual(self.d, 0.0)
		self.p, self.d = spiSeMe_distribution_test.spiSeMe_distribution_test(self.data_A, self.data_B)
		self.assertAlmostEqual(self.p, 0.5, delta=0.1)
		self.assertAlmostEqual(self.d, 0.04, delta=0.01)

	def test_aux_autocorr_arg(self):
		print("Autocorrelation: Intentionally providing bad arguments...")
		A, lags = spiSeMe_event_autocorrelation.spiSeMe_event_autocorrelation(self.data_Z, 0.1, 1.0)
		self.assertFalse(A)
		self.assertFalse(lags)
		A, lags = spiSeMe_event_autocorrelation.spiSeMe_event_autocorrelation(self.data_A, -0.1, 1.0)
		self.assertFalse(A)
		self.assertFalse(lags)
		A, lags = spiSeMe_event_autocorrelation.spiSeMe_event_autocorrelation(self.data_A, 0.1, -1.0)
		self.assertFalse(A)
		self.assertFalse(lags)
		A, lags = spiSeMe_event_autocorrelation.spiSeMe_event_autocorrelation(self.data_A, 1.0, 0.1)
		self.assertFalse(A)
		self.assertFalse(lags)
		A, lags = spiSeMe_event_autocorrelation.spiSeMe_event_autocorrelation(self.data_A, 'x', 1.0)
		self.assertFalse(A)
		self.assertFalse(lags)
		A, lags = spiSeMe_event_autocorrelation.spiSeMe_event_autocorrelation(self.data_A, 0.1, 'x')
		self.assertFalse(A)
		self.assertFalse(lags)
	def test_aux_autocorr_run(self):
		print("Autocorrelation: Checking output formats...")
		A, lags = spiSeMe_event_autocorrelation.spiSeMe_event_autocorrelation(self.data_A, 0.1, 1.0)
		self.assertTrue(isinstance(A,np.ndarray))
		self.assertTrue(isinstance(lags,np.ndarray))
		self.assertEqual(A.ndim, 1)
		self.assertEqual(lags.ndim, 1)
		self.assertEqual(A.size, lags.size)
		self.assertEqual(A.size, 10)

	def test_aux_crosscorr_arg(self):
		print("Cross Correlation: Intentionally providing bad arguments...")
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_Z, self.data_A, 0.1, 1.0)
		self.assertFalse(C)
		self.assertFalse(lags)
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_A, self.data_Z, 0.1, 1.0)
		self.assertFalse(C)
		self.assertFalse(lags)
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_A, self.data_B, 'x', 1.0)
		self.assertFalse(C)
		self.assertFalse(lags)
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_A, self.data_B, 0.1, 'x')
		self.assertFalse(C)
		self.assertFalse(lags)
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_A, self.data_B, -0.1, 1.0)
		self.assertFalse(C)
		self.assertFalse(lags)
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_A, self.data_B, 0.1, -1.0)
		self.assertFalse(C)
		self.assertFalse(lags)
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_A, self.data_B, 1.0, 0.1)
		self.assertFalse(C)
		self.assertFalse(lags)
	def test_aux_crosscorr_run(self):
		print("Cross Correlation: Checking output formats...")
		C, lags = spiSeMe_event_cross_correlation.spiSeMe_event_cross_correlation(self.data_A, self.data_B, 0.1, 1.0)
		self.assertTrue(isinstance(C,np.ndarray))
		self.assertTrue(isinstance(lags,np.ndarray))
		self.assertEqual(C.ndim, 1)
		self.assertEqual(lags.ndim, 1)
		self.assertEqual(C.size, lags.size)
		self.assertEqual(C.size, 21)


if __name__ == '__main__':
    unittest.main()
