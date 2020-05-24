import unittest
import numpy as np
import sys
sys.path.insert(0, "../spiSeMe")
import spiSeMe_distribution_test
import spiSeMe_surrogate_dither
import spiSeMe_surrogate_iaaft
import spiSeMe_surrogate_jodi
import spiSeMe_surrogate_sa

class TestStringMethods(unittest.TestCase):
	def setUp(self):
		self.data_A = np.loadtxt('test_data_A.dat')
		self.data_Z = [1, 2, 3, 4, 5, 6]

	def test_surrogate_jodi_arg(self):
		print("JODI: Intentionally providing bad arguments...")
		self.assertFalse(spiSeMe_surrogate_jodi.spiSeMe_surrogate_jodi(self.data_Z))
		self.assertFalse(spiSeMe_surrogate_jodi.spiSeMe_surrogate_jodi(self.data_A, 'x'))
		self.assertFalse(spiSeMe_surrogate_jodi.spiSeMe_surrogate_jodi(self.data_A, 1, 'x'))
	def test_surrogate_jodi_run(self):
		print("JODI: Checking output formats...")
		self.M = 3
		self.data_out = spiSeMe_surrogate_jodi.spiSeMe_surrogate_jodi(self.data_A, self.M, False)
		self.assertTrue(isinstance(self.data_out,np.ndarray))
		self.assertEqual(self.data_out.ndim, 2)
		self.assertEqual(self.data_out.size, self.M * self.data_A.size)
		self.assertEqual(self.data_out[0].size, self.data_A.size)
		self.p, self.d = spiSeMe_distribution_test.spiSeMe_distribution_test(self.data_A, self.data_out[0])
		self.assertEqual(self.p, 1)

	def test_surrogate_iaaft_arg(self):
		print("IAAFT: Intentionally providing bad arguments...")
		self.assertFalse(spiSeMe_surrogate_iaaft.spiSeMe_surrogate_iaaft(self.data_Z))
		self.assertFalse(spiSeMe_surrogate_iaaft.spiSeMe_surrogate_iaaft(self.data_A, 3.1))
		self.assertFalse(spiSeMe_surrogate_iaaft.spiSeMe_surrogate_iaaft(self.data_A, 'invalid'))
		self.assertFalse(spiSeMe_surrogate_iaaft.spiSeMe_surrogate_iaaft(self.data_A, 'distribution', 'x'))
		self.assertFalse(spiSeMe_surrogate_iaaft.spiSeMe_surrogate_iaaft(self.data_A, 'distribution', 1, 'x'))
	def test_surrogate_iaaft_run(self):
		print("IAAFT: Checking output formats...")
		self.M = 3
		self.data_out = spiSeMe_surrogate_iaaft.spiSeMe_surrogate_iaaft(self.data_A, 'distribution', self.M, False)
		self.assertTrue(isinstance(self.data_out,np.ndarray))
		self.assertEqual(self.data_out.ndim, 2)
		self.assertEqual(self.data_out.size, self.M * self.data_A.size)
		self.assertEqual(self.data_out[0].size, self.data_A.size)
		self.p, self.d = spiSeMe_distribution_test.spiSeMe_distribution_test(self.data_A, self.data_out[0])
		self.assertEqual(self.p, 1)

	def test_surrogate_sa_arg(self):
		print("SA: Intentionally providing bad arguments...")
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_Z, 0.1, 1.0))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 'x', 1.0))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 'x'))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 1.0, 0.1))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, -0.1, 1.0))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, -1.0))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 'x'))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 1.1))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, 'x'))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, -1.0))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, 0.1, 1.0, 'invalid'))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, 0.1, 1.0, 'max', 'x'))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, 0.1, 1.0, 'max', 100, 'x'))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, 0.1, 1.0, 'max', 100, 50, 'x'))
		self.assertFalse(spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, 0.1, 1.0, 'max', 100, 50, 1, 'x'))
	def test_surrogate_sa_run(self):
		print("SA: Checking output formats...")
		self.M = 3
		self.data_out = spiSeMe_surrogate_sa.spiSeMe_surrogate_sa(self.data_A, 0.1, 1.0, 0.9, 0.1, 1.0, 'max', 100, 50, self.M, False)
		self.assertTrue(isinstance(self.data_out,np.ndarray))
		self.assertEqual(self.data_out.ndim, 2)
		self.assertEqual(self.data_out.size, self.M * self.data_A.size)
		self.assertEqual(self.data_out[0].size, self.data_A.size)
		self.p, self.d = spiSeMe_distribution_test.spiSeMe_distribution_test(self.data_A, self.data_out[0])
		self.assertEqual(self.p, 1)

	def test_surrogate_dither_arg(self):
		print("Dither: Intentionally providing bad arguments...")
		self.assertFalse(spiSeMe_surrogate_dither.spiSeMe_surrogate_dither(self.data_Z))
		self.assertFalse(spiSeMe_surrogate_dither.spiSeMe_surrogate_dither(self.data_A, 3.1))
		self.assertFalse(spiSeMe_surrogate_dither.spiSeMe_surrogate_dither(self.data_A, 'invalid'))
		self.assertFalse(spiSeMe_surrogate_dither.spiSeMe_surrogate_dither(self.data_A, 'uniform', 'x'))
		self.assertFalse(spiSeMe_surrogate_dither.spiSeMe_surrogate_dither(self.data_A, 'uniform', 0.1, 'x'))
		self.assertFalse(spiSeMe_surrogate_dither.spiSeMe_surrogate_dither(self.data_A, 'uniform', 0.1, 1, 'x'))
	def test_surrogate_dither_run(self):
		print("Dither: Checking output formats...")
		self.M = 3
		self.data_out = spiSeMe_surrogate_dither.spiSeMe_surrogate_dither(self.data_A, 'uniform', 0.01, self.M, False)
		self.assertTrue(isinstance(self.data_out,np.ndarray))
		self.assertEqual(self.data_out.ndim, 2)
		self.assertEqual(self.data_out.size, self.M * self.data_A.size)
		self.assertEqual(self.data_out[0].size, self.data_A.size)

if __name__ == '__main__':
    unittest.main()
