import unittest
import sys
sys.path.append("../src")
from galois import *
import random

class TestPolynomial(unittest.TestCase):
    
    '''Tests for correct initalization of a polynomial without optional arguments'''
    def test_polynomial_creation_basic(self):
        num = random.randint(0,1000)
        p = polynomial(num)
        self.assertEqual(p.size, num)
        self.assertEqual(len(p.coeffs), 0)
    
    ''' Tests the optional poly arugment of the creation of a poly'''
    def test_polynomial_creation_with_field(self):
        num = random.randint(0,1000)
        p = polynomial(num)
        self.assertEqual(p.size, num)
        self.assertEqual(len(p.coeffs), 0)

        num2 = random.randint(0,1000)
        p2 = polynomial(num2, p)
        self.assertEqual(p2.size, num)
        self.assertEqual(len(p2.coeffs), 0)
    
    '''Tests the set size method for polynomial'''
    def test_set_size_method(self):
        size = random.randint(0,100000)
        p = polynomial(2)
        p.set_size(size)
        self.assertEqual(p.size, size)

    
    '''Tests if the coeffs are set correctly '''
    def test_set_coeffs(self):
        #1 + 2x + 3x^2
        coeffs = [1,2,3]
        p = polynomial(3)
        p.set_coeffs(coeffs)
        self.assertEqual([1,2,3], p.coeffs)

        coeffs_bad = [1,2,3,4]
        p.set_coeffs(coeffs)
        self.assertEqual([1,2,3], p.coeffs)

    '''Test the resize method for a polynomial object'''
    def test_resize_method(self):

        coeffs = [1,2,3,0,0,0]
        p = polynomial(6)
        p.set_coeffs(coeffs)
        self.assertEqual([1,2,3,0,0,0], p.coeffs)
        self.assertEqual(p.size, 6)

        p.resize()
        self.assertEqual([1,2,3], p.coeffs)
        self.assertEqual(p.size, 3)
    
    '''Tests for a correct evaluation of the zeros
        The array the method returns is a list of powers of the generator where the polynomial is zero.'''
    def test_find_all_zeros(self):
        # The prime polynomial is a prime within a field that we choose.
        prime_poly = 0b100011101
        # The generator element is a element that when multiplied to itself will generate all elements within the field
        generator = 0b10000011
        #This is the default field
        GF = GaloisField(8, prime_poly, generator)
        
        #1 + 2x + 3x^2
        coeffs = [1,1]
        p = polynomial(3)
        p.set_coeffs(coeffs)
        zeros = p.find_zeros(GF)
        zeros_check =[255]
        self.assertEqual(zeros, zeros_check)



if __name__ == "__main__":
    unittest.main()