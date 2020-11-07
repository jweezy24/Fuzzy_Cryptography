import unittest
import sys
sys.path.append("../src")
from galois import *
import random

# The prime polynomial is a prime within a field that we choose.
prime_poly = 0b100011101
# The generator element is a element that when multiplied to itself will generate all elements within the field
generator = 0b10000011
#This is the default field
GF = GaloisField(8, prime_poly, generator)
class TestPolynomial(unittest.TestCase):
    ''' This test verifies the arithmetic object's creation. '''
    def test_creation(self):
        PA = polynomial_arithmetic(GF)
        self.assertTrue(PA.field == GF)

    ''' This test examines polynomial addition.
        Working in a Galois field requires that when we add two polynomials the result be in the same set.
        '''
    def test_addition(self):
        PA = polynomial_arithmetic(GF)
        
        # 1 + x + 3x^2
        co1 = [1,1,3]
        # 10 + 2x + 6x^2 + 200x^3
        co2 = [10, 2, 6, 200]

        cores = [11, 3, 5, 200]

        p1 = polynomial(3)
        p2 = polynomial(4)

        p1.set_coeffs(co1)
        p2.set_coeffs(co2)

        result = PA.add(p1, p2)

        self.assertEqual(result.coeffs, cores)

    def test_addition_complex(self):
        PA = polynomial_arithmetic(GF)
        
        # 177 + 24x + 76x^2 + 104x^3 + 147x^4
        co1 = [177,24,76,104,147]
        # 212 + 225x + 63x^2 + 173x^3 + 11x^4
        co2 = [212, 225, 63, 173, 11]

        cores = [101, 249, 115, 197, 152]

        p1 = polynomial(len(co1))
        p2 = polynomial(len(co2))

        p1.set_coeffs(co1)
        p2.set_coeffs(co2)

        result = PA.add(p1, p2)

        self.assertEqual(result.coeffs, cores)
    
    ''' These tests evaluate multiplication for a simple case and a complex case.
    '''

    def test_multiplcation(self):
        PA = polynomial_arithmetic(GF)
        
        # 1 + x + 3x^2
        co1 = [1,1,3]
        # 10 + 2x + 6x^2 + 200x^3
        co2 = [10, 2, 6, 200]

        cores = [10, 8, 26, 200, 194, 69]

        p1 = polynomial(3)
        p2 = polynomial(4)

        p1.set_coeffs(co1)
        p2.set_coeffs(co2)

        result = PA.mult(p1, p2)

        self.assertEqual(result.coeffs, cores)

    def test_multiplication_complex(self):
        PA = polynomial_arithmetic(GF)
        
        # 177 + 24x + 76x^2 + 104x^3 + 147x^4
        co1 = [177,24,76,104,147]
        # 212 + 225x + 63x^2 + 173x^3 + 11x^4
        co2 = [212, 225, 63, 173, 11]

        cores = [204,208,250,193,96,19,227,39,68]

        p1 = polynomial(len(co1))
        p2 = polynomial(len(co2))

        p1.set_coeffs(co1)
        p2.set_coeffs(co2)

        result = PA.mult(p1, p2)

        self.assertEqual(result.coeffs, cores)

    ''' These tests evaluate division for a simple case and a complex case.
        We also verify the remainders.
    '''

    def test_multiplcation(self):
        PA = polynomial_arithmetic(GF)
        
        # 1 + x + 3x^2
        co1 = [1,1,3]
        # 10 + 2x + 6x^2 + 200x^3
        co2 = [10, 2, 6, 200]

        cores = [0]
        coremainder = [1, 1, 3]

        p1 = polynomial(3)
        p2 = polynomial(4)

        p1.set_coeffs(co1)
        p2.set_coeffs(co2)

        result = PA.div(p1, p2)
        remainder = PA.div(p1, p2, 1)

        self.assertEqual(result.coeffs, cores)
        self.assertEqual(remainder.coeffs, coremainder)

    def test_multiplication_complex(self):
        PA = polynomial_arithmetic(GF)
        
        # 177 + 24x + 76x^2 + 104x^3 + 147x^4
        co1 = [177,24,76,104,147]
        # 212 + 225x + 63x^2 + 173x^3 + 11x^4
        co2 = [212,225,63,173,11]

        cores = [79]
        coremainder = [106, 60, 111, 5]
        
        p1 = polynomial(len(co1))
        p2 = polynomial(len(co2))

        p1.set_coeffs(co1)
        p2.set_coeffs(co2)

        result = PA.div(p1, p2)
        remainder = PA.div(p1, p2, 1)

        self.assertEqual(result.coeffs, cores)
        self.assertEqual(remainder.coeffs, coremainder)

        result = PA.div(p2, p1)
        remainder = PA.div(p2, p1, 1)
    

        cores = [147]
        coremainder = [115,200,150,229]
        
        self.assertEqual(result.coeffs, cores)
        self.assertEqual(remainder.coeffs, coremainder)




if __name__ == "__main__":
    unittest.main()