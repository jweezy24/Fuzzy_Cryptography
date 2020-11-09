import unittest
import sys
sys.path.append("../src")
from galois import *
import random

# The prime polynomial is a prime within a field that we choose.
prime_poly = 0b100011101
# The generator element is a element that when multiplied to itself will generate all elements within the field
generator = 0b10000011
class TestFields(unittest.TestCase):

    '''Tests field creation'''
    def test_field_creation(self):
       GF = GaloisField(8, prime_poly, generator)
       self.assertEqual(GF.max_num, 2**8)
       self.assertEqual(GF.prime_poly, prime_poly)
       self.assertEqual(GF.generator, generator)

    ''' Tests the setup_tables piece of initalization.
        The setup_tables method can't have any repeats besides 1 within the inverse log.
        This is supposed to be a 1-1 mapping.
        There are two 1s because its also supposed to be a cycle and the 1s denote the beginning and end.'''
    def test_setup_tables(self):
        GF = GaloisField(8, prime_poly, generator)
        
        for i in range(0, GF.max_num):
            test_ele = GF.gflog[i]
            for j in range(0, GF.max_num):
                if i != j:
                    self.assertFalse(test_ele == GF.gflog[j])

        for i in range(0, GF.max_num):
            test_ele = GF.gfilog[i]
            for j in range(0, GF.max_num):
                if i != j  and test_ele != 1:
                    self.assertFalse(test_ele == GF.gfilog[j])

    ''' Because we are working in a field, we must test addition, multiplication, and division.
        Subtraction in our Galois field is uneccisary because we are working in a base 2 field meaning that we are looking at integers as polynomials.
        Base 2 Galois fields do not have negative number representations.
        Addition in a base 2 Galois field is a XOR operation.'''
    def test_addition(self):
        GF = GaloisField(8, prime_poly, generator)
        
        for i in range(0, GF.max_num):
            for j in range(0, GF.max_num):
                result = GF.add(i,j)
                if i == j:
                    self.assertEqual(result, 0)
                else:
                    self.assertTrue(result < GF.max_num)
                    self.assertTrue(j^i == result)

        result = GF.add(1000, 2)
        self.assertEqual(result, 0)

    '''Multiplication needs to be closed (no number greater than 256).
       We also need to check if it is able to catch edge cases. 
       Being that there are 256 elements within out field we need to tests if they, when multiplied together give a nmumber within the set.
       '''
    def test_multiplication(self):
        GF = GaloisField(8, prime_poly, generator)
        
        for i in range(0, GF.max_num):
            for j in range(0, GF.max_num):
                result = GF.mult(i,j)
                self.assertTrue(result < GF.max_num)

        
        result = GF.mult(1000, 2)
        self.assertEqual(result, 0)

    ''' Our division algorithm is an inverse of the multiplication.
        '''
    def test_division(self):
        GF = GaloisField(8, prime_poly, generator)
        
        for i in range(0, GF.max_num):
            for j in range(0, GF.max_num):
                result = GF.div(i,j)
                self.assertTrue(result < GF.max_num)

        
        result = GF.mult(1000, 2)
        self.assertEqual(result, 0)

    ''' Tests the exponential function'''
    def test_exponential(self):
        GF = GaloisField(8, prime_poly, generator)
        res = 1
        test_arr = []

        for i in range(0, 255):
            res = GF.pow(generator, i)
            self.assertTrue(res not in test_arr)
            test_arr.append(res)

    '''This test will check to see if the polynomial objects can evaluate correctly '''
    def test_polynomial_evaluation(self):
        GF = GaloisField(8, prime_poly, generator)

        co2 = [212, 225, 63, 173, 11]
        
        p = polynomial(len(co2))
        p.set_coeffs(co2)

        for i in range(0,10):
            val = GF.eval_poly(p, i)
            if i == 0:
                self.assertEqual(val, 0)
            else:
                self.assertTrue(val < 256)

    ''' Here we test if the inverse exists for all elements in the set.
        0 Does not have a mutiplicative inverse so we do not test for it'''
    def test_inverse_value_finder(self):
        GF = GaloisField(8, prime_poly, generator)

        for i in range(1,256):
            val = GF.get_inverse(i)
            self.assertTrue(val != -1)
        






if __name__ == "__main__":
    unittest.main()