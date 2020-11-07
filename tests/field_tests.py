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




if __name__ == "__main__":
    unittest.main()