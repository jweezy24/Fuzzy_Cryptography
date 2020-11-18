import unittest
import sys
sys.path.append("../src")
from galois import *
import random

# The prime polynomial is a prime within a field that we choose.
prime_poly = 0b10011
# The generator element is a element that when multiplied to itself will generate all elements within the field
generator = 0b0010
#This is the default field
GF = GaloisField(4, prime_poly, generator)
class TestReedSolomon(unittest.TestCase):
    ''' Tests the creation of a ReedSolomonObj'''
    def test_init(self):
        n = 10
        k = 1
        RS = ReedSolomonObj(GF, n, k)

        self.assertEqual(RS.n, n)
        self.assertEqual(RS.k, k)
        self.assertTrue(RS.field)
        self.assertTrue(RS.PA)

    ''' This tests if the generator polynomial is created correctly. '''
    def test_generator_polynomial(self):
        n = 15
        k = 9
        RS = ReedSolomonObj(GF, n, k)

        cores = [GF.pow(generator, 6), GF.pow(generator, 9),GF.pow(generator, 6), GF.pow(generator, 4),
        GF.pow(generator, 14), GF.pow(generator, 10), 1]

        g = RS.get_generator_poly()
        
        self.assertEqual(cores, g.coeffs)

    '''Here we test the syndrome calculation from a scrambled message 
        The polynomials that we tests come from https://ntrs.nasa.gov/citations/19900019023'''
    def test_syndrome_calculation(self):
        n = 15
        k = 9
        RS = ReedSolomonObj(GF, n, k)


        g = RS.get_generator_poly() 
        C = polynomial(9)
        C_Coeffs = [GF.pow(generator,12),GF.pow(generator,8), GF.pow(generator,3), 
        GF.pow(generator,4), GF.pow(generator,10), GF.pow(generator,8), 0, GF.pow(generator,11), 1]
        C.set_coeffs(C_Coeffs)

        correct_syndromes = [1,1,GF.pow(generator,5),1,0,GF.pow(generator,10)]
        poly,syndromes = RS.calculate_syndrome(C, g)
        self.assertEqual(correct_syndromes, poly.coeffs)
    
    '''Here we test the error locator polynomial creation to verify it's creation 
        The polynomials that we tests come from https://ntrs.nasa.gov/citations/19900019023'''
    def test_error_locator_poly(self):
        n = 15
        k = 9
        RS = ReedSolomonObj(GF, n, k)
        
        g = RS.get_generator_poly() 
        C = polynomial(9)
        C_Coeffs = [GF.pow(generator,12),GF.pow(generator,8), GF.pow(generator,3), 
        GF.pow(generator,4), GF.pow(generator,10), GF.pow(generator,8), 0, GF.pow(generator,11), 1]
        C.set_coeffs(C_Coeffs)

        correct_syndromes = [1,1,GF.pow(generator,5),1,0,GF.pow(generator,10)]
        poly,syndromes = RS.calculate_syndrome(C, g)
        self.assertEqual(correct_syndromes, poly.coeffs)

        sig_r = RS.get_sigma_r(poly)

        terms = [1,1,GF.pow(generator,5),1,0,GF.pow(generator,10)]
        terms.reverse()

        self.assertEqual(sig_r.coeffs, terms)
    
    '''Here we test the berlecamp alg from a scrambled message 
        The polynomials that we tests come from https://ntrs.nasa.gov/citations/19900019023'''
    def test_berlecamp_massy(self):
        n = 15
        k = 9
        RS = ReedSolomonObj(GF, n, k)
        
        g = RS.get_generator_poly() 
        C = polynomial(9)
        C_Coeffs = [GF.pow(generator,12),GF.pow(generator,8), GF.pow(generator,3), 
        GF.pow(generator,4), GF.pow(generator,10), GF.pow(generator,8), 0, GF.pow(generator,11), 1]
        C.set_coeffs(C_Coeffs)

        correct_syndromes = [1,1,GF.pow(generator,5),1,0,GF.pow(generator,10)]
        poly,syndromes = RS.calculate_syndrome(C, g)
        self.assertEqual(correct_syndromes, poly.coeffs)

        correct_sigma = [1, 1, GF.pow(generator, 10)]
        sig = RS.berlecamp_alg(poly)
        self.assertEqual(sig.coeffs, correct_sigma)
        
    '''Here we test if our code is able to properly decode the example in the Nasa tech report.
        The polynomials that we tests come from https://ntrs.nasa.gov/citations/19900019023'''
    def test_error_correction(self):
        n = 15
        k = 9
        RS = ReedSolomonObj(GF, n, k)
        
        
        g = RS.get_generator_poly() 
        C = polynomial(9)
        C_Coeffs = [GF.pow(generator,12),GF.pow(generator,8), GF.pow(generator,3), 
        GF.pow(generator,4), GF.pow(generator,10), GF.pow(generator,8), 0, GF.pow(generator,11), 1]
        C.set_coeffs(C_Coeffs)

        correct_syndromes = [1,1,GF.pow(generator,5),1,0,GF.pow(generator,10)]
        poly,syndromes = RS.calculate_syndrome(C, g)
        self.assertEqual(correct_syndromes, poly.coeffs)

        correct_sigma = [1, 1, GF.pow(generator, 10)]
        sig = RS.berlecamp_alg(poly)
        self.assertEqual(sig.coeffs, correct_sigma)
        
        correct_sigma_r = [1, 1, GF.pow(generator, 10)]
        correct_sigma_r.reverse()
        s_r = RS.get_sigma_r(sig)

        self.assertEqual(s_r.coeffs, correct_sigma_r)

        zeros = s_r.find_zeros(GF)
        zeros_check =[2,8]
        self.assertEqual(zeros, zeros_check)

        found_errors = RS.find_error_values(poly, zeros)
        actual_errors = [1,1]
        self.assertEqual(actual_errors, found_errors)


        corrected_C =  [GF.pow(generator,12),GF.pow(generator,8), GF.pow(generator,3)^1, 
        GF.pow(generator,4), GF.pow(generator,10), GF.pow(generator,8), 0, GF.pow(generator,11), 1^1]
        C2 = polynomial(9)
        C2.set_coeffs(corrected_C)
        C2.resize()
        RS.correct_found_errors(C,zeros,found_errors)
        self.assertEqual(C.coeffs,C2.coeffs)

    '''Here we test the algorithm against random errors. UNDER CONSTRUCTION
        The polynomials that we tests come from https://ntrs.nasa.gov/citations/19900019023'''
    def test_error_correction_2(self):
        n = 15
        k = 9
        RS = ReedSolomonObj(GF, n, k)
        
        
        g = RS.get_generator_poly() 
        C = polynomial(9)
        C_Coeffs = [GF.pow(generator,12),GF.pow(generator,8), GF.pow(generator,3), 
        GF.pow(generator,4), GF.pow(generator,10), GF.pow(generator,8), 0, GF.pow(generator,11), 1]
        C.set_coeffs(C_Coeffs)

        correct_syndromes = [1,1,GF.pow(generator,5),1,0,GF.pow(generator,10)]
        poly,syndromes = RS.calculate_syndrome(C, g)
        self.assertEqual(correct_syndromes, poly.coeffs)

        correct_sigma = [1, 1, GF.pow(generator, 10)]
        sig = RS.berlecamp_alg(poly)
        self.assertEqual(sig.coeffs, correct_sigma)
        
        correct_sigma_r = [1, 1, GF.pow(generator, 10)]
        correct_sigma_r.reverse()
        s_r = RS.get_sigma_r(sig)

        self.assertEqual(s_r.coeffs, correct_sigma_r)

        zeros = s_r.find_zeros(GF)
        zeros_check =[2,8]
        self.assertEqual(zeros, zeros_check)

        found_errors = RS.find_error_values(poly, zeros)
        actual_errors = [1,1]
        self.assertEqual(actual_errors, found_errors)


        corrected_C =  [GF.pow(generator,12),GF.pow(generator,8), GF.pow(generator,3)^1, 
        GF.pow(generator,4), GF.pow(generator,10), GF.pow(generator,8), 0, GF.pow(generator,11), 1^1]
        C2 = polynomial(9)
        C2.set_coeffs(corrected_C)
        C2.resize()
        RS.correct_found_errors(C,zeros,found_errors)
        self.assertEqual(C.coeffs,C2.coeffs)

        


        

if __name__ == "__main__":
    unittest.main()