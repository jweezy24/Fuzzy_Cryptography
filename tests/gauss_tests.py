import unittest
import sys
sys.path.append("../src")
from galois import *
from gauss import *
import random


# The prime polynomial is a prime within a field that we choose.
prime_poly = 0b100011101
# The generator element is a element that when multiplied to itself will generate all elements within the field
generator = 0b10000011
#This is the default field
GF = GaloisField(8, prime_poly, generator)
Gauss = GaussianObj(GF)

class TestReedSolomon(unittest.TestCase):

    def test_find_zero_constant(self):

        row1 = [1,2,3]
        row2 = [2,3,4]

        t1 =Gauss.find_constant_for_row_addition(row1, row2, 0)
        t2 =Gauss.find_constant_for_row_addition(row1, row2, 1)
        t3 =Gauss.find_constant_for_row_addition(row1, row2, 2)

        self.assertEqual(GF.mult(row1[0], t1) ^ row2[0], 0)
        self.assertEqual(GF.mult(row1[1], t2) ^ row2[1], 0)
        self.assertEqual(GF.mult(row1[2], t3) ^ row2[2], 0)

    def test_triangular_form_simple(self):
        matrix = [[13,32,73], [2,3,4]]
        Gauss.triangular_form(matrix)
        Gauss.print_matrix(matrix)

    def test_triangular_form_complex(self):
        matrix = [
            [210 ,214 ,130 ,125 , 20 ,  6],
            [153 ,215 ,108 ,119 ,125 , 92],
            [234 ,122 , 25 , 98 , 24 , 45],
            [186 ,179 ,191 ,218 , 58 , 85],
            [209 ,201 , 26 ,225 ,138 ,129]]
        Gauss.triangular_form(matrix)

    def test_substitution_simple(self):
        matrix = [[13,32,73], [2,3,4]]
        Gauss.triangular_form(matrix)
        #Gauss.print_matrix(matrix)
        mat =Gauss.substitute_values(matrix)


    def test_substitution_complex(self):
        matrix = [
            [210 ,214 ,130 ,125 , 20 ,  6],
            [153 ,215 ,108 ,119 ,125 , 92],
            [234 ,122 , 25 , 98 , 24 , 45],
            [186 ,179 ,191 ,218 , 58 , 85],
            [209 ,201 , 26 ,225 ,138 ,129]]
        Gauss.triangular_form(matrix)
        #Gauss.print_matrix(matrix)
        Gauss.substitute_values(matrix)

    def test_system_solving_simple(self):
        matrix = [[13,32,73], [2,3,4]]
        mat = Gauss.solve_system(matrix)
        print(mat)
        #Gauss.print_matrix(mat)



if __name__ == "__main__":
    unittest.main()