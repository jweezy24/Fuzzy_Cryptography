# Example Prime poly and generator
# unsigned int prim_poly = (unsigned int) 0b100011101;
# int generator = (unsigned int)0b10000011;

class GaloisField():
    def __init__(self, size, prime_poly, generator):
        self.size = size
        self.max_num = 2**size
        self.prime_poly = prime_poly
        self.gflog = [0 for i in range(0, 2**size)]
        self.gfilog = [ 0 for i in range(0, 2**size)]

        self.setup_tables()
    
    def setup_tables(self):
        b = 1
        log = 0
        for i in range(0, 2**self.size):
            self.gflog[b] = log
            self.gfilog[log] = b
            b = b << 1
            if(b & self.max_num):
                b = (b ^ self.prime_poly)
            
if __name__ == "__main__":
    prime_poly = 0b100011101
    generator = 0b10000011
    test = GaloisField(8, prime_poly, generator)

