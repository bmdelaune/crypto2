import numpy #You may have to install this library
import sys


def isprime(number):  
    if number<=1:  
        return 0  
    check=2  
    maxneeded=number  
    while check<maxneeded+1:  
        maxneeded=number/check  
        if number%check==0:  
            return 0  
        check+=1  
        return 1 

if __name__ == '__main__':
    sz = long(10**30)
    prime_list = [i for i in xrange(1,sz) if isprime(i)]
