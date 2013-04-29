#import numpy #You may have to install this library
import sys
import cmath as cm
import math as m
import threading
#from multiprocessing import Process
"""import scipy
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import rand
#import scipy
from scipy import linalg, matrix
import numbthy"""
#17196283094983 #

n = 87463 #36443380917683970574932716576131
np = 917641 #3976418276500331
nq = 18739663 #9164876123081801
B = []

class myThread (threading.Thread):
    def __init__(self, _n, _p):
        threading.Thread.__init__(self)
        self._n = _n
        self._p = _p
    def run(self):
        legendreManager(self._n,self._p)


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
def lengthOfB(x):
    return m.floor((cm.e**((1/cm.sqrt(2))*cm.sqrt(cm.log(x)*cm.log(cm.log(x))))).real)
        
def isBSmooth(x,xB=[0]):
    global B
    if(xB == [0]):
         xB = [0]*len(B)
    for i in range(len(B)):
        if(x % B[i] == 0):
            x = x/B[i]
            xB[i] = (xB[i] + 1) % 2
            if(x == 1):
                return xB
            if(isBSmooth(x,xB)):
                return xB
    xB[0] = -1
    return xB

def legendre(_n,_p):
    return _n**((_p-1)/2)%_p
    
def legendreManager(_n,_p):
    global B
    if(legendre(_n,_p) == 1):
        B.append(p)
    
def prime4(upto=100):
    primes=[2]
    num = 3
    for num in range(3,upto,2):
        isprime=True
        for factor in range(3,1+int(m.sqrt(num)),2):
            if not num % factor: isprime=False; break
        if isprime: primes.append(num)
    return primes
    '''while(len(primes)<upto):
        isprime=True
        for factor in range(3,1+int(m.sqrt(num)),2):
            if not num % factor: isprime=False; break
        if isprime: primes.append(num)
        num = num + 2'''
    

"""def null(A, eps=1e-15):
    u, s, vh = scipy.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)"""
    
def myGauss(m):
    #eliminate columns
    for col in range(len(m[0])):
        for row in range(col+1, len(m)):
            r = [(rowValue * (-(m[row][col] / m[col][col]))) for rowValue in m[col]]
            m[row] = [sum(pair) for pair in zip(m[row], r)]
    #now backsolve by substitution
    ans = []
    m.reverse() #makes it easier to backsolve
    for sol in range(len(m)):
            if sol == 0:
                ans.append(m[sol][-1] / m[sol][-2])
            else:
                inner = 0
                #substitute in all known coefficients
                for x in range(sol):
                    inner += (ans[x]*m[sol][-2-x])
                #the equation is now reduced to ax + b = c form
                #solve with (c - b) / a
                ans.append((m[sol][-1]-inner)/m[sol][-sol-2])
    ans.reverse()
    return ans
    
if __name__ == '__main__':#
    #fake_primes = [2,3,5,7,11,13,17,19,23,29,31,37]
    Blength = long(lengthOfB(n))
    primes = prime4(Blength)
    print "Length of B: ",Blength
    #print primes
    #prime_list = [i for i in xrange(1,Blength) if isprime(i)]
    threads = []
    for p in primes:
        #curr = thread.start_new_thread(legendreManager,(n,p))
        curr = myThread(n,p)
        threads.append(curr)
        curr.start()
    for t in threads:
        t.join()
    print "Exiting Main Thread"
    root_n = int(m.floor(m.sqrt(n)))
    #print root_n
    XY = []
    sparse = []
    global B
    #print B
    print "len(B):",len(B)
    M=2
    N=0
    while(len(XY) < len(B)+1):
        for x in xrange(root_n-M,root_n-N):
            #print x
            y = (x)**2 % n
            ysmooth = isBSmooth(y)
            if(ysmooth[0] <> -1):
                XY.append([x,y])
                sparse.append(ysmooth)
                #print x,y,ysmooth[0]
        for x in xrange(root_n+N+1,root_n+M):
            #print x
            y = (x)**2 % n
            ysmooth = isBSmooth(y)
            if(ysmooth[0] <> -1):
                XY.append([x,y])
                sparse.append(ysmooth)
                #print x,y,ysmooth[0]
        N = M
        M = M+1
    #sparse = numpy.matrix(scipy.transpose(sparse))
    #nulls = null(sparse)
    #nulls = numpy.hsplit(nulls,nulls.shape[1])
    #print nulls
    xt = 1
    yt = 1
    print "XY length:",len(XY)
    #print XY
    print "Sparse size:",len(sparse[0]),len(sparse)
    """print "sparse size:",sparse.shape
    print "nulls size:",len(nulls),",",len(nulls[0])
    for j in range(0,len(nulls)):
        for i in range(0,nulls[j].size):
            if(nulls[j][i] <> 0):
                #print i," ",XY[i][0]," ",XY[i][1]
                xt = xt * XY[i][0]
                yt = yt * XY[i][1]
        xt = xt % n
        yt = yt % n
        #print "gcd(",xt,"+",yt,",",n,") = ",numbthy.gcd((xt+yt),n)
        #print "gcd(-,",n,") = ",numbthy.gcd((xt-yt),n)
        gcd = numbthy.gcd(xt+yt,n)
        if(gcd <> 1 and gcd <> n):
            print "WE HAVE A WINNER!: ",gcd,"*",n/gcd
        #print ""
        xt = 1
        yt = 1
    #print nulls[1]"""
    
    
        
        
