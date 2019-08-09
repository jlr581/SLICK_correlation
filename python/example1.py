'''
Created on 24 July 2019

Example code for calling the cslick correlation method.

@author: Lenneke Jong
'''
from correlation import correlation

import numpy as np


if __name__ == '__main__':
 
    series1=np.loadtxt('test_series1.txt')
    series2=np.loadtxt('test_series2.txt')

    print 'Call the cslick correlation function as follows:'
    print 'cslick=correlation(x1,y1,x2,y2,hc)'
    print '''where (x1,y1) and (x2,y2) are two differently and unevenly sampled series to be correlated,
in this implementation these are numpy arrays. Optional argument hc (default=0.4) is a coefficient to
tune  how closely data from the two series must be to be included in the calculation.\n'''
    

    print "test_series1.txt"
    for row in series1[:]:
        print "{:1.6f} {:1.6f}".format(row[0],row[1]) 
    print "\ntest_series2.txt"
    for row in series2[:]:
        print "{:1.6f} {:1.6f}".format(row[0],row[1]) 
    
    cslick=correlation(series1[:,1],series2[:,1],series1[:,0],series2[:,0],hc=0.4)     
    
    print  '\nCorrelation: {}'.format(cslick)
