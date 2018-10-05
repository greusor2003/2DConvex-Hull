# -*- coding: utf-8 -*-
"""
Algorithm for finding the convex hull of a list of points.
The set of points is implemented as a (2,) ndarray (numpy type), each column
representing the (x,y) coordinates of a point.

Step 0: The functions scatter_points() and no_dupli() generate L, the list of 
        distinct points to find the convex hull for. CH = the list of the 
        corner points in the boundary of the convex hull.
        L = the algorithm's input
        CH = the algorithms output

Step 1: Among the points of the lowest y, find the point with the highest x and
        include it in CH. Remove it from L. 

Step 2: If L nonempty, find a point q in L so that the line from the last point
        in CH to q slopes the least through positive values from the horizontal. 
        Include q in CH. Remove q from L.

Step 3: As long as P is nonempty, find the list P of points on the left side of
        the vector from the beginning to the last point of CH. If P empty, stop
        and return CH. If P is nonempty, the next point in CH is searched for 
        in P as q was at Step 2, but instead of the horizontal, the direction
        of the vector from the second last to the last point in CH is
        used instead. Remove the point thus found from P. Repeat.
        When P is empty, return CH.
        
Created on Sat Dec  2 13:19:24 2017

@author: mikos
"""

import numpy as np
import matplotlib.pyplot as plt # for graphics output
import time # for performance tests only


def no_dupli(L):
    """ removes duplicates from the input list of (1,2) ndarrays.
    """
    N = [L.pop(0)]
    while L != []:
        k = 0
        flag = True
        while k < len(N) and L != [] and flag:
            if (N[k] == L[0]).all():
                L.pop(0)
                flag = False
            else:
                k = k + 1
                if k == len(N):
                    N.append(L.pop(0))
                    flag = False
    
    return N
    

def plot_points(L, color):
    """
    scatterplot of color 'color' for the points of L
    """
    
    X = list()
    Y = list()
    for p in L:
        X.append(p[0])
        Y.append(p[1])
    plt.scatter(X, Y, c=color)

def scatter_points(n):
    """
    returns a list of n random points whose (x,y) coordinates make up a 
    shape (n,2) ndarray;
    the n-random points must be distinct in order to assign L as input
    in the algorithm
    """
    P1 = np.random.randn(int(np.ceil(n/2)), 2) - 4
    P2 = 3 * np.random.rand(int(np.ceil(n/4)), 2) - np.array([10, 0])
    P3 = np.random.randn(int(np.ceil(n/4)), 2) + 3
    """
    P1=np.floor(P1)
    P2=np.floor(P2)
    P3=np.floor(P3)
    """
    L = list(np.concatenate((P1,P2,P3), axis=0))
    
    return L 
    #return  no_dupli(L)


def find_start(L):
    """
    returns the index of the rightmost lowest y-coordinate point  :
    """
    y_min = L[0][1] # initialization
    mY=[0] #holds the indices of the min y-coordinate points from min to max x-coord 
    
    for k in range(len(L)): #so designed as L is traversed once
        
        lx = L[k][0]
        ly = L[k][1]
        
        if ly < y_min:
            y_min = ly
            mY.clear()
            mY.append(k)
        else:
            if ly == y_min and mY != []: 
                if lx > L[mY[len(mY)-1]][0]:
                    mY.append(k)
                else:
                    if lx < L[mY[0]][0]:
                       mY.insert(0,k)
                        
    return (mY[-1])

"""
function returns the ndarray of unit vectors starting from the point
p in the directions of the points in the list L
"""
def normalize(p,L):
    """
    returns the list of unit vectors starting from the point
    p in the directions of the points in the list L
    """    
    Norm = []
    for k in L:
        x = k - p
        norm =  x.dot(x)**0.5
        if norm != 0:
            Norm.append(x/norm)
        else:
            print('warning: attempt to normalize 0')
        
    return Norm
        

def side_points(p, v, L):
    """
    function returning the points from the list L situated on the left side of 
    the line passing through p of direction v
    p = point on the line, v=directional vector of the line
    p and v are (2,) ndarrays
    """ 
    u = np.array([-v[1], v[0]]) # positive normal of v:
    N = list() # list of points on one side of the line p,v:
    for k in range(len(L)):
        if (L[k] - p).dot(u) >= 0:
                N.append(L[k])
    
    return N

def next_in_hull(p, v, L):
    """
    finds the next point of the convex hull;
    p = curent point; v= direction of support line;
    (p, v) defines the support line trough p of direction v;
    L = list of points on the left side of the support line; 
    where the new point of the convex hull is searched for;
    """ 
    N = normalize(p, L)
    if N != []:
        q = N[0]
        index = 0
        for k in range(1, len(N)):
            if (N[k] - q).dot(v) >= 0: # points on support line included
                 q = N[k]
                 index = k
        
        return index


def convex_hull(L):
    """
    the main function returning corner points of the convex hull of L
    """
    CH=list()
    if L != []:
        P = list(L)
        # find the starting point of the algorithm and add it to the convex hull:
        ind0 = find_start(P)
        CH.append(P.pop(ind0))
        # find the next point and add it to the convex hull list CH:
        if P != []:
              ind1 = next_in_hull(CH[0], np.array([1,0]), P)
              CH.append(P.pop(ind1))
              # use the hyperplane criterion as function side_points to complete CH:
              while P != []:
                  p = CH[-2]
                  q = CH[-1]
                  v = q - p 
                  P = side_points(CH[0], CH[-1] - CH[0], P)
                  ind = next_in_hull(q, v, P)
                  if P != []:
                      CH.append(P.pop(ind))
        return CH

"running an example:"
if __name__ == '__main__':
    tic = time.time()
    L = list(scatter_points(10000))
    toc = time.time()
    print('time to generate the initial list of points: ', toc - tic)
    
    plot_points(L, 'g')
    tic = time.time()
    CH = convex_hull(L)
    toc = time.time()
    print("time to generate the convex hull's boundary vertices: ", toc - tic)
    plot_points(CH, 'r')
    CH.append(CH[0])
    D = np.array(CH)
    plt.plot(D[:, 0], D[:, 1])

    plt.axis('equal')
    plt.show()
