# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 00:10:16 2017
Andrew's monotone chain 2D convex hull algorithm
Source: https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain

"""

import time  # to test execution time
import numpy as np
import matplotlib.pyplot as plt
 
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


def convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
    starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list. 
    return lower[:-1] + upper[:-1]

#driven program
    
tic = time.time()
Q = list(scatter_points(10000))
L = []
for q in Q:
    L.append((q[0], q[1]))
toc = time.time()
print('time to generate the initial list of points: ', toc - tic)
plot_points(L,'b')
tic = time.time()
CH=convex_hull(L)
toc = time.time()
print("time to generate the convex hull's boundary vertices: ", toc - tic)
plot_points(CH, 'r')
CH.append(CH[0])
D=np.array(CH)
plt.plot(D[:,0],D[:,1])

plt.axis('equal')
plt.show()




