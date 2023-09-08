#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
# 
import numpy as np

def rotate_matrix(matrix, n=1):
    
    transposed = []
    
    for i in np.arange(n):
        for row in zip(*matrix):
            transposed.append(list(row))
        rotated = []
        for row in transposed:
            rotated.append(row[::-1])
        
    return rotated

# Function to rotate the matrix 90 degree clockwise
def rotate90Clockwise(arr, n) :


    # printing the matrix on the basis of
    # observations made on indices.
    for j in range(n) :
        for i in range(n - 1, -1, -1) :
            print(arr[i][j], end = " ")
            
    return arr