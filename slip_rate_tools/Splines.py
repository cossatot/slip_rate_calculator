# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 10:26:52 2012

@author: itchy
"""

import numpy as nmp



def spline1d(x_out, x_data, y_data, x_slope, y_slope, *args):
    """
    SPLINE1D	1-D interpolation using Green's function for a spline in tension
    SPLINE1D will find a spline-based curve using continuous curvature splines
    in tension (if set).  The algorithm uses the Green's function for the spline.
    You can supply data constrains, slope constrains, or a mix of both.
    Solution can be evaluated at arbitrary locations
	
    Use one of the following 3 call formats:
        y = spline1d (x_out, x_data, y_data, x_slope, y_slope)
        y = spline1d (x_out, x_data, y_data, x_slope, y_slope, t)
        y = spline1d (x_out, x_data, y_data, x_slope, y_slope, t, cutoff)
        
    The input parameters are:
            
     x_out -  Desired output x positions
     x_data  -   coordinates of points with data constraints
     y_data  -   data constraints at the above points
     x_slope  -   coordinates of points with slope constraints
     y_slope  -   slope constraints at the above points
     t  - tension to use, 0 <= t <= 1
	  if t is a vector of length 2 the second value is taken as the lengthscale
     cutoff  - if set, eigenvalues whose ratio to the maximum eigenvalue are smaller
     than cutoff are zeroed out before the curve is evaluated.
	
     One of (x_data, y_data) and (x_slope, y_slope) can be ([], [])
     t, if not set, defaults to 0 (cubic spline).  t = 1 gives linear interpolation

      The loutput values are:
          y  - the interpolation
          l  - optionally, the eigenvalues of the linear system
      
      See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A 
      Green function Approach, Math. Geol., 30, 77-93.
     
    Ported to Python by Richard Styron, June 2012
    """
        
    
    # Pick a reasonable(?) lengthscale 
    length_scale = (nmp.amax(x_out) - nmp.amin(x_out)) / 50.
    
    if len(args) == 0:  # no tension selected, set default
        t = 0.
    elif len(args) == 1:    # cutoff not set
        t = args
        cutoff = 0.
    elif len(args) == 2:    
        t = args[0]
        cutoff = args[1]

    t = nmp.array([t])
    cutoff = 0.


        
    nt = len(t)
    if nt == 2:     # user gave both tension and lengthscale
        length_scale = t[1]
        t = t[0]
    
    # TODO:  Add error/exception for values of t outside of [0,1]
    
    
    # Misc initializations
    if t < 1:
        p = nmp.sqrt(t / (1 - t))
        p = p / length_scale
        
    n0 = 0
    n1 = 1
        
        
        # First we must enforce the use of column vectors for the data constraints
        # FIX THIS: some_vector.shape unpacks to (m,) instead of (m,1)
        
    x_out = nmp.matrix(x_out)
    [m,n] = x_out.shape
    if m < n:
        #x_out = x_out.reshape(x_out.shape[0], -1)
        x_out = x_out.T
        
    x_data = nmp.matrix(x_data)            
    [m,n] = x_data.shape
    if m < n:
        #x_data = x_data.reshape(x_data.shape[0], -1)
        x_data = x_data.T      
        
        y_data = nmp.matrix(y_data)
    [m,n] = y_data.shape
    if m < n:
        #y_data = y_data.reshape(y_data.shape[0], -1)        
        y_data = y_data.T    
    
    [n0, m0] = y_data.shape

    x_slope = nmp.mat(x_slope)
    [m,n] = x_slope.shape
    if m < n:
        #x_slope = x_slope.reshape(x_slope.shape[0], -1)
        x_slope = x_slope.T
    
    y_slope = nmp.mat(y_slope)
    [m,n] = y_slope.shape
    if m < n:
    #y_slope = y_slope.reshape(y_slope.shape[0], -1)
        y_slope = y_slope.T
    
    #n1 = len(x_slope)
    n1 = 0
    
    # Assembly final xp, yp vectors (possibly combination of data and slopes)
    
    #xp = nmp.array([[x_data] , [x_slope]])
    #yp = nmp.array([[y_data] , [y_slope]])
    
    # TODO: fix slope constrain problems by putting an 'if' statement here:
    #   if slopes exist, add to the vectors
    
    xp = nmp.matrix(x_data)
    yp = nmp.matrix(y_data)
    
    # Now build the square n x n linear system that must be solved for the alpha's
    
    n = n0 + n1
    
    A = nmp.zeros((n, n))
    
    for i in nmp.arange(0,n0):  # First add equations for data constraints
    
        r = xp[i] - xp
        ar = nmp.abs(r)
        
        if t == 0:
            B = (ar ** 3)            
            A[i,:] = B.T
        
        elif t == 1:
            B = (ar)
            A[i,:] = B.T
            
        else: 
            B = nmp.exp(nmp.multiply(-p, ar)) + nmp.multiply(p, ar)
            A[i,:] = B.T
            
    if n1 > 0:            
        for i in nmp.arange(0,n1):  # Then add equations for slope constraints
            
            j = i + n0
            r = xp[j] - xp
            ar = nmp.abs(r)
            
            if t == 0:
                B = 3.0 * (r * ar)
                A[j,:] = B.T
            
            elif t == 1:
                B = nmp.sign(r)
                A[j,:] = B.T
            
    # Done building square linear system, now solve it
        
    # TODO: fix for cutoff > 0.0 -- deal with nargout for SVD        
            
        if cutoff > 0.0:     # solve using SVD
        #    U, S, V = svd(A)
        #    s = nmp.diag(S)
        #    if 
            pass
        
    else:
        alpha = nmp.linalg.solve(A, yp)
        
        
        # Now evaluate final solution at output locations
        
    y = nmp.zeros((len(x_out),m0))
        
    for i in nmp.arange(0,n):
        r = xp[i] - x_out
        ar = nmp.abs(r)
            
        if t == 0:
            y = y + (ar ** 3) * alpha[i,:]
                
        elif t == 1:
            y = y + ar * alpha[i,:]
                
        else:
            B = nmp.exp(nmp.multiply(-p, ar)) + nmp.multiply(p, ar)
            y = y + nmp.multiply(B, alpha[i,:])

    y = nmp.array(y.T)

    y = y[0,:]
    
    return y

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
