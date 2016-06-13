from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport round

ctypedef np.float64_t dtype_t
ctypedef np.int_t int_t
ctypedef np.uint8_t bool_t

cdef inline dtype_t dydtMS(dtype_t y, dtype_t t, dtype_t p0, dtype_t p1) nogil:
    return -p0*y + p1

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def euler_integrate(np.ndarray[dtype_t,ndim=2] lengthAdj, np.ndarray[dtype_t,ndim=1] p,
        np.ndarray[dtype_t,ndim=2] y0, dtype_t yth, dtype_t delta, dtype_t eps,
        int_t T,int_t M, bool_t fo, bool_t sos):
    '''
    Cythonized integrator.
    '''
    # some data and loop variables
    cdef dtype_t dt = (1.0*T)/M
    cdef dtype_t ampDev
    cdef int_t delay_ij
    cdef Py_ssize_t nNodes = y0.shape[0]
    cdef Py_ssize_t ycs = 1

    # arrays
    cdef np.ndarray[dtype_t,ndim=2] y = np.zeros((nNodes,M+1),dtype=np.float64)
    cdef np.ndarray[int_t,ndim=2] s = np.zeros((nNodes,M+1),dtype=np.int)
    cdef np.ndarray[dtype_t,ndim=2] pulses = np.zeros((nNodes,M+1),dtype=np.float64)
    cdef np.ndarray[int_t,ndim=1] nodesToReset = np.zeros((nNodes,),dtype=np.int)

    cdef Py_ssize_t i = 0
    cdef Py_ssize_t j,k,n,nn
    cdef int_t quitStep = M

    with nogil:
      for i in xrange(nNodes):
          y[i,0] = y0[i,0]

      # start stepping (M+1 ensures we store 0,dt,2*dt,...,M*dt=T)
      for i in xrange(1,M+1):

          # --- Check threshold ---
          for j in xrange(nNodes):
              if y[j,i-1] > yth:
                nodesToReset[j] = 1

          # --- Fire to neighbors ---
          for n in xrange(nNodes):
              if nodesToReset[n] > 0:
                  y[n,i-1] = 0
                  s[n,i-1] = 1
                  for nn in xrange(nNodes):
                      if lengthAdj[n,nn] > 0:
                          # here's where we use the distances
                          delay_ij = int(delta*round(lengthAdj[n,nn]/dt))
                          # if a pulse is to be added after T, ignore it
                          #   (it will never fire)
                          if i-1+delay_ij < M+1:
                              pulses[nn,i-1+delay_ij] += eps


          # --- Check for synchronization ---
          ampDev = 0.0
          for k in xrange(nNodes):
              ampDev = ampDev + (y[0,i-1] - y[k,i-1])*(y[0,i-1] - y[k,i-1])

          if ampDev < 1.0e-10 and sos > 0:
              quitStep = i-1
              break

          # --- Resolve pulses and clear nodesToReset ---
          for j in xrange(nNodes):
              y[j,i-1] += pulses[j,i-1]
              y[j,i] = y[j,i-1] + dt*dydtMS(y[j,i-1],(i-1)*dt,p[0],p[1])
              nodesToReset[j] = 0

    # out of the for loop
    if fo < 1:
        return np.reshape(y[:,quitStep],(nNodes,1)),np.reshape(y[:,quitStep],(nNodes,1))
    return y[:,:(quitStep+1)],s[:,:(quitStep+1)]
