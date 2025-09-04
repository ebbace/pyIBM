CB----------------------------------------------------------------------
C                             ***************
C                             ***  D9JR3  ***
C                             ***************
C ----------------------------------------------------------------------
C Update:        ... original code by J.P. Draayer
C          05/01 ... modified according to SU3CGVCS by C. Bahri
C ----------------------------------------------------------------------
C     9j symbol for SO(3) 
C ----------------------------------------------------------------------
C References:
C  1. D.M. Brink and G.R. Satchler, Angular Momentum (Oxford)
C ----------------------------------------------------------------------
      function d9jr3( j1t, j2t, j3t, j4t, j5t, j6t, j7t, j8t, j9t )
C     ------------------------------------------------------------------
      implicit none
      double precision d9jr3
      integer j1t, j2t, j3t, j4t, j5t, j6t, j7t, j8t, j9t
C
      integer iabs, min0, max0
      double precision dfloat, drr3
C
C     Temporary variables
      integer it, itmin, itmax
C
      double precision dlogf
C
      common / BKDF / dlogf(0:2000)
C
      d9jr3 = 0.d0
      itmin = max0( iabs(j1t-j9t), iabs(j2t-j6t), iabs(j4t-j8t) )
      itmax = min0( j1t+j9t, j2t+j6t, j4t+j8t )
      if( itmin .gt. itmax ) return
      do it=itmin,itmax,2
        d9jr3 = d9jr3 + dfloat(it+1) * drr3(j1t,j9t,j4t,j8t,it,j7t) *
     1      drr3(j2t,j6t,j8t,j4t,it,j5t)*drr3(j1t,j9t,j2t,j6t,it,j3t)
      end do
C
      return
C ---*end of d9jr3*-----------------------------------------------------
      end
