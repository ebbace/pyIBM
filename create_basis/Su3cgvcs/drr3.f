CB----------------------------------------------------------------------
C                             ***************
C                             ***  D6JR3  ***
C                             ***************
C ----------------------------------------------------------------------
C Update:        ... original code by J.P. Draayer
C          05/01 ... modified according to SU3CGVCS by C. Bahri
C ----------------------------------------------------------------------
C     6j symbols for SO(3) -- triangle relations checked in delta
C ----------------------------------------------------------------------
C References:
C  1. M. Rotenberg, R. Bivins, N. Metropolis, and J.K. Wooten,
C     The 3j- and 6j-symbols (MIT Press)
C ----------------------------------------------------------------------
      function drr3( j1t, j2t, j3t, l1t, l2t, l3t )
C     ------------------------------------------------------------------
      implicit none
      double precision drr3
      integer j1t, j2t, j3t, l1t, l2t, l3t
C
C     Temporary variables
      integer i1, i2, i3, i4, i5, i6, i7, it, itmin, itmax
      double precision dx, dc, dsum
C
      double precision dexp, dlogf
      logical          btest
C
      common / BKDF / dlogf(0:2000)
C
      double precision delta
      delta( j1t, j2t, j3t ) = dlogf(j1t+j2t-j3t) + dlogf(j2t+j3t-j1t)
     1        + dlogf(j3t+j1t-j2t) - dlogf(j1t+j2t+j3t+2)
C
      drr3 = 0.d0
        dc = delta( j1t, j2t, j3t )
        dc = dc + delta( l1t, l2t, j3t )
        dc = dc + delta( l1t, j2t, l3t )
        dc = dc + delta( j1t, l2t, l3t )
        dc = dc/2.d0
C
      i1 = j3t + l3t - j1t - l1t
      i2 = j3t + l3t - j2t - l2t
      i3 = j1t + j2t + l1t + l2t + 2
      i4 = j1t + j2t - j3t
      i5 = l1t + l2t - j3t
      i6 = j1t + l2t - l3t
      i7 = l1t + j2t - l3t
C
      itmin = max0(  0,-i1,-i2 )
      itmax = min0( i3, i4, i5, i6, i7 )
        if( itmin .gt. itmax ) return
C
      do it = itmin, itmax, 2
        dsum = dexp( dc + dlogf(i3-it) - ( dlogf(i4-it) + dlogf(i5-it) +
     1               dlogf(i6-it) + dlogf(i7-it) + dlogf(it) +
     2               dlogf(i1+it) + dlogf(i2+it) ) )
        if( btest( it, 1 ) )then
          drr3 = drr3 - dsum
        else
          drr3 = drr3 + dsum
        end if
      end do
C
      return
C ---*end of drr3*------------------------------------------------------
      end
