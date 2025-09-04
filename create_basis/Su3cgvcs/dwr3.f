CB----------------------------------------------------------------------
C                             ****************
C                             ***   DWR3   ***
C                             ****************
C ----------------------------------------------------------------------
C Update:        ... original code by J.P. Draayer
C          05/01 ... modified according to SU3CGVCS by C. Bahri
C ----------------------------------------------------------------------
C     Wigner coefficients for SO(3) -- triangle relations checked in
C     delta.
C ----------------------------------------------------------------------
C References:
C  1. M.E. Rose, Elementary Theory of Angular Momentum (Wiley)
C ----------------------------------------------------------------------
      function dwr3( j1t, j2t, j3t, m1t, m2t, m3t )
C     ------------------------------------------------------------------
      implicit none
      double precision dwr3
      integer j1t, j2t, j3t, m1t, m2t, m3t
C
C     Temporary variables
      integer i1, i2, i3, i4, i5, it, itmin, itmax
      double precision dc, dtop, dbot, dsum
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
      dwr3 = 0.d0
        if( m1t+m2t-m3t .ne. 0 ) return
        dc = delta( j1t, j2t, j3t )
C
      i1 = j3t - j2t + m1t
      i2 = j3t - j1t - m2t
      i3 = j1t + j2t - j3t
      i4 = j1t - m1t
        if( btest( i4, 0 ) ) return
      i5 = j2t + m2t
        if( btest( i5, 0 ) ) return
C
      itmin = max0( 0, -i1, -i2 )
      itmax = min0( i3, i4, i5 )
        if( itmin .gt. itmax ) return
C
        dtop = ( dlog(dfloat(j3t+1)) + dc + dlogf(j1t+m1t) +
     1         dlogf(j1t-m1t) + dlogf(j2t+m2t) + dlogf(j2t-m2t) +
     2         dlogf(j3t+m3t) + dlogf(j3t-m3t) ) / dfloat(2)
C
      do it = itmin, itmax, 2
        dbot = dlogf(i3-it) + dlogf(i4-it) + dlogf(i5-it) +
     1         dlogf(it) + dlogf(i1+it) + dlogf(i2+it)
        dsum = dexp(dtop-dbot)
        if( btest( it, 1 ) ) then
          dwr3 = dwr3 - dsum
        else
          dwr3 = dwr3 + dsum
        end if
      end do
C
      return
C ---*end of dwr3*------------------------------------------------------
      end
