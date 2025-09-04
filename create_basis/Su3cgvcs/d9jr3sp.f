CB----------------------------------------------------------------------
C                             ***************
C                             *** D9JR3SP ***
C                             ***************
C ----------------------------------------------------------------------
C Update:        ... original code by C. Bahri
C          05/01 ... modified according to SU3CGVCS by C. Bahri
C ----------------------------------------------------------------------
C     Special 9j symbol for SO(3) -- triangle relations checked in
C     delta
C ----------------------------------------------------------------------
C References:
C  1. D.M. Brink and G.R. Satchler, Angular Momentum (Oxford)
C ----------------------------------------------------------------------
      function d9jr3sp( j1t, j2t, j3t, j4t, j5t, j6t, j7t, j8t, j9t )
C     ------------------------------------------------------------------
      implicit none
      double precision d9jr3sp
      integer j1t, j2t, j3t, j4t, j5t, j6t, j7t, j8t, j9t
C
C     Temporary variables
      integer i1, i2, i3, i4, i5, i6, i7, i8, i9
      integer k1, k2, k3, k4, l1, l2, ix, ixmin, ixmax, iy, iymin, iymax
      double precision dx, dc, dtop, dbot, dsum
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
      d9jr3sp = 0.d0
      if( j7t .eq. j1t+j4t ) then
        dc = delta( j8t, j9t, j7t )
        dc = dc - delta( j1t, j2t, j3t )
        dc = dc - delta( j2t, j5t, j8t )
        dc = dc - delta( j4t, j5t, j6t )
        dc = dc - delta( j3t, j6t, j9t )
        dc = (dc + dlogf(2*j1t)+dlogf(2*j4t)-dlogf(2*j7t+2))/2.d0
          i1 = - j1t + j2t + j3t
          i2 =   j2t + j5t - j8t
          i3 = - j4t + j5t + j6t
          i4 =   j3t + j6t - j9t
          i5 =   j7t + j8t + j9t + 2
          i6 =   j1t + j2t + j3t + 2
          i7 =   j2t + j5t + j8t + 2
          i8 =   j4t + j5t + j6t + 2
          i9 =   j3t + j6t + j9t + 2
          k1 = - j3t + j6t + j9t
          k2 =   j2t - j5t + j8t
          k3 = - j1t + j3t - j5t + j8t
          k4 = - j3t - j4t + j5t + j9t
          dc = dc + dlogf(i1)+dlogf(i2)+dlogf(i3)+dlogf(i4)+dlogf(i5)
     1            -(dlogf(i6)+dlogf(i7)+dlogf(i8)+dlogf(i9))
          ixmin = max0( 0, -k4 )
          ixmax = min0( 2*j3t, i4 )
          do ix = ixmin, ixmax, 2
            dx = dc + dlogf(2*j3t-ix)+dlogf(k1+ix) - (dlogf(ix)
     1           +dlogf(i4-ix))
            l1 = k3 - ix
            l2 = k4 + ix
            iymin = max0( 0, -l1 )
            iymax = min0( 2*j5t, i2, l2 )
            do iy = iymin, iymax, 2
              dsum = dexp(dx + dlogf(2*j5t-iy)+dlogf(k2+iy)
     1           - (dlogf(iy)+dlogf(i2-iy)+dlogf(l1+iy)+dlogf(l2-iy)))
              if( btest( ix+iy, 1 ) ) then
                d9jr3sp = d9jr3sp - dsum
              else
                d9jr3sp = d9jr3sp + dsum
              end if
            end do
          end do
      else
        d9jr3sp = 1.d+38
      end if
C
      return
C ---*end of d9jr3sp*---------------------------------------------------
      end
