CB---------------------------------------------------------------------
C                             ************
C                             *** CGU3 ***
C                             ************
C ---------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 11/97)
C ---------------------------------------------------------------------
C Update:  11/97 ... original code by C. Bahri
C          10/99 ... cleanup
C ---------------------------------------------------------------------
C     Function to calculate the extremal CG of SU(3) > SU(2)xU(1)
C ---------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C  3. D.J. Rowe and J. Repka, JMP 38, 4363-4388 (1997)
C ---------------------------------------------------------------------
      subroutine cgu3( lm1,mu1, lm2,mu2, lm3,mu3, j1t,I1t, j2t,I2t,
     1                 j3t,I3t, nhw, ifu3e, dcgu3e, dcgu3, nmax )
      implicit none
      integer lm1,mu1, lm2,mu2, lm3,mu3         ! SU(3) irreps
      integer j1t,I1t, j2t,I2t, j3t,I3t         ! U(1)xSU(2) irreps
      integer nhw
      integer nmax
      integer          ifu3e( * )
      double precision dcgu3e( nmax, * )
      double precision dcgu3( * )
C
      integer kromax
      integer ial                               ! outer multiplicity
C
C Indices and temporary variables
      integer njI
      integer nbit8
C
      double precision factor, dtemp
      integer s1t,s2t,s3t
      integer lkey, itemp
      integer j1tp,I1tp, j2tp,I2tp
      integer j1td, j2td
      integer a, b, c, d, e, f
C
C Extrinsic functions
      double precision d9jr3sp
      double precision klw
      double precision dsqrt, dfloat, dexp
      integer ishft, iand, iabs                 ! intrinsic functions
      logical btest
C
C Common blocks from file # 4
      integer NFACT, NBINO, NTWOS, NMBIN
      double precision dlogf, dlnbin, twos
      parameter( NFACT = 2000, NBINO = 400, NTWOS = 128,
     1           NMBIN  = NBINO*(NBINO+3)/2 )
      common / BKDF / dlogf(0:NFACT)
      common / BBIN / dlnbin(0:NMBIN), twos(-NTWOS:NTWOS)
C
C Intrinsic functions
      integer iupac
      double precision dbino2, drr3sp
C
C     Unpacking
      data nbit8  / 255 /            ! zff
      iupac( j1t, j2t, j3t ) = iand( ishft( j1t, j2t ), j3t )
C
C     Binomials
      dbino2( a, b ) = dlnbin( a*(2*a+1) + 2*b )
C
C     Racah coefficients for special arguments (c = a + b)
      drr3sp( a, b, e, d, c, f ) = dexp( 0.5d0*( (dlogf(2*a)+dlogf(2*b)+
     1   dlogf(a+b+d+e+2)+dlogf(a+b-d+e)+dlogf(a+b+d-e)+dlogf(-a+e+f)+
     2   dlogf(-b+d+f)) -
     3  (dlogf(2*a+2*b+2)+dlogf(-a-b+d+e)+dlogf(a+e-f)+dlogf(a-e+f)+
     4   dlogf(a+e+f+2)+dlogf(b+d-f)+dlogf(b-d+f)+dlogf(b+d+f+2)) ) )
C
      s1t = mu1
      s2t = mu2
      s3t = mu3
C     if( nmax .lt. nhw(-10 ) ) call error(' Impossible')
      call u3mult( lm1,mu1, lm2,mu2, lm3,mu3, kromax, *900 )
C ----------------------------------------------------------------------
              do ial = 1, kromax
                dcgu3( ial ) = 0.d0
              end do ! ial
C ----------------------------------------------------------------------
      if( I1t .ge. iabs(s1t-j1t) .and. I1t .le. s1t+j1t .and.
     2    I2t .ge. iabs(s2t-j2t) .and. I2t .le. s2t+j2t .and.
     3    I3t .ge. iabs(s3t-j3t) .and. I3t .le. s3t+j3t .and.
     4    I3t .ge. iabs(I1t-I2t) .and. I3t .le. I1t+I2t ) then
C
          factor = klw(lm1,mu1,j1t,I1t) * klw(lm2,mu2,j2t,I2t)
     1        / klw(lm3,mu3,j3t,I3t) * dsqrt( dfloat( (j3t+1)*(s3t+1)
     2        * (I1t+1)*(I2t+1) ) ) * dexp( 0.5d0
     3        * (dlogf(2*j1t+2)+dlogf(2*j2t+2)) )
C
C         putting the correct phase
          if( btest( I3t-s3t-j3t, 1 ) ) factor = -factor
C
C         retrieving the extremal CG
C
          do njI = 1, nhw
            lkey = ifu3e( njI )
            j1tp = iupac( lkey, -24, nbit8 )
            I1tp = iupac( lkey, -16, nbit8 )
            j2tp = iupac( lkey,  -8, nbit8 )
            I2tp = iupac( lkey,   0, nbit8 )
            j1td = j1t - j1tp
            j2td = j2t - j2tp
C
            if( j2td .ge. 0 .and.
     1          I2tp .ge. iabs(I2t-j2td) .and. I2tp .le. I2t+j2td .and.
     2          j1td .ge. 0 .and.
     3          I1tp .ge. iabs(I1t-j1td) .and. I1tp .le. I1t+j1td .and.
     4          s3t .ge. iabs(I1tp-I2tp) .and. s3t .le. I1tp+I2tp ) then
	      itemp = j3t * (2*j3t+1) + 2*j1td
              dtemp = dsqrt( dfloat( (I1tp+1)*(I2tp+1) ) ) * dexp(-0.5d0
     1            * (dlogf(2*j1tp)+dlogf(2*j1td)+dlogf(2*j2tp)
     2               +dlogf(2*j2td) - dbino2( j3t, j1td )))
     3            * drr3sp(j1td,j1tp,I1t,s1t,j1t,I1tp)
     4            * drr3sp(j2td,j2tp,I2t,s2t,j2t,I2tp)
     5            * d9jr3sp( j1td, I1t, I1tp,j2td,I2t,I2tp,j3t,I3t,s3t )
     6            / ( klw(lm1,mu1,j1tp,I1tp) * klw(lm2,mu2,j2tp,I2tp) )
C
C ----------------------------------------------------------------------
              do ial = 1, kromax
                dcgu3( ial ) = dcgu3( ial ) + dtemp*dcgu3e( njI,ial )
              end do ! ial
C ----------------------------------------------------------------------
            end if ! possible recouplings
          end do ! njI
C ----------------------------------------------------------------------
              do ial = 1, kromax
                dcgu3( ial ) = dcgu3( ial ) * factor
              end do ! ial
C ----------------------------------------------------------------------
      end if ! I3t
C
900   continue
      return
C ---*end of cgu3*------------------------------------------------------
      end
