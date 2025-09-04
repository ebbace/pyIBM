CB---------------------------------------------------------------------
C                             **************
C                             *** CGU3HW ***
C                             **************
C ---------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 11/98)
C ---------------------------------------------------------------------
C Update:  11/98 ... original code by C. Bahri
C ---------------------------------------------------------------------
C     Function to calculate CG of SU(3) > U(1) x SU(2) with the highest
C     weight in the third entry.
C ---------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C ---------------------------------------------------------------------
      subroutine cgu3hw( lm1,mu1, lm2,mu2, lm3,mu3, nhw, ! krothy,
     1                   ifu3e, dcgu3e, nmax )
C
      implicit none
      integer lm1,mu1, lm2,mu2, lm3,mu3         ! SU(3) irreps
      integer nmax                              ! physical dimension
      integer nhw                               ! # of 3-hw CGs
      integer krothy                            ! rho-theory
      integer          ifu3e( * )               ! F-labels
      double precision dcgu3e( nmax, * )        ! F-functions/CGs of hw
C
      integer kromax                            ! outer multiplicity
      integer j1t,I1t, j2t,I2t, s1t, s2t, s3t   ! U(1)xSU(2) irreps
      integer jt, st, stmin, stmax              ! coupled index irreps
C
C Indices and temporary variables
      integer njI, njIF                         ! allowed and forbidden
      integer ncstr
      integer kro
C     integer jtot
C     integer nbit8
      integer NASI
      double precision a, b, dtemp
C
      parameter ( NASI = 128 )
C     integer          ifu3eF( NASI )
      double precision dmtx( NASI, NASI ), asi( NASI, NASI )
      double precision dtempa( NASI )
C
C Extrinsic functions
      double precision dfloat, dsqrt            ! intrinsic functions
      double precision klw, dexp, d9jr3sp       !, djhr3, d9jr3
      integer iabs, ishft, ior !, iand          ! intrinsic functions
      logical btest
C
C Common blocks from file # 4
      double precision dlogf
      common / BKDF / dlogf(0:2000)
C
C Intrinsic functions
      integer ipack !, iupac
C
C     Packing and unpacking
C     data nbit8  / 255 /            ! zff
      ipack( j1t, j2t, st, jt ) = ior(jt, ishft( ior(st, ishft(
     1        ior(j2t, ishft( j1t, 8 )),8 )),8 ))
C     iupac( j1t, j2t, j3t ) = iand( ishft( j1t, j2t ), j3t )
C
C     G(lm1,om1;lm2,lm2||om al lm3,lm3)bt = (-)^n1/K1 K2 sqrt(n!/n1!n2!)
C                * U9lm(lm1 lm2 om *; n1 n2 n -; om1 om2 om3 *; - - -)
C
C     Checking the allowed couplings ...
C
      call u3mult( lm1,mu1, lm2,mu2, lm3,mu3, kromax, *900 )
C
C     Initializations
C
      s1t = mu1
      s2t = mu2
      s3t = mu3
      j1t = ( 2*lm1+mu1 + 2*lm2+mu2 - (2*lm3+mu3) ) / 3
      j2t = 0
      jt  = j1t
      stmin = max0( iabs(s1t-s2t), iabs(s3t-jt) )
      stmax = min0( s1t+s2t, s3t+jt )
C
C     Non-orthonormalized CG
C
      njI  = 0
      njIF = 0
      do while ( j1t .ge. 0 .and. j2t .le. lm2+mu2)
C
        do I1t  = iabs(s1t-j1t), s1t+j1t, 2
        do I2t  = iabs(s2t-j2t), s2t+j2t, 2
        if( s3t .ge. iabs(I1t-I2t) .and. s3t .le. I1t+I2t ) then
C
C         Allowed states
          if( I1t .le. s1t-j1t+2*lm1 .and. I2t .le. s2t-j2t+2*lm2 ) then
C
            njI = njI + 1
            ifu3e( njI ) = ipack( j1t, I1t, J2t, I2t )
C
            kro = 0
            do st = stmax, stmin, -2
C!!!!         jtot = s1t+I1t + s2t+I2t + st+s3t
              kro = kro + 1
C!!!!         if( btest(j1t,0) .eor. btest(jtot,1) ) then
C!!!!         if( btest(j1t+jtot/2,0) ) then
              if( btest(j1t,0) ) then
                dcgu3e( njI, kro ) = - dexp( 0.5d0 * (dlogf(2*jt)
     1            - dlogf(2*j1t) - dlogf(2*j2t)) )
     2            * dsqrt(dfloat((st+1)*(jt+1)*(I1t+1)*(I2t+1)))
     3            * d9jr3sp(j1t,I1t,s1t, j2t,I2t,s2t, jt,s3t,st)
     4            / (klw(lm1,mu1,j1t,I1t)*klw(lm2,mu2,j2t,I2t))
              else
                dcgu3e( njI, kro ) =   dexp( 0.5d0 * (dlogf(2*jt)
     1            - dlogf(2*j1t) - dlogf(2*j2t)) )
     2            * dsqrt(dfloat((st+1)*(jt+1)*(I1t+1)*(I2t+1)))
     3            * d9jr3sp(j1t,I1t,s1t, j2t,I2t,s2t, jt,s3t,st)
     4            / (klw(lm1,mu1,j1t,I1t)*klw(lm2,mu2,j2t,I2t))
              end if ! phase
            end do ! st
C
C         Forbidden states -> constraint equations
          else
C
            njIF = njIF + 1
            if( njIF .gt. NASI ) call error(' CGU3HW: Increase NASI.')
            kro  = 0
C           ifu3eF( njIF ) = ipack( j1t, I1t, J2t, I2t )
            do st = stmax, stmin, -2
              kro = kro + 1
              dmtx( njIF, kro ) = 
     1              dsqrt(dfloat((st+1)*(jt+1)*(I1t+1)*(I2t+1)))
     2            * d9jr3sp(j1t,I1t,s1t, j2t,I2t,s2t, jt,s3t,st)
            end do ! st
C
          end if ! I1t and I2t
C
        end if ! s3t
        end do ! I2t
        end do ! I1t
C
        j1t = j1t - 1
        j2t = j2t + 1
C
      end do ! j1t & j2t
      nhw    = njI
      ncstr  = njIF
      krothy = kro
CW    write(6,120) nhw, krothy, kromax, ncstr
C     if( krothy .ne. kromax ) then
C       write(6,120) nhw, krothy, kromax, ncstr
C       write(6,130) krothy - kromax
C       call attn(' CGU3HW: index multiplicity does not match',*110)
C110    continue
C120    format(' number of coefficients for <--;--||HW> ', i4/
C    1         ' number of allowed s-coupling           ', i4/
C    2         ' number of multiplicity                 ', i4/
C    3         ' number of constraints                  ', i4)
C130    format(' number of independent constraints      ', i4)
C     end if
C
      if( krothy .ne. kromax ) then
C
C     Solving the constraint equations
C
      call homgauss( ncstr, krothy, NASI, dmtx, asi )
C
      do njI = 1, nhw
          do st = 1, krothy
            dtempa( st ) = dcgu3e( njI, st )
          end do ! st
        do kro = 1, kromax ! = krothy - ncstr
          dcgu3e( njI, kro ) = 0.d0
          do st = 1, krothy
            dcgu3e( njI, kro ) = dcgu3e( njI, kro )
     1                           + asi( st, kro ) * dtempa( st )
          end do ! st
        end do ! kro
      end do ! njI
      if( krothy-ncstr .ne. kromax ) then
        call attn(' CGU3HW: index multiplicity does not match',*140)
 140    continue
      end if 
C
      end if 
C
C     Resolution of outer rho-multiplicity so the element of the algebra
C     lies on rho=1.
C
      if( kromax .gt. 1 ) then !.and. lm3 .gt. mu3 ) then
      dtemp = dsqrt( dfloat(lm3) / dfloat(lm3+mu3+1) )
      a = dsqrt(dfloat( mu3+2 )) * ( dtemp + dfloat( mu3 )/dtemp )
      b = dsqrt(dfloat( mu3 )) * (-dfloat( mu3+2 )*dtemp + 1.d0/dtemp )
      do njI = 1, nhw
        dcgu3e( njI, 1 ) = a * dcgu3e( njI, 1 ) + b * dcgu3e( njI, 2 )
      end do ! njI
      end if ! kromax
C
C     Gram-Schmidt orthonormalization
C
      call orthon( nhw, kromax, nmax, dcgu3e )
C
900   continue
C
      return
C ---*end of cgu3hw*----------------------------------------------------
      end
