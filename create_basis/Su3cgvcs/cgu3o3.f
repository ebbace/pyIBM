CB----------------------------------------------------------------------
C                       **************************
C                       ***  CG SU(3) > SO(3)  ***
C                       **************************
C ----------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 05/99)
C ----------------------------------------------------------------------
C Update: 05/99 ... original code
C         08/99 ... change for binomial coefficients (CB)
C ----------------------------------------------------------------------
C     Subroutine to calculate the Clebsch Gordan coefficients for SU(3)
C     > SO(3).
C
C     <(lm1 mu1) al1 L1; (lm2 mu2) al2 L2 || rho (lm3 mu3) al3 L3 >
C       = sqrt(8 pi^2/2*L+1) Sum_{M;1,2} Sum_{K;1,2,3>=0} Sum_{jIN;1,2}
C         * ( L1 M1; L2 M2 | L3 K3 ) Kbar^(lm3 mu3)_K3,al3(L3)
C         * < 1 | 1 > Kbar^(lm1 mu1)_K1,al1(L1)
C         * < 2 | 2 > Kbar^(lm2 mu2)_K2,al2(L2)( I1 N1; I2 N2 | s3 s3 )
C         * <(lm1 mu1) j1 I1; (lm2 mu2) j2 I2 || rho (lm3 mu3) 0 s3 >
C     with
C     < 1 | 1 > = < K1 L1 M1 | j1 I1 N1 >
C
C ----------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C ----------------------------------------------------------------------
      subroutine cgu3o3( lm1,mu1, lm2,mu2, lm3,mu3, L1,L2,L3,
     1                   mfmax, ifu3e, dcgu3e, dcgu3o3, nmax )
      implicit none
C
      integer NRO
C
      parameter( NRO = 10 )
C
      integer lm1,mu1, lm2,mu2, lm3,mu3         ! SU(3) irrep
      integer L1, L2, L3                        ! SO(3) sub-irrep
      integer nmax                              ! physical dimension
      integer mfmax                             ! # of 3-hw CGs
      integer ifu3e( * )                        !   their SU(2) sub-IR
      double precision dcgu3e( nmax, * )        !   their values
      double precision dcgu3o3( NRO,NRO,NRO,NRO )! CGs of SU(3) > SO(3)
C
      integer kromax                            ! outer multiplicity
      integer kp1max, kp2max, kp3max            ! inner multiplicity
      integer K1start, K2start, K3start         ! K values
      integer K3max
      integer M1min, M1max, M2max               ! limits of M1, M2
C
      integer s3t
      integer j1t, I1t, N1t, j2t, I2t, N2t      ! SU(2) sub-irrep
C
C Parameter constants
      integer MXK, NDIM
      double precision PIV8, ZERO
      parameter( NDIM = 150, MXK = 10 )
      parameter( PIV8 = 8.885765876317d0, ZERO = 1.d-12 )
C
C Indices and temporary variables
      integer njI, kro, icount !, ie1, ie2
      integer ial1, ial2, ial3
      integer kp1, kp2, kp3
      integer K1, K2, K3
      integer M1, M2
      integer lkey                              ! temporary key
      integer iphase1, iphase2                  ! imaginary overlap phi
CT    integer i, j
      double precision dwigIN, dwigLM, factor
C
      double precision xkbar1( MXK, MXK ), xkbar2( MXK, MXK ),
     1                 xkbar3( MXK, MXK )
      double precision xk1( NDIM, 2 )
      double precision ovrlap1, ovrlap2 ! ovrkbar1, ovrkbar2, ovrkbar3
C
      integer nbit8
      integer index, iover, ibits
C
C Extrinsic functions
      double precision ovrlap
      double precision dwr3
      double precision dsqrt, dfloat, dabs
      integer ishft, iand
      integer max0, min0, mod, iabs
      logical btest
C
C Intrinsic functions
C
      integer kmult, kstart
      integer iupac
C
C     K-multiplicity
      kmult( lm3, mu3, L3 ) = max0( 0, (lm3 + mu3 + 2 - L3)/2 )
     1    - max0( 0, (lm3 + 1 - L3)/2 ) - max0( 0, (mu3 + 1 - L3)/2 )
C
      kstart( lm3, mu3, L3 ) = mod( mu3, 2 ) + 2 * ( max0( 0,
     1    (L3 - lm3)/2 ) + mod( mu3+1, 2 ) * mod( iabs(L3 - lm3), 2 ) )
C
C     Packing and unpacking
      data nbit8  / 255 /            ! zff
      iupac( index, iover, ibits ) = iand( ishft( index,iover ), ibits )
C
C     Initializations
C
      s3t = mu3
      kp1max = kmult( lm1, mu1, L1 )
      if( kp1max .eq. 0 ) return
      kp2max = kmult( lm2, mu2, L2 )
      if( kp2max .eq. 0 ) return
      kp3max = kmult( lm3, mu3, L3 )
      if( kp3max .eq. 0 ) return
      K1start = kstart( lm1, mu1, L1 )
      K2start = kstart( lm2, mu2, L2 )
      K3start = kstart( lm3, mu3, L3 )
      K3max   = K3start + 2 * ( kp3max - 1 )
      call u3mult( lm1,mu1, lm2,mu2, lm3,mu3, kromax, *100 )
      do kro  = 1, kromax
      do ial3 = 1, kp3max
      do ial2 = 1, kp2max
      do ial1 = 1, kp1max
        dcgu3o3( ial1, ial2, ial3, kro ) = 0.d0
      end do ! 3
      end do ! 2
      end do ! 1
      end do ! rho
C
C     ... K-matrix
C
      call  kmat( lm1, mu1, L1, 150, xk1 )
        index = 0
        do kp1 = 1, kp1max
CU        do ial1 = kp1, kp1max           ! Upper row
          do ial1 = 1, kp1                ! Lower row
            index = index + 1
            xkbar1( kp1, ial1 ) = xk1( index, 2 )
          end do ! ial1
        end do                            ! Lower
        do kp1 = 1, kp1max                ! Lower
CU        do ial1 = 1, kp1 - 1            ! Upper row
          do ial1 = kp1+1, kp1max         ! Lower row
            xkbar1( kp1, ial1 ) = xkbar1( ial1, kp1 )
          end do ! ial1
        end do ! kp1
C
      call  kmat( lm2, mu2, L2, 150, xk1 )
        index = 0
        do kp2 = 1, kp2max
CU        do ial2 = kp2, kp2max           ! Upper row
          do ial2 = 1, kp2                ! Lower row
            index = index + 1
            xkbar2( kp2, ial2 ) = xk1( index, 2 )
          end do ! ial2
        end do                            ! Lower
        do kp2 = 1, kp2max                ! Lower
CU        do ial2 = 1, kp2 - 1            ! Upper row
          do ial2 = kp2+1, kp2max         ! Lower row
            xkbar2( kp2, ial2 ) = xkbar2( ial2, kp2 )
          end do ! ial2
        end do ! kp2
C
      call  kmat( lm3, mu3, L3, 150, xk1 )
        index = 0
        do kp3 = 1, kp3max
CU        do ial3 = kp3, kp3max           ! Upper row
          do ial3 = 1, kp3                ! Lower row
            index = index + 1
            xkbar3( kp3, ial3 ) = xk1( index, 2 )
          end do ! ial3
        end do                            ! Lower
        do kp3 = 1, kp3max                ! Lower
CU        do ial3 = 1, kp3 - 1            ! Upper row
          do ial3 = kp3+1, kp3max         ! Lower row
            xkbar3( kp3, ial3 ) = xkbar3( ial3, kp3 )
          end do ! ial3
        end do ! kp3
C
C     Calculations begin ...
C
Cr    do kro = 1, kromax
C
C       Retrieve the highest weight state of 3rd entry
C
        factor =  PIV8 / dsqrt( dfloat(2*L3 + 1) )
C
        do njI = 1, mfmax
C
          lkey = ifu3e( njI )
          j1t  = iupac( lkey, -24, nbit8 )
          I1t  = iupac( lkey, -16, nbit8 )
          j2t  = iupac( lkey,  -8, nbit8 )
          I2t  = iupac( lkey,   0, nbit8 )
C
C         M1max = min0( L1, I1t )
C         if( btest( I1t+M1max, 0 ) ) M1max = M1max - 1
          M2max = min0( L2, I2t )
          if( btest( I2t+M2max, 0 ) ) M2max = M2max - 1
C
C         Run over Ns
C
          N2t = s3t + I1t + 2
          do N1t = -I1t, I1t, 2
          N2t = N2t - 2
 1200     format(' (',2(i2,'/2'),';',2(i2,'/2'),a,2(i2,'/2'),' ) =',
     1           f16.5)
C
          if( iabs( N2t ) .le. I2t ) then
C
            dwigIN = dwr3( I1t, I2t, mu3, N1t, N2t, mu3 )
C
            iphase1 = (mu1 + j1t - N1t)/2
            iphase2 = (mu2 + j2t - N2t)/2
C           write(6,1010) 2, mu2, j2t, N2t, iphase2
C1010       format(' *phase ', i1, 3(i3,'/2'), ' :', i3)
C
C           Run over Ms
C
            do M2 = -M2max, M2max, 2
C
            M1min = max0( K3start - M2, -L1, -I1t )
            M1max = min0( K3max - M2, L1, I1t )
            if( btest( I1t+M1min, 0 ) ) M1min = M1min + 1
            if( btest( I1t+M1max, 0 ) ) M1max = M1max - 1
            do M1 = M1min, M1max, 2
C
C1020         format(' M1 =', i3, ' M2 =', i3)
C
C             Run over Ks
C
              K2 = K2start - 2
              do kp2 = 1, kp2max
              K2 = K2 + 2
                ovrlap2 = ovrlap( lm2,mu2, K2,L2,M2, j2t,I2t,N2t )
                icount = icount + 1
CW              ie2 = 2*lm2 + mu2 - 3*j2t
C               write(6,9000) K2,L2,M2,ie2,I2t,N2t,icount,2
C    1             , ovrlap2
C9000 format(' Calling for OVRLAP <',3i4,'|',3i4,'>', i10,'(',i1,')',
CW   2               f12.5 )
                if( dabs( ovrlap2 ) .gt. ZERO ) then
              do ial2 = 1, kp2max
C
              K1 = K1start - 2
              do kp1 = 1, kp1max
              K1 = K1 + 2
                ovrlap1 = ovrlap( lm1,mu1, K1,L1,M1, j1t,I1t,N1t )
                icount = icount + 1
CW              ie1 = 2*lm1 + mu1 - 3*j1t
C               write(6,9000) K1,L1,M1,ie1,I1t,N1t,icount,1
CW   1             , ovrlap1
                if( dabs( ovrlap1 ) .gt. ZERO ) then
              do ial1 = 1, kp1max
C
              K3 = K3start - 2
              do kp3 = 1, kp3max
              K3 = K3 + 2
                if( K3 .eq. M1+M2 ) then
C
                dwigLM = dwr3( 2*L1,2*L2,2*L3, 2*M1,2*M2,2*K3 ) 
C
              do ial3 = 1, kp3max
C
CF              factor = PIV8 / dsqrt( dfloat(2*L3 + 1) )
C    1            * xkbar3(kp3,ial3)
C    2            * dwr3( 2*L1,2*L2,2*L3, 2*M1,2*M2,2*K3 ) * dwigIN
C
                  do kro = 1, kromax
C
                  if( btest( iphase1+iphase2, 1 ) ) then
                    dcgu3o3( ial1, ial2, ial3, kro )
     1                = dcgu3o3( ial1, ial2, ial3, kro ) - ovrlap1
     2                  * xkbar1(kp1,ial1) * ovrlap2 * xkbar2(kp2,ial2)
     3                  * factor * dwigIN * dwigLM * dcgu3e( njI, kro )
     4                  * xkbar3(kp3,ial3)
                  else
                    dcgu3o3( ial1, ial2, ial3, kro )
     1                = dcgu3o3( ial1, ial2, ial3, kro ) + ovrlap1
     2                  * xkbar1(kp1,ial1) * ovrlap2 * xkbar2(kp2,ial2)
     3                  * factor * dwigIN * dwigLM * dcgu3e( njI, kro )
     4                  * xkbar3(kp3,ial3)
C    3                  * factor * dcgu3e( njI, kro )
                  end if ! phase
C
                end do ! kro
C
              end do ! ial3
                end if ! K3
              end do ! K3
C
              end do ! ial1
                end if
              end do ! K1
C
              end do ! ial2
                end if
              end do ! K2
C
            end do ! M1
C
            end do ! M2
C
          end if ! N2t
C
          end do ! N1t
C
        end do ! njI
C     write(6,*) ' icount ', icount
C
Cr    end do ! kro
C       write(6,*) ' j1t,I1t,j2t,I2t,0,s3t', j1t,I1t,j2t,I2t,0,s3t
C       write(6,*) njI, kro, dcgu3e( njI, kro )
C             write(6,*) ' >> N1t N2t', N1t, N2t
100   return
C
C ---*end of cgu3o3*----------------------------------------------------
      end
