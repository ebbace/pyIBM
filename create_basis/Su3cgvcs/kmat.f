CB----------------------------------------------------------------------
C                       **************************
C                       ***     K - MATRIX     ***
C                       **************************
C ----------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 03/99)
C ----------------------------------------------------------------------
C Update: 03/99 ... original code
C         08/99 ... change for binomial coefficients (CB)
C         05/02 ... change for lookup table (CB)
C ----------------------------------------------------------------------
C     Subroutine to calculate the K-matrix and its inverse for SU(3) 
C     > SO(3) directly.
C
C     S^L(K,K') = 4 Sum_nu bino(mu,nu) I(mu-nu,nu,K) I(mu-nu,nu,K')
C                 x Sum_n^(L-K) d(K,l,K',n) Sum_r (-)^r (lm+nu,r)
C                 x    (lm+nu+L,(K-K'+2n+2r)/2)^(-1)
C
C ----------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C ----------------------------------------------------------------------
      subroutine kmat( lm, mu, L, idimk, xk )
C
      implicit none
C
      integer lm, mu, L
      integer idimk
      double precision xk( idimk, * )
C
C Parameter constants
      integer MXK
      double precision PI, TWOPI, ZERO, ONE, TWO
      parameter ( ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0 )
      parameter ( MXK = 10, PI = 3.14159265359D0, TWOPI = TWO*PI )
C
C Indices and temporary variables
      integer kpmax, mli
      integer ikl, ikr, Kl, Kr
      integer i1, i2, i3, index
      integer nu, n, p, q
      double precision sumnu, sumn, factor
      double precision phmu, phKl, phKr           ! phase
      double precision x, xi, z
      double precision xk1( MXK, MXK ), xk2( MXK, MXK )
      double precision e( MXK ), y( MXK )
C
C Extrinsic functions
      logical btest
      integer max0, mod
      double precision dsqrt, dfloat, dexp
      double precision dmat
C
C Common blocks from file # 4
      integer NFACT, NBINO, NTWOS, N2BIN, NMBIN, NSQR
      double precision dlogf, dlnbin, twos, icoef, scoef
      parameter( NFACT = 2000, NBINO = 128, NTWOS = 128,
     1           N2BIN  = NBINO*(2*NBINO+3), NMBIN  = NBINO*(NBINO+3)/2,
     2           NSQR = NBINO*(NBINO+2) )
      common / BKDF / dlogf(0:NFACT)
      common / BBIN / dlnbin(0:N2BIN), twos(-NTWOS:NTWOS)
      common / DINT / icoef( 0:NMBIN, 0:2*NBINO )
      common / DSNT / scoef( 0:NSQR,  0:2*NBINO )
C     
C Intrinsic functions
C
      integer ind, snd
      integer kmult, kstart
C
C     Indexing for I- and S-functions
      ind( p, q ) = p*(p+1)/2 + q
      snd( p, q ) = p*(NBINO+1) + q
C
C     d-function
c	dmat( a, b, c, d ) = dexp( 
c	1   0.5d0*( dlnbin(ind(2*b,b+a)) - dlnbin(ind(2*b,b+c)) ) 
c    2   + dlnbin(ind(b-a,d)) + dlnbin(ind(b+a,b+c-d)) )
C
C     K-multiplicity
      kmult( lm, mu, L ) = max0( 0, (lm + mu + 2 - L)/2 )
     1    - max0( 0, (lm + 1 - L)/2 ) - max0( 0, (mu + 1 - L)/2 )
C
      kstart( lm, mu, L ) = mod( mu, 2 ) + 2 * ( max0( 0, (L - lm)/2 )
     1    + mod( mu+1, 2 ) * mod( iabs(L - lm), 2 ) )
C
C Initializations
      mli = 1
      kpmax = kmult( lm, mu, L )
C
C Checking the range of I- and S- functions
      if( lm+mu .gt. NBINO ) call error(' KMAT: increase NBINO.' )
      if( L     .gt. NBINO ) call error(' KMAT: increase NBINO.' )
C
C ---------------------------------------------------------------------
C
C Calculating S-matrix
C
C     icount = 0
      if( btest( mu, 0 ) ) then
	  phmu = -ONE
      else
	  phmu = +ONE
      end if
      Kl   = kstart( lm, mu, L ) - 2
      phKl = -phmu
      do ikl = 1, kpmax
        Kl = Kl + 2
	  phKl = -phKl
      Kr = kstart( lm, mu, L ) - 2
      phKr = -phmu
      do ikr = 1, kpmax
        Kr = Kr + 2
	  phKr = -phKr
C
        sumnu = ZERO
C
        do nu = 0, mu/2
C
          sumn = ZERO
          do n = max0( 0, Kr-Kl ), L-Kl
            q = n + (Kl - Kr)/2
	    sumn = sumn + dmat(2*Kl,2*L,2*Kr,2*n)*scoef(snd(lm+nu,L),q)
          end do ! n
C
	  sumnu = sumnu + sumn * icoef( ind(mu-nu,nu), (mu+Kl)/2 ) 
     1            * icoef( ind(mu-nu,nu), (mu+Kr)/2 )/dfloat(lm+nu+L+1)
     2            * dexp( dlnbin( ind(mu,nu) ) )
	  end do
C
        do nu = mu/2+1, mu
C
          sumn = ZERO
          do n = max0( 0, Kr-Kl ), L-Kl
            q = n + (Kl - Kr)/2
	    sumn = sumn + dmat(2*Kl,2*L,2*Kr,2*n)*scoef(snd(lm+nu,L),q)
          end do ! n
C
	  sumnu = sumnu + sumn * icoef( ind(nu,mu-nu), (mu+Kl)/2 ) 
     1            * icoef( ind(nu,mu-nu), (mu+Kr)/2 )/dfloat(lm+nu+L+1)
     2            * phKl * phKr * dexp( dlnbin( ind(mu,nu) ) )
	  end do
C
        xk1( ikl, ikr ) = TWO * sumnu
C
      end do ! Kr
      end do ! Kl
C
C ---------------------------------------------------------------------
C
C Calculating K-matrix and its inverse
C
C     Taking the square root of S-matrix
C
      factor = TWOPI !* twos(-mu )
      if( kpmax .eq. 1 ) then
        xk1( 1, 1 )  = factor * dsqrt( xk1( 1, 1 ) )
        xk( mli, 1 ) = xk1( 1, 1 )
        xk( mli, 2 ) = ONE / xk1( 1, 1 )
      else
C
C ...   square root
C
C     -----------------------------------------------------------------
          do i1 = 1, kpmax
C         do i2 = 1, kpmax
          do i2 = 1, i1
            xk2( i1, i2 ) = xk1( i1, i2 )
          end do ! i2
          end do ! i1
          call tred2( xk2, kpmax, MXK, e, y )
          call tqli( e, y, kpmax, MXK, xk2 )
C     -----------------------------------------------------------------
          do i1 = 1, kpmax
            e( i1 ) = factor * dsqrt( e( i1 ) )            ! eig. K
          end do ! i1
 1100     format(/' ** Eigenvalues of K ** '/5e15.4)
C
C         Upper triangle of matrix by row
          index = -1
          do i1 = 1, kpmax
CU        do i2 = i1, kpmax                            ! Upper row
          do i2 = 1, i1                                ! Lower row
            index = index + 1
            x  = ZERO
            xi = ZERO
            do i3 = 1, kpmax
              z  = xk2(i1,i3) * xk2(i2,i3)                      
              x  = x  + z * e(i3)
              xi = xi + z / e(i3)
            end do ! i3
            xk( mli + index, 1 ) = x                   ! new K
            xk( mli + index, 2 ) = xi                  ! new K^-1
          end do ! i2
          end do ! i1
      end if ! kpmax
CW    write(6,*) ' Number of sum:', icount
CW    write(6,'(8f10.4)')(xk(index,1), index=1,kpmax*(kpmax+1))
CW    write(6,'(8f10.4)')(xk(index,2), index=1,kpmax*(kpmax+1))
C
C ---*end of kmat*------------------------------------------------------
      end
