CB----------------------------------------------------------------------
C                       **************************
C                       ***      OVERLAP       ***
C                       **************************
C ----------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 03/99)
C ----------------------------------------------------------------------
C Update: 03/99 ... original code
C         08/99 ... change for binomial coefficients (CB)
C         05/02 ... change for lookup table (CB)
C ----------------------------------------------------------------------
C     Subroutine to calculate the overlap of SO(3) and SU(2) h.w. of
C     SU(3) irrep.
C
C**   < al L M | j I N > = sqrt(2*L+1 / 8 pi^2) Sum_K K^_K,al(L)
C     <  K L M | j I N > = sqrt(2*L+1 / 8 pi^2) * (-)^(I+s-j+2*lm+mu)
C          * Sum_m U( lm+mu-2j, 2j+s+m, lm, s-m, lm+s-m, I )
C          *   d( s, s, m, 0 ) * I( s+m, s-m, K )
C          * Sum_mu d( j+m, I, N, mu ) * I( 2I+N-j-m-2mu, j+m-N+2mu, M )
C          * - (-2)^(2j+1)/(2J+L+1) * ( 2J,2j )^(1/2)
C          *   Sum_nu d( K, L, M, nu ) * Sum_r (-)^r ( 2J-2j,r )
C          *   ( 2J+L,(K-M+2r+2nu+2j)/2 )^(-1)
C
C ----------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C ----------------------------------------------------------------------
      function ovrlap( lm, mu, K, L, M, jt, It, Nt )
C
      implicit none
C
      integer lm, mu, K, L, M, jt, It, Nt
      double precision ovrlap
C
C Parameter constants
      integer MXK
      double precision PI, TWOPI, ZERO, ONE, TWO
      parameter ( ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0 )
      parameter ( MXK = 10, PI = 3.14159265359D0, TWOPI = TWO*PI )
C
C Indices and temporary variables
C     integer kpmax
      integer nu, phnu, p, q, r
      integer st, mt
      integer Jbigt !, icount
CPH   integer Jbigt, iphase
      integer a, b, c, d, e, f
      double precision summt, sumnu, sumr, factor
C
C Extrinsic functions
      logical btest
C     integer min0, max0, mod
      integer min0, max0, iabs
      double precision dsqrt, dfloat, dexp
      double precision dmat
C6j   double precision dmat, d6jr3
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
      double precision drr3sp
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
C     Racah coefficients for special arguments
      drr3sp( a, b, e, d, c, f ) = dexp( 0.5d0*( (dlogf(2*a)+dlogf(2*b)+
     1   dlogf(a+b+d+e+2)+dlogf(a+b-d+e)+dlogf(a+b+d-e)+dlogf(-a+e+f)+
     2   dlogf(-b+d+f)) -
     3  (dlogf(2*a+2*b+2)+dlogf(-a-b+d+e)+dlogf(a+e-f)+dlogf(a-e+f)+
     4   dlogf(a+e+f+2)+dlogf(b+d-f)+dlogf(b-d+f)+dlogf(b+d+f+2)) ) )
C
C Initializations
      st = mu
C
C Checking couplings
      if( btest( st+K, 0 ) .or. btest( It+M, 0 ) ) then
C
        ovrlap = ZERO
C
      else
C
C       allowed couplings
        summt = ZERO
c       a = 0
C --*----*----*----*----*----*----*----*----*----*----*----*----*----*--
        do mt = -st, st, 2
C
          if( iabs( jt + mt ) .le. It ) then
C
          Jbigt = lm + (st - mt)/2
          if( Jbigt .ge. iabs(lm+mu-2*jt-(st+mt)/2)
     1        .or. Jbigt .le. lm+mu+(st+mt)/2 ) then
C
            sumnu = ZERO
c	      icount = 0
C ----------------------------------------------------------------------
 1          format(a, ' max=', i6, ' min=', i6 )
            do nu = max0( 0, M-K ), min0( L-K, L+M )
C
c             icount = icount + 1
              sumnu = sumnu + dmat( 2*K, 2*L, 2*M, 2*nu ) * scoef( snd(
     1                Jbigt-jt, jt+L), (K-M+2*nu+jt)/2 )
c	        write(6,*) 'd',2*K,2*L,2*M,2*nu,dmat(2*K,2*L,2*M,2*nu ),
c    1            '  S', Jbigt-jt, jt+L, (K-M+2*nu+jt)/2,
c    2            scoef( snd(Jbigt-jt, jt+L), (K-M+2*nu+jt)/2 ), sumnu
C
            end do ! nu 1
            factor = sumnu
c	    a = a + icount
c	    write(6,*) ' icount', icount, a
C ----------------------------------------------------------------------
C
            q = (jt + mt - Nt)/2
            sumnu = ZERO
c	    icount = 0
	    if( btest( max0( 0,-q ), 0 ) ) then
	      phnu = +ONE
          else
	      phnu = -ONE
          end if
C ----------------------------------------------------------------------
            do nu = max0( 0,-q ), min0( (It-jt-mt)/2, (It+Nt)/2 )
C
c              icount = icount + 1
              phnu = -phnu
	        if( It .ge. 2*q+4*nu ) then
                sumnu = sumnu + phnu * dmat( jt+mt, It, Nt, 2*nu )
     1                  * icoef( ind(It-q-2*nu, q+2*nu), (It+M)/2 )
              else
                if( btest( It+M, 1 ) ) then
                sumnu = sumnu - phnu * dmat( jt+mt, It, Nt, 2*nu )
     1                  * icoef( ind(q+2*nu, It-q-2*nu), (It+M)/2 )
	          else
                sumnu = sumnu + phnu * dmat( jt+mt, It, Nt, 2*nu )
     1                  * icoef( ind(q+2*nu, It-q-2*nu), (It+M)/2 )
	          end if
	        end if
c	      write(6,*) 'd', jt+mt,It,Nt,2*nu,dmat(jt+mt,It,Nt,2*nu),
c    1            '  I1', It-q-2*nu, q+2*nu, (It+M)/2,
c    2            icoef( ind(It-q-2*nu, q+2*nu), (It+M)/2 ), sumnu
C
            end do ! nu 2
            factor = factor * sumnu
c	    a = a + icount
c	    write(6,*) ' icount', icount, a
C ----------------------------------------------------------------------
C
c              icount = icount + 1
	      if( mt .ge. 0 ) then
		    sumnu = + dmat( st, st, mt, 0 ) 
     1                  * icoef( ind( (st+mt)/2,(st-mt)/2 ), (st+K)/2 )
	      else
	        if( btest( st+K, 1 ) ) then
		      sumnu = - dmat( st, st, mt, 0 ) 
     1                  * icoef( ind( (st-mt)/2,(st+mt)/2 ), (st+K)/2 )
	        else
		      sumnu = + dmat( st, st, mt, 0 ) 
     1                  * icoef( ind( (st-mt)/2,(st+mt)/2 ), (st+K)/2 )
	        end if
	      end if
c	      write(6,*) 'd', st,st,mt,0, dmat(st,st,mt,0),
c    1            '  I2', (st+mt)/2, (st-mt)/2, (st+K)/2,
c    2            icoef( ind( (st+mt)/2,(st-mt)/2 ), (st+K)/2 ), sumnu
            factor = factor * sumnu
c	    a = a + icount
c	    write(6,*) ' icount', icount, a
C ----------------------------------------------------------------------
C 
            factor = factor*dexp( 0.5d0 * dlnbin( ind(Jbigt,jt) ) )
     1               / dfloat(Jbigt+L+1)
C
C The overall factor and phase are put outside the sum over mt
C     factor = factor * TWOPI * dsqrt( 0.5d0*dfloat(2*L+1) )
C                     * twos(jt+1)
C     factor = -(-)^(jt+1) * i^((st+jt-Nt)/2) * factor
C
C --*----*----*----*----*----*----*----*----*----*----*----*----*----*--
            summt = summt + factor * dsqrt( dfloat( (Jbigt+1)*(It+1) ) )
     1          * drr3sp( (st-mt)/2,lm,jt+(st+mt)/2,lm+mu-jt, Jbigt,It )
C --*----*----*----*----*----*----*----*----*----*----*----*----*----*--
C
C     summt = (-)^(lm+(It+st-jt)/2) * summt
C
          end if ! Jbigt
C
          end if ! It
C
        end do ! mt
C
        if( btest( (It+st-jt)/2+ 2*lm+mu + jt, 0 ) ) then
          ovrlap = -summt * TWOPI * dsqrt( 0.5d0*dfloat(2*L+1) )
     1                * twos(jt+1)
        else
          ovrlap =  summt * TWOPI * dsqrt( 0.5d0*dfloat(2*L+1) )
     1                * twos(jt+1)
        end if ! phase factor
C
C       imaginary phase factor is put outside the subroutine
CPH     if( btest( st+jt-Nt, 0 ) ) then
C         write(6,*) ' Warning ... coupling ...'
C       else
C         iphase = (st + jt - Nt)/2
C         if( mod( iphase, 4 ) .eq. 0 ) then
C           write(6,1010) K,L,M,jt,It,Nt,ovrlap, '(+)',lm,mu
C         else if( mod( iphase, 4 ) .eq. 1 ) then
C           write(6,1010) K,L,M,jt,It,Nt,ovrlap, '(+i)',lm,mu
C         else if( mod( iphase, 4 ) .eq. 2 ) then
C           write(6,1010) K,L,M,jt,It,Nt,ovrlap, '(-)',lm,mu
C         else if( mod( iphase, 4 ) .eq. 3 ) then
C           write(6,1010) K,L,M,jt,It,Nt,ovrlap, '(-i)',lm,mu
C         end if
CPH     end if ! iphase
C
      end if ! st, It
 1010 format(': <',3i3,'|',3(i2,'/'),'>',t25,f10.5,a,' for (',2i3,')')
C
C ---*end of ovrlap*----------------------------------------------------
      end
