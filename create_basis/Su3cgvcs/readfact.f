      subroutine readfact
C     -----------------------------------------------------------------
C     Binomials (inverse beta) and factorials (gamma).
C     Gamma for integer and half-integer values.
C     -----------------------------------------------------------------
C     Author: C. Bahri     (UofT, 3/99) - inclusion of half-integers
C             based on su3genbk by J.P. Draayer
C     -----------------------------------------------------------------
      implicit double precision(d) 
      double precision PI, SQRTPI, LOGPI, LOGHFFACT
      parameter( NFACT = 2000, NBINO = 400, NTWOS = 128,
     1           NMBIN = NBINO*(NBINO+3)/2 )
C
      common / BKDF / dlogf( 0:NFACT )
      common / BBIN / dlnbin( 0:NMBIN ), dtwos(-NTWOS:NTWOS )
C
      parameter ( PI = 3.141592653589793d0, 
     1            SQRTPI = 1.772453850905516d0, 
     2            LOGPI = 1.144729885849400d0, 
     3            LOGHFFACT = -0.1207822376352453d0 )
C
C     Factorials (Gamma function) - logarithm
C
      dlogf( 0 ) = 0.d0
      dlogf( 1 ) = LOGHFFACT
C
      do i = 2, 2000, 2
	dlogf( i ) = dlogf( i-2 ) + dlog( dfloat(i/2) )
      end do
C
      do i = 3, 1999, 2
	dlogf( i ) = dlogf( i-2 ) + dlog( 0.5 * dfloat(i) )
      end do 
C
C     Binomials (Inverse Beta function)
C
      do i = 0, 2
        dlnbin( i ) = 0.d0
      end do ! i
C
      do i = 2, NBINO
	do j = 0, i/2
	  index = i * (i+1) / 2 + j
	  dgen  = dlogf( i ) - dlogf( j ) - dlogf( i-j ) 
	  dlnbin( index ) = dgen
	  index = i * (i+3) / 2 - j
	  dlnbin( index ) = dgen
	end do ! j
      end do ! i
C
C     Power of 2 (2^N)
C
      dtwos( 0 ) = 1.d0
      do i = 1, NTWOS
	dtwos( +i ) = dtwos( i-1 ) * 2.d0
	dtwos( -i ) = dtwos( 1-i ) / 2.d0
      end do ! i
C
      return
C----*end of readfact*-------------------------------------------------
      end
