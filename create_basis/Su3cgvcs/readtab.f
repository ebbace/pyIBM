CB----------------------------------------------------------------------
C                       **************************
C                       ***       Tables       ***
C                       **************************
C ----------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto/LSU, 05/02)
C ----------------------------------------------------------------------
C Update: 05/02 ... original code
C ----------------------------------------------------------------------
C     Subroutine to generate lookup tables for SU(3) CG packages:
C
C         - gamma function (factorial)
C         - inverse beta function (binomial)
C         - I(p,q,sigma)
C         - S(p,q,sigma)
C
C     I(p,q,sigma) = 1/2^(p+q) Sum_n (-)^n binom(p, sigma-n) binom(q,n)
C     S(p,q,sigma) = Sum_n (-)^n binom(p, n) binom(p+q,sigma-n)^(-1)
C
C     Extended precision version (may not work for some compilers).
C ----------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C ----------------------------------------------------------------------
      subroutine readtab
C
      implicit none
C
      integer NFACT, NBINO, NTWOS, N2BIN, NMBIN, NSQR
      parameter ( NFACT = 2000, NBINO = 128, NTWOS = 128,
     1            N2BIN = NBINO*(2*NBINO+3),
     2            NMBIN = NBINO*(NBINO+3)/2, NSQR = NBINO*(NBINO+2) )
C    1    NMBIN = NBINO*(NBINO+3)/2, NBINI = NMBIN*(NMBIN+1)/2+NBINO,
C    2    NSQR  = NBINO*NBINO,       NBINS = NSQR*(NSQR+1)/2+NBINO )
C
      double precision PI, SQRTPI, LOGPI, LOGHFFACT
C                             1 2345678901234567890123456789012345
      parameter ( PI        = 3.1415926535897932384626433832795259d0,
     1            SQRTPI    = 1.7724538509055160272981674833411452d0,
     2            LOGPI     = 1.1447298858494001741434273513530587d0,
     3            LOGHFFACT =-0.12078223763524522234551844578164721d0 )
      double precision ZERO, ONE, TWO, HALF, QUARTER
      parameter ( ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0,
     1            HALF = 0.5d0, QUARTER = 0.25d0 )
C
      double precision qlogf, qlnbin, qtwos
c     double precision dlogf, dlnbin, dtwos
C
      common / BKDF / qlogf( 0:NFACT )
      common / BBIN / qlnbin( 0:N2BIN ), qtwos(-NTWOS:NTWOS )
c     common / BKDF / dlogf( 0:NFACT )
c     common / BBIN / dlnbin( 0:N2BIN ), dtwos(-NTWOS:NTWOS )
C
      double precision qintg, qsntg
c     double precision dintg, dsntg
      common / DINT / qintg( 0:NMBIN, 0:2*NBINO )
      common / DSNT / qsntg( 0:NSQR,  0:2*NBINO )
c     common / DINT / dintg( 0:NMBIN, 0:2*NBINO )
c     common / DSNT / dsntg( 0:NSQR,  0:2*NBINO )
C
      integer          p, q, sigma, n, icount
C     integer          nmin, nmax
      double precision gen
      double precision phase, sums, tums
C
C     Intrinsic functions
      double precision dfloat, dlog, dexp, dnint
C
C Clock
C
c     integer ibeg, iend
c     integer mclock
c     double precision xbeg, xend, xtime, DSCALE
c     parameter( DSCALE = 1.d-6 )
C
C     Indexing and binomials
      integer ind, snd
C     double precision binom
      ind( p, q ) = p*(p+1)/2 + q
      snd( p, q ) = p*(NBINO+1) + q
C
C ----------------------------------------------------------------------
      icount = 1000000
c     write(6,*) NBINO, NMBIN, NSQR, NBINI, NBINS
c     write(6,*) ' EXTENDED PRECISION '
C
C ----------------------------------------------------------------------
C
C     Factorials (Gamma function) - logarithm
C
C ----------------------------------------------------------------------
c     ibeg = mclock()
c     xbeg = dfloat( ibeg )
      qlogf( 0 ) = ZERO
      qlogf( 1 ) = LOGHFFACT
C
      do p = 2, NFACT, 2
        qlogf( p ) = qlogf( p-2 ) + dlog( dfloat(p/2) )
      end do
C
      do p = 3, NFACT-1, 2
        qlogf( p ) = qlogf( p-2 ) + dlog( HALF * dfloat(p) )
      end do
C
c     iend = mclock()
c     xend = dfloat( iend )
c         write(6,'(a)') ' Factorial '
c         xtime = DSCALE * ( xend - xbeg )
c         write(6,300) 1000000, xtime, icount, xtime/dfloat(icount)
C ----------------------------------------------------------------------
C
C     Binomials (Inverse Beta function) - logarithm
C     integer only
C
C ----------------------------------------------------------------------
c     ibeg = mclock()
c     xbeg = dfloat( ibeg )
      qlnbin( 0 ) = ZERO
      qlnbin( 1 ) = ZERO
C
      do p = 0, 2*NBINO
        do q = 0, p/2
          gen  = qlogf( 2*p ) - qlogf( 2*q ) - qlogf( 2*p-2*q )
          icount = p * (p+1) / 2 + q
          qlnbin( icount ) = gen
          icount = p * (p+3) / 2 - q
          qlnbin( icount ) = gen
        end do ! q
      end do ! p
C
c     iend = mclock()
c     xend = dfloat( iend )
c         write(6,'(a)') ' Binomial '
c         xtime = DSCALE * ( xend - xbeg )
c         write(6,300) 1000000, xtime, icount, xtime/dfloat(icount)
C ----------------------------------------------------------------------
C
C     Power of 2 (2^N)
C
C ----------------------------------------------------------------------
c     ibeg = mclock()
c     xbeg = dfloat( ibeg )
      qtwos( 0 ) = ONE
      do p = 1, NTWOS
        qtwos( +p ) = qtwos( p-1 ) * TWO
        qtwos( -p ) = qtwos( 1-p ) / TWO
      end do ! i
C
c     iend = mclock()
c     xend = dfloat( iend )
c         write(6,'(a)') ' Power of 2 '
c         xtime = DSCALE * ( xend - xbeg )
c         write(6,300) 1000000, xtime, icount, xtime/dfloat(icount)
C ----------------------------------------------------------------------
C
C     I( p, q, sigma ) - recurrence of half sum
C
C ----------------------------------------------------------------------
c     ibeg = mclock()
c     xbeg = dfloat( ibeg )
C
      do p = 0, NBINO
        icount = ind( p, 0 )
        do sigma = 0, p
          qintg( icount, sigma ) = dnint(dexp(qlnbin(ind( p, sigma ))))
     1                             * qtwos(-p )
        end do
      end do
      icount = -1
      do p = 0, NBINO
        icount = icount + 1
        if( icount .ne. ind(p,0) ) write(6,*) ' NOT OK '
      do q = 1, p
        icount = icount + 1
        if( icount-p-1 .ne. ind(p-1,q-1) ) write(6,*) ' NOT OK FURTHER'
        qintg( icount, 0 ) = qintg( icount-p-1, 0 )
     1                             * QUARTER
        qintg( icount, 1 ) = qintg( icount-p-1, 1 )
     1                             * QUARTER
        do sigma = 2, p+q
          qintg( icount, sigma ) = ( qintg( icount-p-1, sigma )
     1                             - qintg( icount-p-1, sigma-2 ) )
     2                             * QUARTER
        end do
      end do
      end do
C
c     iend = mclock()
c     xend = dfloat( iend )
c         write(6,'(a)') ' I(p,q,sigma)'
c         xtime = DSCALE * ( xend - xbeg )
c         write(6,300) 1000000, xtime, icount, xtime/dfloat(icount)
C ----------------------------------------------------------------------
C
C     S( p, q, sigma )
C
C ----------------------------------------------------------------------
c     ibeg = mclock()
c     xbeg = dfloat( ibeg )
C
      icount = -1
      do q = 0, NBINO
        icount = icount + 1
        do sigma = 0, q
          qsntg( icount, sigma ) = ONE/dnint(dexp(qlnbin(ind(q,sigma))))
        end do
      end do
      do p = 1, NBINO
      do q = 0, NBINO-1
        icount = icount + 1
	do sigma = 0, q-1           ! product
	  qsntg( icount, sigma ) = (dfloat(q - sigma) * qsntg( icount-1, 
     1           sigma ) + dfloat(p) * qsntg( icount-NBINO-1, sigma ))
     2           / dfloat(p + q)
	end do
        do sigma = q, p+q-1       ! division
          qsntg( icount, sigma ) = qsntg( icount-NBINO, sigma ) 
     1                             - qsntg( icount-NBINO, sigma+1 )
        end do
          qsntg( icount, p+q ) = ONE
      end do
	q = NBINO
        icount = icount + 1
	do sigma = 0, q-1           ! product
	  qsntg( icount, sigma ) = (dfloat(q - sigma) * qsntg( icount-1, 
     1           sigma ) + dfloat(p) * qsntg( icount-NBINO-1, sigma ))
     2           / dfloat(p + q)
	end do
	  phase = ONE
	  sums  = ZERO
	  do n = 0, p
	    sums  = sums + phase*dexp( qlnbin(ind( sigma+n, n )) )
	    phase = -phase
	  end do
	  tums  = dexp( qlnbin(ind( p+sigma, p )) )
	  qsntg( icount, q ) = dnint(sums) / dnint(tums)
	do sigma = q+1, p+q         ! division
	  phase = ONE
	  sums  = ZERO
	  do n = 0, p+q-sigma
	    sums = sums + phase*dexp( qlnbin(ind( p, n ))
     1                               -qlnbin(ind( p+q, sigma+n )) )
	    phase = -phase
	  end do
	  qsntg( icount, sigma ) = sums
	end do
          qsntg( icount, p+q ) = ONE
      end do
c     iend = mclock()
c     xend = dfloat( iend )
c         write(6,'(a)') ' S(p,q,sigma)'
c         xtime = DSCALE * ( xend - xbeg )
c         write(6,300) 1000000, xtime, icount, xtime/dfloat(icount)
C ----------------------------------------------------------------------
C
c     icount = -1
c     do p = 0, NBINO
c     do q = 0, p
c        icount = icount + 1
c        if( icount .ne. ind(p,q) ) write(51,*) ' NOT OK'
c     do sigma = 0, p+q
c        write(51,*) p, q, sigma, qintg( icount, sigma )
c     end do
c     end do
c     end do
C
c     icount = -1
c     do p = 0, NBINO
c     do q = 0, NBINO
c        icount = icount + 1
c        if( icount .ne. snd(p,q) ) write(53,*) ' NOT OK'
c     do sigma = 0, p+q
c        write(53,*) p, q, sigma, icount, qsntg( icount, sigma )
c     end do
c     end do
c     end do
100   format( i6, 2(1pd22.12), 1pd22.12 )
110   format(2i6, 2(1pd22.12), 1pd22.12 )
120   format(3i6, 2(1pd22.12), 1pd22.12 )
300       format(/1x, i4, 20x, ' T o t a l :',e12.5,' seconds'/
     1    i15, ' calc. ',  3x, ' A v e .   :',e12.5,' seconds'/)
310       format(a,'time:', 2i5, e12.5)
C
      return
c     stop
C ---*end of readtab*---------------------------------------------------
      end
