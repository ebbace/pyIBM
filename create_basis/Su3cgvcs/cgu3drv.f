      program cgu3drv
CB---------------------------------------------------------------------
C                        ********************
C                        ***  CGU3 DRIVER ***
C                        ********************
C ---------------------------------------------------------------------
C Update:  11/98 ... original code by C. Bahri
C ---------------------------------------------------------------------
C     Driver to calculate the Clebsch-Gordan coefficients of SU(3)
C     in SU(2) bases.
C ---------------------------------------------------------------------
C References:
C  1. D.J. Rowe and C. Bahri, J. Math. Phys. (2000)
C  2. C. Bahri, D.J. Rowe, and J.P. Draayer, Comp. Phys. Comm. (in
C     press)
C ---------------------------------------------------------------------
C
      implicit none
      integer NMAX, NRO, NLBLMAX
      character*57 HEADER
      parameter ( NMAX = 1000, NRO = 10, NLBLMAX = 10000 )
      parameter ( HEADER = 
     5  ' ((lm1 mu1)j1tI1t; (lm2 mu2)j2tI2t || rh (lm3 mu3)j3tI3t)' )
C                 1         2         3         4         5
C        123456789012345678901234567890123456789012345678901234567
C
C Arrays
      integer          ifu3e( NMAX )            ! packed labels
      double precision dcgu3e( NMAX * NRO )     ! seed coefficients
      double precision dcgu3( NRO )             ! CG coefficients
C
C Indices and temporary variables
      logical gosu3, gosu2, gorite
      integer iopt                              ! option number
      integer idx, icount                       ! index and counter
      integer kro
      integer nbit8
      character*2 ans, no, yes
      data yes/'Y'/, no/'N'/
C
C Irreps and weights
      integer lm1,mu1, lm2,mu2, lm3,mu3, kromax ! irreps
      integer j1t,I1t, j2t,I2t, j3t,I3t         ! weights
      integer ie1,     ie2,     ie3             ! alt weights
      integer nhw                               ! # of seed coefficients
C
C Arrays of weights
      integer iw1a( NLBLMAX ), iw2a( NLBLMAX ), iw3a( NLBLMAX )
      integer lbl1, lbl2, lbl3
      integer lbl1max, lbl2max, lbl3max
C
C Clock
C
cc    integer ibeg, iend
cc    integer mclock
cc    double precision xbeg, xend, xtime, dfloat, DSCALE
cc    parameter( DSCALE = 1.d-6 )
C
C External functions
      logical chkwght
      integer iabs
C
C Intrinsic functions
      integer ipack, iupac
C
C     Packing and unpacking
      data nbit8  / 255 /            ! zff
      ipack( j1t, j2t, j3t, kro ) = ior(kro, ishft( ior(j3t, ishft(
     1        ior(j2t, ishft( j1t, 8 )),8 )),8 ))
      iupac( j1t, j2t, j3t ) = iand( ishft( j1t, j2t ), j3t )
C
C     Initializations
C     call blocks                               ! load factorials
      call readfact                             ! load factorials
      call intro                                ! introductory page
      call option( iopt )                       ! choose option
C
C     Proceed ...
      gosu3 = .true.
C
      do while ( gosu3 )
        write(6,50) ' Enter the irreps: lm1,mu1, lm2,mu2, lm3,mu3::'
        read(5,*)   lm1,mu1, lm2,mu2, lm3,mu3
        call u3mult( lm1,mu1, lm2,mu2, lm3,mu3, kromax, *5 )
5       if( kromax .eq. 0 ) then
          write(6,100) HEADER, lm1,mu1, lm2,mu2, lm3,mu3
          go to 40
        end if ! kromax
C
cc      ibeg = mclock()
cc      xbeg = dfloat( ibeg )
C
C ----*------*------*------*------*------*------*------*------*------*-
        call cgu3hw( lm1,mu1, lm2,mu2, lm3,mu3, nhw,
     1                   ifu3e, dcgu3e, NMAX )
C ----*------*------*------*------*------*------*------*------*------*-
C
cc      iend = mclock()
cc      xend = dfloat( iend )
cc      xtime = DSCALE * ( xend - xbeg )
C       write(6,200) 8888, xtime, nhw, xtime/dfloat(nhw)
C200    format(/1x, i4, 20x, '>T o t a l :',e12.5,' seconds'/
C    1  i15, ' nhwc. ',  3x, ' A v e .   :',e12.5,' seconds'/)
C
        gosu2 = .true.
C
        do while ( gosu2 )
          gorite = .true.
C
cc        ibeg = mclock()
cc        xbeg = dfloat( ibeg )
C         write(6,310) ' Initial ', ibeg, mclock(), xbeg
C
C Options
C ---------------------------------------------------------------------
          if( iopt .eq. 1 ) then
C ---------------------------------------------------------------------
C
            write(6,50) ' Enter: 2*j1,2*I1, 2*j2,2*I2, 2*j3,2*I3::'
            read(5,*)   j1t,I1t, j2t,I2t, j3t,I3t
C
            ie1 = 2*lm1 + mu1 - 3*j1t
            ie2 = 2*lm2 + mu2 - 3*j2t
            ie3 = 2*lm3 + mu3 - 3*j3t
C
            if( ie3 .ne. ie1+ie2 .or. .not. chkwght(lm1,mu1,j1t,I1t)
     2          .or. .not. chkwght(lm2,mu2,j2t,I2t)
     3          .or. .not. chkwght(lm3,mu3,j3t,I3t)
     4          .or. I3t .lt. iabs(I1t-I2t) .or. I3t .gt. I1t+I2t ) then
              write(6,140) HEADER, lm1, mu1, j1t, I1t,
     1                             lm2, mu2, j2t, I2t,
     2                             lm3, mu3, j3t, I3t
              go to 30
            end if
C
            lbl1max = 1
            lbl2max = 1
            lbl3max = 1
C
            iw1a( 1 ) = ipack( 0, 0, j1t, I1t )
            iw2a( 1 ) = ipack( 0, 0, j2t, I2t )
            iw3a( 1 ) = ipack( 0, 0, j3t, I3t )
C
C ---------------------------------------------------------------------
          else if( iopt .eq. 2 ) then
C ---------------------------------------------------------------------
C
            write(6,50) ' Enter: 2*j1,2*I1, 2*j2,2*I2::'
            read(5,*)   j1t,I1t, j2t,I2t
C
            if( .not. chkwght(lm1,mu1,j1t,I1t)
     1          .or. .not. chkwght(lm2,mu2,j2t,I2t) ) then
              write(6,140) HEADER, lm1, mu1, j1t, I1t,
     1                             lm2, mu2, j2t, I2t,
     2                             lm3, mu3, NMAX, NMAX
              go to 30
            end if
C
            lbl1max = 1
            lbl2max = 1
C
            iw1a( 1 ) = ipack( 0, 0, j1t, I1t )
            iw2a( 1 ) = ipack( 0, 0, j2t, I2t )
C
            call genwght( lm3, mu3, iw3a, lbl3max )
C
C ---------------------------------------------------------------------
          else if( iopt .eq. 3 ) then
C ---------------------------------------------------------------------
C
            write(6,50) ' Enter: 2*j3,2*I3::'
            read(5,*)   j3t,I3t
C
            if( .not. chkwght(lm3,mu3,j3t,I3t) ) then
              write(6,140) HEADER, lm1, mu1, NMAX, NMAX,
     1                             lm2, mu2, NMAX, NMAX,
     2                             lm3, mu3, j3t, I3t
              go to 30
            end if
C
            lbl3max = 1
            iw3a( 1 ) = ipack( 0, 0, j3t, I3t )
C
            call genwght( lm1, mu1, iw1a, lbl1max )
            call genwght( lm2, mu2, iw2a, lbl2max )
C
C ---------------------------------------------------------------------
          else if( iopt .eq. 4 ) then
C ---------------------------------------------------------------------
C
            call genwght( lm1, mu1, iw1a, lbl1max )
            call genwght( lm2, mu2, iw2a, lbl2max )
            call genwght( lm3, mu3, iw3a, lbl3max )
C
          end if
          if( lbl1max .gt. NLBLMAX .or. lbl2max .gt. NLBLMAX
     1                             .or. lbl3max .gt. NLBLMAX )
     2      call error(' CGU3DRV: Increase NLBLMAX.')
C
C         Calculating CGs ...
C
          icount = 0
          do lbl1 = 1, lbl1max
            idx = iw1a( lbl1 )
            j1t = iupac( idx, -8, nbit8 )
            I1t = iupac( idx,  0, nbit8 )
            ie1 = 2*lm1 + mu1 - 3*j1t
          do lbl2 = 1, lbl2max
            idx = iw2a( lbl2 )
            j2t = iupac( idx, -8, nbit8 )
            I2t = iupac( idx,  0, nbit8 )
            ie2 = 2*lm2 + mu2 - 3*j2t
          do lbl3 = 1, lbl3max
            idx = iw3a( lbl3 )
            j3t = iupac( idx, -8, nbit8 )
            I3t = iupac( idx,  0, nbit8 )
            ie3 = 2*lm3 + mu3 - 3*j3t
C
          if( ie3 .eq. ie1+ie2 .and. I3t .ge. iabs(I1t-I2t)
     1       .and. I3t .le. I1t+I2t ) then
C
C ----*------*------*------*------*------*------*------*------*------*-
            call cgu3( lm1,mu1, lm2,mu2, lm3,mu3, j1t,I1t, j2t,I2t,
     1                 j3t,I3t, nhw, ifu3e, dcgu3e, dcgu3, NMAX )
C ----*------*------*------*------*------*------*------*------*------*-
C
            icount = icount + 1
C
            if( gorite ) then
                write(6,120) HEADER, lm1, mu1, j1t, I1t,
     1                       lm2, mu2, j2t, I2t,
     2                       1, lm3, mu3, j3t, I3t, dcgu3(1),
     3                       ( kro, dcgu3(kro), kro=2,kromax )
C               if( logfile.eq.yes ) write(ifl,120) HEADER,
C    1                       lm1, mu1, j1t, I1t,  lm2, mu2, j2t, I2t,
C    2                       1, lm3, mu3, j3t, I3t, dcgu3(1),
C    3                       ( kro, dcgu3(kro), kro=2,kromax )
                gorite = .false.
            else
                write(6,130) j1t, I1t, j2t, I2t, 1, j3t, I3t, dcgu3(1),
     1                       ( kro, dcgu3(kro), kro=2,kromax )
C               if( logfile.eq.yes ) write(ifl,130) j1t,I1t, j2t,I2t,
C    1              1,j3t,I3t,dcgu3(1), ( kro,dcgu3(kro), kro=2,kromax )
            end if ! writing
C
          end if ! allowed coupling
C
          end do ! lbl3
          end do ! lbl2
          end do ! lbl1
C
 30       if( iopt .lt. 4 ) then
 35         write(6,50) ' More j,I''s?'
            call readin( ans, *35 )
            if( ans .ne. yes ) gosu2 = .false.
          else
            gosu2 = .false.
          end if
C
cc        iend = mclock()
cc        xend = dfloat( iend )
C         write(6,310) ' Final ', iend, mclock(), xend
cc        xtime = DSCALE * ( xend - xbeg )
cc        write(6,300) 1000000, xtime, icount, xtime/dfloat(icount)
	    write(6,300) icount
cc300     format(/1x, i4, 20x, ' T o t a l :',e12.5,' seconds'/
cc   1    i15, ' calc. ',  3x, ' A v e .   :',e12.5,' seconds'/)
300		format(' Total: ', i15, ' nodes.')
cc310     format(a,'time:', 2i5, e12.5)
C
        end do ! su2
C
 40     write(6,50) ' More (lm,mu)''s?'
        call readin( ans, *40 )
        if( ans .ne. yes ) gosu3 = .false.
      end do ! su3
C
50    format(/a:a)
      continue
100   format(/a/2(2x,2i4,7x),7x,2i4,9x,' does not exist')
110   format(/a/2(2x,2i4,7x),7x,2i4,9x,' rho>9')
120   format(/a/2(2x,3i4,i3),3x,i3,1x,3i4,i3,f12.7/9(37x,i3,f28.7/))
130   format(   2(10x,i4,i3),3x,i3,9x, i4,i3,f12.7/9(37x,i3,f28.7/))
140   format(/a/2(2x,3i4,i3),7x,3i4,i3,3x,' does not exist')
C100   format(/' (lm1mu1j1tI1t;lm2mu2j1tI2t|| rhlm3mu3j1tI3t)'/
C    1       ' ',2(1x,2i3,6x),5x,2i3,9x,' does not exist')
C110   format(/' (lm1mu1j1tI1t;lm2mu2j1tI2t|| rhlm3mu3j1tI3t)'/
C    1       ' ',2(1x,2i3,6x),5x,2i3,9x,' rho>9')
C120   format(/' (lm1mu1j1tI1t;lm2mu2j1tI2t|| rhlm3mu3j1tI3t)'/
C    1       ' ',2(1x,4i3),2x,5i3,f12.7/9(29x,i3,f24.7/))
C130   format(' ',2(7x,2i3),2x,i3,6x,2i3,f12.7/9(29x,i3,f24.7/))
C140   format(/' (lm1mu1j1tI1t;lm2mu2j1tI2t|| rhlm3mu3j1tI3t)'/
C    1       ' ',2(1x,4i3),5x,4i3,3x,' does not exist')
C150  format(' ',2(7x,i3,3x),8x,i3,6x,'     epsilon equivalence')
      stop
      end
C
CB---------------------------------------------------------------------
C                        ***  INTRODUCTION ***
C ---------------------------------------------------------------------
C Update:  05/00 ... original code by C. Bahri
C ---------------------------------------------------------------------
      subroutine intro
C
      implicit logical (a-z)
C
      write(6,100)
 100  format(
     1 /' *** ------------------------------------------------ ***'/
     2  ' ***        Clebsch Gordan coefficient of SU(3)       ***'/
     3  ' ***                   in SU(2) bases                 ***'/
     4  ' *** ------------------------------------------------ ***'/
     5  ' ((lm1 mu1) j1 I1; (lm2 mu2) j2 I2 || rh (lm3 mu3) j3 I3)')
      return
      end
