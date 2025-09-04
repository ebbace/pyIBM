      program cgu3o3drv
CB---------------------------------------------------------------------
C                        **********************
C                        ***  CGU3O3 DRIVER ***
C                        **********************
C ---------------------------------------------------------------------
C Update:  11/98 ... original code by C. Bahri
C          05/02 ... look up tables for I- and S-functions (CB)
C ---------------------------------------------------------------------
C     Driver to calculate the Clebsch-Gordan coefficients of SU(3)
C     in SO(3) bases.
C ---------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C ---------------------------------------------------------------------
C
      implicit none
      integer NMAX, NRO, NLBLMAX
      character*57 HEADER
      parameter ( NMAX = 1000, NRO = 10, NLBLMAX = 10000 )
      parameter ( HEADER = 
     1  ' ((lm1 mu1)al1 L1; (lm2 mu2)al2 L2 || rh (lm3 mu3)al3 L3)' )
C                 1         2         3         4         5
C        123456789012345678901234567890123456789012345678901234567
C
C Arrays
      integer          ifu3e( NMAX )                 ! packed labels
      double precision dcgu3e( NMAX * NRO )          ! seed coefficients
      double precision dcgu3o3( NRO, NRO, NRO, NRO ) ! CG coefficients
      double precision dtemp
C
C Indices and temporary variables
      logical gosu3, goso3, gorite
      integer iopt,iovalue                              ! option number
      integer icount                            ! index and counter
      integer kro, kp1, kp2, kp3
      character*2 ans, no, yes
      data yes/'Y'/, no/'N'/
C
C Irreps and weights
      integer lm1,mu1, lm2,mu2, lm3,mu3, kromax ! irreps
      integer L1, L2, L3
      integer nhw                               ! # of seed coefficients
C
C Arrays of weights
      integer iw1a( NLBLMAX ), iw2a( NLBLMAX ), iw3a( NLBLMAX )
      integer lbl1, lbl2, lbl3
      integer lbl1max, lbl2max, lbl3max
C
C Clock
C
c     integer ibeg, iend
c     integer mclock
c     double precision xbeg, xend, xtime, dfloat, DSCALE
c     parameter( DSCALE = 1.d-6 )
C
C Extrinsic functions
      integer iabs, max0 !, min0, mod
C
C Intrinsic functions
      integer kmult !, kstart
C
C     K-multiplicity
C
      kmult( lm3, mu3, L3 ) = max0( 0, (lm3 + mu3 + 2 - L3)/2 )
     1    - max0( 0, (lm3 + 1 - L3)/2 ) - max0( 0, (mu3 + 1 - L3)/2 )
C
CK    kstart( lm3, mu3, L3 ) = mod( lm3, 2 ) + 2*( max0(0, (L3 - mu3)/2)
CK   1    + mod( lm3+1, 2 ) * mod( iabs(L3 - mu3), 2 ) )
C
C     Initializations
C
C     Initializations
      call readfact				! load factorials
      call readtab                              ! load I- and S-functions
      call intro                                ! introductory page
cq      call option( iopt )                       ! choose option
      iopt=1
C
C     Proceed ...
      gosu3 = .true.
C
        write(6,50) ' Enter: lm1,mu1, lm2,mu2, lm3,mu3::'
        read (5,*)   lm1,mu1, lm2,mu2, lm3,mu3
        call u3mult( lm1,mu1, lm2,mu2, lm3,mu3, kromax, *5 )
5       if( kromax .eq. 0 ) then
          write(6,100) HEADER, lm1,mu1, lm2,mu2, lm3,mu3
          go to 40
        end if ! kromax
C
C       write(6,*) ' Entering CG-HW'
C
C ----*------*------*------*------*------*------*------*------*------*-
        call cgu3hw( lm1,mu1, lm2,mu2, lm3,mu3, nhw,
     1                   ifu3e, dcgu3e, 1000 )
C ----*------*------*------*------*------*------*------*------*------*-
C
        goso3 = .true.
C        
        gorite=.true.
        do while ( goso3 )
C
c         ibeg = mclock()
c         xbeg = dfloat( ibeg )
C         write(6,210) ' Initial ', ibeg, mclock(), xbeg
C
C Options
C ---------------------------------------------------------------------
          if( iopt .eq. 1 ) then
C ---------------------------------------------------------------------
C
c            write(6,50) ' Enter: L1, L2, L3::'
           read(5,*,iostat=iovalue)   L1, L2, L3
           if (iovalue<0)goto 40
C
            if(      kmult( lm1, mu1, L1 ) .eq. 0
     2          .or. kmult( lm2, mu2, L2 ) .eq. 0
     3          .or. kmult( lm3, mu3, L3 ) .eq. 0
     4          .or. L3 .lt. iabs(L1-L2) .or. L3 .gt. L1+L2 ) then
              write(6,140) HEADER, lm1, mu1, 0, L1,
     1                             lm2, mu2, 0, L2,
     2                             lm3, mu3, 0, L3
              go to 30
            end if
C
            lbl1max = 1
            lbl2max = 1
            lbl3max = 1
C
            iw1a( 1 ) = L1
            iw2a( 1 ) = L2
            iw3a( 1 ) = L3
C
C ---------------------------------------------------------------------
          else if( iopt .eq. 2 ) then
C ---------------------------------------------------------------------
C
            write(6,50) ' Enter: L1, L2::'
            read(5,*)   L1, L2
C
            if(      kmult( lm1, mu1, L1 ) .eq. 0
     1          .or. kmult( lm2, mu2, L2 ) .eq. 0 ) then
              write(6,140) HEADER, lm1, mu1, 0, L1,
     1                             lm2, mu2, 0, L2,
     2                             lm3, mu3, NMAX, NMAX
              go to 30
            end if
C
            lbl1max = 1
            lbl2max = 1
C
            iw1a( 1 ) = L1
            iw2a( 1 ) = L2
C
            call genL( lm3, mu3, iw3a, lbl3max )
C
C ---------------------------------------------------------------------
          else if( iopt .eq. 3 ) then
C ---------------------------------------------------------------------
C
            write(6,50) ' Enter: L3::'
            read(5,*)   L3
C
            if( kmult( lm3, mu3, L3 ) .eq. 0 ) then
              write(6,140) HEADER, lm1, mu1, NMAX, NMAX,
     1                             lm2, mu2, NMAX, NMAX,
     2                             lm3, mu3, 0, L3
              go to 30
            end if
C
            lbl3max = 1
            iw3a( 1 ) = L3
C
            call genL( lm1, mu1, iw1a, lbl1max )
            call genL( lm2, mu2, iw2a, lbl2max )
C
C ---------------------------------------------------------------------
          else if( iopt .eq. 4 ) then
C ---------------------------------------------------------------------
C
            call genL( lm1, mu1, iw1a, lbl1max )
            call genL( lm2, mu2, iw2a, lbl2max )
            call genL( lm3, mu3, iw3a, lbl3max )
C
          end if
          if( lbl1max .gt. NLBLMAX .or. lbl2max .gt. NLBLMAX 
     1                             .or. lbl3max .gt. NLBLMAX ) 
     2      call error(' CGU3O3DRV: Increase NLBLMAX.')
C
C         Calculating CGs ...
C
          icount = 0
          do lbl1 = 1, lbl1max
            L1 = iw1a( lbl1 )
          do lbl2 = 1, lbl2max
            L2 = iw2a( lbl2 )
          do lbl3 = 1, lbl3max
            L3 = iw3a( lbl3 )
C
          if( L3 .ge. iabs(L1-L2) .and. L3 .le. L1+L2 ) then
C
C ----*------*------*------*------*------*------*------*------*------*-
            call cgu3o3( lm1,mu1, lm2,mu2, lm3,mu3, L1,L2,L3,
     1                   nhw, ifu3e, dcgu3e, dcgu3o3, NMAX )
C ----*------*------*------*------*------*------*------*------*------*-
C
            dtemp = 0.d0
              do kp1 = 1, kmult( lm1, mu1, L1 )
              do kp2 = 1, kmult( lm2, mu2, L2 )
              do kp3 = 1, kmult( lm3, mu3, L3 )
C
            icount = icount + 1
		do kro = 1, kromax
	      dtemp  = dtemp + dcgu3o3( kp1, kp2, kp3, kro )**2
		end do
C
            if( gorite ) then
                write(6,120) HEADER, lm1, mu1, 1, L1,
     1                       lm2, mu2, 1, L2,
     2                       1, lm3, mu3, 1, L3, dcgu3o3(1,1,1,1),
     3                       ( kro, dcgu3o3(1,1,1,kro), kro=2, kromax )



C               if( logfile.eq.yes ) write(ifl,120) HEADER,
C    1                       lm1, mu1, 1, L1,  lm2, mu2, 1, L2,
C    2                       1, lm3, mu3, 1, L3, dcgu3o3(1,1,1,1),
C    3                       ( kro, dcgu3o3(1,1,1,kro), kro=2, kromax )
              gorite = .false.
            else
C               if( kp1 .ne. 1 .and. kp2 .ne. 1 .and. kp3 .ne. 1 ) then
                write(6,130) kp1,L1, kp2,L2, 1, kp3,L3,
     1                   dcgu3o3( kp1, kp2, kp3, 1 ), 
     2                   ( kro, dcgu3o3(kp1,kp2,kp3,kro), kro=2,kromax )
C               if( logfile.eq.yes ) write(ifl,130) kp1, L1, kp2, L2,
C    1                   1, kp3, L3, dcgu3o3( kp1, kp2, kp3, 1 ), 
C    2                   ( kro, dcgu3o3(kp1,kp2,kp3,kro), kro=2,kromax )
C               end if
            end if ! writing
C
              end do
              end do
              end do
C
            if( dtemp .gt. 1.d0 ) write( 6, * ) 'WARNING'
          end if ! allowed coupling
C
          end do ! lbl3
          end do ! lbl2
          end do ! lbl1
C
c         iend = mclock()
c         xend = dfloat( iend )
c         write(6,'(1pq32.26)') dcgu3o3(1,1,1,1)
C         write(6,210) ' Final ', iend, mclock(), xend
c         xtime = DSCALE * ( xend - xbeg )
c         write(6,200) 1000000, xtime, icount, xtime/dfloat(icount)
c		write(6,200) icount
c200      format(/1x, i4, 20x, ' T o t a l :',e12.5,' seconds'/
c    1    i15, ' calc. ',  3x, ' A v e .   :',e12.5,' seconds'/)
 200		format(' Total: ', i15, ' nodes.')
c210      format(a,'time:', 2i5, e12.5)
C
c 30       if( iopt .lt. 4 ) then
c 35         write(6,50) ' More L''s?'
c            call readin( ans, *35 )
c            if( ans .ne. yes ) goso3 = .false.
c          else
c            goso3 = .false.
c          end if
C
30      continue
        end do ! so3
C
c40      write(6,50) ' More (lm,mu) s?'
c        call readin( ans, *40 )
c        if( ans .ne. yes ) gosu3 = .false.
c      end do ! su3
C
40    continue
50    format(/a:a)
      continue
100   format(/a/2(2x,2i4,7x),7x,2i4,9x,' does not exist')
110   format(/a/2(2x,2i4,7x),7x,2i4,9x,' rho>9')
120   format(/a/2(2x,3i4,i3),3x,i3,1x,3i4,i3,f12.7/9(37x,i3,f28.7/))
130   format(   2(10x,i4,i3),3x,i3,9x, i4,i3,f12.7/9(37x,i3,f28.7/))
140   format(/a/2(2x,3i4,i3),7x,3i4,i3,3x,' does not exist')
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
     3  ' ***                   in SO(3) bases                 ***'/
     4  ' *** ------------------------------------------------ ***'/
     5  ' ((lm1 mu1)al1 L1; (lm2 mu2)al2 L2 || rh (lm3 mu3)al3 L3)')
      return
      end
