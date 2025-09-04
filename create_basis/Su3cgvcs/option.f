CB---------------------------------------------------------------------
C                          ***  OPTION ***
C ---------------------------------------------------------------------
C Update:  05/00 ... original code by C. Bahri
C ---------------------------------------------------------------------
      subroutine option( iopt )
C
      implicit none
      integer iopt
C
      write(6,100)
      read (5,*)   iopt
      write(6,'(a,i4)') ' ***> option ', iopt
 100  format(/' *** Choose the following options to calculate:'/
     1        '     1) specific CG'/
     2        '     2) all 3rd weights'/
     3        '     3) all 1st and 2nd weights'/
     4        '     4) all weights')
      return
      end
