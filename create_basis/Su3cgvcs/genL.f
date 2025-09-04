CB---------------------------------------------------------------------
C                            ***************
C                            ***   GENL  ***
C                            ***************
C ---------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 11/98)
C ---------------------------------------------------------------------
C Update:  11/98 ... original code by C. Bahri
C ---------------------------------------------------------------------
C     Function to generate all Ls in the irrep.
C ---------------------------------------------------------------------
      subroutine genL( lm, mu, iwa, idx )
C
      implicit none
C
      integer lm, mu
      integer iwa( * )
      integer idx
C
      integer L
C
C Intrinsic functions
      integer kmult
C
C     K-multiplicity
C
      kmult( lm, mu, L ) = max0( 0, (lm + mu + 2 - L)/2 )
     1    - max0( 0, (lm + 1 - L)/2 ) - max0( 0, (mu + 1 - L)/2 )
C
      idx = 0
C
      do L = 0, lm+mu
        if( kmult( lm, mu, L ) .ne. 0 ) then
          idx = idx + 1
          iwa( idx ) = L
        end if
      end do ! jt
C
      return
C ---*end of genwght*---------------------------------------------------
      end
