CB---------------------------------------------------------------------
C                            ***************
C                            *** GENWGHT ***
C                            ***************
C ---------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 11/98)
C ---------------------------------------------------------------------
C Update:  11/98 ... original code by C. Bahri
C ---------------------------------------------------------------------
C     Function to generate all weights in the irrep.
C ---------------------------------------------------------------------
      subroutine genwght( lm, mu, iwa, idx )
C
      implicit none
C
      integer lm, mu
      integer iwa( * )
      integer idx
C
      integer st, jt, It
C
C Intrinsic functions
      integer ipack
C
C     Packing
      ipack( lm, mu, st, jt ) = ior(jt, ishft( ior(st, ishft(
     1        ior(mu, ishft( lm, 8 )),8 )),8 ))
C
      st  = mu
      idx = 0
C
      do jt = 0, lm+mu
        do It = iabs( st-jt ), min0( st+jt, st-jt+2*lm ), 2
          idx = idx + 1
          iwa( idx ) = ipack( 0, 0, jt, It )
C         write(6,'(2x,z8,2i6)') iwa( idx ), jt, It
        end do ! It
      end do ! jt
C
      return
C ---*end of genwght*---------------------------------------------------
      end
