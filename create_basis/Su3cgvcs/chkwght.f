CB---------------------------------------------------------------------
C                            ***************
C                            *** CHKWGHT ***
C                            ***************
C ---------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 11/98)
C ---------------------------------------------------------------------
C Update:  11/98 ... original code by C. Bahri
C ---------------------------------------------------------------------
C     Function to check whether the weight is in the irrep.
C ---------------------------------------------------------------------
      logical function chkwght( lm, mu, jt, It )
C
      implicit logical( a-z )
C
      integer lm, mu
      integer st, jt, It
C
      st = mu
C
      chkwght = .false.
C
      if( jt .ge. 0 .and. jt .le. lm+mu .and. It .ge. iabs( st-jt ) 
     1    .and. It .le. min0( st+jt, st-jt+2*lm ) ) chkwght = .true.
C      
      return
C ---*end of chkwght*---------------------------------------------------
      end
