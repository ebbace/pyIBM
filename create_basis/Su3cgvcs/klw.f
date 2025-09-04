CB----------------------------------------------------------------------
C                              *************      
C                              *** K l w ***
C                              *************      
C ----------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 11/97)
C ----------------------------------------------------------------------
C Update:  11/97 ... original code by C. Bahri
C ----------------------------------------------------------------------
C     Function to calculate the K-matrix of SU(2)xU(1).
C ----------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C  3. D.J. Rowe and J. Repka, JMP 38, 4363-4388 (1997)
C ----------------------------------------------------------------------
      function klw( lm, mu, jt, It )
      implicit none
      integer lm, mu, jt, It, st         ! SU(3) > SU(2)xU(1) irrep
      double precision klw               ! .t = 2 .
      double precision dlogf
      double precision dtop, dbot        ! temporary variables
      double precision dexp              ! intrinsic functions
C
      common / BKDF / dlogf(0:2000)
C
C     K(lm,om) = dsqrt((lm1+1)! lm2! / (om1+1)! om2! )
C
C Input/Output Variables
C     lm*    > 
C
      st   = mu
      dtop = dlogf(2*(lm+st+1)) + dlogf(2*lm) 
      dbot = dlogf(2*lm+st+It-jt+2) + dlogf(2*lm+st-It-jt)
      klw  = dexp(0.5d0 * (dtop - dbot))
C
      return
C ---*end of klw*-------------------------------------------------------
      end
