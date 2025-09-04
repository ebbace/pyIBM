CB----------------------------------------------------------------------
C                             ****************    
C                             *** d matrix ***
C                             ****************    
C ----------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 03/99)
C ----------------------------------------------------------------------
C Update:  03/99 ... original code by C. Bahri
C ----------------------------------------------------------------------
C     Function to calculate the Wigner d-matrix.
C ----------------------------------------------------------------------
C References:
C  1. C. Bahri, D.J. Rowe, and J.P. Draayer, CPC (in press)
C  2. D.J. Rowe and C. Bahri, JMP 41, 6544-6565 (2000)
C ----------------------------------------------------------------------
      function dmat( Mt, Jt, Nt, mut )
      implicit none
      integer Mt, Jt, Nt, mut
      double precision dmat
      double precision dlogf
      double precision dexp              ! intrinsic functions
      logical          btest
C
      common / BKDF / dlogf(0:2000)
C
C     d(M,J,N,mu) = (-)^(M-N+mu) sqrt((J+M)!(J-M)!(J+N)!(J-N)!)/
C                   (J+N-mu)!mu!(J-M-mu)!(M-N+mu)!
C
      dmat = 0.5d0 * ( dlogf( Jt + Mt ) + dlogf( Jt - Mt )
     1               + dlogf( Jt + Nt ) + dlogf( Jt - Nt ) )
     2       - ( dlogf( Jt + Nt - mut ) + dlogf( mut ) 
     3         + dlogf( Jt - Mt - mut ) + dlogf( Mt - Nt + mut ) )
      if( btest( Mt - Nt + mut, 1 ) ) then
        dmat = - dexp( dmat )
      else
        dmat =   dexp( dmat )
      end if 
C
      return
C ---*end of dmat*------------------------------------------------------
      end
