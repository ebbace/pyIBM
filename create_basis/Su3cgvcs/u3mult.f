CB----------------------------------------------------------------------
C                             ****************
C                             ***  U3MULT  ***
C                             ****************
C ----------------------------------------------------------------------
C Update:  11/89 ... original code by J.P. Draayer
C          05/01 ... modified according to SU3CGVCS by C. Bahri
C ----------------------------------------------------------------------
C     Multiplicity in outer product of two SU(3) irreps
C ----------------------------------------------------------------------
      subroutine u3mult( lx1,mx1, lx2,mx2, lx3,mx3, multu3, * )               
C     ------------------------------------------------------------------
      implicit logical( a-z )
	integer lx1, mx1, lx2, mx2, lx3, mx3, multu3
	integer mx, nx, mu, nu, l1, m1, l2, m2, l3, m3
C
      multu3 = 0                                                          
      nx = lx1 + lx2 - lx3 - mx1 - mx2 + mx3                                        
      mx = nx / 3                                                           
      if( 3*mx .ne. nx ) return 1                                            
      if( mx .ge. 0) then                                                   
        l1 = lx1                                                         
        l2 = lx2                                                         
        l3 = lx3                                                         
        m1 = mx1                                                         
        m2 = mx2                                                         
        m3 = mx3                                                         
      else                                                              
        l1 = mx1                                                         
        l2 = mx2                                                         
        l3 = mx3                                                         
        m1 = lx1                                                         
        m2 = lx2                                                         
        m3 = lx3                                                         
        mx =-mx                                                         
      end if                                                             
      nx = mx + m1 + m2 - m3                                                    
      mu = min0( l1-mx, m2 )                                                 
      if( mu .lt. 0 ) return 1                                               
      nu = min0( l2-mx, m1 )                                                 
      if( nu .lt. 0 ) return 1                                               
      multu3 = max0( min0( nx, nu ) -max0( nx-mu, 0 )+1, 0 )                        
      if( multu3 .ne. 0 ) return                                             
      return 1                                                          
C ---*end of u3mult*----------------------------------------------------
      end                                                               
