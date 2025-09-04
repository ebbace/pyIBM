      subroutine readin(ans,*)                                                  
C     -----------------------------------------------------------------
C     Subprogram to limit response to Yes or No answer.
C     -----------------------------------------------------------------
      parameter(loopmax=25)                                                     
      character*2 ans,yes,yea,no,nop                                            
      logical ok                                                                
      save loop                                                                 
      data yes/'y'/,yea/'y'/,no/'n'/,nop/'n'/                                   
C                                                                               
      ok=.false.                                                                
      read(5,'(a)',end=10) ans                                                  
      if(ans .eq. yea) then                                                     
         ans = yes                                                              
         ok=.true.                                                              
      else if (ans .eq. nop) then                                               
         ans = no                                                               
         ok=.true.                                                              
      else if (ans .eq. yes) then                                               
         ok=.true.                                                              
      else if (ans .eq. no) then                                                
         ok=.true.                                                              
      end if                                                                    
10    if(.not. ok) then                                                         
        write(6,'(a)') ' **Use Y (y) for yes and N (n) for no!'                 
        loop=loop + 1                                                           
        if(loop .lt. loopmax) then                                              
           return 1                                                              
        else                                                                    
           stop '@READIN: enough attempts!'                                     
        end if                                                                  
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
