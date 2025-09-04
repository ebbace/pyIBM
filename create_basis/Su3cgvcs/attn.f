      subroutine attn( liter, * )
C     -----------------------------------------------------------------
C     Warning because of gross error.
C     -----------------------------------------------------------------
      character*(*) liter
      write(6,'(//a,a)') ' ***** ATTENTION ', liter
      return 1
c     return
C
      entry error( liter )
C     -----------------------------------------------------------------
C     Exit program because of gross error.
C     -----------------------------------------------------------------
      write(6,'(//a,a)') ' ***** FATAL ERROR ', liter
      stop
c ---*end of attn*------------------------------------------------------
      end
