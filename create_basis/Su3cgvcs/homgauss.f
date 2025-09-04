CB---------------------------------------------------------------------
C                            ****************
C                            *** HOMGAUSS ***
C                            ****************
C ---------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 06/00)
C ---------------------------------------------------------------------
C Update:  06/00 ... original code by C. Bahri
C ---------------------------------------------------------------------
C     Function to obtain the solution of the set of homogeneous linear
C     equation using Gauss elimination/row reduction.
C ---------------------------------------------------------------------
      subroutine homgauss( nrow, ncol, nmax, dmtx, dsol )
      implicit none
      integer          nrow                     ! number of equations
      integer          ncol                     ! number of variables
      integer          nmax                     ! physical dimension
      double precision dmtx( nmax, * )          ! the matrix
      double precision dsol( nmax, * )          ! the solution
      integer          nmult                    ! number of multiplicity
      integer          irow, jcol, ipick, kmult !, jrun
      double precision dabs                     ! intrinsic function
      double precision big, sum, temp
      double precision ZERO !, TOL
      parameter ( ZERO = 1.d-5 ) !, TOL = 1.d-5 )
C
C     The largest elements for the main diagonal of the matrix.
C     -> pivot
C
      do irow = 1, nrow
C
C       pivot
C
        big = dabs( dmtx( irow, irow ) )
        do ipick = irow + 1, nrow
          if( dabs( dmtx( ipick, irow ) ) .gt. big ) then
            big = dabs( dmtx( ipick, irow ) )
C
C       row interchange
C
            do jcol = irow, ncol
              temp = dmtx( irow, jcol )
              dmtx( irow, jcol )  = dmtx( ipick, jcol )
              dmtx( ipick, jcol ) = temp
            end do ! jcol
          end if ! big
        end do ! ipick
C
C       row normalization
C
          if( dabs( dmtx( irow, irow ) ) .gt. ZERO ) then
            temp = 1.d0/dmtx( irow, irow )
            do jcol = irow+1, ncol
              dmtx( irow, jcol ) = dmtx( irow, jcol ) * temp
            end do ! jcol
              dmtx( irow, irow ) = 1.d0
          end if ! irow
C
C       row reduction
C
        do ipick = irow+1, nrow
          temp = dmtx( ipick, irow )
C	  write(6,*) irow, ipick, nrow, temp
          do jcol = irow+1, ncol
            dmtx( ipick, jcol ) = dmtx( ipick, jcol )
     1                          - temp * dmtx( irow, jcol )
          end do ! jcol
            dmtx( ipick, irow ) = 0.d0
        end do ! ipick
C
      end do ! irow
CW    do irow = 1, nrow
C	write(6,'(4d16.4)') (dmtx(irow,jcol),jcol=1,ncol)
CW    end do
C
C     rank checking
C
      do irow = nrow, 1, -1
C       sum = 0.d0
C       do jcol = 1, ncol
C         sum = sum + dabs( dmtx( irow, jcol ) )
C       end do ! jcol
C       sum = dabs( dmtx( irow, irow ) )
        if( dabs( dmtx(irow,irow) ) .gt. ZERO ) then
CW        write(6,*) ' In row', irow, dmtx( irow, irow )
	  goto 100
        end if ! ZERO
      end do ! irow
 100  nrow = irow
CW    write(6,1000) nrow, ncol
C1000 format(' There are ', i4, ' linearly independent equations'
CW   1      /' with ', 5x, i4, ' variables.')
      if( nrow .ge. ncol ) then
        if( nrow .eq. ncol ) write(6,1010) ' has trivial solutions.'
        if( nrow .gt. ncol ) write(6,1010) ' is overdetermined.'
        do irow = 1, nrow
        do jcol = 1, ncol
          dsol( irow, jcol ) = 0.d0
        end do
        end do
        return
C     else
C       if( nrow .lt. ncol ) write(6,1010) ' is underdetermined.'
      end if
 1010 format(' The system of equations', a)
C
C     back substitution ... different from regular Gauss elimination
C
      nmult = ncol - nrow
      do irow = nrow, 1, -1
C
C       independent solutions
C
        do jcol = 1, nmult
          do kmult = 1, nmult
            if( kmult .eq. jcol ) then
              dsol( nrow+jcol, kmult ) = 1.d0
            else
              dsol( nrow+jcol, kmult ) = 0.d0
            end if ! kmult
          end do ! kmult
        end do ! jcol
C
C       dependent solutions
C
C       temp = 1.d0 / dmtx( irow, irow )
          do kmult = 1, nmult
C         -------------------------------------------------------
            sum = 0.d0
            do jcol = irow+1, ncol
              sum = sum + dmtx( irow, jcol ) * dsol( jcol, kmult )
            end do ! jcol
C           dsol( irow, kmult ) = - sum * temp
            dsol( irow, kmult ) = - sum !* temp
C         -------------------------------------------------------
          end do ! kmult
          
*DEB    do jcol = nrow, irow, -1
C         temp = 1.d0 / dmtx( irow, jcol )
*         do kmult = 1, nmult
C         -------------------------------------------------------
*           sum = 0.d0
*           do jrun = jcol+1, ncol
*             sum = sum + dmtx( irow, jrun ) * dsol( jrun, kmult )
*           end do ! jrun
*           dsol( jcol, kmult ) = - sum * temp
C         -------------------------------------------------------
*         end do ! kmult
*DEB    end do ! jcol
      end do ! irow
      return
C ---*end of homgauss*--------------------------------------------------
      end
