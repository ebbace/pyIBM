CB---------------------------------------------------------------------
C                             **************
C                             *** ORTHON ***
C                             **************
C ---------------------------------------------------------------------
C Author: Chairul Bahri (University of Toronto, 11/97)
C ---------------------------------------------------------------------
C Update:  11/97 ... original code by C. Bahri
C ---------------------------------------------------------------------
C     Function to orthonormalize vectors.
C ---------------------------------------------------------------------
      subroutine orthon( nsize, nvec, nmax, dvec )
      implicit none
      integer          nsize                    ! dimension of vectors
      integer          nvec                     ! number of vectors
      integer          nmax                     ! physical dimension
      double precision dvec( nmax, * )          ! the vectors
      integer          ivec, jvec, isize
      double precision dsqrt                    ! intrinsic function
      double precision sum
      double precision ZERO
      parameter ( ZERO = 1.d-12 )
C
C     Gramm-Schmidt orthonormalization
C
      do ivec = 1, nvec
        do jvec = 1, ivec - 1
C         -------------------------------------------------------
          sum = 0.d0
          do isize = 1, nsize
            sum = sum + dvec( isize, ivec ) * dvec( isize, jvec )
          end do ! isize
C         -------------------------------------------------------
C         orthogonalization
C         -----------------
          do isize = 1, nsize
            dvec( isize, ivec ) = dvec( isize, ivec )
     1                           - sum * dvec( isize, jvec )
          end do ! isize
        end do ! jvec
C         -------------------------------------------------------
          sum = 0.d0
          do isize = 1, nsize
            sum = sum + dvec( isize, ivec ) * dvec( isize, ivec )
          end do ! isize
C         -------------------------------------------------------
C         normalization
C         -------------
        if ( sum .gt. ZERO ) then
          sum = 1.d0 / dsqrt( sum )  ! could be a phase problem
          do isize = 1, nsize
            dvec( isize, ivec ) = dvec( isize, ivec ) * sum
          end do ! isize
        end if ! sum
      end do ! ivec
      return
C ---*end of orthon*----------------------------------------------------
      end
