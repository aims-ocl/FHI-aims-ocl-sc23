      subroutine invert_simple_complex (matr , inv_matr ,n_matrix)

! PURPOSE  
! 
! Invert the Green Function  G(w)_ij -> [G(w)^-1]_ij  at each frequency points through
! LU decomposition (this method could not be the fastest but it is applicable to
! almost every square matrix)
!

     use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit

      implicit none
      integer n_matrix
      complex*16  matr (n_matrix,n_matrix)
      complex*16  inv_matr (n_matrix,n_matrix)
!      complex*16  inv_green_fn_freq(n_matrix,n_matrix,nomega)
     
! Counters
      integer i_freq 
      integer i_basis

!auxiliary stuff
!      complex*16  product_of_the_two(n_matrix,n_matrix,nomega)

! quantities to compute the matrix inversion     
      integer info
      integer :: ipiv(n_matrix)
      complex*16 :: work (n_matrix)

! first calculate the LU factorization

! to avoid that the green function is overwritten by computing it's inverse 
      inv_matr(:,:) = matr(:,:)

! at each frequency point evaluate the LU factorization, and check for errors
!      do i_freq = 1, nomega, 1
        call zgetrf( n_matrix, n_matrix, inv_matr, &
                     n_matrix, ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
            stop
          endif
        endif
! this invert the matrix in the LU factorized form
        call zgetri(n_matrix, inv_matr, n_matrix, &
                 ipiv, work, n_matrix, info)

       if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
            stop
          endif
        endif


!      enddo

! test if  green_fn_freq*inv_green_fn_freq = 1

!      product_of_the_two(:,:,:) = 0.d0
!
!      do i_freq = 1, nomega, 1 
!        call zgemm ('N','N',n_matrix, n_matrix, n_matrix, &
!           1.d0, green_fn_freq(1,1,i_freq), n_matrix, &
!          inv_green_fn_freq(1,1,i_freq), n_matrix, 0.d0,&
!         product_of_the_two(1,1,i_freq), n_matrix  )
!      enddo

!      do i_freq=1, nomega, 5
!        write(use_unit,*) " "
!        write(use_unit,*) i_freq
!        write(use_unit,*) " "
!        do i_basis = 1, n_matrix, 1
!          write(use_unit,*) product_of_the_two(i_basis,i_basis,1), product_of_the_two(i_basis,1,1)
!        enddo
!      enddo
!      if(myid.eq.0)print* , product_of_the_two(:,:,:)
!
!      if(myid.eq.0)then
!        write(use_unit,*) "  --- Inversion of the Green's Function: DONE ---"
!        write(use_unit,*) " "
!      endif

      end subroutine invert_simple_complex

!-----------------------------------------------------------------
      subroutine invert_green (green_fn_freq,nomega, inv_green_fn_freq,n_matrix)

! PURPOSE  
! 
! Invert the Green Function  G(w)_ij -> [G(w)^-1]_ij  at each frequency points through
! LU decomposition (this method could not be the fastest but it is applicable to
! almost every square matrix)
!

     use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit

      implicit none

      integer n_matrix
      integer nomega
      complex*16  green_fn_freq(n_matrix,n_matrix,nomega)
      complex*16  inv_green_fn_freq(n_matrix,n_matrix,nomega)
     
! Counters
      integer i_freq 
      integer i_basis

!auxiliary stuff
!      complex*16  product_of_the_two(n_matrix,n_matrix,nomega)

! quantities to compute the matrix inversion     
      integer info
      integer :: ipiv(n_matrix)
      complex*16 :: work (n_matrix)

! first calculate the LU factorization

! to avoid that the green function is overwritten by computing it's inverse 
      inv_green_fn_freq(:,:,:) = green_fn_freq(:,:,:)

! at each frequency point evaluate the LU factorization, and check for errors
      do i_freq = 1, nomega, 1
        call zgetrf( n_matrix, n_matrix, inv_green_fn_freq(1,1,i_freq), &
                     n_matrix, ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
            stop
          endif
        endif
! this invert the matrix in the LU factorized form
        call zgetri(n_matrix, inv_green_fn_freq(1,1,i_freq), n_matrix, &
                 ipiv, work, n_matrix, info)

       if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
            stop
          endif
        endif


      enddo

! test if  green_fn_freq*inv_green_fn_freq = 1

!      product_of_the_two(:,:,:) = 0.d0
!
!      do i_freq = 1, nomega, 1 
!        call zgemm ('N','N',n_matrix, n_matrix, n_matrix, &
!           1.d0, green_fn_freq(1,1,i_freq), n_matrix, &
!          inv_green_fn_freq(1,1,i_freq), n_matrix, 0.d0,&
!         product_of_the_two(1,1,i_freq), n_matrix  )
!      enddo

!      do i_freq=1, nomega, 5
!        write(use_unit,*) " "
!        write(use_unit,*) i_freq
!        write(use_unit,*) " "
!        do i_basis = 1, n_matrix, 1
!          write(use_unit,*) product_of_the_two(i_basis,i_basis,1), product_of_the_two(i_basis,1,1)
!        enddo
!      enddo
!      if(myid.eq.0)print* , product_of_the_two(:,:,:)
!
!      if(myid.eq.0)then
!        write(use_unit,*) "  --- Inversion of the Green's Function: DONE ---"
!        write(use_unit,*) " "
!      endif

      end subroutine invert_green
