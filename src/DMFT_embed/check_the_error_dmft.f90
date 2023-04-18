      subroutine check_the_error_dmft (polar_freq, polar_green_freq, &
      !    n_matrix, omega, nomega, womega,&
                 name_of_quantity,&
                 average_error, max_error, printout )


      use dimensions
      use runtime_choices
      use species_data
      use pbc_lists
      use physics
      use prodbas
      use scgw_grid
      use constants
      use mpi_tasks
!      use synchronize_mpi
      use hartree_fock
      use localized_basbas
      use gw_para
      use dmft_para
      use poles_fit
      use gt
      use timing

      !n_matrix = n_basis
!ARGUMENTS
!      integer n_matrix, nomega,i_basis_1
      !real*8 omega(-nomega:nomega)
!      real*8 omega(nomega)
!      real*8 womega(nomega)
      !real*8 polar_freq       (n_basbas, n_loc_prodbas,1:nomega)
      !complex*16 polar_freq (n_matrix, n_matrix,nomega)
      complex*16, dimension(n_basis, n_basis,nomega) :: polar_freq 
      !real*8 polar_green_freq (n_basbas, n_loc_prodbas,1:nomega)
      !complex*16 polar_green_freq (n_matrix, n_matrix,nomega)
      complex*16, dimension(n_basis, n_basis,nomega) ::  polar_green_freq
      logical printout

! statistical stuff      
      !real*8 error (n_matrix, n_matrix)
      real*8, allocatable:: error (:, :)
      real*8, allocatable:: error_RE (:, :,:)
      real*8, allocatable:: error_IM (:, :,:)

      real*8, allocatable:: average_polar_matr (:, :)
!      real*8 error_row(n_matrix)
      real*8 average_polar
      real*8 average_error
!      real*8 average_relative_error
      real*8 max_error
!      real*8 hystogram_interval
!      integer n_hystogram_columns
!      real*8, dimension(:,:), allocatable :: hystogram

!      real*8 reference_value (n_matrix, n_matrix)
 
!counters
      integer i_freq
      integer i_basis
      integer j_basis
!      integer i_interval
!     integer counter
      character*30 name_of_quantity
!      character*40 file_name
!      character*4 extension     
!      logical end_of_hystogram_loop          
      real*8 threshold 

!     if(myid.eq.0) write(use_unit,*) , nomega, n_basis

      !if(.not. allocated(error)) allocate(error (n_matrix, n_matrix))
      !if(.not. allocated(error_RE)) allocate(error_RE (n_matrix, n_matrix,nomega))
      !if(.not. allocated(error_IM)) allocate(error_IM (n_matrix, n_matrix,nomega))
      !if(.not. allocated(average_polar_matr)) allocate(average_polar_matr (n_matrix, n_matrix))
      if(.not. allocated(error)) allocate(error (n_basis, n_basis))
      if(.not. allocated(error_RE)) allocate(error_RE (n_basis, n_basis,nomega))
      if(.not. allocated(error_IM)) allocate(error_IM (n_basis, n_basis,nomega))
      if(.not. allocated(average_polar_matr)) allocate(average_polar_matr (n_basis, n_basis))
!      do i_freq = 1, nomega/2, 1
!      do i_basis_1 = 1, n_matrix, 1
!            write(use_unit,*)'embed_gf from check_the_norm!!!!!!!!!',  i_basis_1, &
!             polar_freq(i_basis_1, i_basis_1,i_freq), polar_green_freq(i_basis_1, i_basis_1,i_freq)
!      enddo
!      enddo
!stop


      error(:,:) = 0.d0
      error_RE(:,:,:) = 0.d0
      error_IM(:,:,:) = 0.d0
      max_error = 0.d0
      !reference_value (:,:) = 0.d0
      average_polar = 0.d0
      average_polar_matr (:,:)= 0.d0
      do i_freq=1,nomega, 1
        !do i_basis = 1, n_matrix,1 
        !  do j_basis = 1, n_matrix,1
        do i_basis = 1, n_basis,1 
          do j_basis = 1, n_basis,1

!            write(use_unit,*) i_freq
!            error(i_basbas,j_basbas) = error(i_basbas,j_basbas)+(polar_freq(i_basbas,j_basbas,i_freq)-&
!            polar_green_freq (i_basbas,j_basbas,i_freq))**2 *womega(i_freq)

           error_RE(i_basis,j_basis,i_freq) =(real(polar_freq(i_basis,j_basis,i_freq)-&
           polar_green_freq (i_basis,j_basis,i_freq)))**2 
           error_IM(i_basis,j_basis,i_freq) =(aimag(polar_freq(i_basis,j_basis,i_freq)-&
           polar_green_freq (i_basis,j_basis,i_freq)))**2 
!           error(:,:) = error(:,:)+ (error_RE(:,:)+&
!           error_IM(:,:)) *womega(i_freq)
          
!           average_polar_matr (:,:) =  average_polar_matr (:,:)+&
!           polar_freq (:,:, i_freq)*womega(i_freq) 
!           reference_value (:,:)= reference_value (:,:)+ &
!           polar_freq(:,:,i_freq)**2 *womega(i_freq)

          enddo
        enddo
      enddo
!      do i_basis_1 = 1, n_matrix, 1
!            write(use_unit,*)'error_RE from check_the_norm!!!!!!!!!',  i_basis_1, &
!             error_RE(i_basis_1, i_basis_1,100), error_IM(i_basis_1, i_basis_1,100)
!      enddo

      do i_freq=1,nomega, 1
      !do i_basis = 1, n_matrix,1
      !  do j_basis = 1, n_matrix,1
      do i_basis = 1, n_basis,1
        do j_basis = 1, n_basis,1

           error(i_basis,j_basis) = error(i_basis,j_basis)+ (error_RE(i_basis,j_basis,i_freq)+&
           error_IM(i_basis,j_basis,i_freq)) *womega(i_freq)

           average_polar_matr (i_basis,j_basis) =  average_polar_matr (i_basis,j_basis)+&
           polar_freq (i_basis,j_basis, i_freq)*womega(i_freq)
         enddo
       enddo
      enddo

!      do i_basis_1 = 1, n_matrix, n1
!            write(use_unit,*)'error from check_the_norm!!!!!!!!!',  i_basis_1, &
!             error(i_basis_1, i_basis_1)
!      enddo


!  some statistics on the error

! average error----------------------------------------------------------
      average_error = 0.d0
      !naverage_relative_error = 0.d0
      threshold = 1.d-6
      !do i_basis = 1, n_matrix,1 
      !  do j_basis = 1, n_matrix,1
      do i_basis = 1, n_basis,1 
        do j_basis = 1, n_basis,1

          if (printout)then
!           if (error(i_basbas,j_basbas).gt.threshold) then
!            if(myid.eq.0)then

!             write(use_unit,*) " *** matrix element ", i_basbas,j_basbas
!
!             if(i_basbas.lt.10)then
!                if(j_basbas .lt.10)then 
!                   write(extension,'(A,I1,A,I1)') "0",i_basbas,"0",j_basbas
!                else 
!                   write(extension,'(A,I1,I2)') "0",i_basbas,j_basbas
!                endif
!             else
!                if(j_basbas .lt.10)then
!                   write(extension,'(I2,A,I1)') i_basbas,"0",j_basbas
!                else 
!                   write(extension,'(I2,I2)') i_basbas,j_basbas
!                endif
!             endif

!             if(.false.)then
!                file_name = 'green'//extension//'.dat'
!                open(44,file=file_name)
!                  do i_freq=1,nomega, 1
!                    write(44,*) omega(i_freq),&
!                polar_green_freq (i_basbas,j_basbas,i_freq),&
!                polar_freq (i_basbas,j_basbas,i_freq) 
!                  enddo 
!                close(44)
!             endif
! 
!            endif
!           endif
          endif

          average_error = average_error + &
           error(i_basis,j_basis)

          average_polar = average_polar +&
             average_polar_matr(i_basis,j_basis)

         

          if (error(i_basis,j_basis).gt.max_error)then
             max_error = error(i_basis,j_basis)
          endif

        enddo
      enddo    
 
      max_error = sqrt(max_error)
      !average_polar = average_polar / (n_matrix*n_matrix)
      average_polar = average_polar / (n_basis*n_basis)
      !average_error = sqrt(average_error /(n_matrix*n_matrix))
      average_error = sqrt(average_error /(n_basis*n_basis))
!      average_relative_error = average_error / average_polar

      if(myid.eq.0)then
        write(use_unit,'(A,A,E14.3)') "          | Average Deviation on  ",&
!        write(use_unit,*) "          | Average Deviation on  ",&
                   name_of_quantity, average_error
!        write(use_unit,'(A,A,E14.3)') "          | Relative Deviation on  ",&
!                   name_of_quantity , average_relative_error
        !write(use_unit,'(A,A,E14.3)') "          | Maximum Deviation on  ",&
        !        name_of_quantity , max_error
      endif
!------------------------------------------------------------------------
      if( allocated(error)) deallocate(error)
      if(allocated(error_RE)) deallocate(error_RE )
      if(allocated(error_IM)) deallocate(error_IM )
      if(allocated(average_polar_matr)) deallocate(average_polar_matr )
 
      end subroutine check_the_error_dmft
