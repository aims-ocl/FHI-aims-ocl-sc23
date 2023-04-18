
  subroutine get_pbc_and_cluster_quantities(full_hamiltonian, &
                                            full_ovlp_matrix, &
                                            !inv_full_ovlp_matrix,&
                                            free_cluster_ham, &
                                            on_site_xc_matr,&
                                            full_ovlp_matrix_sqrt, &
                                            !inv_full_ovlp_matrix_sqrt,&
                                            on_site_PBE_xc_matr,&
                                            on_site_PBE_x_matr,&
                                            on_site_PBE_c_matr,&
                                            full_k_PBE_xc_matr,&
                                            k_summed_overlap_matr)!
!on_site_xc_matr_no_ovlp) ! (KS_eigenvalue)

  use dimensions
  use runtime_choices
  use species_data
  use pbc_lists 
  use physics
  use prodbas
  use scgw_grid
  use constants
  use mpi_tasks
  use synchronize_mpi 
  use hartree_fock
  use localized_basbas
  use gw_para
  use dmft_para
  use gt
  use timing

   implicit none


! local parameters
   integer l,a,b,n, i_a, j_a 
   integer i_basis, i_basis_1, i_basis_2
   integer j_basis
   integer k_basis
   integer i_count
   integer i_state
   integer i_freq
   integer i_k_point, i_q_point
   integer ham_size
 
! subroutine actual parameters
   complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix
   complex*16, dimension(n_basis,n_basis, n_spin,n_k_points) :: &
   full_hamiltonian
   real*8, dimension(n_basis,n_basis)  :: free_cluster_ham
   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr
   real*8, dimension(n_basis,n_basis) :: on_site_PBE_xc_matr
   real*8, dimension(n_basis,n_basis) :: on_site_PBE_x_matr
   real*8, dimension(n_basis,n_basis) :: on_site_PBE_c_matr
   real*8, dimension(:,:), allocatable :: inv_k_summed_overlap_matr
   real*8, dimension(n_basis,n_basis):: k_summed_overlap_matr


   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 
   real*8, dimension(:), allocatable :: real_ipiv
   real*8, dimension(:), allocatable :: real_work 
   complex*16, dimension(:,:), allocatable :: full_inv_ovlp_matrix
   real*8, dimension(:), allocatable:: overlap_matrix_loc_w
   complex*16, dimension(:), allocatable:: overlap_matrix_loc_w_complex
   real*8, dimension(:,:), allocatable ::     hamiltonian_loc_w
   real*8, dimension(:,:), allocatable ::     xc_pot_w
   real*8, dimension(:,:), allocatable ::     x_pot_w
   real*8, dimension(:,:), allocatable ::     c_pot_w
   real*8, dimension(:,:), allocatable ::     on_site_xc_matr_temp
   real*8, dimension(:,:), allocatable ::     free_cluster_ham_temp
   complex*16, dimension(:,:), allocatable:: hamiltonian_loc_w_complex
   complex*16, dimension(:,:), allocatable:: xc_pot_w_complex
   complex*16, dimension(:,:), allocatable:: x_pot_w_complex
   complex*16, dimension(:,:), allocatable:: c_pot_w_complex
!   complex*16, dimension(n_basis,n_basis,n_spin,n_k_points) :: full_k_xc_matr
   complex*16, dimension(:,:,:), allocatable :: full_k_xc_matr
   complex*16, dimension(n_basis,n_basis,n_spin,n_k_points) :: &
   full_k_PBE_xc_matr
   complex*16, dimension(:,:,:),allocatable :: full_k_x_matr
   complex*16, dimension(:,:,:),allocatable :: full_k_c_matr
   complex*16, dimension(:,:),allocatable :: inv_full_ovlp_matrix_sqrt
   complex*16, dimension(:,:),allocatable  :: inv_full_ovlp_matrix
   complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix_sqrt
   complex*16, dimension(:,:), allocatable :: sum_full_ovlp_matrix_sqrt

   !complex*16, dimension(n_basis,n_basis,n_k_points) ::full_free_hamiltonian
   complex*16, dimension(:,:), allocatable ::full_free_hamiltonian
!   real*8, dimension(:,:), allocatable :: work_ovl
!   real*8, dimension(:,:), allocatable :: work_ham
   complex*16, dimension(:,:,:), allocatable:: ovlp_full_k_free_cluster_ham

      real*8, allocatable :: xc_realspace(:,:)
      real*8, allocatable :: x_realspace(:,:)
      real*8, allocatable :: c_realspace(:,:)
      real*8 :: alpha_PBE0


!------------------------------------------------------   

       if(.not.allocated(xc_realspace)) then
         allocate(xc_realspace(n_hamiltonian_matrix_size,n_spin))
       endif

       if(.not.allocated(x_realspace)) then
         allocate(x_realspace(n_hamiltonian_matrix_size,n_spin))
       endif

       if(.not.allocated(c_realspace)) then
         allocate(c_realspace(n_hamiltonian_matrix_size,n_spin))
       endif

                  call integrate_xc_realspace_p2 &
                                      (hartree_potential, &
                                      rho, rho_gradient,  &
                                      kinetic_density, &
                                      partition_tab, &
                                      l_shell_max, en_xc,&
                                      en_pot_xc, &
                                      xc_realspace ,& 
                                      x_realspace ,&
                                      c_realspace )
 
!       if (.not. allocated(work_ovl))then
!          allocate(work_ovl(n_centers_basis_I , n_centers_basis_I))
!       endif

!       if (.not. allocated(work_ham))then
!          allocate(work_ham(n_centers_basis_I , n_centers_basis_I))
!       endif

       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif



       if(.not.allocated(inv_k_summed_overlap_matr)) then
         allocate(inv_k_summed_overlap_matr(n_basis,n_basis))
       endif

free_cluster_ham(:,:)=0.d0
on_site_xc_matr(:,:)=0.d0
on_site_PBE_xc_matr(:,:)=0.d0
on_site_PBE_x_matr(:,:)=0.d0
on_site_PBE_c_matr(:,:)=0.d0
inv_k_summed_overlap_matr(:,:)=0.d0


     do i_k_point=1, n_k_points
                

       if (.not. allocated(hamiltonian_loc_w))then
          allocate(hamiltonian_loc_w(n_basis*(n_basis+1)/2,n_spin))
       endif

       if (.not. allocated(overlap_matrix_loc_w))then
          allocate(overlap_matrix_loc_w(n_basis*(n_basis+1)/2))
       endif

       if (.not. allocated(xc_pot_w))then
          allocate(xc_pot_w(n_basis*(n_basis+1)/2,n_spin))
       endif
       if (.not. allocated(xc_pot_w_complex))then
          allocate(xc_pot_w_complex(n_basis*(n_basis+1)/2,n_spin))
       endif
       if (.not. allocated(x_pot_w))then
          allocate(x_pot_w(n_basis*(n_basis+1)/2,n_spin))
       endif
       if (.not. allocated(x_pot_w_complex))then
          allocate(x_pot_w_complex(n_basis*(n_basis+1)/2,n_spin))
       endif
       if (.not. allocated(c_pot_w))then
          allocate(c_pot_w(n_basis*(n_basis+1)/2,n_spin))
       endif
       if (.not. allocated(c_pot_w_complex))then
          allocate(c_pot_w_complex(n_basis*(n_basis+1)/2,n_spin))
       endif
       if (.not. allocated(hamiltonian_loc_w_complex))then
          allocate(hamiltonian_loc_w_complex(n_basis*(n_basis+1)/2,n_spin))
       endif
       if (.not. allocated(overlap_matrix_loc_w_complex))then
          allocate(overlap_matrix_loc_w_complex(n_basis*(n_basis+1)/2))
       endif



                  call  construct_hamiltonian(hamiltonian, &
                                              hamiltonian_loc_w,&
                                              hamiltonian_loc_w_complex,&
                                              i_k_point)!,work_ham )

                  call  construct_overlap( overlap_matrix, &
                                           overlap_matrix_loc_w, &
                                           overlap_matrix_loc_w_complex, &
                                           i_k_point)!, work_ovl)


                  call  construct_xc_matr_kspace(xc_realspace, &
                                                 xc_pot_w, & 
                                                 xc_pot_w_complex, &
                                                 i_k_point)

                  call  construct_xc_matr_kspace(x_realspace,&
                                                 x_pot_w, &
                                                 x_pot_w_complex, &
                                                 i_k_point)

                  call  construct_xc_matr_kspace(c_realspace, &
                                                 c_pot_w, &
                                                 c_pot_w_complex, &
                                                 i_k_point)



       if (.not. allocated(full_k_x_matr))then
          allocate(full_k_x_matr(n_basis,n_basis,n_spin))
       endif

       if (.not. allocated(full_k_c_matr))then
          allocate(full_k_c_matr(n_basis,n_basis,n_spin))
       endif

       if (.not. allocated(full_k_xc_matr))then
          allocate(full_k_xc_matr(n_basis,n_basis,n_spin))
       endif


 if(real_eigenvectors) then


                    n = 0
                        do l = 1, n_basis
                           full_ovlp_matrix(1:l,l,i_k_point) = &
                           overlap_matrix_loc_w(n+1:n+l)   
                           ! filling the upper triangle
                           full_ovlp_matrix(l,1:l-1,i_k_point) = &
                           overlap_matrix_loc_w(n+1:n+l-1) 
                           ! filling the lower triangle

                           full_hamiltonian(1:l,l,1,i_k_point) = &
                           hamiltonian_loc_w (n+1:n+l,1)
                           full_hamiltonian(l,1:l-1,1,i_k_point) = &
                           hamiltonian_loc_w (n+1:n+l-1,1)

                           full_k_xc_matr(1:l,l,1) =   xc_pot_w (n+1:n+l,1)
                           full_k_xc_matr(l,1:l-1,1) = xc_pot_w (n+1:n+l-1,1)

                           full_k_x_matr(1:l,l,1) =   x_pot_w (n+1:n+l,1)
                           full_k_x_matr(l,1:l-1,1) = x_pot_w (n+1:n+l-1,1)

                           full_k_c_matr(1:l,l,1) =   c_pot_w (n+1:n+l,1)
                           full_k_c_matr(l,1:l-1,1) = c_pot_w (n+1:n+l-1,1)

                           n = n+l
                        enddo

       if (allocated(hamiltonian_loc_w))then
          deallocate(hamiltonian_loc_w)
       endif

       if ( allocated(overlap_matrix_loc_w))then
          deallocate(overlap_matrix_loc_w)
       endif

       if ( allocated(xc_pot_w))then
          deallocate(xc_pot_w)
       endif
       if ( allocated(x_pot_w))then
          deallocate(x_pot_w)
       endif
       if ( allocated(c_pot_w))then
          deallocate(c_pot_w)
       endif


 else



                    n = 0
                        do l = 1, n_basis

                           full_ovlp_matrix(1:l,l,i_k_point) = &
                           overlap_matrix_loc_w_complex(n+1:n+l)  
                           ! filling the upper triangle
                           full_ovlp_matrix(l,1:l-1,i_k_point) = &
                           dconjg(overlap_matrix_loc_w_complex(n+1:n+l-1))
                           ! filling the lower triangle


                           full_hamiltonian(1:l,l,1,i_k_point) = &
                           hamiltonian_loc_w_complex (n+1:n+l,1)
                           full_hamiltonian(l,1:l-1,1,i_k_point) = &
                           dconjg(hamiltonian_loc_w_complex (n+1:n+l-1,1))

                           full_k_xc_matr(1:l,l,1) = &
                           xc_pot_w_complex (n+1:n+l,1)
                           full_k_xc_matr(l,1:l-1,1) = &
                           dconjg(xc_pot_w_complex (n+1:n+l-1,1))

                           full_k_x_matr(1:l,l,1) = &
                           x_pot_w_complex (n+1:n+l,1)
                           full_k_x_matr(l,1:l-1,1) = &
                           dconjg(x_pot_w_complex (n+1:n+l-1,1))

                           full_k_c_matr(1:l,l,1) = &
                           c_pot_w_complex (n+1:n+l,1)
                           full_k_c_matr(l,1:l-1,1) = &
                           dconjg(c_pot_w_complex (n+1:n+l-1,1))


                           n = n+l
                        enddo

       if ( allocated(xc_pot_w_complex))then
          deallocate(xc_pot_w_complex)
       endif
       if ( allocated(hamiltonian_loc_w_complex))then
          deallocate(hamiltonian_loc_w_complex)
       endif
       if ( allocated(overlap_matrix_loc_w_complex))then
          deallocate(overlap_matrix_loc_w_complex)
       endif
       if ( allocated(x_pot_w_complex))then
          deallocate(x_pot_w_complex)
       endif
       if ( allocated(c_pot_w_complex))then
          deallocate(c_pot_w_complex)
       endif

 endif
!enddo

       if (.not. allocated(inv_full_ovlp_matrix))then
          allocate(inv_full_ovlp_matrix(n_basis,n_basis))
       endif


       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


      !inv_full_ovlp_matrix(:,:,i_k_point) = full_ovlp_matrix(:,:,i_k_point)
      inv_full_ovlp_matrix(:,:) = full_ovlp_matrix(:,:,i_k_point)

! at each frequency point evaluate the LU factorization, and check for errors
!    do i_k_point = 1, n_k_points, 1

         call zgetrf( n_basis, n_basis, inv_full_ovlp_matrix(:,:), &
                    n_basis , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_k_point
            write(use_unit,*) " * Error info = ", info
  !          stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call zgetri(n_basis, inv_full_ovlp_matrix(:,:), n_basis, &
                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_k_point
            write(use_unit,*) " * Error info = ", info
          endif
        endif


       if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

        full_ovlp_matrix_sqrt(:,:,i_k_point) = full_ovlp_matrix(:,:,i_k_point)   
        !full_ovlp_matrix_sqrt(:,:) = full_ovlp_matrix(:,:,i_k_point)   

     call power_genmat_lapack_complex(n_basis, &
                                 full_ovlp_matrix_sqrt(:,:,i_k_point), &
                                 0.5d0, &
                                 safe_minimum, 1.d-5, '')


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!---- BUILDING THE inverse of the SQRT OVERLAP MATRIX at each k-point -------

       if (.not. allocated(inv_full_ovlp_matrix_sqrt))then
          allocate(inv_full_ovlp_matrix_sqrt(n_basis,n_basis))
       endif

      inv_full_ovlp_matrix_sqrt(:,:) = full_ovlp_matrix_sqrt(:,:,i_k_point)

       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


! at each frequency point evaluate the LU factorization, and check for errors

         call zgetrf( n_basis, n_basis, inv_full_ovlp_matrix_sqrt(:,:), &
                    n_basis , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_k_point
            write(use_unit,*) " * Error info = ", info
  !          stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call zgetri(n_basis, inv_full_ovlp_matrix_sqrt(:,:), n_basis, &
                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_k_point
            write(use_unit,*) " * Error info = ", info
          endif
        endif


       if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif


!-----------------------------------------------------------------------------------------


       if (.not. allocated(full_free_hamiltonian))then
          allocate(full_free_hamiltonian(n_basis,n_basis))
       endif


    full_free_hamiltonian(:,:)=(0.d0,0.d0)
       do i_basis = 1, n_basis, 1
        do j_basis = 1, n_basis, 1

         full_free_hamiltonian(i_basis,j_basis)= &
         full_free_hamiltonian(i_basis,j_basis)+ &
         (full_hamiltonian(i_basis,j_basis, 1,i_k_point) &
         -full_k_xc_matr(i_basis,j_basis,1) )

        enddo
      enddo



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------






!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

 full_k_PBE_xc_matr(:,:,1,i_k_point) = (0.d0,0.d0)

        full_k_PBE_xc_matr(:,:,1,i_k_point) =  (1-hybrid_alpha)*&
                                               full_k_x_matr(:,:,1) &
                                               + 1*full_k_c_matr(:,:,1)


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!do i_k_point = 1, n_k_points

      free_cluster_ham(:,:) = free_cluster_ham(:,:) + &
                              full_k_PBE_xc_matr(:,:,1,i_k_point)

      on_site_xc_matr(:,:) = on_site_xc_matr(:,:) + &
                             full_k_xc_matr(:,:,1) 

      on_site_PBE_xc_matr(:,:) = on_site_PBE_xc_matr(:,:) + &
                                 full_k_PBE_xc_matr(:,:,1,i_k_point) 

      on_site_PBE_x_matr(:,:) = on_site_PBE_x_matr(:,:) + &
                                full_k_x_matr(:,:,1) 

      on_site_PBE_c_matr(:,:) = on_site_PBE_c_matr(:,:) + &
                                full_k_c_matr(:,:,1)

      inv_k_summed_overlap_matr(:,:) = inv_k_summed_overlap_matr(:,:) + &
                                       inv_full_ovlp_matrix(:,:)

       if ( allocated(full_free_hamiltonian))then
          deallocate(full_free_hamiltonian)
       endif

       if (allocated(full_k_x_matr))then
          deallocate(full_k_x_matr)
       endif

       if ( allocated(full_k_c_matr))then
          deallocate(full_k_c_matr)
       endif

       if ( allocated(full_k_xc_matr))then
          deallocate(full_k_xc_matr)
       endif

       if ( allocated(full_free_hamiltonian))then
          deallocate(full_free_hamiltonian)
       endif

       if ( allocated(inv_full_ovlp_matrix_sqrt))then
          deallocate(inv_full_ovlp_matrix_sqrt)
       endif

       if ( allocated(inv_full_ovlp_matrix))then
          deallocate(inv_full_ovlp_matrix)
       endif

enddo

       if(allocated(xc_realspace)) then
         deallocate(xc_realspace)
       endif
       if(allocated(x_realspace)) then
         deallocate(x_realspace)
       endif
       if(allocated(c_realspace)) then
         deallocate(c_realspace)
       endif


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

    free_cluster_ham(:,:)=(1./(n_k_points))*free_cluster_ham(:,:)
    on_site_xc_matr(:,:)=(1./(n_k_points))*on_site_xc_matr(:,:)
    on_site_PBE_xc_matr(:,:)=(1./(n_k_points))*on_site_PBE_xc_matr(:,:)
    on_site_PBE_x_matr(:,:)=(1./(n_k_points))*on_site_PBE_x_matr(:,:)
    on_site_PBE_c_matr(:,:)=(1./(n_k_points))*on_site_PBE_c_matr(:,:)
    inv_k_summed_overlap_matr(:,:)=(1./(n_k_points))*&
                                   inv_k_summed_overlap_matr(:,:)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


         k_summed_overlap_matr = inv_k_summed_overlap_matr


       if (.not. allocated(real_ipiv))then
          allocate(real_ipiv(n_basis))
       endif

       if (.not. allocated(real_work))then
          allocate(real_work (n_basis))
       endif


if(.true.) then
         call dgetrf( n_basis, n_basis, k_summed_overlap_matr, &
                    n_basis , real_ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Error info = ", info
  !          stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call dgetri(n_basis, k_summed_overlap_matr, n_basis, &
                 real_ipiv, real_work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Error info = ", info
          endif
        endif
endif

       if (allocated(real_ipiv))then
          deallocate(real_ipiv)
       endif

       if (allocated(real_work))then
          deallocate(real_work)
       endif

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------



!------------------------------------------------------   

!       if (allocated(work_ovl))then
!          deallocate(work_ovl)
!       endif
!       if (allocated(work_ham))then
!          deallocate(work_ham)
!       endif
       if (allocated(hamiltonian_loc_w))then
          deallocate(hamiltonian_loc_w)
       endif

       if ( allocated(overlap_matrix_loc_w))then
          deallocate(overlap_matrix_loc_w)
       endif

       if ( allocated(xc_pot_w))then
          deallocate(xc_pot_w)
       endif
       if ( allocated(xc_pot_w_complex))then
          deallocate(xc_pot_w_complex)
       endif
       if ( allocated(x_pot_w))then
          deallocate(x_pot_w)
       endif
       if ( allocated(x_pot_w_complex))then
          deallocate(x_pot_w_complex)
       endif
       if ( allocated(c_pot_w))then
          deallocate(c_pot_w)
       endif
       if ( allocated(c_pot_w_complex))then
          deallocate(c_pot_w_complex)
       endif

       if ( allocated(hamiltonian_loc_w_complex))then
          deallocate(hamiltonian_loc_w_complex)
       endif
       if ( allocated(overlap_matrix_loc_w_complex))then
          deallocate(overlap_matrix_loc_w_complex)
       endif

       if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work )
       endif





  end subroutine get_pbc_and_cluster_quantities 
