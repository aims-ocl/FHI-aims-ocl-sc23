!****s* FHI-aims/add_HF_to_first_order_H_polarizability
!  NAME
!   add_HF_to_first_order_H_polarizability
!  SYNOPSIS

subroutine add_HF_to_first_order_H_polarizability &
       ( first_order_density_matrix, first_order_H)

! PURPOSE
! add first_order HF (fock matrix) to first_order_H_polarizability. 
! This RI-V version of HF, ref: evaluate_exchange_matr_v0.f90
! 
! USES

 use dimensions
 use prodbas
 use hartree_fock
 use runtime_choices
 use mpi_tasks
 use synchronize_mpi
 use constants
 use localorb_io, only: use_unit

 implicit none

! ARGUMENTS

 real*8, dimension(3, n_basis,n_basis, n_spin), intent(in) :: first_order_density_matrix
 real*8, dimension(3, n_basis,n_basis, n_spin), intent(inout) :: first_order_H


!  INPUTS
!  o  occ_numbers -- real array,
!       the occupation number of the electrons for each eigenstate and each spin
!  o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle (KS/HF) self-consistent calculation
!  OUTPUTS
!  o  first_order_H
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!     auxiliary matrices for Level 3 Blas matrix multiplications

 real*8, dimension(:,:), allocatable ::  ovlp_denmat
 real*8, dimension(:,:), allocatable ::  temp_ovlp_matr
 real*8, dimension(:,:), allocatable ::  fock_matrix

!     counters

 integer :: i_coord

 integer i_spin
 integer i_basis_1
 integer i_basis_2
 integer i_basis_3
 integer i_index

!     begin work


 if(myid.eq.0) then
   write(use_unit,'(2X,A)')"Evaluating first_order_exchange matrix ... "
 endif

 allocate(temp_ovlp_matr(n_basis,n_basis))
 allocate(ovlp_denmat(n_basis, n_basis))
 allocate(fock_matrix(n_basis, n_basis))

 do i_coord = 1, 3
 do i_spin = 1, n_spin
 
 fock_matrix(:,:) = 0.d0
   do i_basis_1 = 1, n_loc_prodbas, 1
     i_index = 0
     temp_ovlp_matr(:,:) = 0.d0
     do i_basis_2 = 1, n_basis, 1
       do i_basis_3 = 1, i_basis_2, 1

         i_index = basis_nghr(i_basis_3, i_basis_2)
         if(i_index .gt. 0) then
           temp_ovlp_matr(i_basis_3, i_basis_2) = &
            ovlp_3fn(i_index, i_basis_1)
         endif


        enddo
      enddo
      ovlp_denmat(:,:) = 0.d0
      call dsymm('L', 'U', n_basis, n_basis, &
                  1.d0, &
                  temp_ovlp_matr, n_basis, &
                  first_order_density_matrix(i_coord,:,:,i_spin), n_basis, &
                  0.d0, &
                  ovlp_denmat, n_basis &
                 )

      call dsymm('R', 'U', n_basis, n_basis, &
                  1.d0, &
                  temp_ovlp_matr, n_basis, &
                  ovlp_denmat, n_basis, &
                  1.d0, &
                  fock_matrix(:,:), n_basis &
                 )

   enddo ! end of loop over i_basis_1

   if(use_mpi) then
     call sync_matrix(fock_matrix(:,:), &
               n_basis, n_basis)
   endif


   if(n_spin.eq.1) then 
     first_order_H(i_coord,:,:,i_spin) = first_order_H(i_coord,:,:,i_spin) - 0.50d0*hybrid_coeff*fock_matrix(:,:)
   else ! n_spin =2 
     first_order_H(i_coord,:,:,i_spin) = first_order_H(i_coord,:,:,i_spin) - 1.0d0*hybrid_coeff*fock_matrix(:,:)
   endif 


 enddo ! i_spin

 enddo ! i_coord

 if(allocated(ovlp_denmat)) then
    deallocate(ovlp_denmat)
 endif
 if(allocated(temp_ovlp_matr)) then
    deallocate(temp_ovlp_matr)
 endif

 end subroutine add_HF_to_first_order_H_polarizability
!---------------------------------------------------------------------
!******
