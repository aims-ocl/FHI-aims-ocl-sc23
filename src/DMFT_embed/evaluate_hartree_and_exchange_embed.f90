
      subroutine evaluate_hartree_and_exchange_embed(density_matrix, &
                                                     loc_self_enrg,&
                                                     main_loop_reiterations,&
                                                     inner_loop_reiterations, &
                                                     hartree_pot_LDA,&
                                                     hartree_pot_HF,& 
                                                     exchange_self_energy, &
                                                     on_site_xc_matr,&
                                                     full_ovlp_matrix_sqrt,& 
                                                     on_site_PBE_xc_matr,& 
                                                     k_summed_overlap_matr) 



! PURPOSE
! construct the exchange matrix, in this version of the construction density matrix
! is used. This reduces the memory consumption, but slows down the calculation 
! in particular we the occupied states are few.
! 
! USES
      use physics
      use dimensions
      use prodbas
      use hartree_fock
      use runtime_choices
      use mpi_tasks
      use synchronize_mpi_basic
      use constants
      use localized_basbas
      use dmft_para

      implicit none

! ARGUMENTS

      integer, intent(in) :: main_loop_reiterations 
      integer, intent(in) :: inner_loop_reiterations 

       real*8, dimension(:,:), allocatable ::  dens_matr_ovlp
       real*8, dimension(n_basis,n_basis), intent(in) ::  density_matrix
       real*8, dimension(n_basis,n_basis) :: hartree_pot_HF 
       real*8, dimension(n_basis,n_basis), intent(out):: hartree_pot_LDA 
       real*8, dimension(n_basis,n_basis) :: exchange_self_energy
       real*8, dimension (:,:,:) , allocatable  ::  aux_ovlp_3fn
       real*8, dimension(n_basis,n_basis), intent(out) :: loc_self_enrg 
       complex*16, dimension(n_basis,n_basis), intent(in) :: &
       k_summed_overlap_matr 
       real*8, dimension(n_basis,n_basis), intent(in) :: on_site_xc_matr 
       real*8, dimension(n_basis,n_basis), intent(in) :: on_site_PBE_xc_matr
       complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: &
       full_ovlp_matrix_sqrt


!     counters
   integer info

      integer :: i_state

      integer i_a,i_k_point,n,l
      integer i_basbas, j_basbas
      integer i_spin
      integer i_basis
      integer j_basis
      integer i_basis_1
      integer i_basis_2
      integer i_basis_3
      integer i_index, i_count
      logical :: output
      logical :: inner_loop_converged
      real*8 trace
   character*2 iter
   character*27 filename
       complex*16, dimension (:,:) , allocatable  ::  coeff_product 

   real*8, dimension(:), allocatable :: real_ipiv
   real*8, dimension(:), allocatable :: real_work
   real*8 :: symmetry_check


!     begin work

      if(myid.eq.0) then
        write(use_unit,'(2X,A)')"Evaluating the Hartree and &
        &Exchange matrix for the local DMFT self-energy... "
      endif

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      if(.not. allocated( aux_ovlp_3fn ))then
         allocate (aux_ovlp_3fn (n_basis,n_basis, n_loc_prodbas))!local
      endif
 
      if(.not. allocated (dens_matr_ovlp) )then
         allocate (dens_matr_ovlp(n_basis, n_basis))
      endif


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      aux_ovlp_3fn (:,:,:) = 0.d0
      do i_basis =1 , n_basis, 1
        do j_basis =1, n_basis, 1
          if(basis_nghr(i_basis,j_basis).gt.0)then
            aux_ovlp_3fn (i_basis, j_basis,:) = &
             ovlp_3fn(basis_nghr(i_basis,j_basis),:)
          else
            aux_ovlp_3fn (i_basis, j_basis,:) = 0.d0
          endif
        enddo
      enddo

exchange_self_energy(:,:) = 0.d0
hartree_pot_HF(:,:) = 0.d0
      do i_basbas = 1, n_loc_prodbas, 1
dens_matr_ovlp(:,:) = 0.d0


         call dgemm ('N','N', n_basis,&
            n_basis,n_basis,1.d0,&
            density_matrix , n_basis, &
            aux_ovlp_3fn (:,:,i_basbas), &
            n_basis,0.d0, &
            dens_matr_ovlp , n_basis )
!

         call dgemm ( 'N','N', n_basis,&
            n_basis,n_basis,1.d0,&
            aux_ovlp_3fn (:,:,i_basbas), &
            n_basis, dens_matr_ovlp ,&
            n_basis,1.d0, &
            exchange_self_energy , n_basis )



         trace = 0.d0
         do i_basis = 1, n_basis , 1
            trace =  trace + dens_matr_ovlp (i_basis,i_basis)
         enddo

         hartree_pot_HF(:,:) = hartree_pot_HF(:,:)+&
           2*trace*aux_ovlp_3fn (:,:,i_basbas)


enddo


       if(use_mpi) then
         call sync_matrix(exchange_self_energy(:,:), &
                   n_basis, n_basis)
         call sync_matrix(hartree_pot_HF(:,:), &
                   n_basis, n_basis)
       endif


if(main_loop_reiterations.eq.0)then
if(inner_loop_reiterations.eq.0)then

       hartree_pot_LDA (:,:) = hartree_pot_HF (:,:)

endif
endif

loc_self_enrg(:,:) = 0.d0
   do i_basis = 1, n_basis, 1
       do j_basis = 1, n_basis, 1

         loc_self_enrg(i_basis,j_basis) =loc_self_enrg(i_basis,j_basis) - &
                         1.d0*exchange_self_energy(i_basis,j_basis)


       enddo
   enddo





   if(myid.eq.0) then

     output = .false.

      if(output) then


          if( inner_loop_reiterations.lt.10 ) then
             write(iter,'(A,I1)') "0", inner_loop_reiterations
          else
             write(iter,'(I2)') inner_loop_reiterations
          endif
         
          filename = "diag_self_enrg_iter_"//iter//".dat"


          open(57, file=filename)

        do i_basis_2 = 1, i_basis_1, 1
!            write(57,*)'loc_self_energ, EXX, Hartree_LDA, Hartree_HF'
!       do i_basis_1 = 1, n_basis, 1
!               write(57,*)'Delta-Hartree, Delta-EXX', i_basis_1, &
!               hartree_pot_LDA(i_basis_1, i_basis_1)-hartree_pot_HF_loc(i_basis_1, i_basis_1),&
!               -exchange_self_energy(i_basis_1, i_basis_1)-on_site_xc_matr(i_basis_1, i_basis_1)

!         enddo
       enddo

 


          close(57)

!          open(65, file='self_enrg_VS_iter')

!          write(65,*) inner_loop_reiterations, &
!               loc_self_enrg(1,1), &
!               -aux_loc_self_enrg(1,1),&
!               hartree_pot_LDA(1,1),&
!               hartree_pot_HF_loc(1,1)


!          close(65)
      endif
      endif


      if(allocated(aux_ovlp_3fn)) then
         deallocate(aux_ovlp_3fn)
      endif
      if(allocated(dens_matr_ovlp)) then
         deallocate(dens_matr_ovlp)
      endif

      end subroutine evaluate_hartree_and_exchange_embed
!---------------------------------------------------------------------
!******
