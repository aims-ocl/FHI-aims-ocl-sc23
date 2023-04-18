!****s* FHI-aims/initialize_cpt2_para
!  NAME
!   initialize_cpt2_para
!  SYNOPSIS

      subroutine initialize_cpt2_para(cpt2_para, weight_total,p_index)

!  PURPOSE
!  Subroutine initializes the relevant matricies required for PT2 calculation
!  in the "p_index" round.
!

! USES
      use dimensions
      use prodbas
      use pbc_lists
      use hartree_fock
      use mpi_tasks
      use synchronize_mpi
      use timing
      !
      use hartree_fock_p0
      use runtime_choices

      implicit none

! ARGUMENTS 

      integer :: cpt2_para(6,n_tasks)
      integer :: p_index 
      real*8  :: weight_total(n_tasks)

! INPUTS
!
! OUTPUT
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

!  local variables

      

      character(*), parameter :: func = 'initialize_cpt2_para.f90'
      integer :: info, mpierr

!     counters

      integer :: i_index, i_index_local, i_thread
      integer :: i_k_point, i_kp_point, i_q_point, i_qp_point
      integer :: i_qk_point, i_qpk_point, i_irqk_point
      integer :: i_index_block, n_block
      integer :: i_block_1, i_block_1_local

!     begin work

      cpt2_para     = 0
      i_index       = 0
      i_index_block = 0
      n_block       = n_k_points * n_k_points
      do i_qk_point = 1, n_k_points, 1

         if(.not. irk_point_included(i_qk_point) ) cycle

         i_index_block = i_index_block + 1

         i_block_1 = i_index_block * n_block
         i_block_1_local = (i_block_1-1)/n_tasks+1

         if (i_block_1_local .lt. p_index) then
             i_index = i_index + n_block
             cycle
         endif

         do i_k_point = 1, n_k_points, 1

             ! Determine the q grid
             i_q_point = kpq_point_list(i_k_point, i_qk_point)

            do i_qp_point = 1 , n_k_points, 1

                i_index = i_index + 1

                i_index_local = (i_index-1)/n_tasks+1

                if (i_index_local .lt. p_index) then
                    cycle
                elseif (i_index_local .gt. p_index) then
                    return
                endif

                ! distribute the task to the right thread
                i_thread = mod(i_index, n_tasks)

                if (i_thread .eq. 0) i_thread = n_tasks

                ! Determine the k' grid
                i_kp_point = kpq_point_list(i_qp_point,i_qk_point)

                ! Determine the q'-k grid
                i_qpk_point = kq_point_list(i_qp_point,i_k_point)

                ! Determine the k-mesh weight: 
                i_irqk_point = irk_point_mapping(i_qk_point)
                weight_total(i_thread) = irk_weight(i_irqk_point)*k_weights(i_k_point) &
                         * k_weights(i_qp_point)

                ! Determine the distribution of memory in this round "p_index"
                cpt2_para(1,i_thread) = i_k_point
                cpt2_para(2,i_thread) = i_kp_point
                cpt2_para(3,i_thread) = i_q_point
                cpt2_para(4,i_thread) = i_qp_point
                cpt2_para(5,i_thread) = i_qk_point
                cpt2_para(6,i_thread) = i_qpk_point

            enddo ! for i_qp_point
         enddo ! for i_k_point
      enddo ! for i_qk_point

      return

      end subroutine initialize_cpt2_para
