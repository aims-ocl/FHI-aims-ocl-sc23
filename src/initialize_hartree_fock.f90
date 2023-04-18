!****s* FHI-aims/initialize_hartree_fock
!  NAME
!   initialize_hartree_fock
!  SYNOPSIS

      subroutine initialize_hartree_fock &
             ()

!  PURPOSE  
!  prepare all the necessary ingredients for Hartree-Fock and post-Hartree-Fock
!  calcualtions. These includes
!  * construct the auxiliary basis, 
!  * calculate the 3-center O integrals 
!  * calcualte the 2-center Coulomb interacton integrals
!  * multiply O to the square root of the Coulomb integral.

!  USES
      use dimensions
      use runtime_choices
      use prodbas
      use physics
      use hartree_fock
      use hartree_fock_p0
      use lvl_triples
      use tight_binding_auxmat
      use timing
      use synchronize_mpi_basic
      use mpi_tasks
      use lc_wpbeh
      use lrc_pt2
      use species_data
      use localorb_io, only: use_unit
      implicit none 

!  INPUTS
!  none
!  OUTPUTS
!  none
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


!    for debug
      integer :: i, j

      integer :: i_spin, i_state
      integer :: info
      real*8  :: maxR

      integer*8 :: array_dimension

      character(*), parameter :: func = 'initialize_hartree_fock'

      call get_timestamps(time_prodbas_total, clock_time_prodbas_total)

      ! --- Initialize

      ! First, clean HF quantities which depend on pairs
      ! (might change with relaxation)
      call cleanup_basbas()

      n_homo_max = n_electrons/2.d0 + 2  ! Is needed for allocations.
      n_homo_max = min(n_homo_max, n_states)

      ! Construct product basis & basis pairs
      call initialize_prodbas()

      ! Allocate Hartree-Fock related arrays
      call allocate_hartree_fock()
      fock_matr(:,:,:) = 0.d0

      ! --- n_homo & occ_numbers

      ! Needs to be done after allocate_hartree_fock().
      n_homo = 0
      do i_state = 1, n_states
       do i_spin = 1, n_spin
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
         endif
       enddo
      enddo
      if (n_spin .eq. 2 .and. use_hf_multiplicity) then
        n_homo (1) =  int((n_electrons+1)/2.d0) + &
                      (hf_multiplicity - 1)/2
        n_homo (2) =  int((n_electrons)/2.d0) - &
                      (hf_multiplicity - 1)/2
        occ_numbers(1:n_homo(1),1,1) = 1.d0
        occ_numbers(n_homo(1)+1:n_states,1,1) = 0.d0
        occ_numbers(1:n_homo(2),2,1) = 1.d0
        occ_numbers(n_homo(2)+1:n_states,2,1) = 0.d0
      endif

      ! --- Calculate either ovlp_3fn or  coeff_3fn_ten&coulomb_matr_lvl

      if (.not. use_hf_kspace) then
        ! For LC-wPBEh the ovlp_3fn matrix has to be calculated twice,
        ! since we need the range-separation.
        !   ovlp_3fn will store the SR-Part
        !   ovlp_3fn_LR will store the LR-Part
        ! To keep mostyl everything concerning the LC-wPBEh functional together
        ! all new subroutines are kept in the lc_wpbeh module.
        if (use_lc_wpbeh) then
         ! Allocation of SR ovlp_3fn only if hybrid_coeff .ne. 0
          if (hybrid_coeff .ne. 0.0d0) then
		      allocate(ovlp_3fn_SR(n_basis_pairs,n_loc_prodbas), stat=info)
		      if (info.ne.0) then
		        array_dimension = n_basis_pairs*n_loc_prodbas*8
		        write(use_unit,*)
		        write(use_unit,*)
		        write(use_unit,'(1X,A,I8,A,A,A)') '*** Process ', myid, & 
		          ': Allocation of the ovlp_3fn_SR array failed in subroutine ', func, '().'
		        write(use_unit,'(1X,A,I8,A,I18,A)') '*** Process ', myid, & 
		          ': The memory requirement was ', array_dimension, ' Bytes.'
		        write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
		          ': Number of basis pairs = ', n_basis_pairs
		        write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
		          ': Number of product basis functions = ', n_loc_prodbas
		        write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
		          ': Perhaps retry with more memory (more CPUs to parallelize over).'
		        write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
		          ': Aborting the calculation - sorry.'
		        call aims_stop
		      end if
		  end if
          ! Allocation of LR ovlp_3fn
          if (allocated(ovlp_3fn)) deallocate(ovlp_3fn)
          allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=info)
          if (info.ne.0) then
            array_dimension = n_basis_pairs*n_loc_prodbas*8
            write(use_unit,*)
            write(use_unit,*)
            write(use_unit,'(1X,A,I8,A,A,A)') '*** Process ', myid, & 
              ': Allocation of the ovlp_3fn array failed in subroutine ', func, '().'
            write(use_unit,'(1X,A,I8,A,I18,A)') '*** Process ', myid, & 
              ': The memory requirement was ', array_dimension, ' Bytes.'
            write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
              ': Number of basis pairs = ', n_basis_pairs
            write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
              ': Number of product basis functions = ', n_loc_prodbas
            write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
              ': Perhaps retry with more memory (more CPUs to parallelize over).'
            write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
              ': Aborting the calculation - sorry.'
            call aims_stop
          end if
          ! Now call the function which will calculate all ovlp_3fn matrices
          call get_ovlp_3fn_for_lc_wpbeh()
        else if (lrc_pt2_started) then
          ! For lrc-PT2, the ovlp_3fn matrix should be the long-range
          ! corrected one
          if (allocated(ovlp_3fn)) deallocate(ovlp_3fn)
          allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=info)
          if (info.ne.0) then
            array_dimension = n_basis_pairs*n_loc_prodbas*8
            write(use_unit,*)
            write(use_unit,*)
            write(use_unit,'(1X,A,I8,A,A,A)') '*** Process ', myid, & 
              ': Allocation of the ovlp_3fn array failed in subroutine ' ,& 
              func, '().'
            write(use_unit,'(1X,A,I8,A,I18,A)') '*** Process ', myid, & 
              ': The memory requirement was ', array_dimension, &
              ' Bytes.'
            write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
              ': Number of basis pairs = ', n_basis_pairs
            write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
              ': Number of product basis functions = ', n_loc_prodbas
            write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
              ': Perhaps retry with more memory ', &
              '(more CPUs to parallelize over).'
            write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
              ': Aborting the calculation - sorry.'
            call aims_stop
          end if
          ! Now call the function which will calculate all ovlp_3fn matricies
          call get_ovlp_3fn_for_lrc_pt2()
        else
          if (sparse_o3fn) then
             if (.not. allocated(coulomb_matr_lvl)) then
               allocate(coulomb_matr_lvl(n_basbas, n_loc_prodbas), stat=info)
               call check_allocation(info, 'coulomb_matr_lvl', func)
             end if
          end if
          if (.not. sparse_o3fn .or. RI_type == RI_LVL_2nd) then
               if (allocated(ovlp_3fn)) deallocate(ovlp_3fn)
               allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=info)
               if (info.ne.0) then
                 ! in this case, the allocation failed. check_allocation would
                 ! report this, but we add a bit more information here because this is
                 ! the singular most likely array to fail.
                 array_dimension = n_basis_pairs*n_loc_prodbas*8
                 write(use_unit,*)
                 write(use_unit,*)
                 write(use_unit,'(1X,A,I8,A,A,A)') '*** Process ', myid, & 
                   ': Allocation of the ovlp_3fn array failed in subroutine ', func, '().'
                 write(use_unit,'(1X,A,I8,A,I18,A)') '*** Process ', myid, & 
                   ': The memory requirement was ', array_dimension, ' Bytes.'
                 write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
                   ': Number of basis pairs = ', n_basis_pairs
                 write(use_unit,'(1X,A,I8,A,I18)') '*** Process ', myid, & 
                   ': Number of product basis functions = ', n_loc_prodbas
                 write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
                   ': Perhaps retry with more memory (more CPUs to parallelize over).'
                 write(use_unit,'(1X,A,I8,A)') '*** Process ', myid, & 
                   ': Aborting the calculation - sorry.'
                 call aims_stop
               end if
          end if

          select case (RI_type)
          case(RI_SVS)
             call get_coeff_3fn_svs(ovlp_3fn)
          case(RI_V)
             if (use_2d_corr) then
                call get_coeff_3fn_v_2d(ovlp_3fn)
             else
                call get_coeff_3fn_v_1d(ovlp_3fn)
             end if
          case (RI_LVL)
             call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr_lvl)
          case (RI_LVL_full)
             call get_coeff_3fn_lvl_full(ovlp_3fn)
          case (RI_LVL_2nd)
             call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr_lvl)
           ! Bare Coulomb integrals to ovlp_3fn
             call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, OVLP_TYPE_COULOMB)
          case default
           call aims_stop("Invalid version of RI (resolution of identity)", func)
          end select        
        end if
        
!        write(use_unit,*) "ovlp_3fn matrix:"
!        do i= 1, n_basis_pairs
!            do j = 1, n_loc_prodbas
!                write(use_unit,*) i,j, ovlp_3fn(i,j)
!            end do
!        end do
!        stop
      else
! for periodic system
! Compute the LVL triple coefficients in real space

          call initialize_lvl_triples(OVLP_TYPE_COULOMB)
          call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)

          call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,lvl_tricoeff_recip2)

          call cleanup_lvl_triples()
          call deallocate_tb_auxmat()

!          write(use_unit,*) "cutCb_rcut:", cutCb_rcut
          if(use_hse .and. hse_omega_hf /= 0.d0 .and. .not. use_gw_and_hse &
              .and. .not. use_dftpt2_and_hse) then
            call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
          else
            call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
!            call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
!            call initialize_periodic_tb_auxmat(1, 1.d0)
          endif

          call determine_k_minus_q_list(kq_point_list,kpq_point_list)
          call get_coulomb_matr_recip(coulomb_matr_recip,1)
          call deallocate_tb_auxmat()

          call allocate_hartree_fock_p0()
! end if (.not. use_hf_kspace)
       endif
       call get_times(time_prodbas_total, clock_time_prodbas_total,& 
                   tot_time_prodbas_total, tot_clock_time_prodbas_total)

      return
      end subroutine initialize_hartree_fock

!******
