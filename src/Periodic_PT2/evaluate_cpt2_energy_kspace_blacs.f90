module cpt2_kspace_mod

  !  USES
  use dimensions
  use prodbas, only : atom2basbas_off, sp2n_basbas_sp, max_n_basbas_sp, basbas_atom
  !use hartree_fock, only : lvl_tricoeff_mod_r, n_homo_max, n_homo, kq_point_list, kpq_point_list, coulomb_matr_blacs
  use hartree_fock, only : n_homo_max, n_homo, kq_point_list, kpq_point_list, coulomb_matr_blacs
  !use physics, only: KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex
  use physics, only: KS_eigenvalue
  use hartree_fock_p0
  use mpi_tasks
  use synchronize_mpi
  use pbc_lists, only : Cbasis_to_atom
  use runtime_choices
  use timing
  use constants
  use pbc_lists
  use geometry, only : species
  use cpt2_blacs
  use exchange_trico_cpt2
  use exchange_ev_cpt2
  use exchange_coulomb_cpt2
  use restart_pt2
  implicit none

  save
  private
  public evaluate_cpt2_energy_kspace_blacs


  character(*), parameter :: func = 'evaluate_cpt2_energy_kspace_blacs.f90'

  integer, dimension(:), allocatable:: lb_atom, ub_atom

  !
  ! for the mp2 evaluation in each pattern
  !
  ! for calculation in COMPLEX
  complex(kind=8), dimension(:,:,:,:), allocatable, &
      target :: lvl_tricoeff_k, lvl_tricoeff_q ! LVL triple coefficients in k space
  complex(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r ! LVL triple coefficients in k space
  complex(kind=8), dimension(:,:,:,:), allocatable, &
      target :: lvl_tricoeff_kp, lvl_tricoeff_qp ! LVL triple coefficients in k space
  complex(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_recip_kp_r, lvl_tricoeff_recip_qp_r ! LVL triple coefficients in k space

  complex(kind=8), allocatable, dimension(:,:,:), &
      target :: KS_eigenvector_k, KS_eigenvector_kp
  complex(kind=8), allocatable, dimension(:,:,:), &
      target :: KS_eigenvector_q, KS_eigenvector_qp
  complex(kind=8), dimension(:,:,:), &
      pointer :: KS_eigenvector_q_r, KS_eigenvector_qp_r
  complex(kind=8), dimension(:,:,:), &
      pointer :: KS_eigenvector_k_r, KS_eigenvector_kp_r

  complex(kind=8), allocatable, dimension(:,:), &
      target :: aux_eri_kqkpqp, aux_eri_kqpkpq
  complex(kind=8), dimension(:,:), &
      pointer :: aux_eri_kqkpqp_r, aux_eri_kqpkpq_r

  complex(kind=8), allocatable, dimension(:,:,:,:), &
      target :: lvl_tricoeff_kq, lvl_tricoeff_kpqp
  complex(kind=8), allocatable, dimension(:,:,:,:), &
      target :: lvl_tricoeff_kpq, lvl_tricoeff_kqp
  complex(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_kq_r, lvl_tricoeff_kpqp_r
  complex(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_kpq_r, lvl_tricoeff_kqp_r

  complex(kind=8), allocatable, dimension(:,:), &
      target :: coulomb_matr_qk, coulomb_matr_qpk
  complex(kind=8), dimension(:,:), &
      pointer :: coulomb_matr_qk_r, coulomb_matr_qpk_r

  ! for calculations in REAL
  real(kind=8), dimension(:,:,:,:), allocatable, &
      target :: lvl_tricoeff_k_real, lvl_tricoeff_q_real
  real(kind=8), dimension(:,:,:,:), allocatable, &
      target :: lvl_tricoeff_kp_real, lvl_tricoeff_qp_real
  real(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_recip_k_r_real, lvl_tricoeff_recip_q_r_real
  real(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_recip_kp_r_real, lvl_tricoeff_recip_qp_r_real

  real(kind=8), allocatable, dimension(:,:,:), &
      target :: KS_eigenvector_k_real, KS_eigenvector_kp_real
  real(kind=8), allocatable, dimension(:,:,:), &
      target :: KS_eigenvector_q_real, KS_eigenvector_qp_real
  real(kind=8), dimension(:,:,:), &
      pointer :: KS_eigenvector_q_r_real, KS_eigenvector_k_r_real  
  real(kind=8), dimension(:,:,:), &
      pointer :: KS_eigenvector_kp_r_real, KS_eigenvector_qp_r_real

  real(kind=8), allocatable, dimension(:,:), target :: aux_eri_kqkpqp_real 
  real(kind=8), allocatable, dimension(:,:), target :: aux_eri_kqpkpq_real 
  real(kind=8), dimension(:,:), &
      pointer :: aux_eri_kqkpqp_r_real, aux_eri_kqpkpq_r_real

  real(kind=8), allocatable, dimension(:,:,:,:), &
      target :: lvl_tricoeff_kq_real, lvl_tricoeff_kpqp_real
  real(kind=8), allocatable, dimension(:,:,:,:), &
      target :: lvl_tricoeff_kpq_real, lvl_tricoeff_kqp_real
  real(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_kq_r_real, lvl_tricoeff_kpqp_r_real
  real(kind=8), dimension(:,:,:,:), &
      pointer :: lvl_tricoeff_kpq_r_real, lvl_tricoeff_kqp_r_real

  real(kind=8), allocatable, dimension(:,:), &
      target :: coulomb_matr_qk_real, coulomb_matr_qpk_real
  real(kind=8), dimension(:,:), &
      pointer :: coulomb_matr_qk_r_real, coulomb_matr_qpk_r_real

  !==============================================================
  ! about the generation of the k3-pattern list
  integer, allocatable, dimension(:,:) :: kq_pair_1, kq_pair_2, &
      kq_pair_e, kq_pair_ne
  real(kind=8), allocatable, dimension(:) :: total_weight, &
      total_weight_e, total_weight_ne
  integer :: max_pair_1, max_pair_2, s_block
  ! for kq_ne
  integer :: paircount_ne, paircount_local_ne, &
      lbpair_ne, ubpair_ne, max_pair_local_ne
  ! for kq_e
  integer :: paircount_e, paircount_local_e, &
      lbpair_e, ubpair_e, max_pair_local_e
  ! for kq_pair_2
  integer :: paircount, paircount_local, &
      lbpair, ubpair, max_pair_local
  !--------------------------------------------------------------
  !FIXME::   abandon at this stage, will think of usability in the future
  logical :: irk_point_symm = .false.
  integer, allocatable, dimension(:,:) :: irkq_mapping_cpt2
  integer, allocatable, dimension(:) :: irkq_mapping_num_cpt2
  !==============================================================

  ! mpi window related variables
  integer:: win_tri_occ, win_ev_occ, win_coul
  integer:: win_tri_unocc, win_ev_unocc
  integer :: mpierr

  ! record time consumption
  real(kind=8)  time_pt2(8)
  real(kind=8)  time_distribution(6)

  contains

!****s* FHI-aims/evaluate_cpt2_energy_kspace
!  NAME
!   evaluate_cpt2_energy_kspace
!  SYNOPSIS

  subroutine evaluate_cpt2_energy_kspace_blacs &
           ( occ_numbers, &
            pt2_c_energy, &
            pt2_c_energy_os, &
            pt2_c_energy_ss &
           )
      !KS_eigenvalue, &
      !KS_eigenvector, &
      !KS_eigenvector_complex, &

!  PURPOSE
!  Subroutine evaluate_cpt2_energy_kspace evaluates the correlation
!  energy at the PT2 level, which is the simplest post-SCF correlation.
!

! USES
      !use physics, only: KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue

      implicit none

! ARGUMENTS 

      real(kind=8)  :: occ_numbers(n_states,n_spin,n_k_points)

!     output
      real(kind=8)  :: pt2_c_energy
      real(kind=8)  :: pt2_c_energy_os
      real(kind=8)  :: pt2_c_energy_ss

! INPUTS
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the polarisability calcua            ltions
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate. In the present case, n_high_state >= n_homo
!            should be fine. 
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
! o  chemical_potential -- real number, the chemical potential of the system
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  KS_eigenvector_complex -- complex array,
!            the complex eigenvector of the single-particle calculation,
!            used when "real_eigenvectors == .false."
!           
!
! OUTPUT
! o  pt2_c_energy -- real number, the calculated PT2 correlation energy
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

    integer :: n_low_state
    integer :: n_high_state

    real(kind=8)  :: pt2_c_energy_local
    real(kind=8)  :: pt2_c_energy_start
    real(kind=8)  :: pt2_c_energy_local_os
    real(kind=8)  :: pt2_c_energy_start_os
    real(kind=8)  :: pt2_c_energy_local_ss
    real(kind=8)  :: pt2_c_energy_start_ss
      

!    local timing

    ! for parallelization
    integer       :: i_index, j_index, n_index, m_index
    real(kind=8)  :: pt2_term
    real(kind=8)  :: pt2_term_os
    real(kind=8)  :: pt2_term_ss

    character(*), parameter :: func = 'evaluate_cpt2_energy_kspace.f90'
    integer :: info
    character*150 :: info_str

    integer :: max_k_points_task

!   timing

    ! Three loops of k-mesh: 1) q-k (qk); 2) k; 3) q'(qp)
    integer :: i_qk_point
    integer :: i_qk_point_local
    integer :: i_k_point
    integer :: i_k_point_local
    integer :: i_qp_point
    integer :: i_qp_point_local
    ! We also need other three independent k-mesh arguments
    ! 1) q; 2) k'(kp); 3) q'-k (qpk)
    integer :: i_q_point
    integer :: i_q_point_local
    integer :: i_kp_point
    integer :: i_kp_point_local
    integer :: i_qpk_point
    integer :: i_qpk_point_local
    ! index for basis func, MO(occ, vir), spin, prodbas
    integer :: i_spin, a_state



    ! new counter
    integer :: oldatom, i_basis, i_atom

    ! for the loop of kq_pair_1
    integer :: lbsp,ubsp,n_sp_curr

    ! tmp irkq_mapping
    integer, allocatable, dimension(:,:) :: irkq_mapping_tmp
    integer :: irk_index, i_irk_point, i_irkp_point, max_ir_num

!   begin work

    if(myid.eq.0) then
      write(use_unit,*)
      write(use_unit,*)" -----------------------------------------------------------------------"
      write(use_unit,'(2X,A)') ' '
      write(use_unit,'(2X,A)') &
            "Start to calculate the periodic PT2 correlation energy  ... "
      write(use_unit,'(2X,A)') ' '
      write(use_unit,'(2X,A)') &
            "Parallel in terms of k-point, basis, and aux-basis ... "
      write(use_unit,'(2X,A)') ' '
      write(use_unit,"(2X,'n_tasks_kq : ',I3,', n_tasks_bl : ',I3,'X', I3)") &
            n_tasks_kq_cpt2,n_tasks_row_cpt2,n_tasks_col_cpt2
    endif
    call localorb_info('')

    allocate(lb_atom(n_atoms),ub_atom(n_atoms))
    oldatom=-1
    do i_basis = 1, n_basis
       i_atom = Cbasis_to_atom(i_basis)
       if(i_atom.ne.oldatom) then
          lb_atom(i_atom)=i_basis
          if (i_atom.gt.1) ub_atom(i_atom-1)=i_basis-1
          oldatom=i_atom
       end if
    end do
    ub_atom(n_atoms)=n_basis

    n_low_state = n_low_state_cpt2

    !because of the global mpi_fence operations, all mpi tasks have to run all iteration so of the kq loop
    if (mod(n_k_points,n_tasks_kq_cpt2).eq.0) then    
       max_k_points_task=n_k_points/n_tasks_kq_cpt2
    else
       max_k_points_task=n_k_points/n_tasks_kq_cpt2+1
    end if

    call get_timestamps(time_pt2(1), time_pt2(2) )

    if (real_eigenvectors) then
      if (gamma_only_cpt2) then
        call evaluate_memory_consumption_cpt2 (n_low_state,n_homo_max,n_lumo_min,3)
      else
        call evaluate_memory_consumption_cpt2 (n_low_state,n_homo_max,n_lumo_min,2)
      end if
    else
        call evaluate_memory_consumption_cpt2 (n_low_state,n_homo_max,n_lumo_min,1)
    end if

    ! =====================================================
    ! calculate the square root of coulomb_matr_blacs
    ! and replace the original array coulomb_matr_blacs
    do i_k_point_local = 1, n_kq_points_task
       if (real_eigenvectors) then
         call power_auxmat_scalapack_real_cpt2(n_basbas,0.5d0,&
             coulomb_matr_blacs_real(:,:,i_k_point_local),'')
       else
         call power_auxmat_scalapack_complex_cpt2(n_basbas,0.5d0,coulomb_matr_blacs(:,:,i_k_point_local),'')
       end if
    enddo

    call get_timestamps(time_pt2(3), time_pt2(4) )

    if (myid .eq. 0) then
        write(use_unit, "(2X,'Timing for V^{0.5} of coulomb_matr_recip     :',1X,F12.3,' s',4X,F12.3,' s')") &
            time_pt2(3)-time_pt2(1), &
            time_pt2(4)-time_pt2(2)
        write(use_unit, "(2X,' ')")
    endif

    !cpt2_total_grid = n_irk_points * n_k_points * n_k_points
    !cpt2_max_per_thread = (cpt2_total_grid-1) / n_tasks + 1

    time_distribution       = 0.d0

    ! =====================================================
    ! pt2_c_energy       :: the final PT2 correlation
    ! os :: opposite-spin pair
    ! ss :: parallel-spin pair
    pt2_c_energy_local      = 0.d0 
    pt2_c_energy            = 0.d0 
    pt2_c_energy_local_os   = 0.d0 
    pt2_c_energy_os         = 0.d0 
    pt2_c_energy_local_ss   = 0.d0 
    pt2_c_energy_ss         = 0.d0 
    ! =====================================================
    ! reading the restart info., Igor
    ! the story behind the following five lines is:
    ! 1) in the 0th thread, loading the pt2_c_energy_start for the restart file, 
    !       if restart_pt2_read = .true.
    ! 2) broadcasting pt2_c_energy_start to all threads from the 0th thread
    ! 3) broadcasting pt2_finish to all threads from the 0th thread
    pt2_finish              = 0
    current_pt2             = 0.0d0
    current_pt2_os          = 0.0d0
    current_pt2_ss          = 0.0d0
    call read_restart_pt2_info_blacs()
    pt2_c_energy_start      = current_pt2    ! only make sence for myid.eq.0 if restart_pt2_read = .true.
    pt2_c_energy_start_os   = current_pt2_os
    pt2_c_energy_start_ss   = current_pt2_ss
    call bcast_real(pt2_c_energy_start,0)
    call bcast_real(pt2_c_energy_start_os,0)
    call bcast_real(pt2_c_energy_start_ss,0)
    call mp_bcast_int(pt2_finish)
    ! =====================================================

    ! generate the kq list for the (k,k') patterns
    if (irk_point_symm) then
        max_pair_1 = n_irk_points*(n_irk_points+1)/2
        allocate(irkq_mapping_tmp(n_k_points,n_irk_points), stat=info)
        call check_allocation(info, 'irkq_mapping_tmp', func)
        allocate(irkq_mapping_num_cpt2(n_irk_points), stat=info)
        call check_allocation(info, 'irkq_mapping_num_cpt2', func)
        allocate(kq_pair_1(2,max_pair_1), stat=info)
        call check_allocation(info, 'kq_pair_1', func)
        paircount = 0
        irkq_mapping_num_cpt2=0
        irkq_mapping_tmp=0
        do i_k_point = 1, n_k_points, 1
          i_irk_point = irk_point_mapping(i_k_point)
          irkq_mapping_num_cpt2(i_irk_point) = &
              irkq_mapping_num_cpt2(i_irk_point) + 1
          irk_index = irkq_mapping_num_cpt2(i_irk_point)
          irkq_mapping_tmp(irk_index,i_irk_point) = i_k_point
          if (.not. irk_point_included(i_k_point)) cycle
          do i_kp_point = i_k_point, n_k_points, 1
             if (.not. irk_point_included(i_kp_point)) cycle
             i_irkp_point = irk_point_mapping(i_kp_point)
             paircount = paircount + 1
             kq_pair_1(1,paircount) = i_irk_point
             kq_pair_1(2,paircount) = i_irkp_point
          end do
        end do
        ! shrinking the irkq_mapping
        max_ir_num = maxval(irkq_mapping_num_cpt2)
        write(info_str, "(2X,'Max degeneracy in reducible k-grids : ',I6)") &
              max_ir_num
        allocate(irkq_mapping_cpt2(max_ir_num,n_irk_points), stat=info)
        call check_allocation(info, 'irkq_mapping_cpt2', func)
        do i_irk_point = 1, n_irk_points, 1
          irkq_mapping_cpt2(1:max_ir_num,i_irk_point) = &
              irkq_mapping_tmp(1:max_ir_num, i_irk_point)
        end do
        deallocate(irkq_mapping_tmp)
    else
        max_pair_1 = n_k_points*(n_k_points+1)/2
        allocate(kq_pair_1(2,max_pair_1), stat=info)
        call check_allocation(info, 'kq_pair_1', func)
        paircount = 0
        do i_k_point = 1, n_k_points, 1
          do i_kp_point = i_k_point, n_k_points, 1
             paircount = paircount + 1
             kq_pair_1(1,paircount) = i_k_point
             kq_pair_1(2,paircount) = i_kp_point
          end do
        end do
    end if

    if (real_eigenvectors) then
      if (gamma_only_cpt2) then
        call allocate_cpt2_kspace_gamma_only(n_low_state, n_homo_max)
      else
        call allocate_cpt2_kspace_real(n_low_state, n_homo_max)
      end if
    else
      call allocate_cpt2_kspace_complex(n_low_state, n_homo_max)
    end if


    if (set_pair_block_cpt2) then
        s_block = min(pair_block_cpt2,max_pair_1-pt2_finish)
    else
        s_block = 1
    end if

    lbsp=pt2_finish+1
    ubsp=min(lbsp+s_block-1,max_pair_1)
    n_sp_curr = ubsp-lbsp+1

    do while (lbsp.le.max_pair_1)

      call get_timestamps(time_pt2(7),time_pt2(8))
      call get_timestamps(time_pt2(3),time_pt2(4))

      if (irk_point_symm) then
        call generate_cpt2_kq_list_symm(lbsp, ubsp)
      else
        call generate_cpt2_kq_list(lbsp, ubsp)
      end if

      ! now calculate the mp2 energy according to kq_pair_2
      if (real_eigenvectors) then
        if (gamma_only_cpt2) then
          call cpt2_kq_list_gamma_only(paircount, kq_pair_2, total_weight, &
              n_low_state, n_homo_max, &
              n_homo_k, n_lumo_k, &
              pt2_c_energy_local, pt2_c_energy_local_ss, &
              pt2_c_energy_local_os)
        else
          call cpt2_kq_list_real(paircount, kq_pair_2, total_weight, &
              n_low_state, n_homo_max, &
              n_homo_k, n_lumo_k, &
              pt2_c_energy_local, pt2_c_energy_local_ss, &
              pt2_c_energy_local_os)
        end if
      else
        call cpt2_kq_list(paircount, kq_pair_2, total_weight, &
            n_low_state, n_homo_max, &
            n_homo_k, n_lumo_k, &
            pt2_c_energy_local, pt2_c_energy_local_ss, &
            pt2_c_energy_local_os)
      end if
            !n_low_state, occ_numbers, n_homo_max, &
            !n_homo_k, n_lumo_k, coulomb_matr_blacs, &
            !KS_eigenvector, KS_eigenvector_complex, &

      ! calculate the timing
      call get_timestamps(time_pt2(5),time_pt2(6))
      time_pt2(3) = time_pt2(5) - time_pt2(3)
      time_pt2(4) = time_pt2(6) - time_pt2(4)
      time_pt2(7) = time_pt2(3)
      time_pt2(8) = time_pt2(4)
      call sync_timing(time_pt2(3))
      call sync_timing(time_pt2(4))

      pt2_c_energy_ss = pt2_c_energy_local_ss
      pt2_c_energy_os = pt2_c_energy_local_os

      call sync_real_number(pt2_c_energy_ss)
      call sync_real_number(pt2_c_energy_os)

      !pt2_c_energy = pt2_c_energy_ss + pt2_c_energy_os + &
      !    pt2_c_energy_start

      ! Reflash the restarting file, Igor
      ! pt2_c_energy_local :: the PT2 correlation in each thread
      ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
      ! pt2_c_energy       :: the collected PT2 correlation up to now
      pt2_finish      = ubsp
      current_pt2_os  = pt2_c_energy_os + pt2_c_energy_start_os
      current_pt2_ss  = pt2_c_energy_ss + pt2_c_energy_start_ss
      current_pt2     = current_pt2_os + current_pt2_ss
      call write_restart_pt2_info_blacs()

      write(info_str, "(2X,'The progress of PT2 corr. : ',F8.3,'% (',F16.8,' Ha)')") &
           & real(ubsp)/real(max_pair_1)*100.d0, &
           current_pt2
      call localorb_info(info_str)
      write(info_str, "('Evaluate ',I8,' out of ',I8,' k3-tasks')") &
           & paircount, max_pair_2
      call output_timeheader('2X', info_str)
      call output_timer('Total cost', time_pt2(3:4))
      write(info_str, &
          "(2X,'| Cut down ',I6,' k3-tasks with permutation symmetry in (k,kp) and (q,qp)')") &
          & max_pair_2-paircount
      call localorb_info(info_str)
      write(info_str, &
          "(2X,'| Speed up ',I6,' k3-tasks with k=kp and q=qp')") &
          & paircount_e
      call localorb_info(info_str)
      call localorb_info('')

      lbsp = lbsp + s_block
      ubsp = min(lbsp+s_block-1,max_pair_1)
      n_sp_curr = ubsp-lbsp+1
    end do !do while (lbsp.le.max_pair_1)

    ! Note again:
    ! pt2_c_energy_local :: the PT2 correlation in each process
    ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
    ! pt2_c_energy       :: the final total PT2 correlation
    !call sync_real_number(pt2_c_energy_local)
    call sync_real_number(pt2_c_energy_local_os)
    call sync_real_number(pt2_c_energy_local_ss)
    !pt2_c_energy    = pt2_c_energy_start + pt2_c_energy_local
    pt2_c_energy_os = pt2_c_energy_start_os + pt2_c_energy_local_os
    pt2_c_energy_ss = pt2_c_energy_start_ss + pt2_c_energy_local_ss
    pt2_c_energy    = pt2_c_energy_os + pt2_c_energy_ss

    call get_timestamps(time_pt2(5),time_pt2(6))
    time_pt2(3) = time_pt2(5) - time_pt2(1)
    time_pt2(4) = time_pt2(6) - time_pt2(2)
    call sync_timing(time_pt2(3))
    call sync_timing(time_pt2(4))
    call sync_timing(time_distribution(1))
    call sync_timing(time_distribution(2))

    call localorb_info('')
    write(info_str,*)"----------------------------------------------------", &
                "-------------------------"
    call localorb_info(info_str)
    write(info_str,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
        " PT2 correlation        :", pt2_c_energy, "Ha,", &
         pt2_c_energy*hartree, "eV"
    call localorb_info(info_str)
    write(info_str,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
        " PT2 correlation (os)   :", pt2_c_energy_os, "Ha,", &
         pt2_c_energy_os*hartree, "eV"
    call localorb_info(info_str)
    write(info_str,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
        " PT2 correlation (ss)   :", pt2_c_energy_ss, "Ha,", &
         pt2_c_energy_ss*hartree, "eV"
    call localorb_info(info_str)
    call localorb_info('')
    call output_timeheader('2X', 'PT2 corr. calculation')
    call output_timer('Total cost', time_pt2(3:4))
    call output_timer('Communication', time_distribution(1:2))


    if (real_eigenvectors) then
      call deallocate_cpt2_kspace_real()
    else
      call deallocate_cpt2_kspace_complex()
    end if

  end subroutine evaluate_cpt2_energy_kspace_blacs

  subroutine generate_cpt2_kq_list(lbsp, ubsp)

    integer :: lbsp, ubsp

    ! counter
    integer :: i_pair_1, i_irqk_point
    integer :: i_k_point, i_kp_point
    integer :: i_q_point, i_qp_point
    integer :: i_qk_point, i_qpk_point
    integer :: i_k_point_tmp, i_kp_point_tmp

    integer :: lkq, lkq_eri_symm

    integer :: info

    integer :: i,j

    integer :: n_rest_tasks, n_rest_tasks_e, n_rest_tasks_ne
    integer :: start_rest, start_rest_e, start_rest_ne
    integer :: end_rest, end_rest_e, end_rest_ne

    if (allocated(kq_pair_2)) deallocate(kq_pair_2)
    if (allocated(kq_pair_e)) deallocate(kq_pair_e)
    if (allocated(kq_pair_ne)) deallocate(kq_pair_ne)

    if (allocated(total_weight)) deallocate(total_weight)
    if (allocated(total_weight_e)) deallocate(total_weight_e)
    if (allocated(total_weight_ne)) deallocate(total_weight_ne)


    max_pair_2 = 0
    do i_pair_1 = lbsp, ubsp, 1
      i_k_point = kq_pair_1(1,i_pair_1)
      i_kp_point = kq_pair_1(2,i_pair_1)
      if (i_k_point .eq. i_kp_point) then
        max_pair_2 = max_pair_2+n_irk_points
      else
        max_pair_2 = max_pair_2+2*n_irk_points
      end if
    end do

    allocate(kq_pair_e(6,max_pair_2), stat=info)
    call check_allocation(info, 'kq_pair_e', func)
    allocate(total_weight_e(max_pair_2), stat=info)
    call check_allocation(info, 'total_weight_e', func)

    allocate(kq_pair_ne(6,max_pair_2), stat=info)
    call check_allocation(info, 'kq_pair_ne', func)
    allocate(total_weight_ne(max_pair_2), stat=info)
    call check_allocation(info, 'total_weight_ne', func)

    paircount_e = 0
    total_weight_e = 0
    kq_pair_e = 0

    paircount_ne = 0
    total_weight_ne = 0
    kq_pair_ne = 0

    do i_pair_1 = lbsp, ubsp, 1

        i_k_point = kq_pair_1(1,i_pair_1)
        i_kp_point = kq_pair_1(2,i_pair_1)

        lkq = paircount_ne + 1
        lkq_eri_symm = paircount_e + 1

        do i_irqk_point = 1, n_irk_points, 1
          i_qk_point = inv_irk_point_mapping(i_irqk_point)

          i_k_point_tmp = i_k_point
          i_kp_point_tmp = i_kp_point

          i_q_point = kpq_point_list(i_k_point_tmp, i_qk_point)
          i_qp_point = kq_point_list(i_kp_point_tmp, i_qk_point)
          i_qpk_point = kq_point_list(i_qp_point, i_k_point_tmp)

          call check_kqkpqp_pattern(paircount_ne, max_pair_2, lkq, kq_pair_ne, total_weight_ne, &
                            paircount_e, lkq_eri_symm, kq_pair_e, total_weight_e, &
                            i_k_point_tmp, i_q_point, i_kp_point_tmp, i_qp_point, &
                            i_irqk_point, i_qk_point, i_qpk_point)

          if (i_k_point .ne. i_kp_point) then
            i_k_point_tmp = i_kp_point
            i_kp_point_tmp = i_k_point

            i_q_point = kpq_point_list(i_k_point_tmp, i_qk_point)
            i_qp_point = kq_point_list(i_kp_point_tmp, i_qk_point)
            i_qpk_point = kq_point_list(i_qp_point, i_k_point_tmp)

            call check_kqkpqp_pattern(paircount_ne, max_pair_2, lkq, kq_pair_ne, total_weight_ne, &
                              paircount_e, lkq_eri_symm, kq_pair_e, total_weight_e, &
                              i_k_point_tmp, i_q_point, i_kp_point_tmp, i_qp_point, &
                              i_irqk_point, i_qk_point, i_qpk_point)
          end if
        end do ! i_irqk_point = 1, n_irk_points, 1
    end do !

    ! reorder the pair list to minimize communication
    !   sorting the pair list in terms of i_qk_point
    call sorting_by_shell(paircount_ne,5,kq_pair_ne(:,1:paircount_ne),&
                          total_weight_ne(1:paircount_ne))
    !   further sorting the pair list in terms of i_k_point
    call fine_sorting_by_shell(paircount_ne,5,1,kq_pair_ne(:,1:paircount_ne),&
                          total_weight_ne(1:paircount_ne))
    !   sorting the pair list in terms of i_qk_point
    call sorting_by_shell(paircount_e,5,kq_pair_e(:,1:paircount_e),&
                          total_weight_e(1:paircount_e))
    !   further sorting the pair list in terms of i_k_point
    call fine_sorting_by_shell(paircount_e,5,1,kq_pair_e(:,1:paircount_e),&
                          total_weight_e(1:paircount_e))

    ! initialize the local tasks for kq_pair_ne
    n_rest_tasks_ne = mod(paircount_ne, n_tasks_kq_cpt2)
    if (n_rest_tasks_ne.eq.0) then
        start_rest_ne = paircount_ne
        end_rest_ne = paircount_ne-1
        max_pair_local_ne = paircount_ne/n_tasks_kq_cpt2
        paircount_local_ne = max_pair_local_ne
        lbpair_ne = max_pair_local_ne*myid_kq_cpt2 + 1
        ubpair_ne = lbpair_ne + paircount_local_ne - 1
    else
        start_rest_ne = 0
        end_rest_ne = n_rest_tasks_ne-1
        max_pair_local_ne = paircount_ne/n_tasks_kq_cpt2 + 1
        if (myid_kq_cpt2 .lt. n_rest_tasks_ne) then
            paircount_local_ne = max_pair_local_ne
            lbpair_ne = max_pair_local_ne*myid_kq_cpt2+1
            ubpair_ne = lbpair_ne + paircount_local_ne - 1
        else
            paircount_local_ne = max_pair_local_ne - 1
            lbpair_ne = n_rest_tasks_ne*max_pair_local_ne + &
                (myid_kq_cpt2-n_rest_tasks_ne)*paircount_local_ne + 1
            ubpair_ne = lbpair_ne + paircount_local_ne - 1
        end if
    end if

    ! initialize the local tasks for kq_eri_symm
    n_rest_tasks_e = mod(paircount_e, n_tasks_kq_cpt2)
    if (n_rest_tasks_e.eq.0) then
        start_rest_e = 0
        end_rest_e = -1
        max_pair_local_e = paircount_e/n_tasks_kq_cpt2
        paircount_local_e = max_pair_local_e
        lbpair_e = max_pair_local_e*myid_kq_cpt2 + 1
        ubpair_e = lbpair_e + paircount_local_e - 1
    else
        start_rest_e = n_tasks_kq_cpt2-n_rest_tasks_e
        end_rest_e = n_tasks_kq_cpt2-1
        max_pair_local_e = paircount_e/n_tasks_kq_cpt2 + 1
        if (myid_kq_cpt2 .lt. start_rest_e) then
            paircount_local_e = max_pair_local_e - 1
            lbpair_e = paircount_local_e*myid_kq_cpt2+1
            ubpair_e = lbpair_e + paircount_local_e - 1
        else
            paircount_local_e = max_pair_local_e
            lbpair_e = start_rest_e*(max_pair_local_e-1) + &
                (myid_kq_cpt2-start_rest_e)*max_pair_local_e + 1
            ubpair_e = lbpair_e + paircount_local_e - 1
        end if
    end if

    ! initialize the local indices for the final kq_pair_2 and total_weight
    paircount = paircount_e + paircount_ne

    n_rest_tasks = mod(paircount, n_tasks_kq_cpt2)
    if (n_rest_tasks.eq.0) then
        max_pair_local = paircount/n_tasks_kq_cpt2
        paircount_local = max_pair_local
        lbpair = max_pair_local*myid_kq_cpt2 + 1
        ubpair = lbpair + paircount_local - 1
    else
        max_pair_local = paircount/n_tasks_kq_cpt2 + 1
        paircount_local = paircount_local_e + paircount_local_ne
        if (end_rest_ne .lt. start_rest_e) then
            start_rest = end_rest_ne
            end_rest = start_rest_e-1
            if (myid_kq_cpt2 .lt. start_rest) then
                lbpair = max_pair_local*myid_kq_cpt2+1
                ubpair = lbpair + paircount_local - 1
            else if (myid_kq_cpt2 .ge. start_rest .and. &
                    myid_kq_cpt2 .le. end_rest) then
                lbpair = start_rest*max_pair_local + &
                    (myid_kq_cpt2-start_rest)*(max_pair_local-1) + 1
                ubpair = lbpair + paircount_local - 1
            else
                lbpair = start_rest*max_pair_local + &
                    (end_rest-start_rest+1)*(max_pair_local-1) + &
                    (myid_kq_cpt2-end_rest-1)*max_pair_local + 1
                ubpair = lbpair + paircount_local - 1
            end if
        else 
            start_rest = start_rest_e
            end_rest = end_rest_ne
            if (myid_kq_cpt2 .lt. start_rest) then
                lbpair = (max_pair_local-1)*myid_kq_cpt2+1
                ubpair = lbpair + paircount_local - 1
            else if (myid_kq_cpt2 .ge. start_rest .and. &
                    myid_kq_cpt2 .le. end_rest) then
                lbpair = start_rest*(max_pair_local-1) + &
                    (myid_kq_cpt2-start_rest)*max_pair_local + 1
                ubpair = lbpair + paircount_local - 1
            else
                lbpair = start_rest*(max_pair_local-1) + &
                    (end_rest-start_rest+1)*max_pair_local + &
                    (myid_kq_cpt2-end_rest-1)*(max_pair_local-1) + 1
                ubpair = lbpair + paircount_local - 1
            end if
        end if
    end if
    ! loading kq_pair_2 and total_weight and
    allocate(kq_pair_2(6,paircount), stat=info)
    call check_allocation(info, 'kq_pair_2', func)
    allocate(total_weight(paircount), stat=info)
    call check_allocation(info, 'total_weight', func)

    kq_pair_2 = 0
    total_weight = 0
    if (paircount_local_ne.gt.0) then
      kq_pair_2(:,lbpair:lbpair+paircount_local_ne-1) = kq_pair_ne(:,lbpair_ne:ubpair_ne)
      total_weight(lbpair:lbpair+paircount_local_ne-1) = total_weight_ne(lbpair_ne:ubpair_ne)
    end if
    if (paircount_local_e.gt.0) then
      kq_pair_2(:,ubpair-paircount_local_e+1:ubpair) = kq_pair_e(:,lbpair_e:ubpair_e)
      total_weight(ubpair-paircount_local_e+1:ubpair) = total_weight_e(lbpair_e:ubpair_e)
    end if

    !if (myid_bl_cpt2.eq.0) then
    !  write(use_unit,"('igor debug a2',I2,', paircount',3I5,', local',3I5,', lbpair',3I6)") myid_kq_cpt2, &
    !        paircount, paircount_e, paircount_ne, &
    !        paircount_local, paircount_local_e, paircount_local_ne, &
    !        lbpair, lbpair_e, lbpair_ne
    !end if

    deallocate(kq_pair_e, kq_pair_ne, total_weight_e, total_weight_ne)

  end subroutine generate_cpt2_kq_list

  subroutine check_kqkpqp_pattern(paircount, max_pair, lkq, kq_pair, total_weight, &
                              paircount_eri, lkq_eri, kq_eri, tw_eri, &
                              i_k_point, i_q_point, i_kp_point, i_qp_point, &
                              i_irqk_point,i_qk_point,i_qpk_point)
    implicit none

    integer, intent(in) :: max_pair, lkq, lkq_eri
    integer, intent(in) :: i_k_point, i_q_point, i_kp_point, i_qp_point
    integer, intent(in) :: i_irqk_point, i_qk_point, i_qpk_point
    integer, intent(inout) :: paircount, paircount_eri
    integer, dimension(6, max_pair), intent(inout) :: kq_pair, kq_eri
    real(kind=8), dimension(max_pair) :: total_weight, tw_eri
    ! tmp
    real(kind=8) :: k3_weight
    logical :: found, eri_symm
    integer :: i_pair

    ! determine if the symmetry in ERIs can be used or not.
    if (i_k_point.eq.i_kp_point .and. i_q_point.eq.i_qp_point) then
        eri_symm = .true.
    else
        eri_symm = .false.
    end if

    k3_weight = irk_weight(i_irqk_point)*k_weights(i_k_point) &
             * k_weights(i_kp_point)

    found=.false.
    if (.not. eri_symm) then
      do i_pair = lkq, paircount, 1
          if ((i_k_point.eq.kq_pair(3,i_pair)).and.&
              (i_kp_point.eq.kq_pair(1,i_pair)).and. &
              (i_q_point.eq.kq_pair(4,i_pair)).and. &
              (i_qp_point.eq.kq_pair(2,i_pair))) then
            found=.true.
            exit
          endif
      end do
      if (.not. found) then
          paircount = paircount + 1
          kq_pair(1,paircount) = i_k_point
          kq_pair(2,paircount) = i_q_point
          kq_pair(3,paircount) = i_kp_point
          kq_pair(4,paircount) = i_qp_point
          kq_pair(5,paircount) = i_qk_point
          kq_pair(6,paircount) = i_qpk_point
          total_weight(paircount) = k3_weight
      else
          total_weight(i_pair) = total_weight(i_pair) + k3_weight
      end if
    else
      do i_pair = lkq_eri, paircount_eri, 1
          if ((i_k_point.eq.kq_eri(3,i_pair)).and.&
              (i_kp_point.eq.kq_eri(1,i_pair)).and. &
              (i_q_point.eq.kq_eri(4,i_pair)).and. &
              (i_qp_point.eq.kq_eri(2,i_pair))) then
            found=.true.
            exit
          endif
      end do
      if (.not. found) then
          paircount_eri = paircount_eri + 1
          kq_eri(1,paircount_eri) = i_k_point
          kq_eri(2,paircount_eri) = i_q_point
          kq_eri(3,paircount_eri) = i_kp_point
          kq_eri(4,paircount_eri) = i_qp_point
          kq_eri(5,paircount_eri) = i_qk_point
          kq_eri(6,paircount_eri) = i_qpk_point
          tw_eri(paircount_eri) = k3_weight
      else
          tw_eri(i_pair) = tw_eri(i_pair) + k3_weight
      end if
    end if ! if (eri_symm)
  end subroutine check_kqkpqp_pattern


  subroutine generate_cpt2_kq_list_symm(lbsp, ubsp)

    integer :: lbsp, ubsp

    ! counter
    integer :: i_pair_1, i_irqk_point
    integer :: i_k_point, i_kp_point
    integer :: i_q_point, i_qp_point
    integer :: i_qk_point, i_qpk_point
    integer :: i_k_point_tmp, i_kp_point_tmp
    integer :: i,j
    integer :: lkq, ukq

    integer :: i_irk_point, i_irkp_point

    ! kq_pair_1 without irreducible symmetry
    integer, allocatable, dimension(:,:) :: decoded_kq_pair_1
    integer :: kq1_counter, i_kq1
    logical :: found

    integer :: info

    if (allocated(kq_pair_2)) deallocate(kq_pair_2)
    if (allocated(total_weight)) deallocate(total_weight)

    allocate(decoded_kq_pair_1(2,n_k_points*(n_k_points+1)/2), stat=info)

    ! At first, decode the given (k,k') pairs for this round of calculations.
    decoded_kq_pair_1 = 0
    kq1_counter = 0
    max_pair_2 = 0
    do i_pair_1 = lbsp, ubsp, 1
      i_irk_point = kq_pair_1(1,i_pair_1)
      i_irkp_point = kq_pair_1(2,i_pair_1)
      ! go through all possible (k,k') patterns with the same symmetry
      do i = 1, irkq_mapping_num_cpt2(i_irk_point), 1
        do j =1, irkq_mapping_num_cpt2(i_irkp_point), 1
          i_k_point = irkq_mapping_cpt2(i,i_irk_point)
          i_kp_point = irkq_mapping_cpt2(j,i_irkp_point)
          if (i_kp_point.eq.i_k_point) then
            kq1_counter = kq1_counter + 1
            decoded_kq_pair_1(1,kq1_counter) = i_k_point
            decoded_kq_pair_1(2,kq1_counter) = i_kp_point
            max_pair_2 = max_pair_2 + n_irk_points
          else if (i_kp_point.gt.i_k_point) then
            kq1_counter = kq1_counter + 1
            decoded_kq_pair_1(1,kq1_counter) = i_k_point
            decoded_kq_pair_1(2,kq1_counter) = i_kp_point
            max_pair_2 = max_pair_2 + 2*n_irk_points
          else
            found = .false.
            ! check if (k',k) has been included or not
            do i_kq1 = 1, kq1_counter, 1
              if ((decoded_kq_pair_1(1,i_kq1).eq.i_kp_point) .and. &
                  (decoded_kq_pair_1(2,i_kq1).eq.i_k_point)) then
                found = .true.
                exit
              end if
            end do
            if (.not. found) then
              kq1_counter = kq1_counter+1
              decoded_kq_pair_1(1,kq1_counter)=i_kp_point
              decoded_kq_pair_1(2,kq1_counter)=i_k_point
              max_pair_2 = max_pair_2 + 2*n_irk_points
            end if
          end if
        end do
      end do
    end do

    allocate(kq_pair_2(6,max_pair_2), stat=info)
    call check_allocation(info, 'kq_pair_2', func)
    allocate(total_weight(max_pair_2), stat=info)
    call check_allocation(info, 'total_weight', func)

    paircount = 0
    kq_pair_2 = 0
    do i_pair_1 = 1, kq1_counter, 1

      i_k_point = decoded_kq_pair_1(1,i_pair_1)
      i_kp_point = decoded_kq_pair_1(2,i_pair_1)

      lkq = 1

      do i_irqk_point = 1, n_irk_points, 1

        i_qk_point = inv_irk_point_mapping(i_irqk_point)

        i_k_point_tmp = i_k_point
        i_kp_point_tmp = i_kp_point

        i_q_point = kpq_point_list(i_k_point_tmp, i_qk_point)
        i_qp_point = kq_point_list(i_kp_point_tmp, i_qk_point)
        i_qpk_point = kq_point_list(i_qp_point, i_k_point_tmp)

        call check_kqkpqp_pattern_symm(paircount, max_pair_2, lkq,  &
                          kq_pair_2, total_weight, &
                          i_k_point_tmp, i_q_point, i_kp_point_tmp, i_qp_point, &
                          i_irqk_point, i_qk_point, i_qpk_point)

        if (i_k_point .ne. i_kp_point) then

          i_k_point_tmp = i_kp_point
          i_kp_point_tmp = i_k_point

          i_q_point = kpq_point_list(i_k_point_tmp, i_qk_point)
          i_qp_point = kq_point_list(i_kp_point_tmp, i_qk_point)
          i_qpk_point = kq_point_list(i_qp_point, i_k_point_tmp)

          call check_kqkpqp_pattern_symm(paircount, max_pair_2, lkq, &
                            kq_pair_2, total_weight, &
                            i_k_point_tmp, i_q_point, i_kp_point_tmp, i_qp_point, &
                            i_irqk_point,i_qk_point, i_qpk_point)
        end if

      end do ! i_irqk_point
    end do ! i_pair_1

    ! reorder the pair list to minimize communication
    !   sorting the pair list in terms of i_qk_point
    call sorting_by_shell(paircount,5,kq_pair_2(:,1:paircount),total_weight(1:paircount))
    !   further sorting the pair list in terms of i_k_point
    call fine_sorting_by_shell(paircount,5,1,kq_pair_2(:,1:paircount),total_weight(1:paircount))

    ! initialize the local tasks
    if (mod(paircount, n_tasks_kq_cpt2).eq.0) then
        max_pair_local = paircount/n_tasks_kq_cpt2
        paircount_local = max_pair_local
        lbpair = max_pair_local*myid_kq_cpt2 + 1
        ubpair = lbpair + paircount_local - 1
    else
        max_pair_local = paircount/n_tasks_kq_cpt2 + 1
        if (myid_kq_cpt2 .lt. mod(paircount, n_tasks_kq_cpt2)) then
            paircount_local = max_pair_local
            lbpair = max_pair_local*myid_kq_cpt2+1
            ubpair = lbpair + paircount_local - 1
        else
            paircount_local = max_pair_local - 1
            lbpair = mod(paircount, n_tasks_kq_cpt2)*max_pair_local + &
                (myid_kq_cpt2-mod(paircount, n_tasks_kq_cpt2))*paircount_local + 1
            ubpair = lbpair + paircount_local - 1
        end if
    end if
  end subroutine generate_cpt2_kq_list_symm

  subroutine check_kqkpqp_pattern_symm(paircount, max_pair, lkq, &
                              kq_pair, total_weight, &
                              i_k_point, i_q_point, i_kp_point, i_qp_point, &
                              i_irqk_point,i_qk_point,i_qpk_point)
    implicit none

    integer, intent(in) :: max_pair, lkq
    integer, intent(in) :: i_k_point, i_q_point, i_kp_point, i_qp_point
    integer, intent(in) :: i_irqk_point, i_qk_point, i_qpk_point
    integer, intent(inout) :: paircount
    integer, dimension(6, max_pair), intent(inout) :: kq_pair
    real(kind=8), dimension(max_pair) :: total_weight
    ! tmp
    real(kind=8) :: k3_weight
    logical :: found
    integer :: i_pair
    logical :: ir_symm_tmp
    integer :: i_irk_point, i_irkp_point, i_irq_point, i_irqp_point
    integer :: i_irk_point_ex, i_irkp_point_ex, i_irq_point_ex, i_irqp_point_ex

    k3_weight = irk_weight(i_irqk_point)*k_weights(i_k_point) &
             * k_weights(i_kp_point)

    found=.false.
    do i_pair = lkq, paircount, 1
      ! check the symmetry in permuting (k,k') and permuting (q,q')
      if ((i_k_point.eq.kq_pair(3,i_pair)).and.&
          (i_kp_point.eq.kq_pair(1,i_pair)).and. &
          (i_q_point.eq.kq_pair(4,i_pair)).and. &
          (i_qp_point.eq.kq_pair(2,i_pair))) then
        found=.true.
        exit
      endif
      ! check the symmetry accoding to irreducible k-points
      if ((i_qk_point.eq.kq_pair(5,i_pair)).and.&
          (i_qpk_point.eq.kq_pair(6,i_pair))) then

        i_irk_point=irk_point_mapping(i_k_point)
        i_irkp_point=irk_point_mapping(i_kp_point)
        i_irq_point=irk_point_mapping(i_q_point)
        i_irqp_point=irk_point_mapping(i_qp_point)

        i_irk_point_ex=irk_point_mapping(kq_pair(1,i_pair))
        i_irkp_point_ex=irk_point_mapping(kq_pair(3,i_pair))
        i_irq_point_ex=irk_point_mapping(kq_pair(2,i_pair))
        i_irqp_point_ex=irk_point_mapping(kq_pair(4,i_pair))

        ir_symm_tmp = (i_irk_point.eq.i_irk_point_ex) .and. &
                      (i_irkp_point.eq.i_irkp_point_ex) .and. &
                      (i_irq_point.eq.i_irq_point_ex) .and. &
                      (i_irqp_point.eq.i_irqp_point_ex)

        ir_symm_tmp = ir_symm_tmp .or. &
                      (i_irk_point.eq.i_irk_point_ex) .and. &
                      (i_irkp_point.eq.i_irkp_point_ex) .and. &
                      (i_irq_point.eq.i_irqp_point_ex) .and. &
                      (i_irqp_point.eq.i_irq_point_ex)

        ir_symm_tmp = ir_symm_tmp .or. &
                      (i_irk_point.eq.i_irkp_point_ex) .and. &
                      (i_irkp_point.eq.i_irk_point_ex) .and. &
                      (i_irq_point.eq.i_irq_point_ex) .and. &
                      (i_irqp_point.eq.i_irqp_point_ex)

        ir_symm_tmp = ir_symm_tmp .or. &
                      (i_irk_point.eq.i_irkp_point_ex) .and. &
                      (i_irkp_point.eq.i_irk_point_ex) .and. &
                      (i_irq_point.eq.i_irqp_point_ex) .and. &
                      (i_irqp_point.eq.i_irq_point_ex)

        if (ir_symm_tmp) then
            found = .true.
            exit
        end if
      end if
    end do
    if (.not. found) then
      paircount = paircount + 1
      kq_pair(1,paircount) = i_k_point
      kq_pair(2,paircount) = i_q_point
      kq_pair(3,paircount) = i_kp_point
      kq_pair(4,paircount) = i_qp_point
      kq_pair(5,paircount) = i_qk_point
      kq_pair(6,paircount) = i_qpk_point
      total_weight(paircount) = k3_weight
    else
      total_weight(i_pair) = total_weight(i_pair) + k3_weight
    end if
  end subroutine check_kqkpqp_pattern_symm

  subroutine sorting_by_shell(paircount,order_index, kq_pair,total_weight)
    ! FIXME :: not sure if the efficiency could be improved
    real(kind=8), parameter :: ALN2I=1./0.69314718,TINY=1.E-5
    integer,intent(in) :: order_index
    integer,intent(in) :: paircount
    integer,intent(inout) :: kq_pair(6,paircount)
    real(kind=8), intent(inout) :: total_weight(paircount)

    ! tmp variables
    integer :: m, nn, LOGNB2
    integer :: i,j,k,l
    integer :: t(6)
    real(kind=8)  :: tt

    LOGNB2=INT(ALOG(FLOAT(paircount))*ALN2I+TINY)
    m=paircount
    do nn=1,LOGNB2
      m=m/2; k=paircount-m
      do j=1,k
        i=j
 10     continue
        l=i+m
        if(kq_pair(order_index,l).LT.kq_pair(order_index,i)) then

          t=kq_pair(:,i)
          kq_pair(:,i)=kq_pair(:,l)
          kq_pair(:,l)=t

          tt=total_weight(i)
          total_weight(i)=total_weight(l)
          total_weight(l)=tt

          i=i-m
          if(i.GE.1) GOTO 10
        end if
      end do
    end do
    return
  end subroutine sorting_by_shell

  subroutine fine_sorting_by_shell(paircount,order_index, order_index_2, kq_pair,total_weight)
    ! FIXME :: not sure if the efficiency could be improved
    integer,intent(in) :: order_index, order_index_2
    integer,intent(in) :: paircount
    integer,intent(inout) :: kq_pair(6,paircount)
    real(kind=8), intent(inout) :: total_weight(paircount)

    ! tmp variables
    integer :: i,j
    logical :: found
    integer :: n_pair, lb_pairs, ub_pairs
    integer,allocatable,dimension(:,:) :: tabular_kq
    integer :: info

    if (allocated(tabular_kq)) deallocate(tabular_kq)
    allocate(tabular_kq(2,paircount), stat=info)
    

    tabular_kq = 0
    n_pair = 0
    ! determine the size and the indices of the block with the same kq_pair(order_index,:)
    do i = 1, paircount,1
      found = .false.
      do j = 1, n_pair, 1
        if (kq_pair(order_index,i).eq.tabular_kq(1,j)) then
            found=.true.
            exit
        end if
      end do
      if (found) then
          tabular_kq(2, j) = tabular_kq(2, j) + 1
      else
          n_pair = n_pair + 1
          tabular_kq(1, n_pair) = kq_pair(order_index,i)
          tabular_kq(2, n_pair) = 1
      end if
    end do
    
    ub_pairs = 0 
    do i=1,n_pair,1
      lb_pairs = ub_pairs + 1
      ub_pairs = lb_pairs + tabular_kq(2,i) - 1
      call sorting_by_shell(tabular_kq(2,i),order_index_2, &
                            kq_pair(:,lb_pairs:ub_pairs),&
                            total_weight(lb_pairs:ub_pairs))
    end do

  end subroutine fine_sorting_by_shell

  subroutine cpt2_kq_list(paircount, kq_pair, total_weight, &
                   n_low_state, n_homo_max, &
                   n_homo_k, n_lumo_k, &
                   pt2_c_energy_local, pt2_c_energy_local_ss, &
                   pt2_c_energy_local_os)
                   !n_homo_k, n_lumo_k, coulomb_matr_blacs_target, &
                   !KS_eigenvector, KS_eigenvector_complex, &
                   !KS_eigenvalue, pt2_c_energy_local, &

    integer, intent(in) :: paircount, n_low_state, n_homo_max
    integer, intent(in), dimension(6,paircount) :: kq_pair
    real(kind=8),  intent(in), dimension(paircount) :: total_weight
    
    integer, dimension(n_k_points,n_spin),intent(in) :: n_homo_k
    integer, dimension(n_k_points,n_spin),intent(in) :: n_lumo_k

    real(kind=8), intent(inout) :: pt2_c_energy_local, pt2_c_energy_local_os, pt2_c_energy_local_ss

    ! counter
    integer :: i_spin, i_spin_2, i_pair
    integer :: a_state, b_state, n_state, m_state
    integer :: m_state_shifted, n_state_shifted
    integer :: n_state_p, i_task_row

    real(kind=8)  :: k3_weight

    ! for k_point counter
    integer :: i_k_point, i_q_point, i_kp_point, i_qp_point
    integer :: i_qk_point,i_qpk_point
    integer :: i_k_point_old,i_q_point_old,i_kp_point_old,i_qp_point_old
    integer :: i_qk_point_old,i_qpk_point_old
    integer :: i_k_point_next,i_q_point_next,i_kp_point_next,i_qp_point_next
    integer :: i_qk_point_next,i_qpk_point_next



    ! record the PT2 correlation terms locally
    real(kind=8)      :: pt2_term, pt2_term_os, pt2_term_ss, e_diff
    complex(kind=8)  :: E_mp2_tmp_os, E_mp2_tmp_ss
    complex(kind=8)  :: E_mp2_a, E_mp2_b

    integer, dimension(6) :: kq_pair_new, kq_pair_old, kq_pair_next
    logical :: curr_run, next_run

    complex(kind=8), dimension(n_unocc_max,lpb_col:upb_col), &
        target :: aux_eri_tmp_d, aux_eri_tmp_x
    complex(kind=8), dimension(:,:), &
        pointer :: aux_eri_tmp_d_r, aux_eri_tmp_x_r

    complex(kind=8), pointer :: ptr_aux_eri_in_1(:,:)
    complex(kind=8), pointer :: ptr_aux_eri_in_2(:,:)
    complex(kind=8), pointer :: ptr_aux_eri_out(:,:)

    real*8   :: tmp_pt2_os, tmp_pt2_ss

    ! about the symmetry of the ab pair
    logical :: eri_symm
    integer :: b_state_start

    kq_pair_old = 0
    ! ================================================
    ! loading the k-point pattern for the first calculation in this loop
    ! ================================================
    if (((max_pair_local.eq.1) .and. &
        (paircount_local.lt.max_pair_local)) .or. &
        .not. kblacs_member) then
      kq_pair_new = 0
      curr_run = .false.
    else
      kq_pair_new = kq_pair(:,lbpair)
      curr_run = .true.
    end if
    i_k_point = kq_pair_new(1)
    i_qp_point = kq_pair_new(4)
    i_qk_point = kq_pair_new(5)
    i_qpk_point = kq_pair_new(6)
    !
    i_qk_point_old = 0
    i_qpk_point_old = 0
    ! load four lvl_tricoeff vectors
    call saccess_4_tricoeffs_cpt2_complex&
        (n_spin,win_tri_occ,win_tri_unocc,lvl_tricoeff_recip_k_r, &
         lvl_tricoeff_recip_q_r,lvl_tricoeff_recip_kp_r,lvl_tricoeff_recip_qp_r, &
         kq_pair_new,kq_pair_old,curr_run)
    ! load four KS_eigenvectors
    call saccess_4_evs_cpt2_complex&
        (n_spin,win_ev_occ,win_ev_unocc,KS_eigenvector_k_r,&
        KS_eigenvector_q_r,KS_eigenvector_kp_r,KS_eigenvector_qp_r,&
        kq_pair_new,kq_pair_old,curr_run)
    ! load two coulomb matrices
    call saccess_2_coulomb_cpt2_complex(n_spin,win_coul,&
        i_qk_point,i_qk_point_old,coulomb_matr_qk_r,&
        i_qpk_point,i_qpk_point_old,coulomb_matr_qpk_r)

    ! ================================================
    ! start the set of k-point patterns
    ! ================================================
    do i_pair = 1, max_pair_local, 1
      ! loading the k-point pattern for the this term
      if (((i_pair.eq.max_pair_local) .and. &
          (paircount_local.lt.max_pair_local)) .or. &
          .not. kblacs_member) then
        kq_pair_new = 0
        k3_weight = 0.0d0
        curr_run = .false.
      else
        kq_pair_new = kq_pair(:,lbpair+i_pair-1)
        k3_weight = total_weight(lbpair+i_pair-1)
        curr_run = .true.
      end if

      i_k_point = kq_pair_new(1)
      i_q_point = kq_pair_new(2)
      i_kp_point = kq_pair_new(3)
      i_qp_point = kq_pair_new(4)
      i_qk_point = kq_pair_new(5)
      i_qpk_point = kq_pair_new(6)

      ! loading the k-point pattern for the next term
      if ((i_pair+1.gt.max_pair_local) .or. &
          ((i_pair+1.eq.max_pair_local) .and. &
          (paircount_local.lt.max_pair_local)) .or. &
          .not. kblacs_member) then
        kq_pair_next = 0
        next_run = .false.
      else
        kq_pair_next = kq_pair(:,lbpair+i_pair)
        next_run = .true.
      end if

      i_k_point_next = kq_pair_next(1)
      i_q_point_next = kq_pair_next(2)
      i_kp_point_next = kq_pair_next(3)
      i_qp_point_next = kq_pair_next(4)
      i_qk_point_next = kq_pair_next(5)
      i_qpk_point_next = kq_pair_next(6)

      ! loading the k-point pattern for the previous term
      i_k_point_old = kq_pair_old(1)
      i_q_point_old = kq_pair_old(2)
      i_kp_point_old = kq_pair_old(3)
      i_qp_point_old = kq_pair_old(4)
      i_qk_point_old = kq_pair_old(5)
      i_qpk_point_old = kq_pair_old(6)

      call get_timestamps(time_distribution(3),time_distribution(4))
      ! ================================================
      ! load four lvl_tricoeff vectors
      ! ================================================
      call ssync_trico_cpt2_complex&
          (win_tri_occ,win_tri_unocc,lvl_tricoeff_recip_k_r, &
           lvl_tricoeff_recip_kp_r,&
           kq_pair_new,kq_pair_old,curr_run)
      ! ================================================
      ! load four KS_eigenvectors
      ! ================================================
      call ssync_ev_cpt2_complex&
          (win_ev_occ,win_ev_unocc,KS_eigenvector_k_r,&
           KS_eigenvector_kp_r,&
           kq_pair_new,kq_pair_old,curr_run)
      ! ================================================
      ! load two coulomb matrices
      ! ================================================
      call ssync_coulomb_cpt2(win_coul)

      !============
      ! Analysis the timing for the memory distribution.
      !============
      call get_timestamps(time_distribution(5), time_distribution(6))
      time_distribution(3) = time_distribution(5) - time_distribution(3)
      time_distribution(4) = time_distribution(6) - time_distribution(4)
      time_distribution(1) = time_distribution(1) + time_distribution(3)
      time_distribution(2) = time_distribution(2) + time_distribution(4)

      if (curr_run) then
        ! now prepare the M(k,q)
        if (.not. (i_k_point .eq. i_k_point_old .and. i_q_point .eq. i_q_point_old)) then
          if (i_k_point .eq. i_k_point_old .and. i_q_point .eq. i_qp_point_old) then
            lvl_tricoeff_kq_r = lvl_tricoeff_kqp_r
          else
            !lvl_tricoeff_kq_r => lvl_tricoeff_kq
            do i_spin = 1, n_spin
              !! FIXME :: for real eigenvectors, it should be an efficient treatment.
              call compute_m_kq_complex(n_low_state, i_spin, 1, i_k_point,&
                   i_q_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_k_r, KS_eigenvector_q_r, &
                   lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &  
                   coulomb_matr_qk_r, lvl_tricoeff_kq_r)
              ! end loop over i_spin
            enddo
          end if
        end if
        ! prepare the M(kp,qp)
        if (.not. (i_kp_point .eq. i_kp_point_old .and. i_qp_point .eq. i_qp_point_old)) then
          if (i_kp_point .eq. i_kp_point_old .and. i_qp_point .eq. i_q_point_old) then
            lvl_tricoeff_kpqp_r = lvl_tricoeff_kpq_r
          else
            !lvl_tricoeff_kpqp_r => lvl_tricoeff_kpqp
            do i_spin = 1, n_spin
              !! FIXME :: for real eigenvectors, it should be an efficient treatment.
              call compute_m_kq_complex(n_low_state, i_spin, 2, i_kp_point,&
                   i_qp_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_kp_r, KS_eigenvector_qp_r, &
                   lvl_tricoeff_recip_kp_r, lvl_tricoeff_recip_qp_r, &
                   coulomb_matr_qk_r, lvl_tricoeff_kpqp_r)
            enddo ! end loop over i_spin
          end if
        end if
        ! prepare the M(k,qp)
        if (.not. (i_k_point .eq. i_k_point_old .and. i_qp_point .eq. i_qp_point_old)) then
          if (i_qp_point .eq. i_q_point) then
            lvl_tricoeff_kqp_r = lvl_tricoeff_kq_r
          else
            do i_spin = 1, n_spin
              !! FIXME :: for real eigenvectors, it should be an efficient treatment.
              call compute_m_kq_complex(n_low_state, i_spin, 1, i_k_point,&
                   i_qp_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_k_r, KS_eigenvector_qp_r, &
                   lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_qp_r, &  
                   coulomb_matr_qpk_r, lvl_tricoeff_kqp_r)
            enddo ! end loop over i_spin
          end if
        end if
        ! prepare the M(kp,q)
        if (.not. (i_kp_point .eq. i_kp_point_old .and. i_q_point .eq. i_q_point_old)) then
          if (i_qp_point .eq. i_q_point) then
            lvl_tricoeff_kpq_r = lvl_tricoeff_kpqp_r
          else
            do i_spin = 1, n_spin
              !! FIXME :: for real eigenvectors, it should be another treatment.
              call compute_m_kq_complex(n_low_state, i_spin, 2, i_kp_point,&
                   i_q_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_kp_r, KS_eigenvector_q_r, &
                   lvl_tricoeff_recip_kp_r, lvl_tricoeff_recip_q_r, &  
                   coulomb_matr_qpk_r, lvl_tricoeff_kpq_r)
            enddo ! end loop over i_spin
          end if
        end if
      end if !if (curr_run) then
      call get_timestamps(time_pt2(7), time_pt2(8))

      ! load four lvl_tricoeff vectors for the next term
      call saccess_4_tricoeffs_cpt2_complex&
          (n_spin,win_tri_occ,win_tri_unocc,lvl_tricoeff_recip_k_r, &
           lvl_tricoeff_recip_q_r,lvl_tricoeff_recip_kp_r,lvl_tricoeff_recip_qp_r, &
           kq_pair_next,kq_pair_new,next_run)
      ! load four KS_eigenvectors for the next term
      call saccess_4_evs_cpt2_complex&
          (n_spin,win_ev_occ,win_ev_unocc,KS_eigenvector_k_r,&
          KS_eigenvector_q_r,KS_eigenvector_kp_r,KS_eigenvector_qp_r,&
          kq_pair_next,kq_pair_new,next_run)
      ! load two coulomb matrices for the next term
      call saccess_2_coulomb_cpt2_complex&
          (n_spin,win_coul,&
          i_qk_point_next,i_qk_point,coulomb_matr_qk_r,&
          i_qpk_point_next,i_qpk_point,coulomb_matr_qpk_r)

      ! if the current task is fake, jump to the next term after comunication.
      if (.not. curr_run) then
        kq_pair_old = 0
        cycle
      end if

      ! determine if the symmetry in ERIs can be used or not.
      if (i_k_point.eq.i_kp_point .and. i_q_point.eq.i_qp_point) then
          eri_symm = .true.
      else
          eri_symm = .false.
      end if

      ! Now do PT2 calculation
      ! FIXME :: at this moment, it supports close-shell cases only.
      pt2_term_os = 0.0d0
      pt2_term_ss = 0.0d0
      do i_spin = 1, n_spin, 1
      do i_spin_2 = 1, n_spin, 1
      do a_state = n_low_state, n_homo_k(i_k_point, i_spin), 1
        
        if (i_spin .eq. i_spin_2 .and. eri_symm) then
            b_state_start = a_state
        else
            b_state_start = n_low_state
        end if
        do b_state = b_state_start, n_homo_k(i_kp_point,i_spin_2), 1
          ptr_aux_eri_in_1 => lvl_tricoeff_kq_r(lbb_row:ubb_row,lpb_col:upb_col,a_state,i_spin)
          ptr_aux_eri_in_2 => lvl_tricoeff_kpqp_r(lbb_row:ubb_row,lpb_col:upb_col,b_state,i_spin_2)
          ptr_aux_eri_out => aux_eri_kqkpqp_r(lpb_row:upb_row,lpb_col:upb_col)
          call pzgemm('T', 'N', n_unocc_max, n_unocc_max, &
             n_basbas, (1.0d0,0.0d0), &
             ptr_aux_eri_in_1,1,1, pb2desc, &
             ptr_aux_eri_in_2,1,1, pb2desc, (0.d0,0.d0), &
             ptr_aux_eri_out, 1,1, pp2desc &
             )
          aux_eri_tmp_d = (0.0d0,0.0d0)
          aux_eri_tmp_d(lpb_row:upb_row,lpb_col:upb_col) =  &
              aux_eri_kqkpqp_r(lpb_row:upb_row,lpb_col:upb_col)
          call sync_vector_complex(aux_eri_tmp_d,size(aux_eri_tmp_d),comm_blacs_row_cpt2)

          if (i_spin .eq. i_spin_2) then
            if (eri_symm) then
              call pztranu( n_unocc_max, n_unocc_max, &
                   (1.0d0,0.0d0), aux_eri_kqkpqp_r, 1, 1, pp2desc, &
                   (0.0d0,0.0d0), aux_eri_kqpkpq_r, 1, 1, pp2desc )
            else
              aux_eri_kqpkpq_r = (0.0d0,0.0d0)
              ptr_aux_eri_in_1 => lvl_tricoeff_kpq_r(lbb_row:ubb_row,lpb_col:upb_col,b_state,i_spin)
              ptr_aux_eri_in_2 => lvl_tricoeff_kqp_r(lbb_row:ubb_row,lpb_col:upb_col,a_state,i_spin)
              ptr_aux_eri_out => aux_eri_kqpkpq_r(lpb_row:upb_row,lpb_col:upb_col)
              call pzgemm('T', 'N', n_unocc_max, n_unocc_max, &
                     n_basbas, (1.0d0,0.0d0), &
                     ptr_aux_eri_in_1,1,1, pb2desc, &
                     ptr_aux_eri_in_2,1,1, pb2desc, (0.d0,0.d0), &
                     ptr_aux_eri_out, 1,1, pp2desc &
                     )
            end if
            aux_eri_tmp_x = (0.0d0,0.0d0)
            aux_eri_tmp_x(lpb_row:upb_row,lpb_col:upb_col) =  &
                aux_eri_kqpkpq_r(lpb_row:upb_row,lpb_col:upb_col)
            call sync_vector_complex(aux_eri_tmp_x,size(aux_eri_tmp_x),comm_blacs_row_cpt2)
          end if

          do m_state = n_lumo_k(i_qp_point,i_spin_2), n_states_k(i_qp_point),1

            !! parallization with respect to m_state with i_qp_point in column
            !if (.not. (m_state .ge. lpb_col .and. m_state .le. upb_col)) cycle
            ! modify for shifted index
            m_state_shifted = m_state - n_lumo_min + 1
            if (.not. (m_state_shifted .ge. lpb_col .and. m_state_shifted .le. upb_col)) cycle
            ! end of modification

            do n_state = n_lumo_k(i_q_point,i_spin), n_states_k(i_q_point), 1

              !if (myid_row_cpt2 .ne. mod(n_state-n_lumo_min,n_tasks_row_cpt2)) cycle
              ! modify for shifted index
              n_state_shifted = n_state - n_lumo_min + 1
              if (myid_row_cpt2 .ne. mod(n_state_shifted-1,n_tasks_row_cpt2)) cycle
              ! end of modification

              e_diff = KS_eigenvalue(a_state,i_spin,i_k_point)  &
                     + KS_eigenvalue(b_state,i_spin_2,i_kp_point) &
                     - KS_eigenvalue(n_state,i_spin,i_q_point)  &
                     - KS_eigenvalue(m_state,i_spin_2,i_qp_point)
              if (abs(e_diff).lt.1e-6)then
                 write(use_unit,'(10X,A)') &
                      "****************************************"
                 write(use_unit,'(10X,2A)') "| Warning :", &
                   " too close to degeneracy"
                 write(use_unit,'(10X,A)') &
                      "****************************************"
              endif

              E_mp2_a = aux_eri_tmp_d(n_state_shifted,m_state_shifted)
              if (n_spin .eq. 1) then
                E_mp2_b = aux_eri_tmp_x(n_state_shifted,m_state_shifted)
                E_mp2_tmp_os = conjg(E_mp2_a)*E_mp2_a/e_diff
                E_mp2_tmp_ss = conjg(E_mp2_a)*(E_mp2_a - E_mp2_b)/e_diff
                !pt2_term_os = pt2_term_os + real(E_mp2_tmp_os)
                !pt2_term_ss = pt2_term_ss + real(E_mp2_tmp_ss)
              else
                if (i_spin .eq. i_spin_2) then
                  E_mp2_b = aux_eri_tmp_x(n_state_shifted,m_state_shifted)
                  E_mp2_tmp_ss = conjg(E_mp2_a)*(E_mp2_a - E_mp2_b)/e_diff
                  E_mp2_tmp_os = (0.0d0, 0.0d0)
                  !pt2_term_ss = pt2_term_ss + real(E_mp2_tmp_ss)
                else
                  E_mp2_tmp_os = conjg(E_mp2_a)*E_mp2_a/e_diff
                  E_mp2_tmp_ss = (0.0d0, 0.0d0)
                  !pt2_term_os = pt2_term_os + real(E_mp2_tmp_os)
                end if
              end if

              ! utilize the symmetry of ERIs
              if (a_state.ne.b_state .and. eri_symm) then
                  E_mp2_tmp_os = E_mp2_tmp_os*2.0d0
                  E_mp2_tmp_ss = E_mp2_tmp_ss*2.0d0
              end if

              pt2_term_os = pt2_term_os + real(E_mp2_tmp_os)
              pt2_term_ss = pt2_term_ss + real(E_mp2_tmp_ss)

            enddo ! n_state for the virtual states of i_kp_point
          enddo ! m_state for the virtual states of i_qp_point
        enddo ! b_state for the occupied states of i_q_point
      enddo ! a_state for the occupied states of i_k_point
      enddo ! i_spin
      enddo ! i_spin_2

      pt2_c_energy_local_os = pt2_c_energy_local_os + pt2_term_os * k3_weight

      pt2_c_energy_local_ss = pt2_c_energy_local_ss + pt2_term_ss * k3_weight

      kq_pair_old = kq_pair_new

    end do !do i_pair = 1, max_pair_local, 1

    call ssync_trico_cpt2_complex&
        (win_tri_occ,win_tri_unocc,lvl_tricoeff_recip_k_r, &
         lvl_tricoeff_recip_kp_r, &
         kq_pair_new,kq_pair_old,.false.)
    call ssync_ev_cpt2_complex&
        (win_ev_occ,win_ev_unocc,KS_eigenvector_k_r,&
         KS_eigenvector_kp_r, &
         kq_pair_new,kq_pair_old,.false.)
    call ssync_coulomb_cpt2(win_coul)

  end subroutine cpt2_kq_list

  subroutine compute_m_kq_complex(n_low_state, i_spin, iop, &
       i_k_point, i_q_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
       KS_eigenvector_k, KS_eigenvector_q, &
       lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &  
       coulomb_matr, lvl_tricoeff_kq)
  
    implicit none
  
    integer, intent(in) :: i_spin, i_k_point, i_q_point, iop, n_low_state
  
    integer, intent(in) :: n_lumo_k(n_k_points,n_spin), n_homo_k(n_k_points,n_spin)
    real(kind=8), dimension(n_states,n_spin,n_k_points), &
        intent(in) :: KS_eigenvalue
    complex(kind=8), dimension(n_basis,n_low_state:n_homo_max,n_spin), &
        intent(in) :: KS_eigenvector_k
    complex(kind=8), dimension(n_basis,lpb_col:upb_col,n_spin), &
        intent(in) :: KS_eigenvector_q
    complex(kind=8), dimension(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin), &
        intent(in) :: lvl_tricoeff_recip_k_r
    complex(kind=8), dimension(lbb_row:ubb_row,max_n_basis_sp,lpb_col:upb_col,n_spin), &
        intent(in) :: lvl_tricoeff_recip_q_r  
    complex(kind=8), dimension(lbb_row:ubb_row,lbb_col:ubb_col), &
        intent(in) :: coulomb_matr 
  
    !real*8, intent(in), dimension(n_states,n_spin,n_k_points) :: occ_numbers

    complex(kind=8), intent(inout), &
        dimension(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin), &
        target :: lvl_tricoeff_kq
  
    ! counter
    integer:: lb, ub, i_state, i_state_2, i_state_3, i_atom, i_species
    integer:: lbb_atom, ubb_atom, ind, brange

    ! a tmp pointer
    complex(kind=8), pointer :: ptr_lvl_tricoeff_kq_in(:,:,:)
  
!    call perfon('polfreq')

    do i_state=n_low_state, n_homo_k(i_k_point,i_spin)
      !do i_state_2=lpb_col, upb_col
      ! modify for shifted index
      do i_state_3=lpb_col, upb_col
        i_state_2 = n_lumo_min-1+i_state_3
      ! end of modification
        do i_atom=basbas_atom(lbb_row), basbas_atom(ubb_row)
          ! range in basis dimension
          lb=lb_atom(i_atom)
          ub=ub_atom(i_atom)
          brange=ub-lb+1
          ! range in basbas dimension
          i_species = species(i_atom)
          lbb_atom=max(lbb_row,atom2basbas_off(i_atom)+1)
          ubb_atom=min(ubb_row,atom2basbas_off(i_atom)+sp2n_basbas_sp(i_species))

          lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_3,i_state,i_spin)=0.
          do ind=1,brange
             lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_3,i_state,i_spin)=&
                  lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_3,i_state,i_spin) + &
                  conjg(KS_eigenvector_q(lb+ind-1,i_state_3,i_spin)) * &
                  lvl_tricoeff_recip_k_r(lbb_atom:ubb_atom,ind,i_state,i_spin) + &
                  KS_eigenvector_k(lb+ind-1,i_state,i_spin) * &
                  conjg(lvl_tricoeff_recip_q_r(lbb_atom:ubb_atom,ind,i_state_3,i_spin))
          end do
        end do
      end do
    end do

    ptr_lvl_tricoeff_kq_in => lvl_tricoeff_kq(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,i_spin)
    call get_v_multi_ovlp3fn_complex_blacs &
                 ( n_low_state, n_homo_max, n_homo_k(i_k_point,i_spin), coulomb_matr, &
                   ptr_lvl_tricoeff_kq_in, iop)
    !call get_v_multi_ovlp3fn_complex_blacs &
    !             ( n_low_state, n_homo_max, n_homo_k(i_k_point,i_spin), coulomb_matr, &
    !               lvl_tricoeff_kq(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,i_spin), iop)

  end subroutine compute_m_kq_complex



  subroutine allocate_cpt2_kspace_complex(n_low_state, n_homo_max)

    implicit none
    integer, intent(in) :: n_low_state, n_homo_max

    character(*), parameter :: func = 'evaluate_cpt2_energy_kspace_blacs.f90'
    integer :: info
    character*150 :: info_str

    ! create the windows for KS_eigenvectors, lvl_tricoeffs, and coulomb matrix
    call sinit_access_ev_cpt2(n_k_points,n_ks_points_task, &
                              KS_eigenvectors_occ,win_ev_occ, &
                              KS_eigenvectors_unocc, win_ev_unocc)

    call sinit_access_trico_cpt2(n_k_points,n_ks_points_task,&
                                   lvl_tricoeff_occ,win_tri_occ, &
                                   lvl_tricoeff_unocc,win_tri_unocc)

    call sinit_access_coulomb_cpt2(n_k_points,n_kq_points_task,coulomb_matr_blacs,win_coul)

    ! allocate relevant arrays

    if (n_tasks_kq_cpt2.gt.1) then
      allocate(KS_eigenvector_k(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_k', func)
      KS_eigenvector_k_r => KS_eigenvector_k
      allocate(KS_eigenvector_kp(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_kp', func)
      KS_eigenvector_kp_r => KS_eigenvector_kp

      allocate(lvl_tricoeff_k(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_k', func)
      lvl_tricoeff_recip_k_r => lvl_tricoeff_k
      allocate(lvl_tricoeff_kp(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_kp', func)
      lvl_tricoeff_recip_kp_r => lvl_tricoeff_kp

      allocate(KS_eigenvector_q(n_basis,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_q', func)
      KS_eigenvector_q_r => KS_eigenvector_q
      allocate(KS_eigenvector_qp(n_basis,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_qp', func)
      KS_eigenvector_qp_r => KS_eigenvector_qp

      allocate(lvl_tricoeff_q(lbb_row:ubb_row,max_n_basis_sp,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_q', func)
      lvl_tricoeff_recip_q_r => lvl_tricoeff_q
      allocate(lvl_tricoeff_qp(lbb_row:ubb_row,max_n_basis_sp,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_qp', func)
      lvl_tricoeff_recip_qp_r => lvl_tricoeff_qp

      allocate(coulomb_matr_qk(lbb_row:ubb_row,lbb_col:ubb_col),stat=info) 
      call check_allocation(info, 'coulomb_matr_qk', func)
      coulomb_matr_qk_r => coulomb_matr_qk
      allocate(coulomb_matr_qpk(lbb_row:ubb_row,lbb_col:ubb_col),stat=info) 
      call check_allocation(info, 'coulomb_matr_qpk', func)
      coulomb_matr_qpk_r => coulomb_matr_qpk
    else
      if (myid_col_cpt2.ne.0) then
        allocate(KS_eigenvector_k(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'KS_eigenvector_k', func)
        KS_eigenvector_k_r => KS_eigenvector_k
        allocate(KS_eigenvector_kp(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'KS_eigenvector_kp', func)
        KS_eigenvector_kp_r => KS_eigenvector_kp

        allocate(lvl_tricoeff_k(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'lvl_tricoeff_k', func)
        lvl_tricoeff_recip_k_r => lvl_tricoeff_k
        allocate(lvl_tricoeff_kp(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'lvl_tricoeff_kp', func)
        lvl_tricoeff_recip_kp_r => lvl_tricoeff_kp
      end if
    end if

    allocate(lvl_tricoeff_kq(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kq', func)
    lvl_tricoeff_kq_r=>lvl_tricoeff_kq
    allocate(lvl_tricoeff_kqp(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kqp', func)
    lvl_tricoeff_kqp_r=>lvl_tricoeff_kqp
    allocate(lvl_tricoeff_kpq(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kpq', func)
    lvl_tricoeff_kpq_r=>lvl_tricoeff_kpq
    allocate(lvl_tricoeff_kpqp(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kpqp', func)
    lvl_tricoeff_kpqp_r=>lvl_tricoeff_kpqp

    allocate(aux_eri_kqkpqp(lpb_row:upb_row,lpb_col:upb_col),stat=info)
    call check_allocation(info, 'aux_eri_kq', func)
    aux_eri_kqkpqp_r=>aux_eri_kqkpqp
    allocate(aux_eri_kqpkpq(lpb_row:upb_row,lpb_col:upb_col),stat=info)
    call check_allocation(info, 'aux_eri_kq', func)
    aux_eri_kqpkpq_r=>aux_eri_kqpkpq

  end subroutine allocate_cpt2_kspace_complex


  subroutine deallocate_cpt2_kspace_complex()

    call sfinalize_access_trico_cpt2(win_tri_occ,win_tri_unocc)
    call sfinalize_access_ev_cpt2(win_ev_occ,win_ev_unocc)
    call sfinalize_access_coulomb_cpt2(win_coul)

    if (allocated(KS_eigenvector_k)) deallocate(KS_eigenvector_k)
    if (allocated(KS_eigenvector_q)) deallocate(KS_eigenvector_q)
    if (allocated(KS_eigenvector_kp)) deallocate(KS_eigenvector_kp)
    if (allocated(KS_eigenvector_qp)) deallocate(KS_eigenvector_qp)

    if (allocated(lvl_tricoeff_k)) deallocate(lvl_tricoeff_k)
    if (allocated(lvl_tricoeff_q)) deallocate(lvl_tricoeff_q)
    if (allocated(lvl_tricoeff_kp)) deallocate(lvl_tricoeff_kp)
    if (allocated(lvl_tricoeff_qp)) deallocate(lvl_tricoeff_qp)

    if (allocated(coulomb_matr_qk)) deallocate(coulomb_matr_qk)
    if (allocated(coulomb_matr_qpk)) deallocate(coulomb_matr_qpk)

    if (allocated(lvl_tricoeff_kq)) deallocate(lvl_tricoeff_kq)
    if (allocated(lvl_tricoeff_kqp)) deallocate(lvl_tricoeff_kqp)
    if (allocated(lvl_tricoeff_kpq)) deallocate(lvl_tricoeff_kpq)
    if (allocated(lvl_tricoeff_kpqp)) deallocate(lvl_tricoeff_kpqp)

    if (allocated(aux_eri_kqkpqp)) deallocate(aux_eri_kqkpqp)
    if (allocated(aux_eri_kqpkpq)) deallocate(aux_eri_kqpkpq)

  end subroutine deallocate_cpt2_kspace_complex


  subroutine cpt2_kq_list_real(paircount, kq_pair, total_weight, &
                   n_low_state, n_homo_max, &
                   n_homo_k, n_lumo_k, &
                   pt2_c_energy_local, pt2_c_energy_local_ss, &
                   pt2_c_energy_local_os)

    integer, intent(in) :: paircount, n_low_state, n_homo_max
    integer, intent(in), dimension(6,paircount) :: kq_pair
    real(kind=8),  intent(in), dimension(paircount) :: total_weight
    
    integer, dimension(n_k_points,n_spin),intent(in) :: n_homo_k
    integer, dimension(n_k_points,n_spin),intent(in) :: n_lumo_k

    real(kind=8), intent(inout) :: pt2_c_energy_local, &
        pt2_c_energy_local_os, pt2_c_energy_local_ss

    ! counter
    integer :: i_spin, i_spin_2, i_pair
    integer :: a_state, b_state, n_state, m_state
    integer :: m_state_shifted, n_state_shifted
    integer :: n_state_p, i_task_row

    real(kind=8)  :: k3_weight

    ! for k_point counter
    integer :: i_k_point, i_q_point, i_kp_point, i_qp_point
    integer :: i_qk_point,i_qpk_point
    integer :: i_k_point_old,i_q_point_old,i_kp_point_old,i_qp_point_old
    integer :: i_qk_point_old,i_qpk_point_old
    integer :: i_k_point_next,i_q_point_next,i_kp_point_next,i_qp_point_next
    integer :: i_qk_point_next,i_qpk_point_next

    ! record the PT2 correlation terms locally
    real(kind=8)  :: pt2_term, pt2_term_os, pt2_term_ss, e_diff
    real(kind=8)  :: E_mp2_tmp_os, E_mp2_tmp_ss
    real(kind=8)  :: E_mp2_a, E_mp2_b

    integer, dimension(6) :: kq_pair_new, kq_pair_old, kq_pair_next
    logical :: curr_run, next_run

    real(kind=8), dimension(n_unocc_max,lpb_col:upb_col), &
        target :: aux_eri_tmp_d, aux_eri_tmp_x
    real(kind=8), dimension(:,:), &
        pointer :: aux_eri_tmp_d_r, aux_eri_tmp_x_r
    real(kind=8), pointer :: ptr_aux_eri_in_1(:,:)
    real(kind=8), pointer :: ptr_aux_eri_in_2(:,:)
    real(kind=8), pointer :: ptr_aux_eri_out(:,:)

    real(kind=8)   :: tmp_pt2_os, tmp_pt2_ss

    ! about the symmetry of the ab pair
    logical :: eri_symm
    integer :: b_state_start, n_symm_pattern

    kq_pair_old = 0
    n_symm_pattern = 0
    ! ================================================
    ! loading the k-point pattern for the first calculation in this loop
    ! ================================================
    if (((max_pair_local.eq.1) .and. &
        (paircount_local.lt.max_pair_local)) .or. &
        .not. kblacs_member) then
      kq_pair_new = 0
      curr_run = .false.
    else
      kq_pair_new = kq_pair(:,lbpair)
      curr_run = .true.
    end if
    i_k_point = kq_pair_new(1)
    i_qp_point = kq_pair_new(4)
    i_qk_point = kq_pair_new(5)
    i_qpk_point = kq_pair_new(6)
    !
    i_qk_point_old = 0
    i_qpk_point_old = 0
    ! load four lvl_tricoeff vectors
    call saccess_4_tricoeffs_cpt2_real(n_spin,win_tri_occ,win_tri_unocc,&
        lvl_tricoeff_recip_k_r_real, lvl_tricoeff_recip_q_r_real,&
        lvl_tricoeff_recip_kp_r_real,lvl_tricoeff_recip_qp_r_real, &
        kq_pair_new,kq_pair_old,curr_run)
    ! load four KS_eigenvectors
    call saccess_4_evs_cpt2_real(n_spin,win_ev_occ,win_ev_unocc,&
        KS_eigenvector_k_r_real,KS_eigenvector_q_r_real,&
        KS_eigenvector_kp_r_real,KS_eigenvector_qp_r_real,&
        kq_pair_new,kq_pair_old,curr_run)
    ! load two coulomb matrices
    call saccess_2_coulomb_cpt2_real(n_spin,win_coul,&
        i_qk_point,i_qk_point_old,coulomb_matr_qk_r_real,&
        i_qpk_point,i_qpk_point_old,coulomb_matr_qpk_r_real)

    ! ================================================
    ! start the set of k-point patterns
    ! ================================================
    do i_pair = 1, max_pair_local, 1
      ! loading the k-point pattern for the this term
      if (((i_pair.eq.max_pair_local) .and. &
          (paircount_local.lt.max_pair_local)) .or. &
          .not. kblacs_member) then
        kq_pair_new = 0
        k3_weight = 0.0d0
        curr_run = .false.
      else
        kq_pair_new = kq_pair(:,lbpair+i_pair-1)
        k3_weight = total_weight(lbpair+i_pair-1)
        curr_run = .true.
      end if

      i_k_point = kq_pair_new(1)
      i_q_point = kq_pair_new(2)
      i_kp_point = kq_pair_new(3)
      i_qp_point = kq_pair_new(4)
      i_qk_point = kq_pair_new(5)
      i_qpk_point = kq_pair_new(6)

      ! loading the k-point pattern for the next term
      if ((i_pair+1.gt.max_pair_local) .or. &
          ((i_pair+1.eq.max_pair_local) .and. &
          (paircount_local.lt.max_pair_local)) .or. &
          .not. kblacs_member) then
        kq_pair_next = 0
        next_run = .false.
      else
        kq_pair_next = kq_pair(:,lbpair+i_pair)
        next_run = .true.
      end if

      i_k_point_next = kq_pair_next(1)
      i_q_point_next = kq_pair_next(2)
      i_kp_point_next = kq_pair_next(3)
      i_qp_point_next = kq_pair_next(4)
      i_qk_point_next = kq_pair_next(5)
      i_qpk_point_next = kq_pair_next(6)

      ! loading the k-point pattern for the previous term
      i_k_point_old = kq_pair_old(1)
      i_q_point_old = kq_pair_old(2)
      i_kp_point_old = kq_pair_old(3)
      i_qp_point_old = kq_pair_old(4)
      i_qk_point_old = kq_pair_old(5)
      i_qpk_point_old = kq_pair_old(6)

      call get_timestamps(time_distribution(3),time_distribution(4))
      ! ================================================
      ! load four lvl_tricoeff vectors
      ! ================================================
      call ssync_trico_cpt2_real(win_tri_occ,win_tri_unocc,lvl_tricoeff_recip_k_r_real, &
           lvl_tricoeff_recip_kp_r_real, &
           kq_pair_new,kq_pair_old,curr_run)
      ! ================================================
      ! load four KS_eigenvectors
      ! ================================================
      call ssync_ev_cpt2_real(win_ev_occ,win_ev_unocc,KS_eigenvector_k_r_real,&
           KS_eigenvector_kp_r_real,&
           kq_pair_new,kq_pair_old,curr_run)
      ! ================================================
      ! load two coulomb matrices
      ! ================================================
      call ssync_coulomb_cpt2(win_coul)

      !============
      ! Analysis the timing for the memory distribution.
      !============
      call get_timestamps(time_distribution(5), time_distribution(6))
      time_distribution(3) = time_distribution(5) - time_distribution(3)
      time_distribution(4) = time_distribution(6) - time_distribution(4)
      time_distribution(1) = time_distribution(1) + time_distribution(3)
      time_distribution(2) = time_distribution(2) + time_distribution(4)

      if (curr_run) then
        ! now prepare the M(k,q)
        if (.not. (i_k_point .eq. i_k_point_old .and. i_q_point .eq. i_q_point_old)) then
          if (i_k_point .eq. i_k_point_old .and. i_q_point .eq. i_qp_point_old) then
            lvl_tricoeff_kq_r_real = lvl_tricoeff_kqp_r_real
          else
            do i_spin = 1, n_spin
              call compute_m_kq_real(n_low_state, i_spin, 1, i_k_point,&
                   i_q_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_k_r_real, KS_eigenvector_q_r_real, &
                   lvl_tricoeff_recip_k_r_real, lvl_tricoeff_recip_q_r_real, &  
                   coulomb_matr_qk_r_real, lvl_tricoeff_kq_r_real)
            enddo ! end loop over i_spin
          end if
        end if
        ! prepare the M(kp,qp)
        if (.not. (i_kp_point .eq. i_kp_point_old .and. i_qp_point .eq. i_qp_point_old)) then
          if (i_kp_point .eq. i_kp_point_old .and. i_qp_point .eq. i_q_point_old) then
            lvl_tricoeff_kpqp_r_real = lvl_tricoeff_kpq_r_real
          else
            do i_spin = 1, n_spin
              call compute_m_kq_real(n_low_state, i_spin, 2, i_kp_point,&
                   i_qp_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_kp_r_real, KS_eigenvector_qp_r_real, &
                   lvl_tricoeff_recip_kp_r_real, lvl_tricoeff_recip_qp_r_real, &
                   coulomb_matr_qk_r_real, lvl_tricoeff_kpqp_r_real)
            enddo ! end loop over i_spin
          end if
        end if
        ! prepare the M(k,qp)
        if (.not. (i_k_point .eq. i_k_point_old .and. i_qp_point .eq. i_qp_point_old)) then
          if (i_qp_point .eq. i_q_point) then
            lvl_tricoeff_kqp_r_real = lvl_tricoeff_kq_r_real
          else
            do i_spin = 1, n_spin
              call compute_m_kq_real(n_low_state, i_spin, 1, i_k_point,&
                   i_qp_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_k_r_real, KS_eigenvector_qp_r_real, &
                   lvl_tricoeff_recip_k_r_real, lvl_tricoeff_recip_qp_r_real, &  
                   coulomb_matr_qpk_r_real, lvl_tricoeff_kqp_r_real)
            enddo ! end loop over i_spin
          end if
        end if
        ! prepare the M(kp,q)
        if (.not. (i_kp_point .eq. i_kp_point_old .and. i_q_point .eq. i_q_point_old)) then
          if (i_qp_point .eq. i_q_point) then
            lvl_tricoeff_kpq_r_real = lvl_tricoeff_kpqp_r_real
          else
            do i_spin = 1, n_spin
              call compute_m_kq_real(n_low_state, i_spin, 2, i_kp_point,&
                   i_q_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
                   KS_eigenvector_kp_r_real, KS_eigenvector_q_r_real, &
                   lvl_tricoeff_recip_kp_r_real, lvl_tricoeff_recip_q_r_real, &
                   coulomb_matr_qpk_r_real, lvl_tricoeff_kpq_r_real)
            enddo ! end loop over i_spin
          end if
        end if
      end if ! if (curr_run) then
      call get_timestamps(time_pt2(7), time_pt2(8))

      ! load four lvl_tricoeff vectors for the next term
      call saccess_4_tricoeffs_cpt2_real(n_spin,win_tri_occ,win_tri_unocc,&
          lvl_tricoeff_recip_k_r_real, lvl_tricoeff_recip_q_r_real, &
          lvl_tricoeff_recip_kp_r_real,lvl_tricoeff_recip_qp_r_real, &
          kq_pair_next,kq_pair_new,next_run)
      ! load four KS_eigenvectors for the next term
      call saccess_4_evs_cpt2_real(n_spin,win_ev_occ,win_ev_unocc,&
          KS_eigenvector_k_r_real,KS_eigenvector_q_r_real,&
          KS_eigenvector_kp_r_real,KS_eigenvector_qp_r_real,&
          kq_pair_next,kq_pair_new,next_run)
      ! load two coulomb matrices for the next term
      call saccess_2_coulomb_cpt2_real(n_spin,win_coul,&
          i_qk_point_next,i_qk_point,coulomb_matr_qk_r_real,&
          i_qpk_point_next,i_qpk_point,coulomb_matr_qpk_r_real)

      ! if the current task is fake, jump to the next term.
      if (.not. curr_run) then
        kq_pair_old = 0
        cycle
      end if

      ! determine if the symmetry in ERIs can be used or not.
      if (i_k_point.eq.i_kp_point .and. i_q_point.eq.i_qp_point) then
          eri_symm = .true.
      else
          eri_symm = .false.
      end if

      ! Now do PT2 calculation
      ! FIXME :: at this moment, it supports close-shell cases only.
      pt2_term_os = 0.0d0
      pt2_term_ss = 0.0d0
      do i_spin = 1, n_spin, 1
      do i_spin_2 = 1, n_spin, 1
      do a_state = n_low_state, n_homo_k(i_k_point, i_spin), 1

        ! explore the symmetry in ERIs
        if (i_spin.eq.i_spin_2 .and. eri_symm) then
            b_state_start = a_state
        else
            b_state_start = n_low_state
        end if

        do b_state = b_state_start, n_homo_k(i_kp_point,i_spin_2), 1

          !aux_eri_kqkpqp_r_real = 0.0d0
          ptr_aux_eri_in_1 => lvl_tricoeff_kq_r_real(lbb_row:ubb_row,lpb_col:upb_col,a_state,i_spin)
          ptr_aux_eri_in_2 => lvl_tricoeff_kpqp_r_real(lbb_row:ubb_row,lpb_col:upb_col,b_state,i_spin_2)
          ptr_aux_eri_out => aux_eri_kqkpqp_r_real(lpb_row:upb_row,lpb_col:upb_col)
          call pdgemm('T', 'N', n_unocc_max, n_unocc_max, &
             n_basbas, 1.0d0, &
             ptr_aux_eri_in_1,1,1, pb2desc, &
             ptr_aux_eri_in_2,1,1, pb2desc, 0.d0, &
             ptr_aux_eri_out, 1,1, pp2desc &
             )
          aux_eri_tmp_d = 0.0d0
          aux_eri_tmp_d(lpb_row:upb_row,lpb_col:upb_col) =  &
              aux_eri_kqkpqp_r_real(lpb_row:upb_row,lpb_col:upb_col)
          call sync_vector(aux_eri_tmp_d,size(aux_eri_tmp_d),comm_blacs_row_cpt2)

          if (i_spin .eq. i_spin_2) then
            if (eri_symm) then
              call pdtran( n_unocc_max, n_unocc_max, &
                   1.0d0, aux_eri_kqkpqp_r_real, 1, 1, pp2desc, &
                   0.0d0, aux_eri_kqpkpq_r_real, 1, 1, pp2desc )
            else
              !aux_eri_kqpkpq_r_real = 0.0d0
              ptr_aux_eri_in_1 => lvl_tricoeff_kpq_r_real(lbb_row:ubb_row,lpb_col:upb_col,b_state,i_spin)
              ptr_aux_eri_in_2 => lvl_tricoeff_kqp_r_real(lbb_row:ubb_row,lpb_col:upb_col,a_state,i_spin)
              ptr_aux_eri_out => aux_eri_kqpkpq_r_real(lpb_row:upb_row,lpb_col:upb_col)
              call pdgemm('T', 'N', n_unocc_max, n_unocc_max, &
                     n_basbas, 1.0d0, &
                     ptr_aux_eri_in_1,1,1, pb2desc, &
                     ptr_aux_eri_in_2,1,1, pb2desc, 0.d0, &
                     ptr_aux_eri_out, 1,1, pp2desc &
                     )
            end if
            aux_eri_tmp_x = 0.0d0
            aux_eri_tmp_x(lpb_row:upb_row,lpb_col:upb_col) =  &
                aux_eri_kqpkpq_r_real(lpb_row:upb_row,lpb_col:upb_col)
            call sync_vector(aux_eri_tmp_x,size(aux_eri_tmp_x),comm_blacs_row_cpt2)
          end if

          do m_state = n_lumo_k(i_qp_point,i_spin_2), n_states_k(i_qp_point),1

            !! parallization with respect to m_state with i_qp_point in column
            !if (.not. (m_state .ge. lpb_col .and. m_state .le. upb_col)) cycle
            ! modify for shifted index
            m_state_shifted = m_state - n_lumo_min + 1
            if (.not. (m_state_shifted .ge. lpb_col .and. m_state_shifted .le. upb_col)) cycle
            ! end of modification

            do n_state = n_lumo_k(i_q_point,i_spin), n_states_k(i_q_point), 1

              !if (myid_row_cpt2 .ne. mod(n_state-n_lumo_min,n_tasks_row_cpt2)) cycle
              ! modify for shifted index
              n_state_shifted = n_state - n_lumo_min + 1
              if (myid_row_cpt2 .ne. mod(n_state_shifted-1,n_tasks_row_cpt2)) cycle
              ! end of modification

              e_diff = KS_eigenvalue(a_state,i_spin,i_k_point)  &
                     + KS_eigenvalue(b_state,i_spin_2,i_kp_point) &
                     - KS_eigenvalue(n_state,i_spin,i_q_point)  &
                     - KS_eigenvalue(m_state,i_spin_2,i_qp_point)
              if (abs(e_diff).lt.1e-6)then
                 write(use_unit,'(10X,A)') &
                      "****************************************"
                 write(use_unit,'(10X,2A)') "| Warning :", &
                   " too close to degeneracy"
                 write(use_unit,'(10X,A)') &
                      "****************************************"
              endif

              E_mp2_a = aux_eri_tmp_d(n_state_shifted,m_state_shifted)

              if (n_spin .eq. 1) then
                E_mp2_b = aux_eri_tmp_x(n_state_shifted,m_state_shifted)
                E_mp2_tmp_os = E_mp2_a*E_mp2_a/e_diff
                E_mp2_tmp_ss = E_mp2_a*(E_mp2_a - E_mp2_b)/e_diff
                !pt2_term_os = pt2_term_os + E_mp2_tmp_os
                !pt2_term_ss = pt2_term_ss + E_mp2_tmp_ss
              else
                if (i_spin .eq. i_spin_2) then
                  E_mp2_b = aux_eri_tmp_x(n_state_shifted,m_state_shifted)
                  E_mp2_tmp_ss = E_mp2_a*(E_mp2_a - E_mp2_b)/e_diff
                  E_mp2_tmp_os = 0.0d0
                  !pt2_term_ss = pt2_term_ss + E_mp2_tmp_ss
                else
                  E_mp2_tmp_os = E_mp2_a*E_mp2_a/e_diff
                  E_mp2_tmp_ss = 0.0d0
                  !pt2_term_os = pt2_term_os + E_mp2_tmp_os
                end if
              end if

              ! utilize the symmetry of ERIs
              if (a_state.ne.b_state .and. eri_symm) then
                  E_mp2_tmp_os = E_mp2_tmp_os*2.0d0
                  E_mp2_tmp_ss = E_mp2_tmp_ss*2.0d0
              end if

              pt2_term_os = pt2_term_os + E_mp2_tmp_os
              pt2_term_ss = pt2_term_ss + E_mp2_tmp_ss

            enddo ! n_state for the virtual states of i_kp_point
          enddo ! m_state for the virtual states of i_qp_point
        enddo ! b_state for the occupied states of i_q_point
      enddo ! a_state for the occupied states of i_k_point
      enddo ! i_spin
      enddo ! i_spin_2

      pt2_c_energy_local_os = pt2_c_energy_local_os + pt2_term_os * k3_weight

      pt2_c_energy_local_ss = pt2_c_energy_local_ss + pt2_term_ss * k3_weight

      kq_pair_old = kq_pair_new

    end do !do i_pair = 1, max_pair_local, 1

    ! stop a fake communication
    call ssync_trico_cpt2_real(win_tri_occ,win_tri_unocc,lvl_tricoeff_recip_k_r_real, &
         lvl_tricoeff_recip_kp_r_real, &
         kq_pair_new,kq_pair_old,.false.)
    call ssync_ev_cpt2_real(win_ev_occ,win_ev_unocc,KS_eigenvector_k_r_real,&
         KS_eigenvector_kp_r_real, &
         kq_pair_new,kq_pair_old,.false.)
    call ssync_coulomb_cpt2(win_coul)

  end subroutine cpt2_kq_list_real

  subroutine compute_m_kq_real(n_low_state, i_spin, iop, &
       i_k_point, i_q_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
       KS_eigenvector_k, KS_eigenvector_q, &
       lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &  
       coulomb_matr, lvl_tricoeff_kq)
  
    implicit none
  
    integer, intent(in) :: i_spin, i_k_point, i_q_point, iop, n_low_state
  
    integer, intent(in) :: n_lumo_k(n_k_points,n_spin), n_homo_k(n_k_points,n_spin)
    real(kind=8), dimension(n_states,n_spin,n_k_points), &
        intent(in) :: KS_eigenvalue
    real(kind=8), dimension(n_basis,n_low_state:n_homo_max,n_spin), &
        intent(in) :: KS_eigenvector_k
    real(kind=8), dimension(n_basis,lpb_col:upb_col,n_spin), &
        intent(in) :: KS_eigenvector_q
    real(kind=8), dimension(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin), &
        intent(in) :: lvl_tricoeff_recip_k_r
    real(kind=8), dimension(lbb_row:ubb_row,max_n_basis_sp,lpb_col:upb_col,n_spin), &
        intent(in) :: lvl_tricoeff_recip_q_r  
    real(kind=8), dimension(lbb_row:ubb_row,lbb_col:ubb_col), &
       intent(in) :: coulomb_matr 
  
    !rea8, intent(in), dimension(n_states,n_spin,n_k_points) :: occ_numbers

    real(kind=8), intent(inout), &
        dimension(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin), &
        target :: lvl_tricoeff_kq
  
    ! counter
    integer:: lb, ub, i_state, i_state_2, i_state_3, i_atom, i_species
    integer:: lbb_atom, ubb_atom, ind, brange
    real(kind=8), pointer :: ptr_lvl_tricoeff_kq_in(:,:,:)
  
!    call perfon('polfreq')

    do i_state=n_low_state, n_homo_k(i_k_point,i_spin)
      !do i_state_2=lpb_col, upb_col
      ! modify for shifted index
      do i_state_3=lpb_col, upb_col
        i_state_2 = n_lumo_min-1+i_state_3
      ! end of modification
        do i_atom=basbas_atom(lbb_row), basbas_atom(ubb_row)
          ! range in basis dimension
          lb=lb_atom(i_atom)
          ub=ub_atom(i_atom)
          brange=ub-lb+1
          ! range in basbas dimension
          i_species = species(i_atom)
          lbb_atom=max(lbb_row,atom2basbas_off(i_atom)+1)
          ubb_atom=min(ubb_row,atom2basbas_off(i_atom)+sp2n_basbas_sp(i_species))

          lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_3,i_state,i_spin)=0.
          do ind=1,brange
             lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_3,i_state,i_spin)=&
                  lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_3,i_state,i_spin) + &
                  KS_eigenvector_q(lb+ind-1,i_state_3,i_spin) * &
                  lvl_tricoeff_recip_k_r(lbb_atom:ubb_atom,ind,i_state,i_spin) + &
                  KS_eigenvector_k(lb+ind-1,i_state,i_spin) * &
                  lvl_tricoeff_recip_q_r(lbb_atom:ubb_atom,ind,i_state_3,i_spin)
          end do
        end do
      end do
    end do

    ptr_lvl_tricoeff_kq_in => lvl_tricoeff_kq(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,i_spin)
    call get_v_multi_ovlp3fn_real_blacs &
                 ( n_low_state, n_homo_max, n_homo_k(i_k_point,i_spin), coulomb_matr, &
                   ptr_lvl_tricoeff_kq_in, iop)

  end subroutine compute_m_kq_real

  subroutine allocate_cpt2_kspace_real(n_low_state, n_homo_max)

    implicit none
    integer, intent(in) :: n_low_state, n_homo_max

    character(*), parameter :: func = 'evaluate_cpt2_energy_kspace_blacs.f90'
    integer :: info
    character*150 :: info_str

    ! create the windows for KS_eigenvectors, lvl_tricoeffs, and coulomb matrix
    call sinit_access_ev_real(n_k_points,n_ks_points_task, &
                              KS_eigenvectors_occ_real,win_ev_occ, &
                              KS_eigenvectors_unocc_real, win_ev_unocc)

    call sinit_access_trico_real(n_k_points,n_ks_points_task,&
                                   lvl_tricoeff_occ_real,win_tri_occ, &
                                   lvl_tricoeff_unocc_real,win_tri_unocc)

    call sinit_access_coulomb_real(n_k_points,n_kq_points_task,coulomb_matr_blacs_real,win_coul)

    ! allocate relevant arrays

    if (n_tasks_kq_cpt2.gt.1) then
      allocate(KS_eigenvector_k_real(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_k', func)
      KS_eigenvector_k_r_real => KS_eigenvector_k_real
      allocate(KS_eigenvector_kp_real(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_kp', func)
      KS_eigenvector_kp_r_real => KS_eigenvector_kp_real

      allocate(lvl_tricoeff_k_real(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_k', func)
      lvl_tricoeff_recip_k_r_real => lvl_tricoeff_k_real
      allocate(lvl_tricoeff_kp_real(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_kp', func)
      lvl_tricoeff_recip_kp_r_real => lvl_tricoeff_kp_real

      allocate(KS_eigenvector_q_real(n_basis,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_q', func)
      KS_eigenvector_q_r_real => KS_eigenvector_q_real
      allocate(KS_eigenvector_qp_real(n_basis,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_qp', func)
      KS_eigenvector_qp_r_real => KS_eigenvector_qp_real

      allocate(lvl_tricoeff_q_real(lbb_row:ubb_row,max_n_basis_sp,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_q', func)
      lvl_tricoeff_recip_q_r_real => lvl_tricoeff_q_real
      allocate(lvl_tricoeff_qp_real(lbb_row:ubb_row,max_n_basis_sp,lpb_col:upb_col,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_qp', func)
      lvl_tricoeff_recip_qp_r_real => lvl_tricoeff_qp_real

      allocate(coulomb_matr_qk_real(lbb_row:ubb_row,lbb_col:ubb_col),stat=info) 
      call check_allocation(info, 'coulomb_matr_qk', func)
      coulomb_matr_qk_r_real => coulomb_matr_qk_real
      allocate(coulomb_matr_qpk_real(lbb_row:ubb_row,lbb_col:ubb_col),stat=info) 
      call check_allocation(info, 'coulomb_matr_qpk', func)
      coulomb_matr_qpk_r_real => coulomb_matr_qpk_real
    else
      if (myid_col_cpt2.ne.0) then
        allocate(KS_eigenvector_k_real(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'KS_eigenvector_k', func)
        KS_eigenvector_k_r_real => KS_eigenvector_k_real
        allocate(KS_eigenvector_kp_real(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'KS_eigenvector_kp', func)
        KS_eigenvector_kp_r_real => KS_eigenvector_kp_real

        allocate(lvl_tricoeff_k_real(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'lvl_tricoeff_k', func)
        lvl_tricoeff_recip_k_r_real => lvl_tricoeff_k_real
        allocate(lvl_tricoeff_kp_real(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'lvl_tricoeff_kp', func)
        lvl_tricoeff_recip_kp_r_real => lvl_tricoeff_kp_real
      end if
    end if

    allocate(lvl_tricoeff_kq_real(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kq', func)
    lvl_tricoeff_kq_r_real=>lvl_tricoeff_kq_real
    allocate(lvl_tricoeff_kqp_real(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kqp', func)
    lvl_tricoeff_kqp_r_real=>lvl_tricoeff_kqp_real
    allocate(lvl_tricoeff_kpq_real(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kpq', func)
    lvl_tricoeff_kpq_r_real=>lvl_tricoeff_kpq_real
    allocate(lvl_tricoeff_kpqp_real(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kpqp', func)
    lvl_tricoeff_kpqp_r_real=>lvl_tricoeff_kpqp_real

    allocate(aux_eri_kqkpqp_real(lpb_row:upb_row,lpb_col:upb_col),stat=info)
    call check_allocation(info, 'aux_eri_kq', func)
    aux_eri_kqkpqp_r_real=>aux_eri_kqkpqp_real
    allocate(aux_eri_kqpkpq_real(lpb_row:upb_row,lpb_col:upb_col),stat=info)
    call check_allocation(info, 'aux_eri_kq', func)
    aux_eri_kqpkpq_r_real=>aux_eri_kqpkpq_real

  end subroutine allocate_cpt2_kspace_real

  subroutine deallocate_cpt2_kspace_real()

    if (.not. gamma_only_cpt2) then
       call sfinalize_access_trico_cpt2(win_tri_occ,win_tri_unocc)
       call sfinalize_access_ev_cpt2(win_ev_occ,win_ev_unocc)
       call sfinalize_access_coulomb_cpt2(win_coul)
   end if

    if (allocated(KS_eigenvector_k_real)) deallocate(KS_eigenvector_k_real)
    if (allocated(KS_eigenvector_q_real)) deallocate(KS_eigenvector_q_real)
    if (allocated(KS_eigenvector_kp_real)) deallocate(KS_eigenvector_kp_real)
    if (allocated(KS_eigenvector_qp_real)) deallocate(KS_eigenvector_qp_real)

    if (allocated(lvl_tricoeff_k_real)) deallocate(lvl_tricoeff_k_real)
    if (allocated(lvl_tricoeff_q_real)) deallocate(lvl_tricoeff_q_real)
    if (allocated(lvl_tricoeff_kp_real)) deallocate(lvl_tricoeff_kp_real)
    if (allocated(lvl_tricoeff_qp_real)) deallocate(lvl_tricoeff_qp_real)

    if (allocated(coulomb_matr_qk_real)) deallocate(coulomb_matr_qk_real)
    if (allocated(coulomb_matr_qpk_real)) deallocate(coulomb_matr_qpk_real)

    if (allocated(lvl_tricoeff_kq_real)) deallocate(lvl_tricoeff_kq_real)
    if (allocated(lvl_tricoeff_kqp_real)) deallocate(lvl_tricoeff_kqp_real)
    if (allocated(lvl_tricoeff_kpq_real)) deallocate(lvl_tricoeff_kpq_real)
    if (allocated(lvl_tricoeff_kpqp_real)) deallocate(lvl_tricoeff_kpqp_real)

    if (allocated(aux_eri_kqkpqp_real)) deallocate(aux_eri_kqkpqp_real)
    if (allocated(aux_eri_kqpkpq_real)) deallocate(aux_eri_kqpkpq_real)

  end subroutine deallocate_cpt2_kspace_real

  subroutine evaluate_memory_consumption_cpt2 (n_low_state,n_homo_max,n_lumo_min,iop)
    implicit none
    integer, intent(in) :: n_low_state, n_homo_max, n_lumo_min, iop
    ! iop = 1 : multiple k-grids simulation in COMPLEX
    ! iop = 2 : multiple k-grids simulation in REAL
    ! iop = 3 : gamma-only in REAL

    real(kind=8)  :: memory_consumption(3,10)
    real(kind=8) :: dcsize, munit, tmp_memory
    character*150 :: info_str
    
    munit=2.0**30

    if (iop.eq.2 .or. iop.eq.3) then
      dcsize=8.0d0
      write(info_str, &
          "(2X,'It is a canonical PT2 calculation with',I5,' k-grid(s) in REAL.')") &
          n_k_points
      call localorb_info(info_str)
    else if (iop.eq.1) then
      dcsize=16.0d0
      write(info_str, &
          "(2X,'It is a canonical PT2 calculation with',I5,' k-grids in COMPLEX.')") &
          n_k_points
      call localorb_info(info_str)
    end if

    memory_consumption = 0.0d0

    ! for lvl_tri_coeff
    memory_consumption(1,1) = real(int(n_k_points,8)*int(n_spin,8)*int(n_states,8)*&
                               int(max_n_basis_sp,8)*int(n_basbas,8),8)/munit*dcsize
    memory_consumption(2,1) = memory_consumption(1,1)/real(n_tasks,8)
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,1)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,1)
    ! for coulomb_matr 
    memory_consumption(1,2) = real(int(n_k_points,8)*int(n_basbas,8)*&
                               int(n_basbas,8),8)/munit*dcsize
    memory_consumption(2,2) = memory_consumption(1,2)/n_tasks
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,2)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,2)
    ! for KS_eigenvectors and KS_eigenvalue
    memory_consumption(1,3) = n_spin*(n_basis+1)*n_states/munit*n_k_points*dcsize
    memory_consumption(2,3) = n_spin*(n_basis/n_tasks+1)*n_states/munit*n_k_points*dcsize
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,3)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,3)

    ! for KS_eigenvectors
    if (iop.eq.1 .or. iop.eq.2) then
      tmp_memory=n_spin*n_basis*(n_homo_max-n_low_state+1)
      if (n_tasks_kq_cpt2.gt.1) then
        tmp_memory=tmp_memory + &
            int(n_spin*n_basis,8)*int((n_states-n_lumo_min+1)/n_tasks_col_cpt2+1,8)
      end if
      memory_consumption(2,4) = 2.0d0*tmp_memory*dcsize/munit
    else if (iop.eq.3) then
      tmp_memory=n_spin*n_basis*(n_homo_max-n_low_state+1)
      memory_consumption(2,4) = 1.0d0*tmp_memory*dcsize/munit
    end if
    memory_consumption(1,4) = memory_consumption(2,4)*n_tasks
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,4)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,4)

    ! for lvl_tricoeff_recip_k
    if (iop.eq.1 .or. iop.eq.2) then
      tmp_memory=&
          n_spin*max_n_basis_sp*(n_homo_max-n_low_state+1)*int(n_basbas/n_tasks_row_cpt2+1)
      if (n_tasks_kq_cpt2.gt.1) then
        tmp_memory=tmp_memory+&
          n_spin*max_n_basis_sp*int((n_states-n_lumo_min+1)/n_tasks_col_cpt2+1)*&
          int(n_basbas/n_tasks_row_cpt2+1)
      end if
      memory_consumption(2,5) = 2.0d0*tmp_memory*dcsize/munit
    else if (iop.eq.3) then
      tmp_memory=&
          n_spin*max_n_basis_sp*(n_homo_max-n_low_state+1)*int(n_basbas/n_tasks_row_cpt2+1)
      memory_consumption(2,5) = 1.0d0*tmp_memory*dcsize/munit
    end if
    memory_consumption(1,5) = memory_consumption(2,5)*n_tasks
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,5)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,5)

    ! for coulomb_matr_k
    tmp_memory=0.0d0
    if (n_tasks_kq_cpt2.gt.1) then
        tmp_memory=int(n_basbas/n_tasks_row_cpt2+1)**2
    end if
    if (iop.eq.1 .or. iop.eq.2) then
      memory_consumption(2,7) = 2.0d0*tmp_memory*dcsize/munit
    else if (iop.eq.3) then
      memory_consumption(2,7) = 1.0d0*tmp_memory*dcsize/munit
    end if
    memory_consumption(1,7) = memory_consumption(2,7)*n_tasks
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,7)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,7)

    ! for lvl_full_kq
    tmp_memory=n_spin*(n_homo_max-n_low_state+1)*&
               int(n_basbas/n_tasks_row_cpt2+1)*&
               int((n_states-n_lumo_min+1)/n_tasks_col_cpt2+1)
    if (iop.eq.1 .or. iop.eq.2) then
      memory_consumption(2,6) = 4.0d0*tmp_memory*dcsize/munit
    else if (iop.eq.3) then
      memory_consumption(2,6) = 1.0d0*tmp_memory*dcsize/munit
    end if
    memory_consumption(1,6) = memory_consumption(2,6)*n_tasks
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,6)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,6)

    ! for aux_eri
    tmp_memory=int((n_states-n_lumo_min+1)/n_tasks_row_cpt2+1)*&
        int((n_states-n_lumo_min+1)/n_tasks_col_cpt2+1)
    tmp_memory= tmp_memory + tmp_memory * n_tasks_row_cpt2
    memory_consumption(2,8) = 2.0d0*tmp_memory*dcsize/munit
    memory_consumption(1,8) = memory_consumption(2,8)*n_tasks
    memory_consumption(1,10) = memory_consumption(1,10)+memory_consumption(1,8)
    memory_consumption(2,10) = memory_consumption(2,10)+memory_consumption(2,8)

    ! Print out the memory needed in this PT2 calculation
    write(info_str, *) '   -----------------------------------------',&
        '-------------------------------------------'
    call localorb_info(info_str)
    write(info_str, "(4X,A40,' : ',A16,',',A16,' (in Gb)')") &
        'Detailed memory counting                ','per-task','global'
    call localorb_info(info_str)
    write(info_str, *) '   -----------------------------------------',&
        '-------------------------------------------'
    call localorb_info(info_str) 
    write(info_str, "(4X,A40,' : ',F16.3,',',F16.3)") &
        '|- Minimal memory request               ', &
        memory_consumption(2,10), memory_consumption(1,10)
    call localorb_info(info_str) 
    write(info_str, *) '   -----------------------------------------',&
        '-------------------------------------------'
    call localorb_info(info_str) 

    write(info_str, "(4X,'|-Data distributed among n_tasks    (', I5, ') --------')") &
        n_tasks
    call localorb_info(info_str)
    write(info_str, "(4X, A40,' : ',F16.3,',',F16.3)") &
        '|- LVL 3-center overlap                 ', &
        memory_consumption(2,1), memory_consumption(1,1)
    call localorb_info(info_str) 
    write(info_str, "(4X,A40,' : ',F16.3,',',F16.3)") &
        '|- Coulomb matrix                       ', &
        memory_consumption(2,2), memory_consumption(1,2)
    call localorb_info(info_str) 
    write(info_str, "(4X,A40,' : ',F16.3,',',F16.3)") &
        '|- KS eigenvectors+eigenvalues          ', &
        memory_consumption(2,3), memory_consumption(1,3)
    call localorb_info(info_str) 

    write(info_str, *) '   -----------------------------------------',&
        '-------------------------------------------'
    call localorb_info(info_str) 
    write(info_str, "(4X,'|-Data distributed among n_tasks_bl (', I5, ') --------')") &
        n_tasks_bl_cpt2
    call localorb_info(info_str)

    if (n_tasks_kq_cpt2.gt.1) then
      write(info_str, "(4X,A40,' : ',F16.3,',',F16.3)") &
        '|- Coulomb matrix                       ', &
        memory_consumption(2,7), memory_consumption(1,7)
      call localorb_info(info_str) 
    end if
    write(info_str, "(4X, A40,' : ',F16.3,',',F16.3)") &
        '|- LVL 3-center overlap                 ', &
        memory_consumption(2,5), memory_consumption(1,5)
    call localorb_info(info_str) 
    write(info_str, "(4X,A40,' : ',F16.3,',',F16.3)") &
        '|- KS eigenvectors                      ', &
        memory_consumption(2,4), memory_consumption(1,4)
    call localorb_info(info_str) 

    write(info_str, "(4X,A40,' : ',F16.3,',',F16.3)") &
        '|- LVL-full 3-center overlap            ', &
        memory_consumption(2,6), memory_consumption(1,6)
    call localorb_info(info_str) 

    write(info_str, "(4X,A40,' : ',F16.3,',',F16.3)") &
        '|- Virtual-orbital part of ERIs         ', &
        memory_consumption(2,8), memory_consumption(1,8)
    call localorb_info(info_str) 

    write(info_str, *) '   -----------------------------------------',&
        '-------------------------------------------'
    call localorb_info(info_str) 
    write(info_str, *) ''
    call localorb_info(info_str)

  end subroutine evaluate_memory_consumption_cpt2

  subroutine allocate_cpt2_kspace_gamma_only(n_low_state, n_homo_max)

    implicit none
    integer, intent(in) :: n_low_state, n_homo_max

    character(*), parameter :: func = 'evaluate_cpt2_energy_kspace_blacs.f90'
    integer :: info
    character*150 :: info_str

    ! allocate relevant arrays

    if (myid_col_cpt2.eq.0) then
        KS_eigenvector_k_r_real => KS_eigenvectors_occ_real(:,:,:,1)
        lvl_tricoeff_recip_k_r_real => lvl_tricoeff_occ_real(:,:,:,:,1)
    else
        allocate(KS_eigenvector_k_real(n_basis,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'KS_eigenvector_k', func)
        KS_eigenvector_k_r_real => KS_eigenvector_k_real
        allocate(lvl_tricoeff_k_real(lbb_row:ubb_row,max_n_basis_sp,n_low_state:n_homo_max,n_spin),stat=info) 
        call check_allocation(info, 'lvl_tricoeff_k', func)
        lvl_tricoeff_recip_k_r_real => lvl_tricoeff_k_real
    end if
    KS_eigenvector_q_r_real => KS_eigenvectors_unocc_real(:,:,:,1)
    lvl_tricoeff_recip_q_r_real => lvl_tricoeff_unocc_real(:,:,:,:,1)
    coulomb_matr_qk_r_real => coulomb_matr_blacs_real(:,:,1)

    allocate(lvl_tricoeff_kq_real(lbb_row:ubb_row,lpb_col:upb_col,n_low_state:n_homo_max,n_spin),stat=info)
    call check_allocation(info, 'lvl_tricoeff_kq', func)
    lvl_tricoeff_kq_r_real=>lvl_tricoeff_kq_real

    allocate(aux_eri_kqkpqp_real(lpb_row:upb_row,lpb_col:upb_col),stat=info)
    call check_allocation(info, 'aux_eri_kq', func)
    aux_eri_kqkpqp_r_real=>aux_eri_kqkpqp_real
    allocate(aux_eri_kqpkpq_real(lpb_row:upb_row,lpb_col:upb_col),stat=info)
    call check_allocation(info, 'aux_eri_kq', func)
    aux_eri_kqpkpq_r_real=>aux_eri_kqpkpq_real

  end subroutine allocate_cpt2_kspace_gamma_only

  subroutine cpt2_kq_list_gamma_only(paircount, kq_pair, total_weight, &
                   n_low_state, n_homo_max, &
                   n_homo_k, n_lumo_k, &
                   pt2_c_energy_local, pt2_c_energy_local_ss, &
                   pt2_c_energy_local_os)

    integer, intent(in) :: paircount, n_low_state, n_homo_max
    integer, intent(in), dimension(6,paircount) :: kq_pair
    real(kind=8),  intent(in), dimension(paircount) :: total_weight
    
    integer, dimension(n_k_points,n_spin),intent(in) :: n_homo_k
    integer, dimension(n_k_points,n_spin),intent(in) :: n_lumo_k

    real(kind=8), intent(inout) :: pt2_c_energy_local, &
        pt2_c_energy_local_os, pt2_c_energy_local_ss

    ! counter
    integer :: i_spin, i_spin_2, i_pair
    integer :: a_state, b_state, n_state, m_state
    integer :: m_state_shifted, n_state_shifted
    integer :: n_state_p, i_task_row

    real(kind=8)  :: k3_weight

    ! for k_point counter
    integer :: i_k_point, i_q_point, i_kp_point, i_qp_point
    integer :: i_qk_point,i_qpk_point
    integer :: i_k_point_old,i_q_point_old,i_kp_point_old,i_qp_point_old
    integer :: i_qk_point_old,i_qpk_point_old
    integer :: i_k_point_next,i_q_point_next,i_kp_point_next,i_qp_point_next
    integer :: i_qk_point_next,i_qpk_point_next

    integer :: b_state_start

    ! record the PT2 correlation terms locally
    real(kind=8)  :: pt2_term, pt2_term_os, pt2_term_ss, e_diff
    real(kind=8)  :: E_mp2_tmp_os, E_mp2_tmp_ss
    real(kind=8)  :: E_mp2_a, E_mp2_b

    integer, dimension(6) :: kq_pair_new, kq_pair_old, kq_pair_next
    logical :: curr_run, next_run

    real(kind=8), dimension(n_unocc_max,lpb_col:upb_col), &
        target :: aux_eri_tmp_d, aux_eri_tmp_x
    real(kind=8), dimension(:,:), &
        pointer :: aux_eri_tmp_d_r, aux_eri_tmp_x_r
    real(kind=8), pointer :: ptr_aux_eri_in_1(:,:)
    real(kind=8), pointer :: ptr_aux_eri_in_2(:,:)
    real(kind=8), pointer :: ptr_aux_eri_out(:,:)

    real(kind=8)   :: tmp_pt2_os, tmp_pt2_ss

    kq_pair_old = 0
    i_pair = 1
    ! loading the k-point pattern for the this term
    if (((i_pair.eq.max_pair_local) .and. &
        (paircount_local.lt.max_pair_local)) .or. &
        .not. kblacs_member) then
      kq_pair_new = 0
      k3_weight = 0.0d0
      curr_run = .false.
    else
      kq_pair_new = kq_pair(:,lbpair+i_pair-1)
      k3_weight = total_weight(lbpair+i_pair-1)
      curr_run = .true.
    end if

    i_k_point = kq_pair_new(1)
    i_q_point = kq_pair_new(2)
    i_kp_point = kq_pair_new(3)
    i_qp_point = kq_pair_new(4)
    i_qk_point = kq_pair_new(5)
    i_qpk_point = kq_pair_new(6)

    call get_timestamps(time_distribution(3),time_distribution(4))
    ! ================================================
    ! sync lvl_tricoeff_recip_k_r_real
    ! ================================================
    call ssync_trico_cpt2_gamma_only(lvl_tricoeff_recip_k_r_real,curr_run)
    ! ================================================
    ! sync KS_eigenvector_k_r_real
    ! ================================================
    call ssync_ev_cpt2_gamma_only(KS_eigenvector_k_r_real,curr_run)

    !============
    ! Analysis the timing for the memory distribution.
    !============
    call get_timestamps(time_distribution(5), time_distribution(6))
    time_distribution(3) = time_distribution(5) - time_distribution(3)
    time_distribution(4) = time_distribution(6) - time_distribution(4)
    time_distribution(1) = time_distribution(1) + time_distribution(3)
    time_distribution(2) = time_distribution(2) + time_distribution(4)

    !============
    ! Calculate the RI-V coefficient
    !============
    if (curr_run) then
      do i_spin = 1, n_spin
        call compute_m_kq_real(n_low_state, i_spin, 1, i_k_point,&
             i_q_point, n_lumo_k, n_homo_k, KS_eigenvalue, &
             KS_eigenvector_k_r_real, KS_eigenvector_q_r_real, &
             lvl_tricoeff_recip_k_r_real, lvl_tricoeff_recip_q_r_real, &  
             coulomb_matr_qk_r_real, lvl_tricoeff_kq_r_real)
      enddo ! end loop over i_spin
    end if

    call get_timestamps(time_pt2(7), time_pt2(8))

    ! if the current task is fake, jump to the next term after comunication.
    if (.not. curr_run) then
      kq_pair_old = 0
      return
    end if

    !============
    ! Now do PT2 calculation
    !============
    pt2_term_os = 0.0d0
    pt2_term_ss = 0.0d0
    do i_spin = 1, n_spin, 1
    do i_spin_2 = i_spin, n_spin, 1  ! only terms that are not covered by symmetry
    do a_state = n_low_state, n_homo_k(i_k_point, i_spin), 1

      if (i_spin.eq.i_spin_2) then
          ! Symmetry in ERIs may be exploited
          b_state_start = a_state
      else
          b_state_start = n_low_state
      end if

      do b_state = b_state_start, n_homo_k(i_kp_point,i_spin_2), 1

        !aux_eri_kqkpqp_r_real = 0.0d0
        ptr_aux_eri_in_1 => lvl_tricoeff_kq_r_real(lbb_row:ubb_row,lpb_col:upb_col,a_state,i_spin)
        ptr_aux_eri_in_2 => lvl_tricoeff_kq_r_real(lbb_row:ubb_row,lpb_col:upb_col,b_state,i_spin_2)
        ptr_aux_eri_out => aux_eri_kqkpqp_r_real(lpb_row:upb_row,lpb_col:upb_col)
        call pdgemm('T', 'N', n_unocc_max, n_unocc_max, &
           n_basbas, 1.0d0, &
           ptr_aux_eri_in_1,1,1, pb2desc, &
           ptr_aux_eri_in_2,1,1, pb2desc, 0.d0, &
           ptr_aux_eri_out, 1,1, pp2desc &
           )
        aux_eri_tmp_d = 0.0d0
        aux_eri_tmp_d(lpb_row:upb_row,lpb_col:upb_col) =  &
            aux_eri_kqkpqp_r_real(lpb_row:upb_row,lpb_col:upb_col)
        call sync_vector(aux_eri_tmp_d,size(aux_eri_tmp_d),comm_blacs_row_cpt2)

        if (i_spin .eq. i_spin_2) then
            call pdtran( n_unocc_max, n_unocc_max, &
                 1.0d0, aux_eri_kqkpqp_r_real, 1, 1, pp2desc, &
                 0.0d0, aux_eri_kqpkpq_r_real, 1, 1, pp2desc )
            aux_eri_tmp_x = 0.0d0
            aux_eri_tmp_x(lpb_row:upb_row,lpb_col:upb_col) =  &
                aux_eri_kqpkpq_r_real(lpb_row:upb_row,lpb_col:upb_col)
            call sync_vector(aux_eri_tmp_x,size(aux_eri_tmp_x),comm_blacs_row_cpt2)
        end if

        do m_state = n_lumo_k(i_qp_point,i_spin_2), n_states_k(i_qp_point),1
          m_state_shifted = m_state - n_lumo_min + 1
          ! parallization with respect to m_state with i_qp_point in column
          if (.not. (m_state_shifted .ge. lpb_col .and. m_state_shifted .le. upb_col)) cycle

          do n_state = n_lumo_k(i_q_point,i_spin), n_states_k(i_q_point), 1
            n_state_shifted = n_state - n_lumo_min + 1
            if (myid_row_cpt2 .ne. mod(n_state_shifted-1,n_tasks_row_cpt2)) cycle

            e_diff = KS_eigenvalue(a_state,i_spin,i_k_point)  &
                   + KS_eigenvalue(b_state,i_spin_2,i_kp_point) &
                   - KS_eigenvalue(n_state,i_spin,i_q_point)  &
                   - KS_eigenvalue(m_state,i_spin_2,i_qp_point)
            if (abs(e_diff).lt.1e-6)then
               write(use_unit,'(10X,A)') &
                    "****************************************"
               write(use_unit,'(10X,2A)') "| Warning :", &
                 " too close to degeneracy"
               write(use_unit,'(10X,A)') &
                    "****************************************"
            endif

            E_mp2_a = aux_eri_tmp_d(n_state_shifted,m_state_shifted)
            if (n_spin .eq. 1) then
              E_mp2_b = aux_eri_tmp_x(n_state_shifted,m_state_shifted)
              E_mp2_tmp_os = E_mp2_a*E_mp2_a/e_diff
              E_mp2_tmp_ss = E_mp2_a*(E_mp2_a - E_mp2_b)/e_diff
            else
              if (i_spin .eq. i_spin_2) then
                E_mp2_b = aux_eri_tmp_x(n_state_shifted,m_state_shifted)
                E_mp2_tmp_ss = E_mp2_a*(E_mp2_a - E_mp2_b)/e_diff
                E_mp2_tmp_os = 0.0d0
              else
                E_mp2_tmp_os = E_mp2_a*E_mp2_a/e_diff
                E_mp2_tmp_ss = 0.0d0
              end if
            end if

            ! If symmetry is exploited, take symetric term int account
            if (i_spin.ne.i_spin_2 .or. a_state.ne.b_state) then
                E_mp2_tmp_os = E_mp2_tmp_os*2.0d0
                E_mp2_tmp_ss = E_mp2_tmp_ss*2.0d0
            end if

            pt2_term_os = pt2_term_os + E_mp2_tmp_os
            pt2_term_ss = pt2_term_ss + E_mp2_tmp_ss

          enddo ! n_state for the virtual states of i_kp_point
        enddo ! m_state for the virtual states of i_qp_point
      enddo ! b_state for the occupied states of i_q_point
    enddo ! a_state for the occupied states of i_k_point
    enddo ! i_spin
    enddo ! i_spin_2

    pt2_c_energy_local_os = pt2_c_energy_local_os + pt2_term_os * k3_weight

    pt2_c_energy_local_ss = pt2_c_energy_local_ss + pt2_term_ss * k3_weight

  end subroutine cpt2_kq_list_gamma_only


end module cpt2_kspace_mod
