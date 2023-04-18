module evaluate_polarisability_kspace_mod

  !  USES
  use dimensions
  use prodbas, only : atom2basbas_off, sp2n_basbas_sp, max_n_basbas_sp, basbas_atom
  use hartree_fock, only : lvl_tricoeff_mod_r, n_homo_max, n_homo, kq_point_list
  use hartree_fock_p0
  use mpi_tasks
  use synchronize_mpi
  use pbc_lists
  use runtime_choices
  use timing
  use constants
  use geometry, only : species
  use crpa_blacs
  use exchange_trico
  use exchange_ev
  implicit none

  private
  public evaluate_polarisability_kspace, evaluate_polarisability_kspace_list

  !MODULE VARIABLES
  integer, allocatable :: n_homo_k(:,:) ! the HOMO level at a given k-point q
  integer, allocatable :: n_lumo_k(:,:) ! # of non-fully occupied states at a given k-point q
  complex*16, dimension(:,:,:,:), allocatable, target :: lvl_tricoeff_k, lvl_tricoeff_q ! LVL triple coefficients in k space
  complex*16, dimension(:,:,:,:), pointer :: lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r ! LVL triple coefficients in k space

  character(*), parameter :: func = 'evaluate_polarisability_kspace.f90'

  integer, dimension(:), allocatable:: lb_atom, ub_atom

!  logical:: use_adjoint=.true.
  logical:: use_adjoint=.false.
  logical:: doublereal

  contains

  !  NAME evaluate_polarisability_kspace
  !   
  !  SYNOPSIS
  
  subroutine evaluate_polarisability_kspace &
       ( n_low_state, omega_n, KS_eigenvalue, KS_eigenvector, &
       KS_eigenvector_complex, &
       occ_numbers, polar_kspace)

    !  PURPOSE
    !  Subroutine evaluate_polarisability_kspace.f90 is a wrapper for evaluate_polarisability_kspace_list 
    !  for a single frequency point. 
    !

    !  ARGUMENTS

    real*8 :: omega_n
    real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers
    complex*16, dimension(n_basbas,n_basbas,n_irkq_points_task) :: polar_kspace
    real*8, dimension(:,:), allocatable :: polar_kspace_real
    integer :: max_irk_points_task, i_irkq_point_local, n_low_state

    ! local variable
    ! BL: evaulate_polarisability_kspace_list needs omega_n to be of type
    ! real*8, dimension(n_freq)
    ! if n_freq = 1 this would be a simple real*8
    ! However, the ibm compiler does not allow the "identity" real*8, dim(1) =
    ! real*8
    ! And I am also not sure that this is actually an identity.
    real*8, dimension(1) :: omega_n_vec

    omega_n_vec(1) = omega_n 
    if (n_k_points.eq.1) allocate(polar_kspace_real(n_basbas,n_basbas))

    !because of the global mpi_fence operations, all mpi tasks have to run all iteration so of the irkq loop
    if (mod(n_irk_points,n_tasks_irkq).eq.0) then    
       max_irk_points_task=n_irk_points/n_tasks_irkq
    else
       max_irk_points_task=n_irk_points/n_tasks_irkq+1
    end if
    do i_irkq_point_local=1,max_irk_points_task
       call evaluate_polarisability_kspace_list &
            (n_low_state, i_irkq_point_local, 1, omega_n_vec, KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            occ_numbers, polar_kspace_real, polar_kspace(:,:,i_irkq_point_local))
    end do
    if (n_k_points.eq.1) deallocate(polar_kspace_real)

  end subroutine evaluate_polarisability_kspace 
 
  subroutine evaluate_polarisability_kspace_list &
       (n_low_state, i_irkq_point_local,n_freq, omega_n, KS_eigenvalue, KS_eigenvector, &
       KS_eigenvector_complex, occ_numbers, polar_kspace_real, polar_kspace_complex)

    !  PURPOSE
    !  Subroutine evaluate_polarisability_kspace.f90 evaluates the exact-exchange part of
    !  Hartree-Fock hamiltonian in a periodic system. The algorithm used here are based
    !  on the reciprocal space and localized resolution of identity (RI-LVL)
    !  USES
    use basis, only: basis_atom
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_freq, n_low_state, i_irkq_point_local
    real*8, dimension(n_freq)  :: omega_n
    real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers
    real*8, dimension(lbb_row:ubb_row,lbb_col:ubb_col,n_freq) :: polar_kspace_real
    complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col,n_freq) :: polar_kspace_complex



    !  INPUTS
    !  o  n_freq -- number of frequency points
    !  o  omega_n -- frequency points (in the imaginary axis)
    !  o  occ_numbers -- real array,
    !       the occupation number of the electrons for each eigenstate and each spin
    !  o  KS_eigenvalue -- real array,
    !            the eigenvalue of the single-particle (KS/HF) self-consistent calculation
    !  o  KS_eigenvector -- real array,
    !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
    !  o  KS_eigenvector_complex -- complex array,
    !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
    !  OUTPUTS
    !  o  polar_kspace -- the polarisability matrix in k-space
    !  the exact exchange matrix (the "hf_exchange_matr_complex" in the source code) is defined in MODULE
    !  hartree_fock_p0
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

    ! FIXME: need to change the data structure of lvl_tricoeff_kq_2
    ! FIXME: need to distribute the memory storage of polar_kspace

    character(8) :: myid_string, call_counter_string
    integer :: k, j, i

    complex*16, allocatable :: KS_eigenvector_q(:,:,:) ! KS eigenvector at a given q point in the BZ
    complex*16, allocatable :: KS_eigenvector_k(:,:,:) ! KS eigenvector at a given k point in the BZ

    complex*16, dimension(:,:),allocatable :: adjoint
    integer :: info, mpierr

    integer(kind=MPI_ADDRESS_KIND):: nbytes, offset

    !  integer :: n_homo_max 
    integer :: n_unocc_max
    integer :: n_lumo(n_spin)
    integer :: n_lumo_min
    integer :: n_lumo_eff
    integer :: n_unocc_eff
    !  real*8 :: k_minus_q(3), k_minus_q_lattvec(3)
    
    ! counter 
    integer i_state
    integer i_spin
    integer i_k_point
    integer i_q_point
    integer i_kq_point
    integer i_irkq_point
    integer i_irkq_point_1

    integer i_k_point_old, i_q_point_old
    integer i_basis, i_atom, oldatom
    integer i_basis_1, i_basis_2, count, i_freq
    integer:: lb, ub, lbb, ubb, bbrange

    ! mpi related variables
    integer:: win_tri, win_ev

    integer, dimension(2,n_k_points):: kq_val
    integer, dimension(n_k_points):: adjoint_exists
    integer:: paircount
    integer:: k_task, q_task, i_pair, i_pair2, max_irk_points_task, i_irk
    integer, dimension(2):: kq_val_sav
    logical:: found, newq
    complex*16:: zero, one
    integer:: lpair, upair, pairset, n_adjoint, aesav
    call perfon('ev_pol')

    !initialization of module variables (belong to kq_pair routine)
    allocate(n_homo_k(n_k_points,n_spin),stat=info) 
    call check_allocation(info, 'n_homo_k', func)
    allocate(n_lumo_k(n_k_points,n_spin),stat=info) 
    call check_allocation(info, 'n_unocc_k', func)

    n_homo_k(:,:) = 0
    n_lumo_k(:,:)= n_states
    do i_spin = 1, n_spin, 1
       do i_k_point = 1, n_k_points, 1
          do i_state = 1, n_states
             if(occ_numbers(i_state,i_spin,i_k_point) .gt. 1.e-12) then
                n_homo_k(i_k_point,i_spin) = i_state
             endif
          enddo

          do i_state = n_states, 1, -1
             if(occ_numbers(i_state,i_spin,i_k_point) .lt. 2.d0/dble(n_spin)) then
                n_lumo_k(i_k_point,i_spin) = i_state 
             endif
          enddo
       enddo
    enddo
    n_homo(1)=maxval(n_homo_k(:,1), n_k_points)
    n_homo(n_spin)=maxval(n_homo_k(:,n_spin), n_k_points)
    n_homo_max = max(n_homo(1), n_homo(n_spin)) 

    n_lumo(1)=minval(n_lumo_k(:,1), n_k_points)
    n_lumo(n_spin)=minval(n_lumo_k(:,n_spin), n_k_points) 
    n_lumo_min = min(n_lumo(1), n_lumo(n_spin)) 
    n_unocc_max = n_states - n_lumo_min + 1

    allocate(lb_atom(n_atoms),ub_atom(n_atoms))
    oldatom=-1
    do i_basis = 1, n_basis
!       i_atom = Cbasis_to_atom(i_basis)
       i_atom = basis_atom(i_basis)

       if(i_atom.ne.oldatom) then
          lb_atom(i_atom)=i_basis
          if (i_atom.gt.1) ub_atom(i_atom-1)=i_basis-1
          oldatom=i_atom
       end if
    end do
    ub_atom(n_atoms)=n_basis

    allocate(KS_eigenvector_k(n_basis,n_states,n_spin),stat=info) 
    call check_allocation(info, 'KS_eigenvector_k', func)

    !create windows to access distributed arrays

    if (real_eigenvectors) then
       call init_access_ev_real(n_k_points,n_k_points_task,KS_eigenvector,win_ev)
    else
       call init_access_ev_complex(n_k_points,n_k_points_task,KS_eigenvector_complex,win_ev)
    end if

    if ((n_tasks_irkq.gt.1).and.(myid_col.eq.0)) then
       allocate(lvl_tricoeff_q(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin),stat=info) 
       call check_allocation(info, 'lvl_tricoeff_q', func)
       lvl_tricoeff_recip_q_r => lvl_tricoeff_q
    end if

    !routine variables
    allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
    call check_allocation(info, 'KS_eigenvector_q', func)
    
    call sinit_access_trico(n_k_points,n_ks_points_task,lvl_tricoeff_mod_r,win_tri)
    
    !determine the k and q required
    i_irkq_point=n_tasks_irkq*(i_irkq_point_local-1) + myid_irkq + 1

    n_adjoint=0
    if(n_k_points.gt.1) then        

       if ((n_tasks_irkq.gt.1).and.(myid_col.eq.0)) then
          allocate(lvl_tricoeff_k(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin),stat=info) 
          call check_allocation(info, 'lvl_tricoeff_recip_k_r', func)
          lvl_tricoeff_recip_k_r => lvl_tricoeff_k
       end if

       paircount=0
       kq_val=-1
       adjoint_exists=0
       do i_q_point = 1, n_k_points
          do i_k_point = 1, n_k_points
             i_kq_point = kq_point_list(i_k_point,i_q_point)
             if(i_kq_point .eq. 0) cycle
             
             ! only handle the irreducible k points
             if(.not. irk_point_included(i_kq_point) ) cycle
             
             ! accept only values for the right i_irkq_point
             i_irkq_point_1 = irk_point_mapping(i_kq_point)
             if (i_irkq_point_1.ne.i_irkq_point) cycle

             found=.false.
             if(use_adjoint) then
                !check if adjoint pair is already in the list
                do i_pair=1,paircount             
                   if ((i_k_point.eq.kq_point_list(1,kq_val(2,i_pair))).and.(i_q_point.eq.kq_point_list(1,kq_val(1,i_pair)))) then
                      found=.true.
                      exit
                   end if
                end do
             end if
             if (.not.found) then
                !add new kq pair
                paircount=paircount+1
                kq_val(1,paircount) = i_k_point
                kq_val(2,paircount) = i_q_point
             else
                !the adjoint pair has been found
                adjoint_exists(i_pair)=adjoint_exists(i_pair)+1
             end if
          end do
       end do


       if(use_adjoint) then
          !sort pairs, put entries with adjoint to the front
          i_pair=1
          n_adjoint=paircount
          do while (i_pair.le.n_adjoint)
             if (adjoint_exists(i_pair).ne.1) then
                !find position at the back of the list
                i_pair2=n_adjoint
!                do while ((.not.adjoint_exists(i_pair2)).and.(i_pair2.gt.i_pair))
                do while ((adjoint_exists(i_pair2).eq.0).and.(i_pair2.gt.i_pair))
                   i_pair2=i_pair2-1
                   n_adjoint=n_adjoint-1
                end do
                !reorder
                kq_val_sav=kq_val(:,i_pair)
                kq_val(:,i_pair)=kq_val(:,i_pair2)
                kq_val(:,i_pair2)=kq_val_sav
                aesav=adjoint_exists(i_pair)
                adjoint_exists(i_pair)=adjoint_exists(i_pair2)
                adjoint_exists(i_pair2)=aesav
                n_adjoint=n_adjoint-1
             end if
             i_pair=i_pair+1
          end do

          !check if all entries without adjoint are hermitian
          do i_pair=n_adjoint+1,paircount
             if (kq_val(2,i_pair).ne.kq_point_list(1,kq_val(1,i_pair))) then
                print*,'problem with symmetries in evaluate_polarisability_kspace'
                print*,'please contact RPA/GW developers'
                print*,kq_val(:,i_pair)
                stop
             end if
          end do

          !check for additional symmetries
          adjoint_exists=1
          i_pair=2
          do while (i_pair.le.paircount-1)
             i_pair2=i_pair+1
             do while (i_pair2.lt.paircount)
                if ((kq_val(1,i_pair).eq.kq_val(2,i_pair2)).and.(kq_val(2,i_pair).eq.kq_val(1,i_pair2))) then
                   !i_pair2 can be replaced by complex conjugate of i_pair
                   adjoint_exists(i_pair)=2
                   if (i_pair2.lt.paircount) kq_val(:,i_pair2:paircount-1)=kq_val(:,i_pair2+1:paircount)
                   if (i_pair2.le.n_adjoint) n_adjoint=n_adjoint-1
                   paircount=paircount-1
                   exit
                end if
                i_pair2=i_pair2+1
             end do
             i_pair=i_pair+1
          end do
       end if

       !reorder the pair list to minimize communication
       do i_pair=2, paircount
          found=.false.
          if(.not.use_adjoint) then
             do i_pair2=i_pair,n_k_points
                if ((kq_val(2,i_pair2).eq.kq_point_list(1,kq_val(1,i_pair-1))).and.(kq_val(1,i_pair2).eq.kq_point_list(1,kq_val(2,i_pair-1)))) then
                   found=.true.
                   exit
                end if
             end do
          else
             if(i_pair.ge.n_adjoint) exit
             do i_pair2=i_pair,n_adjoint
                if (kq_val(2,i_pair2).eq.kq_val(1,i_pair-1)) then
                   found=.true.
                   exit
                end if
             end do
          end if
          !reorder
          if (found) then
             kq_val_sav=kq_val(:,i_pair)
             kq_val(:,i_pair)=kq_val(:,i_pair2)
             kq_val(:,i_pair2)=kq_val_sav             
          end if
       end do
    else
       !only one k point!
       paircount=1
       kq_val(1,paircount) = 1
       kq_val(2,paircount) = 1
    end if

    !set paircount to zero if nothing has to be computed on this k task
    if (i_irkq_point_local.gt.n_irkq_points_task) then
       paircount=0
       n_adjoint=0
    end if

    if(n_k_points.gt.1) polar_kspace_complex = (0.d0, 0.d0)

    !because of the global mpi_fence operations, all mpi tasks have to run all iterations of the irkq loop
    i_k_point_old=-1
    i_q_point_old=-1

    !split the loop over the pairs into the pairs that need a +h.c. and the pairs that don't
    do pairset=1,2
       if (pairset.eq.1) then
          !first do the entries where the adjoint has to be added
          lpair=1
          upair=n_adjoint
       else
          !do the other entries
          lpair=n_adjoint+1
          upair=n_k_points
       end if

       do i_pair=lpair,upair
          if (i_pair.gt.paircount) then
             call ssync_trico(win_tri)
             call sync_ev(win_ev)
             cycle
          end if
          
          if(adjoint_exists(i_pair).eq.2) then
             doublereal=.true.
          else
             doublereal=.false.
          end if

          i_k_point=kq_val(1,i_pair)
          i_q_point=kq_val(2,i_pair)

          if (i_k_point.eq.i_q_point) then
             !only i_q_point data needed if not up to date
             call saccess_tricoeff(n_spin,win_tri,i_q_point, &
                  lvl_tricoeff_recip_q_r)
             call access_ev(win_ev,k_point_loc(:,i_q_point),KS_eigenvector_q)

             call kq_pair(KS_eigenvector_q, i_k_point, KS_eigenvector_q, i_q_point, KS_eigenvalue, &
                  occ_numbers, omega_n(1:n_freq), n_freq, n_low_state, n_unocc_max, n_basis, &
                  max_n_basbas_sp, n_k_points_task, lvl_tricoeff_recip_q_r, lvl_tricoeff_recip_q_r, &
                  polar_kspace_real, polar_kspace_complex)
          else
             !both i_q_point and i_k_point_data is needed
             if ((i_q_point.ne.i_k_point_old).or.(n_tasks_irkq.eq.1)) then
                call saccess_2_tricoeffs(n_spin,win_tri,i_q_point,lvl_tricoeff_recip_q_r, &
                     i_k_point,lvl_tricoeff_recip_k_r)
                call access_2_evs(win_ev,k_point_loc(:,i_q_point),KS_eigenvector_q,k_point_loc(:,i_k_point),KS_eigenvector_k)
             else
                !reuse the old data
                call scopy_tricoeff(lvl_tricoeff_recip_k_r,lvl_tricoeff_recip_q_r)
                KS_eigenvector_q=KS_eigenvector_k
                call saccess_tricoeff(n_spin,win_tri,i_k_point, &
                     lvl_tricoeff_recip_k_r)
                call access_ev(win_ev,k_point_loc(:,i_k_point),KS_eigenvector_k)
             end if
             call kq_pair(KS_eigenvector_k, i_k_point, KS_eigenvector_q, i_q_point, KS_eigenvalue, &
                  occ_numbers, omega_n(1:n_freq), n_freq, n_low_state, n_unocc_max, n_basis, &
                  max_n_basbas_sp, n_k_points_task, lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &
                  polar_kspace_real, polar_kspace_complex)
             i_k_point_old=i_k_point
          end if

          i_q_point_old=i_q_point
          
          ! end pair loop
       enddo

       if (use_adjoint.and.(pairset.eq.1).and.(n_adjoint.gt.0)) then       
          !add adjoint contribution (i.e. the -i_q_point,-i_k_point contribution to i_irkq_point that was dropped 
          ! in the computation of kq_pair )
          if (i_irkq_point_local.le.n_irkq_points_task) then             
             allocate(adjoint(n_bb_row,n_bb_col))
             one=1.
             zero=0.
             do i_freq=1,n_freq
                call pztranc( n_basbas, n_basbas, &
                     one, polar_kspace_complex(lbb_row,lbb_col,i_freq), 1, 1, bb2desc, &
                     zero, adjoint, 1, 1, bb2desc )
                polar_kspace_complex(:,:,i_freq) = polar_kspace_complex(:,:,i_freq) + adjoint
             end do
             deallocate(adjoint)
          end if
       end if
    end do

    if ((n_k_points.gt.1).and.(i_irkq_point_local.le.n_irkq_points_task)) then
       i_irkq_point=n_tasks_irkq*(i_irkq_point_local-1) + myid_irkq + 1
       i_k_point = inv_irk_point_mapping(i_irkq_point)
       polar_kspace_complex = polar_kspace_complex * k_weights(i_k_point)    
    end if

    if ((n_tasks_irkq.gt.1).and.(myid_col.eq.0)) then
       deallocate(lvl_tricoeff_q)
       if (n_k_points.gt.1) deallocate(lvl_tricoeff_k)
    end if
    deallocate(KS_eigenvector_q)

    call sfinalize_access_trico(win_tri)
    call finalize_access_ev(win_ev)

    deallocate(n_homo_k)
    deallocate(n_lumo_k)
    deallocate(KS_eigenvector_k)

    deallocate(lb_atom, ub_atom)

    call perfoff

    return

  end subroutine evaluate_polarisability_kspace_list
  !---------------------------------------------------------------------
  !******

subroutine kq_pair(KS_eigenvector_k, i_k_point, KS_eigenvector_q, i_q_point, KS_eigenvalue, occ_numbers,&
     omega_n, n_freqs, n_low_state, n_unocc_max, n_basis, max_n_basbas_sp, n_k_points_task, &
     lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, polar_kspace_real, polar_kspace_complex)

  implicit none
  integer, intent(in) :: n_freqs, n_unocc_max, n_basis, n_k_points_task, max_n_basbas_sp, n_low_state

  complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_states, &
       n_spin), intent(in) :: lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r
  complex*16, dimension(n_basis,n_states,n_spin), intent(in) :: KS_eigenvector_k, KS_eigenvector_q
  real*8, dimension(n_states,n_spin,n_k_points), intent(in) :: KS_eigenvalue
  integer, intent(in) :: i_k_point, i_q_point
  real*8, dimension(n_states,n_spin,n_k_points), intent(in) :: occ_numbers
  real*8, intent(in) :: omega_n(1:n_freqs)
  real*8, dimension(lbb_row:ubb_row,lbb_col:ubb_col,1:n_freqs), intent(inout) :: polar_kspace_real
  complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col,1:n_freqs), intent(inout) :: polar_kspace_complex

! counter 
  integer i_state, lb1, ub1
  integer i_spin, maxpairs

  call perfon('kq_pair')

!!$omp parallel do                                                          &
!!$omp  DEFAULT(NONE)                                                       &
!!$omp private(i_spin)                   &
!!$omp shared(n_basis, lbb_row, ubb_row,    &
!!$omp        n_states, n_spin, lvl_tricoeff_kq_row, lvl_tricoeff_recip_k,lvl_tricoeff_recip_q,               &
!!$omp        lvl_tricoeff_kq_row_2, KS_eigenvector_q, KS_eigenvector_k)        &
!!$omp   SCHEDULE(STATIC)
!!$omp end parallel do

  do i_spin = 1, n_spin
     maxpairs = (n_states_k(i_k_point)-n_low_state+1)*(n_states_k(i_q_point)-n_low_state+1) &
          - (n_lumo_k(i_k_point,i_spin)-n_low_state)*(n_lumo_k(i_q_point,i_spin)-n_low_state) &
          - (n_states_k(i_k_point)-n_homo_k(i_k_point,i_spin))*(n_states_k(i_q_point)-n_homo_k(i_q_point,i_spin))
     if (real_eigenvectors.and.(i_k_point.eq.i_q_point)) then
        polar_kspace_real=0.
        call compute_pol_kq_real(n_low_state, i_spin, n_freqs, i_k_point, &
             i_q_point, n_lumo_k, n_homo_k, occ_numbers, KS_eigenvalue,   &
             KS_eigenvector_q, lvl_tricoeff_recip_q_r, &  
             omega_n, polar_kspace_real, maxpairs)
        if(n_k_points.gt.1) polar_kspace_complex = polar_kspace_complex + polar_kspace_real
     else
        call compute_pol_kq_complex(n_low_state, i_spin, n_freqs, i_k_point, &
             i_q_point, n_lumo_k, n_homo_k, occ_numbers, KS_eigenvalue,   &
             KS_eigenvector_q, KS_eigenvector_k, &
             lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &  
             omega_n, polar_kspace_complex, maxpairs)
     end if
     ! end loop over i_spin
  enddo

  call perfoff

end subroutine kq_pair


subroutine compute_lvl_kq_complex(i_spin, &
     lvl_tricoeff_kq, KS_eigenvector_q, KS_eigenvector_k, lvl_tricoeff_recip_k, &
     lvl_tricoeff_recip_q, lbb,ubb, state_pairs, n_pairs)

  implicit none
  integer :: lbb, ubb, i_state, lb1, ub1
  integer , intent(in) :: i_spin
  complex*16, dimension(lbb:ubb,max_n_basis_sp,n_states,n_spin), intent(in) :: lvl_tricoeff_recip_k, lvl_tricoeff_recip_q  
  complex*16, intent(inout) :: lvl_tricoeff_kq(lbb:ubb,n_pairs) ! Full LVL triple coefficient, one basis index is transformed to

  complex*16, dimension(n_basis,n_states,n_spin), intent(in) :: KS_eigenvector_q, KS_eigenvector_k

  integer:: i_atom, i_species, lb, ub, i_state_2, brange, lbb_atom, ubb_atom, ind, n_pairs, i_pair
  integer, dimension(2,n_pairs):: state_pairs

  call perfon('polind')

!  do i_state_2=lb1,ub1
  do i_pair=1, n_pairs
     i_state=state_pairs(1,i_pair)
     i_state_2=state_pairs(2,i_pair)
     do i_atom=basbas_atom(lbb),basbas_atom(ubb)        
        !range in basis dimension
        lb=lb_atom(i_atom)
        ub=ub_atom(i_atom)
        brange=ub_atom(i_atom)-lb_atom(i_atom)+1
        
        !range in basbas dimension
        i_species = species(i_atom)
        lbb_atom=max(lbb,atom2basbas_off(i_atom)+1)
        ubb_atom=min(ubb,atom2basbas_off(i_atom)+sp2n_basbas_sp(i_species))

        lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)=0.
        do ind=1,brange
           lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)=lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)+&
                conjg(KS_eigenvector_q(lb+ind-1,i_state_2,i_spin)) * lvl_tricoeff_recip_k(lbb_atom:ubb_atom,ind,i_state,i_spin) +&
                KS_eigenvector_k(lb+ind-1,i_state,i_spin) * conjg(lvl_tricoeff_recip_q(lbb_atom:ubb_atom,ind,i_state_2,i_spin))
        end do

     end do
  end do

  call perfoff

end subroutine compute_lvl_kq_complex

subroutine compute_lvl_kq_real(i_spin, &
     lvl_tricoeff_kq, KS_eigenvector, lvl_tricoeff_recip, &
     lbb,ubb, state_pairs, n_pairs)

  implicit none
  integer :: lbb, ubb, i_state, lb1, ub1
  integer , intent(in) :: i_spin
  complex*16, dimension(lbb:ubb,max_n_basis_sp,n_states,n_spin), intent(in) :: lvl_tricoeff_recip
  real*8, intent(inout) :: lvl_tricoeff_kq(lbb:ubb,n_pairs) ! Full LVL triple coefficient, one basis index is transformed to

  complex*16, dimension(n_basis,n_states,n_spin), intent(in) :: KS_eigenvector

  integer:: i_atom, i_species, lb, ub, i_state_2, brange, lbb_atom, ubb_atom, ind, n_pairs, i_pair
  integer, dimension(2,n_pairs):: state_pairs

  call perfon('polind')

  do i_pair=1, n_pairs
     i_state=state_pairs(1,i_pair)
     i_state_2=state_pairs(2,i_pair)
     do i_atom=basbas_atom(lbb),basbas_atom(ubb)        
        !range in basis dimension
        lb=lb_atom(i_atom)
        ub=ub_atom(i_atom)
        brange=ub_atom(i_atom)-lb_atom(i_atom)+1
        
        !range in basbas dimension
        i_species = species(i_atom)
        lbb_atom=max(lbb,atom2basbas_off(i_atom)+1)
        ubb_atom=min(ubb,atom2basbas_off(i_atom)+sp2n_basbas_sp(i_species))

        lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)=0.
        do ind=1,brange
           lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)=lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair) + &
                real(KS_eigenvector(lb+ind-1,i_state_2,i_spin)) * real(lvl_tricoeff_recip(lbb_atom:ubb_atom,ind,i_state,i_spin)) + &
                real(lvl_tricoeff_recip(lbb_atom:ubb_atom,ind,i_state_2,i_spin)) * real(KS_eigenvector(lb+ind-1,i_state,i_spin))
        end do

     end do
  end do

  call perfoff

end subroutine compute_lvl_kq_real
   
subroutine compute_pol_kq_complex(n_low_state, i_spin, n_freqs,&
     i_k_point, i_q_point, n_lumo_k, n_homo_k, occ_numbers, KS_eigenvalue, &
     KS_eigenvector_q, KS_eigenvector_k, &
     lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &  
     omega, polar_kspace_out, maxpairs)

  use dimensions, only : n_states_k

  implicit none

  integer, intent(in) :: i_spin, i_k_point, i_q_point, n_freqs, n_low_state, maxpairs

  integer, intent(in) :: n_lumo_k(n_k_points,n_spin), n_homo_k(n_k_points,n_spin)
  real*8, dimension(n_states,n_spin,n_k_points), intent(in) :: KS_eigenvalue
  complex*16, dimension(n_basis,n_states,n_spin), intent(in) :: KS_eigenvector_q, KS_eigenvector_k
  complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin), intent(in) :: lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r  

  real*8, intent(in), dimension(n_states,n_spin,n_k_points) :: occ_numbers
  real*8, intent(in), dimension(n_freqs) :: omega
  complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col,n_freqs), intent(inout) :: polar_kspace_out
  complex*16, dimension(:,:,:),allocatable :: polar_kspace_tmp

  complex*16, dimension(:,:), allocatable, target:: lvl_tricoeff_kq_row
  complex*16, dimension(:,:), allocatable, target:: lvl_tricoeff_kq_col_arr
  complex*16, dimension(:,:), pointer:: lvl_tricoeff_kq_col

  integer             :: i_state_1, i_state_2, i_freq
  integer             :: n_bbf_task
  complex*16, dimension(:,:,:), allocatable:: lvl_tricoeff_tmp

  real*8              :: occ_diff, e_diff
  complex*16:: pre
  integer:: lb1, ub1, i_state, i_pair, n_pairs, lbsp, ubsp, s_block, n_sp_curr
  integer:: count_row, count_col, mpierr
  integer, dimension(2, maxpairs):: state_pairs

  call perfon('polfreq')

  n_bbf_task=n_bb_col*n_freqs
  i_pair=0

  do i_state = n_low_state, n_states_k(i_k_point)

     if (i_k_point.eq.i_q_point) then
        if (i_state.gt.n_homo_k(i_k_point,i_spin)) exit
        lb1=max(i_state+1,n_lumo_k(i_q_point,i_spin))
        ub1=n_states_k(i_q_point)
        doublereal=.true.
     else
        if(i_state.lt.n_lumo_k(i_k_point,i_spin)) then
           lb1=n_lumo_k(i_q_point,i_spin)
           ub1=n_states_k(i_q_point)
        elseif(i_state.le.n_homo_k(i_k_point,i_spin)) then
           lb1=n_low_state
           ub1=n_states_k(i_q_point)
        else
           lb1=n_low_state
           ub1=n_homo_k(i_q_point,i_spin)
        end if
     end if

     do i_state_2=lb1,ub1
        i_pair=i_pair+1
        state_pairs(1,i_pair)=i_state
        state_pairs(2,i_pair)=i_state_2
     end do
  end do
  n_pairs=i_pair

  !limit the memory consumption by the lvl_tricoeff_tmp array to 800 MB per task
  s_block=min(40000000/(bb_bl_row*(n_freqs+2)),maxpairs)

  allocate(lvl_tricoeff_kq_row(lbb_row:ubb_row,s_block))

  !use new array for lvl_tricoeff_kq_col if needed
  if (myid_row.eq.myid_col) then
     lvl_tricoeff_kq_col => lvl_tricoeff_kq_row
  else
     allocate(lvl_tricoeff_kq_col_arr(lbb_col:ubb_col,s_block))
     lvl_tricoeff_kq_col => lvl_tricoeff_kq_col_arr
  end if  
  allocate(lvl_tricoeff_tmp(lbb_col:ubb_col,n_freqs,s_block))

  count_row=size(lvl_tricoeff_kq_row)
  count_col=size(lvl_tricoeff_kq_col)

  lbsp=1
  ubsp = min(lbsp+s_block-1,n_pairs)
  n_sp_curr=ubsp-lbsp+1 

  if (doublereal) then
     allocate(polar_kspace_tmp(lbb_row:ubb_row,lbb_col:ubb_col,n_freqs))
     polar_kspace_tmp=0.
  end if

  do while (lbsp.le.n_pairs)  
     if (myid_col.eq.0) then
        call compute_lvl_kq_complex(i_spin, lvl_tricoeff_kq_row, KS_eigenvector_q, &
             KS_eigenvector_k, lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, lbb_row, ubb_row, state_pairs(:,lbsp:ubsp), n_sp_curr)
     end if
     call mpi_bcast(lvl_tricoeff_kq_row,count_row,MPI_DOUBLE_COMPLEX,0,comm_blacs_col,mpierr)
     call mpi_bcast(lvl_tricoeff_kq_col,count_col,MPI_DOUBLE_COMPLEX,myid_col,comm_blacs_row,mpierr)

     do i_pair = lbsp, ubsp
        i_state=state_pairs(1,i_pair)
        i_state_1=state_pairs(2,i_pair)
        
        occ_diff = occ_numbers(i_state,i_spin,i_k_point) - occ_numbers(i_state_1,i_spin,i_q_point)     
        e_diff = KS_eigenvalue(i_state,i_spin,i_k_point)- KS_eigenvalue(i_state_1,i_spin,i_q_point)
        
        i_state_2 = i_state_1 - lb1 + 1     
        do i_freq=1, n_freqs
!begin test cohsex
!gw
           pre=occ_diff/(e_diff - (0.d0,1.d0)*omega(i_freq))
!cohsex
!           pre=occ_diff/(e_diff - (0.d0,1.d0)*omega(1))
!end test cohsex
           lvl_tricoeff_tmp(:,i_freq,i_pair-lbsp+1) = pre * lvl_tricoeff_kq_col(:,i_pair-lbsp+1)
        end do
     end do
     
     
     ! comment lvl_tricoeff_kq(i_prodbas_1,1:n_homo_max,1:n_homo_max,i_spin) = lvl_tricoeff_kq_2(i_prodbas_1,1:n_homo_max,1:n_homo_max,i_spin) = 

     if (doublereal) then
        call zgemm('N','C', n_bb_row, n_bbf_task, n_sp_curr, (1.d0,0.d0), &
             lvl_tricoeff_kq_row, n_bb_row, &
             lvl_tricoeff_tmp, n_bbf_task, (1.d0,0.d0), &
             polar_kspace_tmp, n_bb_row )
     else
        call zgemm('N','C', n_bb_row, n_bbf_task, n_sp_curr, (1.d0,0.d0), &
             lvl_tricoeff_kq_row, n_bb_row, &
             lvl_tricoeff_tmp, n_bbf_task, (1.d0,0.d0), &
             polar_kspace_out, n_bb_row )
     end if

     lbsp = lbsp + s_block
     ubsp = min(lbsp+s_block-1,n_pairs)
     n_sp_curr=ubsp-lbsp+1
  end do

  if (doublereal) then
     polar_kspace_out = polar_kspace_out + polar_kspace_tmp + conjg(polar_kspace_tmp)
     deallocate(polar_kspace_tmp)
  end if

  deallocate(lvl_tricoeff_kq_row)
  if (myid_row.ne.myid_col) deallocate(lvl_tricoeff_kq_col_arr)
  deallocate(lvl_tricoeff_tmp)
  call perfoff
end subroutine compute_pol_kq_complex
   
subroutine compute_pol_kq_real(n_low_state, i_spin, n_freqs,&
     i_k_point, i_q_point, n_lumo_k, n_homo_k, occ_numbers, KS_eigenvalue, &
     KS_eigenvector_q, lvl_tricoeff_recip_q_r, &  
     omega, polar_kspace_real, maxpairs)

  use dimensions, only : n_states_k

  implicit none

  integer, intent(in) :: i_spin, i_k_point, i_q_point, n_freqs, n_low_state, maxpairs

  integer, intent(in) :: n_lumo_k(n_k_points,n_spin), n_homo_k(n_k_points,n_spin)
  real*8, dimension(n_states,n_spin,n_k_points), intent(in) :: KS_eigenvalue
  complex*16, dimension(n_basis,n_states,n_spin), intent(in) :: KS_eigenvector_q
  complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin), intent(in) :: lvl_tricoeff_recip_q_r  

  real*8, intent(in), dimension(n_states,n_spin,n_k_points) :: occ_numbers
  real*8, intent(in), dimension(n_freqs) :: omega

  real*8, dimension(lbb_row:ubb_row,lbb_col:ubb_col,n_freqs),intent(inout) :: polar_kspace_real

  real*8, dimension(:,:), allocatable, target:: lvl_tricoeff_kq_row
  real*8, dimension(:,:), allocatable, target:: lvl_tricoeff_kq_col_arr
  real*8, dimension(:,:), pointer:: lvl_tricoeff_kq_col

  integer             :: i_state_1, i_state_2, i_freq
  integer             :: n_bbf_task

  real*8, dimension(:,:,:), allocatable:: lvl_tricoeff_tmp

  real*8              :: occ_diff, e_diff
  complex*16:: pre
  integer:: lb1, ub1, i_state, i_pair, n_pairs, lbsp, ubsp, s_block, n_sp_curr
  integer:: count_row, count_col, mpierr
  integer, dimension(2, maxpairs):: state_pairs

  call perfon('polfreq')

  n_bbf_task=n_bb_col*n_freqs
  i_pair=0

  do i_state = n_low_state, n_homo_k(i_k_point,i_spin)
     lb1=max(i_state+1,n_lumo_k(i_q_point,i_spin))
     ub1=n_states_k(i_q_point)
     do i_state_2=lb1,ub1
        i_pair=i_pair+1
        state_pairs(1,i_pair)=i_state
        state_pairs(2,i_pair)=i_state_2
     end do
  end do
  n_pairs=i_pair

  !limit the memory consumption by the lvl_tricoeff_tmp array to 800 MB per task
  s_block=min(80000000/(bb_bl_row*(n_freqs+2)),maxpairs)
  if(myid_bl.eq.0) print*,'lvlkqsize:',8*s_block*bb_bl_row*(n_freqs+2)/1000000.
  allocate(lvl_tricoeff_kq_row(lbb_row:ubb_row,s_block))

  !use new array for lvl_tricoeff_kq_col if needed
  if (myid_row.eq.myid_col) then
     lvl_tricoeff_kq_col => lvl_tricoeff_kq_row
  else
     allocate(lvl_tricoeff_kq_col_arr(lbb_col:ubb_col,s_block))
     lvl_tricoeff_kq_col => lvl_tricoeff_kq_col_arr
  end if  
  allocate(lvl_tricoeff_tmp(lbb_col:ubb_col,n_freqs,s_block))

  count_row=size(lvl_tricoeff_kq_row)
  count_col=size(lvl_tricoeff_kq_col)

  lbsp=1
  ubsp = min(lbsp+s_block-1,n_pairs)
  n_sp_curr=ubsp-lbsp+1 

  do while (lbsp.le.n_pairs)
     if (myid_col.eq.0) then
        call compute_lvl_kq_real(i_spin, lvl_tricoeff_kq_row, KS_eigenvector_q, &
             lvl_tricoeff_recip_q_r, lbb_row, ubb_row, state_pairs(:,lbsp:ubsp), n_sp_curr)
     end if
     call perfon('bc1')
     call mpi_bcast(lvl_tricoeff_kq_row,count_row,MPI_DOUBLE_PRECISION,0,comm_blacs_col,mpierr)
     !the processes with myid_row.eq.myid_col have the right data
     call mpi_bcast(lvl_tricoeff_kq_col,count_col,MPI_DOUBLE_PRECISION,myid_col,comm_blacs_row,mpierr)
     call perfoff

     do i_pair = lbsp, ubsp
        i_state=state_pairs(1,i_pair)
        i_state_1=state_pairs(2,i_pair)
        
        occ_diff = occ_numbers(i_state,i_spin,i_k_point) - occ_numbers(i_state_1,i_spin,i_q_point)     
        e_diff = KS_eigenvalue(i_state,i_spin,i_k_point)- KS_eigenvalue(i_state_1,i_spin,i_q_point)
        
        i_state_2 = i_state_1 - lb1 + 1     
        do i_freq=1, n_freqs
           !the complex conjugate due to the (s1,s2) (s2,s1) symmetry can be applied directly here
!begin cohsex test
! gw
           pre=2*real(occ_diff/(e_diff - (0.d0,1.d0)*omega(i_freq)))
! cohsex 
!           pre=2*real(occ_diff/(e_diff - (0.d0,1.d0)*omega(1)))
!end cohsex test
           lvl_tricoeff_tmp(:,i_freq,i_pair-lbsp+1) = pre * lvl_tricoeff_kq_col(:,i_pair-lbsp+1)
        end do
     end do
     
     
     ! comment lvl_tricoeff_kq(i_prodbas_1,1:n_homo_max,1:n_homo_max,i_spin) = lvl_tricoeff_kq_2(i_prodbas_1,1:n_homo_max,1:n_homo_max,i_spin) = 
     
     call dgemm('N','C', n_bb_row, n_bbf_task, n_sp_curr, (1.d0,0.d0), &
          lvl_tricoeff_kq_row, n_bb_row, &
          lvl_tricoeff_tmp, n_bbf_task, (1.d0,0.d0), &
          polar_kspace_real, n_bb_row )

     lbsp = lbsp + s_block
     ubsp = min(lbsp+s_block-1,n_pairs)
     n_sp_curr=ubsp-lbsp+1
  end do

  deallocate(lvl_tricoeff_kq_row)
  if (myid_row.ne.myid_col) deallocate(lvl_tricoeff_kq_col_arr)
  deallocate(lvl_tricoeff_tmp)
  call perfoff
end subroutine compute_pol_kq_real

end module evaluate_polarisability_kspace_mod
