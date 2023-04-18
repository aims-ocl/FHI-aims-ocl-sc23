module g_times_w_p0
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
  use crpa_blacs
  use exchange_trico
  use exchange_ev
  implicit none

  save
  private
  public compute_g_times_w_p0, kq_pair, init_compute_g_times_w_p0, finalize_compute_g_times_w_p0

  !MODULE VARIABLES
  complex*16, dimension(:,:,:,:), allocatable, target :: lvl_tricoeff_k, lvl_tricoeff_q ! LVL triple coefficients in k space
  complex*16, dimension(:,:,:,:), pointer :: lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r ! LVL triple coefficients in k space

  complex*16, allocatable :: KS_eigenvector_q(:,:,:) ! KS eigenvector at a given q point in the BZ
  complex*16, allocatable :: KS_eigenvector_k(:,:,:) ! KS eigenvector at a given k point in the BZ

  integer, dimension(:), allocatable:: lb_atom, ub_atom
  real*8, dimension(:), allocatable :: omega_full, womega_full, omega
  integer:: max_freq, n_freq, n_full_freq, n_low_state, n_KS_states
  integer:: dcsize, blocksize, n_blocks

contains

  subroutine init_compute_g_times_w_p0(in_n_low_state, in_n_KS_states,&
       in_max_freq, in_n_freq, in_n_full_freq, in_omega_full, in_womega_full, in_omega)

    use basis, only: basis_atom
    implicit none

    integer :: in_n_low_state
    integer :: in_n_KS_states
    integer :: in_n_freq
    integer :: in_n_full_freq
    integer :: in_max_freq
    real*8 :: in_omega_full(in_n_full_freq)
    real*8 :: in_womega_full(in_n_full_freq)
    real*8 :: in_omega(in_n_freq)

    character(*), parameter :: func = 'compute_g_times_w_p0.f90'
    integer:: info, mpierr
    integer:: i_basis, i_atom, oldatom, totalstates, rest

    max_freq=in_max_freq
    n_freq=in_n_freq
    n_full_freq=in_n_full_freq
    n_low_state=in_n_low_state
    n_KS_states=in_n_KS_states
    allocate(omega_full(n_full_freq), womega_full(n_full_freq), omega(n_freq))
    omega_full=in_omega_full
    womega_full=in_womega_full
    omega=in_omega
    
    allocate(lb_atom(n_atoms),ub_atom(n_atoms))
! start test
!     write(use_unit,*) "myid, i_basis, Cbasis_to_atom(i_basis), basis_atom(i_basis)"
!    do i_basis=1, n_basis, 1
!      write(use_unit,*) myid, i_basis, Cbasis_to_atom(i_basis), basis_atom(i_basis)
!    enddo
! end test
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
    allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
    call check_allocation(info, 'KS_eigenvector_q', func)


    totalstates=n_KS_states - n_low_state + 1
    !limit the memory consumption by the lvl_tricoeff_kq arrays to 800 MB per task
    blocksize=min(50000000/(3*bb_bl_row*n_states),totalstates)
!    blocksize=4

    n_blocks=totalstates/blocksize
    rest=mod(totalstates,blocksize)
    if (rest.gt.0) n_blocks=n_blocks+1
    if (myid.eq.0) print*,'ojo gwblocks',blocksize,totalstates, n_blocks



    call mpi_type_size(MPI_DOUBLE_COMPLEX,dcsize,mpierr) 


  end subroutine init_compute_g_times_w_p0

  subroutine finalize_compute_g_times_w_p0
    deallocate(KS_eigenvector_k)
    deallocate(KS_eigenvector_q)



    deallocate(lb_atom, ub_atom)
    deallocate(omega_full, womega_full, omega)
  end subroutine finalize_compute_g_times_w_p0

  subroutine compute_g_times_w_p0 &
       ( n_homo_k, i_irkq_point_local, KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
       KS_eigenvector_irk, KS_eigenvector_complex_irk, &
       occ_numbers, chemical_potential_spin, &
       screened_coulomb, gw_selfenergy1, delta_gw_selfenergy1)
    
    !  PURPOSE
    !  Subroutine compute_g_times_w_p0.f90 evaluates the exact-exchange part of
    !  Hartree-Fock hamiltonian in a periodic system. The algorithm used here are based
    !  on the reciprocal space and localized resolution of identity (RI-LVL)
    !
    !  USES
    
    use dimensions
    use prodbas
    use hartree_fock
    use hartree_fock_p0
    use mpi_tasks
    use synchronize_mpi
    use pbc_lists
    use runtime_choices
    use timing
    use constants
    implicit none
    
    !  ARGUMENTS
    
    integer :: n_homo_k(n_spin,n_k_points)
    integer i_irkq_point_local
    real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_basis,n_states,n_spin,n_irk_points_task) :: KS_eigenvector_irk
    complex*16, dimension(n_basis,n_states,n_spin,n_irk_points_task) :: KS_eigenvector_complex_irk
    real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers
    real*8 :: chemical_potential_spin(n_spin)
    complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col,n_full_freq) :: screened_coulomb
    complex*16, dimension(n_freq,n_low_state:n_KS_states,n_spin,n_irk_points) :: gw_selfenergy1
    complex*16, dimension(n_freq,n_low_state:n_KS_states,n_spin,n_irk_points) :: delta_gw_selfenergy1
    
    !  INPUTS
    !  o  n_low_state -- the lowest single-particle states whose GW self-energy will be calculated
    !  o  n_KS_states -- number of single-particle states whose GW self-energy will be calculated
    !  o  n_freq -- number of frequency points for the self-energy
    !  o  n_full_freq -- number of frequency points for the screened Coulomb matrix
    !  o  max_freq -- the number of frequency point to consider
    !  o  n_homo_k -- the k-point- and spin-depdenent HOMO levels 
    !  o  omega_full -- the frequency grid for W (in the imaginary axis)
    !  o  womega_full -- the integration weight for the frequency grid for W
    !  o  omega -- full frequency grid for the self-energy
    !  o  occ_numbers -- real array,
    !       the occupation number of the electrons for each eigenstate and each spin
    !  o  KS_eigenvalue -- real array,
    !            the eigenvalue of the single-particle (KS/HF) self-consistent calculation
    !  o  KS_eigenvector -- real array,
    !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
    !  o  KS_eigenvector_complex -- complex array,
    !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
    !  o  KS_eigenvector_irk -- the real eigenvector of the single-particle (KS/HF) self-consistent calculation
    !            on an irreducible set of k vectors
    !  o  KS_eigenvector_complex_irk -- the complex eigenvector of the single-particle (KS/HF) self-consistent 
    !            calculation on an irreducible set of k vectors
    !  o  screened_coulomb -- the polarisability matrix in k-space
    !  the exact exchange matrix (the "hf_exchange_matr_complex" in the source code) is defined in MODULE
    !  hartree_fock_p0
    !
    !  INPUTS/OUTPUTS
    !  o gw_selfenergy1 -- complex, periodic GW self-energy, each calling to this routine adds one more
    !                   frequency point in the convulation
    !
    !  INPUTS/OUTPUTS
    !  o delta_gw_selfenergy1 -- a correction term to the GW self-energy at the frequency point omega_n
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
    
    
    integer :: info, mpierr
    character(*), parameter :: func = 'compute_g_times_w_p0.f90'
    character*150 :: info_str
    
    
    ! counter 
    integer i_k_point
    integer i_irk_point
    integer i_q_point
    integer i_kq_point
    integer i_k_point_local
    integer i_irk_point_local
    integer i_q_point_local
    integer i_kq_point_local
    integer i_irkq_point

    integer win_ev, win_tri
    integer paircount, count

    integer(kind=MPI_ADDRESS_KIND):: nbytes, offset
    integer, dimension(:,:), allocatable:: kq_val
    integer, dimension(3):: kq_val_sav
    logical:: found
    integer i_pair, i_pair2, k_task, q_task, i_k_point_old, irkq_task, conjugate
    integer:: max_paircount
    
    integer:: oldatom, i_basis, i_atom

    call perfon('gtimesw')

    if (real_eigenvectors) then
       call init_access_ev_real(n_k_points,n_k_points_task,KS_eigenvector,win_ev)
    else
       call init_access_ev_complex(n_k_points,n_k_points_task,KS_eigenvector_complex,win_ev)
    end if
    call sinit_access_trico(n_k_points,n_ks_points_task,lvl_tricoeff_mod_r,win_tri)

    !routine variables
    if ((n_tasks_irkq.gt.1).and.(myid_col.eq.0)) then
       allocate(lvl_tricoeff_q(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin),stat=info) 
       call check_allocation(info, 'lvl_tricoeff_q', func)
       lvl_tricoeff_recip_q_r => lvl_tricoeff_q
    end if

    if(n_k_points.gt.1) then 
       if ((n_tasks_irkq.gt.1).and.(myid_col.eq.0)) then
          allocate(lvl_tricoeff_k(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin),stat=info) 
          call check_allocation(info, 'lvl_tricoeff_k', func)
          lvl_tricoeff_recip_k_r => lvl_tricoeff_k
       end if

       if (i_irkq_point_local.gt.n_irkq_points_task) then
          i_irkq_point=-1
       else
          i_irkq_point=n_tasks_irkq*(i_irkq_point_local-1)+myid_irkq+1
       end if

       !find number of all k-q points that map to irkq_point_local
       paircount=0
       do i_kq_point=1,n_k_points
          if(irk_point_mapping(i_kq_point).eq.i_irkq_point) then
             !find all relevant k,q pairs
             do i_irk_point=1,n_irk_points
                i_k_point = inv_irk_point_mapping(i_irk_point)
                do i_q_point=1,n_k_points
                   if(kq_point_list(i_k_point,i_q_point).eq.i_kq_point) then
                      paircount=paircount+1
                   end if
                end do
             end do
          end if
       end do

       allocate(kq_val(3,paircount))

       !save all k-q points that map to irkq_point_local
       paircount=0
       do i_kq_point=1,n_k_points
          if(irk_point_mapping(i_kq_point).eq.i_irkq_point) then
             !find all relevant k,q pairs
             do i_irk_point=1,n_irk_points
                i_k_point = inv_irk_point_mapping(i_irk_point)
                do i_q_point=1,n_k_points
                   if(kq_point_list(i_k_point,i_q_point).eq.i_kq_point) then
                      paircount=paircount+1

                      kq_val(1,paircount) = i_k_point
                      kq_val(2,paircount) = i_q_point
                      if (irk_point_included(i_kq_point)) then
                         kq_val(3,paircount) = 0
                      else
                         kq_val(3,paircount) = 1
                      end if
                   end if
                end do
             end do
          end if
       end do

       call mpi_allreduce(paircount, max_paircount, 1, &
            MPI_INTEGER, MPI_MAX, mpi_comm_world, mpierr)       
       
       !reorder the pair list to minimize communication
       do i_pair=2,paircount
          found=.false.
          do i_pair2=i_pair,paircount
             if (kq_val(2,i_pair2).eq.kq_val(1,i_pair-1)) then
                found=.true.
                exit
             end if
          end do
          if (found) then
             kq_val_sav=kq_val(:,i_pair)
             kq_val(:,i_pair)=kq_val(:,i_pair2)
             kq_val(:,i_pair2)=kq_val_sav             
          end if
       end do
    else
       if(irkblacs_member) then
          paircount=1
       else
          paircount=0
       end if
       max_paircount=1
       allocate(kq_val(3,1))
       kq_val(:,1)=1
    end if

    !because of the global mpi_fence operations, all mpi tasks have to run all iterations of the irkq loop
    i_k_point_old=-1
    do i_pair=1, max_paircount

       if (i_pair.gt.paircount) then
          call ssync_trico(win_tri)
          call sync_ev(win_ev)
          cycle
       end if
       
       i_k_point=kq_val(1,i_pair)
       i_q_point=kq_val(2,i_pair)
       conjugate=kq_val(3,i_pair)
       
       i_irk_point = irk_point_mapping(i_k_point)
       if (i_k_point.eq.i_q_point) then
          !only i_q_point data needed
          call saccess_tricoeff(n_spin,win_tri,i_q_point,lvl_tricoeff_recip_q_r)
          call access_ev(win_ev,k_point_loc(:,i_q_point),KS_eigenvector_q)
          call kq_pair(KS_eigenvector_q, KS_eigenvector_q, KS_eigenvalue(:,:,i_q_point), k_weights(i_q_point), &
               lvl_tricoeff_recip_q_r, lvl_tricoeff_recip_q_r, &
               conjugate, n_homo_k(:,i_q_point), chemical_potential_spin, screened_coulomb, &
               gw_selfenergy1(:,:,:,i_irk_point), delta_gw_selfenergy1(:,:,:,i_irk_point))          
       else
          !both i_q_point and i_k_point_data is needed
          if (i_q_point.ne.i_k_point_old) then
             call saccess_2_tricoeffs(n_spin,win_tri,i_q_point,lvl_tricoeff_recip_q_r, &
                  i_k_point,lvl_tricoeff_recip_k_r)
             call access_2_evs(win_ev,k_point_loc(:,i_q_point),KS_eigenvector_q,k_point_loc(:,i_k_point),KS_eigenvector_k)
          else
             call scopy_tricoeff(lvl_tricoeff_recip_k_r,lvl_tricoeff_recip_q_r)
             KS_eigenvector_q=KS_eigenvector_k
             call saccess_tricoeff(n_spin,win_tri,i_k_point,lvl_tricoeff_recip_k_r)
             call access_ev(win_ev,k_point_loc(:,i_k_point),KS_eigenvector_k)
          end if
          call kq_pair(KS_eigenvector_k, KS_eigenvector_q, KS_eigenvalue(:,:,i_q_point), k_weights(i_q_point), &
               lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &
               conjugate, n_homo_k(:,i_q_point), chemical_potential_spin, screened_coulomb, &
               gw_selfenergy1(:,:,:,i_irk_point), delta_gw_selfenergy1(:,:,:,i_irk_point))   
       end if
       i_k_point_old=i_k_point
       ! end pair loop
    enddo

    if ((n_tasks_irkq.gt.1).and.(myid_col.eq.0)) then
       deallocate(lvl_tricoeff_q)
       if(n_k_points.gt.1) deallocate(lvl_tricoeff_k)
    end if

    deallocate(kq_val)



    call sfinalize_access_trico(win_tri)    
    call finalize_access_ev(win_ev)
    !  if(myid.eq.0) then
    !   write(use_unit,'(2X,A,f16.3,A,f16.3,A)') "RPA response function part 1: ", tot_time_polar_1, " s ", tot_clock_time_polar_1, " s. "
    !   write(use_unit,'(2X,A,f16.3,A,f16.3,A)') "RPA response function part 2: ", tot_time_polar_2, " s ", tot_clock_time_polar_2, " s. "
    !   write(use_unit,'(2X,A,f16.3,A,f16.3,A)') "RPA response function part 3: ", tot_time_polar_3, " s ", tot_clock_time_polar_3, " s. "
    !  endif

    call perfoff

    return

  end subroutine compute_g_times_w_p0
  !---------------------------------------------------------------------
  !******

  subroutine kq_pair(KS_eigenvector_k, KS_eigenvector_q, KS_eigenvalue, k_weight, &
       lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r, &
       conjugate, n_homo_k, chemical_potential_spin, screened_coulomb, gw_selfenergy, delta_gw_selfenergy)
    
    implicit none
    
    complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_states, &
         n_spin), intent(in) :: lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r
    complex*16, dimension(n_basis,n_states,n_spin), intent(in) :: KS_eigenvector_k, KS_eigenvector_q
    real*8, dimension(n_states,n_spin), intent(in) :: KS_eigenvalue
    integer, intent(in) :: conjugate
    integer :: n_homo_k(n_spin)
    complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col,n_full_freq):: screened_coulomb
    real*8 :: chemical_potential_spin(n_spin)
    real*8:: k_weight
    complex*16, dimension(n_freq,n_low_state:n_KS_states,n_spin) :: gw_selfenergy, delta_gw_selfenergy
    
    !temporary arrays
    complex*16, dimension(:,:,:), allocatable, target :: lvl_tricoeff_kq_row, lvl_tricoeff_kq_col_arr
    complex*16, dimension(:,:,:), pointer :: lvl_tricoeff_kq_col 


    ! counter 
    integer:: i_state, i_spin, mpierr, bind, i_block, i_state_bl, n_thisblock, count
    
    call perfon('gwkq_pair')
    

    allocate(lvl_tricoeff_kq_row(lbb_row:ubb_row,n_states,blocksize)) 

    !use new array for lvl_tricoeff_kq_col for the tasks where it is not identical with lvl_tricoeff_kq_row
    if (myid_row.eq.myid_col) then
       lvl_tricoeff_kq_col => lvl_tricoeff_kq_row
    else
       allocate(lvl_tricoeff_kq_col_arr(lbb_col:ubb_col,n_states,blocksize))
       lvl_tricoeff_kq_col => lvl_tricoeff_kq_col_arr
    end if

    do i_spin = 1, n_spin
       do i_block = 1, n_blocks
          n_thisblock=min(blocksize,(n_KS_states-n_low_state)+1-(i_block-1)*blocksize)
          if(myid_col.eq.0) then
             do i_state_bl = 1, n_thisblock
                i_state=(i_block-1)*blocksize+i_state_bl-1+n_low_state
                call compute_frequency_independent_part(lvl_tricoeff_kq_row, KS_eigenvector_q(:,:,i_spin), &
                     KS_eigenvector_k(:,:,i_spin), lvl_tricoeff_recip_k_r(:,:,:,i_spin), &
                     lvl_tricoeff_recip_q_r(:,:,:,i_spin), i_state, i_state_bl)
             end do
          end if
          
          if(n_tasks_bl.gt.1) then
             count=n_bb_row*n_states*n_thisblock
             call mpi_bcast(lvl_tricoeff_kq_row,count,MPI_DOUBLE_COMPLEX,0,comm_blacs_col,mpierr)
             count=n_bb_col*n_states*n_thisblock
             call mpi_bcast(lvl_tricoeff_kq_col,count,MPI_DOUBLE_COMPLEX,myid_col,comm_blacs_row,mpierr)
          end if
          
          i_state=(i_block-1)*blocksize+n_low_state
          call compute_frequency_dependent_part(n_homo_k(i_spin), KS_eigenvalue(:,i_spin), k_weight,  &
               lvl_tricoeff_kq_row, lvl_tricoeff_kq_col, conjugate, screened_coulomb, chemical_potential_spin(i_spin), &
               gw_selfenergy(:,i_state:i_state+n_thisblock-1,i_spin), delta_gw_selfenergy(:,i_state:i_state+n_thisblock-1,i_spin), &
               n_thisblock)             
       end do
       ! end loop over i_spin
    enddo

    deallocate(lvl_tricoeff_kq_row)
    if (myid_row.ne.myid_col) then
       deallocate(lvl_tricoeff_kq_col_arr)
    end if
    call perfoff
    
  end subroutine kq_pair


  subroutine compute_frequency_independent_part(lvl_tricoeff_kq, KS_eigenvector_q, KS_eigenvector_k, &
       lvl_tricoeff_recip_k, lvl_tricoeff_recip_q, i_state, i_state_bl)
    use geometry, only: species
    implicit none
    integer :: i_state, i_state_bl

    complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_states), intent(in) :: lvl_tricoeff_recip_k, lvl_tricoeff_recip_q  
    complex*16, intent(inout) :: lvl_tricoeff_kq(lbb_row:ubb_row,n_states,blocksize) ! Full LVL triple coefficient, one basis index is transformed to
    complex*16, intent(in) :: KS_eigenvector_q(n_basis,n_states), KS_eigenvector_k(n_basis,n_states)
    
    integer:: i_atom, i_species, lb, ub, lbb_atom, ubb_atom, i_state_2, brange, ind

    call perfon('gwind')
    
    do i_state_2=1,n_states
       do i_atom=basbas_atom(lbb_row),basbas_atom(ubb_row)        
          !range in basis dimension
          lb=lb_atom(i_atom)
          ub=ub_atom(i_atom)
          brange=ub_atom(i_atom)-lb_atom(i_atom)+1
          
          !range in basbas dimension
          i_species = species(i_atom)
          lbb_atom=max(lbb_row,atom2basbas_off(i_atom)+1)
          ubb_atom=min(ubb_row,atom2basbas_off(i_atom)+sp2n_basbas_sp(i_species))

          lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_2,i_state_bl)=0. 
          
          do ind=1,brange
             lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_2,i_state_bl) = lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_2,i_state_bl) &
                  + conjg(lvl_tricoeff_recip_q(lbb_atom:ubb_atom, ind, i_state_2)) * KS_eigenvector_k(lb-1+ind,i_state)
          end do
          
          do ind=1,brange
             lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_2,i_state_bl) = lvl_tricoeff_kq(lbb_atom:ubb_atom,i_state_2,i_state_bl) + &
                  conjg(KS_eigenvector_q(lb-1+ind,i_state_2)) * lvl_tricoeff_recip_k(lbb_atom:ubb_atom, ind, i_state)
          end do
       end do
    end do
    call perfoff

  end subroutine compute_frequency_independent_part
  
  
  subroutine compute_frequency_dependent_part(n_homo_k, KS_eigenvalue, k_weight,&
       lvl_tricoeff_kq_row, lvl_tricoeff_kq_col, conjugate, screened_coulomb, chemical_potential_spin, &
       gw_selfenergy, delta_gw_selfenergy, n_thisblock)
    
    implicit none
    complex*16, dimension(n_freq,n_thisblock), intent(inout) :: gw_selfenergy, delta_gw_selfenergy
    integer:: n_homo_k, conjugate, n_thisblock
    real*8:: k_weight
    real*8, dimension(n_states), intent(in) :: KS_eigenvalue
    complex*16, dimension(lbb_row:ubb_row,n_states,blocksize),intent(in) :: lvl_tricoeff_kq_row
    complex*16, dimension(lbb_col:ubb_col,n_states,blocksize),intent(in) :: lvl_tricoeff_kq_col
    complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col,n_full_freq):: screened_coulomb
    real*8:: chemical_potential_spin
    
    complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col):: screened_coulomb_conjg
    complex*16, dimension(:,:,:), allocatable:: screened_coulomb_tmp
    complex*16, dimension(n_states,n_thisblock,n_full_freq):: aux_screened_coulomb
    real*8::  omega_n
    real*8::  womega_n
    complex*16:: omega_minus_en, omega_minus_en2
    complex*16:: delta_gw_tmp
    
    integer:: i_freq, i_freq_1, j_state, i_irk_point, mpierr, i_state_bl
    
    ! intrinsic function
    complex*16 zdotc

    call perfon('gwfreq')
    allocate(screened_coulomb_tmp(lbb_row:ubb_row,n_states,n_thisblock)) 

    call perfon('gwgemm')
    do i_freq=1,max_freq
       if(conjugate.eq.0) then
          call zgemm('N','N', n_bb_row, n_states*n_thisblock, n_bb_col, (1.d0,0.d0), &
               screened_coulomb(lbb_row,lbb_col,i_freq), n_bb_row, &
               lvl_tricoeff_kq_col(lbb_col,1,1), n_bb_col, (0.d0,0.d0), &
               screened_coulomb_tmp(lbb_row,1,1),n_bb_row)
       else
          screened_coulomb_conjg=conjg(screened_coulomb(:,:,i_freq))
          call zgemm('N','N', n_bb_row, n_states*n_thisblock, n_bb_col, (1.d0,0.d0), &
               screened_coulomb_conjg(lbb_row,lbb_col), n_bb_row, &
               lvl_tricoeff_kq_col(lbb_col,1,1), n_bb_col, (0.d0,0.d0), &
               screened_coulomb_tmp(lbb_row,1,1),n_bb_row)
       end if
       call perfon('gwdot')
       do i_state_bl=1,n_thisblock
          do j_state = 1, n_states, 1
             aux_screened_coulomb(j_state,i_state_bl,i_freq) =  &
                  zdotc(n_bb_row, lvl_tricoeff_kq_row(lbb_row,j_state,i_state_bl), 1, &
                  screened_coulomb_tmp(lbb_row,j_state,i_state_bl), 1)
          end do
       end do
       call perfoff
    end do
    call perfoff
    deallocate(screened_coulomb_tmp)

    call mpi_allreduce(MPI_IN_PLACE, aux_screened_coulomb, n_states*max_freq*n_thisblock, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, comm_blacs, mpierr)

    do i_freq=1,max_freq
       omega_n = omega_full(i_freq)
       womega_n = womega_full(i_freq)

       do j_state = 1, n_states
          do i_freq_1 = 1, n_freq, 1
             omega_minus_en = (0.d0,1.d0)*omega(i_freq_1) + chemical_potential_spin - &
                  KS_eigenvalue(j_state)
             do i_state_bl=1,n_thisblock             
                ! also including the negative frequency point -omega_n here
                gw_selfenergy(i_freq_1,i_state_bl) = &
                     gw_selfenergy(i_freq_1,i_state_bl) - &
                     2.d0*omega_minus_en/(omega_minus_en*omega_minus_en + omega_n*omega_n) * &
                     aux_screened_coulomb(j_state,i_state_bl,i_freq)*womega_n*k_weight
!                 if(i_freq_1 .eq. 1 .and. i_freq .eq. max_freq) then
!                   write(use_unit,'(2I4,4f18.8)') i_state_bl, j_state, gw_selfenergy(i_freq_1,i_state_bl),aux_screened_coulomb(j_state,i_state_bl,i_freq) 
!                 endif
             end do
          end do

!  end of loop over j_state
      enddo
!  end of loop over i_freq
    enddo
          
    do i_freq=1, n_freq
       omega_n = omega(i_freq)

       do j_state = 1, n_states

          delta_gw_tmp=(0.d0,0.d0)
          do i_freq_1 = 1, n_full_freq, 1
             !  correction term
             omega_minus_en = (0.d0,1.d0)*omega_n + chemical_potential_spin - &
                  KS_eigenvalue(j_state)
             
             delta_gw_tmp = delta_gw_tmp + &
                  2.d0*omega_minus_en / &
                  (omega_minus_en*omega_minus_en + omega_full(i_freq_1)*omega_full(i_freq_1)) * &
                  womega_full(i_freq_1)
          enddo

          if (j_state .le. n_homo_k) then
             do i_state_bl=1,n_thisblock
                delta_gw_selfenergy(i_freq,i_state_bl) =  &
                     delta_gw_selfenergy(i_freq,i_state_bl) + &
                     (delta_gw_tmp/2.d0/pi - 0.5d0) * &
                     aux_screened_coulomb(j_state,i_state_bl,i_freq)*k_weight
             end do
          else
             do i_state_bl=1,n_thisblock
                delta_gw_selfenergy(i_freq,i_state_bl) = &
                     delta_gw_selfenergy(i_freq,i_state_bl) + &
                     (delta_gw_tmp/2.d0/pi + 0.5d0) * &
                     aux_screened_coulomb(j_state,i_state_bl,i_freq)*k_weight
             end do
          endif
          ! end of loop over j_state
       enddo

       ! end of loop over i_freq
    enddo

    call perfoff
    
  end subroutine compute_frequency_dependent_part
  
end module g_times_w_p0

