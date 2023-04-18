  module cpt2_blacs
    use mpi_tasks
    use dimensions
    use localorb_io
    use runtime_choices, only : flag_frozen_core_postSCF, flag_frozen_core, i_start_mp2
    use hartree_fock, only : n_homo_max, n_homo

    implicit none

    integer               :: blacsdim_cpt2, pair_block_cpt2

    !communicators
    integer:: comm_blacs_row_cpt2, comm_blacs_col_cpt2, comm_cart_cpt2, comm_blacs_cpt2, comm_kq_cpt2
    integer:: n_tasks_row_cpt2, n_tasks_col_cpt2, n_tasks_bl_cpt2, n_tasks_kq_cpt2
    integer:: myid_row_cpt2, myid_col_cpt2, myid_col_bl_cpt2, myid_bl_cpt2, myid_kq_cpt2
    integer:: pb_bl_col, lpb_col, upb_col, n_pb_col
    integer:: pb_bl_row, lpb_row, upb_row, n_pb_row
    !BLACS descriptors
    integer, dimension(9):: b1desc, bf1desc, bb2desc, bbf2desc, pb2desc, pp2desc
    !integer, dimension(:,:) :: blacs_ids

    logical               :: kblacs_member

    !module variables, move them here for the determine of pp2desc
    integer, allocatable :: n_homo_k(:,:) ! the HOMO level at a given k-point q
    integer, allocatable :: n_lumo_k(:,:) ! # of non-fully occupied states at a given k-point q
    integer, dimension(2) :: n_lumo
    integer :: n_unocc_max, n_lumo_min, n_occ_max, n_low_state_cpt2
    !new distribution of KS_eigenvectors and lvl_tricoeff for MP2 implementation
    complex(kind=8), allocatable, dimension(:,:,:,:,:), target :: lvl_tricoeff_occ, lvl_tricoeff_unocc
    complex(kind=8), allocatable, dimension(:,:,:,:), target ::KS_eigenvectors_occ, KS_eigenvectors_unocc

    real(kind=8), allocatable, dimension(:,:,:,:,:), target :: lvl_tricoeff_occ_real, lvl_tricoeff_unocc_real
    real(kind=8), allocatable, dimension(:,:,:,:), target :: KS_eigenvectors_occ_real, KS_eigenvectors_unocc_real
    real(kind=8), allocatable, dimension(:,:,:), target :: coulomb_matr_blacs_real

    real*8  :: eex_energy_cpt2

    ! flag for gamma-only calculation
    logical :: gamma_only_cpt2 = .false.

  contains

    subroutine distribute_tasks_cpt2(n_tasks,n_k_points,n_tasks_bl_cpt2,n_tasks_kq_cpt2)

      integer, intent(in):: n_tasks, n_k_points
      integer, intent(out):: n_tasks_bl_cpt2,n_tasks_kq_cpt2
      integer:: max_sqrt_bl, i, bldim, kdim, eff_k_points
      real:: total, total_max, imbalance

      max_sqrt_bl=floor(sqrt(real(n_tasks)))
      total_max=0.
      
      do i=max_sqrt_bl,1,-1
         bldim=i*i
         kdim=min(n_tasks/bldim,n_k_points)
         if (mod(n_k_points,kdim).gt.0) then
            !compute the load imbalance due to different number of k points on the different tasks
            eff_k_points=(n_k_points/kdim+1)*kdim
            imbalance=1.0*eff_k_points/n_k_points
         else
            imbalance=1.
         endif

         !efficiency of the parallelisation (fraction of tasks used divided by imbalance)
         total=1.0*bldim*kdim/n_tasks/imbalance
         if (total.ge.total_max) then
            n_tasks_bl_cpt2=bldim
            n_tasks_kq_cpt2=kdim
            total_max=total
         end if
      end do

    end subroutine distribute_tasks_cpt2

    subroutine init_cpt2_blacs_distribution(occ_numbers)
 
!      use synchronize_mpi
      implicit none

      real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers 

      integer, dimension(3):: dims
      logical, dimension(3):: periods, remain
      integer:: ld_bb, ld_pp
      integer:: cntxt
      integer, external:: numroc
      integer:: i_k_point, i_spin, a_state
      integer:: info, mpierr
      character(*), parameter :: func = 'mod_cpt2_blacs.f90'



      ! homo and lumo information (in particular n_lumo_min) is necessary for the initialization of pp2decs
      allocate(n_homo_k(n_k_points,n_spin),stat=info) 
      call check_allocation(info, 'n_homo_k', func)
      allocate(n_lumo_k(n_k_points,n_spin),stat=info) 
      call check_allocation(info, 'n_unocc_k', func)

      if (flag_frozen_core_postSCF) then ! count the frozen core states
          call count_frozen_core_states(n_low_state_cpt2)
      else if (flag_frozen_core) then
          n_low_state_cpt2 = i_start_mp2
      else
          n_low_state_cpt2 = 1
      endif


      n_homo_k(:,:) = 0
      n_lumo_k(:,:)= n_states
      do i_spin = 1, n_spin, 1
        do i_k_point = 1, n_k_points, 1
          do a_state = 1, n_states
            if(occ_numbers(a_state,i_spin,i_k_point) .gt. 1.e-12) then
             n_homo_k(i_k_point,i_spin) = a_state
            endif
          enddo
         
          do a_state = n_states, 1, -1
            if(occ_numbers(a_state,i_spin,i_k_point) .lt. 1.d0) then
             n_lumo_k(i_k_point,i_spin) = a_state 
            endif
          enddo
        enddo
      enddo
      n_homo(1)=maxval(n_homo_k(:,1), n_k_points)
      n_homo(n_spin)=maxval(n_homo_k(:,n_spin), n_k_points)
      n_homo_max = max(n_homo(1), n_homo(n_spin)) 
      n_occ_max = n_homo_max - n_low_state_cpt2 + 1

      n_lumo(1)=minval(n_lumo_k(:,1), n_k_points)
      n_lumo(n_spin)=minval(n_lumo_k(:,n_spin), n_k_points)
      n_lumo_min = min(n_lumo(1), n_lumo(n_spin)) 
      n_unocc_max = n_states - n_lumo_min + 1

      !split MPI_COMM_WORLD into communicators for parallelization over kq points and parallelization of the matrices for one kq
      if ((blacsdim_cpt2.ge.1).and.(blacsdim_cpt2.le.n_tasks)) then
         !keep parallelisation if it has been set in control.in
         n_tasks_bl_cpt2=blacsdim_cpt2
         n_tasks_kq_cpt2=n_tasks/blacsdim_cpt2
      else
         if (blacsdim_cpt2.gt.n_tasks) write(use_unit,*) &
             'blacsdim from control.in file is larger than number of MPI tasks and is therefore ignored'
         !get best distribution
         call distribute_tasks_cpt2(n_tasks,n_k_points,n_tasks_bl_cpt2,n_tasks_kq_cpt2)
      end if
      

      !allocate(blacs_ids(10,n_tasks_bl_cpt2),stat=info)
      !call check_allocation(info, 'blacs_ids', func)
      !blacs_ids=0

      n_tasks_row_cpt2 = INT( SQRT( REAL(n_tasks_bl_cpt2) ) )
      n_tasks_col_cpt2 = n_tasks_bl_cpt2 / n_tasks_row_cpt2
      if (n_tasks_row_cpt2 .ne. n_tasks_col_cpt2) &
          call aims_stop_coll("processes in row should be the same as those in column", func)

      dims=(/n_tasks_kq_cpt2,n_tasks_row_cpt2,n_tasks_col_cpt2/)
      periods=(/.false., .false., .false./)
      call mpi_cart_create(MPI_COMM_WORLD,3,dims,periods,.false.,comm_cart_cpt2,mpierr)

      !block sizes
      if (mod(n_basbas,n_tasks_row_cpt2).eq.0) then
         bb_bl_row=n_basbas/n_tasks_row_cpt2
      else
         bb_bl_row=n_basbas/n_tasks_row_cpt2+1
      end if
      if (mod(n_basbas,n_tasks_col_cpt2).eq.0) then
         bb_bl_col=n_basbas/n_tasks_col_cpt2
      else
         bb_bl_col=n_basbas/n_tasks_col_cpt2+1
      end if
      ! block sizes for n_basis
      if (mod(n_unocc_max,n_tasks_row_cpt2).eq.0) then
         pb_bl_row=n_unocc_max/n_tasks_row_cpt2
      else
         pb_bl_row=n_unocc_max/n_tasks_row_cpt2+1
      end if
      if (mod(n_unocc_max,n_tasks_col_cpt2).eq.0) then
         pb_bl_col=n_unocc_max/n_tasks_col_cpt2
      else
         pb_bl_col=n_unocc_max/n_tasks_col_cpt2+1
      end if

      if (myid.ge.n_tasks_kq_cpt2*n_tasks_bl_cpt2) then
         !tasks that do not participate in the (k x n_basbas²) parallelisation
         kblacs_member=.false.
         n_kq_points_task = 0
         n_ks_points_task = 0
         lbb_row=0
         ubb_row=0
         lbb_col=0
         ubb_col=0
      else
         !cartesian communicators have a c-like ordering, i.e. the last index is the 'fast' one and should be used for the BLACS context, the first one is used for the kq parallelization
         kblacs_member=.true.
         remain=(/.true., .false., .false./)
         call mpi_cart_sub(comm_cart_cpt2, remain, comm_kq_cpt2, mpierr)
!         call mpi_comm_size(comm_kq_cpt2,n_tasks_kq_cpt2,mpierr)
         call mpi_comm_rank(comm_kq_cpt2,myid_kq_cpt2,mpierr)

         remain=(/.false., .true., .true./)
         call mpi_cart_sub(comm_cart_cpt2, remain, comm_blacs_cpt2, mpierr)
!         call mpi_comm_size(comm_blacs_cpt2,n_tasks_bl_cpt2,mpierr)
         call mpi_comm_rank(comm_blacs_cpt2,myid_bl_cpt2,mpierr)

         remain=(/.false., .true., .false./)
         call mpi_cart_sub(comm_cart_cpt2, remain, comm_blacs_row_cpt2, mpierr)

         remain=(/.false., .false., .true./)
         call mpi_cart_sub(comm_cart_cpt2, remain, comm_blacs_col_cpt2, mpierr)

         n_kq_points_task = 0
         do i_k_point = 1, n_k_points
            if(myid_kq_cpt2.eq.mod(i_k_point-1,n_tasks_kq_cpt2) .and. myid_kq_cpt2 .lt. n_k_points) then
               n_kq_points_task = n_kq_points_task + 1
            endif
         enddo
         
         !new local k range for arrays that are (blacsdim)x(k) parallel
         n_ks_points_task = 0
         do i_k_point = 1, n_k_points
            if(myid_kq_cpt2.eq.mod(i_k_point-1,n_tasks_kq_cpt2) .and. myid_kq_cpt2 .lt. n_k_points) then
               n_ks_points_task = n_ks_points_task + 1
            endif
         enddo         

         !create BLACS context for each kq process
         cntxt=comm_blacs_cpt2         
         if (myid.eq.0) print*,'communicator split kq, row, col: ', n_tasks_kq_cpt2, n_tasks_row_cpt2, n_tasks_col_cpt2
         call blacs_gridinit(cntxt,'r',n_tasks_row_cpt2, n_tasks_col_cpt2)
         call blacs_gridinfo(cntxt, n_tasks_row_cpt2, n_tasks_col_cpt2, myid_row_cpt2, myid_col_cpt2)
                  
         n_bb_row = NUMROC( n_basbas, bb_bl_row, myid_row_cpt2, 0, n_tasks_row_cpt2 )
         n_bb_col = NUMROC( n_basbas, bb_bl_col, myid_col_cpt2, 0, n_tasks_col_cpt2 )

         n_pb_row = NUMROC( n_unocc_max, pb_bl_row, myid_row_cpt2, 0, n_tasks_row_cpt2 )
         n_pb_col = NUMROC( n_unocc_max, pb_bl_col, myid_col_cpt2, 0, n_tasks_col_cpt2 )
         
         !leading dimension for descinit call
         ld_bb=max(1,n_bb_row)
         ld_pp=max(1,n_pb_row)
         call descinit(bb2desc, n_basbas, n_basbas, bb_bl_row, bb_bl_col, 0, 0, cntxt, ld_bb, info)
         ! for left and right ovlp3fn matrices
         call descinit(pb2desc, n_basbas, n_unocc_max, bb_bl_row, pb_bl_col, 0, 0, cntxt, ld_bb, info)
         call descinit(pp2desc, n_unocc_max, n_unocc_max, pb_bl_row, pb_bl_col, 0, 0, cntxt, ld_pp, info)
         
         
         !lower and upper bound for n_basbas dimensions
         lbb_row=myid_row_cpt2*bb_bl_row+1
         ubb_row=min(myid_row_cpt2*bb_bl_row+n_bb_row,n_basbas)
         lbb_col=myid_col_cpt2*bb_bl_col+1
         ubb_col=min(myid_col_cpt2*bb_bl_col+n_bb_col,n_basbas)
         !lower and upper bound for n_basis dimensions
         !lpb_row=(n_lumo_min-1)+myid_row_cpt2*pb_bl_row+1
         !upb_row=(n_lumo_min-1)+min(myid_row_cpt2*pb_bl_row+n_pb_row,n_unocc_max)
         !lpb_col=(n_lumo_min-1)+myid_col_cpt2*pb_bl_col+1
         !upb_col=(n_lumo_min-1)+min(myid_col_cpt2*pb_bl_col+n_pb_col,n_unocc_max)
         lpb_row=myid_row_cpt2*pb_bl_row+1
         upb_row=min(myid_row_cpt2*pb_bl_row+n_pb_row,n_unocc_max)
         lpb_col=myid_col_cpt2*pb_bl_col+1
         upb_col=min(myid_col_cpt2*pb_bl_col+n_pb_col,n_unocc_max)

      end if

    end subroutine init_cpt2_blacs_distribution

  end module cpt2_blacs
