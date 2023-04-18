
  module crpa_blacs
    use mpi_tasks
    use dimensions
    use localorb_io

    implicit none

    integer               :: blacsdim

    !communicators
    integer:: comm_blacs_row, comm_blacs_col, comm_cart, comm_blacs, comm_irkq
    integer:: n_tasks_row, n_tasks_col, n_tasks_bl, n_tasks_irkq
    integer:: myid_row, myid_col, myid_bl, myid_irkq

    !BLACS descriptors
    integer, dimension(9):: b1desc, bf1desc, bb2desc, bbf2desc

    logical               :: irkblacs_member

  contains

    subroutine distribute_tasks(n_tasks,n_irk_points,n_tasks_bl,n_tasks_irkq)

      integer, intent(in):: n_tasks, n_irk_points
      integer, intent(out):: n_tasks_bl,n_tasks_irkq
      integer:: max_sqrt_bl, i, bldim, kdim, eff_k_points
      real:: total, total_max, imbalance

      max_sqrt_bl=floor(sqrt(real(n_tasks)))
      total_max=0.
      
      do i=max_sqrt_bl,1,-1
         bldim=i*i
         kdim=min(n_tasks/bldim,n_irk_points)
         if (mod(n_irk_points,kdim).gt.0) then
            !compute the load imbalance due to different number of k points on the different tasks
            eff_k_points=(n_irk_points/kdim+1)*kdim
            imbalance=1.0*eff_k_points/n_irk_points
         else
            imbalance=1.
         endif

         !efficiency of the parallelisation (fraction of tasks used divided by imbalance)
         total=1.0*bldim*kdim/n_tasks/imbalance
         if (total.ge.total_max) then
            n_tasks_bl=bldim
            n_tasks_irkq=kdim
            total_max=total
         end if
      end do

    end subroutine distribute_tasks

    subroutine init_crpa_blacs_distribution(n_full_freq)
 
!      use synchronize_mpi
      implicit none
      integer, intent(in):: n_full_freq
      integer, dimension(3):: dims
      logical, dimension(3):: periods, remain
      integer:: ld_bb
      integer:: cntxt
      integer, external:: numroc
      integer:: i_irk_point, i_k_point
      integer:: info, mpierr

      !split MPI_COMM_WORLD into communicators for parallelization over irkq points and parallelization of the matrices for one irkq

      if ((blacsdim.ge.1).and.(blacsdim.le.n_tasks)) then
         !keep parallelisation if it has been set in control.in
         n_tasks_bl=blacsdim
         n_tasks_irkq=n_tasks/blacsdim
      else
         if (blacsdim.gt.n_tasks) write(use_unit,*) 'blacsdim from control.in file is larger than number of MPI tasks and is therefore ignored'
         !get best distribution
         call distribute_tasks(n_tasks,n_irk_points,n_tasks_bl,n_tasks_irkq)
      end if

      n_tasks_row = INT( SQRT( REAL(n_tasks_bl) ) )
      n_tasks_col = n_tasks_bl / n_tasks_row

      dims=(/n_tasks_irkq,n_tasks_row,n_tasks_col/)
      periods=(/.false., .false., .false./)
      call mpi_cart_create(MPI_COMM_WORLD,3,dims,periods,.false.,comm_cart,mpierr)

      !block sizes
      if (mod(n_basbas,n_tasks_row).eq.0) then
         bb_bl_row=n_basbas/n_tasks_row
      else
         bb_bl_row=n_basbas/n_tasks_row+1
      end if
      if (mod(n_basbas,n_tasks_col).eq.0) then
         bb_bl_col=n_basbas/n_tasks_col
      else
         bb_bl_col=n_basbas/n_tasks_col+1
      end if

      if (myid.ge.n_tasks_irkq*n_tasks_bl) then
         !tasks that do not participate in the (k x n_basbas²) parallelisation
         irkblacs_member=.false.
         n_irkq_points_task = 0
         n_ks_points_task = 0
         lbb_row=0
         ubb_row=0
         lbb_col=0
         ubb_col=0
      else
         !cartesian communicators have a c-like ordering, i.e. the last index is the 'fast' one and should be used for the BLACS context, the first one is used for the irkq parallelization
         irkblacs_member=.true.
         remain=(/.true., .false., .false./)
         call mpi_cart_sub(comm_cart, remain, comm_irkq, mpierr)
!         call mpi_comm_size(comm_irkq,n_tasks_irkq,mpierr)
         call mpi_comm_rank(comm_irkq,myid_irkq,mpierr)

         remain=(/.false., .true., .true./)
         call mpi_cart_sub(comm_cart, remain, comm_blacs, mpierr)
!         call mpi_comm_size(comm_blacs,n_tasks_bl,mpierr)
         call mpi_comm_rank(comm_blacs,myid_bl,mpierr)

         remain=(/.false., .true., .false./)
         call mpi_cart_sub(comm_cart, remain, comm_blacs_row, mpierr)

         remain=(/.false., .false., .true./)
         call mpi_cart_sub(comm_cart, remain, comm_blacs_col, mpierr)

         n_irkq_points_task = 0
         do i_irk_point = 1, n_irk_points
            if(myid_irkq.eq.mod(i_irk_point-1,n_tasks_irkq) .and. myid_irkq .lt. n_irk_points) then
               n_irkq_points_task = n_irkq_points_task + 1
            endif
         enddo
         
         !new local k range for arrays that are (blacsdim)x(k) parallel
         n_ks_points_task = 0
         do i_k_point = 1, n_k_points
            if(myid_irkq.eq.mod(i_k_point-1,n_tasks_irkq) .and. myid_irkq .lt. n_k_points) then
               n_ks_points_task = n_ks_points_task + 1
            endif
         enddo         

         !create BLACS context for each irkq process
         cntxt=comm_blacs         
         if (myid.eq.0) print*,'communicator split irkq, row, col: ', n_tasks_irkq, n_tasks_row, n_tasks_col
         call blacs_gridinit(cntxt,'r',n_tasks_row, n_tasks_col)
         call blacs_gridinfo(cntxt, n_tasks_row, n_tasks_col, myid_row, myid_col)
                  
         n_bb_row = NUMROC( n_basbas, bb_bl_row, myid_row, 0, n_tasks_row )
         n_bb_col = NUMROC( n_basbas, bb_bl_col, myid_col, 0, n_tasks_col )
         
         !leading dimension for descinit call
         ld_bb=max(1,n_bb_row)
         call descinit(bb2desc, n_basbas, n_basbas, bb_bl_row, bb_bl_col, 0, 0, cntxt, ld_bb, info)
         call descinit(bbf2desc, n_basbas, n_basbas*n_full_freq, bb_bl_row, bb_bl_col*n_full_freq, 0, 0, cntxt, ld_bb, info)
         
         
         !lower and upper bound for n_basbas dimensions
         lbb_row=myid_row*bb_bl_row+1
         ubb_row=min(myid_row*bb_bl_row+n_bb_row,n_basbas)
         lbb_col=myid_col*bb_bl_col+1
         ubb_col=min(myid_col*bb_bl_col+n_bb_col,n_basbas)
      end if
    end subroutine init_crpa_blacs_distribution

  end module crpa_blacs
