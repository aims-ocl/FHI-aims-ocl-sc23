  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/distribute_irreducible_eigenvectors
  !  NAME
  !    distribute_irreducible_eigenvectors
  !  SYNOPSIS

  subroutine distribute_irreducible_eigenvectors &
             ( KS_eigenvector, KS_eigenvector_complex, &
               KS_eigenvector_irk, KS_eigenvector_complex_irk )

    !  PURPOSE
    !
    !   distribute the KS eigenectors on a set of irreducible k points
    !   over the processors
    !   
    !
    !  USES

    use dimensions
    use runtime_choices
    use pbc_lists
    use scalapack_wrapper
    use mpi_tasks
    use synchronize_mpi_basic
    use aims_memory_tracking
    
    implicit none

    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_basis,n_states,n_spin,n_irk_points_task) :: KS_eigenvector_irk
    complex*16, dimension(n_basis,n_states,n_spin,n_irk_points_task) :: KS_eigenvector_complex_irk

    !  ARGUMENTS

    !  INPUTS
    !  o  KS_eigenvector -- real eigenvectors from the single-particle (KS/HF) self-consistent calculation
    !                    on full set of k points
    !  o  KS_eigenvector_complex --  complex eigenvectors from the single-particle (KS/HF) self-consistent 
    !                    calculation on full set of k points
    !  OUTPUTS
    !  o  KS_eigenvector_irk -- real eigenvectors from the single-particle (KS/HF) self-consistent calculation
    !                    on a irreducible set of k points
    !  o  KS_eigenvector_complex_irk --  complex eigenvectors from the single-particle (KS/HF) self-consistent 
    !                    calculation on a irreducible set of k points
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales

    integer, dimension(:), allocatable :: kpoint_task
    integer, dimension(:), allocatable :: id_kpoint
    real*8, allocatable :: KS_eigenvector_tmp(:,:,:)
    complex*16, allocatable :: KS_eigenvector_complex_tmp(:,:,:)
    integer :: info

    ! counter
    integer :: i_k_point
    integer :: i_k_point_local
    integer :: i_irk_point
    integer :: i_irk_point_local
    integer :: i_basis_1
    integer :: i_task
    integer :: id_send
    integer :: id_recv

    character(*), parameter :: func = 'distribute_irreducible_eigenvectors'

    if(real_eigenvectors) then
        allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin),stat=info) 
        call check_allocation(info, 'KS_eigenvector_tmp', func)
        allocate(KS_eigenvector_complex_tmp(1,1,1),stat=info)
        call check_allocation(info, 'KS_eigenvector_complex_tmp', func)
    else
        allocate(KS_eigenvector_tmp(1,1,1),stat=info)
        call check_allocation(info, 'KS_eigenvector_tmp', func)
        allocate(KS_eigenvector_complex_tmp(n_basis,n_states,n_spin),stat=info) 
        call check_allocation(info, 'KS_eigenvector_complex_tmp', func)
    endif
  
    if(use_scalapack) then

        call aims_allocate(kpoint_task, n_tasks, 'kpoint_task')
        call aims_allocate(id_kpoint, n_k_points, 'id_kpoint')
!     determine the first MPI task that a give k point lives on 

          kpoint_task(:) = 0
          do i_task = 1, n_tasks, 1
            if(myid .eq. i_task -1) then
              kpoint_task (i_task) = my_k_point
            endif
          enddo
          call sync_integer_vector(kpoint_task, n_tasks)
          do i_k_point = 1, n_k_points
            do i_task = 1, n_tasks, 1
              if(i_k_point .eq. kpoint_task(i_task)) then
                id_kpoint(i_k_point) = i_task-1
                if(i_k_point .eq. i_task - 1) then
                  exit
                endif
              endif
            enddo
          enddo
     endif


   do i_k_point = 1, n_k_points, 1

      if( .not. irk_point_included(i_k_point)) cycle

      i_irk_point = irk_point_mapping(i_k_point)
      i_k_point_local = (i_k_point-1)/n_tasks + 1
      i_irk_point_local = (i_irk_point-1)/n_tasks + 1
    
      if(use_scalapack) then
         id_send=mod(id_kpoint(i_k_point),n_tasks)
      else
         id_send = mod(i_k_point, n_tasks)
      endif
      id_recv = mod(i_irk_point, n_tasks)

      if(myid.eq.id_send) then
        if(real_eigenvectors) then
           KS_eigenvector_tmp(:,:,:) = KS_eigenvector(:,:,:,i_k_point_local)
        else
           KS_eigenvector_complex_tmp(:,:,:) = KS_eigenvector_complex(:,:,:,i_k_point_local)
        endif
      endif

      if(id_send .ne. id_recv) then
        if(real_eigenvectors) then
          if(myid .eq. id_send) then
             call send_real_vector(KS_eigenvector_tmp,n_basis*n_states*n_spin,id_recv)
          elseif(myid .eq. id_recv) then
             call receive_real_vector(KS_eigenvector_tmp,n_basis*n_states*n_spin,id_send)
          endif
        else
          if(myid .eq. id_send) then
             call send_complex_vector(KS_eigenvector_complex_tmp,n_basis*n_states*n_spin,id_recv)
          elseif(myid .eq. id_recv) then
             call receive_complex_vector(KS_eigenvector_complex_tmp,n_basis*n_states*n_spin,id_send)
          endif
        endif
      endif

      if(myid.eq.id_recv) then
        if(real_eigenvectors) then
          KS_eigenvector_irk(:,:,:,i_irk_point_local) = KS_eigenvector_tmp(:,:,:)
        else
          KS_eigenvector_complex_irk(:,:,:,i_irk_point_local) = KS_eigenvector_complex_tmp(:,:,:)
        endif
      endif

   enddo

   if(allocated(kpoint_task))then
      call aims_deallocate(kpoint_task, 'kpoint_task')
   endif
   if(allocated(id_kpoint))then
      call aims_deallocate(id_kpoint, 'id_kpoint')
   endif
   deallocate(KS_eigenvector_tmp)
   deallocate(KS_eigenvector_complex_tmp)

   end subroutine distribute_irreducible_eigenvectors

