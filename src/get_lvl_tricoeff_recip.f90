
  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/get_lvl_tricoeff_recip
  !  NAME
  !    get_lvl_tricoeff_recip
  !  SYNOPSIS

  subroutine get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,lvl_tricoeff_recip2)

    !  PURPOSE
    !
    !    Compute the LVL triple expansion coefficients in real space (i.e., for a set
    !    of Bravais vectors which separate the two unit cells where the basis pairs
    !    live) and Fourier transfrom the LVL triple coefficients from real space
    !    (on a Bravais lattice) to reciprocal space.
    !
    !
    !  USES

    use dimensions
    use prodbas
    use pbc_lists
    use geometry
    use localorb_io
    use timing
    use basis
    use mpi_tasks
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_cells_task
    complex*16, intent(OUT) :: lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_task)
    real*8, intent(OUT)     :: lvl_tricoeff_recip2(max_n_basbas_sp,n_basis,n_basis)

    !  INPUTS
    !    o n_cells_task :: the number of unit cells (within the Born-von Karmen supercell) per task
    !  OUTPUTS
    !    o lvl_tricoeff_recip1: LVL triple coefficents in k space for first atom
    !    o lvl_tricoeff_recip2: LVL triple coefficents in k space for second atom
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    integer :: i_atom_1, i_atom_2, i_atom_aux
    integer :: i_species_1, i_species_2, i_species_aux, i_species_other
    integer :: basis_off_1, basis_off_2, n_sp_basis_1, n_sp_basis_2
    integer :: basis_off_own, basis_off_other, n_sp_basis_own, n_sp_basis_other
    integer :: bboff_1, n_spbb_1, bboff_2, n_spbb_2
    integer :: i_sp_basis_1, i_sp_basis_2
    integer :: i_cell_1, i_cell_2, i_cell_3
    integer :: i_cell, i_cell_local, n_cell_local
    integer :: i_basis_1,i_basis_2,i_basis_3,i_basis_4,i_prodbas_1,i_prodbas_2
    real*8 :: Dvec(3), Cvec(3)
    real*8 :: Rvecs(3,n_cells_task), dummy(1,1,1,1,1,1)
    real*8, allocatable :: coeff_3fn(:,:,:,:,:), coeff_3fn_tmp(:,:,:,:)
    real*8, allocatable :: Cvec_length(:)
    integer :: id_root, mpierr, i_k_point, i_k_point_local
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip'

    call get_timestamps(time_lvl_triple_recip, clock_time_lvl_triple_recip)

    if(output_priority .le. OL_norm ) then
      write(info_str,'(2X,A)') "Computing the triple LVL expansion coefficents in real space and"
      call localorb_info(info_str)
      write(info_str,'(2X,A)') "Fourier transform from real space to reciprocal space ..."
      call localorb_info(info_str)
    endif

    if(.not. allocated(coeff_3fn)) then
      allocate(coeff_3fn(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2, n_cells_task), stat=info)
      call check_allocation(info, 'coeff_3fn', func)
    endif
    if(.not. allocated(coeff_3fn_tmp)) then
      allocate(coeff_3fn_tmp(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2), stat=info)
      call check_allocation(info, 'coeff_3fn_tmp', func)
    endif
    if(.not. allocated(Cvec_length)) then
      allocate(Cvec_length(n_cells), stat=info)
      call check_allocation(info, 'Cvec_length', func)
    endif

    lvl_tricoeff_recip1(:,:,:,:) = (0.d0,0.d0)
    lvl_tricoeff_recip2(:,:,:)   = 0.d0

    do i_atom_1 = 1, n_atoms, 1
       i_species_1 = species(i_atom_1)
       basis_off_1 = atom2basis_off(i_atom_1)
       n_sp_basis_1 = sp2n_basis_sp(i_species_1)
       bboff_1 = atom2basbas_off(i_atom_1)
       n_spbb_1 = sp2n_basbas_sp(i_species_1)

       do i_atom_2 = 1, n_atoms, 1
          i_species_2 = species(i_atom_2)
          basis_off_2 = atom2basis_off(i_atom_2)
          n_sp_basis_2 = sp2n_basis_sp(i_species_2)
          bboff_2 = atom2basbas_off(i_atom_2)
          n_spbb_2 = sp2n_basbas_sp(i_species_2)

          Dvec = coords(:, i_atom_2) - coords(:, i_atom_1)

          i_cell_local=0
          do i_cell = 1, n_cells, 1

             if(myid .ne. mod(i_cell, n_tasks)) cycle
             i_cell_local = (i_cell-1)/n_tasks + 1

             i_cell_1 = cell_index(i_cell, 1)
             i_cell_2 = cell_index(i_cell, 2)
             i_cell_3 = cell_index(i_cell, 3)
             ! distance between two unit cells
             Cvec = matmul(lattice_vector, (/i_cell_1, i_cell_2, i_cell_3/))
             Cvec_length(i_cell)=sqrt(Cvec(1)*Cvec(1) + Cvec(2)*Cvec(2) + Cvec(3)*Cvec(3))

             Rvecs(:, i_cell_local) = Dvec + Cvec
          end do
          n_cell_local = i_cell_local   ! Last value

          ! Get coeffs in real space

          coeff_3fn = 0 ! Safety only
          if(n_cell_local.gt.0) then
            call get_pairwise_coeff_3fn(i_species_1, i_species_2, n_cell_local, Rvecs, coeff_3fn, dummy, .false.)
          endif

          ! Fourier transform

          do i_cell = 1, n_cells, 1
             id_root = mod(i_cell,n_tasks)
             i_cell_local = (i_cell-1)/n_tasks + 1

!             if(Cvec_length(i_cell) .gt. cutCB_rcut) cycle

             if(myid.eq.id_root) then
               coeff_3fn_tmp(:,:,:,:) = coeff_3fn(:,:,:,:,i_cell_local)
             endif

             call mpi_bcast(coeff_3fn_tmp, size(coeff_3fn_tmp), &
                             MPI_REAL8, id_root, mpi_comm_global, mpierr)

             ! The on-site coefficients have to be halved
             if(i_cell==1) coeff_3fn_tmp(:,:,:,:) = 0.5*coeff_3fn_tmp(:,:,:,:)

             ! The coeffs for the first atom need a regular Fourier transform

             do i_k_point = 1, n_k_points, 1

               if(myid.eq.mod(i_k_point,n_tasks)) then
                 i_k_point_local = (i_k_point-1)/n_tasks + 1

                 do i_sp_basis_2 = 1, n_sp_basis_2
                 do i_sp_basis_1 = 1, n_sp_basis_1
                   lvl_tricoeff_recip1(1:n_spbb_1, basis_off_1+i_sp_basis_1,         &
                   &                   basis_off_2+i_sp_basis_2, i_k_point_local) =  &
                   lvl_tricoeff_recip1(1:n_spbb_1, basis_off_1+i_sp_basis_1,         &
                   &                   basis_off_2+i_sp_basis_2, i_k_point_local) +  &
                   & coeff_3fn_tmp(i_sp_basis_1,i_sp_basis_2,1:n_spbb_1, 1)*k_phase(i_cell,i_k_point)
!                   & coeff_3fn_tmp(i_sp_basis_1,i_sp_basis_2,1:n_spbb_1, 1)*conjg(k_phase(i_cell,i_k_point))
!                    if(i_k_point .eq. 1 .and. abs(coeff_3fn_tmp(i_sp_basis_1,i_sp_basis_2,1, 1)).gt.1.e-5) then
!                       write(use_unit,'(6I4,4f16.8)') i_cell, cell_index(i_cell,:), &
!                              basis_off_1+i_sp_basis_1, basis_off_2+i_sp_basis_2, &
!                              lvl_tricoeff_recip1(1:n_spbb_1, basis_off_1+i_sp_basis_1,         &
!                                    basis_off_2+i_sp_basis_2, i_k_point_local), &
!                                  coeff_3fn_tmp(i_sp_basis_1,i_sp_basis_2,1:n_spbb_1, 1)
!                    endif
                 enddo
                 enddo
               endif
            enddo

             ! The coeffs for the second atom are nonzero only for i_cell == 1
             ! The Fourier transform thus is the identity for all k-points

            if((i_cell == 1) .and. (i_atom_2 .ne. i_atom_1)) then
               do i_sp_basis_2 = 1, n_sp_basis_2
               do i_sp_basis_1 = 1, n_sp_basis_1
                 lvl_tricoeff_recip2(1:n_spbb_2, basis_off_1+i_sp_basis_1, &
                 &                   basis_off_2+i_sp_basis_2) =           &
                 & coeff_3fn_tmp(i_sp_basis_1,i_sp_basis_2,1:n_spbb_2, 2)
               enddo
               enddo
            endif

! end loop over i_cell
          enddo

! end loop over i_atom_2
       enddo
! end loop over i_atom_1
    enddo

    if(allocated(coeff_3fn)) then
      deallocate(coeff_3fn)
    endif
    if(allocated(coeff_3fn_tmp)) then
      deallocate(coeff_3fn_tmp)
    endif
    if(allocated(Cvec_length)) then
      deallocate(Cvec_length)
    endif

    call get_times(time_lvl_triple_recip, clock_time_lvl_triple_recip, tot_time_lvl_triple_recip, tot_clock_time_lvl_triple_recip)

  end subroutine get_lvl_tricoeff_recip

module lvl_tricoeff
  use basis
  use prodbas
  use dimensions
  use pbc_lists, only: cell_index, k_point_loc
  use geometry
  use localorb_io
  use timing
  use crpa_blacs
  
  implicit none
  
  real*8,allocatable, dimension(:,:,:,:) :: lvl_tricoeff_cell, bcell_fac
  real*8, allocatable :: trico_cell(:,:,:,:)
  complex*16, allocatable :: trico_k(:,:,:,:,:)
  integer, allocatable:: ap_pos(:,:), sizearr(:), disparr(:)
  
  integer:: n_c3fn, n_c3fn_task, n_c3fn_k_task
  integer:: myoffset, rest, lb_ap, ub_ap, n_ap, lb_b2, ub_b2
  logical:: save_cell
  complex*16, allocatable:: k_phase_save(:,:)    
contains

  subroutine my_fill_coeff_3fn(i_atom_1, i_atom_2, i_species_1, i_species_2, i_cell_in, &
       n_cells, n_atoms, sizeLatticeVector, max_n_basis_sp,max_n_basbas_sp, &
         coords, cell_index, lattice_vector, coeff_3fn)
      
    implicit none
    integer, intent(in) :: i_atom_1, i_atom_2, i_species_1, i_species_2, i_cell_in, n_cells, n_atoms
    integer, intent(in) :: max_n_basis_sp, max_n_basbas_sp, &
         cell_index(n_cells,3), sizeLatticeVector
    real*8, intent(in)   :: coords(3,n_atoms), lattice_vector(3,sizeLatticeVector)
    real*8              :: Dvec(3), Cvec(3), Rvec(3), dummy(1,1,1,1,1,1)
    
    integer             :: i_cell, i_cell_1, i_cell_2, i_cell_3
    
    real*8, intent(out) :: coeff_3fn(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2)
    call perfon('myfill')
    
    if(i_cell_in.gt.n_cells) then
       i_cell=i_cell_in-n_cells
    else
       i_cell=i_cell_in
    end if
    
    Dvec = coords(:, i_atom_2) - coords(:, i_atom_1)
    
    i_cell_1 = cell_index(i_cell, 1)
    i_cell_2 = cell_index(i_cell, 2)
    i_cell_3 = cell_index(i_cell, 3)
    ! distance between two unit cells
    Cvec = matmul(lattice_vector, (/i_cell_1, i_cell_2, i_cell_3/))
    
    Rvec(:) = Dvec + Cvec
    
    
    ! Get coeffs in real space
    
    coeff_3fn = 0.d0 ! Safety only
    call my_get_pairwise_coeff_3fn(i_species_1, i_species_2, 1, Rvec, coeff_3fn, dummy, .false.)
    call perfoff
  end subroutine my_fill_coeff_3fn


  subroutine gw_init_lvl_tricoeff_recip(n_cells, n_cells_task, n_k_points_special, n_k_points_special_task, n_ks_points_special_task,&
       k_phase_special)

    !  PURPOSE
    !
    !    Compute the LVL triple expansion coefficients in real space (i.e., for a set
    !    of Bravais vectors which separate the two unit cells where the basis pairs
    !    live) and Fourier transfrom the LVL triple coefficients from real space
    !    (on a Bravais lattice) to reciprocal space.


    implicit none

    !  ARGUMENTS
    integer, intent(IN) :: n_cells, n_cells_task, n_k_points_special
    integer, intent(IN) :: n_k_points_special_task, n_ks_points_special_task
    complex*16, intent(IN)  :: k_phase_special(n_cells,n_k_points_special)

    !  INPUTS
    !    o n_cells_task :: the number of unit cells (within the Born-von Karmen supercell) per task
    !  OUTPUTS
    !    o lvl_tricoeff_mod: combined LVL triple coefficents in k space for first and second atom
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    integer :: i_atom_1, i_atom_2
    integer :: i_species_1, i_species_2
    integer :: basis_off_1, basis_off_2, n_sp_basis_1, n_sp_basis_2
    integer :: bboff_1, n_spbb_1, bboff_2, n_spbb_2
    integer :: i_sp_basis_1, i_sp_basis_2
    integer :: i_cell, i_cell_map, i, lbi, ubi, lbj, ubj
    integer:: i_c3fn, i_ap, per_task, i_task, i_k_point, i_ks_point_local

    real*8, allocatable :: coeff_3fn(:,:,:,:,:)
    real*8, allocatable:: coeff(:,:,:)
    integer, allocatable:: find_cell(:,:,:,:)

    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip'
    integer:: mpierr
    integer:: sizeLatticeVector

    integer :: omp_get_max_threads,  threadId, numThreads, omp_get_thread_num
    integer:: brest, n_b_task, myboffset, mycount, task
    integer(kind=MPI_ADDRESS_KIND):: nbytes, offset
    integer:: fencecount, maxfence, win_tri


    if (.not.irkblacs_member) then
       call get_timestamps(time_lvl_triple_recip, clock_time_lvl_triple_recip)
       call get_times(time_lvl_triple_recip, clock_time_lvl_triple_recip, tot_time_lvl_triple_recip, tot_clock_time_lvl_triple_recip)
       return
    end if

    call perfon('initlvl')
    call get_timestamps(time_lvl_triple_recip, clock_time_lvl_triple_recip)


!    if(use_gw.and.out_band) then
!       save_cell=.true.
!    else
!       save_cell=.false.
!    end if
    save_cell=.true.


    !compute number of triple coefficients for this task
    n_c3fn=n_atoms*n_atoms*(n_cells+1)
    n_c3fn_task=n_c3fn/n_tasks_bl
    rest=n_c3fn-n_tasks_bl*n_c3fn_task
    if(myid_bl.lt.rest) then
       n_c3fn_task=n_c3fn_task+1
       myoffset=myid_bl*n_c3fn_task
    else
       myoffset=myid_bl*n_c3fn_task+rest
    end if

    !corresponding lower and upper bound for atom pairs for this task
    lb_ap=myoffset/(n_cells+1)
    ub_ap=(myoffset+n_c3fn_task-1)/(n_cells+1)
    n_ap=ub_ap-lb_ap+1
    allocate(ap_pos(3,0:n_atoms*n_atoms-1))

    i_task=0
    if(i_task.lt.rest) then
       per_task=n_c3fn/n_tasks_bl+1
    else
       per_task=n_c3fn/n_tasks_bl
    end if
    i_cell=0
    do i_ap=0,n_atoms*n_atoms-1
       do while(i_cell.ge.per_task)
          i_cell=i_cell-per_task
          i_task=i_task+1
          if(i_task.lt.rest) then
             per_task=n_c3fn/n_tasks_bl+1
          else
             per_task=n_c3fn/n_tasks_bl
          end if
       end do
       !first task and offset that contains the data for the i_ap pair
       ap_pos(1,i_ap)=i_task
       ap_pos(2,i_ap)=i_cell/(n_cells+1)
       if (ap_pos(2,i_ap)*(n_cells+1).lt.i_cell) ap_pos(2,i_ap)=ap_pos(2,i_ap)+1
       if (i_ap.gt.0) then
          !last task that contains data for the i_ap pair
          if (ap_pos(2,i_ap).gt.0) then
             ap_pos(3,i_ap-1)=i_task
          else
             ap_pos(3,i_ap-1)=i_task-1
          end if
       end if
       i_cell=i_cell+n_cells+1
    end do
    ap_pos(3,n_atoms*n_atoms-1)=n_tasks_bl-1

    numThreads = 1
    threadId=1
    sizeLatticeVector = size(lattice_vector,dim=2)
!$  numThreads = omp_get_max_threads()
    allocate(coeff_3fn(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2, numThreads), stat=info)
    call check_allocation(info, 'coeff_3fn', func)

    if (save_cell) then
       allocate(k_phase_save(n_cells,n_k_points_special))
       k_phase_save=k_phase_special
       allocate(trico_cell(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,0:n_c3fn_task-1), stat=info)
       call check_allocation(info, 'trico_cell', func)
       trico_cell=0.
    else
       allocate(trico_k(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,n_ks_points_special_task,n_ap), stat=info)
       call check_allocation(info, 'trico_k', func)
       trico_k=0.
    end if
    allocate(find_cell(2,n_cells+1,n_atoms,n_atoms))
    find_cell=0

    do i_c3fn=0,n_c3fn_task-1
       i_ap=(i_c3fn+myoffset)/(n_cells+1)
       i_cell=i_c3fn+myoffset-i_ap*(n_cells+1)+1
       i_cell_map=mod(i_cell-1,n_cells)+1
       i_atom_2=i_ap/n_atoms+1
       i_atom_1=i_ap-(i_atom_2-1)*n_atoms+1
       i_species_1 = species(i_atom_1)
       basis_off_1 = atom2basis_off(i_atom_1)
       n_sp_basis_1 = sp2n_basis_sp(i_species_1)
       bboff_1 = atom2basbas_off(i_atom_1)
       n_spbb_1 = sp2n_basbas_sp(i_species_1)
       
       i_species_2 = species(i_atom_2)
       basis_off_2 = atom2basis_off(i_atom_2)
       n_sp_basis_2 = sp2n_basis_sp(i_species_2)
       bboff_2 = atom2basbas_off(i_atom_2)
       n_spbb_2 = sp2n_basbas_sp(i_species_2)

       find_cell(1,i_cell,i_atom_1,i_atom_2)=myid_bl
       find_cell(2,i_cell,i_atom_1,i_atom_2)=i_c3fn

       if (i_cell.le.n_cells) then
          call my_fill_coeff_3fn(i_atom_1, i_atom_2, i_species_1, i_species_2, i_cell, &
               n_cells, n_atoms, sizeLatticeVector, max_n_basis_sp, &
               max_n_basbas_sp, coords, cell_index, lattice_vector, coeff_3fn(:,:,:,:,threadId))
          
          ! The on-site coefficients have to be halved
          if (i_cell.eq.1) coeff_3fn(:,:,:,1,threadId)=0.5*coeff_3fn(:,:,:,1,threadId) 

          if(save_cell) then          
             trico_cell(:,:,:,i_c3fn)=coeff_3fn(:,:,:,1,threadId)
          else
             !fourier transform directly
             do i_ks_point_local = 1,n_ks_points_special_task 
                i_k_point=(i_ks_point_local-1)*n_tasks_irkq+myid_irkq+1
                trico_k(:,:,:,i_ks_point_local,i_ap-lb_ap+1) = trico_k(:,:,:,i_ks_point_local,i_ap-lb_ap+1) + &
                     coeff_3fn(:,:,:,1,threadId) * k_phase_special(i_cell_map,i_k_point)
             end do
          end if
       else
          call my_fill_coeff_3fn(i_atom_2, i_atom_1, i_species_2, i_species_1, 1, &
               n_cells, n_atoms, sizeLatticeVector, max_n_basis_sp, &
               max_n_basbas_sp, coords, cell_index, lattice_vector, coeff_3fn(:,:,:,:,threadId))
          
          if(save_cell) then
             ! The on-site coefficients have to be halved, 
             ! transpose to compensate for atom1 <-> atom2 exchange
             do i=1, max_n_basbas_sp
                trico_cell(:,:,i,i_c3fn)=0.5*transpose(coeff_3fn(:,:,i,2,threadId))
             end do

          else
             !fourier transform directly                
             do i_ks_point_local = 1,n_ks_points_special_task 
                i_k_point=(i_ks_point_local-1)*n_tasks_irkq+myid_irkq+1
                do i=1, max_n_basbas_sp    
                   trico_k(:,:,i,i_ks_point_local,i_ap-lb_ap+1) = trico_k(:,:,i,i_ks_point_local,i_ap-lb_ap+1) + &
                        0.5*transpose(coeff_3fn(:,:,i,2,threadId))* k_phase_special(i_cell_map,i_k_point)
                end do
             end do
          end if
       end if
    end do

    deallocate(coeff_3fn)
    call mpi_allreduce(MPI_IN_PLACE,find_cell,size(find_cell),MPI_INTEGER,MPI_SUM,comm_blacs,mpierr)

    if (save_cell) then

       call perfon('redistlvl')

!     !compute number of elements to store for this task
!     n_bcell=n_basis*(n_cells+1)
!     n_bcell_task=n_bcell/n_tasks_col
!     brest=n_bcell-n_bcell_task*n_tasks_col
!     if(myid_col.lt.brest) then
!        n_bcell_task=n_bcell_task+1
!        bmyoffset=myid_col*n_bcell_task
!     else
!        bmyoffset=myid_col*n_bcell_task+brest
!     end if
!     !corresponding lower and upper bound for basis index for this task
!     lb_b2=bmyoffset/(n_cells+1)+1
!     ub_b2=(bmyoffset+n_bcell_task-1)/(n_cells+1)+1
!     n_b2=ub_b2-lb_b2+1

       n_b_task=n_basis/n_tasks_col
       brest=n_basis-n_b_task*n_tasks_col
       if(myid_col.lt.brest) then
          n_b_task=n_b_task+1
          myboffset=myid_col*n_b_task
       else
          myboffset=myid_col*n_b_task+brest
       end if
       lb_b2=myboffset+1
       ub_b2=myboffset+n_b_task
       
       allocate(coeff(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp), stat=info)
       allocate(bcell_fac(lbb_row:ubb_row,max_n_basis_sp,lb_b2:ub_b2,n_cells+1), stat=info)
       call check_allocation(info, 'bcell_fac', func)
       bcell_fac=0.

       maxfence=(basbas_atom(ubb_row)-basbas_atom(lbb_row)+1)*(basis_atom(ub_b2)-basis_atom(lb_b2)+1)*(n_cells+1)
       call mpi_allreduce(MPI_IN_PLACE,maxfence,1,MPI_INTEGER,MPI_MAX,comm_blacs,mpierr)       
       nbytes=int(max_n_basis_sp,MPI_ADDRESS_KIND)*max_n_basis_sp*max_n_basbas_sp*n_c3fn_task
       call mpi_win_create(trico_cell,nbytes,8,MPI_INFO_NULL,comm_blacs,win_tri,mpierr)
       call mpi_win_fence(0,win_tri,mpierr)


       fencecount=0
       do i_atom_1=basbas_atom(lbb_row),basbas_atom(ubb_row)    
          i_species_1 = species(i_atom_1)
          n_sp_basis_1 = sp2n_basis_sp(i_species_1)
          bboff_1 = atom2basbas_off(i_atom_1)
          n_spbb_1 = sp2n_basbas_sp(i_species_1)
          lbi=max(lbb_row-bboff_1,1)
          ubi=min(ubb_row-bboff_1,n_spbb_1)
          nbytes=int(max_n_basis_sp,MPI_ADDRESS_KIND)*max_n_basis_sp*(ubi-lbi+1)
          
          do i_atom_2=basis_atom(lb_b2),basis_atom(ub_b2)
             i_species_2 = species(i_atom_2)
             basis_off_2 = atom2basis_off(i_atom_2)
             n_sp_basis_2 = sp2n_basis_sp(i_species_2)
             lbj=max(lb_b2-basis_off_2,1)
             ubj=min(ub_b2-basis_off_2,n_sp_basis_2)
             nbytes=int(max_n_basis_sp,MPI_ADDRESS_KIND)*max_n_basis_sp*(ubi-lbi+1)
             
             do i_cell=1,n_cells+1
                task=find_cell(1,i_cell,i_atom_1,i_atom_2)
                offset=find_cell(2,i_cell,i_atom_1,i_atom_2)*max_n_basis_sp*max_n_basis_sp*max_n_basbas_sp+ &
                     max_n_basis_sp*max_n_basis_sp*(lbi-1)
                
                call mpi_get(coeff, nbytes, MPI_DOUBLE_PRECISION, task, &
                     offset, nbytes, MPI_DOUBLE_PRECISION, win_tri, mpierr)
                call mpi_win_fence(0,win_tri,mpierr)
                fencecount=fencecount+1
                
                if(i_cell.le.n_cells) then 
                   do i_sp_basis_2 = lbj, ubj
                      do i_sp_basis_1 = 1, n_sp_basis_1
                         do i=lbi,ubi
                            bcell_fac(bboff_1+i, i_sp_basis_1, i_sp_basis_2+basis_off_2,i_cell) =  &
                                 coeff(i_sp_basis_1,i_sp_basis_2,i-lbi+1)
                         end do
                      end do
                   end do
                else

                   !n_cells+1 can be added to 1st element
                   do i_sp_basis_2 = lbj, ubj
                      do i_sp_basis_1 = 1, n_sp_basis_1
                         do i=lbi,ubi
                            bcell_fac(bboff_1+i, i_sp_basis_1, i_sp_basis_2+basis_off_2,1) =  &
                                 bcell_fac(bboff_1+i, i_sp_basis_1, i_sp_basis_2+basis_off_2,1) +  &
                                 coeff(i_sp_basis_1,i_sp_basis_2,i-lbi+1)
                         end do
                      end do
                   end do
                end if
             end do
          end do
       end do
       
       do while (fencecount.lt.maxfence) 
          call mpi_win_fence(0,win_tri,mpierr)   
          fencecount=fencecount+1
       end do

       deallocate(coeff)
       call mpi_win_free(win_tri,mpierr)

       mycount=(ubb_row-lbb_row+1)*max_n_basis_sp*(ub_b2-lb_b2+1)    
       allocate(sizearr(0:n_tasks_col-1))
       allocate(disparr(0:n_tasks_col-1))

       call mpi_allgather(mycount,1,MPI_INTEGER,&
            sizearr,1,MPI_INTEGER,&
            comm_blacs_col,mpierr)
       disparr=0
       do i_task=1,n_tasks_col-1
          disparr(i_task)=disparr(i_task-1)+sizearr(i_task-1)
       end do
       
       call perfoff
       deallocate(trico_cell)
    end if
    deallocate(find_cell)

    call get_times(time_lvl_triple_recip, clock_time_lvl_triple_recip, tot_time_lvl_triple_recip, tot_clock_time_lvl_triple_recip)
    call perfoff

  end subroutine gw_init_lvl_tricoeff_recip

  
  subroutine gw_cleanup_lvl_tricoeff_recip
    if(irkblacs_member) then
       if (save_cell) then
          deallocate(bcell_fac)
          deallocate(k_phase_save)
          deallocate(sizearr)
          deallocate(disparr)
       else
          deallocate(trico_k)
       end if
       deallocate(ap_pos)
    end if

  end subroutine gw_cleanup_lvl_tricoeff_recip


  subroutine gw_get_lvl_tricoeff_recip(n_cells, n_k_points_special, n_k_points_special_task, n_ks_points_special_task,&
       k_phase_special, KS_eigenvector, KS_eigenvector_complex, lvl_tricoeff_mod_r)

    !  PURPOSE
    !
    !    Compute the LVL triple expansion coefficients in real space (i.e., for a set
    !    of Bravais vectors which separate the two unit cells where the basis pairs
    !    live) and Fourier transfrom the LVL triple coefficients from real space
    !    (on a Bravais lattice) to reciprocal space.

    use runtime_choices, only: real_eigenvectors
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_cells, n_k_points_special
    integer, intent(IN) :: n_k_points_special_task, n_ks_points_special_task
    complex*16, intent(IN)  :: k_phase_special(n_cells,n_k_points_special)
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_special_task), intent(IN) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_special_task), intent(IN) :: KS_eigenvector_complex

    complex*16, intent(OUT) :: lvl_tricoeff_mod_r(lbb_row:,:,:,:,:)

    !  INPUTS

    !  OUTPUTS
    !    o lvl_tricoeff_mod: combined LVL triple coefficents in k space for first and second atom
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    integer :: i_cell

    complex*16, allocatable, dimension(:,:,:,:) :: lvl_tricoeff_tmp
    complex*16, allocatable, dimension(:,:,:,:) :: trico_tmp
    complex*16, allocatable, dimension(:,:,:) :: trico_tmp2

    integer :: mpierr, i_k_point, i_k_point_local, i_ks_point_local, i_spin, count, k_task, win_ev, count_new
    integer(kind=MPI_ADDRESS_KIND):: nbytes, offset
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip'
    complex*16, dimension(n_basis,n_states,n_spin) :: KS_eigenvector_k
    real*8, dimension(:,:,:), allocatable :: KS_eigenvector_r
    integer:: i_atom_1, i_atom_2, win_tri, i_exc, n_exc, ppos
    integer:: i_ap, c0pos, fencecount
    integer:: lbi, ubi, i, i_sp_basis_1, i_sp_basis_2, n_spbb_1, bboff_1, i_species_1, i_species_2
    integer:: n_sp_basis_1, n_sp_basis_2, basis_off_2
    integer:: atoms, maxatoms
    logical:: i_collect

    call perfon('getlvl')
    if (irkblacs_member) then
!       if(myid_col.eq.0) then
          allocate(lvl_tricoeff_tmp(lbb_row:ubb_row,max_n_basis_sp,n_basis,n_ks_points_special_task),stat=info)
          call check_allocation(info, 'lvl_tricoeff_tmp', func)
          lvl_tricoeff_tmp=0.
!       end if
       
       if(save_cell) then
          !fourier transform
          allocate(trico_tmp2(lbb_row:ubb_row,max_n_basis_sp,lb_b2:ub_b2),stat=info)

          do i_ks_point_local = 1,n_ks_points_special_task 
!             call perfon('nlvlfour')
             trico_tmp2=0.
             i_k_point=(i_ks_point_local-1)*n_tasks_irkq+myid_irkq+1
             do i_cell=1,n_cells
                trico_tmp2(:,:,lb_b2:ub_b2) = trico_tmp2(:,:,lb_b2:ub_b2) + &
                     bcell_fac(:,:,:,i_cell) * k_phase_special(i_cell,i_k_point)
             end do
!             call perfoff
             
!             call perfon('nlvlaggr')   
             
             call mpi_gatherv(trico_tmp2,sizearr(myid_col),MPI_DOUBLE_COMPLEX,&
                  lvl_tricoeff_tmp(lbb_row,1,1,i_ks_point_local), sizearr, disparr,MPI_DOUBLE_COMPLEX,&
                  0,comm_blacs_col,mpierr)
             
!             call perfoff
          end do
          deallocate(trico_tmp2)

       else

          !aggregate data
          call perfon('lvlagg1')


          nbytes=int(max_n_basis_sp,MPI_ADDRESS_KIND)*max_n_basis_sp*max_n_basbas_sp*n_ks_points_special_task*n_ap
          call mpi_win_create(trico_k,nbytes,16,MPI_INFO_NULL,comm_blacs,win_tri,mpierr)
          call mpi_win_fence(0,win_tri,mpierr)
          
          !check if data for last atom pair is complete and collect if necessary
          c0pos=ub_ap*(n_cells+1)
          i_collect=(c0pos.ge.myoffset).and.(c0pos.lt.myoffset+n_c3fn_task)          
          if(i_collect) allocate(trico_tmp(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,n_ks_points_special_task))

          nbytes=int(max_n_basis_sp,MPI_ADDRESS_KIND)*max_n_basis_sp*max_n_basbas_sp*n_ks_points_special_task
          offset=0
          ppos=myid_bl
          n_exc=(n_cells+1)/(n_c3fn/n_tasks_bl)+1
          i_cell=n_c3fn_task+myoffset-ub_ap*(n_cells+1)

          do i_exc=1,n_exc
             if (i_collect.and.(i_cell.lt.(n_cells+1))) then
                ppos=ppos+1
                i_cell=i_cell+n_c3fn/n_tasks_bl
                if(ppos.lt.rest) i_cell=i_cell+1
                call mpi_get(trico_tmp, nbytes, MPI_DOUBLE_COMPLEX, ppos, &
                     offset, nbytes, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
                call mpi_win_fence(0,win_tri,mpierr)
                trico_k(:,:,:,:,n_ap)=trico_k(:,:,:,:,n_ap)+trico_tmp
             else
                call mpi_win_fence(0,win_tri,mpierr)
             end if
          end do

          if(i_collect) deallocate(trico_tmp)

          call perfoff
          
          call perfon('lvlaggr2')

          
          ! the number of calls to mpi_win_fence has to be the same for all tasks
          ! get the maximum number of atoms in a lbb_row:ubb_row range
          maxatoms=0
          do i = 0, n_tasks_row-1
             lbi = i*bb_bl_row+1
             ubi = min((i+1)*bb_bl_row,n_basbas)
             atoms = basbas_atom(ubi) - basbas_atom(lbi) + 1
             if (atoms.gt.maxatoms) maxatoms=atoms
          end do
          
          !collect data from all pairs that contribute to the lbb_row:ubb_row section
          fencecount=0       
          if(myid_col.eq.0) then   
             allocate(trico_tmp(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,n_ks_points_special_task))
             do i_atom_1=basbas_atom(lbb_row),basbas_atom(ubb_row) 
                i_species_1 = species(i_atom_1)
                n_sp_basis_1 = sp2n_basis_sp(i_species_1)
                bboff_1 = atom2basbas_off(i_atom_1)
                n_spbb_1 = sp2n_basbas_sp(i_species_1)
                
                lbi=max(lbb_row-bboff_1,1)
                ubi=min(ubb_row-bboff_1,n_spbb_1)
                
                do i_atom_2=1,n_atoms
                   i_species_2 = species(i_atom_2)
                   basis_off_2 = atom2basis_off(i_atom_2)
                   n_sp_basis_2 = sp2n_basis_sp(i_species_2)
                   
                   ! compute data position
                   i_ap=(i_atom_2-1)*n_atoms+i_atom_1-1
                   offset=ap_pos(2,i_ap)*nbytes
                   
                   call mpi_get(trico_tmp, nbytes, MPI_DOUBLE_COMPLEX, ap_pos(1,i_ap), &
                        offset, nbytes, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
                   call mpi_win_fence(0,win_tri,mpierr)
                   
                   fencecount=fencecount+1
                   do i_k_point = 1, n_ks_points_special_task
                      do i_sp_basis_2 = 1, n_sp_basis_2
                         do i_sp_basis_1 = 1, n_sp_basis_1
                            do i=lbi,ubi
                               lvl_tricoeff_tmp(bboff_1+i, i_sp_basis_1,         &
                                    basis_off_2+i_sp_basis_2,i_k_point) =  &
                                    lvl_tricoeff_tmp(bboff_1+i, i_sp_basis_1,         &
                                    basis_off_2+i_sp_basis_2,i_k_point) +  &
                                    trico_tmp(i_sp_basis_1,i_sp_basis_2,i,i_k_point)
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             deallocate(trico_tmp)
          end if
          
          do while (fencecount.lt.n_atoms*maxatoms)
             call mpi_win_fence(0,win_tri,mpierr)
             fencecount=fencecount+1
          end do

!          deallocate(trico_tmp)

          call mpi_win_free(win_tri,mpierr)
       
          call perfoff
       end if
    end if

    call perfon('lvlmult')       
    !create windows to access distributed arrays
    !to multiply lvl_tricoeff_k with eigenvector
    if (mod(myid+n_tasks-1,n_tasks).lt.n_k_points_special) then
       nbytes=int(n_basis,MPI_ADDRESS_KIND)*n_states*n_spin*n_k_points_special_task*8
    else
       nbytes=0 
    end if
    
    if(real_eigenvectors) then
       call mpi_win_create(KS_eigenvector,nbytes,8,MPI_INFO_NULL,mpi_comm_world,win_ev,mpierr)
       if (myid_bl.eq.0) allocate(KS_eigenvector_r(n_basis,n_states,n_spin))  
    else
       nbytes=nbytes*2
       call mpi_win_create(KS_eigenvector_complex,nbytes,16,MPI_INFO_NULL,mpi_comm_world,win_ev,mpierr)
    end if

    call mpi_win_fence(0,win_ev,mpierr)

    do i_ks_point_local = 1, n_k_points_special/n_tasks_irkq+1
       i_k_point=(i_ks_point_local-1)*n_tasks_irkq+myid_irkq+1
       if ((i_k_point.gt.n_k_points_special).or.(.not.irkblacs_member)) then
          call mpi_win_fence(0,win_ev,mpierr)
          cycle
       end if

       count=n_basbas*max_n_basis_sp
       !only task 0 of the group with the same k uses mpi_get and then does a bcast to avoid too many p2p messages
       if (myid_bl.eq.0) then
          !get the mapping of the eigenvector parallelization
          k_task=k_point_loc(1,i_k_point)
          i_k_point_local=k_point_loc(2,i_k_point)
          offset=(i_k_point_local-1)*int(n_basis,MPI_ADDRESS_KIND)*n_states*n_spin
          if(real_eigenvectors) then
             call mpi_get(KS_eigenvector_r, n_basis*n_states*n_spin, MPI_DOUBLE_PRECISION, k_task, &
                  offset, n_basis*n_states*n_spin, MPI_DOUBLE_PRECISION, win_ev, mpierr)
          else
             call mpi_get(KS_eigenvector_k, n_basis*n_states*n_spin, MPI_DOUBLE_COMPLEX, k_task, &
                  offset, n_basis*n_states*n_spin, MPI_DOUBLE_COMPLEX, win_ev, mpierr)
          end if
       end if
       call mpi_win_fence(0,win_ev,mpierr)
       if (real_eigenvectors.and.(myid_bl.eq.0)) KS_eigenvector_k=KS_eigenvector_r

       if (.not.irkblacs_member) cycle

       if(myid_col.eq.0) then
          call mpi_bcast(KS_eigenvector_k, n_basis*n_states*n_spin, MPI_DOUBLE_COMPLEX, 0, &
               comm_blacs_row, mpierr)
!          call perfon('lvlgemm')
          do i_spin = 1 , n_spin
             ! new count
             count_new = (ubb_row-lbb_row+1)*max_n_basis_sp
             call zgemm('N','N',count_new,n_states,n_basis,(1.d0,0.d0),&
                  lvl_tricoeff_tmp(lbb_row,1,1,i_ks_point_local),count_new,&
                  KS_eigenvector_k(1,1,i_spin),n_basis,(0.d0,0.d0),&
                  lvl_tricoeff_mod_r(lbb_row,1,1,i_spin,i_ks_point_local),count_new)
          end do
!          call perfoff
       end if
    end do

    call perfoff

    call mpi_win_free(win_ev,mpierr)

    if (allocated(lvl_tricoeff_tmp)) deallocate(lvl_tricoeff_tmp)
    if (allocated(KS_eigenvector_r)) deallocate(KS_eigenvector_r)


    call perfoff
  end subroutine gw_get_lvl_tricoeff_recip



 subroutine gw_get_single_lvl_tricoeff(n_cells, &
       k_phase_special, KS_eigenvector_complex, lvl_tricoeff_mod_r)

    !  PURPOSE
    !
    !    Compute the (myid_row parallel) LVL triple expansion coefficients in real space (i.e., for a set
    !    of Bravais vectors which separate the two unit cells where the basis pairs
    !    live), Fourier transfrom the LVL triple coefficients from real space
    !    (on a Bravais lattice) to reciprocal space, multiplies with the eigenvector
    !    and broadcasts to all myid_col. 

    implicit none

    !  ARGUMENTS

    complex*16, intent(IN)  :: k_phase_special(n_cells)
    integer:: n_cells
    complex*16, dimension(n_basis,n_states,n_spin), intent(IN) :: KS_eigenvector_complex
    complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin), &
         intent(OUT) :: lvl_tricoeff_mod_r


    !  OUTPUTS
    !    o lvl_tricoeff_mod: combined LVL triple coefficents in k space for first and second atom
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    complex*16, allocatable, dimension(:,:,:) :: lvl_tricoeff_tmp
    integer :: mpierr, i_spin, count, i_cell
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip'

    call perfon('getslvl')

    if (.not.irkblacs_member) return

    allocate(lvl_tricoeff_tmp(lbb_row:ubb_row,max_n_basis_sp,lb_b2:ub_b2),stat=info)
    call check_allocation(info, 'lvl_tricoeff_tmp', func)


    !in place
!    call perfon('nlvlfour')
    lvl_tricoeff_tmp = 0.
    do i_cell=1,n_cells
       lvl_tricoeff_tmp = lvl_tricoeff_tmp + &
            bcell_fac(:,:,:,i_cell) * k_phase_special(i_cell)
    end do

!    call perfoff

!    call perfon('lvlgemm')
    do i_spin = 1 , n_spin
       count = n_bb_row*max_n_basis_sp
       call zgemm('N','N',count,n_states,ub_b2-lb_b2+1,(1.d0,0.d0),&
            lvl_tricoeff_tmp(:,:,lb_b2),count,&
            KS_eigenvector_complex(lb_b2,1,i_spin),n_basis,(0.d0,0.d0),&
            lvl_tricoeff_mod_r(lbb_row,1,1,i_spin),count)
    end do
!    call perfoff

    deallocate(lvl_tricoeff_tmp)    

!    call perfon('nlvlaggr')       
    count=size(lvl_tricoeff_mod_r)
    call mpi_allreduce(MPI_IN_PLACE,lvl_tricoeff_mod_r,count,MPI_DOUBLE_COMPLEX,MPI_SUM,&
         comm_blacs_col,mpierr)
    
!    call perfoff

    call perfoff

  end subroutine gw_get_single_lvl_tricoeff

 subroutine gw_get_single_lvl_tricoeff_0(n_cells, &
       k_phase_special, KS_eigenvector_complex, lvl_tricoeff_mod_r)

    !  PURPOSE
    !
    !    Compute the (myid_row parallel) LVL triple expansion coefficients in real space (i.e., for a set
    !    of Bravais vectors which separate the two unit cells where the basis pairs
    !    live), Fourier transfrom the LVL triple coefficients from real space
    !    (on a Bravais lattice) to reciprocal space, and multiplies with the eigenvector.
    !    The result is only stored on myid_col=0.

    implicit none

    !  ARGUMENTS

    complex*16, intent(IN)  :: k_phase_special(n_cells)
    integer:: n_cells
    complex*16, dimension(n_basis,n_states,n_spin), intent(IN) :: KS_eigenvector_complex
    complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin), &
         intent(OUT) :: lvl_tricoeff_mod_r


    !  OUTPUTS
    !    o lvl_tricoeff_mod: combined LVL triple coefficents in k space for first and second atom
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    complex*16, allocatable, dimension(:,:,:) :: lvl_tricoeff_tmp, trico_tmp
    integer :: mpierr, i_spin, count, i_cell
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip'

    call perfon('getslvl0')

    if (.not.irkblacs_member) return

    allocate(lvl_tricoeff_tmp(lbb_row:ubb_row,max_n_basis_sp,lb_b2:ub_b2),stat=info)
    call check_allocation(info, 'lvl_tricoeff_tmp', func)


    !in place
!    call perfon('nlvlfour')
    lvl_tricoeff_tmp = 0.
    do i_cell=1,n_cells
       lvl_tricoeff_tmp = lvl_tricoeff_tmp + &
            bcell_fac(:,:,:,i_cell) * k_phase_special(i_cell)
    end do

!    call perfoff

    if(myid_col.eq.0) allocate(trico_tmp(lbb_row:ubb_row,max_n_basis_sp,n_basis),stat=info)
       
!    call perfon('nlvlaggr')   
    
    call mpi_gatherv(lvl_tricoeff_tmp,sizearr(myid_col),MPI_DOUBLE_COMPLEX,&
         trico_tmp, sizearr, disparr,MPI_DOUBLE_COMPLEX,&
         0,comm_blacs_col,mpierr)
    
!    call perfoff
    
    if(myid_col.eq.0) then       
!       call perfon('lvlgemm')
       do i_spin = 1 , n_spin
          count = n_bb_row*max_n_basis_sp
          call zgemm('N','N',count,n_states,n_basis,(1.d0,0.d0),&
               trico_tmp,count,&
               KS_eigenvector_complex(1,1,i_spin),n_basis,(0.d0,0.d0),&
               lvl_tricoeff_mod_r(lbb_row,1,1,i_spin),count)
       end do
!       call perfoff
       deallocate(trico_tmp)
    end if

    deallocate(lvl_tricoeff_tmp)

    call perfoff
  end subroutine gw_get_single_lvl_tricoeff_0


  subroutine gw_get_single_lvl_tricoeff_noev(n_cells, &
       k_phase_special, lvl_tricoeff)

    !  PURPOSE
    !
    !    Compute the LVL triple expansion coefficients in real space (i.e., for a set
    !    of Bravais vectors which separate the two unit cells where the basis pairs
    !    live) and Fourier transfrom the LVL triple coefficients from real space
    !    (on a Bravais lattice) to reciprocal space.

    implicit none

    !  ARGUMENTS

    complex*16, intent(IN)  :: k_phase_special(n_cells)
    integer:: n_cells
    complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_basis), intent(OUT) :: lvl_tricoeff


    !  INPUTS

    !  OUTPUTS
    !    o lvl_tricoeff: combined LVL triple coefficents in k space for first and second atom
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales

    integer :: i_cell, mpierr

    call perfon('ngetsnlvl')
    if (.not.irkblacs_member) return

    !in place
!    call perfon('nlvlfour')
    lvl_tricoeff = 0.
    do i_cell=1,n_cells
       lvl_tricoeff(:,:,lb_b2:ub_b2) = lvl_tricoeff(:,:,lb_b2:ub_b2) + &
            bcell_fac(:,:,:,i_cell) * k_phase_special(i_cell)
    end do
!    call perfoff

!    call perfon('nlvlaggr')   
    call mpi_allgatherv(MPI_IN_PLACE,sizearr(myid_col),MPI_DOUBLE_COMPLEX,&
         lvl_tricoeff, sizearr, disparr,MPI_DOUBLE_COMPLEX,&
         comm_blacs_col,mpierr)
!    call perfoff

    call perfoff
  end subroutine gw_get_single_lvl_tricoeff_noev


  subroutine get_lvl_col_from_row(lvl_col, lvl_row, n_spin_loc)

    complex*16, dimension(lbb_col:ubb_col,max_n_basis_sp,n_basis,n_spin_loc):: lvl_col
    complex*16, dimension(lbb_row:ubb_row,max_n_basis_sp,n_basis,n_spin_loc):: lvl_row
    integer:: n_spin_loc
    integer:: count, mirror, mpierr
    integer:: status(MPI_STATUS_SIZE)
 
    if (.not.irkblacs_member) return
call perfon('colfromrow')
    count=n_bb_row*max_n_basis_sp*n_basis*n_spin_loc
    mirror=myid_row + myid_col*n_tasks_row
    if (myid_bl.eq.0) then
       lvl_col=lvl_row
    elseif(myid_col.eq.0) then
       call mpi_send(lvl_row, count, MPI_DOUBLE_COMPLEX, mirror, 333, &
            comm_blacs, mpierr)
    elseif(myid_row.eq.0) then
       count=n_bb_col*max_n_basis_sp*n_basis*n_spin_loc
       call mpi_recv(lvl_col, count, MPI_DOUBLE_COMPLEX, mirror, 333, &
            comm_blacs, status, mpierr)
    end if
    count=n_bb_col*max_n_basis_sp*n_basis*n_spin_loc
    call mpi_bcast(lvl_col, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_row, mpierr)
call perfoff
  end subroutine get_lvl_col_from_row
end module lvl_tricoeff
