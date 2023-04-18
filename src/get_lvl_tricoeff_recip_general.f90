  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/get_lvl_tricoeff_recip_general
  !  NAME
  !    get_lvl_tricoeff_recip_general
  !  SYNOPSIS

  subroutine get_lvl_tricoeff_recip_general &
             ( n_cells, n_cells_task, n_kpoints_special, n_kpoints_special_task, &
               kphase_special, lvl_tricoeff_recip)

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
    use basis
    use geometry
    use pbc_lists, only : cell_index
    use localorb_io
    use mpi_tasks
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_cells
    integer, intent(IN) :: n_cells_task
    integer, intent(IN) :: n_kpoints_special
    integer, intent(IN) :: n_kpoints_special_task
    complex*16, intent(IN)  :: kphase_special(n_cells,n_kpoints_special)
    complex*16, intent(OUT) :: lvl_tricoeff_recip(max_n_basbas_sp,n_basis,n_basis,n_kpoints_special_task)

    !  INPUTS
    !    o n_cells_task :: the number of unit cells per task
    !    o n_cells      :: the number of realspace unit cells in the calculations (note this typically 
    !                      goes beyond the Born-von-Karmen supercell)
    !    o n_kpoints_special      :: the number of k points in the 1st BZ, not nessessarily same as the regualr
    !                                k grid.
    !    o n_kpoints_special_task :: the number of k points in the 1st BZ per task
    !    o kphase_special :: the K phase factor for a special set of k points (not necessarily the
    !                         regular k grid)
    !  OUTPUTS
    !    o lvl_tricoeff_recip: LVL triple coefficents in k space for first atom
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
    integer :: id_root, mpierr, i_k_point, i_k_point_local
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip_general'

!    if(output_priority == OL_high) then
!      write(info_str,'(2X,A)') "Computing the triple LVL expansion coefficents in real space and"
!      call localorb_info(info_str)
!      write(info_str,'(2X,A)') "Fourier transform from real space to reciprocal space ..."
!      call localorb_info(info_str)
!    endif

    if(.not. allocated(coeff_3fn)) then
      allocate(coeff_3fn(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2, n_cells_task), stat=info)
      call check_allocation(info, 'coeff_3fn', func)
    endif
    if(.not. allocated(coeff_3fn_tmp)) then
      allocate(coeff_3fn_tmp(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2), stat=info)
      call check_allocation(info, 'coeff_3fn_tmp', func)
    endif

    lvl_tricoeff_recip(:,:,:,:) = (0.d0,0.d0)

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
             if(myid.eq.id_root) then
               i_cell_local = (i_cell-1)/n_tasks + 1
               coeff_3fn_tmp(:,:,:,:) = coeff_3fn(:,:,:,:,i_cell_local)
             endif

             call mpi_bcast(coeff_3fn_tmp, size(coeff_3fn_tmp), &
                             MPI_REAL8, id_root, mpi_comm_global, mpierr)

             ! The on-site coefficients have to be halved
             if(i_cell==1) coeff_3fn_tmp(:,:,:,:) = 0.5*coeff_3fn_tmp(:,:,:,:)

             ! The coeffs for the first atom need a regular Fourier transform

!             if(n_kpoints_special .eq. 1) then
!               do i_sp_basis_2 = 1, n_sp_basis_2
!                  do i_sp_basis_1 = 1, n_sp_basis_1
!                      lvl_tricoeff_recip(1:n_spbb_1, basis_off_1+i_sp_basis_1,         &
!                       &                   basis_off_2+i_sp_basis_2, 1) =  &
!                      lvl_tricoeff_recip(1:n_spbb_1, basis_off_1+i_sp_basis_1,         &
!                       &                   basis_off_2+i_sp_basis_2, 1) +  &
!                       & coeff_3fn_tmp(i_sp_basis_1,i_sp_basis_2,1:n_spbb_1, 1)*kphase_special(i_cell,1)
!                  enddo
!               enddo
!             else
               do i_k_point = 1, n_kpoints_special, 1

                  if(myid.eq.mod(i_k_point,n_tasks)) then
                     i_k_point_local = (i_k_point-1)/n_tasks + 1

                     do i_sp_basis_2 = 1, n_sp_basis_2
                       do i_sp_basis_1 = 1, n_sp_basis_1
                          lvl_tricoeff_recip(1:n_spbb_1, basis_off_1+i_sp_basis_1,         &
                        &                   basis_off_2+i_sp_basis_2, i_k_point_local) =  &
                          lvl_tricoeff_recip(1:n_spbb_1, basis_off_1+i_sp_basis_1,         &
                        &                   basis_off_2+i_sp_basis_2, i_k_point_local) +  &
                        & coeff_3fn_tmp(i_sp_basis_1,i_sp_basis_2,1:n_spbb_1, 1)*kphase_special(i_cell,i_k_point)
                       enddo
                     enddo
                  endif
                enddo
!   end of if n_k_points_speical .eq. 1

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

  end subroutine get_lvl_tricoeff_recip_general
