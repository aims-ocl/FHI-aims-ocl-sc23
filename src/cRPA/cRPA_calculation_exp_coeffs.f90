!Calculates the RI-LVL expansion coefficients
MODULE cRPA_calculation_exp_coeffs
    use localized_basbas
    use tight_binding_auxmat

    USE cRPA_view
    USE cRPA_storage
    USE cRPA_parallelism

    IMPLICIT NONE

    type basis_and_species_offsets
       INTEGER :: i_species, &
                  offset_basis, &
                  offset_basbas, &
                  extend_basis, &
                  extend_basbas,&
                  start_basis, &
                  start_basbas, &
                  end_basis, &
                  end_basbas
    end type basis_and_species_offsets

CONTAINS
    subroutine calc_lvl_tricoeff(exp_coeffs,my_basbas_dist)

    use dimensions
    use prodbas
    use pbc_lists
    use geometry
    use localorb_io
    implicit none


    type(expansion_coefficients) :: exp_coeffs
    type(basbas_distribution),INTENT(IN) :: my_basbas_dist


    ! Local variabales
    integer :: i_atom_1, i_atom_2, i_basbas
    integer :: i_cell_1, i_cell_2, i_cell_3, i_on_or_offsite
    integer :: i_cell
    integer :: mpierr
    integer :: info, i_atom_pair_my
    INTEGER :: i_cell_block_start, i_cell_block_end, n_cell_block_size
    INTEGER :: i_Rvecs

    real*8 :: Dvec(3), Cvec(3), sum_radius
    real*8 :: dummy(1,1,1,1,1,1)
    real*8, allocatable :: coeff_3fn(:,:,:,:,:,:), coeff_one_basbas(:,:), coeff_one_basbas_part(:,:), &
                           Rvecs(:,:)
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip'


    type(basis_and_species_offsets) :: atom_1, atom_2
    type(exp_coeff_calculation_distribution) :: int_dist

    write(info_str,'(2X,A)') "Computing the triple LVL expansion coefficents in real space"
    CALL write_stdout(info_str)

!    CALL initialize_localized_basbas()
    CALL initialize_lvl_triples(OVLP_TYPE_COULOMB)
    CALL initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)

    CALL create_exp_coeff_calc_dist(n_atoms,n_cells, comm_exp_coeffs,int_dist)
    n_cell_block_size = n_cells + 1

    ALLOCATE(coeff_one_basbas(n_basis,n_basis), stat=info)
    CALL check_allocation(info, 'coeff_one_basbas', func)
    ALLOCATE(coeff_3fn(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2, &
                       n_cell_block_size,int_dist%n_atom_pairs_my), stat=info)
    CALL check_allocation(info, 'coeff_3fn', func)
    ALLOCATE(coeff_one_basbas_part(n_basis,n_basis), stat=info)
    CALL check_allocation(info, 'coeff_one_basbas_part', func)

    ALLOCATE(Rvecs(3,n_cell_block_size), stat=info)
    CALL check_allocation(info, 'Rvecs', func)

    do i_cell_block_start = 1, n_cells + 1, n_cell_block_size

        i_cell_block_end = min(i_cell_block_start + n_cell_block_size - 1, &
                               n_cells +1)
        coeff_3fn = 0

       do i_atom_pair_my = 1, int_dist%n_atom_pairs_my
            i_Rvecs = 0

            i_atom_1 = int_dist%i_atom_pair_my_to_atom1(i_atom_pair_my)
            i_atom_2 = int_dist%i_atom_pair_my_to_atom2(i_atom_pair_my)

            CALL set_offsets(i_atom_1, atom_1)
            CALL set_offsets(i_atom_2, atom_2)

            do i_cell = i_cell_block_start, i_cell_block_end

               Dvec = coords(:, i_atom_2) - coords(:, i_atom_1)
               if(i_cell <= n_cells) then
                   i_cell_1 = cell_index(i_cell, 1)
                   i_cell_2 = cell_index(i_cell, 2)
                   i_cell_3 = cell_index(i_cell, 3)
               elseif(i_cell == n_cells + 1) then
                   i_cell_1 = cell_index(1, 1)
                   i_cell_2 = cell_index(1, 2)
                   i_cell_3 = cell_index(1, 3)
               endif
              ! distance between two unit cells

               Cvec = matmul(lattice_vector, (/i_cell_1, i_cell_2, i_cell_3/))
               i_Rvecs = i_Rvecs + 1
               Rvecs(:, i_Rvecs) = Dvec + Cvec
!               if(myid==0) write(use_unit,*) "rvecs", i_Rvecs, Rvecs(:,i_RVecs)
          enddo
CALL write_debug("cell block"//num2str(i_cell_block_start)//num2str(n_cell_block_size)&
                //num2str(i_atom_1)//num2str(i_atom_2))

               call get_pairwise_coeff_3fn(atom_1%i_species, atom_2%i_species, &
                                       i_Rvecs, &
                                       Rvecs, coeff_3fn(:,:,:,:,:,i_atom_pair_my), &
                                       dummy, .false.)
!CALL write_debug("chksum"//num2str(chksum_3d_matrix(coeff_3fn(:,:,:,1,1,i_atom_pair_my))))

       enddo



!           sum_radius = atom_radius(species(i_atom_1)) &
!                           + atom_radius(species(i_atom_2))
!           if(sqrt(SUM(Rvecs(:,1)**2)) >=1.3d0* sum_radius) CYCLE


       do i_cell = i_cell_block_start, i_cell_block_end

       i_Rvecs = i_cell - i_cell_block_start + 1

       if(i_cell <=n_cells) then
         i_on_or_offsite = 1
       elseif(i_cell == n_cells+1) then
         i_on_or_offsite = 2
       endif

       do i_basbas = 1,n_basbas
         coeff_one_basbas(:,:) = 0
         coeff_one_basbas_part(:,:) = 0

         do i_atom_pair_my = 1, int_dist%n_atom_pairs_my
           i_atom_1 = int_dist%i_atom_pair_my_to_atom1(i_atom_pair_my)
           i_atom_2 = int_dist%i_atom_pair_my_to_atom2(i_atom_pair_my)

           CALL set_offsets(i_atom_1, atom_1)
           CALL set_offsets(i_atom_2, atom_2)

            if(i_on_or_offsite==1) then
               if(i_basbas < atom_1%start_basbas &
                  .OR. &
                  atom_1%end_basbas < i_basbas) CYCLE
               coeff_one_basbas_part(atom_1%start_basis:atom_1%end_basis, &
                                     atom_2%start_basis:atom_2%end_basis) =  &
               coeff_one_basbas_part(atom_1%start_basis:atom_1%end_basis, &
                                atom_2%start_basis:atom_2%end_basis) + &
                coeff_3fn(1:atom_1%extend_basis, &
                          1:atom_2%extend_basis, &
                          i_basbas-atom_1%offset_basbas, 1,i_Rvecs,i_atom_pair_my)

!CALL write_debug("basbas"//num2str(i_basbas)//num2str(i_basbas-atom_1%offset_basbas)&
!                 //num2str(i_Rvecs)&
!                 //num2str(chksum_2d_matrix(coeff_3fn(:,:, &
!                          i_basbas-atom_1%offset_basbas, 1,i_Rvecs,i_atom_pair_my)))&
!                                  //num2str(chksum_2d_matrix(coeff_one_basbas_part)) )


            elseif((i_on_or_offsite==2) .and. (i_atom_2 .ne. i_atom_1)) then
               if(i_basbas < atom_2%start_basbas &
                  .OR. &
                  atom_2%end_basbas < i_basbas) CYCLE
               coeff_one_basbas_part(atom_1%start_basis:atom_1%end_basis, &
                                     atom_2%start_basis:atom_2%end_basis) =  &
               coeff_one_basbas_part(atom_1%start_basis:atom_1%end_basis, &
                                     atom_2%start_basis:atom_2%end_basis) + &
                coeff_3fn(1:atom_1%extend_basis, &
                          1:atom_2%extend_basis, &
                          i_basbas-atom_2%offset_basbas,2, i_Rvecs,i_atom_pair_my)

!CALL write_debug("offsite basbas"//num2str(i_basbas)//num2str(i_basbas-atom_2%offset_basbas)&
!                 //num2str(i_Rvecs)&
!                 //num2str(chksum_2d_matrix(coeff_3fn(:,:, &
!                          i_basbas-atom_2%offset_basbas, 2,i_Rvecs,i_atom_pair_my)))&
!                                  //num2str(chksum_2d_matrix(coeff_one_basbas_part)) )

           endif
         enddo

         ! The on-site coefficients have to be halved
         if(i_cell==1 .OR. i_cell == n_cells +1) coeff_one_basbas_part(:,:) = 0.5d0*coeff_one_basbas_part(:,:)

         CALL MPI_ALLREDUCE(coeff_one_basbas_part, &
                            coeff_one_basbas, size(coeff_one_basbas), &
                            MPI_REAL8, MPI_SUM,&
                            comm_exp_coeffs%mpi_comm , mpierr)

         if(my_basbas_dist%n_offset_rows < i_basbas &
            .AND. &
            i_basbas <= my_basbas_dist%n_offset_rows+my_basbas_dist%n_rows) then
!CALL write_debug("Storing"//num2str(i_basbas)//num2str(chksum_2d_matrix(coeff_one_basbas)) )
                CALL store_cell_compressed(exp_coeffs, i_basbas, i_cell, &
                                           coeff_one_basbas(:,:))
         endif
       enddo
       enddo
    enddo

    DEALLOCATE(Rvecs)
    DEALLOCATE(coeff_one_basbas)
    if(allocated(coeff_3fn)) then
      deallocate(coeff_3fn)
    endif
    if(allocated(coeff_one_basbas_part)) then
      deallocate(coeff_one_basbas_part)
    endif

    call cleanup_lvl_triples()
    call deallocate_tb_auxmat()

    CALL free_exp_coeff_calc_dist(int_dist)
  END SUBROUTINE calc_lvl_tricoeff

  LOGICAL FUNCTION is_my_atom(i_atom_1,i_atom_2)
    INTEGER, INTENT(IN) :: i_atom_1,i_atom_2

    INTEGER :: i

    i=i_atom_1 * n_atoms + i_atom_2

    if(mod(i,n_tasks) == myid) then
        is_my_atom = .TRUE.
        return
    endif

    is_my_atom = .FALSE.
  END FUNCTION

   SUBROUTINE set_offsets(i_atom, offsets)
      INTEGER, INTENT(IN) :: i_atom
      type(basis_and_species_offsets) :: offsets

      offsets%i_species = species(i_atom)

      offsets%offset_basis = atom2basis_off(i_atom)
      offsets%extend_basis = sp2n_basis_sp(offsets%i_species)
      offsets%start_basis = atom2basis_off(i_atom) + 1
      offsets%end_basis = atom2basis_off(i_atom) + sp2n_basis_sp(offsets%i_species)

      offsets%offset_basbas = atom2basbas_off(i_atom)
      offsets%extend_basbas = sp2n_basbas_sp(offsets%i_species)
      offsets%start_basbas = atom2basbas_off(i_atom) + 1
      offsets%end_basbas = atom2basbas_off(i_atom) + sp2n_basbas_sp(offsets%i_species)

   END SUBROUTINE set_offsets

    subroutine calc_lvl_tricoeff_old(exp_coeffs,my_basbas_dist)

    use dimensions
    use prodbas
    use pbc_lists
    use geometry
    use localorb_io
    implicit none


    type(expansion_coefficients) :: exp_coeffs
    type(basbas_distribution),INTENT(IN) :: my_basbas_dist


    ! Local variabales
    integer :: i_atom_1, i_atom_2, i_basbas
    integer :: i_cell_1, i_cell_2, i_cell_3
    integer :: i_cell
    integer :: mpierr
    integer :: info

    real*8 :: Dvec(3), Cvec(3), sum_radius
    real*8 :: Rvecs(3,n_cells_task), dummy(1,1,1,1,1,1)
    real*8, allocatable :: coeff_3fn(:,:,:,:,:), coeff_3fn_tmp(:,:,:,:)
    REAL*8, DIMENSION(n_basis,n_basis,max_n_basbas_sp) :: lvl_tricoeff_cell, &
                                                          lvl_tricoeff_cell_part, &
                                                          lvl_tricoeff_offsite, &
                                                          lvl_tricoeff_offsite_part
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_recip'


    type(basis_and_species_offsets) :: atom_1, atom_2

    write(info_str,'(2X,A)') "Computing the triple LVL expansion coefficents in real space"
    call write_stdout(info_str)

!    CALL initialize_localized_basbas()
    call initialize_lvl_triples(OVLP_TYPE_COULOMB)
    call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)


    if(.not. allocated(coeff_3fn)) then
      allocate(coeff_3fn(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2, 1), stat=info)
      call check_allocation(info, 'coeff_3fn', func)
    endif
    if(.not. allocated(coeff_3fn_tmp)) then
      allocate(coeff_3fn_tmp(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2), stat=info)
      call check_allocation(info, 'coeff_3fn_tmp', func)
    endif


    do i_cell = 1, n_cells, 1

        lvl_tricoeff_cell(:,:,:) = 0.d0
        lvl_tricoeff_cell_part(:,:,:) = 0.d0
        lvl_tricoeff_offsite(:,:,:) = 0.d0
        lvl_tricoeff_offsite_part(:,:,:) = 0.d0

        do i_atom_1 = 1, n_atoms, 1
           CALL set_offsets(i_atom_1, atom_1)

           do i_atom_2 = 1, n_atoms, 1
              CALL set_offsets(i_atom_2, atom_2)

              if(.NOT. is_my_atom(i_atom_1,i_atom_2)) CYCLE

                Dvec = coords(:, i_atom_2) - coords(:, i_atom_1)
                i_cell_1 = cell_index(i_cell, 1)
                i_cell_2 = cell_index(i_cell, 2)
                i_cell_3 = cell_index(i_cell, 3)
                ! distance between two unit cells
                Cvec = matmul(lattice_vector, (/i_cell_1, i_cell_2, i_cell_3/))
                Rvecs(:, 1) = Dvec + Cvec


                sum_radius = atom_radius(species(i_atom_1)) &
                             + atom_radius(species(i_atom_2))
               if(sqrt(SUM(Rvecs(:,1)**2)) >=1.3d0* sum_radius) CYCLE


               coeff_3fn = 0 ! Safety only
               call get_pairwise_coeff_3fn(atom_1%i_species, atom_2%i_species, &
                                           1, &
                                           Rvecs, coeff_3fn, dummy, .false.)

               coeff_3fn_tmp(:,:,:,:) = coeff_3fn(:,:,:,:,1)

               ! The on-site coefficients have to be halved
               if(i_cell==1) coeff_3fn_tmp(:,:,:,:) = 0.5*coeff_3fn_tmp(:,:,:,:)

                  lvl_tricoeff_cell_part(atom_1%start_basis:atom_1%end_basis, &
                                         atom_2%start_basis:atom_2%end_basis, &
                                         1:atom_1%extend_basbas) =  &
                  lvl_tricoeff_cell_part(atom_1%start_basis:atom_1%end_basis, &
                                         atom_2%start_basis:atom_2%end_basis, &
                                         1:atom_1%extend_basbas) + &
                      coeff_3fn_tmp(atom_1%start_basis:atom_1%end_basis, &
                                    atom_2%start_basis:atom_2%end_basis, &
                                    1:atom_1%extend_basbas, 1)

            if((i_cell == 1) .and. (i_atom_2 .ne. i_atom_1)) then
               lvl_tricoeff_offsite_part(atom_1%start_basis:atom_1%end_basis, &
                                         atom_2%start_basis:atom_2%end_basis, &
                                         1:atom_2%extend_basbas) =  &
                             coeff_3fn_tmp(:,:,1:atom_2%extend_basbas, 2)

            endif
          enddo
       enddo

       CALL MPI_ALLREDUCE(lvl_tricoeff_cell_part, &
                          lvl_tricoeff_cell, size(lvl_tricoeff_cell_part), &
                          MPI_REAL8, MPI_SUM, mpi_comm_global , mpierr)

       do i_basbas = 1,n_basbas! max_n_basis_sp
          if(my_basbas_dist%n_offset_rows < i_basbas &
             .AND. &
             i_basbas <= my_basbas_dist%n_offset_rows+my_basbas_dist%n_rows) then
                CALL store_cell_compressed(exp_coeffs, i_basbas, i_cell, &
                                           lvl_tricoeff_cell(:,:,i_basbas))
          endif
       enddo

       if(i_cell == 1) then
         CALL MPI_ALLREDUCE(lvl_tricoeff_offsite_part, &
                            lvl_tricoeff_offsite, size(lvl_tricoeff_offsite_part), &
                            MPI_REAL8, MPI_SUM, mpi_comm_global , mpierr)
       endif
    enddo

    do i_basbas = 1,n_basbas! max_n_basis_sp
          if(my_basbas_dist%n_offset_rows < i_basbas &
             .AND. &
             i_basbas <= my_basbas_dist%n_offset_rows+my_basbas_dist%n_rows) then
               CALL store_cell_compressed(exp_coeffs, i_basbas, i_cell, &
                                          lvl_tricoeff_offsite(:,:,i_basbas))
          endif
    enddo


    if(allocated(coeff_3fn)) then
      deallocate(coeff_3fn)
    endif
    if(allocated(coeff_3fn_tmp)) then
      deallocate(coeff_3fn_tmp)
    endif

    call cleanup_lvl_triples()
    call deallocate_tb_auxmat()

  END SUBROUTINE calc_lvl_tricoeff_old

  INTEGER FUNCTION basbas_sp_offset_by_species(i_species)
     INTEGER, INTENT(IN) :: i_species
     if(i_species > 1) then
        basbas_sp_offset_by_species = SUM(sp2n_basbas_sp(1:i_species))
     elseif(i_species == 1) then
        basbas_sp_offset_by_species = 0
     endif
  END FUNCTION

  INTEGER FUNCTION basbas_to_basbas_sp(i_basbas)
     INTEGER, INTENT(IN) :: i_basbas
     INTEGER :: i_species, &
                i_atom, &
                i_basbas_sp_offset, &
                i_basbas_sp_pos

     do i_atom = 1, n_atoms
        if(atom2basbas_off(i_atom)+1>=i_basbas) EXIT
     enddo

     i_species = species(i_atom)

     i_basbas_sp_offset = basbas_sp_offset_by_species(i_species)
     i_basbas_sp_pos = i_basbas - atom2basbas_off(i_atom)

     basbas_to_basbas_sp = i_basbas_sp_offset + &
                           i_basbas_sp_pos

  END FUNCTION

END MODULE cRPA_calculation_exp_coeffs
