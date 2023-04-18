!****h* FHI-aims/localized_basbas
!  NAME
!    localized_basbas
!  SYNOPSIS

module localized_basbas

  !  PURPOSE
  !
  !     Provide utilities for a expansion of basis function products
  !     in auxiliary (product) basis functions which are localized close
  !     to the corresponding atoms.
  !
  !  DESIGN
  !
  !     Make use of locality as much as possible.
  !
  !     To this end, sparse tensors (type(sp_ten)) are used.  Some quantities
  !     do have some kind of symmetries (ovlp3fn is symmetric with respect to
  !     i_basis_1 and i_basis_2, coulomb_matr is symmetric with respect to
  !     i_basbas_1 and i_basbas_2).  In these cases, only the upper right
  !     triangle is saved and to be used.
  !
  !  USES

  use sparse_tensor
  use dimensions
  use basis
  use geometry
  use prodbas
  use mpi_tasks
  use synchronize_mpi_basic
  use scalapack_utils
  use runtime_choices
  use localorb_io
  use lvl_triples
  implicit none

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
  !    Release version, FHI-aims (2010).
  !
  !  SOURCE

  ! Atom pairs this proc is responsible for:
  integer :: n_local_atom_pairs
  integer, allocatable :: local_atom_pairs(:,:)

  ! Maximum number of pairs of basis function for a given pair of atoms
  integer :: n_max_basis_pair
  ! Maximum number of atomic blocks of aux functions for a given pair of atoms
  integer :: n_max_basbas_atom
  ! Maximum number of aux functions for a given pair of atoms
  integer :: n_max_basbas

  ! Parameter for joining or splitting rows of basis pairs
  real*8, parameter :: frac_split_glob = 0.9d0

contains

  !----------------------------------------------------------------------------
  !****s* localized_basbas/initialize_localized_basbas
  !  NAME
  !    initialize_localized_basbas
  !  SYNOPSIS

  subroutine initialize_localized_basbas()

    !  PURPOSE
    !     Allocate and initialize module data structures
    !  USES

    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: info, i_atom
    character(*), parameter :: func = 'initialize_localized_basbas'

    ! Count local_atom_pairs
    call loop_over_atom_pairs(n_local_atom_pairs)
    ! add by igor for use_dftpt2_and_hse or use_gw_and_hse
    if (allocated(local_atom_pairs)) deallocate(local_atom_pairs) 
    ! Get local_atom_pairs
    allocate(local_atom_pairs(2, n_local_atom_pairs), stat=info)
    call check_allocation(info, 'local_atom_pairs')
    call loop_over_atom_pairs(n_local_atom_pairs, local_atom_pairs)

    ! Size hints
    n_max_basbas_atom = 2
    n_max_basis_pair = maxval(sp2n_basis_sp)**2
    n_max_basbas = n_max_basbas_atom * maxval(sp2n_basbas_sp)
    
    return

  contains

    subroutine loop_over_atom_pairs(n_local_atom_pairs, local_atom_pairs)
      implicit none
      integer, intent(INOUT) :: n_local_atom_pairs
      integer, intent(INOUT), optional :: local_atom_pairs(2,n_local_atom_pairs)
      ! Figure out atom pairs to be stored on this proc
      ! if (n_local_atom_pairs == 0) then
      !    Only count
      ! else
      !    Fill local_atom_pairs data structure
      ! end if

      integer :: i_species
      integer :: i_atom_1, i_species_1, i_atom_2, i_species_2
      integer :: n_local_atom_pairs_uptonow, n_atom_pairs_uptonow
      real*8 :: diff_vec(3), atom_distance, atom_radius_sum

      n_local_atom_pairs_uptonow = 0
      n_atom_pairs_uptonow = 0
      do i_atom_1 = 1, n_atoms
         i_species_1 = species(i_atom_1)
         do i_atom_2 = i_atom_1, n_atoms
            i_species_2 = species(i_atom_2)
            atom_radius_sum = atom_radius(i_species_1) + &
            &                 atom_radius(i_species_2)
            diff_vec = coords(:, i_atom_2) - coords(:, i_atom_1)
            atom_distance = sqrt(dot_product(diff_vec, diff_vec))

            if (atom_distance < atom_radius_sum) then
               n_atom_pairs_uptonow = n_atom_pairs_uptonow + 1
               if (myid == modulo(n_atom_pairs_uptonow, n_tasks)) then
                  n_local_atom_pairs_uptonow = n_local_atom_pairs_uptonow + 1
                  if (present(local_atom_pairs)) then
                     if (n_local_atom_pairs_uptonow > n_local_atom_pairs) then
                        call aims_stop('n_local_atom_pairs too small', func)
                     end if
                     local_atom_pairs(1, n_local_atom_pairs_uptonow) = i_atom_1
                     local_atom_pairs(2, n_local_atom_pairs_uptonow) = i_atom_2
                  end if
               end if
            end if

         end do   ! i_atom_2
      end do   ! i_atom_1

      if (present(local_atom_pairs)) then
         if (n_local_atom_pairs /= n_local_atom_pairs_uptonow) then
            call aims_stop('n_local_atom_pairs too large', func)
         end if
      else
         n_local_atom_pairs = n_local_atom_pairs_uptonow
      end if

    end subroutine loop_over_atom_pairs

  end subroutine initialize_localized_basbas
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/cleanup_localized_basbas
  !  NAME
  !    cleanup_localized_basbas
  !  SYNOPSIS

  subroutine cleanup_localized_basbas()

    !  PURPOSE
    !    Deallocate module data structures
    !  USES

    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'cleanup_localized_basbas'

    n_local_atom_pairs = 0
    if (allocated(local_atom_pairs)) deallocate(local_atom_pairs)

  end subroutine cleanup_localized_basbas
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/ovlp3fn_to_coeff3fn_ten
  !  NAME
  !    ovlp3fn_to_coeff3fn_ten
  !  SYNOPSIS

  subroutine ovlp3fn_to_coeff3fn_ten(ovlp3fn_ten, coulomb_ten)

    !  PURPOSE
    !
    !     Loop over the atom pairs and construct the LVL auxiliary expansion
    !     coefficients from the triple- and double-Coulombs.
    !
    !  USES

    use localorb_io, only: use_unit
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(INOUT) :: ovlp3fn_ten
    type(sp_ten), intent(IN) :: coulomb_ten

    !  INPUTS
    !    o ovlp3fn_ten -- [full] Sparse tensor of all 3-centers needed
    !                     for the atom pairs of this proc
    !    o coulomb_ten -- [full] Sparse tensor of all 2-centers needed
    !                     for the atom pairs of this proc
    !  OUTPUTS
    !    o ovlp3fn_ten -- [full] Sparse tensor of all auxiliary expansion
    !                     coefficients of the atom pairs on this proc
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_atom_pair, i_atom_1, i_atom_2, i_basbas_atom, i_atom_3
    integer :: n_basbas_atoms, basbas_atoms(n_atoms)
    integer, allocatable :: row2basbas(:), col2basis_pair(:,:)
    integer :: n_row_uptonow, n_col_uptonow, n_row, n_col
    integer :: basis_off_1, basis_off_2, n_basis_1, n_basis_2
    integer :: i_basis_1, i_basis_2, i_basbas, i_row, i_col
    real*8, allocatable :: loc_ovlp3fn(:,:), loc_coulomb(:,:)
    real*8, allocatable :: loc_ovlp3fn_tmp(:,:)
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'ovlp3fn_to_coeff3fn_ten'

    allocate(row2basbas(n_max_basbas), stat=info)
    call check_allocation(info, 'row2basbas', func)
    allocate(col2basis_pair(2, n_max_basis_pair), stat=info)
    call check_allocation(info, 'col2basis_pair', func)
    allocate(loc_ovlp3fn(n_max_basbas, n_max_basis_pair), stat=info)
    call check_allocation(info, 'loc_ovlp3fn', func)
    allocate(loc_coulomb(n_max_basbas, n_max_basbas), stat=info)
    call check_allocation(info, 'loc_ovlp3fn', func)

    do i_atom_pair = 1, n_local_atom_pairs
       i_atom_1 = local_atom_pairs(1, i_atom_pair)
       i_atom_2 = local_atom_pairs(2, i_atom_pair)

       ! Get row2basbas
       call get_basbas_atoms_for_pair(i_atom_1, i_atom_2, &
       &                              n_basbas_atoms, basbas_atoms)
       n_row_uptonow = 0
       do i_basbas_atom = 1, n_basbas_atoms
          i_atom_3 = basbas_atoms(i_basbas_atom)
          do i_basbas = 1, sp2n_basbas_sp(species(i_atom_3))
             n_row_uptonow = n_row_uptonow + 1
             if (n_row_uptonow > n_max_basbas) then
                call aims_stop('row mismatch (basbas)', func)
             end if
             row2basbas(n_row_uptonow) = atom2basbas_off(i_atom_3) + i_basbas
          end do
       end do
       n_row = n_row_uptonow

       ! Fill upper (right) triangle of loc_coulomb(:,:)
       call extract_matricization(coulomb_ten, coulomb_ten%val, loc_coulomb, &
       &                          1, (/1/), n_row, row2basbas, &
       &                          1, (/2/), n_row, row2basbas)

       ! Get col2basis_pair
       basis_off_1 = atom2basis_off(i_atom_1)
       basis_off_2 = atom2basis_off(i_atom_2)
       n_basis_1 = sp2n_basis_sp(species(i_atom_1))
       n_basis_2 = sp2n_basis_sp(species(i_atom_2))
       n_col_uptonow = 0
       do i_basis_2 = basis_off_2 + 1, basis_off_2 + n_basis_2
          do i_basis_1 = basis_off_1 + 1, basis_off_1 + n_basis_1
             if (i_basis_1 > i_basis_2) cycle
             if (basis_nghr(i_basis_1, i_basis_2) > 0) then
                n_col_uptonow = n_col_uptonow + 1
                if (n_col_uptonow > n_max_basis_pair) then
                   call aims_stop('col mismatch (basis_pair)', func)
                end if
                col2basis_pair(1, n_col_uptonow) = i_basis_1
                col2basis_pair(2, n_col_uptonow) = i_basis_2
             end if
          end do
       end do
       n_col = n_col_uptonow

       ! Fill loc_ovlp3fn(:,:)
       call extract_matricization(ovlp3fn_ten, ovlp3fn_ten%val, loc_ovlp3fn,&
       &                          1, (/3/),    n_row, row2basbas, &
       &                          2, (/1, 2/), n_col, col2basis_pair)

       ! Solve!
       call dposv('U', n_row, n_col, &
       &          loc_coulomb, n_max_basbas, &
       &          loc_ovlp3fn, n_max_basbas, info)
       if (info /= 0) then
          write(info_str, "('* ',A,': dposv: info =',I7,'; N =',I7,A,2I5)") &
          & trim(func), info, n_row, '; for atom pair', i_atom_1, i_atom_2
          write(use_unit, "(A)") trim(info_str)

          ! Fill upper (right) triangle of loc_coulomb(:,:)
          call extract_matricization(coulomb_ten, coulomb_ten%val, loc_coulomb,&
          &                          1, (/1/), n_row, row2basbas, &
          &                          1, (/2/), n_row, row2basbas)
          do i_row = 1, n_row
             do i_col = i_row+1, n_row  ! sic
                loc_coulomb(i_col, i_row) = loc_coulomb(i_row, i_col)
             end do
          end do
          ! Fill loc_ovlp3fn(:,:) again
          call extract_matricization(ovlp3fn_ten, ovlp3fn_ten%val, &
          &                                                    loc_ovlp3fn,&
          &                          1, (/3/),    n_row, row2basbas, &
          &                          2, (/1, 2/), n_col, col2basis_pair)

          ! Least squares:
          call power_genmat_lapack(n_row, loc_coulomb(1:n_row, 1:n_row), &
          &                        -1.d0, safe_minimum, &
          &                        prodbas_threshold, 'local Coulomb')

          allocate(loc_ovlp3fn_tmp(n_row, n_col), stat=info)
          call check_allocation(info, 'loc_ovlp3fn_tmp', func)
          call dgemm('N', 'N', n_row, n_col, n_row, &
          &          1.d0, loc_coulomb, n_max_basbas, &
          &                loc_ovlp3fn, n_max_basbas, &
          &          0.d0, loc_ovlp3fn_tmp, n_row)
          loc_ovlp3fn(1:n_row, 1:n_col) = loc_ovlp3fn_tmp
          deallocate(loc_ovlp3fn_tmp)
       end if

       ! Put back loc_ovlp3fn into ovlp3fn_ten
       call update_matricization(ovlp3fn_ten, ovlp3fn_ten%val, loc_ovlp3fn, &
       &                         1, (/3/),    n_row, row2basbas, &
       &                         2, (/1, 2/), n_col, col2basis_pair)

    end do

    deallocate(row2basbas, col2basis_pair, loc_ovlp3fn, loc_coulomb)

  end subroutine ovlp3fn_to_coeff3fn_ten
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/get_coeff3fn_ten
  !  NAME
  !    get_coeff3fn_ten
  !  SYNOPSIS

  subroutine get_coeff3fn_ten(coeff3fn_ten)

    !  PURPOSE
    !
    !     Loop over the atom pairs and construct the LVL auxiliary expansion
    !     coefficients by calls to get_pairwise_coeff_3fn().
    !
    !  USES

    use tight_binding_auxmat
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(INOUT) :: coeff3fn_ten

    !  INPUTS
    !    o coeff3fn_ten -- [>=blocking] Blocking of auxiliary expansion
    !                                   coefficients
    !  OUTPUTS
    !    o coeff3fn_ten -- [full] Sparse tensor of all auxiliary expansion
    !                             coefficients of the atom pairs on this proc
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: ovlp_type
    integer :: i_species_1, i_species_2, i_atom_1, i_atom_2
    integer :: i_atom_pair, i_pair_last, i_atom_a, i_atom_b, i_atom_bb
    real*8, allocatable :: coeff3fn_calc(:,:,:,:,:)
    integer :: i_Rvec, n_Rvec, n_Rvec_uptonow, this_n_Rvec
    integer, parameter :: max_n_Rvec = 20
    integer, parameter :: max_buf_size = 2**27   ! 128 MiB
    real*8 :: size_per_Rvec
    integer :: Rvec2atom(2, max_n_Rvec)
    logical :: direct_hit, transposed_hit, is_onsite
    logical :: is_transposed(max_n_Rvec)
    real*8 :: Rvec(3, max_n_Rvec), dummy(1,1,1,1,1,1)
    integer :: n_side, i_side
    integer :: n_arr, arr_top(3), arr_shp(3), arr_str(3)
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_coeff3fn_ten'

    if (use_hse .and. (hse_omega_hf /= 0.d0).and. .not. use_gw_and_hse &
        .and. .not. use_dftpt2_and_hse) then
       ovlp_type = OVLP_TYPE_HSE
    else
       ovlp_type = OVLP_TYPE_COULOMB
    end if
    call initialize_lvl_triples(ovlp_type)
    call initialize_tb_auxmat(1, ovlp_type)

    ! figure out n_Rvec
    size_per_Rvec = 8.d0 * max_n_basis_sp**2 * max_n_basbas_sp * 2
    n_Rvec = max(1, min(max_n_Rvec, ceiling(max_buf_size / size_per_Rvec) ))

    call upgrade_blocking(coeff3fn_ten, SPTEN_FULL, 'coeff3fn')

    write(info_str, "(2X,A,': ',A,I7,A)") &
    & func, 'Temporary array needs', &
    & ceiling(size_per_Rvec * n_Rvec / 2.d0**20), ' MiB buffer space on CPU 0.'
    call localorb_info(info_str)
    allocate(coeff3fn_calc(max_n_basis_sp, max_n_basis_sp, &
    &                      max_n_basbas_sp, 2, n_Rvec), stat=info)
    call check_allocation(info, 'coeff3fn_calc', func)
    arr_str(3) = max_n_basis_sp**2
    n_arr = max_n_basis_sp**2 * max_n_basbas_sp 

    do i_species_1 = 1, n_species
       do i_species_2 = 1, i_species_1

          i_pair_last = 0

          PAIR_BLOCKS_L: do

             ! Search for pairs in main list
             n_Rvec_uptonow = 0
             do i_atom_pair = i_pair_last+1, n_local_atom_pairs
                i_atom_a = local_atom_pairs(1, i_atom_pair)
                i_atom_b = local_atom_pairs(2, i_atom_pair)
                direct_hit = (species(i_atom_a) == i_species_1 .and. &
                &             species(i_atom_b) == i_species_2)
                transposed_hit = (species(i_atom_a) == i_species_2 .and. &
                &                 species(i_atom_b) == i_species_1)
                if (direct_hit .or. transposed_hit) then
                   i_pair_last = i_atom_pair
                   n_Rvec_uptonow = n_Rvec_uptonow + 1
                   Rvec(:, n_Rvec_uptonow) = coords(:, i_atom_b) &
                   &                         - coords(:, i_atom_a)
                   if (direct_hit) then
                      is_transposed(n_Rvec_uptonow) = .false.
                      Rvec2atom(1, n_Rvec_uptonow) = i_atom_a
                      Rvec2atom(2, n_Rvec_uptonow) = i_atom_b
                      ! Could also be onsite (-> also transposed_hit == .true.)
                   else
                      is_transposed(n_Rvec_uptonow) = .true.
                      Rvec(:, n_Rvec_uptonow) = - Rvec(:, n_Rvec_uptonow)
                      Rvec2atom(1, n_Rvec_uptonow) = i_atom_b
                      Rvec2atom(2, n_Rvec_uptonow) = i_atom_a
                   end if
                   ! Rvec2atom stores as calculated /= local_atom_pairs
                end if
                if (n_Rvec_uptonow >= n_Rvec) exit
             end do
             if (n_Rvec_uptonow == 0) exit PAIR_BLOCKS_L
             this_n_Rvec = n_Rvec_uptonow

             ! Calculate a bunch of pairs
             call get_pairwise_coeff_3fn(i_species_1, i_species_2, &
             &                           this_n_Rvec, Rvec, coeff3fn_calc, dummy, .false.)

             ! Store them into coeff3fn_ten
             do i_Rvec = 1, this_n_Rvec
                is_onsite = (all(abs(Rvec(:, i_Rvec)) < 1d-10))
                if (is_onsite) then
                   n_side = 1
                else
                   n_side = 2
                end if
                i_atom_1 = Rvec2atom(1, i_Rvec)
                i_atom_2 = Rvec2atom(2, i_Rvec)

                ! Now the difficult part.
                ! The indexing in coeff3fn_ten is according to atom_pair.
                ! In case of (is_transposed(i_Rvec)), this does not match
                ! the order in coeff3fn_calc.  Fix this by using non-default
                ! strides for the array coeff3fn_calc.
                if (is_transposed(i_Rvec)) then
                   i_atom_a = i_atom_2
                   i_atom_b = i_atom_1
                   arr_str(1) = max_n_basis_sp
                   arr_str(2) = 1
                else
                   i_atom_a = i_atom_1
                   i_atom_b = i_atom_2
                   arr_str(1) = 1
                   arr_str(2) = max_n_basis_sp
                end if
                arr_top(1) = atom2basis_off(i_atom_a) + 1
                arr_top(2) = atom2basis_off(i_atom_b) + 1
                arr_shp(1) = sp2n_basis_sp(species(i_atom_a))
                arr_shp(2) = sp2n_basis_sp(species(i_atom_b))

                do i_side = 1, n_side
                   ! Use atom numbering for coeff3fn_calc to get i_atom_bb.
                   if (i_side == 1) then
                      i_atom_bb = i_atom_1
                   else
                      i_atom_bb = i_atom_2
                   end if
                   arr_top(3) = atom2basbas_off(i_atom_bb) + 1
                   arr_shp(3) = sp2n_basbas_sp(species(i_atom_bb))
                   call update_block(coeff3fn_ten, coeff3fn_ten%val, n_arr, &
                   &                 coeff3fn_calc(:,:,:, i_side, i_Rvec), &
                   &                 arr_top, arr_shp, arr_str)
                end do

             end do
          end do PAIR_BLOCKS_L
       end do
    end do

    deallocate(coeff3fn_calc)
    call cleanup_lvl_triples()
    call deallocate_tb_auxmat()

  end subroutine get_coeff3fn_ten
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/construct_LVL_3fn_want
  !  NAME
  !    construct_LVL_3fn_want
  !  SYNOPSIS

  subroutine construct_LVL_3fn_want(want)

    !  PURPOSE
    !
    !     Construct (local) want of ovlp_3fn_ten / coeff_3fn_ten.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(OUT) :: want

    !  INPUTS
    !    none
    !  OUTPUTS
    !    o want -- Sparse blocking naming the parts of ovlp_3fn needed on this
    !              processor.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_atom_pair, i_atom_1, i_atom_2, i_atom_3
    integer :: n_basbas_atoms, basbas_atoms(n_atoms)
    integer :: n_max_nzb, n_nzb_pair
    integer :: i_nzb_pair, i_nzb_basbas, top_basbas, shp_basbas
    integer :: n_nzb_uptonow
    integer, allocatable :: top(:,:), shp(:,:)
    integer :: glb_shp(3)
    integer :: info
    character(*), parameter :: func = 'construct_LVL_3fn_want'

    n_max_nzb = maxval(n_basis_fn_species)
    allocate(top(2, n_max_nzb), shp(2, n_max_nzb), stat=info)
    call check_allocation(info, 'top/shp')
    glb_shp(1:2) = n_basis
    glb_shp(3) = n_basbas

    ! Count number of non-zero blocks
    n_nzb_uptonow = 0
    do i_atom_pair = 1, n_local_atom_pairs
       i_atom_1 = local_atom_pairs(1, i_atom_pair)
       i_atom_2 = local_atom_pairs(2, i_atom_pair)
       call get_basbas_atoms_for_pair(i_atom_1, i_atom_2, &
       &                              n_basbas_atoms, basbas_atoms)
       call get_basis_pairs_for_pair(i_atom_1, i_atom_2, n_max_nzb, &
       &                             top, shp, n_nzb_pair, frac_split_glob)

       n_nzb_uptonow = n_nzb_uptonow + n_nzb_pair * n_basbas_atoms
    end do

    ! Get non-zero blocks
    call alloc_sp_ten(want, 3, glb_shp, SPTEN_BLOCKING, n_nzb_uptonow, 0, &
    &                 'LVL-3fn-want')
    n_nzb_uptonow = 0
    do i_atom_pair = 1, n_local_atom_pairs
       i_atom_1 = local_atom_pairs(1, i_atom_pair)
       i_atom_2 = local_atom_pairs(2, i_atom_pair)
       call get_basbas_atoms_for_pair(i_atom_1, i_atom_2, &
       &                              n_basbas_atoms, basbas_atoms)
       call get_basis_pairs_for_pair(i_atom_1, i_atom_2, n_max_nzb, &
       &                             top, shp, n_nzb_pair, frac_split_glob)

       do i_nzb_basbas = 1, n_basbas_atoms
          i_atom_3 = basbas_atoms(i_nzb_basbas)
          top_basbas = atom2basbas_off(i_atom_3) + 1
          shp_basbas = sp2n_basbas_sp(species(i_atom_3))
          do i_nzb_pair = 1, n_nzb_pair
             n_nzb_uptonow = n_nzb_uptonow + 1
             want%top(1:2, n_nzb_uptonow) = top(:, i_nzb_pair)
             want%shp(1:2, n_nzb_uptonow) = shp(:, i_nzb_pair)
             want%top(3, n_nzb_uptonow) = top_basbas
             want%shp(3, n_nzb_uptonow) = shp_basbas
          end do
       end do

    end do   ! i_atom_pair

  end subroutine construct_LVL_3fn_want
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/construct_LVL_coulomb_want
  !  NAME
  !    construct_LVL_coulomb_want
  !  SYNOPSIS

  subroutine construct_LVL_coulomb_want(want)

    !  PURPOSE
    !
    !     Construct (local) want of coulomb_matr
    !
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(OUT) :: want

    !  INPUTS
    !    none
    !  OUTPUTS
    !    o want -- Sparse blocking naming the parts of coulomb_matr needed
    !              on this processor.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_atom_pair, i_atom_1, i_atom_2
    logical :: pair_needed(n_atoms, n_atoms)
    integer :: n_basbas_atoms, basbas_atoms(n_atoms)
    integer :: i, j
    integer :: glb_shp(2), top(2), shp(2)
    integer :: n_nzb, n_nzb_uptonow
    character(*), parameter :: func = 'construct_LVL_coulomb_want'

    ! Get pair_needed
    pair_needed = .false.
    do i_atom_pair = 1, n_local_atom_pairs
       i_atom_1 = local_atom_pairs(1, i_atom_pair)
       i_atom_2 = local_atom_pairs(2, i_atom_pair)

       call get_basbas_atoms_for_pair(i_atom_1, i_atom_2, &
       &                              n_basbas_atoms, basbas_atoms)
       do i = 1, n_basbas_atoms
          do j = i, n_basbas_atoms
             if (basbas_atoms(i) <= basbas_atoms(j)) then
                pair_needed(basbas_atoms(i), basbas_atoms(j)) = .true.
             else
                pair_needed(basbas_atoms(j), basbas_atoms(i)) = .true.
             end if
          end do
       end do
    end do

    n_nzb = count(pair_needed)
    glb_shp(1:2) = n_basbas
    call alloc_sp_ten(want, 2, glb_shp, SPTEN_BLOCKING, n_nzb, 0, &
    &                 'LVL-2fn-want')

    n_nzb_uptonow = 0
    do i_atom_1 = 1, n_atoms
       top(1) = atom2basbas_off(i_atom_1) + 1
       shp(1) = sp2n_basbas_sp(species(i_atom_1))
       do i_atom_2 = i_atom_1, n_atoms
          if (pair_needed(i_atom_1, i_atom_2)) then
             top(2) = atom2basbas_off(i_atom_2) + 1
             shp(2) = sp2n_basbas_sp(species(i_atom_2))
             n_nzb_uptonow = n_nzb_uptonow + 1
             want%top(:, n_nzb_uptonow) = top
             want%shp(:, n_nzb_uptonow) = shp
          end if
       end do
    end do

  end subroutine construct_LVL_coulomb_want
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/not_so_sparse_ovlp3fn
  !  NAME
  !    not_so_sparse_ovlp3fn
  !  SYNOPSIS

  subroutine not_so_sparse_ovlp3fn(ovlp_3fn, ten)

    !  PURPOSE
    !
    !    Copy the prodbas-distributed ovlp_3fn to a corresponding tensor.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)
    type(sp_ten), intent(OUT) :: ten

    !  INPUTS
    !    o ovlp_3fn -- prodbas distributed triple-overlap
    !  OUTPUTS
    !    o ten -- Same content, same distribution, different storage scheme
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer, allocatable :: basis_pairs(:,:)
    integer :: info
    character(*), parameter :: func = 'not_so_sparse_ovlp3fn'

    ! --- Get blockings for pairs and auxs

    call get_ovlp3fn_basbas_blocking(ten)
    call upgrade_blocking(ten, SPTEN_FULL, '3fn-1D')
    call check_local_aliasing(ten, 'ten', func)

    ! --- Get basis_pairs

    allocate(basis_pairs(2, n_basis_pairs), stat=info)
    call check_allocation(info, 'basis_pairs', func)
    call get_basis_pair_list(basis_pairs)

    ! --- Get values

    call update_matricization(ten, ten%val, ovlp_3fn, &
    &                         2, (/1, 2/), n_basis_pairs, basis_pairs, &
    &                         1, (/3/), n_loc_prodbas, map_prodbas(:, myid+1))

    deallocate(basis_pairs)

  end subroutine not_so_sparse_ovlp3fn
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/back_to_ovlp3fn
  !  NAME
  !    back_to_ovlp3fn
  !  SYNOPSIS

  subroutine back_to_ovlp3fn(ovlp_3fn, ten)

    !  PURPOSE
    !
    !    Copy the contents of the tensor ten to ovlp_3fn.  Only the upper
    !    triangle ('U') of ten is used (see get_basis_pair_list()).
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(OUT) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)
    type(sp_ten), intent(IN) :: ten

    !  INPUTS
    !    o ovlp_3fn -- prodbas distributed triple-overlap
    !  OUTPUTS
    !    o ten -- Same content, same distribution, different storage scheme
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer, allocatable :: basis_pairs(:,:)
    integer :: info
    character(*), parameter :: func = 'back_to_ovlp3fn'

    ! --- Get basis_pairs

    allocate(basis_pairs(2, n_basis_pairs), stat=info)
    call check_allocation(info, 'basis_pairs', func)
    call get_basis_pair_list(basis_pairs)

    ! --- Get values

    call extract_matricization(ten, ten%val, ovlp_3fn, &
    &                          2, (/1, 2/), n_basis_pairs, basis_pairs, &
    &                          1, (/3/), n_loc_prodbas, map_prodbas(:, myid+1))

    deallocate(basis_pairs)

  end subroutine back_to_ovlp3fn
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/global_basis_pairs
  !  NAME
  !    global_basis_pairs
  !  SYNOPSIS

  subroutine global_basis_pairs(basis_pair_ten)

    !  PURPOSE
    !     Get global blocking of basis pairs.
    !     This could also be interesting for ordinary DFT.
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(OUT) :: basis_pair_ten

    !  INPUTS
    !    none
    !  OUTPUTS
    !    o basis_pair_ten -- [blocking] Blocking of overlapping basis funcs
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: n_max_nzb, glb_shp(2)
    integer, allocatable :: top(:,:), shp(:,:)
    integer :: i_atom_pair, i_atom_1, i_atom_2
    integer :: n_local_nzb_uptonow, n_local_nzb, n_nzb, off, i_task
    integer :: i_nzb_pair, n_nzb_pair
    integer :: n_nzb_list(0:n_tasks-1)
    integer :: info
    character(*), parameter :: func = 'global_basis_pairs'

    glb_shp(:) = n_basis
    n_max_nzb = maxval(n_basis_fn_species)
    allocate(top(2, n_max_nzb), shp(2, n_max_nzb), stat=info)
    call check_allocation(info, 'top/shp')

    ! Get n_local_nzb
    n_local_nzb_uptonow = 0
    do i_atom_pair = 1, n_local_atom_pairs
       i_atom_1 = local_atom_pairs(1, i_atom_pair)
       i_atom_2 = local_atom_pairs(2, i_atom_pair)
       call get_basis_pairs_for_pair(i_atom_1, i_atom_2, n_max_nzb, &
       &                             top, shp, n_nzb_pair, frac_split_glob)
       n_local_nzb_uptonow = n_local_nzb_uptonow + n_nzb_pair
    end do
    n_local_nzb = n_local_nzb_uptonow

    ! Synchronize number of blocks
    n_nzb_list = 0
    n_nzb_list(myid) = n_local_nzb
    call sync_int_vector(n_nzb_list, n_tasks)
    n_nzb = sum(n_nzb_list)
    off = 0
    do i_task = 0, myid-1
       off = off + n_nzb_list(i_task)
    end do

    ! Fill locals ...
    call alloc_sp_ten(basis_pair_ten, 2, glb_shp, SPTEN_BLOCKING, n_nzb, 0, &
    &                 'global_basis_pair')
    basis_pair_ten%top = 0
    basis_pair_ten%shp = 0
    do i_atom_pair = 1, n_local_atom_pairs
       i_atom_1 = local_atom_pairs(1, i_atom_pair)
       i_atom_2 = local_atom_pairs(2, i_atom_pair)
       call get_basis_pairs_for_pair(i_atom_1, i_atom_2, n_max_nzb, &
       &                             top, shp, n_nzb_pair, frac_split_glob)
       do i_nzb_pair = 1, n_nzb_pair
          basis_pair_ten%top(:, off+i_nzb_pair) = top(:, i_nzb_pair)
          basis_pair_ten%shp(:, off+i_nzb_pair) = shp(:, i_nzb_pair)
       end do
       off = off + n_nzb_pair
    end do

    ! ... and synchronize them
    call sync_int_vector(basis_pair_ten%top, 2*n_nzb)
    call sync_int_vector(basis_pair_ten%shp, 2*n_nzb)

    deallocate(top, shp)

  end subroutine global_basis_pairs
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/get_basbas_atoms_for_pair
  !  NAME
  !    get_basbas_atoms_for_pair
  !  SYNOPSIS

  subroutine get_basbas_atoms_for_pair(i_atom_1, i_atom_2, &
  &                                    n_basbas_atoms, basbas_atoms)

    !  PURPOSE
    !     For a given pair of atoms, specify whose product basis functions
    !     should be used
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_atom_1, i_atom_2
    integer, intent(OUT) :: n_basbas_atoms
    integer, intent(OUT) :: basbas_atoms(n_atoms)

    !  INPUTS
    !    o atom_1, atom_2 -- Atom numbers of pair
    !  OUTPUTS
    !    o n_basbas_atoms -- Number of atoms whose aux functions we need
    !    o basbas_atoms -- Actual atom numbers
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'get_basbas_atoms_for_pair'


    if (i_atom_1 == i_atom_2) then
       ! As the aux basis is constructed for onsite, that should suffice.
       n_basbas_atoms = 1
       basbas_atoms(1) = i_atom_1
    else
       ! Simple two-center approximation:
       n_basbas_atoms = 2
       basbas_atoms(1) = i_atom_1
       basbas_atoms(2) = i_atom_2
       ! JW: If we are to add the aux funcs of additional atoms, it is
       ! probably best to define a cutoff and take all atoms which are nearer
       ! to *both* corresponding atoms than that cutoff.
    end if
    if (n_basbas_atoms > n_max_basbas_atom) then
       call aims_stop('n_max_basbas_atom too small', func)
    end if

  end subroutine get_basbas_atoms_for_pair
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/get_basis_pairs_for_pair
  !  NAME
  !    get_basis_pairs_for_pair
  !  SYNOPSIS

  subroutine get_basis_pairs_for_pair(i_atom_1, i_atom_2, n_max_nzb, &
  &                                   top, shp, n_nzb, frac_split)

    !  PURPOSE
    !
    !     Figure out a blocking structure for the offsite pairs of basis
    !     functions for a given atom pair.
    !
    !     The current strategy assumes that the basis functions are sorted by
    !     increasing extent for each atom.  For each column, it is figured out
    !     how many rows are needed; if the 
    !     
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_atom_1, i_atom_2
    integer, intent(IN) :: n_max_nzb
    integer, intent(OUT) :: top(2, n_max_nzb), shp(2, n_max_nzb)
    integer, intent(OUT) :: n_nzb
    real*8, intent(IN) :: frac_split

    !  INPUTS
    !    o i_atom_1, i_atom_2 -- Atom numbers of pair
    !    o n_max_nzb -- Maximum number of blocks
    !    o frac_split -- Fraction (between 0 and 1).
    !              Two columns are kept in one block if the shorter one
    !              contains more rows than frac_split times the longer one.
    !  OUTPUTS
    !    o top, shp -- First (top) element and shape (shp) of blocks
    !    o n_nzb -- Number of blocks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: atom_diff(3), atom_dist, radius_1, radius_2
    integer :: basis_off_1, n_basis_1, basis_off_2, n_basis_2
    integer :: i_basis_1, i_basis_2
    integer :: i_fn_1, i_fn_2, n_nzb_uptonow
    logical :: can_be_skipped, add_to_current
    integer :: this_n_basis_1, current_n_basis_1, min_n_basis_1
    character(*), parameter :: func = 'get_basis_pairs_for_pair'

    atom_diff = coords(:, i_atom_2) - coords(:, i_atom_1)
    atom_dist = sqrt(dot_product(atom_diff, atom_diff))

    basis_off_1 = atom2basis_off(i_atom_1)
    n_basis_1 = sp2n_basis_sp(species(i_atom_1))
    basis_off_2 = atom2basis_off(i_atom_2)
    n_basis_2 = sp2n_basis_sp(species(i_atom_2))

    n_nzb_uptonow = 0
    current_n_basis_1 = 0
    min_n_basis_1 = 0
    ! Go backwards (from most extended to most localized) in columns:
    do i_basis_2 = basis_off_2 + n_basis_2, basis_off_2 + 1, -1
       i_fn_2 = basis_fn(i_basis_2)
       radius_2 = outer_radius(i_fn_2)

       ! Skip basis functions which are short-ranged enough
       do i_basis_1 = basis_off_1 + 1, basis_off_1 + n_basis_1
          i_fn_1 = basis_fn(i_basis_1)
          radius_1 = outer_radius(i_fn_1)
          can_be_skipped = (radius_1 + radius_2 < atom_dist)
          if (.not. can_be_skipped) exit
       end do   ! i_basis_1
       this_n_basis_1 = basis_off_1 + n_basis_1 - i_basis_1 + 1

       ! Save
       if (this_n_basis_1 > 0) then
          add_to_current = ((this_n_basis_1 >= min_n_basis_1) .and.&
          &                 (this_n_basis_1 <= current_n_basis_1))
          if (add_to_current) then
             shp(2, n_nzb_uptonow) = shp(2, n_nzb_uptonow) + 1
             top(2, n_nzb_uptonow) = top(2, n_nzb_uptonow) - 1
             if (top(2, n_nzb_uptonow) /= i_basis_2) then
                call aims_stop('i_basis_2 mismatch', func)
             end if
          else
             ! Add new block
             n_nzb_uptonow = n_nzb_uptonow + 1
             if (n_nzb_uptonow > n_max_nzb) then
                call aims_stop('n_max_nzb exceeded', func)
             end if
             top(1, n_nzb_uptonow) = i_basis_1
             top(2, n_nzb_uptonow) = i_basis_2
             shp(1, n_nzb_uptonow) = this_n_basis_1
             shp(2, n_nzb_uptonow) = 1      ! Start with this column only
             current_n_basis_1 = this_n_basis_1
             min_n_basis_1 = floor(frac_split * this_n_basis_1)
             ! JW: Maybe a constant offset would be better...
          end if
       else
          ! Empty column; reset!
          current_n_basis_1 = 0
          min_n_basis_1 = 0
       end if
    end do   ! i_basis_2
    n_nzb = n_nzb_uptonow

  end subroutine get_basis_pairs_for_pair
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/get_basis_pair_list
  !  NAME
  !    get_basis_pair_list
  !  SYNOPSIS

  subroutine get_basis_pair_list(basis_pairs)

    !  PURPOSE
    !    Get an explicit list of the basis pairs from basis_nghr.
    !    Keep invariant i_basis_1 <= i_basis_2 (== lapack 'U').
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(OUT) :: basis_pairs(2, n_basis_pairs)

    !  INPUTS
    !    none
    !  OUTPUTS
    !    basis_pairs -- i_basis_pair -> i_basis_1 <= i_basis_2
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_basis_1, i_basis_2, i_basis_pair
    character(*), parameter :: func = 'get_basis_pair_list'

    basis_pairs = 0
    do i_basis_2 = 1, n_basis
       do i_basis_1 = 1, i_basis_2
          i_basis_pair = basis_nghr(i_basis_1, i_basis_2)
          if (i_basis_pair > 0) then
             if (any(basis_pairs(:, i_basis_pair) > 0)) then
                call aims_stop('Doubled basis pair', func)
             end if
             basis_pairs(1, i_basis_pair) = i_basis_1
             basis_pairs(2, i_basis_pair) = i_basis_2
          end if
       end do
    end do
    if (any(basis_pairs <= 0)) call aims_stop('Missing basis pair', func)

  end subroutine get_basis_pair_list
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/get_ovlp3fn_basbas_blocking
  !  NAME
  !    get_ovlp3fn_basbas_blocking
  !  SYNOPSIS

  subroutine get_ovlp3fn_basbas_blocking(ovlp3fn_bbb)

    !  PURPOSE
    !
    !    Get the blocking of ovlp_3fn roughly(*) corresponding to the
    !    distribution over auxiliary basis functions.
    !
    !    (*) "roughly" because to get a more efficient blocking, the number
    !        of basis pairs is slightly increased.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(OUT) :: ovlp3fn_bbb

    !  INPUTS
    !    none
    !  OUTPUTS
    !    o ovlp3fn_bbb -- [blocking] Description of the local structure of
    !                                ovlp_3fn
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    type(sp_ten) :: basis_pair_ten
    integer :: n_nzb_basbas, n_nzb, n_nzb_uptonow
    integer :: i_nzb_basbas, i_nzb_pair
    integer :: glb_shp(3), this_top(3), this_shp(3)
    integer, allocatable :: basbas2row(:), basbas2col(:)
    integer, allocatable :: col_top(:), col_shp(:), col_pos(:)
    integer :: info
    character(*), parameter :: func = 'get_ovlp3fn_basbas_blocking'

    ! Get blocking of basis pairs
    call global_basis_pairs(basis_pair_ten)

    ! Get (local) blocking of product (auxiliary) basis functions
    allocate(basbas2row(n_basbas), basbas2col(n_basbas), stat=info)
    call check_allocation(info, 'basbas2xxx', func)
    allocate(col_top(n_basbas), col_shp(n_basbas), col_pos(n_basbas), stat=info)
    call check_allocation(info, 'col_xxx', func)
    call get_basbas_to_rowcol(.false., basbas2row, basbas2col)
    call get_topshppos_from_basbas2rowcol(basbas2col, n_nzb_basbas, &
    &                                     col_top, col_shp, col_pos)
    deallocate(basbas2row, basbas2col)
    
    ! Prepare
    glb_shp(1:2) = n_basis
    glb_shp(3) = n_basbas
    n_nzb = basis_pair_ten%n_nzb * n_nzb_basbas
    call alloc_sp_ten(ovlp3fn_bbb, 3, glb_shp, SPTEN_BLOCKING, n_nzb, 0, &
    &                 '3fn-1D-blocking')

    ! Get outer product blocking
    n_nzb_uptonow = 0
    do i_nzb_basbas = 1, n_nzb_basbas
       this_top(3) = col_top(i_nzb_basbas)
       this_shp(3) = col_shp(i_nzb_basbas)
       do i_nzb_pair = 1, basis_pair_ten%n_nzb
          this_top(1:2) = basis_pair_ten%top(:, i_nzb_pair)
          this_shp(1:2) = basis_pair_ten%shp(:, i_nzb_pair)
          n_nzb_uptonow = n_nzb_uptonow + 1
          if (n_nzb_uptonow > ovlp3fn_bbb%n_nzb) then
             call aims_stop('n_nzb mismatch', func)
          end if
          ovlp3fn_bbb%top(:, n_nzb_uptonow) = this_top
          ovlp3fn_bbb%shp(:, n_nzb_uptonow) = this_shp
       end do
    end do
    if (n_nzb_uptonow /= ovlp3fn_bbb%n_nzb) then
       call aims_stop('n_nzb too large', func)
    end if
    call dealloc_sp_ten(basis_pair_ten)

  end subroutine get_ovlp3fn_basbas_blocking
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/get_auxmat_descriptor
  !  NAME
  !    get_auxmat_descriptor
  !  SYNOPSIS

  subroutine get_auxmat_descriptor(distr_2d, auxmat_desc, name)

    !  PURPOSE
    !    Get a sp_ten descriptor for a distributed auxiliary basis matrix.
    !  USES

    implicit none

    !  ARGUMENTS

    logical, intent(IN) :: distr_2d
    type(sp_ten), intent(OUT) :: auxmat_desc
    character(*), intent(IN) :: name

    !  INPUTS
    !    o distr_2d -- In case of use_scalapack, use 1D or 2D distribution?
    !  OUTPUTS
    !    o auxmat_desc -- Descriptor
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer, allocatable :: basbas2row(:), basbas2col(:)
    integer :: n_row_nzb, n_col_nzb, i_row_nzb, i_col_nzb, n_nzb_uptonow
    integer, allocatable :: row_top(:), row_shp(:), row_pos(:)
    integer, allocatable :: col_top(:), col_shp(:), col_pos(:)
    integer :: glb_shp(2)
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_auxmat_descriptor'

    if (use_scalapack) then
       if (distr_2d) then
          call get_scalapack_descriptor(auxmat_desc, &
          & n_basbas, n_basbas, n_basbas, &
          & nb_aux_2d, nb_aux_2d, &
          & nprow_aux_2d, npcol_aux_2d, myprow_aux_2d, mypcol_aux_2d, name)
       else
          call get_scalapack_descriptor(auxmat_desc, &
          & n_basbas, n_basbas, n_basbas, &
          & mb_aux, nb_aux, nprow_aux, npcol_aux, myprow_aux, mypcol_aux, name)
       end if
    else
       allocate(basbas2row(n_basbas), basbas2col(n_basbas), stat=info)
       call check_allocation(info, 'basbas2xxx', func)

       allocate(row_top(n_basbas), row_shp(n_basbas), stat=info)
       call check_allocation(info, 'row_xxx', func)
       allocate(col_top(n_basbas), col_shp(n_basbas), stat=info)
       call check_allocation(info, 'col_xxx', func)
       allocate(row_pos(n_basbas), col_pos(n_basbas), stat=info)
       call check_allocation(info, 'xxx_pos', func)
       
       call get_basbas_to_rowcol(distr_2d, basbas2row, basbas2col)
       call get_topshppos_from_basbas2rowcol(basbas2row, n_row_nzb, &
       &                                     row_top, row_shp, row_pos)
       call get_topshppos_from_basbas2rowcol(basbas2col, n_col_nzb, &
       &                                     col_top, col_shp, col_pos)

       glb_shp(1:2) = n_basbas
       call alloc_sp_ten(auxmat_desc, 2, glb_shp, SPTEN_DESCRIPTOR, &
       &                 n_row_nzb*n_col_nzb, n_basbas*n_loc_prodbas, &
       &                 name, .false.)

       n_nzb_uptonow = 0
       do i_col_nzb = 1, n_col_nzb
          do i_row_nzb = 1, n_row_nzb
             n_nzb_uptonow = n_nzb_uptonow + 1
             auxmat_desc%top(1, n_nzb_uptonow) = row_top(i_row_nzb)
             auxmat_desc%top(2, n_nzb_uptonow) = col_top(i_row_nzb)
             auxmat_desc%shp(1, n_nzb_uptonow) = row_shp(i_col_nzb)
             auxmat_desc%shp(2, n_nzb_uptonow) = col_shp(i_col_nzb)
             auxmat_desc%off(n_nzb_uptonow) = (row_pos(i_row_nzb)-1) &
             &                     + n_basbas*(col_pos(i_col_nzb)-1)
          end do
       end do
       deallocate(row_top, row_shp, row_pos, col_top, col_shp, col_pos)
       deallocate(basbas2row, basbas2col)
    end if
!!$    write(info_str, "(2X,'Having',I6,' non-zero-blocks in auxmat')") &
!!$    & auxmat_desc%n_nzb
!!$    call localorb_allinfo(info_str)

  end subroutine get_auxmat_descriptor
  !******
  !----------------------------------------------------------------------------
  !****s* localized_basbas/get_topshppos_from_basbas2rowcol
  !  NAME
  !    get_topshppos_from_basbas2rowcol
  !  SYNOPSIS

  subroutine get_topshppos_from_basbas2rowcol(basbas2rowcol, &
  &                                           n_nzb, top, shp, pos)

    !  PURPOSE
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: basbas2rowcol(n_basbas)
    integer, intent(OUT) :: n_nzb
    integer, intent(OUT) :: top(n_basbas), shp(n_basbas), pos(n_basbas)

    !  INPUTS
    !    o basbas2rowcol -- Distribution by i_basbas -> i_loc_prodbas
    !  OUTPUTS
    !    o n_nzb -- Number of blocks
    !    o top -- First global index in block
    !    o shp -- Number of entries in block
    !    o pos -- First local index in block
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_basbas, i_rowcol, i_next_rowcol, n_nzb_uptonow
    character(*), parameter :: func = 'get_topshppos_from_basbas2rowcol'

    
    n_nzb_uptonow = 0
    i_next_rowcol = 0
    do i_basbas = 1, n_basbas
       i_rowcol = basbas2rowcol(i_basbas)
       if (i_next_rowcol > 0) then
          ! Within block
          if (i_rowcol == i_next_rowcol) then
             ! Continued block
             i_next_rowcol = i_rowcol + 1
          else
             ! End of block
             shp(n_nzb_uptonow) = i_basbas - top(n_nzb_uptonow)
             i_next_rowcol = 0
          end if
       end if
       if (i_next_rowcol == 0 .and. i_rowcol > 0) then
          ! New block
          n_nzb_uptonow = n_nzb_uptonow + 1
          top(n_nzb_uptonow) = i_basbas
          pos(n_nzb_uptonow) = i_rowcol
          i_next_rowcol = i_rowcol + 1
       end if
    end do
    if (i_next_rowcol > 0) then
       shp(n_nzb_uptonow) = n_basbas - top(n_nzb_uptonow) + 1
       i_next_rowcol = 0
    end if
    n_nzb = n_nzb_uptonow

  end subroutine get_topshppos_from_basbas2rowcol
  !******
end module localized_basbas
!******
