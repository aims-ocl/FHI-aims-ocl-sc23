program test_lvlpairs
  use species_data
  use sparse_tensor
  use basis
  use prodbas
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use localized_basbas
  implicit none
  real*8, allocatable :: coul3fn(:,:)
  real*8, allocatable :: coeff3fn_lvl(:,:)
  real*8, allocatable :: coul2fn(:,:)
  real*8, allocatable :: coul3fn_lvl(:,:)
  type(sp_ten) :: coeff3fn_ten
  integer :: i_basis_1, i_basis_2, i_basis_pair, i_basbas
  integer :: i_atom_1, i_atom_2, i_atom_bb
  integer :: basis_off_1, basis_off_2, basbas_off
  character*150 :: fnname, basis_name_1, basis_name_2, basbas_name
  integer :: n_basis_sp_1, n_basis_sp_2, n_basbas_sp
  integer :: i_basis_sp_1, i_basis_sp_2, i_basbas_sp
  integer :: i_loc_prodbas
  real*8 :: acc_errsq, acc_valsq
  real*8, allocatable :: max_errsq(:,:,:), max_valsq(:,:,:)
  character(*), parameter :: func = 'test_lvlpairs'

  ! --- Initialize

  call init_aims()
  if (.not. use_prodbas) call aims_stop('No prodbas defined', func)
  if (n_tasks > 1) call aims_stop('Only serial run', func)

  ! --- Calculate direct integrals

  allocate(coul3fn(n_basis_pairs, n_loc_prodbas))
  call integrate_ovlp3fn(l_shell_max, coul3fn, OVLP_TYPE_COULOMB)

  ! --- Calculate LVL approximation to integrals

  allocate(coeff3fn_lvl(n_basis_pairs, n_loc_prodbas))
  allocate(coul2fn(n_loc_prodbas, n_basbas))
  allocate(coul3fn_lvl(n_basis_pairs, n_loc_prodbas))
  call get_coeff_3fn_lvl(coeff3fn_ten, coul2fn)
  call back_to_ovlp3fn(coeff3fn_lvl, coeff3fn_ten)
  coul3fn_lvl = coeff3fn_lvl; call get_v_multi_ovlp3fn(coul2fn, coul3fn_lvl)

  ! --- Analysis

  if (n_tasks > 1) call aims_stop('Only serial run', func)

  allocate(max_errsq(n_atoms, n_atoms, n_atoms)); max_errsq = 0.d0
  allocate(max_valsq(n_atoms, n_atoms, n_atoms)); max_valsq = 0.d0

  if (myid == 0) open(22, file="errors.dat")
  ATOM_1: do i_atom_1 = 1, n_atoms
     basis_off_1 = atom2basis_off(i_atom_1)
     n_basis_sp_1 = sp2n_basis_sp(species(i_atom_1))
     ATOM_2: do i_atom_2 = 1, n_atoms
        basis_off_2 = atom2basis_off(i_atom_2)
        n_basis_sp_2 = sp2n_basis_sp(species(i_atom_2))

        BASIS_1: do i_basis_sp_1 = 1, n_basis_sp_1
           i_basis_1 = basis_off_1 + i_basis_sp_1
           fnname = basis_fnname(basis_fn(i_basis_1))
           write(basis_name_1, "(A,I3)") trim(fnname), basis_m(i_basis_1)
           BASIS_2: do i_basis_sp_2 = 1, n_basis_sp_2
              i_basis_2 = basis_off_2 + i_basis_sp_2
              fnname = basis_fnname(basis_fn(i_basis_2))
              write(basis_name_2, "(A,I3)") trim(fnname), basis_m(i_basis_2)

              i_basis_pair = basis_nghr(i_basis_1, i_basis_2)
              if (i_basis_pair <= 0) cycle

              ATOM_BB: do i_atom_bb = 1, n_atoms
                 basbas_off = atom2basbas_off(i_atom_bb)
                 n_basbas_sp = sp2n_basbas_sp(species(i_atom_bb))

                 acc_errsq = 0.d0
                 acc_valsq = 0.d0
                 BASBAS: do i_basbas_sp = 1, n_basbas_sp
                    i_basbas = basbas_off + i_basbas_sp
                    FIND_BB: do i_loc_prodbas = 1, n_loc_prodbas
                       if (map_prodbas(i_loc_prodbas, myid+1) == i_basbas) then
                          acc_valsq = acc_valsq + &
                          &           coul3fn(i_basis_pair, i_loc_prodbas)**2
                          acc_errsq = acc_errsq + &
                          &    (coul3fn_lvl(i_basis_pair, i_loc_prodbas) - &
                          &     coul3fn(i_basis_pair, i_loc_prodbas))**2
                          exit FIND_BB
                       end if
                    end do FIND_BB
                 end do BASBAS
                 call sync_real_number(acc_valsq)
                 call sync_real_number(acc_errsq)
                 
                 if (myid == 0) then
                    write(22, "(3I3,2ES10.2,'   ',2A20)") &
                    & i_atom_1, i_atom_2, i_atom_bb, acc_valsq, acc_errsq, &
                    & trim(basis_name_1), trim(basis_name_2)
                 end if
                 max_errsq(i_atom_1, i_atom_2, i_atom_bb) = &
                 & max(max_errsq(i_atom_1, i_atom_2, i_atom_bb), acc_errsq)
                 max_valsq(i_atom_1, i_atom_2, i_atom_bb) = &
                 & max(max_valsq(i_atom_1, i_atom_2, i_atom_bb), acc_valsq)
              end do ATOM_BB
           end do BASIS_2
        end do BASIS_1
        write(22,*)
        write(22,*)
     end do ATOM_2
     write(22,*)
  end do ATOM_1
  write(22,*)

  if (myid == 0) then
     do i_atom_1 = 1, n_atoms
        do i_atom_2 = 1, i_atom_1
           do i_atom_bb = 1, n_atoms
              write(22,"(3I4,2ES10.2)") &
              & i_atom_1, i_atom_2, i_atom_bb, &
              & max_valsq(i_atom_1, i_atom_2, i_atom_bb), &
              & max_errsq(i_atom_1, i_atom_2, i_atom_bb)
           end do
        end do
     end do
     close(22)
  end if

contains

  ! --------------------------------------------------------------------------

  character*100 function basis_fnname(i_fn)
    implicit none
    integer, intent(IN) :: i_fn
    integer :: i_species
    character :: l_to_str ! function

    i_species = basisfn_species(i_fn)
    write(basis_fnname, "(A,' ',A,' ',I0,' ',A1)") &
    & trim(species_name(i_species)), trim(basisfn_type(i_fn)), &
    & basisfn_n(i_fn), l_to_str(basisfn_l(i_fn))
  end function basis_fnname

  ! --------------------------------------------------------------------------

  character*100 function basbas_fnname(i_fn)
    implicit none
    integer, intent(IN) :: i_fn
    integer :: i_species
    character :: l_to_str ! function

    i_species = basbasfn_species(i_fn)
    write(basbas_fnname, "(A,' ',A1)") &
    & trim(species_name(i_species)), l_to_str(basbasfn_l(i_fn))
  end function basbas_fnname

end program test_lvlpairs
