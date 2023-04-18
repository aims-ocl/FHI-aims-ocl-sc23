!****s* FHI-aims/density_matrix_mixing
!  NAME
!   density_matrix_mixing
!  SYNOPSIS

subroutine density_matrix_mixing_dmft(i_spin, number_of_loops, density_matrix)

  !  PURPOSE
  !
  !    Mixing density matrix for self-consistent Hartree-Fock or hybrid
  !    functional calculations.
  !
  !    Please note that this procedure is also used for Fock matrix mixing for
  !    hf_version == HF_EIGEN because the structures of fock_matr and agree.
  !    In principle, this could always be done (would be cleaner), but the
  !    Fock matrix is much more expensive to set up and is not needed
  !    otherwise in the initialization step.
  !
  !  USES

  use mpi_tasks
  use dimensions
  use runtime_choices
  use hartree_fock
  use mixing, only: pulay_saved_iter_denmat, mixing_factor

  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: i_spin
  integer, intent(IN) :: number_of_loops
  real*8, intent(INOUT) :: density_matrix(n_basis, n_basis)

  !   INPUTS 
  !     o i_spin -- the current spin      
  !     o number_of_loops -- the current self-consistent iteration
  !   INPUTS/OUPUTS
  !     o density_matrix -- upon input, the density_matrix of current iteration,
  !                         upon output, the mixed density_matrix to be used
  !                                      for evaluting the Fock matrix
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

  real*8, dimension(:), allocatable :: delta_denmat ! output - input  (error)
  real*8, dimension(:), allocatable :: denmat_diff  ! mixed - output (update)
  real*8, dimension(:,:), allocatable:: prev_denmat_temp  ! mixed - output (update)
  integer i_store
  character*150 :: info_str
  integer :: info
  character(*), parameter :: func = 'density_matrix_mixing'


  allocate(delta_denmat(n_local_denmat), stat=info)
  call check_allocation(info, 'delta_denmat', func)
  allocate(denmat_diff(n_local_denmat), stat=info)
  call check_allocation(info, 'denmat_diff', func)
if(.not. allocated(prev_denmat_temp))then
allocate( prev_denmat_temp(n_basis*(n_basis+1)/2,n_spin))
endif
write(use_unit,*) '1 line!!!!!!!!!!!!!!!!!', prev_denmat_temp, n_local_denmat

  if (number_of_loops <= 1) then
     ! On first entry, set D^0_in to D^0_out, because this is still better
     ! than using zero.
     prev_denmat_temp(:, i_spin) = 0.d0
write(use_unit,*) '1 line!!!!!!!!!!!!!!!!!', prev_denmat_temp
     call glob2loc(prev_denmat_temp(:, i_spin), density_matrix, 1.d0)
  end if
write(use_unit,*) '2 line!!!!!!!!!!!!!!!!!'

  ! delta_denmat := density_matrix - prev_denmat
  delta_denmat = - prev_denmat_temp(:, i_spin)
  call glob2loc(delta_denmat, density_matrix, 1.d0)

  !    linear mixing
  if(mixer.eq.MIX_LINEAR) then

     denmat_diff = linear_mix_param(i_spin) * delta_denmat

  else if (mixer.eq.MIX_PULAY .and. number_of_loops.le.1) then

     ! I am not completely sure if this should better read linear_mix_param to
     ! match the corresponding charge density Pulay mixer.  But first, that
     ! variable is not guaranteed to be initialized, and second, delta_denmat
     ! is zero anyway.
     denmat_diff = charge_mix_param(i_spin) * delta_denmat

  else if (mixer.eq.MIX_PULAY) then

     !   Update the density matrix difference

     denmat_diff(:) = charge_mix_param(i_spin) * delta_denmat(:)

     if(pulay_saved_iter_denmat .gt. 0) then

        denmat_diff (:) = denmat_diff(:) + &
        mixing_factor(1) * &
        ( prev_denmat_diff(:,1,i_spin) + &
        charge_mix_param(i_spin) * ( &
        delta_denmat(:) - prev_denmat_error(:,1,i_spin) ) &
        )
     endif

     do i_store = 2, pulay_saved_iter_denmat, 1
        denmat_diff(:) = denmat_diff(:) + &
        mixing_factor(i_store) * &
        ( prev_denmat_diff(:,i_store, i_spin) + &
        charge_mix_param(i_spin) * ( &
        prev_denmat_error(:,i_store-1,i_spin) - &
        prev_denmat_error(:,i_store,i_spin) ) &
        )
     enddo

  end if
write(use_unit,*) '3 line!!!!!!!!!!!!!!!!!'

  ! Store result for next cycle
  if (mixer .eq. MIX_PULAY) then

     ! Shift old Pulay stores
     do i_store = pulay_saved_iter_denmat, 2, -1

        prev_denmat_error(:,i_store,i_spin) = &
        prev_denmat_error(:,i_store-1,i_spin)

        prev_denmat_diff(:,i_store,i_spin) = &
        prev_denmat_diff(:,i_store-1,i_spin)

     enddo

     ! Save this Pulay iteration
     prev_denmat_error(:,1,i_spin) = delta_denmat
     if (number_of_loops <= 1) then
        prev_denmat_diff(:,1,i_spin) = delta_denmat * charge_mix_param(i_spin)
     else
        prev_denmat_diff(:,1,i_spin) = denmat_diff(:)
     end if

  endif

  ! Update density_matrix
  density_matrix = 0.d0
  prev_denmat_temp(:, i_spin) = prev_denmat_temp(:, i_spin) + denmat_diff
  call loc2glob(density_matrix, prev_denmat_temp(:, i_spin), 1.d0)
  call sync_matrix(density_matrix, n_basis, n_basis)

  deallocate(delta_denmat, denmat_diff)

contains

  !----------------------------------------------------------------------------
  !****s* FHI-aims/density_matrix_mixing/glob2loc
  !  NAME
  !    glob2loc
  !  SYNOPSIS

  subroutine glob2loc(loc, glob, fac)

    !  PURPOSE
    !
    !    Add fac*glob to loc.  Assume glob to be symmetric.
    !
    !  USES

    implicit none

    !  ARGUMENTS
!    real*8, intent(INOUT) :: loc(n_local_denmat)
    real*8, intent(INOUT) :: loc(n_local_denmat)
    real*8, intent(IN) :: glob(n_basis, n_basis)
    real*8, intent(IN) :: fac

    !  INPUTS
    !    o loc -- local matrix on input
    !    o glob -- global matrix
    !    o fac -- factor for local matrix
    !  OUTPUTS
    !    o loc -- loc + fac*glob
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_glob, i_loc, i_basis_1, i_basis_2
    character(*), parameter :: func = 'density_matrix_mixing:glob2loc'

    i_glob = 0
    i_loc = 0
    do i_basis_1 = 1, n_basis
       do i_basis_2 = 1, i_basis_1
          i_glob = i_glob + 1
          if (modulo(i_glob, n_tasks) == myid) then
             i_loc = i_loc + 1
             if (i_loc > n_local_denmat) call aims_stop('Invalid i_loc', func)
             loc(i_loc) = loc(i_loc) + fac * glob(i_basis_2, i_basis_1)
          end if
       end do
    end do

  end subroutine glob2loc
  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/density_matrix_mixing/loc2glob
  !  NAME
  !    loc2glob
  !  SYNOPSIS

  subroutine loc2glob(glob, loc, fac)

    !  PURPOSE
    !
    !    Add fac*loc to glob.  Assume glob to be symmetric.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT) :: glob(n_basis, n_basis)
    real*8, intent(IN) :: loc(n_local_denmat)
    real*8, intent(IN) :: fac

    !  INPUTS
    !    o glob -- global matrix on input
    !    o loc -- local matrix
    !    o fac -- factor for global matrix
    !  OUTPUTS
    !    o glob -- glob + fac*loc
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_glob, i_loc, i_basis_1, i_basis_2
    character(*), parameter :: func = 'density_matrix_mixing:loc2glob'

    i_glob = 0
    i_loc = 0
    do i_basis_1 = 1, n_basis
       do i_basis_2 = 1, i_basis_1
          i_glob = i_glob + 1
          if (modulo(i_glob, n_tasks) == myid) then
             i_loc = i_loc + 1
             glob(i_basis_2, i_basis_1) = glob(i_basis_2, i_basis_1) &
             &                            + fac * loc(i_loc)
             if (i_basis_1 /= i_basis_2) then
                glob(i_basis_1, i_basis_2) = glob(i_basis_1, i_basis_2) &
                &                            + fac * loc(i_loc)
             end if
          end if
       end do
    end do

  end subroutine loc2glob
  !******

end subroutine density_matrix_mixing_dmft
!******
