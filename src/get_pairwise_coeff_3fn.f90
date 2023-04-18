!****s* FHI-aims/get_pairwise_coeff_3fn
!  NAME
!    get_pairwise_coeff_3fn
!  SYNOPSIS

subroutine get_pairwise_coeff_3fn(i_species_1, i_species_2, n_Rvec, Rvec, &
&                                 coeff_3fn, d_coeff_3fn, calc_deriv)

  !  PURPOSE
  !
  !    Calculate the LVL expansion coefficients for two given species and
  !    different connection vectors
  !
  !    The modules localized_basbas and tight_binding_auxmat must already
  !    be initialized.
  !
  !  USES

  use runtime_choices, only : prodbas_threshold
  use localized_basbas
  use lvl_triples
  use tight_binding_auxmat
  use numerical_utilities
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: i_species_1, i_species_2
  integer, intent(IN) :: n_Rvec
  real*8, intent(IN) :: Rvec(3, n_Rvec)
  real*8, intent(INOUT) :: coeff_3fn(max_n_basis_sp, max_n_basis_sp, &
  &                                  max_n_basbas_sp, 2, n_Rvec)
  real*8, intent(INOUT) :: d_coeff_3fn(max_n_basis_sp, max_n_basis_sp, &
  &                                    max_n_basbas_sp, 2, 3, n_Rvec)
  logical, intent(IN) :: calc_deriv

  !  INPUTS
  !    o i_species_1, i_species_2 -- Atomic species of the two atoms
  !    o n_Rvec -- Number of distinct pairs of these sorts of atoms
  !    o Rvec -- Connecting vectors (coords2 - coords1)
  !    o calc_deriv -- flag if derivative should be calculated
  !  OUTPUTS
  !    o coeff_3fn -- Expansion coefficients in basis_sp^2 and basbas_sp order
  !    o d_coeff_3fn -- Derivative of coeff_3fn with respect to Rvec
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: n_basbas_sp_1, n_basbas_sp_2, n_basbas_sp_both, n1, n2, nb, i
  integer :: i_Rvec, n_side
  logical :: is_onsite
  real*8, allocatable :: this_o3fn(:,:,:), this_c3fn(:,:,:), this_x3fn(:,:,:)
  real*8, allocatable :: onsite_coulomb(:,:,:), offsite_coulomb(:,:), d_offsite_coulomb(:,:,:)
  real*8, allocatable :: inv_coulomb(:,:), d_coulomb(:,:)
  real*8, parameter :: nullvec(3) = (/0.d0, 0.d0, 0.d0/)
  real*8 :: dummy(1,1,1)
  integer :: info
  character(*), parameter :: func = 'get_pairwise_coeff_3fn'

  ! --- Allocate

  allocate(this_o3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_o3fn', func)
  allocate(this_c3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_c3fn', func)
  if(calc_deriv) then
     ! We need another auxiliary array
     allocate(this_x3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
     &        stat=info)
     call check_allocation(info, 'this_c3fn', func)
  endif

  allocate(onsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 2), stat=info)
  call check_allocation(info, 'onsite_coulomb', func)
  allocate(offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp), stat=info)
  call check_allocation(info, 'offsite_coulomb', func)
  allocate(inv_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
  call check_allocation(info, 'inv_coulomb', func)
  if(calc_deriv) then
     allocate(d_offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 3), stat=info)
     call check_allocation(info, 'd_offsite_coulomb', func)
     allocate(d_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
     call check_allocation(info, 'd_coulomb', func)
  else
     allocate(d_offsite_coulomb(1,1,1))
  endif

  ! --- Get expansion coefficients

  coeff_3fn = 0
  if(calc_deriv) d_coeff_3fn = 0

  call calculate_lvl_triples(i_species_1, i_species_2, n_Rvec, Rvec, coeff_3fn, d_coeff_3fn, calc_deriv)

  ! --- Get onsite Coulomb matrix elements

  n_basbas_sp_1 = sp2n_basbas_sp(i_species_1)
  n_basbas_sp_2 = sp2n_basbas_sp(i_species_2)
  n1 = n_basbas_sp_1
  n2 = n_basbas_sp_2

  ! onsite-1 (upper left)
  call fast_calculate_tb_auxmat(i_species_1, i_species_1, nullvec, &
  &                             onsite_coulomb(:,:, 1), dummy, .false.)
  ! onsite-2 (lower right)
  call fast_calculate_tb_auxmat(i_species_2, i_species_2, nullvec, &
  &                             onsite_coulomb(:,:, 2), dummy, .false.)

!DB exlcude ghost atoms without basis fctns
  if((n_basbas_sp_1.eq.0).and.(n_basbas_sp_2.eq.0)) return

  do i_Rvec = 1, n_Rvec

     this_o3fn = 0.d0

     ! 1st onsite term
     inv_coulomb(1:n1, 1:n1) = onsite_coulomb(1:n1, 1:n1, 1)
     this_o3fn(:,:, 1:n1) = coeff_3fn(:,:,1:n1, 1, i_Rvec)

     is_onsite = all(abs(Rvec(:, i_Rvec)) < 1d-10)
     if (is_onsite) then
        n_side = 1
        n_basbas_sp_both = n_basbas_sp_1
        nb = n_basbas_sp_both
     else
        n_side = 2
        n_basbas_sp_both = n_basbas_sp_1 + n_basbas_sp_2
        nb = n_basbas_sp_both

        ! 2nd onsite term
        inv_coulomb(n1+1:nb, n1+1:nb) = onsite_coulomb(1:n2, 1:n2, 2)
        this_o3fn(:,:, n1+1:nb) = coeff_3fn(:,:,1:n2, 2, i_Rvec)

        ! offsite terms
        call fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec(:, i_Rvec), &
                                      offsite_coulomb, d_offsite_coulomb, calc_deriv)
        inv_coulomb(1:n1, n1+1:nb) = offsite_coulomb(1:n1, 1:n2)
        inv_coulomb(n1+1:nb, 1:n1) = transpose(offsite_coulomb(1:n1, 1:n2))
     end if

     ! invert
     call power_genmat_lapack(n_basbas_sp_both, inv_coulomb(1:nb, 1:nb), &
     &                        -1.d0, safe_minimum, prodbas_threshold, '')

     ! apply
     call dgemm('N', 'N',max_n_basis_sp**2, nb, nb, &
                 1.d0, this_o3fn,max_n_basis_sp**2, &
                     inv_coulomb, 2*max_n_basbas_sp,&
                 0.d0, this_c3fn, max_n_basis_sp**2)

     ! put
     coeff_3fn(:,:, 1:n1, 1, i_Rvec) = this_c3fn(:,:,1:n1)
     if(is_onsite) then
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = 0.
        if(calc_deriv) d_coeff_3fn(:,:,:,:,:,i_Rvec) = 0.
     else
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = this_c3fn(:,:,n1+1:nb)
     endif

     if(calc_deriv .and. .not. is_onsite) then

        ! Care about derivative
        do i = 1, 3

           ! Part 1: inv_coulomb * (d coeff_3fn / d Rvec)
           this_o3fn(:,:,1:n1)    = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec)
           this_o3fn(:,:,n1+1:nb) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec)
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb, &
            &          1.d0, this_o3fn, max_n_basis_sp**2, &
            &                inv_coulomb, 2*max_n_basbas_sp, &
            &          0.d0, this_x3fn, max_n_basis_sp**2)

           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = this_x3fn(:,:,n1+1:n1+n2)

           ! Part 2: (d inv_coulomb / d Rvec) * coeff_3fn
           !       = - inv_coulomb * (d coulomb / d Rvec) * inv_coulomb * coeff_3fn
           ! Please note: this_c3fn still contains inv_coulomb * coeff_3fn
           ! The derivative of the onsite Coulomb matrices is 0

           d_coulomb(:,:) = 0
           d_coulomb(1:n1, n1+1:nb) = d_offsite_coulomb(1:n1, 1:n2, i)
           d_coulomb(n1+1:nb, 1:n1) = transpose(d_offsite_coulomb(1:n1, 1:n2, i))
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb, &
                      1.d0, this_c3fn, max_n_basis_sp**2,  &
                            d_coulomb, 2*max_n_basbas_sp,  &
                      0.d0, this_o3fn, max_n_basis_sp**2)
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb,  &
                       1.d0, this_o3fn, max_n_basis_sp**2,  &
                             inv_coulomb, 2*max_n_basbas_sp,&
                       0.d0, this_x3fn, max_n_basis_sp**2)

           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) - this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) - this_x3fn(:,:,n1+1:n1+n2)
        enddo
     endif

  end do

  deallocate(this_o3fn)
  deallocate(this_c3fn)

  deallocate(onsite_coulomb)
  deallocate(offsite_coulomb)
  deallocate(inv_coulomb)

  if(calc_deriv) then
     deallocate(this_x3fn)
     deallocate(d_coulomb)
  endif
  deallocate(d_offsite_coulomb)

end subroutine get_pairwise_coeff_3fn


!****s* FHI-aims/get_pairwise_coeff_3fn_vb
!  NAME
!    get_pairwise_coeff_3fn_vb
!  SYNOPSIS

subroutine get_pairwise_coeff_3fn_vb(i_atom_1, i_atom_2, i_species_1, i_species_2, atom2basis_len2, atom2vb_basis_off, n_atoms2, n_Rvec, Rvec, &
&                                 coeff_3fn, d_coeff_3fn, calc_deriv)

  !  PURPOSE
  !
  !    Calculate the LVL expansion coefficients for two given species and
  !    different connection vectors
  !
  !    The modules localized_basbas and tight_binding_auxmat must already
  !    be initialized.
  !
  !  USES

  use runtime_choices, only : prodbas_threshold
  use localized_basbas
  use lvl_triples
  use tight_binding_auxmat
  use numerical_utilities
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: i_atom_1, i_atom_2, i_species_1, i_species_2,n_atoms2
  integer, intent(IN) :: n_Rvec
  integer, intent(IN) :: atom2basis_len2(n_atoms2),atom2vb_basis_off(n_atoms2)
  real*8, intent(IN) :: Rvec(3, n_Rvec)
  real*8, intent(INOUT) :: coeff_3fn(max_n_basis_sp2, max_n_basis_sp2, &
  &                                  max_n_basbas_sp, 2, n_Rvec)
  real*8, intent(INOUT) :: d_coeff_3fn(max_n_basis_sp2, max_n_basis_sp2, &
  &                                    max_n_basbas_sp, 2, 3, n_Rvec)
  logical, intent(IN) :: calc_deriv

  !  INPUTS
  !    o i_species_1, i_species_2 -- Atomic species of the two atoms
  !    o n_Rvec -- Number of distinct pairs of these sorts of atoms
  !    o Rvec -- Connecting vectors (coords2 - coords1)
  !    o calc_deriv -- flag if derivative should be calculated
  !  OUTPUTS
  !    o coeff_3fn -- Expansion coefficients in basis_sp^2 and basbas_sp order
  !    o d_coeff_3fn -- Derivative of coeff_3fn with respect to Rvec
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: n_basbas_sp_1, n_basbas_sp_2, n_basbas_sp_both, n1, n2, nb, i
  integer :: i_Rvec, n_side
  logical :: is_onsite
  real*8, allocatable :: this_o3fn(:,:,:), this_c3fn(:,:,:), this_x3fn(:,:,:)
  real*8, allocatable :: onsite_coulomb(:,:,:), offsite_coulomb(:,:), d_offsite_coulomb(:,:,:)
  real*8, allocatable :: inv_coulomb(:,:), d_coulomb(:,:)
  real*8, parameter :: nullvec(3) = (/0.d0, 0.d0, 0.d0/)
  real*8 :: dummy(1,1,1)
  integer :: info
  character(*), parameter :: func = 'get_pairwise_coeff_3fn'

  ! --- Allocate

  allocate(this_o3fn(max_n_basis_sp2, max_n_basis_sp2, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_o3fn', func)
  allocate(this_c3fn(max_n_basis_sp2, max_n_basis_sp2, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_c3fn', func)
  if(calc_deriv) then
     ! We need another auxiliary array
     allocate(this_x3fn(max_n_basis_sp2, max_n_basis_sp2, 2*max_n_basbas_sp), &
     &        stat=info)
     call check_allocation(info, 'this_c3fn', func)
  endif

  allocate(onsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 2), stat=info)
  call check_allocation(info, 'onsite_coulomb', func)
  allocate(offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp), stat=info)
  call check_allocation(info, 'offsite_coulomb', func)
  allocate(inv_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
  call check_allocation(info, 'inv_coulomb', func)
  if(calc_deriv) then
     allocate(d_offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 3), stat=info)
     call check_allocation(info, 'd_offsite_coulomb', func)
     allocate(d_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
     call check_allocation(info, 'd_coulomb', func)
  else
     allocate(d_offsite_coulomb(1,1,1))
  endif

  ! --- Get expansion coefficients

  coeff_3fn = 0
  if(calc_deriv) d_coeff_3fn = 0

  call calculate_lvl_triples_vb(i_atom_1, i_atom_2, i_species_1, i_species_2, &
              atom2basis_len2, atom2vb_basis_off, n_atoms2, n_Rvec, Rvec, coeff_3fn, d_coeff_3fn, calc_deriv)

  ! --- Get onsite Coulomb matrix elements

  n_basbas_sp_1 = sp2n_basbas_sp(i_species_1)
  n_basbas_sp_2 = sp2n_basbas_sp(i_species_2)
  n1 = n_basbas_sp_1
  n2 = n_basbas_sp_2

  ! onsite-1 (upper left)
  call fast_calculate_tb_auxmat(i_species_1, i_species_1, nullvec, &
  &                             onsite_coulomb(:,:, 1), dummy, .false.)
  ! onsite-2 (lower right)
  call fast_calculate_tb_auxmat(i_species_2, i_species_2, nullvec, &
  &                             onsite_coulomb(:,:, 2), dummy, .false.)

!DB exlcude ghost atoms without basis fctns
  if((n_basbas_sp_1.eq.0).and.(n_basbas_sp_2.eq.0)) return

  do i_Rvec = 1, n_Rvec

     this_o3fn = 0.d0

     ! 1st onsite term
     inv_coulomb(1:n1, 1:n1) = onsite_coulomb(1:n1, 1:n1, 1)
     this_o3fn(:,:, 1:n1) = coeff_3fn(:,:,1:n1, 1, i_Rvec)

     is_onsite = all(abs(Rvec(:, i_Rvec)) < 1d-10)
     if (is_onsite) then
        n_side = 1
        n_basbas_sp_both = n_basbas_sp_1
        nb = n_basbas_sp_both
     else
        n_side = 2
        n_basbas_sp_both = n_basbas_sp_1 + n_basbas_sp_2
        nb = n_basbas_sp_both

        ! 2nd onsite term
        inv_coulomb(n1+1:nb, n1+1:nb) = onsite_coulomb(1:n2, 1:n2, 2)
        this_o3fn(:,:, n1+1:nb) = coeff_3fn(:,:,1:n2, 2, i_Rvec)

        ! offsite terms
        call fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec(:, i_Rvec), &
                                      offsite_coulomb, d_offsite_coulomb, calc_deriv)
        inv_coulomb(1:n1, n1+1:nb) = offsite_coulomb(1:n1, 1:n2)
        inv_coulomb(n1+1:nb, 1:n1) = transpose(offsite_coulomb(1:n1, 1:n2))
     end if

     ! invert
     call power_genmat_lapack(n_basbas_sp_both, inv_coulomb(1:nb, 1:nb), &
     &                        -1.d0, safe_minimum, prodbas_threshold, '')

     ! apply
     call dgemm('N', 'N',max_n_basis_sp2**2, nb, nb, &
                 1.d0, this_o3fn,max_n_basis_sp2**2, &
                     inv_coulomb, 2*max_n_basbas_sp,&
                 0.d0, this_c3fn, max_n_basis_sp2**2)

     ! put
     coeff_3fn(:,:, 1:n1, 1, i_Rvec) = this_c3fn(:,:,1:n1)
     if(is_onsite) then
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = 0.
        if(calc_deriv) d_coeff_3fn(:,:,:,:,:,i_Rvec) = 0.
     else
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = this_c3fn(:,:,n1+1:nb)
     endif

     if(calc_deriv .and. .not. is_onsite) then

        ! Care about derivative
        do i = 1, 3

           ! Part 1: inv_coulomb * (d coeff_3fn / d Rvec)
           this_o3fn(:,:,1:n1)    = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec)
           this_o3fn(:,:,n1+1:nb) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec)
           call dgemm('N', 'N', max_n_basis_sp2**2, nb, nb, &
            &          1.d0, this_o3fn, max_n_basis_sp2**2, &
            &                inv_coulomb, 2*max_n_basbas_sp, &
            &          0.d0, this_x3fn, max_n_basis_sp2**2)

           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = this_x3fn(:,:,n1+1:n1+n2)

           ! Part 2: (d inv_coulomb / d Rvec) * coeff_3fn
           !       = - inv_coulomb * (d coulomb / d Rvec) * inv_coulomb * coeff_3fn
           ! Please note: this_c3fn still contains inv_coulomb * coeff_3fn
           ! The derivative of the onsite Coulomb matrices is 0

           d_coulomb(:,:) = 0
           d_coulomb(1:n1, n1+1:nb) = d_offsite_coulomb(1:n1, 1:n2, i)
           d_coulomb(n1+1:nb, 1:n1) = transpose(d_offsite_coulomb(1:n1, 1:n2, i))
           call dgemm('N', 'N', max_n_basis_sp2**2, nb, nb, &
                      1.d0, this_c3fn, max_n_basis_sp2**2,  &
                            d_coulomb, 2*max_n_basbas_sp,  &
                      0.d0, this_o3fn, max_n_basis_sp2**2)
           call dgemm('N', 'N', max_n_basis_sp2**2, nb, nb,  &
                       1.d0, this_o3fn, max_n_basis_sp2**2,  &
                             inv_coulomb, 2*max_n_basbas_sp,&
                       0.d0, this_x3fn, max_n_basis_sp2**2)

           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) - this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) - this_x3fn(:,:,n1+1:n1+n2)
        enddo
     endif

  end do

  deallocate(this_o3fn)
  deallocate(this_c3fn)

  deallocate(onsite_coulomb)
  deallocate(offsite_coulomb)
  deallocate(inv_coulomb)

  if(calc_deriv) then
     deallocate(this_x3fn)
     deallocate(d_coulomb)
  endif
  deallocate(d_offsite_coulomb)

end subroutine get_pairwise_coeff_3fn_vb

subroutine my_get_pairwise_coeff_3fn(i_species_1, i_species_2, n_Rvec, Rvec, &
&                                 coeff_3fn, d_coeff_3fn, calc_deriv)

  !  PURPOSE
  !
  !    Calculate the LVL expansion coefficients for two given species and
  !    different connection vectors
  !
  !    The modules localized_basbas and tight_binding_auxmat must already
  !    be initialized.
  !
  !  USES

  use runtime_choices, only : prodbas_threshold
  use localized_basbas, only : SP2N_BASBAS_SP, SAFE_MINIMUM, MAX_N_BASIS_SP, MAX_N_BASBAS_SP
  use my_lvl_triples, only : my_calculate_lvl_triples
  use lvl_triples, only : calculate_lvl_triples
  use tight_binding_auxmat, only : fast_calculate_tb_auxmat, my_fast_calculate_tb_auxmat
  use mpi_tasks, only : myid, check_allocation
  use numerical_utilities
  use dimensions, only : use_threadsafe_gwinit
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: i_species_1, i_species_2
  integer, intent(IN) :: n_Rvec
  real*8, intent(IN) :: Rvec(3, n_Rvec)
  real*8, intent(INOUT) :: coeff_3fn(max_n_basis_sp, max_n_basis_sp, &
  &                                  max_n_basbas_sp, 2, n_Rvec)
  real*8, intent(INOUT) :: d_coeff_3fn(max_n_basis_sp, max_n_basis_sp, &
  &                                    max_n_basbas_sp, 2, 3, n_Rvec)
  logical, intent(IN) :: calc_deriv

  !  INPUTS
  !    o i_species_1, i_species_2 -- Atomic species of the two atoms
  !    o n_Rvec -- Number of distinct pairs of these sorts of atoms
  !    o Rvec -- Connecting vectors (coords2 - coords1)
  !    o calc_deriv -- flag if derivative should be calculated
  !  OUTPUTS
  !    o coeff_3fn -- Expansion coefficients in basis_sp^2 and basbas_sp order
  !    o d_coeff_3fn -- Derivative of coeff_3fn with respect to Rvec
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: n_basbas_sp_1, n_basbas_sp_2, n_basbas_sp_both, n1, n2, nb, i
  integer :: i_Rvec, n_side
  logical :: is_onsite
  real*8, allocatable :: this_o3fn(:,:,:), this_c3fn(:,:,:), this_x3fn(:,:,:)
  real*8, allocatable :: onsite_coulomb(:,:,:), offsite_coulomb(:,:), d_offsite_coulomb(:,:,:)
  real*8, allocatable :: inv_coulomb(:,:), d_coulomb(:,:)
  real*8, parameter :: nullvec(3) = (/0.d0, 0.d0, 0.d0/)
  real*8 :: dummy(1,1,1)
  integer :: info
  character(*), parameter :: func = 'get_pairwise_coeff_3fn'

  ! --- Allocate

  allocate(this_o3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_o3fn', func)
  allocate(this_c3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_c3fn', func)
  if(calc_deriv) then
     ! We need another auxiliary array
     allocate(this_x3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
     &        stat=info)
     call check_allocation(info, 'this_c3fn', func)
  endif

  allocate(onsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 2), stat=info)
  call check_allocation(info, 'onsite_coulomb', func)
  allocate(offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp), stat=info)
  call check_allocation(info, 'offsite_coulomb', func)
  allocate(inv_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
  call check_allocation(info, 'inv_coulomb', func)
  if(calc_deriv) then
     allocate(d_offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 3), stat=info)
     call check_allocation(info, 'd_offsite_coulomb', func)
     allocate(d_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
     call check_allocation(info, 'd_coulomb', func)
  else
     allocate(d_offsite_coulomb(1,1,1))
  endif

  ! --- Get expansion coefficients

  coeff_3fn = 0
  if(calc_deriv) d_coeff_3fn = 0

  if (use_threadsafe_gwinit) then
    call my_calculate_lvl_triples(i_species_1, i_species_2, n_Rvec, Rvec, coeff_3fn, d_coeff_3fn, calc_deriv)
  else
    call calculate_lvl_triples(i_species_1, i_species_2, n_Rvec, Rvec, coeff_3fn, d_coeff_3fn, calc_deriv)
  endif

  ! --- Get onsite Coulomb matrix elements

  n_basbas_sp_1 = sp2n_basbas_sp(i_species_1)
  n_basbas_sp_2 = sp2n_basbas_sp(i_species_2)
  n1 = n_basbas_sp_1
  n2 = n_basbas_sp_2

  if (use_threadsafe_gwinit) then
    ! onsite-1 (upper left)

    call my_fast_calculate_tb_auxmat(i_species_1, i_species_1, nullvec, 1, max_n_basbas_sp, 1, max_n_basbas_sp, &
                                     onsite_coulomb(:,:, 1), dummy, .false., .false.)
        ! onsite-2 (lower right)
    call my_fast_calculate_tb_auxmat(i_species_2, i_species_2, nullvec, 1, max_n_basbas_sp, 1, max_n_basbas_sp, &
                                     onsite_coulomb(:,:, 2), dummy, .false., .false.)
  else
    ! onsite-1 (upper left)

    call fast_calculate_tb_auxmat(i_species_1, i_species_1, nullvec, &
                                  onsite_coulomb(:,:, 1), dummy, .false.)

    ! onsite-2 (lower right)
     call fast_calculate_tb_auxmat(i_species_2, i_species_2, nullvec, &
                                   onsite_coulomb(:,:, 2), dummy, .false.)
  endif


!DB exlcude ghost atoms without basis fctns
  if((n_basbas_sp_1.eq.0).and.(n_basbas_sp_2.eq.0)) return

  do i_Rvec = 1, n_Rvec

     this_o3fn = 0.d0

     ! 1st onsite term
     inv_coulomb(1:n1, 1:n1) = onsite_coulomb(1:n1, 1:n1, 1)
     this_o3fn(:,:, 1:n1) = coeff_3fn(:,:,1:n1, 1, i_Rvec)

     is_onsite = all(abs(Rvec(:, i_Rvec)) < 1d-10)
     if (is_onsite) then
        n_side = 1
        n_basbas_sp_both = n_basbas_sp_1
        nb = n_basbas_sp_both
     else
        n_side = 2
        n_basbas_sp_both = n_basbas_sp_1 + n_basbas_sp_2
        nb = n_basbas_sp_both

        ! 2nd onsite term
        inv_coulomb(n1+1:nb, n1+1:nb) = onsite_coulomb(1:n2, 1:n2, 2)
        this_o3fn(:,:, n1+1:nb) = coeff_3fn(:,:,1:n2, 2, i_Rvec)

        ! offsite terms
        if (use_threadsafe_gwinit) then
          call my_fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec(:, i_Rvec), 1, max_n_basbas_sp, 1, max_n_basbas_sp, &
                                      offsite_coulomb, d_offsite_coulomb, calc_deriv, .false.)
        else
          call fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec(:, i_Rvec), &
                                        offsite_coulomb, d_offsite_coulomb, calc_deriv)
        endif

        inv_coulomb(1:n1, n1+1:nb) = offsite_coulomb(1:n1, 1:n2)
        inv_coulomb(n1+1:nb, 1:n1) = transpose(offsite_coulomb(1:n1, 1:n2))
     end if

     ! invert
     call perfon('gpc_inv')
     call power_genmat_lapack(n_basbas_sp_both, inv_coulomb(1:nb, 1:nb), &
     &                        -1.d0, safe_minimum, prodbas_threshold, '')
     call perfoff

     ! apply
     call dgemm('N', 'N',max_n_basis_sp**2, nb, nb, &
                 1.d0, this_o3fn,max_n_basis_sp**2, &
                     inv_coulomb, 2*max_n_basbas_sp,&
                 0.d0, this_c3fn, max_n_basis_sp**2)


     ! put
     coeff_3fn(:,:, 1:n1, 1, i_Rvec) = this_c3fn(:,:,1:n1)
     if(is_onsite) then
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = 0.
        if(calc_deriv) d_coeff_3fn(:,:,:,:,:,i_Rvec) = 0.
     else
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = this_c3fn(:,:,n1+1:nb)
     endif

     if(calc_deriv .and. .not. is_onsite) then

        ! Care about derivative
        do i = 1, 3

           ! Part 1: inv_coulomb * (d coeff_3fn / d Rvec)
           this_o3fn(:,:,1:n1)    = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec)
           this_o3fn(:,:,n1+1:nb) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec)
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb, &
            &          1.d0, this_o3fn, max_n_basis_sp**2, &
            &                inv_coulomb, 2*max_n_basbas_sp, &
            &          0.d0, this_x3fn, max_n_basis_sp**2)

           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = this_x3fn(:,:,n1+1:n1+n2)

           ! Part 2: (d inv_coulomb / d Rvec) * coeff_3fn
           !       = - inv_coulomb * (d coulomb / d Rvec) * inv_coulomb * coeff_3fn
           ! Please note: this_c3fn still contains inv_coulomb * coeff_3fn
           ! The derivative of the onsite Coulomb matrices is 0

           d_coulomb(:,:) = 0
           d_coulomb(1:n1, n1+1:nb) = d_offsite_coulomb(1:n1, 1:n2, i)
           d_coulomb(n1+1:nb, 1:n1) = transpose(d_offsite_coulomb(1:n1, 1:n2, i))
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb, &
                      1.d0, this_c3fn, max_n_basis_sp**2,  &
                            d_coulomb, 2*max_n_basbas_sp,  &
                      0.d0, this_o3fn, max_n_basis_sp**2)
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb,  &
                       1.d0, this_o3fn, max_n_basis_sp**2,  &
                             inv_coulomb, 2*max_n_basbas_sp,&
                       0.d0, this_x3fn, max_n_basis_sp**2)


           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) - this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) - this_x3fn(:,:,n1+1:n1+n2)
        enddo
     endif

  end do

  deallocate(this_o3fn)
  deallocate(this_c3fn)

  deallocate(onsite_coulomb)
  deallocate(offsite_coulomb)
  deallocate(inv_coulomb)

  if(calc_deriv) then
     deallocate(this_x3fn)
     deallocate(d_coulomb)
  endif
  deallocate(d_offsite_coulomb)

end subroutine my_get_pairwise_coeff_3fn

subroutine my_get_pairwise_coeff_3fn_old(i_species_1, i_species_2, n_Rvec, Rvec, &
&                                 coeff_3fn, d_coeff_3fn, calc_deriv)

  !  PURPOSE
  !
  !    Calculate the LVL expansion coefficients for two given species and
  !    different connection vectors
  !
  !    The modules localized_basbas and tight_binding_auxmat must already
  !    be initialized.
  !
  !  USES

  use runtime_choices, only : prodbas_threshold
  use localized_basbas
  use lvl_triples
  use tight_binding_auxmat
  use numerical_utilities
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: i_species_1, i_species_2
  integer, intent(IN) :: n_Rvec
  real*8, intent(IN) :: Rvec(3, n_Rvec)
  real*8, intent(INOUT) :: coeff_3fn(max_n_basis_sp, max_n_basis_sp, &
  &                                  max_n_basbas_sp, 2, n_Rvec)
  real*8, intent(INOUT) :: d_coeff_3fn(max_n_basis_sp, max_n_basis_sp, &
  &                                    max_n_basbas_sp, 2, 3, n_Rvec)
  logical, intent(IN) :: calc_deriv

  !  INPUTS
  !    o i_species_1, i_species_2 -- Atomic species of the two atoms
  !    o n_Rvec -- Number of distinct pairs of these sorts of atoms
  !    o Rvec -- Connecting vectors (coords2 - coords1)
  !    o calc_deriv -- flag if derivative should be calculated
  !  OUTPUTS
  !    o coeff_3fn -- Expansion coefficients in basis_sp^2 and basbas_sp order
  !    o d_coeff_3fn -- Derivative of coeff_3fn with respect to Rvec
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: n_basbas_sp_1, n_basbas_sp_2, n_basbas_sp_both, n1, n2, nb, i
  integer :: i_Rvec, n_side
  logical :: is_onsite
  real*8, allocatable :: this_o3fn(:,:,:), this_c3fn(:,:,:), this_x3fn(:,:,:)
  real*8, allocatable :: onsite_coulomb(:,:,:), offsite_coulomb(:,:), d_offsite_coulomb(:,:,:)
  real*8, allocatable :: inv_coulomb(:,:), d_coulomb(:,:)
  real*8, parameter :: nullvec(3) = (/0.d0, 0.d0, 0.d0/)
  real*8 :: dummy(1,1,1)
  integer :: info
  character(*), parameter :: func = 'get_pairwise_coeff_3fn'

  ! --- Allocate

  allocate(this_o3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_o3fn', func)
  allocate(this_c3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
  &        stat=info)
  call check_allocation(info, 'this_c3fn', func)
  if(calc_deriv) then
     ! We need another auxiliary array
     allocate(this_x3fn(max_n_basis_sp, max_n_basis_sp, 2*max_n_basbas_sp), &
     &        stat=info)
     call check_allocation(info, 'this_c3fn', func)
  endif

  allocate(onsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 2), stat=info)
  call check_allocation(info, 'onsite_coulomb', func)
  allocate(offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp), stat=info)
  call check_allocation(info, 'offsite_coulomb', func)
  allocate(inv_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
  call check_allocation(info, 'inv_coulomb', func)
  if(calc_deriv) then
     allocate(d_offsite_coulomb(max_n_basbas_sp, max_n_basbas_sp, 3), stat=info)
     call check_allocation(info, 'd_offsite_coulomb', func)
     allocate(d_coulomb(2*max_n_basbas_sp, 2*max_n_basbas_sp), stat=info)
     call check_allocation(info, 'd_coulomb', func)
  else
     allocate(d_offsite_coulomb(1,1,1))
  endif

  ! --- Get expansion coefficients

  coeff_3fn = 0
  if(calc_deriv) d_coeff_3fn = 0

  call calculate_lvl_triples(i_species_1, i_species_2, n_Rvec, Rvec, coeff_3fn, d_coeff_3fn, calc_deriv)

  ! --- Get onsite Coulomb matrix elements

  n_basbas_sp_1 = sp2n_basbas_sp(i_species_1)
  n_basbas_sp_2 = sp2n_basbas_sp(i_species_2)
  n1 = n_basbas_sp_1
  n2 = n_basbas_sp_2

  ! onsite-1 (upper left)
  call fast_calculate_tb_auxmat(i_species_1, i_species_1, nullvec, &
  &                             onsite_coulomb(:,:, 1), dummy, .false.)
  ! onsite-2 (lower right)
  call fast_calculate_tb_auxmat(i_species_2, i_species_2, nullvec, &
  &                             onsite_coulomb(:,:, 2), dummy, .false.)

!DB exlcude ghost atoms without basis fctns
  if((n_basbas_sp_1.eq.0).and.(n_basbas_sp_2.eq.0)) return

  do i_Rvec = 1, n_Rvec

     this_o3fn = 0.d0

     ! 1st onsite term
     inv_coulomb(1:n1, 1:n1) = onsite_coulomb(1:n1, 1:n1, 1)
     this_o3fn(:,:, 1:n1) = coeff_3fn(:,:,1:n1, 1, i_Rvec)

     is_onsite = all(abs(Rvec(:, i_Rvec)) < 1d-10)
     if (is_onsite) then
        n_side = 1
        n_basbas_sp_both = n_basbas_sp_1
        nb = n_basbas_sp_both
     else
        n_side = 2
        n_basbas_sp_both = n_basbas_sp_1 + n_basbas_sp_2
        nb = n_basbas_sp_both

        ! 2nd onsite term
        inv_coulomb(n1+1:nb, n1+1:nb) = onsite_coulomb(1:n2, 1:n2, 2)
        this_o3fn(:,:, n1+1:nb) = coeff_3fn(:,:,1:n2, 2, i_Rvec)

        ! offsite terms
        call fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec(:, i_Rvec), &
                                      offsite_coulomb, d_offsite_coulomb, calc_deriv)
        inv_coulomb(1:n1, n1+1:nb) = offsite_coulomb(1:n1, 1:n2)
        inv_coulomb(n1+1:nb, 1:n1) = transpose(offsite_coulomb(1:n1, 1:n2))
     end if

     ! invert
     call power_genmat_lapack(n_basbas_sp_both, inv_coulomb(1:nb, 1:nb), &
     &                        -1.d0, safe_minimum, prodbas_threshold, '')

     ! apply
     call dgemm('N', 'N',max_n_basis_sp**2, nb, nb, &
                 1.d0, this_o3fn,max_n_basis_sp**2, &
                     inv_coulomb, 2*max_n_basbas_sp,&
                 0.d0, this_c3fn, max_n_basis_sp**2)

     ! put
     coeff_3fn(:,:, 1:n1, 1, i_Rvec) = this_c3fn(:,:,1:n1)
     if(is_onsite) then
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = 0.
        if(calc_deriv) d_coeff_3fn(:,:,:,:,:,i_Rvec) = 0.
     else
        coeff_3fn(:,:, 1:n2, 2, i_Rvec) = this_c3fn(:,:,n1+1:nb)
     endif

     if(calc_deriv .and. .not. is_onsite) then

        ! Care about derivative
        do i = 1, 3

           ! Part 1: inv_coulomb * (d coeff_3fn / d Rvec)
           this_o3fn(:,:,1:n1)    = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec)
           this_o3fn(:,:,n1+1:nb) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec)
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb, &
            &          1.d0, this_o3fn, max_n_basis_sp**2, &
            &                inv_coulomb, 2*max_n_basbas_sp, &
            &          0.d0, this_x3fn, max_n_basis_sp**2)

           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = this_x3fn(:,:,n1+1:n1+n2)

           ! Part 2: (d inv_coulomb / d Rvec) * coeff_3fn
           !       = - inv_coulomb * (d coulomb / d Rvec) * inv_coulomb * coeff_3fn
           ! Please note: this_c3fn still contains inv_coulomb * coeff_3fn
           ! The derivative of the onsite Coulomb matrices is 0

           d_coulomb(:,:) = 0
           d_coulomb(1:n1, n1+1:nb) = d_offsite_coulomb(1:n1, 1:n2, i)
           d_coulomb(n1+1:nb, 1:n1) = transpose(d_offsite_coulomb(1:n1, 1:n2, i))
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb, &
                      1.d0, this_c3fn, max_n_basis_sp**2,  &
                            d_coulomb, 2*max_n_basbas_sp,  &
                      0.d0, this_o3fn, max_n_basis_sp**2)
           call dgemm('N', 'N', max_n_basis_sp**2, nb, nb,  &
                       1.d0, this_o3fn, max_n_basis_sp**2,  &
                             inv_coulomb, 2*max_n_basbas_sp,&
                       0.d0, this_x3fn, max_n_basis_sp**2)

           d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) = d_coeff_3fn(:,:,1:n1,1,i,i_Rvec) - this_x3fn(:,:,1:n1)
           d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) = d_coeff_3fn(:,:,1:n2,2,i,i_Rvec) - this_x3fn(:,:,n1+1:n1+n2)
        enddo
     endif

  end do

  deallocate(this_o3fn)
  deallocate(this_c3fn)

  deallocate(onsite_coulomb)
  deallocate(offsite_coulomb)
  deallocate(inv_coulomb)

  if(calc_deriv) then
     deallocate(this_x3fn)
     deallocate(d_coulomb)
  endif
  deallocate(d_offsite_coulomb)

end subroutine my_get_pairwise_coeff_3fn_old
!******
