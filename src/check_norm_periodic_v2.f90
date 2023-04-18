!******
!-------------------------------------------------------------------------------
!****s* force_occupation/check_norm_periodic_v2
!  NAME
!    check_norm_periodic_v2
!  SYNOPSIS

subroutine check_norm_periodic_v2(chemical_potential, KS_eigenvalue, n_electrons,&
            occ_numbers, diff_electrons, i_counter,current_spin)

!  PURPOSE
!  checks and returns the charge density norm (i.e. the electron count) for
!  a given Fermi level, and for whichever smearing type was chosen ...!
!
!  USES
! 
  use dimensions,       only: spin_degeneracy, force_occupation_basis, &
                              n_periodic, n_states, n_k_points, n_force_occ  
  use runtime_choices,  only: occupation_width, occupation_type, &
                              n_methfessel_paxton
  use arch_specific
  use constants
  use pbc_lists
  use synchronize_mpi
  use physics,          only: KS_eigenvector_complex
  use force_occupation, only: force_occ_state_periodic, &
                              force_occ_max_KS_state, which_basis_function, &
                              force_occ_spin, force_occ_basis_occupation
  use localorB_io,      only: use_unit
  implicit none

!  ARGUMENTS

  real*8, intent(inout) :: chemical_potential
  real*8, dimension(n_states,n_k_points), intent(in) :: KS_eigenvalue
  real*8, intent(in) :: n_electrons
  real*8, dimension(n_states,n_k_points), intent(out) :: occ_numbers
  real*8, intent(out) :: diff_electrons
  integer, intent(inout) :: i_counter 

  integer :: current_spin


!  INPUT
!   o chemical_potential -- chemical potential
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o n_electrons -- number of electrons
!    
!  OUTPUT
!   o occ_numbers -- occupations of the states
!   o diff_electrons -- Electron count with current (not yet converged) Fermi level vs.
!                    the actual electron count if we had the correct Fermi level already
!   o i_counter -- number of iterations that have already passed in the "zeroin" scheme used
!                  to find the Fermi level iteratively 
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
!

  
  !  local variables
  real*8 :: temp_n_electrons
  real*8 :: one_over_ow
  real*8 :: H_odd
  real*8 :: H_even
  real*8 :: hermite_arg
  real*8 :: gauss_weight
  real*8 :: A
  real*8 :: max_exponent
  real*8 :: exp_arg
  integer:: my_basis_fn
     
  !  counter
  integer :: i_state,i_force_occ,i_k_point
  integer :: i_spin,i_basis,i_k,j_state
  integer :: i_mp

  real*8 :: temp,temp_2
  integer :: change(n_states,n_force_occ,n_k_points)
  integer :: change_i
  integer :: change_j
  integer :: force_spin

!  write(use_unit,*) "check_norm... n_spin=", n_spin, "n_states=", n_states
!  do i_spin = 1, n_spin, 1
!     do i_state = 1, n_states, 1
!        write(use_unit,*) i_state, i_spin, KS_eigenvalue(i_state, i_spin)
!     end do
!  end do
  i_counter = i_counter + 1
  one_over_ow   = 1. / occupation_width
  temp_n_electrons   = 0.d0
  select case(occupation_type)
  case (0)
    !  gaussian smearing (0 <= occ_numbers <= spin_degenracy)
    do i_k_point = 1, n_k_points, 1
      do i_state = 1, n_states, 1
        occ_numbers(i_state,i_k_point) = spin_degeneracy * 0.5d0 * & 
             (1 - arch_erf((KS_eigenvalue(i_state,i_k_point) - chemical_potential) * one_over_ow))
        temp_n_electrons = temp_n_electrons + occ_numbers(i_state,i_k_point)*k_weights(i_k_point)
      end do
    end do

  case (1)
    !  fermi smearing (0 <= occ_numbers <= spin_degeneracy)
    max_exponent = maxexponent(chemical_potential) * log(2.d0)
    do i_k_point = 1, n_k_points, 1
      do i_state = 1, n_states, 1
        exp_arg = (KS_eigenvalue(i_state,i_k_point) - chemical_potential) * one_over_ow
        if (exp_arg .lt. max_exponent) then
          occ_numbers(i_state,i_k_point) = spin_degeneracy / (1 + exp(exp_arg))
          temp_n_electrons = temp_n_electrons + occ_numbers(i_state,i_k_point)*k_weights(i_k_point)
        else
          occ_numbers(i_state,i_k_point) = 0.d0
        end if
      end do
    end do
    
  case (2)
    ! Methfessel-Paxton smearing (0 <= occ_numbers <= spin_degeneracy)
    do i_k_point = 1, n_k_points, 1
      do i_state = 1, n_states, 1
        hermite_arg  = (KS_eigenvalue(i_state,i_k_point) - chemical_potential) * one_over_ow 
        gauss_weight = exp(- hermite_arg * hermite_arg)
        !     zero order contribution
        occ_numbers(i_state,i_k_point) = 0.5 * (1 - arch_erf(hermite_arg))
        if (n_methfessel_paxton .gt. 0) then
          !    write(use_unit,*) A, H_even, H_odd, hermite_arg
          !     first order contribution
          A = - 0.25 * pisqrt_inv
          !     H_even = H_0 = 1
          H_even = 1.d0
          !     H_odd =  H_1 = 2 * x
          H_odd  = 2 * hermite_arg
          occ_numbers(i_state,i_k_point) = occ_numbers(i_state,i_k_point) + A * H_odd * gauss_weight
        end if
        if (n_methfessel_paxton .gt. 1) then
          do i_mp = 2, n_methfessel_paxton, 1
            A = - 1. / dble(4 * i_mp) * A
            H_even = 2 * hermite_arg * H_odd  - 2 *  i_mp      * H_even
            H_odd  = 2 * hermite_arg * H_even - 2 * (i_mp + 1) * H_odd
            occ_numbers(i_state,i_k_point) = occ_numbers(i_state,i_k_point) + A * H_odd * gauss_weight
          end do
        end if
        !write(use_unit,*) occ_numbers(i_state, i_spin)
        occ_numbers(i_state,i_k_point) = occ_numbers(i_state,i_k_point) * spin_degeneracy
        !write(use_unit,*) occ_numbers(i_state, i_spin)
        temp_n_electrons = temp_n_electrons + occ_numbers(i_state,i_k_point)*k_weights(i_k_point)
      end do
      !write(use_unit,*) "temp_n_electrons= ", temp_n_electrons, " chemical_potential= ", chemical_potential
    end do

  case (3)
    ! integer occupation (occ_numbers = 0, or spin_degenracy) test by igor
    if (n_k_points.gt.1) then
      write(use_unit,*) "* At present, the integer occupation scheme can not be used to the periodic model."
      write(use_unit,*) "* Abort."
      stop
    end if
    chemical_potential = 0.5d0*(KS_eigenvalue(1,1)+KS_eigenvalue(2,1))
    do i_k_point = 1, n_k_points, 1
      temp_n_electrons = 0.d0
      i_state = 1
      diff_electrons = n_electrons
      do while (diff_electrons .gt. 1.0d-8)
        occ_numbers(i_state,i_k_point) = spin_degeneracy
        temp_n_electrons = temp_n_electrons + occ_numbers(i_state,i_k_point)
        diff_electrons = n_electrons - temp_n_electrons
        i_state = i_state + 1
      end do
      if (KS_eigenvalue(i_state-1,i_k_point).gt.chemical_potential) then
        chemical_potential = 0.5d0*(KS_eigenvalue(i_state-1,i_k_point)+KS_eigenvalue(i_state,i_k_point))
      end if
    end do

  case default
    write(use_unit,*) "* Unknown smearing type in subroutine check_norm."
    write(use_unit,*) "* Abort."
    stop
     
  end select

  ! Occupation of basis function is requested
  if (force_occupation_basis) then

    force_occ_state_periodic = 0

    do i_k_point = 1, n_k_points, 1
      do i_force_occ = 1, n_force_occ, 1
        do i_state = 1, force_occ_max_KS_state(i_force_occ), 1
          change(i_state,i_force_occ,i_k_point) = i_state
        end do
      end do
    end do


    i_k = 0
    do i_k_point = 1, n_k_points, 1
      
      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
      
        i_k = i_k + 1
      
        do i_force_occ = 1, n_force_occ, 1
          my_basis_fn = which_basis_function(i_force_occ)
          force_spin = force_occ_spin(i_force_occ)
          do j_state = 1, force_occ_max_KS_state(i_force_occ), 1
            change_j = change(j_state,i_force_occ,i_k_point)
            do i_state = 1, force_occ_max_KS_state(i_force_occ), 1
              change_i = change(i_state,i_force_occ,i_k_point)
              if (abs(KS_eigenvector_complex(my_basis_fn, change_j, force_spin, i_k))&
                .gt.abs(KS_eigenvector_complex(my_basis_fn, change_i, force_spin, i_k)))&
              then
                change(j_state,i_force_occ,i_k_point) = change_i 
                change(i_state,i_force_occ,i_k_point) = change_j
              end if
            end do
          end do
        end do

        do i_force_occ = 1, n_force_occ, 1
           force_occ_state_periodic(i_force_occ,i_k_point) = change(1,i_force_occ,i_k_point)
        end do
      
      end if ! i_k ...
      
    end do ! k_points
    ! write(use_unit,*)">>>>>>>>>>>>>>>>>>>"
    ! do i_state = 1, force_occ_max_KS_state(1), 1
    ! write(use_unit,*) change(i_state,1),  KS_eigenvector_complex(10,change(i_state,1),1,1)
    ! enddo
    ! write(use_unit,*)"<<<<<<<<<<<<<<<<<<<"
    call sync_integer_vector(force_occ_state_periodic,n_force_occ*n_k_points)
      
      
      
    if (n_periodic .eq. 0) then
      write(use_unit,*) "* Using periodic version of check_norm for non-periodic case."
      write(use_unit,*) "* Abort."
      stop
    else

      do i_force_occ = 1, n_force_occ, 1
        if (current_spin.eq.force_occ_spin(i_force_occ)) then
          do i_k_point = 1, n_k_points, 1                                     

            temp_n_electrons = temp_n_electrons &
              - occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point),i_k_point)*k_weights(i_k_point)
            occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point),i_k_point) = &
              force_occ_basis_occupation(i_force_occ)
            temp_n_electrons = temp_n_electrons &
              + occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point), i_k_point)*k_weights(i_k_point)

          end do
        end if
      end do
    end if 

  end if


  diff_electrons = temp_n_electrons - n_electrons
  
end subroutine check_norm_periodic_v2






