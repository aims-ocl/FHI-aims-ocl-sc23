!****h* FHI-aims/force_occupation
!  NAME
!    force_occupation
!  SYNOPSIS

module force_occupation

!  PURPOSE

! Module force_occupations contains all infrastructure and subroutines needed
! for forcing the occupation numbers of any KS state using a restart file based projection
! 
! M.Gramzow
!
!  Subroutines and functions:
! 
!  list of all to come, e.g. check_norm*, zeroin*,...
!  allocation and deallocation routines

! allocate_force_occupation
! allocate_previous_eigenvector
! deallocate_force_occupation
! deallocate_previous_eigenvector
! check_norm_force_occ_p0
! check_norm_force_occ_v2
! check_norm_constraint_occ
! determine_corehole_projection
! zeroin_occ_p0
! zeroin_occ_v2
! zeroin_constraint_occ
! get_constraint_fermi_occ
! get_occupation_numbers_occ_p0
! get_occupation_numbers_occ_v2

! evaluate_gradient_constraint_occ

! global variable declarations


! force_occ_pr_state    :   number of KS state provided in restart file, read in control.in
! force_occ_spin        :   which spin channel is force_occ_pr_state belonging to
! forced_occ_number     :   desired occupation number for force_occ_pr_state
! force_occ_state       :   number of KS state whose occupation number has been changed during previous scf-cycle
! previous_eigenvector  :   KS eigenvector from previous iteration needed for projection

  use dimensions
!   use geometry
  use basis
  use runtime_choices
  use constraint
  use localorb_io
  use bfgs

  implicit none

! these are used for cluster calculations

  integer,dimension(:),allocatable :: force_occ_pr_state
  integer,dimension(:),allocatable :: force_occ_pr_state2
  integer,dimension(:),allocatable :: force_occ_spin
  real*8,dimension(:),allocatable :: forced_occ_number
  integer,dimension(:),allocatable :: force_occ_state
  real*8,dimension(:,:,:,:),allocatable :: previous_eigenvector

  integer, dimension(:),allocatable :: force_occ_max_KS_state
  integer, dimension(:),allocatable :: force_occ_min_KS_state
  integer, dimension(:),allocatable :: force_occ_autorep
  integer, dimension(:),allocatable :: force_occ_autored
  integer, dimension(:),allocatable :: fop_detect_auto

  integer, dimension(:),allocatable :: force_occ_maxred_KS_state
  integer, dimension(:),allocatable :: force_occ_step_KS_state

! extra dimension for k_points needed or not? I just do not understand right now

  integer,dimension(:,:),allocatable :: force_occ_state_periodic
  integer,dimension(:,:),allocatable :: force_occ_pr_state_periodic
  complex*16,dimension(:,:,:,:),allocatable :: previous_eigenvector_complex

! these are needed for basis projection - no corehole_projection needed

      character*20,dimension(:),allocatable :: force_occ_basis_type
      integer,dimension(:),allocatable :: force_occ_atom
      integer,dimension(:),allocatable :: force_occ_basis_n
      integer,dimension(:),allocatable :: force_occ_basis_l
      integer,dimension(:),allocatable :: force_occ_basis_m
      integer,dimension(:),allocatable :: which_basis_function
      real*8,dimension(:),allocatable :: force_occ_basis_occupation

  real*8   :: force_occupation_smearing_width = 0.01

  contains


!******
!-------------------------------------------------------------------------------
!****s* force_occupation/allocate_force_occ_state
!  NAME
!   allocate_force_occ_state
!  SYNOPSIS

  subroutine allocate_force_occ_state()

  use mpi_tasks, only: aims_stop
  implicit none

  integer :: i_force_occ,i_basis, i_fn

    if (.not. allocated(force_occ_state_periodic)) then

      allocate(force_occ_state_periodic(n_force_occ,n_k_points))

    end if

    which_basis_function = 0

    do  i_force_occ = 1, n_force_occ, 1
       look_basis: do i_basis = 1, n_basis, 1
          i_fn = basis_fn(i_basis)


          if ((basisfn_n(i_fn).eq.force_occ_basis_n(i_force_occ)) .and. &
               (basis_l(i_basis).eq.force_occ_basis_l(i_force_occ))  .and. &
               (basis_m(i_basis).eq.force_occ_basis_m(i_force_occ))  .and. &
               (basis_atom(i_basis).eq.force_occ_atom(i_force_occ))  .and. &
               (basisfn_type(i_fn).eq.force_occ_basis_type(i_force_occ)) ) then
             which_basis_function(i_force_occ) = i_basis
             exit look_basis

     
          end if
       end do look_basis

      if(which_basis_function(i_force_occ) .eq. 0)then
         write(use_unit,*)'Error: force occ basis ', i_force_occ,' not found'
         call aims_stop
      end if

   end do



  end subroutine allocate_force_occ_state
!******
!-------------------------------------------------------------------------------
!****s* force_occupation/allocate_force_occupation_basis
!  NAME
!   allocate_force_occupation_basis
!  SYNOPSIS

  subroutine allocate_force_occupation_basis()


  if (n_periodic .eq. 0) then
    if (.not. allocated(force_occ_state)) then
      allocate(force_occ_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_type)) then
      allocate(force_occ_basis_type(n_force_occ))
    end if
    if (.not. allocated(force_occ_atom)) then
      allocate(force_occ_atom(n_force_occ))
    end if
    if (.not. allocated(force_occ_spin)) then
      allocate(force_occ_spin(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_n)) then
      allocate(force_occ_basis_n(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_l)) then
      allocate(force_occ_basis_l(n_force_occ))
    end if
    if (.not. allocated(force_occ_max_KS_state)) then
      allocate(force_occ_max_KS_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_m)) then
      allocate(force_occ_basis_m(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_occupation)) then
      allocate(force_occ_basis_occupation(n_force_occ))
    end if
    if (.not. allocated(which_basis_function)) then
      allocate(which_basis_function(n_force_occ))
    end if
  else
!     if (.not. allocated(force_occ_state_periodic)) then
!       allocate(force_occ_state_periodic(n_force_occ,n_k_points))
!     end if
    if (.not. allocated(force_occ_basis_type)) then
      allocate(force_occ_basis_type(n_force_occ))
    end if
    if (.not. allocated(force_occ_atom)) then
      allocate(force_occ_atom(n_force_occ))
    end if
    if (.not. allocated(force_occ_spin)) then
      allocate(force_occ_spin(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_n)) then
      allocate(force_occ_basis_n(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_l)) then
      allocate(force_occ_basis_l(n_force_occ))
    end if
    if (.not. allocated(force_occ_max_KS_state)) then
      allocate(force_occ_max_KS_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_m)) then
      allocate(force_occ_basis_m(n_force_occ))
    end if
    if (.not. allocated(force_occ_basis_occupation)) then
      allocate(force_occ_basis_occupation(n_force_occ))
    end if
    if (.not. allocated(which_basis_function)) then
      allocate(which_basis_function(n_force_occ))
    end if
  end if

  end subroutine allocate_force_occupation_basis





!******
!-------------------------------------------------------------------------------
!****s* force_occupation/deallocate_force_occupation_basis
!  NAME
!   deallocate_force_occupation_basis
!  SYNOPSIS

  subroutine deallocate_force_occupation_basis()


  if (n_periodic .eq. 0) then
    if (allocated(force_occ_state)) then
      deallocate(force_occ_state)
    end if
    if (allocated(force_occ_basis_type)) then
      deallocate(force_occ_basis_type)
    end if
    if (allocated(force_occ_atom)) then
      deallocate(force_occ_atom)
    end if
    if (allocated(force_occ_spin)) then
      deallocate(force_occ_spin)
    end if
    if (allocated(force_occ_basis_n)) then
      deallocate(force_occ_basis_n)
    end if
    if (allocated(force_occ_basis_l)) then
      deallocate(force_occ_basis_l)
    end if
    if (allocated(force_occ_max_KS_state)) then
      deallocate(force_occ_max_KS_state)
    end if
    if (allocated(force_occ_basis_m)) then
      deallocate(force_occ_basis_m)
    end if
    if (allocated(force_occ_basis_occupation)) then
      deallocate(force_occ_basis_occupation)
    end if
  else
    if (allocated(which_basis_function)) then
      deallocate(which_basis_function)
    end if
    if (allocated(force_occ_state_periodic)) then
      deallocate(force_occ_state_periodic)
    end if
    if (allocated(force_occ_basis_type)) then
      deallocate(force_occ_basis_type)
    end if
    if (allocated(force_occ_atom)) then
      deallocate(force_occ_atom)
    end if
    if (allocated(force_occ_spin)) then
      deallocate(force_occ_spin)
    end if
    if (allocated(force_occ_basis_n)) then
      deallocate(force_occ_basis_n)
    end if
    if (allocated(force_occ_basis_l)) then
      deallocate(force_occ_basis_l)
    end if
    if (allocated(force_occ_max_KS_state)) then
      deallocate(force_occ_max_KS_state)
    end if
    if (allocated(force_occ_basis_m)) then
      deallocate(force_occ_basis_m)
    end if
    if (allocated(force_occ_basis_occupation)) then
      deallocate(force_occ_basis_occupation)
    end if
  end if

  end subroutine deallocate_force_occupation_basis



!******
!-------------------------------------------------------------------------------
!****s* force_occupation/allocate_force_occupation
!  NAME
!   allocate_force_occupation
!  SYNOPSIS

  subroutine allocate_force_occupation()


  if (n_periodic .eq. 0) then
    if (.not. allocated(force_occ_pr_state)) then
      allocate(force_occ_pr_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_pr_state2)) then
      allocate(force_occ_pr_state2(n_force_occ))
    end if
    if (.not. allocated(force_occ_spin)) then
      allocate(force_occ_spin(n_force_occ))
    end if
    if (.not. allocated(forced_occ_number)) then
      allocate(forced_occ_number(n_force_occ))
    end if
    if (.not. allocated(force_occ_state)) then
      allocate(force_occ_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_max_KS_state)) then
      allocate(force_occ_max_KS_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_min_KS_state)) then
      allocate(force_occ_min_KS_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_maxred_KS_state)) then
      allocate(force_occ_maxred_KS_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_step_KS_state)) then
      allocate(force_occ_step_KS_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_autored)) then
      allocate(force_occ_autored(n_force_occ))
    end if
    if (.not. allocated(force_occ_autorep)) then
      allocate(force_occ_autorep(n_force_occ))
    end if
    if (.not. allocated(fop_detect_auto)) then
      allocate(fop_detect_auto(n_force_occ))
    end if
  else
    if (.not. allocated(force_occ_pr_state)) then
      allocate(force_occ_pr_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_spin)) then
      allocate(force_occ_spin(n_force_occ))
    end if
    if (.not. allocated(forced_occ_number)) then
      allocate(forced_occ_number(n_force_occ))
    end if
    if (.not. allocated(force_occ_max_KS_state)) then
      allocate(force_occ_max_KS_state(n_force_occ))
    end if
    if (.not. allocated(force_occ_min_KS_state)) then
      allocate(force_occ_min_KS_state(n_force_occ))
    end if
  end if

  end subroutine allocate_force_occupation

!******
!-------------------------------------------------------------------------------
!****s* force_occupation/allocate_previous_eigenvector
!  NAME
!   allocate_previous_eigenvector
!  SYNOPSIS

  subroutine allocate_previous_eigenvector()

  if (n_periodic .eq. 0) then
    if (.not.allocated(previous_eigenvector)) then
      allocate(previous_eigenvector(n_basis,n_states,n_spin,n_k_points_task))
    end if
  else
    if (real_eigenvectors)then
      if (.not.allocated(previous_eigenvector)) then
        allocate(previous_eigenvector(n_basis,n_states,n_spin,n_k_points_task))
        previous_eigenvector = 0.d0
      end if
    else
      if (.not.allocated(previous_eigenvector_complex)) then
        allocate(previous_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task))
        previous_eigenvector_complex = 0.d0
      end if

    end if
    if (.not. allocated(force_occ_pr_state_periodic)) then
      allocate(force_occ_pr_state_periodic(n_force_occ,n_k_points))
    end if
    if (.not. allocated(force_occ_state_periodic)) then
      allocate(force_occ_state_periodic(n_force_occ,n_k_points))
    end if
  end if

  end subroutine allocate_previous_eigenvector

!******
!-------------------------------------------------------------------------------
!****s* force_occupation/deallocate_force_occupation
!  NAME
!   deallocate_force_occupation
!  SYNOPSIS

  subroutine deallocate_force_occupation()


  if (n_periodic .eq. 0) then
    if (allocated(force_occ_pr_state)) then
      deallocate(force_occ_pr_state)
    end if
    if (allocated(force_occ_pr_state2)) then
      deallocate(force_occ_pr_state2)
    end if
    if (allocated(force_occ_spin)) then
      deallocate(force_occ_spin)
    end if
    if (allocated(forced_occ_number)) then
      deallocate(forced_occ_number)
    end if
    if (allocated(force_occ_state)) then
      deallocate(force_occ_state)
    end if
    if (allocated(force_occ_min_KS_state)) then
      deallocate(force_occ_min_KS_state)
    end if
    if (allocated(force_occ_max_KS_state)) then
      deallocate(force_occ_max_KS_state)
    end if
    if (allocated(force_occ_maxred_KS_state)) then
      deallocate(force_occ_maxred_KS_state)
    end if
    if (allocated(force_occ_step_KS_state)) then
      deallocate(force_occ_step_KS_state)
    end if
    if (allocated(force_occ_autored)) then
      deallocate(force_occ_autored)
    end if
    if (allocated(force_occ_autorep)) then
      deallocate(force_occ_autorep)
    end if
    if (allocated(fop_detect_auto)) then
      deallocate(fop_detect_auto)
    end if
  else
    if (allocated(force_occ_pr_state)) then
      deallocate(force_occ_pr_state)
    end if
    if (allocated(force_occ_pr_state_periodic)) then
      deallocate(force_occ_pr_state_periodic)
    end if
    if (allocated(force_occ_spin)) then
      deallocate(force_occ_spin)
    end if
    if (allocated(forced_occ_number)) then
      deallocate(forced_occ_number)
    end if
    if (allocated(force_occ_state_periodic)) then
      deallocate(force_occ_state_periodic)
    end if
    if (allocated(force_occ_min_KS_state)) then
      deallocate(force_occ_min_KS_state)
    end if
    if (allocated(force_occ_max_KS_state)) then
      deallocate(force_occ_max_KS_state)
    end if
  end if

  end subroutine deallocate_force_occupation

!******
!-------------------------------------------------------------------------------
!****s* force_occupation/deallocate_previous_eigenvector
!  NAME
!   deallocate_previous_eigenvector
!  SYNOPSIS

  subroutine deallocate_previous_eigenvector()


  if (allocated(previous_eigenvector)) then
    deallocate(previous_eigenvector)
  end if
  if (allocated(previous_eigenvector_complex)) then
    deallocate(previous_eigenvector_complex)
  end if 



  end subroutine deallocate_previous_eigenvector
!******
!-------------------------------------------------------------------------------
!****s*  force_occupation/check_norm_force_occ_p0
!  NAME
!    check_norm_force_occ_p0
!  SYNOPSIS

subroutine check_norm_force_occ_p0(chemical_potential, KS_eigenvalue, n_electrons,&
   occ_numbers, diff_electrons, i_counter)

!  PURPOSE
!  checks and returns the charge density norm (i.e. the electron count) for
!  a given Fermi level, and for whichever smearing type was chosen ...
!
!  USES

!   use dimensions
!   use runtime_choices
  use arch_specific
  use constants
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, intent(in) :: chemical_potential
  real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: KS_eigenvalue
  real*8, intent(in) :: n_electrons
  real*8, dimension(n_states, n_spin,n_k_points), intent(out) :: occ_numbers
  real*8, intent(out) :: diff_electrons
  integer, intent(inout) :: i_counter 
 

!  INPUTS
!   o chemical_potential -- chemical potential
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o n_electrons -- number of electrons
!    
!  OUTPUT
!   o occ_numbers -- occupations of states
!   o diff_electrons -- ?????????
!   o i_counter -- ???????????
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
  real*8,dimension(n_force_occ,n_states,2) :: temp_n_electrons_forced_state
  real*8 :: one_over_ow
  real*8 :: H_odd
  real*8 :: H_even
  real*8 :: hermite_arg
  real*8 :: gauss_weight
  real*8 :: A
  real*8 :: max_exponent
  real*8 :: exp_arg

  !  counter
  integer :: i_state,i_force_occ
  integer :: i_spin
  integer :: i_mp
  integer :: i_k_point,i_k

  
  one_over_ow   = 1.d0 / occupation_width
!  write(use_unit,*) "check_norm... n_spin=", n_spin, "n_states=", n_states
!  do i_spin = 1, n_spin, 1
!     do i_state = 1, n_states, 1
!        write(use_unit,*) i_ifstate, i_spin, KS_eigenvalue(i_state, i_spin)
!        write(use_unit,*) (KS_eigenvalue(i_state, i_spin) - chemical_potential) * one_over_ow
!     end do
  !  end do
  i_counter = i_counter + 1
  temp_n_electrons   = 0.d0
  temp_n_electrons_forced_state = 0.d0
  select case(occupation_type)
  case (0)
     !  gaussian smearing (0 <= occ_numbers <= spin_degenracy)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

             occ_numbers(i_state, i_spin, i_k_point) = spin_degeneracy *  0.5d0 * & 
               (1.d0 - arch_erf((KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential) * one_over_ow))
    
             temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin, i_k_point) * k_weights(i_k_point)
  
           end do
        end do
     end do


  case (1)
     !  fermi smearing (0 <= occ_numbers <= spin_degeneracy)
     max_exponent = maxexponent(chemical_potential) * log(2.d0)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              exp_arg = (KS_eigenvalue(i_state, i_spin,i_k_point) - chemical_potential) * one_over_ow
              if (exp_arg .lt. max_exponent) then

                 occ_numbers(i_state, i_spin, i_k_point) = spin_degeneracy / (1.d0 + exp(exp_arg))

                 temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin, i_k_point) * k_weights(i_k_point)

              else
                 occ_numbers(i_state, i_spin, i_k_point) = 0.d0
              end if
           end do
        end do
     end do

  case (2)
     ! Methfessel-Paxton smearing (0 <= occ_numbers <= spin_degeneracy)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              hermite_arg  = (KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential) * one_over_ow 
              gauss_weight = exp(- hermite_arg * hermite_arg)

              !     zero order contribution
              occ_numbers(i_state, i_spin,i_k_point) = 0.5d0 * (1.d0 - arch_erf(hermite_arg))

              if (n_methfessel_paxton .gt. 0) then
             
                 !     first order contribution
                 A = - 0.25d0 * pisqrt_inv
                 !     H_even = H_0 = 1
                 H_even = 1.d0
                 !     H_odd =  H_1 = 2 * x
                 H_odd  = 2d0 * hermite_arg
                 occ_numbers(i_state, i_spin,i_k_point) = occ_numbers(i_state, i_spin, i_k_point) + A * H_odd * gauss_weight
              end if

              if (n_methfessel_paxton .gt. 1) then
                 do i_mp = 2, n_methfessel_paxton, 1

                    A = - 1.d0 / dble(4d0 * i_mp) * A
                    H_even = 2d0 * hermite_arg * H_odd  - 2d0 *  dble(i_mp)      * H_even
                    H_odd  = 2d0 * hermite_arg * H_even - 2d0 * (dble(i_mp) + 1d0) * H_odd
                    occ_numbers(i_state, i_spin, i_k_point) = occ_numbers(i_state, i_spin,i_k_point) + A * H_odd * gauss_weight
                 end do
              end if


           occ_numbers(i_state, i_spin,i_k_point) = occ_numbers(i_state, i_spin,i_k_point) &
                * spin_degeneracy 

           temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin,i_k_point) * k_weights(i_k_point)
        end do
     end do
  end do

  case default
     write(use_unit,*) "* Unknown smearing type in subroutine check_norm."
     write(use_unit,*) "* Abort."
     stop
     
  end select

!   write(use_unit,*)"IN CHECK NORM HALLO"

! calculating n_electrons in force_occ_state_periodic before and after changing occ_numbers

  do i_force_occ = 1, n_force_occ, 1
    do i_state = 1, n_states, 1
      do i_k_point = 1, n_k_points, 1

        temp_n_electrons_forced_state(i_force_occ,i_state,1) = temp_n_electrons_forced_state(i_force_occ,i_state,1) + &
          occ_numbers(i_state,1,i_k_point) * &
          k_weights(i_k_point)

        if (n_spin > 1) then 
          temp_n_electrons_forced_state(i_force_occ,i_state,2) = temp_n_electrons_forced_state(i_force_occ,i_state,2) + &
            occ_numbers(i_state,2,i_k_point) * &
            k_weights(i_k_point)
        end if
      end do
    end do
  end do


  if (n_periodic .eq. 0) then
    do i_force_occ = 1, n_force_occ, 1
      temp_n_electrons = temp_n_electrons - occ_numbers(force_occ_state(i_force_occ),&
          force_occ_spin(i_force_occ),1)*k_weights(1)
      occ_numbers(force_occ_state(i_force_occ),force_occ_spin(i_force_occ),1) = &
          forced_occ_number(i_force_occ)
      temp_n_electrons = temp_n_electrons + occ_numbers(force_occ_state(i_force_occ), &
          force_occ_spin(i_force_occ),1)*k_weights(1)
    end do
  else

    do i_force_occ = 1, n_force_occ, 1
      do i_k_point = 1, n_k_points, 1                                     

        temp_n_electrons = temp_n_electrons - occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point),&
            force_occ_spin(i_force_occ),i_k_point)*k_weights(i_k_point)
        occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point),force_occ_spin(i_force_occ),i_k_point) = &
            forced_occ_number(i_force_occ)
        temp_n_electrons = temp_n_electrons + occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point), &
            force_occ_spin(i_force_occ),i_k_point)*k_weights(i_k_point)

      end do
      
    end do
  end if 


!   temp_n_electrons_forced_state = 0.d0
!   do i_force_occ = 1, n_force_occ, 1
!     do i_state = 1, n_states, 1
!       do i_k_point = 1, n_k_points, 1
! 
!     temp_n_electrons_forced_state(i_force_occ,i_state,1) = temp_n_electrons_forced_state(i_force_occ,i_state,1) + &
!       occ_numbers(i_state,1,i_k_point) * &
!       k_weights(i_k_point)
! 
!     temp_n_electrons_forced_state(i_force_occ,i_state,2) = temp_n_electrons_forced_state(i_force_occ,i_state,2) + &
!       occ_numbers(i_state,2,i_k_point) * &
!       k_weights(i_k_point)
! 
!       end do
!     end do
!   end do






  diff_electrons = temp_n_electrons - n_electrons
  
end subroutine check_norm_force_occ_p0
!******
!-------------------------------------------------------------------------------
!****s*  force_occupation/check_norm_force_occ_p0_smearing
!  NAME
!    check_norm_force_occ_p0_smearing
!  SYNOPSIS

subroutine check_norm_force_occ_p0_smearing(chemical_potential, KS_eigenvalue, n_electrons,&
   occ_numbers, diff_electrons, i_counter)

!  PURPOSE
!  checks and returns the charge density norm (i.e. the electron count) for
!  a given Fermi level, and for whichever smearing type was chosen and more 
!  importantly, applies constraints plus a smearing which allows to deal with 
!  degenerate KS states
!
!  USES

!   use dimensions
!   use runtime_choices
  use arch_specific
  use constants
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, intent(in) :: chemical_potential
  real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: KS_eigenvalue
  real*8, intent(in) :: n_electrons
  real*8, dimension(n_states, n_spin,n_k_points), intent(out) :: occ_numbers
  real*8, intent(out) :: diff_electrons
  integer, intent(inout) :: i_counter 
 

!  INPUTS
!   o chemical_potential -- chemical potential
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o n_electrons -- number of electrons
!    
!  OUTPUT
!   o occ_numbers -- occupations of states
!   o diff_electrons -- ?????????
!   o i_counter -- ???????????
! 
!  AUTHOR
!    Reinhard J. Maurer, Technical University Munich
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2013).
!  SOURCE
!

  
  !  local variables
  real*8 :: temp_n_electrons
  real*8,dimension(n_force_occ,n_states,2) :: temp_n_electrons_forced_state
  real*8 :: one_over_ow
  real*8 :: H_odd
  real*8 :: H_even
  real*8 :: hermite_arg
  real*8 :: gauss_weight
  real*8 :: A
  real*8 :: max_exponent
  real*8 :: exp_arg

  !  counter
  integer :: i_state,i_force_occ
  integer :: i_spin
  integer :: i_mp
  integer :: i_k_point,i_k

  ! Number of electrons in constraint_regime
  real*8  :: n_electrons_constr_regime
  real*8  :: delta_n_electrons
  integer :: n_state_regime
  real*8  :: regime_part_func
  real*8  :: constr_regime_width
  real*8  :: regime_l
  real*8  :: regime_r

  one_over_ow   = 1.d0 / occupation_width
  i_counter = i_counter + 1
  temp_n_electrons   = 0.d0
  temp_n_electrons_forced_state = 0.d0
  select case(occupation_type)
  case (0)
     !  gaussian smearing (0 <= occ_numbers <= spin_degenracy)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              occ_numbers(i_state, i_spin, i_k_point) = spin_degeneracy *  0.5d0 * & 
                  (1.d0 - arch_erf((KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential) * one_over_ow))
    
              temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin, i_k_point) * k_weights(i_k_point)
  
           end do
        end do
     end do


  case (1)
     !  fermi smearing (0 <= occ_numbers <= spin_degeneracy)
     max_exponent = maxexponent(chemical_potential) * log(2.d0)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              exp_arg = (KS_eigenvalue(i_state, i_spin,i_k_point) - chemical_potential) * one_over_ow
              if (exp_arg .lt. max_exponent) then

                 occ_numbers(i_state, i_spin, i_k_point) = spin_degeneracy / (1.d0 + exp(exp_arg))

                 temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin, i_k_point) * k_weights(i_k_point)

              else
                 occ_numbers(i_state, i_spin, i_k_point) = 0.d0
              end if
           end do
        end do
     end do

  case (2)
     ! Methfessel-Paxton smearing (0 <= occ_numbers <= spin_degeneracy)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              hermite_arg  = (KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential) * one_over_ow 
              gauss_weight = exp(- hermite_arg * hermite_arg)

              !     zero order contribution
              occ_numbers(i_state, i_spin,i_k_point) = 0.5d0 * (1.d0 - arch_erf(hermite_arg))

              if (n_methfessel_paxton .gt. 0) then
             
                 !     first order contribution
                 A = - 0.25d0 * pisqrt_inv
                 !     H_even = H_0 = 1
                 H_even = 1.d0
                 !     H_odd =  H_1 = 2 * x
                 H_odd  = 2d0 * hermite_arg
                 occ_numbers(i_state, i_spin,i_k_point) = occ_numbers(i_state, i_spin, i_k_point) + A * H_odd * gauss_weight
              end if

              if (n_methfessel_paxton .gt. 1) then
                 do i_mp = 2, n_methfessel_paxton, 1

                    A = - 1.d0 / dble(4d0 * i_mp) * A
                    H_even = 2d0 * hermite_arg * H_odd  - 2d0 *  dble(i_mp)      * H_even
                    H_odd  = 2d0 * hermite_arg * H_even - 2d0 * (dble(i_mp) + 1d0) * H_odd
                    occ_numbers(i_state, i_spin, i_k_point) = occ_numbers(i_state, i_spin,i_k_point) + A * H_odd * gauss_weight
                 end do
              end if


           occ_numbers(i_state, i_spin,i_k_point) = occ_numbers(i_state, i_spin,i_k_point) &
                * spin_degeneracy 

           temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin,i_k_point) * k_weights(i_k_point)
        end do
     end do
  end do

  case default
     write(use_unit,*) "* Unknown smearing type in subroutine check_norm."
     write(use_unit,*) "* Abort."
     stop
     
  end select

!Measure the amount of electrons in the constraint_regime (+-force_occupation_smearing_width*3)
  do i_force_occ = 1, n_force_occ, 1
    regime_part_func=0d0
    n_electrons_constr_regime=0d0
    delta_n_electrons=0d0
    n_state_regime=0
    regime_l=KS_eigenvalue(force_occ_state(i_force_occ),force_occ_spin(i_force_occ) ,1)&
             -3d0*force_occupation_smearing_width
    regime_r=KS_eigenvalue(force_occ_state(i_force_occ),force_occ_spin(i_force_occ) ,1)&
             +3d0*force_occupation_smearing_width
    do i_state = 1, n_states, 1
      if (KS_eigenvalue(i_state, force_occ_spin(i_force_occ),1)>regime_l .and. KS_eigenvalue(i_state,&
          force_occ_spin(i_force_occ), 1)<regime_r) then
        n_electrons_constr_regime = n_electrons_constr_regime + occ_numbers(i_state,&
          force_occ_spin(i_force_occ),1)
        n_state_regime = n_state_regime + 1
        regime_part_func= regime_part_func + exp(-(KS_eigenvalue(i_state, force_occ_spin(i_force_occ),1)-&
          KS_eigenvalue(force_occ_state(i_force_occ), force_occ_spin(i_force_occ),1))**2/&
          (2*(force_occupation_smearing_width**2)))
      end if
    end do
    !We know the number of electrons, now we have to know how many electrons to erase
    delta_n_electrons = forced_occ_number(i_force_occ)&
                        - occ_numbers(force_occ_state(i_force_occ),force_occ_spin(i_force_occ),1)

     !write(use_unit,*) "--------------Force_OCC_SMEAR_DEBUG-------------------"
     !write(use_unit,*) "N_CONSTRAINT ", i_force_occ
     !write(use_unit,*) "N_ELECTRONS_CONSTR_REGIME ", n_electrons_constr_regime
     !write(use_unit,*) "DELTA_N_ELECTRONS ", delta_n_electrons
     !write(use_unit,*) "N_STATE_REGIME ", n_state_regime

    ! Applying constraint
    temp_n_electrons = temp_n_electrons - n_electrons_constr_regime
    do i_state = 1,n_states, 1
      if (KS_eigenvalue(i_state, force_occ_spin(i_force_occ), 1)>regime_l .and. KS_eigenvalue(i_state,&
          force_occ_spin(i_force_occ), 1)<regime_r) then
        occ_numbers(i_state, force_occ_spin(i_force_occ),1) = occ_numbers(i_state, force_occ_spin(i_force_occ),1)&
           + delta_n_electrons * ( exp(-(KS_eigenvalue(i_state, force_occ_spin(i_force_occ),1)-&
           KS_eigenvalue(force_occ_state(i_force_occ), force_occ_spin(i_force_occ),1))**2/&
          (2*(force_occupation_smearing_width**2))) / regime_part_func)

        if (occ_numbers(i_state, force_occ_spin(i_force_occ),1) > spin_degeneracy) then
          occ_numbers(i_state,force_occ_spin(i_force_occ),1) = spin_degeneracy
          temp_n_electrons = temp_n_electrons + spin_degeneracy
        elseif (occ_numbers(i_state, force_occ_spin(i_force_occ),1) < 0d0) then
          occ_numbers(i_state,force_occ_spin(i_force_occ),1) = 0d0
        else
          temp_n_electrons = temp_n_electrons + occ_numbers(i_state,force_occ_spin(i_force_occ),1)
        end if
!          write(use_unit,*) "N_state ", i_state
!          write(use_unit,*) "energy ", KS_eigenvalue(i_state, force_occ_spin(i_force_occ),1)
!          write(use_unit,*) "occupation ", occ_numbers(i_state,force_occ_spin(i_force_occ),1)
      end if
    end do
    !temp_n_electrons = temp_n_electrons + delta_n_electrons
!      write(use_unit,*) "temp_n_electrons ", temp_n_electrons, " ", n_electrons
!      write(use_unit,*) "---------------Force_OCC_SMEAR_DEBUG------------"
  end do

  diff_electrons = temp_n_electrons - n_electrons

end subroutine check_norm_force_occ_p0_smearing
!******
!-------------------------------------------------------------------------------
!****s* force_occupation/check_norm_force_occ_v2
!  NAME
!    check_norm_force_occ_v2
!  SYNOPSIS

subroutine check_norm_force_occ_v2(chemical_potential, KS_eigenvalue, n_electrons,&
            occ_numbers, diff_electrons, i_counter,current_spin)

!  PURPOSE
!  checks and returns the charge density norm (i.e. the electron count) for
!  a given Fermi level, and for whichever smearing type was chosen ...!
!
!  USES
! 
!   use dimensions
!   use runtime_choices
  use arch_specific
  use constants
  implicit none

!  ARGUMENTS

  real*8, intent(in) :: chemical_potential
  real*8, dimension(n_states), intent(in) :: KS_eigenvalue
  real*8, intent(in) :: n_electrons
  real*8, dimension(n_states), intent(out) :: occ_numbers
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
!   o diff_electrons -- ??????
!   o i_counter -- ???????
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
     
  !  counter
  integer :: i_state,i_force_occ
  integer :: i_spin
  integer :: i_mp
  
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
     do i_state = 1, n_states, 1
        occ_numbers(i_state) = spin_degeneracy * 0.5d0 * & 
             (1 - arch_erf((KS_eigenvalue(i_state) - chemical_potential) * one_over_ow))
        temp_n_electrons = temp_n_electrons + occ_numbers(i_state)
     end do

  case (1)
     !  fermi smearing (0 <= occ_numbers <= spin_degeneracy)
     max_exponent = maxexponent(chemical_potential) * log(2.d0)
        do i_state = 1, n_states, 1
           exp_arg = (KS_eigenvalue(i_state) - chemical_potential) * one_over_ow
           if (exp_arg .lt. max_exponent) then
              occ_numbers(i_state) = spin_degeneracy / (1 + exp(exp_arg))
              temp_n_electrons = temp_n_electrons + occ_numbers(i_state)
           else
              occ_numbers(i_state) = 0.d0
           end if
        end do
     
  case (2)
     ! Methfessel-Paxton smearing (0 <= occ_numbers <= spin_degeneracy)
     do i_state = 1, n_states, 1
        hermite_arg  = (KS_eigenvalue(i_state) - chemical_potential) * one_over_ow 
        gauss_weight = exp(- hermite_arg * hermite_arg)
        !     zero order contribution
        occ_numbers(i_state) = 0.5 * (1 - arch_erf(hermite_arg))
        if (n_methfessel_paxton .gt. 0) then
           !    write(use_unit,*) A, H_even, H_odd, hermite_arg
           !     first order contribution
              A = - 0.25 * pisqrt_inv
              !     H_even = H_0 = 1
              H_even = 1.d0
              !     H_odd =  H_1 = 2 * x
              H_odd  = 2 * hermite_arg
              occ_numbers(i_state) = occ_numbers(i_state) + A * H_odd * gauss_weight
           end if
           if (n_methfessel_paxton .gt. 1) then
              do i_mp = 2, n_methfessel_paxton, 1
                 A = - 1. / dble(4 * i_mp) * A
                 H_even = 2 * hermite_arg * H_odd  - 2 *  i_mp      * H_even
                 H_odd  = 2 * hermite_arg * H_even - 2 * (i_mp + 1) * H_odd
                 occ_numbers(i_state) = occ_numbers(i_state) + A * H_odd * gauss_weight
              end do
           end if
!           write(use_unit,*) occ_numbers(i_state, i_spin)
           occ_numbers(i_state) = occ_numbers(i_state) * spin_degeneracy
!           write(use_unit,*) occ_numbers(i_state, i_spin)
           temp_n_electrons = temp_n_electrons + occ_numbers(i_state)
        end do
!     write(use_unit,*) "temp_n_electrons= ", temp_n_electrons, " chemical_potential= ", chemical_potential
  case default
     write(use_unit,*) "* Unknown smearing type in subroutine check_norm."
     write(use_unit,*) "* Abort."
     stop
     
  end select

  do i_force_occ = 1, n_force_occ, 1
    if (current_spin .eq. force_occ_spin(i_force_occ)) then
      temp_n_electrons = temp_n_electrons - occ_numbers(force_occ_state(i_force_occ))
      occ_numbers(force_occ_state(i_force_occ)) = &
          forced_occ_number(i_force_occ)
      temp_n_electrons = temp_n_electrons + occ_numbers(force_occ_state(i_force_occ)) 
    end if
  end do

  diff_electrons = temp_n_electrons - n_electrons
  
end subroutine check_norm_force_occ_v2





!******
!-------------------------------------------------------------------------------
!****s* force_occupation/check_norm_force_occ_periodic_v2
!  NAME
!    check_norm_force_occ_periodic_v2
!  SYNOPSIS

subroutine check_norm_force_occ_periodic_v2(chemical_potential, KS_eigenvalue, n_electrons,&
            occ_numbers, diff_electrons, i_counter,current_spin)

!  PURPOSE
!  checks and returns the charge density norm (i.e. the electron count) for
!  a given Fermi level, and for whichever smearing type was chosen ...!
!
!  USES
! 
!   use dimensions
!   use runtime_choices
  use arch_specific
  use constants
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, intent(in) :: chemical_potential
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
!   o diff_electrons -- ??????
!   o i_counter -- ???????
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
     
  !  counter
  integer :: i_state,i_force_occ,i_k_point
  integer :: i_spin
  integer :: i_mp
  
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
!           write(use_unit,*) occ_numbers(i_state, i_spin)
           occ_numbers(i_state,i_k_point) = occ_numbers(i_state,i_k_point) * spin_degeneracy
!           write(use_unit,*) occ_numbers(i_state, i_spin)
           temp_n_electrons = temp_n_electrons + occ_numbers(i_state,i_k_point)*k_weights(i_k_point)
        end do
!     write(use_unit,*) "temp_n_electrons= ", temp_n_electrons, " chemical_potential= ", chemical_potential
      end do

  case default
     write(use_unit,*) "* Unknown smearing type in subroutine check_norm."
     write(use_unit,*) "* Abort."
     stop
     
  end select

  do i_force_occ = 1, n_force_occ, 1
    if (current_spin .eq. force_occ_spin(i_force_occ)) then
      do i_k_point = 1, n_k_points, 1                                     
        temp_n_electrons = temp_n_electrons - occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point),&
            i_k_point)*k_weights(i_k_point)
        occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point),i_k_point) = &
            forced_occ_number(i_force_occ)
        temp_n_electrons = temp_n_electrons + occ_numbers(force_occ_state_periodic(i_force_occ,i_k_point), &
            i_k_point)*k_weights(i_k_point)
      end do
    end if
  end do







  diff_electrons = temp_n_electrons - n_electrons
  
end subroutine check_norm_force_occ_periodic_v2















!******
!-------------------------------------------------------------------------------
!****s* force_occupation/check_norm_constraint_occ
!  NAME
!   check_norm_constraint_occ
!  SYNOPSIS

subroutine check_norm_constraint_occ(chemical_potential, KS_eigenvalue, constraint_electrons, &
                                 constraint_proj, occ_numbers, diff_electrons, i_counter,&
                                current_spin)

!  PURPOSE
!   checks and returns the charge density norm (i.e. the electron count) for
!   a given Fermi level, and for whichever smearing type was chosen ...
!  USES

!   use dimensions
!   use runtime_choices
  use arch_specific
  use constants
  implicit none

!  ARGUMENTS

  real*8, intent(in) :: chemical_potential
  real*8, dimension(n_states), intent(in) :: KS_eigenvalue
  real*8, dimension(n_states), intent(in) :: constraint_proj  
  real*8, intent(in) :: constraint_electrons
  real*8, dimension(n_states), intent(out) :: occ_numbers
  real*8, intent(out) :: diff_electrons
  integer, intent(inout) :: i_counter 

  integer :: current_spin


!  INPUTS
!   o chemical_potential -- chemical potential
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o constraint_proj -- ?????????
!   o constraint_electrons -- number of nostrain electrons
!
!  OUTPUT
!   o occ_numbers -- occupation numbers
!   o diff_electrons -- ???????????
!   o i_counter -- ?????????????
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

     
  !  counter
  integer :: i_state,i_force_occ
  integer :: i_spin

  integer :: i_mp
  
  i_counter = i_counter + 1
  one_over_ow   = 1. / occupation_width
  temp_n_electrons   = 0.d0
  
  select case(occupation_type)
  case (0)
     !  gaussian smearing (0 <= occ_numbers )
!     do i_spin =1, n_spin, 1
        do i_state = 1, n_states, 1
            occ_numbers(i_state) = spin_degeneracy * 0.5d0 * &
           (1 - arch_erf((KS_eigenvalue(i_state) - chemical_potential) * one_over_ow))
           temp_n_electrons = temp_n_electrons + constraint_proj(i_state)*occ_numbers(i_state)
        end do
!     end do

  case (1)
     !  fermi smearing (0 <= occ_numbers )
!     do i_spin =1, n_spin, 1
        do i_state = 1, n_states, 1
           occ_numbers(i_state) = spin_degeneracy / &
           (1 + exp((KS_eigenvalue(i_state) - chemical_potential) * one_over_ow))
           temp_n_electrons = temp_n_electrons + occ_numbers(i_state) * constraint_proj(i_state)
        end do
!     end do
  case (2)
     ! Methfessel-Paxton smearing (0 <= occ_numbers <= pin_degeneracy)
!     do i_spin =1, n_spin, 1
        do i_state = 1, n_states, 1
           hermite_arg  = (KS_eigenvalue(i_state) - chemical_potential) * one_over_ow 
           gauss_weight = exp(- hermite_arg * hermite_arg)
           !     zero order contribution
           occ_numbers(i_state) = 0.5 * (1 - arch_erf(hermite_arg))
           if (n_methfessel_paxton .gt. 0) then
!              write(use_unit,*) A, H_even, H_odd, hermite_arg
              !     first order contribution
              A = - 0.25 * pisqrt_inv
              !     H_even = H_0 = 1
              H_even = 1.d0
              !     H_odd =  H_1 = 2 * x
              H_odd  = 2 * hermite_arg
              occ_numbers(i_state) = occ_numbers(i_state) + A * H_odd * gauss_weight
           end if
           if (n_methfessel_paxton .gt. 1) then
              do i_mp = 2, n_methfessel_paxton, 1
                 A = - 1. / dble(4 * i_mp) * A
                 H_even = 2 * hermite_arg * H_odd  - 2 *  i_mp      * H_even
                 H_odd  = 2 * hermite_arg * H_even - 2 * (i_mp + 1) * H_odd
                 occ_numbers(i_state) = occ_numbers(i_state) + A * H_odd * gauss_weight
              end do
           end if
!     factor 2 due to spin
!            occ_numbers(i_state,i_spin) = occ_numbers(i_state,i_spin) * spin_degeneracy
!            write(use_unit,*) "occ = ", occ_numbers(i_state,i_spin)
            temp_n_electrons = temp_n_electrons + occ_numbers(i_state) * constraint_proj(i_state)
         end do
!     end do
 
  case default
         write(use_unit,*) "* Unknown smearing type in subroutine check_norm_constraint."
         write(use_unit,*) "* Abort."
         stop

  end select


  do i_force_occ = 1, n_force_occ, 1
    if (current_spin .eq. force_occ_spin(i_force_occ)) then
      temp_n_electrons = temp_n_electrons - occ_numbers(force_occ_state(i_force_occ))
      occ_numbers(force_occ_state(i_force_occ)) = &
          forced_occ_number(i_force_occ)
      temp_n_electrons = temp_n_electrons + occ_numbers(force_occ_state(i_force_occ)) 
    end if
  end do



    diff_electrons = temp_n_electrons - constraint_electrons

  

end subroutine check_norm_constraint_occ
!******
!-------------------------------------------------------------------------------
!****s* force_occupation/f_external_constraint_occ
!  NAME
!    f_external_constraint_occ
!  SYNOPSIS 
      subroutine f_external_constraint_occ(n_arguments, arguments, &
           f_value,ifn, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )
!  PURPOSE
!    Provides the target function for the constrained DFT calculations as the
!    l2-norm of the difference between the actual number of electrons in the regions and
!    given target number of electrons in the regions. This routine is passed externally 
!    to the BFGS-optimizer.
!  USES
      use constants
      use constraint
!       use dimensions
      use localorb_io
      use mpi_tasks, only: aims_stop

      implicit none
!  ARGUMENTS
      real*8 :: f_value
      integer :: n_arguments
      real*8, dimension(3*n_arguments) :: arguments
      integer :: ifn
      real*8 hamiltonian( n_basis*(n_basis+1)/2, n_spin )
      real*8 overlap_matrix( n_basis*(n_basis+1)/2 )
      real*8, dimension(n_states, n_spin) :: KS_eigenvalue
      real*8, dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
      real*8, dimension(n_states, n_spin) :: occ_numbers
      real*8 :: n_electrons
      real*8 :: chemical_potential
      real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
      real*8, dimension(n_region, n_spin) :: electrons_in_region
      integer :: current_spin


!  INPUTS
!    o n_arguments -- number of arguments to f
!    o arguments -- the array of arguments to f
!    o ifn -- function call counter
!    o hamiltonian -- the Hamilton matrix
!    o overlap_matrix -- the overlap matrix
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o occ_numbers -- occupation numbers
!    o n_electrons -- number of electrons in the spin channel
!    o chemical_potential -- the chemical potential
!    o constraint_proj -- projection to the constrained regions
!    o electrons_in_region -- number of electrons in the regions
!    o current_spin -- current spin channel
!  OUTPUT
!    o f_value -- value of the target function f
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


!     locals

      character*100 :: info_str

      real*8 :: avg_zero
      integer :: i_region, i_hamiltonian
      logical :: t_out = .false.

      real*8 :: avg_of_arguments
      character(*), parameter :: func = 'f_external_constraint_occ'

!     begin work

      if (constraint_debug) then
        call localorb_info('',use_unit,'(A)')
        write(info_str,'(2X,A,I5,A)') &
        "Constrained DFT: Eigenvector update number ", &
        ifn,"."
        call localorb_info(info_str,use_unit,'(A)')
        call localorb_info('',use_unit,'(A)')
      end if

      avg_of_arguments = 0.0d0

      ! FIXME EITHER I DO NOT UNDERSTAND THIS AVERAGING OR IT IS SIMPLY INCORRECT.

      do i_region = 1, n_arguments, 1
         avg_of_arguments = avg_of_arguments + arguments(i_region) * &
                            constraint_electrons(i_region,current_spin)
      enddo
      avg_of_arguments = avg_of_arguments / n_electrons

      do i_region = 1, n_active_regions - 1, 1
         constraint_potential(i_region,current_spin) = &
              arguments(i_region) - avg_of_arguments
      enddo

      constraint_potential(n_active_regions,current_spin) = &
           - avg_of_arguments

      !  Restore original hamiltonian for current spin component
      do i_hamiltonian=1, n_basis*(n_basis+1)/2, 1
            hamiltonian_work(i_hamiltonian,current_spin) = &
                 hamiltonian(i_hamiltonian,current_spin)
      enddo

      ! ... and add the current constraint potentials.
      call add_constraint_potentials_v2 &
           ( current_spin, overlap_matrix, hamiltonian_work &
           )

      if (constraint_debug) then
                   write(info_str,'(2X,A)') &
                        ''
                   write(info_str,'(2X,A)') &
                        'Region                    Constraint potential'
                   call localorb_info(info_str,use_unit,'(A)')
                   do i_region = 1, n_active_regions, 1
                     write(info_str,'(2X,I5,20X,F20.15)') &
                     i_region, &
                     constraint_potential(i_region,current_spin)*hartree
                     call localorb_info(info_str,use_unit,'(A)')
                   enddo
                   write(info_str,'(2X,A)') &
                        ''
      end if

      ! solve for new KS eigenstates and -energies ...
      if (n_k_points > 1) call aims_stop('Only implemented for 1 k-point', func)
      call improve_real_eigenfunctions &
           ( overlap_matrix, hamiltonian_work(:,current_spin), &
           t_out, &
           KS_eigenvalue(:,current_spin), &
           KS_eigenvector(:,:,current_spin), 1 &
           )

      ! ... and obtain the updated global Fermi level and
      ! and global occupation numbers

      call determine_corehole_projection(overlap_matrix,KS_eigenvector)

      ! Notice n_electrons is here the number of electrons in the current spin channel
      call get_occupation_numbers_occ_v2 &
           ( KS_eigenvalue(:,current_spin), n_electrons, &
           t_out, &
           occ_numbers(:,current_spin), chemical_potential, &
           current_spin)

!test
!      write(use_unit,*) "occ_numbers: "
!      write(use_unit,*) occ_numbers
!      write(use_unit,*)
!test end

      ! update constraint projectors and count number of electrons in each region, each spin
      call get_electrons_per_region_v2 &
           ( current_spin, KS_eigenvector, overlap_matrix, occ_numbers, &
           constraint_proj, &
           electrons_in_region &
           )

      f_value = 0.0d0
      do i_region = 1, n_active_regions, 1
         f_value = f_value + &
              ( electrons_in_region(i_region,current_spin) - &
              constraint_electrons(i_region,current_spin) )**2.d0
      enddo

      f_value = dsqrt(f_value)

      end subroutine f_external_constraint_occ
!******
!-------------------------------------------------------------------------------
!****s* force_occupation/determine_corehole_projection
!  NAME
!    determine_corehole_projection
!  SYNOPSIS 


  subroutine determine_corehole_projection(overlap_matrix,KS_eigenvector)
  
!   use dimensions
  use timing,                  only: number_of_loops
  use mpi_utilities
  
  !  ARGUMENTS
    real*8 overlap_matrix( n_hamiltonian_matrix_size )
    real*8,     dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
    real*8,     dimension(n_states, n_spin) :: KS_eigenvalue
  
  ! local variables
  
    integer :: i_force_occ, reduction
    integer :: i_basis_1,i_basis_2,i_index_real,i_state
    
    real*8 :: biggest,temp

    character(*), parameter :: func = 'determine_corehole_projection'
  ! output
  
 
    do i_force_occ = 1, n_force_occ, 1
      do while (biggest .lt. 0.1d0)
        biggest = 0d0
        do i_state = force_occ_min_KS_state(i_force_occ), force_occ_max_KS_state(i_force_occ), 1
          temp = 0.d0
          do i_basis_1 = 1, n_basis, 1
            do i_basis_2 = i_basis_1+1, n_basis, 1
              i_index_real = (i_basis_2-1)*i_basis_2/2 +  i_basis_1
              temp = temp + &
                previous_eigenvector(i_basis_1,force_occ_pr_state(i_force_occ),&
                force_occ_spin(i_force_occ),1) * &
                KS_eigenvector(i_basis_2,i_state,force_occ_spin(i_force_occ)) * &
                overlap_matrix(i_index_real)
            end do
          end do
          temp = 2.0*temp
          do i_basis_1 = 1, n_basis, 1
            i_index_real = i_basis_1*(1+i_basis_1)/2
            temp = temp + &
              previous_eigenvector(i_basis_1,force_occ_pr_state(i_force_occ),&
              force_occ_spin(i_force_occ),1) * &
              KS_eigenvector(i_basis_1,i_state,force_occ_spin(i_force_occ)) * &
              overlap_matrix(i_index_real)
          end do
          if (abs(temp).gt.biggest) then
  !           write(use_unit,*)"IN PR-ROUTINE",abs(temp),biggest,i_state,force_occ_pr_state
            force_occ_state(i_force_occ) = i_state
            biggest = abs(temp)
          end if
        end do
        ! do this test for every single force_occupation_projector
        if (myid.eq.0) then
          write(use_unit,*) "BIG", biggest
          write(use_unit,*) "FORCED TO: now -- previous" ,force_occ_state(i_force_occ), force_occ_pr_state(i_force_occ)
        end if
        if (biggest .lt. 0.1d0) then
          !call aims_stop('MOM overlap is less than 0.1d0, aborting since this will be garbage anyway.', func)
          if (myid.eq.0) then
            write(use_unit,*) "WARNING: MOM overlap .lt. threshold (0.1d0): force_occ_max_KS_state += 1"
          end if
          force_occ_max_KS_state(i_force_occ) = force_occ_max_KS_state(i_force_occ) + 1
        end if
      end do
      biggest = 0.d0
    end do

    if (force_occupation_projector_auto) then
      ! automatic oscillation detection
      do i_force_occ = 1, n_force_occ, 1
        if (force_occ_autored(i_force_occ) .lt. number_of_loops) then
          if ((force_occ_state(i_force_occ) .eq. force_occ_pr_state2(i_force_occ)) .and. &
             (force_occ_state(i_force_occ) .ne. force_occ_pr_state(i_force_occ))) then
            fop_detect_auto(i_force_occ) = fop_detect_auto(i_force_occ) + 1
            if (force_occ_autorep(i_force_occ) .lt. fop_detect_auto(i_force_occ)) then
              ! we have detected an oscillation, now decrease to force_occ_max_KS_state by one and reset counter
              force_occ_maxred_KS_state(i_force_occ) = force_occ_max_KS_state(i_force_occ) - 1
              force_occ_step_KS_state(i_force_occ) = force_occ_autored(i_force_occ)
              if (myid.eq.0) then
                write (use_unit,'(2X,A,I5,A,I3)') "| SCF [ ", number_of_loops, " ] Detected oscillation in projector ", i_force_occ
                ! Reset the counter
                fop_detect_auto(i_force_occ) = 0
              end if
              ! Decrease the MOM-window upper bound until we're not oscillating anymore
              reduction = force_occ_max_KS_state(i_force_occ) - MAX(force_occ_pr_state(i_force_occ), force_occ_state(i_force_occ)) + 1
              call force_occ_reduce_subspace(number_of_loops, i_force_occ, reduction)
            end if
          else
            fop_detect_auto(i_force_occ) = 0
          end if
        end if
      end do
    end if
!     

  end subroutine determine_corehole_projection



!******
!-------------------------------------------------------------------------------
!****s* force_occupation/determine_corehole_projection_periodic
!  NAME
!    determine_corehole_projection_periodic
!  SYNOPSIS 


  subroutine determine_corehole_projection_periodic(overlap_matrix_w_complex,overlap_matrix_w,&
                                                   KS_eigenvector_complex,KS_eigenvector,i_k,i_k_point)
  
!   use dimensions
!   use runtime_choices
  use mpi_utilities
  
  !  ARGUMENTS
    complex*16 :: overlap_matrix_w_complex( n_basis*(n_basis+1)/2 )
    real*8     :: overlap_matrix_w( n_basis*(n_basis+1)/2 )
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task) ::  KS_eigenvector_complex
    real*8, dimension(n_basis, n_states, n_spin,n_k_points_task) ::  KS_eigenvector
    integer :: i_k
    integer :: i_k_point
  
  ! local variables
  
    complex*16, dimension(n_basis,n_basis):: ovlp_work_full
    integer :: i_force_occ,i_index
    integer :: i_basis_1,i_basis_2,i_index_real,i_state
    
    real*8 :: biggest
    complex*16 :: temp

    character*140 :: info_str
  ! output
  
  
!     i_index = 0
!     do i_basis_2 = 1,n_basis, 1
!     do i_basis_1 = 1,i_basis_2,1
!       
!       i_index = i_index + 1
!       ovlp_work_full(i_basis_1, i_basis_2) = overlap_matrix_w_complex(i_index)
!       
!     end do
!     end do

    do i_force_occ = 1, n_force_occ, 1
      biggest = 0d0
      do i_state = force_occ_min_KS_state(i_force_occ), force_occ_max_KS_state(i_force_occ), 1
        temp = 0.d0
        i_index = 0
        do i_basis_2 = 1, n_basis, 1
          do i_basis_1 = 1, i_basis_2, 1
            i_index = i_index + 1
            if (real_eigenvectors) then
              temp = temp + &
                previous_eigenvector(i_basis_2,&
                force_occ_pr_state_periodic(i_force_occ,i_k_point),&
                force_occ_spin(i_force_occ),i_k) * &
                KS_eigenvector(i_basis_1,i_state,force_occ_spin(i_force_occ),i_k) * &
                overlap_matrix_w(i_index) 
            else
              temp = temp + &
                previous_eigenvector_complex(i_basis_2,&
                force_occ_pr_state_periodic(i_force_occ,i_k_point),&
                force_occ_spin(i_force_occ),i_k) * &
                conjg(KS_eigenvector_complex(i_basis_1,i_state,force_occ_spin(i_force_occ),i_k)) * &
                conjg(overlap_matrix_w_complex(i_index))
            end if
          end do
        end do
        i_index = 0
        do i_basis_2 = 1, n_basis, 1
          do i_basis_1 = 1, i_basis_2-1, 1
            i_index = i_index + 1
            if (real_eigenvectors) then
              temp = temp + &
                previous_eigenvector(i_basis_1,&
                force_occ_pr_state_periodic(i_force_occ,i_k_point),&
                force_occ_spin(i_force_occ),i_k) * &
                KS_eigenvector(i_basis_2,i_state,force_occ_spin(i_force_occ),i_k) * &
                overlap_matrix_w(i_index)
            else
              temp = temp + &
                previous_eigenvector_complex(i_basis_1,&
                force_occ_pr_state_periodic(i_force_occ,i_k_point),&
                force_occ_spin(i_force_occ),i_k) * &
                conjg(KS_eigenvector_complex(i_basis_2,i_state,force_occ_spin(i_force_occ),i_k)) * &
                conjg(overlap_matrix_w_complex(i_index)) 
            end if
          end do
          i_index = i_index + 1
        end do
!         temp = temp + conjg(temp)
!     do i_basis_1 = 1, n_basis, 1
!           i_index = i_basis_1*(1+i_basis_1)/2
!           temp = temp + &
!         previous_eigenvector_complex(i_basis_1,&
!             force_occ_pr_state_periodic(i_force_occ,i_k_point),&
!         force_occ_spin(i_force_occ),i_k) * &
!         conjg(KS_eigenvector(i_basis_1,i_state,force_occ_spin(i_force_occ),i_k)) * &
!         overlap_matrix_w_complex(i_index) 
!         end do
        if (abs(temp).gt.biggest) then
!           write(use_unit,*)"IN PR-ROUTINE",abs(temp),biggest,i_state,force_occ_pr_state
          force_occ_state_periodic(i_force_occ,i_k_point) = i_state
          biggest = abs(temp)
        end if
      end do

      if (myid.eq.0) then
        write(use_unit,*)"BIG",biggest
        write(use_unit,*)"FORCED TO: now -- previous" ,force_occ_state_periodic(i_force_occ,:), &
                                                       force_occ_pr_state_periodic(i_force_occ,:)
      end if

    end do

    
  end subroutine determine_corehole_projection_periodic

!******
!-------------------------------------------------------------------------------
!****s* force_occupation/zeroin_occ_p0
!  NAME
!    zeroin_occ_p0
!  SYNOPSIS
!   VB 2005:
!   Routine taken from http://www.netlib.org/go/
!
!   Copyright: My impression is that these routines come with a liberal (if nonexistant)
!   copyright provision, but netlib states that if in doubt, ask the original
!   author. In any case, it is publicly available for anyone to use.
!
!   Adapted by Ralf Gehrke 2005 to work for our purposes
!
!   Description from netlib
!   lib      ../go
!   for      Golden Oldies:  widely used,  but not in standard libraries.
!   #      Nominations welcome!
!   rel      excellent
!   age      old
!   editor      Eric Grosse
!   master      netlib.bell-labs.com
!
!  To get d1mach, mail netlib
!       send d1mach from core
!  VB: In my opinion, we already have the machine accuracy in question, but we can always get it
!      from netlib by downloading the dependencies also.
!
      subroutine zeroin_occ_p0(ax, bx, fax, fbx, tol, &
           occupation_acc, max_zeroin, &
           n_states, KS_eigenvalue, &
           n_electrons, chemical_potential, diff_electrons, &
           occ_numbers, i_counter)



      implicit none


! imported variables

! input
      real*8, intent(in) :: ax
      real*8, intent(in) :: bx
      real*8, intent(in) :: fax
      real*8, intent(in) :: fbx
      real*8, intent(in) :: tol

      real*8, intent(in) :: occupation_acc
      integer, intent(in) :: max_zeroin

      integer, intent(in) :: n_states
!      integer, intent(in) :: n_spin
      real*8, dimension(n_states), intent(in) :: &
           KS_eigenvalue
      real*8, intent(in) :: n_electrons

! output
      real*8, dimension(n_states, n_spin, n_k_points) :: &
           occ_numbers
!      real*8, dimension(n_states,n_spin,n_k_points), intent(out) ::
!     +     occ_numbers
!      real*8, dimension(:,:,:), intent(out) ::
!     +     occ_numbers
      real*8, intent(out) :: chemical_potential
      real*8, intent(out) :: diff_electrons
      integer, intent(inout) :: i_counter



!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
      real*8  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, &
           q, r, s
      real*8  dabs, d1mach

 10   eps = d1mach(4)
      tol1 = eps+1.0d0
!
      a  = ax
      b  = bx
      fa = fax
      fb = fbx

!test
!      write(use_unit,*) "After 10: "
!      write(use_unit,*) " a = ", a
!      write(use_unit,*) " b = ", b
!      write(use_unit,*) "fa = ", fa
!      write(use_unit,*) "fb = ", fb
!test end

!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(use_unit,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' zeroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if (abs(fb) .le. occupation_acc) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge. &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
!  140 fb=f(b)
 140  continue

      if (force_occupation_smearing) then
        call check_norm_force_occ_p0_smearing(b, KS_eigenvalue, n_electrons, &
           occ_numbers, fb, i_counter)
      else
        call check_norm_force_occ_p0(b, KS_eigenvalue, n_electrons, &
           occ_numbers, fb, i_counter)
      end if
!test
!      write(use_unit,*) "i_counter = ", i_counter
!      write(use_unit,*) "n_electrons: ", n_electrons
!      write(use_unit,*) "occ_numbers: ", occ_numbers
!test end

      if (i_counter .gt. max_zeroin) then
!       force end despite non-convergence
        go to 150
      end if
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 chemical_potential = b
      diff_electrons     = fb
      return

      end subroutine zeroin_occ_p0

!******
!-------------------------------------------------------------------------------
!****s* force_occupation/zeroin_occ_v2
!  NAME
!    zeroin_occ_v2
!  SYNOPSIS

!   VB 2005:
!   Routine taken from http://www.netlib.org/go/
!
!   Copyright: My impression is that these routines come with a liberal (if nonexistant)
!   copyright provision, but netlib states that if in doubt, ask the original
!   author. In any case, it is publicly available for anyone to use.
!
!   Adapted by Ralf Gehrke 2005 to work for our purposes
!
!   Description from netlib
!   lib      ../go
!   for      Golden Oldies:  widely used,  but not in standard libraries.
!   #      Nominations welcome!
!   rel      excellent
!   age      old
!   editor      Eric Grosse
!   master      netlib.bell-labs.com
!
!  To get d1mach, mail netlib
!       send d1mach from core
!  VB: In my opinion, we already have the machine accuracy in question, but we can always get it
!      from netlib by downloading the dependencies also.
!
      subroutine zeroin_occ_v2(ax, bx, fax, fbx, tol, &
           occupation_acc, max_zeroin, &
           n_states, KS_eigenvalue, &
           n_electrons, chemical_potential, diff_electrons, &
           occ_numbers, i_counter,current_spin)


      implicit none

! imported variables

! input
      real*8, intent(in) :: ax
      real*8, intent(in) :: bx
      real*8, intent(in) :: fax
      real*8, intent(in) :: fbx
      real*8, intent(in) :: tol

      real*8, intent(in) :: occupation_acc
      integer, intent(in) :: max_zeroin

      integer, intent(in) :: n_states
!      integer, intent(in) :: n_spin
      real*8, dimension(n_states), intent(in) :: KS_eigenvalue
      real*8, intent(in) :: n_electrons

! output
      real*8, dimension(n_states), intent(out) :: occ_numbers
      real*8, intent(out) :: chemical_potential
      real*8, intent(out) :: diff_electrons
      integer, intent(inout) :: i_counter

      integer :: current_spin



!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
      real*8  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, &
           q, r, s
      real*8  dabs, d1mach

 10   eps = d1mach(4)
      tol1 = eps+1.0d0
!
      a  = ax
      b  = bx
      fa = fax
      fb = fbx

!test
!      write(use_unit,*) "After 10: "
!      write(use_unit,*) " a = ", a
!      write(use_unit,*) " b = ", b
!      write(use_unit,*) "fa = ", fa
!      write(use_unit,*) "fb = ", fb
!test end

!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(use_unit,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' zeroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if (abs(fb) .le. occupation_acc) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge. &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
!  140 fb=f(b)
 140  continue
      call check_norm_force_occ_v2(b, KS_eigenvalue, n_electrons, &
           occ_numbers, fb, i_counter,current_spin)
!test
!      write(use_unit,*) "i_counter = ", i_counter
!      write(use_unit,*) "n_electrons: ", n_electrons
!      write(use_unit,*) "occ_numbers: ", occ_numbers
!test end

      if (i_counter .gt. max_zeroin) then
!       force end despite non-convergence
        go to 150
      end if
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 chemical_potential = b
      diff_electrons     = fb
      return

      end subroutine zeroin_occ_v2



!******
!-------------------------------------------------------------------------------
!****s* force_occupation/zeroin_occ_periodic_v2
!  NAME
!    zeroin_occ_periodic_v2
!  SYNOPSIS

!   VB 2005:
!   Routine taken from http://www.netlib.org/go/
!
!   Copyright: My impression is that these routines come with a liberal (if nonexistant)
!   copyright provision, but netlib states that if in doubt, ask the original
!   author. In any case, it is publicly available for anyone to use.
!
!   Adapted by Ralf Gehrke 2005 to work for our purposes
!
!   Description from netlib
!   lib      ../go
!   for      Golden Oldies:  widely used,  but not in standard libraries.
!   #      Nominations welcome!
!   rel      excellent
!   age      old
!   editor      Eric Grosse
!   master      netlib.bell-labs.com
!
!  To get d1mach, mail netlib
!       send d1mach from core
!  VB: In my opinion, we already have the machine accuracy in question, but we can always get it
!      from netlib by downloading the dependencies also.
!
      subroutine zeroin_occ_periodic_v2(ax, bx, fax, fbx, tol, &
           occupation_acc, max_zeroin, &
           n_states, KS_eigenvalue, &
           n_electrons, chemical_potential, diff_electrons, &
           occ_numbers, i_counter,current_spin)


      implicit none

! imported variables

! input
      real*8, intent(in) :: ax
      real*8, intent(in) :: bx
      real*8, intent(in) :: fax
      real*8, intent(in) :: fbx
      real*8, intent(in) :: tol

      real*8, intent(in) :: occupation_acc
      integer, intent(in) :: max_zeroin

      integer, intent(in) :: n_states
!      integer, intent(in) :: n_spin
      real*8, dimension(n_states,n_k_points), intent(in) :: KS_eigenvalue
      real*8, intent(in) :: n_electrons

! output
      real*8, dimension(n_states,n_k_points), intent(out) :: occ_numbers
      real*8, intent(out) :: chemical_potential
      real*8, intent(out) :: diff_electrons
      integer, intent(inout) :: i_counter

      integer :: current_spin



!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
      real*8  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, &
           q, r, s
      real*8  dabs, d1mach

 10   eps = d1mach(4)
      tol1 = eps+1.0d0
!
      a  = ax
      b  = bx
      fa = fax
      fb = fbx

!test
!      write(use_unit,*) "After 10: "
!      write(use_unit,*) " a = ", a
!      write(use_unit,*) " b = ", b
!      write(use_unit,*) "fa = ", fa
!      write(use_unit,*) "fb = ", fb
!test end

!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(use_unit,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' zeroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if (abs(fb) .le. occupation_acc) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge. &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
!  140 fb=f(b)
 140  continue
      call check_norm_force_occ_periodic_v2(b, KS_eigenvalue, n_electrons, &
           occ_numbers, fb, i_counter,current_spin)
!test
!      write(use_unit,*) "i_counter = ", i_counter
!      write(use_unit,*) "n_electrons: ", n_electrons
!      write(use_unit,*) "occ_numbers: ", occ_numbers
!test end

      if (i_counter .gt. max_zeroin) then
!       force end despite non-convergence
        go to 150
      end if
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 chemical_potential = b
      diff_electrons     = fb
      return

      end subroutine zeroin_occ_periodic_v2








!******
!-------------------------------------------------------------------------------
!****s* force_occupation/zeroin_constraint_occ
!  NAME
!    zeroin_constraint_occ
!  SYNOPSIS

!   VB 2005:
!   Routine taken from http://www.netlib.org/go/
!
!   Copyright: My impression is that these routines come with a liberal (if nonexistant)
!   copyright provision, but netlib states that if in doubt, ask the original
!   author. In any case, it is publicly available for anyone to use.
!
!   Adapted by Ralf Gehrke 2005 to work for our purposes
!
!   Description from netlib
!   lib      ../go
!   for      Golden Oldies:  widely used,  but not in standard libraries.
!   #      Nominations welcome!
!   rel      excellent
!   age      old
!   editor      Eric Grosse
!   master      netlib.bell-labs.com
!
!  To get d1mach, mail netlib
!       send d1mach from core
!  VB: In my opinion, we already have the machine accuracy in question, but we can always get it
!      from netlib by downloading the dependencies also.
!
      subroutine zeroin_constraint_occ(ax, bx, fax, fbx, tol, &
           occupation_acc, max_zeroin, &
           n_states, KS_eigenvalue, &
           constraint_electrons,chemical_potential, diff_electrons, &
           constraint_proj, occ_numbers, i_counter,&
           current_spin)

      implicit none

! imported variables

! input
      real*8, intent(in) :: ax
      real*8, intent(in) :: bx
      real*8, intent(in) :: fax
      real*8, intent(in) :: fbx
      real*8, intent(in) :: tol

      real*8, intent(in) :: occupation_acc
      integer, intent(in) :: max_zeroin

      integer, intent(in) :: n_states
      real*8, dimension(n_states), intent(in) :: KS_eigenvalue
      real*8, dimension(n_states), intent(in) :: constraint_proj
      real*8, intent(in) :: constraint_electrons

! output
      real*8, dimension(n_states), intent(out) :: occ_numbers
      real*8, intent(out) :: chemical_potential
      real*8, intent(out) :: diff_electrons
      integer, intent(inout) :: i_counter

      integer :: current_spin


!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
      double precision  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, &
           q, r, s
      double precision  dabs, d1mach

 10   eps = d1mach(4)
      tol1 = eps+1.0d0
!
      a  = ax
      b  = bx
      fa = fax
      fb = fbx

!test
!      write(use_unit,*) "After 10: "
!      write(use_unit,*) " a = ", a
!      write(use_unit,*) " b = ", b
!      write(use_unit,*) "fa = ", fa
!      write(use_unit,*) "fb = ", fb
!test end

!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(use_unit,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' zeroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if (abs(fb) .le. occupation_acc) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge. &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
!  140 fb=f(b)
 140  continue
      call check_norm_constraint_occ(b, KS_eigenvalue, constraint_electrons, &
           constraint_proj, occ_numbers, fb, i_counter,current_spin)
!test
!      write(use_unit,*) "i_counter = ", i_counter
!      write(use_unit,*) "constraint_electrons: ", constraint_electrons
!      write(use_unit,*) "occ_numbers: ", occ_numbers
!test end

      if (i_counter .gt. max_zeroin) then
!       force end despite non-convergence
        go to 150
      end if
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 chemical_potential = b
      diff_electrons     = fb
      return

      end subroutine zeroin_constraint_occ

!******
!-------------------------------------------------------------------------------
!****s* force_occupation/get_constraint_fermi_occ
!  NAME
!    get_constraint_fermi_occ
!  SYNOPSIS
!
! determines occupation numbers according to two available schemes:
! 1) gaussian smearing according to Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
! 2) fermi smearing
! 3) methfessel-paxton (methfessel et al., PRB 40, 6, 40 (1989))
! R. Gehrke (2005)
!---------------------------------------------------------------------- 

subroutine get_constraint_fermi_occ(KS_eigenvalue,constraint_electrons,constraint_proj,constraint_fermi,&
    current_spin)
  
!   use dimensions
!   use runtime_choices
  use constants
  use mpi_tasks, only: aims_stop
  use localorb_io


  implicit none
  
  
  !  imported variables

  !  input 
  
  real*8, dimension(n_states), intent(in) :: KS_eigenvalue
  real*8, dimension(n_states), intent(in) :: constraint_proj
  real*8, intent(in) :: constraint_electrons

  !  output
  
  real*8, intent(out) :: constraint_fermi
  
  integer :: current_spin
  integer,dimension(n_force_occ) :: force_occ_states

!  local variables
  
  real*8, dimension(n_states) :: occ_numbers
  real*8 :: chemical_potential

  real*8 :: degeneracy_threshold
  integer :: n_degenerate
  integer :: highest_full_state
  integer :: lowest_empty_state
  real*8 :: shared_electrons
  logical :: degenerate
  real*8  :: chemical_potential_l
  real*8  :: chemical_potential_r
  real*8 :: diff_l
  real*8 :: diff_r
  real*8 :: diff_electrons
  real*8 :: diff_electrons_thr
  real*8 :: lowest_eigenvalue
  real*8 :: highest_eigenvalue
  
  real*8 :: electron_count
  
  character*120 :: info_str

  !  counters
  
  integer :: i_state
  integer :: i_state_2
!  integer :: i_spin
  integer :: i_counter
  
!  write(use_unit,'(2X,A)') "Determining occupation numbers for Kohn-Sham eigenstates."

  i_counter = 0

!  write(use_unit,*) "get_occ_numbers.... n_spin=", n_spin, "n_states=", n_states
!  write(use_unit,*)"CONSTRAINT ELECTRONS: ",constraint_electrons

  lowest_eigenvalue = KS_eigenvalue(1)
  do i_state = 1, n_states, 1
     if (KS_eigenvalue(i_state) .lt. lowest_eigenvalue) then
        lowest_eigenvalue = KS_eigenvalue(i_state)
     end if
  end do

  highest_eigenvalue = KS_eigenvalue(n_states)
  do i_state = 1, n_states, 1
     if (KS_eigenvalue(i_state) .gt. highest_eigenvalue) then
        highest_eigenvalue = KS_eigenvalue(i_state)
     end if
  end do
 



!  initialize zeroin algorithm
!  to avoid problems if only one state is calculated (so for the H-atom)
!  do not initialize mit KS_eigenvalue(n_states), because it is then
!  identical to KS_eigenvalue(1)

  chemical_potential_r = highest_eigenvalue
  chemical_potential_l = lowest_eigenvalue
  call check_norm_constraint_occ(chemical_potential_l, KS_eigenvalue, constraint_electrons, &
                             constraint_proj, occ_numbers, diff_l, i_counter,current_spin)
  call check_norm_constraint_occ(chemical_potential_r, KS_eigenvalue, constraint_electrons, & 
                             constraint_proj, occ_numbers, diff_r, i_counter,current_spin)
  do while (diff_l * diff_r .gt. 0.d0)

!     interval for chemical potential still not found

     chemical_potential_l = chemical_potential_r
     chemical_potential_r = chemical_potential_r + 0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
     diff_l = diff_r
     call check_norm_constraint_occ(chemical_potential_r, KS_eigenvalue, constraint_electrons, &
                                constraint_proj, occ_numbers, diff_r, i_counter,current_spin)

!     write(use_unit,*) chemical_potential_l, chemical_potential_r
     !     If the right interval cannot be found, stop potentially
     !     infinite loop
     
     if (i_counter .gt. max_zeroin) then
        write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* Too many iterations needed to find good interval around E_F."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* According to this check, you may not be able to find one."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "Please look carefully at your settings."
        call localorb_info(info_str,use_unit,'(A)')
        call aims_stop("*** Error in finding electronic chemical potential", &
              "force_occupation")
     end if
  end do
      
  !test
  !     write(use_unit,*) "Before zeroin."
  !     write(use_unit,*) "left: ", chemical_potential_l
  !         write(use_unit,*) "right: ", chemical_potential_r
  !         write(use_unit,*) "occ_numbers : ", occ_numbers
  !test end 
  
  call zeroin_constraint_occ(chemical_potential_l, chemical_potential_r, diff_l, &
       diff_r, fermi_acc, occupation_acc, max_zeroin, &
       n_states, KS_eigenvalue, constraint_electrons, chemical_potential, &
       diff_electrons, constraint_proj, occ_numbers, i_counter,current_spin)
      
  !     call check_norm one more time, just to make sure we really have the 
  !     occupation numbers for the final fermi level [this is critical for H]

  call check_norm_constraint_occ( chemical_potential, KS_eigenvalue, constraint_electrons, &
                              constraint_proj, occ_numbers, diff_electrons, i_counter,current_spin)

  !     and decrement i_counter which was falsely incremented ...

  i_counter = i_counter - 1
         
  !       If zeroin did not converge to the desired accuracy, issue
  !       warning and stop, for now; we can create nicer behavior after
  !       we gain more experience.

  if (i_counter .gt. max_zeroin) then
     write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
              call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') "* Too many iterations needed to find correct occupation of levels."
              call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') &
       "* Check variables occupation_acc, max_zeroin before continuing your calculation."
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* Status of get_occ_numbers when aborted: "
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* State #         Eigenvalue       Occ. number"
              call localorb_info(info_str,use_unit,'(A)')

!        write(use_unit,*) "* spin # ", i_spin
        write (info_str,*) "----------"
        do i_state = 1, n_states, 1
           write (info_str,*) i_state, KS_eigenvalue(i_state), occ_numbers(i_state)
        end do
      
     stop
  end if
  
  !     finally, occupation thresholding - if needed

  if (occupation_thr.gt.0.d0) then
     call threshold_occ_numbers_constraint( constraint_electrons, occ_numbers, diff_electrons_thr)
  end if
      
  constraint_fermi = chemical_potential

!  write(use_unit,*)"PROJ_OP, OCC_NUMBERS, OCC_NUMBERS*PROJ_OP, KS_EIGENVALUES"
!  do i_state = 1, n_states, 1
!     write(use_unit,*)"STATE:", i_state
!     write(use_unit,*) constraint_proj(i_state),occ_numbers(i_state),occ_numbers(i_state)*constraint_proj(i_state),KS_eigenvalue(i_state)
!  end do
!  write(use_unit,*)"CONSTRAINT FERMI:", constraint_fermi
!  write(use_unit, '(2X, A, E10.4)') "| Chemical potential (Fermi level) in eV                 : ", chemical_potential * hartree
!  write(use_unit, '(2X, A, E10.4)') "| Error in electron count due to remaining E_F inaccuracy: ", diff_electrons
!  if (occupation_thr.gt.0.d0) then
!     write(use_unit, '(2X, A, E10.4)') "| Error in electron count after thresholding : ", diff_electrons_thr
!  end if
!  write(use_unit,*)
  
end subroutine get_constraint_fermi_occ
      


!******
!-------------------------------------------------------------------------------
!****s* force_occupation/get_occupation_numbers_occ_p0
!  NAME
!    get_occupation_numbers_occ_p0
!  SYNOPSIS

  subroutine get_occupation_numbers_occ_p0(KS_eigenvalue, n_electrons, t_out, occ_numbers,&
      chemical_potential)
  
  !  PURPOSE
  ! Subroutine get_occupation_numbers
  !
  ! determines occupation numbers according to two available schemes:
  ! o 1) gaussian smearing according to Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
  ! o 2) fermi smearing
  ! o 3) methfessel-paxton (methfessel et al., PRB 40, 6, 40 (1989))
  !
  !  USES
  
  !   use dimensions
  !   use runtime_choices
    use localorb_io
    use constants
    implicit none
  
  !  ARGUMENTS
  
    real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: KS_eigenvalue
    real*8, intent(in) :: n_electrons
  
  
    logical :: t_out
  
    real*8, dimension(n_states, n_spin, n_k_points), intent(out) :: occ_numbers
    real*8, intent(out) :: chemical_potential
  
  !  INPUTS
  !  o KS_eigenvalue -- Kohn-Sham eigenvalues
  !  o n_electrons -- number of electrons
  !  o t_out -- is the information printed out or not ?
  !
  !  OUTPUT
  !  o occ_numbers -- occupation weights of different KS states
  !  o chemical_potential -- chemical potential
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
  
  
  
  
    !  local variables
  !  real*8 :: degeneracy_threshold
  !  integer :: n_degenerate
  !  integer :: highest_full_state
  !  integer :: lowest_empty_state
  !  real*8 :: shared_electrons
  !  logical :: degenerate
    real*8 :: chemical_potential_l
    real*8 :: chemical_potential_r
    real*8 :: diff_l
    real*8 :: diff_r
    real*8 :: diff_electrons
    real*8 :: diff_electrons_thr
    real*8 :: lowest_eigenvalue
    real*8 :: highest_eigenvalue
  
  !  real*8 :: electron_count
  
    character*120 :: info_str
  
    !  counters
    integer :: i_state
  !  integer :: i_state_2
    integer :: i_spin
    integer :: i_k_points
    integer :: i_counter
  
  
    if (t_out) then
      write(info_str,'(2X,A)') &
        "Determining occupation numbers for Kohn-Sham eigenstates."
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
  
    i_counter = 0
    !  write(use_unit,*) "get_occ_numbers.... n_spin=", n_spin, "n_states=", n_states
  
    lowest_eigenvalue = KS_eigenvalue(1,1,1)
  
    do i_k_points = 1, n_k_points,1
      do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1
          if (KS_eigenvalue(i_state, i_spin,i_k_points) .lt. lowest_eigenvalue) then
            lowest_eigenvalue = KS_eigenvalue(i_state, i_spin, i_k_points)
          end if
        end do
      end do
    end do
  
    highest_eigenvalue = KS_eigenvalue(n_states,1,1)
    do i_k_points = 1,n_k_points,1
      do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1
          if (KS_eigenvalue(i_state, i_spin, i_k_points) .gt. highest_eigenvalue) then
            highest_eigenvalue = KS_eigenvalue(i_state, i_spin, i_k_points)
          end if
        end do
      end do
    end do
  
    !  initialize zeroin algorithm
    !  to avoid problems if only one state is calculated (so for the H-atom)
    !  do not initialize mit KS_eigenvalue(n_states), because it is then
    !  identical to KS_eigenvalue(1)
  
    if (lowest_eigenvalue.ne.highest_eigenvalue) then
      chemical_potential_r = highest_eigenvalue
    else
      chemical_potential_r = 0.d0
    end if
    chemical_potential_l = lowest_eigenvalue
 
    if (force_occupation_smearing.and.n_periodic.eq.0) then
        call check_norm_force_occ_p0_smearing(chemical_potential_l, KS_eigenvalue, n_electrons,&
             occ_numbers, diff_l, i_counter)
        call check_norm_force_occ_p0_smearing(chemical_potential_r, KS_eigenvalue, n_electrons,&
             occ_numbers, diff_r, i_counter)
    else
        call check_norm_force_occ_p0(chemical_potential_l, KS_eigenvalue, n_electrons,&
             occ_numbers, diff_l, i_counter)
        call check_norm_force_occ_p0(chemical_potential_r, KS_eigenvalue, n_electrons,&
             occ_numbers, diff_r, i_counter)
    end if  

    do while (diff_l * diff_r .gt. 0.d0)
      !     interval for chemical potential still not found
      !     Must extend intervals both ways: 
      !     * For zero electrons, need Fermi level below lowest state
      !     * For all electrons in one channel, need Fermi level above highest state
  
      chemical_potential_l = chemical_potential_l - 0.5d0 * dabs(highest_eigenvalue-lowest_eigenvalue)
      chemical_potential_r = chemical_potential_r + 0.5d0 * dabs(highest_eigenvalue-lowest_eigenvalue)
      diff_l = diff_r

      if (force_occupation_smearing.and.n_periodic.eq.0) then
          call check_norm_force_occ_p0_smearing(chemical_potential_l, KS_eigenvalue, n_electrons, &
               occ_numbers, diff_l, i_counter)
          call check_norm_force_occ_p0_smearing(chemical_potential_r, KS_eigenvalue, n_electrons, &
               occ_numbers, diff_r, i_counter)
      else
          call check_norm_force_occ_p0(chemical_potential_l, KS_eigenvalue, n_electrons, &
               occ_numbers, diff_l, i_counter)
          call check_norm_force_occ_p0(chemical_potential_r, KS_eigenvalue, n_electrons, &
               occ_numbers, diff_r, i_counter)
      end if

      !     If the right interval cannot be found, stop potentially
      !     infinite loop
  
      if (i_counter .gt. max_zeroin) then
        write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* Too many iterations needed to find good interval around E_F."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
          "* According to this check, you may not be able to find one. "
         call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
          "Please look carefully at your settings."
        call localorb_info(info_str,use_unit,'(A)')
        stop
      end if
    end do
  
    call zeroin_occ_p0(chemical_potential_l, chemical_potential_r, diff_l, diff_r, fermi_acc,&
    occupation_acc, max_zeroin, n_states, KS_eigenvalue, n_electrons, chemical_potential, &
    diff_electrons, occ_numbers, i_counter)
  
  
    !     call check_norm one more time, just to make sure we really have the 
    !     occupation numbers for the final fermi level [this is critical for H]

    if (force_occupation_smearing) then
        call check_norm_force_occ_p0_smearing( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, &
            diff_electrons, i_counter)
    else
        call check_norm_force_occ_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, &
             diff_electrons, i_counter)
    end if
    !     and decrement i_counter which was falsely incremented ...
    i_counter = i_counter - 1
  
    !       If zeroin did not converge to the desired accuracy, issue
    !       warning and stop, for now; we can create nicer behavior after
    !       we gain more experience.
    if (i_counter .gt. max_zeroin) then
      write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,'(1X,A)') "* Too many iterations needed to find correct occupation of levels."
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,'(1X,A)') "* Check variables occupation_acc, max_zeroin before continuing your calculation."
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* Status of get_occ_numbers when aborted: "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* State #         Eigenvalue       Occ. number"
      call localorb_info(info_str,use_unit,'(A)')
      do i_k_points = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
          write (info_str,*) "* spin # ", i_spin
          call localorb_info(info_str,use_unit,'(A)')
          write (info_str,*) "----------"
          call localorb_info(info_str,use_unit,'(A)')
          do i_state = 1, n_states, 1
            write (info_str,*) i_state, KS_eigenvalue(i_state, i_spin,i_k_points), occ_numbers(i_state, i_spin,i_k_points)
            call localorb_info(info_str,use_unit,'(A)')
          end do
        end do
      end do
      stop
    end if
  
    !     finally, occupation thresholding - if needed
    if (occupation_thr.gt.0.d0) then
      call threshold_occ_numbers( n_electrons, occ_numbers, diff_electrons_thr)
    end if
  
    if (t_out) then
      write (info_str, '(2X, A, E10.4)') "| Chemical potential (Fermi level) in eV                 : ", chemical_potential * hartree
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      write (info_str, '(2X, A, E10.4)') "| Error in electron count due to remaining E_F inaccuracy: ", diff_electrons
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      if (occupation_thr.gt.0.d0) then
        write (info_str, '(2X, A, E10.4)') "| Error in electron count after thresholding : ", diff_electrons_thr
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
      end if
      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
  
  
  end subroutine get_occupation_numbers_occ_p0


!******
!-------------------------------------------------------------------------------
!****s* force_occupation/get_occupation_numbers_occ_v2
!  NAME
!   get_occupation_numbers_occ_v2
!  SYNOPSIS

  subroutine get_occupation_numbers_occ_v2(KS_eigenvalue, n_electrons, t_out, occ_numbers,&
             chemical_potential,current_spin)
  
  !  PURPOSE
  ! Subroutine get_occupation_numbers determines occupation numbers according to two available schemes:
  ! o 1) gaussian smearing according to Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
  ! o 2) fermi smearing
  ! o 3) methfessel-paxton (methfessel et al., PRB 40, 6, 40 (1989))
  !
  !  USES
  
    use dimensions
    use runtime_choices
    use localorb_io
    use constants
    implicit none
  
  !  ARGUMENTS
  
    real*8, dimension(n_states), intent(in) :: KS_eigenvalue
    real*8, intent(in) :: n_electrons
    logical :: t_out
  
    real*8, dimension(n_states), intent(out) :: occ_numbers
    real*8, intent(out) :: chemical_potential
  
    integer :: current_spin

  
  !  INPUTS
  !  o KS_eigenvalue -- Kohn-Sham eigenvalues
  !  o n_electrons -- number of electrons
  !  o t_out -- is the information printed out or not ?
  !
  !  OUTPUT
  !  o occ_numbers -- occupation weights of different KS states
  !  o chemical_potential -- chemical potential
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


    !  local variables
    real*8 :: chemical_potential_l
    real*8 :: chemical_potential_r
    real*8 :: diff_l
    real*8 :: diff_r
    real*8 :: diff_electrons
    real*8 :: diff_electrons_thr
    real*8 :: lowest_eigenvalue
    real*8 :: highest_eigenvalue
  
  
    character*120 :: info_str
  
    !  counters
    integer :: i_state
    integer :: i_counter
  
  
    if (t_out) then
      write(info_str,'(2X,A)') &
         "Determining occupation numbers for Kohn-Sham eigenstates."
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if
  
    i_counter = 0
    !  write(use_unit,*) "get_occ_numbers.... n_spin=", n_spin, "n_states=", n_states
  
    lowest_eigenvalue = KS_eigenvalue(1)
    do i_state = 1, n_states, 1
      if (KS_eigenvalue(i_state) .lt. lowest_eigenvalue) then
        lowest_eigenvalue = KS_eigenvalue(i_state)
      end if
    end do
  
    highest_eigenvalue = KS_eigenvalue(n_states)
    do i_state = 1, n_states, 1
      if (KS_eigenvalue(i_state) .gt. highest_eigenvalue) then
        highest_eigenvalue = KS_eigenvalue(i_state)
      end if
    end do
  
    !  initialize zeroin algorithm
    !  to avoid problems if only one state is calculated (so for the H-atom)
    !  do not initialize mit KS_eigenvalue(n_states), because it is then
    !  identical to KS_eigenvalue(1)
  
    if (lowest_eigenvalue.ne.highest_eigenvalue) then
      chemical_potential_r = highest_eigenvalue
    else
      chemical_potential_r = 0.d0
    end if
    chemical_potential_l = lowest_eigenvalue
  
    call check_norm_force_occ_v2(chemical_potential_l, KS_eigenvalue, n_electrons, occ_numbers, diff_l,&
      i_counter,current_spin)
    call check_norm_force_occ_v2(chemical_potential_r, KS_eigenvalue, n_electrons, occ_numbers, diff_r,&
      i_counter,current_spin)
    do while (diff_l * diff_r .gt. 0.d0)
      !     interval for chemical potential still not found
      !     Must extend intervals both ways: 
      !     * For zero electrons, need Fermi level below lowest state
      !     * For all electrons in one channel, need Fermi level above highest state
      chemical_potential_l = chemical_potential_l - 0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
      chemical_potential_r = chemical_potential_r + 0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
      diff_l = diff_r
      call check_norm_force_occ_v2(chemical_potential_l, KS_eigenvalue, n_electrons, occ_numbers, diff_l,&
        i_counter,current_spin)
      call check_norm_force_occ_v2(chemical_potential_r, KS_eigenvalue, n_electrons, occ_numbers, diff_r,&
        i_counter,current_spin)
      !     write(use_unit,*) chemical_potential_l, chemical_potential_r
      !     If the right interval cannot be found, stop potentially
      !     infinite loop
      if (i_counter .gt. max_zeroin) then
          write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
          call localorb_info(info_str,use_unit,'(A)')
          write (info_str,'(1X,A)') "* Too many iterations needed to find good interval around E_F."
          call localorb_info(info_str,use_unit,'(A)')
          write (info_str,'(1X,A)') &
              "* According to this check, you may not be able to find one."
          call localorb_info(info_str,use_unit,'(A)')
          write (info_str,'(1X,A)') &
              "Please look carefully at your settings."
          call localorb_info(info_str,use_unit,'(A)')
          stop
      end if
    end do
  
    call zeroin_occ_v2(chemical_potential_l, chemical_potential_r, diff_l, diff_r,&
         fermi_acc, occupation_acc, max_zeroin,n_states, KS_eigenvalue, n_electrons, &
         chemical_potential, diff_electrons, occ_numbers, i_counter,current_spin)

    !     call check_norm one more time, just to make sure we really have the 
    !     occupation numbers for the final fermi level [this is critical for H]
    call check_norm_force_occ_v2( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers,&
         diff_electrons, i_counter,current_spin)
    !     and decrement i_counter which was falsely incremented ...
    i_counter = i_counter - 1
  
    !       If zeroin did not converge to the desired accuracy, issue
    !       warning and stop, for now; we can create nicer behavior after
    !       we gain more experience.
    if (i_counter .gt. max_zeroin) then
      write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,'(1X,A)') "* Too many iterations needed to find correct occupation of levels."
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,'(1X,A)') "* Check variables occupation_acc, max_zeroin before continuing your calculation."
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* Status of get_occ_numbers when aborted: "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* State #         Eigenvalue       Occ. number"
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,*) "----------"
      call localorb_info(info_str,use_unit,'(A)')
      do i_state = 1, n_states, 1
        write (info_str,*) i_state, KS_eigenvalue(i_state), occ_numbers(i_state)
        call localorb_info(info_str,use_unit,'(A)')
      end do
      stop
    end if
  
    !     finally, occupation thresholding - if needed
    if (occupation_thr.gt.0.d0) then
      call threshold_occ_numbers( n_electrons, occ_numbers, diff_electrons_thr)
    end if

    if (t_out) then
      write (info_str, '(2X, A, E10.4)') "| Chemical potential (Fermi level) in eV                 : ", &
        chemical_potential * hartree
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str, '(2X, A, E10.4)') "| Error in electron count due to remaining E_F inaccuracy: ", diff_electrons
      call localorb_info(info_str,use_unit,'(A)')
      if (occupation_thr.gt.0.d0) then
        write (info_str, '(2X, A, E10.4)') "| Error in electron count after thresholding : ", diff_electrons_thr
        call localorb_info(info_str,use_unit,'(A)')
      end if
      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)')
    end if
  
  
  end subroutine get_occupation_numbers_occ_v2

!******
!-------------------------------------------------------------------------------
!****s* force_occupation/get_occupation_numbers_occ_periodic_v2
!  NAME
!   get_occupation_numbers_occ_periodic_v2
!  SYNOPSIS

  subroutine get_occupation_numbers_occ_periodic_v2(KS_eigenvalue, n_electrons, t_out, occ_numbers,&
             chemical_potential,current_spin)
  
  !  PURPOSE
  ! Subroutine get_occupation_numbers determines occupation numbers according to two available schemes:
  ! o 1) gaussian smearing according to Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
  ! o 2) fermi smearing
  ! o 3) methfessel-paxton (methfessel et al., PRB 40, 6, 40 (1989))
  !
  ! needed for fixed-spin moment calculations in periodic case
  !
  !  USES
  
    use dimensions
    use runtime_choices
    use localorb_io
    use constants
    implicit none
  
  !  ARGUMENTS
  
    real*8, dimension(n_states,n_k_points), intent(in) :: KS_eigenvalue
    real*8, intent(in) :: n_electrons
    logical :: t_out
  
    real*8, dimension(n_states,n_k_points), intent(out) :: occ_numbers
    real*8, intent(out) :: chemical_potential
  
    integer :: current_spin

  
  !  INPUTS
  !  o KS_eigenvalue -- Kohn-Sham eigenvalues
  !  o n_electrons -- number of electrons
  !  o t_out -- is the information printed out or not ?
  !
  !  OUTPUT
  !  o occ_numbers -- occupation weights of different KS states
  !  o chemical_potential -- chemical potential
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


    !  local variables
    real*8 :: chemical_potential_l
    real*8 :: chemical_potential_r
    real*8 :: diff_l
    real*8 :: diff_r
    real*8 :: diff_electrons
    real*8 :: diff_electrons_thr
    real*8 :: lowest_eigenvalue
    real*8 :: highest_eigenvalue
  
  
    character*140 :: info_str
  
    !  counters
    integer :: i_state
    integer :: i_counter
    integer :: i_k_points
  
    if (t_out) then
      write(info_str,'(2X,A)') &
       "Determining occupation numbers for Kohn-Sham eigenstates."
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
    end if

    i_counter = 0
    !  write(use_unit,*) "get_occ_numbers.... n_spin=", n_spin, "n_states=", n_states
  
    lowest_eigenvalue = KS_eigenvalue(1,1)

    do i_k_points = 1, n_k_points, 1
      do i_state = 1, n_states, 1
        if (KS_eigenvalue(i_state,i_k_points) .lt. lowest_eigenvalue) then
          lowest_eigenvalue = KS_eigenvalue(i_state,i_k_points)
        end if
      end do
    end do

  
    highest_eigenvalue = KS_eigenvalue(n_states,1)

    do i_k_points = 1, n_k_points, 1
      do i_state = 1, n_states, 1
        if (KS_eigenvalue(i_state,i_k_points) .gt. highest_eigenvalue) then
          highest_eigenvalue = KS_eigenvalue(i_state,i_k_points)
        end if
      end do
    end do
  
    !  initialize zeroin algorithm
    !  to avoid problems if only one state is calculated (so for the H-atom)
    !  do not initialize mit KS_eigenvalue(n_states), because it is then
    !  identical to KS_eigenvalue(1)
  
    if (lowest_eigenvalue.ne.highest_eigenvalue) then
      chemical_potential_r = highest_eigenvalue
    else
      chemical_potential_r = 0.d0
    end if
    chemical_potential_l = lowest_eigenvalue
  
    call check_norm_force_occ_periodic_v2(chemical_potential_l, KS_eigenvalue, n_electrons, occ_numbers, diff_l,&
      i_counter,current_spin)
    call check_norm_force_occ_periodic_v2(chemical_potential_r, KS_eigenvalue, n_electrons, occ_numbers, diff_r,&
      i_counter,current_spin)
    do while (diff_l * diff_r .gt. 0.d0)
      !     interval for chemical potential still not found
      !     Must extend intervals both ways: 
      !     * For zero electrons, need Fermi level below lowest state
      !     * For all electrons in one channel, need Fermi level above highest state
      chemical_potential_l = chemical_potential_l - 0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
      chemical_potential_r = chemical_potential_r + 0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
      diff_l = diff_r
      call check_norm_force_occ_periodic_v2(chemical_potential_l, KS_eigenvalue, n_electrons, occ_numbers, diff_l,&
        i_counter,current_spin)
      call check_norm_force_occ_periodic_v2(chemical_potential_r, KS_eigenvalue, n_electrons, occ_numbers, diff_r,&
        i_counter,current_spin)
      !     write(use_unit,*) chemical_potential_l, chemical_potential_r
      !     If the right interval cannot be found, stop potentially
      !     infinite loop
      if (i_counter .gt. max_zeroin) then
        write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* Too many iterations needed to find good interval around E_F."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
           "* According to this check, you may not be able to find one."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
           "Please look carefully at your settings."
        call localorb_info(info_str,use_unit,'(A)')
        stop
      end if
    end do
  
    call zeroin_occ_periodic_v2(chemical_potential_l, chemical_potential_r, diff_l, diff_r,&
         fermi_acc, occupation_acc, max_zeroin,n_states, KS_eigenvalue, n_electrons, &
         chemical_potential, diff_electrons, occ_numbers, i_counter,current_spin)

    !     call check_norm one more time, just to make sure we really have the 
    !     occupation numbers for the final fermi level [this is critical for H]
    call check_norm_force_occ_periodic_v2( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers,&
         diff_electrons, i_counter,current_spin)
    !     and decrement i_counter which was falsely incremented ...
    i_counter = i_counter - 1
  
    !       If zeroin did not converge to the desired accuracy, issue
    !       warning and stop, for now; we can create nicer behavior after
    !       we gain more experience.
    if (i_counter .gt. max_zeroin) then
      write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,'(1X,A)') "* Too many iterations needed to find correct occupation of levels."
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,'(1X,A)') "* Check variables occupation_acc, max_zeroin before continuing your calculation."
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* Status of get_occ_numbers when aborted, 1st k-point: "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* "
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,*) "* State #         Eigenvalue       Occ. number"
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str,*) "----------"
      call localorb_info(info_str,use_unit,'(A)')
      do i_state = 1, n_states, 1
        write (info_str,*) i_state, KS_eigenvalue(i_state,1), occ_numbers(i_state,1)
        call localorb_info(info_str,use_unit,'(A)')
      end do
      stop
    end if
  
    !     finally, occupation thresholding - if needed
    if (occupation_thr.gt.0.d0) then
      call threshold_occ_numbers( n_electrons, occ_numbers, diff_electrons_thr)
    end if
  
    if (t_out) then
      write (info_str, '(2X, A,I1,A,F15.10)') "| Chemical potential (Fermi level) for spin channel "&
      ,current_spin, " in eV               : ", &
      chemical_potential * hartree
      call localorb_info(info_str,use_unit,'(A)')
      write (info_str, '(2X, A, E10.4)') "| Error in electron count due to remaining E_F inaccuracy: ", diff_electrons
      call localorb_info(info_str,use_unit,'(A)')
      if (occupation_thr.gt.0.d0) then
        write (info_str, '(2X, A, E10.4)') "| Error in electron count after thresholding : ", diff_electrons_thr
        call localorb_info(info_str,use_unit,'(A)')
      end if
      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)')
    end if
  
  
  end subroutine get_occupation_numbers_occ_periodic_v2



!******
!-------------------------------------------------------------------------------
!****s* force_occupation/bfgs_constraint_occ_v2
!  NAME
!    bfgs_constraint_occ_v2
!  SYNOPSIS
      subroutine bfgs_constraint_occ_v2(n_arguments, arguments, &
           max_fn_evals, eps, change_in_f, arg_mag, delta0, &
           n_fn_evals, hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin  &
           )
!  PURPOSE
!    Find a solution to the equations f(x) = 0 using the BFGS-method
!    and finite-differencing for computing the derivatives.
!
!  ARGUMENTS
      integer :: n_arguments
      real*8 :: arguments(n_arguments)
      integer :: max_fn_evals
      real*8 :: eps
      real*8 :: change_in_f
      real*8 :: delta0
!       external :: f_external
      real*8 hamiltonian( n_basis*(n_basis+1)/2, n_spin )
      real*8 overlap_matrix( n_basis*(n_basis+1)/2 )
      real*8, dimension(n_states, n_spin) :: KS_eigenvalue
      real*8, dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
      real*8, dimension(n_states, n_spin) :: occ_numbers
      real*8 :: n_electrons
      real*8 :: chemical_potential
      real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
      real*8, dimension(n_region, n_spin) :: electrons_in_region
      integer :: current_spin
      real*8 :: f_value
      integer :: n_fn_evals
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o n_arguments -- number of arguments for the objective function
!    o arguments -- arguments for the objective function
!    o max_fn_evals -- maximum number of evaluations of the objective function
!    o eps -- tolerance of the objective function, i.e., |f(x)| < eps
!    o change_in_f -- inital value of f, the objective function
!    o delta -- basic differencing parameter for the evaluation of the gradient of f
!    o f_external_constraint_occ -- external routine to evaluate the objective function f
!    o hamiltonian -- the Hamilton matrix
!    o overlap_matrix -- the overlap matrix
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o occ_numbers -- occupation numbers
!    o n_electrons -- total number of electrons
!    o chemical_potential -- current chemical potential
!    o constraint_proj -- projector to the constrained region
!    o electrons_in_region -- number of electrons in the constrained regions
!    o current_spin -- current spin channel
!  OUTPUT
!    o f_value -- final value obtained by the objective function f
!    o n_fn_evals -- number of evaluations of f taken
!  SOURCE

!     locals
      real*8 :: gradient(n_arguments)
      real*8 :: gradient_work(n_arguments)

      real*8 :: hessian_d(n_arguments)
      real*8 :: hessian_d_work(n_arguments)

      real*8 :: hessian_l(n_arguments,n_arguments)
      real*8 :: hessian_l_work(n_arguments,n_arguments)

      real*8 :: s(n_arguments)
      real*8 :: y(n_arguments)

      real*8 :: sHs
      real*8 :: try_arguments(n_arguments)
      real*8 :: arg_mag(n_arguments)

      logical :: converged
      logical :: scaled_update
      logical :: global_improvement
      logical :: longer_step
      logical :: shorter_step

!      integer :: info
      integer :: method

      real*8 :: f_value_1
      real*8 :: f_value_2
!      real*8 :: f_previous
      real*8 :: alpha
      real*8 :: alpha_scale
!      real*8 :: y_scale
      real*8 :: sigma
      real*8 :: s_times_y
      real*8 :: s_times_g
      real*8 :: total_alpha
      real*8 :: deps
      real*8 :: aeps
 !     real*8 :: temp
      real*8 :: diff_in_f
      real*8 :: norm_s

      real*8 :: delta

!      character*120 :: info_str

!      integer :: cause_of_exit

!     counters
      integer :: i_index

!     minimal values for alpha and delta
      aeps = 0.1d0*eps
      deps = 0.1d0*aeps
!      deps = min(0.1d0*aeps,0.001d0*delta0)

!     initialize hessian to be the identity matrix
      do i_index = 1, n_arguments, 1
         hessian_d(i_index) = 1.0d0
      enddo
      hessian_l = 0.0d0

      method = 1
      diff_in_f = change_in_f
      delta = delta0

      n_fn_evals = 1
      call f_external_constraint_occ(n_arguments, arguments, &
           f_value, n_fn_evals, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )

      call evaluate_gradient_constraint_occ(n_arguments, arguments, f_value, &
           n_fn_evals, delta, arg_mag, gradient, method, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )

      converged = .false.

!     this is the main BFGS loop
      do while ((.not.converged).and.(n_fn_evals.lt.max_fn_evals))

         hessian_d_work = hessian_d
         hessian_l_work = hessian_l
         s = -gradient
         call solve_factored_problem(hessian_d_work, &
              hessian_l_work, n_arguments, s)

         sHs = 0.0d0
         norm_s = 0.0d0
         do i_index = 1, n_arguments, 1
            sHs = sHs + gradient(i_index)*s(i_index)
            norm_s = norm_s + s(i_index)**2
         enddo
         norm_s = sqrt(norm_s)
!         alpha = -2.0d0*diff_in_f/sHs
!         alpha = 0.5d0

         alpha = -diff_in_f/(norm_s*sHs)
         if (alpha.gt.1.0d0) then
            alpha = 1.0d0
         end if

         global_improvement = .false.
         longer_step = .false.
         shorter_step = .false.
         total_alpha = 0.0d0

!         print *,'Spin:', current_spin
!         print *,'Fn_evals:', n_fn_evals
!         print *,'differece_method: ', method
!         print *,'eps: ', eps
!         print *,'aeps: ', aeps
!         print *,'deps: ', deps
!         print *,'sHs: ', sHs
!         print *,'norm_s: ', norm_s
!         print *,'alpha: ', alpha
!         print *,'delta: ', delta
!         print *,'f_value: ', f_value

!     first line search step, try = current + alpha*s
         try_arguments = arguments + alpha*s
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, try_arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         if (f_value_1.lt.f_value) then
            f_value = f_value_1
            arguments = try_arguments
            total_alpha = total_alpha + alpha
            global_improvement = .true.
            longer_step = .true.
         else
            shorter_step = .true.
         end if

!     if first step was successful then set alpha <- 2.0d0*alpha as long as
!     there is improvement
         do while (longer_step)

            alpha = 2.0d0*alpha
            try_arguments = arguments + alpha*s
            n_fn_evals = n_fn_evals + 1
            call f_external_constraint_occ(n_arguments, try_arguments, &
                 f_value_1, n_fn_evals, &
                 hamiltonian, &
                 overlap_matrix,KS_eigenvalue, KS_eigenvector, &
                 occ_numbers, &
                 n_electrons, chemical_potential, constraint_proj, &
                 electrons_in_region, &
                 current_spin &
                 )

            if (f_value_1.lt.f_value) then
               f_value = f_value_1
               arguments = try_arguments
               total_alpha = total_alpha + alpha
            else
               longer_step = .false.
            end if

         enddo

!     if there has not been any improvement so far, try setting
!     alpha <- 0.5d0*alpha as long as improvement is found
         do while ((shorter_step).and.(.not.global_improvement))

            alpha = 0.5d0*alpha
            try_arguments = arguments + alpha*s
            n_fn_evals = n_fn_evals + 1
            call f_external_constraint_occ(n_arguments, try_arguments, &
                 f_value_2, n_fn_evals, &
                 hamiltonian, &
                 overlap_matrix,KS_eigenvalue, KS_eigenvector, &
                 occ_numbers, &
                 n_electrons, chemical_potential, constraint_proj, &
                 electrons_in_region, &
                 current_spin &
                 )

            if (f_value_2.lt.f_value) then
               f_value = f_value_2
               arguments = try_arguments
               total_alpha = total_alpha + alpha
               global_improvement = .true.

            else

               alpha_scale = 0.1d0
               if (f_value_1 + f_value .gt. f_value_2 + f_value_2) then
                  alpha_scale = 1.0d0 + 0.5d0 * (f_value - f_value_1) / &
                       (f_value + f_value_1 - f_value_2 - f_value_2)
               end if
               if (alpha_scale .lt. 0.1d0) then
                  alpha_scale = 0.1d0
               end if
               alpha = alpha_scale*alpha

               try_arguments = arguments + alpha*s
               n_fn_evals = n_fn_evals + 1
               call f_external_constraint_occ(n_arguments, try_arguments, &
                    f_value_1, n_fn_evals, &
                    hamiltonian, &
                    overlap_matrix,KS_eigenvalue, KS_eigenvector, &
                    occ_numbers, &
                    n_electrons, chemical_potential, constraint_proj, &
                    electrons_in_region, &
                    current_spin &
                    )
               if (f_value_1.lt.f_value) then
                  f_value = f_value_1
                  arguments = try_arguments
                  total_alpha = total_alpha + alpha
                  global_improvement = .true.
               end if

            end if

            if (alpha.lt.aeps) then
               shorter_step = .false.
            end if

         enddo

!     next, update the Hessian according to the BFGS formula
         if (global_improvement) then

            call evaluate_gradient_constraint_occ(n_arguments, arguments, &
                 f_value, &
                 n_fn_evals, delta, arg_mag, gradient_work, &
                 method, &
                 hamiltonian, &
                 overlap_matrix, KS_eigenvalue, KS_eigenvector, &
                 occ_numbers, &
                 n_electrons, chemical_potential, constraint_proj, &
                 electrons_in_region, &
                 current_spin &
                 )

            s_times_g = 0.0d0
            do i_index = 1, n_arguments, 1
               s_times_g = s_times_g + &
                    s(i_index)*gradient_work(i_index)
            enddo

            s_times_y = s_times_g - sHs

            diff_in_f = f_value

!     update the Hessian only if (s_times_y.gt.0.0d0).and.(sHs.lt.0.0d0)
            if ((s_times_y.gt.0.0d0).and.(sHs.lt.0.0d0)) then

               if (s_times_y + total_alpha*sHs.le.0.0d0) then

!     normal update of the Hessian
                  scaled_update = .false.

                  sigma = -1.0d0/sHs

!                  print *, 'non-scaled update: ',
!     +                 s_times_y + total_alpha*sHs
!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, gradient, scaled_update)


!                  print *,'total_alpha', total_alpha
                  sigma = 1.0d0/(total_alpha*s_times_y)
                  y = gradient_work - gradient

!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, y, scaled_update)


               else

                  scaled_update = .false.
!                  scaled_update = .true.

!                  sigma = total_alpha/(s_times_y - total_alpha*sHs)

                  sigma = -1.0d0/sHs

!                  print *, 'scaled update: ',
!     +                 s_times_y + total_alpha*sHs
!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, gradient, scaled_update)

                  scaled_update = .false.

!                  y = gradient_work + (s_times_y*sigma - 1.0d0)*gradient
!                  sigma = 1.0d0 / (sigma*(s_times_y**2))

                  sigma = 1.0d0/s_times_y
                  y = gradient_work - gradient

!                  print *,'total_alpha', total_alpha
!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, y, scaled_update)

               end if

            end if

         else

            if (method.eq.1) then
               method = 2
!               delta = delta0
               call evaluate_gradient_constraint_occ(n_arguments, &
                    arguments, &
                    f_value, &
                    n_fn_evals, delta, arg_mag, &
                    gradient, &
                    method, hamiltonian, &
                    overlap_matrix, KS_eigenvalue, KS_eigenvector, &
                    occ_numbers, &
                    n_electrons, chemical_potential, &
                    constraint_proj, &
                    electrons_in_region, &
                    current_spin &
                    )

!$$$            elseif (method.eq.2) then
!$$$               method = 3
!$$$               delta = delta0
!$$$               call evaluate_gradient_constraint(n_arguments,
!$$$     +              arguments,
!$$$     +              f_value,
!$$$     +              n_fn_evals, delta, arg_mag, f_external_constraint_occ,
!$$$     +              gradient,
!$$$     +              method, hamiltonian,
!$$$     +              overlap_matrix, KS_eigenvalue, KS_eigenvector,
!$$$     +              occ_numbers,
!$$$     +              n_electrons, chemical_potential,
!$$$     +              constraint_proj,
!$$$     +              electrons_in_region,
!$$$     +              current_spin
!$$$     +              )

            elseif (delta.gt.deps) then
               delta = delta/2.0d0
               method = 1
               call evaluate_gradient_constraint_occ(n_arguments, &
                    arguments, &
                    f_value, &
                    n_fn_evals, delta, arg_mag, &
                    gradient, &
                    method, hamiltonian, &
                    overlap_matrix, KS_eigenvalue, KS_eigenvector, &
                    occ_numbers, &
                    n_electrons, chemical_potential, &
                    constraint_proj, &
                    electrons_in_region, &
                    current_spin &
                    )
            end if

!     end update hessian
         end if

!     test for convergence
         converged = (f_value.lt.eps)

!     main BFGS loop
      enddo

!     finally, one more evaluation of the target function is needed
!     to update all the variables - like eigenvalues, chemical potentials
!     and constraint potentials
      n_fn_evals = n_fn_evals + 1
      call f_external_constraint_occ(n_arguments, arguments, &
           f_value, n_fn_evals, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )

!      print *, 'On exit, f_value: ', f_value
!      print *, 'total_alpha: ', total_alpha
!      print *, 'alpha: ', alpha
!      print *, 'delta: ', delta

      end subroutine bfgs_constraint_occ_v2





!******
!-------------------------------------------------------------------------------
!****s* force_occupation/evaluate_gradient_constraint_occ
!  NAME
!    evaluate_gradient_constraint_occ
!  SYNOPSIS
      subroutine evaluate_gradient_constraint_occ(n_arguments, arguments, &
           f_value, n_fn_evals, delta, arg_mag, gradient, &
           method, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )
!  PURPOSE
!    Evaluate an approximation to the gradient of f at the points given by
!    the array arguments using finite differences.
!
!  ARGUMENTS
      integer :: n_arguments
      real*8 :: arguments(n_arguments)
      real*8 hamiltonian( n_basis*(n_basis+1)/2, n_spin )
      real*8 overlap_matrix( n_basis*(n_basis+1)/2 )
      real*8, dimension(n_states, n_spin) :: KS_eigenvalue
      real*8, dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
      real*8, dimension(n_states, n_spin) :: occ_numbers
      real*8 :: n_electrons
      real*8 :: chemical_potential
      real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
      real*8, dimension(n_region, n_spin) :: electrons_in_region
      integer :: current_spin
!       external :: f_external_constraint_occ
      real*8 :: f_value
      real*8 :: delta
      real*8, dimension(n_arguments) :: arg_mag
      integer :: n_fn_evals
      integer :: method
      real*8 :: gradient(n_arguments)
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o n_arguments -- number of arguments for the objective function
!    o arguments -- arguments for the objective function
!    o max_fn_evals -- maximum number of evaluations of the objective function
!    o eps -- tolerance of the objective function, i.e., |f(x)| < eps
!    o f_value -- value of the objective function
!    o delta -- basic differencing parameter for the evaluation of the gradient of f
!    o f_external_constraint_occ -- external routine to evaluate the objective function f
!    o hamiltonian -- the Hamilton matrix
!    o overlap_matrix -- the overlap matrix
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o occ_numbers -- occupation numbers
!    o n_electrons -- total number of electrons
!    o chemical_potential -- current chemical potential
!    o constraint_proj -- projector to the constrained region
!    o electrons_in_region -- number of electrons in the constrained regions
!    o current_spin -- current spin channel
!    o method -- selector for the finite difference scheme to use
!    o arg_mag -- magnitude of the arguments
!  OUTPUT
!    o gradient -- difference approximation to gradient of f
!  SOURCE

!     locals
      real*8 :: f_value_1, f_value_2, f_value_3, f_value_4
      integer :: i_index
      real*8 :: delta_scaled
      real*8 :: temp_arg

      do i_index = 1, n_arguments, 1

         delta_scaled = delta*arg_mag(i_index)
         temp_arg = arguments(i_index)

         select case(method)

      case(1)

         arguments(i_index) = temp_arg + delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         gradient(i_index) = (f_value_1 - f_value) / delta_scaled

      case(2)

         arguments(i_index) = temp_arg + delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg - delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, arguments, &
              f_value_2, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         gradient(i_index) = (f_value_1 - f_value_2) / &
              (2.0d0*delta_scaled)

      case(3)

         arguments(i_index) = temp_arg + 2.0d0*delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg + delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, arguments, &
              f_value_2, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg - delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, arguments, &
              f_value_3, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg - 2.0d0*delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external_constraint_occ(n_arguments, arguments, &
              f_value_4, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )


         gradient(i_index) = &
            (f_value_4 - 8.0d0*f_value_3 + 8.0d0*f_value_2 - f_value_1)/ &
            (12.0d0*delta_scaled)

      end select

      arguments(i_index) = temp_arg

      enddo

      end subroutine evaluate_gradient_constraint_occ

!*****
!-------------------------------------------------------------------------------
!****s* force_occupation/get_KS_orbitals_forced_occ_p0
!  NAME
!   get_KS_orbitals_forced_occ_p0
!  SYNOPSIS

subroutine get_KS_orbitals_forced_occ_p0(overlap_matrix,hamiltonian,n_electrons, &
     KS_eigenvalue,KS_eigenvector, KS_eigenvector_complex,  occ_numbers,  &
     chemical_potential)

!  PURPOSE
!  The subroutine calculates Kohn-Sham eigenvalues and eigenvectors.
!  It also calculates the occupation weights for eigenstates and 
!  value of the chemical potential.
!  USES

!   use dimensions
  use geometry
  use basis
!   use runtime_choices
  use constraint
  use localorb_io
  use mixing_constraint
  use bfgs
  use physics, only: rho_change
  use constants
  use synchronize_mpi
  use scalapack_wrapper
  use separate_core_states


  implicit none

!  ARGUMENTS

  real*8 :: hamiltonian( n_hamiltonian_matrix_size, n_spin )
  real*8 :: overlap_matrix( n_hamiltonian_matrix_size ) 

  real*8,     dimension(n_states, n_spin,n_k_points) :: occ_numbers
  real*8,     dimension(n_states, n_spin,n_k_points) :: KS_eigenvalue

  real*8,     dimension(n_basis, n_states, n_spin, n_k_points_task) ::  KS_eigenvector

  complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task) ::  KS_eigenvector_complex

  real*8 :: n_electrons
  real*8 :: chemical_potential    


! INPUTS
! o overlap_matrix -- overlap matrix
! o hamiltonian -- Hamiltonian matrix
! o n_electrons -- number of electrons in system
!
! OUTPUT
! o KS_eigenvalue -- Kohn-Sham eigenvalues
! o KS_eigenvector -- Kohn-Sham eigenvectors if real number eigenvectors are in use.
! o KS_eigenvector_complex -- Kohn-Sham eigenvectors if complex number eigenvectors are in use.
! o occ_number  -- occupation weights of the eigenstates
! o chemical_potential -- value of chemical potential
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





  real*8,    dimension(:,:),allocatable :: hamiltonian_w
  real*8,    dimension(:),  allocatable :: overlap_matrix_w
  complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
  complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex



  real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
  real*8, dimension(n_region, n_spin) :: constraint_fermi
  real*8, dimension(n_region, n_spin) :: electrons_in_region
  real*8, dimension(n_spin) :: constraint_fermitot

  integer, dimension(n_spin) :: n_eval
  real*8, dimension(:,:,:),allocatable :: work_ham
  real*8, dimension(:,:),allocatable :: work_ovl
  real*8, allocatable,dimension(:,:,:):: occ_temp
  real*8, allocatable,dimension(:,:,:):: eigen_value_temp
  integer:: info
  integer, dimension(n_spin) :: n_homo




  ! local variables

  logical :: constraint_converged

  logical :: t_out = .true.

  real*8 :: avg_zero
  real*8 :: avg_potential

  character*120 :: info_str

  !  counters

  integer :: i_spin,j_spin, i_k_point,i_force_occ
  integer :: i_basis,j_basis,k_basis,l_basis,i_basis_1,i_basis_2,i_index_real
  integer :: i_region
  integer :: i_state
  integer :: i_hamiltonian
  integer :: number_constraint_iter

  !     variables

  integer :: n_arguments
  real*8, dimension(:), allocatable :: arguments
  real*8 :: eps
  real*8 :: fmin

  real*8, dimension(:), allocatable :: gradient
  real*8, dimension(:), allocatable :: hessian
  real*8, dimension(:), allocatable :: work
  real*8 :: dfn
  real*8, dimension(:), allocatable :: xm
  real*8 :: hh
  real*8 :: min_rho_change

  integer :: mode
  integer :: maxfn
  integer :: iexit, i_k

  real*8, dimension(n_spin) :: n_electrons_in_spin
  real*8, dimension(n_spin) :: chemical_potential_spin
  real*8 :: diff_in_chem_pot

  real*8,dimension(n_states,2):: temp_n_electrons_forced_state

  real*8 :: biggest,temp

       real*8 :: spin_shift_fsm (n_spin)

  !  begin work

  write(info_str,'(2X,A)') ''
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str,'(2X,A,A)') &
       "Updating Kohn-Sham eigenvalues, eigenvectors,", &
       " occupation numbers, and Fermi level."
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  if (.not.use_cg) then
     occ_numbers = 0.0d0
  end if
  
!   if (n_periodic .eq. 0) then
!     call determine_corehole_projection(overlap_matrix,KS_eigenvector)
!   else
! 
!   end if

  ! if locally constrained DFT requested, prepare for it
  if (use_constraint) then
     !  add constraint potential(s) from previous iterations to Hamiltonian

     !  In case of a constraint calculation we need to store the original hamiltonian
     !  because we need it for any inner iteration 
     !  We also allocate a full work version of the ovlp matrix here.

     call allocate_hamiltonian_work  &
          ( overlap_matrix)

     !  copy original hamiltonian to hamiltonian_work, to save it from 
     !  destruction by LAPACK.
     do i_spin=1, n_spin, 1
        do i_hamiltonian=1, n_basis*(n_basis+1)/2
           hamiltonian_work(i_hamiltonian,i_spin) &
                =hamiltonian(i_hamiltonian,i_spin)
        enddo
     enddo

     !       If a constraint was requested, we may already know constraint potentials for different regions. 
     !       These must be added to hamiltonian_work. 

     if (n_active_regions.gt.1) then
        call add_constraint_potentials &
             ( overlap_matrix, hamiltonian &
             )
     end if

     t_out = .false.

  end if ! use_constraint

  call solve_KS_eigen(overlap_matrix, hamiltonian, &
  &                   KS_eigenvalue,KS_eigenvector, KS_eigenvector_complex)

  ! and obtain the usual occupation numbers

! the following is only valid for nonperiodic calculations

  if (n_periodic.eq.0) then

    call determine_corehole_projection(overlap_matrix,KS_eigenvector)
  
    call get_occupation_numbers_occ_p0 &
      ( KS_eigenvalue, n_electrons, t_out, &
      occ_numbers, chemical_potential &
      )

    !  At this point, if an electron number constraint was requested, attempt to enforce
    !  it using auxiliary potentials      

    if (use_constraint) then

      ! count number of electrons in each region, each spin
      do i_spin = 1, n_spin, 1
        call get_electrons_per_region_v2 &
          ( i_spin, KS_eigenvector, overlap_matrix, occ_numbers, &
          constraint_proj,  &
          electrons_in_region &
          )
      enddo
  
      if (constraint_it_lim.gt.1) then
      ! if iterative determination of constraint_potential requested, check whether
      ! the correct electron numbers have already been reached
  
        constraint_converged = .true.
  
        do i_spin = 1, n_spin, 1
          do i_region = 1, n_active_regions, 1    
            constraint_converged = constraint_converged .and. &
            ( dabs(electrons_in_region(i_region,i_spin)- &
            constraint_electrons(i_region,i_spin)).lt. &
            constraint_precision )
          enddo
        enddo
  
        if (constraint_converged) then
          write(info_str,'(2X,A,A)') &
          "Locally constrained DFT: ", &
          "Electron counting constraint fulfilled on entry."
          call localorb_info(info_str,use_unit,'(A)')
        else
          write(info_str,'(2X,A,A)') &
          "Locally constrained DFT: ", &
          "Constraint on electron numbers not yet fulfilled." 
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
          "Entering inner self-consistency loop."
          call localorb_info(info_str,use_unit,'(A)')
  
          ! restore hamiltonian
          hamiltonian = hamiltonian_work
  
          end if
  
      end if
  
      ! begin inner self-consistency loop to determine 
      ! constraint potentials
      if( (mixer_constraint.eq.1) ) then
        call allocate_pulay_constraint ( )
      end if
  
      number_constraint_iter = 0
  
      if (mixer_constraint.le.1) then
  
      ! If constraint already converged, still do a final adjustment of 
      ! the eigenvalue levels etc, just for consistency's sake!!
      !
      ! This is a patch for a problem where significant shifts of the 
      ! eigenvalues of spin-up vs spin-down against each other can arise
      ! even though the overall spin-up <-> spin-down occupancies 
      ! do not change at all. 
  
        call determine_corehole_projection(overlap_matrix,KS_eigenvector)
  
        if (constraint_converged) then
          n_electrons_in_spin = 0.d0
          do i_spin = 1, n_spin, 1
            do i_region = 1, n_active_regions, 1 
              n_electrons_in_spin(i_spin) =  &
              n_electrons_in_spin(i_spin) + &
              constraint_electrons(i_region,i_spin)
            enddo
  
            call get_occupation_numbers_occ_v2 &
            ( KS_eigenvalue(:,i_spin,1),  &
            n_electrons_in_spin(i_spin), &
            t_out, &
            occ_numbers(:,i_spin,1),  &
            chemical_potential_spin(i_spin), &
            i_spin&
            )
  
            ! The number of electrons in region 1 is now by definition correct.
            do i_region = 1, n_active_regions, 1
              electrons_in_region(i_region,i_spin) =  &
                constraint_electrons(i_region,i_spin)
            enddo
  
          enddo
  
          if (n_spin.eq.2) then
  
            diff_in_chem_pot = &
              chemical_potential_spin(1) -  &
              chemical_potential_spin(2)
  
            if (constraint_debug) then
              call localorb_info('',use_unit,'(A)') 
              write (info_str, '(2X,A,A,F15.10,A)') &
                'Difference in chemical potential between ', &
                'spin channels:', diff_in_chem_pot*hartree,  &
                " eV."
              call localorb_info(info_str,use_unit,'(A)')
            end if
  
            spin_shift(1) = n_electrons_in_spin(2)/n_electrons *  &
              diff_in_chem_pot
            spin_shift(2) = - n_electrons_in_spin(1)/n_electrons *  &
              diff_in_chem_pot
  
            do i_spin = 1, n_spin, 1
  
            ! adjust constraint potentials to reflect the spin shift
            do i_region = 1, n_active_regions, 1
  
              constraint_potential(i_region,i_spin) =  &
              constraint_potential(i_region,i_spin) - &
              spin_shift(i_spin)
  
          enddo
  
          ! adjust the final eigenvalues to reflect the spin shift
          do i_state = 1, n_states, 1
              KS_eigenvalue(i_state,i_spin,1) =  &
              KS_eigenvalue(i_state,i_spin,1) - &
              spin_shift(i_spin)
          enddo
  
          chemical_potential_spin(i_spin) =  &
            chemical_potential_spin(i_spin) - spin_shift(i_spin)
  
        enddo
  
        ! Both spin channel chemical potentials must now be equal
        ! simply by construction!
        chemical_potential = chemical_potential_spin(1)
  
        ! end adjustment of different spin channels
        end if
  
        ! end the patch if the constraint was already converged.
      end if
  
      do while ( (.not.constraint_converged) .and. &
          (number_constraint_iter.lt.constraint_it_lim) )
  
        number_constraint_iter = number_constraint_iter+1
  
        if (constraint_debug) then
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A)')  &
            "------------------------------------------------------"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,I5)') &
            "Locally constrained DFT: Inner SCF iteration no. ", &
            number_constraint_iter
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A)')  &
            "------------------------------------------------------"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        end if
  
        ! store present constraint potential for later mixing
        constraint_pot_old = constraint_potential
  
        ! Determine region-specific Fermi levels which would be needed to
        ! satisfy the constraint in each region, based on the present
        ! KS eigenstates alone
  
        call determine_corehole_projection(overlap_matrix,KS_eigenvector)
        do i_region = 1, n_active_regions, 1
        do i_spin = 1, n_spin, 1
          call get_constraint_fermi_occ &
            ( KS_eigenvalue(:,i_spin,1),  &
            constraint_electrons(i_region,i_spin), &
            constraint_proj(:,i_region,i_spin),  &
            constraint_fermi(i_region,i_spin), &
            i_spin&
            )
        enddo
        enddo
  
        ! write results       
        if (constraint_debug) then
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A)')  &
            "Required Fermi levels in different subsystems:"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,5X,A,I2,5X,A,I5)')  &
            "| Region  ", ("Spin ", i_spin, i_spin=1,n_spin,1)
        call localorb_info(info_str,use_unit,'(A)')
        do i_region = 1, n_active_regions, 1
          write(info_str,'(2X,A,I6,2X,F10.5,2X,F10.5)') &
            "| ", &
            constraint_region_label(i_region), &
            (constraint_fermi(i_region,i_spin), i_spin=1,n_spin,1) 
          call localorb_info(info_str,use_unit,'(A)')
        enddo
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        end if
  
        ! compute average zero level to which all potentials will be referenced
        ! variant 1: avg_zero = constraint_fermi(1,1)
        ! variant 2: avg_zero = chemical_potential
        ! variant 3:
        avg_zero = 0.d0
        do i_spin = 1, n_spin, 1
        do i_region = 1, n_active_regions, 1
          avg_zero = avg_zero +  &
            constraint_electrons(i_region,i_spin) * &
            constraint_fermi(i_region,i_spin)
        enddo
        enddo
        avg_zero = avg_zero / n_electrons
  
        ! To each region, we apply a constraint potential that would shift
        ! the local Fermi level to the current overall Fermi level.
        ! We arbitrarily set the constraint potential of the first region 
        ! to zero ... 
        do i_region=1, n_active_regions, 1
        do i_spin=1, n_spin, 1
          delta_constraint_potential(i_region,i_spin)= &
            avg_zero - &
            constraint_fermi(i_region,i_spin)
        enddo
        enddo
  
        ! write results       
        if (constraint_debug) then
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,A)')  &
            "Suggested local potential change", &
            " in different subsystems:"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,5X,A,I2,5X,A,I5)')  &
            "| Region  ", ("Spin ", i_spin, i_spin=1,n_spin,1)
        call localorb_info(info_str,use_unit,'(A)')
        do i_region = 1, n_active_regions, 1
          write(info_str,'(2X,A,I6,2X,F10.5,2X,F10.5)') &
            "| ", &
            constraint_region_label(i_region), &
            ( delta_constraint_potential(i_region,i_spin),  &
            i_spin=1,n_spin,1)
          call localorb_info(info_str,use_unit,'(A)')
        enddo
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        end if
  
        ! Damp additionally needed constraint potential by simple linear mixing factor
  
        if ( (mixer_constraint.eq.0).or.(number_constraint_iter.le. &
          ini_linear_mixing_constraint) ) then
  
        do i_region = 1, n_active_regions, 1
          do i_spin = 1, n_spin, 1
              constraint_potential(i_region,i_spin) =   &
              constraint_pot_old(i_region,i_spin) +  &
              constraint_mix(1) *  &
              delta_constraint_potential(i_region,i_spin)  
          enddo
        enddo
  
        end if
  
        if (mixer_constraint.eq.1) then
  
        if (number_constraint_iter.le.  &
            ini_linear_mixing_constraint) then
  
          do i_spin = 1, n_spin, 1
              do i_region = 1, n_active_regions, 1
            delta_potential_mixing(i_region,i_spin) =  &
                  delta_constraint_potential(i_region,i_spin)  &
                  - delta_constraint_potential(1,i_spin)
              enddo
          enddo
  
          call prepare_pulay_mixing_constraint   &
            ( delta_potential_mixing )
  
        else
  
          do i_spin = 1, n_spin, 1
              do i_region = 1, n_active_regions, 1
            delta_potential_mixing(i_region,i_spin) =  &
                  delta_constraint_potential(i_region,i_spin)  &
                  - delta_constraint_potential(1,1)
              enddo
          enddo
  
          call execute_pulay_mixing_constraint(  &
            number_constraint_iter, constraint_potential,  &
            delta_potential_mixing )
  
        end if
  
        end if
  
        ! write results       
        if (constraint_debug) then
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,A)') &
            "Local constraint potentials in different subsystems after mixing:"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,5X,A,I2,5X,A,I5)') "| Region  ", ("Spin ", i_spin, i_spin=1,n_spin,1)
        call localorb_info(info_str,use_unit,'(A)')
        do i_region = 1, n_active_regions, 1
          write(info_str,'(2X,A,I6,2X,F10.5,2X,F10.5)')  "| ",  &
            constraint_region_label(i_region),  &
            ( constraint_potential(i_region,i_spin), i_spin=1,n_spin,1)
          call localorb_info(info_str,use_unit,'(A)')
        enddo
        write(info_str,*)
        call localorb_info(info_str,use_unit,'(A)')
        end if
  
        !  Restore original hamiltonian ...
        do i_spin=1, n_spin, 1
        do i_hamiltonian=1, n_basis*(n_basis+1)/2
          hamiltonian_work(i_hamiltonian,i_spin)  &
            =hamiltonian(i_hamiltonian,i_spin)
        enddo
        enddo
  
        ! ... and add the current constraint potentials.      
        call add_constraint_potentials  &
          ( overlap_matrix, hamiltonian_work  &
          )
  
        ! solve for new KS eigenstates and -energies ...
        do i_spin = 1, n_spin, 1
  
        call improve_real_eigenfunctions  &
            ( overlap_matrix, hamiltonian_work(:,i_spin),   &
            t_out,   &
            KS_eigenvalue(:,i_spin,1), KS_eigenvector(:,:,i_spin,1), 1  &
            )
  
        enddo
        ! ... and obtain the updated global Fermi level and
        ! and global occupation numbers
  
        call determine_corehole_projection(overlap_matrix,KS_eigenvector)
  
        call get_occupation_numbers_occ_p0  &
          ( KS_eigenvalue, n_electrons, t_out,  &
          occ_numbers, chemical_potential &
          )
  
        ! update constraint projectors and count number of electrons in each region, each spin
        do i_spin = 1, n_spin, 1
        call get_electrons_per_region_v2  &
            ( i_spin, KS_eigenvector, overlap_matrix, occ_numbers,  &
            constraint_proj,   &
            electrons_in_region  &
            )
        enddo
  
        ! Finally, check for convergence:
        constraint_converged = .true.
  
        do i_spin = 1, n_spin, 1
        do i_region = 1, n_active_regions, 1   
          constraint_converged = constraint_converged .and.  &
            ( dabs(electrons_in_region(i_region,i_spin)-  &
            constraint_electrons(i_region,i_spin)).lt.  &
            constraint_precision )
        enddo
        enddo
  
        if (constraint_debug) then
        if (constraint_converged) then
          write(info_str,'(2X,A,A)')  &
            "Electron counting constraint fulfilled - ",  &
            "inner scf cycle converged."
          call localorb_info(info_str,use_unit,'(A)')
        else
          write(info_str,'(2X,A,A,A)')  &
            "Constraint on electron numbers not yet fulfilled."
          call localorb_info(info_str,use_unit,'(A)')
        end if
        end if
  
        ! end of inner SCF loop.
      enddo
  
      !     end if (mixer_constraint.le.1)
      end if
  
      if (mixer_constraint.eq.2) then
      !           call BFGS implementation
  
      !            first, allocations for BFGS
      n_arguments = n_active_regions - 1
  
      if ((n_arguments.gt.0).and.(.not.constraint_converged))   &
          then
        ! This is where we actually run through the optimizer again ...
  
        if (.not.allocated(arguments)) then
        allocate(arguments(n_arguments),stat=info)
        call check_allocation(info, 'overlap_matrix_w              ')
        end if
        if (.not.allocated(xm)) then
        allocate( xm(3*n_arguments),stat=info) 
        call check_allocation(info, 'xm                            ')
        end if
  
        ! parameters for BFGS
        ! (consult the module bfgs.f)
  
        ! magnitude of the solution, xm > 0
        xm = 0.1d0
  
        ! step lenght for calculating the gradient
        hh = 1.0e-6
  
        ! tolerance of convergence with respect to the argumetns, x
        eps = constraint_precision
  
        ! mode of utilisation, mode=1 corresponds to no Hessian given
        mode = 1
  
        ! maximum number of function evaluations premitted
        maxfn = constraint_it_lim
  
        ! loop over spins, so that the spin channels don't mix here
        do  i_spin = 1, n_spin, 1
  
        ! initialize number of required electrons per spin channel
        n_electrons_in_spin(i_spin) = 0.0d0
        do i_region = 1, n_active_regions, 1
          n_electrons_in_spin(i_spin) =   &
            n_electrons_in_spin(i_spin) +  &
            constraint_electrons(i_region,i_spin)
        enddo
  
        if (constraint_debug) then
          write(info_str,'(2X,A,A,I5,A)')  &
            '------------------',  &
            ' Iteration for spin channel ', i_spin,  &
            ' ------------------'
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)')  &
            'Region               Constraint potential'
          call localorb_info(info_str,use_unit,'(A)')
          do i_region = 1, n_active_regions, 1
              write(info_str,'(2X,I5,20X,F15.10)')  &
              i_region,   &
              constraint_potential(i_region,i_spin)*hartree
              call localorb_info(info_str,use_unit,'(A)')
          enddo
        end if
  
        ! preload the chemical potential
        chemical_potential_spin(i_spin) =   &
            chemical_potential
  
        ! prepare the arguments according to the solution last iteration
        do i_region = 1, n_arguments, 1
          arguments(i_region) =   &
            constraint_potential(i_region,i_spin) -  &
            constraint_potential(n_active_regions,i_spin)
        enddo
  
        ! initialize fmin
        fmin = 0.d0
  
        ! dfn is the expected change in the target functional
        dfn = 0.0d0
        do i_region = 1, n_active_regions, 1
          dfn = dfn + &
            ( electrons_in_region(i_region,i_spin) - &
            constraint_electrons(i_region,i_spin) )**2
        enddo
  
        ! call bfgs_constraint_v2 from the bfgs.f module
        call bfgs_constraint_occ_v2(  &
            n_arguments, arguments, maxfn, eps, dfn, xm,  &
            hh, number_constraint_iter, hamiltonian, &
            overlap_matrix, KS_eigenvalue, KS_eigenvector, &
            occ_numbers, n_electrons_in_spin(i_spin),  &
            chemical_potential_spin(i_spin), &
            constraint_proj, electrons_in_region,  &
            i_spin &
            )
        n_eval(i_spin) = number_constraint_iter 
  
        ! end loop over spin channels
        enddo
  
      else if ((n_arguments.eq.0).or.constraint_converged) then
        ! Either: Only one region, but with a spin constraint. 
        ! Or: The constraint was already converged upon entry.
        ! In either case, simply redetermine needed Fermi level 
        ! shifts for each spin channel in a single shot once, but 
        ! do not affect the balance within each individual spin 
        ! channel.
  
        call determine_corehole_projection(overlap_matrix,KS_eigenvector)
  
        do i_spin = 1, n_spin, 1
  
        n_electrons_in_spin(i_spin) = 0.d0
        do i_region = 1, n_active_regions, 1
          n_electrons_in_spin(i_spin) =  &
            n_electrons_in_spin(i_spin) + &
            constraint_electrons(i_region,i_spin)
        enddo
  
        call get_occupation_numbers_occ_v2 &
            ( KS_eigenvalue(:,i_spin,1),  &
            n_electrons_in_spin(i_spin), &
            t_out, &
            occ_numbers(:,i_spin,1),  &
            chemical_potential_spin(i_spin), &
            i_spin&
            )
  
        ! Do one final count of the electrons in each region,
        ! just to catch any possible residual deviations within 
        ! the constraint convergence criterion. 
        !
        ! THIS IS IMPORTANT TO ENSURE THAT THE ENERGY CORRECTION
        ! DUE TO EIGENVALUE SHIFTS IS EXACTLY ZERO ON AVERAGE, AS IT SHOULD BE.
        !
        if (n_active_regions.gt.1) then
          call get_electrons_per_region_v2 &
            ( i_spin, KS_eigenvector, overlap_matrix, &
            occ_numbers, constraint_proj,  &
            electrons_in_region  &
            ) 
        else
          electrons_in_region(1,i_spin) =   &
            constraint_electrons(1,i_spin)
        end if
  
        enddo
  
        n_eval = 0
  
      else
        write(use_unit,*) "* No constraint regions defined?"
        stop
      end if
  
      if (n_spin.eq.2) then
  
        !              Balance the fermi levels of the spin channels
  
        diff_in_chem_pot =  &
          chemical_potential_spin(1) -   &
          chemical_potential_spin(2)
  
        if (constraint_debug) then
        call localorb_info('',use_unit,'(A)')
        write (info_str, '(2X,A,A,F15.10,A)')  &
            'Difference in chemical potential between ',  &
            'spin channels:', diff_in_chem_pot*hartree,   &
            " eV."
        call localorb_info(info_str,use_unit,'(A)')
        end if
  
        ! Shift all constraint potentials between the two spin channels
        ! to create one uniform Fermi level for the overall system
  
        ! Yes, the weighted average requires that 1 <-> 2 are exchanged
        ! in the equations. 
        spin_shift(1) = n_electrons_in_spin(2)/n_electrons *   &
          diff_in_chem_pot
        spin_shift(2) = - n_electrons_in_spin(1)/n_electrons *   &
          diff_in_chem_pot
  
        do i_spin = 1, n_spin, 1
  
        ! adjust constraint potentials to reflect the spin shift
        if (n_active_regions.gt.1) then
          ! in this case, an initial shift makes sense and must be readded here.
          do i_region = 1, n_active_regions, 1
              constraint_potential(i_region,i_spin) =   &
              constraint_potential(i_region,i_spin) -  &
              spin_shift(i_spin)
          enddo
        else
          ! plain fixed-spin-moment constraint - initial shift was not applied.
          constraint_potential(1,i_spin) =   &
            - spin_shift(i_spin)
        end if
  
        ! adjust the final eigenvalues to reflect the spin shift
        do i_state = 1, n_states, 1
          KS_eigenvalue(i_state,i_spin,1) =   &
            KS_eigenvalue(i_state,i_spin,1) - &
            spin_shift(i_spin)
        enddo
  
        chemical_potential_spin(i_spin) =  &
            chemical_potential_spin(i_spin) - spin_shift(i_spin)
  
        enddo
  
        ! Both spin channel chemical potentials must now be equal
        ! simply by construction!
        chemical_potential = chemical_potential_spin(1)
  
        ! end adjustment of spin channels if needed
      end if
  
      !     deallocations
      if (allocated(arguments)) then
        deallocate(arguments)
      end if
      if (allocated(gradient)) then
        deallocate(gradient)
      end if
      if (allocated(hessian)) then
        deallocate(hessian)
      end if
      if (allocated(work)) then
        deallocate(work)
      end if
      if (allocated(xm)) then
        deallocate(xm)
      end if
  
      ! Finally, check for convergence:
      constraint_converged = .true.
  
      do i_spin = 1, n_spin, 1
        do i_region = 1, n_active_regions, 1   
        constraint_converged = constraint_converged .and. &
            ( dabs(electrons_in_region(i_region,i_spin)- &
            constraint_electrons(i_region,i_spin)).lt. &
            constraint_precision )
        enddo
      enddo
  
      !     end if (mixer_constraint.eq.3)
      end if
  
      ! after end of loop, compute average auxiliary potential for all subsystems and shift all eigenvalues,
      ! constraint_potentials uniformly so as to reduce the energy correction to zero.
  
      avg_potential = 0.d0
      do i_spin = 1, n_spin, 1
      do i_region = 1, n_active_regions, 1
        avg_potential = avg_potential +  &
          electrons_in_region(i_region,i_spin) * &
          constraint_potential(i_region,i_spin)
      enddo
      enddo
      avg_potential = avg_potential / n_electrons
  
      do i_region = 1, n_active_regions, 1
      do i_spin = 1, n_spin, 1
        constraint_potential(i_region,i_spin) =  &
          constraint_potential(i_region,i_spin) - &
          avg_potential
      enddo
      enddo
  
      do i_spin = 1, n_spin, 1
      do i_state = 1, n_states, 1
        KS_eigenvalue(i_state,i_spin,1) =  &
          KS_eigenvalue(i_state,i_spin,1) - &
          avg_potential
      enddo
      enddo
  
      chemical_potential = chemical_potential - avg_potential
  
      ! Write all important results ...
      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)')
  
      ! issue warning if the inner scf loop was exited without convergence.
      if(constraint_converged.and.(mixer_constraint.le.1)) then
  
      write(info_str,'(2X,A,I5,A)') "Locally constrained DFT: Converged after ",  &
          number_constraint_iter, " inner scf iterations."
      call localorb_info(info_str,use_unit,'(A)')
  
      else if (constraint_converged) then
      do i_spin = 1, n_spin, 1
  
        write(info_str,'(2X,A,I5,A,I5,A)')"Locally constrained DFT: Spin ", i_spin,  &
          " converged after ", n_eval(i_spin), " iterations."
        call localorb_info(info_str,use_unit,'(A)')
  
      enddo
      else if (mixer_constraint.eq.2) then
  
      write(info_str,'(2X,A)')"Locally constrained DFT: NO CONVERGENCE." 
      call localorb_info(info_str,use_unit,'(A)')
  
      do i_spin = 1, n_spin, 1
        write(info_str,'(2X,A,I5,A,I5,A)') &
          "Spin ", i_spin,  &
          ": ", &
          n_eval(i_spin), " iterations."
        call localorb_info(info_str,use_unit,'(A)')
      enddo
      else
      write(info_str,'(2X,A,I5,A)') &
          "Locally constrained DFT: NO CONVERGENCE after ",  &
          number_constraint_iter, " inner scf iterations."
      call localorb_info(info_str,use_unit,'(A)')
      endif
  
      write(info_str,'(2X,A)') "| "
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A,A)') "| Final constraint potentials in different subsystems:"
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A,10X,A,I2,12X,A,I2)')"|   Region  ", (" Spin ", i_spin, i_spin=1,n_spin,1)
      call localorb_info(info_str,use_unit,'(A)')
  
      do i_region = 1, n_active_regions, 1
  
      write(info_str,'(2X,A,I6,2X,F15.10,A,2X,F15.10,A)')"|   ", constraint_region_label(i_region), &
          ( constraint_potential(i_region,i_spin)*hartree, " eV", i_spin=1,n_spin,1)
      call localorb_info(info_str,use_unit,'(A)')
      enddo
  
      write(info_str,'(2X,A)') "| "
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A)') "| Electron count in different subsystems:"
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A,10X,A,I2,13X,A,I2)')"|   Region   ", ("Spin ", i_spin, i_spin=1,n_spin,1)
      call localorb_info(info_str,use_unit,'(A)')
  
      do i_region = 1, n_active_regions, 1
  
      write(info_str,'(2X,A,I6,5X,F15.10,5X,F15.10)')"|   ", constraint_region_label(i_region), &
          (electrons_in_region(i_region,i_spin), i_spin=1,n_spin,1)
      call localorb_info(info_str,use_unit,'(A)')
  
      enddo
  
      write(info_str,'(2X,A)') "| "
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A,A)')  "| Final average constraint potential:"
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A,F22.10,A,F22.10,A)') "|   ", avg_potential, " Ha   ", avg_potential*hartree, " eV."
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A)') "| "
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A,A)')  "| Overall Fermi level:"
      call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,'(2X,A,F22.10,A,F22.10,A)')"|   ", chemical_potential, " Ha   ", chemical_potential*hartree, " eV."
      call localorb_info(info_str,use_unit,'(A)')
  
      !          write(info_str,'(2X,A)') "| "
      !              call localorb_info(info_str,use_unit,'(A)')
      !          write(info_str,'(2X,A,A)') 
      !     +    "| Energy correction due to ",
      !     +    "auxiliary potentials:"
      !              call localorb_info(info_str,use_unit,'(A)')
      !          write(info_str,'(2X,A,F17.5,A,F17.5,A)')
      !     +    "|   ",
      !     +    constraint_energy_correction, " Ha   ",
      !     +    constraint_energy_correction*hartree, " eV."
      !              call localorb_info(info_str,use_unit,'(A)')
      !          write(info_str,'(2X,A)') "|"
      !              call localorb_info(info_str,use_unit,'(A)')
  
      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)')
  
      if (mixer_constraint.eq.1) then
        call cleanup_pulay_constraint ( )
      end if
  
      do i_spin = 1, n_spin, 1
      do i_region = 1, n_active_regions, 1
        do i_state = 1, n_states, 1
        constraint_proj_final(i_state, i_region, i_spin) = &
            constraint_proj(i_state, i_region, i_spin)
        end do
      end do
      end do
  
      ! end if relates to everything to do with constrained DFT. 
    end if
  
  else ! n_periodic != 0, have to project corehole state for periodic case
! in principle, also constraint for periodic case goes here 
! at first, only consider complex eigenvectors

! here comes everything without constraint, constraint follows later
! new: fixed_spin_moment
! do it like in constraint case using check_norm_v2 type of routine

    !call sync_integer_vector(force_occ_pr_state_periodic,n_force_occ*n_k_points)

    force_occ_state_periodic = 0
    i_k = 0
    do i_k_point = 1, n_k_points,1
      if (myid.eq. MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k = i_k + 1
        if(i_k == 1)then
            if (real_eigenvectors) then 
              allocate(overlap_matrix_w_complex(n_basis*(n_basis+1)/2),stat=info)
              call check_allocation(info, 'overlap_matrix_w_complex      ')
              allocate(overlap_matrix_w(n_basis*(n_basis+1)/2),stat=info)
              call check_allocation(info, 'overlap_matrix_w_complex      ')
              allocate(work_ovl(1,1),stat=info)
              call check_allocation(info, 'work_ovl      ')
            else
              allocate(overlap_matrix_w_complex(n_basis*(n_basis+1)/2),stat=info)
              call check_allocation(info, 'overlap_matrix_w_complex      ')
              allocate(overlap_matrix_w(1),stat=info)
              call check_allocation(info, 'overlap_matrix_w_complex      ')
              allocate(work_ovl(1,1),stat=info)
              call check_allocation(info, 'work_ovl      ')
            end if
        end if
!       here comes the construct_ovl for k_points
!       but I do not know yet how to synchronize the mpi-distributed 
!       results for different k_points to one occ_numbers value
        call construct_overlap(overlap_matrix,overlap_matrix_w,&
          overlap_matrix_w_complex,i_k_point,work_ovl)
!         overlap matrix exists for k_point, but what about about occ_numbers -- synchronising?

        call determine_corehole_projection_periodic &
          (overlap_matrix_w_complex,overlap_matrix_w,KS_eigenvector_complex,KS_eigenvector,i_k,i_k_point)
      end if
    end do

!     write(use_unit,*)"HALLO"
!     do i_k_point = 1, n_k_points,1
!       if (myid .eq.0) then
!         write(use_unit,*) i_k_point,force_occ_state_periodic(:,i_k_point)
!       end if
!     end do
! 
!     write(use_unit,*)"HALLO 2"
!     do i_k_point = 1, n_k_points,1
!       if (myid .eq.1) then
!         write(use_unit,*) i_k_point,force_occ_state_periodic(:,i_k_point)
!       end if
!     end do

    call sync_integer_vector(force_occ_state_periodic,n_force_occ*n_k_points)

!     if (myid .eq. 0) then
!       write(use_unit,*) "force_occ_state_periodic: cpu 0"
!       do i_k_point = 1, n_k_points, 1
!         write(use_unit,*) force_occ_state_periodic(1,i_k_point),i_k_point
!       end do
!       do i_state = 1, n_states, 1
!         write(use_unit,*) "occ_numbers 1st k_point before get_occ"
!         write(use_unit,*) occ_numbers(i_state,1,1),occ_numbers(i_state,2,1),i_state
!       end do
!     end if
! 
!     if (myid .eq. 1) then
!       write(use_unit,*) "force_occ_state_periodic : cpu 1"
!       do i_k_point = 1, n_k_points, 1
!         write(use_unit,*) force_occ_state_periodic(1,i_k_point),i_k_point
!       end do
!       do i_state = 1, n_states, 1
!         write(use_unit,*) "occ_numbers 1st k_point before get_occ"
!         write(use_unit,*) occ_numbers(i_state,1,1),occ_numbers(i_state,2,1),i_state
!       end do
!     end if

    if (.not. fixed_spin_moment) then  
  
      call get_occupation_numbers_occ_p0  &
        ( KS_eigenvalue, n_electrons, t_out,  &
        occ_numbers, chemical_potential &
        )
  
  !     if (myid .eq. 0) then
  !       do i_state = 1, n_states, 1
  !         write(use_unit,*) "occ_numbers 1st k_point after get_occ"
  !         write(use_unit,*) occ_numbers(i_state,1,1),occ_numbers(i_state,2,1),i_state
  !       end do
  !     end if
  
  !     write(use_unit,*) force_occ_state_periodic(:,1)
  !     works ;)
    else ! fixed_spin_moment

      do i_spin = 1, n_spin, 1

        call get_occupation_numbers_occ_periodic_v2(KS_eigenvalue(:,i_spin,:),&
          fixed_spin_moment_electrons(i_spin),.true.,occ_numbers(:,i_spin,:),&
          chemical_potential_spin(i_spin),i_spin)

      end do

! now calculate overall chemical potential
      diff_in_chem_pot = chemical_potential_spin(1) - chemical_potential_spin(2)
      if (.true.) then
!   call localorb_info('',use_unit,'(A)') 
    write (info_str, '(2X,A,A,F15.10,A)') &
          'Difference in chemical potential between ', &
          'spin channels:', diff_in_chem_pot*hartree,  &
          " eV."
    call localorb_info(info_str,use_unit,'(A)')
      end if
 
      spin_shift_fsm (1) = fixed_spin_moment_electrons(2)/n_electrons *  &
        diff_in_chem_pot
      spin_shift_fsm (2) = - fixed_spin_moment_electrons(1)/n_electrons *  &
        diff_in_chem_pot

      chemical_potential_spin(1) =  &
      chemical_potential_spin(1) - spin_shift_fsm (1)
      chemical_potential_spin(2) =  &
      chemical_potential_spin(2) - spin_shift_fsm (2)

      write (info_str, '(2X,A,F15.10,A)') &
          'Overall chemical potential: ', &
           chemical_potential_spin(1)*hartree,  &
          " eV."
    call localorb_info(info_str,use_unit,'(A)')
    call localorb_info('',use_unit,'(A)') 


!       if (extra_forced_occ_debug_out)
!         temp_n_electrons_forced_state = 0.d0
! 
!     do i_state = 1, n_states, 1
!       do i_k_point = 1, n_k_points, 1
!     
!         temp_n_electrons_forced_state(i_state,1) = temp_n_electrons_forced_state(i_state,1) + &
!           occ_numbers(i_state,1,i_k_point) * &
!           k_weights(i_k_point)
!     
!         temp_n_electrons_forced_state(i_state,2) = temp_n_electrons_forced_state(i_state,2) + &
!           occ_numbers(i_state,2,i_k_point) * &
!           k_weights(i_k_point)
!     
!       end do
!     end do
! 
! 
!       if (myid .eq. 0) then
!     write(use_unit,*)"ELECTRONS IN STATES SUMMED OVER K-POINTS"
!     do i_state = 1, n_states, 1
!         write(use_unit,*) i_state, temp_n_electrons_forced_state(i_state,1),temp_n_electrons_forced_state(i_state,2)
!       end do
!       end if



    end if ! fixed_spin_moment

  end if ! n_periodic == 0


end subroutine get_KS_orbitals_forced_occ_p0
!******    


!*****
!-------------------------------------------------------------------------------
!****s* force_occupation/adjust_force_occ_sr
!  NAME
!   adjust_force_occ_sr
!  SYNOPSIS

subroutine adjust_force_occ_sr(sr_to_soc_idxmap, i_k_point)

!  PURPOSE
!  In the case of spin=none and perturbative SOC we need to adjust the forced
!  indices to the new dimensions. The constraint will be applied to i_spin 1
!  and the remainder distributed to i_spin 2
!  USES

  use localorb_io,                 only : localorb_info
  use mpi_utilities,               only : myid
  use synchronize_mpi,             only : sync_force_occupation
  use mpi_tasks,                   only : n_tasks
  use dimensions,                  only : n_states, n_k_points
  use soc_utilities,               only : revert_soc_to_sr_environment, &
                                          convert_sr_to_soc_environment

  implicit none

!  ARGUMENTS


! INPUTS
!
!  sr_to_soc_idx - indexing array mapping sr-ks_indices to soc-indices
!  i_k_point - (optional) k-point index
!
! OUTPUT
!
!  NONE - just reindexing of module variables
!
!  AUTHOR
!    Georg Michelitsch, TU Muenchen
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
!    Feb 2018 - created.
!  SOURCE


  ! local variables

  character*120                                 :: info_str

  integer                                       :: n_force_occ_soc
  integer                                       :: n_spin_sr
  integer,dimension(:),allocatable              :: force_occ_pr_state_soc
  integer,dimension(:),allocatable              :: force_occ_spin_soc
  real*8,dimension(:),allocatable               :: forced_occ_number_soc
  integer,dimension(:),allocatable              :: force_occ_min_KS_state_soc
  integer,dimension(:),allocatable              :: force_occ_max_KS_state_soc

  ! counters

  integer                                       :: i_force_occ, i_force_occ_soc, &
                                                   i_state, i_spin, shift

  ! variables

  integer,dimension(n_states)                   :: sr_to_soc_idxmap
  integer                                       :: i_k_point

  ! begin work

  call revert_soc_to_sr_environment ()
    n_spin_sr = n_spin
  call convert_sr_to_soc_environment ()

  if (myid .eq. 0) then
    write(info_str,'(2X,A)') ''
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A,A)') &
         "Updating force_occupation constraints", &
         " for correct report in perturbative SOC run."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if


  if (n_spin .eq. n_spin_sr) then

    n_force_occ_soc = n_force_occ * 2

    !write(use_unit,*) n_force_occ_soc, myid

!    if (myid .eq. 0) then
    if (myid .eq. MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then

      allocate(force_occ_pr_state_soc(n_force_occ_soc))
      allocate(force_occ_spin_soc(n_force_occ_soc))
      allocate(forced_occ_number_soc(n_force_occ_soc))
      allocate(force_occ_min_KS_state_soc(n_force_occ_soc))
      allocate(force_occ_max_KS_state_soc(n_force_occ_soc))

      preloop : do i_force_occ = 1, n_force_occ
        socloop : do i_force_occ_soc = 2 * i_force_occ - 1, 2 * i_force_occ

          if (n_periodic .eq. 0) then
            force_occ_pr_state_soc(i_force_occ_soc) = 2 * force_occ_pr_state(i_force_occ) - (2 * i_force_occ - i_force_occ_soc)
          else
            force_occ_pr_state_soc(i_force_occ_soc) = 2 * force_occ_pr_state_periodic(i_force_occ, i_k_point) - (2 * i_force_occ - i_force_occ_soc)
          end if

          force_occ_spin_soc(i_force_occ_soc) = force_occ_spin(i_force_occ)

          force_occ_min_KS_state_soc(i_force_occ_soc) = 2 * force_occ_min_KS_state(i_force_occ) - 1
          force_occ_max_KS_state_soc(i_force_occ_soc) = 2 * force_occ_max_KS_state(i_force_occ)

        end do socloop

        ! Decrease by 2 since loop variable is at (bounds + 1) after socloop
        if (forced_occ_number(i_force_occ) .lt. 1.0d0) then
          forced_occ_number_soc(i_force_occ_soc - 2) = forced_occ_number(i_force_occ)
          forced_occ_number_soc(i_force_occ_soc - 1) = 0.d0
        else
          forced_occ_number_soc(i_force_occ_soc - 2) = 1.0d0
          forced_occ_number_soc(i_force_occ_soc - 1) = forced_occ_number(i_force_occ) - 1.0d0
        end if

      end do preloop
    end if ! each MPI task / kpt

    if (n_periodic .eq. 0) then
      deallocate(force_occ_pr_state)
      deallocate(force_occ_state)
      allocate(force_occ_pr_state(n_force_occ_soc))
      allocate(force_occ_state(n_force_occ_soc))
    else
      deallocate(force_occ_pr_state_periodic)
      deallocate(force_occ_state_periodic)
      allocate(force_occ_pr_state_periodic(n_force_occ_soc, n_k_points))
      force_occ_pr_state_periodic = 0
      allocate(force_occ_state_periodic(n_force_occ_soc, n_k_points))
      force_occ_state_periodic = 0
    end if

    deallocate(force_occ_spin)
    deallocate(forced_occ_number)
    deallocate(force_occ_min_KS_state)
    deallocate(force_occ_max_KS_state)

    allocate(force_occ_spin(n_force_occ_soc))
    allocate(forced_occ_number(n_force_occ_soc))
    allocate(force_occ_min_KS_state(n_force_occ_soc))
    allocate(force_occ_max_KS_state(n_force_occ_soc))

    n_force_occ = n_force_occ_soc

!    if (myid .eq. 0) then
    if (myid .eq. MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
      if (n_periodic .eq. 0) then
        force_occ_pr_state = force_occ_pr_state_soc
        force_occ_state = force_occ_pr_state_soc
      else
        force_occ_pr_state_periodic(:, i_k_point) = force_occ_pr_state_soc
        force_occ_state_periodic(:, i_k_point) = force_occ_pr_state_soc
      end if

      force_occ_spin = force_occ_spin_soc
      forced_occ_number = forced_occ_number_soc
      force_occ_min_KS_state = force_occ_min_KS_state_soc
      force_occ_max_KS_state = force_occ_max_KS_state_soc

      deallocate(force_occ_pr_state_soc)
      deallocate(force_occ_spin_soc)
      deallocate(forced_occ_number_soc)
      deallocate(force_occ_min_KS_state_soc)
      deallocate(force_occ_max_KS_state_soc)
    end if ! each MPI task / kpt

  else ! adapt spin collinear force_occupations

    do i_force_occ = 1, n_force_occ
      do i_spin = 1, n_spin_sr
        if (i_spin == 2) then
          shift = n_states / 2
        else
          shift = 0
        end if
        do i_state = 1, n_states
          if (force_occ_spin(i_force_occ) .eq. i_spin) then
            if (n_periodic .eq. 0) then
              if (force_occ_state(i_force_occ) .eq. sr_to_soc_idxmap(i_state) - shift) then
                force_occ_state(i_force_occ) = i_state
                force_occ_pr_state(i_force_occ) = i_state
                force_occ_spin(i_force_occ) = 1
              end if
            else ! periodic case
              if (force_occ_state_periodic(i_force_occ, i_k_point) .eq. sr_to_soc_idxmap(i_state) - shift) then
                force_occ_state_periodic(i_force_occ, i_k_point) = i_state
                force_occ_pr_state_periodic(i_force_occ, i_k_point) = i_state
                force_occ_spin(i_force_occ) = 1
              end if
            end if
            if (force_occ_min_KS_state(i_force_occ) .eq. sr_to_soc_idxmap(i_state) - shift) then
              force_occ_min_KS_state(i_force_occ) = i_state
            end if
            if (force_occ_max_KS_state(i_force_occ) .eq. sr_to_soc_idxmap(i_state) - shift) then
              force_occ_max_KS_state(i_force_occ) = i_state
            end if
          end if
        end do
      end do
    end do

  end if

  if (n_periodic .eq. 0) then
    call sync_force_occupation(force_occ_pr_state, force_occ_state, &
            force_occ_spin, forced_occ_number, force_occ_min_KS_state, &
            force_occ_max_KS_state, n_force_occ, i_k_point)
  else
    call sync_force_occupation(force_occ_pr_state_periodic(:, i_k_point), &
            force_occ_state_periodic(:, i_k_point), force_occ_spin, &
            forced_occ_number, force_occ_min_KS_state, force_occ_max_KS_state, &
            n_force_occ, i_k_point)
  end if

  if (myid .eq. 0) then
!  if (myid .eq. MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
    if (n_periodic .eq. 0) then
      do i_force_occ = 1, n_force_occ
        write (use_unit,'(2X,A,I3,1X,I8,1X,F12.8,1X,I8,1X,I8)') &
        "| Updated projector to ", force_occ_state(i_force_occ), force_occ_spin(i_force_occ), &
        forced_occ_number(i_force_occ), force_occ_min_KS_state(i_force_occ), force_occ_max_KS_state(i_force_occ)
      end do
    else
      do i_force_occ = 1, n_force_occ
        write (use_unit,'(2X,A,I3,1X,I8,1X,F12.8,1X,I8,1X,I8)') &
        "| Updated projector to ", force_occ_state_periodic(i_force_occ, i_k_point), force_occ_spin(i_force_occ), &
        forced_occ_number(i_force_occ), force_occ_min_KS_state(i_force_occ), force_occ_max_KS_state(i_force_occ)
      end do
    end if
  end if ! each MPI task / kpt

end subroutine adjust_force_occ_sr
!******


!*****
!-------------------------------------------------------------------------------
!****s* force_occupation/force_occ_reduce_subspace
!  NAME
!   force_occ_reduce_subspace
!  SYNOPSIS

subroutine force_occ_reduce_subspace(n_scf, i_force_occ, reduction)

!  PURPOSE
!  This subroutine reduces the MOM-subspace in the force_occupation_projector
!  case to overcome oscillations in the SCF, when there are two symmetry-equivalent
!  KS states to be constrained
!  USES

  use localorb_io,                 only : localorb_info
  use mpi_utilities,               only : myid
  use synchronize_mpi,             only : sync_force_occupation
  use mpi_tasks,                   only : n_tasks
  use dimensions,                  only : n_force_occ

  implicit none

!  ARGUMENTS


! INPUTS
!
!  n_scf - number of scf-cycles
!
! OUTPUT
!
!  NONE - just modifying module variables
!
!  AUTHOR
!    Georg Michelitsch, TU Muenchen
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
!    Apr 2018 - created.
!  SOURCE


  ! local variables

  character*120                                 :: info_str

  ! counters

  integer, intent(IN)                           :: i_force_occ
  integer, intent(IN)                           :: n_scf
  integer, intent(IN)                           :: reduction

  ! variables

  ! begin work

  !do i_force_occ=1, n_force_occ
  if ((force_occ_max_KS_state(i_force_occ) .gt. force_occ_maxred_KS_state(i_force_occ)) .and. &
  (n_scf .gt. force_occ_step_KS_state(i_force_occ))) then
    force_occ_max_KS_state(i_force_occ) = force_occ_max_KS_state(i_force_occ) - reduction
    if (myid.eq.0) then
      write (use_unit,'(2X,A,I5,A,I3,A,I3)') &
      "| SCF [ ", n_scf , " ]  Updated projector ", i_force_occ, " to new upper boundary  ", force_occ_max_KS_state(i_force_occ)
    end if
  end if
  !end do

end subroutine force_occ_reduce_subspace
!******

end module force_occupation

