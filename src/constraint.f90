!****h* FHI-aims/constraint
!  NAME
!    constraint
!  SYNOPSIS

module constraint

  !  PURPOSE
  !  Module constraint handles a possible basis-function dependent constraint on occupation numbers.
  !
  !  This module will be used in subroutine get_KS_orbitals if needed.
  !
  !  Subroutines inside the module:
  !  o allocate_constraint
  !  o add_constraint_potentials
  !  o cleanup_constraint
  !
  ! USES

  use dimensions
  use runtime_choices
  use localorb_io

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





  real*8, dimension(:,:), allocatable :: constraint_potential

  integer, dimension(:), allocatable :: constraint_region
  real*8, dimension(:,:), allocatable :: constraint_electrons
  integer, dimension(:), allocatable :: constraint_region_label
  integer, dimension(:), allocatable :: constraint_region_number
  real*8, dimension(:,:,:), allocatable :: constraint_proj_final

  integer :: n_active_regions

  ! for mixing of constraint potentials
  real*8, dimension(:, :), allocatable :: constraint_pot_old
  real*8, dimension(:, :), allocatable :: delta_constraint_potential
  real*8, dimension(:, :), allocatable :: delta_potential_mixing

  ! energy correction for sum of eigenvalues due to artificial constraint potentials
  real*8, private :: constraint_energy_correction

  logical :: constraint_debug = .false.

  ! work Hamiltonian for later use
  real*8, dimension(:,:), allocatable :: hamiltonian_work

  ! work overlap matrix and ovlp * eigenvector matrix to
  ! speed up evaluation of constraint_projectors

  real*8, dimension(:,:), allocatable :: full_ovlp
  real*8, dimension(:,:), allocatable :: ovlp_times_eigenvec

  real*8, dimension(:), allocatable :: spin_shift

  ! mulliken_decomp_total was moved to mulliken module

  character*100, private :: info_str

  !******
contains
  !---------------------------------------------------------------------


  !****s* constraint/allocate_constraint
  !  NAME
  !    allocate_constraint
  !  SYNOPSIS

  subroutine allocate_constraint( )

    !  PURPOSE
    !  Allocates module variables
    !
    !  ARGUMENTS
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE




    use mpi_tasks, only: check_allocation
    implicit none
    integer:: info


    !       begin work

    !        write(use_unit,*) "allocate constraint"

    if (.not.allocated(constraint_potential)) then
       allocate( constraint_potential(n_region, n_spin),stat=info)
       call check_allocation(info, 'constraint_potential          ')
    end if

    if (.not.allocated (constraint_region)) then
       allocate ( constraint_region(n_atoms) ,stat=info)
       call check_allocation(info, 'constraint_region             ')
    end if
    if (.not.allocated(constraint_electrons)) then
       allocate ( constraint_electrons(n_region,2) ,stat=info)
       call check_allocation(info, 'constraint_electrons          ')
    end if
    if (.not.allocated(constraint_region_label)) then
       allocate ( constraint_region_label(n_region) ,stat=info)
       call check_allocation(info, 'constraint_region_label       ')
    end if
    if (.not.allocated(constraint_region_number)) then
       ! The dimension must be n_atoms because the region label,
       ! in principle, is allowed to reach up to n_atoms.
       allocate ( constraint_region_number(n_atoms) ,stat=info)
       call check_allocation(info, 'constraint_region_number      ')
    end if

    if (.not.allocated(constraint_pot_old)) then
       allocate( constraint_pot_old (n_region, n_spin),stat=info)
       call check_allocation(info, 'constraint_pot_old            ')
    end if
    if (.not.allocated(delta_constraint_potential)) then
       allocate( delta_constraint_potential (n_region, n_spin))
       call check_allocation(info, 'delta_constraint_potential    ')
    end if
    !              if (mixer_constraint.eq.1) then
    if (.not.allocated(delta_potential_mixing)) then
       allocate( delta_potential_mixing(n_region, n_spin) ,stat=info)
       call check_allocation(info, 'delta_potential_mixing        ')
    end if
    if (.not.allocated(spin_shift)) then
       allocate( spin_shift(n_spin),stat=info)
       call check_allocation(info, 'spin_shift                    ')
    end if
    !              end if


  end subroutine allocate_constraint
  !******
  !---------------------------------------------------------------------
  !****s* constraint/allocate_constraint_projectors
  !  NAME
  !    allocate_constraint_projectors
  !  SYNOPSIS

  subroutine allocate_constraint_projectors()

    !  PURPOSE
    !  Subroutine allocate_pulay allocates the necessary storage arrays
    !  for constraint.  
    !  USES
    use mpi_tasks, only: check_allocation
    implicit none
    !  ARGUMENTS
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE




    integer:: info

    if (.not.allocated(constraint_proj_final)) then
       allocate( constraint_proj_final(n_states, n_region, n_spin),stat=info)
       call check_allocation(info, 'constraint_proj_final         ')
    end if

  end subroutine allocate_constraint_projectors
  !******
  !------------------------------------------------------------------------------
  !****s* constraint/allocate_hamiltonian_work
  !  NAME
  !    allocate_hamiltonian_work
  !  SYNOPSIS

  subroutine allocate_hamiltonian_work ( overlap_matrix  )

    !  PURPOSE
    !  ???????????
    !
    use mpi_tasks, only: check_allocation
    implicit none
    !  ARGUMENTS

    real*8, dimension(n_basis*(n_basis+1)/2) :: overlap_matrix

    !  INPUTS
    !   o overlap_matrix -- overlap matrix
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    integer info
    integer :: i_basis_1,i_basis_2,i_index


    !       begin work

    if (.not.allocated(hamiltonian_work)) then
       allocate( hamiltonian_work(n_basis*(n_basis+1)/2,n_spin) ,stat=info)
       call check_allocation(info, 'hamiltonian_work              ')
    end if

    if (.not.allocated(full_ovlp)) then
       allocate(full_ovlp(n_basis,n_basis),stat=info)
       call check_allocation(info, 'full_ovlp                     ')
    end if

    if (.not.allocated(ovlp_times_eigenvec)) then
       allocate(ovlp_times_eigenvec(n_basis,n_states),stat=info)
       call check_allocation(info, 'ovlp_times_eigenvec           ')
    end if

    ! expand packed overlap matrix
    ! FIXME - THIS CAN GO AWAY WHEN WE EXPAND THE OVLP / HAM MATRIX TO
    ! THEIR FULLL SIZES IN THE FIRST PLACE
    i_index = 0
    do i_basis_2 = 1, n_basis, 1
       do i_basis_1 = 1, i_basis_2, 1
          i_index = i_index+1
          full_ovlp(i_basis_1,i_basis_2) = overlap_matrix(i_index)
          full_ovlp(i_basis_2,i_basis_1) = overlap_matrix(i_index)
       enddo
    enddo

  end subroutine allocate_hamiltonian_work

  !******
  !---------------------------------------------------------------------
  !****s* constraint/add_constraint_potentials
  !  NAME
  !    add_constraint_potentials
  !  SYNOPSIS

  subroutine add_constraint_potentials ( overlap_matrix, hamiltonian  )

    !  PURPOSE
    !    Subroutine add_constraint_potentials adds auxiliary potentials to
    !    each region in order to enforce a possible constraint  

    ! USES
    use basis
    implicit none
    !  ARGUMENTS

    real*8 overlap_matrix ( n_basis*(n_basis+1)/2 )
    real*8 hamiltonian ( n_basis*(n_basis+1)/2, n_spin )

    !  INPUTS
    !   o overlap_matrix - overlap matrix
    !   o hamiltonian  - ???????????????
    !
    !  OUTPUT
    !   o hamiltonian  - ???????????????
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer :: i_spin
    integer :: i_basis_1
    integer :: i_basis_2

    integer :: i_hamiltonian, info

    !  begin work

    do i_spin = 1, n_spin, 1

       i_hamiltonian = 0
       do i_basis_2 = 1, n_basis, 1
          do i_basis_1 = 1, i_basis_2, 1
             i_hamiltonian = i_hamiltonian + 1

             hamiltonian(i_hamiltonian, i_spin) = &
                  hamiltonian(i_hamiltonian, i_spin) + &
                  0.5d0 * &
                  (   constraint_potential &
                  (constraint_region(basis_atom(i_basis_1)), i_spin ) &
                  + constraint_potential &
                  (constraint_region(basis_atom(i_basis_2)), i_spin ) &
                  ) * overlap_matrix(i_hamiltonian)

          enddo
       enddo

    enddo

  end subroutine add_constraint_potentials
  !******
  !---------------------------------------------------------------------
  !****s* constraint/add_constraint_potentials_v2
  !  NAME
  !    add_constraint_potentials_v2
  !  SYNOPSIS

  subroutine add_constraint_potentials_v2 &
       ( current_spin, overlap_matrix, hamiltonian  )


    !  PURPOSE
    !  Subroutine add_constraint_potentials adds auxiliary potentials to
    !  each region in order to enforce a possible constraint
    !  Here, only one spin component is affected at a time    

    !  USES

    use basis
    implicit none
    !  ARGUMENTS

    real*8 overlap_matrix ( n_basis*(n_basis+1)/2 )
    real*8 hamiltonian ( n_basis*(n_basis+1)/2, n_spin )

    !  INPUTS
    !    o overlap_matrix - overlap matrix
    !    o hamiltonian - ????????
    !  OUTPUT
    !    o hamiltonian - ????????
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE






    integer :: current_spin

    ! counters

    integer :: i_basis_1
    integer :: i_basis_2

    integer :: i_hamiltonian

    !  begin work

    i_hamiltonian = 0
    do i_basis_2 = 1, n_basis, 1
       do i_basis_1 = 1, i_basis_2, 1
          i_hamiltonian = i_hamiltonian + 1

          hamiltonian(i_hamiltonian, current_spin) = &
               hamiltonian(i_hamiltonian, current_spin) + &
               0.5d0 * &
               (   constraint_potential &
               (constraint_region(basis_atom(i_basis_1)), &
               current_spin ) &
               + constraint_potential &
               (constraint_region(basis_atom(i_basis_2)), &
               current_spin ) &
               ) * overlap_matrix(i_hamiltonian)

       enddo
    enddo

  end subroutine add_constraint_potentials_v2
  !******
  !---------------------------------------------------------------------
  !****s* constraint/get_electrons_per_region
  !  NAME
  !    get_electrons_per_region
  !  SYNOPSIS

  subroutine get_electrons_per_region &
       ( KS_eigenvector, overlap_matrix, occ_numbers, &
       constraint_proj,  electrons_in_region )

    !  PURPOSE
    !    Subroutine get_electrons_per_region computes the projection of the
    !    current Kohn-Sham states onto the different subsystems, and from there
    !    evaluates the total number of electrons found in each subsystem  


    !  USES
    use basis
    implicit none
    !  ARGUMENTS

    real*8, dimension(n_basis, n_states, n_spin)  ::  KS_eigenvector
    real*8, dimension(n_basis*(n_basis+1)/2 )     :: overlap_matrix( n_basis*(n_basis+1)/2 )
    real*8, dimension(n_states, n_spin)           :: occ_numbers
    real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
    real*8, dimension(n_region, n_spin)           :: electrons_in_region


    !  INPUTS
    !  o   KS_eigenvector - ??????????
    !  o   overlap_matrix - overlap matrix
    !  o   occ_numbers - ??????????
    !
    !  OUTPUT
    !
    ! o  constraint_proj  - ??????????
    ! o  electrons_in_region  - ??????????
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE






    !  Local variables

    real*8, external :: ddot
    ! counters

    integer :: i_spin
    integer :: i_state
    integer :: i_region
    integer :: i_basis_1
    integer :: i_basis_2

    integer :: i_hamiltonian

    !  begin work

    ! First, compute projection operators onto subsystems, as defined in
    ! Eqs. (7a) and (7b) of Behler et al. PRB
    constraint_proj = 0.d0
    do i_spin = 1, n_spin, 1
       do i_state = 1, n_states, 1

          !test
          !              if ((i_spin.eq.1).and.(i_state.eq.4)) then
          !                write(use_unit,'(A)')
          !     +          "Spin 1, state 4: Assembling constraint_proj:"
          !              end if
          !test end
          i_hamiltonian = 0
          do i_basis_2 = 1, n_basis, 1
             do i_basis_1 = 1, i_basis_2, 1
                i_hamiltonian = i_hamiltonian+1

                constraint_proj (i_state, &
                     constraint_region(basis_atom(i_basis_1)), i_spin ) = &
                     constraint_proj (i_state, &
                     constraint_region(basis_atom(i_basis_1)), i_spin) + &
                     KS_eigenvector(i_basis_1, i_state, i_spin) * &
                     overlap_matrix(i_hamiltonian) * &
                     KS_eigenvector(i_basis_2, i_state, i_spin)

                if (i_basis_1.ne.i_basis_2) then
                   constraint_proj (i_state, &
                        constraint_region(basis_atom(i_basis_2)), i_spin ) = &
                        constraint_proj (i_state, &
                        constraint_region(basis_atom(i_basis_2)), i_spin ) + &
                        KS_eigenvector(i_basis_2, i_state, i_spin) * &
                        overlap_matrix(i_hamiltonian) * &
                        KS_eigenvector(i_basis_1, i_state, i_spin)
                end if

                !test
                !              if ((i_spin.eq.1).and.(i_state.eq.4)) then
                !                  write(use_unit,'(2X,I5,1X,I5,A,F10.5,A,F10.5)')
                !     +            i_basis_1, i_basis_2, ": constraint_proj(1) =",
                !     +            constraint_proj(i_state,1,i_spin),
                !     +            ", constraint_proj(2) =",
                !     +            constraint_proj(i_state,2,i_spin)
                !              end if
                !test end

             enddo
          enddo

       enddo
    enddo

    !test
    !        write(use_unit,*)
    !test end

    ! Next, for each Kohn-Sham eigenstate, compute the contribution of each subsystem
    ! to the overall occupation numbers, and sum up the total number of electrons
    electrons_in_region = 0.d0
    do i_spin = 1, n_spin, 1
       !test
       !          write(use_unit,'(2X,A,I3)') "Spin ", i_spin
       !test end
       do i_region = 1, n_active_regions, 1
          !test
          !            write(use_unit,'(4X,A,I3)') "Region ", i_region
          !test end
          !            do i_state = 1, n_states, 1
          !
          !              constraint_occupation(i_state, i_region, i_spin) =
          !     +          constraint_proj(i_state,i_region,i_spin) *
          !     +          occ_numbers(i_state,i_spin)

          electrons_in_region(i_region,i_spin)= &
               ddot(n_states,constraint_proj(:,i_region,i_spin),1, &
               occ_numbers(:,i_spin),1)

          !test
          !            if ((i_region.eq.1).and.(i_spin.eq.1)) then
          !              write(use_unit,'(6X,A,I5,A,F10.5,A,F10.5)')
          !     +        "State ", i_state, ": Occupation ",
          !     +        constraint_occupation(i_state, i_region, i_spin),
          !     +        ", total electron count ",
          !     +        electrons_in_region(i_region,i_spin)
          !              write(use_unit,'(8X,A,F10.5,A,F10.5)')
          !     +        "constraint_proj: ",
          !     +        constraint_proj(i_state,i_region,i_spin),
          !     +        ", occupation: ", occ_numbers(i_state,i_spin)
          !            end if
          !test end
          !            enddo
       enddo
    enddo

    if (constraint_debug) then
       write(info_str,*)
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A)') &
            "Electron count in different subsystems:"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,10X,A,I2,10X,A,I5)') &
            "| Region  ", ("Spin ", i_spin, i_spin=1,n_spin,1)
       call localorb_info(info_str,use_unit,'(A)')
       do i_region = 1, n_active_regions, 1
          write(info_str,'(2X,A,I6,2X,F15.10,2X,F15.10)') &
               "| ", &
               constraint_region_label(i_region), &
               (electrons_in_region(i_region,i_spin), i_spin=1,n_spin,1)
          call localorb_info(info_str,use_unit,'(A)')
       enddo
       write(info_str,*)
       call localorb_info(info_str,use_unit,'(A)')
    end if

  end subroutine get_electrons_per_region
  !******
  !---------------------------------------------------------------------
  !****s* constraint/get_electrons_per_region_v2
  !  NAME
  !    get_electrons_per_region_v2
  !  SYNOPSIS

  subroutine get_electrons_per_region_v2 &
       ( current_spin, &
       KS_eigenvector, overlap_matrix, occ_numbers, &
       constraint_proj, &
       electrons_in_region &
       )

    !  PURPOSE
    !
    !  Subroutine get_electrons_per_region computes the projection of the
    !  current Kohn-Sham states onto the different subsystems, and from there
    !  evaluates the total number of electrons found in each subsystem
    !
    !  This version only treats one spin channel at a time, current_spin
    !   

    !  USES

    use basis
    implicit none

    !  ARGUMENTS

    integer :: current_spin
    real*8, dimension(n_basis, n_states, n_spin) :: KS_eigenvector
    real*8, dimension( n_basis*(n_basis+1)/2 )   :: overlap_matrix
    real*8, dimension(n_states, n_spin)          :: occ_numbers
    real*8, dimension(n_states, n_region, n_spin):: constraint_proj
    real*8, dimension(n_region, n_spin)          :: electrons_in_region

    !  INPUTS
    !    o KS_eigenvector - ?????
    !    o overlap_matrix - overlap matrix
    !    o occ_numbers - ?????
    !  OUTPUT
    !    o constraint_proj - ?????
    !    o electrons_in_region - ?????
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE




    integer :: max_occ_number

    real*8, external :: ddot

    ! counters

    integer :: i_state
    integer :: i_region
    integer :: i_basis_1
    integer :: i_basis_2

    integer :: i_hamiltonian

    !  begin work

    ! First, compute projection operators onto subsystems, as defined in
    ! Eqs. (7a) and (7b) of Behler et al. PRB
    constraint_proj(:,:,current_spin) = 0.d0

    ! General initialization of max_occ_number
    ! in some cases, we evaluate constraint_fermi at some point; cannot
    ! limit the number of states for which we evaluate constraint_proj
    max_occ_number = n_states

    ! if we do not need the constraint projector outside this subroutine,
    ! only evaluate it for occupied states
    if (mixer_constraint.eq.2) then
       do i_state = n_states, 1, -1
          if (dabs(occ_numbers(i_state,current_spin)).gt.0.d0) then
             max_occ_number = i_state
             exit
          end if
       enddo
    end if

    ! prepare the full product of overlap_matrix * KS_eigenvectors :

    call dgemm &
         ( 'N','N',n_basis,max_occ_number,n_basis, &
         1.d0,full_ovlp,n_basis, &
         KS_eigenvector(:,1:max_occ_number,current_spin),n_basis, &
         0.d0,ovlp_times_eigenvec,n_basis &
         )

    do i_state = 1, max_occ_number, 1

       do i_basis_1 = 1, n_basis, 1

          constraint_proj (i_state, &
               constraint_region(basis_atom(i_basis_1)), &
               current_spin ) = &
               constraint_proj (i_state, &
               constraint_region(basis_atom(i_basis_1)), &
               current_spin) + &
               KS_eigenvector(i_basis_1, i_state, current_spin)* &
               ovlp_times_eigenvec(i_basis_1,i_state)

       enddo

    enddo

    ! Next, for each Kohn-Sham eigenstate, compute the contribution of each subsystem
    ! to the overall occupation numbers, and sum up the total number of electrons
    ! electrons_in_region(:,current_spin) = 0.d0

    do i_region = 1, n_active_regions, 1

       electrons_in_region(i_region,current_spin)= &
            ddot( max_occ_number, &
            constraint_proj(1:max_occ_number,i_region,current_spin),1, &
            occ_numbers(1:max_occ_number,current_spin),1 )

    enddo

    if (constraint_debug) then
       write(info_str,*)
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A)') &
            "Electron count in different subsystems:"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,15X,A,I2)') &
            "| Region  ", "Spin ", current_spin
       call localorb_info(info_str,use_unit,'(A)')
       do i_region = 1, n_active_regions, 1
          write(info_str,'(2X,A,I6,2X,F20.15)') &
               "| ", &
               constraint_region_label(i_region), &
               electrons_in_region(i_region,current_spin)
          call localorb_info(info_str,use_unit,'(A)')
       enddo
       write(info_str,*)
       call localorb_info(info_str,use_unit,'(A)')
    end if

  end subroutine get_electrons_per_region_v2
  !******
  !---------------------------------------------------------------------
  !****s* constraint/get_energy_correction
  !  NAME
  !    get_energy_correction
  !  SYNOPSIS

  subroutine get_energy_correction &
       ( KS_eigenvector, occ_numbers, overlap_matrix &
       )

    !  PURPOSE
    !  Subroutine get_energy_correction computes the correction term to
    !  the sum of eigenvalues which was artificially added through the
    !  constraint Hamiltonian.    

    !  USES

    use basis
    implicit none

    !  ARGUMENTS
    real*8, dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
    real*8, dimension( n_basis*(n_basis+1)/2 )   ::  overlap_matrix
    real*8, dimension(n_states, n_spin)          :: occ_numbers

    !  INPUTS
    !    o  KS_eigenvector  - ???????
    !    o  overlap_matrix  - overlap matrix
    !    o  occ_numbers - ???????
    !
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    !  Local variables

    real*8 :: aux_term_2
    real*8, dimension( n_basis*(n_basis+1)/2,n_spin ) ::  aux_matrix

    integer, dimension(n_spin) :: max_occ_number

    ! counters

    integer :: i_spin
    integer :: i_state
    integer :: i_basis_1
    integer :: i_basis_2

    integer :: i_hamiltonian

    !  begin work

    ! find the maximal occupation number
    do i_spin = 1, n_spin, 1
       do i_state = n_states, 1, -1
          if (dabs(occ_numbers(i_state,i_spin)).gt.0.d0) then
             max_occ_number(i_spin) = i_state
             exit
          end if
       enddo
    enddo
    do i_spin = 1, n_spin, 1

       i_hamiltonian = 0
       do i_basis_2 = 1, n_basis, 1
          do i_basis_1 = 1, i_basis_2, 1
             i_hamiltonian= i_hamiltonian+1

             ! this is really the constrained Hamiltonian
             aux_matrix(i_hamiltonian,i_spin) = &
                  0.5d0 * &
                  (   constraint_potential &
                  (constraint_region(basis_atom(i_basis_1)), i_spin ) &
                  + constraint_potential &
                  (constraint_region(basis_atom(i_basis_2)), i_spin ) &
                  ) * overlap_matrix(i_hamiltonian)

          enddo
       enddo

    enddo

    constraint_energy_correction = 0.d0

    do i_spin = 1, n_spin, 1

       do i_state = 1, max_occ_number(i_spin), 1

          i_hamiltonian = 0
          do i_basis_2 = 1, n_basis, 1

             aux_term_2 = occ_numbers(i_state,i_spin) * &
                  KS_eigenvector(i_basis_2,i_state,i_spin)

             do i_basis_1 = 1, i_basis_2, 1
                i_hamiltonian= i_hamiltonian+1

                ! we are using packed matrices; in the matrix products, must account
                ! also for the lower triangle of the auxiliary Hamiltonian
                if (i_basis_2.eq.i_basis_1) then
                   constraint_energy_correction = &
                        constraint_energy_correction + &
                        aux_term_2 * &
                        KS_eigenvector(i_basis_1,i_state,i_spin) * &
                        aux_matrix(i_hamiltonian,i_spin)
                else
                   constraint_energy_correction = &
                        constraint_energy_correction + &
                        2.d0* aux_term_2 * &
                        KS_eigenvector(i_basis_1,i_state,i_spin) * &
                        aux_matrix(i_hamiltonian,i_spin)
                end if

             enddo
          enddo

       enddo

    enddo


  end subroutine get_energy_correction

  !---------------------------------------------------------------------
  !****s* constraint/get_energy_backshift
  !  NAME
  !    get_energy_backshift
  !  SYNOPSIS

  subroutine get_energy_backshift(KS_eigenvalue)

    !  PURPOSE
    !  Subroutine get_energy_backshift computes the correction terms 
    !  which needs to be added to KS energies, in order to make constraint
    !  spectra comparible to unconstraint caclulated ones    

    !  USES
    use mulliken, only : mulliken_decomp_total
!    use physics
    implicit none

    !  ARGUMENTS


    real*8, dimension(n_states,n_spin,n_k_points):: KS_eigenvalue

    !  INPUTS
    !    o  KS_eigenvector  - ???????
    !    o  occ_numbers - ???????
    !
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    !  Local variables


    ! counters

    integer :: i_spin
    integer :: i_state
    integer :: i_region
    integer :: i_atom
    integer :: i_k_point


    !  begin work

    do i_state=1,n_states,1
    do i_atom = 1, n_atoms
      i_region = constraint_region(i_atom)
    do i_spin = 1, n_spin
    do i_k_point = 1, n_k_points


      KS_eigenvalue(i_state,i_spin,i_k_point) = &
         KS_eigenvalue(i_state,i_spin,i_k_point) - &
         mulliken_decomp_total(i_atom, i_state, i_spin, i_k_point) * & 
         constraint_potential(i_region, i_spin)

    enddo
    enddo
    enddo
    enddo



  end subroutine get_energy_backshift
  !******



  !******
  !---------------------------------------------------------------------
  !****s* constraint/cleanup_constraint
  !  NAME
  !    cleanup_constraint
  !  SYNOPSIS

  subroutine cleanup_constraint()

    !  PURPOSE
    !  Subroutine cleanup_constraint deallocates all storage to do with constraints    


    !  USES
    implicit none
    !  ARGUMENTS
    !  INPUTS
    !    none
    !  OUTPUT
    !    none 
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE




    if (allocated(constraint_potential)) then
       deallocate( constraint_potential )
    end if

    if ( allocated (constraint_region) ) then
       deallocate ( constraint_region )
    end if
    if ( allocated (constraint_region_label) ) then
       deallocate ( constraint_region_label )
    end if
    if ( allocated (constraint_region_number) ) then
       deallocate ( constraint_region_number )
    end if
    if ( allocated (constraint_electrons) ) then
       deallocate ( constraint_electrons )
    end if
    if (allocated(constraint_pot_old)) then
       deallocate( constraint_pot_old )
    end if
    if (allocated(delta_constraint_potential)) then
       deallocate( delta_constraint_potential )
    end if
    if (allocated(delta_potential_mixing)) then
       deallocate( delta_potential_mixing )
    end if
    if (allocated(constraint_proj_final)) then
       deallocate( constraint_proj_final )
    end if
    if (allocated(hamiltonian_work)) then
       deallocate( hamiltonian_work )
    end if
    if (allocated(spin_shift)) then
       deallocate( spin_shift )
    end if
    if (allocated(full_ovlp)) then
       deallocate(full_ovlp)
    end if

    if (allocated(ovlp_times_eigenvec)) then
       deallocate(ovlp_times_eigenvec)
    end if

  end subroutine cleanup_constraint
  !******
  !---------------------------------------------------------------------
end module constraint
