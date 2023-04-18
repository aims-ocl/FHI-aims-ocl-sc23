!****h*  FHI-aims/transport
!  NAME
!    transport
!  SYNOPSIS
!  module transport

module transport

 !  PURPOSE
 !  Calculate electron transmission trought a nanostructure with semi-infinite leads
 !  using Landauer-Buttiker formalism (= Green's functions).
 !
 !  VB,WK,DN: Moved over Paula's transport module into git mainline. No
 !           comments as yet - this functionality is experimental and will
 !           evolve. Note that it is production ready for some things, Paula
 !           has used it for "real" physics after all.
 ! USES
 ! o use pbc_lists
 ! o use species_data
 ! o use synchronize_mpi

  use dimensions
  use runtime_choices
  use pbc_lists
  use species_data
  use synchronize_mpi

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
  !
  ! SOURCE


  implicit none



  complex*16, allocatable, dimension(:,:), private:: lead_self_energy, lead_self_energy_1,lead_self_energy_2, &
       lead_self_energy_3,  lead_self_energy_4
  integer, private:: lead_atoms(0:4), lead_basis(4), lead_atoms_start(4), lead_basis_start(4), lead_basis_end(4)
  integer,private:: n_leads, n_energy_steps
  character*30,private:: lead_file(4), tunne_file
  real*8,private:: energy_start, energy_end
  logical,private:: initialized =.false.
  integer,private::n_boundary_iterations
  real*8,private:: epsilon_end, boundary_treshold, boundary_mix, element_treshold
  complex*16,private:: epsilon_start

  real*8,private:: state_energy_max_r(0:4),state_energy_max_lead(0:4), average_pot(0:4)
  real*8,private:: state_energy_min_r(0:4), state_energy_min_lead(0:4), average_pot_lead(0:4)
  real*8,private:: delta_pot_level = 0.0d0
  integer, private,parameter:: n_spin_lead = 1
  logical,private :: out_min_eigenvalue = .false.
  integer, private:: gate_level_lead = 0
  logical, private:: fermi_level_fix = .false.

  integer,  dimension(3), private::  n_k_points_xyz_tr = 1
  integer, private :: n_k_points_tr = 1
  logical, private :: transport_k_points = .false.

  ! This is the data structure where the information of the leads is stored. Basically, it is one matrix
  ! entry in the periodic structure of aims. It is not the most compact format of storage but retains
  ! required information to recreate the matrices and the k-phase information. The fields are:
  !
  ! index -- absolute index of the entry in the compressed storage format matrix
  ! cell_1, cell_2, cell_3 -- supercell indeces associated with the entry
  ! i_basis_1, i_basis_2 -- row and column index of the entry in the 2d-matrix format
  ! ham -- values of the hamiltonian element, possibly for two spin channels
  ! ovlp -- values for the overlap matrix element
  type matrix_entry
     integer :: index
     integer :: cell_1
     integer :: cell_2
     integer :: cell_3
     integer :: i_basis_1
     integer :: i_basis_2
     real*8, dimension(2) :: ham
     real*8 :: ovlp
  end type matrix_entry

  complex*16, allocatable, dimension(:,:), private :: k_phase_base_tr
  real*8, allocatable, dimension(:), private :: k_weights_tr

  !****
contains
  !----------------------------------------------------------------
  !****s* transport/transport_read_control_in
  !  NAME
  !    transport_read_control_in
  !  SYNOPSIS
  !  subroutine transport_read_control_in(inputline)

  subroutine transport_read_control_in(inputline)

    !  PURPOSE
    !    Read the input settings from control.in for the transport calculations 
    !  INPUTS
    !    o inputline -- character line from the control.in, starting with "transport"
    !  OUTPUT
    !    none
    !  SOURCE

    use localorb_io, only: use_unit
    implicit none
    character(*),intent(in) :: inputline
    character*40 desc_str, desc_str2
    integer:: i_code, i_lead, i_atom, i_basis
    real*8:: r_temp


    ! Default values:
    if(.not. initialized)then
       write(tunne_file,'(A)') 'tunneling_probability'
       initialized = .true.
       n_leads = 0
       lead_basis = 0
       lead_atoms_start = 0
       lead_atoms = 0
       lead_basis = 0
       n_energy_steps = 0
       energy_start = 0.d0
       energy_end = 0.d0
       n_boundary_iterations = 10
       epsilon_end = 0.0001
       epsilon_start = (0,1)*0.02
       boundary_treshold = 1.0
       boundary_mix = 0.7d0
       element_treshold = 1.0E-10
    end if

    read (inputline,*,iostat=i_code) desc_str, desc_str2

    select case(desc_str2)

    case("lead_calculation")

       ! Calculate the boundary condition for the semi-infinite lead
       transport_lead_calculation = .true.
       if(transport_calculation)then
          if(myid==0) write(use_unit,*) 'Error: Both lead_calculation and transport_calculation on.'
          stop
       end if


    case("transport_calculation")
       ! Calculate the Landauer transport curve.
       transport_calculation  = .true.
       if(transport_lead_calculation)then
          if(myid==0) write(use_unit,*) 'Error: Both lead_calculation and transport_calculation on.'
          stop
       end if


    case('fermi_level_fix')
       ! Set the Fermi level of the semi-infinite leads and the center region to be
       ! at the same potential level
       fermi_level_fix = .true.


    case('output_min_eigenvalues')
       out_min_eigenvalue = .true.
       if(myid==0) write(use_unit,*) '| Smallest eigenvalues per atoms written out'


    case("epsilon_end")
       ! Small complex epsilon variable
       read (inputline,*) desc_str, desc_str2, epsilon_end
       if(myid==0) write(use_unit,*) '|End value for transport epsilon: ', epsilon_end 


    case("epsilon_start")
       ! Small complex epsilon variable
       read (inputline,*) desc_str, desc_str2, r_temp
       if(myid==0) write(use_unit,*) '|Begin value for transport epsilon: ', r_temp
       epsilon_start = r_temp*(0,1)


    case("delta_pot_level")
       ! Extra potential shift between the semi-infinite leads and the center region
       read (inputline,*) desc_str, desc_str2, delta_pot_level
       delta_pot_level = delta_pot_level/Hartree



    case("tunneling_file_name")
       ! The tunneling curve file name
       read (inputline,*) desc_str, desc_str2, tunne_file
       if(myid==0) write(use_unit,'(2X,A,A)') '|Tunneling propability is printed out to file: ', trim(tunne_file)


    case("energy_range")
       ! Energy range and number of steps for the tunneling curve.
       read (inputline,*) desc_str, desc_str2,   energy_start, energy_end,n_energy_steps
       if(myid==0) write(use_unit,'(2X,A,F8.3,A,F8.3,A)') '|Tunneling propability in range:', &
            energy_start,'eV - ',energy_end, 'eV'
       if(myid==0) write(use_unit,'(2X,A,I5)') '|Tunneling propability steps:',n_energy_steps
       energy_start = energy_start/Hartree
       energy_end = energy_end/Hartree


    case("number_of_boundary_iterations")
       ! The number of boundary iterations for the boundary conditions.
       read (inputline,*) desc_str, desc_str2, n_boundary_iterations
       if(myid==0) write(use_unit,*) '| Number of boundary iterations set to :',n_boundary_iterations

    case("boundary_treshold")
       ! The number of boundary iterations for the boundary conditions.
       read (inputline,*) desc_str, desc_str2, boundary_treshold
       if(myid==0) write(use_unit,*) '| Boundary treshold set to             :', boundary_treshold

    case("boundary_mix")
       ! The number of boundary iterations for the boundary conditions.
       read (inputline,*) desc_str, desc_str2, boundary_mix
       if(myid==0) write(use_unit,*) '| Boundary mix parameter set to       :', boundary_mix


    case("lead_1")
       ! Boundary conditions set up for the lead 1
       if(lead_atoms_start(1) >0)then
          write(use_unit,*) 'Error lead 1 defined 2 times'
          stop
       end if

       read (inputline,*) desc_str, desc_str2,  lead_atoms_start(1), lead_file(1)
       if(myid==0) write(use_unit,*) '| Lead 1: start atom:', lead_atoms_start(1),',file: ', trim(lead_file(1))
       n_leads = n_leads + 1
       

    case("lead_2")
       ! Boundary conditions set up for the lead 2
       if(lead_atoms_start(2) >0)then
          write(use_unit,*) 'Error lead 2 defined 2 times'
          stop
       end if

       read (inputline,*) desc_str,desc_str2,  lead_atoms_start(2), lead_file(2)
       if(myid==0) write(use_unit,*) '| Lead 2: start atom:', lead_atoms_start(2),',file: ', trim(lead_file(2))
       n_leads = n_leads + 1


    case("lead_3")
       ! Boundary conditions set up for the lead 3
       if(lead_atoms_start(3) >0)then
          write(use_unit,*) 'Error lead 3 defined 2 times'
          stop
       end if

       read (inputline,*) desc_str,desc_str2,  lead_atoms_start(3), lead_file(3)
       if(myid==0) write(use_unit,*) '| Lead 3: start atom:', lead_atoms_start(3),',file: ', trim(lead_file(3))
       n_leads = n_leads + 1


    case("lead_4")
       ! Boundary conditions set up for the lead 4
       if(lead_atoms_start(4) >0)then
          write(use_unit,*) 'Error lead 4 defined 2 times'
          stop
       end if

       read (inputline,*) desc_str,desc_str2,  lead_atoms_start(4), lead_file(4)
       if(myid==0) write(use_unit,*) '| Lead 4: start atom:', lead_atoms_start(4),',file: ', trim(lead_file(4))
       n_leads = n_leads + 1

    case("element_treshold")
       ! The number of boundary iterations for the boundary conditions.
       read (inputline,*) desc_str, desc_str2, element_treshold
       if(myid==0) write(use_unit,*) '| Element treshold for lead matrices set to:', element_treshold

    case("k_grid")
       read(inputline,*) desc_str, desc_str2, n_k_points_xyz_tr(1), n_k_points_xyz_tr(2), n_k_points_xyz_tr(3)
       transport_k_points = .true.
       n_k_points_tr = n_k_points_xyz_tr(1)*n_k_points_xyz_tr(2)*n_k_points_xyz_tr(3)

    case default
       if(myid==0) write(use_unit,*) '*** Error: The following transport subkeyword does not exist:', trim(desc_str2)
       call aims_stop

    end select





  end subroutine transport_read_control_in
  !******----------------------------------------------------------
  !****s* transport/transport_check_control_in
  !  NAME
  !    transport_check_control_in
  !  SYNOPSIS
  !    subroutine transport_check_control_in

  subroutine transport_check_control_in
 
    !  PURPOSE
    !    Make some checks for the user defined transport input parameters.
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  SOURCE


   use localorb_io, only: use_unit
   implicit none


    if(transport_calculation)then

       ! Check are all the leads defined

       if(lead_atoms_start(4) /=0 )then
          if(lead_atoms_start(3) ==0)then
             write(use_unit,*) 'Error: lead 4 defined, but not lead 3!'
             call Aims_stop
          end if
       end if
       if(lead_atoms_start(3) /=0 )then
          if(lead_atoms_start(2) ==0)then
             write(use_unit,*) 'Error: lead 3 defined, but not lead 2!'
             call Aims_stop
          end if
       end if
       if(lead_atoms_start(2) ==0 )then
          write(use_unit,*) 'Error: lead 2 not defined!'
          call Aims_stop
       end if
       if(lead_atoms_start(1) ==0 )then
          write(use_unit,*) 'Error: lead 1 not defined!'
          call Aims_stop
       end if
   
       if(maxval(lead_atoms_start)> n_atoms)then
          write(use_unit,*) 'Error: leads starting point larger than the total number of atoms!'
          call aims_stop
       end if

       if (transport_k_points) then
          if (n_k_points_tr < 1) then
             if (myid==0) write(use_unit,*) 'Error: Less than one k-point for transport calculation!'
             call aims_stop
          end if

          if (n_k_points_xyz_tr(3) > 1) then
             if (myid==0) write(use_unit,*) 'Error: Only one k-point for the third lattice vector allowed!'
             call aims_stop
          end if

          if (lead_atoms_start(3) > 0 .or. lead_atoms_start(4) > 0) then
             if (myid==0) write(use_unit,*) 'Error: Only two leads allowed when using k-points!'
             call aims_stop
          end if

       end if

    end if


  end subroutine transport_check_control_in
  !******----------------------------------------------------------
  !****s* transport/transport_write_boundary_conditions
  !  NAME
  !    transport_write_boundary_conditions
  !  SYNOPSIS
  !      subroutine transport_write_boundary_conditions( overlap_matrix, hamiltonian, chemical_potential)

  subroutine transport_write_boundary_conditions( overlap_matrix, hamiltonian, chemical_potential)
 
    !  PURPOSE
    !    Transport semi-infinite lead calculation: write out the lead_self_energy file.
    !  INPUTS
    !    o overlap_matrix
    !    o hamiltonian
    !    o chemical_potential
    !  OUTPUT
    !    none
    !  SOURCE


    ! Write boundary conditions to the transport calculations to the file "lead_self_energy"

    use localorb_io, only: use_unit
    implicit none
    real*8:: overlap_matrix( n_hamiltonian_matrix_size )
    real*8:: hamiltonian ( n_hamiltonian_matrix_size,n_spin )
    real*8:: chemical_potential

    integer:: max_overlap_distance, i_center, i_center_L, i_layer, i_atom, n_lead_basis, i_basis_1, i_basis_2
    integer:: i_matrix_element_non_zero
    real*8,    allocatable, dimension(:,:)::   lead_ovlp_matrix, lead_connection_matrix_o
    real*8,    allocatable, dimension(:,:,:):: lead_connection_matrix_h, lead_hamiltonian
    logical,   allocatable, dimension(:,:):: matrix_element_non_zero

    integer:: i_cell, i_size, i_spin, i_index_real

    if (transport_k_points) then
       call transport_write_boundary_conditions_k_points( overlap_matrix, hamiltonian,chemical_potential )
       return
    end if

    if(myid==0) then


       write(use_unit,*) '----------------------------------------------------------------'
       write(use_unit,*) 'Writing the boundary conditions for the semi-infinite leads'
       write(use_unit,*) 'to the file: lead_self_energy.'



       ! How many times the lead have to be multiplied?
       ! The atoms at "left" should not tought to "right" atoms.
       ! This is why the lead have to be multiplied to be large enough.
       ! However, this property currently do not work, the geometry of the
       ! system have to be just large enough.

       max_overlap_distance = 0

       do i_center_L = 1, n_centers_integrals

          i_center= centers_basis_integrals(i_center_L)
          max_overlap_distance = max(max_overlap_distance, abs(cell_index(center_to_cell(i_center),3)))

       end do
       write(use_unit,*) '| The number of layers needed to the lead:',max_overlap_distance

       n_lead_basis = max_overlap_distance *n_basis
       write(use_unit,*) '| The number of lead basis:',n_lead_basis


       ! These variables are for the plan to make the region multiplication to work:
       allocate(lead_connection_matrix_o(n_lead_basis, n_lead_basis))
       allocate(lead_connection_matrix_h(n_lead_basis, n_lead_basis,n_spin))
       allocate(lead_ovlp_matrix(n_lead_basis, n_lead_basis))
       allocate(lead_hamiltonian(n_lead_basis, n_lead_basis,n_spin))
       allocate(matrix_element_non_zero(n_lead_basis, n_lead_basis))


       lead_connection_matrix_h = 0.d0
       lead_connection_matrix_o = 0.d0
       lead_ovlp_matrix = 0.d0
       lead_hamiltonian = 0.d0
       matrix_element_non_zero = .false.

       ! Evaluate first non-zero elements.

       do i_cell = 1,n_cells_in_hamiltonian-1
          do i_basis_2 = 1, n_basis
             if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1
                do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
                   i_index_real = i_index_real + 1
                   i_basis_1 =  column_index_hamiltonian(i_index_real)

                   if(cell_index(i_cell,3) == 0)then

                      do i_spin = 1, n_spin
                         if( abs(hamiltonian(  i_index_real, i_spin)) > element_treshold) then
                            matrix_element_non_zero (i_basis_1, i_basis_2) = .true.
                            matrix_element_non_zero (i_basis_2, i_basis_1) = .true.
                         end if
                      end do
                      if( abs(overlap_matrix(  i_index_real)) > element_treshold) then
                         matrix_element_non_zero (i_basis_1, i_basis_2) = .true.
                         matrix_element_non_zero (i_basis_2, i_basis_1) = .true.
                      end if

                   else if(cell_index(i_cell,3) ==  1)then

                      do i_spin = 1, n_spin
                         if( abs(hamiltonian(  i_index_real, i_spin)) > element_treshold)then
                            matrix_element_non_zero (i_basis_2, i_basis_1) = .true.
                         end if
                      end do
                      if( abs(overlap_matrix(  i_index_real)) > element_treshold) then
                         matrix_element_non_zero (i_basis_2, i_basis_1) = .true.
                      end if

                   else if(cell_index(i_cell,3) ==  -1)then

                      do i_spin = 1, n_spin
                         if( abs(hamiltonian(  i_index_real, i_spin)) > element_treshold)then
                            matrix_element_non_zero (i_basis_1, i_basis_2) = .true.
                         end if
                      end do
                      if( abs(overlap_matrix(  i_index_real)) > element_treshold) then
                         matrix_element_non_zero (i_basis_1, i_basis_2) = .true.
                      end if

                   end if
                end do
             end if
          end do
       end do


       ! How many non-zero elements ?
       i_matrix_element_non_zero = 0

       do i_basis_1 = 1, n_lead_basis
          do i_basis_2 = 1, n_lead_basis
             if(matrix_element_non_zero(i_basis_1, i_basis_2))then
                i_matrix_element_non_zero = i_matrix_element_non_zero + 1
             end if
          end do
       end do



       ! Collect matrix elements of Hamiltonian and overlap-matrix
       do i_cell = 1,n_cells_in_hamiltonian-1

          do i_basis_2 = 1, n_basis

             if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

                i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1

                do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                   i_index_real = i_index_real + 1
                   i_basis_1 =  column_index_hamiltonian(i_index_real)

                   if(cell_index(i_cell,3) == 0)then

                      do i_spin = 1, n_spin

                         lead_hamiltonian( i_basis_1,i_basis_2 ,i_spin) =  &
                              hamiltonian(  i_index_real, i_spin)
                      end do

                      lead_ovlp_matrix(i_basis_1, i_basis_2) = overlap_matrix(  i_index_real)

                      if(i_basis_1 /= i_basis_2)then

                         do i_spin = 1, n_spin

                            lead_hamiltonian( i_basis_2,i_basis_1 ,i_spin) =  &
                                 hamiltonian(  i_index_real, i_spin)
                         end do

                         lead_ovlp_matrix(i_basis_2, i_basis_1) = overlap_matrix(  i_index_real)
                      end if

                   else if(cell_index(i_cell,3) ==  1)then

                      do i_spin = 1, n_spin

                         lead_connection_matrix_h(i_basis_2, i_basis_1, i_spin) =   hamiltonian(  i_index_real, i_spin)

                      end do

                      lead_connection_matrix_o(i_basis_2, i_basis_1) =  overlap_matrix(  i_index_real)


                   else if(cell_index(i_cell,3) ==  -1)then

                      do i_spin = 1, n_spin

                         lead_connection_matrix_h(i_basis_1, i_basis_2, i_spin) =   hamiltonian(  i_index_real, i_spin)

                      end do

                      lead_connection_matrix_o(i_basis_1, i_basis_2) =  overlap_matrix(  i_index_real)

                   end if
                end do
             end if
          end do
       end do




       ! Write out the lead information: atom coordinates and self energy matrix.
       open(88, file='lead_self_energy')

       write(88,*) n_atoms *  max_overlap_distance,   n_lead_basis
       write(88,*) average_pot(0)
       write(88,*) chemical_potential
       write(88,*) i_matrix_element_non_zero, n_spin


       do i_basis_1 = 1, n_lead_basis
          do i_basis_2 = 1, n_lead_basis

             if(matrix_element_non_zero(i_basis_1, i_basis_2))then

                write(88,'(2I5)') i_basis_1, i_basis_2

                write(88,'(4E16.8)') lead_hamiltonian(i_basis_1, i_basis_2,1), lead_ovlp_matrix(i_basis_1, i_basis_2), &
                     lead_connection_matrix_h(i_basis_1, i_basis_2,1), lead_connection_matrix_o(i_basis_1, i_basis_2)

             end if
          end do
       end do

       ! If spin polarized calculation, we need to write more.
       if(n_spin > 1)then

          do i_basis_1 = 1, n_lead_basis
             do i_basis_2 = 1, n_lead_basis

                if(matrix_element_non_zero(i_basis_1, i_basis_2))then

                   write(88,'(2I5)') i_basis_1, i_basis_2
                   write(88,'(4E16.8)') lead_hamiltonian(i_basis_1, i_basis_2,2), &
                        lead_connection_matrix_h(i_basis_1, i_basis_2,2)

                end if
             end do
          end do
       end if


       close(88)

       write(use_unit,*) '----------------------------------------------------------------'

    end if ! myid==0


    deallocate(lead_connection_matrix_o)
    deallocate(lead_connection_matrix_h)
    deallocate(lead_ovlp_matrix)
    deallocate(lead_hamiltonian)
    deallocate(matrix_element_non_zero)


  end subroutine transport_write_boundary_conditions
  !******----------------------------------------------------------
  !****s* transport/transport_write_boundary_conditions_k_points
  !  NAME
  !    transport_write_boundary_conditions_k_points
  !  SYNOPSIS
  !      subroutine transport_write_boundary_conditions_k_points( overlap_matrix, hamiltonian, chemical_potential)

  subroutine transport_write_boundary_conditions_k_points( overlap_matrix, hamiltonian, chemical_potential)
 
    !  PURPOSE
    !    Transport semi-infinite lead calculation: write out the lead_self_energy file when k-points are used
    !    i.e. when the transport is across an interface that is periodic in two dimensions.
    !  INPUTS
    !    o overlap_matrix
    !    o hamiltonian
    !    o chemical_potential
    !  OUTPUT
    !    none
    !  SOURCE


    ! Write boundary conditions to the transport calculations to the file "lead_self_energy"

    use localorb_io, only: use_unit
    implicit none
    real*8:: overlap_matrix( n_hamiltonian_matrix_size )
    real*8:: hamiltonian ( n_hamiltonian_matrix_size,n_spin )
    real*8:: chemical_potential

    integer:: max_overlap_distance, i_center, i_center_L, i_layer, i_atom, n_lead_basis, i_basis_1, i_basis_2
    integer:: i_matrix_element_non_zero
    logical, allocatable, dimension(:) :: matrix_element_non_zero

    integer:: i_cell, i_size, i_spin, i_index_real, i_index

    type(matrix_entry), pointer, dimension(:) :: matrix_element

    if(myid==0) then

       write(use_unit,*) '----------------------------------------------------------------'
       write(use_unit,*) 'Writing the boundary conditions for the semi-infinite leads'
       write(use_unit,*) 'for interfacial transport to the file: lead_self_energy_k_points.'



       ! How many times the lead have to be multiplied?
       ! The atoms at "left" should not tought to "right" atoms.
       ! This is why the lead have to be multiplied to be large enough.
       ! However, this property currently do not work, the geometry of the
       ! system have to be just large enough.

       max_overlap_distance = 0

       do i_center_L = 1, n_centers_integrals

          i_center= centers_basis_integrals(i_center_L)
          max_overlap_distance = max(max_overlap_distance, abs(cell_index(center_to_cell(i_center),3)))

       end do
       write(use_unit,*) '| The number of layers needed to the lead:',max_overlap_distance

       n_lead_basis = max_overlap_distance *n_basis
       write(use_unit,*) '| The number of lead basis:',n_lead_basis


       ! These variables are for the plan to make the region multiplication to work:
       allocate(matrix_element_non_zero(n_hamiltonian_matrix_size))

       matrix_element_non_zero = .false.

       ! Evaluate first non-zero elements.

       do i_cell = 1,n_cells_in_hamiltonian-1
          do i_basis_2 = 1, n_basis
             if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1
                do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
                   i_index_real = i_index_real + 1
                   i_basis_1 =  column_index_hamiltonian(i_index_real)

                   if ((cell_index(i_cell,3) == -1).or.(cell_index(i_cell,3) == 0).or. &
                        (cell_index(i_cell,3) == 1)) then

                      do i_spin = 1, n_spin
                         if( abs(hamiltonian(  i_index_real, i_spin)) > element_treshold) then
                            matrix_element_non_zero (i_index_real) = .true.
                         end if
                      end do

                      if( abs(overlap_matrix(  i_index_real)) > element_treshold) then
                         matrix_element_non_zero (i_index_real) = .true.
                      end if

                   end if
                end do
             end if
          end do
       end do

       ! How many non-zero elements ?
       i_matrix_element_non_zero = 0
       do i_basis_1 = 1, n_hamiltonian_matrix_size
          if(matrix_element_non_zero(i_basis_1)) i_matrix_element_non_zero = i_matrix_element_non_zero + 1
       end do

       allocate(matrix_element(i_matrix_element_non_zero))

       ! Collect matrix elements of Hamiltonian and overlap-matrix
       i_index = 0
       do i_cell = 1,n_cells_in_hamiltonian-1

          do i_basis_2 = 1, n_basis

             if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

                i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1

                do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                   i_index_real = i_index_real + 1
                   i_basis_1 =  column_index_hamiltonian(i_index_real)

                   if (matrix_element_non_zero(i_index_real)) then
                      i_index = i_index + 1
                      matrix_element(i_index)%index = i_index_real
                      matrix_element(i_index)%cell_1 = cell_index(i_cell,1)
                      matrix_element(i_index)%cell_2 = cell_index(i_cell,2)
                      matrix_element(i_index)%cell_3 = cell_index(i_cell,3)
                      matrix_element(i_index)%i_basis_1 = i_basis_1
                      matrix_element(i_index)%i_basis_2 = i_basis_2
                      do i_spin = 1, n_spin
                         matrix_element(i_index)%ham(i_spin) = hamiltonian(i_index_real, i_spin)
                      end do
                      matrix_element(i_index)%ovlp = overlap_matrix(i_index_real)
                   end if
                   
                end do
             end if
          end do
       end do


       ! Write out the lead information: atom coordinates and self energy matrix.
       open(88, file='lead_self_energy_k_points')

       write(88,*) n_atoms *  max_overlap_distance,   n_lead_basis
       write(88,*) average_pot(0)
       write(88,*) chemical_potential
       write(88,*) i_matrix_element_non_zero, n_spin

       do i_index = 1, i_matrix_element_non_zero
          write(88,'(1I10)') matrix_element(i_index)%index
          write(88,'(5I7)') matrix_element(i_index)%i_basis_1, matrix_element(i_index)%i_basis_2, &
               matrix_element(i_index)%cell_1, matrix_element(i_index)%cell_2, matrix_element(i_index)%cell_3
          write(88,'(2E16.8)') matrix_element(i_index)%ham(1), matrix_element(i_index)%ovlp
          if(n_spin > 1)then
             write(88,'(1E16.8)') matrix_element(i_index)%ham(2)
          end if
       end do

       close(88)

       write(use_unit,*) '----------------------------------------------------------------'

    end if ! myid==0

    deallocate(matrix_element_non_zero)
    deallocate(matrix_element)

  end subroutine transport_write_boundary_conditions_k_points
  !******----------------------------------------------------------
  !****s* transport/read_lead_information
  !  NAME
  !    read_lead_information 
  !  SYNOPSIS
  !    subroutine read_lead_information( overlap_matrix, hamiltonian, chemical_potential)

  subroutine read_lead_information( overlap_matrix, hamiltonian, chemical_potential)
 
    !  PURPOSE
    !    Evaluation electron tunneling propability. 
    !    Note that for the scalapack calculations there is different subroutine.
    !  
    !  INPUTS
    !    o overlap_matrix
    !    o hamiltonian
    !    o chemical_potential
    !  OUTPUT
    !    none
    !  SOURCE

    use basis, only: basis_atom, n_basis_atom
    use localorb_io, only: use_unit
    implicit none

    real*8:: overlap_matrix( n_hamiltonian_matrix_size )
    real*8:: hamiltonian ( n_hamiltonian_matrix_size,n_spin )

    integer:: i_lead, i_basis_1, i_basis_2, n_lead_basis

    real*8, allocatable, dimension(:,:,:,:) :: lead_hamiltonian
    real*8, allocatable, dimension(:,:,:) :: lead_ovlp_matrix
    
    real*8, allocatable, dimension(:,:,:,:) :: lead_connection_matrix_h
    real*8, allocatable, dimension(:,:,:) :: lead_connection_matrix_o

    complex*16,    allocatable, dimension(:,:,:,:) :: lead_hamiltonian_cmplx
    complex*16,    allocatable, dimension(:,:,:) :: lead_ovlp_matrix_cmplx

    complex*16,    allocatable, dimension(:,:,:,:) :: lead_connection_matrix_h_cmplx
    complex*16,    allocatable, dimension(:,:,:) :: lead_connection_matrix_o_cmplx

    character*31 :: file_name

    integer:: ipiv(n_basis), n_boundary_iterations_max
    integer:: info, i_1, i_2, i_spin, i_cell, i_size, i_index_real, i_energy, i_index
    integer:: i_L1_1, i_L1_2, i_L2_1, i_L2_2, n_lead_basis_max, i_atom, i_basis

    complex*16, allocatable, dimension(:,:):: green, greenL,work, workL, ttt,ttt2, gamma1, gamma2
    complex*16:: tun(4,n_tasks), tun_prev
    real*8:: chemical_potential,  energy_print(n_tasks)
    complex*16:: energy, energy_lead, delta_energy, epsilon, epsilon_start_max
    real*8:: prev_max, lead_chem_pot(4), delta_fermi_level
    integer:: i_task, i_k_point
    complex*16 :: k_phase_tr
    
    integer :: lead_matrix_elements(4)
    
    type(matrix_entry), pointer, dimension(:) :: lead_1, lead_2

    if(myid==0) write(use_unit,*) 'Evaluating tunneling probability'

    n_boundary_iterations_max = max(30, n_boundary_iterations)
    epsilon_start_max = 0.1

    ! Where the boundary conditions are put.
    do i_lead = 1, n_leads

       do i_basis_1 = 1,n_basis
          if(Cbasis_to_atom(i_basis_1)== lead_atoms_start(i_lead))then
             lead_basis_start(i_lead) = i_basis_1
             exit
          end if
       end do
    end do



    ! The 2-4 leads and boundary files.

    do i_lead = 1, n_leads

       ! Read the boundary information from the file.
       ! Every lead have own boundary information file.
       ! First, the headers

       open(88, file=lead_file(i_lead))
       read(88,*) lead_atoms( i_lead ), lead_basis( i_lead)
       read(88,*) average_pot_lead(i_lead)
       read(88,*) lead_chem_pot(i_lead)
       read(88,*) lead_matrix_elements(i_lead)

       if(myid==0) write(use_unit,'(2X,A)') ' '
       if(myid==0) write(use_unit,'(2X,A,I1,A)') 'Lead ',i_lead,   ' :'
       if(myid==0) write(use_unit,'(2X,A,I5,A,I5)') '| Including atoms:   ', lead_atoms_start( i_lead), & 
                   '  to',  lead_atoms_start( i_lead)+ lead_atoms( i_lead )-1
       if(myid==0) write(use_unit,'(2X,A,I5)') '| Starting basis:  ', lead_basis_start(i_lead)
       if(myid==0) write(use_unit,'(2X,A,I5)') '| Number of atoms: ', lead_atoms( i_lead )
       if(myid==0) write(use_unit,'(2X,A,I5)') '| Number of basis: ', lead_basis( i_lead)
       if(myid==0) write(use_unit,'(2X,A,F16.8)') '| Pot ave  : ',  average_pot_lead(i_lead)*Hartree
       if(myid==0) write(use_unit,'(2X,A,F16.8)') '| Lead chem. pot.  : ',  lead_chem_pot(i_lead)*Hartree

       lead_basis_end(i_lead) = lead_basis_start(i_lead) + lead_basis( i_lead) -1

       close(88)

    end do

    n_lead_basis = maxval(lead_basis)

    n_basis_atom = 0
    do i_basis = 1, n_basis
       n_basis_atom(basis_atom(i_basis)) =  n_basis_atom(basis_atom(i_basis)) + 1
    end do
    
    if (transport_k_points) then
       allocate(lead_hamiltonian_cmplx(n_lead_basis, n_lead_basis,n_spin_lead,n_leads))
       lead_hamiltonian_cmplx = 0.d0
       allocate(lead_ovlp_matrix_cmplx(n_lead_basis, n_lead_basis,n_leads))
       lead_ovlp_matrix_cmplx = 0.d0
       allocate(lead_connection_matrix_h_cmplx(n_lead_basis, n_lead_basis, n_spin_lead, n_leads))
       lead_connection_matrix_h_cmplx = 0.0d0
       allocate(lead_connection_matrix_o_cmplx(n_lead_basis, n_lead_basis, n_leads))
       lead_connection_matrix_o_cmplx = 0.0d0
       ! Dummy allocations:
       allocate(lead_hamiltonian(1,1,1,1))
       lead_hamiltonian = 0.d0
       allocate(lead_ovlp_matrix(1,1,1))
       lead_ovlp_matrix = 0.d0
       allocate(lead_connection_matrix_h(1,1,1,1))
       lead_connection_matrix_h = 0.0d0
       allocate(lead_connection_matrix_o(1,1,1))
       lead_connection_matrix_o = 0.0d0
       !----------------------------------------------------------------
    else
       allocate(lead_hamiltonian(n_lead_basis, n_lead_basis, n_spin_lead, n_leads))
       lead_hamiltonian = 0.d0
       allocate(lead_ovlp_matrix(n_lead_basis, n_lead_basis, n_leads))
       lead_ovlp_matrix = 0.d0
       allocate(lead_connection_matrix_h(n_lead_basis, n_lead_basis, n_spin_lead, n_leads))
       lead_connection_matrix_h = 0.0d0
       allocate(lead_connection_matrix_o(n_lead_basis, n_lead_basis, n_leads))
       lead_connection_matrix_o = 0.0d0

       ! Dummy allocations:
       allocate(lead_hamiltonian_cmplx(1,1,1,1))
       lead_hamiltonian_cmplx = 0.d0
       allocate(lead_ovlp_matrix_cmplx(1,1,1))
       lead_ovlp_matrix_cmplx = 0.d0
       allocate(lead_connection_matrix_h_cmplx(1,1,1,1))
       lead_connection_matrix_h_cmplx = 0.0d0
       allocate(lead_connection_matrix_o_cmplx(1,1,1))
       lead_connection_matrix_o_cmplx = 0.0d0
       !----------------------------------------------------------------
    end if

    do i_lead = 1, n_leads

       ! Re-read the header, could also be just dummies
       open(88, file=lead_file(i_lead))
       read(88,*) lead_atoms( i_lead ), lead_basis( i_lead)
       read(88,*) average_pot_lead(i_lead)
       read(88,*) lead_chem_pot(i_lead)
       read(88,*) lead_matrix_elements(i_lead)

       ! Next, read the entries of the matrices for the leads    
       if (transport_k_points) then

          select case(i_lead)

          case(1)

             allocate(lead_1(lead_matrix_elements(i_lead)))
             do i_index = 1, lead_matrix_elements(i_lead)
                read(88,*) lead_1(i_index)%index
                read(88,*) lead_1(i_index)%i_basis_1, lead_1(i_index)%i_basis_2, &
                     lead_1(i_index)%cell_1, lead_1(i_index)%cell_2, lead_1(i_index)%cell_3
                read(88,*) lead_1(i_index)%ham(1), lead_1(i_index)%ovlp
                if(n_spin > 1)then
                   read(88,*) lead_1(i_index)%ham(2)
                end if

             end do

             allocate(lead_self_energy_1(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_1 = 0.d0
             
          case(2)

             allocate(lead_2(lead_matrix_elements(i_lead)))
             do i_index = 1, lead_matrix_elements(i_lead)
                read(88,*) lead_2(i_index)%index
                read(88,*) lead_2(i_index)%i_basis_1, lead_2(i_index)%i_basis_2, &
                     lead_2(i_index)%cell_1, lead_2(i_index)%cell_2, lead_2(i_index)%cell_3
                read(88,*) lead_2(i_index)%ham(1), lead_2(i_index)%ovlp
                if(n_spin > 1)then
                   read(88,*) lead_2(i_index)%ham(2)
                end if

             end do

             allocate(lead_self_energy_2(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_2 = 0.d0

             
          end select

       else ! we have no k-points

          ! More information from the boundary condition files.

          select case(i_lead)

          case(1)

             allocate(lead_self_energy_1(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_1 = 0.d0
                          
          case(2)
             
             allocate(lead_self_energy_2(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_2 = 0.d0
          case(3)
             
             allocate(lead_self_energy_3(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_3 = 0.d0
             
          case(4)
             
             allocate(lead_self_energy_4(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_4 = 0.d0
                          
          end select

          do i_1 = 1, lead_matrix_elements(i_lead)
                
             read(88,*)  i_basis_1,  i_basis_2
             read(88,*) lead_hamiltonian(i_basis_1, i_basis_2,1,i_lead), lead_ovlp_matrix(i_basis_1, i_basis_2, i_lead), &
                  lead_connection_matrix_h(i_basis_1, i_basis_2, 1, i_lead), lead_connection_matrix_o(i_basis_1, i_basis_2, i_lead)
                
          end do

       end if ! k-points or not

       close(88)

       ! Some checks. There is easily errors in control.in

       if(lead_atoms_start(i_lead)+ lead_atoms(i_lead)-1 >n_atoms)then
          write(use_unit,*) 'Error: lead',i_lead, &
             ' start atom must be wrong the final atom is larger than ', &
             'number of atoms in the system!' 
          stop
       end if

       i_basis = 0
       do i_atom =   lead_atoms_start(i_lead), lead_atoms_start(i_lead)+ lead_atoms(i_lead)-1, 1 
          i_basis = i_basis +  n_basis_atom(i_atom)
       end do

!       write(use_unit,*) 'file:',lead_basis(i_lead), 'system:', i_basis 

       if(i_basis /= lead_basis(i_lead))then
          write(use_unit,*) 'Error: the boundary file for lead',i_lead, 'do not match to the system'
          write(use_unit,*) 'file:',lead_basis(i_lead), 'system:', i_basis 
          stop
       end if

    end do ! i_leads

    if (transport_k_points) call transport_init_k_points()

    allocate(work(n_basis, n_basis))
    work = 0.d0
    allocate(green(n_basis, max(lead_basis(1), lead_basis(3))))
    green = 0.d0

    allocate(workL(n_lead_basis, n_lead_basis))
    workL = 0.d0
    allocate(greenL(n_lead_basis, n_lead_basis))
    greenL = 0.d0
    allocate(ttt(n_lead_basis, n_lead_basis))
    ttt = 0.d0
    allocate(ttt2(n_lead_basis, n_lead_basis))
    ttt2 = 0.d0

    ! Results are saved to "tunne_file"

    if(n_spin ==1)then
       open(22,file=tunne_file)
    else

       write(file_name,'(A,I1)') trim(tunne_file), 1
       open(22,file=file_name)

       write(file_name,'(A,I1)') trim(tunne_file), 2
       open(23,file=file_name)
    end if

    do i_lead = 1, n_leads
       if(myid==0) &
          write(use_unit,'(A,I2,A,F16.8)') ' Fermi level - Reference potential for lead', &
             i_lead,':', (chemical_potential- average_pot(i_lead)) * hartree
       if(myid==0) &
          write(use_unit,'(A,I2,A,F16.8)') ' Reference gate voltage for lead: ',i_lead,':', &
             (-(-average_pot(i_lead) + chemical_potential) &
             + (- average_pot_lead(i_lead) + lead_chem_pot(i_lead)))*Hartree

    end do

    ! Possibility to gate.
    if(gate_level_lead /=0)then
       delta_fermi_level = (-( -average_pot(gate_level_lead) + chemical_potential)  &
            + (- average_pot_lead(gate_level_lead)+  lead_chem_pot(gate_level_lead))) 

       if(myid==0) write(use_unit,*) '| delta Fermi level: ', delta_fermi_level*Hartree
    else
       delta_fermi_level = 0.0d0
    end if

    do i_lead = 1, n_leads
        if(myid==0) write(use_unit,'(A,I2,A,F16.8)') ' Reference gate voltage after fix for lead: ',i_lead,':', ( &
            -( -average_pot(i_lead) + chemical_potential + delta_fermi_level )  &
          + (- average_pot_lead(i_lead)+  lead_chem_pot(i_lead)))*Hartree

        if(myid==0) write(22,'(A,I2,2F16.8)') '# Reference gate for lead: ', i_lead, &
             ( -( -average_pot(i_lead) + chemical_potential + delta_fermi_level )  &
             + (- average_pot_lead(i_lead)+  lead_chem_pot(i_lead)))*Hartree, delta_fermi_level
    end do

    if(fermi_level_fix)then
       do i_lead = 1, n_leads

          average_pot_lead(i_lead)= lead_chem_pot(i_lead)
          average_pot(i_lead) = chemical_potential

       end do
    end if

    if(0.02 < epsilon_end) epsilon_start =  (0,1)*epsilon_end

    ! Loop over k-points
    do i_k_point = 1, n_k_points_tr

       if ((n_k_points_tr > 1).and.(myid==0)) then
          write(use_unit,'(A,I4)') '# K-point: ', i_k_point
          write(22,'(A,I4)') '# K-point: ', i_k_point
       end if

       if (transport_k_points) then
          call construct_lead_matrices(lead_1, lead_2, lead_hamiltonian_cmplx, lead_ovlp_matrix_cmplx, &
               lead_connection_matrix_h_cmplx, lead_connection_matrix_o_cmplx, &
               lead_matrix_elements(1), lead_matrix_elements(2), n_lead_basis, i_k_point)
       end if

       ! The loop over energy
       do i_energy = myid,n_energy_steps-1, n_tasks
          
          energy = chemical_potential+delta_fermi_level + energy_start  + (energy_end-energy_start)/ (n_energy_steps-1) * (i_energy)

          call evaluate_lead_self_energies(lead_hamiltonian, lead_ovlp_matrix, lead_connection_matrix_h, lead_connection_matrix_o, &
               lead_hamiltonian_cmplx, lead_ovlp_matrix_cmplx, lead_connection_matrix_h_cmplx, lead_connection_matrix_o_cmplx, &
               n_lead_basis, ipiv, n_boundary_iterations_max, energy, epsilon_start_max)

          ! Construct the matrix for the Green's function of the center region

          do i_spin = 1, n_spin

             energy = real(energy) 
             work = 0.d0

             do i_cell = 1,n_cells_in_hamiltonian-1
             
                do i_basis_2 = 1, n_basis

                   if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                   
                      i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1

                      do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                         i_index_real = i_index_real + 1
                         i_basis_1 =  column_index_hamiltonian(i_index_real)

                         if (cell_index(i_cell,3) == 0) then
                            
                            if (transport_k_points) then

                               k_phase_tr = (k_phase_base_tr(1,i_k_point)**cell_index(i_cell,1))* &
                                    (k_phase_base_tr(2,i_k_point)**cell_index(i_cell,2))
                               
                               work( i_basis_1,i_basis_2) = work( i_basis_1,i_basis_2) + &
                                    k_phase_tr*( energy*overlap_matrix(i_index_real) - &
                                    hamiltonian(i_index_real,i_spin) )

                               if (i_basis_1 /= i_basis_2) then

                                  work( i_basis_2,i_basis_1 ) =  work( i_basis_2,i_basis_1) + &
                                       dconjg(k_phase_tr*( energy*overlap_matrix(i_index_real) - &
                                       hamiltonian(i_index_real,i_spin) ))
                               
                               end if

                            else if(cell_index(i_cell,1) == 0 .and. cell_index(i_cell,2) == 0) then
                               
                               work( i_basis_1,i_basis_2) = energy*overlap_matrix(  i_index_real) -   &
                                    hamiltonian(  i_index_real, i_spin)

                               if(i_basis_1 /= i_basis_2)then

                                  work( i_basis_2,i_basis_1 ) = energy*overlap_matrix(  i_index_real) -   &
                                       hamiltonian(  i_index_real, i_spin)                      

                               end if
                            end if
                         end if
                      end do
                   end if
                end do
             end do

             ! Add the boundary conditions to the matrix of center cell

             do i_lead = 1, n_leads

                i_1 = lead_basis_start(i_lead)
                i_2 = lead_basis_start(i_lead)+lead_basis(i_lead)-1

                select case(i_lead)
                case(1)
                   work( i_1:i_2, i_1:i_2 ) =  work( i_1:i_2, i_1:i_2)  - lead_self_energy_1(1:lead_basis(1), 1:lead_basis(1))

                case(2)
                   work( i_1:i_2, i_1:i_2 ) =  work( i_1:i_2, i_1:i_2)  - lead_self_energy_2(1:lead_basis(2), 1:lead_basis(2))

                case(3)
                   work( i_1:i_2, i_1:i_2 ) =  work( i_1:i_2, i_1:i_2)  - lead_self_energy_3(1:lead_basis(3), 1:lead_basis(3))

                case(4)
                   work( i_1:i_2, i_1:i_2 ) =  work( i_1:i_2, i_1:i_2)  - lead_self_energy_4(1:lead_basis(4), 1:lead_basis(4))
                   
                end select
             end do


             ! Now we can solve the Green's function of the system.


             i_L1_1 = lead_basis_start(1)
             i_L1_2 = lead_basis_start(1)+lead_basis(1)-1

             i_L2_1 = lead_basis_start(2)
             i_L2_2 = lead_basis_start(2)+lead_basis(2)-1

             
             call ZGETRF(n_basis ,n_basis, work, n_basis, ipiv, info )


             green = 0.d0
             do i_basis_1 =  1, lead_basis(1)
                green(lead_basis_start(1)+i_basis_1-1, i_basis_1) = 1.d0
             end do


             call  ZGETRS( 'N',n_basis ,  lead_basis(1), work, n_basis, ipiv, green(1:n_basis, 1:lead_basis(1)), n_basis, info )

             greenL(1:lead_basis(2),1:lead_basis(1)) = green(i_L2_1:i_L2_2,1:lead_basis(1))

             if (transport_k_points) then
                workL(1:lead_basis(2),1:lead_basis(2)) =  (0.0d0,1.0d0)*(lead_self_energy_2(1:lead_basis(2),1:lead_basis(2)) - &
                     dconjg(transpose(lead_self_energy_2(1:lead_basis(2),1:lead_basis(2)))))
             else
                workL(1:lead_basis(2),1:lead_basis(2)) =  2* aimag(lead_self_energy_2(1:lead_basis(2),1:lead_basis(2)))
             end if
                
          
             call ZGEMM('N','N',lead_basis(2),lead_basis(1),lead_basis(2),(1.d0,0.d0), workL, n_lead_basis, &
                  greenL, n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

             if (transport_k_points) then
                workL(1:lead_basis(1),1:lead_basis(1)) =  (0.0d0,1.0d0)*(lead_self_energy_1(1:lead_basis(1),1:lead_basis(1)) - &
                     dconjg(transpose(lead_self_energy_1(1:lead_basis(1),1:lead_basis(1)))))
             else
                workL(1:lead_basis(1),1:lead_basis(1))  = 2* aimag(lead_self_energy_1(1:lead_basis(1),1:lead_basis(1)))
             end if

             call ZGEMM('N','C',lead_basis(1),lead_basis(2),lead_basis(1),(1.d0,0.d0), workL, n_lead_basis, &
                  greenL, n_lead_basis, (0.d0,0.d0), ttt2,n_lead_basis)

             call ZGEMM('N','N',lead_basis(2),lead_basis(2),lead_basis(1),(1.d0,0.d0), ttt, n_lead_basis, ttt2, n_lead_basis, &
                  (0.d0,0.d0), workL,n_lead_basis)

             tun = 0.d0

             do i_basis_1 = 1, lead_basis(2)
                tun(1,myid+1) = tun(1,myid+1) + workL(i_basis_1, i_basis_1)
             end do

             ! In the case we have 4 leads we need to solve more.

             if(n_leads == 4)then

                i_L2_1 = lead_basis_start(3)
                i_L2_2 = lead_basis_start(3)+lead_basis(3)-1

                greenL(1:lead_basis(3),1:lead_basis(1)) = green(i_L2_1:i_L2_2,1:lead_basis(1))

                workL(1:lead_basis(3),1:lead_basis(3)) =  2* aimag(lead_self_energy_3(1:lead_basis(3),1:lead_basis(3)))
                
                call ZGEMM('N','N',lead_basis(3),lead_basis(1),lead_basis(3),(1.d0,0.d0), workL, n_lead_basis, &
                     greenL, n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

                workL(1:lead_basis(1),1:lead_basis(1))  = 2* aimag(lead_self_energy_1(1:lead_basis(1),1:lead_basis(1)))
                
                call ZGEMM('N','C',lead_basis(1),lead_basis(3),lead_basis(1),(1.d0,0.d0), workL, n_lead_basis, &
                  greenL, n_lead_basis, (0.d0,0.d0), ttt2,n_lead_basis)
                
                
                call ZGEMM('N','N',lead_basis(3),lead_basis(3),lead_basis(1),(1.d0,0.d0), ttt, n_lead_basis, ttt2, n_lead_basis, &
                     (0.d0,0.d0), workL,n_lead_basis)


                do i_basis_1 = 1, lead_basis(3)
                   tun(2,myid+1) = tun(2,myid+1) + workL(i_basis_1, i_basis_1)
                end do

                i_L2_1 = lead_basis_start(4)
                i_L2_2 = lead_basis_start(4)+lead_basis(4)-1

                greenL(1:lead_basis(4),1:lead_basis(1)) = green(i_L2_1:i_L2_2,1:lead_basis(1))

                workL(1:lead_basis(4),1:lead_basis(4))  = 2* aimag(lead_self_energy_4(1:lead_basis(4),1:lead_basis(4)))

                call ZGEMM('N','N',lead_basis(4),lead_basis(1),lead_basis(4),(1.d0,0.d0), workL, n_lead_basis, &
                     greenL, n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

                workL(1:lead_basis(1),1:lead_basis(1))  = 2* aimag(lead_self_energy_1(1:lead_basis(1),1:lead_basis(1)))

                call ZGEMM('N','C',lead_basis(1),lead_basis(4),lead_basis(1),(1.d0,0.d0), workL, n_lead_basis, &
                     greenL, n_lead_basis, (0.d0,0.d0), &
                  ttt2,n_lead_basis)

                call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(1),(1.d0,0.d0), ttt, n_lead_basis, ttt2, n_lead_basis, &
                     (0.d0,0.d0), workL,n_lead_basis)


                do i_basis_1 = 1, lead_basis(4)
                   tun(3,myid+1) = tun(3,myid+1) + workL(i_basis_1, i_basis_1)
                end do

             !----------------------------------

                i_L2_1 = lead_basis_start(4)
                i_L2_2 = lead_basis_start(4)+lead_basis(4)-1

                green = 0.d0
                do i_basis_1 = 1, lead_basis(3)
                   green( lead_basis_start(3)+i_basis_1-1, i_basis_1) = 1.d0
                end do

                call  ZGETRS( 'N',n_basis ,  lead_basis(3), work, n_basis, ipiv, green(1:n_basis, 1:lead_basis(3)), n_basis, info )
                
                greenL(1:lead_basis(3),1:lead_basis(3)) = green(i_L2_1:i_L2_2,1:lead_basis(3))                
                workL(1:lead_basis(4),1:lead_basis(4))  = 2* aimag(lead_self_energy_4(1:lead_basis(4),1:lead_basis(4)))

                call ZGEMM('N','N',lead_basis(4),lead_basis(3),lead_basis(4),(1.d0,0.d0), workL, n_lead_basis, greenL, &
                     n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

                workL(1:lead_basis(3),1:lead_basis(3))  = 2* aimag(lead_self_energy_3(1:lead_basis(3),1:lead_basis(3)))

                call ZGEMM('N','C',lead_basis(3),lead_basis(4),lead_basis(3),(1.d0,0.d0), workL, n_lead_basis, greenL, &
                     n_lead_basis, (0.d0,0.d0), ttt2,n_lead_basis)

                call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(3),(1.d0,0.d0), ttt, n_lead_basis, ttt2, n_lead_basis, &
                     (0.d0,0.d0), workL,n_lead_basis)

                do i_basis_1 = 1, lead_basis(4)
                   tun(4,myid+1) = tun(4,myid+1) + workL(i_basis_1, i_basis_1)
                end do

             end if ! n_leads==4

             energy_print = 0.d0
             energy_print(myid+1) = energy

             call sync_vector_complex(tun(1,1), n_tasks*4)
             call sync_vector(energy_print(1), n_tasks)

             ! Print out the results.

             if(myid==0)then
                do i_task = 1,n_tasks
                   if(n_leads == 2)then

                      if(i_spin == 1)then
                         write(22,*) & 
                           real(energy_print(i_task)-chemical_potential-delta_fermi_level )*Hartree, & 
                           real(tun(1,i_task))
                      else
                         write(23,*) & 
                           real(energy_print(i_task)-chemical_potential-delta_fermi_level )*Hartree, & 
                           real(tun(1,i_task))
                      end if
                   
                      write(use_unit,*) & 
                           real(energy_print(i_task)-chemical_potential-delta_fermi_level )*Hartree, & 
                           real(tun(1,i_task))

                   else if(n_leads==4)then
                   
                      if(i_spin == 1)then
                         write(22,'(5F14.6)') &
                              real(energy_print(i_task) -chemical_potential &
                              - delta_fermi_level) * Hartree, &
                              real(tun(1, i_task)), real(tun(2, i_task)), &
                              real(tun(3, i_task)), real(tun(4, i_task))
                      else
                         write(23,'(5F14.6)') &
                              real(energy_print(i_task) - chemical_potential &
                              - delta_fermi_level) * Hartree, & 
                              real(tun(1, i_task)), real(tun(2, i_task)), &
                              real(tun(3, i_task)), real(tun(4, i_task))
                      end if
                      write(use_unit,'(5F14.6)') &
                           real(energy_print(i_task) - chemical_potential &
                           - delta_fermi_level) * Hartree, &
                           real(tun(1, i_task)), real(tun(2, i_task)), &
                           real(tun(3, i_task)), real(tun(4, i_task))
                      
                   end if
                end do
             end if
          end do

          tun_prev = tun(1, 1)

       end do ! i_energy

    end do ! i_k_point

    close(22)
    deallocate(workL)
    deallocate(work)
    deallocate(green)
    deallocate(greenL)
    deallocate(ttt)       
    deallocate(ttt2)
    

  end subroutine read_lead_information
  !******-------------------------------------------------------------------------------------
  !****s* transport/read_lead_information_scalapack
  !  NAME
  !     read_lead_information_scalapack
  !  SYNOPSIS
  !    subroutine read_lead_information_scalapack( chemical_potential)
      
  subroutine read_lead_information_scalapack( chemical_potential)
     
    !  PURPOSE
    !   Evaluation electron tunneling propability. 
    !   Note that this is the scalapack vertion of the "read_lead_information"
    !  
    !  INPUTS
    !    o chemical_potential
    !  OUTPUT
    !    none
    !  SOURCE

    use basis, only: basis_atom, n_basis_atom
    use scalapack_wrapper
    use timing
    use localorb_io
    implicit none
    integer:: i_lead, i_basis_1, i_basis_2, n_lead_basis
    real*8,    allocatable, dimension(:,:):: lead_ovlp_matrix_1, lead_connection_matrix_o_1
    real*8,    allocatable, dimension(:,:,:):: lead_connection_matrix_h_1, lead_hamiltonian_1

    real*8,    allocatable, dimension(:,:):: lead_ovlp_matrix_2, lead_connection_matrix_o_2
    real*8,    allocatable, dimension(:,:,:):: lead_connection_matrix_h_2, lead_hamiltonian_2

    real*8,    allocatable, dimension(:,:):: lead_ovlp_matrix_3, lead_connection_matrix_o_3
    real*8,    allocatable, dimension(:,:,:):: lead_connection_matrix_h_3, lead_hamiltonian_3

    real*8,    allocatable, dimension(:,:):: lead_ovlp_matrix_4, lead_connection_matrix_o_4
    real*8,    allocatable, dimension(:,:,:):: lead_connection_matrix_h_4, lead_hamiltonian_4

    integer:: ipiv(n_basis)
    integer:: info, i_1, i_2, i_spin, i_cell, i_size, i_index_real, i_energy
    integer:: i_L1_1, i_L1_2, i_L2_1, i_L2_2, n_lead_basis_max, i_task, i_atom, i_basis,i_l


    complex*16,    allocatable, dimension(:,:):: workL, ttt, ttt2, greenL
    complex*16:: tun(0:(n_tasks-1),4)
    real*8:: chemical_potential, delta_fermi_level

    complex*16:: epsilon
    real*8:: energy(0:(n_tasks-1)),energy_lead,  lead_chem_pot(4)
    integer:: task_energy(0:n_tasks), task_lead(0:n_tasks), n_task_energy_steps
    integer:: send_12(1:n_tasks),send_13(1:n_tasks),send_14(1:n_tasks),send_34(1:n_tasks)
    integer:: receive_12(1:n_tasks), receive_13(1:n_tasks), receive_14(1:n_tasks), receive_34(1:n_tasks)

!    logical:: ldos = .false.
    real*8, dimension(:,:,:), allocatable :: mulliken_decomp
    character*40 :: proj_dos_filename
    character*100 :: info_str, outputformat
    character*100 :: outputformat_header
    real*8:: time_1a, time_1b, time_2a, time_2b, new_max
    integer:: i_matrix_element_non_zero


    call get_timestamps (time_1a, time_1b)



    if(myid==0) write(use_unit,*) 'Evaluating tunneling probability using scalapack'


    if( mod(n_tasks,n_leads) /=0)then
       write(use_unit,*) 'Error: scalapack parallerization only with n_task= N*n_leads'
       stop
    end if


    n_task_energy_steps = n_tasks/n_leads


    ! Parallerization:

    i_task = -1
    do i_energy = 1,n_task_energy_steps

       send_12(i_energy) = i_task + 1
       receive_12(i_energy) = i_task +2

       send_13(i_energy) = i_task + 1
       receive_13(i_energy) = i_task + 3

       send_14(i_energy) = i_task + 4
       receive_14(i_energy) = i_task + 1

       send_34(i_energy) = i_task + 3
       receive_34(i_energy) = i_task + 4

       do i_lead = 1, n_leads

          i_task = i_task + 1
          task_lead  ( i_task )  = i_lead
          task_energy( i_task ) = i_energy

       end do
    end do

    ! Where to put boundary conditons.
    do i_lead = 1, n_leads

       basis_check: do i_basis_1 = 1,n_basis
          if(Cbasis_to_atom(i_basis_1)== lead_atoms_start(i_lead))then
             lead_basis_start(i_lead) = i_basis_1
             exit basis_check
          end if
       end do basis_check
    end do

    n_lead_basis = maxval(lead_basis)


    ! Read in the boundary condition information from the file.

    do i_lead = 1, n_leads


       open(88, file=lead_file(i_lead))
       read(88,*) lead_atoms( i_lead ), lead_basis( i_lead)
       read(88,*) average_pot_lead(i_lead)
       read(88,*) lead_chem_pot(i_lead)
       read(88,*) i_matrix_element_non_zero



       if(myid==0) write(use_unit,'(2X,A)') ' '
       if(myid==0) write(use_unit,'(2X,A,I1,A)') 'Lead ',i_lead,   ' :'
       if(myid==0) write(use_unit,'(2X,A,I5,A,I5)') '| Including atoms:   ', & 
                        lead_atoms_start( i_lead), '  to',  lead_atoms_start( i_lead)+ lead_atoms( i_lead )-1
       if(myid==0) write(use_unit,'(2X,A,I5)') '| Starting basis:  ', lead_basis_start(i_lead)
       if(myid==0) write(use_unit,'(2X,A,I5)') '| Number of atoms: ', lead_atoms( i_lead )
       if(myid==0) write(use_unit,'(2X,A,I5)') '| Number of basis: ', lead_basis( i_lead)
       if(myid==0) write(use_unit,'(2X,A,F16.8)') '| Pot ave  : ',  average_pot_lead(i_lead)*Hartree
       if(myid==0) write(use_unit,'(2X,A,F16.8)') '| Lead chem. pot.  : ',  lead_chem_pot(i_lead)*Hartree


       lead_basis_end(i_lead) = lead_basis_start(i_lead) + lead_basis( i_lead) -1

       n_basis_atom = 0
       do i_basis = 1, n_basis
          n_basis_atom(basis_atom(i_basis)) =  n_basis_atom(basis_atom(i_basis)) + 1
       end do


       if(i_lead == task_lead(myid))then


          select case(i_lead)

          case(1)

             allocate(lead_connection_matrix_o_1(n_lead_basis, n_lead_basis))
             lead_connection_matrix_o_1 = 0.d0
             allocate(lead_connection_matrix_h_1(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_connection_matrix_h_1 = 0.d0
             allocate(lead_ovlp_matrix_1(n_lead_basis, n_lead_basis))
             lead_ovlp_matrix_1 = 0.d0
             allocate(lead_hamiltonian_1(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_hamiltonian_1 = 0.d0
             allocate(lead_self_energy_1(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_1 = 0.d0

             do i_1 = 1, i_matrix_element_non_zero
             
                read(88,*)  i_basis_1,  i_basis_2
                read(88,*) lead_hamiltonian_1(i_basis_1, i_basis_2,1), lead_ovlp_matrix_1(i_basis_1, i_basis_2), &
                     lead_connection_matrix_h_1(i_basis_1, i_basis_2,1), lead_connection_matrix_o_1(i_basis_1, i_basis_2)
             
             end do


          case(2)

             allocate(lead_connection_matrix_o_2(n_lead_basis, n_lead_basis))
             lead_connection_matrix_o_2 = 0.d0
             allocate(lead_connection_matrix_h_2(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_connection_matrix_h_2 = 0.d0
             allocate(lead_ovlp_matrix_2(n_lead_basis, n_lead_basis))
             lead_ovlp_matrix_2 = 0.d0
             allocate(lead_hamiltonian_2(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_hamiltonian_2 = 0.d0
             allocate(lead_self_energy_2(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_2  = 0.d0


             do i_1 = 1, i_matrix_element_non_zero
             
                read(88,*)  i_basis_1,  i_basis_2
                read(88,*) lead_hamiltonian_2(i_basis_1, i_basis_2,1), lead_ovlp_matrix_2(i_basis_1, i_basis_2), &
                     lead_connection_matrix_h_2(i_basis_1, i_basis_2,1), lead_connection_matrix_o_2(i_basis_1, i_basis_2)

             end do

          case(3)

             allocate(lead_connection_matrix_o_3(n_lead_basis, n_lead_basis))
             lead_connection_matrix_o_3 = 0.d0
             allocate(lead_connection_matrix_h_3(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_connection_matrix_h_3 = 0.d0
             allocate(lead_ovlp_matrix_3(n_lead_basis, n_lead_basis))
             lead_ovlp_matrix_3 = 0.d0
             allocate(lead_hamiltonian_3(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_hamiltonian_3 = 0.d0
             allocate(lead_self_energy_3(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_3 = 0.d0


             do i_1 = 1, i_matrix_element_non_zero
             
                read(88,*)  i_basis_1,  i_basis_2
                read(88,*) lead_hamiltonian_3(i_basis_1, i_basis_2,1), lead_ovlp_matrix_3(i_basis_1, i_basis_2), &
                     lead_connection_matrix_h_3(i_basis_1, i_basis_2,1), lead_connection_matrix_o_3(i_basis_1, i_basis_2)
             end do



          case(4)


             allocate(lead_connection_matrix_o_4(n_lead_basis, n_lead_basis))
             lead_connection_matrix_o_4 = 0.d0
             allocate(lead_connection_matrix_h_4(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_connection_matrix_h_4 = 0.d0
             allocate(lead_ovlp_matrix_4(n_lead_basis, n_lead_basis))
             lead_ovlp_matrix_4 = 0.d0
             allocate(lead_hamiltonian_4(n_lead_basis, n_lead_basis,n_spin_lead))
             lead_hamiltonian_4 = 0.d0
             allocate(lead_self_energy_4(lead_basis( i_lead),lead_basis( i_lead)))
             lead_self_energy_4  = 0.d0



             do i_1 = 1, i_matrix_element_non_zero
             
                read(88,*)  i_basis_1,  i_basis_2
                read(88,*) lead_hamiltonian_4(i_basis_1, i_basis_2,1), lead_ovlp_matrix_4(i_basis_1, i_basis_2), &
                     lead_connection_matrix_h_4(i_basis_1, i_basis_2,1), lead_connection_matrix_o_4(i_basis_1, i_basis_2)

             end do


          end select
       end if
       close(88)


       if(lead_atoms_start(i_lead)+ lead_atoms(i_lead)-1 >n_atoms)then
          write(use_unit,*) 'Error: lead',i_lead, &
             ' start atom must be wrong the final atom is larger than ', &
             'number of atoms in the system!' 
          stop
       end if

       i_basis = 0
       do i_atom =   lead_atoms_start(i_lead), lead_atoms_start(i_lead)+ lead_atoms(i_lead)-1, 1 
          i_basis = i_basis +  n_basis_atom(i_atom)
       end do

       if(i_basis /= lead_basis(i_lead))then
          write(use_unit,*) 'Error: the boundary file for lead',i_lead, 'do not match to the system'
          write(use_unit,*) 'file:',lead_basis(i_lead), 'system:', i_basis 
          stop
       end if

    end do ! n_leads


    do i_lead = 1, n_leads
       if(myid==0) &
          write(use_unit,*) ' Fermi level - Reference potential',i_lead, &
             ':', (chemical_potential- average_pot(i_lead))*hartree
       if(myid==0) &
          write(use_unit,*) ' Reference gate voltage: ',i_lead,':', &
             (-( -average_pot(i_lead) + chemical_potential)  &
             + (- average_pot_lead(i_lead)+  lead_chem_pot(i_lead)))*Hartree

    end do


    ! The file for the results
    open(22,file=tunne_file)


    ! Possibility to use gate.
    if(gate_level_lead /=0)then

       delta_fermi_level = (-( -average_pot(gate_level_lead) + chemical_potential)  &
            + (- average_pot_lead(gate_level_lead)+  lead_chem_pot(gate_level_lead))) 

       if(myid==0) write(use_unit,*) '| delta Fermi level: ', delta_fermi_level*Hartree

    else
       delta_fermi_level = 0.0d0
    end if


    ! The "gate" potential of the system
    do i_lead = 1, n_leads
        if(myid==0) write(use_unit,*) ' Reference gate voltage after fix: ',i_lead,':', ( &
            -( -average_pot(i_lead) + chemical_potential + delta_fermi_level )  &
          + (- average_pot_lead(i_lead)+  lead_chem_pot(i_lead)))*Hartree

        if(myid==0) write(22,'(A,I2,2F16.8)') '# Reference gate: ', i_lead,& 
          ( -( -average_pot(i_lead) + chemical_potential + delta_fermi_level )  &
          + (- average_pot_lead(i_lead)+  lead_chem_pot(i_lead)))*Hartree, delta_fermi_level

    end do

    ! If the Fermi levels fixed.
    if(fermi_level_fix)then
       do i_lead = 1, n_leads

          average_pot_lead(i_lead)= lead_chem_pot(i_lead)
          average_pot(i_lead) = chemical_potential
          
       end do
    end if

    allocate(workL(n_lead_basis, n_lead_basis))
    workL = 0.d0
    allocate(ttt(n_lead_basis, n_lead_basis))
    ttt = 0.d0
    allocate(ttt2(n_lead_basis, n_lead_basis))
    ttt2 = 0.d0
    allocate(greenL(n_lead_basis, n_lead_basis))
    greenL = 0.d0
 
    if(0.02 < epsilon_end) epsilon_start =  (0,1)*epsilon_end

    ! The loop over energy.

    do i_energy = 0,n_energy_steps-1, n_task_energy_steps

       do i_task  = 1, n_task_energy_steps
          energy(i_task) = chemical_potential+delta_fermi_level + & 
                           energy_start  + (energy_end-energy_start)/ (n_energy_steps-1) * (i_energy + i_task-1 )
       end do

       ! Initialize the Green's function of the boundaries.

       do i_lead = 1, n_leads
          if(i_lead == task_lead(myid))then

             energy_lead = energy(task_energy(myid)) + average_pot_lead(i_lead) - average_pot(i_lead) 

             epsilon = epsilon_start

             select case(i_lead)
             case(1)
                workL = ( energy_lead + epsilon)* lead_ovlp_matrix_1 - lead_hamiltonian_1(:,:,1)
             case(2)
                workL = ( energy_lead + epsilon)* lead_ovlp_matrix_2 - lead_hamiltonian_2(:,:,1)
             case(3)
                workL = ( energy_lead + epsilon)* lead_ovlp_matrix_3 - lead_hamiltonian_3(:,:,1)
             case(4)
                workL = ( energy_lead + epsilon)* lead_ovlp_matrix_4 - lead_hamiltonian_4(:,:,1)
             end select


             call ZGETRF(n_lead_basis ,n_lead_basis, workL, n_lead_basis, ipiv, info )

             greenL = 0.d0
             do i_basis_1 = 1, lead_basis(i_lead)
                greenL(i_basis_1, i_basis_1) = 1.d0
             end do

             call  ZGETRS( 'N',lead_basis(i_lead) , lead_basis(i_lead), &
                           workL, n_lead_basis, ipiv, greenL, n_lead_basis,&
                           info )



             select case(i_lead)
             case(1)

                ttt(1:lead_basis(1),1:lead_basis(1)) = &
                   real(energy_lead) * &
                   lead_connection_matrix_o_1(1:lead_basis(1),1:lead_basis(1)) &
                 - lead_connection_matrix_h_1(1:lead_basis(1),1:lead_basis(1),1)

             case(2)

                ttt = real(energy_lead) * lead_connection_matrix_o_2 &
                   - lead_connection_matrix_h_2(:,:,1)

             case(3)

                ttt(1:lead_basis(3),1:lead_basis(3)) = & 
                   real(energy_lead) * &
                   lead_connection_matrix_o_3(1:lead_basis(3),1:lead_basis(3)) &
                 - lead_connection_matrix_h_3(1:lead_basis(3),1:lead_basis(3),1)

             case(4)

                ttt(1:lead_basis(4),1:lead_basis(4)) = & 
                   real(energy_lead) * &
                   lead_connection_matrix_o_4(1:lead_basis(4),1:lead_basis(4)) &
                 - lead_connection_matrix_h_4(1:lead_basis(4),1:lead_basis(4),1)

             end select

             ! Iterate the Green's function of the boundaries

             boundary_scalapack: do i_1 = 1, n_boundary_iterations

                epsilon = epsilon_end*(0,1) &
                   + (epsilon_start - epsilon_end*(0,1)) &
                   * (n_boundary_iterations-i_1) &
                   / real(n_boundary_iterations-1)**2
                ttt2 = greenL

                select case(i_lead)
                case(1)

                   call ZGEMM('T','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),ttt, n_lead_basis,  &
                        greenL,n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                   call ZGEMM('N','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),workL, n_lead_basis, &
                        ttt,n_lead_basis, (0.d0,0.d0), greenL,n_lead_basis)

                   workL(1:lead_basis(1),1:lead_basis(1)) = & 
                        (energy_lead + epsilon) * lead_ovlp_matrix_1 - lead_hamiltonian_1(:,:,1) &
                        - greenL(1:lead_basis(1),1:lead_basis(1))

                case(2)

                   call ZGEMM('T','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),ttt, n_lead_basis, &
                        greenL,n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                   call ZGEMM('N','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),workL, n_lead_basis, &
                        ttt,n_lead_basis, (0.d0,0.d0), greenL,n_lead_basis)

                   workL(1:lead_basis(2),1:lead_basis(2)) = & 
                        (energy_lead + epsilon) * lead_ovlp_matrix_2 - lead_hamiltonian_2(:,:,1) &
                        - greenL(1:lead_basis(2),1:lead_basis(2))

                case(3)

                   call ZGEMM('T','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),ttt, n_lead_basis, &
                        greenL,n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                   call ZGEMM('N','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),workL, n_lead_basis, &
                        ttt,n_lead_basis, (0.d0,0.d0), greenL,n_lead_basis)

                   workL(1:lead_basis(3),1:lead_basis(3)) = & 
                        (energy_lead + epsilon) * lead_ovlp_matrix_3 - lead_hamiltonian_3(:,:,1) &
                        - greenL(1:lead_basis(3),1:lead_basis(3))

                case(4)

                   call ZGEMM('T','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),ttt, n_lead_basis, &
                        greenL,n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                   call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),workL, n_lead_basis, &
                        ttt,n_lead_basis, (0.d0,0.d0), greenL,n_lead_basis)

                   workL(1:lead_basis(4),1:lead_basis(4)) = & 
                        (energy_lead + epsilon) * lead_ovlp_matrix_4 - lead_hamiltonian_4(:,:,1) &
                        - greenL(1:lead_basis(4),1:lead_basis(4))

                end select

                call ZGETRF(lead_basis(i_lead) ,lead_basis(i_lead), workL, n_lead_basis, ipiv, info )

                greenL = 0.d0
                do i_basis_1 = 1, lead_basis(i_lead)
                   greenL(i_basis_1, i_basis_1) = 1.d0
                end do

                call  ZGETRS( 'N',lead_basis(i_lead) , lead_basis(i_lead), & 
                              workL, n_lead_basis, ipiv, greenL, n_lead_basis, info )

                new_max = &
                  maxval(abs(greenL(1:lead_basis(i_lead),1:lead_basis(i_lead))&
                             - ttt2(1:lead_basis(i_lead),1:lead_basis(i_lead))))

                ! Is the boundary condition converged ?
                if(abs(new_max)< boundary_treshold)then
                   exit boundary_scalapack
                end if
                ! Mixing to the previous solution
                greenL = greenL*boundary_mix + ttt2* ( 1.d0 - boundary_mix)

             end do boundary_scalapack

             ! Final steps to boundary condition calculation.

             select case(i_lead)
             case(1)

                call ZGEMM('T','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),ttt, n_lead_basis, greenL, &
                     n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                call ZGEMM('N','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),workL, n_lead_basis, ttt, &
                     n_lead_basis, (0.d0,0.d0), lead_self_energy_1,lead_basis(1))

             case(2)

                call ZGEMM('T','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),ttt, n_lead_basis, greenL, &
                     n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                call ZGEMM('N','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),workL, n_lead_basis, ttt, &
                     n_lead_basis, (0.d0,0.d0), lead_self_energy_2,lead_basis(2))

             case(3)

                call ZGEMM('T','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),ttt, n_lead_basis, greenL, &
                     n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                call ZGEMM('N','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),workL, n_lead_basis, ttt, &
                     n_lead_basis, (0.d0,0.d0), lead_self_energy_3,lead_basis(3))

             case(4)

                call ZGEMM('T','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),ttt, n_lead_basis, greenL, &
                     n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),workL, n_lead_basis, ttt, &
                     n_lead_basis, (0.d0,0.d0), lead_self_energy_4,lead_basis(4))

             end select

          end if
       end do ! i_lead

       ! Initialize and solve the Green's function calculation to center region.
       do i_spin = 1, n_spin

          do i_task = 1, n_task_energy_steps,1
             
             call construct_greenfunction_scalapack( energy(i_task) )

             do i_lead = 1, n_leads

                ttt = 0.d0
                
                if(i_lead == task_lead(myid) .and. i_task==task_energy(myid))then
                   select case(i_lead)
                   case(1)
                      ttt = lead_self_energy_1
                   case(2)
                      ttt = lead_self_energy_2
                   case(3)
                      ttt = lead_self_energy_3
                   case(4)
                      ttt = lead_self_energy_4
                   end select
                end if

                call bcast_complex_vector( ttt,n_lead_basis*n_lead_basis, (i_task-1)*n_leads + (i_lead-1) )

                i_1 = lead_basis_start(i_lead)
                i_2 = lead_basis_start(i_lead)+lead_basis(i_lead)-1

                call add_self_energy_to_greenfunction_scalapack & 
                ( ttt, lead_basis_start(i_lead), lead_basis_end(i_lead), n_lead_basis)

             end do

             call solve_greens_functions(  greenL,& 
                  lead_basis_start(1), lead_basis_end(1), lead_basis_start(2), lead_basis_end(2), &
                  lead_basis_start(3), lead_basis_end(3), lead_basis_start(4), lead_basis_end(4), &
                  lead_basis(1), lead_basis(2), lead_basis(3), lead_basis(4),n_lead_basis,  &
                  n_leads , ttt, n_task_energy_steps, i_task, &
                  receive_12, receive_13, receive_14, receive_34)
             
          end do ! task

          ! Move to the electron transport calculation.
          ! Unfortunately the previous boundary condition can no be used as a starting point to the next one
          ! because it is more important to save memory.

          do i_lead = 1, n_leads

             if(i_lead == task_lead(myid))then

                i_1 = lead_basis_start(i_lead)
                i_2 = lead_basis_start(i_lead)+lead_basis(i_lead)-1

                select case(i_lead)
                case(1)
                   
                   do i_basis_1 = 1,lead_basis(1)
                      do i_basis_2 = 1,lead_basis(1)
                         
                         lead_self_energy_1(i_basis_1, i_basis_2)  = 2 * aimag(lead_self_energy_1(i_basis_1, i_basis_2))
                         
                      end do
                   end do
                   
                case(2)

                   do i_basis_1 = 1,lead_basis(2)
                      do i_basis_2 = 1,lead_basis(2)
                         
                         lead_self_energy_2(i_basis_1, i_basis_2)  = 2 * aimag(lead_self_energy_2(i_basis_1, i_basis_2))

                      end do
                   end do
                   
                case(3)

                   do i_basis_1 = 1,lead_basis(3)
                      do i_basis_2 = 1,lead_basis(3)

                         lead_self_energy_3(i_basis_1, i_basis_2)  = 2 * aimag(lead_self_energy_3(i_basis_1, i_basis_2))

                      end do
                   end do

                case(4)

                   do i_basis_1 = 1,lead_basis(4)
                      do i_basis_2 = 1,lead_basis(4)
                         
                         lead_self_energy_4(i_basis_1, i_basis_2)  = 2 * aimag(lead_self_energy_4(i_basis_1, i_basis_2))

                      end do
                   end do

                end select
             end if
          end do

          tun = 0.d0
          
          if(myid== send_12(task_energy(myid)))then
             call send_complex_vector(lead_self_energy_1, n_lead_basis*n_lead_basis, receive_12(task_energy(myid))) 
          end if
          if(myid==receive_12(task_energy(myid)))then
             call receive_complex_vector(ttt,  n_lead_basis*n_lead_basis,  send_12(task_energy(myid))) 
          end if
       
          if(n_leads >=3 )then
             if(myid== send_13(task_energy(myid)))then
                call send_complex_vector(lead_self_energy_1, n_lead_basis*n_lead_basis, receive_13(task_energy(myid))) 
             end if
             if(myid==receive_13(task_energy(myid)))then
                call receive_complex_vector(ttt,  n_lead_basis*n_lead_basis,  send_13(task_energy(myid))) 
             end if
          end if

          if(n_leads >=4 )then
             if(myid== send_14(task_energy(myid)))then
                call send_complex_vector(lead_self_energy_4, n_lead_basis*n_lead_basis, receive_14(task_energy(myid))) 
             end if
             if(myid==receive_14(task_energy(myid)))then
                call receive_complex_vector(workL,  n_lead_basis*n_lead_basis,  send_14(task_energy(myid))) 
             end if

             if(myid== send_34(task_energy(myid)))then
                call send_complex_vector(lead_self_energy_3, n_lead_basis*n_lead_basis, receive_34(task_energy(myid))) 
             end if
             if(myid==receive_34(task_energy(myid)))then
                call receive_complex_vector(ttt,  n_lead_basis*n_lead_basis,  send_34(task_energy(myid))) 
             end if
          end if

          ! The final tunneling evaluation
          
          if(myid==receive_12(task_energy(myid)))then

             call ZGEMM('N','C',lead_basis(1),lead_basis(2),lead_basis(1),(1.d0,0.d0), ttt, lead_basis(1), greenL,  &
                  n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

             call ZGEMM('N','N',lead_basis(2),lead_basis(1),lead_basis(2),(1.d0,0.d0), lead_self_energy_2, lead_basis(2), &
                  greenL, n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

             call ZGEMM('N','N',lead_basis(2),lead_basis(2),lead_basis(1),(1.d0,0.d0), ttt, n_lead_basis, workL, &
                  n_lead_basis, (0.d0,0.d0), greenL,n_lead_basis)

             do i_basis_1 = 1, lead_basis(2)
                tun(task_energy(myid),1) = tun(task_energy(myid),1) + greenL(i_basis_1, i_basis_1)
             end do
          end if

          if(n_leads >=3 )then

             if(myid==receive_13(task_energy(myid)))then
                
                call ZGEMM('N','C',lead_basis(1),lead_basis(3),lead_basis(1),(1.d0,0.d0), ttt, lead_basis(1), greenL, &
                     n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                call ZGEMM('N','N',lead_basis(3),lead_basis(1),lead_basis(3),(1.d0,0.d0), lead_self_energy_3, lead_basis(3), &
                     greenL, n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

                call ZGEMM('N','N',lead_basis(3),lead_basis(3),lead_basis(1),(1.d0,0.d0), ttt, n_lead_basis, workL,  &
                     n_lead_basis, (0.d0,0.d0), greenL,n_lead_basis)

                do i_basis_1 = 1, lead_basis(3)
                   tun(task_energy(myid),2) = tun(task_energy(myid),2) + greenL(i_basis_1, i_basis_1)
                end do
             end if
          end if

          if(n_leads >=4 )then

             if(myid==receive_14(task_energy(myid)))then

                call ZGEMM('N','N',lead_basis(4),lead_basis(1),lead_basis(4),(1.d0,0.d0), workL, lead_basis(4), greenL, &
                     n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

                call ZGEMM('N','C',lead_basis(1),lead_basis(4),lead_basis(1),(1.d0,0.d0), lead_self_energy_1, lead_basis(1), &
                     greenL, n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(1),(1.d0,0.d0), ttt, n_lead_basis, workL, &
                     n_lead_basis, (0.d0,0.d0), greenL, n_lead_basis)

                do i_basis_1 = 1, lead_basis(4)
                   tun(task_energy(myid),3) = tun(task_energy(myid),3) + greenL(i_basis_1, i_basis_1)
                end do

             end if

             if(myid==receive_34(task_energy(myid)))then

                call ZGEMM('N','C',lead_basis(3),lead_basis(4),lead_basis(3),(1.d0,0.d0), ttt, lead_basis(3), &
                     greenL, n_lead_basis, (0.d0,0.d0), workL,n_lead_basis)

                call ZGEMM('N','N',lead_basis(4),lead_basis(3),lead_basis(4),(1.d0,0.d0), lead_self_energy_4, &
                     lead_basis(4), greenL, n_lead_basis, (0.d0,0.d0), ttt,n_lead_basis)

                call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(3),(1.d0,0.d0), ttt, n_lead_basis, &
                     workL, n_lead_basis, (0.d0,0.d0), greenL,n_lead_basis)

                do i_basis_1 = 1, lead_basis(4)
                   tun(task_energy(myid),4) = tun(task_energy(myid),4) + greenL(i_basis_1, i_basis_1)
                end do
             
             end if
          end if

          call sync_vector_complex(tun(0,1), n_tasks*4)
          
          ! Printout the results

          if(myid==0)then
             if(n_leads == 2)then
                do i_task = 1, n_task_energy_steps
                   write(22,*) (energy(i_task)-chemical_potential-delta_fermi_level )*Hartree, real(tun(i_task,1))
                   write(use_unit,*)  (energy(i_task)-chemical_potential-delta_fermi_level )*Hartree, real(tun(i_task,1))
                end do
             else if(n_leads==4)then
                do i_task = 1, n_task_energy_steps
                   
                   write(22,'(5F14.6)') (energy(i_task)-chemical_potential-delta_fermi_level)*Hartree, &
                        real(tun(i_task,1)), real(tun(i_task,2)), real(tun(i_task,3)),real(tun(i_task,4))
                   
                   write(use_unit,'(5F14.6)') (energy(i_task)-chemical_potential-delta_fermi_level)*Hartree, &
                        real(tun(i_task,1)), real(tun(i_task,2)), real(tun(i_task,3)),real(tun(i_task,4))

                end do
             end if
          end if

       end do ! spin
    end do ! energy

    close(22)

    if(allocated( lead_ovlp_matrix_1          )) deallocate( lead_ovlp_matrix_1         )
    if(allocated( lead_connection_matrix_o_1  )) deallocate( lead_connection_matrix_o_1 )
    if(allocated( lead_connection_matrix_h_1  )) deallocate( lead_connection_matrix_h_1 )
    if(allocated( lead_hamiltonian_1          )) deallocate( lead_hamiltonian_1         )

    if(allocated( lead_ovlp_matrix_2          )) deallocate( lead_ovlp_matrix_2         )
    if(allocated( lead_connection_matrix_o_2  )) deallocate( lead_connection_matrix_o_2 )
    if(allocated( lead_connection_matrix_h_2  )) deallocate( lead_connection_matrix_h_2 )
    if(allocated( lead_hamiltonian_2          )) deallocate( lead_hamiltonian_2         )

    if(allocated( lead_ovlp_matrix_3          )) deallocate( lead_ovlp_matrix_3         )
    if(allocated( lead_connection_matrix_o_3  )) deallocate( lead_connection_matrix_o_3 )
    if(allocated( lead_connection_matrix_h_3  )) deallocate( lead_connection_matrix_h_3 )
    if(allocated( lead_hamiltonian_3          )) deallocate( lead_hamiltonian_3         )

    if(allocated( lead_ovlp_matrix_4          )) deallocate( lead_ovlp_matrix_4         )
    if(allocated( lead_connection_matrix_o_4  )) deallocate( lead_connection_matrix_o_4 )
    if(allocated( lead_connection_matrix_h_4  )) deallocate( lead_connection_matrix_h_4 )
    if(allocated( lead_hamiltonian_4          )) deallocate( lead_hamiltonian_4         )

    if(allocated( workL                       )) deallocate( workL                      )
    if(allocated( ttt                         )) deallocate( ttt                        )

    if(allocated( mulliken_decomp             )) deallocate( mulliken_decomp            )

    if(allocated( greenL                      )) deallocate( greenL                     )


    call get_timestamps (time_2a, time_2b)
    if(myid==0) write(use_unit,*) 'Time for transport:',time_2a-time_1a , time_2b-time_1b  


  end subroutine read_lead_information_scalapack

  !******-------------------------------------------------------------------------------------
  !****s* transport/transport_find_local_min_occupated_state 
  !  NAME
  !     transport_find_local_min_occupated_state
  !  SYNOPSIS
  !      subroutine transport_find_local_min_occupated_state(occ_numbers, KS_eigenvalue,KS_eigenvector, &
  !                      KS_eigenvector_complex, overlap_matrix )

  subroutine transport_find_local_min_occupated_state(occ_numbers, KS_eigenvalue,KS_eigenvector, &
       KS_eigenvector_complex, overlap_matrix )
           
    !  PURPOSE
    !    Evatuate the reference potential, so that the leads match to the center region without potential jumps.
    !    The reference potential is the average of smallest eigenstates projected to atoms in different regions 
    !    (leads and the center region).
    !  
    !  INPUTS
    !   o occ_numbers 
    !   o KS_eigenvalue
    !   o KS_eigenvector
    !   o KS_eigenvector_complex,
    !   o overlap_matrix 
    !
    !  OUTPUT
    !    none
    !  SOURCE


    use scalapack_wrapper
    use localorb_io
    implicit none
    real*8:: overlap_matrix( n_hamiltonian_matrix_size )
    real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector_complex
    real*8,     dimension(n_states, n_spin, n_k_points)             :: occ_numbers
    real*8,     dimension(n_states, n_spin,n_k_points) :: KS_eigenvalue

    real*8:: projected_weight_r(n_atoms), projected_weight, min_eigenvalue(n_atoms,n_k_points)
    integer:: region(n_atoms),  atoms_in_region(0:4)

    integer:: i_lead, i_basis_1, i_basis_2, i_cell, i_state, i_k, i_k_point, i_index, i_spin, i_atom
    real*8:: min_states(n_atoms,n_k_points)
  

    if(myid==0)  write(use_unit,*) '---------------------------------------------------'
    if(myid==0)  write(use_unit,*) '  Electron transport calculation begin               '
    if(myid==0)  write(use_unit,*) 
    if(myid==0)  write(use_unit,*) 'Evaluating the reference potential level'


    if(transport_lead_calculation)then

       ! In the lead calculatin the whole system is center region
       region  = 0
       
    else 
       
       ! The center region is 0.
       region = 0


       do i_lead = 1, n_leads

          open(88, file=lead_file(i_lead))
          read(88,*) lead_atoms( i_lead ), lead_basis( i_lead)
          close(88)


          do i_atom = 1,n_atoms

             if( i_atom >= lead_atoms_start(i_lead) .and. i_atom <=  lead_atoms_start(i_lead)+ lead_atoms(i_lead)-1 )then

                region(i_atom) = i_lead
             end if

          end do
       end do
       lead_atoms(0) = n_atoms - sum(lead_atoms(1:4))

    end if

    min_eigenvalue = 0.0
    min_states = 0
 
    if(use_scalapack)then

       ! The scalapack version of project weights is in scalapack_wrapper

       call transport_proj_weight(  min_eigenvalue,  min_states, overlap_matrix, KS_eigenvalue )

    else

       do i_spin = 1, n_spin
          do i_state = 1, n_states , 1
             
             i_k = 0
             do i_k_point = 1, n_k_points,1

                projected_weight_r = 0.d0

                if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid<= n_k_points ) then
                   i_k = i_k + 1

                   do i_cell = 1,n_cells_in_hamiltonian-1
                      
                      do i_basis_2 = 1, n_basis, 1

                         do i_index = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
                            
                            i_basis_1 =  column_index_hamiltonian(i_index)
                            
                            if(real_eigenvectors)then

                               projected_weight_r(Cbasis_to_atom(i_basis_1)) = projected_weight_r(Cbasis_to_atom(i_basis_1)) &
                                    + KS_eigenvector(i_basis_1,i_state,i_spin,i_k)  &
                                    * KS_eigenvector(i_basis_2,i_state,i_spin,i_k)  &
                                    * overlap_matrix(i_index) &
                                    * dble(k_phase( i_cell,i_k_point))

                               
                               if(i_basis_1.ne.i_basis_2) then
                                  
                                  projected_weight_r(Cbasis_to_atom(i_basis_2)) = projected_weight_r(Cbasis_to_atom(i_basis_2)) &
                                       + KS_eigenvector(i_basis_1,i_state,i_spin,i_k)  &
                                       * KS_eigenvector(i_basis_2,i_state,i_spin,i_k)  &
                                       * overlap_matrix(i_index) &
                                       * dble(k_phase( i_cell,i_k_point))
                               end if
                               
                            else

                               projected_weight_r(Cbasis_to_atom(i_basis_1)) = projected_weight_r(Cbasis_to_atom(i_basis_1)) &
                                    + real(KS_eigenvector_complex(i_basis_1,i_state,i_spin,i_k)  &
                                    * conjg(KS_eigenvector_complex(i_basis_2,i_state,i_spin,i_k))  &
                                    * overlap_matrix(i_index) &
                                    * conjg(k_phase( i_cell,i_k_point) ))

                               if(i_basis_1.ne.i_basis_2) then
                                  
                                  projected_weight_r(Cbasis_to_atom(i_basis_2)) = projected_weight_r(Cbasis_to_atom(i_basis_2)) &
                                       + real(KS_eigenvector_complex(i_basis_2,i_state,i_spin,i_k)  &
                                       * conjg(KS_eigenvector_complex(i_basis_1,i_state,i_spin,i_k))  &
                                       * overlap_matrix(i_index) &
                                       * k_phase( i_cell,i_k_point) ) 
                                  
                               end if
                               
                            end if ! real_eigenvalues
                         end do ! i_index
                      end do ! i_basis_2
                   end do! i_cell


                   ! Next check is the eigenstate smallest and which atoms it belongs.                
                   ! Note that the eigenstate can belong to the several atoms. Here we use
                   ! weighted average over the lowest eigenvalues.
                   projected_weight_r = abs(projected_weight_r)*k_weights(i_k_point)

                   do i_atom = 1, n_atoms          

                      if( projected_weight_r(i_atom)  > 1e-4)then
                      
                         if( min_states(i_atom,i_k_point) < 5e-5) then 
                            ! The first eigenstate

                            min_states(i_atom,i_k_point) = projected_weight_r(i_atom)
                            min_eigenvalue(i_atom,i_k_point) = KS_eigenvalue(i_state,i_spin,i_k_point)*projected_weight_r(i_atom)

                      
                         else if( abs(min_eigenvalue(i_atom,i_k_point)/ min_states(i_atom,i_k_point) &
                              -  KS_eigenvalue(i_state,i_spin,i_k_point)) < 2)then
                            ! The eigenstate is close enough to the lowest one. Let's add it on.


                            min_states(i_atom,i_k_point) =  min_states(i_atom,i_k_point) + projected_weight_r(i_atom)
                            min_eigenvalue(i_atom,i_k_point) =  min_eigenvalue(i_atom,i_k_point) &
                                 +  KS_eigenvalue(i_state,i_spin,i_k_point)*projected_weight_r(i_atom)

                         else if( min_eigenvalue(i_atom,i_k_point)/ min_states(i_atom,i_k_point) &
                              >  KS_eigenvalue(i_state,i_spin,i_k_point)) then
                            ! We found the lower eigenstate.

                            min_states(i_atom,i_k_point) = projected_weight_r(i_atom)
                            min_eigenvalue(i_atom,i_k_point) = KS_eigenvalue(i_state,i_spin,i_k_point)*projected_weight_r(i_atom)
                         end if

                      end if
                   end do
                end if ! myid- k_point
             end do ! i_k_point
          end do ! i_state
       end do ! spin
    end if ! use_scalapack
    
    call sync_vector( min_eigenvalue, n_atoms*n_k_points)
    call sync_vector( min_states, n_atoms*n_k_points)

    do i_atom = 1, n_atoms
       do i_k_point = 1, n_k_points
          
          min_eigenvalue(i_atom,i_k_point) =  min_eigenvalue(i_atom,i_k_point)/ min_states(i_atom,i_k_point)

       end do
    end do
 
    do i_atom = 1, n_atoms
       min_eigenvalue(i_atom,1) =   sum(min_eigenvalue(i_atom,1:n_k_points))/n_k_points
    end do

    ! If the user wants, print out the minimum eigenvalues.
    if(out_min_eigenvalue)then
       if(myid==0)then
          open(889,file='min_eigenvalue.dat')
          do i_atom = 1, n_atoms
             write(889,*) i_atom,  min_eigenvalue(i_atom,1)*hartree
          end do
          close(889)
       end if
    end if

    ! Finally average over regions: leads and the center region.

    atoms_in_region = 0
    average_pot = 0.d0
    do i_atom = 1, n_atoms

       atoms_in_region(region(i_atom))   = atoms_in_region(region(i_atom)) + 1
       average_pot(region(i_atom)) =   average_pot(region(i_atom)) +  min_eigenvalue(i_atom,1)

    end do


    do i_lead = 0, n_leads

       average_pot(i_lead) = average_pot(i_lead) / atoms_in_region(i_lead)

       if(myid==0) write(use_unit,'(A,I3,F16.8)') ' | Average potential in lead:',i_lead,   average_pot(i_lead)*hartree

    end do


    average_pot(:) =   average_pot(:) + delta_pot_level
    if(myid==0) write(use_unit,*) '  | Average extra potential level ',  delta_pot_level
    if(myid==0) write(use_unit,*) 

  end subroutine transport_find_local_min_occupated_state
  !******-------------------------------------------------------------------------------------
  !****s* transport/transport_init_k_points
  !  NAME
  !      transport_init_k_points
  !  SYNOPSIS
  !      subroutine transport_init_k_points()

  subroutine transport_init_k_points()
           
    !  PURPOSE
    !    Initializes the k-point strucutres for transport calculations over intefaces. No k-points in the
    !    transport direction is employed.
    !  USES
    use constants, only: pi
    implicit none
    !  INPUTS
    !    none
    !
    !  OUTPUT
    !    none
    !
    !  SOURCE

    integer :: i_basis, i_k_point, i_x, i_y
    real*8 :: r_x, r_y
    
    allocate(k_phase_base_tr(3, n_k_points_tr))
    allocate(k_weights_tr(n_k_points_tr))

    i_k_point = 0
    do i_x = 1, n_k_points_xyz_tr(1)
       do i_y = 1, n_k_points_xyz_tr(2)
          
          i_k_point = i_k_point + 1

          r_x = dble(i_x-1) / dble(n_k_points_xyz_tr(1))
          r_y = dble(i_y-1) / dble(n_k_points_xyz_tr(2))

          k_phase_base_tr(1,i_k_point) = exp((0.0d0,1.0d0)*2*pi*r_x)
          k_phase_base_tr(2,i_k_point) = exp((0.0d0,1.0d0)*2*pi*r_y)
          k_phase_base_tr(3,i_k_point) = 1.0d0
       end do
    end do

    k_weights_tr(:) = 1.0d0/(n_k_points_xyz_tr(1)*n_k_points_xyz_tr(2))

  end subroutine transport_init_k_points

  !****s* transport/construct_lead_matrices
  !  NAME
  !     construct_lead_matrices
  !  SYNOPSIS
  !      subroutine construct_lead_matrices(lead_1, lead_2, lead_hamiltonian_cmplx, lead_ovlp_matrix_cmplx, &
  !             lead_connection_matrix_h_cmplx, lead_connection_matrix_o_cmplx, &
  !             lme1, lme2, n_lead_basis, i_k)

  subroutine construct_lead_matrices(lead_1, lead_2, lead_hamiltonian_cmplx, lead_ovlp_matrix_cmplx, &
       lead_connection_matrix_h_cmplx, lead_connection_matrix_o_cmplx, &
       lme1, lme2, n_lead_basis, i_k)
           
    !  PURPOSE
    !    Constructs the matrices for calculating the self energies of the leads when k-points are used in
    !    an interface transport calculation
    !  
    !  INPUTS
    !    o lead_1, lead_2 -- data structures containing the information on the leads
    !    o lme1, lme2 -- dimensions of lead_1 and lead_2
    !    o n_lead_basis -- number of basis functions in the lead regions
    !    o i_k -- the current k-point
    !
    !  OUTPUT
    !    o lead_***_matrix_n_cmplx -- complex matrices describing the lead, used to calculate the self-energies
    !                                 of the leads
    !
    !  SOURCE

    complex*16, dimension(n_lead_basis,n_lead_basis,n_spin_lead,n_leads) :: lead_hamiltonian_cmplx
    complex*16, dimension(n_lead_basis,n_lead_basis,n_leads) :: lead_ovlp_matrix_cmplx

    complex*16, dimension(n_lead_basis,n_lead_basis,n_spin_lead,n_leads):: lead_connection_matrix_h_cmplx
    complex*16, dimension(n_lead_basis,n_lead_basis,n_leads) :: lead_connection_matrix_o_cmplx
    
    type(matrix_entry), dimension(lme1) :: lead_1
    type(matrix_entry), dimension(lme2) :: lead_2

    integer :: lme1, lme2, n_lead_basis, i_k, i_lead
    integer :: i_entry, i_basis_1, i_basis_2
    integer :: i_cell_1, i_cell_2

    complex*16 :: k_phase_tr

    ! First, the first lead
    i_lead = 1
    do i_entry = 1, lme1

       i_basis_1 = lead_1(i_entry)%i_basis_1
       i_basis_2 = lead_1(i_entry)%i_basis_2
       i_cell_1 = lead_1(i_entry)%cell_1
       i_cell_2 = lead_1(i_entry)%cell_2

       k_phase_tr = (k_phase_base_tr(1,i_k)**i_cell_1)*(k_phase_base_tr(2,i_k)**i_cell_2)

       select case(lead_1(i_entry)%cell_3)

       case(0)

          lead_hamiltonian_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) = &
               lead_hamiltonian_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) + &
               k_phase_tr*lead_1(i_entry)%ham(n_spin_lead)

          lead_ovlp_matrix_cmplx(i_basis_1,i_basis_2,i_lead) = &
               lead_ovlp_matrix_cmplx(i_basis_1,i_basis_2,i_lead) + &
               k_phase_tr*lead_1(i_entry)%ovlp

          if (i_basis_1 /= i_basis_2) then

             lead_hamiltonian_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) = &
                  lead_hamiltonian_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) + &
                  dconjg(k_phase_tr*lead_1(i_entry)%ham(n_spin_lead))
             
             lead_ovlp_matrix_cmplx(i_basis_2,i_basis_1,i_lead) = &
                  lead_ovlp_matrix_cmplx(i_basis_2,i_basis_1,i_lead) + &
                  dconjg(k_phase_tr*lead_1(i_entry)%ovlp)

          end if

       case(1)

          lead_connection_matrix_h_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) = &
               lead_connection_matrix_h_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) + &
               dconjg(k_phase_tr*lead_1(i_entry)%ham(n_spin_lead))

          lead_connection_matrix_o_cmplx(i_basis_2,i_basis_1,i_lead) = &
               lead_connection_matrix_o_cmplx(i_basis_2,i_basis_1,i_lead) + &
               dconjg(k_phase_tr*lead_1(i_entry)%ovlp)

       case(-1)

          lead_connection_matrix_h_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) = &
               lead_connection_matrix_h_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) + &
               k_phase_tr*lead_1(i_entry)%ham(n_spin_lead)

          lead_connection_matrix_o_cmplx(i_basis_1,i_basis_2,i_lead) = &
               lead_connection_matrix_o_cmplx(i_basis_1,i_basis_2,i_lead) + &
               k_phase_tr*lead_1(i_entry)%ovlp

       end select

    end do

    ! The second lead
    i_lead = 2
    do i_entry = 1, lme2

       i_basis_1 = lead_2(i_entry)%i_basis_1
       i_basis_2 = lead_2(i_entry)%i_basis_2
       i_cell_1 = lead_2(i_entry)%cell_1
       i_cell_2 = lead_2(i_entry)%cell_2

       k_phase_tr = (k_phase_base_tr(1,i_k)**i_cell_1)*(k_phase_base_tr(2,i_k)**i_cell_2)
       
       select case(lead_2(i_entry)%cell_3)

       case(0)

          lead_hamiltonian_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) = &
               lead_hamiltonian_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) + &
               k_phase_tr*lead_2(i_entry)%ham(n_spin_lead)

          lead_ovlp_matrix_cmplx(i_basis_1,i_basis_2,i_lead) = &
               lead_ovlp_matrix_cmplx(i_basis_1,i_basis_2,i_lead) + &
               k_phase_tr*lead_2(i_entry)%ovlp

          if (i_basis_1 /= i_basis_2) then

             lead_hamiltonian_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) = &
                  lead_hamiltonian_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) + &
                  dconjg(k_phase_tr*lead_2(i_entry)%ham(n_spin_lead))
             
             lead_ovlp_matrix_cmplx(i_basis_2,i_basis_1,i_lead) = &
                  lead_ovlp_matrix_cmplx(i_basis_2,i_basis_1,i_lead) + &
                  dconjg(k_phase_tr*lead_2(i_entry)%ovlp)

          end if

       case(1)

          lead_connection_matrix_h_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) = &
               lead_connection_matrix_h_cmplx(i_basis_2,i_basis_1,n_spin_lead,i_lead) + &
               dconjg(k_phase_tr*lead_2(i_entry)%ham(n_spin_lead))

          lead_connection_matrix_o_cmplx(i_basis_2,i_basis_1,i_lead) = &
               lead_connection_matrix_o_cmplx(i_basis_2,i_basis_1,i_lead) + &
               dconjg(k_phase_tr*lead_2(i_entry)%ovlp)

       case(-1)

          lead_connection_matrix_h_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) = &
               lead_connection_matrix_h_cmplx(i_basis_1,i_basis_2,n_spin_lead,i_lead) + &
               k_phase_tr*lead_2(i_entry)%ham(n_spin_lead)

          lead_connection_matrix_o_cmplx(i_basis_1,i_basis_2,i_lead) = &
               lead_connection_matrix_o_cmplx(i_basis_1,i_basis_2,i_lead) + &
               k_phase_tr*lead_2(i_entry)%ovlp

       end select

    end do

  end subroutine construct_lead_matrices

  !****s* transport/evaluate_lead_self_energies
  !  NAME
  !     evaluate_lead_self_energies
  !  SYNOPSIS
  !      subroutine evaluate_lead_self_energies(lead_hamiltonian_1, lead_hamiltonian_2, lead_hamiltonian_3, lead_hamiltonian_4, &
  !       lead_ovlp_matrix_1, lead_ovlp_matrix_2, lead_ovlp_matrix_3, lead_ovlp_matrix_4, &
  !       lead_connection_matrix_h_1, lead_connection_matrix_h_2, lead_connection_matrix_h_3, lead_connection_matrix_h_4, & 
  !       lead_connection_matrix_o_1, lead_connection_matrix_o_2, lead_connection_matrix_o_3, lead_connection_matrix_o_4, &
  !       lead_hamiltonian_cmplx, &
  !       lead_ovlp_matrix_1_cmplx, lead_ovlp_matrix_2_cmplx, &
  !       lead_connection_matrix_h_1_cmplx, lead_connection_matrix_h_2_cmplx, &
  !       lead_connection_matrix_o_1_cmplx, lead_connection_matrix_o_2_cmplx, &
  !       n_lead_basis, ipiv, n_boundary_iterations_max, energy, epsilon_start_max)

  subroutine evaluate_lead_self_energies(lead_hamiltonian, lead_ovlp_matrix, lead_connection_matrix_h, lead_connection_matrix_o, &
       lead_hamiltonian_cmplx, lead_ovlp_matrix_cmplx, lead_connection_matrix_h_cmplx, lead_connection_matrix_o_cmplx, &
       n_lead_basis, ipiv, n_boundary_iterations_max, energy, epsilon_start_max)

    use localorb_io, only: use_unit
    implicit none

    real*8, dimension(n_lead_basis,n_lead_basis,n_spin_lead,n_leads):: lead_hamiltonian
    real*8, dimension(n_lead_basis,n_lead_basis,n_leads):: lead_ovlp_matrix
    real*8, dimension(n_lead_basis,n_lead_basis,n_spin_lead,n_leads):: lead_connection_matrix_h
    real*8, dimension(n_lead_basis,n_lead_basis,n_leads):: lead_connection_matrix_o

    complex*16, dimension(n_lead_basis,n_lead_basis,n_spin_lead,n_leads) :: lead_hamiltonian_cmplx
    complex*16, dimension(n_lead_basis,n_lead_basis,n_leads) :: lead_ovlp_matrix_cmplx
    complex*16, dimension(n_lead_basis,n_lead_basis,n_spin_lead,n_leads):: lead_connection_matrix_h_cmplx
    complex*16, dimension(n_lead_basis,n_lead_basis,n_leads):: lead_connection_matrix_o_cmplx

    complex*16 :: energy, epsilon_start_max
    integer :: n_lead_basis
    integer:: ipiv(n_basis), n_boundary_iterations_max

    complex*16 :: epsilon, energy_lead
    integer :: i_lead, i_basis_1, i_1, apu, info
    real*8 :: new_max
    
    complex*16, allocatable, dimension(:,:):: greenL, workL, ttt, ttt2

    allocate(workL(n_lead_basis, n_lead_basis))
    workL = 0.d0
    allocate(greenL(n_lead_basis, n_lead_basis))
    greenL = 0.d0
    allocate(ttt(n_lead_basis, n_lead_basis))
    ttt = 0.d0
    allocate(ttt2(n_lead_basis, n_lead_basis))
    ttt2 = 0.d0

    apu = 0

    ! Initialize boundary information for the given energy
    do i_lead = 1, n_leads

       energy_lead = energy  + average_pot_lead(i_lead) - average_pot(i_lead)
       epsilon = epsilon_start_max
       
       select case(i_lead)
       case(1)
          if (transport_k_points) then
             workL(1:lead_basis(1),1:lead_basis(1))  = (energy_lead+ epsilon) * lead_ovlp_matrix_cmplx(:,:,i_lead) - &
                  lead_hamiltonian_cmplx(:,:,1,i_lead) - lead_self_energy_1(1:lead_basis(1),1:lead_basis(1))
          else
             workL(1:lead_basis(1),1:lead_basis(1))  = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
                  lead_hamiltonian(:,:,1,i_lead) - lead_self_energy_1(1:lead_basis(1),1:lead_basis(1))
          end if

       case(2)
          if (transport_k_points) then
             workL(1:lead_basis(2),1:lead_basis(2))  = (energy_lead+ epsilon) * lead_ovlp_matrix_cmplx(:,:,i_lead) - &
                  lead_hamiltonian_cmplx(:,:,1,i_lead) - lead_self_energy_2(1:lead_basis(2),1:lead_basis(2))
          else
             workL(1:lead_basis(2),1:lead_basis(2)) = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
                  lead_hamiltonian(:,:,1,i_lead) - lead_self_energy_2(1:lead_basis(2),1:lead_basis(2))
          end if
          
       case(3)
          workL(1:lead_basis(3),1:lead_basis(3)) = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
               lead_hamiltonian(:,:,1,i_lead) - lead_self_energy_3(1:lead_basis(3),1:lead_basis(3))
          
       case(4)
          workL(1:lead_basis(4),1:lead_basis(4)) = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
               lead_hamiltonian(:,:,1,i_lead) - lead_self_energy_4(1:lead_basis(4),1:lead_basis(4))
          
       end select
       
       call ZGETRF(n_lead_basis ,n_lead_basis, workL, n_lead_basis, ipiv, info )
       
       greenL = 0.d0
       do i_basis_1 = 1, n_lead_basis
          greenL(i_basis_1, i_basis_1) = 1.d0
       end do
       
       call  ZGETRS( 'N',n_lead_basis , n_lead_basis, workL, n_lead_basis, ipiv, greenL, n_lead_basis, info )
       
       select case(i_lead)
       case(1)

          if (transport_k_points) then
             ttt(1:lead_basis(1),1:lead_basis(1))  =  & 
                  (real(energy_lead)) *lead_connection_matrix_o_cmplx(1:lead_basis(1),1:lead_basis(1),i_lead) &
                  - lead_connection_matrix_h_cmplx(1:lead_basis(1),1:lead_basis(1),1,i_lead)
          else
             ttt(1:lead_basis(1),1:lead_basis(1))  =  & 
                  (real(energy_lead)) *lead_connection_matrix_o(1:lead_basis(1),1:lead_basis(1),i_lead) &
                  - lead_connection_matrix_h(1:lead_basis(1),1:lead_basis(1),1,i_lead)
          end if
          
       case(2)

          if (transport_k_points) then
             ttt(1:lead_basis(2),1:lead_basis(2))  =  & 
                  (real(energy_lead)) *lead_connection_matrix_o_cmplx(1:lead_basis(2),1:lead_basis(2),i_lead) &
                  - lead_connection_matrix_h_cmplx(1:lead_basis(2),1:lead_basis(2),1,i_lead)
          else          
             ttt(1:lead_basis(2),1:lead_basis(2)) =  & 
                  (real(energy_lead)) *lead_connection_matrix_o(1:lead_basis(2),1:lead_basis(2),i_lead) &
                  - lead_connection_matrix_h(1:lead_basis(2),1:lead_basis(2),1,i_lead)
          end if
          
       case(3)
          
          ttt(1:lead_basis(3),1:lead_basis(3))  =  & 
               (real(energy_lead)) *lead_connection_matrix_o(1:lead_basis(3),1:lead_basis(3),i_lead) &
               - lead_connection_matrix_h(1:lead_basis(3),1:lead_basis(3),1,i_lead)
          
       case(4)
          
          ttt(1:lead_basis(4),1:lead_basis(4))  =  & 
               (real(energy_lead)) * lead_connection_matrix_o(1:lead_basis(4),1:lead_basis(4),i_lead) &
               - lead_connection_matrix_h(1:lead_basis(4),1:lead_basis(4),1,i_lead)
          
       end select

       ! Iterate the Green's function of boundary
       boundary: do i_1 = 1, n_boundary_iterations_max
          
          epsilon = epsilon_end*(0,1) &
               + (epsilon_start_max - epsilon_end*(0,1)) * (n_boundary_iterations_max-i_1) &
               / real(n_boundary_iterations_max-1)**2
          ttt2 = greenL
          
          select case(i_lead)
          case(1)
             
             if (transport_k_points) then
                call ZGEMM('C','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis,&
                     (0.d0,0.d0), workL,n_lead_basis)
             else
                call ZGEMM('T','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis,&
                     (0.d0,0.d0), workL,n_lead_basis)
             end if
             
             call ZGEMM('N','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
                  (0.d0,0.d0), greenL,n_lead_basis)
             
             if(transport_k_points) then
                workL(1:lead_basis(1),1:lead_basis(1)) = (energy_lead+ epsilon) * lead_ovlp_matrix_cmplx(:,:,i_lead) - &
                     lead_hamiltonian_cmplx(:,:,1,i_lead) - greenL(1:lead_basis(1),1:lead_basis(1))
             else
                workL(1:lead_basis(1),1:lead_basis(1)) = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
                     lead_hamiltonian(:,:,1,i_lead) - greenL(1:lead_basis(1),1:lead_basis(1))
             end if
             
          case(2)
             
             if (transport_k_points) then
                call ZGEMM('C','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis,&
                     (0.d0,0.d0), workL,n_lead_basis)
             else
                call ZGEMM('T','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis,&
                     (0.d0,0.d0), workL,n_lead_basis)
             end if
             
             call ZGEMM('N','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
                  (0.d0,0.d0), greenL,n_lead_basis)

             if(transport_k_points) then
                workL(1:lead_basis(2),1:lead_basis(2)) = (energy_lead+ epsilon) * lead_ovlp_matrix_cmplx(:,:,i_lead) - &
                     lead_hamiltonian_cmplx(:,:,1,i_lead) - greenL(1:lead_basis(2),1:lead_basis(2))
             else
                workL(1:lead_basis(2),1:lead_basis(2)) = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
                     lead_hamiltonian(:,:,1,i_lead) - greenL(1:lead_basis(2),1:lead_basis(2))
             end if
             
          case(3)
             
             call ZGEMM('T','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis,&
                  (0.d0,0.d0), workL,n_lead_basis)
             
             call ZGEMM('N','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
                  (0.d0,0.d0), greenL,n_lead_basis)
             
             workL(1:lead_basis(3),1:lead_basis(3)) = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
                  lead_hamiltonian(:,:,1,i_lead) - greenL(1:lead_basis(3),1:lead_basis(3))
             
          case(4)
             
             call ZGEMM('T','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis,&
                  (0.d0,0.d0), workL,n_lead_basis)
             
             call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
                  (0.d0,0.d0), greenL,n_lead_basis)
             
             workL(1:lead_basis(4),1:lead_basis(4)) = (energy_lead+ epsilon) * lead_ovlp_matrix(:,:,i_lead) - &
                  lead_hamiltonian(:,:,1,i_lead) - greenL(1:lead_basis(4),1:lead_basis(4))
             
          end select

          call ZGETRF(lead_basis(i_lead) ,lead_basis(i_lead), workL, &
               n_lead_basis, ipiv, info )
          
          greenL = 0.d0
          do i_basis_1 = 1, lead_basis(i_lead)
             greenL(i_basis_1, i_basis_1) = 1.d0
          end do
          
          call ZGETRS( 'N',lead_basis(i_lead) , lead_basis(i_lead), &
               workL, n_lead_basis, ipiv, greenL, n_lead_basis, info )
          
          apu = apu+1             
          
          new_max = maxval(abs(greenL(1:lead_basis(i_lead),1:lead_basis(i_lead)) &
               - ttt2(1:lead_basis(i_lead),1:lead_basis(i_lead))))
          
          ! Is the boundary condition converged?
          if(abs(new_max)< boundary_treshold)then
             exit boundary
          end if
          
          ! Mixing to previous solution
          greenL = greenL*boundary_mix + ttt2* ( 1.d0 - boundary_mix)
          
       end do boundary
       
       if(i_1 >= n_boundary_iterations_max)then
          write(use_unit,*) 'Warning: ', new_max, ' boundary iterations did not converge.'
       end if
       n_boundary_iterations_max = n_boundary_iterations
       epsilon_start_max = epsilon_start
       
       ! Final steps of Green's functions of the boundary

       select case(i_lead)
       case(1)

          if (transport_k_points) then
             call ZGEMM('C','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis, &
                  (0.d0,0.d0), workL,n_lead_basis)
          else
             call ZGEMM('T','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis, &
                     (0.d0,0.d0), workL,n_lead_basis)
          end if

          call ZGEMM('N','N',lead_basis(1),lead_basis(1),lead_basis(1),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
               (0.d0,0.d0), lead_self_energy_1,lead_basis(1))
          
       case(2)

          if (transport_k_points) then
             call ZGEMM('C','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis, &
                  (0.d0,0.d0), workL,n_lead_basis)
          else
             call ZGEMM('T','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis, &
                  (0.d0,0.d0), workL,n_lead_basis)
          end if

          call ZGEMM('N','N',lead_basis(2),lead_basis(2),lead_basis(2),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
                     (0.d0,0.d0), lead_self_energy_2,lead_basis(2))
          
       case(3)
          
          call ZGEMM('T','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis, &
               (0.d0,0.d0), workL,n_lead_basis)
          
          call ZGEMM('N','N',lead_basis(3),lead_basis(3),lead_basis(3),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
               (0.d0,0.d0), lead_self_energy_3,lead_basis(3))
          
       case(4)
          
          call ZGEMM('T','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),ttt, n_lead_basis, greenL,n_lead_basis, &
               (0.d0,0.d0), workL,n_lead_basis)
          
          call ZGEMM('N','N',lead_basis(4),lead_basis(4),lead_basis(4),(1.d0,0.d0),workL, n_lead_basis, ttt,n_lead_basis, &
               (0.d0,0.d0), lead_self_energy_4,lead_basis(4))
          
       end select
       
    end do ! i_lead

  end subroutine evaluate_lead_self_energies

End module transport
!******
