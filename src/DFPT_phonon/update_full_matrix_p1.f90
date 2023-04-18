!****s* FHI-aims/update_full_matrix_p1
!  NAME
!   update_full_matrix_p1
!  SYNOPSIS

subroutine update_full_matrix_p1 &
     ( n_compute_c, n_compute_a, i_basis, &
     matrix_shell, matrix &
     )

!  PURPOSE

!  shanghui extened to first-order-X version, which need complex matrix

!  Subroutine update_full_matrix adds a part of the integrals in a
!  matrix (overlap or Hamiltonian matrix)
!  (only for the n_compute basis functions that are nonzero at the
!  current integration shell) to the full matrix (which contains
!  all n_basis basis functions). The link between i_compute = 1 ... n_compute
!  and i_basis = 1 ... n_basis is provided by the index array 
!  i_basis(i_compute).
!
!  USES

  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: aims_stop
  use localorb_io, only: localorb_info, use_unit
  implicit none

!  ARGUMENTS

  integer n_compute_c
  integer n_compute_a
  integer i_basis(n_compute_c)
  real*8 matrix_shell(3,n_compute_c,n_compute_c,n_compute_a)
  real*8 matrix( 3,n_centers_in_hamiltonian,n_hamiltonian_matrix_size)

!  INPUTS
!   o n_compute_c == n_compute_a -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
!   o matrix_shell -- hamiltonian / overlap_matrix of relevant basis functions
!
!  OUTPUT
!   o matrix -- date from matrix_shell is added here
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

  !     auxiliary matrices for Level 3 Blas matrix multiplications

  !     counters

  integer :: i_compute_1,i_compute_2
  integer :: i_cell_1,i_cell_2
  integer :: i_place, i_basis_2, i_basis_1, i_cell
!  integer :: direct_citing(n_compute_a,n_compute_a )
!  real*8 :: data(n_cells_in_hamiltonian, n_basis)
!  real*8,dimension(:,:), allocatable :: data
  integer:: i_cell_old
  integer :: i_coord


!NEC_CB
  integer::help1(n_compute_a)
  integer::help2(n_compute_a)

  !     begin work

  !      now add the aux. ham. matrix to the actual ham. matrix
  !      this requires translating between the actually computed matrix elements and the
  !      full ham. matrix ...

  ! When basis is smaller than n_basis

!-----------------shanghui debug if inv_centers_in_hamiltonian right---------------
!  write(use_unit,*) 'i_center, real_center, i_center'
! do i_compute_1=1, n_centers_in_hamiltonian
!    write(use_unit,*) i_compute_1, centers_in_hamiltonian(i_compute_1),  & 
!    inv_centers_in_hamiltonian(centers_in_hamiltonian(i_compute_1))
! enddo
! stop
!-----------------shanghui end debug if inv_centers_in_hamiltonian right---------------


  select case(packed_matrix_format) 

  case(PM_none)
      write(use_unit,*) 'stop! we need always use PM_index'
      stop

  case(PM_index) !------------------------------------------------------------------------------------

        do i_compute_1 = 1, n_compute_a, 1
          help1(i_compute_1)=Cbasis_to_basis(i_basis(i_compute_1))
          help2(i_compute_1)=center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
        end do


    do i_coord=1,3

        do i_compute_1 = 1, n_compute_a

           i_basis_1 = help1(i_compute_1)!Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_1  = help2(i_compute_1)!center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))

           do i_compute_2 = 1, n_compute_a
              
              i_basis_2 = help1(i_compute_2)!Cbasis_to_basis(i_basis(i_compute_2))

              if(i_basis_2 <= i_basis_1) then

                 i_cell_2    =  help2(i_compute_2)!center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))


                 place_2: do i_place = &  
                 index_hamiltonian(1,position_in_hamiltonian(i_cell_1,i_cell_2), i_basis_1),  & 
                 index_hamiltonian(2,position_in_hamiltonian(i_cell_1,i_cell_2), i_basis_1)  
                 

                    if( column_index_hamiltonian( i_place) == i_basis_2)then


                matrix( i_coord,  &
                        inv_centers_in_hamiltonian(  & 
                        cell_and_atom_to_center(position_in_hamiltonian(i_cell_1,i_cell_2), &
                                                Cbasis_to_atom(i_basis(i_compute_1)))),     & 
                        i_place) = & 
                matrix( i_coord,  &
                        inv_centers_in_hamiltonian(  & 
                        cell_and_atom_to_center(position_in_hamiltonian(i_cell_1,i_cell_2), &
                                                Cbasis_to_atom(i_basis(i_compute_1)))),     & 
                        i_place) + & 
                matrix_shell(i_coord,i_compute_1 ,i_compute_1, i_compute_2)



                matrix( i_coord,  &
                        inv_centers_in_hamiltonian(Cbasis_to_atom(i_basis(i_compute_2))),     & 
                        i_place) = & 
                matrix( i_coord,  &
                        inv_centers_in_hamiltonian(Cbasis_to_atom(i_basis(i_compute_2))),     & 
                        i_place) + & 
                matrix_shell(i_coord,i_compute_2 ,i_compute_1, i_compute_2)

                       exit place_2
                 
                    else if(column_index_hamiltonian( i_place) > i_basis_2)then
                       exit place_2  

                    end if
                 end do place_2

              end if
           end do
        end do

    enddo


!-------------begin no_symmetry version----------------------
!    do i_coord=1,3
! 
!        do i_compute_1 = 1, n_compute_a
! 
!           i_basis_1 = help1(i_compute_1)!Cbasis_to_basis(i_basis(i_compute_1))
!           i_cell_1  = help2(i_compute_1)!center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
! 
!           do i_compute_2 = 1, n_compute_a
!              
!              i_basis_2 = help1(i_compute_2)!Cbasis_to_basis(i_basis(i_compute_2))
! 
!              !if(i_basis_2 <= i_basis_1) then
! 
!                 i_cell_2    =  help2(i_compute_2)!center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))
! 
! 
!                 place_2: do i_place = &  
!                 index_hamiltonian_no_symmetry(1,position_in_hamiltonian(i_cell_1,i_cell_2), i_basis_1),  & 
!                 index_hamiltonian_no_symmetry(2,position_in_hamiltonian(i_cell_1,i_cell_2), i_basis_1)  
!                 
! 
!                    if( column_index_hamiltonian_no_symmetry( i_place) == i_basis_2)then
! 
! 
!                matrix( i_coord,  &
!                        inv_centers_in_hamiltonian(  & 
!                        cell_and_atom_to_center(position_in_hamiltonian(i_cell_1,i_cell_2), &
!                                                Cbasis_to_atom(i_basis(i_compute_1)))),     & 
!                        i_place) = & 
!                matrix( i_coord,  &
!                        inv_centers_in_hamiltonian(  & 
!                        cell_and_atom_to_center(position_in_hamiltonian(i_cell_1,i_cell_2), &
!                                                Cbasis_to_atom(i_basis(i_compute_1)))),     & 
!                        i_place) + & 
!                matrix_shell(i_coord,i_compute_1 ,i_compute_1, i_compute_2)
! 
! 
! 
!                matrix( i_coord,  &
!                        inv_centers_in_hamiltonian(Cbasis_to_atom(i_basis(i_compute_2))),     & 
!                        i_place) = & 
!                matrix( i_coord,  &
!                        inv_centers_in_hamiltonian(Cbasis_to_atom(i_basis(i_compute_2))),     & 
!                        i_place) + & 
!                matrix_shell(i_coord,i_compute_2 ,i_compute_1, i_compute_2)
! 
!                       exit place_2
!                 
!                    else if(column_index_hamiltonian_no_symmetry( i_place) > i_basis_2)then
!                       exit place_2  
! 
!                    end if
!                 end do place_2
! 
!             ! end if
!           end do
!        end do
! 
!    enddo
!-------------end no_symmetry version----------------------


  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select





end subroutine update_full_matrix_p1
!---------------------------------------------------------------------
!******
