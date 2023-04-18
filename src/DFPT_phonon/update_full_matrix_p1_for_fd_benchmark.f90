!****s* FHI-aims/update_full_matrix_p1_for_fd_benchmark
!  NAME
!   update_full_matrix_p1_for_fd_benchmark
!  SYNOPSIS

subroutine update_full_matrix_p1_for_fd_benchmark &
     ( n_compute_c, n_compute_a, i_basis, &
     matrix_shell, matrix &
     )

!  PURPOSE

!  shanghui test version, which need complex matrix

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
  use localorb_io, only: localorb_info, use_unit
  use mpi_tasks, only: aims_stop
  implicit none

!  ARGUMENTS

  integer n_compute_c
  integer n_compute_a
  integer i_basis(n_compute_c)
  real*8 matrix_shell(n_compute_c,n_compute_a)
  real*8 matrix(n_hamiltonian_matrix_size)

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


    !do i_coord=1,3

        do i_compute_1 = 1, n_compute_a

           i_basis_1 = help1(i_compute_1)!Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_1  = help2(i_compute_1)!center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))

           do i_compute_2 = 1, n_compute_a
              
              i_basis_2 = help1(i_compute_2)!Cbasis_to_basis(i_basis(i_compute_2))

              if(i_basis_2 <= i_basis_1) then

                 i_cell_2    =  help2(i_compute_2)!center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))
 
                 !------we only count unit-cell here------------------
                 if(i_cell_1.eq.1.and.i_cell_2.eq.1) then 

                  !if(i_basis(i_compute_1).gt.26.or.i_basis(i_compute_2).gt.26) then 
                  !write(use_unit,*) 'i_cell_1,i_cell_2:',i_cell_1,i_cell_2
                  !endif

                 place_2: do i_place = &  
                 index_hamiltonian(1,position_in_hamiltonian(i_cell_1,i_cell_2), i_basis_1),  & 
                 index_hamiltonian(2,position_in_hamiltonian(i_cell_1,i_cell_2), i_basis_1)  
                 

                    if( column_index_hamiltonian( i_place) == i_basis_2)then

                       if (i_compute_2.le.i_compute_1) then

                         matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
                       else

                         matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_1, i_compute_2)
                       end if

                       exit place_2
                 
                    else if(column_index_hamiltonian( i_place) > i_basis_2)then
                       exit place_2  

                    end if
                 end do place_2
             
                 endif ! i_cell_1=i_cell_2=1

              end if
           end do
        end do

    !enddo



  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select





end subroutine update_full_matrix_p1_for_fd_benchmark
!---------------------------------------------------------------------
!******
