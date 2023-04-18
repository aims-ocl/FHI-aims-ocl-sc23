!****s* FHI-aims/update_full_matrix_p0
!  NAME
!   update_full_matrix_p0
!  SYNOPSIS

subroutine update_full_matrix_p0 &
     ( n_compute_c, n_compute_a, i_basis, &
     matrix_shell, matrix &
     )

!  PURPOSE
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
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop
  use pbc_lists
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer n_compute_c
  integer n_compute_a
  integer i_basis(n_compute_c)
  real*8 matrix_shell(n_compute_c,n_compute_a)
  real*8 matrix( n_hamiltonian_matrix_size)

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

!  WPH:  Details about how the packing schemes work (which are necessary to
!        understand this subroutine) may be found in the heading of physics.f90,
!        under the sub-heading "hamiltonian".






  !  local variables

  !     auxiliary matrices for Level 3 Blas matrix multiplications

  !     counters

  integer :: i_compute_1
  integer :: i_compute_2
  integer :: i_offset
  integer :: i_index_real
  integer :: i_offset_first_part
  integer :: i_cell_index, i_cell_1
  integer :: i_max_basis, i_min_basis
  integer :: i_one_part, i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: offset(n_cells) !_in_hamiltonian)
  integer :: offset_end(n_cells) !_in_hamiltonian)
  integer :: help
!  integer :: direct_citing(n_compute_a,n_compute_a )
!  real*8 :: data(n_cells_in_hamiltonian, n_basis)
!  real*8,dimension(:,:), allocatable :: data
  integer:: i_cell_old



  !     begin work

  !      now add the aux. ham. matrix to the actual ham. matrix
  !      this requires translating between the actually computed matrix elements and the
  !      full ham. matrix ...

  ! When basis is smaller than n_basis

  !      write(use_unit,*) n_compute_c,n_compute_a

  select case(packed_matrix_format)

  case(PM_none)
     i_index_real = 0
     do i_compute_2 = 1, n_compute_a, 1

        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2

        do i_compute_1 = 1,i_compute_2,1
           i_index_real = i_offset + i_basis(i_compute_1) 


           matrix(i_index_real) = matrix(i_index_real) &
                + matrix_shell(i_compute_1, i_compute_2)

        enddo
     enddo


  case(PM_index) !------------------------------------------------------------------------------------


     if(n_periodic == 0)then

        do i_compute_1 = 1, n_compute_a, 1

           
!           write(use_unit,*) 'i_compute_2', i_compute_2
           
           i_start =  index_hamiltonian(1,1, i_basis(i_compute_1))
           i_end   =  index_hamiltonian(2,1, i_basis(i_compute_1))

           if(i_end<i_start) cycle ! otherways an undefined i_index_real is used!

!           write(use_unit,*) '-'

           do i_compute_2 = 1,i_compute_1,1
              
              i_basis_2 = i_basis(i_compute_2)
              
              place: do i_place = i_start, i_end, 1
                 
                 if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
                    matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
                    i_index_real = i_place
                    exit place 

                 
                 else if(column_index_hamiltonian( i_place) > i_basis_2)then
                    i_index_real = i_place
                    exit place  

                 end if
              end do place
              i_start = i_index_real
           end do
        end do
        


           
        
     else ! Periodic case----------------



        ! Unfortunately the periodic systems can not use the searching routine used now in the clusters.
        ! This is because the peridic systems have extra packing for supercell information.


        do i_compute_1 = 1, n_compute_a, 1
           do i_compute_2 = 1,i_compute_1-1,1

             
              matrix_shell(i_compute_1, i_compute_2) = &
                 matrix_shell(i_compute_2, i_compute_1) 

           end do
        end do


        do i_compute_1 = 1, n_compute_a, 1

           i_basis_1 = Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_old = center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))

           offset_end = -1
           offset = -1

           do i_cell_1 = 1, n_cells 


              i_cell = position_in_hamiltonian( i_cell_old, i_cell_1) 

              offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
              offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)

           end do



           do i_compute_2 = 1, n_compute_a
              
              i_basis_2 = Cbasis_to_basis(i_basis(i_compute_2))


              if(i_basis_2 <= i_basis_1)then

                 i_cell    =  center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))

                 place_2: do i_place = offset(i_cell), offset_end(i_cell),1 

!                    write(use_unit,*) offset(i_cell), offset_end(i_cell)
                 
                    if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
                       matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
                       exit place_2
                 
                    else if(column_index_hamiltonian( i_place) > i_basis_2)then
                       exit place_2  


                    end if
                 end do place_2

                 offset(i_cell) = i_place

              end if
           end do
        end do
        


         ! Here is the old implementation, just in case.
        

!        do i_compute_2 = 1, n_compute_a, 1
!
!           do i_compute_1 = 1,n_compute_a,1
!
!              !               if(Cbasis_to_basis(i_basis(i_compute_1)) <= Cbasis_to_basis(i_basis(i_compute_2)))then
!
!
!              i_index_real =  find_position_in_hamiltonian( i_basis(i_compute_2), i_basis(i_compute_1)) 
!
!              !                  if( i_index_real == n_trash_start)then
!              !                     if(matrix_shell(i_compute_1, i_compute_2) /= 0)then
!              !                        write(use_unit,*) 'Trash:', matrix_shell(i_compute_1, i_compute_2), i_compute_1, i_compute_2
!              !                     end if
!              !                  end if
!
!
!              !                  write(use_unit,*) i_basis(i_compute_1), i_basis(i_compute_2), i_index_real
!
!              if(i_index_real/=n_trash_start)  help = help + 1
!
!              
!
!
!              if(i_compute_1 <= i_compute_2 )then
!
!                 matrix(i_index_real) = matrix(i_index_real) &
!                      + matrix_shell(i_compute_1, i_compute_2)
!
!              else
!
!                 matrix(i_index_real) = matrix(i_index_real) &
!                      + matrix_shell(i_compute_2, i_compute_1)
!
!              end if
!
!
!              !              end if
!           end do
!        end do
!


     end if ! n_periodic == 0 - else


  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select







end subroutine update_full_matrix_p0
!---------------------------------------------------------------------
!******
subroutine update_full_matrix_p0X &
     ( n_compute_c, n_compute_a, i_basis, &
     matrix_shell, matrix &
     )

!  PURPOSE
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
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop
  use pbc_lists
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer n_compute_c
  integer n_compute_a
  integer i_basis(n_compute_c)
  real*8 matrix_shell(n_compute_c,n_compute_a)
  real*8 matrix( n_hamiltonian_matrix_size)

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

  integer :: i_compute_1
  integer :: i_compute_2
  integer :: i_offset
  integer :: i_index_real
  integer :: i_offset_first_part
  integer :: i_cell_index, i_cell_1
  integer :: i_max_basis, i_min_basis
  integer :: i_one_part, i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: offset(n_cells) !_in_hamiltonian)
  integer :: offset_end(n_cells) !_in_hamiltonian)
  integer :: help
!  integer :: direct_citing(n_compute_a,n_compute_a )
!  real*8 :: data(n_cells_in_hamiltonian, n_basis)
!  real*8,dimension(:,:), allocatable :: data
  integer:: i_cell_old

!NEC_CB
  integer::help1(n_compute_a)
  integer::help2(n_compute_a)

  !     begin work

  !      now add the aux. ham. matrix to the actual ham. matrix
  !      this requires translating between the actually 
  !      computed matrix elements and the full ham. matrix ...

  !      When basis is smaller than n_basis

  !      write(use_unit,*) n_compute_c,n_compute_a


  select case(packed_matrix_format)

  case(PM_none)


     i_index_real = 0
     do i_compute_2 = 1, n_compute_a, 1

        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2


        do i_compute_1 = 1,i_compute_2,1


           i_index_real = i_offset + i_basis(i_compute_1)

           matrix(i_index_real) = matrix(i_index_real) &
                + matrix_shell(i_compute_1, i_compute_2)


        enddo
     enddo


  case(PM_index) !------------------------------------------------------------------------------------


     if(n_periodic == 0)then

        do i_compute_1 = 1, n_compute_a, 1

           
!           write(use_unit,*) 'i_compute_2', i_compute_2
           
           i_start =  index_hamiltonian(1,1, i_basis(i_compute_1))
           i_end   =  index_hamiltonian(2,1, i_basis(i_compute_1))

!           write(use_unit,*) '-'

           do i_compute_2 = 1,i_compute_1,1
              
              i_basis_2 = i_basis(i_compute_2)
              
              place: do i_place = i_start, i_end, 1
                 
                 if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
!NEC_CB                    matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
                    if (i_compute_2.le.i_compute_1) then
                       matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
                    else
                       matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_1, i_compute_2)
                    end if

                    i_index_real = i_place
                    exit place 

                 
                 else if(column_index_hamiltonian( i_place) > i_basis_2)then
                    i_index_real = i_place
                    exit place  

                 end if
              end do place
              i_start = i_index_real
           end do
        end do
        


           
        
     else ! Periodic case----------------



        ! Unfortunately the periodic systems can not use the searching routine used now in the clusters.
        ! This is because the peridic systems have extra packing for supercell information.

!NEC_CB matrix_shell(i,j) is used only with if-clause below for i<=j
!NEC_CB        do i_compute_1 = 1, n_compute_a, 1
!NEC_CB           do i_compute_2 = 1,i_compute_1-1,1
!NEC_CB
!NEC_CB             
!NEC_CB              matrix_shell(i_compute_1, i_compute_2) = matrix_shell(i_compute_2, i_compute_1) 
!NEC_CB
!NEC_CB           end do
!NEC_CB        end do

!NEC_CB Do not read multi-indirectly addressed indices multiple times
        do i_compute_1 = 1, n_compute_a, 1
          help1(i_compute_1)=Cbasis_to_basis(i_basis(i_compute_1))
          help2(i_compute_1)=center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
        end do

        do i_compute_1 = 1, n_compute_a, 1

           i_basis_1 = help1(i_compute_1)!Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_old = help2(i_compute_1)!center_to_cell(Cbasis_to_center(i_basis(i_compute_1))

           offset_end = -1
           offset = -1

           do i_cell_1 = 1, n_cells 


             i_cell = position_in_hamiltonian( i_cell_old, i_cell_1) 

              offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
              offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)

           end do



           do i_compute_2 = 1, n_compute_a
              
              i_basis_2 = help1(i_compute_2)!Cbasis_to_basis(i_basis(i_compute_2))



              if(i_basis_2 <= i_basis_1)then

                 i_cell    =  help2(i_compute_2)!center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))


                 place_2: do i_place = offset(i_cell), offset_end(i_cell),1 

!                    write(use_unit,*) offset(i_cell), offset_end(i_cell)
                 
                    if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
!NEC_CB                       matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
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

                 offset(i_cell) = i_place

              end if
           end do
        end do



         ! Here is the old implementation, just in case.
        

!        do i_compute_2 = 1, n_compute_a, 1
!
!           do i_compute_1 = 1,n_compute_a,1

!              !               if(Cbasis_to_basis(i_basis(i_compute_1)) <= Cbasis_to_basis(i_basis(i_compute_2)))then
!              print *, "Basis 1: ", i_compute_1, " iBasis 1: ", i_basis(i_compute_1)
!              print *, "Basis 2: ", i_compute_2, " iBasis 2: ", i_basis(i_compute_2)

!              i_index_real =  find_position_in_hamiltonian( &
!                    i_basis(i_compute_2), &
!                    i_basis(i_compute_1))

!                 print *, "CPU", &
!                          " Target for basis2 ",i_compute_2, &
!                          " basis1 ",i_compute_1, &
!                          " is ", i_index_real

              !                  if( i_index_real == n_trash_start)then
              !                     if(matrix_shell(i_compute_1, i_compute_2) /= 0)then
              !                        write(use_unit,*) 'Trash:', matrix_shell(i_compute_1, i_compute_2), i_compute_1, i_compute_2
              !                     end if
              !                  end if


              !                  write(use_unit,*) i_basis(i_compute_1), i_basis(i_compute_2), i_index_real

!              if(i_index_real/=n_trash_start)  help = help + 1

              


!              if(i_compute_1 <= i_compute_2 )then
!                 print *, " H = ",matrix(i_index_real), &
!                          " H_Shell = ",matrix_shell(i_compute_1, i_compute_2)
!
!                 matrix(i_index_real) = matrix(i_index_real) &
!                      + matrix_shell(i_compute_1, i_compute_2)
!
!              else
!                 print *, " H = ",matrix(i_index_real), &
!                          " H_Shell = ",matrix_shell(i_compute_2, i_compute_1)
!
!                 matrix(i_index_real) = matrix(i_index_real) &
!                      + matrix_shell(i_compute_2, i_compute_1)
!
!              end if


              !              end if
!           end do
!        end do



     end if ! n_periodic == 0 - else


  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select







end subroutine update_full_matrix_p0X
!---------------------------------------------------------------------
!******

!******
subroutine get_update_full_matrix_map ( i_spin, ld_matrix, n_compute_c, i_basis, map )

!  PURPOSE
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
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop
  implicit none

!  ARGUMENTS

  integer, intent(IN) :: i_spin
  integer, intent(IN) :: ld_matrix
  integer, intent(IN) :: n_compute_c
  integer, intent(IN),  dimension(1:n_compute_c) :: i_basis
  integer, intent(OUT), dimension (1:n_compute_c*n_compute_c) :: map

!  INPUTS
!   o i_spin       -- The spin index.  For a matrix without spin indexing, set to 1.
!   o ld_matrix    -- The leading dimension of the matrix (only needed for spin-indexed matrices)
!   o n_compute_c  -- number of relevant basis functions
!   o i_basis      -- list of relevant basis functions
!
!  OUTPUT
!   o map          -- mapping from matrix shell to global matrix
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

  !     counters

  integer :: i_compute_1
  integer :: i_compute_2
  integer :: element
  integer :: i_offset
  integer :: i_index_real
  integer :: i_index_shell
  integer :: i_offset_first_part
  integer :: i_cell_index, i_cell_1
  integer :: i_max_basis, i_min_basis
  integer :: i_one_part, i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: offset(n_cells) !_in_hamiltonian)
  integer :: offset_end(n_cells) !_in_hamiltonian)
  integer :: help
  integer:: i_cell_old

  integer::help1(n_compute_c)
  integer::help2(n_compute_c)

  !     begin work

  map = 0

  select case(packed_matrix_format)

  case(PM_none)

     i_index_real = 0
     do i_compute_2 = 1, n_compute_c, 1

        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2

        do i_compute_1 = 1,i_compute_2,1

           i_index_shell = i_compute_1 + (i_compute_2-1) * n_compute_c
           i_index_real = i_offset + i_basis(i_compute_1)

           map(i_index_shell) = i_index_real + (i_spin-1)*ld_matrix

        enddo
     enddo


  case(PM_index) 

     if(n_periodic == 0)then

        do i_compute_1 = 1, n_compute_c, 1
           
           i_start =  index_hamiltonian(1,1, i_basis(i_compute_1))
           i_end   =  index_hamiltonian(2,1, i_basis(i_compute_1))

           do i_compute_2 = 1, i_compute_1, 1
              
              i_basis_2 = i_basis(i_compute_2)
              
              place: do i_place = i_start, i_end, 1
                 
                 if( column_index_hamiltonian( i_place) == i_basis_2)then

                    i_index_shell = i_compute_2 + (i_compute_1-1) * n_compute_c

                    map(i_index_shell) = i_place + (i_spin-1)*ld_matrix 

                    i_index_real = i_place
                    exit place 
                 
                 else if(column_index_hamiltonian( i_place) > i_basis_2)then
                    i_index_real = i_place
                    exit place  

                 end if
              end do place
              i_start = i_index_real
           end do
        end do
        
     else ! Periodic case----------------

        ! Unfortunately the periodic systems can not use the searching routine
        ! used now in the clusters.
        ! This is because the peridic systems have extra packing for 
        ! supercell information.

        do i_compute_1 = 1, n_compute_c, 1
          help1(i_compute_1) = Cbasis_to_basis(i_basis(i_compute_1))
          help2(i_compute_1) = center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
        end do

        do i_compute_1 = 1, n_compute_c, 1

           i_basis_1 = help1(i_compute_1)!Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_old = help2(i_compute_1)!center_to_cell(Cbasis_to_center(i_basis(i_compute_1))

           offset_end = -1
           offset = -1

           do i_cell_1 = 1, n_cells 

             i_cell = position_in_hamiltonian( i_cell_old, i_cell_1) 

             offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
             offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)

           end do

           do i_compute_2 = 1, n_compute_c, 1
              
              i_basis_2 = help1(i_compute_2)!Cbasis_to_basis(i_basis(i_compute_2))

              if(i_basis_2 <= i_basis_1)then

                 i_cell    =  help2(i_compute_2)!center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))

                 place_2: do i_place = offset(i_cell), offset_end(i_cell),1 

                    if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
                       if (i_compute_2.le.i_compute_1) then
                         i_index_shell = i_compute_2 + (i_compute_1-1) * n_compute_c
                       else
                         i_index_shell = i_compute_1 + (i_compute_2-1) * n_compute_c
                       end if
                    
                       map(i_index_shell) = i_place + (i_spin-1)*ld_matrix 

                       exit place_2
                 
                    elseif (column_index_hamiltonian( i_place) > i_basis_2) then
                       exit place_2  
                    endif
                 enddo place_2

                 offset(i_cell) = i_place

              end if
           end do
        end do
     end if ! n_periodic == 0 - else


  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select

end subroutine get_update_full_matrix_map
!---------------------------------------------------------------------
!******
