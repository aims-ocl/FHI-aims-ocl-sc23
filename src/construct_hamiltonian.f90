!****s* FHI-aims/construct_hamiltonian
!  NAME
!   construct_hamiltonian
!  SYNOPSIS

subroutine construct_hamiltonian(hamiltonian, hamiltonian_w, &
     hamiltonian_w_complex, k_point,work_ham )


  !  PURPOSE
  !   Contruct hamiltonian matrix for the lapack eigenvalue solver.
  !   This is needed in periodic systems and if the packed matrix format is in use.
  !
  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  use mpi_tasks, only: n_tasks, aims_stop
  use localorb_io, only: localorb_info
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  implicit none

  !  ARGUMENTS

  real*8::     hamiltonian ( n_hamiltonian_matrix_size,n_spin )
  real*8::     hamiltonian_w         ( n_basis*(n_basis+1)/2,n_spin)
  complex*16:: hamiltonian_w_complex ( n_basis*(n_basis+1)/2,n_spin)
  real*8, dimension( n_centers_basis_I , n_centers_basis_I,n_spin ),optional :: work_ham
  integer:: k_point

  ! JW: In principle, optional dummy arguments are only allowed with an
  ! explicit interface (e.g. within a module subroutine).  While I do not
  ! really expect trouble, this might be a possible explanation if there is
  ! trouble with some compilers.

  !  INPUTS
  !    o hamiltonian -- hamiltonian matrix with centers terms separately
  !    o k_point -- k-point wanted to be calculated 
  !    o work_ham -- work space needed if for non-packed format is in use
  !  OUTPUT
  !    o hamiltonian_w -- hamiltonian matrix ready to go lapack eigenvalue solver if real eigenvector are used
  !    o hamiltonian_w_complex  -- hamiltonian matrix ready to go lapack eigenvalue solver if complex eigenvector are used
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE




  integer,dimension(:,:),allocatable :: index_b
  integer,dimension(:,:),allocatable :: index_b_old
  integer,dimension(:,:),allocatable :: index_b_new
  complex*16, dimension(:,:),allocatable :: work


  integer:: n_first_part, i_cell, i_basis, cell_shift, i_index_real, i_size
  integer:: i_index_new, i_index_old, i_index, i_basis_2, i_basis_1

  integer:: i_spin



  select case( packed_matrix_format)

  case(PM_index)


     if(real_eigenvectors)then
        hamiltonian_w = 0.0
     else
        hamiltonian_w_complex = 0.0
     end if


     n_first_part =  n_basis*(n_basis+1)/2

     call aims_allocate(index_b, n_basis,n_basis, "index_b" )


     index_b = 0
     i_index = 0
     do i_basis_2 = 1,n_basis, 1
        do i_basis_1 = 1,i_basis_2,1
           i_index = i_index + 1
           index_b(i_basis_1, i_basis_2) = i_index 
        end do
     end do


     i_index_real = 0



     do i_cell = 1,n_cells_in_hamiltonian-1

        do i_basis_2 = 1, n_basis

           if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

              i_index_real = index_hamiltonian(1,i_cell, i_basis_2)-1

              do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                 i_index_real = i_index_real + 1
                 i_basis_1 =  column_index_hamiltonian(i_index_real)


                 if(real_eigenvectors)then

                    do i_spin = 1, n_spin

                       hamiltonian_w(index_b(i_basis_1,i_basis_2) ,i_spin) =  &
                            hamiltonian_w(index_b(i_basis_1,i_basis_2),i_spin) &
                            + dble(k_phase( i_cell,k_point))  &
                            * hamiltonian(  i_index_real , i_spin)
                    end do


                 else

                    do i_spin = 1, n_spin

                       hamiltonian_w_complex(index_b(i_basis_1,i_basis_2) ,i_spin) = &
                            hamiltonian_w_complex(index_b(i_basis_1,i_basis_2),i_spin) &
                            + k_phase( i_cell,k_point)  &
                            * hamiltonian(  i_index_real , i_spin)

                    end do
                 end if
              end do
           end if
        end do
     end do



     call aims_deallocate( index_b, "index_b" )


  case(PM_none) !-----------------------------------------------------------------------------



     call aims_allocate( index_b_old, n_centers_basis_I, n_centers_basis_I, "index_b_old" )
     call aims_allocate( index_b_new, n_basis, n_basis,                     "index_b_new" )


     index_b_old = 0
     index_b_new = 0
     i_index = 0

     do i_basis_2 = 1,n_basis, 1
        do i_basis_1 = 1,i_basis_2,1

           i_index = i_index + 1
           index_b_new(i_basis_1, i_basis_2) = i_index 

        end do
     end do

     i_index = 0
     do i_basis_2 = 1, n_centers_basis_I, 1
        do i_basis_1 = 1,i_basis_2,1

           i_index = i_index + 1
           index_b_old(i_basis_1, i_basis_2) = i_index 

        end do
     end do



     if(real_eigenvectors)then
        hamiltonian_w = 0.0
     else
        hamiltonian_w_complex = 0.0
     end if

     if(real_eigenvectors)then


        do i_basis_2 = 1,n_centers_basis_I,1 
           do i_basis_1 = 1,n_centers_basis_I,1 


              i_index_old =  max( index_b_old( i_basis_2, i_basis_1), index_b_old( i_basis_1, i_basis_2))

              i_index_new =  index_b_new( Cbasis_to_basis(i_basis_1),  Cbasis_to_basis(i_basis_2)) 

              if(i_index_new /= 0 .and.  i_index_old /= 0)then

                 if(real_eigenvectors)then

                    do i_spin = 1, n_spin


                       hamiltonian_w(i_index_new,i_spin) =    hamiltonian_w(i_index_new,i_spin)   &
                            +  dble(k_phase(center_to_cell(Cbasis_to_center( i_basis_2 )),k_point))  &
                            * dble(k_phase(center_to_cell(Cbasis_to_center( i_basis_1 )),k_point))  &
                            * hamiltonian( i_index_old,i_spin)

                    end do
                 end if
              end if
           end do
        end do





     else ! complex eigenvalues ------------------------------------------------------------


        call aims_allocate( work, n_centers_basis_I, n_basis, "work" )


        if(k_point <= n_tasks)then


           work_ham = 0.d0

           do i_spin = 1, n_spin

              do i_basis_1 = 1,n_centers_basis_I,1 

                 work_ham(1:i_basis_1, i_basis_1,i_spin) =  hamiltonian( index_b_old( 1:i_basis_1, i_basis_1),i_spin)

              end do

              work_ham(:,:,i_spin) = work_ham(:,:, i_spin) + transpose(work_ham(:,:,i_spin))

              do i_basis_2 = 1,n_centers_basis_I,1 
                 work_ham(i_basis_2, i_basis_2,i_spin) = work_ham(i_basis_2, i_basis_2,i_spin)/2
              end do
           end do
        end if




        do i_spin = 1,n_spin

           work(1:n_centers_basis_I, 1:n_basis) = work_ham(1:n_centers_basis_I,1:n_basis,i_spin)      


           do i_basis_2 = 1+n_basis,n_centers_basis_I,1 

              work(1:n_centers_basis_I, Cbasis_to_basis(i_basis_2)) =   &
                   work(1:n_centers_basis_I, Cbasis_to_basis(i_basis_2)) &
                   + work_ham(1:n_centers_basis_I, i_basis_2,i_spin) &
                   *  dconjg(k_phase(center_to_cell(Cbasis_to_center( i_basis_2 )),k_point))
           end do


           do i_basis_1 = 1+n_basis,n_centers_basis_I,1 

              work(Cbasis_to_basis(i_basis_1), 1:n_basis) =   &
                   work(Cbasis_to_basis(i_basis_1), 1:n_basis) &
                   + work(i_basis_1, 1:n_basis) &
                   *  k_phase(center_to_cell(Cbasis_to_center( i_basis_1 )),k_point)
           end do


           do i_basis_2 = 1,n_basis

              hamiltonian_w_complex( index_b_new( 1:i_basis_2, i_basis_2) ,i_spin) =  &
                   + work(1:i_basis_2, i_basis_2)

           end do
        end do


        call aims_deallocate( index_b_old, "index_b_old" )
        call aims_deallocate( index_b_new, "index_b_new" )


     end if

  case default
     
     call localorb_info('Invalid packing type')
     call aims_stop

  end select

  ! Allocatable array that are tracked (probably not needed here, but why not)
  if (allocated(index_b))     call aims_deallocate( index_b,         "index_b" )
  if (allocated(index_b_old)) call aims_deallocate( index_b_old, "index_b_old" )
  if (allocated(index_b_new)) call aims_deallocate( index_b_new, "index_b_new" )
  if (allocated(work))        call aims_deallocate( work,               "work" )

end subroutine construct_hamiltonian
!******
