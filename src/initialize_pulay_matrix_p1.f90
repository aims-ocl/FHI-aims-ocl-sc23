!****s* FHI-aims/initialize_pulay_matrix_p1
!  NAME
!   initialize_pulay_matrix_p1
!  SYNOPSIS

subroutine initialize_pulay_matrix_p1( partition_tab, &
     pulay_saved_iter, previous_rho_error, pulay_matrix )

!  PURPOSE
!  Subroutine initialize_pulay_matrix
!
!  constructs the Pulay matrix from the previously stored densities 
!  from linear mixing.
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use species_data
  use xc
  use lapack_wrapper
  use mpi_utilities
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: partition_tab
  integer :: pulay_saved_iter
  real*8, dimension(n_full_points, n_max_pulay, n_spin) :: previous_rho_error
      
  real*8, dimension((n_max_pulay),(n_max_pulay)) :: pulay_matrix

!  INPUTS
!  o partition_tab -- values of partition function
!  o pulay_saved_iter -- number of iterations saved to Pulay matrix previously
!  o previous_rho_error -- error in electron density from previous iterations.
!
!  OUTPUT
!  o pulay_matrix -- Pulay matrix
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





!  counters

  integer i_index
  integer i_my_batch
  integer i_spin

  integer i_int_grid
  integer i_pulay_store, i_pulay_store_2

  integer :: i_offset
  
!  begin work


  if (myid.eq.0) then
     write(use_unit,'(2X,A,A)') &
          "Switching from linear to Pulay mixing: ", &
          "Initializing Pulay matrix."
  end if
      
  pulay_matrix = 0.d0

!     calculate, use, and store density error err(rho_i) for present iteration i
!     between KS input and output densities from previous iterations n-i, n-j ...
  i_offset = 0
  do i_my_batch = 1, n_my_batches, 1

        ! loop over one batch
        do i_index = 1, batches(i_my_batch)%size, 1           
!           execute only if partition_tab.gt.0 here, i.e. if the integration point
!           makes sense
           i_offset = i_offset + 1
           if (partition_tab(i_offset).gt.0.d0) then

!             determine Pulay matrix elements: < err(rho_i)-err(rho_i+1) | err(rho_j)-err(rho_j+1) >
    
              do i_pulay_store = 1, pulay_saved_iter-1, 1
                 do i_pulay_store_2 = i_pulay_store, pulay_saved_iter-1, 1

                    do i_spin = 1, n_spin, 1

                       pulay_matrix(i_pulay_store,i_pulay_store_2) = &
                            pulay_matrix(i_pulay_store,i_pulay_store_2) + &
                            partition_tab(i_offset) * &
                            ( previous_rho_error(i_offset, i_pulay_store, i_spin) - &
                            previous_rho_error(i_offset, i_pulay_store+1, i_spin))* &
                            ( previous_rho_error(i_offset, i_pulay_store_2, i_spin) - &
                            previous_rho_error(i_offset, i_pulay_store_2+1, i_spin))

                    end do

                 enddo
              enddo

            end if

!         end angular integration
         end do
!     end distribution over threads
      ! end if
      !       end radial integration loop
   end do
!     end loop over atoms
 end subroutine initialize_pulay_matrix_p1
 !******		
