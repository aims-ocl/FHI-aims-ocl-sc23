!****s* FHI-aims/output_potential_p1
!  NAME
!    output_potential_p1
!  SYNOPSIS

subroutine output_effective_potential( hartree_potential, rho, rho_gradient, kinetic_density, &
           free_hartree_potential, rho_free, rho_gradient_free, partition_tab )

!  PURPOSE
!  Subroutine output_potential
!
!  Writes the effective potential on the integration grid. In practise it prints out the Hartree potential part.
!
!  Note: This is not the effective potential for  a real fix, we must evaluate the 
!  xc potential right here so that we can actually write out the relevant bits.
!
!  USES

  use constants, only: pi4_inv
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use mpi_tasks
  use mpi_utilities
  use localorb_io, only: use_unit
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: hartree_potential
  real*8, dimension(n_full_points) :: free_hartree_potential
  real*8, dimension(n_spin, n_full_points) :: rho
  real*8, dimension(3, n_spin, n_full_points) :: rho_gradient
  real*8, dimension(n_spin, n_full_points) :: kinetic_density
  real*8, dimension(n_spin, n_full_points) :: rho_free
  real*8, dimension(3, n_spin, n_full_points) :: rho_gradient_free
  real*8, dimension(n_full_points) :: partition_tab

!  INPUTS
!   o local_potential_parts -- potential which is printed out.
!  OUTPUT
!    none
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

  real*8 grid_coord(3)
  real*8:: local_potential_parts(n_spin), en_density_xc, xc_gradient_deriv(3,3)
  real*8:: local_potential_parts_free(n_spin), rhoF, rho_gradient_F(3,n_spin)
  real*8:: xc_tau_deriv(n_spin)

!  counters

!  integer i_atom, i_atom_2
!  integer i_radial
!  integer i_angular
!  integer i_grid
  integer i_coord
  integer i_spin, i_task, mpierr
  
!  integer i_species

  integer :: i_point, i_my_batch, i_index

  !  begin work
  
  call get_my_task()

  if (myid.eq.0) then
     write(use_unit,'(2X,A,A)') &
          "Writing the effective potential at each grid point", &
          " to file v_eff.dat ."

  end if


  do i_task = 0, n_tasks-1

 !    write(use_unit,*) 'task', i_task, myid
     if(myid == i_task)then
        
 !       write(use_unit,*) 'task sis', i_task, myid

        if(i_task == 0)then
           open (50, file="effective_potential.dat")
        else
!           open (50, file="effective_potential.dat", STATUS='OLD', ACCESS='DIRECT' )
           open (50, file="effective_potential.dat" , status='OLD', position='APPEND')
        end if

        !     calculate potential, atom by atom, and point by point 
        i_point = 0
        do i_my_batch = 1, n_my_batches, 1



           do i_index = 1, batches(i_my_batch)%size, 1

              i_point = i_point + 1

              if (partition_tab(i_point).gt.0.d0) then


              !           calculate grid point coordinate
              grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)



              if(n_periodic > 0)then
                 call map_to_center_cell(grid_coord(1:3) )
              end if




              call evaluate_xc  &
                   ( rho(1,i_point),  rho_gradient(1,1,i_point), kinetic_density(1,i_point), &
                   en_density_xc, local_potential_parts(1), xc_gradient_deriv(1,1), xc_tau_deriv(1), &
                   .false., grid_coord(:) )
                   
              rhoF =   pi4_inv * rho_free(1,i_point)
!              rho_gradient_F(1:3, 1:n_spin) = rho_gradient_free(1:3, 1:n_spin,i_point)

              ! I'm not sure if kinetic_density should be replaced with a dummy argument here?
              ! But then what if we are using a meta-gga for the effective potential? 
              ! The perfect solution would be to calculate the kinetic_density for the free atom.... AJL
              call evaluate_xc  &
                   ( rhoF,  rho_gradient_free(1,1,i_point), kinetic_density(1,i_point), &
                   en_density_xc, local_potential_parts_free(1), xc_gradient_deriv(1,1), & 
                   xc_tau_deriv(1), .false.)




              local_potential_parts(1) = local_potential_parts(1) +  hartree_potential(i_point) ! -  free_hartree_potential(i_point) - local_potential_parts_free(1)
              if(n_spin ==2) local_potential_parts(2) = local_potential_parts(2) &
                   +  hartree_potential(i_point) !-  free_hartree_potential(i_point)- local_potential_parts_free(2)


              local_potential_parts_free(1) = local_potential_parts_free(1) +  free_hartree_potential(i_point) 
              if(n_spin ==2) local_potential_parts_free(2) = local_potential_parts_free(2) &
                   +  free_hartree_potential(i_point)




!              if (myid.eq.0) then
                 write(50,'(7(F20.10,1X))') &
                      (grid_coord(i_coord)*bohr, i_coord=1,3,1), &
                      ( (local_potential_parts(i_spin) - local_potential_parts_free(i_spin))*hartree, i_spin = 1, n_spin, 1 ), &
                      ( local_potential_parts(i_spin)*hartree, i_spin = 1, n_spin, 1 ) 

!              end if

              end if
              ! end loop over a batch
              enddo

           !end if
           ! end loop over batche
        enddo

        close(50)

     end if

!     write(use_unit,*) 'valli', i_task, myid
     call MPI_BARRIER(mpi_comm_global, mpierr)
!     write(use_unit,*) 'end', i_task, myid
  end do



!  stop
  return
end subroutine output_effective_potential
!----------------------------------------------------------------------
!******	
