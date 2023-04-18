module plumed_new_interface

  real*8, dimension(:,:), allocatable :: MD_forces
  real*8, dimension(:), allocatable :: masses
  real*8, dimension(:), allocatable :: charges
  real*8, dimension(:), allocatable :: px,py,pz
  real*8, dimension(:), allocatable :: fx,fy,fz
  real*8 :: plumed_energy,plumed_eunit
  integer :: cell_type, i_atom_plumed

contains

subroutine plumed_new_init ( )
  use runtime_choices
  use physics
  use geometry
  use constants
  use species_data
  use dimensions
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
  use timing
  use relaxation
  real*8 :: step
  allocate(MD_forces(3,n_atoms))
  allocate(masses(n_atoms))
  allocate(charges(n_atoms))
  allocate(px(n_atoms))
  allocate(py(n_atoms))
  allocate(pz(n_atoms))
  allocate(fx(n_atoms))
  allocate(fy(n_atoms))
  allocate(fz(n_atoms))

  if ( n_periodic .ge.1) then
    cell_type=1
  else
    cell_type=0
  endif
  plumed_eunit=1.d0
  do i_atom_plumed = 1, n_atoms
    masses(i_atom_plumed) = species_m(species(i_atom_plumed))
    charges(i_atom_plumed) = 0.d0
  enddo
  if (use_geo_relaxation) then
    step = 1.d0 
  else if (use_molecular_dynamics) then
    step = MD_tstep
  else 
    stop
  endif

  if(myid.eq.0) then
    call init_metadyn_(n_atoms, step, masses, charges,  cell_type , plumed_eunit ,  trim(adjustl(plumed_file))//char(0))
  endif
end subroutine plumed_new_init


subroutine plumed_forces ()
  use runtime_choices
  use physics
  use geometry
  use constants
  use species_data
  use dimensions
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
  use timing
  use relaxation
  integer :: counter

  if (use_geo_relaxation) then
    call clean_force_components(total_forces, forces_lv)
  endif
  if(myid.eq.0) then
  do i_atom_plumed = 1, n_atoms
!        masses(i_atom_plumed) = species_m(species(i_atom_plumed))
!        charges(i_atom_plumed) = 0.d0
        px(i_atom_plumed)=coords(1,i_atom_plumed)
        py(i_atom_plumed)=coords(2,i_atom_plumed)
        pz(i_atom_plumed)=coords(3,i_atom_plumed)
        write(use_unit,*) 'P ',px(i_atom_plumed),py(i_atom_plumed),pz(i_atom_plumed)
        fx(i_atom_plumed)=total_forces(1,i_atom_plumed)
        fy(i_atom_plumed)=total_forces(2,i_atom_plumed)
        fz(i_atom_plumed)=total_forces(3,i_atom_plumed)
        write(use_unit,*) 'F ',fx(i_atom_plumed),fy(i_atom_plumed),fz(i_atom_plumed)
  end do
  if (use_geo_relaxation) then
    counter = relaxation_step_number
  else if (use_molecular_dynamics) then
    counter = MD_stepcount
  else 
    stop
  endif

  call meta_force_calculation_( lattice_vector, counter, px, py, pz , fx, fy, fz,  plumed_energy)
  do i_atom_plumed = 1, n_atoms
          total_forces(1,i_atom_plumed)=fx(i_atom_plumed)
          total_forces(2,i_atom_plumed)=fy(i_atom_plumed)
          total_forces(3,i_atom_plumed)=fz(i_atom_plumed)
  end do
  endif
  call broadcast_MD_velocities(total_forces,0)
!        MD_forces(:,:) = MD_forces(:,:)/MD_KE_factor    
  if (use_geo_relaxation) then
    call clean_force_components(total_forces, forces_lv)
  endif
end subroutine plumed_forces

end module plumed_new_interface
