	Program Fit_esp_charges
	
	use mpi_tasks
        implicit none

	real*8, dimension(:,:), allocatable :: atoms
	real*8, dimension(:,:), allocatable :: atoms_bader
	real*8,  dimension(:,:), allocatable :: a
	integer,  dimension(:), allocatable :: work
	real*8,  dimension(:), allocatable :: b
	real*8,  dimension(:), allocatable :: out_charge
	integer, dimension(:), allocatable  :: species
	real*8, dimension(:), allocatable  :: radius_esp_min_new
	real*8, dimension(:), allocatable  :: radius_esp_max_new
	real*8, dimension(:,:), allocatable  :: esp_atom_constraint
	real*8, dimension(:), allocatable  :: esp_charges
        character*2, dimension(:), allocatable :: species_name
        integer, dimension(:), allocatable :: species_z
	real*8, dimension(:,:), allocatable  :: coords
	real*8, dimension(:,:), allocatable  :: coords_geo
	real*8, dimension(:), allocatable  :: potential
        real*8, dimension (:,:,:), allocatable :: DensTotal
	real*8, dimension(:,:), allocatable  :: r_sum, k_sum
	integer, dimension(:), allocatable  :: task_points
	real*8, dimension(:), allocatable  :: esp_cube_output_org
	real*8, dimension(:), allocatable  :: esp_cube_output
	real*8, dimension(:), allocatable  ::cube_output
	logical, dimension(:), allocatable  :: use_point
        real*8, dimension(:,:), allocatable  :: map_to_center_cell_matrix
        logical :: found
	integer :: i_atom, j_atom, i_full_points, i_x ,i_y, i_z, counter_2, &
                   counter,i_points,i_species, i, j, k, m, task_point
	integer :: info, n_x, n_y, n_z, all_points, cells_atoms, n_species, n_atoms, n_points,max_points
	real*8  :: pot_sum, total_diff, diff,&
                   cell_volume, dist
	real*8  :: coord_current(3), d_r(3),coord_new(3),off(3),esp_dipole(3)
	real*8  :: lattice_vector(3,3),lattice_vector_cube(3,3), recip_lattice_vector(3,3)
	real*8  :: coord(3,4)
        character (100) :: input
        character (30) :: filename,filename_two
	character(30) :: tmp
        real*8 :: beta
        Character (len=30) :: Dummy
	
        real*8, parameter  :: bohr = 0.5291772085900001D0
        real*8, parameter  :: hartree = 27.21138386D0
        real*8, parameter  :: pi = 4.D0*DATAN(1.D0)
	real*8, parameter  :: sqrt_pi = sqrt(pi)
	real*8, parameter  :: pi4 = pi*4d0
        real*8  :: alpha
	real*8  :: R_c
	integer :: rm
	integer :: km
	real*8  :: charge
	integer :: n_periodic = 3
	integer :: esp_constraint = 1
        logical :: output_esp = .True.
	real*8  :: esp_min = 1.1
	real*8  :: esp_max = 4.0
        logical :: write_cube  = .False.


	call getarg(1,filename)
        call getarg(2,input)
        input = trim(input)
        read (input,*) esp_min
        call getarg(3,input)
        input = trim(input)
        read (input,*) esp_max
        call getarg(4,input)
        input = trim(input)
        read (input,*) R_c
        R_c = R_c/bohr
        alpha  = ((3.2/R_c)**2)
        call getarg(5,input)
        input = trim(input)
        read (input,*) rm
        call getarg(6,input)
        input = trim(input)
        read (input,*) km
        call getarg(7,input)
        input = trim(input)
        read (input,*) charge
        call getarg(8,filename_two)
	if (filename_two.ne.'                              ')then
	  esp_constraint = 2
	endif

        write(*,*) 'Fitting ESP charges starts (method 2)'
        write(*,*) 'Reading structural data and potential from: ', filename   
	Open (unit=1, file=filename, status='old')

        read (1,*) Dummy
        read (1,*) Dummy

	read (1,*) n_atoms, off
	do i = 1, 3 ,1
	  read(1,*) (coord(i,j),j = 1, 4, 1)
	enddo
	        
	allocate(atoms(n_atoms,5))
	do i = 1, n_atoms ,1
	  read(1,*) (atoms(i,j),j = 1, 5, 1)
	enddo
	
        n_x = coord(1,1)
        n_y = coord(2,1)
        n_z = coord(3,1)

	allocate (DensTotal(n_x,n_y,n_z))
	do i=1,n_x
		do j=1,n_y
			read(1,*) (DensTotal(i,j,k), k=1,n_z)
		enddo
	enddo
	Close(1)
        if (esp_constraint.eq.2)then
          write(*,*) 'Reading fitting contraints for (method 2)'
	  Open (unit=2, file=filename_two, status='old')

	  read (2,*) Dummy
	  read (2,*) Dummy
	  allocate(atoms_bader(n_atoms,7))
	  do i = 1, n_atoms ,1
	    read(2,*) (atoms_bader(i,j),j = 1, 7, 1)
	  enddo
	  Close(2)
          call getarg(5,tmp)
          tmp = trim(tmp)
          read (tmp,*) beta
        else
          allocate(atoms_bader(1,1))
        endif
        do i = 1, 3 ,1
          lattice_vector_cube(:,i) = coord(i,2:4)*coord(i,1)
        enddo
        allocate(coords(3,n_atoms),stat=info)      
        do i = 1, n_atoms ,1 
           coords(1:3,i) = atoms(i,3:5)
        enddo
        n_points = n_x*n_y*n_z
	allocate (potential(n_points))
        potential = 0d0
        m=0
	do i=1,n_x
		do j=1,n_y
		       do k=1,n_z
                         m=m+1
                         potential(m)=DensTotal(i,j,k)
                      enddo
		enddo
	enddo
        if(allocated(DensTotal)) deallocate(DensTotal)

  allocate(species(n_atoms),stat=info) 
  species = 1 
  i_species = 1
  species(1) = atoms(1,1)
  do i = 1, n_atoms ,1
    do j = 1, i_species, 1
      if (atoms(i,1).eq.species(j))then
        found = .true.
      endif
    enddo
    if (.not.found) then
      i_species = i_species + 1
      species(i_species) = atoms(i,1)
    endif
    found = .false.
  enddo
  n_species = i_species
  allocate(species_z(n_species),stat=info) 
  allocate(species_name(n_species),stat=info) 
  do i_species = 1, n_species, 1
    species_z(i_species) = species(i_species)
  enddo
  do i = 1, n_atoms ,1 
    do i_species = 1, n_species, 1
      if (atoms(i,1).eq.species_z(i_species))then
        species(i) = i_species
      endif
    enddo
  enddo
  if (esp_constraint.eq.2)then
    allocate(esp_atom_constraint(2,n_atoms),stat=info)
    do i = 1, n_atoms ,1
      esp_atom_constraint(1,i)=species_z(species(i))-atoms_bader(i,5)
      esp_atom_constraint(2,i)= beta
    enddo
  else
    allocate(esp_atom_constraint(1,1),stat=info)

  endif
  call esp_name(species_name,species_z,n_species)
  allocate(coords_geo(3,n_atoms),stat=info)
  write(*,*) 'Reading geometry.in to obtain lattice vectors.'
  call read_geo(n_atoms,n_species,n_periodic, coords_geo,lattice_vector,&
                species_name,esp_atom_constraint)
  if (n_periodic > 0) then
    cells_atoms = 27*n_atoms-1
    call get_reciprocal_vectors(n_periodic, lattice_vector, &
    &                                 recip_lattice_vector)
    call get_cell_volume(lattice_vector, cell_volume)
  else
    cells_atoms = (n_atoms-1)
  endif


  call initialize_mpi()

  allocate(task_points(n_tasks+1),stat=info)
  task_points(1)=0
  do i=1, n_tasks-1, 1
    task_points(i+1)=i*floor((1d0*n_points)/(1d0*n_tasks))
  enddo
  task_points(n_tasks+1)=n_points
  max_points = floor(1d0*n_points/n_tasks)+(n_points-n_tasks*floor(1d0*n_points/n_tasks))

  if (n_periodic > 0) then
    allocate(a(n_atoms+2,n_atoms+2),stat=info)
    allocate(b(n_atoms+2),stat=info)
    allocate(k_sum(n_atoms,max_points),stat=info) 
  else
    allocate(a(n_atoms+1,n_atoms+1),stat=info)
    allocate(b(n_atoms+1),stat=info)
    allocate(k_sum(1,1),stat=info) 
  endif
  allocate(r_sum(n_atoms,max_points),stat=info) 
  allocate(use_point(max_points),stat=info)
  pot_sum = 0d0
  all_points = 0
  i_full_points = 0
  task_point = 0
  a = 0d0
  b = 0d0
  counter = 0
  allocate(radius_esp_min_new(n_species),stat=info) 
  allocate(radius_esp_max_new(n_species),stat=info) 
  call esp_radius(radius_esp_min_new, radius_esp_max_new, esp_min, esp_max,&
                  bohr,species_z,n_species)
!call sync_vector(potential,n_points,MPI_COMM_WORLD)
  if (myid==0)then  
    write(*,*) 'Total Number of grid points:', n_points
  endif
  r_sum = 0d0
  k_sum = 0d0
  use_point = .False.
  if (n_periodic.ne.0)then
    allocate(map_to_center_cell_matrix(n_periodic,3),stat=info) 
    call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
    &                                      map_to_center_cell_matrix)
  else
    allocate(map_to_center_cell_matrix(1,1),stat=info) 
  endif
  do i_x = 0, n_x-1, 1 
  do i_y = 0, n_y-1, 1
  do i_z = 0, n_z-1, 1
    d_r(1) = dble(i_x)/n_x
    d_r(2) = dble(i_y)/n_y
    d_r(3) = dble(i_z)/n_z
    i_full_points = i_full_points + 1
    if (i_full_points.gt.task_points(myid+1).and.i_full_points.le.task_points(myid+2))then
      task_point = task_point + 1
      if (mod(task_point,50000)==0) then
        write(*,*) 'Task:', myid, ', ', task_point, ' of ', task_points(myid+2)-task_points(myid+1), ' points done.'
      endif
      coord_new(1:3) = matmul(lattice_vector_cube,d_r(1:3)) +off
      call check_point(1,lattice_vector,radius_esp_min_new,&
		      coord_new,counter_2, n_atoms, n_periodic, &
		      n_species,species,coords,map_to_center_cell_matrix)
      if (counter_2.eq.(cells_atoms+1))then
            use_point(task_point) = .True.
	    counter=counter+1
	    all_points = all_points + 1
	    coord_current(:) = coord_new(1:3)
	    pot_sum = pot_sum - potential(i_full_points)
	    do i_atom = 1, n_atoms, 1
              if(n_periodic.ne.0)then
		call rsummation(coord_current,coords(1:3,i_atom),rm,R_c,&
			      lattice_vector,alpha,r_sum(i_atom,task_point))
		call ksummation(coord_current,coords(1:3,i_atom),km,R_c,&
			      recip_lattice_vector,alpha,k_sum(i_atom,task_point))
		b(i_atom) = b(i_atom)+((-potential(i_full_points))&
				*(r_sum(i_atom,task_point)+(pi4/cell_volume)*&
				k_sum(i_atom,task_point)))
		a(n_atoms+1,i_atom) = a(n_atoms+1,i_atom) +&
					(r_sum(i_atom,task_point)+(pi4/cell_volume)*k_sum(i_atom,task_point))
		a(i_atom,n_atoms+1) = a(i_atom,n_atoms+1) + &
					(r_sum(i_atom,task_point)+(pi4/cell_volume)*k_sum(i_atom,task_point))
              else
                call rsummation_cluster(coord_current,coords(1:3,i_atom),r_sum(i_atom,task_point))
                b(i_atom) = b(i_atom)-r_sum(i_atom,task_point)*potential(i_full_points)
              endif
	    enddo
	  
      endif
    endif    
  enddo
  enddo
  enddo

  i_full_points = 0
  task_point = 0
  do i_full_points =  1, n_points, 1
    if (i_full_points.gt.task_points(myid+1).and.i_full_points.le.task_points(myid+2))then
       task_point = task_point + 1 
       if(use_point(task_point))then
	  do i_atom = 1, n_atoms, 1
	      do j_atom = 1, n_atoms, 1
                if(n_periodic.ne.0)then
		    a(i_atom,j_atom) = a(i_atom,j_atom) + &
				      (r_sum(i_atom,task_point)+(pi4/cell_volume)*k_sum(i_atom,task_point))*&
				      (r_sum(j_atom,task_point)+(pi4/cell_volume)*k_sum(j_atom,task_point))
                else
                    a(i_atom,j_atom) = a(i_atom,j_atom) + &
                                   r_sum(i_atom,task_point)*r_sum(j_atom,task_point)
                endif
	      enddo
	  enddo
       endif
    endif    
  enddo
  if(n_periodic.ne.0)then
    a(n_atoms+1,n_atoms+1) = dble(all_points)
    b(n_atoms+1) = pot_sum
    call sync_vector(b,n_atoms+2,MPI_COMM_WORLD)
    call sync_matrix(a,n_atoms+2 ,n_atoms+2, MPI_COMM_WORLD)
    call sync_real_number(pot_sum,MPI_COMM_WORLD)
    if (myid==0)then  
      write(*,*) 'Number of reduced grid points:', int(a(n_atoms+1,n_atoms+1))
    endif
    a(n_atoms+2,1:n_atoms) = 1d0
    a(1:n_atoms,n_atoms+2) = 1d0
    a(n_atoms+1,n_atoms+2) = 0d0
    a(n_atoms+2,n_atoms+1) = 0d0
    a(n_atoms+2,n_atoms+2) = 0d0
    b(n_atoms+2) = charge
    if (esp_constraint.eq.2)then
	  write(*,'(2X,A)') &
			"| Using supplied constraints for method 2"
	  do i_atom = 1, n_atoms
	    a(i_atom,i_atom) = a(i_atom,i_atom) + &
			      esp_atom_constraint(2,i_atom)
	    b(i_atom) = b(i_atom) + &
				esp_atom_constraint(1,i_atom)*&
				esp_atom_constraint(2,i_atom)
	  enddo
    endif

    allocate(work(n_atoms+2),stat=info)
    if (myid==0)then
	call dgesv(n_atoms+2, 1, a, n_atoms+2, work, b, n_atoms+2 ,info)
	if( info.gt.0 ) then
	  write(*,'(2X,A)')'The diagonal element of the triangular factor of A,'
	  write(*,'(2X,A,I8,A,I8,A)')'U(',info,',',info,') is zero, so that'
	  write(*,'(2X,A)')'A is singular; the solution could not be computed.'
	endif
    endif
   else
      call sync_real_number(1d0*all_points,MPI_COMM_WORLD)
      if (myid==0)then  
	write(*,*) 'Number of reduced grid points:', int(all_points)
      endif
      call sync_vector(b,n_atoms+1,MPI_COMM_WORLD)
      call sync_matrix(a,n_atoms+1 ,n_atoms+1,MPI_COMM_WORLD )
      a(n_atoms+1,1:n_atoms) = 1d0
      a(1:n_atoms,n_atoms+1) = 1d0
      b(n_atoms+1) = charge
      allocate(work(n_atoms+1),stat=info)
      if (myid==0)then
	call dgesv(n_atoms+1, 1, a, n_atoms+1, work, b, n_atoms+1 ,info)
	if( info.gt.0 ) then
	  write(*,'(2X,A)')'The diagonal element of the triangular factor of A,'
	  write(*,'(2X,A,I8,A,I8,A)')'U(',info,',',info,') is zero, so that'
	  write(*,'(2X,A)')'A is singular; the solution could not be computed.'
	endif
      endif
   endif
   allocate(out_charge(n_atoms+1),stat=info)
   out_charge = 0d0
   if (myid==0)then
         do i_atom = 1, n_atoms, 1
            out_charge(i_atom)=b(i_atom)
         enddo
         if(n_periodic.ne.0)then
           out_charge(n_atoms+1)=b(n_atoms+1)
         endif
   endif
   call sync_vector(out_charge,n_atoms+1,MPI_COMM_WORLD)
   if (myid==0.and.output_esp)then
            write(*,'(2X,A)') "| "
	    write(*,'(2X,A)') "| ESP charges fitted to the electrostatic potential : "
	    do i_atom = 1, n_atoms, 1
		write(*,'(2X,A)') "| "
		write(*,'(2X,A,1X,I8,A,A,A,F15.8,F15.8,F15.8,F15.8)') '| Atom ', i_atom&
                   , ": ",trim(species_name(species(i_atom))),", ESP charge: ",&
                   out_charge(i_atom)
	    end do
            if(n_periodic.ne.0)then
	      write(*,'(2X,A)') "| "
	      write(*,'(2X,A,1X,F15.8,A)') '| Potential offset: ', &
		    out_charge(n_atoms+1)*hartree, ' eV'
            endif
	    write(*,'(2X,A)') " "
   endif

   if(write_cube)then
     allocate( esp_cube_output(n_z),stat=info)
     allocate( esp_cube_output_org(n_z),stat=info)
     allocate(cube_output(n_points),stat=info)
   endif
   if (output_esp) then
      !Calculate RMS for fitted charges
      total_diff = 0
      i_full_points = 0
      i_points = 0
      pot_sum = 0d0
      task_point = 0
      if(write_cube)then
        cube_output = 0d0
        esp_cube_output = 0d0
        esp_cube_output_org = 0d0
      endif
      do i_x = 0, n_x-1, 1 
      do i_y = 0, n_y-1, 1
      do i_z = 0, n_z-1, 1
	d_r(1) = dble(i_x)/n_x
	d_r(2) = dble(i_y)/n_y
	d_r(3) = dble(i_z)/n_z
        i_full_points = i_full_points + 1
        if (i_full_points.gt.task_points(myid+1).and.i_full_points.le.task_points(myid+2))then
          task_point = task_point + 1
	  if (use_point(task_point))then

              diff = 0d0
              do i_atom = 1, n_atoms, 1
                 if(n_periodic.ne.0)then
                   diff = diff +  (r_sum(i_atom,task_point)+(pi4/cell_volume)*&
				k_sum(i_atom,task_point))*out_charge(i_atom)
                 else
                   diff = diff +  (r_sum(i_atom,task_point))*out_charge(i_atom)
                 endif
              enddo
              if(write_cube)then
                   if(n_periodic.ne.0)then
                     esp_cube_output(i_z+1) = esp_cube_output(i_z+1)-diff - out_charge(n_atoms+1)
                   else
                     esp_cube_output(i_z+1) = esp_cube_output(i_z+1)-diff
                   endif
                   esp_cube_output_org(i_z+1)= esp_cube_output_org(i_z+1) + potential(i_full_points)
                   cube_output(i_full_points) = -diff
              endif
              if(n_periodic.ne.0)then
                total_diff = total_diff + (( diff + out_charge(n_atoms+1)) +  &
                                         potential(i_full_points))**2
              else
                total_diff = total_diff + (( diff ) +  &
                                         potential(i_full_points))**2
              endif
              pot_sum = pot_sum + potential(i_full_points)**2
          endif
	endif
      enddo
      enddo
      enddo

      call sync_real_number(total_diff,MPI_COMM_WORLD)
      call sync_real_number(pot_sum,MPI_COMM_WORLD)
      if(write_cube)then
        call sync_vector(esp_cube_output,n_z,MPI_COMM_WORLD)
        call sync_vector(esp_cube_output_org,n_z,MPI_COMM_WORLD)
        call sync_vector(cube_output,n_points,MPI_COMM_WORLD)
        if(myid==0)then
	   call esp_output_zcube(filename,esp_min,esp_max,esp_cube_output,esp_cube_output_org,n_z)
          call esp_output_cubefile(cube_output,n_x,n_y,n_z,&
                               species_z,n_atoms,n_periodic,n_species,&
                               species,coords,lattice_vector_cube,&
                               off,radius_esp_max_new,2)
        endif
      endif
      endif
      allocate(esp_charges(n_atoms))
      if(myid==0)then
	do i_atom = 1, n_atoms, 1
	  esp_charges(i_atom) = out_charge(i_atom)
	enddo
	esp_dipole = 0
	do i_atom = 1, n_atoms, 1
		esp_dipole(1)= esp_dipole(1) - out_charge(i_atom)*&
			      (coords(1,i_atom))
		esp_dipole(2)= esp_dipole(2) - out_charge(i_atom)*&
			      (coords(2,i_atom))
		esp_dipole(3)= esp_dipole(3) - out_charge(i_atom)*&
			      (coords(3,i_atom))
	enddo
	write(*,'(2X,A,1X,F15.8,F15.8,F15.8)') "Dipole matrix element: "&
				    , esp_dipole(1), esp_dipole(2), esp_dipole(3)
	write(*,'(2X,A,1X,F15.8,F15.8,F15.8)')  "RRMS: ", &
						    sqrt(total_diff/pot_sum)
	! Write one last line to bound output
	write(*,*) ' '
      endif
  if(allocated(a)) deallocate(a)
  if(allocated(b)) deallocate(b)
  if(allocated(work)) deallocate(work)
  if(allocated(out_charge)) deallocate(out_charge)
  if(allocated(potential)) deallocate(potential)
  if(allocated(atoms)) deallocate(atoms)
  if(allocated(atoms_bader)) deallocate(atoms_bader)
  if(allocated(species)) deallocate(species)
  if(allocated(radius_esp_min_new)) deallocate(radius_esp_min_new)
  if(allocated(radius_esp_max_new)) deallocate(radius_esp_max_new)
  if(allocated(esp_atom_constraint)) deallocate(esp_atom_constraint)
  if(allocated(esp_charges)) deallocate(esp_charges)
  if(allocated(species_name)) deallocate(species_name)
  if(allocated(species_z)) deallocate(species_z)
  if(allocated(coords)) deallocate(coords)
  if(allocated(coords_geo)) deallocate(coords_geo)
  if(allocated(task_points)) deallocate(task_points)
  if(write_cube)then
    if(allocated(esp_cube_output)) deallocate(esp_cube_output)
    if(allocated(esp_cube_output_org)) deallocate(esp_cube_output_org)
    if(allocated(cube_output)) deallocate(cube_output)
  endif
  if(n_periodic.ne.0)then
    if(allocated(map_to_center_cell_matrix)) deallocate(map_to_center_cell_matrix)
  endif
  call finalize_mpi
  contains
  
subroutine esp_radius(radius_esp_min, radius_esp_max, esp_min, esp_max,bohr,&
                      species_z,n_species)
!  PURPOSE
!   Generate the points where the potential will be evaluated. 
!    vdw-Radius of atoms +radius < points < vdw-Radius of atoms +radius2
!  USES
!  ARGUMENTS
!  INPUTS
!   o esp_min -- factor for max radius used for ESP fitting
!   o esp_max -- factor for min radius used for ESP fitting
!  OUTPUT
!   o radius_esp_min -- esp_min*vdw_radius
!   o radius_esp_max -- esp_max*vdw_radius
!   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  real*8,  dimension(n_species) :: radius_esp_min
  real*8,  dimension(n_species) :: radius_esp_max
  real*8 :: esp_min
  real*8 :: esp_max
  real*8, dimension(:),       allocatable :: Atom_vdw
  integer :: i_species, info
  real*8 :: bohr
  integer :: species_z(n_species)
  integer :: n_species

  if(.not. allocated( Atom_vdw))then
     allocate( Atom_vdw(103),stat=info)
  end if
  ! vdw radius for all atoms, numbers have been taken from wikipedia
  Atom_vdw(  1) = 1.10/bohr
  Atom_vdw(  2) = 1.40/bohr
  Atom_vdw(  3) = 1.82/bohr
  Atom_vdw(  4) = 1.53/bohr
  Atom_vdw(  5) = 1.92/bohr
  Atom_vdw(  6) = 1.70/bohr
  Atom_vdw(  7) = 1.55/bohr
  Atom_vdw(  8) = 1.52/bohr
  Atom_vdw(  9) = 1.47/bohr
  Atom_vdw( 10) = 1.54/bohr
  Atom_vdw( 11) = 2.27/bohr
  Atom_vdw( 12) = 1.73/bohr
  Atom_vdw( 13) = 1.84/bohr
  Atom_vdw( 14) = 2.10/bohr
  Atom_vdw( 15) = 1.80/bohr
  Atom_vdw( 16) = 1.80/bohr
  Atom_vdw( 17) = 1.75/bohr
  Atom_vdw( 18) = 1.88/bohr
  Atom_vdw( 19) = 2.75/bohr
  Atom_vdw( 20) = 2.31/bohr
  Atom_vdw( 21) = 1.00
  Atom_vdw( 22) = 1.00
  Atom_vdw( 23) = 1.00
  Atom_vdw( 24) = 1.00
  Atom_vdw( 25) = 1.00
  Atom_vdw( 26) = 1.00
  Atom_vdw( 27) = 1.00
  Atom_vdw( 28) = 1.63/bohr
  Atom_vdw( 29) = 1.40/bohr
  Atom_vdw( 30) = 1.39/bohr
  Atom_vdw( 31) = 1.87/bohr
  Atom_vdw( 32) = 2.11/bohr
  Atom_vdw( 33) = 1.85/bohr
  Atom_vdw( 34) = 1.90/bohr
  Atom_vdw( 35) = 1.85/bohr
  Atom_vdw( 36) = 2.02/bohr
  Atom_vdw( 37) = 3.03/bohr
  Atom_vdw( 38) = 2.49/bohr
  Atom_vdw( 39) = 1.00
  Atom_vdw( 40) = 1.00
  Atom_vdw( 41) = 1.00
  Atom_vdw( 42) = 1.00
  Atom_vdw( 43) = 1.00
  Atom_vdw( 44) = 1.00
  Atom_vdw( 45) = 1.00
  Atom_vdw( 46) = 1.63/bohr
  Atom_vdw( 47) = 1.72/bohr
  Atom_vdw( 48) = 1.58/bohr
  Atom_vdw( 49) = 1.93/bohr
  Atom_vdw( 50) = 2.17/bohr
  Atom_vdw( 51) = 2.06/bohr
  Atom_vdw( 52) = 2.06/bohr
  Atom_vdw( 53) = 1.98/bohr
  Atom_vdw( 54) = 2.16/bohr
  Atom_vdw( 55) = 3.43/bohr
  Atom_vdw( 56) = 2.68/bohr
  Atom_vdw( 57) = 1.00
  Atom_vdw( 58) = 1.00
  Atom_vdw( 59) = 1.00
  Atom_vdw( 60) = 1.00
  Atom_vdw( 61) = 1.00
  Atom_vdw( 62) = 1.00
  Atom_vdw( 63) = 1.00
  Atom_vdw( 64) = 1.00
  Atom_vdw( 65) = 1.00
  Atom_vdw( 66) = 1.00
  Atom_vdw( 67) = 1.00
  Atom_vdw( 68) = 1.00
  Atom_vdw( 69) = 1.00
  Atom_vdw( 70) = 1.00
  Atom_vdw( 71) = 1.00
  Atom_vdw( 72) = 1.00
  Atom_vdw( 73) = 1.00
  Atom_vdw( 74) = 1.00
  Atom_vdw( 75) = 1.00
  Atom_vdw( 76) = 1.00
  Atom_vdw( 77) = 1.00
  Atom_vdw( 78) = 1.75/bohr
  Atom_vdw( 79) = 1.66/bohr
  Atom_vdw( 80) = 1.55/bohr
  Atom_vdw( 81) = 1.96/bohr
  Atom_vdw( 82) = 2.02/bohr
  Atom_vdw( 83) = 2.07/bohr
  Atom_vdw( 84) = 1.97/bohr
  Atom_vdw( 85) = 2.02/bohr
  Atom_vdw( 86) = 2.20/bohr
  Atom_vdw( 87) = 3.48/bohr
  Atom_vdw( 88) = 2.83/bohr
  Atom_vdw( 89) = 1.00
  Atom_vdw( 90) = 1.00
  Atom_vdw( 91) = 1.00
  Atom_vdw( 92) = 1.86/bohr
  Atom_vdw( 93) = 1.00
  Atom_vdw( 94) = 1.00
  Atom_vdw( 95) = 1.00
  Atom_vdw( 96) = 1.00
  Atom_vdw( 97) = 1.00
  Atom_vdw( 98) = 1.00
  Atom_vdw( 99) = 1.00
  Atom_vdw(100) = 1.00
  Atom_vdw(101) = 1.00
  Atom_vdw(102) = 1.00
  Atom_vdw(103) = 1.00

  do i_species = 1, n_species, 1
    radius_esp_min(i_species) = &
       Atom_vdw(int(species_z(i_species)))*(esp_min)
    radius_esp_max(i_species) = &
       Atom_vdw(int(species_z(i_species)))*(esp_max)
  enddo


  if (allocated(Atom_vdw)) then
     deallocate(Atom_vdw)
  end if

endsubroutine esp_radius

subroutine check_point(rm,lattice_vector,radius_esp_min_new,coord_new,counter_2,&
                           n_atoms, n_periodic, n_species,species,coords_center,&
                       map_to_center_cell_matrix)


!  PURPOSE
!    Check if point will be included
!    o n_grid_batches_esp - total number of batches in the grid
!    o n_my_batches_esp - the number of batches the thread has
!    o batch_sizes - size of each batch
!    o batch_task_list_esp - distribution of the batches over the threads
!    o grid_partition - partition of the points of the grid into batches
!    In the end the array batches is created that stores all the information 
!    threadwise
!  USES

      implicit none
!  ARGUMENTS
      integer,INTENT(IN)	:: rm,n_atoms, n_periodic, n_species
      real*8,INTENT(IN)		:: lattice_vector(3,3)
      real*8,INTENT(IN)		:: radius_esp_min_new(n_species)
      real*8,INTENT(IN)		:: coord_new(3)
      integer,INTENT(OUT)	:: counter_2
      integer,INTENT(IN)	:: species(n_atoms)
      real*8,INTENT(IN)		:: coords_center(3,n_atoms)
      real*8,INTENT(IN)		:: map_to_center_cell_matrix(n_periodic,3)
!  INPUTS
!    
!  OUTPUT
!    
!  SOURCE


      integer :: rn_x,rn_y,rn_z, j_atom
      real*8 :: r_nvec(3),r_lattvec(3),dir_r(3)
      counter_2 = 0
      if (n_periodic > 0) then
             call map_to_center_cell(n_periodic,coord_new(1:3),lattice_vector,&
                  map_to_center_cell_matrix)
      endif
      do j_atom = 1, n_atoms, 1
         if (n_periodic > 0) then
            do rn_x=-rm,rm,1
            do rn_y=-rm,rm,1
            do rn_z=-rm,rm,1
	        r_nvec(1) = rn_x
	        r_nvec(2) = rn_y
	        r_nvec(3) = rn_z
	        r_lattvec(1:3) = matmul(lattice_vector,&
					r_nvec(1:3))
	        dir_r(:)=coord_new(:)-(coords_center( :,j_atom )+&
				r_lattvec(:))
		dist=sqrt(sum(dir_r**2))
		if (dist>=radius_esp_min_new(species(j_atom)))then
			  counter_2=counter_2+1
		endif
            enddo
            enddo
            enddo
         else
	    dist =  sqrt(sum((coord_new(:)-coords_center( :,j_atom ))**2))
	     if (dist>=radius_esp_min_new(species(j_atom)))then
		  counter_2=counter_2+1
	     endif
         endif
      enddo
    end subroutine check_point

  subroutine get_reciprocal_vectors(n_periodic, lattice_vector, &
  &                                 recip_lattice_vector)

    !  PURPOSE
    !
    !   Calculates reciprocal lattice vectors for periodic systems from lattice
    !   vectors
    !
    !     b1 = 2p (a2 x a3 / a1 . (a2 x a3))
    !    
    !     b2 = 2p (a3 x a1 / a1 . (a2 x a3))
    !    
    !     b3 = 2p (a1 x a2 / a1 . (a2 x a3)) 
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_periodic
    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(OUT) :: recip_lattice_vector(3, n_periodic)

    !  INPUTS
    !    o n_periodic -- Number of periodic dimensions
    !    o lattice_vector -- Bravais lattice
    !  OUTPUT
    !    o recip_lattice_vector -- Reciprocal lattice vectors
    !                              2*pi* transpose(lattice_vector^-1)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8 :: inv_bra(n_periodic, 3)
    integer :: rank
    character(*), parameter :: func = 'get_reciprocal_vectors'

    call pseudo_inverse(func, 3, n_periodic, lattice_vector, inv_bra, &
    &                   1d-10, rank)
    recip_lattice_vector = 2*pi * transpose(inv_bra)

  end subroutine get_reciprocal_vectors

  subroutine pseudo_inverse(caller, M, N, A, Ainv, eps, rank)
    !  PURPOSE
    !    Calculate the pseudoinverse A^+ of A with the properties
    !      A A^+ A == A, A^+ A A^+ == A, (A A^+)^* == A A^+, (A^+ A)^* == A A^+
    !    as described in
    !    http://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse
    !    The Pseudoinverse has the property that
    !      x := A^+ b
    !    gives the 'least squares' solution of the problem
    !      A b == x
    !    in the sense that |Ax - b|_2 and |x|_2 are both minimal.
    !  USES

    implicit none
    character(*), intent(IN)            :: caller
    real*8, intent(IN)                  :: A(M,N)
    real*8, intent(OUT)                 :: Ainv(N,M)
    real*8, intent(IN)                  :: eps
    integer, intent(OUT), optional      :: rank

    !  INPUTS
    !    o caller -- calling procedure (for error messages)
    !    o M, N -- matrix dimensions
    !    o A -- M x N matrix
    !    o eps -- threshold for small eigenvalues
    !  OUTPUTS
    !    o Ainv -- right hand side
    !    o rank (optional) -- number of nonzero singular values
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer                             :: M, N, MNmax, MNmin, LDA, LDU, LDVT
    real*8, allocatable                 :: A_tmp(:,:)
    real*8, allocatable                 :: U(:,:), VT(:,:), S(:)
    integer                             :: INFO
    real*8, allocatable                 :: WORK(:)
    real*8                              :: WORK_tmp(1)
    integer                             :: LWORK
    character(*), parameter             :: JOBU = 'A' ! All columns of U
    character(*), parameter             :: JOBVT = 'A' ! All rows of VT
    character(*), parameter             :: func = 'pseudo_inverse'
    character*150                       :: info_str
    integer                             :: i, j
    external                               dgesvd, dgemm

    ! --- prepare

    MNmax = max(M, N)
    MNmin = min(M, N)
    allocate(A_tmp(M, N), U(M, M), VT(N, N), S(MNmin))
    A_tmp = A
    LDA = M
    LDU = M
    LDVT = N

    ! --- work space query

    LWORK = -1
    call dgesvd(JOBU, JOBVT, M, N, A_tmp, LDA, S, U, LDU, VT, LDVT, &
    &           WORK_tmp, LWORK, INFO)
    if (INFO /= 0) then
       write(*, "('For ',A,': DGESVD workspace info =',I5)") &
       & trim(caller), INFO
    end if
    LWORK = nint(WORK_tmp(1))
    allocate(WORK(LWORK))

    ! --- call

    call dgesvd(JOBU, JOBVT, M, N, A_tmp, LDA, S, U, LDU, VT, LDVT, &
    &           WORK, LWORK, INFO)
    if (INFO /= 0) then
       write(*, "('For ',A,': DGESVD workspace info =',I5)") &
       & trim(caller), INFO
    end if

    ! --- invert

    if (present(rank)) rank = 0
    do i = 1, MNmin
       if (S(i) > eps .and. present(rank)) rank = rank + 1
       ! S(i) = S(i) / (S(i)**2 + eps**2)
       if (S(i) > eps) then
          S(i) = 1.d0 / S(i)
       else
          S(i) = 0.d0
       end if
    end do

    ! --- calculate Ainv

    ! A_tmp = ( V.S^+ )^T = S^+ . VT
    do j = 1, N
       do i = 1, MNmin
          A_tmp(i, j) = S(i) * VT(i, j)
       end do
       do i = MNmin+1, M
          A_tmp(i, j) = 0.d0
       end do
    end do

    ! Ainv = A_tmp^T . U^T
    ! Ainv = matmul(transpose(A_tmp), transpose(U))
    call dgemm('T', 'T', N, M, M, 1.d0, A_tmp, M, U, M, 0.d0, Ainv, N)

    ! --- test

    A_tmp = matmul(A, matmul(Ainv, A))
    if (maxval(A_tmp - A) > 2*MNmax*eps) then
       write(info_str, &
       &     "('*** ',A,', for ',A,': error of',ES10.2,'; eps:',ES10.2)") &
       & trim(func), trim(caller), maxval(A_tmp - A), eps
       write(*,"(2X,'Task',I8,': ',A)") myid, trim(info_str)
       ! call aims_stop(info_str)
    end if

  end subroutine pseudo_inverse
  !******
  subroutine get_cell_volume(lattice_vector, cell_volume)

    !  PURPOSE
    !
    !       Calculates the supercell volume in the periodic systems
    !
    !  USES
    
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: lattice_vector(3, 3)
    real*8, intent(OUT) :: cell_volume

    !  INPUTS
    !    o lattice_vector -- Bravais lattice
    !  OUTPUT
    !    o cell_volume -- (Positive) volume of the unit cell.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    
    character(*), parameter :: func = 'get_cell_volume'

    cell_volume = &
       abs(   lattice_vector(1,1) * lattice_vector(2,2) * lattice_vector(3,3)  &
            + lattice_vector(1,2) * lattice_vector(2,3) * lattice_vector(3,1)  &
            + lattice_vector(1,3) * lattice_vector(2,1) * lattice_vector(3,2)  &
            - lattice_vector(1,1) * lattice_vector(2,3) * lattice_vector(3,2)  &
            - lattice_vector(1,2) * lattice_vector(2,1) * lattice_vector(3,3)  &
            - lattice_vector(1,3) * lattice_vector(2,2) * lattice_vector(3,1)  &
          )

    if( cell_volume < 1d-10 ) then
       write(*, "(A,ES10.2)") &
       & 'Error: linear dependent lattice vectors; cell_volume =', cell_volume
    end if

  end subroutine get_cell_volume

subroutine esp_name(species_name,species_z,n_species)
!  PURPOSE
!   Generate the points where the potential will be evaluated. 
!    vdw-Radius of atoms +radius < points < vdw-Radius of atoms +radius2
!  USES
!  ARGUMENTS
!  INPUTS
!   o esp_min -- factor for max radius used for ESP fitting
!   o esp_max -- factor for min radius used for ESP fitting
!  OUTPUT
!   o radius_esp_min -- esp_min*vdw_radius
!   o radius_esp_max -- esp_max*vdw_radius
!   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  character*2,  dimension(n_species) :: species_name
  integer :: species_z(n_species)
  integer :: n_species
  integer :: i_species
  character(2), dimension(:),       allocatable :: Atom_vdw
  if(.not. allocated( Atom_vdw))then
     allocate( Atom_vdw(103),stat=info)
  end if
  ! vdw radius for all atoms, numbers have been taken from wikipedia
  Atom_vdw(  1) = 'H'
  Atom_vdw(  2) = 'He'
  Atom_vdw(  3) = 'Li'
  Atom_vdw(  4) = 'Be'
  Atom_vdw(  5) = 'B'
  Atom_vdw(  6) = 'C'
  Atom_vdw(  7) = 'N'
  Atom_vdw(  8) = 'O'
  Atom_vdw(  9) = 'F'
  Atom_vdw( 10) = 'Ne'
  Atom_vdw( 11) = 'Na'
  Atom_vdw( 12) = 'Mg'
  Atom_vdw( 13) = 'Al'
  Atom_vdw( 14) = 'Si'
  Atom_vdw( 15) = 'P'
  Atom_vdw( 16) = 'S'
  Atom_vdw( 17) = 'Cl'
  Atom_vdw( 18) = 'Ar'
  Atom_vdw( 19) = 'K'
  Atom_vdw( 20) = 'Ca'
  Atom_vdw( 21) = 'Sc'
  Atom_vdw( 22) = 'Ti'
  Atom_vdw( 23) = 'V'
  Atom_vdw( 24) = 'Cr'
  Atom_vdw( 25) = 'Mn'
  Atom_vdw( 26) = 'Fe'
  Atom_vdw( 27) = 'Co'
  Atom_vdw( 28) = 'Ni'
  Atom_vdw( 29) = 'Cu'
  Atom_vdw( 30) = 'Zn'
  Atom_vdw( 31) = 'Ga'
  Atom_vdw( 32) = 'Ge'
  Atom_vdw( 33) = 'As'
  Atom_vdw( 34) = 'Se'
  Atom_vdw( 35) = 'Br'
  Atom_vdw( 36) = 'Kr'
  Atom_vdw( 37) = 'Rb'
  Atom_vdw( 38) = 'Sr'
  Atom_vdw( 39) = 'Y'
  Atom_vdw( 40) = 'Zr'
  Atom_vdw( 41) = 'Nb'
  Atom_vdw( 42) = 'Mo'
  Atom_vdw( 43) = 'Tc'
  Atom_vdw( 44) = 'Ru'
  Atom_vdw( 45) = 'Rh'
  Atom_vdw( 46) = 'Pd'
  Atom_vdw( 47) = 'Ag'
  Atom_vdw( 48) = 'Cd'
  Atom_vdw( 49) = 'In'
  Atom_vdw( 50) = 'Sn'
  Atom_vdw( 51) = 'Sb'
  Atom_vdw( 52) = 'Te'
  Atom_vdw( 53) = 'I'
  Atom_vdw( 54) = 'Xe'
  Atom_vdw( 55) = 'Cs'
  Atom_vdw( 56) = 'Ba'
  Atom_vdw( 57) = ''
  Atom_vdw( 58) = ''
  Atom_vdw( 59) = ''
  Atom_vdw( 60) = ''
  Atom_vdw( 61) = ''
  Atom_vdw( 62) = ''
  Atom_vdw( 63) = ''
  Atom_vdw( 64) = ''
  Atom_vdw( 65) = ''
  Atom_vdw( 66) = ''
  Atom_vdw( 67) = ''
  Atom_vdw( 68) = ''
  Atom_vdw( 69) = ''
  Atom_vdw( 70) = ''
  Atom_vdw( 71) = ''
  Atom_vdw( 72) = ''
  Atom_vdw( 73) = ''
  Atom_vdw( 74) = ''
  Atom_vdw( 75) = ''
  Atom_vdw( 76) = ''
  Atom_vdw( 77) = ''
  Atom_vdw( 78) = ''
  Atom_vdw( 79) = ''
  Atom_vdw( 80) = ''
  Atom_vdw( 81) = ''
  Atom_vdw( 82) = ''
  Atom_vdw( 83) = ''
  Atom_vdw( 84) = ''
  Atom_vdw( 85) = ''
  Atom_vdw( 86) = ''
  Atom_vdw( 87) = ''
  Atom_vdw( 88) = ''
  Atom_vdw( 89) = ''
  Atom_vdw( 90) = ''
  Atom_vdw( 91) = ''
  Atom_vdw( 92) = ''
  Atom_vdw( 93) = ''
  Atom_vdw( 94) = ''
  Atom_vdw( 95) = ''
  Atom_vdw( 96) = ''
  Atom_vdw( 97) = ''
  Atom_vdw( 98) = ''
  Atom_vdw( 99) = ''
  Atom_vdw(100) = ''
  Atom_vdw(101) = ''
  Atom_vdw(102) = ''
  Atom_vdw(103) = ''

  do i_species = 1, n_species, 1
    species_name(i_species) = Atom_vdw(species_z(i_species))
  enddo
  if (allocated(Atom_vdw)) then
     deallocate(Atom_vdw)
  end if

endsubroutine esp_name

subroutine ksummation(coord_current,coords,km,R_c,recip_lattice_vector,alpha,&
                      k_sum_1)
!  PURPOSE
!   Generate the points where the potential will be evaluated. 
!    vdw-Radius of atoms +radius < points < vdw-Radius of atoms +radius2
!  USES
!  ARGUMENTS
!  INPUTS
!   o esp_min -- factor for max radius used for ESP fitting
!   o esp_max -- factor for min radius used for ESP fitting
!  OUTPUT
!   o radius_esp_min -- esp_min*vdw_radius
!   o radius_esp_max -- esp_max*vdw_radius
!   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  real*8,intent(in) :: coord_current(3)
  real*8,intent(in) :: coords(3)
  integer,intent(in) :: km
  real*8,intent(in) :: R_c 
  real*8,intent(in) :: recip_lattice_vector(3,3)
  real*8,intent(in) :: alpha
  real*8,intent(out) :: k_sum_1
  
  real*8 :: dir_k(3),k_lattvec(3)
  real*8 :: di_k
  integer :: nx, ny, nz
  integer :: k_nvec(3)

  dir_k(1)=coord_current(1)-coords(1)
  dir_k(2)=coord_current(2)-coords(2)
  dir_k(3)=coord_current(3)-coords(3)
  di_k=sqrt(sum(dir_k**2))
  k_sum_1 = 0d0
  if (di_k.le.R_c)then
     do nx=-km,km,1
     do ny=-km,km,1
     do nz=-km,km,1
        k_nvec(1) = nx
        k_nvec(2) = ny
        k_nvec(3) = nz
        if (sqrt(sum(k_nvec**2.)).le.km.and.sqrt(sum(k_nvec**2.))&
                      .ne.0)then
           k_lattvec(1:3) = matmul(recip_lattice_vector,&
                                     k_nvec(1:3))
           k_sum_1 =  k_sum_1 + cos(dot_product(k_lattvec,dir_k))*&
                             !(exp(-(sum(k_lattvec**2))/&
                             !(4.*alpha**2))/sum(k_lattvec**2))
                             (exp(-(sum(k_lattvec**2)*pi**2)/&
                             ((alpha)))/sum(k_lattvec**2))
        endif
     enddo
     enddo
     enddo
  endif

endsubroutine ksummation

subroutine rsummation(coord_current,coords,rm,R_c,lattice_vector,alpha,&
                      r_sum_1)
!  PURPOSE
!   Generate the points where the potential will be evaluated. 
!    vdw-Radius of atoms +radius < points < vdw-Radius of atoms +radius2
!  USES
!  ARGUMENTS
!  INPUTS
!   o esp_min -- factor for max radius used for ESP fitting
!   o esp_max -- factor for min radius used for ESP fitting
!  OUTPUT
!   o radius_esp_min -- esp_min*vdw_radius
!   o radius_esp_max -- esp_max*vdw_radius
!   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  real*8,intent(in) :: coord_current(3)
  real*8,intent(in) :: coords(3)
  integer,intent(in) :: rm
  real*8,intent(in) :: R_c 
  real*8,intent(in) :: lattice_vector(3,3)
  real*8,intent(in) :: alpha
  real*8,intent(out) :: r_sum_1
  
  real*8 :: dir_r(3),r_lattvec(3)
  real*8 :: di_r
  integer :: rn_x, rn_y, rn_z
  integer :: r_nvec(3)

  r_sum_1 = 0d0
  do rn_x=-rm,rm,1
  do rn_y=-rm,rm,1
  do rn_z=-rm,rm,1
     r_nvec(1) = rn_x
     r_nvec(2) = rn_y
     r_nvec(3) = rn_z
     r_lattvec(1:3) = matmul(lattice_vector,&
                                   r_nvec(1:3))
     dir_r(1)=coord_current(1)-(coords(1)+r_lattvec(1))
     dir_r(2)=coord_current(2)-(coords(2)+r_lattvec(2))
     dir_r(3)=coord_current(3)-(coords(3)+r_lattvec(3))
     di_r=sqrt(sum(dir_r**2))
     if (di_r.le.R_c)then
        !r_sum_1  = r_sum_1 + erfc(((alpha))*&
        !                     di_r)/di_r
        r_sum_1  = r_sum_1 + erfc(sqrt(alpha)*&
                                         di_r)/di_r
     endif
  enddo
  enddo
  enddo

endsubroutine rsummation

subroutine rsummation_cluster(coord_current,coords,r_sum_1)
!  PURPOSE
!   Generate the points where the potential will be evaluated. 
!    vdw-Radius of atoms +radius < points < vdw-Radius of atoms +radius2
!  USES
!  ARGUMENTS
!  INPUTS
!   o esp_min -- factor for max radius used for ESP fitting
!   o esp_max -- factor for min radius used for ESP fitting
!  OUTPUT
!   o radius_esp_min -- esp_min*vdw_radius
!   o radius_esp_max -- esp_max*vdw_radius
!   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  real*8,intent(in) :: coord_current(3)
  real*8,intent(in) :: coords(3)
  real*8,intent(out) :: r_sum_1
  
  real*8 :: dir_r(3)
  real*8 :: di_r
  dir_r(1)=coord_current(1)-coords(1)
  dir_r(2)=coord_current(2)-coords(2)
  dir_r(3)=coord_current(3)-coords(3)
  di_r=sqrt(sum(dir_r**2))
  r_sum_1  =  1d0/di_r

endsubroutine rsummation_cluster

subroutine read_geo(n_atoms,n_species,n_periodic,coords,lattice_vector,species_name,esp_atom_constraint)

      implicit none

      integer,intent(in) :: n_atoms
      integer,intent(in) :: n_species
      integer,intent(out) :: n_periodic
      real*8, dimension(3,n_atoms),intent(out)  :: coords
      real*8, dimension(3,3),intent(out)  :: lattice_vector
      character*2, dimension(n_species),intent(in) :: species_name
      real*8, dimension(3,n_atoms),intent(out)  :: esp_atom_constraint
      integer :: i_code,i_atom,i_periodic,i_coord,i_species
      logical :: eof
      character*40 :: desc_str
      character*20 :: species_temp
      logical :: found
      logical :: lv_found_before =.false.
      real*8, parameter  :: bohr = 0.5291772085900001D0

      i_periodic = 0
      eof = .False.
      i_atom = 0
      open (8, FILE="geometry.in" )

!  read first line

      read (8,*,iostat=i_code) desc_str

      if (i_code.ne.0) then

         if (myid.eq.0) then
            write(*,*) "Invalid file geometry.in, or file not found."
            write(*,*) "Aborting."
         end if

         stop

      end if

      do while (.not.eof)

        if (desc_str(1:1).eq."#") then

          continue
        else if (desc_str.eq."lattice_vector") then
          backspace(8)

         i_periodic = i_periodic+1
         lv_found_before =.true.

          read(8,*) desc_str, &
            (lattice_vector(i_coord,i_periodic), i_coord=1,3,1)

          do i_coord = 1,3,1
            lattice_vector(i_coord,i_periodic) = &
            lattice_vector(i_coord,i_periodic) / bohr
          enddo
        else if ((desc_str.eq."atom") .or. (desc_str.eq."atom_frac")) then
          backspace(8)

          i_atom = i_atom+1

          lv_found_before =.false.

          read(8,*) desc_str, (coords(i_coord,i_atom), i_coord=1,3,1), &
            species_temp

!  transform coordinates into atomic units (bohr)

          do i_coord = 1,3,1
            coords(i_coord,i_atom) = coords(i_coord,i_atom)/bohr
          enddo

          found = .false.

          do i_species=1, n_species, 1
            if (species_temp.eq.species_name(i_species)) then
              found = .true.
            end if
          enddo

          if (.not.found) then
             if (myid.eq.0) then
                write(*,'(1X,A,A,A,A,A)') &
                     "* Species ", species_temp, ", listed in ", &
                     "geometry.in, is not described in input file ", &
                     "control.in."
                write(*,*) "* List of ", n_species, " species: "
             end if
             do i_species = 1,n_species,1
                if (myid.eq.0) then
                   write(*,*) "* ",i_species,":",species_name(i_species)
                end if
            enddo
            if (myid.eq.0) then
               write(*,*) "* Aborting."
            end if

            stop
          end if

        else if (desc_str.eq."esp_constraint") then
          !contraints for ESP charge fitting
          if (esp_constraint.eq.1)then
            backspace(8)
            read(8,*) desc_str, &
            esp_atom_constraint(1,i_atom), esp_atom_constraint(2,i_atom),&
            esp_atom_constraint(3,i_atom)
          elseif (esp_constraint.eq.2)then
            backspace(8)
            read(8,*) desc_str, &
            esp_atom_constraint(1,i_atom), esp_atom_constraint(2,i_atom)
          else
            if (myid.eq.0) then
            write(*,*) &
            "No esp constrains (method 1/2) requested - ignoring constraint."
            endif
          endif
        end if

        read (8,*,iostat=i_code) desc_str

        if (i_code.ne.0) then
          eof = .true.
        end if

      end do
      n_periodic = i_periodic
endsubroutine read_geo


subroutine esp_output_zcube(filename, esp_min, esp_max, esp_cube_output,esp_cube_output_org,n_z)
      implicit none

      CHARACTER(*),INTENT(IN):: filename
      real*8,INTENT(IN) ::  esp_min
      real*8,INTENT(IN) ::  esp_max  
      integer,INTENT(IN) :: n_z
      real*8, dimension(n_z) ,INTENT(IN) :: esp_cube_output
      real*8, dimension(n_z) ,INTENT(IN) :: esp_cube_output_org
      CHARACTER(len=len(trim(filename))+12) :: outname
      CHARACTER(len=4) :: value1, value2
      write(value1,'(1F3.1)') esp_min
      write(value2,'(1F3.1)') esp_max
      outname = trim(filename)//'_'//trim(value1)//'_'//trim(value2)//'.dat'
      open (10, FILE=outname )
	do i=1, n_z
	    WRITE (10,fmt = '(2ES17.4)') esp_cube_output(i), esp_cube_output_org(i)
	end do
     close(unit=10)
endsubroutine esp_output_zcube

!---------------------------------------------------------------------------
!****s* esp_cahrges/esp_output_cube
!  NAME
!    esp_output_cube
!  SYNOPSIS
subroutine esp_output_cubefile(cube_output,n_x,n_y,n_z,&
                               species_z,n_atoms,n_periodic,n_species,&
                               species,coords,lattice_vector,&
                               local_offset,radius_esp_max_new,i_esp)
!  PURPOSE
!    Output the potential into a cube file
!
!  USES
  implicit none
!  use plot

  real*8,     dimension(n_x*n_y*n_z) :: cube_output
  integer  :: n_x,n_y,n_z,n_atoms,n_periodic,n_species,i_esp
  integer,  dimension(n_species) :: species_z
  integer,  dimension(n_atoms) :: species
  real*8,   dimension(3,3) :: lattice_vector
  real*8,   dimension(3,n_atoms) :: coords
  real*8,  dimension(n_species) :: radius_esp_max_new
  real*8,   dimension(3) ::   local_offset
!  ARGUMENTS
!  INPUTS
!   o cube_output -- Array with potential data
!   o n_x, n_y, n_z -- number of points in x, y, z
!  OUTPUT
!   o cube file "potential_esp.cube"
!   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  real*8  coords_temp(3,n_atoms)
  real*8, dimension(3,n_atoms) :: periodic_unit_cell_translations_cart
  integer :: descriptor
  character(LEN=15):: file_format
  character(LEN=29):: filename_esp
  integer :: blocksize, i_point, i_atom,i_coord, n_points
  real*8, dimension(3) :: d_r
  real*8, dimension(3) :: cube_coord_min
  real*8, dimension(3) :: cube_coord_max
  CHARACTER(len=1) :: valuek1
  CHARACTER(len=2) :: valuek2
  CHARACTER(len=3) :: valuek3
  CHARACTER(len=4) :: valuek4
  CHARACTER(len=5) :: valuek5
  real*8  :: lattice_vector_new(3,3)
  real*8  :: off
  integer :: i_x,i_y,i_z

     if (myid.eq.0) then
       if (i_esp<=9) then
	  write(valuek1,'(I1)') i_esp
	  filename_esp = 'potential_esp_fit_'//trim(valuek1)//'.cube'
       elseif (i_esp<=99) then
	  write(valuek2,'(I2)') i_esp
	  filename_esp = 'potential_esp_fit_'//trim(valuek2)//'.cube'
       elseif (i_esp<=999) then
	  write(valuek3,'(I3)') i_esp
	  filename_esp = 'potential_esp_fit_'//trim(valuek3)//'.cube'
       elseif (i_esp<=9999) then
	  write(valuek4,'(I4)') i_esp
	  filename_esp = 'potential_esp_fit_'//trim(valuek4)//'.cube'
       elseif (i_esp<=99999) then
	  write(valuek5,'(I5)') i_esp
	  filename_esp = 'potential_esp_fit_'//trim(valuek5)//'.cube'
       endif
       !Open file
       descriptor = 11+i_esp
       open(unit= descriptor, file=filename_esp, ACTION='WRITE')
       n_points = n_x*n_y*n_z
       !Write the header
       coords_temp = coords


       write (descriptor,*) "CUBE FILE written by FHI-AIMS"
       write (descriptor,*) "*****************************"

       write (descriptor,fmt='(1X,I4,3F12.6)') &
              n_atoms, (local_offset(i_coord), i_coord=1,3,1)
       write(descriptor,fmt='(1X,I4,3F12.6)') &
           n_x, lattice_vector(1,1:3)*n_x
       write(descriptor,fmt='(1X,I4,3F12.6)') &
           n_y, lattice_vector(2,1:3)*n_y
       write(descriptor,fmt='(1X,I4,3F12.6)') &
           n_z, lattice_vector(3,1:3)*n_z

       do i_atom = 1,n_atoms,1
         write (descriptor,fmt='(1X,I4,4F12.6)') &
         int(species_z(species(i_atom))), &
         (coords_temp(i_coord,i_atom),i_coord=1,3,1)
       enddo
       write(file_format,'(A)') '(1E13.5, $)'
       blocksize = 1 
       i_point = 0
	do i_x = 0, n_x-1, 1
	do i_y = 0, n_y-1, 1
	do i_z = 0, n_z-1, 1
         i_point = i_point +1
                      write (descriptor,fmt=file_format) &
                           cube_output(i_point)
                       if(mod(blocksize,6).eq.0) then
                          write(descriptor,*) ""
                          blocksize=0
                      elseif(mod(i_point,n_z).eq.0) then
                           write(descriptor,*) ""
                           blocksize=0
                      endif
                      blocksize= blocksize+1
         enddo
         enddo
         enddo
      ! Add a carriage return at the end
      write(descriptor,*) ""
      !close the file
      close(descriptor)
  endif !myid = 0
endsubroutine esp_output_cubefile

subroutine map_to_center_cell(n_periodic,coord,lattice_vector,&
                              map_to_center_cell_matrix)

!  PURPOSE
!   Maps the given coordinate point to the "center cell" (Wigner Seitz)
!   around origin (periodic systems).
!   
!   ****************************************************************************************** 
!   **** If you make ANY changes here, please update "map_to_first_octant.f90" accordingly *** 
!   ****************************************************************************************** 
!
!  USES

  implicit none

!  ARGUMENTS
  integer:: n_periodic
  real*8:: coord(3)
  real*8:: lattice_vector(3,3)
  real*8:: map_to_center_cell_matrix(n_periodic,3)
!  INPUTS
!   o  coord -- coordinate point before mapping
!  OUTPUT
!   o  coord -- coordinate point after mapping
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

  real*8:: length(3)

  if (n_periodic .eq. 0) return

! Atoms right on top of the cell limit will be moved to one
! or the other side depending on a rounding error within
! fortran. In order to have a more controlled shifting of
! the atoms, a very small fraction is substracted from the
! fractional coordinate and before converting it back
! into real-space coordinates the small fraction is added
! again, to ensure that atoms are not moved.

  length = matmul(map_to_center_cell_matrix, coord)
  length = length-1D-8
  length = length - nint(length) + 1D-8
  coord = matmul(lattice_vector, length)

end subroutine map_to_center_cell

 subroutine get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
  &                                        map_to_center_cell_matrix)

    !  PURPOSE
    !
    !  Updates map to center cell matrix. 
    !  Mind that :
    !     (1) map_to_center_cell_matrix == inverse of lattice vector matrix
    !         --> somewhat odd naming
    !     (2) in the case that the unit cell changes, the center cell matrix
    !         changes as well
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN)  :: n_periodic
    real*8,  intent(IN)  :: lattice_vector(3, n_periodic)
    real*8,  intent(OUT) :: map_to_center_cell_matrix(n_periodic, 3)
    
    !  INPUTS
    !    o n_periodic -- Number of periodic dimensions
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    o map_to_center_cell_matrix -- inverse of lattice_vector
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank
    character(*), parameter :: func = 'get_map_to_center_cell_matrix'

    call pseudo_inverse(func, 3, n_periodic, lattice_vector, &
    &                   map_to_center_cell_matrix, 1d-10, rank)
    if(rank /= n_periodic) then
       stop
       write(*,*) 'ERROR: lattice_vector is singular!'
    end if

  end subroutine get_map_to_center_cell_matrix
End Program Fit_esp_charges
        
        