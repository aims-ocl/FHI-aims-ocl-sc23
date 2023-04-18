! CHECK ON "MPI implementation" and parallel DOS

! arguments: hessian, geometry.in.central, control.in.phonon
! should be enough to get all the necessary data, hessian transforms etc

program phonon_hessian_diagonalization
  use MPI_ROUTINES
  use arch_specific, only : arch_erf
  implicit none
  character*100 :: control_name, geometry_name, hessian_name, arg, keyword, band_file_name, dos_file_name, work_dir
  character*5   :: frequency_unit
  integer       :: n_bands, n_atoms, n_species, n_supercells, ioerr, i_point, i_freq, dos_density
  logical       :: get_dos, have_soft_mode, get_cv
  real*8        :: dos_fstart, dos_fend, dos_broad, buf, cell_volume, k_weight, freq, deltaT, ZPE(1), omegaT
  real*8        :: freq1, freq2, deltaf, dos_factor, frequency_factor, cv_start, cv_end, exp_term, omega, T
  integer       :: output_unit, k_index, ik1, ik2, ik3, i_dos, iT
  integer       :: dos_fpoints, sc_max, sc_max_total, cv_points, cv_density, cv_points_start
  integer       :: index, i_coord, is1, is2, is3, i_atom2, i_species, i_band, i_atom, i_lattice, i_coord2
  integer, dimension(3  ) :: supercell, sc_start, sc_end
  real*8,  dimension(3,3) :: lattice_vector, recip_lattice
  real*8,  dimension(3  ) :: k_point_rel, k_point_abs
  character*20, dimension(:            ), allocatable :: species_names, species
  real*8,       dimension(:            ), allocatable :: species_masses, masses
  real*8,       dimension(:,:          ), allocatable :: band_info
  integer,      dimension(:            ), allocatable :: band_points
  real*8,       dimension(:,:          ), allocatable :: dynamic_matrix
  real*8,       dimension(:,:          ), allocatable :: posi_cart, posi_frac
  real*8,       dimension(:,:,:,:,:,:,:), allocatable :: hessian
  real*8,       dimension(:            ), allocatable :: hessian_line
  real*8,       dimension(:            ), allocatable :: eigenvalues
  integer,      dimension(:,:          ), allocatable :: sc_index
  real*8,       dimension(:,:          ), allocatable :: sc_prefactor
  real*8,       dimension(:            ), allocatable :: dos
  real*8,       dimension(:            ), allocatable :: cv
  real*8,       dimension(:            ), allocatable :: free_energy
  real*8,       dimension(:            ), allocatable :: internal_energy
  real*8,       dimension(:            ), allocatable :: ev_use
  ! physical constants - these are in SI units, 
  !       unlike the values in the module constants.f90 in the main FHI-aims source directory
  real*8 :: pi               = 3.14159265358979323846d0
  real*8 :: const_u          = 1.660538782d-27
  real*8 :: const_eV         = 1.602176487d-19
  real*8 :: const_c          = 299792458d0
  real*8 :: const_Angstr     = 1d-10
  real*8 :: const_N_avogadro = 6.02214179d23
  real*8 :: const_h          = 6.62606896d-34
  real*8 :: boltzmann_kB     = 1.3806504d-23
  real*8 :: one_over_sqrt2   = 0.70710678d0
  real*8 :: hessian_factor

  call initialize_mpi ( )

  if (myid.eq.0) then
     write(use_unit,*) "=========================================================="
     write(use_unit,*) 
     write(use_unit,*) "Entering Hessian diagonalization program for FHI-aims"
     write(use_unit,*)
     write(use_unit,*) "                 Version XXXXXX"
     write(use_unit,*) 
     write(use_unit,*) "=========================================================="
  end if

  ! obtain arguments from command line
  call getarg(1,hessian_name)
  call getarg(2,geometry_name)
  call getarg(3,control_name)
  call getarg(4,work_dir)

  ! default for frequency units in absence of user specs
  frequency_unit   = "cm^-1"
  frequency_factor = 1d0/(200d0*pi*const_c) ! set default as inverse cm - as done in the vibrations script 

  ! parse control files to determine necessary dimensions
  if (myid.eq.0) write(use_unit,*)
  if (myid.eq.0) write(use_unit,*) "Parsing control.in ... "
  open(unit = 20, file = control_name, status = "old")
  ioerr        = 0
  n_bands      = 0
  n_species    = 0
  n_supercells = 0
  get_dos      = .false.
  get_cv       = .false.
  cv_density   = 0
  dos_density  = 0
  do while (ioerr.eq.0)
     read(unit = 20, fmt = *, iostat = ioerr) keyword
     if (ioerr.eq.0) then
        if (keyword.eq."phonon") then
           backspace(20) 
           read(20,*) keyword, keyword
           if (keyword.eq."band") then
              n_bands = n_bands + 1 
           else if ((keyword.eq."dos").or.(keyword.eq."DOS")) then
              get_dos = .true.
              backspace(20)
              read(20,*) keyword, keyword, dos_fstart, dos_fend, dos_fpoints, dos_broad, dos_density
           else if (keyword.eq."supercell") then
              backspace(20)
              read(20,*) keyword, keyword, supercell(1), supercell(2), supercell(3)
              n_supercells = supercell(1)*supercell(2)*supercell(3)
              if (myid.eq.0) write(use_unit,'(2X,A,I2,A,I2,A,I2)') '| Found supercell ',supercell(1),' x ',supercell(2),' x ',supercell(3)
              if (myid.eq.0) write(use_unit,'(2X,A,I2,A)') '|   for a total number of ',n_supercells,' unit cells.'
           else if (keyword.eq."frequency_unit") then
              backspace(20)
              read(20,*) keyword, keyword, frequency_unit 
              if (trim(frequency_unit).eq.'cm^-1') then
                 frequency_factor = 1d0/(200d0*pi*const_c)
              else if (trim(frequency_unit).eq.'THz') then                 
                 frequency_factor = 1d0/(2d12*pi)
              else
                 if (myid.eq.0) then
                    write(use_unit,*) "* WARNING: Frequency unit ", trim(frequency_unit), " unknown."
                    write(use_unit,*) "*          Possible values are THz and cm^-1. Please correct."
                    write(use_unit,*) "* Aborting."
                 end if
                 stop
              end if
              if (myid.eq.0) write(use_unit,'(2X,2A)') "| Output using frequency unit : ", trim(frequency_unit)
           else if ((keyword.eq."cv").or.(keyword.eq."free_energy")) then
              backspace(20)
              read(20,*) keyword, keyword, cv_start, cv_end, cv_points, cv_density
              get_cv  = .true.
           end if
        else if (keyword.eq."species") then
           n_species = n_species + 1
        end if
     end if
  end do
  close(unit = 20)

  if (n_supercells.eq.0) then
     if (myid.eq.0) write(use_unit,*) "* WARNING: Number of supercells appears to be zero!"
     if (myid.eq.0) write(use_unit,*) "*          This can not be correct. Please check. Aborting."
     stop
  end if
  if (myid.eq.0) write(use_unit,*)
  if (myid.eq.0) write(use_unit,*) "Parsing geometry.in .... "
  open(unit = 20, file = geometry_name, status = "old")
  ioerr   = 0
  n_atoms = 0
  do while (ioerr.eq.0) 
     read(unit=20,fmt=*,iostat=ioerr) keyword
     if (ioerr.eq.0) then
        if (keyword.eq."atom") n_atoms = n_atoms + 1
     end if
  end do
  close(unit = 20)

  ! allocate data
  if (myid.eq.0) write(use_unit,*)
  if (myid.eq.0) write(use_unit,*) "Allocating data ... "
  allocate(posi_cart(3,n_atoms))  ! cartesian of atoms in the primitive cell
  allocate(posi_frac(3,n_atoms))  ! fractional coord of atoms in primitive cell
  allocate(species_names(n_species))
  allocate(species(n_atoms))
  allocate(species_masses(n_species))
  allocate(masses(n_atoms))
  allocate(hessian(3,n_atoms,3,n_atoms,supercell(1),supercell(2),supercell(3)))
  allocate(dynamic_matrix(3*n_atoms,3*n_atoms))
  allocate(band_info(6,n_bands))
  allocate(band_points(n_bands))
  allocate(hessian_line(supercell(1)*supercell(2)*supercell(3)*n_atoms*3))
  allocate(eigenvalues(3*n_atoms))
  if (myid.eq.0) write(use_unit,*)
  if (myid.eq.0) write(use_unit,*) "Reading control file ... "
  open(unit = 20, file = control_name, status = "old")
  ioerr     = 0
  i_species = 0
  i_band    = 0
  do while (ioerr.eq.0)
     read(unit = 20, fmt = *, iostat = ioerr) keyword
     if (ioerr.eq.0) then
        if (keyword.eq.'phonon') then
           backspace(20)
           read(20,*) keyword, keyword
           if (keyword.eq.'band') then
              backspace(20)
              i_band = i_band + 1 
              read(20,*) keyword, keyword, (band_info(i_coord,i_band), i_coord = 1, 6), band_points(i_band)
              if (myid.eq.0) write(use_unit,'(2X,A)')        "| found phonon band."
              if (myid.eq.0) write(use_unit,'(2X,A,3F10.6)') "|     starting point = ",band_info(1,i_band),band_info(2,i_band),band_info(3,i_band)
              if (myid.eq.0) write(use_unit,'(2X,A,3F10.6)') "|     ending   point = ",band_info(4,i_band),band_info(5,i_band),band_info(6,i_band)
              if (myid.eq.0) write(use_unit,'(2X,A,I4)')     "|   number of points = ", band_points(i_band)
           else if ((keyword.eq.'dos').or.(keyword.eq.'DOS')) then
              backspace(20)
              read(20,*) keyword, keyword, dos_fstart, dos_fend, dos_fpoints, dos_broad, dos_density
              if (myid.eq.0) then
                 write(use_unit,'(2X,A)') "| found phonon DOS specifications : "
                 write(use_unit,'(2X,A,F14.6,2A)') '|     fstart        = ', dos_fstart, ' ',trim(frequency_unit)
                 write(use_unit,'(2X,A,F14.6,2A)') '|     fend          = ', dos_fend,   ' ',trim(frequency_unit)
                 write(use_unit,'(2X,A,I8)')       '|     Npoints       = ', dos_fpoints
                 write(use_unit,'(2X,A,F14.6,2A)') '|     broadening    = ', dos_broad,  ' ',trim(frequency_unit)
                 write(use_unit,'(2X,A,I8)')       '|     point density = ', dos_density
              end if              
           else if (keyword.eq."cv") then
              backspace(20)
              read(20,*) keyword, keyword, cv_start, cv_end, cv_points, cv_density
              if (myid.eq.0) then 
                 write(use_unit,'(2X,A)'      ) '| found phonon specific heat (c_v) specifications : '
                 write(use_unit,'(2X,A,F14.6)') '|     Tstart        = ', cv_start
                 write(use_unit,'(2X,A,F14.6)') '|     Tend          = ', cv_end
                 write(use_unit,'(2X,A,I8)'   ) '|     Npoints       = ', cv_points
                 write(use_unit,'(2X,A,I8)'   ) '|     point density = ', cv_density
              end if
           end if
        else if (keyword.eq.'species') then 
           backspace(20)
           i_species = i_species + 1 
           read(20,*) keyword, species_names(i_species)
           if (myid.eq.0) write(use_unit,'(2X,2A)') '| Found species ', trim(species_names(i_species))
        else if (keyword.eq.'mass') then
           backspace(20)
           read(20,*) keyword, species_masses(i_species)
           if (myid.eq.0) write(use_unit,'(2X,A,F10.5,A)') '|    mass = ', species_masses(i_species), ' amu'
        end if
     end if
  end do
  close(unit = 20)

  if (get_cv.and.get_dos) then
     if (cv_density.ne.dos_density) then
        cv_density = max(cv_density,dos_density)
        if (myid.eq.0) then
           write(use_unit,'(2X,A)')      "* WARNING: number of q-points does not match between DOS and cv calculation! "
           write(use_unit,'(2X,A,I5,A)') "*          using the higher value of ",cv_density," for integration. "
        end if
     end if
  end if
  dos_density = max(dos_density,cv_density)   ! use one single point density for loops, whatever it happens to be in the end. 

  ! read geometry.in 
  if (myid.eq.0) write(use_unit,*)
  if (myid.eq.0) write(use_unit,*) "Reading geometry file ... "
  open(unit = 20, file = geometry_name, status = "old")
  ioerr  = 0
  i_atom = 0
  i_lattice = 0
  do while (ioerr.eq.0)
     read(unit = 20, fmt = *, iostat = ioerr) keyword
     if (ioerr.eq.0) then
        if (keyword.eq.'atom') then
           backspace(20)
           i_atom = i_atom + 1 
       !    read(20,*) keyword, buf, buf, buf, species(i_atom)
         read(20,*) keyword, posi_cart(1,i_atom), posi_cart(2,i_atom),posi_cart(3,i_atom),species(i_atom)
           if (myid.eq.0) write(use_unit,'(2X,2A)') '| Found atom ', trim(species(i_atom))
        else if (keyword.eq.'lattice_vector') then
           backspace(20)
           i_lattice = i_lattice + 1
           read(20,*) keyword, (lattice_vector(i_coord, i_lattice), i_coord = 1, 3)
           if (myid.eq.0) write(use_unit,'(2X,A,I2,A,3F10.6)') '| Found lattice vector ',i_lattice," : ",(lattice_vector(i_coord, i_lattice), i_coord = 1, 3)
        end if
     end if
  end do
  close(20)

  ! transform cartesian atomic coordinates in geometry.in to fractional
  call cart2frac(lattice_vector,posi_cart,posi_frac,n_atoms)

  ! read hessian file
  hessian(:,:,:,:,:,:,:) = 0d0
  if (myid.eq.0) write(use_unit,*)
  if (myid.eq.0) write(use_unit,*) "Reading Hessian input file ... "
  open(unit = 20, file = hessian_name, status = "old")
  do i_atom = 1, n_atoms
     do i_coord = 1, 3
        ! read line into one single long buffer
        read(20,*) (hessian_line(index),index = 1, supercell(1)*supercell(2)*supercell(3)*n_atoms*3)
        index = 1
        ! then disperse it into a more useful fashion ... 
        do is1 = 1, supercell(1)
           do is2 = 1, supercell(2)
              do is3 = 1, supercell(3)
                 do i_atom2 = 1, n_atoms
                    do i_coord2 = 1, 3
                       hessian(i_coord,i_atom,i_coord2,i_atom2,is1,is2,is3) = hessian_line(index)
                       index = index + 1 
                    end do
                 end do
              end do
           end do
        end do        
     end do
  end do
  close(unit = 20)

  ! assign masses to atoms via species
  masses = 0d0
  do i_atom = 1, n_atoms
     do i_atom2 = 1, n_species
        if (species_names(i_atom2).eq.species(i_atom)) masses(i_atom) = species_masses(i_atom2)
     end do
  end do

  ! start core calculations: 
  ! count out  the supercells
  sc_start(:)  = -supercell(:)/2
  sc_end  (:)  =  supercell(:)/2
 ! check even or odd of supercell(:)
  do i_coord = 1, 3
   if (mod(supercell(i_coord),2) .eq. 0) then
     sc_end(i_coord) = sc_end(i_coord) - 1
   end if
  end do
  sc_max       = 0
  sc_max_total = 0
  do i_coord = 1, 3
     sc_max       = max(sc_end(i_coord), sc_max)
     sc_max_total = max(supercell(i_coord),sc_max_total)
  end do
     sc_max = sc_max + 1 

  ! calculate supercell index for easier calculation of dynamic matrix, also get prefactor [i.e. 1/(number of times to count hessian)]
  allocate(sc_index(3,-sc_max:sc_max))
  allocate(sc_prefactor(3,sc_max_total))
  sc_prefactor(:,:) = 0d0
  sc_index(:,:)     = 0
  do i_coord = 1, 3
     do is1 = sc_start(i_coord), sc_end(i_coord)
       if (is1.lt.0)  then
          sc_index(i_coord,is1) = supercell(i_coord) + 1 + is1
       else
          sc_index(i_coord,is1) = is1 + 1
       end if
       sc_prefactor(i_coord, sc_index(i_coord,is1)) = sc_prefactor(i_coord, sc_index(i_coord,is1)) + 1d0 
     end do 
  end do
  
  ! invert prefactor if greater than zero
  do i_coord = 1, 3
     do is1 = 1, sc_max_total
        if (sc_prefactor(i_coord,is1).gt.0d0) sc_prefactor(i_coord,is1) = 1d0/sc_prefactor(i_coord,is1)
     end do
  end do
   
    
  ! insert mass vector into hessian matrix, need 1/sqrt(m_i m_j) in each element H_ij to have a symmetric soln to Newton's equations
  hessian_factor = const_eV/(const_u*const_Angstr*const_Angstr)  ! factor to bring hessian into SI units ... 
  do i_atom = 1, n_atoms
     do i_atom2 = 1, n_atoms
        hessian(:,i_atom,:,i_atom2,:,:,:) = hessian(:,i_atom,:,i_atom2,:,:,:)*hessian_factor/sqrt(masses(i_atom)*masses(i_atom2))
     end do
  end do
  
  ! calculate cell volume (unit cell!) and reciprocal lattice vectors
  cell_volume =   lattice_vector(1,1)*lattice_vector(2,2)*lattice_vector(3,3) &
                + lattice_vector(1,2)*lattice_vector(2,3)*lattice_vector(3,1) &
                + lattice_vector(1,3)*lattice_vector(2,1)*lattice_vector(3,2) &
                - lattice_vector(1,1)*lattice_vector(2,3)*lattice_vector(3,2) &
                - lattice_vector(1,2)*lattice_vector(2,1)*lattice_vector(3,3) &
                - lattice_vector(1,3)*lattice_vector(2,2)*lattice_vector(3,1) 
  recip_lattice(1,1) = 2 * pi * (lattice_vector(2,2) * lattice_vector(3,3) - lattice_vector(3,2) * lattice_vector(2,3)) / cell_volume
  recip_lattice(2,1) = 2 * pi * (lattice_vector(3,2) * lattice_vector(1,3) - lattice_vector(1,2) * lattice_vector(3,3)) / cell_volume
  recip_lattice(3,1) = 2 * pi * (lattice_vector(1,2) * lattice_vector(2,3) - lattice_vector(2,2) * lattice_vector(1,3)) / cell_volume
  recip_lattice(1,2) = 2 * pi * (lattice_vector(2,3) * lattice_vector(3,1) - lattice_vector(3,3) * lattice_vector(2,1)) / cell_volume
  recip_lattice(2,2) = 2 * pi * (lattice_vector(3,3) * lattice_vector(1,1) - lattice_vector(1,3) * lattice_vector(3,1)) / cell_volume
  recip_lattice(3,2) = 2 * pi * (lattice_vector(1,3) * lattice_vector(2,1) - lattice_vector(2,3) * lattice_vector(1,1)) / cell_volume
  recip_lattice(1,3) = 2 * pi * (lattice_vector(2,1) * lattice_vector(3,2) - lattice_vector(3,1) * lattice_vector(2,2)) / cell_volume
  recip_lattice(2,3) = 2 * pi * (lattice_vector(3,1) * lattice_vector(1,2) - lattice_vector(1,1) * lattice_vector(3,2)) / cell_volume
  recip_lattice(3,3) = 2 * pi * (lattice_vector(1,1) * lattice_vector(2,2) - lattice_vector(2,1) * lattice_vector(1,2)) / cell_volume

  if (myid.eq.0) then
     write(use_unit,*)
     write(use_unit,'(2X,A,F10.6)') 'Unit cell volume : ',cell_volume
     write(use_unit,'(2X,A)') 'Reciprocal lattice : '
     do i_coord = 1, 3
        write(use_unit,'(2X,A,3F10.6)')' | ',(recip_lattice(i_coord2,i_coord),i_coord2 = 1, 3)
     end do
  end if

  ! go through the various bands & calculate eigenvalues - this can be parallelized over bands ... (if so, watch for the unit number. better not do this)
  if (n_bands.gt.0) then
     
     if (myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*) "Computing phonon band structure ... "
     end if

     ! hold all operations until everybody has caught up to this point, messed output otherwise
     call wait_for_all_tasks()

     do i_band = 1, n_bands
        if ((mod(i_band,n_tasks).eq.myid).and.(myid.le.n_bands)) then
           if(i_band < 10)then
              write(band_file_name, '(A,A15,I1,A4)') trim(work_dir),'/phonon_band100', i_band,'.out'
           else if(i_band < 100)then
              write(band_file_name, '(A,A14,I2,A4)') trim(work_dir),'/phonon_band10', i_band,'.out'
           else if(i_band < 1000)then
              write(band_file_name, '(A,A13,I3,A4)') trim(work_dir),'/phonon_band1', i_band,'.out'
           end if
           write(use_unit,'(2X,A,I4,2A)') " | Working on band ", i_band, ", writing to file ", trim(band_file_name)
           output_unit = 20 + myid
           open(unit = output_unit, file = band_file_name, status = 'unknown')
           have_soft_mode = .false.
           do i_point = 1, band_points(i_band)
              
              ! calculate relative k-point
              k_point_rel(1) = dble(i_point-1)*(band_info(4,i_band)-band_info(1,i_band))/dble(band_points(i_band)-1)+band_info(1,i_band)
              k_point_rel(2) = dble(i_point-1)*(band_info(5,i_band)-band_info(2,i_band))/dble(band_points(i_band)-1)+band_info(2,i_band)
              k_point_rel(3) = dble(i_point-1)*(band_info(6,i_band)-band_info(3,i_band))/dble(band_points(i_band)-1)+band_info(3,i_band)

              ! calculate absolute k-point
              k_point_abs(:) = 0
              do i_coord = 1, 3
                 do i_coord2 = 1, 3
                    k_point_abs(i_coord) = k_point_abs(i_coord) + recip_lattice(i_coord,i_coord2)*k_point_rel(i_coord2)
                 end do
              end do

              ! get dynamic matrix and diagonalize while at it
              call diagonalize_dynamic_matrix(k_point_abs, lattice_vector, hessian, supercell, n_atoms, sc_max, sc_max_total, &
                   sc_start, sc_end, sc_index, sc_prefactor, eigenvalues, posi_frac)
              
              ! print eigenvalues to file              
              write(output_unit,'(I6,3F14.6,$)') i_point, (k_point_abs(i_coord), i_coord = 1, 3)
              do i_coord = 1, 3*n_atoms
                 write(output_unit,'(f10.3,E14.4E3,$)') sign(1.0d0,eigenvalues(i_coord)),eigenvalues(i_coord)*frequency_factor
              end do
              write(output_unit,*)
              
              ! check on soft modes ... 
              have_soft_mode = (eigenvalues(1).lt.0d0)           
           end do
           if (have_soft_mode) then
              write(use_unit,'(2X,A,I4)') " | * WARNING: suspecting soft mode in band ", i_band
           end if
           close(unit = output_unit)        
        end if
     end do
  end if

  call wait_for_all_tasks()

  ! calculate DOS, preferrably in parallel, using same diagonalization routine
  if (get_dos.or.get_cv) then

     if (get_dos) then
        allocate(dos(dos_fpoints))
        allocate(ev_use(3*n_atoms))
        deltaf = (dos_fend-dos_fstart)/dble(dos_fpoints-1)
        ! normalization inside error functions for DOS calculation
        dos_factor = one_over_sqrt2/dos_broad
     end if
     if (get_cv) then
        allocate(cv(cv_points))
        allocate(free_energy(cv_points))        
        allocate(internal_energy(cv_points))
        deltaT = (cv_end-cv_start)/dble(cv_points-1)
        ! determine the first temperature above zero, the method gives cv=0 for T=0 and is not defined for T < 0 ...
        ! by determining the starting index of the temperature range T > 0, one can get rid of an if statement inside a loop later,
        ! therefore providing vectorizable code.
        cv_points_start = 1
        T = cv_start
        do while (T.le.0)
           cv_points_start = cv_points_start + 1 
           T               = cv_start + dble(cv_points_start-1)*deltaT
        end do
     end if
     if (myid.eq.0) then
        write(use_unit,*)
        if (get_dos.and.(.not.get_cv)) then
           write(use_unit,*) 'Computing phonon density of states'
        else if (get_dos.and.get_cv) then
           write(use_unit,*) 'Computing phonon density of states and free energy components'
        else
           write(use_unit,*) 'Computing phonon free energy components'
        end if
     end if

     ! k-point weight 
     k_weight = 1d0/dble(dos_density**3)

     k_index = 0
     if (get_dos) dos(:) = 0d0
     if (get_cv ) then
        cv             (:) = 0d0
        free_energy    (:) = 0d0
        internal_energy(:) = 0d0
        ZPE                = 0d0
     end if

     ! loop through all k-points
     do ik1 = 1, dos_density
        k_point_rel(1) = (dble(ik1)+5d-1)/dble(dos_density-1)
        do ik2 = 1, dos_density
           k_point_rel(2) = (dble(ik2)+5d-1)/dble(dos_density-1)           
           do ik3 = 1, dos_density
              k_point_rel(3) = (dble(ik3)+5d-1)/dble(dos_density-1)
              k_index = k_index + 1 

              if (mod(k_index,n_tasks).eq.myid) then

                 ! absolute k_points 
                 k_point_abs(:) = 0
                 do i_coord = 1, 3
                    do i_coord2 = 1, 3
                       k_point_abs(i_coord) = k_point_abs(i_coord) + recip_lattice(i_coord,i_coord2)*k_point_rel(i_coord2)
                    end do
                 end do
                 
                 ! diagonalize k-point here ... 
                 call diagonalize_dynamic_matrix(k_point_abs, lattice_vector, hessian, supercell, n_atoms, sc_max, sc_max_total, &
                      sc_start, sc_end, sc_index, sc_prefactor, eigenvalues, posi_frac)
                 
                 if (get_dos) then
                    ! add up dos components
                    do i_dos = 1, dos_fpoints                       
                       freq  = dos_fstart + dble(i_dos-1)*deltaf   
                       freq1 = freq - deltaf/2d0
                       freq2 = freq + deltaf/2d0
                       ev_use(:) = eigenvalues(:)*frequency_factor
                       do i_freq = 1, 3*n_atoms
                          dos(i_dos) = dos(i_dos) + &
                               (arch_erf((freq2-ev_use(i_freq))*dos_factor)-arch_erf((freq1-ev_use(i_freq))*dos_factor))
                       end do
                    end do
                 end if
                 if (get_cv) then
                    ! calculate specific heat components from this particular system
                    do i_freq = 1, 3*n_atoms
                       freq = eigenvalues(i_freq)
                       omega  = const_h*freq/(2d0*pi)
                       ZPE    = ZPE + omega/2
                       ! calculate contribution to each temperature of a given mode and add it into cv ... 
                       do iT = cv_points_start, cv_points
                          T               = cv_start + dble(iT-1)*deltaT
                          omegaT          = omega/(boltzmann_kB*T)
                          if (omegaT.gt.10d0) then     ! numerical stabilizer for integrands - split into two pieces
                             exp_term            = exp(-omegaT)
                             cv             (iT) = cv(iT) +  omegaT*omegaT*exp_term/((1d0-exp_term)**2d0)
                                 free_energy(iT) = free_energy(iT) + T*log(1d0-exp_term)
                             internal_energy(iT) = internal_energy(iT) + omega*exp_term/(1d0-exp_term)
                          else if (omegaT.gt.0d0) then ! don't count soft phonons, numerical or otherwise
                             exp_term            = exp(omegaT)
                             cv             (iT) = cv(iT)          + omegaT*omegaT*exp_term/((exp_term-1d0)**2d0)
                                 free_energy(iT) = free_energy(iT) - T*omegaT + T*log(exp_term-1d0)
                             internal_energy(iT) = internal_energy(iT) + omega/(exp_term-1)
                          end if
                       end do
                    end do
                 end if
              end if
           end do
        end do
     end do

     ! add integration weight and prefactor at end
     if (get_dos) then
        dos(:) = dos(:)*k_weight/deltaf
        call sync_vector(dos,dos_fpoints)
     end if

     if (get_cv ) then
        call sync_vector(ZPE,1)
        call sync_vector(             cv,cv_points)
        call sync_vector(    free_energy,cv_points)
        call sync_vector(internal_energy,cv_points)
        ZPE                = ZPE*k_weight
        cv         (:)     = boltzmann_kB*         cv(:)*k_weight        
        free_energy(:)     = boltzmann_kB*free_energy(:)*k_weight + ZPE(1)
        internal_energy(:) =          internal_energy(:)*k_weight + ZPE(1)
     end if

     ! output dos if on thread zero
     if (myid.eq.0) then
        if (get_dos) then
           write(use_unit,*)
           write(use_unit,*) 'Writing density of states to file phonon_DOS.dat'
           open(unit = 20, file = 'phonon_DOS.dat', status = 'unknown')
           write(20,*) '# Phonon density of states'
           write(20,*) '# Frequency (THz)     DOS'
           do i_dos = 1, dos_fpoints
              freq  = dos_fstart + dble(i_dos-1)*deltaf   
              write(20,'(2E18.8E4)') freq, dos(i_dos) 
           end do
           close(unit = 20)
        end if
        if (get_cv) then
           write(use_unit,*)
           write(use_unit,*) 'Writing phonon specific heat to file phonon_free_energy.dat'
           open(unit = 20, file = 'phonon_free_energy.dat', status = 'unknown')
           write(20,'(2X,A)')           '# Phonon free energy and specific heat'
           write(20,'(2X,A,E18.8E4,A)') '# Phonon zero point energy = ',ZPE/const_eV,' eV'
           write(20,'(2X,A)')           '# Temperature (K)     Free energy (eV/cell)   Internal energy(eV/cell)   cv (kB/unit cell)       -TS_vib(eV/cell)'
           do iT = 1, cv_points
              T = cv_start + dble(iT-1)*deltaT
              write(20,'(F10.2,12X,2E24.14E4,3X,2E24.14E4)') T, free_energy(iT)/const_eV, internal_energy(iT)/const_eV, &
                   cv(iT)/boltzmann_kB, (free_energy(iT)-internal_energy(iT))/const_eV
           end do
           close(unit = 20)
        end if
     end if
     
     if (allocated(dos        )) deallocate(dos   )
     if (allocated(ev_use     )) deallocate(ev_use)
     if (allocated(cv         )) deallocate(cv    )
     if (allocated(free_energy)) deallocate(free_energy)
  end if

  ! deallocate data and quit
  if (myid.eq.0) then
     write(use_unit,*)
     write(use_unit,*) "Deallocating data ... "
  end if
  deallocate(eigenvalues)
  deallocate(sc_index)
  deallocate(hessian_line)
  deallocate(band_points)
  deallocate(band_info)
  deallocate(species_names)
  deallocate(species)
  deallocate(species_masses)
  deallocate(masses)
  deallocate(hessian)
  deallocate(dynamic_matrix)

  if (myid.eq.0) then
     write(use_unit,*) 
     write(use_unit,*) "Leaving Hessian diagonalization program for FHI-aims"
     write(use_unit,*) 
     write(use_unit,*) "=========================================================="
  end if

  call finalize_mpi ( )

end program phonon_hessian_diagonalization

! subroutine to obtain the dynamic matrix from the input data, given some k-point as input
subroutine diagonalize_dynamic_matrix(k_point, lattice_vector, hessian, supercell, n_atoms, sc_max, sc_max_total, sc_start, &
     sc_end, sc_index, sc_prefactor, eigenvalues, posi_frac)
  implicit none
  integer :: n_atoms, sc_max_total, sc_max
  integer,    dimension(3)                   :: supercell
  real*8,     dimension(3)                   :: k_point, R, R_image,lattice_length
  real*8,     dimension(3)                   :: dist
  real*8,     dimension(3,3)                 :: lattice_vector
  real*8,     dimension(3,n_atoms)           :: posi_frac
  real*8,     dimension(3*n_atoms)           :: eigenvalues
  complex*16, dimension(3,n_atoms,3,n_atoms) :: dynamic_matrix
  integer,    dimension(3)                   :: sc_start, sc_end   
  complex*16, dimension(9*n_atoms)           :: workspace
  real*8,     dimension(9*n_atoms)           :: r_workspace 
  real*8,     dimension(3,sc_max_total)      :: sc_prefactor
  integer,    dimension(3,-sc_max:sc_max)    :: sc_index
  real*8,     dimension(3,n_atoms,3,n_atoms,supercell(1),supercell(2),supercell(3)) :: hessian
  integer    :: isc1,isc2,isc3 ,isc1_corr,isc2_corr,isc3_corr
  integer    :: i_coord, info, i_atom, i_coord_2, i_atom_2, n_image
  integer    :: ext_isc1, ext_isc2, ext_isc3,num_ext_sc
  complex*16 :: sc_phase
  real*8     :: bufreal, bufimaginary, R_length, R_image_length,diag_length
  real*8     :: pi = 3.14159265358979323846d0

  ! calculate unit cell lattice (a, b, and c) length
  lattice_length(1)=sqrt(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 + lattice_vector(3,1)**2)
  lattice_length(2)=sqrt(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 + lattice_vector(3,2)**2)
  lattice_length(3)=sqrt(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 + lattice_vector(3,3)**2)

  
  ! calculate dynamic matrix based on the index arrays obtained before
  dynamic_matrix(:,:,:,:) = 0d0
! use subroutin extended_supercell to find how many supercells to be used
	call extended_supercell(lattice_vector,supercell,diag_length,num_ext_sc)

    do i_atom = 1, n_atoms
     do i_atom_2 = 1, n_atoms
! calculate D matrix for atom pair i_atom & i_atom_2 in the primitive cell

      do isc1 = sc_start(1), sc_end(1)
       do isc2 = sc_start(2), sc_end(2)
        do isc3 = sc_start(3), sc_end(3)
! for this specific unit cell, use extended supercells to find all the atomms with equivalent distance to the central atom

           do i_coord = 1, 3  ! here i_coord for x, y, z direction
             R(i_coord) = (dble(isc1)+posi_frac(1,i_atom_2) - posi_frac(1,i_atom))*lattice_vector(i_coord,1) &
                        +(dble(isc2)+posi_frac(2,i_atom_2) - posi_frac(2,i_atom))*lattice_vector(i_coord,2) &
                       +(dble(isc3)+posi_frac(3,i_atom_2) - posi_frac(3,i_atom))*lattice_vector(i_coord,3)
           end do
           R_length = sqrt(dot_product(R,R)) ! initialize minimun distance
           n_image = 0      ! initialize # of equivalent atoms
           sc_phase = (0,0) ! initialize complex phase factor 

           do ext_isc1= -num_ext_sc, num_ext_sc
            do ext_isc2= -num_ext_sc, num_ext_sc
             do ext_isc3= -num_ext_sc, num_ext_sc
              do i_coord = 1, 3
              R_image(i_coord)= &
         (ext_isc1*supercell(1)+dble(isc1)+posi_frac(1,i_atom_2) - posi_frac(1,i_atom))*lattice_vector(i_coord,1) &
        +(ext_isc2*supercell(2)+dble(isc2)+posi_frac(2,i_atom_2) - posi_frac(2,i_atom))*lattice_vector(i_coord,2) &
        +(ext_isc3*supercell(3)+dble(isc3)+posi_frac(3,i_atom_2) - posi_frac(3,i_atom))*lattice_vector(i_coord,3)
              end do  
              R_image_length=sqrt(dot_product(R_image,R_image))
              if (abs(R_image_length - R_length) .lt. 0.0001) then
                n_image = n_image + 1
                sc_phase = sc_phase + exp((0d0,1d0)*dot_product(R_image,k_point))
              else if (R_image_length .lt. R_length) then
                  n_image = 1
                  R_length = R_image_length
                  sc_phase = exp((0d0,1d0)*dot_product(R_image,k_point))
              end if
             end do  ! ext_isc3
            end do
           end do    ! ext_isc1
           
           sc_phase = sc_phase/n_image  ! average sc_phase

           dynamic_matrix(:,i_atom,:,i_atom_2) = dynamic_matrix(:,i_atom,:,i_atom_2) &
           + sc_phase*hessian(:,i_atom,:,i_atom_2,sc_index(1,isc1),sc_index(2,isc2),sc_index(3,isc3))
! Now summation of dynamic matrix is done for one specific unit cell  
           
        end do ! isc(3)
      end do   ! isc(2)
    end do     ! isc(1)
   end do  ! i_atom_2
  end do   ! i_atom

 ! symmetrize dynamic matrix 
  do i_atom = 1, n_atoms
     do i_atom_2 = 1, i_atom
        do i_coord = 1, 3
           do i_coord_2 = 1, i_coord
              bufreal      = ( real(dynamic_matrix(i_coord,i_atom,i_coord_2,i_atom_2)) +  real(dynamic_matrix(i_coord_2,i_atom_2,i_coord,i_atom)))/2d0
              bufimaginary = (aimag(dynamic_matrix(i_coord,i_atom,i_coord_2,i_atom_2)) - aimag(dynamic_matrix(i_coord_2,i_atom_2,i_coord,i_atom)))/2d0
              dynamic_matrix(i_coord,i_atom,i_coord_2,i_atom_2) = (1d0,0d0)*bufreal + (0d0,1d0)*bufimaginary
              dynamic_matrix(i_coord_2,i_atom_2,i_coord,i_atom) = (1d0,0d0)*bufreal - (0d0,1d0)*bufimaginary
           end do
        end do
     end do
  end do

  ! diagonalize dynamic matrix 
  call ZHEEV('V','U',3*n_atoms,dynamic_matrix,3*n_atoms,eigenvalues,workspace,9*n_atoms,r_workspace,info)

  ! turn from omega^2 into frequencies in units of THz
  eigenvalues(:) = sign(sqrt(abs(eigenvalues(:))),eigenvalues(:))

end subroutine diagonalize_dynamic_matrix

! subroutine to transform cartesian to fractional
subroutine cart2frac(lattice_vector,posi_cart,posi_frac,n_atoms)
  implicit none
  integer                       :: n_atoms, i
  real*8,  dimension(3,3)       :: lattice_vector  
  real*8,  dimension(3,n_atoms) :: posi_cart, posi_frac  
  real*8                        :: dert0,dert1,dert2,dert3

  dert0 = lattice_vector(1,1)*lattice_vector(2,2)*lattice_vector(3,3) + lattice_vector(1,3)*lattice_vector(2,1)*lattice_vector(3,2) + &
          lattice_vector(1,2)*lattice_vector(2,3)*lattice_vector(3,1) - lattice_vector(1,3)*lattice_vector(2,2)*lattice_vector(3,1) - &
          lattice_vector(1,1)*lattice_vector(2,3)*lattice_vector(3,2) - lattice_vector(1,2)*lattice_vector(2,1)*lattice_vector(3,3)
  do i = 1, n_atoms
   dert1 = posi_cart(1,i)*lattice_vector(2,2)*lattice_vector(3,3) + lattice_vector(1,3)*posi_cart(2,i)*lattice_vector(3,2) + &
           lattice_vector(1,2)*lattice_vector(2,3)*posi_cart(3,i) - lattice_vector(1,3)*lattice_vector(2,2)*posi_cart(3,i) - &
           posi_cart(1,i)*lattice_vector(2,3)*lattice_vector(3,2) - lattice_vector(1,2)*posi_cart(2,i)*lattice_vector(3,3)
   dert2 = lattice_vector(1,1)*posi_cart(2,i)*lattice_vector(3,3) + lattice_vector(1,3)*lattice_vector(2,1)*posi_cart(3,i) + &
           posi_cart(1,i)*lattice_vector(2,3)*lattice_vector(3,1) - lattice_vector(1,3)*posi_cart(2,i)*lattice_vector(3,1) - &
           lattice_vector(1,1)*lattice_vector(2,3)*posi_cart(3,i) - posi_cart(1,i)*lattice_vector(2,1)*lattice_vector(3,3)
   dert3 = lattice_vector(1,1)*lattice_vector(2,2)*posi_cart(3,i) + posi_cart(1,i)*lattice_vector(2,1)*lattice_vector(3,2) + &
           lattice_vector(1,2)*posi_cart(2,i)*lattice_vector(3,1) - posi_cart(1,i)*lattice_vector(2,2)*lattice_vector(3,1) - &
           lattice_vector(1,1)*posi_cart(2,i)*lattice_vector(3,2) - lattice_vector(1,2)*lattice_vector(2,1)*posi_cart(3,i)
   posi_frac(1,i) = dert1/dert0
   posi_frac(2,i) = dert2/dert0
   posi_frac(3,i) = dert3/dert0
  end do 
end subroutine cart2frac

! subroutine to find how large the extended supercell should be
	subroutine extended_supercell(lattice_vector,supercell,diag_length,num_ext_sc)
	implicit none
	integer, dimension(3) :: supercell
	real*8, dimension(3) :: sc_a,sc_b,sc_c,diagonal
	real*8, dimension(3,3) :: lattice_vector
	real*8 :: sc_a_length,sc_b_length,sc_c_length,diag_length,diag_length_new
	real*8 :: dot_sc_ab,dot_sc_bc,dot_sc_ca,alpha,beta,gamma
	real*8 :: sc_a_height,sc_b_height,sc_c_height,sc_volume
	integer :: num_ext_sc,i,j

	sc_a(:) = lattice_vector(:,1)*supercell(1)  ! define supercell lattice
	sc_b(:) = lattice_vector(:,2)*supercell(2)
	sc_c(:) = lattice_vector(:,3)*supercell(3)
	sc_a_length = sqrt(dot_product(sc_a,sc_a))
        sc_b_length = sqrt(dot_product(sc_b,sc_b))
        sc_c_length = sqrt(dot_product(sc_c,sc_c))

! now find the longest of the four body diagonals of the supercell
	diagonal(:) = sc_a(:) + sc_b(:) + sc_c(:)
        diag_length = sqrt(dot_product(diagonal,diagonal))
	diagonal(:) = -sc_a(:) + sc_b(:) + sc_c(:)
        diag_length_new = sqrt(dot_product(diagonal,diagonal))
        if (diag_length_new .gt. diag_length) then
         diag_length = diag_length_new
        end if
        diagonal(:) = sc_a(:) - sc_b(:) + sc_c(:)
	diag_length_new = sqrt(dot_product(diagonal,diagonal))
        if (diag_length_new .gt. diag_length) then
         diag_length = diag_length_new
        end if
        diagonal(:) = sc_a(:) + sc_b(:) - sc_c(:)
	diag_length_new = sqrt(dot_product(diagonal,diagonal))
        if (diag_length_new .gt. diag_length) then
         diag_length = diag_length_new
        end if
!  now calculate diameter of the inner contagent sphere within the supercell
	dot_sc_ab = sc_a(1)*sc_b(1) + sc_a(2)*sc_b(2) + sc_a(3)*sc_b(3) 
        dot_sc_bc = sc_b(1)*sc_c(1) + sc_b(2)*sc_c(2) + sc_b(3)*sc_c(3)
        dot_sc_ca = sc_c(1)*sc_a(1) + sc_c(2)*sc_a(2) + sc_c(3)*sc_a(3)

	alpha = acos(dot_sc_bc/(sc_b_length*sc_c_length)) 
	beta  = acos(dot_sc_ca/(sc_c_length*sc_a_length)) 
	gamma = acos(dot_sc_ab/(sc_a_length*sc_b_length)) 
        sc_volume= sc_a_length*sc_b_length*sc_c_length &
                      *sqrt(1+2*cos(alpha)*cos(beta)*cos(gamma)-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma))
	sc_a_height=sc_volume/(sc_b_length*sc_c_length*sin(alpha))
        sc_b_height=sc_volume/(sc_c_length*sc_a_length*sin(beta))
	if (sc_b_height .lt. sc_a_height) then
	 sc_a_height = sc_b_height
	end if
	sc_c_height=sc_volume/(sc_a_length*sc_b_length*sin(gamma))
	if (sc_c_height .lt. sc_a_height) then
         sc_a_height = sc_c_height
        end if
! sc_a_height is the shortest height of the supercell, i.e., the sphere diameter

	num_ext_sc = int(2*diag_length/sc_a_height) + 1 ! # of extended supercell
	end subroutine extended_supercell




