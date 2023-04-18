!
!  small code to take the Hessian matrix computed by aims as well as the 
!  masses of all atoms involved and calculate the vibration frequencies of 
!  the molecule/cluster in question
!
!  First version: Felix Hanke, 2007
!
!    mkl compiler options on the thlc's:
! ifort -o hessian_diagonalization_jmol.x hessian_diagonalization_jmol.f90 -L/opt/intel/mkl/8.1/lib/32/ -lmkl_lapack -lmkl_ia32 -lmkl
!
!  The basic procedure is described in R.M.Martin, "Electronic Structure - Basic Theory and Practical Methods", Ch 19.1
!  For details on the matrix symmetrization and the calculation of the eigenvectors, HJ Kreuzer/ZW Gortel "Physisorption Kinetics"
!          contains a similar but much more detailed discussion in their chapter on the 'Mesoscopic Master Equation'
!                       

!
program hessian_diagonalization
use MPI_ROUTINES
implicit none
  character*80 :: name_inputhessian, name_inputmasses, inputbuf, name_grad_dipole, name_xyzbuffer, name_ir_buffer, name_grad_polar, name_raman_buffer
  integer :: i_atoms, i_coords, ix, iy, n_eigenvalues, i_atoms_2, i_coords_2, lwork, info, n_atoms
  real*8, dimension(:,:), allocatable :: inputhessian, hessian, atomic_positions
  real*8, dimension(:,:), allocatable :: grad_dipole, grad_polar
  real*8, dimension(:), allocatable :: atomic_masses, eigenvalues, workspace, mass_vector, frequencies, reduced_mass
  real*8 :: const_u, const_eV, const_Angstr, hessian_factor, const_thresh, const_c, buf, mass, grad_dipole_factor, const_N_avogadro, pi, raman_factor
  real*8 :: const_hbar, const_kB, temp, press
  real*8, dimension(3,3) :: I_matrix
  real*8, dimension(3)   :: I_eigenvalues, CMposition
  character*10, dimension(:), allocatable :: element_names
  character*22, dimension(:), allocatable :: symmetries
  real*8, dimension(:,:), allocatable ::grad_dipole_internal, grad_polar_internal
  real*8, dimension(:), allocatable ::infrared_intensity, raman_intensity, raman_iso, raman_aniso
  real*8 :: ir_factor, norm, ZPE, ZPE_cumulative, ZPE_first_six, ZPE_first_three, vib_free_E, rot_free_E, trans_free_E
  integer :: i_coords_3, i, j
  integer :: periodic

  real*8, external :: ddot

  ! The vibrational and rotational parts are safe and can always be printed.
  logical :: flag_free_energy = .true.
  !
  ! The translational part is the one that depends on the external pressure and that can be quite confusing.
  ! Hence this will only be printed if the user requests it.
  !
  logical :: flag_trans_free_energy = .false.
  !
  ! VB: 298.15 K are standard ambient conditions (textbook, e.g., Wedler.)
  !     Wedler also gives 1.013 bar as the standard pressure; the "extra" 25 are 
  !     presumably the roundoff error; in any case, 1.01325 is used by NIST.
  !
  real*8  :: Tstart = 0.0, Tend = 298.15, pstart = 101325.0, pend = 101325.0
  integer :: Tpoints = 2, ppoints = 1

  ! initialize MPI and only do something if we are on the 0-CPU ... 
  call initialize_mpi ( )
  call get_my_task()
  if (myid.eq.0) then
    ! start work
    write(STDOUT,*) 'Entering program hessian_diagonalization'
    write(STDOUT,*) 'For FHI-AIMS, Version XXXXXX'
  
    ! various physical constants and conversion factors
    const_u          = 1.66053886d-27            ! http://www.physics.nist.gov/cuu/Constants/index.html
    const_eV         = 1.60217653d-19
    const_c          = 299792458d0
    const_Angstr     = 1d-10
    const_hbar       = 1.054571628d-34 
    const_N_avogadro = 6.02214179d23
    pi               = 3.14159265d0
    const_kB         = 8.6173324e-5 ! in eV/K
    hessian_factor   = const_eV/(const_u*const_Angstr*const_Angstr)  ! factor to bring hessian into SI units ... 
    !  grad_dipole_factor = const_eV ! bring grad_dipole to SI
    grad_dipole_factor = 4.80320577161 ! bring grad_dipole to D/Ang
    !ir_factor = const_N_avogadro * pi / (3.d0 * const_c)
    ir_factor = 1
    raman_factor = 0.078415972 ! To obtain raman_intensity in Ang^4/amu
    
    ! extract input parameters: number of atoms, input files
    call getarg(1,inputbuf)
    read(inputbuf,*) n_atoms
    call getarg(2,inputbuf)
    read(inputbuf,*) name_inputhessian
    call getarg(3,inputbuf)
    read(inputbuf,*) name_grad_dipole
    call getarg(4,inputbuf)
    read(inputbuf,*) name_grad_polar
    call getarg(5,inputbuf)
    read(inputbuf,*) name_inputmasses
    call getarg(6,inputbuf)
    read(inputbuf,*) const_thresh
    call getarg(7,inputbuf)
    read(inputbuf,*) name_xyzbuffer
    call getarg(8,inputbuf)
    read(inputbuf,*) name_ir_buffer
    call getarg(9,inputbuf)
    read(inputbuf,*) name_raman_buffer
    call getarg(10,inputbuf)
    read(inputbuf,*) periodic
    write(STDOUT,*) 'Number of atoms                = ', n_atoms
    write(STDOUT,*) 'Name of Hessian input file     = ', trim(name_inputhessian)
    write(STDOUT,*) 'Name of grad dipole input file = ', trim(name_grad_dipole)
    write(STDOUT,*) 'Name of grad polarizability input file = ', trim(name_grad_polar)
    write(STDOUT,*) 'Name of Masses  input file     = ', trim(name_inputmasses)
    write(STDOUT,*) 'Name of XYZ output file        = ', trim(name_xyzbuffer)
    write(STDOUT,*) 'Threshold for Matrix elements  = ', const_thresh
    if (const_thresh < 0d0) write(STDOUT,*) '     All matrix elements are taken into account by default'
    write(STDOUT,*)
    
    ! now know array sizes, allocate input array and working array for 
    ! diagonalization routines ... 
    allocate(atomic_masses(n_atoms))
    allocate(inputhessian(n_atoms*3,n_atoms*3))
    allocate(grad_dipole(3,n_atoms*3))
    !allocate(grad_polar(n_atoms*3, 6))
    allocate(grad_polar(6,n_atoms*3))
    allocate(hessian(n_atoms*3,n_atoms*3))
    allocate(eigenvalues(3*n_atoms))
    allocate(reduced_mass(3*n_atoms))
    allocate(mass_vector(3*n_atoms))
    allocate(symmetries(3*n_atoms))
    allocate(atomic_positions(n_atoms,3))
    allocate(frequencies(n_atoms*3))
    allocate(element_names(3*n_atoms))
    allocate(grad_dipole_internal(3, n_atoms*3))
    allocate(grad_polar_internal(6, n_atoms*3))
    allocate(infrared_intensity(n_atoms*3))
    allocate(raman_intensity(n_atoms*3))
    allocate(raman_iso(n_atoms*3))
    allocate(raman_aniso(n_atoms*3))
    lwork = 9*n_atoms
    allocate(workspace(lwork))

    ! read input from control.in-file
    call read_c
    
    ! read input hessian from file 
    open(20,FILE=name_inputhessian)
    do iy = 1, 3*n_atoms, 1
       read(20,*) (inputhessian(iy,ix),ix = 1, 3*n_atoms, 1)
       do ix = 1, 3*n_atoms, 1
          if (abs(inputhessian(iy,ix)).lt.const_thresh) inputhessian(iy,ix) = 0d0
       end do
    end do
    close(20)
    
    ! read input grad dipole from file 
    open(20,FILE=name_grad_dipole)
    do iy = 1, 3*n_atoms, 1
       read(20,*) (grad_dipole(ix,iy),ix = 1, 3, 1)
    end do
    close(20)

    open(20,FILE=name_grad_polar)
    do iy = 1, 3*n_atoms, 1
       read(20,*) (grad_polar(ix,iy),ix = 1, 6, 1)
    end do
    close(20)
   
 
    ! read input atomic masses and atomic coordinates
    open(20,FILE=name_inputmasses)
    do i_atoms = 1, n_atoms
       read(20,*) atomic_masses(i_atoms), (atomic_positions(i_atoms,i_coords), i_coords = 1, 3, 1), element_names(i_atoms)
       do i_coords = 1, 3
          mass_vector(3*i_atoms+i_coords-3) = 1d0/sqrt(atomic_masses(i_atoms))
       end do
    end do
    close(20)
    
    ! calculate the symmetric matrix for appropriate LAPACK routine
    ! change everything to SI!!!!!!!! - right now, we have 
    ! [HESSIAN] = eV / Angstrom^2
    ! [masses]  = atomic mass units ... 
    ! This matrix is symmetric (the original matrix from the vibration problem is NOT),
    ! at the price that the displacement eigenvectors have the form sqrt(Mass_I)u(atom_I)
    ! - after the diagonalization, one should divide by the square root of the mass to 
    !   get the appropriate eigenvectors...
    hessian(:,:) = inputhessian(:,:)                    ! initialize hessian matrix
    do i_coords = 1, 3*n_atoms, 1                        ! normalize with the inverse square root of the masses .... 
       hessian(:,i_coords) = hessian(:,i_coords)*mass_vector(i_coords)       ! .... and the proper unit conversion factor
       hessian(i_coords,:) = hessian(i_coords,:)*hessian_factor*mass_vector(i_coords) 
   end do
    
    grad_dipole(:,:) = grad_dipole(:,:) * grad_dipole_factor
    
    do i_coords = 1, 3*n_atoms
       do i_coords_2 = 1, i_coords - 1 
          buf = (hessian(i_coords_2,i_coords)+hessian(i_coords,i_coords_2))/2d0
          hessian(i_coords_2,i_coords) = buf
          hessian(i_coords,i_coords_2) = buf
       end do
    end do
    
    ! diagonalize matrix using LAPACK routine for a real, symmetric matrix - including the calculation of eigenvectors etc
    write(STDOUT,*) 'Solving eigenvalue system for Hessian Matrix'
    call DSYEV('V','U',3*n_atoms,hessian,3*n_atoms,eigenvalues,workspace,lwork,info)
    write(STDOUT,*) 'Done ... '
    write(STDOUT,*)
    
    ! calculate the eigenvectors in cartesian coordinates ?
    do i_coords = 1, 3*n_atoms
!       norm = sqrt(ddot(3*n_atoms, hessian(:, i_coords), 1, hessian(:, i_coords), 1))
       hessian(:,i_coords) = hessian(:,i_coords)*mass_vector(:)
    end do
    
    ! transform dipole derivative to internal coordinates via directional derivative
    ! d/dQ = d/dR * Q_normalized, where Q_normalized is nothing but displacement eigenvector
    ! hessian(:,:) not normalized anymore since mass was divided out
    do i_coords = 1, 3
       ! loop over modes
       do i_coords_2 = 1, 3*n_atoms
          grad_dipole_internal(i_coords, i_coords_2) = ddot(3*n_atoms, grad_dipole(i_coords,:), 1, hessian(:, i_coords_2), 1)
       end do
    end do
    
    ! transform polarizability derivative to internal coordinates via directional derivative
    ! d/dQ = d/dR * Q_normalized, where Q_normalized is nothing but displacement eigenvector
    ! hessian(:,:) not normalized anymore since mass was divided out
    do i_coords = 1, 6
       ! loop over modes
       do i_coords_2 = 1, 3*n_atoms
          grad_polar_internal(i_coords, i_coords_2) = ddot(3*n_atoms, grad_polar(i_coords,:), 1, hessian(:, i_coords_2), 1)
       end do
    end do
    
    ! get infrared intensities
    do i_coords = 1, 3*n_atoms, 1
       infrared_intensity(i_coords) = ddot(3, grad_dipole_internal(:, i_coords), 1, grad_dipole_internal(:, i_coords), 1)
    end do
    
    ! get Raman intensities
    do i_coords = 1, 3*n_atoms, 1  !these are modes
       ! 1) Isotropic part
       do i_coords_2 = 1, 3  !these are cartesian coordinates
         !raman_intensity(i_coords) = ddot(3, grad_polar_internal(:, i_coords), 1, grad_polar_internal(:, i_coords), 1)
         raman_iso(i_coords) = raman_iso(i_coords) + grad_polar_internal(i_coords_2, i_coords)
         !raman_intensity(i_coords) = raman_intensity(i_coords) + grad_polar_internal(i_coords_2, i_coords)*grad_polar_internal(i_coords_2, i_coords)
       end do
       raman_iso(i_coords) = (raman_iso(i_coords)*1./3)**2
       ! 2) Anisotropic part
       raman_aniso(i_coords) = (grad_polar_internal(1,i_coords)-grad_polar_internal(2,i_coords))**2 &
                             + (grad_polar_internal(2,i_coords)-grad_polar_internal(3,i_coords))**2 &
                             + (grad_polar_internal(3,i_coords)-grad_polar_internal(1,i_coords))**2 &
                             + 6*grad_polar_internal(4,i_coords)**2 &
                             + 6*grad_polar_internal(5,i_coords)**2 &
                             + 6*grad_polar_internal(6,i_coords)**2
       ! The Raman intensity is calculated as in Eq.(42) in Neugebauer et al., J Comput Chem 23: 895–910, 2002
       raman_intensity(i_coords) = raman_factor*(45*raman_iso(i_coords)+7./2*raman_aniso(i_coords))
    end do

    
    ! scale infrared intensities
    infrared_intensity(:) = infrared_intensity(:) * ir_factor

    ! Renormalize eigenvectors for output - norm has units of 1/sqrt(mass) though 
    do i_coords = 1, 3*n_atoms
!      hessian(:,i_coords) = hessian(:,i_coords)/mass_vector(:)
      reduced_mass(i_coords) = ddot(3*n_atoms, hessian(:, i_coords), 1, hessian(:, i_coords), 1)
      norm = sqrt(reduced_mass(i_coords))
      hessian(:,i_coords) = hessian(:,i_coords)/norm
    end do
 

   
    ! check results according to the value of omega^2 ... 
    write(STDOUT,*) 'Results: '
    write(STDOUT,*)
    write(STDOUT,*) 'List of all frequencies found:'
    !write(STDOUT,'(A13,2A25,A27)') 'Mode number','Frequency [cm^(-1)]','Zero point energy [eV]','IR-intensity [D^2/Ang^2]'
    write(STDOUT,'(A13,2A25,A27,A35)') 'Mode number','Frequency [cm^(-1)]','Zero point energy [eV]','IR-intensity [D^2/Ang^2]', 'Raman intensity [Angstrom**4/amu]'
    ZPE_cumulative = 0d0
    ZPE_first_six  = 0d0
    ZPE_first_three  = 0d0
    do i_coords = 1, 3*n_atoms, 1
       if (eigenvalues(i_coords).gt.0d0) then
          frequencies(i_coords) = sqrt(eigenvalues(i_coords))
          symmetries(i_coords) = 'stable frequency at '
       else if (eigenvalues(i_coords).lt.0d0) then
          frequencies(i_coords) = -sqrt(-eigenvalues(i_coords))
          symmetries(i_coords) = 'unstable frequency at '
       else 
          symmetries(i_coords) = 'translation or rotation '
          frequencies(i_coords) = 0d0
       end if
       ZPE = const_hbar*frequencies(i_coords)/(2d0*const_eV)
       ZPE_cumulative = ZPE_cumulative + ZPE
       if (i_coords.le.6) then
          ZPE_first_six  = ZPE_first_six + ZPE
       end if
       if (i_coords.le.3) then
          ZPE_first_three  = ZPE_first_three + ZPE
       end if
       write(STDOUT,'(I13,2F25.8,F27.8,F27.8)') i_coords, frequencies(i_coords)/(200*pi*const_c),ZPE,infrared_intensity(i_coords), raman_intensity(i_coords)
    end do
    write(STDOUT,*)
    write(STDOUT,*) 'Summary of zero point energy for entire system:'
    write(STDOUT,'(2X,A,F15.8,A)') '| Cumulative ZPE               = ',ZPE_cumulative, ' eV'
    if (.not. periodic) then
       write(STDOUT,'(2X,A,F15.8,A)') '| without first six eigenmodes = ',ZPE_cumulative-ZPE_first_six,' eV'
    else
      write(STDOUT,'(2X,A,F15.8,A)') '| without first three eigenmodes = ',ZPE_cumulative-ZPE_first_three,' eV'
    endif
    write(STDOUT,*)
    write(STDOUT,*)

    write(STDOUT,*) 'Stability checking - eigenvalues should all be positive for a stable structure. '
    if (.not. periodic) then
       write(STDOUT,*) 'The six smallest frequencies should be (almost) zero for a nonperiodic system:'
       write(STDOUT,'(62F25.8)')  frequencies(1)/(200*pi*const_c), &
                          &  frequencies(2)/(200*pi*const_c), &
                          &  frequencies(3)/(200*pi*const_c), &
                          &  frequencies(4)/(200*pi*const_c), &
                          &  frequencies(5)/(200*pi*const_c), &
                          &  frequencies(6)/(200*pi*const_c) 
    else
       write(STDOUT,*) 'The three smallest frequencies should be (almost) zero for a periodic system:'
       write(STDOUT,'(62F25.8)')  frequencies(1)/(200*pi*const_c), &
                          &  frequencies(2)/(200*pi*const_c), &
                          &  frequencies(3)/(200*pi*const_c) 
    endif
    write(STDOUT,*)
    write(STDOUT,*)
    write(STDOUT,*) 'Compare this with the largest eigenvalue, '
!    write(STDOUT,'(E14.4E3)') eigenvalues(3*n_atoms) 
    write(STDOUT,'(F25.8)')  frequencies(3*n_atoms)/(200*pi*const_c)
    
    ! output the results to the XYZ file 
    open(unit = 20, file = name_xyzbuffer)
    do i_coords = 1, 3*n_atoms
       write(20,*) n_atoms
       write(20,'(A,f10.3,A,E10.4,A,f10.3,A,f5.3,A,f5.3,A)') &
        symmetries(i_coords),frequencies(i_coords)/(200*pi*const_c),' 1/cm IR int. is ',&
        infrared_intensity(i_coords), ' D^2/Ang^2; Raman int. is ', raman_intensity(i_coords), &
        ' Ang^4/amu; red. mass is ', 1d0/reduced_mass(i_coords), &
        ' a.m.u.; force const. is ', (frequencies(i_coords))**2*1d0/reduced_mass(i_coords) * const_u * 1d-2, ' mDyne/Ang.'
       do i_atoms = 1, n_atoms
          write(20,'(2X,A10,6f10.4)') element_names(i_atoms),               &
               (atomic_positions(i_atoms,i_coords_2), i_coords_2 = 1, 3, 1), &
               (hessian(i_coords_2,i_coords), i_coords_2 = (i_atoms-1)*3+1,(i_atoms-1)*3+3,1)
       end do
    end do
    close(unit = 20)
    
    ! output the results to the ir file
    open(unit = 20, file = name_ir_buffer)
    
    write (20,'(A)') "[INT]"
    do i_coords = 1, 3*n_atoms
       write (20,'(E10.4)') infrared_intensity(i_coords)
    end do
    close(unit = 20)
    
    open(unit = 21, file = name_raman_buffer)
    
    write (21,'(A)') "[INT]"
    do i_coords = 1, 3*n_atoms
       write (21,'(E10.4)') raman_intensity(i_coords)
    end do
    close(unit = 21)

    ! end diagonalizing hessian
    ! ---------------------------------------------------------------------------------------------------------------
    ! begin diagonalization of Moment of inertia tensor 
    
    ! calculate center of mass position, starting with the total mass ...
    mass = 0d0
    do i_atoms = 1, n_atoms
       atomic_masses(i_atoms) = atomic_masses(i_atoms)*const_u             ! change into SI units of kg .... 
       mass = mass + atomic_masses(i_atoms)
    end do
    
    ! change atomic positions into SI as well - just remember to change it back when the final output comes ... 
    atomic_positions(:,:) = atomic_positions(:,:)*const_Angstr
    CMposition(:) = 0d0
    do i_atoms = 1, n_atoms
       CMposition(:) = CMposition(:) + atomic_masses(i_atoms)*atomic_positions(i_atoms,:)
    end do
    CMposition(:) = CMposition(:)/mass
    ! change coordinate system into CMCS
    do i_atoms = 1, n_atoms
       atomic_positions(i_atoms,:) = atomic_positions(i_atoms,:) - CMposition(:)
    end do
    
    ! now, calculate all the elements of the moment of inertia tensor
    I_matrix(:,:) = 0d0
    do i_atoms = 1, n_atoms
       I_matrix(1,1) = I_matrix(1,1) + &
            atomic_masses(i_atoms)*(atomic_positions(i_atoms,2)**2+atomic_positions(i_atoms,3)**2)
       I_matrix(1,2) = I_matrix(1,2) - atomic_masses(i_atoms)*atomic_positions(i_atoms,1)*atomic_positions(i_atoms,2)
       I_matrix(1,3) = I_matrix(1,3) - atomic_masses(i_atoms)*atomic_positions(i_atoms,1)*atomic_positions(i_atoms,3)
       I_matrix(2,2) = I_matrix(2,2) + &
            atomic_masses(i_atoms)*(atomic_positions(i_atoms,1)**2+atomic_positions(i_atoms,3)**2)
       I_matrix(2,3) = I_matrix(2,3) - atomic_masses(i_atoms)*atomic_positions(i_atoms,2)*atomic_positions(i_atoms,3)
       I_matrix(3,3) = I_matrix(3,3) + &
            atomic_masses(i_atoms)*(atomic_positions(i_atoms,2)**2+atomic_positions(i_atoms,1)**2)
    end do
    I_matrix(2,1) = I_matrix(1,2)
    I_matrix(3,1) = I_matrix(1,3)
    I_matrix(3,2) = I_matrix(2,3)
    
    ! diagonalize moment of inertia tensor
    call DSYEV('V','U',3,I_matrix,3,I_eigenvalues,workspace,lwork,info)
    
    ! output the moments of inertia, their product, and mass information
    write(STDOUT,*)
    write(STDOUT,'(A,2X,F18.8)')      'Molecular mass [g/mol]:                                ', &
           mass*1000*6.02214129E23 !Should Avogadro's constant go to constants?
    write(STDOUT,'(A,2X,E18.8E4)')    'Total mass [kg]                                        ', &
          mass 
    write(STDOUT,'(A,2X,3E18.8E4,A)') 'Principal moments of inertia [kg m^2]                  ', & 
          I_eigenvalues
    write(STDOUT,'(A,2X,E18.8E4)')    'Product of the principal moments of inertia [kg^3 m^6] ',&
         I_eigenvalues(1)*I_eigenvalues(2)*I_eigenvalues(3)
    write(STDOUT,*)
    write(STDOUT,*) 'Eigenvectors of the inertia tensor'
    do i_coords = 1, 3
       write(STDOUT,'(3f10.4)') (I_matrix(i_coords,ix), ix = 1, 3) 
    end do



    ! output the free energy
    if (flag_free_energy) then

       write(STDOUT,*)
       write(STDOUT,'(A)') 'Output of some temperature-dependent free energy contributions.'
       if (.not. flag_trans_free_energy) then
          write(STDOUT,'(A20,2A35)') 'temperature [K]', 'vibrational free energy [eV]' , 'rotational free energy [eV]'
       end if

       do i=1,Tpoints

          if (Tpoints.eq.1) then
            ! legislate that the end point counts
            temp = Tend
          else
            ! We use all points from Tstart to Tend. 
            ! Note the off-by-one, but the number of points is now correct
            temp = Tstart + (Tend-Tstart)*dble(i-1)/dble(Tpoints-1)
          end if

          ! calculate vibrational free energy
          vib_free_E = 0
             do i_coords = 1, 3*n_atoms, 1
                ZPE = const_hbar*frequencies(i_coords)/(2d0*const_eV)
                if (.not. periodic) then
                   if (i_coords.gt.6) then
                      vib_free_E = vib_free_E + ZPE + (const_kB*temp)*log(1-exp(-2*ZPE/(const_kB*temp)))
                   end if
                else
                   if (i_coords.gt.3) then
                      vib_free_E = vib_free_E + ZPE + (const_kB*temp)*log(1-exp(-2*ZPE/(const_kB*temp)))
                   end if
                endif
             end do

          ! calculate free energy of a rigid rotor
          if (temp.gt.0.d0) then
            rot_free_E = -const_kB*temp*(0.5*log(pi) - 1.5*log(const_hbar**2/(2.0*const_kB*const_eV)) + &
              0.5*log(I_eigenvalues(1)*I_eigenvalues(2)*I_eigenvalues(3)) + 1.5*log(temp))
          else
            rot_free_E = 0.d0
          end if

          ! calculate translational contribution to free energy
          if (flag_trans_free_energy) then
            do j=1,ppoints

               if (ppoints.eq.1) then
                 press = pend
               else
                 press = pstart + (pend-pstart) * dble(j-1) / dble(ppoints-1)
               end if

               if (temp.gt.0.d0) then
                 trans_free_E = -const_kB*temp*(3.0/2.0*log((mass*const_kB*const_eV*temp)/((const_hbar**2)*2.0*pi)) + log(const_kB*const_eV*temp/press)) -const_kB*temp
               else
                 trans_free_E = 0.d0
               end if
               if (j .eq. 1) then
                  if (i .ne. 1) then
                     write(STDOUT,'(A)') "-----------------------------------------------------------------------------------------------------------------------------------------------"
                  end if
                  write(STDOUT,'(A20,2A35,A20,A35)') 'temperature [K]', 'vibrational free energy [eV]' , 'rotational free energy [eV]', &
                      'pressure [Pa]', 'translational free energy [eV]'
                  write(STDOUT,'(F20.8,2F35.8,F20.8,F35.8)') temp, vib_free_E, rot_free_E, press, trans_free_E
               else 
                  write(STDOUT,'(90X,F20.8,F35.8)')  press, trans_free_E
               end if
            end do
          else
            write(STDOUT,'(F20.8,2F35.8)') temp, vib_free_E, rot_free_E 
          end if ! trans free energies

       end do ! loop over temperatures

    end if !print free energies



    
    ! clean up variables at the end ...   
    deallocate(frequencies)
    deallocate(atomic_positions)
    deallocate(symmetries)
    deallocate(workspace)
    deallocate(mass_vector)
    deallocate(eigenvalues)
    deallocate(hessian)
    deallocate(inputhessian)
    deallocate(atomic_masses)
    deallocate(reduced_mass)

    write(STDOUT,*) 'Leaving the Hessian diagonalizer.'

 end if ! on thread 0 ??? 

 call finalize_mpi()

 contains

 subroutine read_c

   implicit none

   integer :: i_code 
   character*132 inputline
   character*40 desc_str


   open (20, FILE="control.in")
   do
      read(20,'(A)',iostat=i_code) inputline
      if(i_code<0) exit               ! end of file reached

      read(inputline,*,iostat=i_code) desc_str
      if(i_code/=0) cycle              ! skip empty line
      if (desc_str(1:1).eq.'#') cycle  ! skip comment

      if (desc_str .eq. 'vibrations') then
         read(inputline,*,iostat=i_code) desc_str, desc_str
         if (desc_str .eq. 'free_energy') then
            read(inputline,*,iostat=i_code) desc_str, desc_str, Tstart, Tend, Tpoints
            if (i_code .eq. 0) then
               flag_free_energy = .true.
            else
               write(STDOUT,*) "WARNING: Options for tag vibrations free_energy not correctly specified!"
            end if

         else if (desc_str .eq. 'trans_free_energy') then
            read(inputline,*,iostat=i_code) desc_str, desc_str, pstart, pend, ppoints
            if (i_code .eq. 0) then
               flag_trans_free_energy = .true.
            else
               write(STDOUT,*) "WARNING: Options for tag vibrations flag_trans_free_energy not correctly specified!"
            end if

         else
            write(STDOUT,*) "WARNING: Unknown sub-tag for tag vibrations"
         end if

      end if

   end do
   close(20)
  
   if (flag_trans_free_energy) then
      if (.not. flag_free_energy) then
         write(STDOUT,'(2A)') "WARNING: Translational free energy will only be printed out" , &
         " if a temperature range is chosen using vibrations free_energy Tstart Tend Tpoints"
      end if
   end if
   
   return
 end subroutine read_c

end program hessian_diagonalization

