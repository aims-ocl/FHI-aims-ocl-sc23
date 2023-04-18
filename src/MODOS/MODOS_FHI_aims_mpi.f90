!PURPOSE: MODOS (Molecular Orbital projected Density Of States) calculation (based on FHI-aims output; MPI version)
!Author: Yong Xu, FHI, Berlin, Germany
!Send comments/suggestions/bug-reports to yongxu@fhi-berlin.mpg.de
!Reference: PhD thesis of Romaner Lorenz (TU Graz)
!
!Procedules:
!
!1. FHI-aims calculation for (1) adsorption system and (2) the isolated molecular part
!   control.in: use key word 'output eigenvec_ovlp'
!----------------------------------------------------------------------
!   this will automatically turn on some options in FHI-aims code: 
!                use_full_spectrum .true. (n_states = n_basis) 
!                collect_eigenvectors .true. if use_scalapack
!                symmetry_reduced_k_grid  .false.
!----------------------------------------------------------------------
!   geometry.in: put all the molecular geometry information after that of the substrate.
!
!2. use "comine_eigenvector_files.f90" to combine all the eigenvec.out files into one
!
!3. run MODOS, which requires input files:
!'MODOS_control.in'
!'eigenvec.out', 'ovlp_mat.out' [of calculation (1)] and 'eigenvec.out_mol' [copied from 'eigenvec.out' of calculation (2)]

!Memory requirement:
!There are four large matrices of dimension (n_basis, n_basis) plus one matrix of dimension (nE,n_states,n_spin,n_k_points) if output MODOS(E)


!Possible improvements of this program
!1. here Fermi level is set to be zero and the sigma must be the same as in FHI-aims calculation
!   determine the Fermi level position from eigenvalues for different sigma
!2. the symmetry of k_points is not considered here (symmetry_reduced_k_grid  .false.)
!   if using irreducible k_points, the weight of k_points should be considered.
!3. the number of states is required to be the same as the number of basis
!   Pay attention to the dimensions of arrays when doing corresponding changes.
!4. Loewdin charge analysis
!   Note: This concerns the square root of a complex number!
!5. solve memory problem for very large system
!   Note: (1) don't calculation MODOS
!         (2) change to distribution of the work (array: mpi_id), use more hosts and open less number of arrays for each host

!Changes:
!Jan. 02, 2012: 
!               fix an error: change occ_E(nE) into occ_E(n_states)
!               to save memory: close some large matrice (M, Minv, ctmp) (Note: KS_eigenvector and overlap matrix is not the original one any more)
!               to save memory: allocate(MODOS_m_ik, MODOS_m, Delta_E) only if(output_MODOS)
!Feb. 11, 2012:
!               Speed up the code by optimizing the most time consuming part, now the calculation of MODOS(E) costs little additional time.
!               Output the MODOS occupation into a separate file: MODOS_occ.dat
!Feb. 25, 2013:
!BB:            Replaced mpi_send/mpi_recive by mpi_allreduce, occ_m_ik and 
!               MODOS_m_ik removed

  program MODOS_FHI_aims_mpi
  use reader
  implicit none
  include 'mpif.h'
  integer mpi_rank, mpi_size, mpi_err, status(MPI_STATUS_SIZE)
  integer, allocatable :: mpi_id(:)

  !Control parameters read from the file: MODOS_control.in
  real(8) :: sigma, E0, E1
  integer :: nE
  logical :: output_MODOS, Mulliken, project_substrate_DOS, restart_CM_SM, check_norm, output_MOOP

  !adsorption system-related
  integer :: n_basis, n_states, n_atoms, n_spin, n_k_points
  integer :: n_basis2, n_spin2, n_k_points2
  integer, allocatable :: basis_atom(:)
  real*8, allocatable :: KS_eigenvalue(:)
  complex*16, allocatable :: KS_eigenvector(:,:)
  complex*16, allocatable :: overlap_matrix(:,:)

  !molecule-related
  integer :: n_basis_mol, n_states_mol, n_atoms_mol, n_spin_mol, n_k_points_mol
  integer, allocatable :: basis_atom_mol(:)
  real*8 :: dummy

  !substrate-related
  integer :: n_basis_sub, n_states_sub, n_atoms_sub, n_spin_sub, n_k_points_sub
  integer, allocatable :: basis_atom_sub(:)
  !MODOS related
  complex*16, allocatable :: CM(:,:), M(:,:), MOOP(:), MOOP_molstate(:,:)
  complex*16, allocatable :: Minv(:,:), Minv_mol(:,:), Minv_sub(:,:)
  complex*16, allocatable :: ctmp(:,:), SM(:,:)
  complex*16, allocatable :: MODOS_m(:,:)
  complex*16, allocatable :: occ_m(:)
  real*8, allocatable :: E(:), DOS(:), DOS_mol(:), DOS_sub(:), N_electrons_atom(:)
  real*8, allocatable :: occ_E(:)
  real*8  Delta_E
  real*8, external :: Gaussian, Gauss_occ
  real*8 :: dE, EE, E_Fermi
  real*8 :: N_electrons, N_electrons_mol, N_electrons_sub
  logical :: CM_file_exist, SM_file_exist, files_exist, cal_CM_SM
  
  !MOOP
  integer :: MOOP_project_basisfunc
  real*8 :: MOOP_scale
  real*8 :: MOOP_norm
  logical :: MOOP_return_Overlap, MOOP_mol_basis, MOOP_project_substrate
  real*8, allocatable :: MOOP_overlap(:)

  !files (formatted or unformatted) related
  character*40 :: filename1, filename2, filename3, filename4, filename5, filename6
  integer :: IRECL1, IRECL2, IRECL3, IRECL4, IRECL5, IRECL6
  integer :: iline1, iline2, iline3, iline4, iline5, iline6

  !counters
  integer :: i, j, k, mm, nn, l, ik, ik2
  integer :: i_state, i_spin, i_k_point, i_basis, i_basis_1, i_basis_2

  !tmp variales
  integer :: itmp, LEN
  real*8 :: rtmp
  complex*16 :: ztmp
  
  call MPI_INIT( mpi_err )
  call MPI_COMM_RANK( MPI_COMM_WORLD, mpi_rank, mpi_err )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, mpi_size, mpi_err )
!write(*,*) 'hello world! from process ', mpi_rank, ' of ',mpi_size

  INQUIRE (IOLENGTH=LEN) itmp

  !read parameters from the control file: 'MODOS_control.in'
  call parameter_reader(sigma, output_MODOS, E0, E1, nE, Mulliken, project_substrate_DOS, restart_CM_SM, check_norm, output_MOOP, MOOP_project_basisfunc, MOOP_mol_basis, MOOP_project_substrate)
  
  ! hardcoded for developement
  ! ***************************************************************
  MOOP_return_Overlap = .false.
  ! ***************************************************************

  if(mpi_rank==0)then
     
     write(*,"('sigma                 = ', es12.3)") sigma
     write(*,"('output_MODOS          = ', L)") output_MODOS
     write(*,"('output_MOOP           = ', L)") output_MOOP

     if (output_MODOS) then
        write(*,"('E0                    = ', f8.3)") E0
        write(*,"('E1                    = ', f8.3)") E1
        write(*,"('nE                    = ', I8)") nE
     end if
     
     if (output_MOOP) then
        write(*,"('MOOP_project_basisfunc = ', I)") MOOP_project_basisfunc
        write(*,"('MOOP_mol_basis  = ', L)") MOOP_mol_basis
        write(*,"('MOOP_project_substrate  = ', L)") MOOP_project_substrate
     end if
     
     write(*,"('Mulliken              = ', L)") Mulliken
     write(*,"('project_substrate_DOS = ', L)") project_substrate_DOS
     write(*,"('restart_CM_SM         = ', L)") restart_CM_SM
     write(*,"('check_norm            = ', L)") check_norm
     
     write(*,*)
  
  end if


!**********************************************************************************************
!Part 1: Read the head lines of files including eigenectors and overlap matrice
!**********************************************************************************************

!**********************************************************************************************
!Step 1: the file including eigenvectors of the adsorption system
!**********************************************************************************************
  filename1 = 'eigenvec.out'

  IRECL1 = 5*LEN
  open (UNIT=51, FILE=filename1, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL1)
  read(51, REC=1) n_basis, n_states, n_spin, n_k_points, n_atoms
  close(51)

  !output the information
  if(mpi_rank==0)then
     write(*,"('n_k_points = ', I6)") n_k_points
     write(*,"('n_spin     = ', I6)") n_spin
     write(*,*)
     write(*,"('n_basis    = ', I6)") n_basis
     write(*,"('n_states   = ', I6)") n_states
     write(*,"('n_atoms    = ', I6)") n_atoms
  end if

  !n_basis = n_states is required
  if (n_basis /= n_states) then
     if(mpi_rank==0)then
        write(*,"('Error: n_basis /= n_states')")
     end if
     stop
  end if


!**********************************************************************************************
!Step 2: the file including overlap matrix of the adsorption system
!**********************************************************************************************
  filename2 = 'ovlp_mat.out'

  IRECL2 = 3*LEN
  open (UNIT=52, FILE=filename2, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL2)
  read(52, REC=1) n_basis2, n_spin2, n_k_points2
  close(52)

  if (n_basis2 /= n_basis .or. n_spin2 /= n_spin .or. n_k_points2 /= n_k_points) then
     if(mpi_rank==0)then
        write(*, "('Error: n_basis/n_spin/n_k_points read from the two files differs!')")
     end if
     stop
  end if

  
!**********************************************************************************************
!Step 3: the file including eigenvectors of the molecular part
!**********************************************************************************************
  filename3 = 'eigenvec.out_mol'

  IRECL3 = 5*LEN
  open (UNIT=53, FILE=filename3, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL3)
  read(53, REC=1) n_basis_mol, n_states_mol, n_spin_mol, n_k_points_mol, n_atoms_mol
  close(53)

  !check: n_spin_mol, n_k_points_mol
  if (n_spin_mol /= n_spin .or. n_k_points_mol /= n_k_points) then
     if(mpi_rank==0)then
        write(*, "('Error: n_spin/n_k_points read from the two files differs!')")
     end if
     stop
  end if

  !output the information
  if(mpi_rank==0)then
     write(*,*)
     write(*,"('n_basis_mol    = ', I6)") n_basis_mol
     write(*,"('n_states_mol   = ', I6)") n_states_mol
     write(*,"('n_atoms_mol    = ', I6)") n_atoms_mol
  end if

  !n_basis = n_states is required
  if (n_basis_mol /= n_states_mol) then
     if(mpi_rank==0)then
        write(*,"('Error: n_basis_mol /= n_states_mol')")
     end if
     stop
  end if

!**********************************************************************************************
!Step 4: the file including eigenvectors of the substrate part (read it if necessary)
!**********************************************************************************************
  if (project_substrate_DOS .or. MOOP_mol_basis) then

     filename4 = 'eigenvec.out_sub'
     
     IRECL4 = 5*LEN
     open (UNIT=54, FILE=filename4, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL4)
     read(54, REC=1) n_basis_sub, n_states_sub, n_spin_sub, n_k_points_sub, n_atoms_sub
     close(54)

     !check: n_spin_sub, n_k_points_sub
     if (n_spin_sub /= n_spin .or. n_k_points_sub /= n_k_points) then
        if(mpi_rank==0)then
           write(*, "('Error: n_spin/n_k_points read from the two files differs!')")
        end if
        stop
     end if 

     !output the information
     if(mpi_rank==0)then
        write(*,*)
        write(*,"('n_basis_sub    = ', I6)") n_basis_sub
        write(*,"('n_states_sub   = ', I6)") n_states_sub
        write(*,"('n_atoms_sub    = ', I6)") n_atoms_sub
     end if

     !check: n_atoms, n_states, n_basis
     if (n_atoms_sub + n_atoms_mol /= n_atoms) then
        if(mpi_rank==0)then
           write(*, "('Error: n_atoms_sub + n_atoms_mol /= n_atoms!')")
        end if
        stop
     end if
     if (n_states_sub + n_states_mol /= n_states) then
        if(mpi_rank==0)then
           write(*, "('Error: n_states_sub + n_states_mol /= n_states!')")
        end if
        stop
     end if
     if (n_basis_sub + n_basis_mol /= n_basis) then
        if(mpi_rank==0)then
           write(*, "('Error: n_basis_sub + n_basis_mol /= n_basis!')")
        end if
        stop
     end if

  else

     n_basis_sub = n_basis - n_basis_mol
     n_states_sub = n_states - n_states_mol
     n_atoms_sub = n_atoms - n_atoms_mol
  
  end if



!**********************************************************************************************
!Part 2: Read the bulk data of files including eigenectors and overlap matrice
!**********************************************************************************************

!**********************************************************************************************
!Step 1: allocate arrays
!**********************************************************************************************
!the adsorption system:
  allocate(basis_atom(n_basis))
  allocate(KS_eigenvalue(n_states))
  allocate(overlap_matrix(n_basis, n_basis))

!the molecular system:
  allocate(basis_atom_mol(n_basis_mol))

!the substrate system:
  if (project_substrate_DOS .or. MOOP_mol_basis) then
    allocate(basis_atom_sub(n_basis_sub))
  endif

!**********************************************************************************************
!Step 2: open files
!**********************************************************************************************
  !(1) the file including eigenvectors of the adsorption system
  IRECL1 = (2 + n_basis*4)*LEN
  open (UNIT=51, FILE=filename1, STATUS='Old', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL1)  
  iline1 = 2
  read(51, REC=iline1) basis_atom(1:n_basis)
  
  !(2) the file including overlap matrix of the adsorption system
  IRECL2 = n_basis*4*LEN
  open (UNIT=52, FILE=filename2, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL2)
  iline2 = 1

  !(3) the file including eigenvectors of the molecular part
  IRECL3 = (2 + n_basis_mol*4)*LEN
  open (UNIT=53, FILE=filename3, STATUS='Old', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL3)
  iline3 = 2
  read(53, REC=iline3) basis_atom_mol(1:n_basis_mol)

  !(4) the file including eigenvectors of the substrate part (if necessary)
  if (project_substrate_DOS .or. MOOP_mol_basis) then
     
     IRECL4 = (2 + n_basis_sub*4)*LEN
     open (UNIT=54, FILE=filename4, STATUS='Old', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL4)
     iline4 = 2
     read(54, REC=iline4) basis_atom_sub(1:n_basis_sub)

  end if

!**********************************************************************************************
!Check: in geometry.in of FHI-aims the molecular part information is AFTER that of substrate part
!**********************************************************************************************
  !check: basis_atom(:) for the molecular part
  do i=1,n_basis_mol     
     if (basis_atom_mol(i) + n_atoms_sub - basis_atom(i+n_basis_sub) /= 0) then
        if(mpi_rank==0)then
           write(*,"('Error: basis_atom_mol(i) + n_atoms_sub /= basis_atom(i+n_basis_sub)')")
           write(*,"('Make sure that the molecule-related information is AFTER that of substrate part in geometry.in.')")
        end if
        stop
     end if     
  end do

  !check: basis_atom(:) for the substrate part
  if (project_substrate_DOS .or. MOOP_mol_basis) then          

     do i=1,n_basis_sub     
        if (basis_atom_sub(i) /= basis_atom(i)) then
           if(mpi_rank==0)then
              write(*,"('Error: basis_atom_sub(i) /= basis_atom(i)')")
              write(*,"('Make sure that the molecule-related information is AFTER that of substrate part in geometry.in.')")
           end if
           stop
        end if     
     end do

  end if


  !allocate arrays before calculations
  !Here we assume that n_states=n_basis.
  !Need to be careful if n_states/=n_basis (The dimensions of CM and SM dependent on whether we also project DOS on substrates orbitals.)
  !to save memory we reuse large matrice for different purposes
  allocate(CM(n_states,n_states), SM(n_states,n_states))
  allocate(Minv_mol(n_states_mol,n_basis_mol))
  
  if (project_substrate_DOS .or. MOOP_mol_basis) then   
    allocate(Minv_sub(n_states_sub,n_basis_sub))
  endif
  
  allocate(E(nE))
  allocate(occ_m(n_states))
  
  if(output_MODOS) then
    allocate(MODOS_m(nE,n_states))
  endif

  ! for the orbital overlap potential a second coefficient matrix MOOP is required
  if (output_MOOP) then
    allocate(MOOP(nE))
  endif
  
  if (MOOP_return_Overlap) then
     allocate(MOOP_overlap(n_states))
  endif
  
  !the energy points
  if (nE == 1) then
     E(1) = E0
  elseif (nE > 1) then
     dE = (E1 - E0) / dfloat(nE-1)
     do i=1,nE
        E(i) = E0 + dfloat(i-1)*dE
     end do
  else
     if(mpi_rank==0)then
        write(*,"('Error: nE < 1')")
     end if
     stop
  end if


!**********************************************************************************************
!Main part: read eigenvector and overlap matrix and do MODOS calculation for each spin and k-point
!**********************************************************************************************
  if (Mulliken) then
     if(mpi_rank==0)then
        write(*,"('Do Mulliken analysis instead of MODOS!')")
     end if
  else
     !restart option: read restart files if they exist, otherwise write restart files
     if (restart_CM_SM) then
        filename5 = "restart_CM.out"
        filename6 = "restart_SM.out"
     
        INQUIRE(FILE=filename5, EXIST=CM_file_exist)
        INQUIRE(FILE=filename6, EXIST=SM_file_exist)

        files_exist = .false.
        if (CM_file_exist .and. SM_file_exist) files_exist = .true.
     
        IRECL5 = n_states*4*LEN
        IRECL6 = n_states*4*LEN
        if (files_exist) then
           open (UNIT=55, FILE=filename5, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL5)
           open (UNIT=56, FILE=filename6, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL6)
        else
           open (UNIT=55, FILE=filename5, STATUS='REPLACE', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='WRITE', RECL=IRECL5)
           open (UNIT=56, FILE=filename6, STATUS='REPLACE', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='WRITE', RECL=IRECL6)
        end if

     end if

     !check whether we need to calculate CM/SM or not
     cal_CM_SM = .true.
     if (restart_CM_SM .and. files_exist) cal_CM_SM = .false.

     if (.not. cal_CM_SM) then
        if(mpi_rank==0)then
           write(*,"('Read CM and SM from restart files.')")
        end if
     end if

  end if


  !distribute myid for each (i_spin, i_k_point)
  allocate(mpi_id(n_k_points*n_spin))
  do i_k_point = 1, n_k_points
     do i_spin = 1, n_spin
        ik = i_spin + (i_k_point-1)*n_spin
        mpi_id(ik) = mod(ik,mpi_size)
     end do
  end do
  
  occ_m(:) = dcmplx(0.0d0)
  MODOS_m(:,:) = dcmplx(0.0d0)
  
  if (output_MOOP) then
     MOOP = dcmplx(0.0d0)
  end if
  
  write(*,"('loop over k')")
  
  do i_k_point = 1, n_k_points     
  do i_spin = 1, n_spin
     
     ik = i_spin + (i_k_point-1)*n_spin
     
     !running on different nodes
     if (mpi_rank .eq. mpi_id(ik)) then
     
        write(*,"('The calculation of i_kpoint/i_spin = ', I3, '/', I1, ' is running on process ', I4)") i_k_point, i_spin, mpi_rank
        

        !***************************************************************************************
        !Step 1: read bulk data (even for the case the restart file can be used)
        !***************************************************************************************

        !(2) overlap matrix of the adsorption system
        !***Note that only the upper part of the stored overlap_matrix contains useful data
        !the elements in the lower part: overlap_matrix(i_basis_1, i_basis_2) = 0 when i_basis_1 > i_basis_2
        do i_basis = 1, n_basis
           iline2 = 1 + i_basis + (ik-1)*n_basis
           read(52, REC=iline2) overlap_matrix(1:i_basis, i_basis)
        end do

        !Use the Hermitian property of the overlap matrix to restore the correct values for them:
        do i_basis_1 = 1, n_basis
           do i_basis_2 = 1, i_basis_1-1
              !!Use the Hermitian property of the overlap matrix to restore the correct values for them:
              overlap_matrix(i_basis_1, i_basis_2) = &
              dconjg(overlap_matrix(i_basis_2, i_basis_1))
           end do
        end do


        !***************************************************************************************
        !Step 2: calculate MOOP for different (i_spin, i_k_point)
        !***************************************************************************************
        if (output_MOOP) then
           
           allocate(KS_eigenvector(n_basis, n_states))
           
           do i_state = 1, n_states
              iline1 = 2 + i_state + (ik-1)*n_states
              read(51, REC=iline1) KS_eigenvalue(i_state), KS_eigenvector(1:n_basis, i_state)
           end do
           
           ! calculate MOOP in molecular orbital basis
           if (MOOP_mol_basis) then
              
              SM = dcmplx(0.0d0)

              !(3) eigenvectors of the molecular part
              do i_state = 1, n_states_mol
                 iline3 = 2 + i_state + (ik-1)*n_states_mol
                 ! Minv_mol ... contains eigenvectors of the molecular part
                 read(53, REC=iline3) dummy, Minv_mol(1:n_basis_mol, i_state)
              end do
              
              do i=1,n_states_mol
                 do j=1,n_basis_mol
                    ! save eigenvectors of the molecular part in SM
                    SM(j+n_basis_sub, i+n_states_sub) = Minv_mol(j, i)
                 end do
              end do
              
              !(4) eigenvectors of the substrate part
              do i_state = 1, n_states_sub
                 iline4 = 2 + i_state + (ik-1)*n_states_sub
                 ! Minv_sub ... eigenvectors of the substrate part
                 read(54, REC=iline4) dummy, Minv_sub(1:n_basis_sub, i_state)
              end do
                 
              do i=1,n_states_sub
                 do j=1,n_basis_sub
                    ! save eigenvectors of the substrate part in SM
                    SM(j, i) = Minv_sub(j, i)
                 end do
              end do

              ! #################################################################################################
              ! SM is the overlap matrix in the basis of the states of molecular and substrate part. To attain SM
              ! the overlap_matrix (S) is transfromed.
              ! SM = M^H * S * M
              ! where M is the eigenvector in the basis of the molecular and substrate part and S is the overlap_matrix.
              ! in the code:
              ! CM = M^H * S
              ! SM = M
              call zgemm('C', 'N', n_states, n_basis, n_basis, dcmplx(1.0d0), SM, n_basis, overlap_matrix, n_basis, dcmplx(0.0d0), CM, n_states)
              
              overlap_matrix = SM
              
              call zgemm('N', 'N', n_states, n_states, n_basis, dcmplx(1.0d0), CM, n_states, overlap_matrix, n_basis, dcmplx(0.0d0), SM, n_states)
              
              ! invert eigenectors of molecular and substrate part
              !Minv_mol = M_mol^-1
              call invmat_z(n_basis_mol, Minv_mol)
        
              !Minv_sub = M_sub^-1        
              call invmat_z(n_basis_sub, Minv_sub)

              ! overlap_matrix is reused to contain inverted eigenectors of molecular and substrate part -->  eigenvectors**-1
              overlap_matrix = dcmplx(0.0d0)
              do i=1,n_states_mol
                 do j=1,n_basis_mol
                    overlap_matrix(j+n_basis_sub, i+n_states_sub) = Minv_mol(j,i)
                 end do
              end do
              
              do i=1,n_states_sub
                 do j=1,n_basis_sub
                    overlap_matrix(j, i) = Minv_sub(j, i)
                 end do
              end do
              
              ! CM = M^-1 * c
              ! where M is the eigenvector in the basis of the molecular and substrate part and c the eigenvector in the adsorption basis.
              call zgemm('N', 'N', n_states, n_states, n_states, dcmplx(1.0d0), overlap_matrix, n_states, KS_eigenvector, n_states, dcmplx(0.0d0), CM, n_states)
              
              do j=1,nE
              do i=1,n_states
                 EE = E(j) - KS_eigenvalue(i)
                 !*** the sigma here could be different from the sigma used in occupation function
                 !*** if changing sigma, the Fermi level may shift away from zero.
                 Delta_E = Gaussian(sigma,EE)
                 
                 if (Delta_E > 1.0d-15) then !only consider states that have energy close to E
                    do mm=1, n_basis_sub
                       
                       ! use all eigenvectors of molecule
                       if (MOOP_project_basisfunc .eq. 0) then
                          do nn=1, n_basis_mol
                             MOOP(j) = MOOP(j) + dconjg(CM(mm, i)) * CM(n_basis_sub+nn, i) * SM(mm, n_basis_sub+nn) * Delta_E
                          end do
                       
                       ! use specific eigenvector
                       else
                          MOOP(j) = MOOP(j) + dconjg(CM(mm, i)) * CM(n_basis_sub+MOOP_project_basisfunc, i) * SM(mm, n_basis_sub+MOOP_project_basisfunc) * Delta_E
                       end if
                          
                    end do
                 end if
              end do !i
              end do !j
              
              overlap_matrix = dcmplx(0.0d0)
              
              ! read fresh overlap_matrix since it was "misused" a couple of times
              do i_basis = 1, n_basis
                 iline2 = 1 + i_basis + (ik-1)*n_basis
                 read(52, REC=iline2) overlap_matrix(1:i_basis, i_basis)
              end do

              !Use the Hermitian property of the overlap matrix to restore the correct values for them:
              do i_basis_1 = 1, n_basis
                 do i_basis_2 = 1, i_basis_1-1
                    !!Use the Hermitian property of the overlap matrix to restore the correct values for them:
                    overlap_matrix(i_basis_1, i_basis_2) = &
                    dconjg(overlap_matrix(i_basis_2, i_basis_1))
                 end do
              end do
           
           ! calculate MOOP in atom basisfunction basis
           else
           
           do j=1,nE
              do i=1,n_states
                 EE = E(j) - KS_eigenvalue(i)
                 !*** the sigma here could be different from the sigma used in occupation function
                 !*** if changing sigma, the Fermi level may shift away from zero.
                 Delta_E = Gaussian(sigma,EE)
                 
                 if (Delta_E > 1.0d-15) then !only consider states that have energy close to E

                    ! use all eigenvectors of molecule
                    if (MOOP_project_basisfunc .eq. 0) then
                       do mm=1, n_basis_sub
                          do nn=1, n_basis_mol
                             MOOP(j) = MOOP(j) + dconjg(KS_eigenvector(mm, i)) * KS_eigenvector(n_basis_sub+nn, i) * overlap_matrix(mm, n_basis_sub+nn) * Delta_E
                          end do
                       end do
                    
                    ! use specific substrate basis function
                    elseif (MOOP_project_substrate) then
                       do nn=1, n_basis_mol
                          MOOP(j) = MOOP(j) + dconjg(KS_eigenvector(MOOP_project_basisfunc, i)) * KS_eigenvector(n_basis_sub+nn, i) * overlap_matrix(MOOP_project_basisfunc, n_basis_sub+nn) * Delta_E
                       end do
                    
                    ! use specific molecule basis function
                    else
                       do mm=1, n_basis_sub
                          MOOP(j) = MOOP(j) + dconjg(KS_eigenvector(mm, i)) * KS_eigenvector(n_basis_sub+MOOP_project_basisfunc, i) * overlap_matrix(mm, n_basis_sub+MOOP_project_basisfunc) * Delta_E
                       end do
                    end if
                 end if
              end do !i
           end do !j
           end if
           
           deallocate(KS_eigenvector)
           
        end if
        
        !***************************************************************************************
        !Step 3: calculate CM and SM for different (i_spin, i_k_point)
        !***************************************************************************************
        if (Mulliken) then
           allocate(KS_eigenvector(n_basis, n_states))
           
           do i_state = 1, n_states
              iline1 = 2 + i_state + (ik-1)*n_states
              read(51, REC=iline1) KS_Eigenvalue(i_state), CM(1:n_basis, i_state)
           end do
           
           SM = overlap_matrix
                
        else
           if (cal_CM_SM) then
              !M: eigenvectors of isolated molecular and substrate parts (molecular oribtals)
              !initialization
              !M = dcmplx(0.0d0)
			  !***M --> SM
			  
              SM = dcmplx(0.0d0)

              !(3) eigenvectors of the molecular part
              do i_state = 1, n_states_mol
                 iline3 = 2 + i_state + (ik-1)*n_states_mol
                 ! Minv_mol ... contains eigenvectors of the molecular part
                 read(53, REC=iline3) dummy, Minv_mol(1:n_basis_mol, i_state)
              end do
              
              do i=1,n_states_mol
                 do j=1,n_basis_mol
                    !M(j+n_basis_sub, i+n_states_sub) = KS_eigenvector_mol(j, i)
                    ! save eigenvectors of the molecular part in SM
                    SM(j+n_basis_sub, i+n_states_sub) = Minv_mol(j, i)
                 end do
              end do
              
              !M from substrate part
              !two options: (1) use and (2) don't use orbitals of the substrate
              if (project_substrate_DOS) then
                 !(4) eigenvectors of the substrate part
                 do i_state = 1, n_states_sub
                    iline4 = 2 + i_state + (ik-1)*n_states_sub
                    ! Minv_sub ... eigenvectors of the substrate part
                    read(54, REC=iline4) dummy, Minv_sub(1:n_basis_sub, i_state)
                 end do
                 
                 do i=1,n_states_sub
                    do j=1,n_basis_sub
                       !M(j, i) = KS_eigenvector_sub(j, i)
                       ! save eigenvectors of the substrate part in SM
                       SM(j, i) = Minv_sub(j, i)
                    end do
                 end do
           
              else
                 if (n_states_sub /= n_basis_sub) then
                    write(*,"('Error: n_states_sub /= n_basis_sub')")
                    stop
                 end if

                 do i=1,n_states_sub
                    !M(i, i) = dcmplx(1.0d0)
                    
                    ! If the DOS should not be projected onto the substrate, set eigenvectors of the substrate part to 1.0.
                    SM(i, i) = dcmplx(1.0d0)
                 end do

              end if !project_substrate_DOS


              !SM = M^H * overlap_matrix * M
              !ctmp = (M^H) .x. overlap_matrix   [ctmp(n_states,n_basis), M(n_basis,n_states), overlap_matrix(n_basis,n_basis)]
              !***ctmp --> CM, M --> SM
			  !call zgemm('C', 'N', n_states, n_basis, n_basis, dcmplx(1.0d0), M, n_basis, overlap_matrix, n_basis, dcmplx(0.0d0), ctmp, n_states)
			  
              ! CM = M^H * S --> just use as temporary storage
              ! where M is the eigenvector in the basis of the molecular and substrate part and S is the overlap_matrix.
			  call zgemm('C', 'N', n_states, n_basis, n_basis, dcmplx(1.0d0), SM, n_basis, overlap_matrix, n_basis, dcmplx(0.0d0), CM, n_states)
              !SM = ctmp .x. M  [SM(n_states,n_states), ctmp(n_states,n_basis), M(n_basis,n_states)]       
              !call zgemm('N', 'N', n_states, n_states, n_basis, dcmplx(1.0d0), ctmp, n_states, M, n_basis, dcmplx(0.0d0), SM, n_states)
              !***ctmp --> CM, M --> SM --> overlap_matrix (***overlap_matrix has been changed)
              
			  overlap_matrix = SM
              
              ! #################################################################################################
              ! SM is the overlap matrix in the basis of the states of molecular and substrate part. To attain SM
              ! the overlap_matrix (S) is transfromed.
              ! SM = M^H * S * M
              ! where M is the eigenvector in the basis of the molecular and substrate part and S is the overlap_matrix.
              ! in the code:
              ! CM = M^H * S
              ! SM = M
              call zgemm('N', 'N', n_states, n_states, n_basis, dcmplx(1.0d0), CM, n_states, overlap_matrix, n_basis, dcmplx(0.0d0), SM, n_states)
              
              ! invert eigenectors of molecular and substrate part
              !Minv_mol = M_mol^-1
              call invmat_z(n_basis_mol, Minv_mol)
        
              !Minv_sub = M_sub^-1
              if (project_substrate_DOS) then          
                 call invmat_z(n_basis_sub, Minv_sub)
              end if


              !construct Minv (= M^-1) by Minv_mol and Minv_sub
			  !***Minv --> overlap_matrix
              !Minv = dcmplx(0.0d0)
              ! overlap_matrix is reused to contain inverted eigenectors of molecular and substrate part -->  eigenvectors**-1
              overlap_matrix = dcmplx(0.0d0)
              do i=1,n_states_mol
                 do j=1,n_basis_mol
                    !Minv(j+n_basis_sub, i+n_states_sub) = Minv_mol(j,i)
                    overlap_matrix(j+n_basis_sub, i+n_states_sub) = Minv_mol(j,i)
                 end do
              end do

              if (project_substrate_DOS) then 
                 do i=1,n_states_sub
                    do j=1,n_basis_sub
                       !Minv(j, i) = Minv_sub(j, i)
                       overlap_matrix(j, i) = Minv_sub(j, i)
                    end do
                 end do
                 
              else
                 do i=1,n_states_sub
                    do j=i,n_basis_sub
                       if (i.eq.j)then
                          overlap_matrix(j, i)=dcmplx(1.0d0)
                       else
                          overlap_matrix(j, i)=dcmplx(0.0d0)
                          overlap_matrix(i, j)=dcmplx(0.0d0)
                       end if
                    end do
                 end do               
              end if
              
              !(1) eigenvectors of the adsorption system
              allocate(KS_eigenvector(n_basis, n_states))
              
              do i_state = 1, n_states
                 iline1 = 2 + i_state + (ik-1)*n_states
                 read(51, REC=iline1) KS_eigenvalue(i_state), KS_eigenvector(1:n_basis, i_state)
              end do
              
              !M * CM = KS_eigenvector -> CM = M^-1 * KS_eigenvector
              !CM = Minv * KS_eigenvector   
			  !***Minv --> overlap_matrix
              !call zgemm('N', 'N', n_states, n_states, n_states, dcmplx(1.0d0), Minv, n_states, KS_eigenvector, n_states, dcmplx(0.0d0), CM, n_states)
              
              ! CM = M^-1 * c
              ! where M is the eigenvector in the basis of the molecular and substrate part and c the eigenvector in the adsorption basis.
              call zgemm('N', 'N', n_states, n_states, n_states, dcmplx(1.0d0), overlap_matrix, n_states, KS_eigenvector, n_states, dcmplx(0.0d0), CM, n_states)
              deallocate(KS_eigenvector)
       
              if (restart_CM_SM) then
                 !output CM and SM into restart files
                 do i_state = 1, n_states
                    iline5 = i_state + (ik-1)*n_states
                    write(55, REC=iline5) CM(1:n_states, i_state)
                    iline6 = i_state + (ik-1)*n_states
                    write(56, REC=iline6) SM(1:n_states, i_state)
                 end do
              end if

           else
        
              !read CM and SM from restart files
              do i_state = 1, n_states
                 iline5 = i_state + (ik-1)*n_states
                 read(55, REC=iline5) CM(1:n_states, i_state)
                 iline6 = i_state + (ik-1)*n_states
                 read(56, REC=iline6) SM(1:n_states, i_state)
              end do

           end if !cal_CM_SM

        end if !Mulliken


        !check the normalization condition of eigenvector and overlap matrix (norm = 1) (from Eqn. 3.50)
        if (check_norm) then

           do i=1,n_states
              
              ztmp = 0.0d0
              do mm=1,n_states
                 do l=1,n_states
                    ztmp = ztmp + dconjg(CM(mm,i)) * CM(l,i) * SM(mm,l)
                 end do
              end do
              
              if (abs(ztmp - 1.0d0) > 1.0d-6) then
                 write(*,"('WARNING: the normalization condition may be not satisfied!')")
                 write(*,"(3I6, 2es18.6)") i, i_spin, i_k_point, ztmp
              end if
           
           end do

        end if !check_norm


        !***************************************************************************************
        !Step 4: MODOS for different (i_spin, i_k_point)
        !***************************************************************************************
        !Fermi level is set to be zero when outputing the eigenvalues
        !***FIXME: needs the automatically change Fermi level to give correct number of electrons
        E_Fermi = 0.0d0

        !(1) occupation of mm states for (i_spin, i_k_point)
        allocate(occ_E(n_states))
        !calculate the occupation of all the states
        do i=1,n_states
           occ_E(i) = Gauss_occ(sigma,KS_eigenvalue(i)-E_Fermi)
        end do
        
        !****************The most time consuming part******************
        !use overlap_matrix(mm,i) to store data: sum_{l} (dconjg(CM(mm,i)*CM(l,i) * SM(mm,l))
        !the data can be used for both MODOS occupation and MODOS(E) calculation
        !only occupied states are considered if only calculating MODOS occupation
        
        if (output_MODOS) then !to output MODOS(E), we need the full data
           call zgemm('N', 'N', n_states, n_states, n_states, dcmplx(1.0d0), SM, n_states, CM, n_states, dcmplx(0.0d0), overlap_matrix, n_states)
           overlap_matrix = overlap_matrix * dconjg(CM)
        else
           overlap_matrix = dcmplx(0.0d0)
           do i=1,n_states
              if (occ_E(i) > 1.0d-10) then !only consider occupied states              
                 do mm=1,n_states
                    do l=1,n_states
                       overlap_matrix(mm,i) = overlap_matrix(mm,i) + CM(l,i) * SM(mm,l)      
                    end do
                    overlap_matrix(mm,i) = overlap_matrix(mm,i) * dconjg(CM(mm,i))
                 end do !mm
              end if
           end do !i
        end if !output_MODOS

        !MODOS occupation
        do i=1,n_states
           if (occ_E(i) > 1.0d-10) then !only consider occupied states
              do mm=1,n_states
                 occ_m(mm) = occ_m(mm) + occ_E(i)*overlap_matrix(mm,i) 
              end do

           end if
        end do
        deallocate(occ_E)
        !MODOS(E)
        if (output_MODOS) then
           do j=1,nE
              do i=1,n_states
                 EE = E(j) - KS_eigenvalue(i)
                 !*** the sigma here could be different from the sigma used in occupation function
                 !*** if changing sigma, the Fermi level may shift away from zero.
                 Delta_E = Gaussian(sigma,EE)  
                 if (Delta_E > 1.0d-15) then !only consider states that have energy close to E
                    do mm=1,n_states
                       MODOS_m(j,mm) = MODOS_m(j,mm) + Delta_E*overlap_matrix(mm,i)                        
                    end do
                 end if
              
              end do !i
           end do !j
        end if

     end if !mpi_rank

  end do !i_spin
  end do !i_k_point

  close(51)
  close(52)
  close(53)
  close(54)


  !send the occupation information from node mpi_id(ik) to node 0
  call sync_vector_complex(occ_m, n_states, MPI_COMM_WORLD)
     
  !send the MODOS information from node mpi_id(ik) to node 0
  if (output_MODOS) then
     call sync_vector_complex(MODOS_m, nE*n_states, MPI_COMM_WORLD)
  end if
  
  !send the MOOP information from node mpi_id(ik) to node 0
  if (output_MOOP) then
     call sync_vector_complex(MOOP, nE, MPI_COMM_WORLD)
  end if

  
  !deallocate(M, CM, SM, ctmp)
  !deallocate(Minv, Minv_mol, Minv_sub)
  deallocate(CM, SM, Minv_mol)
  !deallocate arrays
  deallocate(KS_eigenvalue, overlap_matrix)
  deallocate(basis_atom_mol)
  if (project_substrate_DOS .or. MOOP_mol_basis) then 
    deallocate(basis_atom_sub, Minv_sub)
  endif
  
!**********************************************************************************************
!Final part: result output
!**********************************************************************************************

if(mpi_rank==0)then

!**********************************************************************************************
!occupation for different orbitals
!**********************************************************************************************
  !occupation (complex*16)
  !note: divided a factor to give correct number of electrons
  occ_m = occ_m / dfloat(n_k_points*n_spin) * 2.0d0


!**********************************************************************************************
!number of electrons
!**********************************************************************************************
  !molecular part
  N_electrons_mol = 0.0d0
  do mm=1+n_states_sub,n_states_mol+n_states_sub
     N_electrons_mol = N_electrons_mol + DBLE(occ_m(mm)) 
  end do
  
  !substrate part
  N_electrons_sub = 0.0d0
  do mm=1,n_states_sub
     N_electrons_sub = N_electrons_sub + DBLE(occ_m(mm)) 
  end do

  N_electrons = N_electrons_mol + N_electrons_sub

  !output the information
  write(*,*)
  write(*,"('N_electrons_mol   = ', f12.6)") N_electrons_mol
  write(*,"('N_electrons_sub   = ', f12.6)") N_electrons_sub
  write(*,"('N_electrons       = ', f12.6)") N_electrons

  !number of electrons for each atom if doing Mulliken analysis
  if (Mulliken) then

     allocate(N_electrons_atom(n_atoms))
     N_electrons_atom = 0.0d0
     do mm=1,n_states
        N_electrons_atom(basis_atom(mm)) = N_electrons_atom(basis_atom(mm)) + DBLE(occ_m(mm))        
     end do

     write(*,*)
     write(*,*)"Number of electrons for each atom"
     write(*,"(A6, A18)")'i_atom', 'N_electrons'

     do i=1,n_atoms
        write(*,"(I6, f18.6)")i, N_electrons_atom(i)
     end do

     deallocate(N_electrons_atom)

  else if (.not. project_substrate_DOS) then

     !number of electrons for each atom of the substrate
     allocate(N_electrons_atom(n_atoms_sub))
     N_electrons_atom = 0.0d0
     do mm=1,n_states_sub
        N_electrons_atom(basis_atom(mm)) = N_electrons_atom(basis_atom(mm)) + DBLE(occ_m(mm))     
     end do

     write(*,*)
     write(*,*)"Number of electrons for each atom belongs to substrate"
     write(*,"(A6, A18)")'i_atom', 'N_electrons'

     do i=1,n_atoms_sub
        write(*,"(I6, f18.6)")i, N_electrons_atom(i)
     end do
     
     !project DOS only onto molecular orbitals
     write(*,*)
     write(*,*)"Atomic orbital occupation for substrate:"
     write(*,"(A6,A18, A10)")'state', 'occupation', 'i_atom'
     do mm=1,n_states_sub
        write(*,"(I6, f18.6, I10)")mm, DBLE(occ_m(mm)), basis_atom(mm)
     end do

     write(*,*)
     write(*,*)"Molecular orbital occupation:"
     write(*,"(A6, A18)")'state', 'occupation'
     do mm=1,n_states_mol
        write(*,"(I6, f18.6)")mm, DBLE(occ_m(mm+n_states_sub))
     end do

     !output the MODOS occupation into 'MODOS_occ.dat'
     open (UNIT=40, FILE='MODOS_occ.dat', STATUS='REPLACE')
     write(40,"(A6, A18)")'MO', 'MODOS-occupation'
     do mm=1,n_states_mol
        write(40,"(I6, f18.6)")mm, DBLE(occ_m(mm+n_states_sub))
     end do
     close(40)

     deallocate(N_electrons_atom)
  
  else
  
     !project DOS also onto substrate's orbitals
     write(*,*)
     write(*,*)"Substrate-orbital occupation:"
     write(*,"(A6, A18)")'state', 'occupation'
     do mm=1,n_states_sub
        write(*,"(I6, f18.6)")mm, DBLE(occ_m(mm))
     end do     

     write(*,*)
     write(*,*)"Molecular orbital occupation:"
     write(*,"(A6, A18)")'state', 'occupation'
     do mm=1,n_states_mol
        write(*,"(I6, f18.6)")mm, DBLE(occ_m(mm+n_states_sub))
     end do     

  end if
  deallocate(basis_atom)


!**********************************************************************************************
! wirte orbital overlap potential to file
!**********************************************************************************************
  if (output_MOOP) then
    
     MOOP = MOOP / dfloat(n_k_points*n_spin)
    
     open (UNIT=57, FILE='MOOP.dat', STATUS='REPLACE')
     write(57,"(A12, A18)")'Energy', 'MOOP' 
     do j=1,nE
        write(57,"(f12.6, f18.6)") E(j), DBLE(MOOP(j))
     end do
     close(57)
     
     deallocate(MOOP)
  endif

  if (MOOP_return_Overlap) then
  
     open (UNIT=58, FILE='MOOP_overlap.dat', STATUS='REPLACE')
     write(58,"(A18)")'MOOP_overlap' 
     do j=1,n_states
        write(58,"(f18.6)") DBLE( MOOP_overlap(j) )
     end do
     close(58)
  
     deallocate(MOOP_overlap)
  endif
  
!**********************************************************************************************
!MODOS related
!**********************************************************************************************
  if (output_MODOS) then

     MODOS_m = MODOS_m / dfloat(n_k_points*n_spin) * 2.0d0

     allocate(DOS(nE), DOS_mol(nE), DOS_sub(nE))
     !DOS
     DOS_mol = 0.0d0
     DOS_sub = 0.0d0
     do j=1,nE
        do mm=1+n_states_sub,n_states_mol+n_states_sub
           DOS_mol(j) = DOS_mol(j) + DBLE(MODOS_m(j,mm)) 
        end do
     
        do mm=1,n_states_sub
           DOS_sub(j) = DOS_sub(j) + DBLE(MODOS_m(j,mm))
        end do
     end do
  
     DOS = DOS_mol + DOS_sub
     
     !output DOS
     open (UNIT=50, FILE='MODOS.dat', STATUS='REPLACE')
     write(50,"(A12, 3A18)")'Energy', 'DOS_total', 'DOS_mol' 
     do j=1,nE
        write(50,"(f12.6, 2f18.6, 1000f18.6)") E(j), DOS(j), DOS_mol(j), DBLE(MODOS_m(j,1+n_states_sub:n_states_mol+n_states_sub))
     end do
     close(50)
     deallocate(DOS, DOS_mol, DOS_sub)
  end if !output_MODOS

end if

  if(output_MODOS) deallocate( MODOS_m)
  deallocate( occ_m)
  deallocate(E)
  deallocate(mpi_id)

  call MPI_FINALIZE(mpi_err)

  stop
  end program MODOS_FHI_aims_mpi





  !=========================================================================== 
  !Inverse a complex*16 matrix
  !=========================================================================== 
  subroutine invmat_z(n, A)
  implicit none
  integer, intent(in) :: n
  complex*16, intent(inout) :: A(n,n)
  integer :: lwork, info
  integer, allocatable :: ipiv(:)
  complex*16, allocatable :: work(:)
      
  !call zgetrf ( m, n, a, lda, ipiv, info )
  !call zgetri (n, a, lda, ipiv, work, lwork, info)
  lwork = n
  allocate(ipiv(n), work(lwork))
      
  call zgetrf(n, n, A, n, ipiv, info )
  if (info/=0) then
     write(*,"('Error: wrong in zgetrf',I8)")info
     stop
  end if
      
  call zgetri (n, a, n, ipiv, work, lwork, info)
  if (info/=0) then
     write(*,"('Error: wrong in zgetri',I8)")info
     stop
  end if
      
  deallocate(ipiv, work)
      
  return
  end subroutine invmat_z  
    
  !BB: subroutine sync_vector_complex copied from synchronize_mpi_basic.f90 of
  !    FHI-aims
  subroutine sync_vector_complex(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary double complex vector over all tasks in
    !    a given MPI-communicator in manageable parts.
    implicit none
    include 'mpif.h'
    !  ARGUMENTS
    integer:: dim
    complex*16 :: vector(dim)
    integer, optional :: mpi_comm
    !  INPUTS
    !    o dim -- dimension of the array to be synchronized
    !    o vector -- the array itself
    !    o mpi_comm -- the communicator over which the synchronization takes 
    !                  place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks in 
    !                mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    complex*16, allocatable, dimension(:) :: temp_mpi

    integer:: comm, i, len
    integer :: mpierr

    !     adjust if needed
    integer, parameter :: max_len = 500000


     comm = mpi_comm



    allocate(temp_mpi(min(dim,max_len)),stat=i)

    do i=1,dim,max_len

       if(dim-i+1 < max_len) then
          len = dim-i+1
       else
          len = max_len
       endif

       call MPI_ALLREDUCE(vector(i), temp_mpi, len, &
       MPI_DOUBLE_COMPLEX, MPI_SUM, comm, mpierr)
       vector(i:i-1+len) = temp_mpi(1:len)

    enddo

    deallocate(temp_mpi)

  end subroutine sync_vector_complex

  !http://mathworld.wolfram.com/GaussianFunction.html
  !Gaussian broadening function
  real*8 function Gaussian(sigma,E)
  implicit none
  real*8, intent(in) :: sigma, E
  real*8 ::  PI, f0

  PI = datan(1.0d0)*4.0d0
  
  f0 = 1.0d0/sigma/dsqrt(2.0d0*PI)

  Gaussian = f0 * dexp(-(E*E)/(2.0d0*sigma*sigma))

  end function Gaussian



  !Gaussian occupation function: intergration of Gassian function
  real*8 function Gauss_occ(sigma,E)
  implicit none
  real*8, intent(in) :: sigma, E
!  real*8, external :: derf
  real*8 :: derf
  real*8, parameter :: zero = 1.0d-10

  if (sigma > 0.0d0) then
     Gauss_occ = 0.5d0 * (1.0d0 - derf(E/sigma))
  elseif (sigma < zero .and. sigma >= 0.0d0) then
     if (E > 0.0d0) then
        Gauss_occ = 0.0d0
     else
        Gauss_occ = 1.0d0  !Note: without broadening E = 0 may cause problem!!!
     end if  
  else
     write(*,"('Error: sigma < 0')")
     stop
  end if

  end function Gauss_occ


  
