!****** FHI-aims/friction
module friction
!  NAME
!    friction 
!  SYNOPSIS
!  AUTHOR
!   Reinhard J. Maurer, Yale University (2015)

!  PURPOSE
!  calculates nonadiabatic coupling in the 
!  ground state Kohn-Sham single particle picture 
!  and constructs the electron-hole pair Friction tensor
!  USES

    use dimensions
    use runtime_choices
    use pbc_lists
    use timing
    use physics
    use species_data
    use localorb_io
    use geometry
    use constants
    use mpi_tasks
    use scalapack_wrapper, only: my_k_point, eigenvec, eigenvec_complex, &
                           l_col, l_row , mxld, mxcol, &
                           ovlp, ovlp_complex, ham, ham_complex
    use synchronize_mpi_basic
    use debugmanager, only: module_is_debugged, debugprint
    use basis

    implicit none
    private

    !TODO include more check_allocation statements
    !TODO DFPT response
    !TODO get SCALAPACK to work
    !TODO   - 1 define scalapack distributed H1, S1, G1, and G2
    !TODO   - 2 do partial <psi|G|psi> and synchronize over k and bands

    !TODO IN THE CASE OF .not. use_local_index scalapack is only used 
    !TODO to diagonalize, which means we have hamiltonian and overlap e
    !TODO anyway in full, the question is what to do now

    !TODO do we implement scalapack in a way that every node gets a piece 
    !TODO of the matrices and we push this trough the code, or do we 
    !TODO just take the right kpoint parall. the first would be more efficient
    !for Nk< Ncpus

    !TODO currently only gaussian broadening works without explicit friction_spectrum

    !------------------------------------------------------!
    !            PRIVATE VARIABLES                         !
    !------------------------------------------------------!
    !real*8,dimension(:,:,:,:,:,:),allocatable :: first_order_H
    !complex*16,dimension(:,:,:,:,:,:),allocatable :: first_order_H_cmplx
    real*8,dimension(:,:,:,:,:,:),allocatable :: first_order_G
    complex*16,dimension(:,:,:,:,:,:),allocatable :: first_order_G_cmplx
    real*8,dimension(:,:,:,:,:),allocatable :: first_order_S
    complex*16,dimension(:,:,:,:,:),allocatable :: first_order_S_cmplx
    !SPARSE
    !real*8,dimension(:,:,:),allocatable :: first_order_S
    !real*8,dimension(:,:,:,:),allocatable :: first_order_H
    
    real*8, dimension(:,:), allocatable :: coords_backup

    real*8, dimension(:,:,:,:), allocatable :: KS_eigenvector_backup
    complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_complex_backup
    real*8, dimension(:,:,:),   allocatable :: occ_numbers_backup
    real*8, dimension(:,:,:),   allocatable :: KS_eigenvalue_backup
    
    
    real*8, dimension(:), allocatable   :: masses
    
    !------------------------------------------------------!
    !            PRIVATE CONSTANTS                         !
    !------------------------------------------------------!
   
    !real*8, parameter :: second = 2.418884326505d−17
    real*8, parameter :: ps = 2.418884326505d-5
    !------------------------------------------------------!
    !            PUBLIC VARIABLES                          !
    !------------------------------------------------------!
    !TODO calculate friction tensor as a function of q
    complex*16, dimension(:,:,:),allocatable, public :: friction_tensor 
    complex*16, dimension(:,:,:),allocatable, public :: friction_eigvecs 
    real*8, dimension(:,:), allocatable, public :: friction_eigenvalues 
    !complex*16, dimension(:,:),allocatable, public :: friction_tensor 
    !complex*16, dimension(:,:),allocatable, public :: friction_eigvecs 
    !real*8, dimension(:), allocatable, public :: friction_eigenvalues 
    integer, public :: friction_n_active_atoms = 0
   
    integer, public :: friction_iter_limit = 20
    real*8, public :: friction_accuracy_etot=0.00001/hartree
    real*8, public :: friction_accuracy_eev=0.01/hartree
    real*8, public :: friction_accuracy_rho=0.00001
    real*8, public :: friction_accuracy_potjump=1.01/hartree

    real*8, public :: friction_temperature = 300.0

    logical, public :: numerical_friction = .true.
    logical, public :: friction_active_atoms = .false.
    real*8, public :: friction_numeric_disp = 0.0025*bohr !should be 0.0025
    logical, public :: output_first_order_matrices = .false.
    logical, public :: output_friction_eigenvectors = .false.


    logical, dimension(:), allocatable, public :: friction_atoms_list
    integer, dimension(:), allocatable, public :: friction_index_list
    logical, public :: friction_knotk=.false.
    logical, public :: friction_read_matrices=.false.
    
    logical, public ::  friction_output_spectrum = .false. 

    !calculate friction from H1 and S1
    character(len=20), public :: friction_delta_type='gaussian' !is one of square, gaussian, lorentzian, sine
    real*8, public :: friction_window_size = 0.010/hartree
    real*8, public :: friction_broadening_width = 0.30/hartree
    real*8, public :: friction_perturbation = 0.d0
    
    real*8, public :: friction_max_energy = 2.00/hartree
    real*8, public :: friction_discretization_length = 0.01/hartree

    logical, public :: friction_use_complex_matrices = .true.
    !logical, public :: friction_use_complex_matrices = .false.

    integer, public :: friction_coupling_matrix_mode = 0
    integer, public :: friction_n_q_points = 1

    !------------------------------------------------------!
    !            INTERFACES                                !
    !------------------------------------------------------!

    interface friction_read_matrix
      module procedure friction_read_matrices_H_S
      module procedure friction_read_matrix_G
    end interface

    interface friction_read_matrix_p1
      module procedure friction_read_matrices_H_S_p1
      module procedure friction_read_matrix_G_p1
    end interface
    
    interface friction_output_matrix
      module procedure friction_output_matrices_H_S
      module procedure friction_output_matrix_G 
    end interface
    
    interface friction_output_matrix_p1
      module procedure friction_output_matrices_H_S_p1 
      module procedure friction_output_matrix_G_p1

    end interface

    !------------------------------------------------------!
    !            PUBLIC  ROUTINES                          !
    !------------------------------------------------------!
    public :: friction_calculation !main entry point to module 
    public :: allocate_friction
    public :: deallocate_friction

    !PRIVATE ROUTINES

    !allocate_friction
    !deallocate_friction
    !allocate_backup_physics
    !friction_calculation
    !friction_calculate_tensor
    !friction_calculate_tensor_p1
    !friction_calculate_S1_DFPT
    !friction_calculate_S1_DFPT_p1
    !friction_calculate_H1_and_S1_num 
    !friction_calculate_H1_and_S1_DFPT
    !friction_construct_complex_matrix
    !friction_construct_real_matrix
    !friction_output_matrices 
    !friction_output_matrices_p1 
    !friction_read_matrices
    !friction_read_matrices_p1

contains
    subroutine allocate_friction()
    !  NAME
    !    allocate_friction
    !  SYNOPSIS
    !  allocates important dynamical arrays
    !
    !    AJL: If you add a new allocation *PLEASE PLEASE* include a
    !    deallocation in the cleanup routine, otherwise problems are
    !    encountered when treating aims as a library
    !
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)

        implicit none
       
        integer :: info 
        character(*), parameter :: func = 'allocate_friction'

        allocate(friction_atoms_list(n_atoms),stat=info)
        friction_atoms_list = .false.
        call check_allocation(info, "friction_atoms_list",func)

        !IF FRICTION_KNOTK =.true., 
        if (friction_knotk) then
          friction_use_complex_matrices = .false.
        endif

    end subroutine allocate_friction
    
    subroutine allocate_backup_physics()
    !  NAME
    !    allocate_backup_physics
    !  SYNOPSIS
    !  allocates important dynamical arrays
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)
     
        implicit none
       
        integer :: info 
        character(*), parameter :: func = 'allocate_backup_physics'

        !TODO fix all possible allocation options for SCALAPACK
        !allocate H1 and S1
        !if (use_scalapack) then
          !if (n_periodic==0.or.real_eigenvectors) then
            !allocate(first_order_S(3,friction_n_active_atoms,1,mxld,mxcol))
            !allocate(first_order_H(3,friction_n_active_atoms,1,mxld,mxcol,n_spin))
          !else
            !allocate(first_order_S_cmplx(3,friction_n_active_atoms,1,mxld,mxcol))
            !allocate(first_order_H_cmplx(3,friction_n_active_atoms,1,mxld,mxcol,n_spin))
          !endif 
        !else

        !endif 
              
        if (friction_use_complex_matrices) then
          if (n_periodic==0) then
            !allocate(first_order_G(3,friction_n_active_atoms,1,n_basis,n_basis,n_spin))
            if (.not. allocated(first_order_G)) then
                allocate(first_order_G(3,friction_n_active_atoms,1,n_basis*(n_basis+1)/2,1,n_spin),stat=info)
                call check_allocation(info, 'first_order_G        ')
            endif
            if (.not. allocated(first_order_S)) then
                allocate(first_order_S(3,friction_n_active_atoms,1,n_basis*(n_basis+1)/2,1),stat=info)
                call check_allocation(info, 'first_order_S        ')
            endif
          else
            if (real_eigenvectors) then
              !allocate(first_order_G(3,friction_n_active_atoms,n_k_points_task,n_basis,n_basis,n_spin))
              if (.not. allocated(first_order_G)) then
                  allocate(first_order_G(3,friction_n_active_atoms,n_k_points_task,n_basis*(n_basis+1)/2,1,n_spin),stat=info)
                  call check_allocation(info, 'first_order_G        ')
              endif
              if (.not. allocated(first_order_S)) then
                  allocate(first_order_S(3,friction_n_active_atoms,n_k_points_task,n_basis*(n_basis+1)/2,1),stat=info)
                  call check_allocation(info, 'first_order_S        ')
              endif
            else
              if (.not. allocated(first_order_G_cmplx)) then
                  allocate(first_order_G_cmplx(3,friction_n_active_atoms,n_k_points_task,n_basis*(n_basis+1)/2,1,n_spin),stat=info)
                  call check_allocation(info, 'first_order_G_cmplx        ')
              endif
              if (.not. allocated(first_order_S_cmplx)) then
                  allocate(first_order_S_cmplx(3,friction_n_active_atoms,n_k_points_task,n_basis*(n_basis+1)/2,1),stat=info)
                  call check_allocation(info, 'first_order_S_cmplx        ')
              endif
            endif
          endif
        else
          if (n_periodic==0) then
            !allocate(first_order_G(3,friction_n_active_atoms,1,n_basis,n_basis,n_spin))
            if (.not. allocated(first_order_G)) then
                allocate(first_order_G(3,friction_n_active_atoms,1,n_basis*(n_basis+1)/2,1,n_spin),stat=info)
                call check_allocation(info, 'first_order_G        ')
            endif
            if (.not. allocated(first_order_S)) then
                allocate(first_order_S(3,friction_n_active_atoms,1,n_basis*(n_basis+1)/2,1),stat=info)
                call check_allocation(info, 'first_order_S        ')
            endif
          else
              !allocate(first_order_G(3,friction_n_active_atoms,n_cells_in_hamiltonian,n_basis,n_basis,n_spin))
              if (.not. allocated(first_order_G)) then
                  allocate(first_order_G(3,friction_n_active_atoms,n_cells_in_hamiltonian,n_basis*(n_basis+1)/2,1,n_spin),stat=info)
                  call check_allocation(info, 'first_order_G        ')
              endif
              if (.not. allocated(first_order_S)) then
                  allocate(first_order_S(3,friction_n_active_atoms,n_cells_in_hamiltonian,n_basis*(n_basis+1)/2,1),stat=info)
                  call check_allocation(info, 'first_order_S        ')
              endif
          endif
        endif
        
        if (.not.allocated(first_order_G))allocate(first_order_G(1,1,1,1,1,1)) ! allocate dummies
        if (.not.allocated(first_order_G_cmplx))allocate(first_order_G_cmplx(1,1,1,1,1,1)) ! allocate dummies
        if (.not.allocated(first_order_S))allocate(first_order_S(1,1,1,1,1)) ! allocate dummies
        if (.not.allocated(first_order_S_cmplx))allocate(first_order_S_cmplx(1,1,1,1,1)) ! allocate dummies
       
        if (.not. allocated(friction_index_list)) then
            allocate(friction_index_list(friction_n_active_atoms),stat=info)
            call check_allocation(info, "friction_index_list",func)
        endif
        if (.not. allocated(masses)) then
            allocate(masses(friction_n_active_atoms),stat=info)
            call check_allocation(info, "masses",func)
        endif

    end subroutine allocate_backup_physics

    subroutine deallocate_friction()
    !  NAME
    !    deallocate_friction 
    !  SYNOPSIS
    !  deallocates important dynamical arrays
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)

        implicit none
        
        character(*), parameter :: func = 'deallocate_friction'

        if (allocated(first_order_G)) deallocate(first_order_G)
        if (allocated(first_order_G_cmplx)) deallocate(first_order_G_cmplx)
        if (allocated(first_order_S)) deallocate(first_order_S)
        if (allocated(first_order_S_cmplx)) deallocate(first_order_S_cmplx)
       
        if (allocated(friction_index_list)) deallocate(friction_index_list)
        if (allocated(masses)) deallocate(masses)

        if (allocated(friction_atoms_list)) deallocate(friction_atoms_list)

    end subroutine deallocate_friction

    subroutine friction_calculation(converged_scf) 
    !  NAME
    !    friction_calculation 
    !  SYNOPSIS
    !    main entrance routine that manages the calculation of 
    !    the friction tensor. distributes work to subroutines 
    !    for isolated, periodic cases 
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)

        implicit none

        logical, intent(in) :: converged_scf

        character(len=300) :: info_str
        character(len=300) :: tmp_str

        integer :: i_atom, i_atom2, i_coord, j_coord
        integer :: n_k_points2,i_q_point
        integer :: success
        integer :: pstart, pend, n_print, p
        real*8, dimension(:), allocatable :: work
        complex*16, dimension(:), allocatable :: work_cmplx
        
        character(*), parameter :: func = 'friction_calculation'


        if (myid.eq.0) then
            write (info_str,'(2X,A,2X)') &
            "************************FRICTION**********************************"
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "          _____     _      _   _               "
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "         |  ___| __(_) ___| |_(_) ___  _ __    "
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "         | |_ | '__| |/ __| __| |/ _ \| '_ \   "
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "         |  _|| |  | | (__| |_| | (_) | | | |  "
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "         |_|  |_|  |_|\___|\__|_|\___/|_| |_|  "
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "************************FRICTION**********************************"
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "+-+-+-+-+   Please cite Phys. Rev. B  94, 115432 (2016)  +-+-+-+-+"
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A,2X)') &
            "************************FRICTION**********************************"
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(4X,A,2X)') &
            "FRICTION Calculating cartesian friction tensor due &
            & to electron-phonon coupling"
            call localorb_info(info_str,use_unit,'(A)')
        endif
        
        if (.not. converged_scf) then
            if (myid.eq.0) then
                write (info_str,'(4X,A)') &
                    "FRICTION Cannot calculate friction if scf is not converged. &
                    & Skipping Friction!"
                call localorb_info(info_str,use_unit,'(A)')
            endif
            return
        endif 
        if (friction_n_active_atoms <1) then
            if (myid.eq.0) then
                write (info_str,'(4X,A)') &
                    "No atoms were included in the friction calculation. &
                    & You need to set calculate_friction .true. for atoms in geometry.in! Goodbye!"
                call localorb_info(info_str,use_unit,'(A)')
            endif
            return
        endif
        
        !!!AT THE MOMENT q-integration doesnt work, SO WE SKIP IT
        !TODO will only work for the DFPT case
        friction_knotk = .false.
        !!!AT THE MOMENT q-integration doesnt work, SO WE SKIP IT
        if (friction_knotk) then
            !we need to work with real matrices
            friction_use_complex_matrices = .false.
        endif
        !if (n_periodic.eq.0) then
        !   friction_use_complex_matrices = .true.
        !endif

        if (friction_coupling_matrix_mode==2 .and. numerical_friction) then
            write (info_str,'(4X,A)') &
                "Exact coupling matrix elements (mode=2) can only be used in combination with DFPT. &
                 & Please evaluate electron-phonon coupling via DFPT."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop('friction_coupling_matrix_mode=2 needs friction_calculation DFPT!')
        endif

        if (.not. numerical_friction .and. n_periodic>0) then
            write(info_str,'(4x,A)') &
                "DFPT currently only works for the cluster case and currently only calculates first_order_S. &
                & Switching back to numerical friction!"
            call localorb_info(info_str,use_unit,'(A)')
            numerical_friction = .true.
        endif

        call allocate_backup_physics()

        friction_index_list = 0
        i_atom2 = 0
        do i_atom=1, n_atoms
            if (friction_atoms_list(i_atom)) then
                i_atom2 = i_atom2 +1
                masses(i_atom2) = species_m(species(i_atom))*mass_atomic_unit
                friction_index_list(i_atom2) = i_atom
            endif
        enddo

        if (myid.eq.0) then
            write (info_str,'(4X,A,I8,A)') &
                "FRICTION Calculating cartesian friction for ",friction_n_active_atoms," atoms"
            call localorb_info(info_str,use_unit,'(A)')
        endif
        if (myid.eq.0 .and. module_is_debugged('friction')) then
            write (info_str,'(4X,100F16.5)') &
                (masses(i_atom), i_atom=1, friction_n_active_atoms)
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(4X,A,I8)') &
                "Maximum number of SCF or DFPT iterations for response is ",friction_iter_limit
            call localorb_info(info_str,use_unit,'(A)')
            if (numerical_friction) then
                write (info_str,'(4X,A)') &
                    "FRICTION First order H and S are calculated via finite differences"
                call localorb_info(info_str,use_unit,'(A)')
            else
                write (info_str,'(4X,A)') &
                    "FRICTION First order H and S are calculated via DFPT"
                call localorb_info(info_str,use_unit,'(A)')
            end if
        endif
        
        if (friction_knotk) then
            allocate(friction_tensor(3*friction_n_active_atoms,&
                3*friction_n_active_atoms,n_k_points))
            allocate(friction_eigvecs(3*friction_n_active_atoms,&
                3*friction_n_active_atoms,n_k_points))
            allocate(friction_eigenvalues(3*friction_n_active_atoms,n_k_points))
        else
            allocate(friction_tensor(3*friction_n_active_atoms,&
                3*friction_n_active_atoms,1))
            allocate(friction_eigvecs(3*friction_n_active_atoms,&
                3*friction_n_active_atoms,1))
            allocate(friction_eigenvalues(3*friction_n_active_atoms,1))
        endif
        friction_tensor = 0.0
        friction_eigvecs = 0.0
        friction_eigenvalues = 0.0

        if (friction_knotk) then
            !pass
        else
            friction_n_q_points = 1
        endif

        !DFPT only works for n_spin=1 and xc_func = LDA
        if (.not. numerical_friction .and.( flag_xc.ne.3 .or. spin_treatment.ne.0)) then
          write (info_str,'(4X,A,I8,A)') &
              "DFPT Friction currently only works for pz-lda and n_spin=1."
          call localorb_info(info_str,use_unit,'(A)')
          write (info_str,'(4X,A,I8,A)') &
              "Switching to numeric finite-difference based evaluation!"
          call localorb_info(info_str,use_unit,'(A)')
          numerical_friction = .true.
        endif

        !CALCULATING TENSOR
        if (n_periodic.eq.0) then
            if (myid.eq.0) then
                write (info_str,'(4X,A,I8,A)') &
                    "Friction will be calculated in isolated system"
                call localorb_info(info_str,use_unit,'(A)')
            endif
            call friction_calculate_tensor()
        elseif (n_periodic>0) then
                write (info_str,'(4X,A,I8,A)') &
                    "Friction will be calculated in periodic system"
                call localorb_info(info_str,use_unit,'(A)')
            call friction_calculate_tensor_p1()
        endif

        !friction tensor is now set
        !now use the tensor, diagonalize analyze, etc.
        !transform friction to SI units

        !!zero out imaginary parts of friction tensor
        !do i_coord=1, 3*friction_n_active_atoms
          !do j_coord=1, 3*friction_n_active_atoms
            !friction_tensor(i_coord,j_coord,:) = cmplx(real(friction_tensor(i_coord,j_coord,:)),0.d0)  
          !enddo
        !enddo

        do i_q_point =1, friction_n_q_points, 1 !n_k_points2, 1

        if (friction_knotk) then
            !if(myid.eq.0 .and. module_is_debugged('friction')) then
            if(myid.eq.0) then
                write(info_str,'(A,3F12.6)') 'friction tensor at q= ',&
                    k_point_list(i_q_point,1),k_point_list(i_q_point,2),k_point_list(i_q_point,3)
                call localorb_info(info_str,use_unit,'(A)')
            endif
        endif

        n_print = friction_n_active_atoms/2+mod(friction_n_active_atoms,2)
        if(myid.eq.0 ) then
            write (info_str,'(4X,A)') &
                "********Printing Friction Tensor in 1/ps*********"
            call localorb_info(info_str,use_unit,'(A)')
            do i_coord=1, 3*friction_n_active_atoms
                write (info_str,'(4X,I6)') i_coord
                do p=1, n_print
                  pstart = 1+(p-1)*6 
                  pend = pstart+5
                  if (pend>(3*friction_n_active_atoms)) pend=3*friction_n_active_atoms 
                  if (p==1) then
                      write (info_str,'(4X,I6,6E18.9)') i_coord, &
                          (dble(friction_tensor(i_coord,j_coord,i_q_point))/ps, j_coord=pstart, pend)
                      call localorb_info(info_str,use_unit,'(A)')
                  else
                      write (info_str,'(10X,6E18.9)') &
                          (dble(friction_tensor(i_coord,j_coord,i_q_point))/ps, j_coord=pstart, pend)
                      call localorb_info(info_str,use_unit,'(A)')
                  endif
                enddo
            enddo
            write (info_str,'(4X,A)') &
                "**********END Printing Friction Tensor************"
            call localorb_info(info_str,use_unit,'(A)')
            ! Pass tensor to i-PI if asked
            if (use_pimd_wrapper .and. ipi_ftensor) then
                write(tmp_str,'(A, 1X, 300I5)') "Friction indexes:", (friction_index_list(p), p=1, friction_n_active_atoms)
                comm_string=trim(comm_string) // trim(tmp_str)
                write(tmp_str,'(1X, A, 1X)') "Friction:"
                comm_string=trim(comm_string) // trim(tmp_str)
                do i_coord=1, 3*friction_n_active_atoms
                    write(tmp_str,'(300E20.6E3)') (dble(friction_tensor(i_coord,j_coord,i_q_point)), j_coord=1, 3*friction_n_active_atoms)
                    comm_string=trim(comm_string) // trim(tmp_str)
                enddo           
            endif
        endif
        if(myid.eq.0 ) then
            write (info_str,'(4X,A)') &
                "**********Printing Lifetime Tensor in ps**********"
            call localorb_info(info_str,use_unit,'(A)')
            do i_coord=1, 3*friction_n_active_atoms
                write (info_str,'(4X,I6)') i_coord
                do p=1, n_print
                  pstart = 1+(p-1)*6 
                  pend = pstart+5
                  if (pend>3*friction_n_active_atoms) pend=3*friction_n_active_atoms 
                  if (p==1) then
                      write (info_str,'(4X,I6,6E18.9)') i_coord, &
                          (ps/dble(friction_tensor(i_coord,j_coord,i_q_point)), j_coord=pstart, pend)
                      call localorb_info(info_str,use_unit,'(A)')
                  else
                      write (info_str,'(10X,6E18.9)') &
                          (ps/dble(friction_tensor(i_coord,j_coord,i_q_point)), j_coord=pstart, pend)
                      call localorb_info(info_str,use_unit,'(A)')
                  endif
                enddo
            enddo
            write (info_str,'(4X,A)') &
                "**********END Printing Lifetime Tensor************"
            call localorb_info(info_str,use_unit,'(A)')
        endif
       
        allocate(work_cmplx(3*3*friction_n_active_atoms)) 
        allocate(work(3*3*friction_n_active_atoms-2)) 
        !diagonalize friction tensor
        !friction_eigvecs(:,:,i_k_point2) = friction_tensor(:,:,i_k_point2)
        friction_eigvecs(:,:,i_q_point) = real(friction_tensor(:,:,i_q_point))
        !call dsyev('V','U',3*friction_n_active_atoms,friction_eigvecs,&
        !    3*friction_n_active_atoms,friction_eigenvalues,work,&
        !    3*3*friction_n_active_atoms,success)
        call zheev('V','U',3*friction_n_active_atoms,friction_eigvecs(:,:,i_q_point),&
            3*friction_n_active_atoms,friction_eigenvalues(:,i_q_point),work_cmplx,&
            3*3*friction_n_active_atoms,work,success)
        deallocate(work)
        deallocate(work_cmplx)

        !TODO query the success of this diagonalization

        if(myid.eq.0 ) then
            write (info_str,'(4X,A)') &
                "*********Diagonal Lifetime Eigvals in ps**********"
            call localorb_info(info_str,use_unit,'(A)')
            do p=1, n_print
              pstart = 1+(p-1)*6 
              pend = pstart+5
              if (pend>3*friction_n_active_atoms) pend=3*friction_n_active_atoms 
                write (info_str,'(10X,6E18.9)') &
                    (ps/friction_eigenvalues(j_coord,i_q_point), j_coord=pstart, pend)
                call localorb_info(info_str,use_unit,'(A)')
            enddo
            write (info_str,'(4X,A)') &
                "*************Diagonalized Lifetime Eigvecs********"
            call localorb_info(info_str,use_unit,'(A)')
            do i_coord=1, 3*friction_n_active_atoms
                write (info_str,'(4X,I6)') i_coord
                do p=1, n_print
                  pstart = 1+(p-1)*6 
                  pend = pstart+5
                  if (pend>3*friction_n_active_atoms) pend=3*friction_n_active_atoms 
                  if (p==1) then
                      write (info_str,'(4X,I6,6E18.9)') i_coord, &
                          (dble(friction_eigvecs(i_coord,j_coord,i_q_point)), j_coord=pstart, pend)
                      call localorb_info(info_str,use_unit,'(A)')
                  else
                      write (info_str,'(10X,6E18.9)') &
                          (dble(friction_eigvecs(i_coord,j_coord,i_q_point)), j_coord=pstart, pend)
                      call localorb_info(info_str,use_unit,'(A)')
                  endif
                enddo
            enddo
        endif
       
        !TODO WHAT HAPPENS AFTER DIAGONALIZATION, FRICTION FORCES?!
        !TODO WHAT ABOUT HESSIAN and comparison to it?

        call friction_output_tensor(friction_tensor(:,:,i_q_point))

        if (output_friction_eigenvectors .and. myid.eq.0) then
            call friction_output_eigenvectors_jmol(friction_eigenvalues(:,i_q_point), &
                friction_eigvecs(:,:,i_q_point))
        endif
        
        enddo

        if (allocated(friction_tensor)) deallocate(friction_tensor)
        if (allocated(friction_eigvecs)) deallocate(friction_eigvecs)
        if (allocated(friction_eigenvalues)) deallocate(friction_eigenvalues)

    end subroutine friction_calculation
            
    subroutine friction_calculate_tensor()
    !  NAME
    !    friction_calculate_tensor 
    !  SYNOPSIS
    !  calculates the friction tensor for an isolated system
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)
        
        implicit none

        !  ARGUMENTS

        ! imported variables

        ! local variables
        integer :: i_atom, j_atom, i_coord, j_coord, i_spin
        integer :: atom_counter, n_matrix_elements
        integer :: i_cart, j_cart, i, j, f, a, b, c, n
        integer :: i_cell, n_axis
        integer :: i_basis, j_basis
        integer :: orb_min, orb_homo, orb_lumo, orb_max
        real*8 :: e, tmp, delta, e_sum
        real*8 :: spin_factor, norm, energy, max_energy
        real*8 :: friction

        integer :: success
        
        real*8:: nacs1, nacs2
        real*8, dimension(:,:), allocatable :: G1, G2, S1, S2
        real*8, dimension(:), allocatable :: x_axis
        real*8, dimension(:), allocatable :: spectrum
        real*8, dimension(:), allocatable :: spectrum_tmp

        logical :: kpoint_has_states

        character(len=1200) :: info_str

        integer :: i_window, i_index
       
        integer, dimension(:,:), allocatable :: basis_index

        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis,1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo
       
        if (friction_read_matrices) then
            call friction_read_matrix(first_order_S, first_order_G)
            !call friction_read_matrix(first_order_G)
            if(myid.eq.0 .and. module_is_debugged('friction')) then
                write (info_str,'(4X,A)') &
                    "Finished reading of coupling matrix G."
                call localorb_info(info_str,use_unit,'(A)')
            endif
        else 
            if (numerical_friction) then
                call friction_calculate_H1_and_S1_num()
            else
                call friction_calculate_H1_and_S1_DFPT()
            end if
            if(myid.eq.0 .and. module_is_debugged('friction')) then
                write (info_str,'(4X,A)') &
                    "Finished calculation of first order H and S"
                call localorb_info(info_str,use_unit,'(A)')
            endif
        endif
        
        if (output_first_order_matrices) then
          call friction_output_matrix(first_order_S, first_order_G)
          !call friction_output_matrix(first_order_G)
        endif

        !!!!!!!!!!!!! FRICTION TENSOR !!!!!!!!!!!!!!!!!
        if(myid.eq.0 .and. module_is_debugged('friction')) then
            write (info_str,'(4X,A)') &
                "FRICTION Starting construction of Friction Tensor"
            call localorb_info(info_str,use_unit,'(A)')
        endif

        allocate(G1(n_basis,n_basis))
        allocate(G2(n_basis,n_basis))
        allocate(S1(n_basis,n_basis))
        allocate(S2(n_basis,n_basis))

        if (friction_max_energy<=4.0*friction_broadening_width+friction_perturbation) then
            max_energy = friction_perturbation + friction_broadening_width*4.0
            if (myid.eq.0) then
                write(info_str, '(4X,A)') 'Friction max_energy was below one, '
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str, '(4X,A)') 'set to 4 times friction_broadening_width, this is the recommended minimum value'
                call localorb_info(info_str,use_unit,'(A)')
            endif
        else
            max_energy = friction_max_energy
        endif

        if (friction_output_spectrum .and. myid.eq.0) then
            open (40, file="nacs-spectrum.out")
            write (40,'(A,I8)') &
                "No of components ",3*friction_n_active_atoms
            write (40,'(A,F8.4)') &
                "Discretization length in eV",friction_discretization_length*hartree
            write (40,'(A,I8)') &
                "Number of Bins ",n_axis
        endif

        friction_eigenvalues(:,:) = 0.d0
        friction_tensor(:,:,:) = (0.d0,0.d0)
        friction_eigvecs(:,:,:) = (0.d0,0.d0)

        if(myid.eq.0 .and. module_is_debugged('friction')) then
            write (info_str,'(4X,A,F12.5)') &
                "!!!!!!!!!!!!USING WINDOW!!!!!!! in eV ",friction_broadening_width*27.2114
            call localorb_info(info_str,use_unit,'(A)')
        endif

        atom_counter = -1
        n_matrix_elements = friction_n_active_atoms*3*((friction_n_active_atoms*3)+1)/2
        i_atom = 1
        i_cart = 0
        ATOM_LOOP:  do i_coord = 1, friction_n_active_atoms*3
            i_cart = i_cart + 1
            if (i_cart>3) then
                i_cart = 1
                i_atom = i_atom + 1
            endif
            j_atom = 1
            j_cart = 0
            ATOM_LOOP2: do j_coord = 1, friction_n_active_atoms*3
                j_cart = j_cart + 1
                if (j_cart>3) then
                    j_cart = 1
                    j_atom = j_atom + 1
                endif

                !INCLUDE only calculate upper triangle with diagonal
                if (j_coord<i_coord) cycle
                
                atom_counter = atom_counter + 1  
                !TRIVIAL PARALLELIZATION over matrix elements
                if (myid/=MOD(atom_counter,n_tasks) .or. (myid>n_matrix_elements) ) then
                    cycle
                else

                if (friction_output_spectrum) then
                    !build discrete grid
                    n_axis = int(max_energy/friction_discretization_length) 
                    allocate(x_axis(n_axis))
                    allocate(spectrum(n_axis))
                    allocate(spectrum_tmp(n_axis))

                    energy = 0.0
                    do n=1, n_axis
                        x_axis(n) = energy
                        energy = energy + friction_discretization_length
                    enddo
                endif
                if (friction_output_spectrum) then
                    spectrum(:) = 0.d0
                endif
                friction = 0.d0

                if(myid.eq.0 ) then
                    write (info_str,'(4X,A,I4,I4)') &
                        "Calculating Friction component for ",i_coord,j_coord
                    call localorb_info(info_str,use_unit,'(A)')
                endif
                if (friction_output_spectrum .and. myid.eq.0) then
                    write (40,'(A,I4,I4)') &
                        "Friction component for ",i_coord,j_coord
                endif
                !---------------------
                SPIN_LOOP: do i_spin = 1, n_spin 
                    orb_min = 1
                    orb_homo = 1
                    orb_lumo = 1
                    orb_max = 1

                    do i=1, n_states
                      if (KS_eigenvalue(i,i_spin,1)<=chemical_potential-2.0*max_energy) orb_min = i
                      if (fermi_pop(KS_eigenvalue(i,i_spin,1),chemical_potential,&
                          friction_temperature, n_spin)>=0.001) orb_homo = i
                      if (fermi_pop(KS_eigenvalue(i,i_spin,1),chemical_potential,&
                          friction_temperature, n_spin)>=1.000) orb_lumo = i
                      if (KS_eigenvalue(i,i_spin,1)<=chemical_potential+2.0*max_energy) orb_max = i
                    enddo
                    if(myid.eq.0 .and. module_is_debugged('friction')) then
                        write (info_str,'(4X,A,I4)') &
                            "spin  ",i_spin 
                        call localorb_info(info_str,use_unit,'(A)')
                        write (info_str,'(4X,A,I4,I4,I4,I4)') &
                            "Boundaries  ",orb_min,orb_homo,orb_lumo,orb_max
                        call localorb_info(info_str,use_unit,'(A)')
                    endif

                    !now we know, so construct G
                    G1(:,:) = 0.d0
                    G2(:,:) = 0.d0
                    S1(:,:) = 0.d0
                    S2(:,:) = 0.d0
                    c = 0
                    do b=1,n_basis
                      do a=1, b
                        c = c +1
                        G1(a,b) = first_order_G(i_cart,i_atom,1,c,1,i_spin)
                        G2(a,b) = first_order_G(j_cart,j_atom,1,c,1,i_spin)
                        G1(b,a)  = G1(a,b) 
                        G2(b,a)  = G2(a,b)
                        S1(a,b) = first_order_S(i_cart,i_atom,1,c,1)
                        S2(a,b) = first_order_S(j_cart,j_atom,1,c,1)
                        S1(b,a)  = S1(a,b) 
                        S2(b,a)  = S2(a,b) 
                      enddo
                    enddo

                    do i=orb_min, orb_homo
                        do f=orb_lumo, orb_max
                            e = KS_eigenvalue(f,i_spin,1) &
                                - KS_eigenvalue(i,i_spin,1) 
                            if (e<=tiny(e)) cycle
                            if (e>1.00*max_energy) cycle
                            tmp = fermi_pop(KS_eigenvalue(i,i_spin,1),chemical_potential, &
                                friction_temperature, n_spin)-&
                                fermi_pop(KS_eigenvalue(f,i_spin,1),chemical_potential, &
                                friction_temperature, n_spin)
                            tmp = tmp*(2.0/n_spin)
                            if (abs(tmp)<=0.0001) cycle
                            nacs1 = 0.d0
                            nacs2 = 0.d0
                            e_sum = (KS_eigenvalue(f,i_spin,1) + KS_eigenvalue(i,i_spin,1))/2.d0
                            if (friction_coupling_matrix_mode==0) then
                                do b=1,n_basis
                                  do a=1, n_basis
                                        nacs1 = nacs1 + KS_eigenvector(a,i,i_spin,1)&
                                            *(G1(a,b)-chemical_potential*S1(a,b))&
                                            *KS_eigenvector(b,f,i_spin,1)
                                        nacs2 = nacs2 + KS_eigenvector(a,i,i_spin,1)&
                                            *(G2(a,b)-chemical_potential*S2(a,b))&
                                            *KS_eigenvector(b,f,i_spin,1)
                                  enddo
                                enddo
                            else if (friction_coupling_matrix_mode==1) then
                                do b=1,n_basis
                                  do a=1, n_basis
                                        nacs1 = nacs1 + KS_eigenvector(a,i,i_spin,1)&
                                            *(G1(a,b) -e_sum*S1(a,b) )&
                                            *KS_eigenvector(b,f,i_spin,1)
                                        nacs2 = nacs2 + KS_eigenvector(a,i,i_spin,1)&
                                            *(G2(a,b) -e_sum*S2(a,b) )&
                                        *KS_eigenvector(b,f,i_spin,1)
                                  enddo
                                enddo
                            else if (friction_coupling_matrix_mode==2) then
                                do b=1,n_basis
                                  do a=1, n_basis
                                        nacs1 = nacs1 + KS_eigenvector(a,i,i_spin,1)&
                                            *(G1(a,b) -KS_eigenvalue(i,i_spin,1)*S1(b,a)&
                                            -KS_eigenvalue(f,i_spin,1)*S1(a,b) )&
                                            *KS_eigenvector(b,f,i_spin,1)
                                        nacs2 = nacs2 + KS_eigenvector(a,i,i_spin,1)&
                                            *(G2(a,b) -KS_eigenvalue(i,i_spin,1)*S2(b,a)&
                                            -KS_eigenvalue(f,i_spin,1)*S2(a,b) )&
                                        *KS_eigenvector(b,f,i_spin,1)
                                  enddo
                                enddo
                            endif
                            tmp = tmp*nacs1*nacs2
                            tmp = tmp  / (e)
                            if(myid.eq.0 .and. module_is_debugged('friction')) then
                                write (info_str,'(4X,A,F12.4,1X,E19.8,1X,E19.8)') &
                                    "Excitation  ",e*hartree, &
                                    occ_numbers(i,i_spin,1)-occ_numbers(f,i_spin,1), &
                                    tmp*pi/(ps*hartree*sqrt(masses(i_atom)*masses(j_atom)))
                                call localorb_info(info_str,use_unit,'(A)')
                            endif 
                            if (friction_output_spectrum) then
                                norm = 0.d0
                                spectrum_tmp(:) = 0.d0
                                !integrate each state
                                do n=1, n_axis 
                                    delta = delta_function(x_axis(n),e,friction_broadening_width)
                                    norm = norm + delta
                                    spectrum_tmp(n) = spectrum_tmp(n) + delta*tmp
                                enddo
                                norm = norm * friction_discretization_length
                                spectrum_tmp(:) = spectrum_tmp(:) * norm
                                spectrum(:) = spectrum(:) + spectrum_tmp(:)
                            else
                                delta = delta_function(e,friction_perturbation,friction_broadening_width)
                                friction = friction + delta*tmp/gaussian_norm(e,friction_broadening_width) 
                            endif
                        enddo
                    enddo
                enddo SPIN_LOOP
                
                if (friction_output_spectrum) then 
                    spectrum(:) = spectrum(:) * pi / sqrt(masses(i_atom)*masses(j_atom))
                else
                    friction = friction * pi / sqrt(masses(i_atom)*masses(j_atom))
                endif

                !no we have a full spectrum, we can perform delta function integration
                if (myid.eq.0 .and. friction_output_spectrum.and.module_is_debugged('friction')) then
                    write (info_str,'(4X,A)') &
                        "Excitation energy in eV   Coupling element in 1/eV*ps" 
                    call localorb_info(info_str,use_unit,'(A)')
                    write (info_str,'(4X,A)') &
                        "==========================================" 
                    call localorb_info(info_str,use_unit,'(A)')
                endif
                if (myid.eq.0 .and. friction_output_spectrum) then
                    write (40,'(A)') &
                        "Excitation energy in eV   Coupling element in 1/eV*ps" 
                    write (40,'(A)') &
                        "==========================================" 
                endif
                norm=0.d0
                if (friction_output_spectrum) then
                    do n=1, n_axis 
                        delta = delta_function(x_axis(n),friction_perturbation,friction_window_size)
                        norm = norm + delta
                        friction = friction + spectrum(n)*delta
                        if (myid.eq.0 .and. module_is_debugged('friction')) then
                            write (info_str,'(E16.6,4X,E16.6)') &
                                x_axis(n)*hartree, spectrum(n)/(ps)
                            call localorb_info(info_str,use_unit,'(A)')
                        endif
                        if (myid.eq.0 .and. friction_output_spectrum) then
                            write (40,'(E16.6,5X,E16.6)') &
                                x_axis(n)*hartree, spectrum(n)/(ps)
                        endif
                    enddo
                    friction = friction/norm
                endif 
                !EVALUATING FRICTION AT 0, Fermi level
                if (allocated(x_axis)) deallocate(x_axis)
                if (allocated(spectrum)) deallocate(spectrum)
                if (allocated(spectrum_tmp)) deallocate(spectrum_tmp)

                friction_tensor(i_coord,j_coord,1) = cmplx(friction,0.d0)
                friction_tensor(j_coord,i_coord,1) = friction_tensor(i_coord,j_coord,1)
                !END OF INTERNAL LOOPS
                !---------------------
                
                !END OF PARALLELIZATION LOOP OVER ATOMS
                endif
            end do ATOM_LOOP2 
        end do ATOM_LOOP 
        
        !SUM OVER ALL CPUS to collect friction matrix
        call sync_matrix_complex(friction_tensor(:,:,1),3*friction_n_active_atoms,3*friction_n_active_atoms) 

        if (friction_output_spectrum) then
            if(myid.eq.0) then
                close(40)
            endif
        endif
        
        deallocate(basis_index)
        deallocate(G1)
        deallocate(G2)
        deallocate(S1)
        deallocate(S2)

        if(myid.eq.0 .and. module_is_debugged('friction')) then
            write (info_str,'(4X,A)') &
                "Finished calculation of Friction Tensor"
            call localorb_info(info_str,use_unit,'(A)')
        endif 

    end subroutine friction_calculate_tensor
    
    subroutine friction_calculate_tensor_p1()
    !  NAME
    !    friction_calculate_tensor_p1 
    !  SYNOPSIS
    !    calculates the friction tensor for a periodic system
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)
        
        implicit none

        !  ARGUMENTS
        ! imported variables
        ! local variables
        integer :: i_atom, j_atom, i_coord, j_coord, i_spin
        integer :: i_k_point, i_k_point2, i_k_point3
        integer :: i_k_task, i_k_task2
        integer :: i_q_point
        integer :: i_cart, j_cart, i, j, f, a, b,c, n
        integer :: i_basis, j_basis, i_index
        integer :: i_basis_1, i_basis_2
        integer :: i_cell, n_axis, n_axis_neg
        integer :: i_cell1, i_cell2, i_cell3, i_cell_n
        integer :: cell_x, cell_y, cell_z
        integer :: ind_cell1, ind_cell2, ind_cell3
        integer :: orb_min, orb_homo, orb_lumo, orb_max
        integer :: id_send, id_recv, n_k_points2
      
        integer :: success

        real*8 :: states_in_window
        real*8 :: e, tmp, tmp2, delta, e_sum
        real*8 :: spin_factor, norm, energy, max_energy
        complex*16 :: friction, friction_tmp
        complex*16 :: kphase2
        complex*16 :: nacs1, nacs2, tmp_cmplx
        
        real*8, dimension(:,:,:), allocatable :: KS_eigenvector_tmp
        complex*16, dimension(:,:,:), allocatable :: KS_eigenvector_tmp_complex
        real*8, dimension(:,:,:), allocatable :: KS_eigenvector_tmp2
        complex*16, dimension(:,:,:), allocatable :: KS_eigenvector_tmp_complex2

        complex*16, dimension(:,:), allocatable :: G1, G2, S1, S2
        real*8, dimension(:), allocatable :: x_axis
        complex*16, dimension(:), allocatable :: spectrum
        complex*16, dimension(:), allocatable :: spectrum_tmp

        complex*16, dimension(:,:,:,:), allocatable :: k_k_tensor
        real*8 :: k_diff(3),q_vec(3)
        complex*16 :: friction_k_k

        real*8 :: time
        logical :: kpoint_has_states, info, found_qdiff
        character(len=1200) :: info_str
        integer :: i_window
        integer, dimension(:,:), allocatable :: basis_index
 
        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis,1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!READ / WRITE !!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (friction_read_matrices) then
            call friction_read_matrix_p1(first_order_S,first_order_G,first_order_S_cmplx,first_order_G_cmplx)
            !call friction_read_matrix_p1(first_order_G,first_order_G_cmplx)
            if(myid.eq.0 .and. module_is_debugged('friction')) then
                write (info_str,'(4X,A)') &
                    "Finished reading of first order H and S"
                call localorb_info(info_str,use_unit,'(A)')
            endif
        else 
            if (numerical_friction) then
                call friction_calculate_H1_and_S1_num()
            else
                call friction_calculate_H1_and_S1_DFPT()
                !call aims_stop('DFPT Friction not implemented yet!')
            end if
            if(myid.eq.0 .and. module_is_debugged('friction')) then
                write (info_str,'(4X,A)') &
                    "Finished calculation of first order H and S"
                call localorb_info(info_str,use_unit,'(A)')
            endif
        endif
        
        if (output_first_order_matrices) then
          !call friction_output_matrix_p1(first_order_G, first_order_G_cmplx)
          call friction_output_matrix_p1(first_order_S,first_order_G,first_order_S_cmplx,first_order_G_cmplx)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!READ / WRITE !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!! FRICTION TENSOR !!!!!!!!!!!!!!!!!
        if(myid.eq.0 ) then
            write (info_str,'(4X,A)') &
                "FRICTION Starting construction of Friction Tensor"
            call localorb_info(info_str,use_unit,'(A)')
        endif

        allocate(G1(n_basis,n_basis))
        allocate(G2(n_basis,n_basis))
        allocate(S1(n_basis,n_basis))
        allocate(S2(n_basis,n_basis))

        if (friction_max_energy<4.0*friction_broadening_width+friction_perturbation) then
            max_energy = friction_perturbation + friction_broadening_width*4.0
            if (myid.eq.0) then
                write(info_str, '(4X,A)') '**Friction max_energy was below one, **'
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str, '(4X,A)') '**set to 4 times friction_broadening_width, this is the recommended minimum value**'
                call localorb_info(info_str,use_unit,'(A)')
            endif
        else
            max_energy = friction_max_energy
        endif
        if (friction_output_spectrum) then
          !build discrete grid
          n_axis = int(max_energy/friction_discretization_length)
          allocate(x_axis(n_axis))
          allocate(spectrum(n_axis))
          allocate(spectrum_tmp(n_axis))
        
          energy = 0.0
          do n=1, n_axis
              x_axis(n) = energy
              energy = energy + friction_discretization_length
          enddo
          energy = 0.0
        endif

        if (friction_output_spectrum .and. myid.eq.0) then
            open (40, file="nacs-spectrum.out")
            write (40,'(A,I8)') &
                "No of components ",3*friction_n_active_atoms
            write (40,'(A,F8.4)') &
                "Discretization length in eV ",friction_discretization_length*hartree
            write (40,'(A,I8)') &
                "Number of Bins ",n_axis
        endif

        friction_tensor(:,:,:) = (0.d0,0.d0)
        friction_eigvecs(:,:,:) = (0.d0,0.d0)
        friction_eigenvalues(:,:) = (0.d0,0.d0)

        if(myid.eq.0 .and. module_is_debugged('friction')) then
          write (info_str,'(4X,A,F12.5)') &
              "!!!!!!!!!!!!USING Broadening!!!!!!! in eV ",friction_broadening_width*27.2114
          call localorb_info(info_str,use_unit,'(A)')
        endif

        !have not figured out knotk and scalapack yet
        if (myid.eq.0 .and. use_scalapack .and. friction_knotk) then
          write (info_str,'(4X,A)') &
              "SCALAPACK AND friction_q_integration are currently not compatible ... set friction_q_integration to .false."
          call localorb_info(info_str,use_unit,'(A)')
          friction_knotk = .false.
        endif
          
        !ALL STUFF WE NEED FOR KNOTK  
        if (real_eigenvectors) then
            allocate( KS_eigenvector_tmp(n_basis,n_states,n_spin),stat=success)
        else
            allocate( KS_eigenvector_tmp_complex(n_basis,n_states,n_spin),stat=success)
        endif
        if (friction_knotk) then
            n_k_points2 = n_k_points
            allocate(k_k_tensor(3*friction_n_active_atoms,3*friction_n_active_atoms,&
                n_k_points,n_k_points))
        else
            n_k_points2 = 1
        endif
       
        do i_q_point=1, friction_n_q_points
            if (friction_knotk) then 
                k_k_tensor=0.d0
                q_vec = k_point_list(i_q_point,:)    
            else
                q_vec = 0.d0 
            endif
            if(myid.eq.0 .and. friction_knotk .and. module_is_debugged('friction')) then
                write (info_str,'(4X,A,3F12.6)') &
                    "qpoint = ",k_point_list(i_q_point,1),k_point_list(i_q_point,2),k_point_list(i_q_point,3)
                call localorb_info(info_str,use_unit,'(A)')
            endif 

        i_k_task = 0
        KPOINT_LOOP: do i_k_point=1, n_k_points,1
        !if (use_scalapack) then
            !if (my_k_point /= i_k_point) cycle
        !else
            !if (myid/=MOD(i_k_point, n_tasks) .or. (myid>n_k_points) ) then 
              !cycle
            !else
        
        KPOINT_LOOP2: do i_k_point3=1, n_k_points2, 1
            found_qdiff = .false.
            i_k_point2 = i_k_point3 
            k_diff = k_point_list(i_k_point2,:)-k_point_list(i_k_point,:)
            if (abs(k_diff(1)-q_vec(1))<=tiny(1.d0) .and.&
                abs(k_diff(2)-q_vec(2))<=tiny(1.d0) .and.&
                abs(k_diff(3)-q_vec(3))<=tiny(1.d0) ) then
                found_qdiff = .true.
                exit
            endif
        enddo KPOINT_LOOP2
        if (.not. friction_knotk)  i_k_point2 = i_k_point
        if (.not. found_qdiff .and.friction_knotk) cycle

        id_send = MOD(i_k_point2,n_tasks)
        id_recv = MOD(i_k_point,n_tasks)
        
        !BEGIN MPI 
        if (myid.eq.MOD(i_k_point, n_tasks) .and. (myid<=n_k_points) ) then 
        !endif
        i_k_task = i_k_task +1

        !IF WE WANT TO DO k not k' every CPU needs to have 
        !all eigenvectors at all k points
        if (id_send .ne. id_recv) then
          if (myid.eq. id_recv) then
            if (real_eigenvectors) then
               call receive_real_vector(KS_eigenvector_tmp,n_basis*n_states*n_spin,id_send)
            else
               call receive_complex_vector(KS_eigenvector_tmp_complex,n_basis*n_states*n_spin,id_send)
            endif
          endif
        else
          if (real_eigenvectors) then
            KS_eigenvector_tmp = KS_eigenvector(:,:,:,i_k_task)
          else
            KS_eigenvector_tmp_complex =KS_eigenvector_complex(:,:,:,i_k_task) 
          endif
        endif

        !if(myid.eq.0 .and. module_is_debugged('friction')) then
        if(myid.eq.0) then
            write (info_str,'(4X,A,I4,1X,I4)') &
                "kpoint",i_k_point,i_k_point2
            call localorb_info(info_str,use_unit,'(A)')
        endif
        
        SPIN_LOOP: do i_spin = 1, n_spin 
            orb_min = 1
            orb_homo = 1
            orb_lumo = 1
            orb_max = 1
            do i=1, n_states
              if (KS_eigenvalue(i,i_spin,i_k_point)<=chemical_potential-2.0*max_energy) orb_min = i
              if (fermi_pop(KS_eigenvalue(i,i_spin,i_k_point),chemical_potential, &
                        friction_temperature, n_spin)>=0.001) orb_homo = i
              if (fermi_pop(KS_eigenvalue(i,i_spin,i_k_point2),chemical_potential, &
                        friction_temperature, n_spin)>=1.000) orb_lumo = i
              if (KS_eigenvalue(i,i_spin,i_k_point2)<=chemical_potential+2.0*max_energy) orb_max = i
            enddo
            if(myid.eq.0 .and. module_is_debugged('friction')) then
                write (info_str,'(4X,A,I4,1X,A,I4,I4,I4,I4)') &
                    "spin  ", i_spin, "Boundaries  ",orb_min,orb_homo,orb_lumo,orb_max
                call localorb_info(info_str,use_unit,'(A)')
            endif
            !start by checking if there are any important states to consider
            kpoint_has_states = .false.
            do i=orb_min, orb_homo
              do f=orb_lumo, orb_max
                  e = KS_eigenvalue(f,i_spin,i_k_point2) &
                      - KS_eigenvalue(i,i_spin,i_k_point) 
                  if (e<=tiny(e)) cycle
                  if (e>(1.0*max_energy)) cycle
                  !if (abs(e)>1.0*max_energy) cycle
                  tmp = fermi_pop(KS_eigenvalue(i,i_spin,i_k_point),chemical_potential, &
                      friction_temperature, n_spin)-&
                      fermi_pop(KS_eigenvalue(f,i_spin,i_k_point2),chemical_potential, &
                      friction_temperature, n_spin)
                  !tmp = fermi_pop(KS_eigenvalue(i,i_spin,i_k_point),chemical_potential, &
                  !    friction_temperature, n_spin)
                  !tmp2 =fermi_pop(KS_eigenvalue(f,i_spin,i_k_point2),chemical_potential, &
                  !    friction_temperature, n_spin)
                  !if initial state is empty or final state is already full, cycle
                  !if (tmp<=0.0001 .or. tmp2>=(1.9999/n_spin)) cycle
                  !tmp = tmp*(1.0-tmp2)*(2.0/n_spin)
                  tmp = tmp*(2.0/n_spin)
                  !tmp = fermi_pop(KS_eigenvalue(i,i_spin,i_k_point),chemical_potential, &
                      !friction_temperature, n_spin)-&
                      !fermi_pop(KS_eigenvalue(f,i_spin,i_k_point2),chemical_potential, &
                      !friction_temperature, n_spin)
                  if (abs(tmp)<=0.0001) cycle
                  kpoint_has_states = .true.
                  exit
              enddo
            enddo
            if (.not. kpoint_has_states) cycle
        !---------------------

        i_atom = 1
        i_cart = 0
        ATOM_LOOP:  do i_coord = 1, friction_n_active_atoms*3
            i_cart = i_cart + 1
            if (i_cart>3) then
                i_cart = 1
                i_atom = i_atom + 1
            endif
            j_atom = 1
            j_cart = 0
            ATOM_LOOP2: do j_coord = 1, friction_n_active_atoms*3
                j_cart = j_cart + 1
                if (j_cart>3) then
                    j_cart = 1
                    j_atom = j_atom + 1
                endif
                
                friction = (0.d0,0.d0)
                if (friction_output_spectrum) then 
                    spectrum(:) = (0.d0,0.d0)
                endif

                !INCLUDE to only calculate upper triangle with diagonal
                if (j_coord<i_coord) cycle
                if(myid.eq.0 .and. module_is_debugged('friction')) then
                    write (info_str,'(4X,A,I4,I4)') &
                        "Calculating Friction component for ",i_coord,j_coord
                    call localorb_info(info_str,use_unit,'(A)')
                endif
                if (friction_output_spectrum .and. myid.eq.0) then
                    write (40,'(A,I4,I4)') &
                        "Friction component for ",i_coord,j_coord
                endif

                    !now we know, so construct G
                    G1(:,:) = (0.d0,0.d0)
                    G2(:,:) = (0.d0,0.d0)
                    S1(:,:) = (0.d0,0.d0)
                    S2(:,:) = (0.d0,0.d0)
                    if (friction_use_complex_matrices) then
                        if (real_eigenvectors) then
                            c = 0
                            do b=1, n_basis
                              do a=1, b
                                c = c + 1
                                G1(a,b) = cmplx(first_order_G(i_cart,i_atom,i_k_task,c,1,i_spin),0.d0)
                                G2(a,b) = cmplx(first_order_G(j_cart,j_atom,i_k_task,c,1,i_spin),0.d0)
                                G1(b,a) = conjg(G1(a,b))
                                G2(b,a) = conjg(G2(a,b))
                                S1(a,b) = cmplx(first_order_S(i_cart,i_atom,i_k_task,c,1),0.d0)
                                S2(a,b) = cmplx(first_order_S(j_cart,j_atom,i_k_task,c,1),0.d0)
                                S1(b,a) = conjg(S1(a,b))
                                S2(b,a) = conjg(S2(a,b))
                              enddo
                            enddo
                        else
                            c = 0
                            do b=1, n_basis
                              do a=1, b
                                c = c + 1
                                G1(a,b) = first_order_G_cmplx(i_cart,i_atom,i_k_task,c,1,i_spin)
                                G2(a,b) = first_order_G_cmplx(j_cart,j_atom,i_k_task,c,1,i_spin)
                                G1(b,a) = conjg(G1(a,b))
                                G2(b,a) = conjg(G2(a,b))
                                S1(a,b) = first_order_S_cmplx(i_cart,i_atom,i_k_task,c,1)
                                S2(a,b) = first_order_S_cmplx(j_cart,j_atom,i_k_task,c,1)
                                S1(b,a) = conjg(S1(a,b))
                                S2(b,a) = conjg(S2(a,b))
                              enddo
                            enddo
                        endif
                    else
                        if (friction_knotk) then
                          !build the DFPT supercell matrices, S_full, M_full
                          c = 0
                          do b=1, n_basis
                            do a=1, b
                              c = c + 1
                              do i_cell1 = 1, n_cells_in_hamiltonian-1
                                do i_cell2 = 1, n_cells_in_sc_DFPT
                                  cell_x = cell_index_sc_DFPT(i_cell2,1)
                                  cell_y = cell_index_sc_DFPT(i_cell2,2)
                                  cell_z = cell_index_sc_DFPT(i_cell2,3)
                                  kphase2 = k_phase_base(1,i_k_point2) ** cell_x &
                                      & * k_phase_base(2,i_k_point2) ** cell_y &
                                      & * k_phase_base(3,i_k_point2) ** cell_z
                                  !add the normal contribution
                                  G1(a,b) = G1(a,b) +&
                                    first_order_G(i_cart,i_atom,i_cell1,c,1,i_spin)*&
                                    k_phase(i_cell1,i_k_point)*conjg(kphase2)
                                  G2(a,b) = G2(a,b) +&
                                    first_order_G(j_cart,j_atom,i_cell1,c,1,i_spin)*&
                                    k_phase(i_cell1,i_k_point)*conjg(kphase2)
                                  S1(a,b) = S1(a,b) +&
                                    first_order_S(i_cart,i_atom,i_cell1,c,1)*&
                                    k_phase(i_cell1,i_k_point)*conjg(kphase2)
                                  S2(a,b) = S2(a,b) +&
                                    first_order_S(j_cart,j_atom,i_cell1,c,1)*&
                                    k_phase(i_cell1,i_k_point)*conjg(kphase2)
                                enddo
                              enddo
                              G1(b,a) = conjg(G1(a,b))
                              G2(b,a) = conjg(G2(a,b))
                              S1(b,a) = conjg(S1(a,b))
                              S2(b,a) = conjg(S2(a,b))
                            enddo
                          enddo

                          !do b=1, n_basis_sc_DFPT
                            !i_basis_2 = Cbasis_to_basis(b)
                            !i_cell2 = center_to_cell(Cbasis_to_center(b))
                            !!i_cell2 = center_to_cell(Cbasis_to_center(b))
                            !do a=1,n_basis_sc_DFPT
                              !!!do i_cell2=1, n_cells
                              !i_basis_1 = Cbasis_to_basis(a)
                              !if (i_basis_1<=i_basis_2) then
                                !i_cell1 = center_to_cell(Cbasis_to_center(a))
                                !i_cell_n = position_in_hamiltonian_PBC(i_cell1, i_cell2)
                                !!if (i_cell_n<n_cells_in_hamiltonian) then
                                !if (i_cell_n<n_cells_in_sc_DFPT) then
                                  !G1(i_basis_1,i_basis_2) = G1(i_basis_1,i_basis_2) + &
                                  !k_phase(i_cell2,i_k_point2)*conjg(k_phase(i_cell1,i_k_point))*&
                                  !(first_order_G(i_cart,i_atom,i_cell_n,basis_index(i_basis_1,i_basis_2),1,i_spin))
                                  !G2(i_basis_1,i_basis_2) = G2(i_basis_1,i_basis_2) + &
                                  !k_phase(i_cell2,i_k_point2)*conjg(k_phase(i_cell1,i_k_point))*&
                                  !(first_order_G(j_cart,j_atom,i_cell_n,basis_index(i_basis_1,i_basis_2),1,i_spin))
                                  !S1(i_basis_1,i_basis_2) = S1(i_basis_1,i_basis_2) + &
                                  !k_phase(i_cell2,i_k_point2)*conjg(k_phase(i_cell1,i_k_point))*&
                                  !(first_order_S(i_cart,i_atom,i_cell_n,basis_index(i_basis_1,i_basis_2),1))
                                  !S2(i_basis_1,i_basis_2) = S2(i_basis_1,i_basis_2) + &
                                  !k_phase(i_cell2,i_k_point2)*conjg(k_phase(i_cell1,i_k_point))*&
                                  !(first_order_S(j_cart,j_atom,i_cell_n,basis_index(i_basis_1,i_basis_2),1))
                                !endif
                              !endif
                            !enddo
                          !enddo
                          !do b = 1, n_basis
                            !do a = 1, b  
                              !G1(b,a) = conjg(G1(a,b))
                              !G2(b,a) = conjg(G2(a,b))
                              !S1(b,a) = conjg(S1(a,b))
                              !S2(b,a) = conjg(S2(a,b))
                            !end do
                          !end do
                        else
                          c = 0
                          do b=1,n_basis
                            do a=1, b
                              c = c + 1
                              do i_cell=1, n_cells_in_hamiltonian-1
                                G1(a,b) = G1(a,b) +&
                                  first_order_G(i_cart,i_atom,i_cell,c,1,i_spin)*&
                                  k_phase(i_cell,i_k_point)
                                G2(a,b) = G2(a,b) +&
                                  first_order_G(j_cart,j_atom,i_cell,c,1,i_spin)*&
                                  k_phase(i_cell,i_k_point)
                                S1(a,b) = S1(a,b) +&
                                  first_order_S(i_cart,i_atom,i_cell,c,1)*&
                                  k_phase(i_cell,i_k_point)
                                S2(a,b) = S2(a,b) +&
                                  first_order_S(j_cart,j_atom,i_cell,c,1)*&
                                  k_phase(i_cell,i_k_point)
                              enddo
                              G1(b,a) = conjg(G1(a,b))
                              G2(b,a) = conjg(G2(a,b))
                              S1(b,a) = conjg(S1(a,b))
                              S2(b,a) = conjg(S2(a,b))
                            enddo
                          enddo
                        endif
                    endif

                    do i=orb_min, orb_homo
                        do f=orb_lumo, orb_max
                            !THIS IS THE POINT WHERE WE SHOULD DISTINGUISH BETWEEN DIFFERENT 
                            !SCALAPACK tasks
                            e = KS_eigenvalue(f,i_spin,i_k_point2) &
                                - KS_eigenvalue(i,i_spin,i_k_point) 
                            if (e<=tiny(e)) cycle
                            if (e>1.0*max_energy) cycle
                            tmp = fermi_pop(KS_eigenvalue(i,i_spin,i_k_point),chemical_potential, &
                                friction_temperature, n_spin)-&
                                fermi_pop(KS_eigenvalue(f,i_spin,i_k_point2),chemical_potential, &
                                friction_temperature, n_spin)
                            tmp = tmp*(2.0/n_spin)
                            if (abs(tmp)<=0.0001) cycle
                            nacs1 = (0.d0,0.d0)
                            nacs2 = (0.d0,0.d0)
                            e_sum = (KS_eigenvalue(f,i_spin,i_k_point2) + KS_eigenvalue(i,i_spin,i_k_point))/2.d0 
                            !ONLY AT THIS POINT WE ACTUALLY NEED EIGENVECTORS, namely for state i and state f
                            if (friction_coupling_matrix_mode==0) then
                                if (real_eigenvectors) then
                                  do a=1, n_basis
                                    do b=1,n_basis
                                      nacs1 = nacs1 + cmplx(KS_eigenvector(a,i,i_spin,i_k_task),0.d0)&
                                        *(G1(a,b)-chemical_potential*S1(a,b))&
                                        *cmplx(KS_eigenvector_tmp(b,f,i_spin),0.d0)
                                      nacs2 = nacs2 + cmplx(KS_eigenvector(a,i,i_spin,i_k_task),0.d0)&
                                        *(G2(a,b)-chemical_potential*S2(a,b))&
                                        *cmplx(KS_eigenvector_tmp(b,f,i_spin),0.d0)
                                    enddo
                                  enddo
                                else
                                  do a=1, n_basis
                                    do b=1,n_basis
                                      nacs1 = nacs1 + conjg(KS_eigenvector_complex(a,i,i_spin,i_k_task))&
                                        *(G1(a,b)-chemical_potential*S1(a,b))&
                                        *KS_eigenvector_tmp_complex(b,f,i_spin)
                                      nacs2 = nacs2 + conjg(KS_eigenvector_complex(a,i,i_spin,i_k_task))&
                                        *(G2(a,b)-chemical_potential*S2(a,b))&
                                        *KS_eigenvector_tmp_complex(b,f,i_spin)
                                    enddo
                                  enddo
                                endif
                            else if (friction_coupling_matrix_mode==1) then
                                if (real_eigenvectors) then
                                  do a=1, n_basis
                                    do b=1,n_basis
                                      nacs1 = nacs1 + cmplx(KS_eigenvector(a,i,i_spin,i_k_task),0.d0)&
                                        *(G1(a,b) - e_sum*S1(a,b))&
                                        *cmplx(KS_eigenvector_tmp(b,f,i_spin),0.d0)
                                      nacs2 = nacs2 + cmplx(KS_eigenvector(a,i,i_spin,i_k_task),0.d0)&
                                        *(G2(a,b) -e_sum*S2(a,b))&
                                        *cmplx(KS_eigenvector_tmp(b,f,i_spin),0.d0)
                                    enddo
                                  enddo
                                else
                                  do a=1, n_basis
                                    do b=1,n_basis
                                      nacs1 = nacs1 + conjg(KS_eigenvector_complex(a,i,i_spin,i_k_task))&
                                        *(G1(a,b) -e_sum*S1(a,b))&
                                        *KS_eigenvector_tmp_complex(b,f,i_spin)
                                      nacs2 = nacs2 + conjg(KS_eigenvector_complex(a,i,i_spin,i_k_task))&
                                        *(G2(a,b) -e_sum*S2(a,b))&
                                        *KS_eigenvector_tmp_complex(b,f,i_spin)
                                    enddo
                                  enddo
                                endif
                            else if (friction_coupling_matrix_mode==2) then
                                if (real_eigenvectors) then
                                  do a=1, n_basis
                                    do b=1,n_basis
                                      nacs1 = nacs1 + cmplx(KS_eigenvector(a,i,i_spin,i_k_task),0.d0)&
                                        *(G1(a,b) - KS_eigenvalue(i,i_spin,i_k_point)*conjg(S1(b,a)) &
                                        -KS_eigenvalue(f,i_spin,i_k_point2)*S1(a,b))&
                                        *cmplx(KS_eigenvector_tmp(b,f,i_spin),0.d0)
                                      nacs2 = nacs2 + cmplx(KS_eigenvector(a,i,i_spin,i_k_task),0.d0)&
                                        *(G2(a,b) -KS_eigenvalue(i,i_spin,i_k_point)*conjg(S2(b,a)) &
                                        -KS_eigenvalue(f,i_spin,i_k_point2)*S2(a,b))&
                                        *cmplx(KS_eigenvector_tmp(b,f,i_spin),0.d0)
                                    enddo
                                  enddo
                                else
                                  do a=1, n_basis
                                    do b=1,n_basis
                                      nacs1 = nacs1 + conjg(KS_eigenvector_complex(a,i,i_spin,i_k_task)) &
                                        *(G1(a,b) -KS_eigenvalue(i,i_spin,i_k_point)*conjg(S1(b,a)) &
                                        -KS_eigenvalue(f,i_spin,i_k_point2)*S1(a,b))&
                                        *KS_eigenvector_tmp_complex(b,f,i_spin)
                                      nacs2 = nacs2 + conjg(KS_eigenvector_complex(a,i,i_spin,i_k_task))&
                                        *(G2(a,b) -KS_eigenvalue(i,i_spin,i_k_point)*conjg(S2(b,a))- &
                                        KS_eigenvalue(f,i_spin,i_k_point2)*S2(a,b)) &
                                        *KS_eigenvector_tmp_complex(b,f,i_spin)
                                    enddo
                                  enddo
                                endif
                            else
                                !pass
                            endif
                            tmp_cmplx = cmplx(tmp,0.d0)*conjg(nacs1)*nacs2
                            tmp_cmplx = tmp_cmplx / &
                                    cmplx(e,0.d0)
                            if(myid.eq.0 .and. module_is_debugged('friction')) then
                                write (info_str,'(4X,A,I4,I4,F12.4,1X,2E19.8,1X,E19.8)') &
                                    "Excitation  ",i,f,e*hartree, tmp,&
                                        (tmp_cmplx*pi*k_weights(i_k_point))/&
                                            cmplx(ps*hartree*sqrt(masses(i_atom)*masses(j_atom)),0.d0)
                                call localorb_info(info_str,use_unit,'(A)')
                            endif 
                            
                            if (friction_output_spectrum) then
                                spectrum_tmp(:) = (0.d0,0.d0)
                                !integrate each state
                                norm = 0.d0
                                do n=1, n_axis 
                                    delta = delta_function(x_axis(n),e,friction_broadening_width)
                                    norm = norm + delta
                                    spectrum_tmp(n) = spectrum_tmp(n) + delta*tmp_cmplx 
                                enddo
                                norm = norm * friction_discretization_length
                                if (friction_knotk) then
                                    spectrum_tmp(:) = spectrum_tmp(:) * &
                                        cmplx(((k_weights(i_k_point)*k_weights(i_k_point2))/norm),0.d0)
                                else
                                    spectrum_tmp(:) = spectrum_tmp(:) * &
                                        cmplx((k_weights(i_k_point)/norm),0.d0)
                                endif
                                spectrum(:) = spectrum(:) + spectrum_tmp(:)
                            endif
                            delta = delta_function(e,friction_perturbation,friction_broadening_width)
                            friction_tmp = delta*tmp_cmplx/gaussian_norm(e,friction_broadening_width)
                            if (friction_knotk) then
                                friction_tmp = friction_tmp * &
                                    cmplx((k_weights(i_k_point)*k_weights(i_k_point2)),0.d0)
                            else
                                friction_tmp = friction_tmp * &
                                    cmplx(k_weights(i_k_point),0.d0)
                            endif
                            friction = friction + friction_tmp
                        enddo
                    enddo
                if (friction_output_spectrum) then
                    spectrum(:) = spectrum(:) * cmplx(pi/sqrt(masses(i_atom)*masses(j_atom)),0.d0)
                else
                    friction = friction * cmplx(pi/sqrt(masses(i_atom)*masses(j_atom)),0.d0)
                endif
            
                !!now we have a full spectrum, we can perform delta function integration
                !if (myid.eq.0 .and. friction_output_spectrum .and. module_is_debugged('friction')) then
                    !write (info_str,'(4X,A)') &
                        !"Excitation energy in eV   Re(Coupling element)   Im(Coupling element) in 1/ps" 
                    !call localorb_info(info_str,use_unit,'(A)')
                    !write (info_str,'(4X,A)') &
                        !"==========================================" 
                    !call localorb_info(info_str,use_unit,'(A)')
                !endif
                !if (myid.eq.0 .and. friction_output_spectrum) then
                    !write (40,'(A)') &
                        !"Excitation energy in eV   Re(Coupling element)   Im(Coupling element) in 1/ps" 
                    !write (40,'(A)') &
                        !"==========================================" 
                !endif
                norm = 0.d0
                !TODO OUTPUT
                if (friction_output_spectrum) then
                    do n=1, n_axis 
                        delta = delta_function(x_axis(n),friction_perturbation,friction_window_size)
                        norm = norm + delta
                        friction = friction + spectrum(n)*cmplx(delta,0.d0)
                        !if (myid.eq.0 .and. module_is_debugged('friction')) then
                            !write (info_str,'(4X,E16.6,1X,2E16.6)') &
                                !x_axis(n)*hartree, spectrum(n)/cmplx(ps,0.d0)
                            !call localorb_info(info_str,use_unit,'(A)')
                        !endif
                        !if (myid.eq.0) then
                            !write (40,'(E16.6,5X,2E16.6)') &
                                !x_axis(n)*hartree, spectrum(n)/cmplx(ps,0.d0)
                        !endif
                    enddo
                    friction = (friction)/cmplx(norm,0.d0)
                endif

                if (friction_knotk) then
                    k_k_tensor(i_coord,j_coord,i_k_point,i_k_point2) = &
                        k_k_tensor(i_coord,j_coord,i_k_point,i_k_point2) + friction
                    if (i_coord /= j_coord) then
                      k_k_tensor(j_coord,i_coord,i_k_point,i_k_point2) = & 
                        k_k_tensor(j_coord,i_coord,i_k_point,i_k_point2) + conjg(friction)
                    endif
                else
                    friction_tensor(i_coord,j_coord,1) = friction_tensor(i_coord,j_coord,1) +&
                        friction
                    if (i_coord /= j_coord) then
                      friction_tensor(j_coord,i_coord,1) = friction_tensor(j_coord,i_coord,1)+&
                          conjg(friction)
                    endif
                endif
                !END OF INTERNAL LOOPS
                !---------------------
            end do ATOM_LOOP2 
        end do ATOM_LOOP 
        enddo SPIN_LOOP
        !collect all different k k' components
        else
          if (myid.eq. id_send) then
            i_k_task2 = (i_k_point2/n_tasks)+1
            if (real_eigenvectors) then
                call send_real_vector(KS_eigenvector(:,:,:,i_k_task2),n_basis*n_states*n_spin,id_recv)
            else
                call send_complex_vector(KS_eigenvector_complex(:,:,:,i_k_task2),n_basis*n_states*n_spin,id_recv)
            endif
          else 
            !pass 
          endif
        endif
        !END MPI
        enddo KPOINT_LOOP

        if (friction_knotk) then
            call sync_vector_complex(k_k_tensor,9*friction_n_active_atoms*friction_n_active_atoms*&
                n_k_points*n_k_points,MPI_COMM_WORLD)
            !call sync_matrix_complex(k_k_tensor,n_k_points,n_k_points) 
            do i_k_point=1, n_k_points
              do i_k_point2=1, n_k_points
                friction_tensor(:,:,i_q_point) = friction_tensor(:,:,i_q_point) + &
                    k_k_tensor(:,:,i_k_point,i_k_point2)
              enddo
            enddo
        else
            call sync_vector_complex(friction_tensor,9*friction_n_active_atoms*&
                friction_n_active_atoms,MPI_COMM_WORLD)
        endif

        enddo !i_q_point

        if (friction_output_spectrum) then
            if(myid.eq.0) then
                close(40)
            endif
        endif

        if (allocated(KS_eigenvector_tmp)) deallocate(KS_eigenvector_tmp)
        if (allocated(KS_eigenvector_tmp_complex)) deallocate(KS_eigenvector_tmp_complex)
        if (allocated(x_axis)) deallocate(x_axis)
        if (allocated(spectrum)) deallocate(spectrum)
        if (allocated(spectrum_tmp)) deallocate(spectrum_tmp)
        if (allocated(k_k_tensor)) deallocate(k_k_tensor)

        deallocate(G1)
        deallocate(G2)
        deallocate(basis_index)

        if(myid.eq.0 .and. module_is_debugged('friction')) then
            write (info_str,'(4X,A)') &
                "Finished calculation of Friction Tensor"
            call localorb_info(info_str,use_unit,'(A)')
        endif 

    end subroutine friction_calculate_tensor_p1
    
    subroutine friction_calculate_H1_and_S1_num( opt_S)
    !  NAME
    !    friction_calculate_H1_and_S1_num
    !  SYNOPSIS
    !    calculates the H1 and S1 using finite differences.
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)

        implicit none

        ! imported variables
        integer :: info 

        logical, optional, intent(in) :: opt_S
        logical :: calculate_S=.true.
        ! local variables
    
        !finite difference calculation of H and S 
        character(len=300) :: info_str
        integer :: i_atom, i_cart, i_k_point, i_spin, i_basis, j_basis
        integer :: i_k_task
        integer :: i2, i3
        integer :: number_of_loops, atom_counter

        logical :: converged, enough_walltime
        
        integer ::  sc_iter_limit_backup
        character(len=40) :: output_level_backup
        real*8 :: sc_accuracy_potjump_backup, sc_accuracy_etot_backup
        real*8 :: sc_accuracy_eev_backup, sc_accuracy_rho_backup
        real*8 :: fermi_energy
        logical :: use_molecular_dynamics_backup,use_wf_extrapolation_backup
        integer :: postprocess_anyway_backup
    
        real*8,dimension(:,:),allocatable :: first_order_S_tmp_real
        real*8,dimension(:,:,:),allocatable :: first_order_H_tmp_real
        complex*16,dimension(:,:),allocatable :: first_order_S_tmp_cmplx
        complex*16,dimension(:,:,:),allocatable :: first_order_H_tmp_cmplx

        real*8, dimension(3,3) :: frac_coords_tmp
        real*8, dimension(3) :: trans_corr_vec

        character(*), parameter :: func = 'friction_calculate_H1_and_S1_num'
  
        if (present(opt_S)) calculate_S=opt_S

        if (myid.eq.0) then
            write (info_str,'(4X,A,A)') &
                "FRICTION Starting finite difference calculation of H1 and S1"
            call localorb_info(info_str,use_unit,'(A)')
        endif
        
        if(real_eigenvectors)then
           if (.not.allocated(KS_eigenvector_backup)) then
              allocate(KS_eigenvector_backup(n_basis,n_states,n_spin,n_k_points_task),stat=info)
              call check_allocation(info, 'KS_eigenvector_backup                ')
              KS_eigenvector_backup = 0.d0
           end if
        else
           if (.not.allocated(KS_eigenvector_complex_backup)) then
              allocate( KS_eigenvector_complex_backup(n_basis,n_states,n_spin,n_k_points_task),stat=info)
              call check_allocation(info, 'KS_eigenvector_complex_backup        ')
              KS_eigenvector_complex_backup = 0.d0
           end if
        end if
        if (.not.allocated(KS_eigenvector_backup))allocate(KS_eigenvector_backup(1,1,1,1)) ! allocate dummies
        if (.not.allocated(KS_eigenvector_complex_backup))allocate(KS_eigenvector_complex_backup(1,1,1,1))
        
        allocate( KS_eigenvalue_backup(n_states,n_spin,n_k_points),stat=info )
        call check_allocation(info, "KS_eigenvalue_backup",func)
        allocate( occ_numbers_backup(n_states,n_spin,n_k_points),stat=info )
        call check_allocation(info, "occ_numbers_backup",func)
      
        allocate(coords_backup(3,n_atoms),stat=info)
        call check_allocation(info, "coords_backup",func)
              
        !save old SCF data 
        KS_eigenvector_backup(:,:,:,:) = KS_eigenvector(:,:,:,:)
        KS_eigenvector_complex_backup(:,:,:,:) = KS_eigenvector_complex(:,:,:,:)
        occ_numbers_backup(:,:,:) = occ_numbers(:,:,:)
        KS_eigenvalue_backup(:,:,:) = KS_eigenvalue(:,:,:)
        coords_backup(:,:) = coords(:,:) 
        
        !finite difference workaround to solve strange error when 
        !displacement crosses fractional coordinate values of 0.5
        !call cart2frac(lattice_vectors, coords, frac_coords)
        !frac_coords_tmp
        !atom_counter = 0
        !trans_corr_vec = 0.0
        !do i_atom=1, n_atoms 
            !if (friction_atoms_list(i_atom)) then
                !atom_counter = atom_counter + 1
                !do i_cart=1, 3,1
                    !coords(i_cart,i_atom) = coords(i_cart,i_atom)+friction_numeric_disp
                    !call cart2frac(lattice_vector,coords,frac_coords_tmp)
                    !!check if frac_coords +- friction_numeric_disp crossed the 0.5 line 
                    !if (frac_coords_tmp(i_cart,i_atom)>0.50 .and. frac_coords(i_cart,i_atom)<0.50) then
                        !trans_corr_vec(i_cart) = -0.50
                    !endif
                    !coords(i_cart,i_atom) = coords(i_cart,i_atom)-2.0*friction_numeric_disp
                    !call cart2frac(lattice_vector,coords,frac_coords_tmp)
                    !if (frac_coords_tmp(i_cart,i_atom)>0.50 .and. frac_coords(i_cart,i_atom)<0.50) then
                        !trans_corr_vec(i_cart) = 0.50
                    !endif
                    !coords(i_cart,i_atom) = coords_backup(i_cart,i_atom)
                !end do
            !endif
        !enddo
        !do i_atom=1, n_atoms
            !coords(:,i_atom) = coords(:,i_atom) + trans_corr_vec(:)
            !coords_backup(:,i_atom) = coords_backup(:,i_atom) + trans_corr_vec(:)
        !enddo

        !TODO IF FOR COMPLEX MATRICES
        !TODO READ/OUTPUT MATRICES FOR COMPLEX MATRICES 
        !if (use_scalapack) then
            !if (n_periodic==0) then
              !allocate(first_order_S_tmp_real(1,mxld, mxcol))
              !allocate(first_order_H_tmp_real(1,mxld, mxcol,n_spin))
            !else 
              !allocate(first_order_S_tmp_real(n_cells_in_hamiltonian,mxld, mxcol))
              !allocate(first_order_H_tmp_real(n_cells_in_hamiltonian,mxld, mxcol,n_spin))
            !endif
        !else
        
        if (friction_use_complex_matrices) then
            if (n_periodic==0) then
              !allocate(first_order_S_tmp_real(1,n_basis,n_basis))
              !allocate(first_order_H_tmp_real(1,n_basis,n_basis,n_spin))
              allocate(first_order_S_tmp_real(1,n_basis*(n_basis+1)/2),stat=info)
              call check_allocation(info, "first_order_S_tmp_real",func)
              allocate(first_order_H_tmp_real(1,n_basis*(n_basis+1)/2,n_spin),stat=info)
              call check_allocation(info, "first_order_H_tmp_real",func)
            else 
              if (real_eigenvectors) then
                allocate(first_order_S_tmp_real(n_k_points_task,n_basis*(n_basis+1)/2),stat=info)
                call check_allocation(info, "first_order_S_tmp_real",func)
                allocate(first_order_H_tmp_real(n_k_points_task,n_basis*(n_basis+1)/2,n_spin),stat=info)
                call check_allocation(info, "first_order_H_tmp_real",func)
              else
                allocate(first_order_S_tmp_cmplx(n_k_points_task,n_basis*(n_basis+1)/2),stat=info)
                call check_allocation(info, "first_order_S_tmp_cmplx",func)
                allocate(first_order_H_tmp_cmplx(n_k_points_task,n_basis*(n_basis+1)/2,n_spin),stat=info)
                call check_allocation(info, "first_order_H_tmp_cmplx",func)
              endif
            endif
        else
            if (n_periodic==0) then
              allocate(first_order_S_tmp_real(1,n_basis*(n_basis+1)/2),stat=info)
              call check_allocation(info, "first_order_S_tmp_real",func)
              allocate(first_order_H_tmp_real(1,n_basis*(n_basis+1)/2,n_spin),stat=info)
              call check_allocation(info, "first_order_H_tmp_real",func)
            else 
              allocate(first_order_S_tmp_real(n_cells_in_hamiltonian,n_basis*(n_basis+1)/2),stat=info)
              call check_allocation(info, "first_order_S_tmp_real",func)
              allocate(first_order_H_tmp_real(n_cells_in_hamiltonian,n_basis*(n_basis+1)/2,n_spin),stat=info)
              call check_allocation(info, "first_order_H_tmp_real",func)
            endif
        endif
        if (.not.allocated(first_order_S_tmp_real)) allocate(first_order_S_tmp_real(1,1))
        if (.not.allocated(first_order_H_tmp_real)) allocate(first_order_H_tmp_real(1,1,n_spin))
        if (.not.allocated(first_order_S_tmp_cmplx)) allocate(first_order_S_tmp_cmplx(1,1))
        if (.not.allocated(first_order_H_tmp_cmplx)) allocate(first_order_H_tmp_cmplx(1,1,n_spin))
        !OVERWRITE SCF SETTINGS
        ini_linear_mixing = 0
        sc_iter_limit_backup = sc_iter_limit
        sc_iter_limit = friction_iter_limit
        sc_accuracy_potjump_backup = sc_accuracy_potjump
        sc_accuracy_potjump = friction_accuracy_potjump
        sc_accuracy_etot_backup = sc_accuracy_etot
        sc_accuracy_etot = friction_accuracy_etot
        sc_accuracy_rho_backup = sc_accuracy_rho
        sc_accuracy_rho = friction_accuracy_rho
        sc_accuracy_eev_backup = sc_accuracy_eev
        sc_accuracy_eev = friction_accuracy_eev

        output_level_backup = output_level
        output_level = 'MD_light' 
        use_forces = .false.
        postprocess_anyway_backup = postprocess_anyway
        postprocess_anyway = PP_ANYWAY_EVERYTHING
       
        !!+!+!wf extrapolation cannot deal with this
        !use_wf_extrapolation_backup = use_wf_extrapolation
        use_wf_extrapolation = .false.
        !!+!+!wf_extrapolation
        !use_molecular_dynamics_backup = use_molecular_dynamics
        use_molecular_dynamics = .false.
        
        fermi_energy = chemical_potential

        restart_write = .false.

        atom_counter = 0
        ATOM_LOOP: do i_atom=1, n_atoms 
            
            if (friction_atoms_list(i_atom)) then
                
                atom_counter = atom_counter + 1
                if (myid.eq.0) then
                    write (info_str,'(4X,A,I4,A,I4,A,I4)') &
                        "FRICTION Displacing atom ",atom_counter," of ",&
                        friction_n_active_atoms ," active atoms, index: ",i_atom
                    call localorb_info(info_str,use_unit,'(A)')
                endif
            
                CART_LOOP: do i_cart=1, 3,1
                  !MODIFY GEOMETRY HERE
                  if (myid.eq.0) then
                      if (i_cart .eq. 1) then
                          write (info_str,'(6X,A)') &
                              "FRICTION Calculating cartesian displacement in x direction"
                      else if (i_cart .eq. 2) then
                          write (info_str,'(6X,A)') &
                              "FRICTION Calculating cartesian displacement in y direction"
                      else
                          write (info_str,'(6X,A)') &
                              "FRICTION Calculating cartesian displacement in z direction"
                      endif
                      call localorb_info(info_str,use_unit,'(A)')
                  endif

                  !POSITIVE DISPLACEMENT
                  !update geometry
                  coords(i_cart, i_atom) = coords_backup(i_cart, i_atom) +& 
                  friction_numeric_disp
                  !call map_to_center_cell(coords(:,i_atom))
                  !save old KS_eigenvectors
                  !KS_eigenvector(:,:,:,:) = KS_eigenvector_backup(:,:,:,:)
                  !KS_eigenvector_complex(:,:,:,:) = KS_eigenvector_complex_backup(:,:,:,:)
                  !occ_numbers(:,:,:) = occ_numbers_backup(:,:,:)
                  !KS_eigenvalue(:,:,:) = KS_eigenvalue_backup(:,:,:)
                  
                  !SCF FOR FINITE DIFFERENCE

                  call reinitialize_scf(converged)
                  call scf_solver(converged, enough_walltime)
                  if (friction_use_complex_matrices) then
                      first_order_S_tmp_cmplx(:,:) = 0.0
                      first_order_H_tmp_cmplx(:,:,:) = 0.0
                      first_order_S_tmp_real(:,:) = 0.0
                      first_order_H_tmp_real(:,:,:) = 0.0
                      call friction_construct_complex_matrix(overlap_matrix, &
                        first_order_S_tmp_real(:,:),first_order_S_tmp_cmplx(:,:))
                      do i_spin=1, n_spin
                        call friction_construct_complex_matrix(hamiltonian(:,i_spin), &
                            first_order_H_tmp_real(:,:,i_spin),&
                            first_order_H_tmp_cmplx(:,:,i_spin))
                      enddo 
                      if (n_periodic==0 .or. real_eigenvectors) then
                        first_order_G(i_cart,atom_counter,:,:,1,:) = first_order_H_tmp_real(:,:,:)
                        if (calculate_S) first_order_S(i_cart,atom_counter,:,:,1) = first_order_S_tmp_real(:,:)
                      else
                        first_order_G_cmplx(i_cart,atom_counter,:,:,1,:) = first_order_H_tmp_cmplx(:,:,:)
                        if (calculate_S) first_order_S_cmplx(i_cart,atom_counter,:,:,1) = first_order_S_tmp_cmplx(:,:)
                      endif
                  else
                      !if (use_scalapack) then
                          !do i_spin=1, n_spin
                            !first_order_S_tmp_real
                          !enddo 
                          !TODO FIRST ASSURE THAT HAM AND OVLP ARE THERE IN SCALAPACK FORM
                          !TODO GENERATE A FUNCTION THAT builds S1 and H1 in SCALAPACK FORM
                          !call friction_construct_real_matrix_scalapack(overlap_matrix, &
                            !first_order_S_tmp_real(:,:))
                      !else 
                          first_order_S_tmp_real(:,:) = 0.0
                          first_order_H_tmp_real(:,:,:) = 0.0
                          call friction_construct_real_matrix(overlap_matrix, &
                            first_order_S_tmp_real(:,:))
                          do i_spin=1, n_spin
                            call friction_construct_real_matrix(hamiltonian(:,i_spin), &
                                first_order_H_tmp_real(:,:,i_spin))  
                          enddo 
                      !endif
                      if (calculate_S) then
                          first_order_S(i_cart,atom_counter,:,:,1) = first_order_S_tmp_real(:,:)
                      endif
                      first_order_G(i_cart,atom_counter,:,:,1,:) = first_order_H_tmp_real(:,:,:)
                  endif
                  !NEGATIVE DISPLACEMENT
                  !update geometry
                  coords(i_cart, i_atom) = coords_backup(i_cart, i_atom) -&
                     friction_numeric_disp
                  !call map_to_center_cell(coords(:,i_atom))
                  !save old KS_eigenvectors
                  !KS_eigenvector(:,:,:,:) = KS_eigenvector_backup(:,:,:,:)
                  !KS_eigenvector_complex(:,:,:,:) = KS_eigenvector_complex_backup(:,:,:,:)
                  !occ_numbers(:,:,:) = occ_numbers_backup(:,:,:)
                  !KS_eigenvalue(:,:,:) = KS_eigenvalue_backup(:,:,:)

                  !SCF FOR FINITE DIFFERENCE
                  call reinitialize_scf(converged)
                  call scf_solver(converged, enough_walltime)
                  
                  if (friction_use_complex_matrices) then
                      first_order_S_tmp_cmplx(:,:) = 0.0
                      first_order_H_tmp_cmplx(:,:,:) = 0.0
                      first_order_S_tmp_real(:,:) = 0.0
                      first_order_H_tmp_real(:,:,:) = 0.0
                      call friction_construct_complex_matrix(overlap_matrix, &
                        first_order_S_tmp_real(:,:),first_order_S_tmp_cmplx(:,:))
                      do i_spin=1, n_spin
                        call friction_construct_complex_matrix(hamiltonian(:,i_spin), &
                            first_order_H_tmp_real(:,:,i_spin),&
                            first_order_H_tmp_cmplx(:,:,i_spin))
                      enddo 
                      if (n_periodic==0 .or. real_eigenvectors) then
                          if (calculate_S) then
                              first_order_S(i_cart,atom_counter,:,:,1) = &
                                  first_order_S(i_cart,atom_counter,:,:,1) - &
                                  first_order_S_tmp_real(:,:)
                          endif
                          first_order_G(i_cart,atom_counter,:,:,1,:) = &
                              first_order_G(i_cart,atom_counter,:,:,1,:)-&
                              first_order_H_tmp_real(:,:,:)
                      else
                          first_order_G_cmplx(i_cart,atom_counter,:,:,1,:) = & 
                              first_order_G_cmplx(i_cart,atom_counter,:,:,1,:) - &
                              first_order_H_tmp_cmplx(:,:,:)
                          if (calculate_S) then
                              first_order_S_cmplx(i_cart,atom_counter,:,:,1) = & 
                                  first_order_S_cmplx(i_cart,atom_counter,:,:,1) - &
                                  first_order_S_tmp_cmplx(:,:)
                          endif
                      endif
                  else
                      !if (use_scalapack) then
                          !do i_spin=1, n_spin
                            !call friction_construct_real_matrix_scalapack(hamiltonian(:,i_spin), &
                                !first_order_H_tmp_real(:,:,i_spin))  
                          !enddo 
                          !call friction_construct_real_matrix_scalapack(overlap_matrix, &
                            !first_order_S_tmp_real(:,:))
                      !else
                          first_order_S_tmp_real(:,:) = 0.0
                          first_order_H_tmp_real(:,:,:) = 0.0
                          call friction_construct_real_matrix(overlap_matrix, &
                            first_order_S_tmp_real(:,:))
                          do i_spin=1, n_spin
                            call friction_construct_real_matrix(hamiltonian(:,i_spin), &
                                first_order_H_tmp_real(:,:,i_spin))  
                          enddo 
                      !endif
                      first_order_G(i_cart,atom_counter,:,:,1,:) = &
                          first_order_G(i_cart,atom_counter,:,:,1,:) - first_order_H_tmp_real(:,:,:)
                      if (calculate_S) then
                          first_order_S(i_cart,atom_counter,:,:,1) = &
                              first_order_S(i_cart,atom_counter,:,:,1) - first_order_S_tmp_real(:,:)
                      endif
                  endif
                  coords(i_cart, i_atom) = coords_backup(i_cart, i_atom)
                  !call map_to_center_cell(coords(:,i_atom))
                end do CART_LOOP
            endif
        end do ATOM_LOOP
        first_order_G(:,:,:,:,:,:) = first_order_G(:,:,:,:,:,:) &
            / (2.0*friction_numeric_disp)
        first_order_G_cmplx(:,:,:,:,:,:) = first_order_G_cmplx(:,:,:,:,:,:) &
            / (2.0*friction_numeric_disp)
        if (calculate_S) then
            first_order_S_cmplx(:,:,:,:,:) = first_order_S_cmplx(:,:,:,:,:) &
            / (2.0*friction_numeric_disp)
        first_order_S(:,:,:,:,:) = first_order_S(:,:,:,:,:) &
            / (2.0*friction_numeric_disp)
        endif
        if (allocated(first_order_S_tmp_real)) deallocate(first_order_S_tmp_real)
        if (allocated(first_order_H_tmp_real)) deallocate(first_order_H_tmp_real)
        if (allocated(first_order_S_tmp_cmplx)) deallocate(first_order_S_tmp_cmplx)
        if (allocated(first_order_H_tmp_cmplx)) deallocate(first_order_H_tmp_cmplx)
        !RESTORE electronic structure to equilibrium
        coords(:,:) = coords_backup(:,:)
        !call map_to_center_cell(coords(:,i_atom))
        call reinitialize_scf(converged)
        !call scf_solver(converged, enough_walltime)
        KS_eigenvector(:,:,:,:) = KS_eigenvector_backup(:,:,:,:)
        KS_eigenvector_complex(:,:,:,:) = KS_eigenvector_complex_backup(:,:,:,:)
        occ_numbers(:,:,:) = occ_numbers_backup(:,:,:)
        KS_eigenvalue(:,:,:) = KS_eigenvalue_backup(:,:,:)

        sc_iter_limit = sc_iter_limit_backup
        sc_accuracy_potjump = sc_accuracy_potjump_backup
        sc_accuracy_etot = sc_accuracy_etot_backup
        sc_accuracy_rho = sc_accuracy_rho_backup
        sc_accuracy_eev = sc_accuracy_eev_backup
        output_level = output_level_backup
        postprocess_anyway = postprocess_anyway_backup

        !use_wf_extrapolation = use_wf_extrapolation_backup
        !use_molecular_dynamics = use_molecular_dynamics_backup
        
        if (allocated(KS_eigenvector_backup)) deallocate(KS_eigenvector_backup)
        if (allocated(KS_eigenvector_complex_backup)) deallocate(KS_eigenvector_complex_backup)
        if (allocated(KS_eigenvalue_backup)) deallocate(KS_eigenvalue_backup)
        if (allocated(occ_numbers_backup)) deallocate(occ_numbers_backup)
        if (allocated(coords_backup)) deallocate(coords_backup)

    end subroutine friction_calculate_H1_and_S1_num
    
    subroutine friction_calculate_S1_DFPT(i_atom, i_cart, atom_counter, S_real, S_real_right)
    !  NAME
    !    friction_calculate_S1_DFPT
    !  SYNOPSIS
    !    calculates the one-sided nuclear derivative of the overlap matrix
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2016)
        
        implicit none

        !  ARGUMENTS

        ! imported variables
        integer, intent(in) :: i_atom
        integer, intent(in) :: i_cart
        integer, intent(inout) :: atom_counter 
        real*8, intent(out) :: S_real(:,:)
        real*8, intent(out) :: S_real_right(:,:)

        ! local variables
        integer :: i_k_point
        integer :: i_index_real, i_basis_1, i_basis_2
   
        !integer, dimension(:,:),allocatable :: basis_index
        
        character(len=300) :: info_str
        
        !allocate(S_real(n_basis,n_basis))
        !allocate(basis_index(n_basis,n_basis))
        
        !basis_index = 0
        !i_index_real = 0
        !do i_basis_2 = 1, n_basis, 1
          !do i_basis_1 = 1, i_basis_2, 1
            !i_index_real = i_index_real + 1
            !basis_index(i_basis_1,i_basis_2) = i_index_real
          !enddo
        !enddo
        
        S_real(:,:) = 0.d0
        if (friction_coupling_matrix_mode==2) then
            call integrate_first_order_S_right_reduce_memory(partition_tab, l_shell_max, S_real_right,i_atom,i_cart)
            do i_basis_1=1, n_basis
              do i_basis_2=1, n_basis
                S_real(i_basis_1,i_basis_2) = S_real_right(i_basis_1,i_basis_2)+ &
                    & S_real_right(i_basis_2,i_basis_1)
              enddo
            enddo
        else
            call integrate_first_order_S_reduce_memory(partition_tab, l_shell_max, S_real,i_atom,i_cart)
        endif

        !deallocate(basis_index)
        !deallocate(S_real)

    end subroutine friction_calculate_S1_DFPT
    
    subroutine friction_calculate_S1_DFPT_p1(i_atom, i_cart, i_q_point, atom_counter, S_complex, S_complex_right)
    !  NAME
    !    friction_calculate_S1_DFPT_p1
    !  SYNOPSIS
    !    calculates the one-sided nuclear derivative of the overlap matrix
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2016)
        
        implicit none

        !  ARGUMENTS

        ! imported variables
        integer, intent(in) :: i_atom
        integer, intent(in) :: i_cart
        integer, intent(in) :: i_q_point
        integer, intent(in) :: atom_counter 
        complex*16, intent(out) :: S_complex(:,:,:)
        complex*16, intent(out) :: S_complex_right(:,:,:)

        ! local variables
        integer :: i_index_real, i_basis_1, i_basis_2
        
       
        complex*16, allocatable :: S_real(:)
        complex*16, allocatable :: S_real_right(:)
        !integer, dimension(:,:),allocatable :: basis_index
        
        character(len=300) :: info_str
        
        allocate(S_real(n_hamiltonian_matrix_size))
        allocate(S_real_right(n_hamiltonian_matrix_size))
        !allocate(S_complex(n_basis,n_basis,n_k_points_task))
        !allocate(basis_index(n_basis,n_basis))
        
        !basis_index = 0
        !i_index_real = 0
        !do i_basis_2 = 1, n_basis, 1
          !do i_basis_1 = 1, i_basis_2, 1
            !i_index_real = i_index_real + 1
            !basis_index(i_basis_1,i_basis_2) = i_index_real
          !enddo
        !enddo
        S_complex = (0.d0,0.d0) 
        S_complex_right = (0.d0,0.d0) 
        if (friction_coupling_matrix_mode==2) then
            call integrate_first_order_S_right_phonon_reduce_memory(partition_tab,&
                l_shell_max,i_q_point, i_atom, i_cart, S_real_right) 
            call construct_first_order_S_phonon_reduce_memory(&
                S_real_right, S_complex_right) 

            do i_basis_1=1, n_basis
              do i_basis_2=1, n_basis
                S_complex(i_basis_1,i_basis_2,:) = S_complex_right(i_basis_1,i_basis_2,:)+ &
                   & conjg(S_complex_right(i_basis_2,i_basis_1,:))
              enddo
            enddo
        else
            call integrate_first_order_S_phonon_reduce_memory(partition_tab,&
                l_shell_max,i_q_point, i_atom, i_cart, S_real) 
            call construct_first_order_S_phonon_reduce_memory(&
                S_real, S_complex) 
        endif 

        !deallocate(basis_index)
        deallocate(S_real_right)
        deallocate(S_real)
        !deallocate(S_complex)

    end subroutine friction_calculate_S1_DFPT_p1
    
    subroutine friction_calculate_H1_and_S1_DFPT(i_q_point_in)
    !  NAME
    !    friction_calculate_H1_and_S1_DFPT
    !  SYNOPSIS
    !    calculates the H1 and S1 using DFPT 
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)
        
        implicit none

        ! imported variables
        integer, intent(inout), optional :: i_q_point_in
        ! local variables
    
        !finite difference calculation of H and S 
        character(len=300) :: info_str
        integer :: i_atom, i_cart, i_k_point, i_spin, i_basis, j_basis
        integer :: i_k_task, i_stat, i_state, i_q_point
        integer :: i_point
        integer :: i_index_real
        integer :: number_of_loops, atom_counter

        logical :: converged, enough_walltime
        logical :: first_iteration, below_it_limit
        
        integer ::  sc_iter_limit_backup
        character(len=40) :: output_level_backup
        integer :: postprocess_anyway_backup
        
        integer :: info 
        integer ::  max_occ_number(n_spin)
        real*8  ::  change_of_first_order_DM
      
        integer, dimension(:,:), allocatable :: basis_index

        !-----------------cluster--------------! 
        real*8, allocatable :: first_order_U(:,:)
        real*8, allocatable :: first_order_E(:)
        real*8, allocatable :: first_order_H(:,:)
        real*8, allocatable :: S_real(:,:)
        real*8, allocatable :: S_real_right(:,:)
        real*8, allocatable :: density_matrix(:,:) 
        real*8, allocatable :: energy_density_matrix(:,:) 
        real*8, allocatable :: first_order_density_matrix(:,:)
        real*8, allocatable :: first_order_energy_density_matrix(:,:)
        real*8, allocatable :: old_first_order_density_matrix(:,:)

        real*8, allocatable :: first_order_rho(:)
        real*8, allocatable :: first_order_potential(:)
        !-----------------periodic-------------!
        !DFPT_phonon
        real*8, allocatable        :: first_order_S_sparse(:,:,:)

        !reduce_memory
        complex*16, allocatable    :: first_order_U_complex(:,:,:)
        complex*16, allocatable    :: first_order_H_sparse(:)
        complex*16, allocatable    :: first_order_H_complex(:,:,:) 
        complex*16, allocatable    :: S_complex(:,:,:) 
        complex*16, allocatable    :: S_complex_right(:,:,:) 

        real*8, allocatable        ::  density_matrix_sparse(:)
        complex*16, allocatable    ::  first_order_density_matrix_sparse(:)
        complex*16, allocatable    ::  old_first_order_density_matrix_sparse(:)

        real*8, allocatable        ::  energy_density_matrix_sparse(:) 
        complex*16, allocatable    ::  first_order_energy_density_matrix_sparse(:)

        complex*16, allocatable    :: first_order_rho_complex(:)
        complex*16, allocatable    :: first_order_potential_complex(:)

        real*8, allocatable        :: first_order_rho_Re(:)
        real*8, allocatable        :: first_order_rho_Im(:)
        real*8, allocatable        :: first_order_potential_Re(:)
        real*8, allocatable        :: first_order_potential_Im(:)

        complex*16, allocatable    :: rho_free_gradient(:)
        complex*16, allocatable    :: v_free_gradient(:)

        complex*16, allocatable    :: hellman_feynman_dynamical_matrix_delta_part(:,:,:,:)

        character(*), parameter :: func = 'friction_calculate_H1_and_S1_DFPT'

        DFPT_mixing = 0.20d0

        allocate(basis_index(n_basis,n_basis))
        basis_index = 0
        i_index_real = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis, 1
            i_index_real = i_index_real + 1
            basis_index(i_basis,j_basis) = i_index_real
          enddo
        enddo
        
        if (present(i_q_point_in)) then
            i_q_point = i_q_point_in
        else
            i_q_point = 1
        endif

        if (myid.eq.0) then
            write (info_str,'(4X,A,A)') &
                "FRICTION Starting DFPT calculation of H1 and S1"
            call localorb_info(info_str,use_unit,'(A)')
        endif
             
        if (myid.eq.0) then
          write(info_str,'(A)') &
              'This routine currently only calculates first order S for the cluster case! &
              & First order H is still calculated with finite difference'
            call localorb_info(info_str,use_unit,'(A)')
        endif

        if(real_eigenvectors)then
           if (.not.allocated(KS_eigenvector_backup)) then
              allocate(KS_eigenvector_backup(n_basis,n_states,n_spin,n_k_points_task),stat=info)
              call check_allocation(info, 'KS_eigenvector_backup                ')
              KS_eigenvector_backup = 0.d0
           end if
        else
           if (.not.allocated(KS_eigenvector_complex_backup)) then
              allocate( KS_eigenvector_complex_backup(n_basis,n_states,n_spin,n_k_points_task),stat=info)
              call check_allocation(info, 'KS_eigenvector_complex_backup        ')
              KS_eigenvector_complex_backup = 0.d0
           end if
        end if
        if (.not.allocated(KS_eigenvector_backup))allocate(KS_eigenvector_backup(1,1,1,1)) ! allocate dummies
        if (.not.allocated(KS_eigenvector_complex_backup))allocate(KS_eigenvector_complex_backup(1,1,1,1))
        
        allocate( KS_eigenvalue_backup(n_states,n_spin,n_k_points),stat=info )
        call check_allocation(info, "KS_eigenvalue_backup",func)
        allocate( occ_numbers_backup(n_states,n_spin,n_k_points),stat=info )
        call check_allocation(info, "occ_numbers_backup",func)
      
        !save old SCF data 
        KS_eigenvector_backup(:,:,:,:) = KS_eigenvector(:,:,:,:)
        KS_eigenvector_complex_backup(:,:,:,:) = KS_eigenvector_complex(:,:,:,:)
        occ_numbers_backup(:,:,:) = occ_numbers(:,:,:)
        KS_eigenvalue_backup(:,:,:) = KS_eigenvalue(:,:,:)

        !-------------allocate--------------!
        if (n_periodic.eq.0) then
           !cluster
           allocate(first_order_density_matrix(n_basis,n_basis))
           allocate(first_order_energy_density_matrix(n_basis,n_basis))
           allocate(old_first_order_density_matrix(n_basis,n_basis))
           allocate(first_order_rho(n_full_points))
           allocate(first_order_potential(n_full_points))
          
           allocate(first_order_U(n_basis,n_basis))
           allocate(S_real(n_basis,n_basis))
           allocate(S_real_right(n_basis,n_basis))
           allocate(first_order_E(n_basis))
           allocate(first_order_H(n_basis,n_basis))

           allocate(density_matrix(n_basis, n_basis)) 
           allocate(energy_density_matrix(n_basis, n_basis))
           !periodic dummies
           allocate(first_order_rho_complex(1))
           allocate(first_order_potential_complex(1))

           allocate(first_order_rho_Re(1))
           allocate(first_order_rho_Im(1))
           allocate(first_order_potential_Re(1))
           allocate(first_order_potential_Im(1))

           allocate(rho_free_gradient(1))
           allocate(v_free_gradient(1))
           
           allocate(first_order_H_sparse(1))
           allocate(first_order_H_complex(1,1,1))
           allocate(first_order_U_complex(1,1,1))

           allocate(density_matrix_sparse(1))
           allocate(first_order_density_matrix_sparse(1))
           allocate(old_first_order_density_matrix_sparse(1))
           allocate(energy_density_matrix_sparse(1))
           allocate(first_order_energy_density_matrix_sparse(1))
        else
           !cluster dummies
           allocate(first_order_density_matrix(1,1))
           allocate(first_order_energy_density_matrix(1,1))
           allocate(old_first_order_density_matrix(1,1))
           allocate(first_order_rho(1))
           allocate(first_order_potential(1))
           allocate(first_order_U(1,1))
           allocate(first_order_E(1))
           allocate(first_order_H(1,1))
           allocate(density_matrix(1,1)) 
           allocate(energy_density_matrix(1,1))
           !periodic
           allocate(first_order_rho_complex(n_full_points))
           allocate(first_order_potential_complex(n_full_points))

           allocate(first_order_rho_Re(n_full_points))
           allocate(first_order_rho_Im(n_full_points))
           allocate(first_order_potential_Re(n_full_points))
           allocate(first_order_potential_Im(n_full_points))

           allocate(rho_free_gradient(n_full_points))
           allocate(v_free_gradient(n_full_points))
           
           allocate(first_order_H_sparse(n_hamiltonian_matrix_size))
           allocate(first_order_H_complex(n_basis, n_basis,n_k_points_task))
           allocate(S_complex(n_basis, n_basis,n_k_points_task))
           allocate(S_complex_right(n_basis, n_basis,n_k_points_task))
           allocate(first_order_U_complex(n_basis,n_basis,n_k_points_task))

           allocate(density_matrix_sparse(n_hamiltonian_matrix_size))
           allocate(first_order_density_matrix_sparse(n_hamiltonian_matrix_size))
           allocate(old_first_order_density_matrix_sparse(n_hamiltonian_matrix_size))
           allocate(energy_density_matrix_sparse(n_hamiltonian_matrix_size))
           allocate(first_order_energy_density_matrix_sparse(n_hamiltonian_matrix_size))
           allocate(hellman_feynman_dynamical_matrix_delta_part(3,n_atoms,3,n_atoms))

           !phonon
           allocate(first_order_S_sparse(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size))

        endif

        !---------------------------------!
        if (n_periodic .eq. 0) then
            !find the max_occ_number
            do i_spin = 1, n_spin, 1
                max_occ_number(i_spin) = 0
                do i_state = n_states, 1, -1
                 if (dabs(occ_numbers(i_state,i_spin,1)).gt.0.d0) then
                  max_occ_number(i_spin) = i_state
                  exit
                 endif
                enddo
            enddo
            call evaluate_zero_order_DM_reduce_memory( & 
                 KS_eigenvector, occ_numbers, max_occ_number, density_matrix)
        else
           if(packed_matrix_format.ne.PM_index) then
           call aims_stop('shanghui only use sparse matrix for DFPT_phonon_reduce_memory', & 
                          'cpscf_solver_phonon_reduce_memory')
           endif
            !periodic
        endif 
      
        !non_memory_reduce versions 
        !DFPT_supercell phonon first_order_S
        !if (n_periodic>0) then
            !call integrate_first_order_S_p1(partition_tab, l_shell_max, &
            !first_order_S_sparse)
        !else
            !non-sparse
            !call integrate_first_order_S(partition_tab, l_shell_max, first_order_S)
        !endif

        atom_counter = 0
        ATOM_LOOP: do i_atom=1, n_atoms 
            
            if (friction_atoms_list(i_atom)) then
                
                atom_counter = atom_counter + 1
                CART_LOOP: do i_cart=1, 3,1
                    if (myid.eq.0) then
                        write(info_str,'(2X,A,1X,I4,5X,A,1X,I4)') 'FRICTION CPSCF working for i_atom =',&
                            & friction_index_list(atom_counter),'i_cart =',i_cart
                        call localorb_info(info_str, use_unit,'(A)')
                    endif

                    !CPSCF calculate first_order_S
                    if (myid.eq.0) then
                        write (info_str,'(4X,A,A)') &
                            "FRICTION Calculate analytical first order S1"
                        call localorb_info(info_str,use_unit,'(A)')
                    endif
                    if (n_periodic .eq. 0) then
                        call friction_calculate_S1_DFPT(i_atom, i_cart, atom_counter, S_real,&
                            & S_real_right)
                    else
                        call friction_calculate_S1_DFPT_p1(i_atom, i_cart, i_q_point, atom_counter,S_complex,&
                            & S_complex_right)
                    endif

                    !Calculate first order H
                    first_order_H_complex=0.0d0
                    first_order_U_complex=0.0d0 
                    first_order_H=0.0d0
                    first_order_U=0.0d0 
                    first_order_E=0.0d0
                
                    !!evaluate first order DM as start guess
                    !if (n_periodic .eq. 0) then
                        !call evaluate_first_order_DM_reduce_memory(S_real,  &
                            !KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number,   &
                            !first_order_U,old_first_order_density_matrix & 
                            !) 
                    !else
                        !!periodic
                        !call  integrate_free_atom_sum_gradient_phonon_reduce_memory &
                              !(partition_tab, i_q_point, i_atom, i_cart, rho_free_gradient,v_free_gradient) 

                        !call  evaluate_zero_order_DM_phonon_reduce_memory &
                              !(KS_eigenvector, KS_eigenvector_complex, occ_numbers,density_matrix_sparse)
                        !first_order_U_complex = (0.0d0, 0.0d0)
                        !call evaluate_first_order_DM_phonon_reduce_memory(S_complex,  &
                             !KS_eigenvector_complex, KS_eigenvalue, occ_numbers,   &
                             !first_order_U_complex,old_first_order_density_matrix_sparse)
                    !endif

                    converged = .false.
                    below_it_limit = .true.
                    number_of_loops = 0

                    !CPSCF calculate first_order_S
                    !if (myid.eq.0) then
                        !write (info_str,'(4X,A,A)') &
                            !"FRICTION Calculate analytical first order H1"
                        !call localorb_info(info_str,use_unit,'(A)')
                    !endif
                    !!START SCF
                    !SCF_LOOP: do while ( (.not.converged) .and. below_it_limit )
                        !number_of_loops = number_of_loops + 1

                        !write(info_str, '(4X,A,1X,I4)') "Begin CP-self-consistency iteration #", number_of_loops
                        !call localorb_info(info_str, use_unit,'(A)')
                        !if (n_periodic.eq.0) then
                            
                            !!evaluate DM
                            !call evaluate_first_order_DM_reduce_memory(S_real, &
                                !KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number, &
                                !first_order_U,first_order_density_matrix)
                            
                            !!mix DM
                            !change_of_first_order_DM =0.0d0         
                            !do i_basis = 1,n_basis
                                !do j_basis = 1,n_basis
                                    !change_of_first_order_DM =                &
                                    !max( change_of_first_order_DM,             &
                                    !dabs(first_order_density_matrix(i_basis,j_basis)  &
                                    !- old_first_order_density_matrix(i_basis,j_basis)) )
                             
                                    !first_order_density_matrix(i_basis,j_basis) =       &
                                    !(1.0d0-DFPT_mixing)*old_first_order_density_matrix(i_basis,j_basis)+  &
                                    !DFPT_mixing*first_order_density_matrix(i_basis,j_basis)
                         
                                    !old_first_order_density_matrix(i_basis,j_basis) =   &
                                    !first_order_density_matrix(i_basis,j_basis)
                                !enddo
                            !enddo
                            !!integrate first order rho
                            !call integrate_first_order_rho_reduce_memory(partition_tab, l_shell_max,  &
                                !KS_eigenvector(:,:,1,1),density_matrix,first_order_density_matrix,max_occ_number, &
                                !first_order_rho,i_atom,i_cart)

                            !!calculate H1
                            !call update_hartree_potential_p2_shanghui &
                                !( hartree_partition_tab,first_order_rho(1:n_full_points),& 
                                !delta_v_hartree_part_at_zero, &
                                !delta_v_hartree_deriv_l0_at_zero, &
                                !multipole_moments, multipole_radius_sq, &
                                !l_hartree_max_far_distance, &
                                !outer_potential_radius )

                            !call sum_up_whole_potential_p2_shanghui &
                                    !( delta_v_hartree_part_at_zero, &
                                    !delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                                    !partition_tab, first_order_rho(1:n_full_points), &
                                    !first_order_potential(1:n_full_points),  & !<--------get first_order_DM_potential
                                    !.false., multipole_radius_sq, &
                                    !l_hartree_max_far_distance, &
                                    !outer_potential_radius)

                            !call  integrate_first_order_H_reduce_memory &
                                 !(hartree_potential,first_order_potential, rho, rho_gradient,&
                                  !first_order_rho, & 
                                  !partition_tab, l_shell_max, &
                                  !first_order_H, i_atom, i_cart &
                                 !)
                        
                            !!calculate U1
                            !call evaluate_first_order_U_reduce_memory(first_order_H, S_real,  &
                                 !KS_eigenvector, KS_eigenvalue, occ_numbers,  max_occ_number, &
                                 !first_order_U,first_order_E)
                        !else
                            !!periodic
                            
                            !!evaluate DM
                            !call evaluate_first_order_DM_phonon_reduce_memory(S_complex,  &
                                 !KS_eigenvector_complex, KS_eigenvalue, occ_numbers,   &
                                 !first_order_U_complex,first_order_density_matrix_sparse)
                            !!mix DM
                            !change_of_first_order_DM =0.0d0         

                            !do i_basis = 1, n_hamiltonian_matrix_size - 1 
                                !change_of_first_order_DM =                &
                                !max( change_of_first_order_DM,             &
                                !dabs( dble(first_order_density_matrix_sparse(i_basis)  &
                                   !- old_first_order_density_matrix_sparse(i_basis))) )
                                 
                                !first_order_density_matrix_sparse(i_basis) =       &
                                !(1.0d0-DFPT_mixing)*old_first_order_density_matrix_sparse(i_basis)+  &
                                !DFPT_mixing*first_order_density_matrix_sparse(i_basis)
                             
                                !old_first_order_density_matrix_sparse(i_basis) =   &
                                !first_order_density_matrix_sparse(i_basis)
                             !enddo

                            !!integrate first order rho
                            !call integrate_first_order_rho_phonon_reduce_memory(partition_tab, l_shell_max,  &
                                 !density_matrix_sparse,first_order_density_matrix_sparse, &
                                 !i_q_point, i_atom, i_cart, &
                                 !first_order_rho_complex)

                            !first_order_rho_complex(1:n_full_points) = &
                                !first_order_rho_complex(1:n_full_points) +&
                                !rho_free_gradient(1:n_full_points)

                            !do i_point =1 ,n_full_points
                                !first_order_rho_Re(i_point)=  &           
                                !dble(first_order_rho_complex(i_point)) 

                                !first_order_rho_Im(i_point)=  &           
                                !dimag(first_order_rho_complex(i_point))  
                            !enddo  

                            !!calculate H1
                            !call update_hartree_potential_shanghui_phonon_reduce_memory &
                                !(hartree_partition_tab,first_order_rho_Re(1:n_full_points),& 
                                 !delta_v_hartree_part_at_zero, &
                                 !delta_v_hartree_deriv_l0_at_zero, &
                                 !multipole_moments, multipole_radius_sq, &
                                 !l_hartree_max_far_distance, &
                                 !outer_potential_radius )

                            !call sum_up_whole_potential_shanghui_phonon_reduce_memory &
                                !(delta_v_hartree_part_at_zero, &
                                 !delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                                 !partition_tab, first_order_rho_Re(1:n_full_points), &
                                 !first_order_potential_Re(1:n_full_points),  & 
                                 !.false., multipole_radius_sq, &
                                 !l_hartree_max_far_distance, &
                                 !outer_potential_radius, & 
                                 !hellman_feynman_dynamical_matrix_delta_part) 

                            !call update_hartree_potential_shanghui_phonon_reduce_memory &
                                !(hartree_partition_tab,first_order_rho_Im(1:n_full_points),& 
                                 !delta_v_hartree_part_at_zero, &
                                 !delta_v_hartree_deriv_l0_at_zero, &
                                 !multipole_moments, multipole_radius_sq, &
                                 !l_hartree_max_far_distance, &
                                 !outer_potential_radius )

                            !call sum_up_whole_potential_shanghui_phonon_reduce_memory &
                                !(delta_v_hartree_part_at_zero, &
                                 !delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                                 !partition_tab, first_order_rho_Im(1:n_full_points), &
                                 !first_order_potential_Im(1:n_full_points),  & 
                                 !.false., multipole_radius_sq, &
                                 !l_hartree_max_far_distance, &
                                 !outer_potential_radius, & 
                                 !hellman_feynman_dynamical_matrix_delta_part) 

                            !do i_point =1 ,n_full_points
                               !first_order_potential_complex(i_point)=  &           
                               !dcmplx(first_order_potential_Re(i_point),first_order_potential_Im(i_point)) 
                            !enddo  
                     
                            !first_order_potential_complex(1:n_full_points)= &      !@
                                !first_order_potential_complex(1:n_full_points)   &     !@
                                !+(-v_free_gradient(1:n_full_points))           !@

                            !call integrate_first_order_rho_phonon_reduce_memory(partition_tab, l_shell_max,  &
                                   !density_matrix_sparse,first_order_density_matrix_sparse, &
                                   !i_q_point, i_atom, i_cart, &
                                   !first_order_rho_complex)

                            !call integrate_first_order_H_phonon_reduce_memory &
                                   !(hartree_potential,first_order_potential_complex, & 
                                   !rho, rho_gradient, first_order_rho_complex, &
                                   !partition_tab, l_shell_max,    &
                                   !i_q_point, i_atom, i_cart,    &
                                   !first_order_H_sparse)

                            !i_k_task = 0
                            !do i_k_point = 1,n_k_points, 1
                              !if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
                                !i_k_task = i_k_task + 1   
                     
                                !call construct_first_order_matrix_phonon_reduce_memory(first_order_H_sparse, &
                                    !first_order_H_complex(1:n_basis,1:n_basis,i_k_task), i_k_point)
                                !!calculate U1
                                !call evaluate_first_order_U_phonon_reduce_memory(& 
                                  !first_order_H_complex(:,:,i_k_task),& 
                                  !S_complex(:,:,i_k_task), &
                                  !KS_eigenvector_complex(:,:,:,i_k_task), KS_eigenvalue(:,:,i_k_point), &  
                                  !occ_numbers(:,:,i_k_point),  &
                                  !first_order_U_complex(:,:,i_k_task))
                              !endif
                            !enddo
                        !endif

                        !write(info_str, '(4X,A,1X,F16.8)') "CPSCF Change of first_order_DM ", change_of_first_order_DM
                        !call localorb_info(info_str, use_unit,'(A)')
                        !!check convergence
                        !converged = (change_of_first_order_DM.lt.1.0d-4).and.(number_of_loops.gt.1) 

                        !if (converged) then
                          !write(info_str,'(2X,A)') "CP-self-consistency cycle converged."
                          !call localorb_info(info_str,use_unit,'(A)',OL_norm)

                        
                        !else if (number_of_loops.ge.friction_iter_limit) then
                            !below_it_limit = .false.
                        !endif

                    !enddo SCF_LOOP
                      !!copy H to G
                      if (n_periodic.eq.0) then
                        if (friction_coupling_matrix_mode==2) then
                            S_real(:,:) = S_real_right(:,:)
                        endif
                          do j_basis=1, n_basis
                            do i_basis=1, j_basis
                              !first_order_G(i_cart,atom_counter,1,basis_index(i_basis,j_basis):,1,1) = &
                                  !& first_order_H(i_basis,j_basis)
                              first_order_S(i_cart,atom_counter,1,basis_index(i_basis,j_basis),1) = &
                                  & S_real(i_basis,j_basis) 
                            enddo
                          enddo
                      else
                        !periodic
                        if (friction_coupling_matrix_mode==2) then
                            S_complex(:,:,:) = S_complex_right(:,:,:)
                        endif
                        do j_basis = 1, n_basis
                          do i_basis = 1, j_basis
                            first_order_S_cmplx(i_cart,atom_counter,:,&
                                basis_index(i_basis,j_basis),1) = S_complex(i_basis,j_basis,:) 
                            !!n_spin = 1
                            first_order_G_cmplx(i_cart,atom_counter,:,&
                                basis_index(i_basis,j_basis),1,1) = (0.d0,0.d0)
                            !first_order_G_cmplx(i_cart,atom_counter,:,&
                                !basis_index(i_basis,j_basis),1,1) = first_order_H_complex(i_basis,j_basis,:) 
                          enddo
                        enddo
                      endif
                end do CART_LOOP
            endif
        end do ATOM_LOOP

        !DEALLOCATE WHATEVER WE DONT REALLY NEED
        !cluster
        if (allocated(first_order_rho)) deallocate(first_order_rho)
        if (allocated(first_order_potential)) deallocate(first_order_potential)
        if (allocated(first_order_U)) deallocate(first_order_U)
        if (allocated(S_real)) deallocate(S_real)
        if (allocated(S_real_right)) deallocate(S_real_right)
        if (allocated(first_order_E)) deallocate(first_order_E)
        if (allocated(density_matrix)) deallocate(density_matrix)
        if (allocated(energy_density_matrix)) deallocate(energy_density_matrix)
        if (allocated(first_order_density_matrix)) deallocate(first_order_density_matrix)
        if (allocated(first_order_energy_density_matrix)) deallocate(first_order_energy_density_matrix)
        if (allocated(old_first_order_density_matrix)) deallocate(old_first_order_density_matrix)

        !periodic
        if (allocated(first_order_rho_complex)) deallocate(first_order_rho_complex)
        if (allocated(first_order_potential_complex)) deallocate(first_order_potential_complex)
        if (allocated(first_order_rho_Re)) deallocate(first_order_rho_Re)
        if (allocated(first_order_rho_Im)) deallocate(first_order_rho_Im)
        if (allocated(first_order_potential_Re)) deallocate(first_order_potential_Re)
        if (allocated(first_order_potential_Im)) deallocate(first_order_potential_Im)

        if (allocated(rho_free_gradient)) deallocate(rho_free_gradient)
        if (allocated(v_free_gradient)) deallocate(v_free_gradient)

        if (allocated(S_complex_right)) deallocate(S_complex_right)
        if (allocated(S_complex)) deallocate(S_complex)
        if (allocated(first_order_H_sparse)) deallocate(first_order_H_sparse)
        if (allocated(first_order_H_complex)) deallocate(first_order_H_complex)
 
        if (allocated(first_order_U_complex)) deallocate(first_order_U_complex)

        if (allocated(density_matrix_sparse)) deallocate(density_matrix_sparse)
        if (allocated(first_order_density_matrix_sparse)) deallocate(first_order_density_matrix_sparse)
        if (allocated(old_first_order_density_matrix_sparse)) deallocate(old_first_order_density_matrix_sparse)
        if (allocated(energy_density_matrix_sparse)) deallocate(energy_density_matrix_sparse)
        if (allocated(first_order_energy_density_matrix_sparse)) deallocate(first_order_energy_density_matrix_sparse)
        if (allocated(hellman_feynman_dynamical_matrix_delta_part)) deallocate(hellman_feynman_dynamical_matrix_delta_part)

        !RESTORE electronic structure to equilibrium
        !coords(:,:) = coords_backup(:,:)
        !call reinitialize_scf(converged)
        KS_eigenvector(:,:,:,:) = KS_eigenvector_backup(:,:,:,:)
        KS_eigenvector_complex(:,:,:,:) = KS_eigenvector_complex_backup(:,:,:,:)
        occ_numbers(:,:,:) = occ_numbers_backup(:,:,:)
        KS_eigenvalue(:,:,:) = KS_eigenvalue_backup(:,:,:)
        
        if (allocated(KS_eigenvector_backup)) deallocate(KS_eigenvector_backup)
        if (allocated(KS_eigenvector_complex_backup)) deallocate(KS_eigenvector_complex_backup)
        if (allocated(KS_eigenvalue_backup)) deallocate(KS_eigenvalue_backup)
        if (allocated(occ_numbers_backup)) deallocate(occ_numbers_backup)
       
        !phonon
        if (allocated(first_order_S_sparse)) deallocate(first_order_S_sparse)

        deallocate(basis_index)
        
        !H1 doesnt work yet, we still have to calculate with finite-difference
        call friction_calculate_H1_and_S1_num(.false.)

    end subroutine friction_calculate_H1_and_S1_DFPT

    !subroutine friction_construct_complex_matrix(matrix_sparse, matrix_real, matrix_complex)
    !!  NAME
    !!    friction_construct_complex_matrix
    !!  SYNOPSIS
    !!  constructs the complex k-dependent versions of overlap or hamiltonian 
    !!  AUTHOR
    !!   Reinhard J. Maurer, Yale University (2015)
        
        !implicit none

        !real*8, intent(in) :: matrix_sparse(:)
        !real*8, intent(OUT) :: matrix_real(:,:)
        !complex*16, intent(OUT) :: matrix_complex(:,:)
        !! imported variables

        !! local variables
        !integer :: i_k_point, i_k_task,i_index_real
        !integer :: i_cell, i_place, i_basis_1, i_basis_2
        !integer :: k_cell, k_atom, k_cell_new, k_center_new
        
        !real*8, dimension(:,:,:), allocatable :: work
        !complex*16, dimension(:,:,:), allocatable :: work_cmplx
        !integer, dimension(:,:),allocatable :: basis_index
        
        !allocate(basis_index(n_basis,n_basis))
        
        !basis_index = 0
        !i_index_real = 0
        !do i_basis_2 = 1, n_basis, 1
          !do i_basis_1 = 1, i_basis_2, 1
            !i_index_real = i_index_real + 1
            !basis_index(i_basis_1,i_basis_2) = i_index_real
          !enddo
        !enddo
        
        !matrix_real= 0.0d0
        !matrix_complex= (0.0d0,0.0d0)
        !if (packed_matrix_format .eq. PM_none) then
          !do i_basis_2 = 1, n_basis
            !do i_basis_1 = 1, i_basis_2
              !matrix_real(1,basis_index(i_basis_1,i_basis_2)) = matrix_sparse(basis_index(i_basis_1,i_basis_2))
              !!matrix_real(1,i_basis_2,i_basis_1) = matrix_real(1,i_basis_1,i_basis_2)
            !enddo
          !enddo
        !else
          !!finite difference calculation of H and S
          !matrix_real= 0.0d0
          !matrix_complex= (0.0d0,0.0d0)
          !if (real_eigenvectors) then
            !allocate(work(n_k_points_task,n_basis,n_basis)) 
            !i_k_task = 0
            !do i_k_point = 1,n_k_points, 1
                !if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
                    !i_k_task = i_k_task + 1
        
                   !do i_cell = 1,n_cells_in_hamiltonian-1
                     !do i_basis_1 = 1, n_basis
                       !if( index_hamiltonian(1,i_cell, i_basis_1) > 0 )then
                         !do i_place = index_hamiltonian(1,i_cell, i_basis_1), &
                                      !index_hamiltonian(2,i_cell, i_basis_1)
                           !i_basis_2 =  column_index_hamiltonian(i_place)

                           !!matrix_real(i_k_task,i_basis_2, i_basis_1) =  &
                                   !!matrix_real(i_k_task,i_basis_2, i_basis_1) + &
                                   !!dble(k_phase(i_cell,i_k_point))* &
                                   !!matrix_sparse(i_place)
                           !work(i_k_task,i_basis_1, i_basis_2) =  &
                                   !work(i_k_task,i_basis_1, i_basis_2) + &
                                   !dble(k_phase(i_cell,i_k_point))* &
                                   !matrix_sparse(i_place)
                           !if (i_basis_1 .ne.i_basis_2) then
                             !!matrix_real(i_k_task,i_basis_1, i_basis_2) =  &
                                 !!matrix_real(i_k_task,i_basis_1, i_basis_2) + &
                                 !!dble(k_phase(i_cell,i_k_point))* &
                                 !!matrix_sparse(i_place)
                             !work(i_k_task,i_basis_2, i_basis_1) =  &
                                 !work(i_k_task,i_basis_2, i_basis_1) + &
                                 !dble(k_phase(i_cell,i_k_point))* &
                                 !matrix_sparse(i_place)
                           !endif
                         !enddo !i_place
                       !endif !index_hamiltonian
                     !enddo ! i_basis_1
                   !enddo ! i_cell
                    !!build matrix_complex
                !endif
            !enddo
            !do i_basis_2 = 1, n_basis
              !do i_basis_1 = 1, i_basis_2
                !matrix_real(:,basis_index(i_basis_1,i_basis_2)) = work(:,i_basis_1,i_basis_2)
              !enddo
            !enddo
          !else
            !allocate(work_cmplx(n_k_points_task,n_basis,n_basis)) 
            !i_k_task = 0
            !do i_k_point = 1,n_k_points, 1
                !if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
                    !i_k_task = i_k_task + 1   

                   !do i_cell = 1,n_cells_in_hamiltonian-1
                     !do i_basis_1 = 1, n_basis
                       !if( index_hamiltonian(1,i_cell, i_basis_1) > 0 )then
                         !do i_place = index_hamiltonian(1,i_cell, i_basis_1), &
                                      !index_hamiltonian(2,i_cell, i_basis_1)
                           !i_basis_2 =  column_index_hamiltonian(i_place)

                           !!matrix_complex(i_k_task, i_basis_2, i_basis_1) =  &
                                   !!matrix_complex(i_k_task, i_basis_2, i_basis_1) + &
                                   !!k_phase(i_cell,i_k_point)* &
                                   !!matrix_sparse(i_place)
                           !!work_cmplx(i_k_task, i_basis_2, i_basis_1) =  &
                                   !!work_cmplx(i_k_task, i_basis_2, i_basis_1) + &
                                   !!k_phase(i_cell,i_k_point)* &
                                   !!matrix_sparse(i_place)
                           !work_cmplx(i_k_task, i_basis_1, i_basis_2) =  &
                                   !work_cmplx(i_k_task, i_basis_1, i_basis_2) + &
                                   !k_phase(i_cell,i_k_point)* &
                                   !matrix_sparse(i_place)
                           !if (i_basis_1 .ne.i_basis_2) then
                             !!matrix_complex(i_k_task, i_basis_1, i_basis_2) =  &
                                 !!matrix_complex(i_k_task, i_basis_1, i_basis_2) + &
                                 !!dconjg(k_phase(i_cell,i_k_point))* &
                                 !!matrix_sparse(i_place)
                             !!work_cmplx(i_k_task, i_basis_1, i_basis_2) =  &
                                 !!work_cmplx(i_k_task, i_basis_1, i_basis_2) + &
                                 !!dconjg(k_phase(i_cell,i_k_point))* &
                                 !!matrix_sparse(i_place)
                             !work_cmplx(i_k_task, i_basis_2, i_basis_1) =  &
                                 !work_cmplx(i_k_task, i_basis_2, i_basis_1) + &
                                 !dconjg(k_phase(i_cell,i_k_point))* &
                                 !matrix_sparse(i_place)

                           !endif
                         !enddo !i_place
                       !endif !index_hamiltonian
                     !enddo ! i_basis_1
                   !enddo ! i_cell
                    !!build matrix_complex
                !endif
            !enddo
            !do i_basis_2 = 1, n_basis
              !do i_basis_1 = 1, i_basis_2
                !matrix_complex(:,basis_index(i_basis_1,i_basis_2)) = work_cmplx(:,i_basis_1,i_basis_2)
              !enddo
            !enddo
          !endif
        !endif

        !if (allocated(work)) deallocate(work)
        !if (allocated(work_cmplx)) deallocate(work_cmplx)
        !if (allocated(basis_index)) deallocate(basis_index)

    !end subroutine friction_construct_complex_matrix
    
    !subroutine friction_construct_complex_matrix2(matrix_sparse, matrix_real, matrix_complex)
    !!  NAME
    !!    friction_construct_real_matrix
    !!  SYNOPSIS
    !!  constructs the dense real space versions of overlap or hamiltonian 
    !!  AUTHOR
    !!   Reinhard J. Maurer, Yale University (2015)
        
        !implicit none

        !real*8, intent(in) :: matrix_sparse(:)
        !!real*8, intent(OUT) :: matrix_real(n_basis,n_cells_in_hamiltonian,n_basis)
        !real*8, intent(out) :: matrix_real(:,:)
        !complex*16, intent(OUT) :: matrix_complex(:,:)
        !! imported variables

        !! local variables
        !integer :: i_k_point, i_k_task
        !integer :: i_cell, i_place, i_basis_1, i_basis_2
        !integer :: i_index 
        
        !integer, dimension(:,:),allocatable :: basis_index
        !real*8, dimension(:,:,:), allocatable :: work
        
        !integer :: i_index_real

        !allocate(basis_index(n_basis,n_basis))
        !allocate(work(n_cells_in_hamiltonian,n_basis,n_basis)) 

        !basis_index = 0
        !i_index_real = 0
        !do i_basis_2 = 1, n_basis, 1
          !do i_basis_1 = 1, i_basis_2, 1
            !i_index_real = i_index_real + 1
            !basis_index(i_basis_1,i_basis_2) = i_index_real
          !enddo
        !enddo

        !!finite difference calculation of H and S 
        !work = 0.d0
        !matrix_real= 0.0d0
        !if (packed_matrix_format .eq. PM_none) then

          !do i_basis_2 = 1, n_basis
            !do i_basis_1 = 1, i_basis_2
              !work(1,i_basis_1,i_basis_2) = matrix_sparse(basis_index(i_basis_1,i_basis_2))
              !!matrix_real(1,i_basis_1,i_basis_2) = matrix_sparse(basis_index(i_basis_1,i_basis_2))
            !enddo
          !enddo 
          !!matrix_real(1,:) = matrix_sparse(:)
        !else

          !do i_cell = 1,n_cells_in_hamiltonian-1
            !do i_basis_1 = 1, n_basis
                !if( index_hamiltonian(1,i_cell, i_basis_1) > 0 )then
                    !do i_place = index_hamiltonian(1,i_cell, i_basis_1), & 
                                !index_hamiltonian(2,i_cell, i_basis_1)
                        !i_basis_2 =  column_index_hamiltonian(i_place)
        
                        !work(i_cell,i_basis_1,i_basis_2) =  &
                           !work(i_cell,i_basis_1,i_basis_2) + &
                           !matrix_sparse(i_place)
                        !!matrix_real(i_cell,i_basis_1,i_basis_2) =  &
                          !!matrix_real(i_cell,i_basis_1,i_basis_2) + &
                          !!matrix_sparse(i_place)

                        !if (i_basis_1.ne.i_basis_2) then
                          !work(i_cell,i_basis_2,i_basis_1) =  &
                             !work(i_cell,i_basis_2,i_basis_1) + &
                             !matrix_sparse(i_place)
                          !!matrix_real(i_cell,i_basis_2,i_basis_1) =  &
                             !!matrix_real(i_cell,i_basis_2,i_basis_1) + &
                             !!matrix_sparse(i_place)
                        !endif
                    !enddo !i_place
                !endif !index_hamiltonian
            !enddo ! i_basis_1
          !enddo ! i_cell
        !endif
       
        !if (n_periodic.eq.0 .or. real_eigenvectors) then
          !do i_basis_2 = 1, n_basis
            !do i_basis_1 = 1, i_basis_2
              !matrix_real(:,basis_index(i_basis_1,i_basis_2)) = work(:,i_basis_1,i_basis_2)
            !enddo
          !enddo
        !else
          !i_k_task = 0
          !do i_k_point = 1, n_k_points, 1
            !if (myid.eq. MOD(i_k_point,n_tasks) .and. myid <= n_k_points) then
              !i_k_task = i_k_task +1
              !do i_basis_2 = 1, n_basis
                !do i_basis_1 = 1, i_basis_2
                  !do i_cell = 1, n_cells_in_hamiltonian-1
                    !matrix_complex(i_k_task,basis_index(i_basis_1,i_basis_2)) = &
                        !matrix_complex(i_k_task,basis_index(i_basis_1, i_basis_2)) + &
                        !work(i_k_task,i_basis_1,i_basis_2)*k_phase(i_cell,i_k_point)
                  !enddo
                !enddo
              !enddo
            !endif
          !enddo 
        !endif 

        !deallocate(work)
        !deallocate(basis_index)

    !end subroutine friction_construct_complex_matrix2
    
    subroutine friction_construct_complex_matrix(matrix_sparse, matrix_real, matrix_complex)
    !  NAME
    !    friction_construct_complex_matrix
    !  SYNOPSIS
    !  constructs the complex k-dependent versions of overlap or hamiltonian 
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2016)
        
        implicit none

        real*8, intent(in) :: matrix_sparse(:)
        real*8, intent(OUT) :: matrix_real(:,:)
        complex*16, intent(OUT) :: matrix_complex(:,:)
        ! imported variables

        ! local variables
        integer :: i_k_point, i_k_task,i_index_real
        integer :: i_cell, i_place, i_basis_1, i_basis_2
        integer :: k_cell, k_atom, k_cell_new, k_center_new
        
        integer, dimension(:,:),allocatable :: basis_index
        
        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index_real = 0
        do i_basis_2 = 1, n_basis, 1
          do i_basis_1 = 1, i_basis_2, 1
            i_index_real = i_index_real + 1
            basis_index(i_basis_1,i_basis_2) = i_index_real
          enddo
        enddo
        
        matrix_real= 0.0d0
        matrix_complex= (0.0d0,0.0d0)
        if (packed_matrix_format .eq. PM_none) then
          do i_basis_2 = 1, n_basis
            do i_basis_1 = 1, i_basis_2
              matrix_real(1,basis_index(i_basis_1,i_basis_2)) = matrix_sparse(basis_index(i_basis_1,i_basis_2))
            enddo
          enddo
        else
          !finite difference calculation of H and S
          matrix_real= 0.0d0
          matrix_complex= (0.0d0,0.0d0)
          i_k_task = 0
          do i_k_point = 1,n_k_points, 1
                if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
                   i_k_task = i_k_task + 1
        
                   do i_cell = 1,n_cells_in_hamiltonian-1
                     do i_basis_2 = 1, n_basis
                       if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                         i_index_real = index_hamiltonian(1,i_cell,i_basis_2)-1
                         do i_place = index_hamiltonian(1,i_cell, i_basis_2), &
                                      index_hamiltonian(2,i_cell, i_basis_2)
                           i_index_real = i_index_real +1
                           i_basis_1 =  column_index_hamiltonian(i_index_real)

                           if (real_eigenvectors) then
                               matrix_real(i_k_task,basis_index(i_basis_1, i_basis_2)) =  &
                                   matrix_real(i_k_task,basis_index(i_basis_1, i_basis_2)) + &
                                   dble(k_phase(i_cell,i_k_point))* &
                                   matrix_sparse(i_index_real)
                           else
                               matrix_complex(i_k_task,basis_index(i_basis_1, i_basis_2)) =  &
                                   matrix_complex(i_k_task,basis_index(i_basis_1, i_basis_2)) + &
                                   k_phase(i_cell,i_k_point)* &
                                   matrix_sparse(i_index_real)

                           endif
                         enddo !i_place
                       endif !index_hamiltonian
                     enddo ! i_basis_1
                   enddo ! i_cell
                    !build matrix_complex
                endif
          enddo
        endif

        if (allocated(basis_index)) deallocate(basis_index)

    end subroutine friction_construct_complex_matrix
    
    subroutine friction_construct_real_matrix(matrix_sparse, matrix_real)
    !  NAME
    !    friction_construct_real_matrix
    !  SYNOPSIS
    !  constructs the dense real space versions of overlap or hamiltonian 
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)
        
        implicit none

        real*8, intent(in) :: matrix_sparse(:)
        !real*8, intent(OUT) :: matrix_real(n_basis,n_cells_in_hamiltonian,n_basis)
        real*8, intent(out) :: matrix_real(:,:)
        ! imported variables

        ! local variables
        integer :: i_cell, i_place, i_basis_1, i_basis_2
        integer :: i_index 
        
        integer, dimension(:,:),allocatable :: basis_index
        real*8, dimension(:,:,:), allocatable :: work
        
        integer :: i_index_real

        allocate(basis_index(n_basis,n_basis))
        allocate(work(n_cells_in_hamiltonian,n_basis,n_basis)) 

        basis_index = 0
        i_index_real = 0
        do i_basis_2 = 1, n_basis, 1
          do i_basis_1 = 1, i_basis_2, 1
            i_index_real = i_index_real + 1
            basis_index(i_basis_1,i_basis_2) = i_index_real
          enddo
        enddo
        
        !finite difference calculation of H and S 
        work = 0.d0
        matrix_real= 0.0d0
        if (packed_matrix_format .eq. PM_none) then
          !do i_basis_2 = 1, n_basis
            !do i_basis_1 = 1, i_basis_2
              !matrix_real(1,basis_index(i_basis_1,i_basis_2)) = matrix_sparse(basis_index(i_basis_1,i_basis_2))
            !enddo
          !enddo
         matrix_real(1,:) = matrix_sparse(:) 
        else

          do i_cell = 1,n_cells_in_hamiltonian-1
            do i_basis_1 = 1, n_basis
                if( index_hamiltonian(1,i_cell, i_basis_1) > 0 )then
                    do i_place = index_hamiltonian(1,i_cell, i_basis_1), & 
                                index_hamiltonian(2,i_cell, i_basis_1)
                        i_basis_2 =  column_index_hamiltonian(i_place)
        
                        work(i_cell,i_basis_1,i_basis_2) =  &
                           work(i_cell,i_basis_1,i_basis_2) + &
                           matrix_sparse(i_place)
                        !matrix_real(i_cell,i_basis_1,i_basis_2) =  &
                          !matrix_real(i_cell,i_basis_1,i_basis_2) + &
                          !matrix_sparse(i_place)

                        if (i_basis_1.ne.i_basis_2) then
                          work(i_cell,i_basis_2,i_basis_1) =  &
                             work(i_cell,i_basis_2,i_basis_1) + &
                             matrix_sparse(i_place)
                          !matrix_real(i_cell,i_basis_2,i_basis_1) =  &
                             !matrix_real(i_cell,i_basis_2,i_basis_1) + &
                             !matrix_sparse(i_place)
                        endif
                    enddo !i_place
                endif !index_hamiltonian
            enddo ! i_basis_1
          enddo ! i_cell

          do i_basis_2 = 1, n_basis
            do i_basis_1 = 1, i_basis_2
              matrix_real(:,basis_index(i_basis_1,i_basis_2)) = work(:,i_basis_1,i_basis_2)
            enddo
          enddo
        endif

        deallocate(work)
        deallocate(basis_index)

    end subroutine friction_construct_real_matrix
    
    !subroutine friction_construct_real_matrix_scalapack(matrix_real)
    !!  NAME
    !!    friction_construct_real_matrix_scalapack
    !!  SYNOPSIS
    !!  constructs the dense real space versions of overlap or hamiltonian using scalapack 
    !!  AUTHOR
    !!   Reinhard J. Maurer, Yale University (2015)
       
        !implicit none

        !real*8, intent(in) :: matrix_sparse(n_hamiltonian_matrix_size)
        !real*8, intent(out) :: matrix_real(:,:,:)
        !! imported variables

        !! local variables
        !integer :: i_cell, i_place, i_basis_1, i_basis_2
        !integer :: i_index 
        
        !real*8, dimension(:,:,:), allocatable :: work

        !integer, dimension(:,:),allocatable :: basis_index
        
        !allocate(basis_index(n_basis,n_basis))

        !basis_index = 0
        !i_index = 0
        !do i_basis_2 = 1, n_basis, 1
          !do i_basis_1 = 1, i_basis_2
            !i_index = i_index + 1
            !basis_index(i_basis_1,i_basis_2) = i_index
          !enddo
        !enddo

        !!!finite difference calculation of H and S 
        !!matrix_real= 0.0d0
        !!work = 0.0d0
        !!if (packed_matrix_format .eq. PM_none) then

          !!!do i_basis_1 = 1, n_basis
            !!!do i_basis_2 = i_basis_1, n_basis
              !!!!work(i_basis_1,1,i_basis_2) = matrix_sparse(basis_index(i_basis_1,i_basis_2))
              !!!!matrix_real(1,basis_index(i_basis_1,i_basis_2)) = matrix_sparse(i_basis_1 +((i_basis_2-1)*i_basis_2)/2) 
              !!!!work(i_basis_2,1,i_basis_1) = work(i_basis_1,1,i_basis_2)
            !!!enddo
          !!!enddo 
          !!matrix_real(1,:) = matrix_sparse(:)

        !!else

          !!do i_cell = 1,n_cells_in_hamiltonian-1
            !!do i_basis_1 = 1, n_basis
                !!if( index_hamiltonian(1,i_cell, i_basis_1) > 0 )then
                    !!do i_place = index_hamiltonian(1,i_cell, i_basis_1), & 
                                !!index_hamiltonian(2,i_cell, i_basis_1)
                        !!i_basis_2 =  column_index_hamiltonian(i_place)
        
                        !!work(i_basis_1,i_cell,i_basis_2) =  &
                           !!work(i_basis_1,i_cell,i_basis_2) + &
                           !!matrix_sparse(i_place)
                       !!!matrix_real(i_cell,basis_index(i_basis_1,i_basis_2)) =  &
                           !!!matrix_real(i_cell,basis_index(i_basis_1,i_basis_2)) + &
                           !!!matrix_sparse(i_place)

                    !!if (i_basis_1.ne.i_basis_2) then
                        !!work(i_basis_2,i_cell,i_basis_1) =  &
                           !!work(i_basis_2,i_cell,i_basis_1) + &
                           !!matrix_sparse(i_place)
                    !!endif
                    !!enddo !i_place
                !!endif !index_hamiltonian
            !!enddo ! i_basis_1
          !!enddo ! i_cell

          !!do i_basis_2 = 1, n_basis
            !!do i_basis_1 = 1, i_basis_2
              !!matrix_real(:,basis_index(i_basis_1,i_basis_2)) = work(i_basis_2,:,i_basis_1)
            !!enddo
          !!enddo
        !!endif

        !deallocate(work)
        !deallocate(basis_index)

    !end subroutine friction_construct_real_matrix_scalapack
    !----------------------------------------------!
    !!!!!Friction UTILITIES
    !----------------------------------------------!

    subroutine friction_output_matrices_H_S &
            ( first_order_S, first_order_H )

        use dimensions
        use basis
        use localorb_io
        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(in) :: first_order_S(&
            3,friction_n_active_atoms,1,(n_basis*(n_basis+1))/2,1)
        real*8, intent(in) :: first_order_H(&
            3,friction_n_active_atoms,1,(n_basis*(n_basis+1))/2,1,n_spin)

        !  INPUTS
        !   o first_order_S -- overlap matrix
        !   o first_order_H -- Hamiltonian matrix
        !
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

        real*8 output_element(n_basis)

        character(len=30),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        
        integer, dimension(:,:),allocatable :: basis_index
        integer :: i_index

        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis, 1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        !  begin work

        write (info_str,'(2X,A)') &
            "Writing basis function properties and H^(1) and S^(1) matrices."
        call localorb_info(info_str,use_unit,'(A)')

        write (info_str,'(2X,I5,A)') n_basis, " basis functions used."
        call localorb_info(info_str,use_unit,'(A)')

        !  write basis function properties

        open (50, file="basis-indices.out")

        write(info_str,*)
        call localorb_info(info_str,50,'(A)')

        write(info_str,'(A5,1X,A8,1X,A3,1X,A3,1X,A3,1X,A3)') &
            "fn.", "  type  ", "at.", "n", &
            "l", "m"
        call localorb_info(info_str,50,'(A)')

        do i_basis = 1, n_basis, 1
            i_fn = basis_fn(i_basis)

            write(info_str,'(I5,1X,A8,1X,I3,1X,I3,1X,I3,1X,I3)') &
                i_basis, basisfn_type(i_fn), &
                basis_atom(i_basis), basisfn_n(i_fn), basis_l(i_basis), &
                basis_m(i_basis)
            call localorb_info(info_str,50,'(A)')
        enddo
        close (50)

        do i_atom = 1, friction_n_active_atoms, 1
        do i_cart = 1, 3, 1
            !  now, write the overlap matrix
            if (i_cart==1) then 
                write(file_name(1),'(A,I0.3,A)') "first_order_S_",friction_index_list(i_atom),'_x.out'
            else if (i_cart==2) then
                write(file_name(1),'(A,I0.3,A)') "first_order_S_",friction_index_list(i_atom),'_y.out'
            else 
                write(file_name(1),'(A,I0.3,A)') "first_order_S_",friction_index_list(i_atom),'_z.out'
            endif
            open (50, file=file_name(1), form='unformatted')
            do j_basis = 1, n_basis, 1
                write(50) (first_order_S(i_cart,i_atom,1,basis_index(i_basis,j_basis),1),i_basis=1,j_basis )
            enddo
            close(50)
            
            !  now, write the Hamiltonian matrix
            if (n_spin.eq.1) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_",friction_index_list(i_atom),'_z.out'
                endif
            else if (n_spin.eq.2) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_up_",friction_index_list(i_atom),'_x.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_H_dn_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_up_",friction_index_list(i_atom),'_y.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_H_dn_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_up_",friction_index_list(i_atom),'_z.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_H_dn_",friction_index_list(i_atom),'_z.out'
                endif
            end if

            do i_spin = 1, n_spin, 1
                open (50, file=file_name(i_spin), form='unformatted')
                do j_basis = 1, n_basis, 1
                    write(50) (first_order_H(i_cart,i_atom,1,basis_index(i_basis,j_basis),1,i_spin),&
                        i_basis=1,j_basis )
                enddo
                close(50)
            ! spin
            enddo
        enddo
        enddo

        deallocate(basis_index)

    end subroutine friction_output_matrices_H_S
    
    subroutine friction_output_matrix_G &
            ( first_order_G )

        use dimensions
        use basis
        use localorb_io
        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(in) :: first_order_G(&
            3,friction_n_active_atoms,1,(n_basis*(n_basis+1))/2,1,n_spin)

        !  INPUTS
        !   o first_order_G -- coupling matrix
        !
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

        real*8 output_element(n_basis)

        character(len=30),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        
        integer, dimension(:,:),allocatable :: basis_index
        integer :: i_index

        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis, 1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        !  begin work

        write (info_str,'(2X,A)') &
            "Writing basis function properties and G matrix."
        call localorb_info(info_str,use_unit,'(A)')

        write (info_str,'(2X,I5,A)') n_basis, " basis functions used."
        call localorb_info(info_str,use_unit,'(A)')

        !  write basis function properties

        open (50, file="basis-indices.out")

        write(info_str,*)
        call localorb_info(info_str,50,'(A)')

        write(info_str,'(A5,1X,A8,1X,A3,1X,A3,1X,A3,1X,A3)') &
            "fn.", "  type  ", "at.", "n", &
            "l", "m"
        call localorb_info(info_str,50,'(A)')

        do i_basis = 1, n_basis, 1
            i_fn = basis_fn(i_basis)

            write(info_str,'(I5,1X,A8,1X,I3,1X,I3,1X,I3,1X,I3)') &
                i_basis, basisfn_type(i_fn), &
                basis_atom(i_basis), basisfn_n(i_fn), basis_l(i_basis), &
                basis_m(i_basis)
            call localorb_info(info_str,50,'(A)')
        enddo
        close (50)

        do i_atom = 1, friction_n_active_atoms, 1
        do i_cart = 1, 3, 1
            
            !  now, write the Hamiltonian matrix
            if (n_spin.eq.1) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_",friction_index_list(i_atom),'_z.out'
                endif
            else if (n_spin.eq.2) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_up_",friction_index_list(i_atom),'_x.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_G_dn_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_up_",friction_index_list(i_atom),'_y.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_G_dn_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_up_",friction_index_list(i_atom),'_z.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_G_dn_",friction_index_list(i_atom),'_z.out'
                endif
            end if

            do i_spin = 1, n_spin, 1
                open (50, file=file_name(i_spin), form='unformatted')
                do j_basis = 1, n_basis, 1
                    write(50) (first_order_G(i_cart,i_atom,1,basis_index(i_basis,j_basis),1,i_spin),&
                        i_basis=1,j_basis )
                enddo
                close(50)
            ! spin
            enddo
        enddo
        enddo

        deallocate(basis_index)

    end subroutine friction_output_matrix_G
    
    subroutine friction_output_matrices_H_S_p1 &
            ( first_order_S, first_order_H, first_order_S_cmplx, first_order_H_cmplx)

        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(in) :: first_order_S(:,:,:,:,:)
        real*8, intent(in) :: first_order_H(:,:,:,:,:,:)
        complex*16, intent(in) :: first_order_S_cmplx(:,:,:,:,:)
        complex*16, intent(in) :: first_order_H_cmplx(:,:,:,:,:,:)

        !  INPUTS
        !   o first_order_S -- overlap matrix
        !   o first_order_H -- Hamiltonian matrix
        !   o first_order_S_cmplx -- overlap matrix
        !   o first_order_H_cmplx -- Hamiltonian matrix
        !
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

        real*8 output_element(n_basis)

        character(len=50),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer :: i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        integer :: i_cell, i_index, i_k_task,i_k_point 
        
        integer, dimension(:,:),allocatable :: basis_index
        
        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis, 1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        !  begin work

        write (info_str,'(2X,A)') &
            "Writing basis function properties and H^(1) and S^(1) matrices."
        call localorb_info(info_str,use_unit,'(A)')

        write (info_str,'(2X,I5,A)') n_basis, " basis functions used."
        call localorb_info(info_str,use_unit,'(A)')

        !  write basis function properties

        open (50, file="basis-indices.out")
        write(info_str,*)
        call localorb_info(info_str,50,'(A)')

        write(info_str,'(A5,1X,A8,1X,A3,1X,A3,1X,A3,1X,A3)') &
            "fn.", "  type  ", "at.", "n", &
            "l", "m"
        call localorb_info(info_str,50,'(A)')
        do i_basis = 1, n_basis, 1
            i_fn = basis_fn(i_basis)
            write(info_str,'(I5,1X,A8,1X,I3,1X,I3,1X,I3,1X,I3)') &
                i_basis, basisfn_type(i_fn), &
                basis_atom(i_basis), basisfn_n(i_fn), basis_l(i_basis), &
                basis_m(i_basis)
            call localorb_info(info_str,50,'(A)')
        enddo
        close (50)
        
        if(friction_use_complex_matrices) then
            
            i_k_task = 0
            do i_k_point = 1, n_k_points, 1
            if (myid.eq.MOD(i_k_point,n_tasks) .and. myid <= n_k_points) then
                i_k_task = i_k_task + 1
                write (info_str,'(2X,A,I8)') &
                    "Writing H^(1) and S^(1) matrices for k point ", i_k_point
                call localorb_info(info_str,use_unit,'(A)')
                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the overlap matrix
                    if (i_cart==1) then 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_x.out'
                    else if (i_cart==2) then
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_y.out'
                    else 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_z.out'
                    endif
                    open (50, file=file_name(1),form='unformatted')
                    if (real_eigenvectors) then
                        do j_basis = 1, n_basis, 1
                            write(50) (first_order_S(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1),&
                                i_basis=1,j_basis )
                        enddo
                    else
                        do j_basis = 1, n_basis, 1
                            write(50) (first_order_S_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1),&
                                i_basis=1,j_basis )
                        enddo
                    endif
                    close(50)
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if
                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),form='unformatted')
                        if (real_eigenvectors) then
                            do j_basis = 1, n_basis, 1
                                write(50) (first_order_H(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1,j_basis )
                            enddo
                        else
                            do j_basis = 1, n_basis, 1
                                write(50) (first_order_H_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1,j_basis )
                            enddo
                        endif
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            endif
            enddo
        else
            do i_cell=1, n_cells_in_hamiltonian
              if (myid.eq.MOD(i_cell,n_tasks)) then
              !if (myid.eq. 0) then
                do i_atom = 1, friction_n_active_atoms, 1
                  do i_cart = 1, 3, 1
                    !  now, write the overlap matrix
                    if (i_cart==1) then 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_x.out'
                    else if (i_cart==2) then
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_y.out'
                    else 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_z.out'
                    endif
                    open (50, file=file_name(1),form='unformatted')
                    do j_basis = 1, n_basis, 1
                        write(50) (first_order_S(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1),&
                            i_basis=1, j_basis)
                    enddo
                    close(50)
                    
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if

                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),form='unformatted')
                        do j_basis = 1, n_basis, 1
                            write(50) (first_order_H(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1,i_spin),& 
                                i_basis=1, j_basis )
                        enddo
                        close(50)
                    ! spin
                    enddo
                  enddo
                enddo
              endif
            enddo
        endif

        deallocate(basis_index)

    end subroutine friction_output_matrices_H_S_p1
    
    subroutine friction_output_matrix_G_p1 &
            ( first_order_G,first_order_G_cmplx)

        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(in) :: first_order_G(:,:,:,:,:,:)
        complex*16, intent(in) :: first_order_G_cmplx(:,:,:,:,:,:)

        !  INPUTS
        !   o first_order_G -- Hamiltonian matrix
        !   o first_order_G_cmplx -- coupling matrix
        !
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

        real*8 output_element(n_basis)

        character(len=50),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer :: i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        integer :: i_cell, i_index, i_k_task,i_k_point 
        
        integer, dimension(:,:),allocatable :: basis_index
        
        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis, 1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        !  begin work

        write (info_str,'(2X,A)') &
            "Writing basis function properties and H^(1) and S^(1) matrices."
        call localorb_info(info_str,use_unit,'(A)')

        write (info_str,'(2X,I5,A)') n_basis, " basis functions used."
        call localorb_info(info_str,use_unit,'(A)')

        !  write basis function properties

        open (50, file="basis-indices.out")
        write(info_str,*)
        call localorb_info(info_str,50,'(A)')

        write(info_str,'(A5,1X,A8,1X,A3,1X,A3,1X,A3,1X,A3)') &
            "fn.", "  type  ", "at.", "n", &
            "l", "m"
        call localorb_info(info_str,50,'(A)')
        do i_basis = 1, n_basis, 1
            i_fn = basis_fn(i_basis)
            write(info_str,'(I5,1X,A8,1X,I3,1X,I3,1X,I3,1X,I3)') &
                i_basis, basisfn_type(i_fn), &
                basis_atom(i_basis), basisfn_n(i_fn), basis_l(i_basis), &
                basis_m(i_basis)
            call localorb_info(info_str,50,'(A)')
        enddo
        close (50)
        
        if(friction_use_complex_matrices) then
            
            i_k_task = 0
            do i_k_point = 1, n_k_points, 1
            if (myid.eq.MOD(i_k_point,n_tasks) .and. myid <= n_k_points) then
                i_k_task = i_k_task + 1
                write (info_str,'(2X,A,I8)') &
                    "Writing G matrix for k point ", i_k_point
                call localorb_info(info_str,use_unit,'(A)')
                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if
                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),form='unformatted')
                        if (real_eigenvectors) then
                            do j_basis = 1, n_basis, 1
                                write(50) (first_order_G(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1,j_basis )
                            enddo
                        else
                            do j_basis = 1, n_basis, 1
                                write(50) (first_order_G_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1,j_basis )
                            enddo
                        endif
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            endif
            enddo
        else
            do i_cell=1, n_cells_in_hamiltonian
              if (myid.eq.MOD(i_cell,n_tasks)) then
                do i_atom = 1, friction_n_active_atoms, 1
                  do i_cart = 1, 3, 1
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if

                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),form='unformatted')
                        do j_basis = 1, n_basis, 1
                            write(50) (first_order_G(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1,i_spin),& 
                                i_basis=1, j_basis )
                        enddo
                        close(50)
                    ! spin
                    enddo
                  enddo
                enddo
              endif
            enddo
        endif

        deallocate(basis_index)

    end subroutine friction_output_matrix_G_p1
    
    subroutine friction_read_matrices_H_S &
            ( first_order_S, first_order_H)

        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(out) :: first_order_S(&
            3,friction_n_active_atoms,1,(n_basis*(n_basis+1))/2,1)
        real*8, intent(out) :: first_order_H(&
            3,friction_n_active_atoms,1,(n_basis*(n_basis+1))/2,1,n_spin)

        !  INPUTS
        !   o first_order_S -- overlap matrix
        !   o first_order_H -- Hamiltonian matrix
        !
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

        real*8 output_element(n_basis)

        character(len=30),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        
        integer, dimension(:,:),allocatable :: basis_index
        integer :: i_index

        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis, 1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        !  begin work

        write (info_str,'(2X,A)') &
            "Reading H^(1) and S^(1) matrices."
        call localorb_info(info_str,use_unit,'(A)')

        do i_atom = 1, friction_n_active_atoms, 1
        do i_cart = 1, 3, 1
            !  now, write the overlap matrix
            if (i_cart==1) then 
                write(file_name(1),'(A,I0.3,A)') "first_order_S_",friction_index_list(i_atom),'_x.out'
            else if (i_cart==2) then
                write(file_name(1),'(A,I0.3,A)') "first_order_S_",friction_index_list(i_atom),'_y.out'
            else 
                write(file_name(1),'(A,I0.3,A)') "first_order_S_",friction_index_list(i_atom),'_z.out'
            endif
            open (50, file=file_name(1), form='unformatted', status='old')
            do j_basis = 1, n_basis, 1
                read(50) (first_order_S(i_cart,i_atom,1,basis_index(i_basis,j_basis),1),i_basis=1,j_basis )
            enddo
            close(50)
            
            !  now, write the Hamiltonian matrix
            if (n_spin.eq.1) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_",friction_index_list(i_atom),'_z.out'
                endif
            else if (n_spin.eq.2) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_up_",friction_index_list(i_atom),'_x.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_H_dn_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_up_",friction_index_list(i_atom),'_y.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_H_dn_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_H_up_",friction_index_list(i_atom),'_z.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_H_dn_",friction_index_list(i_atom),'_z.out'
                endif
            end if

            do i_spin = 1, n_spin, 1
                open (50, file=file_name(i_spin), form='unformatted', status='old')
                do j_basis = 1, n_basis, 1
                    read(50) (first_order_H(i_cart,i_atom,1,basis_index(i_basis,j_basis),1,i_spin),&
                        i_basis=1,j_basis )
                enddo
                close(50)
            ! spin
            enddo
        enddo
        enddo

        deallocate(basis_index)

    end subroutine friction_read_matrices_H_S
    
    subroutine friction_read_matrix_G &
            ( first_order_G)

        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(out) :: first_order_G(&
            3,friction_n_active_atoms,1,(n_basis*(n_basis+1))/2,1,n_spin)

        !  INPUTS
        !   o first_order_G -- coupling matrix
        !
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

        real*8 output_element(n_basis)

        character(len=30),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        
        integer, dimension(:,:),allocatable :: basis_index
        integer :: i_index

        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis, 1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        !  begin work

        write (info_str,'(2X,A)') &
            "Reading G matrix."
        call localorb_info(info_str,use_unit,'(A)')

        do i_atom = 1, friction_n_active_atoms, 1
        do i_cart = 1, 3, 1
            
            !  now, write the Hamiltonian matrix
            if (n_spin.eq.1) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_",friction_index_list(i_atom),'_z.out'
                endif
            else if (n_spin.eq.2) then
                if (i_cart==1) then 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_up_",friction_index_list(i_atom),'_x.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_G_dn_",friction_index_list(i_atom),'_x.out'
                else if (i_cart==2) then
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_up_",friction_index_list(i_atom),'_y.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_G_dn_",friction_index_list(i_atom),'_y.out'
                else 
                    write(file_name(1),'(A,I0.3,A)') "first_order_G_up_",friction_index_list(i_atom),'_z.out'
                    write(file_name(2),'(A,I0.3,A)') "first_order_G_dn_",friction_index_list(i_atom),'_z.out'
                endif
            end if

            do i_spin = 1, n_spin, 1
                open (50, file=file_name(i_spin), form='unformatted', status='old')
                do j_basis = 1, n_basis, 1
                    read(50) (first_order_G(i_cart,i_atom,1,basis_index(i_basis,j_basis),1,i_spin),&
                        i_basis=1,j_basis )
                enddo
                close(50)
            ! spin
            enddo
        enddo
        enddo

        deallocate(basis_index)

    end subroutine friction_read_matrix_G
    
    subroutine friction_read_matrix_G_p1 &
            ( first_order_G, first_order_G_cmplx)

        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(out) :: first_order_G(:,:,:,:,:,:)
        complex*16, intent(out) :: first_order_G_cmplx(:,:,:,:,:,:)

        !  INPUTS
        !   o first_order_G -- Coupling matrix
        !   o first_order_G_cmplx -- Coupling matrix
        !
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

        real*8 output_element(n_basis)

        character(len=50),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer :: i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        integer :: i_cell, i_index, i_k_task,i_k_point 
        
        integer, dimension(:,:),allocatable :: basis_index
        
        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis,1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        if (friction_use_complex_matrices) then
            if (myid.eq.0) then
                write (info_str,'(2X,A)') &
                    "Reading complex G matrix."
                call localorb_info(info_str,use_unit,'(A)')
            endif
            i_k_task = 0
            do i_k_point = 1, n_k_points, 1
            if (myid.eq.MOD(i_k_point,n_tasks) .and. myid <= n_k_points) then
                i_k_task = i_k_task + 1

                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if
                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),status='old',form='unformatted')
                        if (real_eigenvectors) then
                            do j_basis = 1, n_basis, 1
                                read(50) (first_order_G(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1, j_basis )
                            enddo
                        else
                            do j_basis = 1, n_basis, 1
                                read(50) (first_order_G_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1, j_basis )
                            enddo
                        endif
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            endif
            enddo
        else
            if (myid.eq.0) then
                write (info_str,'(2X,A)') &
                    "Reading real space G matrix."
                call localorb_info(info_str,use_unit,'(A)')
            endif
            do i_cell=1, n_cells_in_hamiltonian
                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_up_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_G_dn_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if

                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),status='old',form='unformatted')
                        do j_basis = 1, n_basis, 1
                            read(50) (first_order_G(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1,i_spin),&
                                i_basis=1, j_basis )
                        enddo
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            enddo
        endif

        deallocate(basis_index)

    end subroutine friction_read_matrix_G_p1

    subroutine friction_read_matrices_H_S_p1 &
            ( first_order_S, first_order_H, first_order_S_cmplx, first_order_H_cmplx)

        implicit none

        !  ARGUMENTS
        !  imported variables

        real*8, intent(out) :: first_order_S(:,:,:,:,:)
        real*8, intent(out) :: first_order_H(:,:,:,:,:,:)
        complex*16, intent(out) :: first_order_S_cmplx(:,:,:,:,:)
        complex*16, intent(out) :: first_order_H_cmplx(:,:,:,:,:,:)

        !  INPUTS
        !   o first_order_S -- overlap matrix
        !   o first_order_H -- Hamiltonian matrix
        !   o first_order_S_cmplx -- overlap matrix
        !   o first_order_H_cmplx -- Hamiltonian matrix
        !
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

        real*8 output_element(n_basis)

        character(len=50),dimension(n_spin) :: file_name

        character(len=100) :: info_str

        !  counters

        integer :: i_cart, i_atom, i_basis, j_basis, i_spin,i_fn
        integer :: i_cell, i_index, i_k_task,i_k_point 
        
        integer, dimension(:,:),allocatable :: basis_index
        
        allocate(basis_index(n_basis,n_basis))
        
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis,1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo
        !  write basis function properties
        if (friction_use_complex_matrices) then
            if (myid.eq.0) then
                write (info_str,'(2X,A)') &
                    "Reading complex H^(1) and S^(1) matrices."
                call localorb_info(info_str,use_unit,'(A)')
            endif
            i_k_task = 0
            do i_k_point = 1, n_k_points, 1
            if (myid.eq.MOD(i_k_point,n_tasks) .and. myid <= n_k_points) then
                i_k_task = i_k_task + 1
                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the overlap matrix
                    if (i_cart==1) then 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_x.out'
                    else if (i_cart==2) then
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_y.out'
                    else 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_z.out'
                    endif
                    open (50, file=file_name(1),status='old',form='unformatted')
                    if (real_eigenvectors) then
                        do j_basis = 1, n_basis, 1
                            read(50) (first_order_S(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1),&
                                i_basis=1, j_basis )
                        enddo
                    else
                        do j_basis = 1, n_basis, 1
                            read(50) (first_order_S_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1),&
                                i_basis=1, j_basis )
                        enddo
                    endif
                    close(50)
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if
                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),status='old',form='unformatted')
                        if (real_eigenvectors) then
                            do j_basis = 1, n_basis, 1
                                read(50) (first_order_H(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1, j_basis )
                            enddo
                        else
                            do j_basis = 1, n_basis, 1
                                read(50) (first_order_H_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1, j_basis )
                            enddo
                        endif
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            endif
            enddo
        else
            if (myid.eq.0) then
                write (info_str,'(2X,A)') &
                    "Reading real space H^(1) and S^(1) matrices."
                call localorb_info(info_str,use_unit,'(A)')
            endif
            do i_cell=1, n_cells_in_hamiltonian
                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the overlap matrix
                    if (i_cart==1) then 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_x.out'
                    else if (i_cart==2) then
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_y.out'
                    else 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_z.out'
                    endif
                    open (50, file=file_name(1),status='old',form='unformatted')
                    do j_basis = 1, n_basis, 1
                        read(50) (first_order_S(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1),&
                            i_basis=1, j_basis )
                    enddo
                    close(50)
                    
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if

                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),status='old',form='unformatted')
                        do j_basis = 1, n_basis, 1
                            read(50) (first_order_H(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1,i_spin),&
                                i_basis=1, j_basis )
                        enddo
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            enddo
        
        endif

        deallocate(basis_index)

    end subroutine friction_read_matrices_H_S_p1
        
    subroutine friction_print_matrix(matrix)

        implicit none

        !  ARGUMENTS
        !  imported variables

        !real*8, intent(in) :: matrix(:,:)
        real*8, intent(in) :: matrix(:)
        !  local variables
        character(len=400) :: info_str
        !  counters

        integer :: i_cart, i_atom, i_basis, j_basis
        integer :: p, pstart, pend, i_index, j_coord,n_print
        
        integer, dimension(:,:),allocatable :: basis_index
        real*8, dimension(:,:),allocatable :: matrix2
        allocate(basis_index(n_basis,n_basis))
        allocate(matrix2(n_basis,n_basis))
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis,1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
            matrix2(i_basis,j_basis) = matrix(basis_index(i_basis,j_basis))
            matrix2(j_basis,i_basis) = matrix(basis_index(i_basis,j_basis))
            !matrix2(i_basis,j_basis) = matrix(i_basis,j_basis)
            !matrix2(j_basis,i_basis) = matrix(j_basis,i_basis)
          enddo
        enddo

        if(myid.eq.0 ) then
            write (info_str,'(4X,A)') &
                "*************Printing S1, 11*************"
            call localorb_info(info_str,use_unit,'(A)')
            n_print = n_basis/6+mod(n_basis,6)
            do i_basis=1, n_basis
                write (info_str,'(4X,I6)') i_basis
                do p=1, n_print
                  pstart = 1+(p-1)*6 
                  pend = pstart+5
                  if (pend>n_basis) pend=n_basis
                  if (p==1) then
                      do j_basis=pstart, pend
                          write (info_str,'(4X,I6,6E18.9)') i_basis, &
                              (matrix2(i_basis,j_coord), j_coord=pstart, pend)
                          call localorb_info(info_str,use_unit,'(A)')
                      enddo
                  else
                      write (info_str,'(10X,6E18.9)') &
                          (matrix2(i_basis,j_coord), j_coord=pstart, pend)
                      call localorb_info(info_str,use_unit,'(A)')
                  endif
                enddo
            enddo
            write (info_str,'(4X,A)') &
                "**********END Printing S1************"
            call localorb_info(info_str,use_unit,'(A)')
        endif
   
        deallocate(matrix2)
        deallocate(basis_index)

    end subroutine friction_print_matrix

    subroutine friction_output_eigenvectors_jmol &
            (eigvals, eigvecs, filename_inp )

        implicit none

        !  ARGUMENTS
        !  imported variables

        !  INPUTS
        !   o first_order_S -- overlap matrix
        !   o first_order_H -- Hamiltonian matrix
        !
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

        real*8, dimension(:), intent(in) :: eigvals
        complex*16, dimension(:,:), intent(in) :: eigvecs
        character(len=100), intent(inout), optional :: filename_inp
        character(len=100) :: filename
        integer :: n, i, j
        real*8, dimension(:),allocatable :: vec
        integer :: v, atom_counter
        real*8 :: norm

        if (.not. present(filename_inp)) then
            filename = 'friction_eigenvectors.jmol'
        else
            filename = filename_inp
        endif
        
        open(78, file=filename)
        allocate(vec(3*friction_n_active_atoms))

        !loop over friction modes
        do n=1, 3*friction_n_active_atoms
            write(78,'(I6)') n_atoms
            write(78,'(A,I3,A,F16.5,A)') 'Mode #',n,' f = ',ps/eigvals(n),' ps'
            j = 0
            norm = 0.0
            do i=1, n_atoms
                if (friction_atoms_list(i)) then
                    do v=1, 3
                      vec(j*3+v) = real(eigvecs(j*3+v,n))/sqrt(masses(j+1))
                      norm = norm + vec(j*3+v)*vec(j*3+v)
                    enddo
                    j = j +1
                endif
            enddo
            vec(:) = vec(:) /sqrt(norm)
            j = 0
            do i=1, n_atoms
                if (friction_atoms_list(i)) then
                    write(78,'(A,2X,6F12.6)') species_name(species(i)), &
                       coords(1,i)*bohr, coords(2,i)*bohr, coords(3,i)*bohr, &
                       vec(j*3+1), vec(j*3+2), vec(j*3+3)
                       !real(friction_eigvecs(n,j*3+1)), real(friction_eigvecs(n,j*3+2)), &
                       !real(friction_eigvecs(n,j*3+3))
                       !friction_eigvecs(n,j*3+1), friction_eigvecs(n,j*3+2), friction_eigvecs(n,j*3+3)
                       !abs(friction_eigvecs(n,j*3+1))*&
                       !sign(atan2(aimag(friction_eigvecs(n,j*3+1)),real(friction_eigvecs(n,j*3+1)))),&
                       !abs(friction_eigvecs(n,j*3+2))*&
                       !sign(atan2(aimag(friction_eigvecs(n,j*3+2)),real(friction_eigvecs(n,j*3+2)))),&
                       !abs(friction_eigvecs(n,j*3+3))*&
                       !sign(atan2(aimag(friction_eigvecs(n,j*3+3)),real(friction_eigvecs(n,j*3+3)))),&
                    j = j +1
                else
                    write(78,'(A,2X,6F12.6)') species_name(species(i)), &
                       coords(1,i)*bohr, coords(2,i)*bohr, coords(3,i)*bohr, &
                       0.000000,  0.000000,  0.000000
                endif
            enddo
        enddo

        deallocate(vec)
        close(78)

    end subroutine friction_output_eigenvectors_jmol

    subroutine friction_output_tensor &
            (tensor, filename_inp )

        implicit none

        !  ARGUMENTS
        !  imported variables

        !  INPUTS
        !   o tensor -- matrix
        !
        !  OUTPUT
        !    none

        !  local variables

        complex*16, dimension(:,:), intent(in) :: tensor 
        character(len=100), intent(inout), optional :: filename_inp
        character(len=100) :: filename
        integer :: n, i_cart, i, j

        if (.not. present(filename_inp)) then
            filename = 'friction_tensor.out'
        else
            filename = filename_inp
        endif
        
        open(79, file=filename)

        !loop over friction modes
        i = 0
        write(79,'(A)') '# friction_tensor in 1/ps'
        do n=1, friction_n_active_atoms
          do i_cart=1, 3
            i = i + 1
            write(79,'(A,I3,A,I3,A,I3)') '# n_atom ',friction_index_list(n),' i_cart ', i_cart, ' index ',i
            
            write(79,'(200F12.6)') (dble(tensor(i,j))/ps, j=1, 3*friction_n_active_atoms)
          enddo
        enddo

        close(79)

    end subroutine friction_output_tensor

    function delta_function(x,x0,s)
        implicit none

        real*8, intent(IN) :: x
        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s
        
        real*8 :: delta_function

        select case (friction_delta_type)
            case('square')
                delta_function = delta_function_square(x,x0,s)
            case('gaussian')
                delta_function = delta_function_gaussian(x,x0,s)
            case('squashed_fermi')
                delta_function = delta_function_squashed_fermi(x,s)
            case('lorentzian')
                delta_function = delta_function_lorentzian(x,x0,s)
            case('sine')
                delta_function = delta_function_sine(x,x0,s)
            case default
                call aims_stop('friction: Unknown friction_delta_type, Aims Stops!')
        end select

    end function delta_function

    !function delta_function_methfessel(e,e0,s,N)
        !implicit none

        !real*8, intent(IN) :: e
        !real*8, intent(IN) :: e0
        !real*8, intent(IN) :: s
        !integer,intent(IN) :: N

        !real*8 :: delta_function_methfessel, tmp,x
        
        !integer :: i

        !x = (e-e0)/s 

        !delta_function_methfessel = 0.d0

        !do i=0, N
            !tmp = 0.d0
            !select case(i) 
                !case(0)
                    !tmp = exp(-x*x)/sqrt(pi)
                !case(1)
                    !tmp = -(4.d0*x*x-2)*exp(-x*x)/(4*sqrt(pi))
                !case(2)
                    !tmp = (16.d0*x**4-48.d0*x*x*12)*exp(-x*x)/(32.d0*sqrt(pi))
                !case(3)
                    !tmp = -(64.d0*x**6-480.d0*x**4+720.d0*x*x-120.d0)*exp(-x*x)/(384.d0*sqrt(pi))
                !case(4)
                    !tmp = (256.d0*x**8-3584.d0*x**6+13440.d0*x**4-13440.d0*x*x+1680.d0)*&
                        !exp(-x*x)/(6144.d0*sqrt(pi))
                !case(5)
                    !tmp = -(1024.d0*x**10-23040.d0*x**8+161280.d0*x**6-403200.d0*x**4+302400.d0*x*x-30240.d0)*&
                        !exp(-x*x)/(122880.d0*sqrt(pi))
                !!case(6)
                    !!tmp = (1024.d0*x**10-23040.d0*x**8+161280.d0*x**6-403200.d0*x**4+302400.d0*x*x-30240.d0)*&
                        !!exp(-x*x)/(122880.d0*sqrt(pi))

                !case default
                    !call aims_stop('delta_function_methfessel: N>6 not allowed')
            !end select
            !delta_function_methfessel = delta_function_methfessel + tmp
        !end do

    !end function delta_function_methfessel

    function delta_function_gaussian(x,x0,s)
        implicit none

        real*8, intent(IN) :: x
        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s
        
        real*8 :: delta_function_gaussian

        delta_function_gaussian = (exp(-0.5*((x-x0)*(x-x0))/(s*s))/(s*sqrt_pi))*one_over_sqrt2

    end function delta_function_gaussian
    
    function delta_function_square(x,x0,s)
        implicit none

        real*8, intent(IN) :: x
        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s
        
        real*8 :: delta_function_square

        if (abs(x-x0)>s) then
            delta_function_square = 0.0
        else
            delta_function_square = 1./s
        endif

    end function delta_function_square
    
    function delta_function_squashed_fermi(x,s)
        implicit none

        real*8, intent(IN) :: x
        real*8, intent(IN) :: s
        
        real*8 :: delta_function_squashed_fermi

        real*8 :: y

        y = (x/s)*sqrt(2/pi) +sqrt(0.5)
        delta_function_squashed_fermi = 2*(y/s)*sqrt(2/pi)*exp(0.5-y*y)

    end function

    function delta_function_lorentzian(x,x0,s)
        implicit none

        real*8, intent(IN) :: x
        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s
        
        real*8 :: delta_function_lorentzian

        delta_function_lorentzian = (1.d0/pi)*((0.5d0*s)/&
            ((x-x0)*(x-x0)+(0.5d0*s)*(0.5d0*s)))

    end function delta_function_lorentzian
    
    function delta_function_sine(x,x0,s)
        implicit none

        real*8, intent(IN) :: x
        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s
        
        real*8 :: delta_function_sine, t

        t = 1/s
        if ((x-x0)<=tiny(x)) then
            delta_function_sine = t/pi
        else
            delta_function_sine = (sin((x-x0)*(t))*sin((x-x0)*(t)))/&
            ((x-x0)*(x-x0)*t*t)
        endif
    end function delta_function_sine

    function fermi_pop(x,x0, T, n_spin)
        implicit none

        real*8, intent(IN) :: x
        real*8, intent(IN) :: x0
        real*8, intent(IN) :: T 
        integer,intent(IN) :: n_spin 

        real*8 :: fermi_pop

        if (T<tiny(1.0)) then
            if (x<x0) then
                !fermi_pop = 2.0/n_spin
                fermi_pop = 1.0
            else
                fermi_pop = 0.0
            end if
        else
            !fermi_pop = (2.0/n_spin)*(1./(exp((x-x0)/(boltzmann_kB*T))+1))
            fermi_pop = (1./(exp((x-x0)/(boltzmann_kB*T))+1))
        end if

    end function fermi_pop
    
    function gaussian_norm(x0, s)
        implicit none

        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s

        real*8 :: gaussian_norm 

        gaussian_norm = 0.5d0 * (1-erf((-x0/s)*(one_over_sqrt2)))

    end function gaussian_norm

    function gaussian_norm2(x0, s)
        implicit none

        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s

        real*8 :: gaussian_norm2

        gaussian_norm2 = 0.5d0 * (1-erf((x0/s)*(one_over_sqrt2)))

    end function gaussian_norm2

    function gaussian_norm3(x0, s)
        implicit none

        real*8, intent(IN) :: x0
        real*8, intent(IN) :: s

        real*8 :: gaussian_norm3

        gaussian_norm3 = 0.5d0 * (1+erf((x0/s)*(one_over_sqrt2)))

    end function gaussian_norm3

end module friction    

!*-----------------------------------------------------*!

subroutine calculate_nonadiabatic_friction(converged_scf)
    !  NAME
    !    calculate_nonadiabatic_friction 
    !  SYNOPSIS
    !  deallocates important dynamical arrays
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)

        use friction

        implicit none
        
        logical, intent(inout) :: converged_scf

        call friction_calculation(converged_scf)

end subroutine calculate_nonadiabatic_friction

subroutine friction_coordinate_workaround()
    !  NAME
    !    calculate_nonadiabatic_friction 
    !  SYNOPSIS
    !  fixes a problem that exists in the finite difference evaluation of friction
    !  AUTHOR
    !   Reinhard J. Maurer, Yale University (2015)

        use runtime_choices
        use dimensions
        use pbc_lists
        use geometry
        use friction
        use localorb_io, only: localorb_info, use_unit

        implicit none
         
        integer :: atom_counter, i_cart, i_atom
        real*8, dimension(3,n_atoms) :: coords_backup, frac_coords_tmp
        real*8, dimension(3) :: trans_corr_vec
        character(len=300) :: info_str 

        coords_backup(:,:) = coords(:,:)

        !atom_counter = 0
        trans_corr_vec = 0.50
        do i_atom=1, n_atoms 
                !atom_counter = atom_counter + 1
                do i_cart=1, 3,1
                    coords(i_cart,i_atom) = coords(i_cart,i_atom)+friction_numeric_disp
                    call cart2frac(lattice_vector,coords,frac_coords_tmp)
                    !check if frac_coords +- friction_numeric_disp crossed the 0.5 line 
                    if (frac_coords_tmp(i_cart,i_atom)>0.50 .and. frac_coords(i_cart,i_atom)<0.50) then
                        trans_corr_vec(i_cart) = -2.*friction_numeric_disp 
                    endif
                    coords(i_cart,i_atom) = coords(i_cart,i_atom)-2.0*friction_numeric_disp
                    call cart2frac(lattice_vector,coords,frac_coords_tmp)
                    if (frac_coords_tmp(i_cart,i_atom)<0.50 .and. frac_coords(i_cart,i_atom)>0.50) then
                        trans_corr_vec(i_cart) = 2.*friction_numeric_disp
                    endif
                    coords(i_cart,i_atom) = coords_backup(i_cart,i_atom)
                end do
        enddo
        do i_atom=1, n_atoms
            coords(:,i_atom) = coords(:,i_atom) + trans_corr_vec(:)
            !write(use_unit,*) coords(:,i_atom)
            !call map_to_center_cell(coords(:,i_atom))
            !coords_backup(:,i_atom) = coords_backup(:,i_atom) + trans_corr_vec(:)
        enddo
        !call cart2frac(lattice_vector,coords,frac_coords) 
        write (info_str,'(4X,A)') &
            "FRICTION workaraound shift of coordinates to avoid finite difference issue " 
        call localorb_info(info_str,use_unit,'(A)')
        
end subroutine friction_coordinate_workaround
