!****s* FHI-aims/get_KS_orbitals_bandplot
!  NAME
!   get_KS_orbitals_bandplot
!  SYNOPSIS

subroutine get_KS_orbitals_bandplot(overlap_matrix,hamiltonian, n_electrons, &
               KS_eigenvalue, KS_eigenvector_complex, occ_numbers, band_par)

!  PURPOSE
!  The subroutine calculates Kohn-Sham eigenvalues and eigenvectors
!  for a given set (nonuniform) of k points, relavant for band plotting.
!  USES

  use dimensions
  use geometry
  use basis
  use runtime_choices
  use localorb_io
  use synchronize_mpi
  use lapack_wrapper
  use crpa_blacs
  use elsi_wrapper, only: aims_elsi_occ
  use pbc_lists, only: k_weights
  use ks_wrapper, only: solve_KS_elsi_serial
  implicit none

!  ARGUMENTS

  real*8 :: n_electrons
  real*8 :: hamiltonian(n_hamiltonian_matrix_size, n_spin )
  real*8 :: overlap_matrix(n_hamiltonian_matrix_size )

  real*8,     dimension(n_states, n_spin,n_k_points) :: KS_eigenvalue
  real*8,     dimension(n_states, n_spin,n_k_points) :: occ_numbers

  complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task) ::  KS_eigenvector_complex

  logical:: band_par

! INPUTS
! o n_electrons -- number of electrons
! o overlap_matrix -- overlap matrix
! o hamiltonian -- Hamiltonian matrix
!
! OUTPUT
! o KS_eigenvalue -- Kohn-Sham eigenvalues for a set of special k points
! o occ_numbers -- occpation numbers for a set of special k points
! o KS_eigenvector_complex -- Kohn-Sham eigenvectors if complex number eigenvectors are in use.
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


  real*8 chemical_potential_tmp
  real*8,    dimension(:,:),allocatable :: hamiltonian_w
  real*8,    dimension(:),  allocatable :: overlap_matrix_w
  complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
  complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex

  real*8, dimension(:,:,:),allocatable :: work_ham
  real*8, dimension(:,:),allocatable :: work_ovl
  complex*16, dimension(:), allocatable, save :: ovlp_work
  complex*16, dimension(:,:), allocatable, save :: ovlp_work_full
  complex*16, dimension(:,:), allocatable, save :: hamiltonian_work

  integer, save :: n_nonsingular
  real*8, dimension(:), allocatable, save :: ovlp_eigenvalues
  complex*16, dimension(:,:), allocatable, save :: ovlp_transform
  complex*16, dimension(:), allocatable, save :: trafo_hamiltonian
  real*8 :: ev_sqrt
  complex*16, dimension(:,:), allocatable, save :: aux_eigenvector

  integer:: info
  integer:: n_tasks_par, myid_par, comm_par

  complex*16 :: one = (1.0d0, 0.0d0)
  complex*16 :: zero = (0.0d0, 0.0d0)

  !  counters

  integer :: i_spin, i_k_point
  integer :: i_k
  integer :: i_trafo_1
  integer :: i_basis_1, i_basis_2
  integer :: dim_work, i_index

  character(*), parameter :: func='get_KS_orbitals_bandplot'

!  begin work

  if (real_eigenvectors) then
    call aims_stop('KS orbitals are expected to be complex in band structure.',func)
  endif

  if (band_par) then
     myid_par=myid_bl
     n_tasks_par=n_tasks_bl
     comm_par=comm_blacs
  else
     myid_par=myid
     n_tasks_par=n_tasks
     comm_par=mpi_comm_world
  end if

  allocate(hamiltonian_w_complex(n_basis*(n_basis+1)/2,n_spin),stat=info)
  call check_allocation(info,'hamiltonian_w_complex',func)
  allocate(overlap_matrix_w_complex(n_basis*(n_basis+1)/2),stat=info)
  call check_allocation(info,'overlap_matrix_w_complex',func)
  allocate(hamiltonian_w(1,1))
  allocate(overlap_matrix_w(1))

  if(packed_matrix_format == PM_none)then
      allocate(work_ham(n_centers_basis_I, n_centers_basis_I, n_spin))
      call check_allocation(info,'work_ham',func)
      allocate(work_ovl(n_centers_basis_I, n_centers_basis_I))
      call check_allocation(info,'work_ovl',func)
  else
     ! dummy only, never touched
      allocate(work_ham( 1, 1, 1))
      allocate(work_ovl( 1, 1))
  end if

!  flag_KS_k_points(:) = BASIS_SINGULAR_NOT_TESTED
  i_k=0
  do i_k_point = 1, n_k_points, 1
!     if(myid.eq.0) then
!       write(use_unit,*) "flag_KS_k_points ", i_k_point, flag_KS_k_points(i_k_point)
!     endif
     if(myid_par .eq. MOD(i_k_point, n_tasks_par)) then
        i_k = i_k + 1
        call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
             hamiltonian_w, overlap_matrix_w, &
             hamiltonian_w_complex, overlap_matrix_w_complex, &
             i_k_point, work_ham, work_ovl)

!        if(i_k_point.eq.2) then
!          write(use_unit,*) i_k_point
!          write(use_unit,'(100(6f16.8,/))') overlap_matrix_w_complex
!          write(use_unit,*)
!          write(use_unit,'(100(6f16.8,/))') hamiltonian_w_complex(:,1)
!        endif
        if(allocated(ovlp_complex)) deallocate(ovlp_complex)
        do i_spin = 1, n_spin, 1

           if(allocated(ovlp_complex)) deallocate(ovlp_complex)

           if(use_elsi .and. .not. use_scalapack) then
              call solve_KS_elsi_serial(hamiltonian_w_complex(:,i_spin),&
                   overlap_matrix_w_complex,KS_eigenvalue(:,i_spin,i_k_point),&
                   KS_eigenvector_complex(:,:,i_spin,i_k),i_spin,i_k_point)
           else
              call improve_complex_eigenfunctions(overlap_matrix_w_complex,&
                   hamiltonian_w_complex(:,i_spin),&
                   KS_eigenvalue(:,i_spin,i_k_point),&
                   KS_eigenvector_complex(:,:,i_spin,i_k),i_k_point)
           end if
        enddo

        if (allocated(ovlp_eigenvalues)) then
              deallocate(ovlp_eigenvalues)
        end if
        if (allocated(ovlp_transform)) then
            deallocate( ovlp_transform )
        end if
     else
       KS_eigenvalue(:,:,i_k_point) = 0.d0
     endif
  enddo

  call sync_vector(KS_eigenvalue,n_states*n_spin*n_k_points,comm_par)


  if(mu_method == 0) then ! zeroin
     call get_occupation_numbers_p0(KS_eigenvalue,n_electrons,.false.,&
             occ_numbers,chemical_potential_tmp)
  else ! bisection
     call aims_elsi_occ(n_electrons,n_states,n_spin,n_k_points,k_weights,&
             KS_eigenvalue,occ_numbers,chemical_potential_tmp)
  endif

!  if(myid.eq.0) then
!    do i_k_point = 1, n_k_points, 1
!      write(use_unit,*) "i_k_point:", i_k_point
!      write(use_unit,*) KS_eigenvalue(:,1,i_k_point)
!      write(use_unit,*)
!    enddo
!  endif

  deallocate(hamiltonian_w_complex)
  deallocate(overlap_matrix_w_complex)
  deallocate(hamiltonian_w)
  deallocate(overlap_matrix_w)
  deallocate(work_ham)
  deallocate(work_ovl)


end subroutine get_KS_orbitals_bandplot
!******
