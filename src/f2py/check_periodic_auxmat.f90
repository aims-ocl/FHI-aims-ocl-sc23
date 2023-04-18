!****s* FHI-aims/check_periodic_auxmat
!  NAME
!    check_periodic_auxmat
!  SYNOPSIS

program check_periodic_auxmat

  !  PURPOSE
  !
  !  USES

  use prodbas
  use sbt_overlap_aims
  use tight_binding_auxmat
  use debug_output
  implicit none

  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  real*8, allocatable :: C1(:,:), C2(:,:), C3(:,:)
  integer :: ovlp_type
  integer :: info
  character*150 :: info_str, filename
  character*2 :: name
  integer :: i
  integer, parameter :: cell_a(3) =  (/ 0, 1, 0/)
  real*8 :: cellvec(3)
  real*8 :: ewald_gamma
  character(*), parameter :: func = 'check_periodic_auxmat'

  call init_aims()
  if (.not. use_prodbas) call aims_stop('No product basis from control.in')
  call initialize_bc_dependent_lists()


  allocate(C1(n_basbas, n_basbas), stat=info)
  call check_allocation(info, 'C1', func)
  allocate(C2(n_basbas, n_basbas), stat=info)
  call check_allocation(info, 'C2', func)
  allocate(C3(n_basbas, n_basbas), stat=info)
  call check_allocation(info, 'C3', func)

  cellvec = matmul(lattice_vector, dble(cell_a))

  write(0,*) 'n_k_points', n_k_points

  do i = 3, 4
     select case (i)
     case(1)
        ovlp_type = OVLP_TYPE_OVERLAP
        name = 'OO'
     case(2)
        ovlp_type = OVLP_TYPE_CUT
        name = 'CC'
     case(3)
        ewald_gamma = 1.d0
        ovlp_type = OVLP_TYPE_COULOMB
        name = 'C1'
     case(4)
        ewald_gamma = 2.d0
        ovlp_type = OVLP_TYPE_COULOMB
        name = 'C2'
     case(5)
        use_hse = .true.
        hse_omega_hf = 0.11d0
        ovlp_type = OVLP_TYPE_HSE
        call cleanup_basbas()
        call initialize_prodbas()
        name = 'HS'
     case default
        cycle
     end select

     write(0,*) '=== ', OVLP_TYPE_NAMES(ovlp_type)
     write(0,*) '* q-space'
     call qspace_and_back(ovlp_type, cellvec, C1)
     write(0,*) '* direct'
     call integrate_auxmat_by_atomic_sbt(C2, ovlp_type, .false., cell_a)
     if (ovlp_type /= OVLP_TYPE_COULOMB) then
        write(0,*) '* BvK'
        call all_BvK_cells(ovlp_type, cell_a, C3)
     end if

     write(0,"('max|q-space - direct|',ES10.2)") maxval(abs(C1-C2))
     if (ovlp_type /= OVLP_TYPE_COULOMB) then
        write(0,"('max|q-space - BvK|',ES10.2)") maxval(abs(C1-C3))
        write(0,"('max|dirct - BvK|',ES10.2)") maxval(abs(C2-C3))
     end if
     write(filename, "(A,'-qspace.dat')") name
     call debug_array(22, C1, tofile=trim(filename))
     write(filename, "(A,'-direct.dat')") name
     call debug_array(22, C2, tofile=trim(filename))
     if (ovlp_type /= OVLP_TYPE_COULOMB) then
        write(filename, "(A,'-BvK.dat')") name
        call debug_array(22, C3, tofile=trim(filename))
     end if
  end do


contains

  subroutine qspace_and_back(ovlp_type, cellvec, C)
    implicit none
    integer, intent(IN) :: ovlp_type
    real*8, intent(IN) :: cellvec(3)
    real*8, intent(OUT) :: C(n_basbas, n_basbas)

    complex*16, allocatable :: c_auxmats(:,:,:)
    complex*16, allocatable :: caux(:,:)
    real*8 :: qvecs(3, n_k_points)
    real*8 :: phi
    complex*16 :: phase
    integer :: info, i_k_point

    allocate(c_auxmats(n_basbas, n_basbas, n_k_points), stat=info)
    call check_allocation(info, 'c_auxmats', func)
    allocate(caux(n_basbas, n_basbas), stat=info)
    call check_allocation(info, 'caux', func)

    if (ovlp_type == OVLP_TYPE_COULOMB) then
       call initialize_periodic_tb_auxmat(1, ewald_gamma)
    else
       call initialize_tb_auxmat(1, ovlp_type)
    end if

    qvecs = matmul(recip_lattice_vector, transpose(k_point_list))

    call get_qspace_auxmat(n_k_points, qvecs, c_auxmats)

    call deallocate_tb_auxmat()

    caux = 0.d0
    do i_k_point = 1, n_k_points
       phi = dot_product(qvecs(:, i_k_point), cellvec)
       phase = cmplx(cos(phi), -sin(phi), kind(0.d0))
       caux = caux + phase * c_auxmats(:,:, i_k_point) / n_k_points
    end do
    if (any(abs(aimag(caux)) > 1d-12)) then
       write(0,*) 'imag(caux):', maxval(abs(aimag(caux)))
    end if


    C = real(caux, kind(0.d0))

  end subroutine qspace_and_back

  subroutine all_BvK_cells(ovlp_type, cell_a, C)
    use bravais
    implicit none
    integer, intent(IN) :: ovlp_type
    integer, intent(IN) :: cell_a(3)
    real*8, intent(OUT) :: C(n_basbas, n_basbas)

    integer :: n_supercell(3), a1, a2, a3, a123(3), i
    real*8 :: BvK(3,3)
    real*8 :: tmp_C(n_basbas, n_basbas)
    real*8 :: Rmax

    do i = 1, 3
       BvK(:, i) = n_k_points_xyz(i) * lattice_vector(:, i)
    end do

    select case (ovlp_type)
    case(OVLP_TYPE_OVERLAP)
       Rmax = 2*maxval(charge_radius_basbas_fn)
    case(OVLP_TYPE_HSE)
       Rmax = maxval(charge_radius_basbas_fn)+maxval(field_radius_basbas_fn)
    case(OVLP_TYPE_CUT)
       Rmax = 2*maxval(charge_radius_basbas_fn) + cutCb_rcut * cutCb_width**6
    case default
       call aims_stop('Invalid ovlp_type', func)
    end select
    call get_n_supercells(3, BvK, Rmax, n_supercell)

    C = 0.d0
    n_supercell = n_supercell
    do a1 = -n_supercell(1), n_supercell(1)
       do a2 = -n_supercell(2), n_supercell(2)
          do a3 = -n_supercell(3), n_supercell(3)
             a123 = (/a1, a2, a3/)
             call integrate_auxmat_by_atomic_sbt(tmp_C, ovlp_type, .false., &
             &                                   cell_a + a123*n_k_points_xyz)
             C = C + tmp_C
          end do
       end do
    end do

  end subroutine all_BvK_cells


end program check_periodic_auxmat
!******
