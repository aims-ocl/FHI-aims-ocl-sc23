!****h* FHI-aims/trust_radius_method
!  NAME
!    trust_radius_method
!  SYNOPSIS

module trust_radius_method

  !  PURPOSE
  !   * trm_trusted_step() predicts a new geometry step DX from the forces F
  !     and the approximate Hessian H within the trust radius Delta.
  !   * trm_BFGS_update() updates an approximate Hessian to fulfill the
  !     secant equation H DX = DF.
  !  USES
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
  !    Release version, FHI-aims (2010).
  !  SOURCE
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use constants
  use distributed_hessian
  use elsi_wrapper, only: aims_elsi_stdevp
  use localorb_io, only: localorb_info
  use numerical_utilities, only: diagonalize_rmatrix, twonorm
  use runtime_choices, only: use_distributed_hessian, relax_geo, RELAX_LTRM_dev
  use symmetry_constrained_relaxation, only: use_symm_const_geo_relaxation
  use synchronize_mpi, only: sync_vector_trm
  use dimensions, only: n_periodic, relax_unit_cell, use_symmetric_forces

  implicit none

  real*8, parameter, private    :: separator = 1.d-10
  character*150, private        :: info_str

  private bum

  ! ------------------------------------------------------------------------ !

contains

  ! ------------------------------------------------------------------------ !
  ! -------------------------------- driver -------------------------------- !
  ! ------------------------------------------------------------------------ !

  subroutine trm_trusted_step(F, H, Delta, DX)
    implicit none
    real*8, intent(in)                  :: F(ncol_hess)
    real*8, intent(in)                  :: H(nrow_hess,ncol_hess)
    real*8, intent(in)                  :: Delta
    real*8, intent(out)                 :: DX(ncol_hess)
    ! Return the minimum of
    !   E~(DX) = -dot_product(F, DX) + 0.5*dot_product(DX, matmul(H, DX))
    ! within a radius of Delta:
    !   dot_product(DX, DX) <= Delta**2.
    !-
    real*8, allocatable                 :: S(:,:), D(:)
    real*8, allocatable                 :: StF(:), StDX(:)
    real*8, allocatable                 :: tmp_real(:,:)
    integer :: i
    character(*), parameter             :: func = 'trm_trusted_step'

    ! --- initialize
    call aims_allocate(S, nrow_hess, ncol_hess, "S")
    call aims_allocate(D, ncol_hess, "D")
    call aims_allocate(StF, ncol_hess, "StF")
    call aims_allocate(StDX, ncol_hess, "StDX")

    S    = H
    D    = 0.d0
    StF  = 0.d0
    StDX = 0.d0
    DX   = 0.d0

    ! --- get eigenbasis and transform
    if (.not. use_distributed_hessian) then
       call diagonalize_rmatrix(ncol_hess, S, D, .true.)
    else
       if (is_worker) then
          if (ncol_hess == 1) then
             D(1) = S(1,1)
             S(1,1) = 1.d0
          else
             call aims_allocate(tmp_real, nrow_hess, ncol_hess, "tmp_real")

             call aims_elsi_stdevp(ncol_hess, nrow_hess, ncol_hess, ncol_hess, &
                     S, D, tmp_real, ctxt_hess, blk_hess, comm_hess)

             S = tmp_real

             call aims_deallocate(tmp_real, "tmp_real")
          end if
       end if
    end if

    ! StF = matmul(transpose(S), F)
    if (.not. use_distributed_hessian) then
       call dgemv('T', ncol_hess, ncol_hess, 1.d0, S, ncol_hess, F, 1, 0.d0, &
               StF, 1)
    else
       if (is_worker) then
          call pdgemv("T", ncol_hess, ncol_hess, 1.d0, S, 1, 1, desc_hess, &
                  F(1+offset), 1, 1, desc_hess, 1, 0.d0, StF(1+offset), 1, 1, &
                  desc_hess, 1)
       end if

       call sync_vector_trm(is_worker, n_worker, blk_hess, nrow_hess, &
               ncol_hess, StF)
    end if
    if (is_worker) then
       ! --- diagnose Hessian
       call trm_examine_Hess(D, StF)

       ! --- TRM
       call trm_TRM(D, StF, StDX, Delta)
       ! DX = matmul(S, StDX)
       call dgemv('N', nrow_hess, ncol_hess, 1.d0, S, nrow_hess, StDX, 1, &
               0.d0, DX(1+offset), 1)
    end if

    if (use_distributed_hessian) then
       call sync_vector_trm(is_worker, n_worker, blk_hess, nrow_hess, &
               ncol_hess, DX)
    end if

    ! --- tidy up
    call aims_deallocate(S, "S")
    call aims_deallocate(D, "D")
    call aims_deallocate(StF, "StF")
    call aims_deallocate(StDX, "StDX")

  end subroutine trm_trusted_step

  ! ------------------------------------------------------------------------ !

  subroutine trm_BFGS_update(DF, DX, H)
    implicit none
    real*8,  intent(in)                 :: DF(ncol_hess)
    real*8,  intent(in)                 :: DX(ncol_hess)
    real*8,  intent(inout)              :: H(nrow_hess,ncol_hess)
    ! see e.g. Bakken and Helgaker, JCP 117, p9160 (2002)
    !-
    real*8                              :: DXDF, DXHDX
    real*8, allocatable                 :: HDX(:)
    integer                             :: i, j
    real*8, external                    :: ddot
    character(*), parameter             :: func = 'trm_BFGS_update'

    call aims_allocate(HDX, ncol_hess, "HDX")

    HDX = 0.d0
    ! HDX = matmul(H, DX)
    if (is_worker) then
       call dgemv('N', nrow_hess, ncol_hess, 1.d0, H, nrow_hess, DX, 1, 0.d0, &
               HDX(1+offset), 1)
    end if

    if (use_distributed_hessian) then
       call sync_vector_trm(is_worker, n_worker, blk_hess, nrow_hess, &
               ncol_hess, HDX)
    end if

    ! DXDF = dot_product(DX, DF)
    DXDF = ddot(ncol_hess, DX, 1, DF, 1)

    ! DXHDX = dot_product(DX, HDX)
    DXHDX = ddot(ncol_hess, DX, 1, HDX, 1)
    do j = 1, ncol_hess
       do i = 1, nrow_hess
          H(i, j) = H(i, j) - DF(i+offset)*DF(j)/DXDF - &
                       HDX(i+offset)*HDX(j)/DXHDX
       end do
    end do

    call aims_deallocate(HDX, "HDX")

  end subroutine trm_BFGS_update

  ! ------------------------------------------------------------------------ !
  ! -------------------------------- helpers ------------------------------- !
  ! ------------------------------------------------------------------------ !

  subroutine trm_examine_Hess(D, StF)
    implicit none
    real*8, intent(in)                  :: D(ncol_hess)
    real*8, intent(inout)               :: StF(ncol_hess)
    ! Check for definiteness and pull out zero-mode StF
    !-
    integer                             :: n_neg, n_zero, n_pos
    real*8                              :: P0F, D0
    character(*), parameter             :: func = 'trm_examine_Hess'

    n_neg  = count(D < - separator)
    n_zero = count(D <   separator) - n_neg
    n_pos  = ncol_hess - n_neg - n_zero

    if (n_pos < ncol_hess) then
       write(info_str, "(2X,'| Hessian has ',I2,' negative and ',I2,' zero eigenvalues.')") &
       & n_neg, n_zero
       call localorb_info(info_str)
       if (n_pos > 0) then
          write(info_str, &
          &     "(2X, '| Positive eigenvalues (eV/A^2):',ES9.2,' ...',ES9.2)") &
          & D(n_neg+n_zero+1)*hartree/bohr**2, D(ncol_hess)*hartree/bohr**2
          call localorb_info(info_str)
       end if
       if (n_neg > 0) then
          write(info_str, &
          &     "(2X, '| Negative eigenvalues (eV/A^2):',ES10.2,' ...',ES10.2)") &
          & D(1)*hartree/bohr**2, D(n_neg)*hartree/bohr**2
          call localorb_info(info_str)
       end if
    else
       write(info_str, &
       &     "(2X, '| Hessian eigenvalues (eV/A^2):',ES9.2,' ...',ES9.2)") &
       & D(n_neg+n_zero+1)*hartree/bohr**2, D(ncol_hess)*hartree/bohr**2
       call localorb_info(info_str)
    end if

    if (n_zero > 0) then
       P0F = twonorm(StF(n_neg+1:n_neg+n_zero))
       D0 = maxval(abs(D(n_neg+1:n_neg+n_zero)))
       if (P0F > separator) then
          write(info_str,"(2X, '| Force in zero mode:',ES10.2,' eV/A')") &
          & P0F*hartree/bohr
          call localorb_info(info_str)
          write(info_str, "(2X, '| where zero is',ES10.2,' eV/A^2')") &
          & D0*hartree/bohr**2
          call localorb_info(info_str)
       end if
       StF(n_neg+1:n_neg+n_zero) = 0.d0
    end if

  end subroutine trm_examine_Hess

  ! ------------------------------------------------------------------------ !

  subroutine trm_TRM(D, StF, StDX, Delta)
    implicit none
    real*8, intent(in)                  :: D(ncol_hess)
    real*8, intent(in)                  :: StF(ncol_hess)
    real*8, intent(out)                 :: StDX(ncol_hess)
    real*8, intent(in)                  :: Delta
    ! find root of
    !   g(lambda) = twonorm(DX)^2 - Delta^2
    !             = sum_i (DF_i / (H_ii + lambda))^2 - Delta^2
    ! using newtons method
    !-
    real*8                              :: lambda, lambda_lower, lambda_upper
    real*8                              :: g, dg
    integer                             :: i_cyc
    real*8                              :: minD, nDX0
    real*8, parameter                   :: lambda_upper_init = 1.d30
    integer, parameter                  :: max_lambda_cyc = 100
    character(*), parameter             :: func = 'trm_TRM'

    StDX = 0.d0
    where(abs(D) > separator) StDX = StF / D
    nDX0 = twonorm(StDX)

    if (Delta <= 0.d0) return ! TRM disabled
    minD = minval(D, mask=abs(StF)>0.d0)
    if (minD > 0.d0) then
       ! positive definite matrix -> use trivial result if within trust radius
       if (twonorm(StDX) < Delta) then
          write(info_str, &
          & "(2X,'| Use Quasi-Newton step of length |H^-1 F| =',ES9.2,' A.')") &
          & twonorm(StDX)*bohr
          call localorb_info(info_str)
          return
       end if
    end if

    if ((minD > 0.0d0) .and. (n_periodic > 0) &
        .and. (relax_unit_cell > 0) .and. (RELAX_LTRM_dev)) then
      ! FlK: in ltrm, just scale the step to lie within the trust radius:
      !         s_reduced = Delta * s / twonorm(s), see
      ! https://aims-git.rz-berlin.mpg.de/aims/FHIaims/merge_requests/14
      ! 14.12.18: Only use when hessian is positive definite (minD > 0)
      ! 03.02.19: Only use when non-periodic and unit cell relaxation is
      !           requested
      ! 04.03.19: Toggle by `RELAX_LTRM_dev` flag

      StDX = Delta * StDX / twonorm(StDX)

      write(info_str, "(2X,'| Step |H^-1 F|=',ES8.2,' A scaled down" // &
                      " to ',ES8.2,' A.')") nDX0*bohr, twonorm(StDX)*bohr
      call localorb_info(info_str)

      ! Information on how the step was reduced.
      write(info_str, "(2X,'| Hessian was scaled to reduce step." // &
                            " Preserves search direction.')")
      call localorb_info(info_str)
      return
    end if

    lambda_lower = max(0.d0, -minD)
    lambda_upper = lambda_upper_init
    lambda = min(lambda_lower + 0.5d0, 0.5d0 * (lambda_upper + lambda_lower))
    do i_cyc = 1, max_lambda_cyc
      StDX = StF / (D + lambda)
      g  = sum(StDX**2) - Delta**2
      dg = - 2.d0 * sum(StDX**2 / (D + lambda))

      if (abs(g/dg) < separator .or. abs(g) < 1.d-13) exit

      ! update valid interval (for nested intervals)
      ! g(lambda) should have negative slope
      if (g < 0.d0) then
         lambda_upper = min(lambda, lambda_upper)
      else
         lambda_lower = max(lambda, lambda_lower)
      end if
      if (dg > 0.d0 .or. lambda_lower > lambda_upper) then
         call bum(func, 'assertion failed')
      end if

      lambda = lambda - g / dg
      if (lambda <= lambda_lower .or. lambda >= lambda_upper) then
         ! if not allowed, use nested intervals
         lambda = 0.5d0 * (lambda_upper + lambda_lower)
      end if
    end do
    if (i_cyc > max_lambda_cyc) then
       write(info_str, &
       &     "('** TRM iteration not converged; last g/g'' =',ES10.2)") g/dg
       call localorb_info(info_str)
       write(info_str, &
       &     "('** Lowest eigenvalue % corr. StF:',2ES10.2)") D(1),StF(1)
       call localorb_info(info_str)
    end if

    ! log
    write(info_str, "(2X,'| Step |H^-1 F|=',ES8.2,' A reduced" // &
    &               " to |(H+',ES8.2,'eV/A^2)^-1 F|=',ES8.2,' A.')") &
    & nDX0*bohr, lambda*hartree/bohr**2, twonorm(StDX)*bohr
    call localorb_info(info_str)

  end subroutine trm_TRM

  ! ------------------------------------------------------------------------ !
  ! -------------------------- technical utilities ------------------------- !
  ! ------------------------------------------------------------------------ !

  subroutine bum(caller, text)
    use mpi_tasks, only: aims_stop
    implicit none
    character(*), intent(in)            :: caller
    character(*), intent(in)            :: text
    !-
    character*200                       :: info_str
    character(*), parameter             :: func = 'bum'

    write(info_str,"('*** bum: ',A,': ',A)") trim(caller), trim(text)
    call localorb_info(info_str)
    call aims_stop

  end subroutine bum

  ! ------------------------------------------------------------------------ !

end module trust_radius_method
