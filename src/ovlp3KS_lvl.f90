!****s* FHI-aims/ovlp3KS_lvl_1d
!  NAME
!    ovlp3KS_lvl_1d
!  SYNOPSIS

subroutine ovlp3KS_lvl_1d(n_KS_states, KS_eigenvector, &
&                         coeff_3fn_ten, coulomb_matr, ovlp_3KS)

  !  PURPOSE
  !
  !    Convenience interface to ovlp3KS_lvl_generic for 1D distribution over
  !    product basis.
  !
  !  USES

  use dimensions
  use prodbas
  use mpi_utilities
  use sparse_tensor
  use timing
  use runtime_choices, only: safe_minimum
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_KS_states
  real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points)
  type(sp_ten), intent(IN) :: coeff_3fn_ten
  real*8, intent(IN) :: coulomb_matr(n_basbas, n_loc_prodbas)
  real*8, intent(OUT) :: ovlp_3KS(n_loc_prodbas, n_states, n_KS_states, n_spin)

  !  INPUTS
  !    o n_KS_states -- the number of  KS states used in the transformation
  !                     (and in later calculations), second state dimension
  !                     of ovlp_3KS.
  !    o KS_eigenvector -- KS eigenvector (should already be broadcasted)
  !    o coeff_3fn_ten -- Sparse matrix of LVL product expansion coefficients
  !                       should correspond to 1D prodbas distribution.
  !    o coulomb_matr -- Coulomb metric of product basis (1d).
  !  OUTPUTS
  !    o ovlp_3KS -- Expansion coefficients of KS eigenstate products in
  !                  orthonormalized product basis (KSKSon).
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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  character(*), parameter :: func = 'ovlp3KS_lvl_1d'

  call ovlp3KS_lvl_generic(n_KS_states, KS_eigenvector, .false., &
  &                        coeff_3fn_ten, coulomb_matr, &
  &                        ovlp_3KS, n_loc_prodbas, n_states, n_KS_states)

end subroutine ovlp3KS_lvl_1d
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/ovlp3KS_lvl_2d
!  NAME
!    ovlp3KS_lvl_2d
!  SYNOPSIS

subroutine ovlp3KS_lvl_2d(n_KS_states, KS_eigenvector, &
&                         coeff_3fn_ten, coulomb_matr, ovlp_3KS)

  !  PURPOSE
  !
  !    Convenience interface to ovlp3KS_lvl_generic for 2D distribution over
  !    states.
  !
  !  USES

  use dimensions
  use prodbas
  use mpi_utilities
  use sparse_tensor
  use timing
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_KS_states
  real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points)
  type(sp_ten), intent(IN) :: coeff_3fn_ten
  real*8, intent(IN) :: coulomb_matr(n_basbas, n_loc_prodbas)
  real*8, intent(OUT) :: ovlp_3KS(n_basbas, ndim1_o3KS, ndim2_o3KS, n_spin)

  !  INPUTS
  !    o n_KS_states -- the number of  KS states used in the transformation
  !                     (and in later calculations), second state dimension
  !                     of ovlp_3KS.
  !    o KS_eigenvector -- KS eigenvector (should already be broadcasted)
  !    o coeff_3fn_ten -- Sparse matrix of LVL product expansion coefficients
  !                       should correspond to 1D prodbas distribution.
  !    o coulomb_matr -- Coulomb metric of product basis (1d).
  !  OUTPUTS
  !    o ovlp_3KS -- Expansion coefficients of KS eigenstate products in
  !                  orthonormalized product basis (KSKSon).
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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  character(*), parameter :: func = 'ovlp3KS_lvl_2d'

  call ovlp3KS_lvl_generic(n_KS_states, KS_eigenvector, .true., &
  &                        coeff_3fn_ten, coulomb_matr, &
  &                        ovlp_3KS, n_basbas, ndim1_o3KS, ndim2_o3KS)

end subroutine ovlp3KS_lvl_2d
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/ovlp3KS_lvl_generic
!  NAME
!    ovlp3KS_lvl_generic
!  SYNOPSIS

subroutine ovlp3KS_lvl_generic(n_KS_states, KS_eigenvector, use_2d_o3KS, &
&                              coeff_3fn_ten, coulomb_matr, &
&                              ovlp_3KS, ld1_o3KS, ld2_o3KS, ld3_o3KS)

  !  PURPOSE
  !
  !     Calculate ovlp_3KS (or coeff_KSKSon, i.e. coeff_3fn contracted twice
  !     with the KS eigenvectors and orthonormalized with respect to the
  !     Coulomb metric).  The rationale orthonormalizing right here is that
  !     sparsity is lost anyway, here, but not before.
  !
  !  USES

  use constants, only: pi
  use dimensions
  use prodbas
  use mpi_utilities
  use sparse_tensor
  use timing
  use runtime_choices, only: safe_minimum, use_scalapack
  use synchronize_mpi_basic, only: sync_vector
  use localorb_io, only: localorb_info, OL_low, OL_norm, use_unit
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_KS_states
  real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points)
  logical, intent(IN) :: use_2d_o3KS
  type(sp_ten), intent(IN) :: coeff_3fn_ten
  real*8, intent(IN) :: coulomb_matr(n_basbas, n_loc_prodbas)
  integer, intent(IN) :: ld1_o3KS, ld2_o3KS, ld3_o3KS
  real*8, intent(OUT) :: ovlp_3KS(ld1_o3KS, ld2_o3KS, ld3_o3KS, n_spin)

  !  INPUTS
  !    o n_KS_states -- the number of  KS states used in the transformation
  !                     (and in later calculations), second state dimension
  !                     of ovlp_3KS.
  !    o KS_eigenvector -- KS eigenvector (should already be broadcasted)
  !    o use_2d_o3KS -- Return ovlp_3KS in 1D prodbas or 2D state distribution?
  !    o coeff_3fn_ten -- Sparse matrix of LVL product expansion coefficients
  !                       should correspond to 1D prodbas distribution.
  !    o coulomb_matr -- Coulomb metric of product basis (1d).
  !    o ld?_o3KS -- Leading dimensions of ovlp_3KS
  !  OUTPUTS
  !    o ovlp_3KS -- Expansion coefficients of KS eigenstate products in
  !                  orthonormalized product basis (KSKSon).
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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer, parameter :: max_n_real = 2**23 ! ~ 8e6 reals, 64 MB
  real*8 :: time_Vsqrt(4), time_ass(4)
  real*8 :: time_fnKSbb(4), time_KSKSbb(4), time_KSKSon(4)
  integer :: n_prefac, n_state_dim, n_block
  real*8, allocatable :: coulomb_matr_sqrt(:,:), KS_vec(:,:)
  integer, allocatable :: basbas2row(:), basbas2col(:)
  real*8, allocatable :: coeff_fnKSbb(:,:,:)
  real*8, allocatable :: coeff_KSKSbb(:,:,:)
  real*8, allocatable :: coeff_KSKSon(:,:,:)
  real*8, allocatable :: glob_KSKSbb(:,:,:), glob_KSKSon(:,:,:)
  real*8, allocatable :: loc_KSKSon(:,:,:)
  integer :: i_spin, i_block, state_off, n_state_block
  integer :: i_loc_state, i_KS_state, i_state
  integer :: i_np1, i_loc_KS_state, i_basbas, n_rows
  integer, parameter :: dlen_ = 9
  integer :: n_KSlbb, n_KSKS, LDcoeff
  integer :: desc_KSKSbb(dlen_), desc_bbbb(dlen_)
  integer :: info, mpierror
  character*150 :: info_str
  real*8 :: fn_charge, max_charge_err
  real*8, allocatable :: charges(:,:)
  integer :: i_loc_prodbas, i_basbas_fn
  character(*), parameter :: func = 'ovlp3KS_lvl_generic'

  call localorb_info('Calculating RI-LVL ovlp_3KS', use_unit, '(2X,A)', OL_norm)

  time_Vsqrt = 0.d0; time_fnKSbb = 0.d0; time_KSKSbb = 0.d0; time_KSKSon = 0.d0
  time_ass = 0.d0

  ! This is crucial for the 2D distribution as ndim?_o3KS might be larger
  ! than the highest state.
  ovlp_3KS = 0.d0

  ! --- checks

  if (n_k_points > 1) call aims_stop('Multiple k-points not implemented', func)

  if (use_2d_o3KS) then
     if (ld1_o3KS /= n_basbas)   call aims_stop('ld1_o3KS [2d]', func)
     if (ld2_o3KS /= ndim1_o3KS) call aims_stop('ld2_o3KS [2d]', func)
     if (ld3_o3KS /= ndim2_o3KS) call aims_stop('ld3_o3KS [2d]', func)
  else
     if (ld1_o3KS /= n_loc_prodbas)   call aims_stop('ld1_o3KS [1d]', func)
     if (ld2_o3KS /= n_states)        call aims_stop('ld2_o3KS [1d]', func)
     if (ld3_o3KS /= n_KS_states)     call aims_stop('ld3_o3KS [1d]', func)
  end if

  ! --- buffer size

  ! Make sure that n_prefac is the same on every process.
  n_prefac = 0
  n_prefac = n_prefac + n_basis                        ! KS_vec
  n_prefac = n_prefac + n_basis * n_max_loc_prodbas    ! coeff_fnKSbb
  n_prefac = n_prefac + n_states * n_max_loc_prodbas   ! coeff_KSKSbb
  n_prefac = n_prefac + n_states * n_max_loc_prodbas   ! coeff_KSKSon
  if (.not. use_scalapack) then
     if (n_tasks > 2 .and. n_states*n_basbas > 2**20) then
        ! Large calculation (> 10^6 numbers) on truely parallel machine.
        call localorb_info('  ** Use scalapack to reduce memory footprint.')
     end if
     n_prefac = n_prefac + n_states * n_basbas         ! glob_KSKSbb
     n_prefac = n_prefac + n_states * n_basbas         ! glob_KSKSon
  end if
  if (use_2d_o3KS) then
     n_prefac = n_prefac + ndim1_o3KS * n_max_loc_prodbas  ! loc_KSKSon
     n_prefac = n_prefac + ndim1_o3KS * n_basbas           ! glob_KSKSon
  end if

  n_state_dim = floor(dble(max_n_real) / dble(n_prefac))
  n_state_dim = min(n_KS_states, n_state_dim)
  n_state_dim = max(n_state_dim, 1)
  n_block = int(ceiling(dble(n_KS_states) / dble(n_state_dim)))

  write(info_str, "(A,I6,A)") &
  & '| Using bunches of', n_state_dim, ' states for ovlp_3KS construction.'
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  write(info_str, "(A,I8,A)") &
  & '| Need', ceiling(dble(8 * n_state_dim * n_prefac) / 2.d0**20), &
  & ' MiB per proc for ovlp_3KS construction buffers.'
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  call localorb_info('', use_unit, '(2X,A)', OL_norm)




  ! --- coulomb_matr_sqrt

  call start_timer(time_Vsqrt)
  if (use_scalapack) then
     allocate(coulomb_matr_sqrt(n_basbas, n_loc_prodbas), stat=info)
     call check_allocation(info, 'coulomb_matr_sqrt', func)
     coulomb_matr_sqrt = coulomb_matr
     ! In principle, a Cholesky decomposition should suffice.  But for now...
     call power_auxmat_scalapack(coulomb_matr_sqrt, 0.5d0, 'Coulomb')
  else
     ! Put global V^0.5 into local array.
     allocate(coulomb_matr_sqrt(n_basbas, n_basbas), stat=info)
     call check_allocation(info, 'coulomb_matr_sqrt', func)
     call gather_auxmat(coulomb_matr_sqrt, coulomb_matr, n_basbas)
     call power_genmat_lapack(n_basbas, coulomb_matr_sqrt, 0.5d0, &
     &                        safe_minimum, 0.d0, 'Coulomb')
  end if
  call stop_timer(time_Vsqrt)

  ! --- allocate other arrays

  allocate(KS_vec(n_basis, n_state_dim), stat=info)
  call check_allocation(info, 'KS_vec', func)
  call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_low)
  allocate(coeff_fnKSbb(n_basis, n_state_dim, n_loc_prodbas), stat=info)
  call check_allocation(info, 'coeff_fnKSbb', func)
  allocate(coeff_KSKSbb(n_states, n_state_dim, n_loc_prodbas),stat=info)
  call check_allocation(info, 'coeff_KSKSbb', func)
  allocate(coeff_KSKSon(n_states, n_state_dim, n_loc_prodbas),stat=info)
  call check_allocation(info, 'coeff_KSKSon', func)
  allocate(basbas2row(n_basbas), basbas2col(n_basbas), stat=info)
  call check_allocation(info, 'basbas2xxx', func)
  call get_basbas_to_rowcol(.false., basbas2row, basbas2col)

  allocate(charges(n_states, n_state_dim), stat=info)
  call check_allocation(info, 'charges', func)

  max_charge_err = 0.d0
  do i_spin = 1, n_spin
     do i_block = 1, n_block   ! block of states
        state_off = (i_block-1) * n_state_dim
        n_state_block = min(n_state_dim, n_KS_states - state_off)

        call start_timer(time_fnKSbb)

        ! --- multiply coeff3fn (& transpose) to KS_eigenvector -> fnKSbb

        do i_loc_state = 1, n_state_block
           i_KS_state = state_off + i_loc_state
           KS_vec(:, i_loc_state) = KS_eigenvector(:, i_KS_state, i_spin, 1)
        end do

        call get_fnKSbb(n_state_block, n_state_dim, basbas2col, &
        &               coeff_3fn_ten, KS_vec, coeff_fnKSbb)

        call stop_timer(time_fnKSbb)

        ! --- fnKSbb -> KSKSbb

        call start_timer(time_KSKSbb)
        n_KSlbb = n_state_dim * n_loc_prodbas
        call dgemm('T', 'N', n_states, n_KSlbb, n_basis, &
        &          1.d0, KS_eigenvector(:,:, i_spin, 1), n_basis, &
        &                coeff_fnKSbb, n_basis, &
        &          0.d0, coeff_KSKSbb, n_states)
        call stop_timer(time_KSKSbb)

        ! --- Check orthonormalization of states

        ! JW: The array coeff_KSKSbb contains the expansion coefficients of
        ! the Kohn-Sham state products in terms of basbas functions.
        ! Therefore, the integral of a product of KS states over whole space
        ! is the weighted sum of integrals over basbas functions, i.e. the
        ! weighted sum over the basbas charges.  The latter are zero for L>0
        ! and proportional to the multipole momement for L=0.
        !
        ! The orthonormalization criterion says that the integral over
        ! products of states should be either one (i_state_1==i_state_2) or
        ! zero (i_state_1/=i_state_2).  Let's check this.
        !
        ! As we do not use charge conservation as a boundary condition when
        ! optimizing the LVL coefficients, the degree of orthonormality is a
        ! real measure for the quality of expansion.

        ! Calculate integrals over products of KS states
        charges = 0.d0
        do i_loc_prodbas = 1, n_loc_prodbas
           i_basbas = map_prodbas(i_loc_prodbas, myid+1)
           if (basbas_l(i_basbas) /= 0) cycle
           i_basbas_fn = basbas_fn(i_basbas)
           ! multipole = 1/(2l+1) * \int_0^\infty dr r^(2+l) f(r)
           ! where l=0 and f(\vec r) = f(r) * Y_{lm}(\hat r).
           ! Y_{00} = 1/sqrt(4pi), and factor 4*pi from angular integration.
           fn_charge = sqrt(4*pi)*multipole_basbas_fn(i_basbas_fn)
           charges = charges + fn_charge * coeff_KSKSbb(:,:, i_loc_prodbas)
        end do
        call sync_vector(charges, n_states*n_state_dim)

        ! Subtract identity matrix
        do i_loc_state = 1, n_state_block
           i_state = state_off + i_loc_state
           if (i_state <= n_states) then
              charges(i_state, i_loc_state) &
              & = charges(i_state, i_loc_state) - 1.d0
           end if
        end do
        max_charge_err = max(max_charge_err, maxval(abs(charges)))

        ! --- KSKSbb * V^0.5 -> KSKSon

        call start_timer(time_KSKSon)
        n_KSKS = n_states * n_state_block
        if (use_scalapack) then
           LDcoeff = n_states * n_state_dim
           call descinit(desc_KSKSbb, n_KSKS, n_basbas, mb_aux, nb_aux, &
           &             0, 0, my_blacs_ctxt_aux, LDcoeff, info)
           if (info /= 0) call aims_stop('descinit(desc_fnKSbb) failed', func)
           call descinit(desc_bbbb, n_basbas, n_basbas, mb_aux, nb_aux, &
           &             0, 0, my_blacs_ctxt_aux, n_basbas, info)
           if (info /= 0) call aims_stop('descinit(desc_bbbb) failed', func)

           call pdgemm('N', 'N', n_KSKS, n_basbas, n_basbas, &
           &           1.d0, coeff_KSKSbb, 1, 1, desc_KSKSbb, &
           &                 coulomb_matr_sqrt, 1, 1, desc_bbbb, &
           &           0.d0, coeff_KSKSon, 1, 1, desc_KSKSbb)

        else
           allocate(glob_KSKSbb(n_states, n_state_block, n_basbas), stat=info)
           call check_allocation(info, 'glob_KSKSbb', func)
           allocate(glob_KSKSon(n_states, n_state_block, n_basbas), stat=info)
           call check_allocation(info, 'glob_KSKSon', func)

           call gather_auxmat(glob_KSKSbb, coeff_KSKSbb(:,1:n_state_block,:), &
           &                  n_KSKS)
           if (myid == 0) then
              call dgemm('N', 'N', n_KSKS, n_basbas, n_basbas, &
              &          1.d0, glob_KSKSbb, n_KSKS, &
              &                coulomb_matr_sqrt, n_basbas, &
              &          0.d0, glob_KSKSon, n_KSKS)
           end if
           call scatter_auxmat(glob_KSKSon, coeff_KSKSon(:,1:n_state_block,:),&
           &                   n_KSKS)
           deallocate(glob_KSKSbb, glob_KSKSon)
        end if
        call stop_timer(time_KSKSon)

        ! --- KSKSon -> ovlp_3KS

        call start_timer(time_ass)
        if (use_2d_o3KS) then
           
           allocate(loc_KSKSon(ndim1_o3KS, n_state_block, n_loc_prodbas), &
           &        stat=info)
           call check_allocation(info, 'loc_KSKSon', func)
           allocate(glob_KSKSon(ndim1_o3KS, n_state_block, n_basbas), stat=info)
           call check_allocation(info, 'glob_KSKSon', func)

           do i_np1 = 0, np1_o3KS-1
              ! Prepare send buffer with all states for row i_np1
              loc_KSKSon = 0.d0
              do i_state = 1, n_states
                 if (own_dim1_o3KS(i_state) /= i_np1) cycle
                 i_loc_state = loc_dim1_o3KS(i_state)
                 loc_KSKSon(i_loc_state, :,:) = coeff_KSKSon(i_state, :,:)
              end do

              ! Allgather states for processor row i_np1
              n_rows = ndim1_o3KS * n_state_block
              call gather_auxmat(glob_KSKSon, loc_KSKSon, n_rows)
              call MPI_Bcast(glob_KSKSon, n_rows*n_basbas, &
              &              MPI_DOUBLE_PRECISION, 0, mpi_comm_global, mpierror)
              if (mpierror /= MPI_SUCCESS) call aims_stop('Bcast error', func)

              ! Save if I am on row i_np1
              if (i_np1 == myp1_o3KS) then
                 do i_loc_state = 1, n_state_block
                    i_KS_state = state_off + i_loc_state
                    if (own_dim2_o3KS(i_KS_state) == myp2_o3KS) then
                       i_loc_KS_state = loc_dim2_o3KS(i_KS_state)
                       do i_basbas = 1, n_basbas
                          ovlp_3KS(i_basbas,:,i_loc_KS_state, i_spin) &
                          & = glob_KSKSon(:, i_loc_state, i_basbas)
                       end do
                    end if
                 end do
              end if
           end do

           deallocate(loc_KSKSon, glob_KSKSon)

        else
           ! .not. use_2d_o3KS:
           do i_loc_state = 1, n_state_block
              i_KS_state = state_off + i_loc_state
              do i_state = 1, n_states
                 ovlp_3KS(:, i_state, i_KS_state, i_spin) &
                 & = coeff_KSKSon(i_state, i_loc_state, :)
              end do
           end do
        end if
        call stop_timer(time_ass)

     end do
  end do


  ! --- Output timing and Deallocation

  write(info_str, "(2X,'| Maximum charge error:',ES10.2)") max_charge_err
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('', use_unit, '(A)', OL_norm)

  call output_timeheader('2X', 'ovlp_3KS timings:', OL_norm)
  call output_timer('Time for V^0.5', time_Vsqrt(3:4), '2X', OL_norm)
  call output_timer('Time for fnfnbb -> fnKSbb', time_fnKSbb(3:4), '2X',OL_norm)
  call output_timer('Time for fnKSbb -> KSKSbb', time_KSKSbb(3:4), '2X',OL_norm)
  call output_timer('Time for KSKSbb -> KSKSon', time_KSKSon(3:4), '2X',OL_norm)
  call output_timer('Time for synchronization', time_ass(3:4), '2X', OL_norm)
  call localorb_info('', use_unit, '(A)', OL_low)

  deallocate(basbas2row, basbas2col)
  deallocate(KS_vec, coeff_fnKSbb, coeff_KSKSbb, coeff_KSKSon)
  deallocate(coulomb_matr_sqrt)

end subroutine ovlp3KS_lvl_generic
!******
