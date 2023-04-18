
! The Exact Two-Component (X2C) method for PBC systems. See J. Chem. Phys., 144, 044105 (2016).
! This file contains the driver subroutines for X2C and the major functions.
! -- Rundong Zhao, Nov. 2018 @ Durham, NC

 subroutine x2c_1st_iteration()
  use rel_x2c_mod
  use dimensions, only : n_k_points
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use constants, only : light_speed
  use physics
  use localorb_io, only: use_unit
  implicit none
  real*8, dimension(:,:,:), allocatable :: x2c_ss
  real*8, dimension(:,:,:), allocatable :: sesc_u 
  integer :: iesc=1 ! =1 NESC; =2 SESC.
  integer :: xdim, jsyml(2), norb(2)
  integer :: i,j,k,l,m,n
  real*8,allocatable :: cc(:,:),ee(:)

  xdim=dim_matrix_rel*dim_matrix_rel*4
  jsyml=1; norb=dim_matrix_rel*2

 ! 1. Transform STVW matrice to LOWER triangular form. Note, its LOWER !!!
  ! temporary arrays:
  call aims_allocate( x2c_ss,    dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_ss")
  call aims_allocate( sesc_u,   xdim, n_k_points, 2 ,"sesc_u" )
  ! globally standing arrays:
  call aims_allocate( x2c_t,     dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_t" )
  call aims_allocate( x2c_w,     dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_w" )
  call aims_allocate( x2c_v,     dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_v" )
  call aims_allocate( x2c_s,     dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_s" )
  call aims_allocate( x2c_v_old, dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_v_old" )
  call aims_allocate( x2c_x,     xdim, n_k_points, 2 ,"x2c_x" )
  call aims_allocate( x2c_fw,    xdim, n_k_points, 2 ,"x2c_fw" )
  x2c_t=0.d0; x2c_w=0.d0; x2c_ss=0.d0; sesc_u=0.d0
  x2c_v=0.d0; x2c_s=0.d0; x2c_x=0.d0; x2c_fw=0.d0

  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_v, x2c_v)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_t, x2c_t)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_w, x2c_w)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_ss,x2c_ss)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_s, x2c_s)

  x2c_v_old = x2c_v


 ! 2. Generate X matrix and NESC form kinetic and potential matrix:
 !    get X matrix X = BA^{-1}, by solving full Dirac equation

  do k=1, n_k_points

    allocate( cc(16*dim_matrix_rel*dim_matrix_rel,2), ee(4*dim_matrix_rel) )
  
   ! Solve full Dirac equation
    call x2c_calc_fdirac_solver(1,1,1, jsyml, norb, light_speed, &
         x2c_v_old(1,k,1), x2c_v_old(1,k,2), x2c_w(1,k,1), x2c_w(1,k,2), x2c_t(1,k,1), x2c_t(1,k,2), &
         x2c_s(1,k,1), x2c_s(1,k,2), x2c_ss(1,k,1), x2c_ss(1,k,2), cc, cc(1,2), ee)

                      ! write(use_unit,*)'Four-component 1st step eigenvalues:'
                      ! write(use_unit,"(10f18.8)")ee

   ! Get X matrix by equation: X = BA^{-1}
    call x2c_calc_xmat_generator(1,1,1, jsyml, norb, cc, cc(:,2), x2c_x(1,k,1), x2c_x(1,k,2) )

    deallocate( cc,ee )

   ! Generate FW transformation matrix and sesc tilde{S}S^{-1}; also NESC form kinetic matrix with new X matrix
    call x2c_gen_sescu_nesct_fw_mat( iesc, 1,1,1, jsyml, norb, light_speed, &
         x2c_s(1,k,1), x2c_s(1,k,2), x2c_t(1,k,1), x2c_t(1,k,2), &
         x2c_x(1,k,1), x2c_x(1,k,2), x2c_fw(1,k,1), x2c_fw(1,k,2), x2c_t(1,k,1), x2c_t(1,k,2), &
         x2c_ss(1,k,1), x2c_ss(1,k,2), sesc_u(1,k,1), sesc_u(1,k,2) )

    ! W matrix ==> XWX
    call x2c_w_to_xwx(1,1,1, jsyml, norb, x2c_x(1,k,1), x2c_x(1,k,2), x2c_w(1,k,1), x2c_w(1,k,2) ) ! STO scheme

  enddo ! end of k

 ! add XWX to V
  x2c_v = x2c_v_old + x2c_w


 ! 3. FW Transformation:


  do k=1, n_k_points

   ! FW Transformation for NESC T matrix:
    call x2c_w_to_xwx( 1,1,1, jsyml, norb, x2c_fw(1,k,1), x2c_fw(1,k,2), x2c_t(1,k,1), x2c_t(1,k,2) )

   ! V:
    call x2c_w_to_xwx( 1,1,1, jsyml, norb, x2c_fw(1,k,1), x2c_fw(1,k,2), x2c_v(1,k,1), x2c_v(1,k,2) )

  enddo

  ! Sum V and T (both FW-Transformed) to generate NESC H matrix:
  l = dim_matrix_rel*(2*dim_matrix_rel+1)*n_k_points*2
  call apbtoc(l, x2c_v, 1, x2c_t, 1, x2c_v, 1)

  call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
       occ_numbers, chemical_potential, chemical_potential_spin)

 !call aims_deallocate( x2c_t,  "x2c_t" )
 !call aims_deallocate( x2c_w,  "x2c_w" )
  call aims_deallocate( x2c_ss, "x2c_ss")
  call aims_deallocate( sesc_u, "sesc_u")

 end subroutine x2c_1st_iteration


 subroutine x2c_scf()
  use rel_x2c_mod
  use dimensions, only : n_k_points
  implicit none
  integer :: xdim, jsyml(2), norb(2)
  integer :: i,j,k,l,m,n

  jsyml=1; norb=dim_matrix_rel*2
  l=dim_matrix_rel*(2*dim_matrix_rel+1)

  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_v,x2c_v)

  x2c_v_old=x2c_v

 ! sum up to get total Fock matrix for X2C, first un-FW-trans Fock
  do k=1,n_k_points
     call apbtoc(l,x2c_v(1,k,1),1,x2c_w(1,k,1),1,x2c_v(1,k,1),1)
     call apbtoc(l,x2c_v(1,k,2),1,x2c_w(1,k,2),1,x2c_v(1,k,2),1)
  enddo


 ! FW-Transformation:
  do k=1,n_k_points
   ! FW-Transformation of V matrix:
    call x2c_w_to_xwx( 1,1,1,jsyml,norb, x2c_fw(1,k,1), x2c_fw(1,k,2), x2c_v(1,k,1), x2c_v(1,k,2) )
   !call x2c_w_to_xwx_nao( norb(1), x2c_fw(1,k,1), x2c_fw(1,k,2), x2c_v_sq(1,1,k,1), x2c_v_sq(1,1,k,2) )

   ! Add kinetic matrix (has been FW-Transformed) to NESC total Fock matrix:
    call apbtoc( l, x2c_v(1,k,1), 1, x2c_t(1,k,1), 1, x2c_v(1,k,1), 1 )
    call apbtoc( l, x2c_v(1,k,2), 1, x2c_t(1,k,2), 1, x2c_v(1,k,2), 1 )
  enddo
              ! if(scf_iteration.eq.2)stop

 end subroutine x2c_scf

 subroutine x2c_elsi_lapack_wrapper( i_k, i_k_point, n_spin, ham_work_c, ovlp_work_c, &
  flag_KS_k_points, n_states_k, KS_eigenvalue, KS_eigenvector_complex )
  use rel_x2c_mod, only : dim_matrix_rel, x2c_v, x2c_s, scf_iteration
  use elsi_wrapper, only : eh_scf, aims_elsi_set_illcond_check, aims_elsi_get_n_illcond, &
                           aims_elsi_set_output, aims_elsi_ev
  use dimensions, only : n_states, n_k_points, n_k_points_task
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use mpi_tasks, only: myid
  implicit none
  integer,intent(in) :: i_k, i_k_point, n_spin
  complex*16,intent(inout) :: ham_work_c(dim_matrix_rel*(2*dim_matrix_rel+1)), &
                              ovlp_work_c(dim_matrix_rel*(2*dim_matrix_rel+1))
  integer,intent(inout) :: flag_KS_k_points, n_states_k
  real*8,intent(out) :: KS_eigenvalue(n_states,n_spin,n_k_points)
  complex*16,intent(out) :: KS_eigenvector_complex(2*dim_matrix_rel,n_states,n_spin,n_k_points_task)

  integer :: i_spin, i,j,n
  complex*16,allocatable :: ovlp_work_cs(:,:), ham_work_cs(:,:), evec_work_cs(:,:)
  real*8,allocatable :: eval_tmp(:)

  character*50 :: file_name

  ! Note, x2c_v and x2c_s are saved in LOWER triangular form!
  ham_work_c(:) = dcmplx( x2c_v(:,i_k_point,1), -x2c_v(:,i_k_point,2) )
  ovlp_work_c(:) = dcmplx( x2c_s(:,i_k_point,1), -x2c_s(:,i_k_point,2) )

  ! Solve KS-equations for both spins
 !do i_spin = 1, n_spin
     call aims_allocate(ovlp_work_cs,2*dim_matrix_rel,2*dim_matrix_rel,'ovlp_work_cs')
     call aims_allocate(ham_work_cs,2*dim_matrix_rel,2*dim_matrix_rel,'ham_work_cs')

     ovlp_work_cs = (0.0d0,0.0d0)
     ham_work_cs = (0.0d0,0.0d0)

     ! Unpack packed matrices
     n = 0
     do i=1, 2*dim_matrix_rel
        do j=1,i
           n=n+1
           ovlp_work_cs(j,i) = ovlp_work_c(n)
           ovlp_work_cs(i,j) = dconjg(ovlp_work_c(n))
           ham_work_cs(j,i) = ham_work_c(n)
           ham_work_cs(i,j) = dconjg(ham_work_c(n))
        enddo
     enddo

     call aims_allocate(evec_work_cs,2*dim_matrix_rel,2*dim_matrix_rel,'evec_work_cs')
     call aims_allocate(eval_tmp,2*dim_matrix_rel,'eval_tmp')

     ! Skip singularity check when possible
     if(flag_KS_k_points == 0) then
        call aims_elsi_set_illcond_check(eh_scf,0)
     else
        call aims_elsi_set_illcond_check(eh_scf,1)
     endif

     call aims_elsi_ev(eh_scf,ham_work_cs,ovlp_work_cs,eval_tmp,evec_work_cs)

  do i_spin=1, n_spin 
     KS_eigenvalue(:,i_spin,i_k_point) = eval_tmp(1:n_states)
     KS_eigenvector_complex(1:2*dim_matrix_rel,1:n_states,i_spin,i_k) = evec_work_cs(1:2*dim_matrix_rel,1:n_states)
  enddo

     call aims_deallocate(ovlp_work_cs,'ovlp_work_cs')
     call aims_deallocate(ham_work_cs,'ham_work_cs')
     call aims_deallocate(evec_work_cs,'evec_work_cs')
     call aims_deallocate(eval_tmp,'eval_tmp')

     call aims_elsi_get_n_illcond(eh_scf,n_states_k)
     n_states_k = 2*dim_matrix_rel - n_states_k

     ! Overlap singular for this k-point?
     call aims_elsi_get_n_illcond(eh_scf,flag_KS_k_points)
     call aims_elsi_set_output(eh_scf,0)
 !enddo

 end subroutine x2c_elsi_lapack_wrapper


subroutine x2c_calc_fdirac_solver(ifrel,nspn,nsym,jsyml,norb,clight, &
  vr,vi,wr,wi,tr,ti,sr,si,ssr,ssi,cr,ci,ee)
!> solve Full Dirac equation from gotten V T S W matrix, output MO
!! coefficients and eigenvalues both negative and positive states
!! gotten here, saving by independant symmetry block order
  implicit none
  integer,intent(in)  :: ifrel,nspn
  integer,intent(in)  :: nsym,jsyml(nsym),norb(nsym)
  real*8, intent(in)  :: clight
  real*8, intent(in)  :: vr(*),vi(*),wr(*),wi(*),tr(*),ti(*),sr(*),si(*),ssr(*),ssi(*)
  real*8, intent(out) :: cr(*),ci(*),ee(*)
  real*8  :: fac
  integer :: k,q,p,nsos,isym

  fac = dsqrt(2.d0)*clight
  k   = 1
  q   = 1
  p   = 1
  do isym=1,nsym
    nsos=norb(isym)
    if(jsyml(isym).eq.0.or.nsos.le.0) cycle ! skip diagonalization for equivalent representations
   !if(ifrel.ne.0)then
      call x2c_fdirac_solveigen_onesym(nsos,fac,tr(k),ti(k),sr(k),si(k),ssr(k),ssi(k),vr(k), &
        vi(k),wr(k),wi(k),cr(q),ci(q),ee(p) )
   !else
   !  call x1c_fdirac_solveigen_onesym(nsos,fac,tr(k),sr(k),vr(k),wr(k),cr(q),ee(p))
   !endif
    k = k+nsos*(nsos+1)/2
    q = q+nsos*nsos*4
    p = p+nsos+nsos
  enddo
end subroutine x2c_calc_fdirac_solver

subroutine x2c_calc_fdirac_solver_nao(ifrel,nspn,nsym,jsyml,norb,clight, &
  vr,vi,wr,wi,tr,ti,tsr,tsi,sr,si,ssr,ssi,cr,ci,ee)
!> solve Full Dirac equation from gotten V T S W matrix, output MO
!! coefficients and eigenvalues both negative and positive states
!! gotten here, saving by independant symmetry block order
  implicit none
  integer,intent(in)  :: ifrel,nspn
  integer,intent(in)  :: nsym,jsyml(nsym),norb(nsym)
  real*8, intent(in)  :: clight
  real*8, intent(in)  :: vr(*),vi(*),wr(*),wi(*),tr(*),ti(*),tsr(*),tsi(*),sr(*),si(*),ssr(*),ssi(*)
  real*8, intent(out) :: cr(*),ci(*),ee(*)
  real*8  :: fac
  integer :: k,q,p,nsos,isym

  fac = dsqrt(2.d0)*clight
  k   = 1
  q   = 1
  p   = 1
  do isym=1,nsym
    nsos=norb(isym)
    if(jsyml(isym).eq.0.or.nsos.le.0) cycle ! skip diagonalization for equivalent representations

    call x2c_fdirac_solveigen_onesym_nao(nsos,fac,tr(k),ti(k),tsr(k),tsi(k),sr(k),si(k),ssr(k),ssi(k),vr(k), &
      vi(k),wr(k),wi(k),cr(q),ci(q),ee(p) )

    k = k+nsos*(nsos+1)/2
    q = q+nsos*nsos*4
    p = p+nsos+nsos
  enddo
end subroutine x2c_calc_fdirac_solver_nao

subroutine x2c_fdirac_solveigen_onesym(n,c2,tr,ti,sr,si,ssr,ssi,vr,vi,wr,wi,cr,ci,e)
!> solve FC=SCE:
! \f{eqnarray*} {
! && FC=SCE ==> S^{-1/2} F S^{-1/2} C'= C' E\\
! && C'=S^{1/2}C ==> C=S^{-1/2} C'
! \f}
  use localorb_io, only: use_unit
  implicit none
  integer,intent(in)  :: n
  real*8, intent(in)  :: tr(*),wr(*),vr(*),sr(*),ssr(*),c2 ! c2 = dsqrt(2)*clight
  real*8, intent(in)  :: ti(*),wi(*),vi(*),si(*),ssi(*)
  real*8, intent(out) :: cr(n*2,n*2),ci(n*2,n*2),e(n*2)
  complex*16, parameter   :: alph=dcmplx(1.d0,0.d0),beta=dcmplx(0.d0,0.d0)
  complex*16, allocatable :: sh(:,:),th(:,:),ff(:,:),s1(:,:),zwork(:)
  real*8,     allocatable :: rwork(:)
  integer :: lwork,info,n2
  integer :: i,j

  n2    = n+n
  allocate( sh(n,n), th(n,n), ff(n2,n2), s1(n2,n2) )
 ! get S^{-1/2} 
  call cvtzmx(n,sr,si,sh)
  call zmatsqrtz(sh,n,0)
 ! get Ss^{-1/2} for aims NAO scheme:
 !call cvtzmx(n,ssr,ssi,th)
 ! get T^{-1/2} for aims STO scheme:
  call cvtzmx(n,tr,ti,th)
  call zmatsqrtz(th,n,0)

  call x2c_calc_xmat_get_fdirac(n,vr,vi,tr,ti,wr,wi,ssr,ssi,sh,th,c2,ff,s1)
  deallocate( th,sh )

  lwork = 8*n2
  info  = 0
  allocate( sh(n2,n2),zwork(lwork), rwork(n2*3) )
  call zgemm('N','N',n2,n2,n2,alph,s1,n2,ff,n2, beta,sh,n2)
  call zgemm('N','N',n2,n2,n2,alph,sh,n2,s1,n2, beta,ff,n2)
  call zheev( 'V', 'U', n2, ff, n2, e, zwork, lwork, rwork, info)
  if( info.ne.0 )then
    write(use_unit,'(a,i5)') ' zheev failed, info=',info
    stop
  endif
  call zgemm('N','N',n2,n2,n2,alph,s1,n2,ff,n2, beta,sh,n2)
  cr = dreal(sh)
  ci = dimag(sh)
  deallocate(ff,s1,sh,zwork,rwork)

  ! (Rundong) It is necessary to check linear dependence for some bases:
  do i=1, 2*n
    if(dabs(e(i)).lt.1.d-5)then
      write(use_unit,"(a,a,a)")&
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',&
    '!!!!! Warning: there may be linear dependent bases, check the eigenvalues !!!!!',&
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    endif
  enddo
end subroutine x2c_fdirac_solveigen_onesym


subroutine x2c_fdirac_solveigen_onesym_nao(n,c2,tr,ti,tsr,tsi,sr,si,ssr,ssi,vr,vi,wr,wi,cr,ci,e)
!> solve FC=SCE:
! \f{eqnarray*} {
! && FC=SCE ==> S^{-1/2} F S^{-1/2} C'= C' E\\
! && C'=S^{1/2}C ==> C=S^{-1/2} C'
! \f}
  use localorb_io, only: use_unit
  implicit none
  integer,intent(in)  :: n
  real*8, intent(in)  :: tr(*),tsr(*),wr(*),vr(*),sr(*),ssr(*),c2 ! c2 = dsqrt(2)*clight
  real*8, intent(in)  :: ti(*),tsi(*),wi(*),vi(*),si(*),ssi(*)
  real*8, intent(out) :: cr(n*2,n*2),ci(n*2,n*2),e(n*2)
  complex*16, parameter   :: alph=dcmplx(1.d0,0.d0),beta=dcmplx(0.d0,0.d0)
  complex*16, allocatable :: sh(:,:),th(:,:),ff(:,:),s1(:,:),ee(:),zvl(:,:),zvr(:,:),zwork(:)
  real*8,     allocatable :: rwork(:)
  integer :: lwork,info,n2

  n2    = n+n
  allocate( sh(n,n), th(n,n), ff(n2,n2), s1(n2,n2), ee(n2), zvl(n2,n2), zvr(n2,n2) )
 ! get S^{-1/2} 
  call cvtzmx(n,sr,si,sh)
  call zmatsqrtz(sh,n,0)
 ! get Ss^{-1/2} 
  call cvtzmx(n,ssr,ssi,th)
  call zmatsqrtz(th,n,0)

  call x2c_calc_xmat_get_fdirac_nao(n,vr,vi,tr,ti,tsr,tsi,wr,wi,ssr,ssi,sh,th,c2,ff,s1)
  deallocate( th,sh )

  lwork = 8*n2
  info  = 0
  allocate( sh(n2,n2),zwork(lwork), rwork(n2*3) )
  call zgemm('N','N',n2,n2,n2,alph,s1,n2,ff,n2, beta,sh,n2)
  call zgemm('N','N',n2,n2,n2,alph,sh,n2,s1,n2, beta,ff,n2)
 !call zheev( 'V', 'U', n2, ff, n2, e, zwork, lwork, rwork, info)
  call zgeev( 'N', 'V', n2, ff, n2, ee, zvl,n2, zvr,n2, zwork, lwork, rwork, info)
  if( info.ne.0 )then
    write(use_unit,'(a,i5)') ' zheev failed, info=',info
    stop
  endif
 !call zgemm('N','N',n2,n2,n2,alph,s1,n2,ff,n2, beta,sh,n2)
  call zgemm('N','N',n2,n2,n2,alph,s1,n2,zvr,n2, beta,sh,n2)
  cr = dreal(sh)
  ci = dimag(sh)
  e  = dreal(ee)
  if(sum(dimag(ee)).gt.1.d-4)&
  write(use_unit,"('Warning: eigenvalues from 4C DKS are complex! Sum of the imaginary part:',f10.5)")sum(dimag(ee))
  call reorder_ascending_eigenvec(n2,n2,e,cr,ci)
  deallocate(ff,s1,sh,ee,zvl,zvr,zwork,rwork)
end subroutine x2c_fdirac_solveigen_onesym_nao


 subroutine q4c_1st_iteration()
  use rel_x2c_mod
  use dimensions, only : n_k_points
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use physics
  implicit none
  real*8, dimension(:,:,:), allocatable :: x2c_ss
  real*8, dimension(:,:,:), allocatable :: x2c_ssc2
  integer :: i,j,k,l,m,n

 ! 1. Transform STVW matrice to LOWER triangular form. Note, its LOWER !!!
  ! globally standing arrays:
  call aims_allocate( x2c_v_old, dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_v_old" )
  call aims_allocate( x2c_v,     dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_v" )
  call aims_allocate( x2c_w,     dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_w" )
  call aims_allocate( x2c_s,     dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_s" )
  ! temporary arrays:
  call aims_allocate( x2c_ss,    dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_ss")
  call aims_allocate( x2c_ssc2,  dim_matrix_rel*(2*dim_matrix_rel+1), n_k_points, 2 ,"x2c_ssc2")
  x2c_v=0.d0; x2c_w=0.d0; x2c_s=0.d0; x2c_ss=0.d0; x2c_ssc2=0.d0; x2c_v_old=0.d0

  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_v, x2c_v)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_w, x2c_w)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_s, x2c_s)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_ss,x2c_ss)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_ssc2,x2c_ssc2)

 ! 2. Construct H and S matrix:

  l = dim_matrix_rel*(2*dim_matrix_rel+1)*n_k_points*2
  call apbtoc(l, x2c_v_old, 1, x2c_ssc2, 1, x2c_v_old, 1)

  if(upw)then
    call apbtoc(l, x2c_v, 1, x2c_v_old, 1, x2c_v, 1)
    call apbtoc(l, x2c_v, 1, x2c_w,     1, x2c_v, 1) ! H matrix for diagonalization
  else
    call apbtoc(l, x2c_v_old, 1, x2c_w, 1, x2c_v_old, 1)
    call apbtoc(l, x2c_v, 1, x2c_v_old, 1, x2c_v,     1) ! H matrix for diagonalization
  endif
  call apbtoc(l, x2c_s, 1, x2c_ss, 1, x2c_s, 1) ! S matrix

 ! 3. Diagonalization:

  call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
       occ_numbers, chemical_potential, chemical_potential_spin)

  if(.not.upw) call aims_deallocate( x2c_w, "x2c_w")
  call aims_deallocate( x2c_ss,   "x2c_ss")
  call aims_deallocate( x2c_ssc2, "x2c_ssc2")

  if(.not.upw) call aims_deallocate( dirac_w, "dirac_w")
  call aims_deallocate( dirac_ss,   "dirac_ss")
  call aims_deallocate( dirac_ssc2, "dirac_ssc2")
  call aims_deallocate( dirac_s,    "dirac_s")

 end subroutine q4c_1st_iteration

 subroutine q4c_scf()
  use rel_x2c_mod
  use dimensions, only : n_k_points
  implicit none
  integer :: i,j,k,l,m,n

  l=dim_matrix_rel*(2*dim_matrix_rel+1)*n_k_points*2

  call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_v,x2c_v)

  call apbtoc(l, x2c_v_old, 1, x2c_v, 1, x2c_v, 1) ! H matrix for diagonalization

  if(upw)then
    call pbc_trans_fabc_to_fflt(dim_matrix_rel,n_k_points,dirac_w, x2c_w)
    call apbtoc(l, x2c_w, 1, x2c_v, 1, x2c_v, 1)
  endif

 end subroutine q4c_scf


 subroutine rel_band_plot_k(i_k_point, KS_eigenvalue)
  use rel_x2c_mod
  use dimensions, only: n_k_points, n_spin, n_states
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  implicit none
  integer,intent(in) :: i_k_point
  real*8,intent(out) :: KS_eigenvalue(n_states,n_spin)
  real*8, allocatable :: smihalf(:,:), es(:)
  real*8, dimension(:,:,:), allocatable :: x2c_ss
  real*8, dimension(:,:,:), allocatable :: x2c_ssc2
  real*8 :: eigenvector(2*dim_matrix_rel*2*dim_matrix_rel,2), eigenval(2*dim_matrix_rel)
  integer :: xdim
  integer :: i,j,k,l,m,n

 ! 1. Transform STVW matrice to LOWER triangular form. Note, its LOWER !!!

  if(allocated(x2c_v))    call aims_deallocate(x2c_v,    "x2c_v")
  if(allocated(x2c_w))    call aims_deallocate(x2c_w,    "x2c_w")
  if(allocated(x2c_s))    call aims_deallocate(x2c_s,    "x2c_s")
  if(allocated(x2c_ss))   call aims_deallocate(x2c_ss,   "x2c_ss")
  if(allocated(x2c_ssc2)) call aims_deallocate(x2c_ssc2, "x2c_ssc2")
  call aims_allocate( x2c_v,     dim_matrix_rel*(2*dim_matrix_rel+1), 1, 2 ,"x2c_v" )
  call aims_allocate( x2c_w,     dim_matrix_rel*(2*dim_matrix_rel+1), 1, 2 ,"x2c_w" )
  call aims_allocate( x2c_s,     dim_matrix_rel*(2*dim_matrix_rel+1), 1, 2 ,"x2c_s" )
  call aims_allocate( x2c_ss,    dim_matrix_rel*(2*dim_matrix_rel+1), 1, 2 ,"x2c_ss")
  call aims_allocate( x2c_ssc2,  dim_matrix_rel*(2*dim_matrix_rel+1), 1, 2 ,"x2c_ssc2")
  x2c_v=0.d0; x2c_w=0.d0; x2c_s=0.d0; x2c_ss=0.d0; x2c_ssc2=0.d0

  call pbc_trans_fabc_to_fflt(dim_matrix_rel,1,dirac_v(1,1,i_k_point,1),    x2c_v)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,1,dirac_w(1,1,i_k_point,1),    x2c_w)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,1,dirac_s(1,1,i_k_point,1),    x2c_s)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,1,dirac_ss(1,1,i_k_point,1),   x2c_ss)
  call pbc_trans_fabc_to_fflt(dim_matrix_rel,1,dirac_ssc2(1,1,i_k_point,1), x2c_ssc2)


 ! 2. Construct H and S matrix:

  l = dim_matrix_rel*(2*dim_matrix_rel+1)*2
  call apbtoc(l, x2c_v, 1, x2c_w,    1, x2c_v, 1)
  call apbtoc(l, x2c_v, 1, x2c_ssc2, 1, x2c_v, 1)  ! H matrix for diagonalization
  call apbtoc(l, x2c_s, 1, x2c_ss,   1, x2c_s, 1)  ! S matrix


 ! 3. Diagonalization:
  allocate( smihalf(2*dim_matrix_rel*2*dim_matrix_rel,2), es(2*dim_matrix_rel) )
  call pbc_diagmat_SX(.false., 2*dim_matrix_rel, 1.d-12, x2c_s, x2c_s(1,1,2), es, smihalf,smihalf(1,2), xdim )
  ! solve FC = SCE
  call pbc_diagmat_F(.false., 2*dim_matrix_rel,xdim, x2c_v, x2c_v(1,1,2), smihalf,smihalf(1,2), eigenvector,eigenvector(1,2),eigenval)
  deallocate( smihalf, es )

  KS_eigenvalue(1:n_states, 1)      = eigenval(1:n_states)
  KS_eigenvalue(1:n_states, n_spin) = eigenval(1:n_states)

  call aims_deallocate( x2c_v,    "x2c_v" )
  call aims_deallocate( x2c_w,    "x2c_w" )
  call aims_deallocate( x2c_s,    "x2c_s" )
  call aims_deallocate( x2c_ss,   "x2c_ss")
  call aims_deallocate( x2c_ssc2, "x2c_ssc2")

 end subroutine rel_band_plot_k



 
