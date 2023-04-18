!****h* FHI-aims/separate_core_states
!  NAME
!    separate_core_states
!  SYNOPSIS
module separate_core_states
  !  PURPOSE
  !
  !  This module includes routines for the frozen core approximation.
  !  Ref: Koepernik and Eschrig, Full-potential nonorthogonal local-orbital
  !  minimum-basis band-structure scheme, Physical Review B 59 (1999) 1743.
  !
  !  USES
  use aims_memory_tracking, only: aims_allocate,aims_deallocate
  use basis, only: n_basis_atom,basis_atom,basis_fn,basisfn_type,basisfn_n,&
      basisfn_l
  use constants, only: hartree
  use dimensions, only: n_states,n_basis,n_core_states,n_spin,n_valence_basis
  use free_atoms, only: free_wave_eigenval
  use geometry, only: species
  use localorb_io, only: OL_norm,use_unit,localorb_info
  use mpi_tasks, only: aims_stop_coll
  use runtime_choices, only: frozen_core_scf,frozen_core_scf_factor,&
      frozen_core_scf_cutoff,use_scalapack,use_elsi_ev,real_eigenvectors,&
      sc_accuracy_rho
  use scalapack_wrapper, only: mxld,mxcol,mb,nb,sc_desc,ham,ovlp,eigenvec,&
      ham_complex,ovlp_complex,eigenvec_complex,my_scalapack_comm_work,&
      my_blacs_ctxt,myprow,mypcol,nprow,npcol,l_row,l_col
  use synchronize_mpi, only: sync_vector
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !    the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
  implicit none

  private

  public :: initialize_separate_core_states
  public :: deallocate_separate_core_states
  public :: permute_ham_scalapack
  public :: permute_ovlp_scalapack
  public :: permute_evec_scalapack
  public :: unpermute_ovlp_scalapack
  public :: fetch_core_eval_scalapack
  public :: permute_ham_lapack
  public :: permute_ovlp_lapack
  public :: permute_evec_lapack
  public :: fetch_core_eval_lapack

  interface permute_ham_lapack
     module procedure permute_ham_lapack_real
     module procedure permute_ham_lapack_cmplx
  end interface

  interface permute_ovlp_lapack
     module procedure permute_ovlp_lapack_real
     module procedure permute_ovlp_lapack_cmplx
  end interface

  interface permute_evec_lapack
     module procedure permute_evec_lapack_real
     module procedure permute_evec_lapack_cmplx
  end interface

  interface fetch_core_eval_lapack
     module procedure fetch_core_eval_lapack_real
     module procedure fetch_core_eval_lapack_cmplx
  end interface

  real*8, allocatable, public :: ham_real_vv(:,:)
  real*8, allocatable, public :: ovlp_real_vv(:,:)
  real*8, allocatable, public :: evec_real_vv(:,:)

  complex*16, allocatable, public :: ham_cmplx_vv(:,:)
  complex*16, allocatable, public :: ovlp_cmplx_vv(:,:)
  complex*16, allocatable, public :: evec_cmplx_vv(:,:)

  integer :: mxld_vv
  integer :: mxcol_vv
  integer :: sc_desc_vv(9)

  integer, allocatable :: perm(:)

  character*200 :: msg

contains

  !--------------------------------------------------------------
  !****s* separate_core_states/initialize_separate_core_states
  !  NAME
  !   initialize_separate_core_states
  subroutine initialize_separate_core_states

    implicit none

    integer :: i_func
    integer :: i_atom
    integer :: i_basis
    integer :: i_state
    integer :: i_core
    integer :: n_save
    integer :: l_save
    integer :: tmp

    integer, allocatable :: core_idx(:)

    integer, external :: numroc

    n_basis_atom = 0

    do i_basis = 1,n_basis
       n_basis_atom(basis_atom(i_basis)) = n_basis_atom(basis_atom(i_basis))+1
    end do

    n_core_states = 0
    n_valence_basis = n_basis
    i_func = 1
    n_save = 1
    l_save = 0

    ! Count core states
    ! TODO: There must be a better way of doing this
    do i_basis = 1,n_basis
       i_atom = basis_atom(i_basis)

       if(basisfn_type(basis_fn(i_basis)) == "atomic") then
          if(basisfn_n(basis_fn(i_basis)) /= n_save&
             .or. basisfn_l(basis_fn(i_basis)) /= l_save) then
             i_func = i_func+1
          end if

          n_save = basisfn_n(basis_fn(i_basis))
          l_save = basisfn_l(basis_fn(i_basis))

          if(free_wave_eigenval(species(i_atom),i_func)&
             <= frozen_core_scf_cutoff) then
             n_core_states = n_core_states+1
          end if
       end if

       if(i_atom /= basis_atom(min(i_basis+1,n_basis))) then
          i_func = 1
          n_save = 1
          l_save = 0
       end if
    end do

    if(n_core_states > 0 .and. frozen_core_scf) then
       if(.not. use_elsi_ev) then
          write(msg,"(X,A)") "***"
          call localorb_info(msg,use_unit,"(A)",OL_norm)
          write(msg,"(X,A)")&
             "*** ERROR: Frozen core does not support density matrix solvers."
          call localorb_info(msg,use_unit,"(A)",OL_norm)
          write(msg,"(X,A)") "***"
          call localorb_info(msg,use_unit,"(A)",OL_norm)
          call aims_stop_coll("Error in frozen core.")
       end if

       n_valence_basis = n_basis-n_core_states

       call aims_allocate(core_idx,n_core_states,"core_idx")
       call aims_allocate(perm,n_basis,"perm")

       i_core = 1
       i_func = 1
       n_save = 1
       l_save = 0
       core_idx = 0

       do i_basis = 1,n_basis
          i_atom = basis_atom(i_basis)

          if(basisfn_type(basis_fn(i_basis)) == "atomic") then
             if(basisfn_n(basis_fn(i_basis)) /= n_save&
                .or. basisfn_l(basis_fn(i_basis)) /= l_save) then
                i_func = i_func+1
             end if

             n_save = basisfn_n(basis_fn(i_basis))
             l_save = basisfn_l(basis_fn(i_basis))

             if(free_wave_eigenval(species(i_atom),i_func)&
                <= frozen_core_scf_cutoff) then
                core_idx(i_core) = i_basis
                i_core = i_core+1
             end if
          end if

          if(i_atom /= basis_atom(min(i_basis+1,n_basis))) then
             i_func = 1
             n_save = 1
             l_save = 0
          end if

          perm(i_basis) = i_basis
       end do

       do i_state = 1,n_core_states
          i_basis = core_idx(i_state)
          tmp = perm(i_state)
          perm(i_state) = perm(i_basis)
          perm(i_basis) = tmp
       end do

       call aims_deallocate(core_idx,"core_idx")

       write(msg,"(2X,A,F10.3,A)") "| Frozen core approximation cutoff: ",&
          frozen_core_scf_cutoff*hartree," eV"
       call localorb_info(msg,use_unit,"(A)",OL_norm)
       write(msg,"(2X,A,I10,A)") "| Found ",n_core_states," core states"
       call localorb_info(msg,use_unit,"(A)",OL_norm)
       if(frozen_core_scf_factor > 1.d0) then
          write(msg,"(2X,A,E13.6)")&
             "| They will be frozen until charge density converges to",&
             sc_accuracy_rho*frozen_core_scf_factor
          call localorb_info(msg,use_unit,"(A)",OL_norm)
       end if

       if(use_scalapack) then
          mxld_vv = numroc(n_valence_basis,mb,myprow,0,nprow)
          mxcol_vv = numroc(n_valence_basis,nb,mypcol,0,npcol)

          call descinit(sc_desc_vv,n_valence_basis,n_valence_basis,mb,nb,0,0,&
               my_blacs_ctxt,max(1,mxld_vv),tmp)

          if(real_eigenvectors) then
             call aims_allocate(ham_real_vv,mxld_vv,mxcol_vv,"ham_real_vv")
             call aims_allocate(ovlp_real_vv,mxld_vv,mxcol_vv,"ovlp_real_vv")
             call aims_allocate(evec_real_vv,mxld_vv,mxcol_vv,"evec_real_vv")
          else
             call aims_allocate(ham_cmplx_vv,mxld_vv,mxcol_vv,"ham_cmplx_vv")
             call aims_allocate(ovlp_cmplx_vv,mxld_vv,mxcol_vv,"ovlp_cmplx_vv")
             call aims_allocate(evec_cmplx_vv,mxld_vv,mxcol_vv,"evec_cmplx_vv")
          end if
       end if
    else
       frozen_core_scf = .false.
       n_core_states = 0
       n_valence_basis = n_basis
    end if

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/deallocate_separate_core_states
  !  NAME
  !   deallocate_separate_core_states
  subroutine deallocate_separate_core_states

    implicit none

    if(allocated(perm)) then
       call aims_deallocate(perm,"perm")
    end if

    if(allocated(ham_real_vv)) then
       call aims_deallocate(ham_real_vv,"ham_real_vv")
    end if

    if(allocated(ovlp_real_vv)) then
       call aims_deallocate(ovlp_real_vv,"ovlp_real_vv")
    end if

    if(allocated(evec_real_vv)) then
       call aims_deallocate(evec_real_vv,"evec_real_vv")
    end if

    if(allocated(ham_cmplx_vv)) then
       call aims_deallocate(ham_cmplx_vv,"ham_cmplx_vv")
    end if

    if(allocated(ovlp_cmplx_vv)) then
       call aims_deallocate(ovlp_cmplx_vv,"ovlp_cmplx_vv")
    end if

    if(allocated(evec_cmplx_vv)) then
       call aims_deallocate(evec_cmplx_vv,"evec_cmplx_vv")
    end if

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_ham_scalapack
  !  NAME
  !   permute_ham_scalapack
  subroutine permute_ham_scalapack(i_spin)

    implicit none

    integer, intent(in) :: i_spin

    integer :: i

    if(real_eigenvectors) then
       eigenvec(:,:,i_spin) = ham(:,:,i_spin)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pdcopy(n_basis,ham(:,:,i_spin),1,perm(i),sc_desc,1,&
                  eigenvec(:,:,i_spin),1,i,sc_desc,1)
          end if
       end do

       ham(:,:,i_spin) = eigenvec(:,:,i_spin)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pdcopy(n_basis,eigenvec(:,:,i_spin),perm(i),1,sc_desc,&
                  n_basis,ham(:,:,i_spin),i,1,sc_desc,n_basis)
          end if
       end do

       call pdtran(n_valence_basis,n_valence_basis,1.d0,ham(:,:,i_spin),&
            n_core_states+1,n_core_states+1,sc_desc,0.d0,ham_real_vv,1,1,&
            sc_desc_vv)

       ! H_vv = H_vv - S_vc * H_cc * S_cv
       ! For accuracy, cannot assume diagonal H_cc and do pdsyrk
       call pdgemm("N","N",n_valence_basis,n_core_states,n_core_states,1.d0,&
            ovlp,n_core_states+1,1,sc_desc,ham(:,:,i_spin),1,1,sc_desc,0.d0,&
            eigenvec(:,:,i_spin),1,1,sc_desc)

       call pdgemm("N","N",n_valence_basis,n_valence_basis,n_core_states,-1.d0,&
            eigenvec(:,:,i_spin),1,1,sc_desc,ovlp,1,n_core_states+1,sc_desc,&
            1.d0,ham_real_vv,1,1,sc_desc_vv)
    else
       eigenvec_complex(:,:,i_spin) = ham_complex(:,:,i_spin)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pzcopy(n_basis,ham_complex(:,:,i_spin),1,perm(i),sc_desc,1,&
                  eigenvec_complex(:,:,i_spin),1,i,sc_desc,1)
          end if
       end do

       ham_complex(:,:,i_spin) = eigenvec_complex(:,:,i_spin)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pzcopy(n_basis,eigenvec_complex(:,:,i_spin),perm(i),1,&
                  sc_desc,n_basis,ham_complex(:,:,i_spin),i,1,sc_desc,n_basis)
          end if
       end do

       call pztranc(n_valence_basis,n_valence_basis,(1.d0,0.d0),&
            ham_complex(:,:,i_spin),n_core_states+1,n_core_states+1,sc_desc,&
            (0.d0,0.d0),ham_cmplx_vv,1,1,sc_desc_vv)

       ! H_vv = H_vv - S_vc * H_cc * S_cv
       ! For accuracy, cannot assume diagonal H_cc and do pzherk
       call pzgemm("N","N",n_valence_basis,n_core_states,n_core_states,&
            (1.d0,0.d0),ovlp_complex,n_core_states+1,1,sc_desc,&
            ham_complex(:,:,i_spin),1,1,sc_desc,(0.d0,0.d0),&
            eigenvec_complex(:,:,i_spin),1,1,sc_desc)

       call pzgemm("N","N",n_valence_basis,n_valence_basis,n_core_states,&
            (-1.d0,0.d0),eigenvec_complex(:,:,i_spin),1,1,sc_desc,ovlp_complex,&
            1,n_core_states+1,sc_desc,(1.d0,0.d0),ham_cmplx_vv,1,1,sc_desc_vv)
    end if

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_ham_lapack_real
  !  NAME
  !   permute_ham_lapack_real
  subroutine permute_ham_lapack_real(ham_work,ovlp_work,work)

    implicit none

    real*8, intent(inout) :: ham_work(n_basis,n_basis)
    real*8, intent(in) :: ovlp_work(n_basis,n_basis)
    real*8, intent(inout) :: work(n_basis,n_basis)

    integer :: i

    work = ham_work

    do i = 1,n_basis
       if(i /= perm(i)) then
          work(:,i) = ham_work(:,perm(i))
       end if
    end do

    ham_work = work

    do i = 1,n_basis
       if(i /= perm(i)) then
          ham_work(i,:) = work(perm(i),:)
       end if
    end do

    ! H_vv = H_vv - S_vc * H_cc * S_cv
    ! For accuracy, cannot assume diagonal H_cc and do dsyrk
    call dgemm("N","N",n_valence_basis,n_core_states,n_core_states,1.d0,&
         ovlp_work(n_core_states+1,1),n_basis,ham_work,n_basis,0.d0,work,&
         n_basis)

    call dgemm("N","N",n_valence_basis,n_valence_basis,n_core_states,-1.d0,&
         work,n_basis,ovlp_work(1,n_core_states+1),n_basis,1.d0,&
         ham_work(n_core_states+1,n_core_states+1),n_basis)

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_ham_lapack_cmplx
  !  NAME
  !   permute_ham_lapack_cmplx
  subroutine permute_ham_lapack_cmplx(ham_work,ovlp_work,work)

    implicit none

    complex*16, intent(inout) :: ham_work(n_basis,n_basis)
    complex*16, intent(in) :: ovlp_work(n_basis,n_basis)
    complex*16, intent(inout) :: work(n_basis,n_basis)

    integer :: i

    work = ham_work

    do i = 1,n_basis
       if(i /= perm(i)) then
          work(:,i) = ham_work(:,perm(i))
       end if
    end do

    ham_work = work

    do i = 1,n_basis
       if(i /= perm(i)) then
          ham_work(i,:) = work(perm(i),:)
       end if
    end do

    ! H_vv = H_vv - S_vc * H_cc * S_cv
    ! For accuracy, cannot assume diagonal H_cc and do zherk
    call zgemm("N","N",n_valence_basis,n_core_states,n_core_states,(1.d0,0.d0),&
         ovlp_work(n_core_states+1,1),n_basis,ham_work,n_basis,(0.d0,0.d0),&
         work,n_basis)

    call zgemm("N","N",n_valence_basis,n_valence_basis,n_core_states,&
         (-1.d0,0.d0),work,n_basis,ovlp_work(1,n_core_states+1),n_basis,&
         (1.d0,0.d0),ham_work(n_core_states+1,n_core_states+1),n_basis)

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_ovlp_scalapack
  !  NAME
  !   permute_ovlp_scalapack
  subroutine permute_ovlp_scalapack()

    implicit none

    integer :: i

    if(real_eigenvectors) then
       eigenvec(:,:,1) = ovlp(:,:)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pdcopy(n_basis,ovlp,1,perm(i),sc_desc,1,eigenvec(:,:,1),1,i,&
                  sc_desc,1)
          end if
       end do

       ovlp(:,:) = eigenvec(:,:,1)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pdcopy(n_basis,eigenvec(:,:,1),perm(i),1,sc_desc,n_basis,&
                  ovlp,i,1,sc_desc,n_basis)
          end if
       end do

       call pdtran(n_valence_basis,n_valence_basis,1.d0,ovlp,n_core_states+1,&
            n_core_states+1,sc_desc,0.d0,ovlp_real_vv,1,1,sc_desc_vv)

       eigenvec(:,:,n_spin) = ovlp(:,:)

       ! S_vv = S_vv - S_vc * S_cv
       ! pdsyrk should be faster, but this is not critical anyway
       call pdgemm("N","N",n_valence_basis,n_valence_basis,n_core_states,-1.d0,&
            eigenvec(:,:,n_spin),n_core_states+1,1,sc_desc,ovlp,1,&
            n_core_states+1,sc_desc,1.d0,ovlp_real_vv,1,1,sc_desc_vv)
    else
       eigenvec_complex(:,:,1) = ovlp_complex(:,:)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pzcopy(n_basis,ovlp_complex,1,perm(i),sc_desc,1,&
                  eigenvec_complex(:,:,1),1,i,sc_desc,1)
          end if
       end do

       ovlp_complex(:,:) = eigenvec_complex(:,:,1)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pzcopy(n_basis,eigenvec_complex(:,:,1),perm(i),1,sc_desc,&
                  n_basis,ovlp_complex,i,1,sc_desc,n_basis)
          end if
       end do

       call pztranc(n_valence_basis,n_valence_basis,(1.d0,0.d0),ovlp_complex,&
            n_core_states+1,n_core_states+1,sc_desc,(0.d0,0.d0),ovlp_cmplx_vv,&
            1,1,sc_desc_vv)

       eigenvec_complex(:,:,1) = ovlp_complex(:,:)

       ! S_vv = S_vv - S_vc * S_cv
       ! pzherk should be faster, but this is not critical anyway
       call pzgemm("N","N",n_valence_basis,n_valence_basis,n_core_states,&
            (-1.d0,0.d0),eigenvec_complex(:,:,1),n_core_states+1,1,sc_desc,&
            ovlp_complex,1,n_core_states+1,sc_desc,(1.d0,0.d0),ovlp_cmplx_vv,1,&
            1,sc_desc_vv)
    end if

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_ovlp_lapack_real
  !  NAME
  !   permute_ovlp_lapack_real
  subroutine permute_ovlp_lapack_real(ovlp_work,work)

    implicit none

    real*8, intent(inout) :: ovlp_work(n_basis,n_basis)
    real*8, intent(inout) :: work(n_basis,n_basis)

    integer :: i

    work = ovlp_work

    do i = 1,n_basis
       if(i /= perm(i)) then
          work(:,i) = ovlp_work(:,perm(i))
       end if
    end do

    ovlp_work = work

    do i = 1,n_basis
       if(i /= perm(i)) then
          ovlp_work(i,:) = work(perm(i),:)
       end if
    end do

    work = ovlp_work

    ! S_vv = S_vv - S_vc * S_cv
    call dsyrk("U","N",n_valence_basis,n_core_states,-1.d0,&
         work(n_core_states+1,1),n_basis,1.d0,&
         ovlp_work(n_core_states+1,n_core_states+1),n_basis)

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_ovlp_lapack_cmplx
  !  NAME
  !   permute_ovlp_lapack_cmplx
  subroutine permute_ovlp_lapack_cmplx(ovlp_work,work)

    implicit none

    complex*16, intent(inout) :: ovlp_work(n_basis,n_basis)
    complex*16, intent(inout) :: work(n_basis,n_basis)

    integer :: i

    work = ovlp_work

    do i = 1,n_basis
       if(i /= perm(i)) then
          work(:,i) = ovlp_work(:,perm(i))
       end if
    end do

    ovlp_work = work

    do i = 1,n_basis
       if(i /= perm(i)) then
          ovlp_work(i,:) = work(perm(i),:)
       end if
    end do

    work = ovlp_work

    ! S_vv = S_vv - S_vc * S_cv
    call zherk("U","N",n_valence_basis,n_core_states,(-1.d0,0.d0),&
         work(n_core_states+1,1),n_basis,(1.d0,0.d0),&
         ovlp_work(n_core_states+1,n_core_states+1),n_basis)

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_evec_scalapack
  !  NAME
  !   permute_evec_scalapack
  subroutine permute_evec_scalapack(i_spin)

    implicit none

    integer, intent(in) :: i_spin

    integer :: i

    if(real_eigenvectors) then
       eigenvec(:,:,i_spin) = 0.d0

       do i = 1,n_core_states
          if(l_col(i) > 0 .and. l_row(i) > 0) then
             eigenvec(l_row(i),l_col(i),i_spin) = 1.d0&
                /sqrt(ovlp(l_row(i),l_col(i)))
          end if
       end do

       call pdgemr2d(n_valence_basis,n_valence_basis,evec_real_vv,1,1,&
            sc_desc_vv,eigenvec(:,:,i_spin),n_core_states+1,n_core_states+1,&
            sc_desc,my_blacs_ctxt)

       ! C_cv = -S_cv * C_vv
       call pdgemm("N","N",n_core_states,n_valence_basis,n_valence_basis,-1.d0,&
            ovlp,1,n_core_states+1,sc_desc,evec_real_vv,1,1,sc_desc_vv,0.d0,&
            eigenvec(:,:,i_spin),1,n_core_states+1,sc_desc)

       ham(:,:,i_spin) = eigenvec(:,:,i_spin)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pdcopy(n_basis,ham(:,:,i_spin),i,1,sc_desc,n_basis,&
                  eigenvec(:,:,i_spin),perm(i),1,sc_desc,n_basis)
          end if
       end do
    else
       eigenvec_complex(:,:,i_spin) = (0.d0,0.d0)

       do i = 1,n_core_states
          if(l_col(i) > 0 .and. l_row(i) > 0) then
             eigenvec_complex(l_row(i),l_col(i),i_spin) = 1.d0&
                /sqrt(real(ovlp_complex(l_row(i),l_col(i)),kind=8))
          end if
       end do

       call pzgemr2d(n_valence_basis,n_valence_basis,evec_cmplx_vv,1,1,&
            sc_desc_vv,eigenvec_complex(:,:,i_spin),n_core_states+1,&
            n_core_states+1,sc_desc,my_blacs_ctxt)

       ! C_cv = -S_cv * C_vv
       call pzgemm("N","N",n_core_states,n_valence_basis,n_valence_basis,&
            (-1.d0,0.d0),ovlp_complex,1,n_core_states+1,sc_desc,evec_cmplx_vv,&
            1,1,sc_desc_vv,(0.d0,0.d0),eigenvec_complex(:,:,i_spin),1,&
            n_core_states+1,sc_desc)

       ham_complex(:,:,i_spin) = eigenvec_complex(:,:,i_spin)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pzcopy(n_basis,ham_complex(:,:,i_spin),i,1,sc_desc,n_basis,&
                  eigenvec_complex(:,:,i_spin),perm(i),1,sc_desc,n_basis)
          end if
       end do
    end if

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_evec_lapack_real
  !  NAME
  !   permute_evec_lapack_real
  subroutine permute_evec_lapack_real(ovlp_work,evec_work)

    implicit none

    real*8, intent(inout) :: ovlp_work(n_basis,n_basis)
    real*8, intent(inout) :: evec_work(n_basis,n_basis)

    integer :: i

    evec_work(:,1:n_core_states) = 0.d0

    do i = 1,n_core_states
       evec_work(i,i) = 1.d0/sqrt(ovlp_work(i,i))
    end do

    ! C_cv = -S_cv * C_vv
    call dgemm("N","N",n_core_states,n_valence_basis,n_valence_basis,-1.d0,&
         ovlp_work(1,n_core_states+1),n_basis,evec_work(n_core_states+1,&
         n_core_states+1),n_basis,0.d0,evec_work(1,n_core_states+1),n_basis)

    ovlp_work = evec_work

    do i = 1,n_basis
       if(i /= perm(i)) then
          evec_work(perm(i),:) = ovlp_work(i,:)
       end if
    end do

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/permute_evec_lapack_cmplx
  !  NAME
  !   permute_evec_lapack_cmplx
  subroutine permute_evec_lapack_cmplx(ovlp_work,evec_work)

    implicit none

    complex*16, intent(inout) :: ovlp_work(n_basis,n_basis)
    complex*16, intent(inout) :: evec_work(n_basis,n_basis)

    integer :: i

    evec_work(:,1:n_core_states) = (0.d0,0.d0)

    do i = 1,n_core_states
       evec_work(i,i) = 1.d0/sqrt(real(ovlp_work(i,i),kind=8))
    end do

    ! C_cv = -S_cv * C_vv
    call zgemm("N","N",n_core_states,n_valence_basis,n_valence_basis,&
         (-1.d0,0.d0),ovlp_work(1,n_core_states+1),n_basis,&
         evec_work(n_core_states+1,n_core_states+1),n_basis,(0.d0,0.d0),&
         evec_work(1,n_core_states+1),n_basis)

    ovlp_work = evec_work

    do i = 1,n_basis
       if(i /= perm(i)) then
          evec_work(perm(i),:) = ovlp_work(i,:)
       end if
    end do

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/unpermute_ovlp_scalapack
  !  NAME
  !   unpermute_ovlp_scalapack
  subroutine unpermute_ovlp_scalapack()

    implicit none

    integer :: i

    if(real_eigenvectors) then
       ham(:,:,1) = ovlp(:,:)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pdcopy(n_basis,ovlp,1,i,sc_desc,1,ham(:,:,1),1,perm(i),&
                  sc_desc,1)
          end if
       end do

       ovlp(:,:) = ham(:,:,1)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pdcopy(n_basis,ham(:,:,1),i,1,sc_desc,n_basis,ovlp,perm(i),1,&
                  sc_desc,n_basis)
          end if
       end do
    else
       ham_complex(:,:,1) = ovlp_complex(:,:)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pzcopy(n_basis,ovlp_complex,1,i,sc_desc,1,ham_complex(:,:,1),&
                  1,perm(i),sc_desc,1)
          end if
       end do

       ovlp_complex(:,:) = ham_complex(:,:,1)

       do i = 1,n_basis
          if(i /= perm(i)) then
             call pzcopy(n_basis,ham_complex(:,:,1),i,1,sc_desc,n_basis,&
                  ovlp_complex,perm(i),1,sc_desc,n_basis)
          end if
       end do
    end if

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/fetch_core_eval_scalapack
  !  NAME
  !   fetch_core_eval
  subroutine fetch_core_eval_scalapack(i_spin,core_eval)

    implicit none

    integer, intent(in) :: i_spin
    real*8, intent(out) :: core_eval(n_core_states)

    integer :: i

    core_eval = 0.d0

    if(real_eigenvectors) then
       do i = 1,n_core_states
          if(l_col(i) > 0 .and. l_row(i) > 0) then
             core_eval(i) = ham(l_row(i),l_col(i),i_spin)&
                /ovlp(l_row(i),l_col(i))
          end if
       end do
    else
       do i = 1,n_core_states
          if(l_col(i) > 0 .and. l_row(i) > 0) then
             core_eval(i) = real(ham_complex(l_row(i),l_col(i),i_spin),kind=8)&
                /real(ovlp_complex(l_row(i),l_col(i)),kind=8)
          end if
       end do
    end if

    call sync_vector(core_eval,n_core_states,my_scalapack_comm_work)
    call heapsort(core_eval,n_core_states)

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/fetch_core_eval_lapack_real
  !  NAME
  !   fetch_core_eval_lapack_real
  subroutine fetch_core_eval_lapack_real(ham_work,ovlp_work,core_eval)

    implicit none

    real*8, intent(in) :: ham_work(n_basis,n_basis)
    real*8, intent(in) :: ovlp_work(n_basis,n_basis)
    real*8, intent(out) :: core_eval(n_core_states)

    integer :: i

    core_eval = 0.d0

    do i = 1,n_core_states
       core_eval(i) = ham_work(i,i)/ovlp_work(i,i)
    end do

    call heapsort(core_eval,n_core_states)

  end subroutine
  !******
  !--------------------------------------------------------------
  !****s* separate_core_states/fetch_core_eval_lapack_cmplx
  !  NAME
  !   fetch_core_eval_lapack_cmplx
  subroutine fetch_core_eval_lapack_cmplx(ham_work,ovlp_work,core_eval)

    implicit none

    complex*16, intent(in) :: ham_work(n_basis,n_basis)
    complex*16, intent(in) :: ovlp_work(n_basis,n_basis)
    real*8, intent(out) :: core_eval(n_core_states)

    integer :: i

    core_eval = 0.d0

    do i = 1,n_core_states
       core_eval(i) = real(ham_work(i,i),kind=8)/real(ovlp_work(i,i),kind=8)
    end do

    call heapsort(core_eval,n_core_states)

  end subroutine
  !******
  !--------------------------------------------------------------
end module
!******
