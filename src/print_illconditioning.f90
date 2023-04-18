!****s* FHI-aims/print_illconditioning
!  NAME
!    print_illconditioning
!  SYNOPSIS
subroutine print_illconditioning(ev_min_kpt,ev_max_kpt,n_good_basis_kpt)
!  PURPOSE
!    Prints a summary of ill-conditioning for all k-points, including lowest
!    and highest eigenvalues of overlap matrices and number of basis functions
!    to be dropped.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement.
!  HISTORY
!    Created: March 2019.
!  USES
  use dimensions, only: n_periodic,n_k_points,n_basis
  use localorb_io, only: OL_NORM,localorb_info,use_unit
  use mpi_tasks, only: aims_stop_coll
  use runtime_choices, only: basis_threshold,override_illconditioning,&
      use_scalapack
  use scalapack_wrapper, only: my_scalapack_id
  use synchronize_mpi, only: sync_vector
! SOURCE
  implicit none

  real*8, intent(inout) :: ev_min_kpt(n_k_points)
  real*8, intent(inout) :: ev_max_kpt(n_k_points)
  integer, intent(in) :: n_good_basis_kpt(n_k_points)

  integer :: i_kpt
  integer :: n_bad_kpt

  logical :: bad
  logical :: very_bad

  character :: aux1
  character :: aux2
  character(200) :: msg

  bad = .false.
  very_bad = .false.

  if(n_periodic > 0) then
     if(use_scalapack .and. my_scalapack_id /= 0) then
        ev_min_kpt = 0.0d0
        ev_max_kpt = 0.0d0
     end if

     call sync_vector(ev_min_kpt,n_k_points)
     call sync_vector(ev_max_kpt,n_k_points)

     n_bad_kpt = 0

     do i_kpt = 1,n_k_points
        if(n_good_basis_kpt(i_kpt) < n_basis) then
           n_bad_kpt = n_bad_kpt+1
        end if
     end do

     if(n_bad_kpt > 0) then
        bad = .true.

        call localorb_info('',use_unit,'(A)',OL_norm)

        write(aux1,"(I1)") int(log(real(n_bad_kpt,8))/log(10.0d0),4)+1
        write(aux2,"(I1)") int(log(real(n_k_points,8))/log(10.0d0),4)+1

        write(msg,'(A,I'//aux1//',A,I'//aux2//',A)') 'The overlap matrix is'//&
           ' singular on ', n_bad_kpt,' out of ',n_k_points,' k-points'
        call localorb_info(msg,use_unit,'(2X,A)',OL_norm)

        write(msg,'(A,E11.4,A)') 'k-points on which the lowest eigenvalue of'//&
           ' the overlap matrix is below the prescribed threshold ',&
           basis_threshold,' include:'
        call localorb_info(msg,use_unit,'(2X,A)',OL_norm)

        call localorb_info('',use_unit,'(A)',OL_norm)

        write(msg,'(A,4X,A,4X,A,4X,A)') 'k-point','Lowest eigenvalue',&
           'Highest eigenvalue','Reduced number of basis functions'
        call localorb_info(msg,use_unit,'(2X,A)',OL_norm)

        do i_kpt = 1,n_k_points
           if(n_good_basis_kpt(i_kpt) < n_basis) then
              write(msg,'(I7,E21.4,E22.4,I37)') i_kpt,ev_min_kpt(i_kpt),&
                 ev_max_kpt(i_kpt),n_good_basis_kpt(i_kpt)
              call localorb_info(msg,use_unit,'(2X,A)',OL_norm)

              if(ev_min_kpt(i_kpt) < 0.0d0) then
                 very_bad = .true.
              end if
           end if
        end do
     end if
  else ! n_periodic > 0
     if(n_good_basis_kpt(1) < n_basis) then
        bad = .true.

        call localorb_info('',use_unit,'(A)',OL_norm)

        write(msg,'(A)') 'The overlap matrix is singular'
        call localorb_info(msg,use_unit,'(2X,A)',OL_norm)

        write(msg,'(A,E11.4,A,E11.4)') 'The lowest eigenvalue of the'//&
           ' overlap matrix ',ev_min_kpt(1),' is below the prescribed'//&
           ' threshold ',basis_threshold
        call localorb_info(msg,use_unit,'(2X,A)',OL_norm)
     end if

     if(ev_min_kpt(1) < 0.0d0) then
        very_bad = .true.
     end if
  end if ! n_periodic > 0

  if(very_bad) then
     call localorb_info('',use_unit,'(A)',OL_norm)

     call aims_stop_coll('The overlap matrix is not positive definite.')
  else if(bad .and. .not. override_illconditioning) then
     call localorb_info('',use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* Attention: Ill-conditioned overlap matrix detected.'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     call localorb_info('',use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* This may mean that you are running with a very large basis set,'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* that you are running an excessively large cutoff radius in a'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* periodic calculaton, or that you are using diffuse Gaussian basis'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* functions with an insufficient integration grid try to increase'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* "radial_multiplier" for all species in that case). You should not'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* usually encounter this behaviour in normal production calculations.'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* In any case, this behavior can often simply be overcome by setting'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* a large enough value of the "basis_threshold" parameter in'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* "control.in". However, even then you MUST NOW ALSO SET the flag'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* "override_illconditioning .true." explicitly - else FHI-aims will'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* stop (as in this case). Be ABSOLUTELY SURE to only ever set the'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* "override_illconditioning" flag if you understand what you are'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* doing, and always check convergence thoroughly. otherwise, running'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* with a near-singular basis set may lead to absolutely'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* uncontrollable, completely wrong numerical results. When using'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* "basis_threshold" in an intelligent way, there is usually no'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* problem - but by forcing you to set the "override_illconditioning"'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* flag first, we are trying to make sure that you are aware of the'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     write(msg,'(A)') '* problem. Apologies for any inconvenience.'
     call localorb_info(msg,use_unit,'(X,A)',OL_norm)

     call localorb_info('',use_unit,'(X,A)',OL_norm)

     call aims_stop_coll('Ill-conditioned overlap matrix found.')
  end if

end subroutine print_illconditioning
!******
