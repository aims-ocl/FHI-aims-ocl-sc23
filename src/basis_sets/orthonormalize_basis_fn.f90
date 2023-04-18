!****s* FHI-aims/orthonormalize_basis_fn
!  NAME
!   orthonormalize_basis_fn
!  SYNOPSIS
 subroutine orthonormalize_basis_fn &
  ( i_function, i_first_fn, basis_wave, n_grid, r_grid, &
  & accuracy, basis_kinetic, basis_deriv, success, &
  & do_wave_at_nucleus, wave_at_nucleus)

!  PURPOSE
!  Subroutine orthonormalize_basis_fn orthonormalizes a present basis function against
!  the past ones of same angular momentum.
!  It is specific to the basis_wave array. It modifies basis_wave(i_function),
!  basis_kinetic(i_function), and also (if needed) basis_deriv(i_function)
!
!  USES

  use dimensions,  only : n_max_grid, n_max_basis_fns, use_basis_gradients
  use runtime_choices, only : flag_rel, REL_x2c, REL_4c_dks
  use mpi_tasks, only : STDERR
  implicit none

!  ARGUMENTS
  integer i_function
  integer i_first_fn
  real*8 basis_wave (n_max_grid, n_max_basis_fns)
  integer n_grid
  real*8 r_grid (n_grid)
  real*8 accuracy
  real*8 basis_kinetic (n_max_grid, n_max_basis_fns)
  real*8 basis_deriv (n_max_grid, n_max_basis_fns)
  integer success
  ! Whether to perform the orthonormalization operation on wave_at_nucleus
  logical, intent(in) :: do_wave_at_nucleus
  ! Wavefunction values at the nucleus
  real*8, intent(in out) :: wave_at_nucleus(n_max_basis_fns)
!
!  INPUTS
!  o i_function
!  o i_first_fn
!  o basis_wave
!  o n_grid
!  o r_grid
!  o accuracy
!  o basis_kinetic
!  o basis_deriv
!  o success
!  o do_wave_at_nucleus
!  o wave_at_nucleus
!  OUTPUTS
!  
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals:
!    FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject
!   to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!  local variables

  real*8 norm
  real*8 scalar_prod(n_max_basis_fns)
  character l_shell_str
  real*8 ortho

!  counters

  integer i_prev, i_grid, i_pass

!  functions

  real*8 int_log_mesh
      
! (Rundong) For relativistic cases, currently I have no idea about how to Gram-Schmidt 
! orthogonalize the large and small components while retaining their RKB condtion. 
! I suppose that a Lowdin orthogonalization before the SCF should be enough.
  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
    success=i_function
    return 
  endif

!  begin work

!      check orthogonality

  ortho = 0.
  norm = 1.
  do i_prev = i_first_fn, i_function-1, 1

    scalar_prod(i_prev) = int_log_mesh ( &
      basis_wave(1,i_function), basis_wave(1,i_prev), n_grid, r_grid )

    ortho = ortho + scalar_prod(i_prev) * scalar_prod(i_prev)

  enddo

  norm = 1. - ortho

  if (norm.gt.accuracy) then

    i_pass = 0

!   if orthogonality is not fulfilled, then re-orthogonalize
    do while (ortho.gt.1.d-08)

!       orthonormalize basis_wave; apply the same operations
!       to kinetic energy term
        do i_grid = 1, n_grid,1

          do i_prev = i_first_fn, i_function-1, 1
            basis_wave (i_grid, i_function) = &
            basis_wave(i_grid, i_function) - &
            scalar_prod(i_prev) * basis_wave(i_grid, i_prev)

            basis_kinetic(i_grid, i_function) = &
            basis_kinetic(i_grid,i_function) - &
            scalar_prod(i_prev) * basis_kinetic(i_grid,i_prev)

          enddo

          basis_wave (i_grid, i_function) = &
          basis_wave (i_grid, i_function)/sqrt(norm)

          basis_kinetic (i_grid, i_function) = &
          basis_kinetic (i_grid, i_function)/sqrt(norm)

        enddo

        if (do_wave_at_nucleus) then
           do i_prev = i_first_fn, i_function-1
              wave_at_nucleus(i_function) = &
                   & wave_at_nucleus(i_function) - &
                   & scalar_prod(i_prev) * wave_at_nucleus(i_prev)
           end do
           wave_at_nucleus(i_function) = &
                & wave_at_nucleus(i_function)/sqrt(norm)
        end if

!       if required, apply exactly the same operations to the stored radial derivative
!       of the basis function
!       notice that use_basis_gradients is conveniently imported from module dimensions
        if (use_basis_gradients) then
          do i_grid = 1, n_grid,1

            do i_prev = i_first_fn, i_function-1, 1
              basis_deriv (i_grid, i_function) = &
              basis_deriv(i_grid, i_function) - &
              scalar_prod(i_prev) * basis_deriv(i_grid, i_prev)
            enddo
            basis_deriv (i_grid, i_function) = &
            basis_deriv (i_grid, i_function)/sqrt(norm)

          enddo
        end if

!       due to numerical inaccuracies, basis_wave may no longer be properly normalized
!       So, normalize what we have left.

        norm = int_log_mesh &
          (basis_wave(1,i_function), basis_wave(1,i_function), n_grid,r_grid )

        do i_grid = 1, n_grid, 1

          basis_wave (i_grid, i_function) = &
          basis_wave (i_grid, i_function)/sqrt(norm)

          basis_kinetic (i_grid, i_function) = &
          basis_kinetic (i_grid, i_function)/sqrt(norm)

        enddo

        if (use_basis_gradients) then
          do i_grid = 1, n_grid, 1
            basis_deriv (i_grid, i_function) = &
            basis_deriv (i_grid, i_function)/sqrt(norm)
          enddo
        end if

!       check if the completely ON function is still orthonormal to all others after first pass:

        ortho = 0.
        do i_prev = i_first_fn, i_function-1, 1

          scalar_prod(i_prev) = int_log_mesh ( &
            basis_wave(1,i_function), basis_wave(1,i_prev), n_grid, r_grid )

          ortho = ortho + scalar_prod(i_prev) * scalar_prod(i_prev)

        enddo

        norm = 1. - ortho

        i_pass = i_pass + 1

        if (i_pass.gt.3) then
          write(STDERR,'(1X,A,A,I5,A)') &
            "* More than three passes needed to ","orthonormalize basis fn ",i_function,"."
          write(STDERR,'(1X,A,A)') &
            "* Aborting - perhaps choose smaller basis size"," or accuracy!"
          stop
        end if

    enddo

    success = i_function

  else

    success = 0

  end if

! thats all folks

  if (do_wave_at_nucleus) wave_at_nucleus(i_function) = &
       & wave_at_nucleus(i_function)/sqrt(norm)
  return

 end subroutine orthonormalize_basis_fn
!******

 subroutine orthonormalize_basis_fn_rel &
  ( i_species, i_function, i_first_fn, i_l, kappa, &
    basis_wave, basis_small, n_grid, r_grid, accuracy, &
    basis_kinetic, basis_deriv, success_last, success )

  use dimensions,  only : n_max_grid, n_species, n_max_basis_fns, use_basis_gradients
  use runtime_choices, only : flag_rel, REL_x2c, REL_4c_dks
  use mpi_tasks, only : STDERR
  implicit none

!  ARGUMENTS
  integer,intent(in) :: i_species
  integer,intent(in) :: i_function   ! current function
  integer,intent(in) :: i_first_fn   ! the starting function of the present angular momentum i_l
  integer,intent(in) :: i_l, kappa
  integer,intent(in) :: n_grid
  real*8,intent(in)  :: r_grid (n_grid)
  real*8,intent(in)  :: accuracy
  integer,intent(in) :: success_last ! the value of the success of the last i_function

  real*8,intent(inout) :: basis_wave (n_max_grid, n_max_basis_fns)
  real*8,intent(inout) :: basis_small (n_max_grid, n_max_basis_fns)
  real*8,intent(inout) :: basis_kinetic (n_max_grid, n_max_basis_fns)
  real*8,intent(inout) :: basis_deriv (n_max_grid, n_max_basis_fns)
  integer,intent(inout) :: success(n_species,n_max_basis_fns)

!  local variables
  real*8 :: norm, ortho, temp
  real*8 :: scalar_prod(n_max_basis_fns)
  character :: l_shell_str
  real*8 :: rkb(n_grid,2)

!  counters
  integer :: i_prev, i_grid, i_pass

!  functions
  real*8 :: int_log_mesh
      
  do i_grid=1, n_grid
     rkb(i_grid,1) = basis_small(i_grid,i_function) / basis_wave(i_grid,i_function)
     if(dabs(basis_wave(i_grid,i_function)).lt.1.d-10 .or. &
        dabs(basis_small(i_grid,i_function)).lt.1.d-10 .or. &
        dabs(rkb(i_grid,1)).gt.1.d3)  rkb(i_grid,1) = rkb(i_grid-1,1)
  enddo

  ortho = 0.d0; norm = 1.d0

  if(i_l.eq.0)then   ! s orbitals

     ! check orthogonality
     do i_prev = i_first_fn, i_function-1, 1

         scalar_prod(i_prev) = &
         int_log_mesh( basis_wave(1,i_function), basis_wave(1,i_prev), n_grid,r_grid ) &
         + int_log_mesh( basis_small(1,i_function), basis_small(1,i_prev), n_grid,r_grid )

         ortho = ortho + scalar_prod(i_prev) * scalar_prod(i_prev)
     enddo

     norm = 1.d0 - ortho
    !write(6,"(50x,'s orbital: i_function=',i3,5x,'ortho=',f18.14)")i_function,ortho

     if (dabs(norm).gt.accuracy) then

       i_pass = 0

       ! if orthogonality is not fulfilled, then re-orthogonalize
       do while (ortho.gt.1.d-08)

           ! orthonormalize basis_wave; apply the same operations
           ! to kinetic energy term
           do i_prev = i_first_fn, i_function-1, 1

              norm = int_log_mesh( basis_wave(1,i_prev), basis_wave(1,i_prev), n_grid,r_grid ) + & 
                     int_log_mesh( basis_small(1,i_prev), basis_small(1,i_prev), n_grid,r_grid )
              temp = scalar_prod(i_prev)/norm

              do i_grid = 1, n_grid
                 basis_wave (i_grid, i_function) = &
                 basis_wave (i_grid, i_function) - temp * basis_wave(i_grid, i_prev)

                 basis_small (i_grid, i_function) = &
                 basis_small (i_grid, i_function) - temp * basis_small(i_grid, i_prev)

                 basis_kinetic(i_grid, i_function) = &
                 basis_kinetic(i_grid, i_function) - temp * basis_kinetic(i_grid, i_prev)
              enddo

           enddo

           ! if required, apply exactly the same operations to the stored radial derivative
           ! of the basis function
           ! notice that use_basis_gradients is conveniently imported from module dimensions
           if (use_basis_gradients) then
           do i_prev = i_first_fn, i_function-1, 1
              norm = int_log_mesh( basis_wave(1,i_prev), basis_wave(1,i_prev), n_grid,r_grid ) + & 
                     int_log_mesh( basis_small(1,i_prev), basis_small(1,i_prev), n_grid,r_grid )
              temp = scalar_prod(i_prev)/norm
              do i_grid = 1, n_grid
                 basis_deriv (i_grid, i_function) = &
                 basis_deriv (i_grid, i_function) - temp * basis_deriv(i_grid, i_prev)
              enddo
           enddo
           end if

           ! due to numerical inaccuracies, basis_wave may no longer be properly normalized
           ! So, normalize what we have left.

           norm = int_log_mesh( basis_wave(1,i_function), basis_wave(1,i_function), n_grid,r_grid ) &
             + int_log_mesh( basis_small(1,i_function), basis_small(1,i_function), n_grid,r_grid )

           do i_grid = 1, n_grid, 1

             basis_wave (i_grid, i_function) = basis_wave (i_grid, i_function)/sqrt(norm)
             basis_small (i_grid, i_function) = basis_small(i_grid, i_function)/sqrt(norm)
             basis_kinetic (i_grid, i_function) = basis_kinetic (i_grid, i_function)/sqrt(norm)

           enddo

           if (use_basis_gradients) then
             do i_grid = 1, n_grid, 1
               basis_deriv (i_grid, i_function) = basis_deriv (i_grid, i_function)/sqrt(norm)
             enddo
           end if

           ! check if the completely ON function is still orthonormal to all others after first pass:

           ortho = 0.d0
           do i_prev = i_first_fn, i_function-1, 1

             scalar_prod(i_prev) = &
               int_log_mesh( basis_wave(1,i_function), basis_wave(1,i_prev), n_grid, r_grid ) &
               + int_log_mesh( basis_small(1,i_function), basis_small(1,i_prev), n_grid, r_grid )

             ortho = ortho + scalar_prod(i_prev) * scalar_prod(i_prev)

           enddo

           norm = 1.d0 - ortho
           i_pass = i_pass + 1
           write(6,"(80x,'ortho=',f18.14)")ortho

           if (i_pass.gt.3) then
             write(STDERR,'(1X,A,A,I5,A)') "* More than three passes needed to ","orthonormalize basis fn ",i_function,"."
             write(STDERR,'(1X,A,A)') "* Aborting - perhaps choose smaller basis size"," or accuracy!"
             stop
           end if

       enddo

       success(i_species,i_function) = i_function

     else ! if norm.le.accuracy:

       success(i_species,i_function) = 0

     end if

  elseif (kappa.eq.i_l) then ! the spin down branch of p, d, f, g, h orbitals

     do i_prev = i_first_fn, i_function-2, 2

         scalar_prod(i_prev) = &
         int_log_mesh( basis_wave(1,i_function), basis_wave(1,i_prev), n_grid,r_grid ) &
         + int_log_mesh( basis_small(1,i_function), basis_small(1,i_prev), n_grid,r_grid ) &
         + int_log_mesh( basis_wave(1,i_function+1), basis_wave(1,i_prev+1), n_grid,r_grid ) &
         + int_log_mesh( basis_small(1,i_function+1), basis_small(1,i_prev+1), n_grid,r_grid )

         ortho = ortho + scalar_prod(i_prev) * scalar_prod(i_prev)
     enddo

     norm = 1.d0 - ortho
    !write(6,"(50x,'p/d/f/g orbital: i_function=',i3,5x,'ortho=',f18.14)")i_function,ortho

     if (dabs(norm).gt.accuracy) then

       i_pass = 0

       ! if orthogonality is not fulfilled, then re-orthogonalize
       do while (ortho.gt.1.d-08)

           ! orthonormalize basis_wave; apply the same operations to kinetic energy term
           do i_prev = i_first_fn, i_function-2, 2

              norm = int_log_mesh( basis_wave(1,i_prev), basis_wave(1,i_prev), n_grid,r_grid ) + & 
                     int_log_mesh( basis_small(1,i_prev), basis_small(1,i_prev), n_grid,r_grid ) + &
                     int_log_mesh( basis_wave(1,i_prev+1), basis_wave(1,i_prev+1), n_grid,r_grid ) + & 
                     int_log_mesh( basis_small(1,i_prev+1), basis_small(1,i_prev+1), n_grid,r_grid )
              temp = scalar_prod(i_prev)/norm

              do i_grid = 1, n_grid
                 basis_wave (i_grid, i_function) = &
                 basis_wave (i_grid, i_function) - temp * basis_wave(i_grid, i_prev)

                 basis_small (i_grid, i_function) = &
                 basis_small (i_grid, i_function) - temp * basis_small(i_grid, i_prev)

                 basis_kinetic(i_grid, i_function) = &
                 basis_kinetic(i_grid, i_function) - temp * basis_kinetic(i_grid, i_prev)

                 basis_wave (i_grid, i_function+1) = &
                 basis_wave (i_grid, i_function+1) - temp * basis_wave(i_grid, i_prev+1)

                 basis_small (i_grid, i_function+1) = &
                 basis_small (i_grid, i_function+1) - temp * basis_small(i_grid, i_prev+1)

                 basis_kinetic(i_grid, i_function+1) = &
                 basis_kinetic(i_grid, i_function+1) - temp * basis_kinetic(i_grid, i_prev+1)
              enddo
           enddo

           ! if required, apply exactly the same operations to the stored radial derivative
           ! of the basis function
           ! notice that use_basis_gradients is conveniently imported from module dimensions
           if (use_basis_gradients) then
           do i_prev = i_first_fn, i_function-2, 2
              norm = int_log_mesh( basis_wave(1,i_prev), basis_wave(1,i_prev), n_grid,r_grid ) + & 
                     int_log_mesh( basis_small(1,i_prev), basis_small(1,i_prev), n_grid,r_grid ) + &
                     int_log_mesh( basis_wave(1,i_prev+1), basis_wave(1,i_prev+1), n_grid,r_grid ) + & 
                     int_log_mesh( basis_small(1,i_prev+1), basis_small(1,i_prev+1), n_grid,r_grid )
              temp = scalar_prod(i_prev)/norm
              do i_grid = 1, n_grid
                 basis_deriv(i_grid, i_function) = &
                 basis_deriv(i_grid, i_function) - temp * basis_deriv(i_grid, i_prev)

                 basis_deriv(i_grid, i_function+1) = &
                 basis_deriv(i_grid, i_function+1) - temp * basis_deriv(i_grid, i_prev+1)
              enddo
           enddo
           end if

           ! due to numerical inaccuracies, basis_wave may no longer be properly normalized
           ! So, normalize what we have left.

           norm = int_log_mesh( basis_wave(1,i_function), basis_wave(1,i_function), n_grid,r_grid ) &
             + int_log_mesh( basis_small(1,i_function), basis_small(1,i_function), n_grid,r_grid ) &
             + int_log_mesh( basis_wave(1,i_function+1), basis_wave(1,i_function+1), n_grid,r_grid ) &
             + int_log_mesh( basis_small(1,i_function+1), basis_small(1,i_function+1), n_grid,r_grid )

           do i_grid = 1, n_grid, 1
              basis_wave(   i_grid, i_function) = basis_wave(   i_grid, i_function)/sqrt(norm)
              basis_small(  i_grid, i_function) = basis_small(  i_grid, i_function)/sqrt(norm)
              basis_kinetic(i_grid, i_function) = basis_kinetic(i_grid, i_function)/sqrt(norm)

              basis_wave(   i_grid, i_function+1) = basis_wave(   i_grid, i_function+1)/sqrt(norm)
              basis_small(  i_grid, i_function+1) = basis_small(  i_grid, i_function+1)/sqrt(norm)
              basis_kinetic(i_grid, i_function+1) = basis_kinetic(i_grid, i_function+1)/sqrt(norm)
           enddo

           if (use_basis_gradients) then
           do i_grid = 1, n_grid, 1
              basis_deriv (i_grid, i_function)   = basis_deriv (i_grid, i_function)/sqrt(norm)
              basis_deriv (i_grid, i_function+1) = basis_deriv (i_grid, i_function+1)/sqrt(norm)
           enddo
           end if

           ! check if the completely ON function is still orthonormal to all others after first pass:

           ortho = 0.d0
           do i_prev = i_first_fn, i_function-2, 2

             scalar_prod(i_prev) = &
               int_log_mesh( basis_wave(1,i_function), basis_wave(1,i_prev), n_grid, r_grid ) &
               + int_log_mesh( basis_small(1,i_function), basis_small(1,i_prev), n_grid, r_grid ) &
               + int_log_mesh( basis_wave(1,i_function+1), basis_wave(1,i_prev+1), n_grid, r_grid ) &
               + int_log_mesh( basis_small(1,i_function+1), basis_small(1,i_prev+1), n_grid, r_grid )

             ortho = ortho + scalar_prod(i_prev) * scalar_prod(i_prev)

           enddo

           norm = 1.d0 - ortho
           i_pass = i_pass + 1
           write(6,"(80x,'ortho=',f18.14)")ortho

           if (i_pass.gt.3) then
             write(STDERR,'(1X,A,A,I5,A)') "* More than three passes needed to ","orthonormalize basis fn ",i_function,"."
             write(STDERR,'(1X,A,A)') "* Aborting - perhaps choose smaller basis size"," or accuracy!"
             stop
           end if

       enddo

       success(i_species,i_function) = i_function

     else ! if norm.le.accuracy:

       success(i_species,i_function) = 0

     end if

  elseif (kappa.eq.-i_l-1) then
     ! for kappa.eq.-l-1 (viz. the spin up branch), we skip it, because the wavefunctions
     ! have already be processed in the previous spin down branch
     if(success_last.eq.0 )then
        success(i_species,i_function) = 0
     else
        success(i_species,i_function) = i_function
     endif
  endif

 end subroutine orthonormalize_basis_fn_rel



