!****s* FHI-aims/orthonormalize_prodbas_fn
!  NAME
!   orthonormalize_prodbas_fn
!  SYNOPSIS
      subroutine orthonormalize_prodbas_fn &
      ( prodbas_metric, i_l, i_function, i_first_fn, n_max_basbas_fns,basis_wave, n_grid, r_grid, &
        accuracy, success &
      )

!  PURPOSE
!  Orthonormalization specific to the basis_wave array. It modifies basis_wave(i_function),
!  basis_kinetic(i_function), and also (if needed) basis_deriv(i_function)
!
!  USES

      use dimensions
      use mpi_tasks, only: STDERR
      implicit none

!  ARGUMENTS
      integer prodbas_metric
      integer i_l
      integer i_function
      integer i_first_fn
      integer n_max_basbas_fns
      real*8  basis_wave (n_max_grid, n_max_basbas_fns)
      integer n_grid
      real*8  r_grid (n_grid)
      real*8  accuracy
      integer success

!  INPUTS
!  o prodbas_metric -- 0, Normal metric, and 1 Coulomb metric
!  o             using delta metric, i.e., the overlap. if false, the
!  o             basis are orthonormalized using the coulomb metric. 
!  o i_l -- current angular momentum channel
!  o i_function -- the current auxiliary radial function that is under
!                   consideration
!  o i_first_fn -- the first auxiliary radial function for current
!               angular momentum channel
!  o basis_wave -- the value auxiliary radial function
!  o n_grid -- number of logarithmic grid points
!  o r_grid -- logrithmic grid
!  o accuracy -- the auxiliary basis accuracy for on-site orthonormalization
!  o success -- logical, if true, tells whether the current radial function 
!              will be accepted.
!  OUTPUTS
!  orthonormalized basis functions
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
!  Subroutine orthonormalize_prod_basis_fn orthonormalizes a present basis function against
!  the past ones of same angular momentum.

!      real*8 basis_kinetic (n_max_grid, n_max_basis_fns)
!      real*8 basis_deriv (n_max_grid, n_max_basis_fns)

!  imported variables

!  input

!  local variables

      real*8 norm
      real*8 scalar_prod(n_max_basbas_fns)
      character l_shell_str
      real*8 ortho

!  counters

      integer i_prev, i_grid, i_pass

!  functions

      real*8 int_log_mesh
      real*8 int_log_coulomb_metric

!  begin work

!      check orthogonality

              ortho = 0.
              norm = 1.
              do i_prev = i_first_fn, i_function-1, 1
               if(prodbas_metric.eq.0) then
                scalar_prod(i_prev) = &
                  int_log_mesh ( &
                    basis_wave(1,i_function), basis_wave(1,i_prev), &
                    n_grid, &
                    r_grid &
                   )
               elseif(prodbas_metric.eq.1) then
                scalar_prod(i_prev) = &
                  int_log_coulomb_metric ( i_l, &
                    basis_wave(1,i_function), basis_wave(1,i_prev), &
                    n_grid, &
                    r_grid &
                  )
               endif
                ortho = ortho + scalar_prod(i_prev)*scalar_prod(i_prev)
!test
!                if (i_function.eq.10) then
!                  write(use_unit,'(A,2I4,2f16.10)') "scalar_prod**2",i_prev,i_function, &
!                    scalar_prod(i_prev)*scalar_prod(i_prev), ortho
!                end if
!test end
              enddo

              norm = 1. - ortho
!test
!              write(use_unit,*) " Function: ", i_function, &
!               ", norm:", norm
!test

!              l_shell_str = l_to_str(i_l)

              if (norm.gt.accuracy) then

                i_pass = 0

!               if orthogonality is not fulfilled, then re-orthogonalize
                do while (ortho.gt.1.d-20)

!                  write(use_unit,*) "i_pass, norm, ortho", i_pass, norm, ortho

!                 orthonormalize basis_wave; apply the same operations
!                 to kinetic energy term
                  do i_grid = 1, n_grid,1

                    do i_prev = i_first_fn, i_function-1, 1
                      basis_wave (i_grid, i_function) = &
                      basis_wave(i_grid, i_function) - &
                      scalar_prod(i_prev)* basis_wave(i_grid, i_prev)

!                      basis_kinetic(i_grid, i_function) =
!     +                basis_kinetic(i_grid,i_function) -
!     +                scalar_prod(i_prev) * basis_kinetic(i_grid,i_prev)

                    enddo
!                    write(use_unit,*)i_grid, basis_wave(i_grid,i_function)

                    basis_wave (i_grid, i_function) = &
                    basis_wave (i_grid, i_function)/sqrt(norm)

!                    basis_kinetic (i_grid, i_function) =
!     +              basis_kinetic (i_grid, i_function)/sqrt(norm)

                  enddo

!                 if required, apply exactly the same operations to the stored radial derivative
!                 of the basis function
!                 notice that use_basis_gradients is conveniently imported from module dimensions
!                  if (use_basis_gradients) then
!                    do i_grid = 1, n_grid,1
!
!                      do i_prev = i_first_fn, i_function-1, 1
!                        basis_deriv (i_grid, i_function) =
!     +                  basis_deriv(i_grid, i_function) -
!     +                  scalar_prod(i_prev)* basis_deriv(i_grid, i_prev)
!                      enddo
!                      basis_deriv (i_grid, i_function) =
!     +                basis_deriv (i_grid, i_function)/sqrt(norm)
!
!                    enddo
!                  end if

!                 due to numerical inaccuracies, basis_wave may no longer be properly normalized
!                 So, normalize what we have left.

                  
                  if(prodbas_metric.eq.0) then
                     norm = &
                       int_log_mesh ( &
                      basis_wave(1,i_function), basis_wave(1,i_function), &
                      n_grid, &
                      r_grid )
                  elseif(prodbas_metric.eq.1) then
                     norm = &
                       int_log_coulomb_metric ( i_l,&
                      basis_wave(1,i_function), basis_wave(1,i_function), &
                      n_grid, &
                      r_grid )
                  endif  

                  do i_grid = 1, n_grid, 1

                    basis_wave (i_grid, i_function) = &
                    basis_wave (i_grid, i_function)/sqrt(norm)

!                    basis_kinetic (i_grid, i_function) =
!     +              basis_kinetic (i_grid, i_function)/sqrt(norm)

                  enddo

!                  if (use_basis_gradients) then
!                    do i_grid = 1, n_grid, 1
!                      basis_deriv (i_grid, i_function) =
!     +                basis_deriv (i_grid, i_function)/sqrt(norm)
!                    enddo
!                  end if

!                 check if the completely ON function is still orthonormal to all others after first pass:

                  ortho = 0.
                  do i_prev = i_first_fn, i_function-1, 1
                   if(prodbas_metric.eq.0) then
                     scalar_prod(i_prev) = &
                     int_log_mesh ( &
                       basis_wave(1,i_function), basis_wave(1,i_prev), &
                       n_grid, &
                       r_grid )
                   elseif(prodbas_metric.eq.1) then
                     scalar_prod(i_prev) = &
                       int_log_coulomb_metric ( i_l,&
                      basis_wave(1,i_function), basis_wave(1,i_prev), &
                      n_grid, &
                      r_grid )
                   endif 
                      ortho = &
                      ortho+scalar_prod(i_prev)*scalar_prod(i_prev)
                  enddo

!test
!                  write(use_unit,*) "  Orthogonality: ", ortho
!test

                  norm = 1. - ortho


                  i_pass = i_pass + 1

!                  write(use_unit,*) "i_pass, norm, ortho 2", i_pass, norm, ortho

                  if (i_pass.gt.3) then
                    write(STDERR,'(1X,A,A,I5,A)') &
                      "* More than three passes needed to ", &
                      "orthonormalize basis fn ", i_function,"."
                    write(STDERR,'(1X,A,A)') &
                      "* Aborting - perhaps choose smaller basis size", &
                      " or accuracy!"
!                    stop
                    exit
                  end if

                enddo

                if (i_pass.gt.3) then
                  success = -1
                else
                  success = 1
                endif

              else

                success = 0

              end if

!  thats all folks

      return

    end subroutine orthonormalize_prodbas_fn
   !******
