
!****h* FHI-aims/bfgs
!  NAME
!    bfgs
!  SYNOPSIS

      module bfgs

!  PURPOSE
!  ????????????????
!
!
!  USES

      use dimensions
      use localorb_io
      implicit none

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!
!******







      contains

!****s* bfgs/bfgs_constraint_v2
!  NAME
!    bfgs_constraint_v2
!  SYNOPSIS
      subroutine bfgs_constraint_v2(f_external, n_arguments, arguments, &
           max_fn_evals, eps, change_in_f, arg_mag, delta0, &
           n_fn_evals, hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )
!  PURPOSE
!    Find a solution to the equations f(x) = 0 using the BFGS-method
!    and finite-differencing for computing the derivatives.
!
!  ARGUMENTS
      integer :: n_arguments
      real*8 :: arguments(n_arguments)
      integer :: max_fn_evals
      real*8 :: eps
      real*8 :: change_in_f
      real*8 :: delta0
      external :: f_external
      real*8 hamiltonian( n_basis*(n_basis+1)/2, n_spin )
      real*8 overlap_matrix( n_basis*(n_basis+1)/2 )
      real*8, dimension(n_states, n_spin) :: KS_eigenvalue
      real*8, dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
      real*8, dimension(n_states, n_spin) :: occ_numbers
      real*8 :: n_electrons
      real*8 :: chemical_potential
      real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
      real*8, dimension(n_region, n_spin) :: electrons_in_region
      integer :: current_spin
      real*8 :: f_value
      integer :: n_fn_evals
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o n_arguments -- number of arguments for the objective function
!    o arguments -- arguments for the objective function
!    o max_fn_evals -- maximum number of evaluations of the objective function
!    o eps -- tolerance of the objective function, i.e., |f(x)| < eps
!    o change_in_f -- inital value of f, the objective function
!    o delta -- basic differencing parameter for the evaluation of the gradient of f
!    o f_external -- external routine to evaluate the objective function f
!    o hamiltonian -- the Hamilton matrix
!    o overlap_matrix -- the overlap matrix
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o occ_numbers -- occupation numbers
!    o n_electrons -- total number of electrons
!    o chemical_potential -- current chemical potential
!    o constraint_proj -- projector to the constrained region
!    o electrons_in_region -- number of electrons in the constrained regions
!    o current_spin -- current spin channel
!  OUTPUT
!    o f_value -- final value obtained by the objective function f
!    o n_fn_evals -- number of evaluations of f taken
!  SOURCE

!     locals
      real*8 :: gradient(n_arguments)
      real*8 :: gradient_work(n_arguments)

      real*8 :: hessian_d(n_arguments)
      real*8 :: hessian_d_work(n_arguments)

      real*8 :: hessian_l(n_arguments,n_arguments)
      real*8 :: hessian_l_work(n_arguments,n_arguments)

      real*8 :: s(n_arguments)
      real*8 :: y(n_arguments)

      real*8 :: sHs
      real*8 :: try_arguments(n_arguments)
      real*8 :: arg_mag(n_arguments)

      logical :: converged
      logical :: scaled_update
      logical :: global_improvement
      logical :: longer_step
      logical :: shorter_step

!      integer :: info
      integer :: method

      real*8 :: f_value_1
      real*8 :: f_value_2
!      real*8 :: f_previous
      real*8 :: alpha
      real*8 :: alpha_scale
!      real*8 :: y_scale
      real*8 :: sigma
      real*8 :: s_times_y
      real*8 :: s_times_g
      real*8 :: total_alpha
      real*8 :: deps
      real*8 :: aeps
 !     real*8 :: temp
      real*8 :: diff_in_f
      real*8 :: norm_s

      real*8 :: delta

!      character*120 :: info_str

!      integer :: cause_of_exit

!     counters
      integer :: i_index

!     minimal values for alpha and delta
      aeps = 0.1d0*eps
      deps = 0.1d0*aeps
!      deps = min(0.1d0*aeps,0.001d0*delta0)

!     initialize hessian to be the identity matrix
      do i_index = 1, n_arguments, 1
         hessian_d(i_index) = 1.0d0
      enddo
      hessian_l = 0.0d0

      method = 1
      diff_in_f = change_in_f
      delta = delta0

      n_fn_evals = 1
      call f_external(n_arguments, arguments, &
           f_value, n_fn_evals, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )

      call evaluate_gradient_constraint(n_arguments, arguments, f_value, &
           n_fn_evals, delta, arg_mag, f_external, gradient, method, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )

      converged = .false.

!     this is the main BFGS loop
      do while ((.not.converged).and.(n_fn_evals.lt.max_fn_evals))

         hessian_d_work = hessian_d
         hessian_l_work = hessian_l
         s = -gradient
         call solve_factored_problem(hessian_d_work, &
              hessian_l_work, n_arguments, s)

         sHs = 0.0d0
         norm_s = 0.0d0
         do i_index = 1, n_arguments, 1
            sHs = sHs + gradient(i_index)*s(i_index)
            norm_s = norm_s + s(i_index)**2
         enddo
         norm_s = sqrt(norm_s)
!         alpha = -2.0d0*diff_in_f/sHs
!         alpha = 0.5d0

         alpha = -diff_in_f/(norm_s*sHs)
         if (alpha.gt.1.0d0) then
            alpha = 1.0d0
         end if

         global_improvement = .false.
         longer_step = .false.
         shorter_step = .false.
         total_alpha = 0.0d0

!         print *,'Spin:', current_spin
!         print *,'Fn_evals:', n_fn_evals
!         print *,'differece_method: ', method
!         print *,'eps: ', eps
!         print *,'aeps: ', aeps
!         print *,'deps: ', deps
!         print *,'sHs: ', sHs
!         print *,'norm_s: ', norm_s
!         print *,'alpha: ', alpha
!         print *,'delta: ', delta
!         print *,'f_value: ', f_value

!     first line search step, try = current + alpha*s
         try_arguments = arguments + alpha*s
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, try_arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         if (f_value_1.lt.f_value) then
            f_value = f_value_1
            arguments = try_arguments
            total_alpha = total_alpha + alpha
            global_improvement = .true.
            longer_step = .true.
         else
            shorter_step = .true.
         end if

!     if first step was successful then set alpha <- 2.0d0*alpha as long as
!     there is improvement
         do while (longer_step)

            alpha = 2.0d0*alpha
            try_arguments = arguments + alpha*s
            n_fn_evals = n_fn_evals + 1
            call f_external(n_arguments, try_arguments, &
                 f_value_1, n_fn_evals, &
                 hamiltonian, &
                 overlap_matrix,KS_eigenvalue, KS_eigenvector, &
                 occ_numbers, &
                 n_electrons, chemical_potential, constraint_proj, &
                 electrons_in_region, &
                 current_spin &
                 )

            if (f_value_1.lt.f_value) then
               f_value = f_value_1
               arguments = try_arguments
               total_alpha = total_alpha + alpha
            else
               longer_step = .false.
            end if

         enddo

!     if there has not been any improvement so far, try setting
!     alpha <- 0.5d0*alpha as long as improvement is found
         do while ((shorter_step).and.(.not.global_improvement))

            alpha = 0.5d0*alpha
            try_arguments = arguments + alpha*s
            n_fn_evals = n_fn_evals + 1
            call f_external(n_arguments, try_arguments, &
                 f_value_2, n_fn_evals, &
                 hamiltonian, &
                 overlap_matrix,KS_eigenvalue, KS_eigenvector, &
                 occ_numbers, &
                 n_electrons, chemical_potential, constraint_proj, &
                 electrons_in_region, &
                 current_spin &
                 )

            if (f_value_2.lt.f_value) then
               f_value = f_value_2
               arguments = try_arguments
               total_alpha = total_alpha + alpha
               global_improvement = .true.

            else

               alpha_scale = 0.1d0
               if (f_value_1 + f_value .gt. f_value_2 + f_value_2) then
                  alpha_scale = 1.0d0 + 0.5d0 * (f_value - f_value_1) / &
                       (f_value + f_value_1 - f_value_2 - f_value_2)
               end if
               if (alpha_scale .lt. 0.1d0) then
                  alpha_scale = 0.1d0
               end if
               alpha = alpha_scale*alpha

               try_arguments = arguments + alpha*s
               n_fn_evals = n_fn_evals + 1
               call f_external(n_arguments, try_arguments, &
                    f_value_1, n_fn_evals, &
                    hamiltonian, &
                    overlap_matrix,KS_eigenvalue, KS_eigenvector, &
                    occ_numbers, &
                    n_electrons, chemical_potential, constraint_proj, &
                    electrons_in_region, &
                    current_spin &
                    )
               if (f_value_1.lt.f_value) then
                  f_value = f_value_1
                  arguments = try_arguments
                  total_alpha = total_alpha + alpha
                  global_improvement = .true.
               end if

            end if

            if (alpha.lt.aeps) then
               shorter_step = .false.
            end if

         enddo

!     next, update the Hessian according to the BFGS formula
         if (global_improvement) then

            call evaluate_gradient_constraint(n_arguments, arguments, &
                 f_value, &
                 n_fn_evals, delta, arg_mag, f_external, gradient_work, &
                 method, &
                 hamiltonian, &
                 overlap_matrix, KS_eigenvalue, KS_eigenvector, &
                 occ_numbers, &
                 n_electrons, chemical_potential, constraint_proj, &
                 electrons_in_region, &
                 current_spin &
                 )

            s_times_g = 0.0d0
            do i_index = 1, n_arguments, 1
               s_times_g = s_times_g + &
                    s(i_index)*gradient_work(i_index)
            enddo

            s_times_y = s_times_g - sHs

            diff_in_f = f_value

!     update the Hessian only if (s_times_y.gt.0.0d0).and.(sHs.lt.0.0d0)
            if ((s_times_y.gt.0.0d0).and.(sHs.lt.0.0d0)) then

               if (s_times_y + total_alpha*sHs.le.0.0d0) then

!     normal update of the Hessian
                  scaled_update = .false.

                  sigma = -1.0d0/sHs

!                  print *, 'non-scaled update: ',
!     +                 s_times_y + total_alpha*sHs
!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, gradient, scaled_update)


!                  print *,'total_alpha', total_alpha
                  sigma = 1.0d0/(total_alpha*s_times_y)
                  y = gradient_work - gradient

!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, y, scaled_update)


               else

                  scaled_update = .false.
!                  scaled_update = .true.

!                  sigma = total_alpha/(s_times_y - total_alpha*sHs)

                  sigma = -1.0d0/sHs

!                  print *, 'scaled update: ',
!     +                 s_times_y + total_alpha*sHs
!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, gradient, scaled_update)

                  scaled_update = .false.

!                  y = gradient_work + (s_times_y*sigma - 1.0d0)*gradient
!                  sigma = 1.0d0 / (sigma*(s_times_y**2))

                  sigma = 1.0d0/s_times_y
                  y = gradient_work - gradient

!                  print *,'total_alpha', total_alpha
!                  print *, 'Sigma: ', sigma

                  call update_cholesky_factors(hessian_d, hessian_l, &
                       n_arguments, sigma, y, scaled_update)

               end if

            end if

         else

            if (method.eq.1) then
               method = 2
!               delta = delta0
               call evaluate_gradient_constraint(n_arguments, &
                    arguments, &
                    f_value, &
                    n_fn_evals, delta, arg_mag, f_external, &
                    gradient, &
                    method, hamiltonian, &
                    overlap_matrix, KS_eigenvalue, KS_eigenvector, &
                    occ_numbers, &
                    n_electrons, chemical_potential, &
                    constraint_proj, &
                    electrons_in_region, &
                    current_spin &
                    )

!$$$            elseif (method.eq.2) then
!$$$               method = 3
!$$$               delta = delta0
!$$$               call evaluate_gradient_constraint(n_arguments,
!$$$     +              arguments,
!$$$     +              f_value,
!$$$     +              n_fn_evals, delta, arg_mag, f_external,
!$$$     +              gradient,
!$$$     +              method, hamiltonian,
!$$$     +              overlap_matrix, KS_eigenvalue, KS_eigenvector,
!$$$     +              occ_numbers,
!$$$     +              n_electrons, chemical_potential,
!$$$     +              constraint_proj,
!$$$     +              electrons_in_region,
!$$$     +              current_spin
!$$$     +              )

            elseif (delta.gt.deps) then
               delta = delta/2.0d0
               method = 1
               call evaluate_gradient_constraint(n_arguments, &
                    arguments, &
                    f_value, &
                    n_fn_evals, delta, arg_mag, f_external, &
                    gradient, &
                    method, hamiltonian, &
                    overlap_matrix, KS_eigenvalue, KS_eigenvector, &
                    occ_numbers, &
                    n_electrons, chemical_potential, &
                    constraint_proj, &
                    electrons_in_region, &
                    current_spin &
                    )
            end if

!     end update hessian
         end if

!     test for convergence
         converged = (f_value.lt.eps)

!     main BFGS loop
      enddo

!     finally, one more evaluation of the target function is needed
!     to update all the variables - like eigenvalues, chemical potentials
!     and constraint potentials
      n_fn_evals = n_fn_evals + 1
      call f_external(n_arguments, arguments, &
           f_value, n_fn_evals, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )

!      print *, 'On exit, f_value: ', f_value
!      print *, 'total_alpha: ', total_alpha
!      print *, 'alpha: ', alpha
!      print *, 'delta: ', delta

      end subroutine bfgs_constraint_v2
!******
!----------------------------------------------------------------------------------------------
!****s* bfgs/evaluate_gradient_constraint
!  NAME
!    evaluate_gradient_constraint
!  SYNOPSIS
      subroutine evaluate_gradient_constraint(n_arguments, arguments, &
           f_value, n_fn_evals, delta, arg_mag, f_external, gradient, &
           method, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )
!  PURPOSE
!    Evaluate an approximation to the gradient of f at the points given by
!    the array arguments using finite differences.
!
!  ARGUMENTS
      integer :: n_arguments
      real*8 :: arguments(n_arguments)
      real*8 hamiltonian( n_basis*(n_basis+1)/2, n_spin )
      real*8 overlap_matrix( n_basis*(n_basis+1)/2 )
      real*8, dimension(n_states, n_spin) :: KS_eigenvalue
      real*8, dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
      real*8, dimension(n_states, n_spin) :: occ_numbers
      real*8 :: n_electrons
      real*8 :: chemical_potential
      real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
      real*8, dimension(n_region, n_spin) :: electrons_in_region
      integer :: current_spin
      external :: f_external
      real*8 :: f_value
      real*8 :: delta
      real*8, dimension(n_arguments) :: arg_mag
      integer :: n_fn_evals
      integer :: method
      real*8 :: gradient(n_arguments)
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o n_arguments -- number of arguments for the objective function
!    o arguments -- arguments for the objective function
!    o max_fn_evals -- maximum number of evaluations of the objective function
!    o eps -- tolerance of the objective function, i.e., |f(x)| < eps
!    o f_value -- value of the objective function
!    o delta -- basic differencing parameter for the evaluation of the gradient of f
!    o f_external -- external routine to evaluate the objective function f
!    o hamiltonian -- the Hamilton matrix
!    o overlap_matrix -- the overlap matrix
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o occ_numbers -- occupation numbers
!    o n_electrons -- total number of electrons
!    o chemical_potential -- current chemical potential
!    o constraint_proj -- projector to the constrained region
!    o electrons_in_region -- number of electrons in the constrained regions
!    o current_spin -- current spin channel
!    o method -- selector for the finite difference scheme to use
!    o arg_mag -- magnitude of the arguments
!  OUTPUT
!    o gradient -- difference approximation to gradient of f
!  SOURCE

!     locals
      real*8 :: f_value_1, f_value_2, f_value_3, f_value_4
      integer :: i_index
      real*8 :: delta_scaled
      real*8 :: temp_arg

      do i_index = 1, n_arguments, 1

         delta_scaled = delta*arg_mag(i_index)
         temp_arg = arguments(i_index)

         select case(method)

      case(1)

         arguments(i_index) = temp_arg + delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         gradient(i_index) = (f_value_1 - f_value) / delta_scaled

      case(2)

         arguments(i_index) = temp_arg + delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg - delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, arguments, &
              f_value_2, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         gradient(i_index) = (f_value_1 - f_value_2) / &
              (2.0d0*delta_scaled)

      case(3)

         arguments(i_index) = temp_arg + 2.0d0*delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, arguments, &
              f_value_1, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg + delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, arguments, &
              f_value_2, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg - delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, arguments, &
              f_value_3, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )

         arguments(i_index) = temp_arg - 2.0d0*delta_scaled
         n_fn_evals = n_fn_evals + 1
         call f_external(n_arguments, arguments, &
              f_value_4, n_fn_evals, &
              hamiltonian, &
              overlap_matrix,KS_eigenvalue, KS_eigenvector, occ_numbers, &
              n_electrons, chemical_potential, constraint_proj, &
              electrons_in_region, &
              current_spin &
              )


         gradient(i_index) = &
            (f_value_4 - 8.0d0*f_value_3 + 8.0d0*f_value_2 - f_value_1)/ &
            (12.0d0*delta_scaled)

      end select

      arguments(i_index) = temp_arg

      enddo

      end subroutine evaluate_gradient_constraint
!******
!-------------------------------------------------------------------------------------------------
!****s* bfgs/update_cholesky_factors
!  NAME
!    update_cholesky_factors
!  SYNOPSIS
      subroutine update_cholesky_factors(d, l, n, sigma, z, scaled)
!  PURPOSE
!    Compute a rank-one update to the Cholesky factors L and D of the matrix
!    factorization LDL^T (for the purposes of BFGS-algorithm)
!  ARGUMENTS
      integer :: n
      real*8 :: z(n)
      real*8 :: sigma
      logical :: scaled
      real*8 :: d(n)
      real*8 :: l(n,n)
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o n -- dimension of the matrices
!    o z -- vector of the rank-one update: + sigma (z z^T)
!    o sigma -- scalar factor of the rank-one update: + sigma (z z^T)
!    o scaled -- flag indicating the the update factor sigma should be regularized
!    o d -- diagonal of the Cholesky factorization LDL^T
!    o l -- lower triangle of the Cholesky factorization LDL^T
!  OUTPUT
!    none, on exit factors l and d are updated
!  SOURCE

!     locals
      real*8 :: d2(n)
      real*8 :: l2(n,n)

      integer :: i,j,k
      real*8 :: sum1, sum2

      real*8 :: macheps
      real*8, external :: d1mach

      d2 = 0.0d0
      l2 = 0.0d0

      macheps = 10.0d0 * d1mach(4)

      if (scaled) then
         sigma = sigma + macheps
      end if

      do j = 1, n, 1

         sum1 = 0.0d0
         sum2 = 0.0d0
         do k = 1, j-1, 1
            d2(j) = d2(j) + l(j,k)*d(k)*l(j,k)
!            sum1 = sum1 + l(j,k)*d(k)*l(j,k)
            d2(j) = d2(j) - l2(j,k)*d2(k)*l2(j,k)
!            sum2 = sum2 + l2(j,k)*d2(k)*l2(j,k)
         enddo

!         d2(j) = sum1 - sum2 + d(j) + sigma*z(j)**2
         d2(j) = d2(j) + d(j) + sigma*z(j)**2
         l2(j,j) = 1.0d0

         if (abs(d2(j)).lt.macheps) then
            d2(j) = d2(j) + macheps
         end if

         do i = j+1, n, 1
            sum1 = 0.0d0
            sum2 = 0.0d0
            do k = 1, j-1, 1
!               sum1 = sum1 + l(j,k)*d(k)*l(j,k)
               l2(i,j) = l2(i,j) + l(j,k)*d(k)*l(j,k)
!               sum2 = sum2 + l2(j,k)*d2(k)*l2(j,k)
               l2(i,j) = l2(i,j) - l2(j,k)*d2(k)*l2(j,k)
            enddo

!            l2(i,j) = ( sum1 - sum2 + d(j)*l(i,j) + sigma*z(i)*z(j) )
!     +           / d2(j)
            l2(i,j) = l2(i,j) + d(j)*l(i,j) + sigma*z(i)*z(j)
            l2(i,j) = l2(i,j)/d2(j)

         enddo

      enddo

      d = d2
      l = l2

      end subroutine update_cholesky_factors
!******
!---------------------------------------------------------------------------------
!****s* bfgs/solve_factored_problem
!  NAME
!    solve_factored_problem
!  SYNOPSIS
      subroutine solve_factored_problem(d, l, n, b)
!  PURPOSE
!     Solves the symmetric linear system of equations LDL^Tx = b
!     when the matrix is given in a factored form. It is assumed that
!     l(i,i) = 1.0d0, and d(i) is not zero
!  ARGUMENTS
      integer :: n
      real*8 :: d(n)
      real*8 :: l(n,n)
      real*8 :: b(n)
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o n -- dimension of the matrices
!    o d -- diagonal of the Cholesky factorization LDL^T
!    o l -- lower triangle of the Cholesky factorization LDL^T
!    o b -- the right hand side
!  OUTPUT
!    none, on exit b is overwritten by the solution x
!  SOURCE

!     locals
      real*8 :: x(n)

      integer :: i,j

!     first, do the forward substitution to solve Lz = b
      do i = 1, n, 1

         x(i) = b(i)

         do j = 1, i-1, 1
            x(i) = x(i) - l(i,j)*x(j)
         enddo

      enddo

!     next, solve the diagonal problem Dy = z
      do i = 1, n, 1
         x(i) = x(i)/d(i)
      end do

!     finalle, do the backward substitution to solve L^T x = y
      do i = n, 1, -1

         do j = i+1, n, 1
            x(i) = x(i) - l(j,i)*x(j)
         enddo

      enddo

      b = x

      end subroutine solve_factored_problem
!******

      end module bfgs
