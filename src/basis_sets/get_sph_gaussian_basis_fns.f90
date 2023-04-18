!****s* FHI-aims/get_sph_gaussian_basis_fns
!  NAME
!   get_sph_gaussian_basis_fns
!  SYNOPSIS
      subroutine get_sph_gaussian_basis_fns( i_species )
!  PURPOSE
!  Subroutine get_gaussian_basis_fns tabulates cartesian Gaussian
!  radial wave functions as polarisation functions for one species
!  USES
      use dimensions,   only : n_max_grid, use_basis_gradients
      use grids,        only : n_grid, r_grid
      use species_data, only : gaussian_wave, gaussian_kinetic, gaussian_wave_deriv, &
                               gaussian_n, gaussian_alpha, gaussian_coeff, &
                               gaussian_n_contr, gaussian_l, n_gaussian, gaussian_l
      use mpi_tasks,    only : myid
      use localorb_io,  only : use_unit
      use psi_at_nucleus_mod, only: psi_at_nucleus_gauss
      implicit none
!  INPUTS
!    o i_species -- species number in question
!  OUTPUTS
!    none
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


!  imported variables

!  input

      integer i_species

!  local variables

      integer cartesian_l
      integer l_shell
      real*8 alpha
      real*8 coeff
      real*8 norm

      real*8 elem_gaussian(n_max_grid)
      real*8 wave(n_max_grid)
      real*8 wave_deriv(n_max_grid)
      real*8 wave_sec_deriv(n_max_grid)

      ! For computing wavefunction at nucleus
      real*8, allocatable :: coeffs(:), alphas(:)
      real*8 :: total_norm
      integer :: ig

!  counters

      integer i_grid, i_gaussian, i_contracted

!  functions

      real*8 int_log_mesh

!  begin work


!     Tabulate Gaussian radial functions and kinetic energy term

      do i_gaussian = 1, n_gaussian(i_species), 1

!       Gaussians for the same cartesian quantum number L

!       Once we have calculated the function for l=L, we only need
!       to copy it for the other instances.
!        if ( gaussian_n(i_species, i_gaussian) .eq.
!     +       gaussian_l(i_species, i_gaussian)      ) then
!         first instance; calculate everything

          do i_grid = 1, n_grid(i_species), 1
            wave (i_grid) = 0.d0
            wave_deriv (i_grid) = 0.d0
            wave_sec_deriv (i_grid) = 0.d0
          enddo

!         tabulate unnormalized Gaussian radial function, second derivative
          cartesian_l = gaussian_n(i_species, i_gaussian)
          do i_contracted = 1,gaussian_n_contr(i_species, i_gaussian),1
            alpha = gaussian_alpha(i_species, i_gaussian, i_contracted)
            coeff = gaussian_coeff(i_species, i_gaussian, i_contracted)

            do i_grid = 1, n_grid(i_species), 1
              elem_gaussian(i_grid) = &
              (r_grid(i_grid,i_species)**(cartesian_l+1)) * &
              exp (- alpha * r_grid(i_grid,i_species)**2.d0)
            enddo

!           normalize elementary gaussian
            norm = &
              int_log_mesh ( &
              elem_gaussian, &
              elem_gaussian, &
              n_grid(i_species), r_grid(1,i_species) &
                         )

            norm = 1.d0/sqrt(norm)

!           add elementary Gaussian to contracted fn.
            do i_grid = 1, n_grid(i_species), 1
              wave ( i_grid ) = &
              wave ( i_grid ) + &
              norm * coeff * elem_gaussian(i_grid)

!             add up the second derivative in the same way
              if (cartesian_l.eq.0) then
                wave_sec_deriv ( i_grid ) = &
                wave_sec_deriv ( i_grid ) + &
                norm * coeff * &
                ( - 2.d0 * alpha * 3.d0 * &
                    r_grid(i_grid, i_species) &
                  + 4.d0 * (alpha**2.d0) * &
                    r_grid(i_grid,i_species)**3.d0 ) &
                  * exp (- alpha * r_grid(i_grid,i_species)**2.d0)
              else
                wave_sec_deriv ( i_grid ) = &
                wave_sec_deriv ( i_grid ) + &
                norm * coeff * &
                (   cartesian_l * (cartesian_l+1) * &
                    (r_grid(i_grid,i_species)**(cartesian_l-1)) &
                  - 2.d0 * alpha * (2*cartesian_l+3) * &
                    (r_grid(i_grid,i_species)**(cartesian_l+1)) &
                  + 4.d0 * (alpha**2.d0) * &
                    (r_grid(i_grid,i_species)**(cartesian_l+3)) &
                ) &
                * exp (- alpha * r_grid(i_grid,i_species)**2.d0)
              end if

            enddo

!           if needed, compute derivative
            if (use_basis_gradients) then
              do i_grid = 1, n_grid(i_species), 1
!               treat L=0 separately to avoid r^0
                if (cartesian_l.eq.0) then
                  wave_deriv ( i_grid ) = &
                  wave_deriv ( i_grid ) + &
                  norm * coeff * &
                  (   1.d0 &
                    - 2.d0 * alpha * &
                      (r_grid(i_grid,i_species)**2.d0) &
                  ) &
                  * exp (- alpha * r_grid(i_grid,i_species)**2.d0)
                else
                  wave_deriv ( i_grid ) = &
                  wave_deriv ( i_grid ) + &
                  norm * coeff * &
                  (   (cartesian_l+1) * &
                      (r_grid(i_grid,i_species)**(cartesian_l)) &
                    - 2.d0 * alpha * &
                      (r_grid(i_grid,i_species)**(cartesian_l+2)) &
                  ) &
                  * exp (- alpha * r_grid(i_grid,i_species)**2.d0)
                end if
              enddo
            end if

          enddo

!         check normalization for entire Gaussian
          norm = &
            int_log_mesh ( &
            wave, &
            wave, &
            n_grid(i_species), r_grid(1,i_species) &
                         )

          norm = 1.d0/sqrt(norm)

!         normalize Gaussian function, second derivative
          do i_grid = 1, n_grid(i_species), 1

            wave (i_grid) = &
            wave (i_grid) * norm

            wave_sec_deriv (i_grid) = &
            wave_sec_deriv (i_grid) * norm

          enddo

          if (use_basis_gradients) then
!           normalize derivative also
            do i_grid = 1, n_grid(i_species), 1
              wave_deriv (i_grid) = &
              wave_deriv (i_grid) * norm
            enddo
          end if

!        end if

!       now store radial function, kinetic energy fragment for all Gaussian partials
!       which belong to the same cartesian quantum number
        l_shell = gaussian_l (i_species, i_gaussian)
        do i_grid = 1, n_grid(i_species),1
          gaussian_wave (i_grid, i_species, i_gaussian) = wave (i_grid)
          gaussian_kinetic (i_grid, i_species, i_gaussian) = &
            - 0.5d0 * wave_sec_deriv (i_grid) &
            + 0.5d0 * l_shell * (l_shell+1) &
              / (r_grid(i_grid,i_species)**2.d0) &
              * wave(i_grid)
        enddo

        if (use_basis_gradients) then
          do i_grid = 1, n_grid(i_species),1
            gaussian_wave_deriv (i_grid, i_species, i_gaussian) = &
              wave_deriv (i_grid)
          enddo
        end if

        if (cartesian_l == 0) then
           ! See also comments in set_gaussian_free_cut
           allocate(alphas(gaussian_n_contr(i_species, i_gaussian)))
           allocate(coeffs(gaussian_n_contr(i_species, i_gaussian)))
           alphas = gaussian_alpha(i_species, i_gaussian, :size(alphas))
           coeffs = gaussian_coeff(i_species, i_gaussian, :size(coeffs))
           coeffs = coeffs/sqrt(norm0(alphas))
           total_norm = sum(coeffs**2*norm0(alphas))
           do ig = 1, size(alphas)
              total_norm = total_norm + sum(2*coeffs(ig)*coeffs(ig+1:)* &
                   & norm0((alphas(ig)+alphas(ig+1:))/2))
           end do
           psi_at_nucleus_gauss(i_gaussian,i_species) = &
                & sum(coeffs)/sqrt(total_norm)
           deallocate(alphas, coeffs)
        else
           psi_at_nucleus_gauss(i_gaussian,i_species) = 0d0
        end if

      enddo

!     output Gaussian basis data
      if (myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*) " List of cartesian Gaussian basis orbitals: "
        write(use_unit,'(4X,A,4X,A)') &
        "L", "l"
        do i_gaussian = 1, n_gaussian(i_species), 1
          write(use_unit,'(2X,I3,2X,I3)') &
          gaussian_n(i_species,i_gaussian), &
          gaussian_l(i_species,i_gaussian)
        enddo
        write(use_unit,*)
      end if

!  that's all folks

      return

    contains
      ! Norm of an elemental Gaussian for l=0
      elemental real*8 function norm0(alpha) result(y)
        use constants, only: pi
        real*8, intent(in) :: alpha
        y = sqrt(pi/(128*alpha**3))
      end function norm0
    end subroutine get_sph_gaussian_basis_fns
!******
