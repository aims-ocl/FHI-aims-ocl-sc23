!****s* FHI-aims/get_aux_gaussian_basis_fns
!  NAME
!   get_aux_gaussian_basis_fns
!  SYNOPSIS

      subroutine get_aux_gaussian_basis_fns &
        ( i_species &
        )

!  PURPOSE
!   Subroutine get_aux_gaussian_basis_fns tabulates cartesian Gaussian
!   radial wave functions as the auxiliary basis functions for HF and
!   HF-based calculations

      use dimensions
      use grids
      use species_data
      use constants

      implicit none

!  ARGUMENTS
      integer i_species

!  INPUTS
!    o i_species -- species for which to calculate basis fns
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

!  local variables

      integer cartesian_l
      integer l_shell
      real*8 alpha
      real*8 coeff
      real*8 norm

      real*8 elem_gaussian(n_max_grid)
      real*8 wave(n_max_grid)

!  counters

      integer i_grid, i_gaussian, i_contracted

!  functions

      real*8 int_log_mesh

!  begin work

!     Tabulate Gaussian radial functions and kinetic energy term

      do i_gaussian = 1, n_aux_gaussian(i_species), 1

!       Gaussians for the same cartesian quantum number L
!       are tabulated in order l = L, L-2, L-4, ...
!       Once we have calculated the function for l=L, we only need
!       to copy it for the other instances.
        if ( aux_gaussian_n(i_species, i_gaussian) .eq. &
             aux_gaussian_l(i_species, i_gaussian)      ) then
!         first instance; calculate everything

          do i_grid = 1, n_grid(i_species), 1
            wave (i_grid) = 0.d0
          enddo

!         tabulate unnormalized Gaussian radial function, second derivative
          cartesian_l = aux_gaussian_n(i_species, i_gaussian)
          do i_contracted = 1, &
                    aux_gaussian_n_contr(i_species, i_gaussian),1
            alpha = &
              aux_gaussian_alpha(i_species, i_gaussian, i_contracted)
            coeff = &
              aux_gaussian_coeff(i_species, i_gaussian, i_contracted)

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
!              if (cartesian_l.eq.0) then
!                wave_sec_deriv ( i_grid ) =
!     +          wave_sec_deriv ( i_grid ) +
!     +          norm * coeff *
!     +          ( - 2.d0 * alpha * 3.d0 *
!     +              r_grid(i_grid, i_species)
!     +            + 4.d0 * (alpha**2.d0) *
!     +              r_grid(i_grid,i_species)**3.d0 )
!     +            * exp (- alpha * r_grid(i_grid,i_species)**2.d0)
!              else
!                wave_sec_deriv ( i_grid ) =
!     +          wave_sec_deriv ( i_grid ) +
!     +          norm * coeff *
!     +          (   cartesian_l * (cartesian_l+1) *
!     +              (r_grid(i_grid,i_species)**(cartesian_l-1))
!     +            - 2.d0 * alpha * (2*cartesian_l+3) *
!     +              (r_grid(i_grid,i_species)**(cartesian_l+1))
!     +            + 4.d0 * (alpha**2.d0) *
!     +              (r_grid(i_grid,i_species)**(cartesian_l+3))
!     +          )
!     +          * exp (- alpha * r_grid(i_grid,i_species)**2.d0)
!              end if

            enddo

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

          enddo

        end if

!       now store radial function, kinetic energy fragment for all Gaussian partials
!       which belong to the same cartesian quantum number
        do i_grid = 1, n_grid(i_species),1
           aux_gaussian_wave (i_grid, i_species, i_gaussian) &
              = wave (i_grid)
        enddo

      enddo

!     output Gaussian basis data
!      write(use_unit,*)
!      write(use_unit,*)"List of auxilliary cartesian Gaussian basis orbitals: "
!      write(use_unit,'(4X,A,4X,A)')
!     +  "L", "l"
!      do i_gaussian = 1, n_aux_gaussian(i_species), 1
!        write(use_unit,'(2X,I3,2X,I3)')
!     +    aux_gaussian_n(i_species,i_gaussian),
!     +    aux_gaussian_l(i_species,i_gaussian)
!      enddo
!      write(use_unit,*)

!  that's all folks

      return
    end subroutine get_aux_gaussian_basis_fns

!----------------------------------------------------------------------
!*****************
