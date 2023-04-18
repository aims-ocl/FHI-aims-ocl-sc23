!---------------------------------------------------------------------
!  Subroutine evaluate_xc_energy_shell
!  adds up the XC energy contributions from several grid points at a time.
!

      subroutine evaluate_xc_energy_shell &
           ( n_points, partition, en_density_xc, &
             local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, &
             local_rho, local_rho_gradient, local_kinetic_density, &
             en_xc, en_pot_xc &
           )

      use dimensions
      use energy_density
      use runtime_choices

      implicit none

!  imported variables

      ! input

      integer :: n_points
      real*8 partition(n_points)

      real*8, dimension(n_points) :: en_density_xc
      real*8, dimension(n_spin, n_points) :: local_xc_derivs
      real*8, dimension(3,n_spin,n_points) :: xc_gradient_deriv
      real*8, dimension(n_spin, n_points) :: xc_tau_deriv

      real*8, dimension(n_spin,n_points) :: local_rho
      real*8, dimension(3,n_spin,n_points) :: local_rho_gradient
      real*8, dimension(n_spin,n_points) :: local_kinetic_density

      ! output

      real*8 :: en_xc
      real*8 :: en_pot_xc

!  local variables

      real*8, dimension(n_points) :: energy_term

      ! counters

      integer :: i_point, i_spin, i_coord

      ! functions

      real*8, external :: ddot

!  begin work

      energy_term = 0.d0
      do i_point = 1, n_points, 1

        do i_spin = 1, n_spin, 1
          energy_term(i_point) = energy_term(i_point) + &
            en_density_xc(i_point) * local_rho(i_spin,i_point)
        enddo

      enddo

      en_xc = en_xc + &
        ddot(n_points,energy_term,1,partition,1)

      energy_term = 0.d0

      do i_point = 1, n_points, 1

        do i_spin = 1, n_spin, 1
          energy_term(i_point) = energy_term(i_point) + &
            local_xc_derivs(i_spin,i_point) * &
            local_rho(i_spin,i_point)
        enddo

      enddo

      if (use_gga) then
        do i_point = 1, n_points, 1

          do i_spin = 1, n_spin, 1
            do i_coord = 1, 3, 1
              energy_term(i_point) = energy_term(i_point) + &
                2.d0 * xc_gradient_deriv(i_coord,i_spin,i_point) * &
                local_rho_gradient(i_coord,i_spin,i_point)
            enddo
          enddo

        enddo
      end if

      ! Needs the first integration catch until we have a way to evaluate meta-GGAs
      ! for free-atoms. AJL/July 2016
      if (use_meta_gga.and..not.first_integration) then
        do i_point = 1, n_points, 1

          do i_spin = 1, n_spin, 1
             energy_term(i_point) = energy_term(i_point) + &
                xc_tau_deriv(i_spin,i_point) * &
                local_kinetic_density(i_spin,i_point)
          enddo

        enddo
      end if

      en_pot_xc = en_pot_xc + &
        ddot(n_points,energy_term,1,partition,1)

      ! CC: Save local (!) v_xc(r) * rho(r) if Harris-Foulkes-like energy density shall be computed
      if (flag_harris_foulkes_energy_density) then

        ! AJL: Added this if statement as free_energy_superpos does not allocate this array yet.
        ! Is it a requirement for me to allocate this array and assess this energy density?
        ! Or can it be that actually we shouldn't be doing this sum anyway at the start of the calculation?
        ! My concern is this shouldn't be dependent on an "allocated" check, and so any suggestions/corrections
        ! are welcomed as I don't know what is best here.

        if (allocated(ed_i_point_to_i_full_point_map)) then
          do i_point = 1, n_points, 1
            ed_xc_pot_energy_density(ed_i_point_to_i_full_point_map(i_point)) = energy_term(i_point)
          end do
        end if
      end if

      end subroutine evaluate_xc_energy_shell
!---------------------------------------------------------------------
