! VB Subroutine limit_rho_change assesses the length of a charge density step,
!    and reduces it by force if the proposed step length is too large.
!
! This must be done in a separate subroutine outside mixing.f because the 
! arguments are different depending on the mixing strategy
! (linear or pulay)

  subroutine limit_rho_change &
  ( partition_tab, mixer_type, rho_diff, rho_gradient_diff, kinetic_density_diff, charge_multiplier )

    use dimensions
    use runtime_choices
    use localorb_io
    use grids
    use mpi_utilities
    use synchronize_mpi

    implicit none

    ! input variables

    real*8, dimension(n_full_points) :: partition_tab
    integer :: mixer_type

    ! input / output variables

    real*8, dimension( n_full_points, n_spin ) :: rho_diff
    real*8, dimension( 3, n_full_points, n_spin ) :: rho_gradient_diff
    real*8, dimension( n_full_points, n_spin ) :: kinetic_density_diff

    ! output variables

    real*8, dimension(n_spin) :: charge_multiplier

    ! local variables

    real*8, dimension(n_spin) :: rho_change

    character*100 :: info_str

    integer :: i_offset, i_my_batch, i_index
    integer :: i_spin

    ! begin work

    charge_multiplier = 0.d0

    ! Determine the charge density norm ...

    rho_change = 0.d0

    i_offset = 0

    do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1
           
          i_offset = i_offset + 1
!         execute only if partition_tab.gt.0 here, i.e. if the integration point
!         makes sense
          if (partition_tab(i_offset).gt.0.d0) then

            do i_spin = 1, n_spin, 1

              rho_change(i_spin) = rho_change(i_spin) +  &
                partition_tab(i_offset) *  &
                rho_diff(i_offset, i_spin)**2
                       
            enddo

          end if

        ! end loop over current batch
        enddo

      ! end if batch in tasklist
      !end if

    ! end loop over all grid batches
    enddo

    ! broadcast the result to all threads
    call sync_density(rho_change)

    do i_spin = 1, n_spin, 1
      rho_change(i_spin) = sqrt(rho_change(i_spin))
    enddo

    if (mixer_type.eq.0) then
      ! Linear mixing - must factor in the mixing factor before
      ! thresholding

      do i_spin = 1, n_spin, 1
        rho_change(i_spin) = rho_change(i_spin) * linear_mix_param(i_spin)
      enddo

    end if

    ! Now check whether charge or spin density hit their respective mixer thresholds ...
    if ( (max_rho_change(1).gt.0.d0) .and. (max_rho_change(1).lt.rho_change(1)) ) then
      ! Must throttle charge density change!

      write(info_str, '(2X,A,E14.6)' ) & 
      "Change of charge density after mixing: ", rho_change(1)
      call localorb_info( info_str )      

      write(info_str, '(2X,A,E14.6,A)' ) & 
      "Throttling charge density change to ", max_rho_change(1), " ."
      call localorb_info( info_str )      

      charge_multiplier(1) = max_rho_change(1) / rho_change(1)

      rho_diff (:,1) = charge_multiplier(1) * rho_diff (:,1)

      if ( use_density_gradient ) then
        rho_gradient_diff (:,:,1) = charge_multiplier(1) * rho_gradient_diff (:,:,1)
      end if

      if (use_meta_gga) then
        kinetic_density_diff(:,1) = charge_multiplier(1) * kinetic_density_diff(:,1)
      end if

    end if

    if (spin_treatment.eq.1) then
      if ( (max_rho_change(2).gt.0.d0) .and. (max_rho_change(2).lt.rho_change(2)) ) then
        ! Must throttle spin density change!

        write(info_str, '(2X,A,E14.6)' ) & 
        "Change of spin density after mixing:  ", rho_change(2)
        call localorb_info( info_str )      

        write(info_str, '(2X,A,E14.6,A)' ) & 
        "Throttling spin density change to  ", max_rho_change(2), " ."
        call localorb_info( info_str )      

        charge_multiplier(2) = max_rho_change(2) / rho_change(2)

        rho_diff (:,2) = charge_multiplier(2) * rho_diff (:,2)

        if ( use_density_gradient ) then
          rho_gradient_diff (:,:,2) = charge_multiplier(2) * rho_gradient_diff (:,:,2)
        end if

        if (use_meta_gga) then
          kinetic_density_diff(:,2) = charge_multiplier(2) * kinetic_density_diff(:,2)
        end if
      
      end if
    end if

  end subroutine limit_rho_change
