!****h* FHI-aims/hartree_fock
!  NAME
!    lrc_pt2
!  SYNOPSIS

      module lrc_pt2

!  PURPOSE  
!  this module contains the subroutines for the long range corrected PT2 functional.  
!
!  Subroutines:
!
!  USES

implicit none

!  AUTHOR
!    Igor Ying Zhang, FHI
!  SEE ALSO
!    
!  COPYRIGHT
!   
!  HISTORY
!    Introduced December 2016
!
!  SOURCE

!   

     

contains
!########################################################   
      
 subroutine get_ovlp_3fn_for_lrc_pt2 &
             ()

    !  PURPOSE
    !  Calculate ovlp_3fn matricies for lrc-PT2 part
    !
    !  AUTHOR
    !  Igor Ying Zhang, FHI
    !
    !  INPUTS
    !  none
    !  
    !  OUTPUTS
    !
    !  USES

    use dimensions
    use runtime_choices
    use prodbas
    use hartree_fock
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    

    character(*), parameter :: func = 'get_ovlp_3fn_for_lrc_pt2'
        
    ! In the futur we may also change the other RI-Types but for now
    ! the RI_V method has to be sufficient
    select case (RI_type)
    !   case(RI_SVS)
    !     call get_coeff_3fn_svs(ovlp_3fn)
      case(RI_V)
        if (use_2d_corr) then
          call get_coeff_3fn_v_2d_lrc_pt2()
        else
          call get_coeff_3fn_v_1d_lrc_pt2()
        end if
    !   case (RI_LVL)
    !     call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr_lvl)
    !   case (RI_LVL_full)
    !     call get_coeff_3fn_lvl_full(ovlp_3fn)
    !   case (RI_LVL_2nd)
    !     call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr_lvl)
        ! Bare Coulomb integrals to ovlp_3fn
    !     call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, OVLP_TYPE_COULOMB)
      case default
        call aims_stop("Invalid version of RI (resolution of identity)", func)
    end select 
        

  end subroutine get_ovlp_3fn_for_lrc_pt2
        
!##################################################
  subroutine get_coeff_3fn_v_1d_lrc_pt2()

    !  PURPOSE
    !
    !    Calculate the orthonormalized product expansion coefficients in RI-V:
    !
    !       C_{ij}^\mu := (\phi_i \phi_j | P_\mu) V_{\mu\nu}^-0.5
    !
    !    Uses 1d distributed matrices internally.
    !
    !  USES
    !  
    !  SOURCE


    use timing
    use prodbas
    use dimensions
    use sbt_overlap_aims
    use species_data
    use runtime_choices
    use hartree_fock
    use mpi_tasks, only: check_allocation
    use localorb_io, only: localorb_info
    
    implicit none

    !  ARGUMENTS

    !  INPUTS
    !    none
    !  OUTPUTS
    !
    !  AUTHOR
    !    Igor Ying Zhang, FHI
    !  SOURCE
    !    get_coeff_3fn_v_1d(ovlp_3fn)  -- FHIaims

    real*8, allocatable :: coulomb_matr(:,:)
    integer :: info
    integer :: i, j 
    integer*8 :: array_dimension
    character(*), parameter :: func = 'get_coeff_3fn_v_1d_lrc_pt2'

    ! Integrate ovlp3fn (T)
    call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
    !for LR
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
    !ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn )


    ! Integrate V
    call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
    !for LR
    allocate(coulomb_matr(n_basbas, n_loc_prodbas), stat=info)
    call check_allocation(info, 'coulomb_matr', func)
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    if (use_logsbt) then
       call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .false.)
    else
       call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr)
    end if
    !ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_coulomb_matr, clock_time_coulomb_matr, tot_time_coulomb_matr, tot_clock_time_coulomb_matr)
        
    ! V^-0.5
    call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
    ! for LR
    if (use_asym_inv_sqrt) then
       ! transposed==.true. because get_v_multi_ovlp3fn() does
       !   C := V^-0.5T.
       call asym_inv_sqrt_of_auxmat_scalapack(coulomb_matr, "Coulomb_LR", .true.)
    else if (use_scalapack) then
       call power_auxmat_scalapack(coulomb_matr, -0.5d0, "Coulomb_LR")
    else
       call power_auxmat_lapack(coulomb_matr, -0.5d0, "Coulomb_LR")
    end if 
    call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr,tot_time_inv_coulomb_matr, tot_clock_time_inv_coulomb_matr)

    ! C := T V^-0.5
    call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
    ! for LR
    call get_v_multi_ovlp3fn(coulomb_matr, ovlp_3fn)
    call get_times(time_ovlp_multi, clock_time_ovlp_multi,tot_time_ovlp_multi, tot_clock_time_ovlp_multi)
    
    deallocate(coulomb_matr)

    call localorb_info('')

  end subroutine get_coeff_3fn_v_1d_lrc_pt2
!##################################################################
  subroutine get_coeff_3fn_v_2d_lrc_pt2()

    !  PURPOSE
    !
    !    Calculate the orthonormalized product expansion coefficients in RI-V:
    !
    !       C_{ij}^\mu := (\phi_i \phi_j | P_\mu) V_{\mu\nu}^-0.5
    !
    !    Uses 2d distributed matrices internally.  The result should be
    !    identical to get_coeff_3fn_v_1d(), but this routine should be faster.
    !    It depends on scalapack, though.
    !
    !  USES 

    use timing
    use dimensions
    use prodbas
    use runtime_choices
    use species_data
    use sbt_overlap_aims
    use hartree_fock
    use mpi_tasks, only: aims_stop, check_allocation
    use localorb_io, only: localorb_info
    
    implicit none

    !  ARGUMENTS

    !  INPUTS
    !    none
    !  OUTPUTS
    !
    !  AUTHOR
    !    Computational Chemistry University of Potsdam Lukas Gallandi
    !  SOURCE
    !    get_coeff_3fn_v_2d(ovlp_3fn)  -- FHIaims

    real*8 :: t0
    real*8, allocatable :: coulomb_matr(:,:), coulomb_matr_1d(:,:)
    integer :: info
    character(*), parameter :: func = 'get_coeff_3fn_v_2d'

    if(.not.use_scalapack) then
      call aims_stop('Cannot use 2d distribution without scalapack', func)
    endif

   ! Integrate ovlp3fn (T)
    call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
    !for LR
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn )


    ! Integrate V
    call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
    !for LR
    allocate(coulomb_matr(max_row_2d,max_col_2d), stat=info)
    call check_allocation(info, 'coulomb_matr', func)
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    if (use_logsbt) then
       call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .true.)
    else
       allocate(coulomb_matr_1d(n_basbas, n_loc_prodbas), stat=info)
         call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr_1d)
         ! t0 = mpi_wtime()
         call dist_1d_2d(n_basbas, coulomb_matr_1d, ubound(coulomb_matr_1d, 1), &
         &                         coulomb_matr, ubound(coulomb_matr,1))
         ! if(myid==0) print *,'dist_1d_2d coulomb_matr:',mpi_wtime()-t0
         deallocate(coulomb_matr_1d)
    end if
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_coulomb_matr, clock_time_coulomb_matr, tot_time_coulomb_matr, tot_clock_time_coulomb_matr)

    ! V^-0.5
    call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
    ! for LR
    use_hse=.false.
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    if (use_asym_inv_sqrt) then
      ! transposed==.false. because get_v_multi_ovlp3fn_2() does
      !   C := TV^-0.5.
      call asym_inv_sqrt_of_auxmat_scalapack_2d(coulomb_matr, "Coulomb_LR",.false.)
    else
      call power_auxmat_scalapack_2d(coulomb_matr, -0.5d0, "Coulomb_LR")
    end if
    use_hse=.true.
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr, tot_time_inv_coulomb_matr, tot_clock_time_inv_coulomb_matr)

    ! C := T V^-0.5
    call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
    ! for LR
    call get_v_multi_ovlp3fn_2(coulomb_matr, ovlp_3fn)
    call get_times(time_ovlp_multi, clock_time_ovlp_multi, tot_time_ovlp_multi, tot_clock_time_ovlp_multi)

    deallocate(coulomb_matr)

    call localorb_info('')

  end subroutine get_coeff_3fn_v_2d_lrc_pt2
!##########################################################  

end module lrc_pt2
!******
