!****h* FHI-aims/hartree_fock
!  NAME
!    lc_wpbeh
!  SYNOPSIS

      module lc_wpbeh

!  PURPOSE  
!  this module contains the subroutines for the long range corrected LC-wPBEh functional.  
!
!  Subroutines:
!
!  USES

implicit none

!  AUTHOR
!    Computational Chemistry University of Potsdam Lukas Gallandi
!  SEE ALSO
!    
!  COPYRIGHT
!   
!  HISTORY
!    Introduced October 2014
!
!  SOURCE

!   

     

contains
  subroutine evaluate_2nd_fock_matr_for_lc_wpbeh (hamiltonian )
    !
    !  evaluate the fock_matr for the LC-wPBEh
    !
    use dimensions
    use hartree_fock
    use prodbas
    use species_data
    use runtime_choices
    use mpi_tasks, only: aims_stop
    implicit none

    integer i_spin 
    integer i_state
    integer i_basis_1
    integer i_basis_2
    integer i_index
    real*8  hamiltonian (n_basis*(n_basis+1)/2,n_spin)
    character(*), parameter :: func = 'evaluate_2nd_fock_matr_for_lc_wpbeh'
    if (packed_matrix_format /= PM_none) then
      call aims_stop('Invalid choice of packed_matrix_format for HF', func)
    end if
    if (hybrid_coeff .ne. 0.0d0) then  
		 do i_spin = 1, n_spin
		   i_index = 0
		   do i_basis_1 = 1, n_basis 
		     do i_basis_2 = 1, i_basis_1 
		        i_index = i_index + 1
		        hamiltonian(i_index, i_spin) = &
		        hamiltonian(i_index, i_spin) &
		        - hybrid_coeff * fock_matr_SR( i_basis_2, i_basis_1, i_spin) - (1/lc_dielectric_constant)*fock_matr( i_basis_2, i_basis_1, i_spin)
		     enddo 
		   enddo
		 enddo 
	 else
	 	do i_spin = 1, n_spin
		   i_index = 0
		   do i_basis_1 = 1, n_basis 
		     do i_basis_2 = 1, i_basis_1 
		        i_index = i_index + 1
		        hamiltonian(i_index, i_spin) = &
		        hamiltonian(i_index, i_spin) &
		        - (1/lc_dielectric_constant)*fock_matr( i_basis_2, i_basis_1, i_spin)
		     enddo 
		   enddo
		 enddo 
	 end if        
  end  subroutine evaluate_2nd_fock_matr_for_lc_wpbeh 
!########################################################       
  subroutine get_fock_energy_for_lc_wpbeh &
             (alpha,KS_eigenvector,occ_numbers,fock_energy,energy_xc)

    !  PURPOSE
    !  Subroutine get_fock_energy evaluates the exact exchange energy
    !  for the LC-wPBEh functional.
    !
    !  USES

        use dimensions
        use runtime_choices, only: lc_dielectric_constant
        use hartree_fock

        implicit none

    !  ARGUMENTS

        real*8 alpha
        real*8, dimension(n_basis,n_states,n_spin,n_k_points) ::  &
                KS_eigenvector
        real*8, dimension(n_states,n_spin,n_k_points) :: &
                occ_numbers

        real*8  fock_energy 
        real*8  energy_xc
    !  INPUTS
    !  o  alpha -- the mixing parameter for HF-based calculations, 1.0 for HF,
    !         0.25 for PBE0, etc.
    !  o  KS_eigenvector -- real array,
    !         the eigenvector of the single-particle calculation
    !  o  occ_numbers -- real array,
    !         the occupation number of the electrons for each eigenstate and each spin
    !  OUTPUTS
    !  o  fock_energy -- real number, here the exact exchange energy
    !  o  en_xc -- real number, the exchange-correlation energy, differs form fock_energy
    !          in cases of hybrid functional calculations
    !
    !  SEE ALSO
    !    Subroutine get_fock_energy

    !  local variables
        real*8 energy_fock_LR
        real*8 energy_fock_SR

    !     counters

        integer :: i_state
        integer :: i_basis_1
        integer :: i_basis_2
        integer :: i_index
        integer :: i_spin

    !     begin work

    energy_fock_LR=0.d0
    energy_fock_SR=0.d0
    if (alpha .ne. 0.0d0) then 
		 do i_spin = 1, n_spin
		   do i_state = 1, n_homo(i_spin)
		     do i_basis_1 = 1, n_basis
		        do i_basis_2 = 1, n_basis
		          energy_fock_LR = energy_fock_LR + &
		            fock_matr(i_basis_1,i_basis_2,i_spin) * &
		            KS_eigenvector(i_basis_1, i_state, i_spin, 1) * &
		            KS_eigenvector(i_basis_2, i_state, i_spin, 1) * &
		             occ_numbers(i_state,i_spin,1)
		          energy_fock_SR = energy_fock_SR + &
		             fock_matr_SR(i_basis_1,i_basis_2,i_spin) * &
		            KS_eigenvector(i_basis_1, i_state, i_spin, 1) * &
		            KS_eigenvector(i_basis_2, i_state, i_spin, 1) * &
		            occ_numbers(i_state,i_spin,1)
		      enddo
		     enddo
		   enddo
		 enddo
	 else
	 	do i_spin = 1, n_spin
		   do i_state = 1, n_homo(i_spin)
		     do i_basis_1 = 1, n_basis
		        do i_basis_2 = 1, n_basis
		          energy_fock_LR = energy_fock_LR + &
		            fock_matr(i_basis_1,i_basis_2,i_spin) * &
		            KS_eigenvector(i_basis_1, i_state, i_spin, 1) * &
		            KS_eigenvector(i_basis_2, i_state, i_spin, 1) * &
		             occ_numbers(i_state,i_spin,1)
		      enddo
		     enddo
		   enddo
		 enddo
	 end if

    fock_energy = -0.5d0*(alpha*energy_fock_SR+(1/lc_dielectric_constant)*energy_fock_LR)
    energy_xc = energy_xc + 0.5d0*(alpha*energy_fock_SR + (1/lc_dielectric_constant)*energy_fock_LR)

  end subroutine get_fock_energy_for_lc_wpbeh
!########################################################   
      
 subroutine get_ovlp_3fn_for_lc_wpbeh &
             ()

    !  PURPOSE
    !  Calculate ovlp_3fn matricies for SR, LR coulomb part
    !
    !  AUTHOR
    !  Computational Chemistry University of Potsdam Lukas Gallandi
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

    

    character(*), parameter :: func = 'get_ovlp_3fn_for_lc_wpbeh'
        
    ! In the futur we may also change the other RI-Types but for now
    ! the RI_V method has to be sufficient
    select case (RI_type)
    !   case(RI_SVS)
    !     call get_coeff_3fn_svs(ovlp_3fn)
      case(RI_V)
        if (use_2d_corr) then
          call get_coeff_3fn_v_2d_lc_wpbeh()
        else
          call get_coeff_3fn_v_1d_lc_wpbeh()
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
        

  end subroutine get_ovlp_3fn_for_lc_wpbeh
        
!##################################################
  subroutine get_coeff_3fn_v_1d_lc_wpbeh()

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
    use dimensions
    use prodbas
    use sbt_overlap_aims
    use species_data
    use runtime_choices
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
    !    get_coeff_3fn_v_1d(ovlp_3fn)  -- FHIaims

    real*8, allocatable :: coulomb_matr(:,:)
    real*8, allocatable :: coulomb_matr_SR(:,:)
    integer :: info
    integer :: i, j 
    integer*8 :: array_dimension
    character(*), parameter :: func = 'get_coeff_3fn_v_1d_lc_wpbeh'

    ! Integrate ovlp3fn (T)
    call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
    !for SR
    if (hybrid_coeff .ne. 0.d0) then
	    call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn_SR, ovlp_type_bare_or_hse_coul)
	 end if
    !for LR
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn )


    ! Integrate V
    call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
    !for SR
    if (hybrid_coeff .ne. 0.d0) then
    	allocate(coulomb_matr_SR(n_basbas, n_loc_prodbas), stat=info)
    	call check_allocation(info, 'coulomb_matr_SR', func)
    	if (use_logsbt) then
	   	call integrate_auxmat_by_atomic_sbt(coulomb_matr_SR, ovlp_type_bare_or_hse_coul, .false.)
	 	else
	   	call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr_SR)
	 	end if
    end if
    allocate(coulomb_matr(n_basbas, n_loc_prodbas), stat=info)
    call check_allocation(info, 'coulomb_matr', func)
    !for LR
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    if (use_logsbt) then
       call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .false.)
    else
       call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr)
    end if
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_coulomb_matr, clock_time_coulomb_matr, tot_time_coulomb_matr, tot_clock_time_coulomb_matr)
        
    ! V^-0.5
    call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
    ! for SR
    if (hybrid_coeff .ne. 0.d0) then
		 if (use_asym_inv_sqrt) then
		    ! transposed==.true. because get_v_multi_ovlp3fn() does
		    !   C := V^-0.5T.
		    call asym_inv_sqrt_of_auxmat_scalapack(coulomb_matr_SR, "Coulomb_SR", .true.)
		 else if (use_scalapack) then
		    call power_auxmat_scalapack(coulomb_matr_SR, -0.5d0, "Coulomb_SR")
		 else
		    call power_auxmat_lapack(coulomb_matr_SR, -0.5d0, "Coulomb_SR")
		 end if   
	 end if
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
    ! for SR
    if (hybrid_coeff .ne. 0.d0) then
    	call get_v_multi_ovlp3fn(coulomb_matr_SR, ovlp_3fn_SR)
    end if
    ! for LR
    call get_v_multi_ovlp3fn(coulomb_matr, ovlp_3fn)
    call get_times(time_ovlp_multi, clock_time_ovlp_multi,tot_time_ovlp_multi, tot_clock_time_ovlp_multi)
    
    if (allocated(coulomb_matr_SR)) then
    	deallocate(coulomb_matr_SR)
    end if
    deallocate(coulomb_matr)
     

    call localorb_info('')

  end subroutine get_coeff_3fn_v_1d_lc_wpbeh
!##################################################################
  subroutine get_coeff_3fn_v_2d_lc_wpbeh()

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
    use prodbas
    use dimensions
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
    real*8, allocatable :: coulomb_matr_SR(:,:), coulomb_matr_SR_1d(:,:)
    integer :: info
    character(*), parameter :: func = 'get_coeff_3fn_v_2d'

    if(.not.use_scalapack) then
      call aims_stop('Cannot use 2d distribution without scalapack', func)
    endif

   ! Integrate ovlp3fn (T)
    call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
    !for SR
    if (hybrid_coeff .ne. 0.d0) then
	    call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn_SR, ovlp_type_bare_or_hse_coul)
	 end if
    !for LR
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_LR
    call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
    call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn )


    ! Integrate V
    call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
    !for SR
    if (hybrid_coeff .ne. 0.d0) then
    	 allocate(coulomb_matr_SR(max_row_2d,max_col_2d), stat=info)
       call check_allocation(info, 'coulomb_matr_SR', func)
		 if (use_logsbt) then
		   call integrate_auxmat_by_atomic_sbt(coulomb_matr_SR, ovlp_type_bare_or_hse_coul, .true.)
		 else
		   allocate(coulomb_matr_SR_1d(n_basbas, n_loc_prodbas), stat=info)
		   call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr_SR_1d)
		   ! t0 = mpi_wtime()
		   call dist_1d_2d(n_basbas, coulomb_matr_SR_1d, ubound(coulomb_matr_SR_1d, 1), &
		   &                         coulomb_matr_SR, ubound(coulomb_matr_SR,1))
		   ! if(myid==0) print *,'dist_1d_2d coulomb_matr:',mpi_wtime()-t0
		   deallocate(coulomb_matr_SR_1d)
		 end if
	 end if
	 allocate(coulomb_matr(max_row_2d,max_col_2d), stat=info)
    call check_allocation(info, 'coulomb_matr', func)
    !for LR
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
    !for SR
    if (hybrid_coeff .ne. 0.d0) then
		 if (use_asym_inv_sqrt) then
		   ! transposed==.false. because get_v_multi_ovlp3fn_2() does
		   !   C := TV^-0.5.
		   call asym_inv_sqrt_of_auxmat_scalapack_2d(coulomb_matr_SR, "Coulomb_SR",.false.)
		 else
		   call power_auxmat_scalapack_2d(coulomb_matr_SR, -0.5d0, "Coulomb_SR")
		 end if
	 end if
    ! for LR
    use_hse=.false.
    ovlp_type_bare_or_hse_coul=OVLP_TYPE_COULOMB
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
    ! for SR
    if (hybrid_coeff .ne. 0.d0) then
	    call get_v_multi_ovlp3fn_2(coulomb_matr_SR, ovlp_3fn_SR)
	 end if
    ! for LR
    call get_v_multi_ovlp3fn_2(coulomb_matr, ovlp_3fn)
    call get_times(time_ovlp_multi, clock_time_ovlp_multi, tot_time_ovlp_multi, tot_clock_time_ovlp_multi)

    if (allocated(coulomb_matr_SR)) then
    	deallocate(coulomb_matr_SR)
    end if
    deallocate(coulomb_matr)

    call localorb_info('')

  end subroutine get_coeff_3fn_v_2d_lc_wpbeh
!##########################################################  

end module lc_wpbeh
!******
