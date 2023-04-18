module get_dielectric_function
!  PURPOSE
!  AUTHOR
!  HISTORY
!    Development version, FHI-aims (2010).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  SOURCE
!  use localorb_io
!  use mpi_tasks
  implicit none
  save

  real*8 :: gk_homo_level
  real*8 :: gk_lumo_level
! homo and lumo level variables are set in output_eigenfunctions.f90
 contains
!
! get_sparse_matrix_greenwood is the wrapper around the Kubo-Greenwood 
! electronic transport spectra (sigma, Seebeck) using the memory-saving sparse matrix momemtum matrices
! from module calculate_mommat_base
!
subroutine get_sparse_matrix_greenwood &
     (KS_eigen, KS_vec, KS_vec_complex, occ_numbers,chemical_potential, &
      partition_tab, l_shell_max, ep1_in, ep2_in)

  use calculate_mommat_base
  use dimensions
  use runtime_choices
  use localorb_io
  use mpi_utilities
  use synchronize_mpi_basic, only: sync_vector, sync_real_number
  use pbc_lists
  use geometry, only: cell_volume
  implicit none
!  ARGUMENTS

  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  &
                                                                 KS_vec_complex
  real*8, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  KS_vec
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  real*8, INTENT(IN) :: chemical_potential
  real*8, target, dimension(n_full_points), INTENT(IN) :: partition_tab
  integer, dimension(n_species), INTENT(IN) :: l_shell_max 
  CHARACTER(len=1), INTENT(IN):: ep1_in
  CHARACTER(len=1), INTENT(IN):: ep2_in


  real*8, dimension(:), allocatable :: dielectric_function
  real*8, dimension(:), allocatable :: seebeck
  real*8, dimension(:), allocatable :: abtew_cond
  real*8, dimension(:), allocatable :: abtew_seebeck
  real*8, dimension(:), allocatable :: fermideriv
  real*8:: omegapl
  real*8:: fermilevel
  character*128 :: info_str
  integer :: info
  integer :: n_state_min_in
  integer :: n_state_max_in
  logical :: use_absorption

  integer :: i_k
  integer :: new_k_point
  integer :: number_directions
  integer :: count1, coord1, coord2
  integer :: fermiindex, fermisteps
  real*8 :: abtewintegral_cond
  real*8 :: abtewintegral_seeb
  integer :: fermiset_revertflag

  fermiset_revertflag = 0
  !  begin work
    use_absorption=.False.
    write(info_str,'(6X,A,1X,I4)') "Kubo-Greenwood transport post processing starts"
    call localorb_info ( info_str )


    if (.not.allocated( dielectric_function)) then
      allocate( dielectric_function ( n_omega),stat=info)
      call check_allocation(info, 'dielectric_function')
    end if

    if (.not.allocated(seebeck)) then
      allocate(seebeck ( n_omega),stat=info)
      call check_allocation(info, 'seebeck')
    end if

    if (.not.allocated(abtew_cond)) then
      allocate( abtew_cond ( n_greenenergy),stat=info)
      call check_allocation(info, 'abtew_cond')
    end if

    if (.not.allocated( abtew_seebeck)) then
      allocate( abtew_seebeck ( n_greenenergy),stat=info)
      call check_allocation(info, 'abtew_seebeck')
    end if

    if (.not.allocated( fermideriv)) then
      allocate( fermideriv ( n_greenenergy),stat=info)
      call check_allocation(info, 'fermideriv')
    end if

    call get_state_minmax_K(KS_eigen, n_state_min_in, n_state_max_in)

fermisteps=fermistatesabove+fermistatesbelow
do fermiindex = 0, fermisteps
   fermilevel = chemical_potential-((fermispacing/hartree)*fermistatesbelow)+(fermiindex*(fermispacing/hartree))


  if(flag_explicit_fermi) then  
     select case(percent_or_ev)
         case(1)
            fermilevel= gk_homo_level + dist_from_vbm/hartree  ! an absolute value above VBM 
             if (fermilevel>gk_lumo_level) then 
                fermilevel=chemical_potential
                fermiset_revertflag=1
             endif
         case(2)
            fermilevel= gk_homo_level + dist_from_vbm*(gk_lumo_level-gk_homo_level)  ! above VBM in percent of the gap
    end select
  endif

    omegapl=0.0d0
    dielectric_function=0.0d0
    seebeck=0.0d0
    abtew_cond=0.0d0
    abtew_seebeck=0.0d0
    fermideriv=0.0d0

    if (flag_out_dclimit) then
      call calc_fermideriv(fermideriv, fermilevel)
    endif 

if (ep1_in==ep2_in) then

    if (ep1_in=='a') then 
    number_directions=3 
    else 
    number_directions=1
    endif

  do count1 = 1, number_directions

  select case (ep1) 
      case('x') 
        coord1=1
      case('y') 
        coord1=2
      case('z') 
        coord1=3
      case('a') 
        coord1=count1 
   endselect

      call allocate_mommat()
      call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up, mommat_full_oned_low,coord1)
      i_k = 0
      do new_k_point = 1,n_k_points
		if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points ) then
		i_k = i_k+1 !new_k_point
       	call allocate_mommat_k()
     	call allocate_moment_one(n_state_min_in, n_state_max_in)


		 call construct_overlap( mommat_full_oned_up, mommat_full_w_up,&
                                        mommat_full_w_complex_up, new_k_point,&
                                        work_ovl_mom )
		call construct_overlap( mommat_full_oned_low, &
                     mommat_full_w_low, mommat_full_w_complex_low, new_k_point,&
                     work_ovl_mom )
		call calc_moment_p0(moment_one,mommat_full_w_up, &
                     mommat_full_w_low, mommat_full_w_complex_up, &
                     mommat_full_w_complex_low, KS_vec(:,:,:,i_k) ,&
                     KS_vec_complex(:,:,:,i_k) , new_k_point,coord1, &
                     n_state_min_in, n_state_max_in) ! coord needs not to be passed here (it's not used) - remains for consistency. 

		call calc_greenwood(moment_one,moment_one,&
                     dielectric_function,seebeck, abtew_cond, abtew_seebeck,&
                     fermideriv ,omegapl,fermilevel,&
                     KS_eigen(:,:,new_k_point),k_weights(new_k_point),&
                     widthone, &
                     widthtwo, n_state_min_in, n_state_max_in)

        call clean_mommat()
		endif ! k-point-parallism condition
      enddo  ! k-point loop
   enddo ! x-y-z loop end (1 or 3 passages according to number_directions) 


dielectric_function=dielectric_function/number_directions
seebeck=seebeck/number_directions
abtew_cond=abtew_cond/number_directions
abtew_seebeck=abtew_seebeck/number_directions ! for "x-y-z-averaged" spectra, the latter are divided by 3;  otherwise by 1...

  else ! distinction between xx,yy,zz,aa OR xy, yz, xz transport properties

      call allocate_mommat_two()
  select case(ep1) 
      case('x') 
coord1=1
      case('y') 
coord1=2
      case('z') 
coord1=3
   endselect

  select case(ep2) 
      case('x') 
coord2=1
      case('y') 
coord2=2
      case('z') 
coord2=3
   endselect

      call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up, mommat_full_oned_low,coord1)
      call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up, mommat_full_oned_low,coord2)

      i_k = 0
      do new_k_point = 1,n_k_points
		if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points ) then
		i_k = i_k+1 !new_k_point
        call allocate_mommat_k()
    	call allocate_moment_one(n_state_min_in, n_state_max_in)
	    call allocate_moment_two(n_state_min_in, n_state_max_in)

		call construct_overlap( mommat_full_oned_up, mommat_full_w_up,&
                     mommat_full_w_complex_up, new_k_point, work_ovl_mom )
		call construct_overlap( mommat_full_oned_low, &
                     mommat_full_w_low, mommat_full_w_complex_low, new_k_point,&
                     work_ovl_mom )
		call calc_moment_p0(moment_one,mommat_full_w_up, &
                     mommat_full_w_low, mommat_full_w_complex_up, &
                     mommat_full_w_complex_low, KS_vec(:,:,:,i_k) ,&
                     KS_vec_complex(:,:,:,i_k) , new_k_point,coord1, &
                     n_state_min_in, n_state_max_in) ! coord not necessary
		call construct_overlap( mommat_full_oned_two_up, &
                     mommat_full_w_up, mommat_full_w_complex_up, new_k_point, &
                     work_ovl_mom )
		call construct_overlap( mommat_full_oned_two_low, &
                     mommat_full_w_low, mommat_full_w_complex_low, new_k_point,&
                     work_ovl_mom )
		call calc_moment_p0(moment_two,mommat_full_w_up, &
                     mommat_full_w_low, mommat_full_w_complex_up, &
                     mommat_full_w_complex_low, KS_vec(:,:,:,i_k) ,&
                     KS_vec_complex(:,:,:,i_k) , new_k_point, coord1, &
                     n_state_min_in, n_state_max_in) ! coord not necessary
		call calc_greenwood(moment_one,moment_two,&
                     dielectric_function, seebeck, abtew_cond, abtew_seebeck,&
                     fermideriv, omegapl,fermilevel,&
                     KS_eigen(:,:,new_k_point),k_weights(new_k_point),&
                     widthone, &
                     widthtwo, n_state_min_in, n_state_max_in)

       call clean_mommat()
		endif ! k-point-parallism condition
      enddo  ! k-point loop

 endif ! ep1!=/==ep2 condition
 
        if(.not. use_local_index) call sync_vector(dielectric_function, n_omega )  
	    if(.not. use_local_index) call sync_vector(seebeck, n_omega )
        if(.not. use_local_index) call sync_real_number(omegapl)

         if(flag_out_dclimit)then 
    		 if(.not. use_local_index) call sync_vector(abtew_cond, n_greenenergy)
		     if(.not. use_local_index) call sync_vector(abtew_seebeck, n_greenenergy)
         endif  

           if(myid == 0) then
                
                call out_greenwood(dielectric_function, seebeck, abtew_cond, &
                                   abtew_seebeck, fermideriv, ep1, ep2, &
                                   cell_volume, fermilevel)
                write(use_unit,'(2X, A)') &
                   'Plasma-frequency computed with sparse momentum-matrix'
                write(use_unit,'(2X, A, 1X, ES12.4,1X, 6A, ES14.6, A)') &
                   'Fermi-level:',fermilevel*hartree,' eV:', &
                   '  Plasmafrequency ',ep1, ', ' ,ep2,' : ', &
                   sqrt(omegapl)*hartree,' eV'

               if(flag_out_dclimit)then 
                 abtewintegral_cond= ergebnis(abtew_cond, abtew_seebeck, cell_volume, 1)
                 abtewintegral_seeb= ergebnis(abtew_cond, abtew_seebeck, cell_volume, 2)

                write(use_unit,'(2X, A)') 'Kubo-Greenwood "omega-->0" values computed with sparse momentum-matrix'

if(fermiindex==fermistatesbelow .and. .not.flag_explicit_fermi  ) then

  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,' eV:',&
'  Integrated L11(Cond): ', abtewintegral_cond, ' (Ohm-1*cm-1) - intrinsic Fermi level'
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,' eV:',&
'  Integrated Seeb:      ', 1E6*0.01*abtewintegral_seeb/(11604.519*widthtwo*abtewintegral_cond),&
 ' (muV/K) - intrinsic Fermi level' 
! 11604.519 is eV--> Kelvin for widthtwo which is a Temperature; 
! 0.01 to correct for (Ohm*cm)-1 unit sigma comes in; 1E6 for muV/K
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,&
' eV:','  Integrated L12:       ', abtewintegral_seeb, ' (arb.units) - intrinsic Fermi level' 


else if(fermiindex==fermistatesbelow .and. .not.flag_explicit_fermi ) then

!  write(use_unit,'(A, 1X, ES12.4,1X, 2A, ES14.6, A)') 'TestOutput for Lumo:',gk_lumo_level*hartree,' &
!  eV:','  TestOutput for Homo ',gk_homo_level*hartree,' eV' 
! test output of homo/lumo

  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,&
' eV:','  Integrated L11(Cond): ', abtewintegral_cond, ' (Ohm-1*cm-1) - explicit Fermi level'
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,&
' eV:','  Integrated Seeb:      ', 1E6*0.01*abtewintegral_seeb/(11604.519*widthtwo*abtewintegral_cond),&
 ' (muV/K) - explicitt Fermi level' 
! 11604.519 is eV--> Kelvin for widthtwo which is a Temperature; 
!0.01 to correct for (Ohm*cm)-1 unit sigma comes in; 1E6 for muV/K
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',&
fermilevel*hartree,' eV:','  Integrated L12:       ', abtewintegral_seeb, ' (arb.units) - explicit Fermi level' 

else 


  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,&
' eV:','  Integrated L11(Cond): ', abtewintegral_cond, ' (Ohm-1*cm-1)'
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,&
' eV:','  Integrated Seeb:      ', 1E6*0.01*abtewintegral_seeb/(11604.519*widthtwo*abtewintegral_cond), ' (muV/K)'
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,&
' eV:','  Integrated L12:       ', abtewintegral_seeb, ' (arb.units)' 


endif  ! output condition to underline intrinsic or artifical fermi level
endif  ! end of abtew-condition dclimit=true
endif  ! if myid==0 condition         
   
enddo ! end of artifical Fermi-level loop

    call clean_mommat_final()
   if (allocated(dielectric_function))  deallocate(dielectric_function)
   if (allocated(seebeck))  deallocate(seebeck)
   if (allocated(abtew_cond))  deallocate(abtew_cond)
   if (allocated(abtew_seebeck))  deallocate(abtew_seebeck)
   if (allocated(fermideriv))  deallocate(fermideriv)
  write(info_str,'(6X,A,1X,I4)') "Momentum Matrix post processing finished"
  call localorb_info (info_str)

end subroutine get_sparse_matrix_greenwood

!
! get_full_matrix_greenwood calcultes the Kubo-Greenwood electronic transport 
! by using the full-matrix momentum matrices from module calculate_dipolemat.f90 - remains for consistency checks
!
subroutine get_full_matrix_greenwood(KS_eigen, KS_vec, KS_vec_complex, occ_numbers, &
           chemical_potential, partition_tab, l_shell_max, ep1_in, ep2_in)
      use calculate_dipolemat
      use dimensions
      use localorb_io
      use mpi_utilities
      use runtime_choices
      use synchronize_mpi_basic, only: sync_vector, sync_real_number
      use pbc_lists
      use geometry, only: cell_volume

  implicit none

  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  &
                                                                 KS_vec_complex
  real*8, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  KS_vec
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  real*8, INTENT(IN) :: chemical_potential
  real*8, target, dimension(n_full_points), INTENT(IN) :: partition_tab
  integer, dimension(n_species), INTENT(IN) :: l_shell_max 
  CHARACTER(len=1), INTENT(IN):: ep1_in
  CHARACTER(len=1), INTENT(IN):: ep2_in


      character*200 :: info_str

      integer :: i_k,new_k_point
      integer :: count1
      integer :: number_directions
      integer :: fermisteps, fermiindex
      real*8 :: fermilevel

          call get_state_minmax(KS_eigen)
          call allocate_spectra()

fermisteps=fermistatesabove+fermistatesbelow
do fermiindex = 0, fermisteps
   fermilevel = chemical_potential-((fermispacing/hartree)*fermistatesbelow)+(fermiindex*(fermispacing/hartree))

    omegapl=0.0d0
    die_el=0.0d0
    seebeck=0.0d0
    abtew=0.0d0
    abtewseebeck=0.0d0
    fermideriv=0.0d0

    if (flag_out_dclimit) then
      call calc_fermideriv(fermideriv, fermilevel)
    endif 

if (ep1_in==ep2_in) then

    if (ep1_in=='a') then 
    number_directions=3 
    else 
    number_directions=1
    endif

  do count1 = 1, number_directions
     call allocate_dipole_mat() ! allocates dipole_mat_full_oned
 
         call calculate_dipolemat_p0 ( partition_tab, l_shell_max, dipole_mat_full_oned, ep1_in, count1) 
      i_k = 0
      do new_k_point = 1,n_k_points
		if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points ) then
		i_k = i_k+1 !new_k_point
		call allocate_dipole_mat_k()
 		call construct_dipolemat_p1(dipole_mat_full_oned, dipolemat_full_w, dipolemat_full_w_complex, k_phase(:,new_k_point) )
		call calc_dipelement_p0 (dipelement_one, dipolemat_full_w, &
                               dipolemat_full_w_complex, KS_vec(:,:,:,i_k), &
                               KS_vec_complex(:,:,:,i_k) , new_k_point)
		call calc_diel(dipelement_one,dipelement_one,die_el,seebeck,abtew,&
                       abtewseebeck,fermideriv,omegapl,fermilevel,KS_eigen(:,:,new_k_point),KS_eigen(:,:,:),k_weights(new_k_point))
        call clean_dipole_mat()
		endif ! k-point-parallism condition
      enddo  ! k-point loop
   enddo ! x-y-z loop end (1 or 3 passages according to number_directions) 


die_el=die_el/number_directions
seebeck=seebeck/number_directions
abtew=abtew/number_directions
abtewseebeck=abtewseebeck/number_directions ! for "x-y-z-averaged" spectra, the latter are divided by 3;  otherwise by 1...


  else ! distinction between xx,yy,zz,aa OR xy, yz, xz transport properties 

      call allocate_dipole_mat_two() ! allocates dipole_mat_full_oned and dipole_mat_full_oned_two
      call calculate_dipolemat_p0 ( partition_tab, l_shell_max, dipole_mat_full_oned, ep1, 1) ! '1' in the last slot serves as a dummy 
      call calculate_dipolemat_p0 ( partition_tab, l_shell_max, dipole_mat_full_oned_two, ep2, 1)

      i_k = 0
      do new_k_point = 1,n_k_points
		if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points ) then
		i_k = i_k+1 !new_k_point
       call allocate_dipole_mat_k()
	    call construct_dipolemat_p1 (dipole_mat_full_oned, dipolemat_full_w, &
                                    dipolemat_full_w_complex, &
                                    k_phase(:,new_k_point))
	    call calc_dipelement_p0 (dipelement_one, dipolemat_full_w, &
                                dipolemat_full_w_complex, KS_vec(:,:,:,i_k), &
                                KS_vec_complex(:,:,:,i_k) , new_k_point)
		call construct_dipolemat_p1(dipole_mat_full_oned_two, dipolemat_full_w, &
                                  dipolemat_full_w_complex, &
                                  k_phase(:,new_k_point))
		call calc_dipelement_p0 (dipelement_two, dipolemat_full_w, &
                               dipolemat_full_w_complex, KS_vec(:,:,:,i_k), &
                               KS_vec_complex(:,:,:,i_k), new_k_point)
		call calc_diel (dipelement_one, dipelement_two, die_el, seebeck, &
                      abtew, abtewseebeck, fermideriv, omegapl, &
                      fermilevel,KS_eigen(:,:,new_k_point),KS_eigen(:,:,:),k_weights(new_k_point))
		endif ! k-point-parallism condition
      enddo  ! k-point loop

 endif ! ep1!=/==ep2 condition
 
        if(.not. use_local_index) call sync_vector(die_el, n_omega )  
	    if(.not. use_local_index) call sync_vector(seebeck, n_omega )
        if(.not. use_local_index) call sync_real_number(omegapl)

         if(flag_out_dclimit)then 
    		 if(.not. use_local_index) call sync_vector(abtew, n_greenenergy)
		     if(.not. use_local_index) call sync_vector(abtewseebeck, n_greenenergy)
         endif  

            if(myid == 0) then
                call out_die_el (die_el, seebeck,abtew, abtewseebeck, &
                                 fermideriv, ep1,ep2, cell_volume, &
                                 fermilevel)  ! General output
                write(use_unit,'(2X, A)') &
                   'Plasma-frequency computed with full momentum-matrix'
                write(use_unit,'(2X, A, 1X, ES12.4,1X, 6A, ES14.6, A)') &
                   'Fermi-level:', fermilevel*hartree, ' eV:', &
                   '  Plasmafrequency ',ep1, ', ' ,ep2,' : ', &
                   sqrt(omegapl)*hartree,' eV'

               if (flag_out_dclimit) then 
                 abtewintegral_cond = ergebnis(abtew, abtewseebeck, &
                                               cell_volume, 1)
                 abtewintegral_seeb = ergebnis(abtew, abtewseebeck, &
                                               cell_volume, 2)

                write(use_unit,'(2X, A)') 'Kubo-Greenwood "omega-->0" values computed with full momentum-matrix'

if(fermiindex==fermistatesbelow) then

  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,' eV:',&
  '  Integrated L11(Cond): ', abtewintegral_cond, ' (Ohm-1*cm-1) - intrinsic Fermi level'
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') &
     'Fermi-level:',fermilevel*hartree,' eV:', &
     '  Integrated Seeb:      ', &
     1E6*0.01*abtewintegral_seeb/(11604.519*widthtwo*abtewintegral_cond), &
     ' (muV/K) - intrinsic Fermi level' 
! 11604.519 is eV--> Kelvin for widthtwo which is a Temperature; 
! 0.01 to correct for (Ohm*cm)-1 unit sigma comes in; 1E6 for muV/K
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') &
     'Fermi-level:',fermilevel*hartree,&
     ' eV:','  Integrated L12:       ', abtewintegral_seeb, &
     ' (arb.units) - intrinsic Fermi level' 


else

  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,' eV:',&
  '  Integrated L11(Cond): ', abtewintegral_cond, ' (Ohm-1*cm-1)'
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,' eV:',&
  '  Integrated Seeb:      ', 1E6*0.01*abtewintegral_seeb/(11604.519*widthtwo*abtewintegral_cond), ' (muV/K)'
  write(use_unit,'(2X, A, 1X, ES12.4, 1X, 2A, ES14.6, A)') 'Fermi-level:',fermilevel*hartree,' eV:',&
  '  Integrated L12:       ', abtewintegral_seeb, ' (arb.units)' 


endif  ! output condition to underline intrinsic or artifical fermi level



               endif ! end of abtew-condition dclimit=true
            endif  ! if myid==0 condition         
   
enddo ! end of artifical Fermi-level loop

            call clean_dipole_mat_final()
            write(info_str,'(6X,A,1X,I4)') "Momentum Matrix post processing finished"
            call localorb_info ( info_str )

end subroutine get_full_matrix_greenwood

end module get_dielectric_function
