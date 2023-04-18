!****s* FHI-aims/prepare_corr_energy_calc
!  NAME
!   prepare_corr_energy_calc
!  SYNOPSIS

      subroutine prepare_corr_energy_calc()

!  PURPOSE
!  This subroutine evaluate the necessary matrix elements (3-center overlap,
!  coulomb matrix) used for accurate correlation energy calculations beyond DFT
!  (e.g. MP2, RPA, RPA+, etc).
!
!  USES

      use dimensions
      use species_data
      use runtime_choices
      use prodbas
      use hartree_fock
      use timing
      use mpi_tasks
      use sbt_overlap_aims
      use localorb_io
      implicit none

!  ARGUMENTS

!  INPUT
!    none
!  OUTPUT
!    none
!
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
      character*20 filename
      logical :: need_o3fn
      real*8 :: time_prodbas_add, clock_time_prodbas_add

! Counter
      integer i_basis_1

      integer :: info
      character*150 :: info_str
      character(*), parameter :: func = 'prepare_corr_energy_calc'

      call get_timestamps(time_prodbas_add, clock_time_prodbas_add)


      ! --- Ensure ovlp_3fn / coeff_3fn_ten&coulomb_matr_lvl

      ! for non-hartree-fock self-consistent calculations, we need to
      ! construct the product (auxiliary) basis functions, and evaluate
      ! the 3-center overlap and coulomb matrix elements here.


      !Reinitialization of ovlp_3fn needed for G0W0@HSE, since in HSE the 
      !overlap include the erf! For now, this is done by setting to zero the HSE 
      !and going through the whole construction of the product basis again.

      !Reinitialization of ovlp_3fn needed for xDH-PBE0@HSE, since in HSE the 
      !overlap include the erf! For now, this is done by setting to zero the HSE 
      !and going through the whole construction of the product basis again.

      if (.not.use_hartree_fock .or. use_gw_and_hse .or. use_dftpt2_and_hse) then 

        call initialize_prodbas()
        
        ! This sub might be called twice, so check need for ovlp_3fn.
        if (sparse_o3fn) then
           if (allocated(coulomb_matr_lvl)) then
              need_o3fn = .false.
           else
              allocate(coulomb_matr_lvl(n_basbas, n_loc_prodbas), stat=info)
              call check_allocation(info, 'coulomb_matr_lvl', func)
              need_o3fn = .true.
           end if
        else
           if (allocated(ovlp_3fn)) then
              need_o3fn = .false.
           else
              allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=info)
              call check_allocation(info, 'ovlp_3fn', func)
              need_o3fn = .true.
           end if
        end if

        !-----------for G0W0@HSE----------------
        if( use_gw_and_hse)then
           deallocate(ovlp_3fn) 
           need_o3fn = .true.
           allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=info)
           call check_allocation(info, 'ovlp_3fn', func)

!           if(myid.eq.0)then
!              write(use_unit,*)"RI logical parameters: "
!              write(use_unit,*)"RI_type    ", RI_type
!              write(use_unit,*)"RI_SVS     ", RI_SVS
!              write(use_unit,*)"RI_V       ", RI_V
!              write(use_unit,*)"RI_LVL_full", RI_LVL_full
!              write(use_unit,*)"RI_LVL     ", RI_LVL
!           endif
        endif 

        !-----------for xDH-PBE0@HSE----------------
        if( use_dftpt2_and_hse)then
           deallocate(ovlp_3fn) 
           need_o3fn = .true.
           allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=info)
           call check_allocation(info, 'ovlp_3fn', func)
        endif 

!test
!        if(myid.eq.0) write(use_unit,*) " Variable need_o3fn : " , need_o3fn
!test end

        if (need_o3fn) then
           select case (RI_type)
           case(RI_SVS)
              call get_coeff_3fn_svs(ovlp_3fn)
           case(RI_V)
              if (use_2d_corr) then
                 call get_coeff_3fn_v_2d(ovlp_3fn)
              else
                 call get_coeff_3fn_v_1d(ovlp_3fn)
              end if
           case(RI_LVL_full)
              call get_coeff_3fn_lvl_full(ovlp_3fn)
           case (RI_LVL)
              call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr_lvl)
           case default
        call aims_stop("Invalid version of RI (resolution of identity)",&
        &              func)
           end select
        endif
     end if
     

     ! --- Optionally swap

     ! if use_ovlp_swap is true, then write ovlp_3fn to disk
     ! and deallocate the array
      if(use_ovlp_swap) then

        if(use_2d_corr) then
           ! The problem is not so much writing; but there is no
           ! ovlp_3KS construction utility for a swapped ovlp_3fn.
           call aims_stop('INTERNAL ERROR: use_ovlp_swap with 2d',func)
        end if

        call get_timestamps(time_ovlp_swap, clock_time_ovlp_swap )
        if(use_mpi) then
          write(filename,'(A,I0)') 'OVLP', myid
        else
          write(filename,'(A)') 'OVLP'
        endif

        open (60,file=filename,form='unformatted', status='unknown', &
                 access ='direct',recl=8*n_basis_pairs)
        do i_basis_1 = 1, n_loc_prodbas, 1
          write(60,rec=i_basis_1) ovlp_3fn(:,i_basis_1)
        enddo

        deallocate (ovlp_3fn)
        call get_times(time_ovlp_swap, clock_time_ovlp_swap )

        close(60)

      endif

      ! --- Timing

      call get_times(time_prodbas_add, clock_time_prodbas_add, &
      &              time_prodbas_total, clock_time_prodbas_total)

      if (.not. use_hartree_fock .or. use_ovlp_swap) then
         call output_timeheader('2X', 'End of correlation preparation', OL_norm)
         if(.not.use_hartree_fock) then
            call ovlp3fn_timings('2X')
         endif
         if(use_ovlp_swap) then
            call output_times('2X', 'Total time for swapping 3-center matrix', &
            & time_ovlp_swap, clock_time_ovlp_swap, OL_norm)
         endif
      end if

      end subroutine prepare_corr_energy_calc
!***************
