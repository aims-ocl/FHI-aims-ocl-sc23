!****s* FHI-aims/initialize_hartree_fock_p0
!  NAME
!   initialize_hartree_fock_p0
!  SYNOPSIS

      subroutine initialize_hartree_fock_p0 &
             ()

!  PURPOSE  
!  prepare all the necessary ingredients for Hartree-Fock and post-Hartree-Fock
!  calcualtions. These includes
!  * construct the auxiliary basis, 
!  * calculate the 3-center O integrals 
!  * calcualte the 2-center Coulomb interacton integrals
!  * multiply O to the square root of the Coulomb integral.

!  USES
      use dimensions
      use runtime_choices
      use prodbas
      use physics
      use hartree_fock_p0
      use timing
      use synchronize_mpi, only: sync_timing
      use species_data
      use calculate_fock_matrix_p0, only: cleanup_fock_matrix_calculations, &
          init_fock_matrix_calculations, evaluate_exchange_matr_realspace_p0
      use localorb_io, only: use_unit
      use mpi_tasks, only: myid
      implicit none

!  INPUTS
!  none
!  OUTPUTS
!  none
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
      real*8 exx_ene, d_exx_ene(3,n_atoms)
      character*40 :: exx_restart
      logical :: exx_restartfile_exists = .false.

! First, clean HF quantities which depend on pairs (might change with
! relaxation)

      ! Total time for product basis setup initialized
      call get_timestamps(time_prodbas_total, clock_time_prodbas_total)

      call cleanup_basbas()

      if (flag_auxil_basis .eq. PRODBAS_FULL) then
         call shrink_full_auxil_basis_v2()
      else
         write(use_unit,*) " The given choice of the auxiliary basis is NOT ",&
         "allowed with RI_method LVL_fast. STOP! "
         stop
      endif

      call allocate_hartree_fock_p0

      ! Initialize fock matrix calculations
      ! Before doing that, call cleanup_fock_matrix_calculations for the case that
      ! the call is from reinitialize_scf
      ! The cleanup routine may be safely called even without previous initialization

      call cleanup_fock_matrix_calculations
      call init_fock_matrix_calculations(use_forces, use_analytical_stress)

      if(restart_file_exists) then ! SVL: restart
         write(exx_restart,'(A,A)') 'exx_',trim(restart_read_file)
         inquire(file=exx_restart,exist=exx_restartfile_exists)
      endif

      if(exx_restartfile_exists) then
         if(myid == 0) then
            write(use_unit,*) "Reading the Hartree-Fock exchange matrix from disk"
         endif

         open(file=exx_restart,unit=7,status='old',form='unformatted')

         if(real_eigenvectors) then
            read(7) hf_exchange_matr_real
         else
            read(7) hf_exchange_matr_complex
         endif

         close(7)
      else ! No restart
         call evaluate_exchange_matr_realspace_p0(KS_eigenvector,&
                 KS_eigenvector_complex,occ_numbers,exx_ene,d_exx_ene,.false.,&
                 .false.)
      endif

      ! Total time for product basis finalized
      call get_times(time_prodbas_total,clock_time_prodbas_total,&
              tot_time_prodbas_total,tot_clock_time_prodbas_total)

    end subroutine initialize_hartree_fock_p0

!******
