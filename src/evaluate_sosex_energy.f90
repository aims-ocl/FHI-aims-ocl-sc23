!****s* FHI-aims/evaluate_sosex_energy
!  NAME
!   evaluate_sosex_energy
!  SYNOPSIS

      subroutine evaluate_sosex_energy &
       (n_low_state, n_KS_states, n_homo,&
        n_freq, omega_full, womega_full,&
        occ_numbers, KS_eigenvalue, ovlp_3KS, e_sosex )

!  PURPOSE
!  Calculate the screended second order exchange energy.
!
!  Note: It amounts a resummation of one type of exchange
!  type diagrams to infinite order.

!  USES

      use constants
      use dimensions
      use prodbas
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use runtime_choices
      use constants
      use evaluate_polarisability_freq, only: evaluate_polarisability_freq_0
      use localorb_io, only: use_unit

      implicit none

      integer :: n_homo(n_spin)
      integer :: n_low_state
      integer :: n_KS_states
      integer :: n_freq
      real*8  :: omega_full(n_freq)
      real*8  :: womega_full(n_freq)
      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: ovlp_3KS(n_loc_prodbas,n_states,n_KS_states,n_spin)
      real*8  :: e_sosex

! INPUTS
! o  n_low_state -- the first KS state to be included in frozen-core
!                   calculation
! o  n_homo -- integer array, the HOMO level for each spin channel
! o  n_KS_states -- integer number,
! o          the highest KS/HF eigenstate, the same as the n_high_state.
! o          In the present case, n_high_state >= n_homo should do fine.
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  ovlp_3KS -- real array
!            this is the transformed 3-cener overlap integration. Now two orbitals of
!            them are KS ones, and one is the auxiliary basis.
!            Note: for parallel calculations, the auxiliary basis are distribuated
!            among the different processors.
! OUTPUT
! o  e_sosex -- real number, the calculated second-order exchange energy
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

      real*8 :: aux

      real*8, dimension(:,:), allocatable :: polar_freq
      real*8, dimension(:,:), allocatable :: polar_lambda
      real*8, dimension(:,:,:), allocatable :: dielec_func
      real*8, dimension(:,:), allocatable :: aux_eri
      real*8, dimension(:), allocatable :: aux_eri_real
      real*8, dimension(:,:), allocatable :: aux_ovlp3KS
      real*8, dimension(:,:), allocatable :: tmp_ovlp3KS
      real*8, dimension(:), allocatable :: aux_vect
      real*8, dimension(:), allocatable :: aux_vect_2
      real*8, dimension(:), allocatable :: aux_spect

      real*8 :: e_diff

      complex*16 :: imag_unit
      parameter (imag_unit = cmplx(0.d0, 1.d0))

      integer ::  n_lumo(n_spin)
      integer ::  n_homo_max
      integer ::  n_lumo_min
      integer ::  n_unocc
      integer ::  n_occ


      integer :: n_lambda
      parameter (n_lambda=5)
      real*8   lambda_grid(n_lambda)
      real*8   lambda_weight(n_lambda)

!     Integer counters
      integer i_spin
      integer i_spin2
      integer i_state
      integer i_state_2
      integer i_state_3
      integer i_state_4
      integer i_task
      integer i_index
      integer i_start_b
      integer i_freq
      integer i_freq_2
      integer i_prodbas_1
      integer i_basbas
      integer :: i_lambda

!     Error variable
      integer mpierr

!     external functions
      real*8, external ::  ddot


      if (myid.eq.0) then
          write(use_unit,'(A)') " -------------------------------------------------------------------------" 
          write(use_unit,'(2A)') " | Evaluating 2nd order screened", &
                " exchange (SOSEX) energy (1D distribution)..."
          write(use_unit,'(A)') 
      endif
!      if (n_KS_states .ne. n_states) then
!         if( myid.eq.0 ) then
!           write(use_unit,'(2X,3A)') "* Evaluating SOSEX energy: ", &
!             "'state_upper_limit (n_KS_states)' must be set", &
!             " to n_states!"
!         endif
!         stop
!      endif

      n_lumo(:) = n_homo(:) + 1
   
      n_homo_max = max(n_homo(1),n_homo(n_spin))
      n_lumo_min = min(n_lumo(1),n_lumo(n_spin))

      n_unocc = n_states-n_lumo_min+1
      n_occ = n_homo_max-n_low_state+1



      if (.not.allocated(dielec_func)) then
          allocate(dielec_func(n_basbas,n_loc_prodbas,n_freq))
      endif
      if (.not.allocated(polar_freq)) then
          allocate(polar_freq(n_basbas,n_loc_prodbas))
      endif
      if (.not.allocated(polar_lambda)) then
          allocate(polar_lambda(n_basbas,n_loc_prodbas))
      endif

      e_sosex = 0.d0
! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

! Gauss-Legendre grid for the coupling constant integral
      call gauleg(0d0,1.d0,lambda_grid,lambda_weight,n_lambda)

      dielec_func = 0.d0
      e_sosex = 0.d0
      do i_freq = 1, n_freq, 1

! compute the symmetrized polarizability (v^(1/2) * \chi_0 * v^(1/2) )
        call evaluate_polarisability_freq_0 &
           ( n_low_state, n_homo, n_lumo, n_KS_states, &
             occ_numbers, omega_full(i_freq), &
             KS_eigenvalue, ovlp_3KS, polar_freq &
           )

        polar_freq(:,:) = polar_freq(:,:) * 2.d0/dble(n_spin)

        do i_lambda = 1, n_lambda, 1
            polar_lambda(:,:) = -polar_freq(:,:) * lambda_grid(i_lambda)

            do i_prodbas_1 = 1, n_loc_prodbas, 1
               i_index = map_prodbas(i_prodbas_1, myid+1) 
               if(i_index .gt. 0 ) then
                polar_lambda(i_index,i_prodbas_1) = 1.0d0 + &
                 polar_lambda(i_index,i_prodbas_1)
               endif
            enddo
            if(use_scalapack) then
              call power_auxmat_scalapack(polar_lambda,-1.d0," ")
            else
              call power_auxmat_lapack(polar_lambda,-1.d0," ")
            endif

            dielec_func(:,:,i_freq) = dielec_func(:,:,i_freq) + &
                 lambda_grid(i_lambda)*polar_lambda(:,:)*lambda_weight(i_lambda)
        enddo

      enddo

      dielec_func(:,:,:) =  dielec_func(:,:,:)*2.d0

      if (allocated(polar_freq)) then
          deallocate(polar_freq)
      endif
      if (allocated(polar_lambda)) then
          deallocate(polar_lambda)
      endif

      if (.not.allocated(aux_eri)) then
          allocate(aux_eri(n_unocc,n_occ))
      endif
      if (.not.allocated(aux_eri_real)) then
          allocate(aux_eri_real(n_unocc*n_occ))
      endif
      if (.not.allocated(aux_ovlp3KS)) then
          allocate(aux_ovlp3KS(n_unocc*n_occ,n_loc_prodbas))
      endif
      if (.not.allocated(tmp_ovlp3KS)) then
          allocate(tmp_ovlp3KS(n_loc_prodbas,n_occ))
      endif
      if (.not.allocated(aux_vect)) then
          allocate(aux_vect(n_loc_prodbas))
      endif
      if (.not.allocated(aux_vect_2)) then
          allocate(aux_vect_2(n_basbas))
      endif
      if (.not.allocated(aux_spect)) then
          allocate(aux_spect(n_freq))
      endif

      do i_spin = 1, n_spin, 1

         i_index = 0 
         aux_ovlp3KS(:,:)=0.d0
         do i_state = n_low_state, n_homo(i_spin), 1
           do i_state_2 = n_lumo(i_spin), n_states, 1
              i_index = i_index + 1
              aux_ovlp3KS(i_index,1:n_loc_prodbas) = &
              ovlp_3KS(1:n_loc_prodbas, i_state_2, i_state, i_spin)
           enddo
         enddo

         do i_state_2 = n_lumo(i_spin), n_states, 1

           if(myid.eq.0 .and. mod(i_state_2,10).eq.0 ) then
               write(use_unit,'(A,2I6)') " | spin, virtial orb : ", i_spin, i_state_2 
           endif

           tmp_ovlp3KS(:,:) = 0.d0
           do i_state_3 = n_low_state, n_homo(i_spin), 1
            tmp_ovlp3KS(:,i_state_3-n_low_state+1) = ovlp_3KS(:,i_state_2,i_state_3,i_spin)
           enddo

           do i_state = n_low_state, n_homo(i_spin), 1

              ! Igor :: fix the bug according to the illegel case of n_loc_prodbas = 0
              if (n_loc_prodbas.gt.0) then
                  call dgemm('T', 'N', n_unocc, n_occ, &
                          n_loc_prodbas, 1.0d0, &
                          ovlp_3KS(1,n_lumo_min,i_state,i_spin), &
                          n_loc_prodbas, &
                          tmp_ovlp3KS(1,1), n_loc_prodbas, 0.d0, &
                          aux_eri(1,1), &
                          n_unocc )
              else
                  aux_eri = 0.d0
              endif


              call sync_matrix(aux_eri(1,1),n_unocc,n_occ) 

              do i_freq = 1, n_freq, 1

                aux_eri_real(:) = 0.d0
                i_index = 0
                do i_state_3 = n_low_state, n_homo(i_spin), 1
                  do i_state_4 = n_lumo(i_spin), n_states, 1


                    i_index = i_index + 1

                    e_diff = KS_eigenvalue(i_state_3,i_spin,1) - KS_eigenvalue(i_state_4,i_spin,1)

                    aux_eri_real(i_index)  = &
                       aux_eri(i_state_4-n_lumo_min+1,i_state_3-n_low_state+1) * 2.d0 *   &
                       e_diff/(e_diff*e_diff + omega_full(i_freq)*omega_full(i_freq) )

                  enddo
                enddo

                call dgemv('T', n_occ*n_unocc, n_loc_prodbas, &
                      1.d0, aux_ovlp3KS, n_occ*n_unocc, &
                      aux_eri_real, 1, 0.d0, aux_vect(1), 1)

                call dgemv('N', n_basbas, n_loc_prodbas, 1.d0, &
                     dielec_func(1,1,i_freq), n_basbas, aux_vect(1), &
                     1, 0.d0, aux_vect_2(1), 1)

                call sync_vector(aux_vect_2,n_basbas)

                do i_prodbas_1 = 1, n_loc_prodbas, 1
                  i_index = map_prodbas(i_prodbas_1, myid+1) 
                  aux_vect(i_prodbas_1) = aux_vect_2(i_index) 
                enddo

                aux_spect(i_freq) = &
                 ddot(n_loc_prodbas, ovlp_3KS(1,i_state_2,i_state,i_spin),  &
                      1, aux_vect(1), 1) 

              enddo
              call sync_vector(aux_spect,n_freq)

              e_diff = KS_eigenvalue(i_state,i_spin,1) - KS_eigenvalue(i_state_2,i_spin,1)
              do i_freq = 1, n_freq, 1
                  e_sosex = e_sosex +  &
                    aux_spect(i_freq) * 2.d0*e_diff/ & 
                    (omega_full(i_freq)*omega_full(i_freq)+e_diff*e_diff) *  &
                    occ_numbers(i_state,i_spin,1)* womega_full(i_freq)
              enddo

!              write(use_unit,'(2X,A,2I4,2f16.8)') " SOSEX energy now : ", i_state, i_state_2, &
!                    e_sosex/4.d0/pi
          
!   End loop i_state_2
            enddo

!   End loop i_state
        enddo
!   End loop i_spin
      enddo
      e_sosex = e_sosex/pi/4.d0

      if (allocated(aux_eri)) then
        deallocate (aux_eri)
      endif

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A,f19.8,A,f19.8,A)') &
         " 2nd-order screened exchange energy : ",  &
          e_sosex, '  Ha. ',  e_sosex*hartree,  ' eV '
        write(use_unit,*)
      endif

      if (allocated(dielec_func)) then
          deallocate(dielec_func)
      endif
      if (allocated(aux_eri)) then
          deallocate(aux_eri)
      endif
      if (allocated(aux_eri_real)) then
          deallocate(aux_eri_real)
      endif
      if (allocated(aux_ovlp3ks)) then
          deallocate(aux_ovlp3ks)
      endif
      if (allocated(tmp_ovlp3ks)) then
          deallocate(tmp_ovlp3ks)
      endif
      if (allocated(aux_vect)) then
          deallocate(aux_vect)
      endif
      if (allocated(aux_vect_2)) then
          deallocate(aux_vect_2)
      endif
      if (allocated(aux_spect)) then
          deallocate(aux_spect)
      endif

      return

      end subroutine evaluate_sosex_energy


!******
!---------------------------------------------------------------------------------------------------
!****s* FHI-aims/evaluate_sosex_energy_2
!  NAME
!   evaluate_sosex_energy_2
!  SYNOPSIS

      subroutine evaluate_sosex_energy_2 &
       (n_low_state, n_KS_states, n_homo,&
        n_freq, omega_full, womega_full,&
        occ_numbers, KS_eigenvalue, ovlp_3KS, e_sosex )

!  PURPOSE
!  Calculate the screended second order exchange energy.
!
!  Note: It amounts a resummation of one type of exchange
!  type diagrams to infinite order.
!
!  This routine does the same as evaluate_polarisability_freq but it uses
!  a differently distributed ovlp_3KS and outputs a 2D distributed polar_freq


!  USES

      use constants
      use dimensions
      use prodbas
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use runtime_choices
      use constants
      use evaluate_polarisability_freq, only: evaluate_polarisability_freq_2
      use localorb_io, only: use_unit

      implicit none

      integer :: n_homo(n_spin)
      integer :: n_low_state
      integer :: n_KS_states
      integer :: n_freq
      real*8  :: omega_full(n_freq)
      real*8  :: womega_full(n_freq)
      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  :: e_sosex

! INPUTS
! o  n_low_state -- the first KS state to be included in frozen-core
!                   calculation
! o  n_homo -- integer array, the HOMO level for each spin channel
! o  n_KS_states -- integer number,
! o          the highest KS/HF eigenstate, the same as the n_high_state.
! o          In the present case, n_high_state >= n_homo should do fine.
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  ovlp_3KS -- real array
!            this is the transformed 3-cener overlap integration. Now two orbitals of
!            them are KS ones, and one is the auxiliary basis.
!            Note: Parallel distribution is among 2nd and 3rd dimension
! OUTPUT
! o  e_sosex -- real number, the calculated second-order exchange energy
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

      real*8, dimension(:,:), allocatable :: polar_freq
      real*8, dimension(:,:), allocatable :: polar_lambda
      real*8, dimension(:,:,:), allocatable :: dielec_func
      real*8, dimension(:,:), allocatable :: aux_eri
      real*8, dimension(:), allocatable :: aux_eri_real
      real*8, dimension(:,:), allocatable :: tmp_ovlp3KS_1
      real*8, dimension(:,:), allocatable :: tmp_ovlp3KS_2
      real*8, dimension(:), allocatable :: aux_vect
      real*8, dimension(:), allocatable :: aux_vect_loc1
      real*8, dimension(:), allocatable :: aux_vect_loc2
      real*8, dimension(:), allocatable :: aux_spect

      real*8 :: e_diff

      integer ::  n_lumo(n_spin)
      integer ::  n_homo_max
      integer ::  n_lumo_min
      integer ::  min_lumo_loc
      integer ::  max_homo_loc
      integer ::  n_low_state_loc
      integer ::  max_local_states1
      integer ::  max_local_states2


      integer :: n_lambda
      parameter (n_lambda=5)
      real*8   lambda_grid(n_lambda)
      real*8   lambda_weight(n_lambda)

!     Integer counters
      integer i_spin
      integer i_spin2
      integer i_state
      integer i_state_2
      integer i_state_3
      integer i_state_4
      integer i_freq
      integer i_lambda
      integer i_loc
      integer i

!     Error variable
      integer mpierr


      if (myid.eq.0) then
          write(use_unit,'(2A)') " ---------------------------------", &
                    "----------------------------------------" 
          write(use_unit,'(2A)') " | Evaluating 2nd order screened", &
                " exchange (SOSEX) energy (2D distribution) ..."
          write(use_unit,'(A)') 
      endif

      n_lumo(:) = n_homo(:) + 1
   
      n_homo_max = max(n_homo(1),n_homo(n_spin))
      n_lumo_min = min(n_lumo(1),n_lumo(n_spin))

      min_lumo_loc = loc_dim1_o3ks(n_lumo_min) ! overall minimal local value corresponding to n_lumo_min
      max_homo_loc = loc_dim2_o3ks(n_homo_max) ! overall maximal local value corresponding to n_homo_max
      n_low_state_loc = loc_dim2_o3ks(n_low_state)

      max_local_states1 = ndim1_o3KS - min_lumo_loc + 1       ! n_unocc locally
      max_local_states2 = max_homo_loc - n_low_state_loc + 1  ! n_occ locally


      allocate(dielec_func(max_row_2d, max_col_2d, n_freq))
      allocate(polar_freq(max_row_2d, max_col_2d))
      allocate(polar_lambda(max_row_2d, max_col_2d))

      e_sosex = 0.d0

! Gauss-Legendre grid for the coupling constant integral
      call gauleg(0d0,1.d0,lambda_grid,lambda_weight,n_lambda)

      dielec_func = 0.d0
      e_sosex = 0.d0
      do i_freq = 1, n_freq, 1

! compute the symmetrized polarizability (v^(1/2) * \chi_0 * v^(1/2) )
        call evaluate_polarisability_freq_2 &
           ( n_low_state, n_homo, n_lumo, &
             occ_numbers, omega_full(i_freq), &
             KS_eigenvalue, ovlp_3KS, polar_freq &
           )

        polar_freq(:,:) = polar_freq(:,:) * 2.d0/dble(n_spin)

        do i_lambda = 1, n_lambda, 1

            call PDLASET( 'Full', n_basbas, n_basbas, 0.d0, 1.d0, polar_lambda, 1, 1, aux_sc_desc_2d )

            polar_lambda(:,:) = polar_lambda(:,:) - polar_freq(:,:) * lambda_grid(i_lambda)

            ! Note: use_2d_corr (2D distribution) only works with scalapack
            call power_auxmat_scalapack_2d(polar_lambda,-1.d0," ")

            dielec_func(:,:,i_freq) = dielec_func(:,:,i_freq) + &
                 lambda_grid(i_lambda)*polar_lambda(:,:)*lambda_weight(i_lambda)
        enddo

      enddo

      dielec_func(:,:,:) =  dielec_func(:,:,:)*2.d0

      deallocate(polar_freq)
      deallocate(polar_lambda)

      allocate(aux_eri(min_lumo_loc:ndim1_o3KS,n_low_state_loc:max_homo_loc))
      allocate(aux_eri_real(min_lumo_loc:ndim1_o3KS))
      allocate(tmp_ovlp3KS_1(n_basbas,min_lumo_loc:ndim1_o3KS))
      allocate(tmp_ovlp3KS_2(n_basbas,n_low_state_loc:max_homo_loc))
      allocate(aux_vect(n_basbas))
      allocate(aux_vect_loc1(max_col_2d))
      allocate(aux_vect_loc2(max_row_2d))
      allocate(aux_spect(n_freq))


      do i_spin = 1, n_spin, 1

         do i_state_2 = n_lumo(i_spin), n_states, 1

           if(myid.eq.0 .and. mod(i_state_2,10).eq.0 ) then
               write(use_unit,'(A,2I6)') " | spin, virtual orb : ", i_spin, i_state_2 
           endif
              if(own_dim1_o3ks(i_state_2) == myp1_o3KS) &
                 tmp_ovlp3KS_2 = ovlp_3KS(:,loc_dim1_o3ks(i_state_2),n_low_state_loc:max_homo_loc,i_spin)
              call mpi_bcast(tmp_ovlp3KS_2,n_basbas*max_local_states2,MPI_REAL8,own_dim1_o3ks(i_state_2), &
                             mpi_comm_o3ks_2, mpierr)

           do i_state = n_low_state, n_homo(i_spin), 1


              if(own_dim2_o3ks(i_state) == myp2_o3KS) &
                 tmp_ovlp3KS_1(:,:) = ovlp_3KS(:,min_lumo_loc:ndim1_o3KS,loc_dim2_o3ks(i_state),i_spin)
              call mpi_bcast(tmp_ovlp3KS_1,n_basbas*max_local_states1,MPI_REAL8,own_dim2_o3ks(i_state), &
                             mpi_comm_o3ks_1, mpierr)

              call dgemm('T', 'N', max_local_states1, max_local_states2, &
                      n_basbas, 1.0d0, &
                      tmp_ovlp3KS_1, n_basbas, &
                      tmp_ovlp3KS_2, n_basbas, 0.d0, &
                      aux_eri, max_local_states1 )

              do i_freq = 1, n_freq, 1

                aux_vect(:) = 0

                do i_state_3 = n_low_state, n_homo(i_spin), 1
                  if(own_dim2_o3ks(i_state_3) /= myp2_o3KS) cycle
                  aux_eri_real(:) = 0
                  do i_state_4 = n_lumo(i_spin), n_states, 1
                    if(own_dim1_o3ks(i_state_4) /= myp1_o3KS) cycle

                    e_diff = KS_eigenvalue(i_state_3,i_spin,1) - KS_eigenvalue(i_state_4,i_spin,1)

!                    aux_eri_real(loc_dim1_o3ks(i_state_4)) = &
!                       aux_eri(loc_dim1_o3ks(i_state_4),loc_dim2_o3ks(i_state_3)-n_low_state_loc+1) * 2.d0 *   &
                    aux_eri_real(loc_dim1_o3ks(i_state_4)) = &
                       aux_eri(loc_dim1_o3ks(i_state_4),loc_dim2_o3ks(i_state_3)) * 2.d0 *   &
                       e_diff/(e_diff*e_diff + omega_full(i_freq)*omega_full(i_freq) )

                  enddo
                  call dgemv('N', n_basbas, max_local_states1, 1.d0, &
                             ovlp_3KS(1,min_lumo_loc,loc_dim2_o3ks(i_state_3),i_spin), n_basbas, &
                             aux_eri_real(min_lumo_loc), 1, &
                             1.d0, aux_vect, 1)
                enddo

                call sync_vector(aux_vect,n_basbas)

                ! dielec_func is distributed in Scalapack manner, aux_vect is full
                ! For multiplying, aux_vect must be distributed
                i_loc = 0
                do i = 1, n_basbas
                  if(MOD((i-1)/nb_aux_2d,npcol_aux_2d)==mypcol_aux_2d) then
                    i_loc = i_loc+1
                    aux_vect_loc1(i_loc) = aux_vect(i)
                  endif
                enddo

                ! Multiply distributed vector with distributed matrix
                call dgemv('N', max_row_2d, max_col_2d, 1.d0, &
                     dielec_func(1,1,i_freq), max_row_2d, aux_vect_loc1, &
                     1, 0.d0, aux_vect_loc2, 1)

                ! Assemble full vector again
                aux_vect = 0
                i_loc = 0
                do i = 1, n_basbas
                  if(MOD((i-1)/nb_aux_2d,nprow_aux_2d)==myprow_aux_2d) then
                    i_loc = i_loc+1
                    aux_vect(i) = aux_vect_loc2(i_loc)
                  endif
                enddo
                call sync_vector(aux_vect,n_basbas)

                if(own_dim1_o3ks(i_state_2) == myp1_o3KS .and. own_dim2_o3ks(i_state) == myp2_o3KS) then
                   aux_spect(i_freq) = &
                      dot_product(ovlp_3KS(:,loc_dim1_o3ks(i_state_2),loc_dim2_o3ks(i_state),i_spin),aux_vect)
                else
                   aux_spect(i_freq) = 0.
                endif

              enddo

              e_diff = KS_eigenvalue(i_state,i_spin,1) - KS_eigenvalue(i_state_2,i_spin,1)
              do i_freq = 1, n_freq, 1
                  e_sosex = e_sosex +  &
                    aux_spect(i_freq) * 2.d0*e_diff/ & 
                    (omega_full(i_freq)*omega_full(i_freq)+e_diff*e_diff) *  &
                    occ_numbers(i_state,i_spin,1)* womega_full(i_freq)
              enddo

!   End loop i_state
            enddo

!   End loop i_state_2
        enddo
!   End loop i_spin
      enddo

      call sync_real_number(e_sosex)
      e_sosex = e_sosex/pi/4.d0

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A,f19.8,A,f19.8,A)') &
         " 2nd-order screened exchange energy : ",  &
          e_sosex, '  Ha. ',  e_sosex*hartree,  ' eV '
        write(use_unit,*)
      endif

      deallocate(dielec_func)
      deallocate(aux_eri)
      deallocate(aux_eri_real)
      deallocate(tmp_ovlp3ks_1)
      deallocate(tmp_ovlp3ks_2)
      deallocate(aux_vect)
      deallocate(aux_vect_loc1)
      deallocate(aux_vect_loc2)
      deallocate(aux_spect)

      return

      end subroutine evaluate_sosex_energy_2

!******
!---------------------------------------------------------------------------------------------------
