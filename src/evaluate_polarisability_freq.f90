!****h* FHI-aims/evaluate_polarisability
!  NAME
!   evaluate_polarisability_freq
!  SYNOPSIS

module evaluate_polarisability_freq
   
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use timing,                       only: get_timestamps,&
                                              get_times

      implicit none
      
      private

      public :: evaluate_polarisability_freq_0, evaluate_polarisability_freq_1, &
                evaluate_polarisability_freq_2, evaluate_polarisability_freq_3
contains 
! **************************************************************************************************

      subroutine evaluate_polarisability_freq_0 &
           ( n_low_state,n_homo,n_first,n_KS_states, &
             occ_numbers, omega_n, &
             KS_eigenvalue, ovlp_3KS, polar_freq &
            )


!  PURPOSE
!  Subroutine evaluate_polarisability_freq  evaluates the non-interacting 
!  polarisability, represented within auxiliary basis.  Here imaginary 
!  frequency domain is used.
!
!  chi^_lk (w) = - \sum_m^unocc \sum_n^occ O_ln^m* O_nk^m
!                {1/(E_m-E_n - iw) + 1/(E_m - E_n +iw)}
!
!  USES

!  ARGUMENTS

      integer :: n_low_state
      integer :: n_homo(n_spin)
      integer :: n_first(n_spin)
      integer :: n_KS_states

      real*8  occ_numbers(n_states,n_spin)
      real*8  omega_n
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  ovlp_3KS(n_loc_prodbas, n_states,n_KS_states,n_spin)

      real*8  polar_freq(n_basbas, n_loc_prodbas)

!  INPUTS
!  o n_low_state -- the lowest KS/HF eigenstate to construct the polarisability
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o n_first -- the first non-fully occupied eigenstate, differs from n_lumo
!    in case of fractional occupation
!  o n_KS_states -- the number of single-particle states to be considered in
!    GW calculation, works here provided n_KS_states >= n_homo
!  o occ_number -- occupation number for each state and each spin
!  o  KS_eigenvalue -- real array,
!          the eigenvalues of the single-particle calculation. For DFT calculation,
!          this is the KS eigenvalue, but for HF calculation, this is then the HF
!          eigenvalue
!  o  ovlp_3KS -- real array
!          this is the transformed 3-cener overlap integration. Now two orbitals of
!          them are KS ones, and one is the auxiliary basis.
!          Note: for parallel calculations, the auxiliary basis are distribuated
!          among the different processors.
!
!  OUTPUTS 
!  o polar_freq -- the calculated non-interacting polarisability in terms of
!         the auxiliary basis
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
      real*8  zeta
      real*8, dimension(:,:), allocatable :: aux_ovlp3KS
      real*8, dimension(:,:), allocatable :: temp_ovlp3KS
      real*8, dimension(:,:), allocatable :: temp
      real*8  time_start, time_end

      integer n_first_min

!     counters


      integer :: j_state
      integer :: k_state

      integer :: i_basbas
      integer :: j_basbas

      integer :: i_index
      integer :: i_spin
      integer :: i_task

      integer :: n_loc_prodbas_remote, mpierr

!     begin work

      n_first_min=n_states
      do  i_spin = 1, n_spin
       n_first_min = min(n_first_min, n_first(i_spin))
      enddo

      allocate (aux_ovlp3KS(n_basbas, n_first_min:n_states),stat=i_index)
      call check_allocation(i_index, 'aux_ovlp3KS                   ')

      allocate (temp_ovlp3KS(n_loc_prodbas, n_first_min:n_states),stat=i_index)
      call check_allocation(i_index, 'temp_ovlp3KS                  ')


      polar_freq=0.d0

      i_index =0
      do i_spin = 1, n_spin
        do j_state = n_low_state, n_homo(i_spin), 1
!         if(myid.eq.0) then
!          write(use_unit,*) i_spin, j_state, n_homo(i_spin)
          do k_state = n_first(i_spin), n_states, 1
            do i_basbas = 1, n_loc_prodbas, 1
              i_index = map_prodbas(i_basbas, myid+1) 
!              if(i_index .eq. 21) then
!               write(use_unit,'(4I4,f20.12)')i_spin, j_state, k_state, i_index, ovlp_3KS(i_basbas,k_state,j_state,i_spin)
!              endif
            enddo
          enddo 
!         endif

          temp_ovlp3KS(:,:) = 0.d0
          aux_ovlp3KS(:,:) = 0.d0
          do k_state =  n_first(i_spin), n_states, 1
           if (k_state .le. j_state) cycle

           zeta = 2.d0*(KS_eigenvalue(j_state,i_spin) - &
                     KS_eigenvalue(k_state,i_spin))/ &
                 (( KS_eigenvalue(j_state,i_spin) - &
                  KS_eigenvalue(k_state,i_spin))**2 + omega_n*omega_n) * &
                 (occ_numbers(j_state,i_spin)-occ_numbers(k_state,i_spin))* &
                  dble(n_spin)/2.d0

            temp_ovlp3KS(:, k_state) = zeta * &
                    ovlp_3KS(:, k_state, j_state,i_spin)
!           write(use_unit,*) j_state, k_state
!           write(use_unit,*) occ_numbers(j_state,i_spin), occ_numbers(k_state,i_spin)
!           write(use_unit,*) (KS_eigenvalue(j_state,i_spin) - KS_eigenvalue(k_state,i_spin))**2 + omega_n*omega_n
!            if(myid.eq.0) then
!              write(use_unit,'(2I6,2f18.8))') k_state, j_state, zeta, ovlp_3KS(1, k_state, j_state,i_spin)
!            endif

          enddo

!   upfolding

!!          do i_task = 1, n_tasks
!!             do i_basbas = 1, n_loc_prodbas
!!
!!               i_index = map_prodbas(i_basbas,i_task)
!!               if(i_index.gt.0 .and. myid.eq.i_task-1) then
!!
!!                 aux_ovlp3KS(i_index, n_first_min:n_states) = &
!!                 temp_ovlp3KS(i_basbas, n_first_min:n_states)
!!               endif
!!
!!             enddo
!!          enddo
!!
!!!   synchronizing the upfolded matrix
!!          call sync_matrix(aux_ovlp3KS(1:n_basbas,n_first_min:n_states), &
!!                       n_basbas, n_states-n_first_min+1)

          do i_task = 1, n_tasks

             n_loc_prodbas_remote = COUNT(map_prodbas(:,i_task)>0)
             allocate(temp(n_loc_prodbas_remote,n_first_min:n_states))

             if(i_task-1 == myid) temp(:,:) = temp_ovlp3KS(:,:)
              if(use_mpi) then
                call mpi_bcast(temp,n_loc_prodbas_remote*(n_states-n_first_min+1), &
                             MPI_REAL8, i_task-1, mpi_comm_global, mpierr)
              endif

             do i_basbas = 1, n_loc_prodbas_remote

               i_index = map_prodbas(i_basbas,i_task)

               aux_ovlp3KS(i_index, n_first_min:n_states) = &
                 temp(i_basbas, n_first_min:n_states)

             enddo

             deallocate(temp)

          enddo

          call cpu_time(time_start)
          call dgemm('N', 'T', n_basbas, n_loc_prodbas, &
                   n_states-n_first(i_spin) + 1, 1.0d0, &
                   aux_ovlp3KS(:,n_first(i_spin):n_states), &
                   n_basbas, &
                   ovlp_3KS(:,n_first(i_spin):n_states,j_state,i_spin), &
                   n_loc_prodbas, 1.d0, polar_freq, n_basbas &
                  )
          call cpu_time(time_end)

       enddo
      enddo


!    printing out
!      if (flag_polar) then
!        do i_state = 1, 100

!         write(filename,'(A,I0,A)')"polar/chichi.n_",i_state,".dat"
!         open(121, file=filename)

!           do i_freq = 1, n_freq

!              write (121, '(2X, 200F18.10)')
!     +              omega_grid(i_freq),
!     +             ( polar_freq(i_state,j_state,i_freq)*2.d0,
!     +               j_state = i_state, 100)
!           enddo
!         close(121)

!   end of do i_state
!        enddo
!   end of if flag_polar
!      endif
!     if(myid.eq.0) then
!      write(use_unit,*) "polar 0"
!      do j_state = 1, n_loc_prodbas , 1
!         i_index = map_prodbas(j_state, myid+1) 
!         write(use_unit,'(3I4,2f18.8)') i_index, j_state, myid, polar_freq(i_index,j_state)
!      enddo
!     endif


      if (allocated (aux_ovlp3KS)) then
        deallocate (aux_ovlp3KS)
      endif
      if (allocated (temp_ovlp3KS)) then
        deallocate (temp_ovlp3KS)
      endif

      return
      end subroutine evaluate_polarisability_freq_0
!---------------------------------------------------------------------
!  NAME
!   evaluate_polarisability_freq
!  SYNOPSIS

      subroutine evaluate_polarisability_freq_1 &
           ( n_low_state,n_homo,n_first,n_KS_states, &
             occ_numbers, omega_n, &
             KS_eigenvalue, ovlp_3KS, polar_freq &
            )


!  PURPOSE
!  Subroutine evaluate_polarisability_freq  evaluates the non-interacting 
!  polarisability, represented within auxiliary basis.  Here imaginary 
!  frequency domain is used.
!
!  chi^_lk (w) = - \sum_m^unocc \sum_n^occ O_ln^m* O_nk^m
!                {1/(E_m-E_n - iw) + 1/(E_m - E_n +iw)}
!
!  RJ: This routine does the same as evaluate_polarisability_freq but it uses
!      a differently distributed ovlp_3KS and outputs a full polar_freq
!
!  USES

!  ARGUMENTS

      integer :: n_low_state
      integer :: n_homo(n_spin)
      integer :: n_first(n_spin)
      integer :: n_KS_states

      real*8  occ_numbers(n_states,n_spin)
      real*8  omega_n
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  ovlp_3KS(n_basbas,n_states,(n_KS_states-1)/n_tasks+1,n_spin)

      real*8  polar_freq(n_basbas, n_basbas)

!  INPUTS
!  o n_low_state -- the lowest KS/HF eigenstate to construct the polarisability
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o n_first -- the first non-fully occupied eigenstate, differs from n_lumo
!    in case of fractional occupation
!  o n_KS_states -- the number of single-particle states to be considered in
!    GW calculation, works here provided n_KS_states >= n_homo
!  o occ_number -- occupation number for each state and each spin
!  o  KS_eigenvalue -- real array,
!          the eigenvalues of the single-particle calculation. For DFT calculation,
!          this is the KS eigenvalue, but for HF calculation, this is then the HF
!          eigenvalue
!  o  ovlp_3KS -- real array
!          this is the transformed 3-cener overlap integration. Now two orbitals of
!          them are KS ones, and one is the auxiliary basis.
!          Note: for parallel calculations, the auxiliary basis are distribuated
!          among the different processors.
!
!  OUTPUTS 
!  o polar_freq -- the calculated non-interacting polarisability in terms of
!         the auxiliary basis
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
      real*8  zeta
      real*8, dimension(:,:), allocatable :: temp_ovlp3KS
      real*8  time_start, time_end

      integer n_first_min

!     counters


      integer :: j_state
      integer :: k_state
      integer :: i_spin
      integer :: n_loc_prodbas_remote
      integer :: j_state_loc

!     begin work

      n_first_min=n_states
      do  i_spin = 1, n_spin
       n_first_min = min(n_first_min, n_first(i_spin))
      enddo

      allocate (temp_ovlp3KS(n_basbas, n_first_min:n_states),stat=j_state)
      call check_allocation(j_state, 'temp_ovlp3KS                  ')


      polar_freq=0.d0

      do i_spin = 1, n_spin
        do j_state = n_low_state, n_homo(i_spin), 1

          if(MOD(j_state-1,n_tasks) /= myid) cycle
          j_state_loc = (j_state-1)/n_tasks + 1

          temp_ovlp3KS(:,:) = 0.d0
          do k_state =  n_first(i_spin), n_states, 1
           if(k_state.le.j_state) cycle

           zeta = 2.d0*(KS_eigenvalue(j_state,i_spin) - &
                     KS_eigenvalue(k_state,i_spin))/ &
                 (( KS_eigenvalue(j_state,i_spin) - &
                  KS_eigenvalue(k_state,i_spin))**2 + omega_n*omega_n) * &
                 (occ_numbers(j_state,i_spin)-occ_numbers(k_state,i_spin))* &
                  dble(n_spin)/2.d0

            temp_ovlp3KS(:, k_state) = zeta * &
                    ovlp_3KS(:, k_state, j_state_loc ,i_spin)
          enddo

          call cpu_time(time_start)
          call dgemm('N', 'T', n_basbas, n_basbas, &
                   n_states-n_first(i_spin) + 1, 1.0d0, &
                   temp_ovlp3KS(:,n_first(i_spin):n_states), &
                   n_basbas, &
                   ovlp_3KS(:,n_first(i_spin):n_states,j_state_loc,i_spin), &
                   n_basbas, 1.d0, polar_freq, n_basbas &
                  )
          call cpu_time(time_end)
       enddo
      enddo

      call sync_matrix(polar_freq, n_basbas, n_basbas)

      if (allocated (temp_ovlp3KS)) then
        deallocate (temp_ovlp3KS)
      endif

!     if(myid.eq.0) then
!      write(use_unit,*) "polar 1"
!      do j_state = 1, 1 , 1
!         write(use_unit,'(3I4,2f18.8)') j_state, j_state, myid, polar_freq(j_state,j_state)
!      enddo
!     endif

      return
      end subroutine evaluate_polarisability_freq_1
!---------------------------------------------------------------------
!**********
!  NAME
!   evaluate_polarisability_freq
!  SYNOPSIS

      subroutine evaluate_polarisability_freq_2 &
           ( n_low_state,n_homo,n_first, &
             occ_numbers, omega_n, &
             KS_eigenvalue, ovlp_3KS, polar_freq, do_real_freq, eta, &
             do_imaginary, do_scalapack, time_polar) 


!  PURPOSE
!  Subroutine evaluate_polarisability_freq  evaluates the non-interacting 
!  polarisability, represented within auxiliary basis.  Here imaginary 
!  frequency domain is used.
!
!  chi^_lk (w) = - \sum_m^unocc \sum_n^occ O_ln^m* O_nk^m
!                {1/(E_m-E_n - iw) + 1/(E_m - E_n +iw)}
!
!  RJ: This routine does the same as evaluate_polarisability_freq but it uses
!      a differently distributed ovlp_3KS and outputs a 2D distributed polar_freq
!
!  ARGUMENTS

      integer :: n_low_state
      integer :: n_homo(n_spin)
      integer :: n_first(n_spin)

      real*8  occ_numbers(n_states,n_spin)
      real*8  omega_n
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)

      real*8  polar_freq(max_row_2d, max_col_2d)
      logical, intent(in), optional   :: do_real_freq
      real*8, intent(in), optional    :: eta
      logical, intent(in), optional   :: do_imaginary
      logical, intent(in), optional   :: do_scalapack
      real*8, intent(inout), optional :: time_polar(2)

!  INPUTS
!  o n_low_state -- the lowest KS/HF eigenstate to construct the polarisability
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o n_first -- the first non-fully occupied eigenstate, differs from n_lumo
!    in case of fractional occupation
!  o occ_number -- occupation number for each state and each spin
!  o  KS_eigenvalue -- real array,
!          the eigenvalues of the single-particle calculation. For DFT calculation,
!          this is the KS eigenvalue, but for HF calculation, this is then the HF
!          eigenvalue
!  o  ovlp_3KS -- real array
!          this is the transformed 3-cener overlap integration. Now two orbitals of
!          them are KS ones, and one is the auxiliary basis.
!          Note: for parallel calculations, the auxiliary basis are distribuated
!          among the different processors.
!  o do_real_freq -- calculate polarisability for real frequencies omega_n
!  o eta -- broadening parameter for real frequency calculations
!  o do_imaginary -- logical, if imaginary or real part of the polarisability is calculated
!                    only meaningful when eta is larger than zero
!  o do_scalapack -- logical, if we run this routine with scalapack or lapack
!  o time_polar -- 2D real array containing cpu and clock time for evaluating 
!                  polarisibility; remember to sync cpu time in the end
!
!  OUTPUTS 
!  o polar_freq -- the calculated non-interacting polarisability in terms of
!         the auxiliary basis
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
      real*8  zeta
      real*8                              :: epsilon_i, epsilon_a
      real*8, dimension(:,:), allocatable :: temp_ovlp3KS_1, temp_ovlp3KS_2
      real*8, dimension(:,:), allocatable :: polar_freq_sub, polar_freq_trans
      real*8, dimension(:,:), allocatable :: tmp1, tmp2
      real*8, dimension(2)                :: tmp_time
      real*8                              :: my_eta
      real*8                              :: occ_corr
      logical                             :: my_do_real_freq
      logical                             :: my_do_scalapack

!     counters

      integer :: ii,jj
      real(kind=8) :: diff
      integer :: j_state
      integer :: k_state
      integer :: i_spin
      integer :: n_loc_prodbas_remote
      integer :: j_state_loc
      integer :: k_state_loc
      integer :: i_c, i_r, len_c, len_r, np_c, np_r, np, i_c_loc, i_r_loc
      integer :: i_cs, i_rs, len_cs, len_rs
      integer :: i, n_tmp, n_blk, mpierr

      if(present(time_polar)) then
        tmp_time(:) = 0.d0
        call get_timestamps(tmp_time(1), tmp_time(2))
      endif

      my_do_real_freq = .FALSE.
      if(present(do_real_freq)) then
        my_do_real_freq = do_real_freq
      endif

      my_eta = 0.0d0
      if(present(eta)) my_eta = eta
      if(my_eta > 1.0d-10.and..not.present(do_imaginary)) then
        call aims_stop('broadening parameter is not zero; set do_imaginary')
      endif
      
      
      my_do_scalapack = .TRUE. 
      if(present(do_scalapack)) then
         my_do_scalapack = do_scalapack
      endif
  
!     begin work

      ! Block size for matrix multiplies, must be a multiple of nb_aux_2d

      n_blk = (255/nb_aux_2d+1)*nb_aux_2d

      allocate (temp_ovlp3KS_1(n_blk,ndim1_o3KS*ndim2_o3KS*n_spin))
      allocate (temp_ovlp3KS_2(n_blk,ndim1_o3KS*ndim2_o3KS*n_spin))
      allocate (polar_freq_sub(n_blk,n_blk))
      allocate (tmp1(nb_aux_2d,nb_aux_2d))
      allocate (tmp2(nb_aux_2d,nb_aux_2d))


      polar_freq=0.d0

      tmp1 = 0
      tmp2 = 0

      do i_r = 0, n_basbas-1, n_blk
      do i_c = i_r, n_basbas-1, n_blk

        ! size of current block
        len_c = min(n_basbas-i_c, n_blk)
        len_r = min(n_basbas-i_r, n_blk)

        polar_freq_sub(:,:) = 0.d0

        n_tmp = 0
        do i_spin = 1, n_spin
          do j_state = n_low_state, n_homo(i_spin), 1 

            if(own_dim2_o3KS(j_state) /= myp2_o3KS) cycle
            j_state_loc = loc_dim2_o3KS(j_state)

            do k_state =  n_first(i_spin), n_states, 1

              if(own_dim1_o3KS(k_state) /= myp1_o3KS .or. k_state .le. j_state) cycle
              k_state_loc = loc_dim1_o3KS(k_state)

              n_tmp = n_tmp+1

              temp_ovlp3KS_2(1:len_c, n_tmp) = &
                      ovlp_3KS(i_c+1:i_c+len_c, k_state_loc, j_state_loc ,i_spin)

              if(i_c>i_r) cycle ! temp_ovlp3KS_1 needs to be set only once per i_r loop

              epsilon_i = KS_eigenvalue(j_state,i_spin)
              epsilon_a = KS_eigenvalue(k_state,i_spin)
              zeta = get_zeta_for_polarisability(omega_n, epsilon_i, epsilon_a,&
                                                 my_do_real_freq, my_eta, do_imaginary) 
              occ_corr = (occ_numbers(j_state,i_spin)-occ_numbers(k_state,i_spin))* &
                         dble(n_spin)/2.d0
              zeta = zeta*occ_corr 

              temp_ovlp3KS_1(1:len_r, n_tmp) = zeta * &
                      ovlp_3KS(i_r+1:i_r+len_r, k_state_loc, j_state_loc ,i_spin)

!              if(myid.eq.0) then
!                write(use_unit,'(2I6,2f16.8)') k_state, j_state, zeta, temp_ovlp3KS_2(1,n_tmp)
!              endif
            enddo

          enddo
        enddo

        call dgemm('N', 'T', len_r, len_c, &
                   n_tmp, 1.0d0, &
                   temp_ovlp3KS_1, &
                   ubound(temp_ovlp3KS_1,1), &
                   temp_ovlp3KS_2, &
                   ubound(temp_ovlp3KS_2,1), 1.d0, &
                   polar_freq_sub, ubound(polar_freq_sub,1) &
                  )

        ! Only the upper half of polar_freq must be set,
        ! the diagonal must be multiplied by 0.5 (for PDTRAN below)
        if(i_r == i_c) then
           do i=1,len_c
              polar_freq_sub(i,i) = 0.5d0*polar_freq_sub(i,i)
              polar_freq_sub(i+1:,i) = 0.d0
           enddo
        endif

        ! Sum up and distribute polar_freq_sub to the owners

        do i_rs = i_r, min(i_r+n_blk-1, n_basbas-1), nb_aux_2d
        do i_cs = i_c, min(i_c+n_blk-1, n_basbas-1), nb_aux_2d

           if(i_cs < i_rs) cycle ! no need to deal with null blocks

           ! size of current block
           len_cs = min(n_basbas-i_cs, nb_aux_2d)
           len_rs = min(n_basbas-i_rs, nb_aux_2d)

           if(my_do_scalapack) then 
             ! processor col/row of current block
             np_c = MOD(i_cs/nb_aux_2d,npcol_aux_2d)
             np_r = MOD(i_rs/nb_aux_2d,nprow_aux_2d)
             np = global_id(np_r,np_c) ! global owner of current block

             ! local offset of current block
             i_c_loc = (i_cs/(nb_aux_2d*npcol_aux_2d))*nb_aux_2d
             i_r_loc = (i_rs/(nb_aux_2d*nprow_aux_2d))*nb_aux_2d

             tmp1(1:len_rs,1:len_cs) = polar_freq_sub(i_rs-i_r+1:i_rs-i_r+len_rs,i_cs-i_c+1:i_cs-i_c+len_cs)

             call mpi_reduce(tmp1,tmp2,nb_aux_2d*len_cs,MPI_REAL8,MPI_SUM,np,mpi_comm_global,mpierr)

             if(np==myid) &
               polar_freq(i_r_loc+1:i_r_loc+len_rs,i_c_loc+1:i_c_loc+len_cs) = &
                 tmp2(1:len_rs,1:len_cs)
           else
             polar_freq(i_rs+1:i_rs+len_rs,i_cs+1:i_cs+len_cs) &
              = polar_freq_sub(i_rs-i_r+1:i_rs-i_r+len_rs,i_cs-i_c+1:i_cs-i_c+len_cs)
           endif
         enddo
         enddo

      enddo
      enddo

      if (allocated (temp_ovlp3KS_1)) deallocate (temp_ovlp3KS_1)
      if (allocated (temp_ovlp3KS_2)) deallocate (temp_ovlp3KS_2)
      if (allocated (polar_freq_sub)) deallocate (polar_freq_sub)
      if (allocated (tmp1)) deallocate (tmp1)
      if (allocated (tmp2)) deallocate (tmp2)

      ! polar_freq has been set only in the upper half, complete it

      allocate(polar_freq_trans(max_row_2d, max_col_2d))
      polar_freq_trans = polar_freq

      if(my_do_scalapack) then
        call pdtran(n_basbas, n_basbas, 1.d0, polar_freq_trans, 1, 1, aux_sc_desc_2d, &
                    1.d0, polar_freq, 1, 1, aux_sc_desc_2d)
      else
       polar_freq = polar_freq + TRANSPOSE(polar_freq_trans)
      endif

      deallocate(polar_freq_trans)
   
      if(present(time_polar)) call get_times(tmp_time(1), tmp_time(2), time_polar(1), &
                                             time_polar(2), unsynced=.TRUE.) 

!      write(use_unit,*) "polar_2"
!      write(use_unit,*) polar_freq(1,1)
      return

      end subroutine evaluate_polarisability_freq_2
!---------------------------------------------------------------------
!**********
!  NAME
!   evaluate_polarisability_freq_3
!  SYNOPSIS

      subroutine evaluate_polarisability_freq_3 &
           ( n_low_state,n_homo,n_first, n_loc_states, &
             occ_numbers, omega_n, &
             KS_eigenvalue, KS_eigenvector, transformed_ovlp3fn, polar_freq &
            )


!  PURPOSE
!  Subroutine evaluate_polarisability_freq  evaluates the non-interacting 
!  polarisability, represented within auxiliary basis.  Here imaginary 
!  frequency domain is used.
!
!  chi^_lk (w) = - \sum_m^unocc \sum_n^occ O_ln^m* O_nk^m
!                {1/(E_m-E_n - iw) + 1/(E_m - E_n +iw)}
!
!  RJ: This routine does the same as evaluate_polarisability_freq but it uses
!      a differently distributed ovlp_3KS and outputs a full polar_freq
!
!  XR: The different versions differ in its data structure. In this particular
!      version, the input 3-index overlap matrix are between the auxiliary basis,
!      atomic basis, and single-particle orbtials. The output "polar_freq" has
!      dimension (n_basbas, n_basbas)
!
!  ARGUMENTS

      integer :: n_low_state
      integer :: n_homo(n_spin)
      integer :: n_first(n_spin)
      integer :: n_loc_states

      real*8  occ_numbers(n_states,n_spin)
      real*8  omega_n
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  KS_eigenvector(n_basis,n_states,n_spin)
      real*8  transformed_ovlp3fn(n_basbas,n_basis,n_loc_states,n_spin)

      real*8  polar_freq(n_basbas, n_basbas)

!  INPUTS
!  o n_low_state -- the lowest KS/HF eigenstate to construct the polarisability
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o n_first -- the first non-fully occupied eigenstate, differs from n_lumo
!    in case of fractional occupation
!  o n_loc_states -- (n_states-1)/n_tasks + 1, number of single-particle states locally
!  o occ_number -- occupation number for each state and each spin
!  o  KS_eigenvalue -- real array,
!          the eigenvalues of the single-particle calculation. For DFT calculation,
!          this is the KS eigenvalue, but for HF calculation, this is then the HF
!          eigenvalue
!  o  ovlp_3KS -- real array
!          this is the transformed 3-cener overlap integration. Now two orbitals of
!          them are KS ones, and one is the auxiliary basis.
!          Note: for parallel calculations, the auxiliary basis are distribuated
!          among the different processors.
!
!  OUTPUTS 
!  o polar_freq -- the calculated non-interacting polarisability in terms of
!         the auxiliary basis
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
      real*8  zeta
      real*8, dimension(:,:), allocatable :: temp_ovlp3KS
      real*8, dimension(:,:), allocatable :: aux_ovlp3KS
      real*8  time_start, time_end

      integer n_first_min

!     counters


      integer :: j_state
      integer :: k_state
      integer :: i_spin
      integer :: n_loc_prodbas_remote
      integer :: j_state_loc

!     begin work

      n_first_min=n_states
      do  i_spin = 1, n_spin
       n_first_min = min(n_first_min, n_first(i_spin))
      enddo

      allocate (temp_ovlp3KS(n_basbas, n_first_min:n_states),stat=j_state)
      call check_allocation(j_state, 'temp_ovlp3KS                  ')
      allocate (aux_ovlp3KS(n_basbas, n_first_min:n_states),stat=j_state)
      call check_allocation(j_state, 'temp_ovlp3KS                  ')


      polar_freq=0.d0

      do i_spin = 1, n_spin
        do j_state = n_low_state, n_homo(i_spin), 1

          if(MOD(j_state-1,n_tasks) /= myid) cycle
          j_state_loc = (j_state-1)/n_tasks + 1

          call dgemm('N', 'N', n_basbas, n_states-n_first_min+1, &
                 n_basis, 1.0d0, transformed_ovlp3fn(1,1,j_state_loc,i_spin), n_basbas, &
                 KS_eigenvector(1,n_first_min,i_spin), n_basis, &
                 0.d0, temp_ovlp3KS(1,n_first_min), &
                 n_basbas)

          do k_state =  n_first(i_spin), n_states, 1
           if(k_state.le.j_state) cycle

           zeta = 2.d0*(KS_eigenvalue(j_state,i_spin) - &
                     KS_eigenvalue(k_state,i_spin))/ &
                 (( KS_eigenvalue(j_state,i_spin) - &
                  KS_eigenvalue(k_state,i_spin))**2 + omega_n*omega_n) * &
                 (occ_numbers(j_state,i_spin)-occ_numbers(k_state,i_spin))* &
                  dble(n_spin)/2.d0

            aux_ovlp3KS(:, k_state) = zeta * &
                    temp_ovlp3KS(:, k_state)
          enddo

          call cpu_time(time_start)
          call dgemm('N', 'T', n_basbas, n_basbas, &
                   n_states-n_first(i_spin) + 1, 1.0d0, &
                   temp_ovlp3KS(:,n_first(i_spin):n_states), &
                   n_basbas, &
                   aux_ovlp3KS(:,n_first(i_spin):n_states), &
                   n_basbas, 1.d0, polar_freq, n_basbas &
                  )
          call cpu_time(time_end)
       enddo
      enddo

      call sync_matrix(polar_freq, n_basbas, n_basbas)

      if (allocated (temp_ovlp3KS)) then
        deallocate (temp_ovlp3KS)
      endif
      if (allocated (aux_ovlp3KS)) then
        deallocate (aux_ovlp3KS)
      endif

      return
      end subroutine evaluate_polarisability_freq_3

! **************************************************************************************************
!> brief calculation of zeta parameter for evaluate_polarisability_freq_2 routine
!        dependence whether we have a real or imaginary frequency or if broadening is used
!  o omega_n -- frequency
!  o epsilon_i --  KS/HF eigenvalue for occupied MO  
!  o epsilon_a -- KS/HF eigenvalue for virtual MO
!  o my_do_real_freq -- whether we have a real or imaginary frequency
!  o my_eta -- broadening parameter
!  o do_imaginary -- whether we take the real or imaginary part of zeta; output is never the full
!                    complex zeta, but always it's real or imaginary part 
! **************************************************************************************************
  function get_zeta_for_polarisability(omega_n, epsilon_i, epsilon_a, my_do_real_freq,&
                                       my_eta, do_imaginary ) result(zeta)

     real(kind=8), intent(in)                                 :: omega_n
     real(kind=8), intent(in)                                 :: epsilon_i, epsilon_a
     real(kind=8), intent(in)                                 :: my_eta
     logical, intent(in)                                      :: my_do_real_freq,&
                                                                 do_imaginary
     real(kind=8)                                             :: zeta

     complex(kind=8)                                          :: im_unit

     im_unit = (0.0d0, 1.0d0)

     if(my_do_real_freq) then
       if(ABS(my_eta) < 1.0d-10) then 
        !***minus sign: we have a real frequency
         zeta = 2.d0*(epsilon_i - epsilon_a)/ &
                (( epsilon_i - epsilon_a)**2 - omega_n*omega_n ) !* occ_corr
       else
         if(do_imaginary) then
           zeta = aimag(1.0d0/(abs(omega_n) + im_unit*2.0d0*my_eta+(epsilon_i-epsilon_a))&
                       +1.0d0/(-abs(omega_n) +(epsilon_i-epsilon_a)))
         else
           zeta = real(1.0d0/(abs(omega_n) + im_unit*2.0d0*my_eta+(epsilon_i-epsilon_a))&
                       +1.0d0/(-abs(omega_n) +(epsilon_i-epsilon_a)),kind=8)
   
         endif
       endif
         
     else
       if(ABS(my_eta) < 1.0d-10) then 
         zeta = 2.d0*(epsilon_i - epsilon_a)/ &
              (( epsilon_i - epsilon_a)**2 + omega_n*omega_n) !* occ_corr
       else
         if(do_imaginary) then
           zeta = aimag(1.0d0/(im_unit*(omega_n + my_eta) + epsilon_i - epsilon_a) &
                  +1.0d0/(im_unit*(-omega_n + my_eta) + epsilon_i - epsilon_a))
         else
           zeta = real(1.0d0/(im_unit*(omega_n + my_eta) + epsilon_i - epsilon_a) &
                  +1.0d0/(im_unit*(-omega_n + my_eta) + epsilon_i - epsilon_a),kind=8)
         endif  
       endif
     endif
  end function get_zeta_for_polarisability

end module evaluate_polarisability_freq
