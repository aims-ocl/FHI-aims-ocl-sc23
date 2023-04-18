!****h* FHI-aims/evaluate_polarizability
!  NAME
!   evaluate_ospt2_polarizability_freq
!  SYNOPSIS

module evaluate_ospt2_polarizability_freq
   
      use dimensions
      use prodbas
      use crpa_blacs
      use mpi_tasks
      use synchronize_mpi
      use arch_specific, only: arch_erf
      !use runtime_choices, only: renormalized_eg

      implicit none
      
      private :: renorm_eg_wrapper

      public :: evaluate_ospt2_polarizability_freq_0, evaluate_ospt2_polarizability_freq_2
contains 
! **************************************************************************************************

      subroutine evaluate_ospt2_polarizability_freq_0 &
           ( n_low_state,n_homo,n_first,n_KS_states, &
             occ_numbers, omega_n, &
             KS_eigenvalue, ovlp_3KS, polar_freq, &
             renorm_eg &
            )


!  PURPOSE
!  Subroutine evaluate_ospt2_polarizability_freq  evaluates the non-interacting 
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

      real*8  polar_freq(n_basbas, n_loc_prodbas,n_spin)
      real*8  renorm_eg(10)

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

      real*8 :: eg, seg, p1, p2, dc

!     begin work

      dc = 0.0d0
      if (omega_n.le.1.0D-6) dc=1.0D-3

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
!          do k_state = n_first(i_spin), n_states, 1
!            do i_basbas = 1, n_loc_prodbas, 1
!              i_index = map_prodbas(i_basbas, myid+1) 
!              if(i_index .eq. 21) then
!               write(use_unit,'(4I4,f20.12)')i_spin, j_state, k_state, i_index, ovlp_3KS(i_basbas,k_state,j_state,i_spin)
!              endif
!            enddo
!          enddo 
!         endif

          temp_ovlp3KS(:,:) = 0.d0
          aux_ovlp3KS(:,:) = 0.d0
          do k_state =  n_first(i_spin), n_states, 1

           ! IYZ:: for fractional spin and charge investigation
           !       we should consider the excitation between the same
           !       orbital with fractional occupation
           !       I introduce a tiny level shift to avoid the sigularity, but
           !       IN THIS CASE, VERY TIGHT FREQUENCE GRID IS NEEDED!!!
           !       As we set dc=1.0d-6, frequency_points should be more than 4000 
           !       (dc=1.0d-6; fp=4000) allows the convergence of osRPA
           !       into several meVs
           if (k_state .lt. j_state) then
               cycle
           else if (k_state.eq.j_state) then
               dc = 1.0d-6
           else
               dc = 0.0d0
           endif

           eg  = KS_eigenvalue(j_state,i_spin) - KS_eigenvalue(k_state,i_spin)
           !!p1 = (renorm_eg(4+i_spin)**(2.0d0))/&
           !!    (renorm_eg(4+i_spin)**(2.0d0)+(abs(eg))**(-abs(renorm_eg(2))))
           !p1 = arch_erf(renorm_eg(2)*log(renorm_eg(4+i_spin)))
           !p2 = renorm_eg(1)*p1*&
           !     exp(-renorm_eg(3)*abs(eg)**(abs(renorm_eg(4))))
           !!p2 = renorm_eg(1)*exp(-renorm_eg(3)*abs(1.0d0-p1)**(abs(renorm_eg(4))))
           !seg = eg*(1.0d0-p2)-dc
           call renorm_eg_wrapper(seg, eg, renorm_eg, i_spin, dc, 0)
           ! previous one by Xinguo for fractional occupation
           !zeta = 2.d0*seg/(seg**2 + omega_n*omega_n) * &
           !      (occ_numbers(j_state,i_spin)-occ_numbers(k_state,i_spin))* &
           !       dble(n_spin)/2.d0
           ! suggested by Weitao Yang for fractional occupation
           zeta = 2.d0*seg/(seg**2 + omega_n*omega_n) * &
                 ( occ_numbers(j_state,i_spin)*dble(n_spin)/2.0d0 * &
                    (1.0d0-occ_numbers(k_state,i_spin)*dble(n_spin)/2.0d0) &                                                                
                 )
           temp_ovlp3KS(:, k_state) = zeta * &
                    ovlp_3KS(:,k_state, j_state,i_spin)

          enddo ! do k_state

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
                   n_loc_prodbas, 1.d0, polar_freq(:,:,i_spin), n_basbas &
                  )
          call cpu_time(time_end)

       enddo ! do j_state
     enddo ! do  i_spin


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
!      do j_state = 1, n_loc_prodbas , 1
!         i_index = map_prodbas(j_state, myid+1) 
!         write(use_unit,'(3I4,2f18.8)') i_index, j_state, myid, polar_freq(i_index,j_state)
!      enddo


      if (allocated (aux_ovlp3KS)) then
        deallocate (aux_ovlp3KS)
      endif
      if (allocated (temp_ovlp3KS)) then
        deallocate (temp_ovlp3KS)
      endif

      return
      end subroutine evaluate_ospt2_polarizability_freq_0
!---------------------------------------------------------------------
!**********
!  NAME
!   evaluate_ospt2_polarizability_freq
!  SYNOPSIS

      subroutine evaluate_ospt2_polarizability_freq_2 &
           ( n_low_state,n_homo,n_first,n_KS_states, &
             occ_numbers, omega_n, &
             KS_eigenvalue, ovlp_3KS, polar_freq, renorm_eg,&
             do_real_freq, eta, do_imaginary &
            )


!  PURPOSE
!  Subroutine evaluate_ospt2_polarizability_freq  evaluates the non-interacting 
!  polarisability, represented within auxiliary basis.  Here imaginary 
!  frequency domain is used.
!
!  chi^_lk (w) = - \sum_m^unocc \sum_n^occ O_ln^m* O_nk^m
!                {1/(E_m-E_n - iw) + 1/(E_m - E_n +iw)}
!
!  RJ: This routine does the same as evaluate_ospt2_polarizability_freq but it uses
!      a differently distributed ovlp_3KS and outputs a 2D distributed polar_freq
!
!  ARGUMENTS

      integer :: n_low_state
      integer :: n_homo(n_spin)
      integer :: n_first(n_spin)
      integer :: n_KS_states

      real*8  occ_numbers(n_states,n_spin)
      real*8  omega_n
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  renorm_eg(10)

      real*8  polar_freq(max_row_2d, max_col_2d, n_spin)
      logical, intent(in), optional :: do_real_freq
      real*8, intent(in), optional  :: eta
      logical, intent(in), optional :: do_imaginary

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
!  o do_real_freq -- calculate polarisability for real frequencies omega_n
!  o eta -- broadening parameter for real frequency calculations
!  o do_imaginary -- logical, if imaginary or real part of the polarisability is calculated
!                    only meaningful when eta is larger than zero
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
      real*8  time_start, time_end
      real*8                              :: my_eta
      complex(kind=8)                     :: im_unit
      logical                             :: my_do_real_freq

!     counters


      integer :: j_state
      integer :: k_state
      integer :: i_spin
      integer :: n_loc_prodbas_remote
      integer :: j_state_loc
      integer :: k_state_loc
      integer :: i_c, i_r, len_c, len_r, np_c, np_r, np, i_c_loc, i_r_loc
      integer :: i_cs, i_rs, len_cs, len_rs
      integer :: i, n_tmp, n_blk, mpierr

      real*8 :: eg, seg, p1, p2, dc

      if(present(do_real_freq)) then
        my_do_real_freq = do_real_freq
        my_eta = 0.0d0
        if(present(eta)) my_eta = eta
        if(my_eta > 1.0d-10.and..not.present(do_imaginary)) then
          call aims_stop('broadening parameter is not zero; set do_imaginary')
        endif
      else
       my_do_real_freq = .FALSE.
      endif

      im_unit = (0.0d0, 1.0d0)

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

      do i_spin = 1, n_spin
        do i_r = 0, n_basbas-1, n_blk
        do i_c = i_r, n_basbas-1, n_blk

          ! size of current block
          len_c = min(n_basbas-i_c, n_blk)
          len_r = min(n_basbas-i_r, n_blk)

          polar_freq_sub(:,:) = 0.d0

          n_tmp = 0
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

              ! IYZ:: for fractional spin and charge investigation
              !       we should consider the excitation between the same
              !       orbital with fractional occupation
              !       I introduce a tiny level shift to avoid the sigularity, but
              !       IN THIS CASE, VERY TIGHT FREQUENCE GRID IS NEEDED!!!
              !       As we set dc=1.0d-6, frequency_points should be more than 4000 
              !       (dc=1.0d-6; fp=4000) allows the convergence of osRPA
              !       into several meVs
              if (k_state .lt. j_state) then
                  cycle
              else if (k_state.eq.j_state) then
                  dc = 1.0d-6
              else
                  dc = 0.0d0
              endif

              eg  = KS_eigenvalue(j_state,i_spin) - KS_eigenvalue(k_state,i_spin)
              !!p1 = (renorm_eg(4+i_spin)**(2.0d0))/&
              !!    (renorm_eg(4+i_spin)**(2.0d0)+(abs(eg))**(-abs(renorm_eg(2))))                                      
              !p1 = arch_erf(renorm_eg(2)*log(renorm_eg(4+i_spin)))
              !p2 = renorm_eg(1)*p1*&
              !     exp(-renorm_eg(3)*abs(eg)**(abs(renorm_eg(4))))
              !!p2 = renorm_eg(1)*exp(-renorm_eg(3)*abs(1.0d0-p1)**(abs(renorm_eg(4))))
              !seg = eg*(1.0d0-p2)-dc
              call renorm_eg_wrapper(seg, eg, renorm_eg, i_spin, dc, 0)
              ! previous one by Xinguo for fractional occupation
              !zeta = 2.d0*seg/(seg**2 + omega_n*omega_n) * &
              !      (occ_numbers(j_state,i_spin)-occ_numbers(k_state,i_spin))* &
              !       dble(n_spin)/2.d0
              ! suggested by Weitao Yang for fractional occupation
              zeta = 2.d0*seg/(seg**2 + omega_n*omega_n) * &
                    ( occ_numbers(j_state,i_spin)*dble(n_spin)/2.0d0 * &
                       (1.0d0-occ_numbers(k_state,i_spin)*dble(n_spin)/2.0d0) &
                    )


              temp_ovlp3KS_1(1:len_r, n_tmp) = zeta * &
                      ovlp_3KS(i_r+1:i_r+len_r, k_state_loc, j_state_loc ,i_spin)

            enddo ! k_state
          enddo ! j_state

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
               polar_freq(i_r_loc+1:i_r_loc+len_rs,i_c_loc+1:i_c_loc+len_cs,i_spin) = &
                 tmp2(1:len_rs,1:len_cs)
           enddo
           enddo

        enddo
        enddo
      enddo ! i_spin

      if (allocated (temp_ovlp3KS_1)) deallocate (temp_ovlp3KS_1)
      if (allocated (temp_ovlp3KS_2)) deallocate (temp_ovlp3KS_2)
      if (allocated (polar_freq_sub)) deallocate (polar_freq_sub)
      if (allocated (tmp1)) deallocate (tmp1)
      if (allocated (tmp2)) deallocate (tmp2)

      ! polar_freq has been set only in the upper half, complete it

      allocate(polar_freq_trans(max_row_2d, max_col_2d))
      polar_freq_trans=0.0d0

      do i_spin = 1, n_spin, 1
        polar_freq_trans = polar_freq(:,:,i_spin)
        call pdtran(n_basbas, n_basbas, 1.d0, polar_freq_trans, 1, 1, aux_sc_desc_2d, &
                    1.d0, polar_freq(:,:,i_spin), 1, 1, aux_sc_desc_2d)
        !call pdtran(n_basbas, n_basbas, 1.d0, polar_freq_trans, 1, 1, bb2desc, &
        !            1.d0, polar_freq(:,:,i_spin), 1, 1, bb2desc)
      enddo

      deallocate(polar_freq_trans)

      return
      end subroutine evaluate_ospt2_polarizability_freq_2

      subroutine renorm_eg_wrapper(seg, eg, renorm_eg, i_spin, dc, iop)
      real*8  seg 
      real*8  eg
      real*8  renorm_eg(10)
      integer i_spin
      real*8  dc
      integer iop
      ! temp
      real*8  :: p0, p1, p2, p3
      select case(iop)
      case(0)
         seg = eg-dc
      case(1)
         p1 = arch_erf(renorm_eg(2)*log(renorm_eg(4+i_spin)))
         p2 = renorm_eg(1)*p1*&
              (exp(-renorm_eg(3)*abs(eg))+exp(-renorm_eg(4)*abs(renorm_eg(4+i_spin))))
         seg = eg*(1.0d0-p2)-dc
      case(2)
         p0 = exp(-renorm_eg(1)*abs(eg)**(renorm_eg(2)))
         p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         p2 = renorm_eg(4)/(exp(eg)+exp(-renorm_eg(4+i_spin)))
         p3 = p0*arch_erf(p1+p2)
         seg = eg*(1.0d0-p3)-dc
      case(3)
         p0 = exp(-renorm_eg(1)*abs(eg)**(renorm_eg(2)))
         p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         p2 = renorm_eg(4)*(exp(renorm_eg(6+i_spin))+exp(-renorm_eg(4+i_spin)))**(-1.0d0)
         p3 = p0*arch_erf(p1+p2)
         seg = eg*(1.0d0-p3)-dc
      case(4)
         p0 = exp(-renorm_eg(1)*abs(eg)**(renorm_eg(2)))
         !p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         p2 = (renorm_eg(3)*exp(renorm_eg(6+i_spin))+renorm_eg(4)*exp(-renorm_eg(4+i_spin)))**(-1.0d0)
         p3 = p0*arch_erf(p2)
         seg = eg*(1.0d0-p3)-dc
      case(5)
         p0 = exp(-renorm_eg(1)*abs(eg)**(3.0d0))
         !p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         p2 = (renorm_eg(3)*exp(renorm_eg(6+i_spin))+renorm_eg(4)*exp(-renorm_eg(4+i_spin)))**(-1.0d0*renorm_eg(2))
         p3 = p0*arch_erf(p2)
         seg = eg*(1.0d0-p3)-dc
      case(6)
         p0 = exp(-renorm_eg(1)*abs(eg)**(2.0d0))
         !p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         p2 = 1.0d0/(renorm_eg(3)*exp(renorm_eg(6+i_spin))**(0.5d0)+renorm_eg(4)*exp(-renorm_eg(4+i_spin)))**(renorm_eg(2))
         p3 = p0*arch_erf(p2)
         seg = eg*(1.0d0-p3)-dc
      case(7)
         p0 = exp(-renorm_eg(1)*abs(eg)**(renorm_eg(2)))
         !p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         !p2 = 1.0d0/(renorm_eg(3)*exp(renorm_eg(6+i_spin))**(0.5d0)+renorm_eg(4)*exp(-renorm_eg(4+i_spin)))**(renorm_eg(2))
         p2 = exp(renorm_eg(3)*renorm_eg(6+i_spin))+exp(-renorm_eg(4)*renorm_eg(4+i_spin))
         p3 = p0*p2
         seg = eg*p3-dc
      case(8)
         p0 = exp(-renorm_eg(1)*abs(eg)**(3.0d0))
         !p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         p2 = (exp(renorm_eg(3)*renorm_eg(6+i_spin))+exp(-renorm_eg(4)*renorm_eg(4+i_spin)))**(-1.0d0*renorm_eg(2))
         p3 = p0*arch_erf(p2)
         seg = eg*(1.0d0-p3)-dc
      case(9)
         p0 = exp(-renorm_eg(1)*abs(eg)**(3.0d0))
         !p1 = renorm_eg(3)*log(renorm_eg(4+i_spin))
         !p2 = 1.0d0/(renorm_eg(3)*exp(renorm_eg(6+i_spin))**(0.5d0)+renorm_eg(4)*exp(-renorm_eg(4+i_spin)))**(renorm_eg(2))
         p2 = (exp(renorm_eg(3)*renorm_eg(6+i_spin))+exp(-renorm_eg(4)*renorm_eg(4+i_spin)))**(renorm_eg(2))
         p3 = p0*p2
         seg = eg*p3-dc
      case(10)
         p0 = exp(-renorm_eg(1)*abs(eg+renorm_eg(3))**(renorm_eg(2)))
         p1 = (1.0d0-p0*arch_erf(renorm_eg(4)*log(renorm_eg(4+i_spin)+1.0d0)))
         seg = eg*p1-dc
      case(11)
         p1 = arch_erf(renorm_eg(2)*log(renorm_eg(4+i_spin)))
         p2 = renorm_eg(1)*p1*&
              exp(-renorm_eg(3)*abs(eg)**(abs(renorm_eg(4))))
         seg = eg*(1.0d0-p2)-dc
      end select
      end subroutine renorm_eg_wrapper

end module evaluate_ospt2_polarizability_freq
