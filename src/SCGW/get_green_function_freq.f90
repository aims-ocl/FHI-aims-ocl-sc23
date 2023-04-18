!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  get_green_function_freq (&
         green_fn_freq, chemical_potential,&
         occ_numbers, hamiltonian, overlap_matrix )

! PURPOSE
! Evaluate the Green's function in NAO basis, in time domain:
! 
!  G_ll'(t) = Sum_n^unocc  O_ln O_nl' Theta(t) e^-(e_n - E_f)t +
!             Sum_n^occ    O_ln O_nl' Theta(-t)e^-(e_n - E_f)t 
!
! and in frequency domain, in the NAO basis, 
!
!  G_ll'(iw) = Sum_n  O_ln O_nl'  (e_n - E_f +i w )/(w^2 + ( e_n - E_f)^2 )
!
! starting from KS eigenvalue and 

      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use scgw_grid
      use poles_fit
      use localorb_io, only: use_unit

! INPUT
! - ntau               number of points of the time grid
! - nomega             number of points of the frequency grids
! - omega          the grid points of the frequency
! - tau                time points
! - n_homo             the homo level
! - KS_eigenvalue    
! - ovelap_matrix       is the ovlp matrix between kohn_Sham orbitals and NAO numerical orbitals

!OUTPUT
! - green_fn_time      is green function in the time domain and NAO basis
! - green_fn_freq      is green function in the frequency domain and NAO basis 
                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPORTANT: only the Green's function in time is multiplied on both side
! with the inverse overlap matrix, the Green's function in frequency is kept
! "naked", since it's used in the Dyson equation as it is! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     
!ARGUMENTS
      implicit none
      complex*16 green_fn_freq(n_basis,n_basis, nomega,n_spin)
      real*8 chemical_potential
      real*8 occ_numbers (n_states, n_spin)       
!EXTRA
      complex*16, dimension(:,:), allocatable ::  green_tmp 
      character*2 iter
      character*17 filename
      logical output

      real*8 hamiltonian (n_basis*(n_basis+1)/2,n_spin)
      real*8 overlap_matrix (n_basis*(n_basis+1)/2)
      real*8 full_hamiltonian (n_basis,n_basis,n_spin)
      real*8 full_overlap_matrix (n_basis,n_basis)

!parallelization on the basis
      integer , dimension (:,:), allocatable :: map_index_basis
      integer n_loc_basis
      integer n_remain_basis

!COUNTERS 
      integer i_tau
      integer i_task
      integer i_spin
      integer i_freq
      integer i_state
      integer j_state
      integer j_basis, i_basis
      integer k_basis
      integer i_index
      integer counter
      integer start
      integer k

      if(myid.eq.0)then
        write(use_unit,*)"  --- Evaluating the  Green's function ---"
      endif  

      green_fn_freq (:,:,:,:) = 0.d0

!get full Hamiltonian and full overlap matrix 
      full_hamiltonian (:,:,:) = 0.d0
      full_overlap_matrix (:,:) = 0.d0
      do i_spin = 1, n_spin
        i_index = 0
        do j_basis = 1, n_basis, 1
           do i_basis = 1, j_basis, 1
              i_index = i_index+1
              full_hamiltonian (i_basis,j_basis,i_spin) = hamiltonian (i_index,i_spin)
              full_hamiltonian (j_basis,i_basis,i_spin) = hamiltonian (i_index,i_spin)
              if(i_spin .eq. 1)then
                full_overlap_matrix (j_basis,i_basis) = overlap_matrix (i_index)
                full_overlap_matrix (i_basis,j_basis) = overlap_matrix (i_index)
              endif
           enddo
        enddo
      enddo

      if(.not. allocated(green_tmp))then
        allocate(green_tmp (n_basis,n_basis))
      endif

      ! CALCULATE THE GREEN FUNCTION IN FREQUENCY AND NAO BASIS
      do i_spin = 1, n_spin
        do i_freq = 1, nomega , 1 
          green_tmp (:,:) = (((0.d0,1.d0)*omega(i_freq) &
               + chemical_potential) * &
               full_overlap_matrix (:,:) -&
               full_hamiltonian (:,:,i_spin) )  
          
          call invert_simple_complex (green_tmp,&
            green_fn_freq(:,:,i_freq, i_spin),n_basis)
        enddo
      enddo

      if(allocated(green_tmp))then
         deallocate(green_tmp)
      endif

!      output = .false.
!      if(myid.eq.0)then
!        if(output) then
!          do i_basis = 1, n_basis, 1
!            if( i_basis.lt.10 ) then
!               write(iter,'(A,I1)') "0",i_basis
!            else
!               write(iter,'(I2)') i_basis
!            endif
!            filename = "gree_fr_VO_"//iter//".dat"
!            open(77, file=filename)
!              do i_freq = 1, nomega, 1
!                write(77,*) omega(i_freq), &
!                         real(green_fn_freq(i_basis,i_basis,i_freq,1)),&
!                        aimag(green_fn_freq(i_basis,i_basis,i_freq,1))
!              enddo
!            close(77)
!          enddo
!        endif
!      endif

      return

      end subroutine get_green_function_freq
!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  get_green_function_freq_v0 (&!nomega,
         n_homo,&
         KS_eigenvalue, green_fn_freq,&! omega, &
         ovlp_NAO_KS, chemical_potential,&! womega, &
         inv_overlap_matrix, occ_numbers)

! PURPOSE
! Evaluate the Green's function in NAO basis, in time domain:
! 
!  G_ll'(t) = Sum_n^unocc  O_ln O_nl' Theta(t) e^-(e_n - E_f)t +
!             Sum_n^occ    O_ln O_nl' Theta(-t)e^-(e_n - E_f)t 
!
! and in frequency domain, in the NAO basis, 
!
!  G_ll'(iw) = Sum_n  O_ln O_nl'  (e_n - E_f +i w )/(w^2 + ( e_n - E_f)^2 )
!
! starting from KS eigenvalue and 

      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use scgw_grid
      use poles_fit
      use localorb_io, only: use_unit

! INPUT
! - ntau               number of points of the time grid
! - nomega             number of points of the frequency grids
! - omega          the grid points of the frequency
! - tau                time points
! - n_homo             the homo level
! - KS_eigenvalue    
! - ovelap_matrix       is the ovlp matrix between kohn_Sham orbitals and NAO numerical orbitals

!OUTPUT
! - green_fn_time      is green function in the time domain and NAO basis
! - green_fn_freq      is green function in the frequency domain and NAO basis 
                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPORTANT: only the Green's function in time is multiplied on both side
! with the inverse overlap matrix, the Green's function in frequency is kept
! "naked", since it's used in the Dyson equation as it is! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     
!ARGUMENTS
      implicit none
      integer  n_homo(n_spin)
      complex*16 green_fn_freq(n_basis,n_basis, nomega,n_spin)
!n      complex*16 g_full(n_basis,n_basis, nomega,n_spin)
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8 ovlp_NAO_KS (n_states, n_basis,n_spin)
      real*8 chemical_potential
      real*8 inv_overlap_matrix (n_basis, n_basis)
      real*8 occ_numbers (n_states, n_spin)       
!EXTRA
      real*8 e_fermi
      real*8 n_part 
      complex*16, dimension(:,:), allocatable ::  green_tmp 
      character*2 iter
      character*17 filename
      logical output
      logical test
      real*8 orthonormality (n_states, n_states)
      complex*16 aux_inv_overlap_matrix (n_basis, n_basis)
      real*8 aux_overlap_matrix (n_basis, n_basis)
      real*8 overlap_matrix (n_basis*(n_basis+1)/2)
!      real*8 sum_ovlp
!      real*8 full_ovlp(n_basis, n_basis)
!      real*8 aux_g (nomega)
!      real*8 aux_g_im (nomega)
!      complex*16 aux_g_tau (ntau) 
!      complex*16 tail_g_tau (ntau) 

!parallelization on the basis
      integer , dimension (:,:), allocatable :: map_index_basis
      integer n_loc_basis
      integer n_remain_basis

      complex*16  ovlp_NAO_KS_tmp (  n_states ,n_basis)
      complex*16  ovlp_NAO_KS_tmp1 ( n_states ,n_basis)
 
!COUNTERS 
      integer i_tau
      integer i_spin
      integer i_freq
      integer i_state
      integer j_state
      integer j_basis, i_basis
      integer k_basis
      integer i_index
      integer counter
      integer start
      integer k
!spline 
!      complex*16 fit
!      real*8 spl_param (4,nomega)
!      real*8 spl_param_im (4,nomega)
!      complex*16 a,b,c,d,a1,b1,c1,d1
!      real*8 h,t
!      real*8 w1
!      complex*16 I
!      integer max_point
    
!      real*8 interval
!      integer n_points
!      real*8 delta_w 

      if(myid.eq.0)then
        write(use_unit,*)"  --- Evaluating the  Green's function ---"
      endif  

      green_fn_freq (:,:,:,:) = 0.d0
      e_fermi = chemical_potential

      !init parallelization------------------------------------ 
      n_remain_basis = MOD(n_basis, n_tasks)
      if (n_remain_basis.eq.0) then
        n_loc_basis = n_basis / n_tasks
      else
        n_loc_basis = n_basis / n_tasks + 1
      endif

      if (.not. allocated(map_index_basis))then 
        allocate(map_index_basis(n_tasks,n_loc_basis ))
      endif

      call distribute_grid (n_basis, map_index_basis,n_loc_basis)
      
      !--------------------------------------------------------

      ! CALCULATE THE GREEN FUNCTION IN FREQUENCY AND NAO BASIS
      do i_spin = 1, n_spin

      !this version is parallelized on the frequency grid_points
!       do k = 1, n_loc_grid, 1
!         i_freq = map_index (myid+1, k)
!         if(i_freq.gt.0)then
!            do i_state = 1, n_states,1
!              do j_basis = 1, n_basis, 1
!                 ovlp_NAO_KS_tmp (i_state, j_basis) = &
!                 ovlp_NAO_KS    (i_state,j_basis, i_spin) * &
!                 1.d0/((0.d0,1.d0)*(omega(i_freq))&
!                 -(KS_eigenvalue(i_state,i_spin)-e_fermi))
!
!                 ovlp_NAO_KS_tmp1 (i_state, j_basis) = &
!                 ovlp_NAO_KS    (i_state,j_basis, i_spin)
!              enddo
!            enddo
!            call zgemm ('T', 'N', n_basis, n_basis, n_states, &
!             1.d0, ovlp_NAO_KS_tmp, n_basis,&
!             ovlp_NAO_KS_tmp1, n_states, 0.d0, &
!             green_fn_freq(:,:,i_freq,i_spin), n_basis)
!         endif
!       enddo
!       print *, first_state_included

       do i_freq = 1, nomega , 1 
!            do j_basis = 1, n_basis, 1
            do k = 1, n_loc_basis, 1
             j_basis = map_index_basis (myid+1, k)

             if(j_basis .gt.0)then
             do k_basis = 1, n_basis, 1
 
              do i_state = 1 , n_states,1
    
                if(i_state .le. n_homo(i_spin))then 
                   green_fn_freq (j_basis,k_basis,i_freq, i_spin) = &
                    green_fn_freq (j_basis,k_basis,i_freq, i_spin) + &
                    occ_numbers(i_state, i_spin)*real(n_spin)/2.d0*&
                    ovlp_NAO_KS (i_state,j_basis, i_spin) *&
                    ovlp_NAO_KS (i_state,k_basis, i_spin) *&
                    1.d0/((0.d0,1.d0)*(omega(i_freq))&
                    -(KS_eigenvalue(i_state,i_spin)-e_fermi))

!                  g_full (j_basis,k_basis,i_freq, i_spin) = &
!                  green_fn_freq (j_basis,k_basis,i_freq, i_spin)
 
                elseif (i_state .gt. n_homo(i_spin))then
                   green_fn_freq (j_basis,k_basis,i_freq, i_spin) = &
                    green_fn_freq (j_basis,k_basis,i_freq, i_spin) + &
                    (1.d0 - occ_numbers(i_state, i_spin)*real(n_spin)/2.d0)*&
                    ovlp_NAO_KS (i_state,j_basis, i_spin) *&
                    ovlp_NAO_KS (i_state,k_basis, i_spin) *&
                    1.d0/((0.d0,1.d0)*(omega(i_freq))&
                    -(KS_eigenvalue(i_state,i_spin)-e_fermi))
                endif  
 
              enddo
              

!              if(first_state_included .gt.1 )then
!                do i_state = 1, first_state_included -1 , 1
!                   g_full (j_basis,k_basis,i_freq, i_spin) = & 
!                   g_full (j_basis,k_basis,i_freq, i_spin)+ &
!                   occ_numbers(i_state, i_spin)*real(n_spin)/2.d0*&
!                   ovlp_NAO_KS (i_state,j_basis, i_spin) *&
!                   ovlp_NAO_KS (i_state,k_basis, i_spin) *&
!                   1.d0/((0.d0,1.d0)*(omega(i_freq))&
!                   -(KS_eigenvalue(i_state,i_spin)-e_fermi))
!
!                enddo
!              endif

             enddo
             endif !j_basis .gt.0
            enddo
         enddo

       do i_freq = 1, nomega, 1
        call sync_matrix_complex (&
         green_fn_freq  (1:n_basis, 1:n_basis, i_freq,i_spin ), &
         n_basis, n_basis)
       enddo

      enddo! i+spin


      aux_inv_overlap_matrix (:,:) = inv_overlap_matrix (:,:)
      if(.not.allocated(green_tmp))then
        allocate(green_tmp (n_basis,n_basis))
      endif

!      if (.true.)then
       do i_spin = 1, n_spin
        do i_freq = 1, nomega, 1
!       do k = 1, n_loc_grid, 1
!         i_freq = map_index (myid+1, k)
!         if(i_freq.gt.0)then
          green_tmp (:,:)= 0.d0
          call zgemm ('N','N',n_basis,n_basis, n_basis,&
            1.d0, green_fn_freq (1,1,i_freq,i_spin), n_basis, &
            aux_inv_overlap_matrix, n_basis, 0.d0,        &
            green_tmp, n_basis)
  
         green_fn_freq (:,:,i_freq, i_spin) = 0.d0
         call zgemm ('N','N',n_basis,n_basis, n_basis,&
            1.d0, aux_inv_overlap_matrix, n_basis, &
            green_tmp, n_basis, 0.d0, &
            green_fn_freq (1,1,i_freq,i_spin), n_basis)
!        endif
  
        enddo
      enddo
!`      endif
 
      if(allocated(green_tmp))then
         deallocate(green_tmp)
      endif

      output = .false.
      if(myid.eq.0)then
      if(output) then
        do i_basis = 1, n_basis, 1
          if( i_basis.lt.10 ) then
             write(iter,'(A,I1)') "0",i_basis
          else
             write(iter,'(I2)') i_basis
          endif
          filename = "gree_fr_VO_"//iter//".dat"
          open(77, file=filename)
            do i_freq = 1, nomega, 1
              write(77,*) omega(i_freq), &
                       real(green_fn_freq(i_basis,i_basis,i_freq,1)),&
                      aimag(green_fn_freq(i_basis,i_basis,i_freq,1))
            enddo
          close(77)
        enddo
      endif
      endif

      !TEST SPLINE START

!      i_index = 0
!      do j_basis = 1, n_basis, 1
!        do k_basis = 1, j_basis, 1
!         i_index = i_index+1
!         full_ovlp (j_basis , k_basis) = overlap_matrix (i_index)
!         full_ovlp (k_basis , j_basis) = full_ovlp (j_basis , k_basis)
!        enddo
!      enddo

!      t_0 = 0.02
!      h = 1./real(nomega-1)*dlog(omegamax/t_0)
!      if(.true.)then
!      aux_g_tau (:) = 0.d0
!     aux_g     (:) = real  (green_fn_freq(1,1,1:nomega))
!     aux_g_im  (:) = aimag (green_fn_freq(1,1,1:nomega))
!     aux_g_tau(:) = 0.d0
!     tail_g_tau(:) = 0.d0
!     spl_param    (:,:)= 0.d0
!     spl_param_im (:,:)= 0.d0
!
!     call cubic_spline (aux_g (1:nomega),&
!                         nomega, spl_param)
!     call cubic_spline (aux_g_im(1:nomega),&
!                         nomega, spl_param_im)
!
!       ! aux omogeneous grid
!         n_points = 2.d0*omega(nomega-1)*tau(ntau)/2.d0/pi
!        interval = (omega(2)-omega(1))/1.1d0
!         interval = omega(nomega-1) / n_points         
!        n_points = int( omega(nomega-1)/interval)
!         print*, n_points
!         n_points = 600
!         interval = omega(nomega-1) / n_points
!       counter = 1        
!       !if(myid.eq.0) !
!       open(46,file='test_spline.dat')
!
!       do i_freq =1, n_points, 1
!         w1 = interval*(i_freq-1)
!       !   w1 = 1.d0/h*log(w1/t_0+1)
!          if(myid.eq.0) print*, w1
!         if(w1.ge. omega(counter+1)) counter = counter + 1
!
!
!          delta_w = w1-1.d0/h*dlog (omega(counter)/t_0+1.d0 )
!          if(counter.le.1)then 
!             delta_w = (w1- omega(counter+1))/t_0
!          else
!             if(counter.lt.)
!            delta_w = (w1- omega(counter))/(omega(counter+1)-omega(counter))
!          endif
!             delta_w = (w1 - omega(counter+1))/(omega(counter+1)-omega(counter))
!
!         a = spl_param(1,counter) + (0.d0,1.d0) * spl_param_im(1,counter)
!         b = spl_param(2,counter) + (0.d0,1.d0) * spl_param_im(2,counter)
!         c = spl_param(3,counter) + (0.d0,1.d0) * spl_param_im(3,counter)
!         d = spl_param(4,counter) + (0.d0,1.d0) * spl_param_im(4,counter)
!
!         fit = a + b*delta_w + c*delta_w**2 +  d*delta_w**3 
!
!         !if(myid.eq.0) !
!         write(46,*) w1, real(fit)
!
!         do i_tau = 1, ntau
!           tail_g_tau(i_tau) = tail_g_tau(i_tau) +&
!           fit * interval * exp((0.d0,1.d0)*w1*tau(i_tau))/pi
!         enddo
!
!       enddo
!
!       ! if(myid.eq.0) 
!        close(46)
!
!      open(33,file= "FT_spline.dat")
!        do i_tau = 1, ntau, 1
!          write(33, *) tau(-i_tau), &
!           real( tail_g_tau(i_tau))
!        enddo
!      close(33)
!      endif
!! CEDA correction for the empty states
!      do i_freq = 1, nomega , 1
!       do j_basis = 1, n_basis, 1
!        do k_basis = 1, n_basis, 1
!
!         sum_ovlp = 0.d0
!         do i_state = 1, n_states, 1
!          sum_ovlp = sum_ovlp +&
!           ovlp_NAO_KS (i_state,j_basis) *&
!           ovlp_NAO_KS (i_state,k_basis)
!         enddo
!
!         green_fn_freq (j_basis,k_basis,i_freq) = &
!         green_fn_freq (j_basis,k_basis,i_freq) + &
!         (full_ovlp(j_basis,k_basis) - sum_ovlp)* &
!         CMPLX((KS_eigenvalue(n_states,1)+1-e_fermi), omega(i_freq))/ &
!            (omega(i_freq)**2 +(KS_eigenvalue(n_states,1)+1-e_fermi)**2)
!
!        enddo
!       enddo
!      enddo

     
      return
      end subroutine get_green_function_freq_v0
