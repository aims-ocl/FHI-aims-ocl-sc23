!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  get_green_function_time ( n_homo,&
         KS_eigenvalue, KS_eigenvector, green_fn_time, &
         chemical_potential, overlap_matrix, occ_numbers)

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
      use numerical_utilities 
      use prodbas
      use poles_fit 
      use scgw_grid 
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
      real*8  :: green_fn_time (n_basis,n_basis,-ntau:ntau, n_spin)
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,1)
      real*8 chemical_potential
      real*8 overlap_matrix(n_basis*(n_basis+1)/2)
      real*8 occ_numbers(n_states, n_spin)
!EXTRA
      real*8 n_part 
      character*2 iter
      character*17 filename
      logical output 
!COUNTERS 
      integer i_tau
      integer i_spin
      integer i_freq
      integer i_state
      integer j_basis, i_basis
      integer k_basis
      integer i_index, k

!parallelization on the basis
      integer , dimension (:,:), allocatable :: map_index_basis
      integer n_loc_basis
      integer n_remain_basis

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
     
      !CALCULATE G IN TIME AND NAO BASIS
      green_fn_time (:,:,:,:) = 0.d0
      do i_spin = 1, n_spin
        do i_tau= 1, ntau, 1
          do k = 1, n_loc_basis, 1
           j_basis = map_index_basis (myid+1, k)
           if(j_basis .gt.0)then
            do k_basis = 1 , n_basis, 1

              ! negative times and occupied states
              do i_state = 1, n_states, 1
               if(i_state.le. n_homo(i_spin))then 
                 green_fn_time (j_basis,k_basis,-i_tau,i_spin) = &
                 green_fn_time (j_basis,k_basis,-i_tau,i_spin) + &
                 occ_numbers (i_state, i_spin)*real(n_spin)/2.d0*&
                 KS_eigenvector (j_basis,i_state,i_spin,1) *&
                 KS_eigenvector (k_basis,i_state,i_spin,1) *&
                 exp(-(KS_eigenvalue(i_state,i_spin)- chemical_potential)*tau(-i_tau))
                 
                 ! here I calculate the tau=0 point separately
                 if (i_tau.eq.1)then
                    green_fn_time (j_basis,k_basis,0,i_spin) = &
                    green_fn_time (j_basis,k_basis,0,i_spin) + &
                    occ_numbers (i_state, i_spin)*real(n_spin)/2.d0*&
                    KS_eigenvector (j_basis,i_state,i_spin,1) *&
                    KS_eigenvector (k_basis,i_state,i_spin,1)
                 endif
                endif
   
              ! positive times and unoccupied states
                if(i_state.ge. n_homo(i_spin)+1)then
                    green_fn_time (j_basis,k_basis,i_tau,i_spin) = &
                    green_fn_time (j_basis,k_basis,i_tau,i_spin) - &
                    KS_eigenvector (j_basis,i_state,i_spin,1) *&
                    KS_eigenvector (k_basis,i_state,i_spin,1) *&
                    exp(-(KS_eigenvalue(i_state,i_spin)- chemical_potential)*tau(i_tau))
                 endif
              enddo
            enddo
           endif
          enddo

          call sync_matrix (green_fn_time(:,:, i_tau, i_spin ), n_basis, n_basis)
          call sync_matrix (green_fn_time(:,:, -i_tau, i_spin ), n_basis, n_basis)

        enddo
        call sync_matrix (green_fn_time(:,:, 0, i_spin ), n_basis, n_basis)
      enddo


      !calculate the number of particle from the Green's function
      output = .true.
      if(output) then
        n_part = 0
        do i_spin = 1, n_spin
        i_index = 0 
         do i_basis=1, n_basis, 1 
          do j_basis = 1, i_basis, 1
      
             i_index = i_index + 1
      
             if(i_basis .eq. j_basis)then     
               n_part = n_part + green_fn_time (i_basis,j_basis,0,i_spin)*&
               overlap_matrix(i_index)
             else 
               n_part = n_part + 2*green_fn_time (i_basis,j_basis,0,i_spin)*&
               overlap_matrix(i_index)
            
             endif
         
          enddo
         enddo
        enddo
        n_part =(2.d0/n_spin)*n_part
     
        if(myid.eq.0)then
          write(use_unit,*) "Number of particle calculated from the",&
              " Green's function= ",  n_part
        endif
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
!            filename = "green_t_VO_"//iter//".dat"
!            open(77, file=filename)
!              do i_tau = -ntau, ntau, 1
!                write(77,*) tau(i_tau), green_fn_time(i_basis,i_basis,i_tau,1)
!              enddo
!            close(77)
!          enddo
!        endif
!      endif


      return
      end subroutine get_green_function_time

!---------------------------------------------------------------------------------
!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  get_green_function_time_v0 ( ntau, n_homo,&
         KS_eigenvalue, green_fn_time, &
         ovlp_NAO_KS, tau, chemical_potential, inv_overlap_matrix,& 
         overlap_matrix, occ_numbers)

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
      use numerical_utilities 
      use prodbas
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
      integer ntau
      integer  n_homo(n_spin)
      real*8  :: green_fn_time (n_basis,n_basis,-ntau:ntau, n_spin)
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8  :: tau(-ntau:ntau)
      real*8 ovlp_NAO_KS (n_states, n_basis, n_spin)
      real*8 chemical_potential
      real*8 inv_overlap_matrix(n_basis,n_basis)
      real*8 overlap_matrix(n_basis*(n_basis+1)/2)
      real*8 occ_numbers(n_states, n_spin)
      real*8 green_fn (n_basis, n_basis) 
      real*8 gev(n_basis)
!EXTRA
      real*8 e_fermi
      real*8 n_part 
      real*8, dimension(:,:), allocatable ::  green_tmp 
      character*2 iter
      character*17 filename
      logical output 
!COUNTERS 
      integer i_tau
      integer i_spin
      integer i_freq
      integer i_state
      integer j_basis, i_basis
      integer k_basis
      integer i_index
      real*8 aux_ovlp (n_basis,n_basis)

!     if(myid.eq.0)print *, occ_numbers

      i_index = 0
      do i_basis = 1, n_basis, 1
        do k_basis = 1 , i_basis , 1
          i_index= i_index+1 
          aux_ovlp (i_basis,k_basis) = overlap_matrix (i_index)
          aux_ovlp (k_basis,i_basis) = overlap_matrix (i_index)
        enddo
      enddo


      e_fermi = chemical_potential!+0.3 
     
      !CALCULATE G IN TIME AND NAO BASIS
      green_fn_time (:,:,:,:) = 0.d0
      do i_spin = 1, n_spin
        do i_tau= 1, ntau, 1
          do j_basis = 1 , n_basis, 1
            do k_basis = 1 , n_basis, 1
   
              ! negative times and occupied states
              !do i_state = 1, n_homo(i_spin), 1
              do i_state = 1, n_states, 1
               if(i_state.le. n_homo(i_spin))then 
                 green_fn_time (j_basis,k_basis,-i_tau,i_spin) = &
                 green_fn_time (j_basis,k_basis,-i_tau,i_spin) + &
                   occ_numbers (i_state, i_spin)*real(n_spin)/2.d0*&
                   ovlp_NAO_KS (i_state,j_basis,i_spin) *&
                   ovlp_NAO_KS (i_state,k_basis,i_spin) *&
                  exp(-(KS_eigenvalue(i_state,i_spin)- chemical_potential)*tau(-i_tau))
   
                 
                 ! here I calculate the tau=0 point separately
                 if (i_tau.eq.1)then
                    green_fn_time (j_basis,k_basis,0,i_spin) = &
                    green_fn_time (j_basis,k_basis,0,i_spin) + &
                   occ_numbers (i_state, i_spin)*real(n_spin)/2.d0*&
                      ovlp_NAO_KS (i_state,j_basis,i_spin) *&
                      ovlp_NAO_KS (i_state,k_basis,i_spin)
                 endif
                endif
   
!              enddo
   
!              do i_state = n_homo(i_spin)+1, n_states, 1  
                
              ! positive times and unoccupied states
                if(i_state.ge. n_homo(i_spin)+1)then
                    green_fn_time (j_basis,k_basis,i_tau,i_spin) = &
                    green_fn_time (j_basis,k_basis,i_tau,i_spin) - &
                   !(1.d0 - real(n_spin)/2.d0*occ_numbers(i_state, i_spin))*&
                      ovlp_NAO_KS (i_state,j_basis,i_spin) *&
                      ovlp_NAO_KS (i_state,k_basis,i_spin) *&
                      exp(-(KS_eigenvalue(i_state,i_spin)- chemical_potential)*tau(i_tau))
                 endif
              enddo
   
!   correction for the empty states
             ! if(.false.)then
             !  do i_state = 1, n_states, 1
             !    green_fn_time (j_basis,k_basis,i_tau) = &
             !    green_fn_time (j_basis,k_basis,i_tau) - &
             !    !aux_ovlp(i_basis , k_basis) * &
             !    ovlp_NAO_KS (i_state,j_basis) *&
             !    ovlp_NAO_KS (i_state,k_basis) *&
             !    exp(-(KS_eigenvalue(n_states,1)+1)*tau(i_tau))
             !  enddo
             ! endif
               
  
  
            enddo
          enddo
        enddo
      enddo


!      green_fn (:,:)  = green_fn_time(:,:,0,1)*2
!      call diagonalize_rmatrix (n_basis, green_fn, gev, .false.)
!          if(myid.eq.0)then
!            do i_basis =1, n_basis, 1
!              print *,  i_basis, gev (i_basis)
!            enddo
!          endif

      ! now I multiply the Green's function on both sides
      ! with the inverse overlap matrix s_ij^-1. With this procedure
      ! we sort of "renormalize " the Green's function
      !  H_ij  Sum_i'j' (s_ii')^-1 G_i'j' (s_j'j)^-1
 
      if(.not.allocated(green_tmp))then
         allocate(green_tmp (n_basis,n_basis))
      endif
      
      do i_spin = 1, n_spin
        do i_tau = -ntau, ntau, 1
          green_tmp = 0.d0
          call dgemm ('N','N',n_basis,n_basis, n_basis,&
             1.d0, green_fn_time (1,1,i_tau,i_spin), n_basis, &
             inv_overlap_matrix, n_basis, 0.d0,& 
             green_tmp, n_basis)
          green_fn_time (:,:,i_tau,i_spin) = 0.d0
          call dgemm ('N','N',n_basis,n_basis, n_basis,&
             1.d0, inv_overlap_matrix, n_basis, &
             green_tmp, n_basis, 0.d0, &
             green_fn_time (1,1,i_tau,i_spin), n_basis)
        enddo
      enddo
      if(allocated(green_tmp))then
         deallocate(green_tmp)
      endif

      !calculate the number of particle from the Green's function

      output = .true.
      if(output) then
        n_part = 0
        do i_spin = 1, n_spin
        i_index = 0 
         do i_basis=1, n_basis, 1 
          do j_basis = 1, i_basis, 1
      
             i_index = i_index + 1
      
             if(i_basis .eq. j_basis)then     
               n_part = n_part + green_fn_time (i_basis,j_basis,0,i_spin)*&
               overlap_matrix(i_index)
             else 
               n_part = n_part + 2*green_fn_time (i_basis,j_basis,0,i_spin)*&
               overlap_matrix(i_index)
            
             endif
         
          enddo
         enddo
        enddo
        n_part =(2.d0/n_spin)*n_part
     
        if(myid.eq.0)then
          write(use_unit,*) "Number of particle calculated from the",&
              " Green's function= ",  n_part
        endif
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
             filename = "green_t_VO_"//iter//".dat"
             open(77, file=filename)
               do i_tau = -ntau, ntau, 1
                 write(77,*) tau(i_tau), green_fn_time(i_basis,i_basis,i_tau,1)
               enddo
             close(77)
        enddo
      endif
      endif


      return
      end subroutine get_green_function_time_v0

!----------------------------------------------------------------
!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  get_green_function_time_v1 ( n_homo,&
         KS_eigenvalue, KS_eigenvector, green_fn_time, &
         chemical_potential, overlap_matrix, occ_numbers)

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
      use numerical_utilities 
      use prodbas
      use poles_fit 
      use scgw_grid 
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
      real*8  :: green_fn_time (n_basis,n_basis,-ntau:ntau, n_spin)
      real*8  :: KS_eigenvalue(n_states,n_spin)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,1)
      real*8 chemical_potential
      real*8 overlap_matrix(n_basis*(n_basis+1)/2)
      real*8 occ_numbers(n_states, n_spin)
!EXTRA
      real*8 n_part 
      character*2 iter
      character*17 filename
      logical output 
!COUNTERS 
      integer i_tau
      integer i_spin
      integer i_freq
      integer i_state
      integer j_basis, i_basis
      integer k_basis
      integer i_index, k

!parallelization on the basis
      integer , dimension (:,:), allocatable :: map_index_basis
      integer n_loc_basis
      integer n_remain_basis

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
     
      !CALCULATE G IN TIME AND NAO BASIS
      green_fn_time (:,:,:,:) = 0.d0
      do i_spin = 1, n_spin
        do i_tau= 1, ntau, 1
          do k = 1, n_loc_basis, 1
           j_basis = map_index_basis (myid+1, k)
           if(j_basis .gt.0)then
            do k_basis = 1 , n_basis, 1

              ! negative times and occupied states
              do i_state = 1, n_states, 1
               if(i_state.le. n_homo(i_spin))then 
                 green_fn_time (j_basis,k_basis,-i_tau,i_spin) = &
                 green_fn_time (j_basis,k_basis,-i_tau,i_spin) + &
                 occ_numbers (i_state, i_spin)*&
!                 real(n_spin)/2.d0*&
                 KS_eigenvector (j_basis,i_state,i_spin,1) *&
                 KS_eigenvector (k_basis,i_state,i_spin,1) *&
                 exp(-(KS_eigenvalue(i_state,i_spin)- chemical_potential)*tau(-i_tau))
                 
                 ! here I calculate the tau=0 point separately
                 if (i_tau.eq.1)then
                    green_fn_time (j_basis,k_basis,0,i_spin) = &
                    green_fn_time (j_basis,k_basis,0,i_spin) + &
                    occ_numbers (i_state, i_spin)*&
!                    real(n_spin)/2.d0*&
                    KS_eigenvector (j_basis,i_state,i_spin,1) *&
                    KS_eigenvector (k_basis,i_state,i_spin,1)
                 endif
                endif
   
              ! positive times and unoccupied states
                if(i_state.ge. n_homo(i_spin)+1)then
                    green_fn_time (j_basis,k_basis,i_tau,i_spin) = &
                    green_fn_time (j_basis,k_basis,i_tau,i_spin) - &
                    KS_eigenvector (j_basis,i_state,i_spin,1) *&
                    KS_eigenvector (k_basis,i_state,i_spin,1) *&
                    exp(-(KS_eigenvalue(i_state,i_spin)- chemical_potential)*tau(i_tau))
                 endif
              enddo
            enddo
           endif
          enddo

          call sync_matrix (green_fn_time(:,:, i_tau, i_spin ), n_basis, n_basis)
          call sync_matrix (green_fn_time(:,:, -i_tau, i_spin ), n_basis, n_basis)

        enddo
        call sync_matrix (green_fn_time(:,:, 0, i_spin ), n_basis, n_basis)
      enddo


      !calculate the number of particle from the Green's function
      output = .true.
      if(output) then
        n_part = 0
        do i_spin = 1, n_spin
        i_index = 0 
         do i_basis=1, n_basis, 1 
          do j_basis = 1, i_basis, 1
      
             i_index = i_index + 1
      
             if(i_basis .eq. j_basis)then     
               n_part = n_part + green_fn_time (i_basis,j_basis,0,i_spin)*&
               overlap_matrix(i_index)
             else 
               n_part = n_part + 2*green_fn_time (i_basis,j_basis,0,i_spin)*&
               overlap_matrix(i_index)
            
             endif
         
          enddo
         enddo
        enddo
        n_part =(2.d0/n_spin)*n_part
     
        if(myid.eq.0)then
          write(use_unit,*) "Number of particle calculated from the",&
              " Green's function= ",  n_part
        endif
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
!            filename = "green_t_VO_"//iter//".dat"
!            open(77, file=filename)
!              do i_tau = -ntau, ntau, 1
!                write(77,*) tau(i_tau), green_fn_time(i_basis,i_basis,i_tau,1)
!              enddo
!            close(77)
!          enddo
!        endif
!      endif


      return
      end subroutine get_green_function_time_v1


