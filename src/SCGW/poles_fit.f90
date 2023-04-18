!****s* FHI-aims/extrapolar

      module poles_fit
      

!      use scgw_grid

      implicit none
      integer :: n_poles = 100!nomega/2
      real*8 ::  pole_min = 0.01
      real*8 ::  pole_max = 100
      integer :: scgw_it_limit = 40 
      real*8 :: scgw_scf_accuracy = 1.d-5
      real*8 ::  first_state_included = 0.d0
!      real*8 ::  fraction_p1 = 1
      integer :: n_poles_W = 100
      integer :: selfe_matr_el = 1
      real*8 :: exchange_energy_0 = 0.d0
      real*8 :: scgw_mix_param = 0.d0
      logical :: scgw_print_all_spectrum = .false.
      real*8 :: scgw_chem_pot = 0.d0
      logical :: write_G = .false.
      logical :: read_G = .false.

      integer :: pole_dist_type = 0
!      logical :: scgw_print_all_spectrum = .false.
      contains


!--------------------------------------------------------------
!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  transform_sigma_diag (&
         self_energy, mat_dim1,mat_dim2 )
          
      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use scgw_grid
      use physics
      use numerical_utilities
      use localorb_io, only: use_unit
     ! use gw_para
!ARGUMENTS
      implicit none

      INTERFACE
         FUNCTION exp_f (p,freq)
           real*8 exp_f
           real*8, INTENT(IN) :: freq
           real*8, INTENT(IN) :: p
         END FUNCTION exp_f
      END INTERFACE

      INTERFACE
         FUNCTION pole_f (p,freq)
           complex*16 pole_f
           real*8, INTENT(IN) :: freq
           real*8, INTENT(IN) :: p
         END FUNCTION pole_f
      END INTERFACE
 
      complex*16 self_energy (mat_dim1,n_freq)

!      integer n_freq
      integer  mat_dim1
      integer  mat_dim2
      complex*16 matrix(mat_dim1,mat_dim2, n_freq)
      real*8 re_matrix(mat_dim1,mat_dim2, n_freq)
      real*8 im_matrix(mat_dim1,mat_dim2, n_freq)
!      complex*16 aux_matrix_FT(mat_dim1,mat_dim2, -ntau:ntau)
      
!n      real*8 matrix_FT(mat_dim1,mat_dim2, -ntau:ntau)
      character*2 iter       
!EXTRA
      character*17 filename
      logical output
 
!COUNTERS 
      integer i_tau
      integer i_freq
      integer i_state
      integer j_state
      integer j_basis, i_basis

      complex*16 green_KS_poles (mat_dim1, mat_dim2, n_freq)
!      integer n_poles 
      integer i_pole, j_pole
      real*8, dimension(:,:,:), allocatable::  re_G_coeff !(mat_dim1, mat_dim2,n_poles)
      real*8, dimension(:,:,:), allocatable::  im_G_coeff !(mat_dim1, mat_dim2,n_poles)
      real*8, dimension(:,:), allocatable::  re_A!(n_freq,n_poles)
      real*8, dimension(:,:), allocatable::  im_A!(n_freq,n_poles)
      real*8, dimension(:), allocatable:: poles!(n_poles)
      
      real*8 largest_error
      real*8 integ_f
      real*8 alpha
      integer i
  
      !for LSQ fit
     
!      complex*16 A(n_freq,n_poles)
      character*(*), parameter :: func = 'W_fit_coeffs'      
      integer rank
      real*8 highest_pole
      integer i_spin 

!definition of the pole basis

      matrix(:,:,:) = 0.d0
      do i_state = 1, mat_dim1
          matrix(i_state,i_state,:) = self_energy(i_state,:)
      enddo

!     if(myid.eq.0)then
!        write(use_unit,*)"Fourier Transform G(iw) -> G(it)"
!      endif

      !n_poles = 50
      if(.not.allocated(re_G_coeff))then
         allocate(re_G_coeff(mat_dim1, mat_dim2,n_poles))
      endif
      if(.not.allocated(im_G_coeff))then
         allocate(im_G_coeff(mat_dim1, mat_dim2,n_poles))
      endif
      if(.not.allocated(re_A))then
         allocate(re_A (n_freq,n_poles))
      endif
      if(.not.allocated(im_A))then
         allocate(im_A (n_freq,n_poles))
      endif
      if(.not.allocated(poles))then
         allocate(poles(n_poles))
      endif
  
      poles(:) = 0.d0
      re_G_coeff(:,:,:) = 0.d0
      im_G_coeff(:,:,:) = 0.d0
      re_matrix (:,:,:) = real(matrix (:,:,:))
      im_matrix (:,:,:) = aimag(matrix (:,:,:))
      

      poles(1) = 0.01
      highest_pole = 50
      alpha = 1.d0/real(n_poles)*log ((highest_pole-poles(1))/poles(1))
      do i_pole = 2 , n_poles
        poles (i_pole) = poles(1)*(exp(i_pole*alpha))
      enddo
       
      re_A(:,:) = 0.d0
      im_A(:,:) = 0.d0
      do i_freq= 1, n_freq 
        do i_pole = 1, n_poles
          re_A(i_freq,i_pole) =  real(pole_f(poles(i_pole),omega_grid(i_freq)))
          im_A(i_freq,i_pole) = aimag(pole_f(poles(i_pole),omega_grid(i_freq)))
        enddo
      enddo  

      largest_error = 0.d0
      do i_state = 1, mat_dim1, 1
          j_state = i_state

           call solve_LSQ (func, n_freq, n_poles, re_A, &
            re_matrix (i_state,j_state,:),re_G_coeff(i_state,j_state,:) , rank)
           call solve_LSQ (func, n_freq, n_poles, im_A, &
            im_matrix (i_state,j_state,:),im_G_coeff(i_state,j_state,:) , rank)

           integ_f = 0.d0
            do i_freq =1, n_freq, 1
              integ_f = integ_f + abs(matrix (i_state,j_state,i_freq))*womega(i_freq)
            enddo

            if(integ_f.gt.1.d-25)then
            
             green_KS_poles(i_state,j_state,:) = (0.d0,0.d0)
             do i_pole = 1, n_poles
                do i_freq = 1, n_freq, 1
                   green_KS_poles (i_state,j_state, i_freq) =& 
                   green_KS_poles (i_state,j_state, i_freq) +&
                   re_G_coeff (i_state,j_state,i_pole) *&
                   real(pole_f(poles(i_pole),omega_grid(i_freq)))+&
                   (0.d0,1.d0)*im_G_coeff (i_state,j_state,i_pole) *&
                   aimag(pole_f(poles(i_pole),omega_grid(i_freq)))
                enddo
              
              
             enddo

             integ_f = 0.d0
             do i_freq = 1, n_freq, 1
               integ_f = integ_f+abs((green_KS_poles(i_state,j_state,i_freq))-&
                             (matrix(i_state,j_state,i_freq)))*womega(i_freq)
             enddo 
             largest_error = max(largest_error,dabs(integ_f))
             
             output = .false.
             if(myid.eq.0)then
               if(output.and.j_state.eq.i_state.and.i_state .lt.min(mat_dim1,mat_dim2)) then
                   if( i_state.lt.10 ) then
                     write(iter,'(A,I1)') "0",i_state
                   else
                     write(iter,'(I2)') i_state
                   endif
                   filename = "poles_GS_"//iter//".dat"
                   open(77, file=filename)
                   do i_freq = 1, n_freq, 1
                     write(77,*) omega_grid(i_freq), &
                             aimag(green_KS_poles(i_state,j_state,i_freq)),&
                             real(matrix(i_state,j_state,i_freq))
                   enddo
                   close(77)
               endif
             endif
            endif !integral > threshold

!       enddo! j_state
      enddo! i_state

      
      if(allocated(re_G_coeff))then
         deallocate(re_G_coeff)
      endif
      if(allocated(im_G_coeff))then
         deallocate(im_G_coeff)
      endif
      if(allocated(re_A))then
         deallocate(re_A )
      endif
      if(allocated(im_A))then
         deallocate(im_A )
      endif
      if(allocated(poles))then
         deallocate(poles)
      endif

      if(myid.eq.0) write(use_unit,'(A,E14.3)') &
         "          | Fit accuracy for G(w)                ", largest_error 

      end subroutine transform_sigma_diag



!--------------------------------------------------------------
!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  transform_G (&
         matrix, mat_dim1, mat_dim2, matrix_FT )
          
      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use scgw_grid
      use physics
      use numerical_utilities
      use localorb_io, only: use_unit
!ARGUMENTS
      implicit none

!------------------------------------------------------------------------------
! THESE FUNCTIONS ARE USED FOR THE FIT (SEE SCGW/POLES.F90)
      INTERFACE
         FUNCTION exp_f (p,freq)
           real*8 exp_f
           real*8, INTENT(IN) :: freq
           real*8, INTENT(IN) :: p
         END FUNCTION exp_f
      END INTERFACE

      INTERFACE
         FUNCTION pole_f (p,freq)
           complex*16 pole_f
           real*8, INTENT(IN) :: freq
           real*8, INTENT(IN) :: p
         END FUNCTION pole_f
      END INTERFACE
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! INPUT-OUTPUT
      integer  mat_dim1
      integer  mat_dim2
      complex*16 matrix(mat_dim1,mat_dim2, nomega)    !input matrix 
      real*8 matrix_FT(mat_dim1,mat_dim2, -ntau:ntau) !output (FT of the input)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!AUXILIARY STUFF
      real*8, dimension(:), allocatable::  re_G_coeff 
      real*8, dimension(:), allocatable::  im_G_coeff 
      real*8, dimension(:,:), allocatable::  re_A
      real*8, dimension(:,:), allocatable::  im_A
      real*8, dimension(:), allocatable:: poles
      complex*16 green_KS_poles (nomega)
      complex*16 aux_matrix_FT(-ntau:ntau)

      character*17 filename
      character*2 iter       
      logical output
 
!COUNTERS 
      integer i_tau
      integer i_freq
      integer i_state
      integer j_state
      integer j_basis, i_basis
      integer i_pole, j_pole
      integer i
      
      real*8 largest_error
      real*8 integ_f
      real*8 alpha
  
!   for LSQ fit----------------------------------------
      character*(*), parameter :: func = 'W_fit_coeffs'      
      integer rank
      real*8 highest_pole
      integer i_spin 

!  parallelization on the basis----------------------------------------
      integer , dimension (:,:), allocatable :: map_index_basis
      integer n_loc_basis
      integer n_remain_basis
      integer k
!----------------------------------------------------------------------

      !initialize parallelization------------------------------------ 
      n_remain_basis = MOD(mat_dim1, n_tasks)
      if (n_remain_basis.eq.0) then
        n_loc_basis = mat_dim1 / n_tasks
      else
        n_loc_basis = mat_dim1 / n_tasks + 1
      endif

      if (.not. allocated(map_index_basis))then
        allocate(map_index_basis(n_tasks,n_loc_basis ))
      endif

      call distribute_grid (mat_dim1, map_index_basis,n_loc_basis)
      !--------------------------------------------------------

      if(.not.allocated(re_G_coeff))then
         allocate(re_G_coeff(n_poles))
      endif
      if(.not.allocated(im_G_coeff))then
         allocate(im_G_coeff(n_poles))
      endif
      if(.not.allocated(re_A))then
         allocate(re_A (nomega,n_poles))
      endif
      if(.not.allocated(im_A))then
         allocate(im_A (nomega,n_poles))
      endif
      if(.not.allocated(poles))then
         allocate(poles(n_poles))
      endif
 
      !---------------------------------------------------------------
      !defines here the pole basis used to compute the FT analytically 
      poles(:) = 0.d0
      poles(1) = pole_min
      highest_pole = pole_max
      alpha = 1.d0/real(n_poles)*log ((highest_pole-poles(1))/poles(1))


      if (pole_dist_type == 0) then
        do i_pole = 2 , n_poles
          poles (i_pole) = poles(1)*(exp(i_pole*alpha))
        enddo
      elseif  (pole_dist_type == 1 .or. pole_dist_type == 2)then
        !testing for Cu_2
        do i_pole = 1 , n_poles-10
          poles (i_pole) = poles(1)*(exp(i_pole*alpha))
        enddo
        do i_pole = n_poles-9, n_poles
          poles (i_pole) = poles(i_pole-1)*1.5
!        if(myid.eq.0) write(use_unit,*) i_pole,  poles (i_pole), poles (i_pole)*hartree
        enddo

      elseif  (pole_dist_type == 2)then       
        do i_pole = 1 , 10
          poles (i_pole) = 0.001 * i_pole
          poles (i_pole + 10) = 0.01 * i_pole
          poles (i_pole + 20) = 0.1 * i_pole
          poles (i_pole + 30) = 1 * i_pole
          poles (i_pole + 40) = 10 * i_pole
          poles (i_pole + 50) = 100 * i_pole
          poles (i_pole + 60) = 1000 * i_pole
          poles (i_pole + 70) = poles (i_pole + 10) + 0.005
          poles (i_pole + 80) = poles (i_pole + 20) + 0.05
          poles (i_pole + 90) = poles (i_pole + 30) + 0.5
        enddo

      endif

           
      !---------------------------------------------------------------
      ! define the matrix of the residuals for the (non-linear) Least square fit
      re_A(:,:) = 0.d0
      im_A(:,:) = 0.d0
      do i_freq= 1, nomega 
        do i_pole = 1, n_poles
          re_A(i_freq,i_pole) =  real(pole_f(poles(i_pole),omega(i_freq)))
          im_A(i_freq,i_pole) = aimag(pole_f(poles(i_pole),omega(i_freq)))
        enddo
      enddo  
      !---------------------------------------------------------------
         

      matrix_FT (:,:,:) = 0.d0
      aux_matrix_FT (:) = 0.d0
      largest_error = 0.d0
      do k = 1, n_loc_basis, 1
        i_state = map_index_basis (myid+1, k)
        if(i_state .gt.0)then
          do j_state = 1, mat_dim2, 1

            !----------------------------------------------------
            !solve the least square fit
            re_G_coeff(:) = 0.d0
            im_G_coeff(:) = 0.d0
            call solve_LSQ (func, nomega, n_poles, re_A, &
             real(matrix (i_state,j_state,:)),re_G_coeff(:) , rank)
            call solve_LSQ (func, nomega, n_poles, im_A, &
             aimag(matrix (i_state,j_state,:)),im_G_coeff(:) , rank)
            !----------------------------------------------------
            
            integ_f = 0.d0
            do i_freq =1, nomega, 1
              integ_f = integ_f + abs(matrix (i_state,j_state,i_freq))*womega1(i_freq)
            enddo
            
            if(integ_f.gt.1.d-25)then
            
              green_KS_poles(:) = (0.d0,0.d0)
              aux_matrix_FT (:) = (0.d0,0.d0)
              do i_pole = 1, n_poles
              !----------------------------------------------------
              ! this is only to check the accuracy of the fit            
                 do i_freq = 1, nomega, 1
                    green_KS_poles (i_freq) =& 
                    green_KS_poles (i_freq) +&
                    re_G_coeff (i_pole) *&
                    real(pole_f(poles(i_pole),omega(i_freq)))+&
                    (0.d0,1.d0)*im_G_coeff (i_pole) *&
                    aimag(pole_f(poles(i_pole),omega(i_freq)))
                 enddo
              !----------------------------------------------------
               
              !----------------------------------------------------
              ! evaluate the FT analytically !!!
                 do i_freq = 0, ntau, 1
                     aux_matrix_FT (-i_freq) =&
                     aux_matrix_FT (-i_freq) + &
                     (re_G_coeff (i_pole) + &
                      im_G_coeff (i_pole) )*&
                      exp_f(abs(poles(i_pole)),tau(i_freq))/2
               
                    if(i_freq.gt.0)then
                     aux_matrix_FT (i_freq) =&
                     aux_matrix_FT (i_freq) + &
                     (re_G_coeff (i_pole) - &
                      im_G_coeff (i_pole) )*&
                      exp_f(abs(poles(i_pole)),tau(i_freq))/2
                    endif
                 enddo
              !----------------------------------------------------
              enddo
              matrix_FT (i_state,j_state,:) = real(aux_matrix_FT(:))
            
              !----------------------------------------------------
              ! estimate the error of the fitted function with the original function
              integ_f = 0.d0
              do i_freq = 1, nomega, 1
                integ_f = integ_f+abs((green_KS_poles(i_freq))-&
                              (matrix(i_state,j_state,i_freq)))*womega1(i_freq)
              enddo 
              largest_error = max(largest_error,dabs(integ_f))
              !----------------------------------------------------
               
!               output = .false.
!               if(myid.eq.0)then
!                 if(output.and.j_state.eq.i_state.and.i_state .lt.min(mat_dim1,mat_dim2)) then
!                     if( i_state.lt.10 ) then
!                       write(iter,'(A,I1)') "0",i_state
!                     else
!                       write(iter,'(I2)') i_state
!                     endif
!                     filename = "polesfitG"//iter//".dat"
!                     open(77, file=filename)
!                     filename = "FT__polG_"//iter//".dat"
!                     open(78, file=filename)
!                     do i_freq = 1, nomega, 1
!                       write(77,*) omega(i_freq), &
!                               aimag(green_KS_poles(i_state,j_state,i_freq)),&
!                               real(matrix(i_state,j_state,i_freq))
!                     enddo
!                     do i_freq = -ntau, ntau, 1
!                       write(78,*) tau(i_freq), &
!                               real(aux_matrix_FT(i_state,j_state,i_freq)),&
!                               aimag(aux_matrix_FT(i_state,j_state,i_freq))
!                     enddo
!                     close(77)
!                     close(78)
!                 endif
!               endif
            endif !integral > threshold

          enddo! j_state
        endif ! i_state .gt.0
      enddo! i_state
     
      do i_tau = -ntau, ntau, 1
       call sync_matrix (matrix_FT(1,1,i_tau), mat_dim1, mat_dim2)
      enddo  
       
      
      if(allocated(re_G_coeff))then
         deallocate(re_G_coeff)
      endif
      if(allocated(im_G_coeff))then
         deallocate(im_G_coeff)
      endif
      if(allocated(re_A))then
         deallocate(re_A )
      endif
      if(allocated(im_A))then
         deallocate(im_A )
      endif
      if(allocated(poles))then
         deallocate(poles)
      endif

      if(myid.eq.0) write(use_unit,'(A,E14.3)') &
         "          | Fit accuracy for G(w)                ", largest_error 

      end subroutine transform_G

!--------------------------------------------------------------
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  transform_sigma (&
         matrix, mat_dim1, mat_dim2,  &
         matrix_FT) 
          
      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use scgw_grid
      use physics
      use numerical_utilities 
      use localorb_io, only: use_unit
!ARGUMENTS
      implicit none

!------------------------------------------------------------------------------
! THESE FUNCTIONS ARE USED FOR THE FIT (SEE SCGW/POLES.F90)
      INTERFACE
         FUNCTION exp_f (p,time)
           real*8 exp_f
           REAL*8, INTENT(IN) :: time
           REAL*8, INTENT(IN) :: p
         END FUNCTION exp_f
      END INTERFACE

      INTERFACE
         FUNCTION pole_f_cc (p,freq)
           complex*16 pole_f_cc
           REAL*8, INTENT(IN) :: freq
           REAL*8, INTENT(IN) :: p
         END FUNCTION pole_f_cc
      END INTERFACE

      INTERFACE
         FUNCTION pole_f (p,freq)
           complex*16 pole_f
           REAL*8, INTENT(IN) :: freq
           REAL*8, INTENT(IN) :: p
         END FUNCTION pole_f
      END INTERFACE
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! INPUT-OUTPUT
      integer  mat_dim1
      integer  mat_dim2
      real*8 matrix(mat_dim1,mat_dim2, -ntau:ntau)
      complex*16 matrix_FT (mat_dim1,mat_dim2, nomega)
!------------------------------------------------------------------------------
! AUXILIARY MATRICES
      complex*16 green_KS_poles (-ntau:ntau)
      real*8, dimension (:) , allocatable :: G_coeff 
      real*8, dimension (:) , allocatable :: G_coeff_min
      real*8, dimension (:) , allocatable :: poles
      real*8 aux_matrix (1:ntau)
 
!EXTRA
      character*2 iter       
      character*17 filename
      logical output
 
!COUNTERS 
      integer i_tau
      integer i_freq, i_freq2
      integer i_state
      integer j_state
      integer i_spin
      integer j_basis, i_basis
      integer i_pole, j_pole

      real*8 largest_error
      real*8 integ_f
      real*8 alpha
      integer i
      real*8 highest_pole 

!---------------------------------------------------------
      !for LSQ fit
!      real*8 aux_matrix (mat_dim1,mat_dim2, 1:ntau)
      real*8, dimension (:,:) , allocatable :: A
      real*8, dimension (:,:) , allocatable :: A_min
      character*(*), parameter :: func = 'W_fit_coeffs'      
      integer rank

!---------------------------------------------------------
!      parallelization on the basis
      integer , dimension (:,:), allocatable :: map_index_basis
      integer n_loc_basis
      integer n_remain_basis
      integer k
!--------------------------------------------------------------------------------

      !init parallelization------------------------------------ 
      n_remain_basis = MOD(mat_dim1, n_tasks)
      if (n_remain_basis.eq.0) then
        n_loc_basis = mat_dim1 / n_tasks
      else
        n_loc_basis = mat_dim1 / n_tasks + 1
      endif

      if (.not. allocated(map_index_basis))then
        allocate(map_index_basis(n_tasks,n_loc_basis ))
      endif

      call distribute_grid (mat_dim1, map_index_basis,n_loc_basis)
      !--------------------------------------------------------
 
      !--------------------------------------------------------
      allocate(poles(n_poles))
      allocate(G_coeff(n_poles))
      allocate(G_coeff_min(n_poles))
      allocate(A(1:ntau,n_poles))

      !--------------------------------------------------------
      !defines here the pole basis used to compute the FT analytically
      poles(:) = 0.d0
      matrix_FT (:,:,:) = 0.d0

      poles(1) = pole_min
      highest_pole = pole_max
      alpha = 1.d0/real(n_poles)*log ((highest_pole-poles(1))/poles(1))
      if (pole_dist_type == 0) then
        do i_pole = 2 , n_poles
          poles (i_pole) = poles(1)*(exp(i_pole*alpha))
        enddo
      elseif  (pole_dist_type == 1)then
        !testing for Cu_2
        do i_pole = 1 , n_poles-10
          poles (i_pole) = poles(1)*(exp(i_pole*alpha))
        enddo
        do i_pole = n_poles-9, n_poles
          poles (i_pole) = poles(i_pole-1)*1.5
!        if(myid.eq.0) write(use_unit,*) i_pole,  poles (i_pole), poles (i_pole)*hartree
        enddo
      elseif  (pole_dist_type == 2)then
        do i_pole = 1 , 10
          poles (i_pole) = 0.001 * i_pole
          poles (i_pole + 10) = 0.01 * i_pole
          poles (i_pole + 20) = 0.1 * i_pole
          poles (i_pole + 30) = 1 * i_pole
          poles (i_pole + 40) = 10 * i_pole
          poles (i_pole + 50) = 100 * i_pole
          poles (i_pole + 60) = 1000 * i_pole
          poles (i_pole + 70) = poles (i_pole + 10) + 0.005
          poles (i_pole + 80) = poles (i_pole + 20) + 0.05
          poles (i_pole + 90) = poles (i_pole + 30) + 0.5
        enddo

      endif


      !--------------------------------------------------------
      !calculate the matrix of the coeff of the linear sistem for LSQ fit
      A(:,:) = 0.d0
      do i_tau = 1, ntau 
        do i_pole = 1, n_poles
          A(i_tau,i_pole) =  exp_f(poles(i_pole),tau(i_tau))
        enddo
      enddo 
      !--------------------------------------------------------

      largest_error = 0.d0
      do k = 1, n_loc_basis, 1
        i_state = map_index_basis (myid+1, k)
        if(i_state .gt.0)then

          do j_state = 1, mat_dim2, 1

          !----------------------------------------------------
          !solve the least square fit           
           G_coeff(:) = 0.d0
           G_coeff_min(:) = 0.d0

           do i_tau =1, ntau, 1
             aux_matrix(i_tau) = matrix(i_state,j_state,-i_tau)
           enddo

           call solve_LSQ(func, ntau, n_poles, A, &
            matrix (i_state,j_state,1:ntau),G_coeff(:) , rank)
           call solve_LSQ(func, ntau, n_poles, A, &
            aux_matrix (1:ntau),G_coeff_min(:) , rank)
          !----------------------------------------------------

           integ_f = 0.d0
           do i_tau =0, ntau, 1
             integ_f = integ_f + abs(matrix (i_state,j_state,i_tau))*wtau(i_tau)
           enddo

           if(integ_f.gt.1.d-25)then
            !----------------------------------------------------
            ! this is only to check the accuracy of the fit 
             green_KS_poles(:) = (0.d0,0.d0)
             do i_pole = 1, n_poles
               do i_freq = 1, ntau, 1
                  green_KS_poles (i_freq) =& 
                  green_KS_poles (i_freq) +&
                  G_coeff (i_pole) *&
                  exp_f(poles(i_pole),tau(i_freq))
             
                  green_KS_poles (-i_freq) =&
                  green_KS_poles (-i_freq) +&
                  G_coeff_min (i_pole) *&
                  exp_f(poles(i_pole),tau(i_freq))
               enddo
             enddo
            !----------------------------------------------------

            !----------------------------------------------------
            ! evaluate the FT analytically !!!        
             do i_pole = 1, n_poles
               do i_freq2 = 1, nomega
                 matrix_FT (i_state,j_state, i_freq2) = &
                 matrix_FT (i_state,j_state, i_freq2) +&
                  (G_coeff (i_pole)* &
                  pole_f_cc(poles(i_pole),omega(i_freq2))+&
                  G_coeff_min (i_pole) *&
                  pole_f(poles(i_pole),omega(i_freq2)))!/pi/2
               enddo
             enddo
            !----------------------------------------------------

            !----------------------------------------------------
            ! estimate the error of the fitted function with the original function
             integ_f = 0.d0
             do i_freq = 1, ntau, 1
               integ_f = integ_f + (real(green_KS_poles(i_freq))-&
                             real(matrix(i_state,j_state,i_freq)))*wtau(i_freq)
             enddo 
             largest_error = max(largest_error,dabs(integ_f))
            !------------------------------------------------------


!         output = .false.
!         if(myid.eq.0)then
!           if(output.and.j_state.eq.i_state.and.i_state .lt.min(mat_dim1,mat_dim2)) then
!               if( i_state.lt.10 ) then
!                 write(iter,'(A,I1)') "0",i_state
!               else
!                 write(iter,'(I2)') i_state
!               endif
!               filename = "sigma__FT"//iter//".dat"
!               open(87, file=filename)
!               filename = "pol_sigma"//iter//".dat"
!               open(77, file=filename)
!               do i_freq = -ntau, ntau, 1
!                 write(77,*) tau(i_freq), &
!                          real(green_KS_poles(i_freq)),&
!                         real(matrix(i_state,j_state,i_freq))
!               enddo
!               do i_freq = 1 ,nomega
!                  write(87,*) omega (i_freq), real(matrix_FT(i_state,j_state,i_freq)),&
!                         aimag(matrix_FT(i_state,j_state,i_freq))
!               enddo
!               close(77)
!               close(87)
!           endif
!         endif
        endif !integral > threshold
       enddo! j_state
       endif !i_state .gt.0
      enddo! i_state

      do i_freq = 1, nomega, 1
        call sync_matrix_complex (matrix_FT(:,:,i_freq), mat_dim1, mat_dim2)
      enddo


!      if(myid.eq.0) write(use_unit,*) "Largest error in the fit of Sigma", largest_error
            if(myid.eq.0) write(use_unit,'(A,E14.3)') &
         "          | Fit accuracy for Self-Energy         ", largest_error

      deallocate(G_coeff)
      deallocate(G_coeff_min)
      deallocate(A)
      deallocate(poles)


      end  subroutine  transform_sigma
!--------------------------------------------------------------
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  transform_polar (&
         matrix, mat_dim1, mat_dim2, &!poles, G_coeff, n_poles , &
         matrix_FT) 
          
      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use scgw_grid
      use physics
      use numerical_utilities 
      use localorb_io, only: use_unit
!ARGUMENTS
      implicit none

      INTERFACE
         FUNCTION exp_f (p,time)
           real*8 exp_f
           REAL*8, INTENT(IN) :: time
           REAL*8, INTENT(IN) :: p
         END FUNCTION exp_f
      END INTERFACE

      INTERFACE
         FUNCTION pole_f_re (p,freq)
           real*8 pole_f_re
           REAL*8, INTENT(IN) :: freq
           REAL*8, INTENT(IN) :: p
         END FUNCTION pole_f_re
      END INTERFACE

      integer  mat_dim1
      integer  mat_dim2
      real*8 matrix(mat_dim1,mat_dim2, 0:ntau)
      real*8 matrix_FT (mat_dim1,mat_dim2, nomega)
      character*2 iter       
!EXTRA
      character*17 filename
      logical output
 
!COUNTERS 
      integer i_tau
      integer i_freq
      integer i_state
      integer j_state
      integer j_basis, i_basis

      !test 
      complex*16 green_KS_poles (0:ntau)
!      integer n_poles
      integer i_pole, j_pole
      real*8, dimension(:), allocatable:: G_coeff !(mat_dim1, mat_dim2,n_poles)
      real*8, dimension(:,:), allocatable:: A!(0:ntau,n_poles)
      real*8, dimension(:), allocatable::  poles!(n_poles)
 
      real*8 largest_error
      real*8 integ_f
      real*8 alpha
      integer i
      real*8 highest_pole 
      !for LSQ fit
     
       character*(*), parameter :: func = 'W_fit_coeffs'      
      integer rank
 
!definition of the pole basis

      n_poles_W = n_poles
      allocate(G_coeff(n_poles_W))
      allocate(A(0:ntau,n_poles_W))
      allocate(poles(n_poles_W))

      poles(:) = 0.d0
      matrix_FT (:,:,:) = 0.d0


      !below E_F
      poles(1) = 0.01
      highest_pole =  100
      alpha = 1.d0/real(n_poles_W/2)*log ((highest_pole-poles(1))/poles(1))

      if (pole_dist_type == 0) then
        do i_pole = 2 , n_poles
          poles (i_pole) = poles(1)*(exp(i_pole*alpha))
        enddo
      elseif  (pole_dist_type == 1)then
        !testing for Cu_2
        do i_pole = 1 , n_poles-10
          poles (i_pole) = poles(1)*(exp(i_pole*alpha))
        enddo
        do i_pole = n_poles-9, n_poles
          poles (i_pole) = poles(i_pole-1)*1.5
!        if(myid.eq.0) write(use_unit,*) i_pole,  poles (i_pole), poles (i_pole)*hartree
        enddo
      elseif  (pole_dist_type == 2)then
        do i_pole = 1 , 10
          poles (i_pole) = 0.001 * i_pole
          poles (i_pole + 10) = 0.01 * i_pole
          poles (i_pole + 20) = 0.1 * i_pole
          poles (i_pole + 30) = 1 * i_pole
          poles (i_pole + 40) = 10 * i_pole
          poles (i_pole + 50) = 100 * i_pole
          poles (i_pole + 60) = 1000 * i_pole
          poles (i_pole + 70) = poles (i_pole + 10) + 0.005
          poles (i_pole + 80) = poles (i_pole + 20) + 0.05
          poles (i_pole + 90) = poles (i_pole + 30) + 0.5
        enddo

      endif


      A(:,:) = 0.d0
      do i_freq= 0, ntau 
        do i_pole = 1, n_poles_W
          A(i_freq,i_pole) =  exp_f(poles(i_pole),tau(i_freq))
        enddo
      enddo          

      largest_error = 0.d0
      do i_state = 1, mat_dim1, 1
        do j_state = 1, mat_dim2, 1

           G_coeff(:) = 0.d0
           call solve_LSQ(func, ntau+1, n_poles_W, A, &
            matrix (i_state,j_state,:),G_coeff(:) , rank)

           integ_f = 0.d0
            do i_freq =0, ntau, 1
              integ_f = integ_f + abs(matrix (i_state,j_state,i_freq))*wtau(i_freq)
            enddo
           if(integ_f.gt.1.d-25)then

      

        green_KS_poles(:) = (0.d0,0.d0)
        do i_pole = 1, n_poles_W
          do i_freq = 0, ntau, 1
             green_KS_poles (i_freq) =& 
             green_KS_poles (i_freq) +&
             G_coeff (i_pole) *&
             exp_f(poles(i_pole),tau(i_freq))
          enddo
          do i_freq = 1, nomega
            matrix_FT (i_state,j_state, i_freq) = &
             matrix_FT (i_state,j_state, i_freq) +&
             G_coeff (i_pole) *&
             pole_f_re(poles(i_pole),omega(i_freq))*2
          enddo
        enddo


         integ_f = 0.d0
         do i_freq = 1, ntau, 1
           integ_f = integ_f + (real(green_KS_poles(i_freq))-&
                         real(matrix(i_state,j_state,i_freq)))*wtau(i_freq)
         enddo 
         largest_error = max(largest_error,dabs(integ_f))

!         output = .false.
!         if(myid.eq.0)then
!           if(output.and.j_state.eq.i_state.and.i_state .lt.min(mat_dim1,mat_dim2)) then
!               if( i_state.lt.10 ) then
!                 write(iter,'(A,I1)') "0",i_state
!               else
!                 write(iter,'(I2)') i_state
!               endif
!               filename = "pol_chiFT"//iter//".dat"
!               open(87, file=filename)
!               filename = "poles_chi"//iter//".dat"
!               open(77, file=filename)
!               do i_freq = 0, ntau, 1
!                 write(77,*) tau(i_freq), &
!                          real(green_KS_poles(i_state,j_state,i_freq)),&
!                         real(matrix(i_state,j_state,i_freq))
!               enddo
!               do i_freq = 1 ,nomega
!                  write(87,*) omega (i_freq), real(matrix_FT(i_state,j_state,i_freq))!,& 
!           !            aimag(matrix_FT(i_state,j_state,i_freq))
!               enddo
!               close(77)
!               close(87)
!           endif
!         endif
        endif !integral > threshold
       enddo! j_state
      enddo! i_state
!      if(myid.eq.0) write(use_unit,*) "Largest error in the fit ", largest_error
      if(myid.eq.0) write(use_unit,'(A,E14.3)') &
         "          | Fit accuracy for polarizability      ", largest_error

      deallocate(G_coeff)
      deallocate(A)
      deallocate(poles)

      end subroutine transform_polar

!--------------------------------------------------------------
!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  transform_W (&
         matrix, mat_dim1, mat_dim2,&! poles, G_coeff, n_poles ,&
         matrix_FT ) 
          
      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use scgw_grid
      use physics
      use numerical_utilities 
      use localorb_io, only: use_unit
!ARGUMENTS
      implicit none

      INTERFACE
         FUNCTION exp_f (p,time)
           real*8 exp_f
           REAL*8, INTENT(IN) :: time
           REAL*8, INTENT(IN) :: p
         END FUNCTION exp_f
      END INTERFACE

      INTERFACE
         FUNCTION pole_f_re (p,freq)
           real*8 pole_f_re
           REAL*8, INTENT(IN) :: freq
           REAL*8, INTENT(IN) :: p
         END FUNCTION pole_f_re
      END INTERFACE

      integer  mat_dim1
      integer  mat_dim2
      real*8 matrix(mat_dim1,mat_dim2, nomega)
      real*8 matrix_FT(mat_dim1,mat_dim2, 0:ntau)
      character*2 iter       
!EXTRA
      character*17 filename
      logical output
 
!COUNTERS 
      integer i_tau
      integer i_freq
      integer i_state
      integer j_state
      integer j_basis, i_basis

      !test 
      complex*16 green_KS_poles (nomega)
!      integer n_poles
      integer i_pole, j_pole

      real*8, dimension(:), allocatable:: G_coeff !(mat_dim1, mat_dim2,n_poles)
      real*8, dimension(:,:), allocatable:: A!(0:ntau,n_poles)
      real*8, dimension(:), allocatable::  poles!(n_poles)

 
      real*8 integ_f
      real*8 alpha
      integer i
      real*8 largest_error
      real*8 highest_pole 
      !for LSQ fit
      real*8 t1,t2,W1,W2     
       character*(*), parameter :: func = 'W_fit_coeffs'      
      integer rank
 
!definition of the pole basis
      !n_poles = 100
      n_poles_W = n_poles
      allocate(G_coeff(n_poles_W))
      allocate(A(nomega,n_poles_W))
      allocate(poles(n_poles_W))

      poles(:) = 0.d0

      matrix_FT(:,:,:) = 0.d0
      poles(1) = 0.01
      highest_pole =  50
      alpha = 1.d0/real(n_poles_W/2)*log ((highest_pole-poles(1))/poles(1))

      if (pole_dist_type == 0) then
        do i_pole = 2 , n_poles
          poles (i_pole) = poles(1)*(exp(i_pole*alpha))
        enddo
      elseif  (pole_dist_type == 1 .or. pole_dist_type == 2)then
        !testing for Cu_2
        do i_pole = 1 , 10
          poles (i_pole) = 0.001 * i_pole
          poles (i_pole + 10) = 0.01 * i_pole
          poles (i_pole + 20) = 0.1 * i_pole
          poles (i_pole + 30) = 1 * i_pole
          poles (i_pole + 40) = 10 * i_pole
          poles (i_pole + 50) = 100 * i_pole
          poles (i_pole + 60) = 1000 * i_pole
          poles (i_pole + 70) = poles (i_pole + 10) + 0.005
          poles (i_pole + 80) = poles (i_pole + 20) + 0.05
          poles (i_pole + 90) = poles (i_pole + 30) + 0.5
        enddo
      endif

!      do i_pole = 1 , n_poles-10
!        poles (i_pole) = poles(1)*(exp(i_pole*alpha))
!        if(myid.eq.0) write(use_unit,*) i_pole,  poles (i_pole), poles (i_pole)*hartree
!      enddo
!      do i_pole = n_poles-9, n_poles
!        poles (i_pole) = poles(i_pole-1)*1.5
!        if(myid.eq.0) write(use_unit,*) i_pole,  poles (i_pole), poles (i_pole)*hartree
!      enddo


      A(:,:) = 0.d0
      do i_freq= 1, nomega 
        do i_pole = 1, n_poles_W
          A(i_freq,i_pole) =  pole_f_re(poles(i_pole),omega(i_freq))
        enddo
      enddo          

      largest_error = 0.d0
      do i_state = 1, mat_dim1, 1
        do j_state = 1, mat_dim2, 1

           G_coeff(:) = 0.d0
           call solve_LSQ(func, nomega, n_poles_W, A, &
            matrix (i_state,j_state,:),G_coeff(:) , rank)

           integ_f = 0.d0
            do i_freq =1, nomega, 1
              integ_f = integ_f + abs(matrix (i_state,j_state,i_freq))*womega1(i_freq)
            enddo
           if(integ_f.gt.1.d-25)then

              green_KS_poles(:) = (0.d0,0.d0)
              do i_pole = 1, n_poles_W
                do i_freq = 1, nomega, 1
                   green_KS_poles (i_freq) =& 
                   green_KS_poles (i_freq) +&
                   G_coeff (i_pole) *&
                   pole_f_re(poles(i_pole),omega(i_freq))
                enddo
                do i_freq = 1, ntau
                  matrix_FT (i_state,j_state, i_freq) =&
                  matrix_FT (i_state,j_state, i_freq) +&
                  G_coeff (i_pole) *&
                  exp_f (poles(i_pole),tau(i_freq))*pi
                enddo
                !point t=0 is treated here
                !testing point 0
           W1 = matrix_FT (i_state,j_state, 1)  
           W2 = matrix_FT (i_state,j_state, 2) 
           t1 = tau(1)
           t2 = tau(2)
           matrix_FT (i_state,j_state, 0) = W1-(W2-W1)/(t2-t1)*t1

!                matrix_FT (i_state,j_state, 0) =&
!                  matrix_FT (i_state,j_state, 0) +&
!                  G_coeff (i_state,j_state,i_pole) *pi
              enddo

         integ_f = 0.d0
         do i_freq = 1, nomega, 1
           integ_f = integ_f + (real(green_KS_poles(i_freq))-&
                         real(matrix(i_state,j_state,i_freq)))*womega1(i_freq)
         enddo
         largest_error = max(largest_error,dabs(integ_f)) 
!         if(abs(integ_f) .gt. 1.d-04 ) then 
!              write(use_unit,*) i_state,j_state, integ_f 
!              do i_freq = nomega-30, nomega, 1
!                 matrix_FT(i_state,j_state,i_freq) = matrix_FT(i_state,j_state,i_freq) &
!                 * 1/(omega(i_freq)- omega(nomega-30)+0.1)**2
!              enddo
!         endif    

         output = .false.
         if(myid.eq.0)then
           if(output.and.i_state.eq.1 ) then
           !if(output.and.j_state.eq.i_state.and.i_state  .lt. min(mat_dim1,mat_dim2)) then
               if( i_state.lt.10 ) then
                 write(iter,'(A,I1)') "0",i_state
               else
                 write(iter,'(I2)') i_state
               endif
               filename = "FT_pol_W_"//iter//".dat"
               open(78, file=filename)
               filename = "gree_polW_"//iter//".dat"
               open(77, file=filename)
               do i_freq = 1, nomega, 1
                 write(77,*) omega(i_freq), &
                          real(green_KS_poles(i_freq)),&
                         real(matrix(i_state,j_state,i_freq))
               enddo
               do i_freq = 0, ntau, 1
                 write(78,*) tau(i_freq), &
                         real(matrix_FT(i_state,j_state,i_freq))
               enddo
               close(78)
               close(78)
           endif
         endif
        endif !integral > threshold
       enddo! j_state
      enddo! i_state
!      if(myid.eq.0) write(use_unit,*) "largest error for W ", largest_error
      if(myid.eq.0) write(use_unit,'(A,E14.3)') &
         "          | Fit accuracy for W                   ", largest_error
      deallocate(G_coeff)
      deallocate(A)
      deallocate(poles)


      end subroutine transform_W
!------------------------------------------------------------------
      subroutine get_convolution_matrix (poles, n_poles, convolution_matr)
      !calculates analytically the convolution between two poles
      use mpi_tasks
      use scgw_grid

      INTERFACE
         FUNCTION pole_f (p,freq)
           complex*16 pole_f
           real*8, INTENT(IN) :: freq
           real*8, INTENT(IN) :: p
         END FUNCTION pole_f
      END INTERFACE

      integer n_poles
      real*8  poles(n_poles)
      complex*16 convolution_matr (n_poles, n_poles, nomega)
      complex*16 convolution_matr0 (n_poles, n_poles, nomega)
 
      integer i_freq1, i_freq2
      integer i_pole, j_pole

      convolution_matr(:,:,:)  = 0.d0

      if(.false.)then
        convolution_matr0(:,:,:)  = 0.d0
        do i_freq1 = 1, nomega
          do i_freq2 = 1, nomega
            do i_pole = 1, n_poles
              do j_pole = 1, n_poles
                convolution_matr0(i_pole, j_pole, i_freq1) = &
                convolution_matr0(i_pole, j_pole, i_freq1) + &
                (pole_f(poles(i_pole),omega(i_freq2))*& 
                pole_f(poles(j_pole),(omega(i_freq1)-omega(i_freq2))) &
                + pole_f(poles(i_pole),-omega(i_freq2))*&
                pole_f(poles(j_pole),(omega(i_freq1)+omega(i_freq2))))*&
                womega1(i_freq2)
              enddo
            enddo
          enddo
        enddo
      endif

      do i_freq1 = 1, nomega
        do i_pole = 1, n_poles
          do j_pole = 1, n_poles
             if(.not. (poles(i_pole)+poles(j_pole)).eq.0)then
              convolution_matr(i_pole, j_pole, i_freq1) =&
              convolution_matr(i_pole, j_pole, i_freq1) +& 
              1.d0/(poles(i_pole)+poles(j_pole)+&
              (0.d0,1.d0)*omega(i_freq1))*4*sqrt(abs(poles(i_pole)*poles(j_pole)))
             endif
          enddo
        enddo
      enddo
         
      if(.false.)then
       open (44, file = 'conv.dat')
        do i_pole = 1, n_poles/4
          do i_freq1 = 1, nomega
            if(myid.eq.0)then
               write(44,*) omega(i_freq1), real (convolution_matr(i_pole, 1, i_freq1)) 
            endif
          enddo
        enddo
       close(44) 
      endif

       
      end subroutine get_convolution_matrix 
!------------------------------------------------------------------
      end module
