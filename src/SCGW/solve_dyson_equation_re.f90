!-----------------------------------------------------------------------
      subroutine solve_dyson_equation_re_v2 ( &
          self_energy_freq,      & 
          xc_matr,               &
          hartree_pot,           &
          hartree_pot0,          &
          exchange_self_energy0, &
          exchange_self_energy,  &
          number_reiterations,   &
          scgw_converged,        &
          i_spin)

      use runtime_choices, only: hybrid_coeff
      use constants, only: pi, hartree
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use physics
      use scgw_grid
      use poles_fit

      implicit none
! INPUT
      complex*16  self_energy_freq      (n_basis,n_basis,nomega)
      real*8      xc_matr               (n_basis,n_basis)
      real*8      hartree_pot           (n_basis,n_basis)
      real*8      hartree_pot0          (n_basis,n_basis)
      real*8      exchange_self_energy  (n_basis,n_basis)
      real*8      exchange_self_energy0 (n_basis,n_basis)

      integer number_reiterations
      logical scgw_converged

      !AUXILIARY 
      logical output
      integer i_freq
      integer k_state
      integer i_spin
      integer k 
      integer i_basis, j_basis
      integer i_state, j_state
      character*2 iter, string
      integer i_index       
      complex*16 self_energy_diag (n_states,n_states,nomega)

      !the local grid
      integer aux_nomega
      real*8  aux_omegamax
      real*8  , dimension (:),   allocatable :: aux_omega
      real*8  , dimension (:),   allocatable :: aux_womega
      integer , dimension (:,:), allocatable :: aux_map_index
      integer n_loc_grid1
      integer n_remain

      complex*16, dimension (:,:,:), allocatable :: aux_G 
      complex*16, dimension (:,:,:), allocatable :: aux_G1 
      real*8,     dimension (:),     allocatable :: spectrum
      real*8,     dimension (:,:),   allocatable :: spectrum_single_state
      real*8 freq
      integer sign_freq
      integer anacon_type1
      complex*16, dimension(:,:,:), allocatable:: self_en_par_KS
      integer n_homo
      integer n_max_par1
      complex*16 selfe
      real*8 aux_matr (n_basis, n_basis)
      real*8 aux_matr1 (n_states, n_states)
      character*18 filename
      character*15 filename1
      character*24 filename2
      real*8 g1  
      integer state_to_print_out 
      integer n_states_to_print
      logical analyze_peaks 

      first_state_included = 1
      anacon_type1         = 1
      n_max_par1           = nomega / 2
      aux_nomega           = 6000
      aux_omegamax         = 2.50
      allocate(aux_omega(aux_nomega))
      allocate(aux_womega(aux_nomega))
      allocate(aux_G (n_states,n_states,1))
      allocate(aux_G1(n_states,n_states,1))
      allocate(spectrum(aux_nomega))
      allocate(self_en_par_KS (n_max_par1,n_states,n_states))
   
      n_homo = 0
      do i_state = 1, n_states
        if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
          n_homo=i_state
        endif
      enddo

     !select_state_to_print----------------
      analyze_peaks = .false. !(produces a projected spectrum on a given KS state)
      if(analyze_peaks)then
        n_states_to_print = 10
        if (n_states_to_print .gt. n_homo ) n_states_to_print = n_homo
        
        allocate(spectrum_single_state(n_states_to_print,aux_nomega))
        do k_state=1, n_states_to_print, 1

          if ( n_homo +1- k_state .lt.10)then
            write(string,'(A,I1)') "0", n_homo +1- k_state
          else
            write(string,'(I2)') n_homo +1- k_state
          endif

          filename2 = "spectrum_KS_state_"//string//".dat"
          if( myid.eq.0 ) open(34+k_state,file=filename2)
        enddo
      else
        n_states_to_print = 0
      endif 
     !---------------------------------------------
      

      !parallelization of the auxiliary grid---------------- 
       n_remain = MOD(aux_nomega, n_tasks)
       if (n_remain.eq.0) then
         n_loc_grid1 = aux_nomega / n_tasks
       else
         n_loc_grid1 = aux_nomega / n_tasks + 1
       endif
       allocate(aux_map_index(n_tasks, n_loc_grid1))
       call distribute_grid (aux_nomega, aux_map_index,n_loc_grid1)

      !transform sigma in the KS basis and do the analytic continuation
      self_en_par_KS   (:,:,:) = (0.d0 ,0.d0)
      self_energy_diag (:,:,:) = (0.d0 ,0.d0)
      call gauleg(0d0, aux_omegamax, aux_omega, aux_womega, aux_nomega)
      call diagonalize_self_en(self_energy_freq, self_energy_diag)
      call analy_continue_green_fn &
          (anacon_type1,&
           nomega, &
           n_max_par1, &
           self_en_par_KS,omega, &
           self_energy_diag, n_states)

     !convert other matrices to the KS basis
      aux_matr  (:,:) = (0.d0,0.d0)
      aux_matr1 (:,:) = (0.d0,0.d0)
      if (.not. use_hartree_fock) then
        aux_matr (:,:) =                 &
        + xc_matr                (:,:)   &
        + hartree_pot0           (:,:)   &
        - exchange_self_energy   (:,:)   &
        - hartree_pot            (:,:)   
      else
        aux_matr (:,:) =                                 &
        + xc_matr                                (:,:)   &
        + hartree_pot0                           (:,:)   &
        + hybrid_coeff * exchange_self_energy0   (:,:)   &
        - exchange_self_energy                   (:,:)   &
        - hartree_pot                            (:,:)
      endif 

      do i_basis=1, n_basis, 1
        do j_basis=1, n_basis, 1
          do i_state=1, n_states, 1
            do j_state=1, n_states, 1
             aux_matr1 (i_state, j_state)=& 
             aux_matr1 (i_state, j_state)+& 
             KS_eigenvector (i_basis,i_state,i_spin,1)*&
             aux_matr (i_basis,j_basis)*& 
             KS_eigenvector (j_basis,j_state,i_spin,1)
            enddo
          enddo
        enddo
      enddo

      !set file name of the output 
      if(scgw_converged) then
         if(n_spin .eq. 1 )then
         filename1 = 'spectrum_sc.dat'
         if( myid.eq.0 )open(33,file=filename1)
         else
             if (i_spin ==1)  filename = "spectrum_sc_up.dat"
             if (i_spin ==2)  filename = "spectrum_sc_do.dat"
             if( myid.eq.0 )open(33,file=filename)
         endif
      else
         if (number_reiterations.lt.10)then
           write(iter,'(A,I1)') "0", number_reiterations
         else
           write(iter,'(I2)') number_reiterations
         endif

         if(n_spin .eq. 1 )then
           filename1 = "sp_ImG"//iter//".dat"
         else
             if (i_spin ==1)  filename1 = "sp_ImG"//iter//"_SU.dat"
             if (i_spin ==2)  filename1 = "sp_ImG"//iter//"_SD.dat"
         endif
         if( myid.eq.0 )open(33,file=filename1)
      endif

      !EVALUATE THE SPECTRAL FUNCTION
      do sign_freq = -1, 1, 2 !loop of sign + and - of the frequency axis
        spectrum(:) = 0.d0
        if(analyze_peaks) spectrum_single_state(:,:) = 0.d0

        do k = 1, n_loc_grid1, 1
          i_freq = aux_map_index (myid+1, k)
          if(i_freq.gt.0)then
            freq = sign_freq*aux_omega(i_freq)
         
            aux_G (:,:,:) = (0.d0,0.d0)
            aux_G1(:,:,:) = (0.d0,0.d0)
          
            do i_state = 1, n_states , 1
              do j_state = 1, n_states, 1
      
                call get_real_selfenergy (anacon_type1, nomega, &
                omega, freq, n_max_par1, self_en_par_KS(:,i_state,j_state), selfe)
      
                aux_G (i_state,j_state,1) =&
                aux_G (i_state,j_state,1) &
                + aux_matr1(i_state,j_state)-selfe 
      
                !add G_0 (only on the diagonal)
                if(i_state .eq. j_state)then
                  aux_G (i_state,i_state,1) =&
                  aux_G (i_state,i_state,1)+& 
                   (freq- (KS_eigenvalue(i_state,i_spin,1)-chemical_potential)+&
                   (0.d0,1.d0)*0.0001*&
                    sign(1.d0,KS_eigenvalue(i_state,i_spin,1)-chemical_potential))
                endif
              enddo
            enddo
      
            call invert_green(& 
            aux_G(1:n_states,1:n_states,1), 1,&
            aux_G1(1:n_states,1:n_states,1),n_states) 
      
            do i_state=1 , n_states,1
              do j_state=1, n_states,1     
                if(i_state .eq.j_state)then 
                 spectrum (i_freq) = spectrum(i_freq) +&
                 aimag(aux_G1(i_state,j_state,1))*2./pi
                endif 
              enddo
            enddo
            if(analyze_peaks) then      
              do k_state=1, n_states_to_print, 1 
                spectrum_single_state (k_state,i_freq) = &
                spectrum_single_state (k_state,i_freq) + &
                aimag(aux_G1(n_homo+1-k_state, n_homo+1-k_state,1))*2./pi 
              enddo
            endif
          endif 
        enddo ! loop over freq points

        call sync_vector(spectrum, aux_nomega)
        if(analyze_peaks) then
          do k_state = 1, n_states_to_print, 1
            call sync_vector(spectrum_single_state(k_state,:),aux_nomega)
          enddo
        endif
 
       !write the spectrum to file----------
        if(myid.eq.0)then
          do i_freq = 1, aux_nomega, 1 
            write(33,*) (sign_freq*aux_omega(i_freq)+chemical_potential)*hartree, abs(spectrum(i_freq))

            if(analyze_peaks) then
              do k_state =1, n_states_to_print, 1 
                write(34+k_state,*) (sign_freq*aux_omega(i_freq)+chemical_potential)*hartree, &
                     abs(spectrum_single_state(k_state,i_freq))
              enddo
            endif
          enddo
        endif  
      enddo !end loop over sign 

      if(myid.eq.0)then
        close (33)
        if(analyze_peaks) then
          do k_state =1, n_states_to_print, 1 
            close (34+k_state)
          enddo
        endif
      endif

      deallocate (spectrum)
      if(analyze_peaks) deallocate (spectrum_single_state)

      end subroutine solve_dyson_equation_re_v2
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine solve_dyson_equation_re_v0 (inv_green_fn_freq, &
          self_energy_freq, new_green_fn_freq,&
          exchange_self_energy, xc_matr, hartree_pot, &
!          inv_overlap_matrix, &
          exchange_self_energy0, &
          new_chem_pot , number_reiterations,&
          scgw_converged,i_spin)


      use runtime_choices, only: hybrid_coeff
      use constants, only: pi, hartree
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use physics
      use scgw_grid
      use poles_fit

      implicit none
! INPUT
      !integer nomega
      complex*16  inv_green_fn_freq(n_basis,n_basis,nomega)
      complex*16  self_energy_freq (n_basis,n_basis,nomega)
      !real*8  omega(nomega)
      real*8  xc_matr(n_basis, n_basis)
      real*8  exchange_self_energy(n_basis,n_basis)
      real*8  hartree_pot(n_basis,n_basis)
      real*8  exchange_self_energy0 (n_basis,n_basis)
!      real*8  inv_overlap_matrix (n_basis, n_basis)
      real*8  new_chem_pot
!      real*8 ovlp_NAO_KS (n_states, n_basis)
      integer number_reiterations
      logical scgw_converged
!OUTPUT
      complex*16  new_green_fn_freq(n_basis,n_basis,nomega)

!AUXILIARY 
      logical output
      integer i_freq
      integer k_state
      integer i_spin
      integer k 
      integer i_basis, j_basis
      integer i_state, j_state
      character*2 iter , string
      character*2 matr_el
      integer i_index       
      complex*16  self_energy_diag (n_states,n_states,nomega)

!the local grid
      integer aux_nomega
      real*8  aux_omegamax
      real*8 , dimension (:), allocatable :: aux_omega
      real*8 , dimension (:), allocatable :: aux_womega
      integer , dimension(:,:) , allocatable :: aux_map_index
      integer n_loc_grid1
      integer n_remain

      complex*16, dimension (:,:,:),allocatable :: aux_G 
      complex*16, dimension (:,:,:),allocatable :: aux_G1 
      real*8, dimension(:), allocatable :: spectrum
      real*8, dimension(:,:), allocatable :: spectrum_single_state
      real*8 freq
      integer sign_freq
      integer anacon_type1
      complex*16, dimension(:,:,:), allocatable:: self_en_par_KS
      integer n_homo
      integer n_max_par1
      complex*16 selfe
      real*8 aux_matr (n_basis, n_basis)
      real*8 aux_matr1 (n_states, n_states)
      real*8 aux_matr2 (n_states, n_basis)
      character*18 filename
      character*15 filename1
      character*24 filename2
      complex*16 spectrum_SE
      real*8 g1  
      integer state_to_print_out 
      integer n_states_to_print
      logical analyze_peaks 

      !print*, 'here 1'
      first_state_included=1
      anacon_type1 = 1
      n_max_par1 = (nomega)/2
      aux_nomega = 6000
      aux_omegamax = 2.50
      allocate(aux_omega(aux_nomega))
      allocate(aux_womega(aux_nomega))
      allocate(aux_G (n_states,n_states,1))
      allocate(aux_G1(n_states,n_states,1))
      allocate(spectrum(aux_nomega))
      allocate(self_en_par_KS (n_max_par1,n_states,n_states))
   
     !calculate the homo level--------------------
      n_homo = 0
      do i_state = 1, n_states
        if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
          n_homo=i_state
        endif
      enddo
!      if(myid.eq.0) print*, "HOMO level :", n_homo
     !---------------------------------------------
 

     !select_state_to_print----------------
      analyze_peaks = .false.
      if(analyze_peaks)then
        n_states_to_print = 10
        if (n_states_to_print .gt. n_homo ) n_states_to_print = n_homo
        
        allocate(spectrum_single_state(n_states_to_print,aux_nomega))
        do k_state=1, n_states_to_print, 1

          if ( n_homo +1- k_state .lt.10)then
            write(string,'(A,I1)') "0", n_homo +1- k_state
          else
            write(string,'(I2)') n_homo +1- k_state
          endif

          filename2 = "spectrum_KS_state_"//string//".dat"
          if( myid.eq.0 ) open(34+k_state,file=filename2)
!          if( myid.eq.0 ) print *, filename2
        enddo
      else
        n_states_to_print = 0
      endif 
     !---------------------------------------------
      

      !parallelization of the auxiliary grid---------------- 
       n_remain = MOD(aux_nomega, n_tasks)
       if (n_remain.eq.0) then
         n_loc_grid1 = aux_nomega / n_tasks
       else
         n_loc_grid1 = aux_nomega / n_tasks + 1
       endif
       allocate(aux_map_index(n_tasks, n_loc_grid1))

       call distribute_grid (aux_nomega, aux_map_index,n_loc_grid1)
      !----------------------------------------------

      !transform sigma in the KS basis where the analytic continuation works better
      self_en_par_KS(:,:,:) = (0.d0 ,0.d0)
      self_energy_diag (:,:,:) = (0.d0 ,0.d0)

      call gauleg(0d0, aux_omegamax, aux_omega, aux_womega, aux_nomega)
      call diagonalize_self_en(self_energy_freq, self_energy_diag)

      call analy_continue_green_fn &
          (anacon_type1,&
           nomega, &
           n_max_par1, &
           self_en_par_KS,omega, &
           self_energy_diag, n_states)
      !---------------------------------------------------------------

      
!      if (.false.)then
!        if( number_reiterations.lt.10 ) then
!           write(iter,'(A,I1)') "0",number_reiterations
!        else
!           write(iter,'(I2)') number_reiterations
!        endif
!        filename = "hys_G_"//iter//".dat"
!        open(555,file = filename)
!        do i_basis = 1, n_states, 1
!          do j_basis = 1, n_states, 1
!            g1 = 0
!            do i_freq = 1, nomega, 1
!              g1 = g1 + abs(self_energy_diag(i_basis,j_basis,i_freq))*&
!                   womega1(i_freq)
!            enddo
!            write(555,'(I4, I4, f9.5)') i_basis, j_basis, g1
!          enddo
!          write(555,*) " "
!        enddo
!        close(555)
!      endif

       !convert everthing in the KS basis--------------------- 
       aux_matr  (:,:) = (0.d0,0.d0)
       aux_matr1 (:,:) = (0.d0,0.d0)
       aux_matr2 (:,:) = (0.d0,0.d0)
 
       if (.not. use_hartree_fock)then
         aux_matr (:,:) =                 &
         + xc_matr                (:,:)   &
         - exchange_self_energy   (:,:)   &
         - 2.d0 * hartree_pot     (:,:)
       else
!         aux_matr (:,:) =                  &
!         + exchange_self_energy0  (:,:)   &
!         - exchange_self_energy    (:,:)   &
!         - 2.d0 * hartree_pot      (:,:)

         aux_matr (:,:) =                                 &
         + xc_matr                                (:,:)   &
         + hybrid_coeff * exchange_self_energy0  (:,:)   &
         - exchange_self_energy                   (:,:)   &
         - 2.d0 * hartree_pot                     (:,:)
       endif 

!       call dgemm ('N','N',n_states,n_basis,n_basis,&
!       1.d0, ovlp_NAO_KS,n_states,aux_matr,n_basis,0.d0,&
!       aux_matr2,n_states)  
!       call dgemm ('N','T',n_states,n_states,n_basis,&
!       1.d0,aux_matr2,n_states,ovlp_NAO_KS,n_states,0.d0,&
!       aux_matr1,n_states)  

        do i_basis=1, n_basis, 1
          do j_basis=1, n_basis, 1
            do i_state=1, n_states, 1
              do j_state=1, n_states, 1
               aux_matr1 (i_state, j_state)=& 
               aux_matr1 (i_state, j_state)+& 
               KS_eigenvector (i_basis,i_state,i_spin,1)*&
               aux_matr (i_basis,j_basis)*& 
               KS_eigenvector (j_basis,j_state,i_spin,1)
              enddo
            enddo
          enddo
        enddo
       !------------------------------------------------------------

       !set file name of the output --------------------------------    
        if(scgw_converged) then
           if(n_spin .eq. 1 )then
           !  filename = "sp_ImG"//iter//".dat"
           filename1 = 'spectrum_sc.dat'
           if( myid.eq.0 )open(33,file=filename1)
           else
               if (i_spin ==1)  filename = "spectrum_sc_up.dat"
               if (i_spin ==2)  filename = "spectrum_sc_do.dat"
               if( myid.eq.0 )open(33,file=filename)
           endif
        else
           if (number_reiterations.lt.10)then
             write(iter,'(A,I1)') "0", number_reiterations
           else
             write(iter,'(I2)') number_reiterations
           endif

           if(n_spin .eq. 1 )then
             filename1 = "sp_ImG"//iter//".dat"
           else
               if (i_spin ==1)  filename1 = "sp_ImG"//iter//"_SU.dat"
               if (i_spin ==2)  filename1 = "sp_ImG"//iter//"_SD.dat"
           endif
           if( myid.eq.0 )open(33,file=filename1)
        endif
       !------------------------------------------------------------


      !if( myid.eq.0 )  print *, "start evaluation of the spectrum"
      do sign_freq = -1, 1, 2 !loop of sign + and - of the frequency axis
        spectrum(:) = 0.d0
        if(analyze_peaks) spectrum_single_state(:,:) = 0.d0

        do k = 1, n_loc_grid1, 1
          i_freq = aux_map_index (myid+1, k)
          if(i_freq.gt.0)then
            freq = sign_freq*aux_omega(i_freq)
         
            aux_G (:,:,:) = (0.d0,0.d0)
            aux_G1(:,:,:) = (0.d0,0.d0)
          
            do i_state = 1, n_states , 1
              do j_state = 1, n_states, 1
      
                call get_real_selfenergy (anacon_type1, nomega, &
                omega, freq, n_max_par1, self_en_par_KS(:,i_state,j_state), selfe)
!                if(myid.eq.0 .and. i_state.eq.j_state .and. i_state.eq.5)then
!                 write(34,*) freq*hartree, real(selfe),aimag(selfe)
!                endif
      
                aux_G (i_state,j_state,1) =&
                aux_G (i_state,j_state,1) &
                + aux_matr1(i_state,j_state)-selfe 
      
                if(i_state .eq. j_state)then
                  aux_G (i_state,i_state,1) =&
                  aux_G (i_state,i_state,1)+& 
                   (freq- (KS_eigenvalue(i_state,i_spin,1)-chemical_potential)+&
                   (0.d0,1.d0)*0.0001*&
                    sign(1.d0,KS_eigenvalue(i_state,i_spin,1)-chemical_potential))
                endif
              enddo
            enddo
      
            call invert_green(& 
            aux_G(1:n_states,1:n_states,1), 1,&
            aux_G1(1:n_states,1:n_states,1),n_states) 
      
            do i_state=1 , n_states,1
              do j_state=1, n_states,1     
                if(i_state .eq.j_state)then 
                 spectrum (i_freq) = spectrum(i_freq) +&
                 aimag(aux_G1(i_state,j_state,1))*2./pi
                endif 
              enddo
            enddo

            if(analyze_peaks) then      
              do k_state=1, n_states_to_print, 1 
                spectrum_single_state (k_state,i_freq) = &
                spectrum_single_state (k_state,i_freq) + &
                aimag(aux_G1(n_homo+1-k_state, n_homo+1-k_state,1))*2./pi 
!                if(i_freq.eq.1) print* , "writing state ", n_homo+1-k_state
              enddo
            endif

          endif 
        enddo ! loop over freq points

       !syncronization----------------------------------------------------       
        call sync_vector(spectrum, aux_nomega)
        if(analyze_peaks) then
          do k_state = 1, n_states_to_print, 1
            call sync_vector(spectrum_single_state(k_state,:),aux_nomega)
          enddo
        endif
       !------------------------------------------------------------------
 
       !write the spectrum to file----------
        if(myid.eq.0)then
          do i_freq = 1, aux_nomega, 1 
            write(33,*) (sign_freq*aux_omega(i_freq)+chemical_potential)*hartree, abs(spectrum(i_freq))

            if(analyze_peaks) then
              do k_state =1, n_states_to_print, 1 
                write(34+k_state,*) (sign_freq*aux_omega(i_freq)+chemical_potential)*hartree, &
                     abs(spectrum_single_state(k_state,i_freq))
              enddo
            endif
          enddo
        endif  
       !--------------------------
      enddo !end loop over sign 

      if(myid.eq.0)then
        close (33)
        if(analyze_peaks) then
          do k_state =1, n_states_to_print, 1 
            close (34+k_state)
          enddo
        endif
      endif

      deallocate (spectrum)
      if(analyze_peaks) deallocate (spectrum_single_state)
     !------------------------------------------


!-----------------------------------------------------------------------------------
!        if(plot_self_energy)then
!
!
!          if (number_reiterations.lt.10)then
!            write(iter,'(A,I1)') "0", number_reiterations
!          else
!            write(iter,'(I2)') number_reiterations
!          endif
!
!          if (selfe_matr_el.lt.10)then
!            write(matr_el,'(A,I1)') "0", selfe_matr_el
!          else
!            write(matr_el,'(I2)') selfe_matr_el 
!          endif
!  
!          filename1 = "SE_el"//matr_el//"_it"//iter//".dat"
!          open(45,file="spectrum_self_energy_"//iter//".dat")
!
!          open(34,file=filename1)
!          do sign_freq = -1, 1, 2
!            do i_freq = 1, aux_nomega, 1
!               freq = sign_freq*aux_omega(i_freq)
!               call get_real_selfenergy (anacon_type1, nomega, &
!               omega, freq, n_max_par1, self_en_par_KS(:,selfe_matr_el,selfe_matr_el), selfe)
!               write(34,*) (freq+chemical_potential)*hartree, real(selfe),aimag(selfe)
!            enddo
!          enddo
!
!          
!          do sign_freq = -1, 1, 2
!            do i_freq = 1, aux_nomega, 1
!              spectrum_SE = 0.d0
!              do i_state =1, n_states, 1
!               freq = sign_freq*aux_omega(i_freq)
!               call get_real_selfenergy (anacon_type1, nomega, &
!               omega, freq, n_max_par1, self_en_par_KS(:,i_state,i_state), selfe)
!               spectrum_SE = spectrum_SE + selfe
!              enddo
!            write(45,*) (freq+chemical_potential)*hartree, &
!               real(spectrum_SE),aimag(spectrum_SE)
!
!             enddo
!  
!          enddo
!
!          close (34) 
!          close (45) 
!        endif
!      endif
!
      end subroutine solve_dyson_equation_re_v0
!------------------------------------------------------------------------

      subroutine solve_dyson_equation_re (inv_green_fn_freq, &
          self_energy_freq, new_green_fn_freq,&
          exchange_self_energy, xc_matr, hartree_pot, &
!          inv_overlap_matrix, &
          exchange_self_energy0, &
          new_chem_pot, number_reiterations,&
          scgw_converged,i_spin, n_homo)


      use constants, only: pi, hartree
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use physics
      use scgw_grid
      use poles_fit


      implicit none
! INPUT
      !integer nomega
      complex*16  inv_green_fn_freq(n_basis,n_basis,nomega)
      complex*16  self_energy_freq (n_basis,n_basis,nomega)
      !real*8  omega(nomega)
      real*8  xc_matr(n_basis, n_basis)
      real*8  exchange_self_energy(n_basis,n_basis)
      real*8  hartree_pot(n_basis,n_basis)
      real*8  exchange_self_energy0 (n_basis,n_basis)
!      real*8  inv_overlap_matrix (n_basis, n_basis)
      real*8  new_chem_pot
!      real*8 ovlp_NAO_KS (n_states, n_basis)
      integer number_reiterations
      integer n_homo(n_spin) 
      logical scgw_converged
!OUTPUT
      complex*16  new_green_fn_freq(n_basis,n_basis,nomega)

!AUXILIARY 
      logical output
      integer i_freq
      integer i_spin
      integer k
      integer i_basis, j_basis
      integer i_state, j_state
      character*2 iter
      character*2 matr_el
      integer i_index       
      complex*16  self_energy_diag (n_states,n_states,nomega)

!the local grid
      integer aux_nomega
      real*8  aux_omegamax
      real*8 , dimension (:), allocatable :: aux_omega
      real*8 , dimension (:), allocatable :: aux_womega
      integer , dimension(:,:) , allocatable :: aux_map_index
      integer n_loc_grid1
      integer n_remain
!#      integer occ_numbers k

      complex*16, dimension (:,:,:),allocatable :: aux_G 
      complex*16, dimension (:,:,:),allocatable :: aux_G1 
      real*8, dimension(:), allocatable :: spectrum
      real*8, dimension(:,:), allocatable :: spectral_analysis
      real*8 freq
      integer sign_freq
      integer anacon_type1
      integer state_to_print_out 
      complex*16, dimension(:,:,:), allocatable:: self_en_par_KS
      integer n_max_par1
      complex*16 selfe
      real*8 aux_matr (n_basis, n_basis)
      real*8 aux_matr1 (n_states, n_states)
      real*8 aux_matr2 (n_states, n_basis)
      character*18 filename
      character*15 filename1
      complex*16 spectrum_SE
      real*8 g1  
      integer j_spin
!      real*8 occ_numbers(n_states, n_spin, 1) 


!      do j_spin = 1, n_spin
!       do i_state = 1, n_states
!         if(occ_numbers(i_state,j_spin,1).gt.1.d-6) then
!           n_homo(i_spin)=i_state
!         endif
!       enddo
!      enddo


      !print*, 'here 1'
      first_state_included=1
      anacon_type1 = 1
      n_max_par1 = (nomega)/2
!      n_max_par1 = n_max_par
      aux_nomega = 6000
      aux_omegamax = 2.50
      allocate(aux_omega(aux_nomega))
      allocate(aux_womega(aux_nomega))
      allocate(aux_G (n_states,n_states,1))
      allocate(aux_G1(n_states,n_states,1))
      allocate(spectrum(aux_nomega))
      allocate(self_en_par_KS (n_max_par1,n_states,n_states))

      !parallelization of the auxiliary grid 
    
       n_remain = MOD(aux_nomega, n_tasks)
       if (n_remain.eq.0) then
         n_loc_grid1 = aux_nomega / n_tasks
       else
         n_loc_grid1 = aux_nomega / n_tasks + 1
       endif
       allocate(aux_map_index(n_tasks, n_loc_grid1))

       call distribute_grid (aux_nomega, aux_map_index,n_loc_grid1)
      !----------------------------------------------
      !transform sigma in the KS basis where the analytic continuation works better

      self_en_par_KS(:,:,:) = (0.d0 ,0.d0)
      self_energy_diag (:,:,:) = (0.d0 ,0.d0)

      call gauleg(0d0, aux_omegamax, aux_omega, aux_womega, aux_nomega)
      call diagonalize_self_en(self_energy_freq, self_energy_diag)
      !print*, 'here 2'

!      if (.false.)then
!        if( number_reiterations.lt.10 ) then
!           write(iter,'(A,I1)') "0",number_reiterations
!        else
!           write(iter,'(I2)') number_reiterations
!        endif
!        filename = "hys_G_"//iter//".dat"
!        open(555,file = filename)
!        do i_basis = 1, n_states, 1
!          do j_basis = 1, n_states, 1
!            g1 = 0
!            do i_freq = 1, nomega, 1
!              g1 = g1 + abs(self_energy_diag(i_basis,j_basis,i_freq))*&
!                   womega1(i_freq)
!            enddo
!            write(555,'(I4, I4, f9.5)') i_basis, j_basis, g1
!          enddo
!          write(555,*) " "
!        enddo
!        close(555)
!      endif

       
      call analy_continue_green_fn &
          (anacon_type1,&
           nomega, &
           n_max_par1, &
           self_en_par_KS,omega, &
           self_energy_diag, n_states)
      
      !print*, 'here 3'
       !convert everthing in the KS basis 
       !(easier to do the fit because Sigma is diagonal)
       aux_matr  (:,:) = (0.d0,0.d0)
       aux_matr1 (:,:) = (0.d0,0.d0)
       aux_matr2 (:,:) = (0.d0,0.d0)
 
       if (.not. use_hartree_fock)then
         aux_matr (:,:) = -  exchange_self_energy(:,:)+  &
         xc_matr (:,:) - hartree_pot(:,:)*2
       else
         aux_matr (:,:) = -  exchange_self_energy(:,:)+  &
         exchange_self_energy0(:,:) - hartree_pot(:,:)*2
       endif 

!       call dgemm ('N','N',n_states,n_basis,n_basis,&
!       1.d0, ovlp_NAO_KS,n_states,aux_matr,n_basis,0.d0,&
!       aux_matr2,n_states)  
!       call dgemm ('N','T',n_states,n_states,n_basis,&
!       1.d0,aux_matr2,n_states,ovlp_NAO_KS,n_states,0.d0,&
!       aux_matr1,n_states)  

        do i_basis=1, n_basis, 1
          do j_basis=1, n_basis, 1
            do i_state=1, n_states, 1
              do j_state=1, n_states, 1
               aux_matr1 (i_state, j_state)=& 
               aux_matr1 (i_state, j_state)+& 
               KS_eigenvector (i_basis,i_state,i_spin,1)*&
               aux_matr (i_basis,j_basis)*& 
               KS_eigenvector (j_basis,j_state,i_spin,1)
              enddo
            enddo
          enddo
        enddo
    
      !print*, 'here 4'
!      if(myid.eq.0)then
!        if (number_reiterations.lt.10)then
!          write(iter,'(A,I1)') "0", number_reiterations
!        else
!          write(iter,'(I2)') number_reiterations
!        endif
!        filename = "sp_reG"//iter//".dat"
!        if(scgw_converged) then 
!          filename = "spectrum_sc.dat"
!        endif
!        open(33,file=filename)
!      endif

      if(scgw_converged) then
         if(n_spin .eq. 1 )then
         !  filename = "sp_ImG"//iter//".dat"
             filename1 = 'spectrum_sc.dat'
             if( myid.eq.0 )open(33,file=filename1)
         else
             if (i_spin ==1)  filename = "spectrum_sc_up.dat"
             if (i_spin ==2)  filename = "spectrum_sc_do.dat"
             if( myid.eq.0 )open(33,file=filename)
         endif
      else
         if (number_reiterations.lt.10)then
           write(iter,'(A,I1)') "0", number_reiterations
         else
           write(iter,'(I2)') number_reiterations
         endif

         if(n_spin .eq. 1 )then
           filename1 = "sp_ImG"//iter//".dat"
         else
             if (i_spin ==1)  filename1 = "sp_ImG"//iter//"_SU.dat"
             if (i_spin ==2)  filename1 = "sp_ImG"//iter//"_SD.dat"
         endif
         if( myid.eq.0 )open(33,file=filename1)
      endif

!-0-------------------------------------
! if required this prints out the diagonal matrix elements of the spectral function

!      if(.false.)then
!         if(myid.eq.0) print*,  "  Entering the spectral analysis" 
!         if(myid.eq.0) print*,  " HOMO = " , n_homo(1) 
!         if(myid.eq.0 .and. n_spin.gt.1) &
!               write(use_unit,*) "  -> Warning, this may explode with spin!"
!         allocate(spectral_analysis(n_states,aux_nomega))
!        
!         do sign_freq = -1, 1, 2
!         !  spectrum(:) = 0.d0
!           spectral_analysis (:,:) = 0.d0
!        
!           do k = 1, n_loc_grid1, 1
!           i_freq = aux_map_index (myid+1, k)
!             if(i_freq.gt.0)then
!!           do i_freq = 1, aux_nomega, 1
!                freq = sign_freq*aux_omega(i_freq)
!               
!                aux_G (:,:,:) = (0.d0,0.d0)
!                aux_G1(:,:,:) = (0.d0,0.d0)
!              
!                do i_state = 1, n_homo(1) , 1
!!                  do j_state = 1, n_states, 1
!               
!                    call get_real_selfenergy (anacon_type1, nomega, &
!                    omega, freq, n_max_par1, self_en_par_KS(:,i_state,i_state), selfe)
!               
!                    aux_G (i_state,i_state,1) =&
!                    aux_G (i_state,i_state,1) &
!                    + aux_matr1(i_state,i_state)-selfe 
!               
!                     aux_G (i_state,i_state,1) =&
!                     aux_G (i_state,i_state,1)+& 
!                      (freq- (KS_eigenvalue(i_state,1,1)-chemical_potential)+&
!                      (0.d0,1.d0)*0.0001*&
!                       sign(1.d0,KS_eigenvalue(i_state,1,1)-chemical_potential))
!                enddo
!               
!                call invert_green(& 
!                aux_G(1:n_homo(1),1:n_homo(1),1), 1,&
!                aux_G1(1:n_homo(1),1:n_homo(1),1),n_homo(1)) 
!               
!                do i_state=1 , n_homo(1),1
!                      spectral_analysis (i_state,i_freq) = &
!                       spectral_analysis (i_state,i_freq) +&
!                       aimag(aux_G1(i_state,i_state,1))*2./pi
!                enddo
!           
!             endif 
!           enddo
! 
!           do i_state = 1, n_homo(1), 1
!             call sync_vector(spectral_analysis (i_state,1),aux_nomega)
!           enddo
!
!!           state_to_print_out = n_homo(1)
!
!!           if (n_homo(1) .lt. 10 ) then
!!             state_to_print_out  = n_homo (1)-1 
!!           endif
!
!           !print out everything
!           !do i_spin = 1, n_spin
!!             do i_state = n_homo(i_spin),  n_homo(i_spin)-state_to_print_out,-1
!             do i_state = 1,  n_homo(1),1
!                if(scgw_converged) then
!                   if (i_state.lt.10)then
!                     write(iter,'(A,I1)') "0", i_state
!                   else
!                     write(iter,'(I2)') i_state
!                   endif
!                   filename1 = "sp_ana"//iter//".dat"
!                   if( myid.eq.0 ) then
!                      open(36,file=filename1)
!                      do i_freq = 1, aux_nomega, 1
!                        write(36,*) (sign_freq*aux_omega(i_freq)+&
!                       chemical_potential)*hartree, &
!                        abs(spectral_analysis(i_state,i_freq))
!                      enddo
!                      close(36)
!                   endif
!                endif
!             enddo
!           !enddo
!         enddo
!         deallocate(spectral_analysis)
!      endif !if (.true.)

!-----------------------------------------------------------
! Here I calculate the intgrated spectral function

      do sign_freq = -1, 1, 2
        spectrum(:) = 0.d0

        do k = 1, n_loc_grid1, 1
        i_freq = aux_map_index (myid+1, k)
        if(i_freq.gt.0)then
!        do i_freq = 1, aux_nomega, 1
          freq = sign_freq*aux_omega(i_freq)
       
          aux_G (:,:,:) = (0.d0,0.d0)
          aux_G1(:,:,:) = (0.d0,0.d0)
        
          do i_state = 1, n_states , 1
            do j_state = 1, n_states, 1

              call get_real_selfenergy (anacon_type1, nomega, &
              omega, freq, n_max_par1, self_en_par_KS(:,i_state,j_state), selfe)
!              if(myid.eq.0 .and. i_state.eq.j_state .and. i_state.eq.5)then
!               write(34,*) freq*hartree, real(selfe),aimag(selfe)
!              endif

              aux_G (i_state,j_state,1) =&
              aux_G (i_state,j_state,1) &
              + aux_matr1(i_state,j_state)-selfe 
 
              if(i_state .eq. j_state)then
               aux_G (i_state,i_state,1) =&
               aux_G (i_state,i_state,1)+& 
                (freq- (KS_eigenvalue(i_state,i_spin,1)-chemical_potential)+&
                (0.d0,1.d0)*0.0001*&
                 sign(1.d0,KS_eigenvalue(i_state,i_spin,1)-chemical_potential))
              endif
            enddo
          enddo
       !     if(myid.eq.0 .and. i_freq .eq.1) print *, aux_G

       call invert_green(& 
       aux_G(1:n_states,1:n_states,1), 1,&
       aux_G1(1:n_states,1:n_states,1),n_states) 
          spectrum = 0.d0
          do i_state=1 , n_states,1
            do j_state=1, n_states,1     
              if(i_state .eq.j_state)then
                spectrum (i_freq) = spectrum(i_freq) +&
                 aimag(aux_G1(i_state,j_state,1))*2./pi
              endif
            enddo
          enddo
        
          endif 
        enddo

        call sync_vector(spectrum,aux_nomega)
!  write the spectrum to file----------
        if(myid.eq.0)then
          do i_freq = 1, aux_nomega, 1 
            write(33,*) (sign_freq*aux_omega(i_freq)+chemical_potential)*hartree, abs(spectrum(i_freq))
          enddo
        endif  
        !--------------------------
      enddo
      


      if(myid.eq.0)then

        close (33)
      endif
!        if(plot_self_energy)then
!
!
!          if (number_reiterations.lt.10)then
!            write(iter,'(A,I1)') "0", number_reiterations
!          else
!            write(iter,'(I2)') number_reiterations
!          endif
!
!          if (selfe_matr_el.lt.10)then
!            write(matr_el,'(A,I1)') "0", selfe_matr_el
!          else
!            write(matr_el,'(I2)') selfe_matr_el 
!          endif
!  
!          filename1 = "SE_el"//matr_el//"_it"//iter//".dat"
!          open(45,file="spectrum_self_energy_"//iter//".dat")
!
!          open(34,file=filename1)
!          do sign_freq = -1, 1, 2
!            do i_freq = 1, aux_nomega, 1
!               freq = sign_freq*aux_omega(i_freq)
!               call get_real_selfenergy (anacon_type1, nomega, &
!               omega, freq, n_max_par1, self_en_par_KS(:,selfe_matr_el,selfe_matr_el), selfe)
!               write(34,*) (freq+chemical_potential)*hartree, real(selfe),aimag(selfe)
!            enddo
!          enddo
!
!          
!          do sign_freq = -1, 1, 2
!            do i_freq = 1, aux_nomega, 1
!              spectrum_SE = 0.d0
!              do i_state =1, n_states, 1
!               freq = sign_freq*aux_omega(i_freq)
!               call get_real_selfenergy (anacon_type1, nomega, &
!               omega, freq, n_max_par1, self_en_par_KS(:,i_state,i_state), selfe)
!               spectrum_SE = spectrum_SE + selfe
!              enddo
!            write(45,*) (freq+chemical_potential)*hartree, &
!               real(spectrum_SE),aimag(spectrum_SE)
!
!             enddo
!  
!          enddo
!
!          close (34) 
!          close (45) 
!        endif
!      endif
!
      end subroutine solve_dyson_equation_re

