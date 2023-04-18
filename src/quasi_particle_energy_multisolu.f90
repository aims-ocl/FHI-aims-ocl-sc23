!****s* FHI-aims/quasi_particle_energy_multisolu
!  NAME
!   quasi_particle_energy_multisolu
!  SYNOPSIS

      subroutine quasi_particle_energy_multisolu &
           (anacon_type,n_max_par, &
            n_low_state,n_high_state, &
            n_homo, n_lumo, n_freq, omega, &
            sigma_par,occ_numbers,KS_eigenvalue, &
            chemical_potential_spin, &
            exchange_self_energy, &
            xc_KS_matr,x_KS_array,c_KS_array, &
            qp_energy )

!  PURPOSE
!  Subroutine quasi_partilce_energy calculation the quasi particle energies
!  for the GW calculation. Note, now the self-energy has been calculated and
!  continuted to the real frequency axis
!
!  USES

      use dimensions
      use runtime_choices
      use mpi_tasks
      use constants
      use gw_para,only :out_self_energy
      use localorb_io,only :use_unit
      implicit none

!  ARGUMENTS

      integer :: anacon_type 
      integer :: n_max_par
      integer :: n_low_state, n_high_state
      integer :: n_freq
      integer :: n_homo(n_spin)
      integer :: n_lumo(n_spin)

      real*8 :: omega(n_freq)
      real*8 :: occ_numbers(n_states,n_spin)
      real*8 :: KS_eigenvalue(n_states,n_spin)
      real*8 :: exchange_self_energy(n_high_state,n_spin)
      real*8 :: xc_KS_matr(n_states,n_states,n_spin)
      real*8 :: x_KS_array(n_states,n_spin)
      real*8 :: c_KS_array(n_states,n_spin)
      real*8 :: chemical_potential_spin(n_spin)

      complex*16 :: sigma_par(n_max_par,n_states,n_spin)

      real*8 :: qp_energy(n_low_state:n_high_state,n_spin)

!  INPUTS
!  o  anacon_type -- integer number, if 0, the two-pole fitting for analytical
!          continuation; if 1, using Pade approximation for ana. cont.        
!  o  n_max_par -- the number of parameters used for analytical continuation 
!          For anacon_type = 0, recommended n_max_par is  4 or 5. If 4, this will 
!          be the normal two-pole fitting, else if 5, it will be two-pole plus a 
!          (small) constant number
!          For anacon_type = 1, recommended n_max_par is the half of n_freq  
!  o  n_low_state  -- integer number,
!          the lowest KS/HF eigenstate for self-energy correction
!  o  n_high_state -- integer number,
!          the highest KS/HF eigenstate for self-energy correction
!  o  n_freq -- integer number, the number of frequency points for the GW self-energy
!  o  n_homo -- integer array, the HOMO level for each spin channel
!  o  n_lumo -- integer array, the LUMO level for each spin channel
!  o  omega(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the self-energy
!  o  occ_numbers -- occupation numbers of single-particle energy levels
!  o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
!  o  chemical_potential -- real number, the chemical potential of the system
!  o  exchange_self_energy -- real array the exact exchange part of the self-energy 
!            of each KS/HF state for each spin channel
!  o  xc_KS_matr -- real array, the matrix elements of the exchange correlation 
!            potential witin KS/HF orbitals
!  o  x_KS_array -- real array, the matrix elements of the exchange
!            potential witin KS/HF orbitals
!  o  c_KS_array -- real array, the matrix elements of the correlation
!            potential witin KS/HF orbitals
!  o  sigma_par -- complex array, the fitting parameters from analytical continuation
! OUTPUTS
!  o  qp_energy -- real array, the final calculated quasiparticle energy for each
!            concerned state (from n_low_state to n_high_state).

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

      real*8, allocatable :: xc_energy(:,:)
      real*8, allocatable :: qp_energy_old(:,:)
      real*8, allocatable :: correl_energy(:,:)
      real*8, allocatable :: qp_energy_multi(:,:,:)
      real(kind=8), dimension(:,:), allocatable :: my_exchange_self_energy

      integer, allocatable :: n_qpe_multi(:,:)

      real*8  e_diff
      real*8  en, en1
      real*8  selfe_real
      real*8  mu
      real*8  delta_mu

!      real*8 :: qpe_sum

      complex*16 selfe, dselfe

      logical :: converged, new_solution_found

      integer, parameter ::  n_steps = 40
      integer, parameter ::  n_solu = 10
      integer, parameter ::  n_real_freq = 2000
      real*8, parameter ::  qp_energy_thr = 1.d-4
      character*80  filename

!     counters


      integer :: i_state,i_state_1
      integer :: i, i_count
      integer :: i_spin
      integer :: i_try
      integer :: i_solu
      integer :: i_freq

!   external function


!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)') "Quasi particle energy calculation starts..."
      endif

      if(myid.eq.0) then
         write(use_unit,'(2X,3A)') "    state  ", "   spin   ",  "solution_no"
      endif

!      write(use_unit,*) "----------------------------------------------"
!      write(use_unit,*) " Exchange term ..."
!      write(use_unit,*)
!      do i_spin = 1, n_spin
!        do i_state = 1, n_high_state
!          write(use_unit,'(I4,f16.4)') i_state, &
!                 exchange_self_energy(i_state,i_spin)*hartree
!        enddo
!        write(use_unit,*)
!      enddo

      allocate(my_exchange_self_energy(n_high_state,n_spin)) 
      if(use_hartree_fock .and. (.not.use_screx) &
            .and. .not. use_gw_and_hse) then
        my_exchange_self_energy(:,:) = exchange_self_energy * & 
                                 (1.0d0-hybrid_coeff)
      else
        my_exchange_self_energy(:,:) = exchange_self_energy
      endif

!  read in DFT exchange-correlation energy for all states
      allocate(xc_energy(n_states,n_spin))
!      open (106, file='xc_energy.dat')
      if(use_split_xc_gw) then
        do i_spin = 1, n_spin, 1
         do i_state = 1, n_states, 1
            xc_energy(i_state,i_spin) = &
            x_KS_array(i_state,i_spin)+c_KS_array(i_state,i_spin)
         enddo
        enddo
      else
        do i_spin = 1, n_spin, 1
         do i_state = 1, n_states, 1
!         read(106, *) i, xc_energy(i_state,i_spin)
            xc_energy(i_state,i_spin) = &
            xc_KS_matr(i_state,i_state,i_spin)

!         if (i.ne.i_state) then
!          write(use_unit,'(2X,A)') "Wrong in reading the xc energy data file"
!          stop
!         endif
          enddo
        enddo
      endif
!      close(106)

!  quasi particle energy calculation

      allocate(qp_energy_old(n_states,n_spin))
      allocate(correl_energy(n_states,n_spin))

      allocate(qp_energy_multi(n_low_state:n_high_state,n_spin,n_solu))
      allocate(n_qpe_multi(n_low_state:n_high_state,n_spin))
!   determine homo level

      qp_energy(n_low_state:n_high_state,1:n_spin)= &
              KS_eigenvalue(n_low_state:n_high_state,1:n_spin)
!      qp_energy_old(:,:)=KS_eigenvalue(:,:)
      n_qpe_multi(:,:) = 0
      qp_energy_multi(:,:,:) = 0.d0
      do i_try = 1, n_steps
        qp_energy_old(:,:)=KS_eigenvalue(:,:)+ (i_try/dble(n_steps/2.d0) - 0.5d0)
      do i_spin = 1, n_spin
        do i_state = n_low_state, n_high_state

          i_state_1 = i_state - n_low_state + 1

          e_diff = 1.d-3
          i_count =0
          converged = .true.
          do while (abs(e_diff).gt.qp_energy_thr)
            i_count = i_count +1
            qp_energy(i_state,i_spin) = &
                  qp_energy_old(i_state,i_spin) + 0.5d0* e_diff
            qp_energy_old(i_state,i_spin) = qp_energy(i_state,i_spin)

            mu =  chemical_potential_spin(i_spin)

            en = qp_energy(i_state,i_spin)-mu

            call get_real_selfenergy(anacon_type,n_freq,omega, &
                      dcmplx(en,0.d0), n_max_par, &
                      sigma_par(1:n_max_par,i_state_1,i_spin), selfe)

!          write(use_unit,'(I4,2f13.6)') i_state,
!     +          hartree*qp_energy(i_state,i_spin),
!     +          real(selfe)*hartree

           qp_energy(i_state,i_spin)= KS_eigenvalue(i_state,i_spin) &
                      + real(selfe) &
                      + my_exchange_self_energy(i_state,i_spin) &
                      - xc_energy(i_state,i_spin)

           e_diff =  qp_energy(i_state,i_spin) &
                     - qp_energy_old(i_state,i_spin)

          if(i_count .gt. 500) then
                if(myid == 0) then
                   write(use_unit,'(2X,2A,I4,A,I4 )') &
                  "QUASI_PARTILCE_ENERGY: self-consistent quasiparticle ", &
                  "solution can not  be found for state: i_state = ", &
                  i_state, "  i_spin = ", i_spin 
                endif
              
               converged = .false.

            exit
          endif

! end of do while
         enddo
!         if(i_state .eq. 13) then
!          write(use_unit,'(I4,2f18.8)')i_try, qp_energy(i_state,i_spin), qp_energy(i_state,i_spin)*hartree
!         endif

         if(converged) then

           new_solution_found = .true.
       
           do i_solu =1, n_qpe_multi(i_state,i_spin)
              if (abs(qp_energy_multi(i_state,i_spin,i_solu) - qp_energy(i_state,i_spin) ) &
                      .le. qp_energy_thr ) then
                    new_solution_found = .false.
                  exit
              endif
           enddo

         
           if(new_solution_found) then
               n_qpe_multi(i_state,i_spin) = n_qpe_multi(i_state,i_spin) + 1
               i_solu = n_qpe_multi(i_state,i_spin)
               qp_energy_multi(i_state,i_spin,i_solu) = qp_energy(i_state,i_spin)
               if(myid.eq.0) then
                  write(use_unit,'(2X,3I8,2f18.8)') i_state, i_spin, i_solu, &
                      qp_energy_multi(i_state,i_spin,i_solu), &
                      qp_energy_multi(i_state,i_spin,i_solu)*hartree
               endif
           endif
!  end of if(converged)
          endif

         correl_energy(i_state,i_spin) = real(selfe)
! end of do i_state
         enddo
! end of do i_spin
        enddo
! end of do i_try
        enddo

        if (myid.eq.0 .and. out_self_energy) then
          do i_spin = 1, n_spin
            do i_state =  n_low_state, n_high_state, 1

             i_state_1 = i_state - n_low_state + 1
             if(n_spin.eq.2) then
               write(filename,'(A,I0,A,I0,A)') &
                          "self_energy/Sigma.real_freq.n_", &
                           i_state,".s_",i_spin,".dat"
             else
                write(filename,'(A,I0,A)') &
                          "self_energy/Sigma.real_freq.n_", &
                           i_state,".dat"
             endif
 
             open(107, file=filename)

             do i_freq = 1, n_real_freq, 1


               en=qp_energy(i_state,i_spin)+dble(i_freq)*2.d0/dble(n_real_freq) - 1.d0

               en1 = en-chemical_potential_spin(i_spin)

               call get_real_selfenergy(anacon_type,n_freq,omega, &
                      dcmplx(en1,0.d0), n_max_par, &
                      sigma_par(1:n_max_par,i_state,i_spin), selfe)

               selfe_real = KS_eigenvalue(i_state,i_spin) &
                      + real(selfe) &
                      + my_exchange_self_energy(i_state,i_spin) &
                      - xc_energy(i_state,i_spin)

               write (107, '(2X, f16.8, 7X, F18.10)')en*hartree, selfe_real*hartree
             enddo

             close(107)
          enddo
         enddo
!  end of if (myid.eq.0 .and. out_self_energy)
       endif

      if(myid.eq.0) then
!        open (107, file= 'qp_energy.dat')
!        write(107, '(2X, A,6X,A,10X,A,10X,A,10X,A, 10X, A)') &
!          "state", "E_GS", "E_X", "E_XC", "E_C", "E_QP"
!        do i_spin = 1, n_spin
!         write(107,*)"-----------------", &
!         "-------------------------------------------------------------"
!        do i_state = n_low_state, n_high_state, 1
!
!             write (107, '(2X, I4, 2X,6F14.4)') &
!                    i_state, KS_eigenvalue(i_state,i_spin)* hartree, &
!                    exchange_self_energy(i_state,i_spin) * hartree, &
!                    xc_energy(i_state,i_spin) * hartree, &
!                    correl_energy(i_state,i_spin)*hartree, &
!                    qp_energy(i_state,i_spin) * hartree
!         enddo
!        enddo

        write(use_unit,*)
        write(use_unit,'(2A)')"-----------------------------------------", &
         "-------------------------------------------------------------"
        write(use_unit,'(15X,A)')"GW quasi-particle energy levels"
        write(use_unit,*)
        if(use_split_xc_gw) then

          write(use_unit,'(15X,A)') &
              "e_qp = e_gs + e_x^ex - e_x^gs - e_c^gs + e_c^nloc"
          write(use_unit,*)
         write(use_unit, '(2X, A, 5X, A,8X, A, 8X,A,8X,A,8X,A,8X,A,8X, A)') &
             "state", "occ_num", "e_gs", "e_x^ex", "e_x^gs", "e_c^gs", &
              "e_c^nloc", "e_qp"
        else
          write(use_unit,'(15X,A)')"e_qp = e_gs + e_x^ex - e_xc^gs + e_c^nloc"
          write(use_unit,*)
          write(use_unit, '(2X, A, 5X, A,8X, A, 8X,A,8X,A,8X,A, 8X, A)') &
             "state", "occ_num", "e_gs", "e_x^ex", "e_xc^gs", &
              "e_c^nloc", "e_qp"
        endif
        write(use_unit,'(2A)')"-----------------------------------------", &
         "-------------------------------------------------------------"
        do i_spin = 1, n_spin
          if(n_spin.eq.2.and.i_spin.eq.1) then
            write(use_unit,'(35X, A)') "Spin Up"
          endif

          if(n_spin.eq.2.and.i_spin.eq.2) then
            write(use_unit,'(35X, A)') "Spin Down"
          endif

        write(use_unit,'(2A)')"-----------------------------------------", &
         "-------------------------------------------------------------"
          do i_state = n_low_state, n_high_state, 1

           if(use_split_xc_gw) then
             write(use_unit, '(2X, I6, 2X,F8.4, 6F14.4)') &
                   i_state, occ_numbers(i_state,i_spin), &
                   KS_eigenvalue(i_state,i_spin)*hartree, &
                   my_exchange_self_energy(i_state,i_spin)*hartree, &
                   x_KS_array(i_state,i_spin)*hartree, &
                   c_KS_array(i_state,i_spin)*hartree, &
                   correl_energy(i_state,i_spin)*hartree, &
                   qp_energy(i_state,i_spin)*hartree
           else
             write(use_unit, '(2X, I6, 2X, F8.4, 5F14.4)') &
                   i_state, occ_numbers(i_state,i_spin), &
                   KS_eigenvalue(i_state,i_spin)*hartree, &
                   my_exchange_self_energy(i_state,i_spin)*hartree, &
                   xc_energy(i_state,i_spin)*hartree, &
                   correl_energy(i_state,i_spin)*hartree, &
                   qp_energy(i_state,i_spin)*hartree
!  end of if use_split_xc_gw
          endif
        enddo


        if(n_spin.eq.2.and.i_spin.eq.1) then
         write(use_unit,'(2A)')"------------------------------------------", &
          "------------------------------------------------------------"
        endif
! end of i_spin
      enddo
      write(use_unit,'(2A)')"------------------------------------------", &
        "------------------------------------------------------------"
      write(use_unit,*)

      if(n_homo(1).gt.0 .and. n_homo(n_spin).gt.0)then ! for H atom
          write(use_unit,'(2X, A, 2f10.4)') &
          " DFT/Hartree-Fock HOMO level (eV): ", &
         KS_eigenvalue(n_homo(1),1)*hartree, &
         KS_eigenvalue(n_homo(n_spin),n_spin)*hartree
         
        write(use_unit,*)
        write(use_unit,'(2X, A, 2f10.4)') &
        " Quasiparticle HOMO level (eV):    ", &
          qp_energy(n_homo(1),1)*hartree, &
          qp_energy(n_homo(n_spin),n_spin)*hartree
      endif
      write(use_unit,*)
      write(use_unit,'(2A)')"------------------------------------------", &
        "------------------------------------------------------------"
      write(use_unit,*)

! end if myid=0
      endif

      deallocate (my_exchange_self_energy)
      deallocate (xc_energy)
      deallocate (qp_energy_old)
      deallocate (correl_energy)
      deallocate (qp_energy_multi)
      deallocate (n_qpe_multi)

      return
      end subroutine quasi_particle_energy_multisolu


!---------------------------------------------------------------------
!******
