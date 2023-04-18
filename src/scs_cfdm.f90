      module scs_cfdm
      use localorb_io
      use constants
      use dimensions
      use geometry
      use vdw_correction
      use species_data
      use runtime_choices
      use synchronize_mpi
      use mpi_utilities
      use relaxation
      implicit none
      PRIVATE ! DEFAULT PRIVATE
      PUBLIC :: self_consistent_screening_CFDM !<--calculate MBD energy as described in Phys. Rev. Lett. 108, 236402 (2012)
      PUBLIC :: get_mbd_enengy_at_TS           !<--  MBD using TS input paramenter
      PUBLIC :: self_consistent_screening_iso_damp !<-- MBD using isotropic damped tensor in SCS
      !! pipe to force module    
      PUBLIC :: gauss_legendre_grid , SCS_TENSOR_MBD_rsSCS ,MBD_TENSOR_MBD_rsSCS, get_alpha_omega_and_Rp

!     <<ENERGY OUTPUTS_FROM_THIS_MODULE>> 
!      real*8 :: ene_pairwise_ts      !<---C6/R^6 using TS alpha and C6
!      real*8 :: ene_pairwise_scs     !<---C6/R^6-using SCS alpha and C6
!      real*8 :: ene_mbd_scs          !<---MBD using CFDM hamiltonian with plasmon frequencies obtained from SCS
!      real*8 :: ene_mbd_ts           !<---MBD using TS plasmon frequencies Ambrosetti et al[1.A] .
!      real*8 :: ene_mbd_rsSCS        !<---MBD using plasmon frequencies obtained from range seperated SCS  Ambrosetti et al

      character*200 :: info_str  
      real*8,dimension(:,:),allocatable:: coupled_atom_pol
      real*8,dimension(:,:),allocatable:: relay_matrix
      real*8,dimension(:,:),allocatable:: relay_matrix_periodic
      real*8,dimension(:),allocatable:: alpha_omega
      real*8,dimension(:),allocatable:: R_p
      real*8,dimension(:),allocatable:: Rvdw_iso
      real*8,dimension(:),allocatable:: alpha_eff
      real*8,dimension(:),allocatable:: C6_eff
      real*8,dimension(:,:,:),allocatable::fd_conf_cords
      real*8,dimension(:),allocatable::fd_conf_ene_mbd
      real*8 :: fd_step_size
      real*8,dimension(3,3):: mol_pol_tensor
      real*8 :: casmir_omega(0:20)
      real*8 :: casmir_omega_weight(0:20)
      contains

      subroutine self_consistent_screening_CFDM(ene_pairwise_ts,ene_pairwise_scs,ene_mbd_scs)
      ! local variable 
      real*8 :: C6_free
      real*8 :: alpha_free
      real*8 :: R_free
      real*8 :: mol_c6_coeff
      real*8 :: C6_atom
      real*8 :: R_vdw_free
      integer :: i_freq
      integer :: i_myatom
      real*8 :: mol_pol_eigen(3)
      real*8 :: WORK(9)
      integer::errorflag,LWORK
      real*8 :: ene_pairwise_ts,ene_pairwise_scs,ene_mbd_scs 
      ! o INPUT
      ! o hirshfeld_volume from hirshfeld partioning 
      ! o relay_matrix is 3N x 3N (imaginary) frequency dependent matrix
      ! o R_p is effective dipole spread of quantum harmonic oscilator defined
      !   with respect to Mayer representation A. Mayer, Phys. Rev. B, 75,045407(2007)
      ! o alpha_omega is single pole polarizability of atom in molecule
      ! o coupled_atom_pol contains atom resolved screened polarizabilities at every
      !   imaginary frequency 
      ! o alpha_eff is atom resolved screened polarizabilty at ZERO frequency 
      ! o C6_eff is atom resolved screened C6 coefficients calculated by
      !   integrating coupled_atom_pol over entire frequency range 



      if(.not.allocated(relay_matrix))     allocate(relay_matrix(3*n_atoms,3*n_atoms))
      if(n_periodic .gt. 0) then
      if(.not.allocated(relay_matrix_periodic))     allocate(relay_matrix_periodic(3*n_atoms,3*n_atoms))
      endif 
      if(.not.allocated(R_p))              allocate(R_p(n_atoms))
      if(.not.allocated(alpha_omega))      allocate(alpha_omega(n_atoms))
      if(.not.allocated(coupled_atom_pol)) allocate(coupled_atom_pol(20,n_atoms))
      if(.not.allocated(alpha_eff))        allocate(alpha_eff(n_atoms))
      if(.not.allocated(C6_eff))           allocate(C6_eff(n_atoms))


      casmir_omega        = 0.d0
      casmir_omega_weight = 0.d0 
!     call gauleg(0.d0,30.d0,casmir_omega(1:20),casmir_omega_weight(1:20),20)
      ! currently I am using   
      call gauss_legendre_grid(casmir_omega,casmir_omega_weight)

      coupled_atom_pol = 0.d0
          mol_c6_coeff = 0.d0
             alpha_eff = 0.d0
                C6_eff = 0.d0 
 

     write(info_str,'(2x,A)')"| Many-Body Dispersion (MBD@SCS) energy "
     call localorb_info(info_str, use_unit,'(A)',OL_norm)
!     write(info_str,'(2x,A)')"| described in Phys. Rev. Lett. 108, 236402 (2012)"
!     call localorb_info(info_str, use_unit,'(A)',OL_norm)


      if(n_periodic .eq. 0) then
        write(info_str,'(2x,A)')"| Dynamic molecular polarizability alpha(iw) (bohr^3)"
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
      else
        write(info_str,'(2x,A)')"| Dynamic polarizability of unit cell alpha(iw) (bohr^3)"
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
      endif
      write(info_str,'(2x,A)')"| ----------------------------------------------------------------------------------"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,A)')"|  omega(Ha)   alpha_iso      alpha_xx       alpha_yy       alpha_zz"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)


    ! Loop over Casimir-Polder frequencies
      do i_freq=0,20,1
              !! zeroout array before getting actual frequency dependent parameters
              R_p = 0.0d0
              alpha_omega = 0.0d0
              relay_matrix= 0.0d0
             mol_pol_tensor= 0.0d0
              ! loop over atoms
              do i_myatom=1,n_atoms,1
                   if (.not.vdw_hirshfeld_data_external(species(i_myatom))) then                  
                   call get_vdw_param(species_element(species(i_myatom)),&
                                 species_z(species(i_myatom)),C6_free,alpha_free,R_vdw_free)
                   call get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),C6_free,&
                              alpha_free,casmir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                   else !! Case when parameter(s) C_6, alpha and R_vdw are specified in control.in 
                   call get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),vdw_hirshfeld_C6(species(i_myatom)),&
                              vdw_hirshfeld_alpha(species(i_myatom)),casmir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                   endif 

              enddo ! end loop over atoms


              ! get screened polarizability matrix(INPUT:-alpha_omega,R_p)
                call get_alpha_omega_relay_matrix()

              ! Contract relay_matrix^-1 for polarizability
                call contract_matrix(mol_pol_tensor)  

              ! get coupled atomic "local" polarizability by summing over row or columns of
              ! relay tensor
              if (i_freq.eq.0) then
                  call get_coupled_atom_polarizability(alpha_eff)
              endif

              if (i_freq.ge.1) then
                  call get_coupled_atom_polarizability(coupled_atom_pol(i_freq,:))
              endif 

              ! Casmir-Polder integral for molecular C6-coefficent 
              LWORK=9
              call DSYEV('V','U',3,mol_pol_tensor,3,mol_pol_eigen,WORK,LWORK,errorflag) 
              call check_info(errorflag,"self_consistent_screening_CFDM","DSYEV")
              if(i_freq.ge.1)  then
                  mol_c6_coeff = mol_c6_coeff + (casmir_omega_weight(i_freq)* (sum(mol_pol_eigen)/3.0)**2)
              endif

              write(info_str,'(2x,"| ",F10.6,4(e15.6))')casmir_omega(i_freq),&
                     sum(mol_pol_eigen)/3.0,mol_pol_eigen(1),mol_pol_eigen(2),mol_pol_eigen(3)
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
      enddo ! end   over freq


              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  

              if(n_periodic .eq. 0) then   
              Write(info_str,'(2x,A,f20.8,2x,A)')&
              "| Molecular C6 coefficient            :  ",0.954929658551372d0*mol_c6_coeff,"    hartree bohr^6"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              else
              Write(info_str,'(2x,A,f20.8,2x,A)')&
              "| Unit cell C6 coefficient            :  ",0.954929658551372d0*mol_c6_coeff,"    hartree bohr^6"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              endif

              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              if(n_periodic .eq. 0) then
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha (bohr^3)  of atom in molecule"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              else
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha (bohr^3)  of atom in unit cell"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
           
              endif  
              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  

              do i_myatom=1,n_atoms,1
                 C6_atom = 0.0d0 
                 do i_freq=1,20,1
                 C6_atom = C6_atom + (casmir_omega_weight(i_freq)*coupled_atom_pol(i_freq,i_myatom)**2)
                 enddo
                 write(info_str,'(2x,A,I4,2x,A,f12.6,6x,f12.6)')"| ATOM ",i_myatom,&
                                      trim(species_name(species(i_myatom))),&
                                      C6_atom*0.954929658551372d0,alpha_eff(i_myatom)
                 call localorb_info(info_str, use_unit,'(A)',OL_norm)  

                 !store effctive C6 for post SCS routines 
                  C6_eff(i_myatom)=C6_atom*0.954929658551372d0
              enddo 


              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)


              ! Compute TS+van der Waals energy 
               call get_ts_vdw(ene_pairwise_ts) 

              ! Compute TS+SCS van der Waals energy
               call get_ts_scs_vdw(ene_pairwise_scs)

              ! Compute many-body-dispersion energy 
               call get_many_body_dispersion_energy_cfdm(ene_mbd_scs)

   

              if(allocated(R_p))                   deallocate(R_p)
              if(allocated(alpha_omega))           deallocate(alpha_omega)
              if(allocated(relay_matrix))          deallocate(relay_matrix)
              if(allocated(relay_matrix_periodic)) deallocate(relay_matrix_periodic)
              if(allocated(coupled_atom_pol))      deallocate(coupled_atom_pol)
              if(allocated(alpha_eff))             deallocate(alpha_eff)
              if(allocated(C6_eff))                deallocate(C6_eff)



      return
      endsubroutine self_consistent_screening_CFDM

subroutine contract_matrix(tensor)
        implicit none

        integer::ir,ic,i,j,i_row,i_col
        real*8,dimension(3,3)::tensor
        tensor(:,:)=0.0d0

        do ir=1,n_atoms,1
             do ic=1,n_atoms,1
                 i_row=0
                    do i=3*ir-2,3*ir,1
                       i_row=i_row+1
                       i_col=0
                          do j=3*ic-2,3*ic,1
                          i_col=i_col+1
                 tensor(i_row,i_col)=tensor(i_row,i_col) + relay_matrix(i,j)
                          enddo
                    enddo
             enddo
        enddo

        return
endsubroutine contract_matrix

subroutine contract_matrix_v1(n_atoms,temp_relay_matrix,tensor)
        implicit none

        integer::ir,ic,i,j,i_row,i_col,n_atoms
        real*8,dimension(3,3)::tensor
        real*8,dimension(3*n_atoms,3*n_atoms)::temp_relay_matrix
        tensor(:,:)=0.0d0

        do ir=1,n_atoms,1
             do ic=1,n_atoms,1
                 i_row=0
                    do i=3*ir-2,3*ir,1
                       i_row=i_row+1
                       i_col=0
                          do j=3*ic-2,3*ic,1
                          i_col=i_col+1
                           tensor(i_row,i_col)=tensor(i_row,i_col) + temp_relay_matrix(i,j)
                          enddo
                    enddo
             enddo
        enddo

        return
endsubroutine contract_matrix_v1


!       subroutine contract_matrix_old(tensor)
!         implicit none
!         integer::i_index,j_index,i_row,i_col
!         real*8,dimension(3,3)::tensor
!         tensor=0.0d0


!         do i_row=1,n_atoms,1
!            do i_col=1,n_atoms,1
!               do i_index=1,3,1
!                  do j_index=1,3,1
                   
!                     tensor(i_index,j_index) = tensor(i_index,j_index) +& 
!                                               relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)
!                  enddo
!               enddo
!            enddo
!         enddo
!         return
!       endsubroutine contract_matrix_old

      subroutine get_coupled_atom_polarizability(iso_polar_coupled)
        implicit none

        integer::i_row,i_col,i_index,j_index
        real*8,dimension(3,3)::matrix
        real*8,dimension(n_atoms),intent(OUT)::iso_polar_coupled
        real*8 :: WORK(9),eigen(3)
        integer ::   LWORK,errorflag
        iso_polar_coupled=0.d0

        ! o polar_coupled is formed by summing either rows 'blocks' or column
        ! 'blocks' of relay_matrix_screened and contains full anisotropic atom
        ! resolved polarizability 

        do i_row=1,n_atoms,1
        matrix=0.0
          do i_col=1,n_atoms,1
             do i_index=1,3,1
                do j_index=1,3,1
                   matrix(i_index,j_index) = matrix(i_index,j_index) + relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)
                enddo
             enddo
          enddo
       !o iso_polar_coupled contains average of trances of atom resolved
       !polarizabilities in polar_coupled
        LWORK=9  
        call DSYEV('V','U',3,matrix,3,eigen,WORK,LWORK,errorflag)  
        call check_info(errorflag,"get_coupled_atom_polarizability","DSYEV")
        iso_polar_coupled(i_row) =(sum(eigen))/3.d0
        enddo
        return
      endsubroutine get_coupled_atom_polarizability



subroutine get_alpha_omega_and_Rp(HF,C6_free,alpha_free,w,a_w,Rpw)
       implicit none
       real*8,intent(in)::HF,C6_free,alpha_free,w
       real*8,intent(out)::a_w,Rpw
       real*8::eff_freq,w_eff_freq_ratio
       ! alpha(iomega)
       eff_freq=0.0d0
       w_eff_freq_ratio=0.0d0
       eff_freq= ((4.d0/3.d0)*C6_free/(alpha_free**2.0))**2
       w_eff_freq_ratio=(w*w)/eff_freq
       a_w=(HF*alpha_free)/(1.0+w_eff_freq_ratio)

       !Rp(iomega) ! Rp_eff definition from   A. Mayer, Phys. Rev. B, 75,045407(2007)
       Rpw= ((a_w/3.d0)*dsqrt(2.0d0/pi))**0.333333333333333333333333333d0

       return
endsubroutine get_alpha_omega_and_Rp



        ! ONLY vdW energy NO forces 
        subroutine get_ts_vdw(ene_pairwise_ts)

        real*8::Sr
        real*8,dimension(n_atoms):: alpha_array
        real*8,dimension(n_atoms):: R_array
        real*8,dimension(n_atoms):: C6_array 
        real*8,dimension(3)::d_xyz
        real*8,dimension(3)::coord_curr
        real*8 ::d,R01,R02,R0,C61,C62,alpha1,alpha2,C6,R_6,dist,R_gau,f_damp
        integer ::  i_atom,j_atom,i_loop,i1,i2,i3
        logical :: periodic_converged
        integer :: periodic_cell
        real*8  :: vdw_energy_change
        real*8  :: ene_pairwise_ts
 
        do i_loop=1,n_atoms,1
           call get_vdw_param(species_element(species(i_loop)),&
                              species_z(species(i_loop)),&
                              C6_array(i_loop),&
                              alpha_array(i_loop),&
                              R_array(i_loop))
        enddo

        d=20.d0
        select case (flag_xc)
          case (1) !PBE0
               Sr=0.96
          case (7) !HSE
               Sr=0.96
          case (6) ! PBE
               Sr=0.94
          case (9) ! BLYP 
               Sr=0.62
          case (10) ! B3LYP 
               Sr=0.84
          case (12) ! revPBE
                Sr=0.60
          case (20) ! AM05
               Sr=0.84
          case default
               Sr=1.0

          select case (flag_post_xc)
               case (1) ! M06L
                    Sr=1.26
               case (2) ! M06
                    Sr=1.16
               case default
          end select
        end select
  

        ene_pairwise_ts=0.d0
        do i_atom = 1,n_atoms,1
           alpha1 = alpha_array(i_atom)*hirshfeld_volume(i_atom)
              C61 = C6_array(i_atom)*hirshfeld_volume(i_atom)*hirshfeld_volume(i_atom)
              R01 = R_array(i_atom)*(hirshfeld_volume(i_atom)**0.333333333333333333333333333d0)
            if (myid.eq.task_list(i_atom)) then
                do j_atom = 1,n_atoms,1
           alpha2 = alpha_array(j_atom)*hirshfeld_volume(j_atom)
              C62 = C6_array(j_atom)*hirshfeld_volume(j_atom)*hirshfeld_volume(j_atom)
              R02 = R_array(j_atom)*(hirshfeld_volume(j_atom)**0.333333333333333333333333333d0)

                 if(C61.gt.0.d0 .or.C62.gt.0.d0) then
                         C6=2.d0*C61*C62/((alpha2/alpha1)*C61 + (alpha1/alpha2)*C62) 
                 else
                         C6=0.d0   
                 endif
                 dist=0.d0  
                 R0=Sr*(R01+R02)
                 d_xyz=coords(:,i_atom)-coords(:,j_atom)
                 dist=d_xyz(1)*d_xyz(1) + d_xyz(2)*d_xyz(2) +d_xyz(3)*d_xyz(3)

                 if(dist.gt.1e-6) then
                         R_6=dist*dist*dist  
                         f_damp=1.d0 + dexp(-d*((dsqrt(dist)/R0)-1.d0))  
                         ene_pairwise_ts = ene_pairwise_ts + (-C6/(R_6*f_damp))
                 endif 

                enddo ! end over i_atom
            endif
        enddo ! end over j_atom
        ene_pairwise_ts =0.5d0*ene_pairwise_ts
        call sync_vdw_correction_in_scs(ene_pairwise_ts) 
   
      if (n_periodic .gt. 0) then 
        periodic_converged = .false.
        periodic_cell      = 0
        do while (.not.periodic_converged)
           periodic_cell = periodic_cell + 1
           vdw_energy_change = 0.d0
           do i1 = -periodic_cell, periodic_cell
              do i2 = -periodic_cell, periodic_cell
                 do i3 = -periodic_cell, periodic_cell
                 if ((abs(i1).eq.periodic_cell).or.(abs(i2).eq.periodic_cell).or.(abs(i3).eq.periodic_cell)) then
                  do i_atom = 1,n_atoms,1
                     alpha1 = alpha_array(i_atom)*hirshfeld_volume(i_atom)
                        C61 = C6_array(i_atom)*hirshfeld_volume(i_atom)*hirshfeld_volume(i_atom)
                        R01 = R_array(i_atom)*(hirshfeld_volume(i_atom)**(1.d0/3.d0))
                     if (myid.eq.task_list(i_atom)) then
                     do j_atom = 1,n_atoms,1
                        alpha2 = alpha_array(j_atom)*hirshfeld_volume(j_atom)
                           C62 = C6_array(j_atom)*hirshfeld_volume(j_atom)*hirshfeld_volume(j_atom)
                           R02 = R_array(j_atom)*(hirshfeld_volume(j_atom)**(1.d0/3.d0))

                        if(C61.gt.0.d0 .or.C62.gt.0.d0) then
                          C6=2.d0*C61*C62/((alpha2/alpha1)*C61 +(alpha1/alpha2)*C62)
                        else
                         C6=0.d0
                        endif
                        dist=0.d0
                        R0=Sr*(R01+R02)
                        coord_curr(:) = coords(:,j_atom)+i1*lattice_vector(:,1)+&
                                        i2*lattice_vector(:,2)+i3*lattice_vector(:,3) 
                        d_xyz=coords(:,i_atom)- coord_curr(:)
                        dist=d_xyz(1)*d_xyz(1) + d_xyz(2)*d_xyz(2) +d_xyz(3)*d_xyz(3)

                        if(dist.gt.1e-6) then
                         R_6=dist*dist*dist
                         f_damp=1.d0 + dexp(-d*((dsqrt(dist)/R0)-1.d0))
                         vdw_energy_change = vdw_energy_change + (-C6/(R_6*f_damp))
                        endif

                   enddo ! end over i_atom
                  endif  ! task
                 enddo ! end over j_atom
                endif  ! shell
               enddo  ! i1
             enddo    !i2
           enddo    !i3

           vdw_energy_change = 0.5*vdw_energy_change ! factor of 0.5 due to complete ij summation
           call sync_vdw_correction_in_scs(vdw_energy_change)
           periodic_converged = (dabs(vdw_energy_change*hartree).lt.sc_accuracy_etot)
           ene_pairwise_ts = ene_pairwise_ts + vdw_energy_change

         enddo ! periodic converge
        endif    ! check periodic



        return 
        endsubroutine get_ts_vdw
        ! ONLY vdW energy NO forces 
        subroutine get_ts_scs_vdw(ene_pairwise_scs)
        ! local var  
        real*8::Sr
        real*8,dimension(n_atoms):: alpha_array
        real*8,dimension(n_atoms):: R_array
        real*8,dimension(3)::d_xyz
        real*8,dimension(3)::coord_curr
        real*8::d,R01,R02,R0,C61,C62,alpha1,alpha2,C6,R_6,dist,R_gau,f_damp,TEMP_VAR
        integer :: i_atom,j_atom,i_loop,i1,i2,i3
        logical :: periodic_converged
        integer :: periodic_cell
        real*8  :: vdw_energy_change
        real*8  :: ene_pairwise_scs
        do i_loop=1,n_atoms,1
        call get_vdw_param(species_element(species(i_loop)),&
                           species_z(species(i_loop)),&
                           TEMP_VAR,alpha_array(i_loop),&
                           R_array(i_loop))
        enddo


        d=20.d0
        select case (flag_xc)
          case (1) !PBE0
               Sr=0.96
          case (7) !HSE
               Sr=0.96
          case (6) ! PBE
               Sr=0.94
          case (9) ! BLYP 
               Sr=0.62
          case (10) ! B3LYP 
               Sr=0.84
          case (12) ! revPBE
                Sr=0.60
          case (20) ! AM05
               Sr=0.84
          case default
               Sr=1.0

          select case (flag_post_xc)
               case (1) ! M06L
                    Sr=1.26
               case (2) ! M06
                    Sr=1.16
               case default
          end select
        end select


        ene_pairwise_scs=0.d0
        do i_atom = 1,n_atoms,1
           alpha1 = alpha_eff(i_atom)
              C61 = C6_eff(i_atom)
              R01 = R_array(i_atom)*(alpha_eff(i_atom)/alpha_array(i_atom))**(1.d0/3.d0)
 
            if (myid.eq.task_list(i_atom)) then
                do j_atom = 1,n_atoms,1
                   alpha2 = alpha_eff(j_atom)
                      C62 = C6_eff(j_atom)
                      R02 = R_array(j_atom)*(alpha_eff(j_atom)/alpha_array(j_atom))**(1.d0/3.d0)
                 if(C61.gt.0.d0 .or.C62.gt.0.d0) then
                         C6=2.d0*C61*C62/((alpha2/alpha1)*C61 +(alpha1/alpha2)*C62)
                 else
                         C6=0.d0
                 endif
                 dist=0.d0
                 R0=Sr*(R01+R02)
                 d_xyz=coords(:,i_atom)-coords(:,j_atom)
                 dist=d_xyz(1)*d_xyz(1) + d_xyz(2)*d_xyz(2) +d_xyz(3)*d_xyz(3)

                 if(dist.gt.1e-6) then
                         R_6=dist*dist*dist
                         f_damp=1.d0 + dexp(-d*((dsqrt(dist)/R0)-1.d0))
                         ene_pairwise_scs = ene_pairwise_scs + (-C6/(R_6*f_damp))
                 endif

                enddo ! end over i_atom
            endif
        enddo ! end over j_atom
        ene_pairwise_scs =0.5d0*ene_pairwise_scs
        call sync_vdw_correction_in_scs(ene_pairwise_scs)

      if (n_periodic .gt. 0) then
        periodic_converged = .false.
        periodic_cell      = 0
        do while (.not.periodic_converged)
           periodic_cell = periodic_cell + 1
           vdw_energy_change = 0.d0

           do i1 = -periodic_cell, periodic_cell
              do i2 = -periodic_cell, periodic_cell
                 do i3 = -periodic_cell, periodic_cell
                 if((abs(i1).eq.periodic_cell).or.(abs(i2).eq.periodic_cell).or.(abs(i3).eq.periodic_cell))then
                  do i_atom = 1,n_atoms,1
                     alpha1 = alpha_eff(i_atom)
                        C61 = C6_eff(i_atom)
                        R01 = R_array(i_atom)*(alpha_eff(i_atom)/alpha_array(i_atom))**(1.d0/3.d0)

                     if (myid.eq.task_list(i_atom)) then
                     do j_atom = 1,n_atoms,1
                        alpha2 = alpha_eff(j_atom)
                           C62 = C6_eff(j_atom)
                           R02 = R_array(j_atom)*(alpha_eff(j_atom)/alpha_array(j_atom))**(1.d0/3.d0)

                        if(C61.gt.0.d0 .or.C62.gt.0.d0) then
                          C6=2.d0*C61*C62/((alpha2/alpha1)*C61+(alpha1/alpha2)*C62)
                        else
                         C6=0.d0
                        endif
                        dist=0.d0
                        R0=Sr*(R01+R02)
                        coord_curr(:) = coords(:,j_atom)+i1*lattice_vector(:,1)+&
                                        i2*lattice_vector(:,2)+i3*lattice_vector(:,3)
                        d_xyz=coords(:,i_atom)- coord_curr(:)
                        dist=d_xyz(1)*d_xyz(1) + d_xyz(2)*d_xyz(2)+d_xyz(3)*d_xyz(3)

                        if(dist.gt.1e-6) then
                         R_6=dist*dist*dist
                         f_damp=1.d0 + dexp(-d*((dsqrt(dist)/R0)-1.d0))
                         vdw_energy_change = vdw_energy_change +(-C6/(R_6*f_damp))
                        endif

                   enddo ! end over j_atom
                  endif  ! tasks
                 enddo ! end over i_atom
                endif  ! shell
               enddo  ! i1
             enddo    !i2
           enddo    !i3

           vdw_energy_change = 0.5*vdw_energy_change 
           call sync_vdw_correction_in_scs(vdw_energy_change)
           periodic_converged =(dabs(vdw_energy_change*hartree).lt.sc_accuracy_etot)
           ene_pairwise_scs = ene_pairwise_scs + vdw_energy_change

         enddo ! periodic converge
        endif  ! check periodic
        return
        endsubroutine get_ts_scs_vdw




      subroutine TPP_TENSOR(dxyz, r_ij, r_pp, TPP)
      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: r_pp
      real*8,dimension(3,3),intent(out) :: TPP
      ! This needs declarding for the PGI Architecture
      real*8 :: derf 

      ! local vars
      real*8,dimension(3,3) :: TPQ,r_tensor
      real*8:: zeta_l,zeta_r,ratio
      integer :: i,j

      !o dipole-dipole interaction tensor between two quantum harmonic
      !  oscilators Tp-p {see.. A. Mayer, Phys. Rev. B, 75,045407(2007)}
      TPP(:,:)=0.d0  
      ratio=r_ij/r_pp
      zeta_l=(derf(ratio)-(dsqrt(4.0d0/pi) * ratio *dexp(-(ratio**2.0d0))))/r_ij**5.0d0
      zeta_r=(dsqrt(16.0d0/pi) * dexp(-(ratio**2.0d0)))/(r_pp * r_pp * r_pp *r_ij * r_ij)

                         !!Tensor product
                          do i=1,3,1
                            do j=1,3,1
                              r_tensor(i,j)=dxyz(i)*dxyz(j)
                            enddo
                          enddo

                          TPQ=r_tensor*3.0d0
                          do i=1,3,1
                             TPQ(i,i)=TPQ(i,i)-(r_ij*r_ij)
                          enddo

                          TPQ=TPQ*zeta_l
                          r_tensor=r_tensor * zeta_r
                          TPP=TPQ-r_tensor
                          TPP=-1.d0*TPP
      return
      endsubroutine



      subroutine SCS_TENSOR_MBD_rsSCS(dxyz, r_ij, r_pp,Rvdw12,beta, TPP)
      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: r_pp
      real*8,intent(in) :: Rvdw12
      real*8,intent(in) :: beta
      real*8,dimension(3,3),intent(out) :: TPP
      ! This needs declarding for the PGI Architecture
      real*8 :: derf, erf

      ! local vars
      real*8,dimension(3,3) :: TPQ,r_tensor
      real*8:: zeta_l,zeta_r,ratio,d_param,fermi_param
      integer :: i,j

      !o dipole-dipole interaction tensor between two quantum harmonic
      !  oscilators Tp-p damped with (1-erf(rij/Rvdwij))
      TPP(:,:)=0.d0
      ! Fermi damping parameter 
      d_param = 6.0
      fermi_param = r_ij/(beta*Rvdw12)
      !!!!!!!!!!!!!!!!!!!!!!!!!  
      ratio=r_ij/r_pp
      zeta_l=(derf(ratio)-(dsqrt(4.0d0/pi) * ratio*dexp(-(ratio**2.0d0))))/r_ij**5.0d0
      zeta_r=(dsqrt(16.0d0/pi) * dexp(-(ratio**2.0d0)))/(r_pp * r_pp * r_pp*r_ij * r_ij)

                         !!Tensor product
                          do i=1,3,1
                            do j=1,3,1
                              r_tensor(i,j)=dxyz(i)*dxyz(j)
                            enddo
                          enddo

                          TPQ=r_tensor*3.0d0
                          do i=1,3,1
                             TPQ(i,i)=TPQ(i,i)-(r_ij*r_ij)
                          enddo

                          TPQ=TPQ*zeta_l
                          r_tensor=r_tensor * zeta_r
                          TPP=TPQ-r_tensor
                          TPP=-1.d0*TPP
                          ! Old erf type damping 
!!!!!!!!!!!!!!!!!!!!!!    TPP=TPP*(1.d0-erf((r_ij/(Rvdw12*beta))**4))
 
                          ! Fermi type damping  
                          TPP=TPP*(1.0-(1.0/(1.0+exp(-d_param*(fermi_param-1.0)))))
      return
      endsubroutine SCS_TENSOR_MBD_rsSCS
      
     subroutine get_many_body_dispersion_energy_cfdm(ene_mbd_scs)

     ! local vars  
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian 
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian_periodic
     real*8,dimension(:,:),allocatable:: coords_SL
     real*8,dimension(:),allocatable:: cfdm_eigenvalues
     real*8,dimension(:),allocatable:: omega_cfdm_SL
     real*8,dimension(:),allocatable:: R_vdw_SL
     real*8,dimension(:),allocatable::alpha_eff_SL
     real*8,dimension(:),allocatable:: WORK
     real*8,dimension(:),allocatable:: task_list_SL  
     real*8,dimension(3,3):: TPP,lattice_vector_SL
     real*8,dimension(3)::dxyz,coord_curr
     real*8:: r_ij
     real*8:: TEMP_VAR
     real*8:: alpha_free
     real*8:: R_vdw_free
     real*8:: Rvdw_12
     real*8:: beta
     real*8:: CFDM_prefactor
     real*8:: E_int
     real*8:: E_nonint
     real*8:: sigma
     real*8:: ene_mbd_scs
     integer :: errorflag,i_atom,j_atom
     integer :: SL_i,SL_j,SL_k ! MBD super cell index's
     integer :: i_index, j_index ! 3x3 block index
     integer :: i_lattice,j_lattice,k_lattice ! lattice increament index
     integer :: periodic_cell_i,periodic_cell_j,periodic_cell_k
     integer :: NIMAG ! Number of imaginary frequencies in CFDM if any
     integer ::n_atoms_SL
  
     ! LAPACK VARs
     integer :: LWORK
!     Begin work       
        select case (flag_xc)
          case (1) !PBE0
               beta=2.53
          case (6) ! PBE
               beta=2.56
          case (7) !HSE 
               beta=2.53
          case default
               beta=0.0
               write(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"WITHOUT short range damping"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
        end select

      if ( n_periodic .eq. 0) then
      ! For Cluster calculation NO super cell needed 
      ! intiate this index so as to avoide any invalid floting point when
      ! normalizing energy in case of cluster/molecule
      SL_i = 1 
      SL_j = 1
      SL_k = 1

!     number of atoms in cluster /molecule
      n_atoms_SL = n_atoms
 
      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK((3*n_atoms_SL)*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))

      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian=0.0d0 
      cfdm_eigenvalues=0.0d0

!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN CLUSTER/MOLECULE
      do i_atom=1,n_atoms_SL
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks)  
      enddo

          do i_atom =1,n_atoms
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0
                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                TEMP_VAR,alpha_free,R_vdw_free)
                  R_vdw_SL(i_atom)= R_vdw_free*((alpha_eff(i_atom)/alpha_free)**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  coords_SL(:,i_atom)   = coords(:,i_atom)
                  alpha_eff_SL(i_atom) = alpha_eff(i_atom)
                  else 
                  !! In case when parameter(s) C_6, alpha, R_vdw are specified via hirshfeld_param keyword in control.in
                  R_vdw_SL(i_atom)= vdw_hirshfeld_R0(species(i_atom))* & 
                     ((alpha_eff(i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  coords_SL(:,i_atom)   = coords(:,i_atom)
                  alpha_eff_SL(i_atom) = alpha_eff(i_atom)

                  endif 
         enddo
      else
      if(.NOT.mbd_scs_vacuum_axis(1)) then
      SL_i = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2)) 
      else
      SL_i=1
      endif  
      if(.NOT.mbd_scs_vacuum_axis(2)) then
      SL_j = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)) 
      else
      SL_j=1
      endif  
      if(.NOT.mbd_scs_vacuum_axis(3)) then
      SL_k = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)) 
      else
      SL_k=1
      endif  
       
!     number of atoms in super cell 
      n_atoms_SL = n_atoms*SL_i*SL_j*SL_k

      write(info_str,'(2X,A,i3,A,i3,A,i3,A,f6.2,A)')&
     "| Creating super cell of dimension",  SL_i," X ", SL_j," X ", SL_k,&
     " in MBD calculation using" ,mbd_supercell_cutoff, "  Angstrom radius"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,"| containing", I6,A)')n_atoms_SL, "  atom"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

 

      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_hamiltonian_periodic)) allocate(cfdm_hamiltonian_periodic(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK(3*n_atoms_SL*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))
      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian = 0.0d0 
      cfdm_eigenvalues = 0.0d0
      cfdm_hamiltonian_periodic = 0.d0 
!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN SUPERCELL/MOLECULE
      do i_atom=1,n_atoms_SL 
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks) 
      enddo

      ! Get all the parameter  required for MBD with super cell calculation
       j_atom = 0
       do i_lattice = 0, SL_i-1,1
         do j_lattice = 0, SL_j-1,1
            do k_lattice = 0, SL_k-1,1
              do i_atom =1,n_atoms
                  j_atom = j_atom +1  
                  coords_SL(:,j_atom) = coords(:,i_atom) + i_lattice*lattice_vector(:,1) + &
                                                           j_lattice*lattice_vector(:,2) + &
                                                           k_lattice*lattice_vector(:,3)   
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0  

                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                TEMP_VAR,alpha_free,R_vdw_free)
                  R_vdw_SL(j_atom)= R_vdw_free*((alpha_eff(i_atom)/alpha_free)**(1.d0/3.d0))  
                  alpha_eff_SL(j_atom) = alpha_eff(i_atom)
                  omega_cfdm_SL(j_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  else
                  !! In case when parameter(s) C_6, alpha, R_vdw are specified via hirshfeld_param keyword in control.in
                  R_vdw_SL(j_atom)= vdw_hirshfeld_R0(species(i_atom))* & 
                      ((alpha_eff(i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  alpha_eff_SL(j_atom) = alpha_eff(i_atom)
                  omega_cfdm_SL(j_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  endif 

              enddo 
            enddo            
         enddo            
       enddo            
      lattice_vector_SL(:,1)  = lattice_vector(:,1)*SL_i 
      lattice_vector_SL(:,2)  = lattice_vector(:,2)*SL_j 
      lattice_vector_SL(:,3)  = lattice_vector(:,3)*SL_k 
      endif 

      ! Construct cfdm hamiltonian     
      do i_atom=1,n_atoms_SL,1 !$1
         if(myid.eq.task_list_SL(i_atom)) then
         do j_atom=i_atom,n_atoms_SL,1 !$2

         if(i_atom.eq.j_atom) then  !#1

               do i_index=1,3,1
                  do j_index=1,3,1
                     if(i_index.eq.j_index) then
                        cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = omega_cfdm_SL(i_atom)**2.0
                     else
                        cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = 0.0d0
                     endif
                  enddo
               enddo

         else
            r_ij=0.0d0
            TPP=0.0d0
            dxyz=0.d0
            dxyz(:)= coords_SL(:,i_atom)-coords_SL(:,j_atom)
            r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
            Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
            sigma=(r_ij/Rvdw_12)**beta
            call CFDM_TENSOR(dxyz,r_ij,sigma,beta,TPP)
            CFDM_prefactor=omega_cfdm_SL(i_atom)*omega_cfdm_SL(j_atom)* &
                           dsqrt(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))

                   ! Transfer each dipole matrix to CFDM hamiltonian i.j accordingly         
               do i_index=1,3,1
                  do j_index=1,3,1
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = TPP(i_index,j_index)*CFDM_prefactor
                   cfdm_hamiltonian(3*j_atom-3+j_index,3*i_atom-3+i_index) = TPP(i_index,j_index)*CFDM_prefactor
                  enddo
               enddo

         endif !#1
         enddo !$2
       endif ! tasks 
      enddo !$1
      call sync_tensors(cfdm_hamiltonian,3*n_atoms_SL)

! Adds dipole field due to image cells based spherical cutoff mbd_cfdm_dip_cutoff
if (n_periodic .gt. 0) then 

      if(.NOT.mbd_scs_vacuum_axis(1)) then
      periodic_cell_i = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,1)**2 +&
                                                                lattice_vector_SL(2,1)**2 +&
                                                                lattice_vector_SL(3,1)**2)) 
      else
      periodic_cell_i=1
      endif

      if(.NOT.mbd_scs_vacuum_axis(2)) then 
      periodic_cell_j = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,2)**2 +&
                                                                lattice_vector_SL(2,2)**2 +&
                                                                lattice_vector_SL(3,2)**2))
      else
      periodic_cell_j=1
      endif
      if(.NOT.mbd_scs_vacuum_axis(3)) then 
      periodic_cell_k = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,3)**2 +&
                                                                lattice_vector_SL(2,3)**2 +&
                                                                lattice_vector_SL(3,3)**2))
      else
      periodic_cell_k=1
      endif

       
      do i_lattice = -periodic_cell_i, periodic_cell_i,1          !$7
         do j_lattice = -periodic_cell_j, periodic_cell_j,1         !$6  
            do k_lattice = -periodic_cell_k, periodic_cell_k,1        !$5 
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then!$4
                do i_atom=1,n_atoms_SL,1 !$3
                 if(myid.eq.task_list_SL(i_atom)) then !$ tasks
                  do j_atom=i_atom,n_atoms_SL,1 !$2

                      r_ij=0.0d0
                      TPP=0.0d0
                      dxyz=0.d0
                      coord_curr(:) = coords_SL(:,i_atom) + i_lattice*lattice_vector_SL(:,1) + &
                                                            j_lattice*lattice_vector_SL(:,2) + &
                                                            k_lattice*lattice_vector_SL(:,3)

                      dxyz(:)= coords_SL(:,j_atom)- coord_curr(:)
                      r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )

                      if(r_ij.le.mbd_cfdm_dip_cutoff/bohr)  then
                        Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
                        sigma=(r_ij/Rvdw_12)**beta
                        call CFDM_TENSOR(dxyz,r_ij,sigma,beta,TPP)
                        CFDM_prefactor=omega_cfdm_SL(i_atom)*&
                                       omega_cfdm_SL(j_atom)*&
                                       sqrt(alpha_eff_SL(i_atom)*&
                                            alpha_eff_SL(j_atom))
                       
                        
                        do i_index=1,3,1
                          do j_index=1,3,1
                          !FILL UPPER BLOCK(s)
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)=&
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          !FILL LOWER BLOCK(s) ALSO MAKE SURE THAT YOU DONT ADD FIELD
                          !CONTRIBUTION TWO TIMES IN DIAGONAL BLOCK below is the
                          !check
                          if(i_atom.NE.j_atom) then
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)=&
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          endif

                         enddo
                        enddo
                      endif 

                  enddo !$2
                 endif !$tasks
                enddo !$3
               endif !$4
            enddo !$5
         enddo !$6
      enddo !$7
   call sync_tensors(cfdm_hamiltonian_periodic,3*n_atoms_SL)    
   cfdm_hamiltonian = cfdm_hamiltonian + cfdm_hamiltonian_periodic
endif

    errorflag=0
    LWORK=3*n_atoms_SL*(3+(3*n_atoms_SL)/2)
    call print_dsyev_info("get_many_body_dispersion_energy_cfdm",3*n_atoms_SL)
    call DSYEV('V','U',3*n_atoms_SL,cfdm_hamiltonian,3*n_atoms_SL,cfdm_eigenvalues,WORK,LWORK,errorflag)
    call check_info(errorflag,"get_many_body_dispersion_energy_cfdm","DSYEV")

    E_int=0.0d0
    E_nonint =0.0d0
    ene_mbd_scs=   0.0
    NIMAG=0 
    do i_atom =1,n_atoms_SL
      E_nonint = E_nonint + omega_cfdm_SL(i_atom)    
    enddo

    do i_atom =1,3*n_atoms_SL
      if(cfdm_eigenvalues(i_atom).ge.0.d0) then
        E_int = E_int + dsqrt(cfdm_eigenvalues(i_atom))
      else
        NIMAG= NIMAG +1  
      endif 
    enddo
    
   !renormalise interation energies based number of cells used for supercell  
   ene_mbd_scs = ((0.5*E_int)-(1.5* E_nonint))/(SL_i*SL_j*SL_k)
   if(NIMAG.gt.0) then
     write(info_str,'(A,I4,A)')"***WARNING: found ",NIMAG," negative eigenvalues in MBD_SCS energy calculation."
     call localorb_info(info_str, use_unit,'(A)',OL_norm)  
   endif
   !
   if(allocated(cfdm_hamiltonian))          deallocate(cfdm_hamiltonian)
   if(allocated(cfdm_hamiltonian_periodic)) deallocate(cfdm_hamiltonian_periodic) 
   if(allocated(cfdm_eigenvalues))          deallocate(cfdm_eigenvalues)
   if(allocated(coords_SL))                 deallocate(coords_SL)
   if(allocated(omega_cfdm_SL))             deallocate(omega_cfdm_SL)
   if(allocated(R_vdw_SL))                  deallocate(R_vdw_SL)
   if(allocated(alpha_eff_SL))              deallocate(alpha_eff_SL)  
   
   return
  
 endsubroutine get_many_body_dispersion_energy_cfdm

      subroutine CFDM_TENSOR(dxyz,r_ij,sigma,beta,Tij)
      ! o Dipole tensor derived from short range screened coulomb potential
      ! as described in Tkatchenko et al, PRL 108, 236402 (2012)

      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: sigma
      real*8,intent(in) :: beta
      real*8,dimension(3,3),intent(out) :: Tij
      ! local vars
      real*8 :: zeta_l
      real*8 :: zeta_r
      integer :: i_index, j_index

      zeta_l=(3.d0 - (3.d0*EXP(-sigma)) - &
             (4.d0*EXP(-sigma)*beta*sigma) + &
             (EXP(-sigma)*beta*beta*sigma) - &
             (EXP(-sigma)*beta*beta*sigma*sigma))/(r_ij**5)
      zeta_r=(EXP(-sigma) + (EXP(-sigma)*beta*sigma) -1.d0)/(r_ij**3)

             do i_index=1,3,1
               do j_index=1,3,1
                Tij(i_index,j_index)=dxyz(i_index)*dxyz(j_index)
               enddo
             enddo
             Tij = Tij*zeta_l

             do i_index=1,3,1
             Tij(i_index,i_index) = Tij(i_index,i_index) + zeta_r
             enddo

             Tij=-1.d0*Tij
      return
      endsubroutine CFDM_TENSOR

    subroutine get_alpha_omega_relay_matrix()
    integer ::i_index,j_index

    real*8,dimension(3,3)::TPP
    real*8,dimension(3) :: dxyz
    real*8,dimension(3) :: coord_curr
    real*8 :: r_ij
    real*8 :: r_pp
    integer :: i_row, i_col
    integer :: i_lattice, j_lattice, k_lattice
    integer :: errorflag,periodic_cell_i,periodic_cell_j,periodic_cell_k

    !For LAPACK
    integer,dimension(3*n_atoms):: IPIV
    real*8,dimension(3*n_atoms):: WORK

    ! initio values
    relay_matrix=0.0d0
     ! compute relay matrix of  cluster or unit cell
       do i_row=1,n_atoms,1 !#1
        if(myid.eq.task_list(i_row)) then
         do i_col=i_row,n_atoms,1 !#2
         TPP=0.d0
            if(i_row.eq.i_col) then  !$1
               do i_index=1,3,1
                  do j_index=1,3,1
                   if(i_index.eq.j_index) then
                   relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=1.d0/alpha_omega(i_row)
                   else                    
                   relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=0.d0
                   endif 
                  enddo
               enddo

            else
               dxyz(:) = coords(:,i_col)-coords(:,i_row)
               r_ij = dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
               r_pp = dSqrt(R_p(i_row)**2 + R_p(i_col)**2)
               call TPP_TENSOR(dxyz,r_ij,r_pp,TPP)
               do i_index=1,3,1
                  do j_index=1,3,1
                   relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=TPP(i_index,j_index)
                   relay_matrix(3*i_col-3+j_index,3*i_row-3+i_index)=TPP(i_index,j_index)
                  enddo
               enddo
            endif !$1
         enddo   !#2
        endif ! task
       enddo  !#1
  call sync_tensors(relay_matrix,3*n_atoms)

  if (n_periodic .gt. 0) then
 
      relay_matrix_periodic=0.d0 
      if(.NOT.mbd_scs_vacuum_axis(1)) then
      periodic_cell_i = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 +& 
                                                               lattice_vector(2,1)**2 +&
                                                               lattice_vector(3,1)**2))
      else
      periodic_cell_i=1
      endif
         
      if(.NOT.mbd_scs_vacuum_axis(2)) then
      periodic_cell_j = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 +& 
                                                               lattice_vector(2,2)**2 +&
                                                               lattice_vector(3,2)**2))
      else
      periodic_cell_j=1
      endif
      if(.NOT.mbd_scs_vacuum_axis(3)) then
      periodic_cell_k = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 +& 
                                                               lattice_vector(2,3)**2 +&
                                                               lattice_vector(3,3)**2))
      else
      periodic_cell_k=1
      endif
          
      do i_lattice = -periodic_cell_i, periodic_cell_i,1
         do j_lattice = -periodic_cell_j, periodic_cell_j,1
            do k_lattice = -periodic_cell_k, periodic_cell_k,1
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then !#1
                  do i_row = 1, n_atoms, 1 ! atom1 loop
                    if(myid.eq.task_list(i_row)) then
                     do i_col = i_row, n_atoms, 1 ! atom2 loop
                          coord_curr = 0.d0
                                dxyz = 0.d0
                                r_ij = 0.d0
                                TPP  = 0.d0

                          ! find the coordinate of images
                              
                          coord_curr(:) = coords(:,i_col) + i_lattice*lattice_vector(:,1) + &
                                          j_lattice*lattice_vector(:,2) + k_lattice*lattice_vector(:,3)

                          dxyz(:) = coords(:,i_row)- coord_curr(:)
                          r_ij    = sqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                          r_pp    = sqrt(R_p(i_row)**2 + R_p(i_col)**2)
                          if(r_ij.le.mbd_scs_dip_cutoff/bohr) then
                            !
                            call TPP_TENSOR(dxyz,r_ij,r_pp,TPP)
                            do i_index=1,3,1
                            do j_index=1,3,1
                            relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index)=&
                            relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index) + TPP(i_index,j_index)
                                   
                                 if(i_row.ne.i_col) then
                                 ! Since loop over atoms goes over upper blocks
                                 ! of matrix so fill lower blocks if necessary
                                 ! as LAPACK call for sysmetric matrix need either upper or lower triangle of matrix
                                   relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index)=&
                                   relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index)+ TPP(i_index,j_index)
                                 endif 
 
                            enddo
                            enddo
                          endif
                     enddo !atom2 loop
                   endif! parallel
                  enddo !atom1 loop
               endif  !#1
            enddo
         enddo
      enddo
   call sync_tensors(relay_matrix_periodic,3*n_atoms)
   relay_matrix = relay_matrix + relay_matrix_periodic  
   endif

   call DGETRF(3*n_atoms, 3*n_atoms, relay_matrix, 3*n_atoms, IPIV, errorflag)
   if(errorflag.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call localorb_info(info_str, use_unit,'(A)',OL_norm)
   endif
   call check_info(errorflag,"get_alpha_omega_relay_matrix","DGETRF")


   call DGETRI(3*n_atoms, relay_matrix, 3*n_atoms, IPIV, WORK,3*n_atoms,errorflag )
   if(errorflag.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call localorb_info(info_str, use_unit,'(A)',OL_norm)
   endif
   call check_info(errorflag,"get_alpha_omega_relay_matrix","DGETRI")
   
   return
 endsubroutine  get_alpha_omega_relay_matrix

subroutine gauss_legendre_grid(casmir_omega,casmir_omega_weight)
         implicit none 
         real*8,dimension(0:20)  :: casmir_omega,casmir_omega_weight 
         casmir_omega(0) = 0.0000000 ;  casmir_omega_weight(0) = 0.0000000
         casmir_omega(1) = 0.0392901 ;  casmir_omega_weight(1) = 0.0786611
         casmir_omega(2) = 0.1183580 ;  casmir_omega_weight(2) = 0.0796400
         casmir_omega(3) = 0.1989120 ;  casmir_omega_weight(3) = 0.0816475
         casmir_omega(4) = 0.2820290 ;  casmir_omega_weight(4) = 0.0847872
         casmir_omega(5) = 0.3689190 ;  casmir_omega_weight(5) = 0.0892294
         casmir_omega(6) = 0.4610060 ;  casmir_omega_weight(6) = 0.0952317
         casmir_omega(7) = 0.5600270 ;  casmir_omega_weight(7) = 0.1031720
         casmir_omega(8) = 0.6681790 ;  casmir_omega_weight(8) = 0.1136050
         casmir_omega(9) = 0.7883360 ;  casmir_omega_weight(9) = 0.1273500
        casmir_omega(10) = 0.9243900 ; casmir_omega_weight(10) = 0.1456520
        casmir_omega(11) = 1.0817900 ; casmir_omega_weight(11) = 0.1704530
        casmir_omega(12) = 1.2684900 ; casmir_omega_weight(12) = 0.2049170
        casmir_omega(13) = 1.4966100 ; casmir_omega_weight(13) = 0.2544560
        casmir_omega(14) = 1.7856300 ; casmir_omega_weight(14) = 0.3289620
        casmir_omega(15) = 2.1691700 ; casmir_omega_weight(15) = 0.4480920
        casmir_omega(16) = 2.7106200 ; casmir_omega_weight(16) = 0.6556060
        casmir_omega(17) = 3.5457300 ; casmir_omega_weight(17) = 1.0659600
        casmir_omega(18) = 5.0273400 ; casmir_omega_weight(18) = 2.0635700
        casmir_omega(19) = 8.4489600 ; casmir_omega_weight(19) = 5.6851000
        casmir_omega(20) = 25.451700 ; casmir_omega_weight(20) = 50.955800
return
endsubroutine gauss_legendre_grid

  subroutine get_mbd_enengy_at_TS(ene_mbd_ts)

     ! local vars  
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian 
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian_periodic
     real*8,dimension(:,:),allocatable:: coords_SL
     real*8,dimension(:),allocatable:: cfdm_eigenvalues
     real*8,dimension(:),allocatable:: omega_cfdm_SL
     real*8,dimension(:),allocatable:: R_vdw_SL
     real*8,dimension(:),allocatable::alpha_eff_SL
     real*8,dimension(:),allocatable:: WORK
     real*8,dimension(:),allocatable:: task_list_SL  
     real*8,dimension(3,3):: TPP,lattice_vector_SL
     real*8,dimension(3)::dxyz,coord_curr
     real*8:: r_ij
     real*8:: C6_free
     real*8:: alpha_free
     real*8:: R_vdw_free
     real*8:: Rvdw_12
     real*8:: beta
     real*8:: CFDM_prefactor
     real*8:: E_int
     real*8:: E_nonint
     real*8:: sigma
     real*8:: ene_mbd_ts
     integer :: errorflag,i_atom,j_atom
     integer :: SL_i,SL_j,SL_k ! MBD super cell index's
     integer :: i_index, j_index ! 3x3 block index
     integer :: i_lattice,j_lattice,k_lattice ! lattice increament index
     integer :: periodic_cell_i,periodic_cell_j,periodic_cell_k
     integer :: NIMAG
     integer ::n_atoms_SL
 
     ! LAPACK VARs
     integer :: LWORK
!     Begin work       
      write(info_str,'(2x,2A)') "|-------------------------------------", &
             "---------------------------------------------"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,A)')"| Many-Body Dispersion (MBD@TS) energy "
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,2A)') "|-------------------------------------", &
             "---------------------------------------------"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

        select case (flag_xc)
          case (1) !PBE0
               beta=1.08
          case (6) ! PBE
               beta=1.07
          case (7) !HSE 
               beta=1.08

          case default
               beta=1.0
               write(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"WITHOUT short range damping"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
        end select

      if ( n_periodic .eq. 0) then
      ! For Cluster calculation NO super cell needed 
      ! intiate this index so as to avoide any invalid floting point when
      ! normalizing energy in case of cluster/molecule
      SL_i = 1 
      SL_j = 1
      SL_k = 1

!     number of atoms in cluster /molecule
      n_atoms_SL = n_atoms
 
      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK((3*n_atoms_SL)*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))

      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian=0.0d0 
      cfdm_eigenvalues=0.0d0

!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN CLUSTER/MOLECULE
      do i_atom=1,n_atoms_SL
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks)  
      enddo

          do i_atom =1,n_atoms_SL
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0
                  C6_free=0.d0
                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then  
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                C6_free,alpha_free,R_vdw_free)
                  R_vdw_SL(i_atom)= R_vdw_free*(hirshfeld_volume(i_atom)**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6_free/(alpha_free**2.0)
                  coords_SL(:,i_atom)   = coords(:,i_atom)
                  alpha_eff_SL(i_atom)  = hirshfeld_volume(i_atom) * alpha_free
                  else
                  !! Case when parameter(s) C_6, alpha and R_vdw are specified in control.in
                  R_vdw_SL(i_atom)= vdw_hirshfeld_R0(species(i_atom))*(hirshfeld_volume(i_atom)**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*vdw_hirshfeld_C6(species(i_atom))/(vdw_hirshfeld_alpha(species(i_atom))**2.0)
                  coords_SL(:,i_atom)   = coords(:,i_atom)
                  alpha_eff_SL(i_atom)  = hirshfeld_volume(i_atom) * vdw_hirshfeld_alpha(species(i_atom))

                  endif 
         enddo
      else
      if(.NOT.mbd_scs_vacuum_axis(1)) then  
      SL_i = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2)) 
      else
      SL_i =1  
      endif
      if(.NOT.mbd_scs_vacuum_axis(2)) then  
      SL_j = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)) 
      else
      SL_j =1  
      endif
      if(.NOT.mbd_scs_vacuum_axis(3)) then  
      SL_k = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)) 
      else
      SL_k =1  
      endif
       
!     number of atoms in super cell 
      n_atoms_SL = n_atoms*SL_i*SL_j*SL_k
      write(info_str,'(2X,A,i3,A,i3,A,i3,A,f6.2,A)')&
     "| Creating super cell of dimension",  SL_i," X ", SL_j," X ", SL_k,&
     " in MBD calculation using" ,mbd_supercell_cutoff, "  Angstrom radius"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,"| containing", I6,A)')n_atoms_SL, "  atom"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)


      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_hamiltonian_periodic)) allocate(cfdm_hamiltonian_periodic(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK(3*n_atoms_SL*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))
      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian = 0.0d0 
      cfdm_eigenvalues = 0.0d0
      cfdm_hamiltonian_periodic = 0.d0 
!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN SUPERCELL/MOLECULE
      do i_atom=1,n_atoms_SL 
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks) 
      enddo

      ! Get all the parameter  required for MBD with super cell calculation
       j_atom = 0
       do i_lattice = 0, SL_i-1,1
         do j_lattice = 0, SL_j-1,1
            do k_lattice = 0, SL_k-1,1
              do i_atom =1,n_atoms
                  j_atom = j_atom +1  
                  coords_SL(:,j_atom) = coords(:,i_atom) + i_lattice*lattice_vector(:,1) + &
                                                           j_lattice*lattice_vector(:,2) + &
                                                           k_lattice*lattice_vector(:,3)   
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0  

                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then 
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                C6_free,alpha_free,R_vdw_free)

                  R_vdw_SL(j_atom)      = R_vdw_free*(hirshfeld_volume(i_atom)**(1.d0/3.d0))
                  alpha_eff_SL(j_atom)  = hirshfeld_volume(i_atom) * alpha_free
                  omega_cfdm_SL(j_atom) = (4.d0/3.d0)*C6_free/(alpha_free**2.0)
                  else
                  !! Case when parameter(s) C_6, alpha and R_vdw are specified in control.in  
                  R_vdw_SL(j_atom)= vdw_hirshfeld_R0(species(i_atom))*(hirshfeld_volume(i_atom)**(1.d0/3.d0))
                  omega_cfdm_SL(j_atom) = (4.d0/3.d0)*vdw_hirshfeld_C6(species(i_atom))/(vdw_hirshfeld_alpha(species(i_atom))**2.0)
                  alpha_eff_SL(j_atom)  = hirshfeld_volume(i_atom) * vdw_hirshfeld_alpha(species(i_atom))
                  endif    

              enddo 
            enddo            
         enddo            
       enddo            
      lattice_vector_SL(:,1)  = lattice_vector(:,1)*SL_i 
      lattice_vector_SL(:,2)  = lattice_vector(:,2)*SL_j 
      lattice_vector_SL(:,3)  = lattice_vector(:,3)*SL_k 
      endif 

      ! Construct cfdm hamiltonian     
      do i_atom=1,n_atoms_SL,1 !$1
         if(myid.eq.task_list_SL(i_atom)) then 
         do j_atom=i_atom,n_atoms_SL,1 !$2

         if(i_atom.eq.j_atom) then  !#1

               do i_index=1,3,1
                  do j_index=1,3,1
                  if(i_index.eq.j_index) then
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = omega_cfdm_SL(i_atom)**2.0
                  else
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = 0.d0 
                  endif  
                  enddo
               enddo

         else
            r_ij=0.0d0
            TPP=0.0d0
            dxyz=0.d0
            dxyz(:)= coords_SL(:,i_atom)-coords_SL(:,j_atom)
            r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
            Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
            sigma=(r_ij/Rvdw_12)**beta
            call MOD_CFDM_TENSOR(dxyz,r_ij,Rvdw_12,beta,TPP)
            CFDM_prefactor=omega_cfdm_SL(i_atom)*omega_cfdm_SL(j_atom)*&
                           dsqrt(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))

                   ! Transfer each dipole matrix to CFDM hamiltonian i.j accordingly         
               do i_index=1,3,1
                  do j_index=1,3,1
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = TPP(i_index,j_index)*CFDM_prefactor
                   cfdm_hamiltonian(3*j_atom-3+j_index,3*i_atom-3+i_index) = TPP(i_index,j_index)*CFDM_prefactor
                  enddo
               enddo

         endif !#1
         enddo !$2
       endif ! tasks 
      enddo !$1
      call sync_tensors(cfdm_hamiltonian,3*n_atoms_SL)

! Adds dipole field due to image cells based spherical cutoff mbd_cfdm_dip_cutoff
if (n_periodic .gt. 0) then 
      if(.NOT.mbd_scs_vacuum_axis(1)) then 
      periodic_cell_i = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,1)**2 +&
                                                                lattice_vector_SL(2,1)**2 +&
                                                                lattice_vector_SL(3,1)**2)) 
      else
      periodic_cell_i =1
      endif

      if(.NOT.mbd_scs_vacuum_axis(2)) then
      periodic_cell_j = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,2)**2 +&
                                                                lattice_vector_SL(2,2)**2 +&
                                                                lattice_vector_SL(3,2)**2))
      else
      periodic_cell_j =1
      endif     

      if(.NOT.mbd_scs_vacuum_axis(3)) then  
      periodic_cell_k = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,3)**2 +&
                                                                lattice_vector_SL(2,3)**2 +&
                                                                lattice_vector_SL(3,3)**2)) 
      else
      periodic_cell_k =1
      endif     

      do i_lattice = -periodic_cell_i, periodic_cell_i,1          !$7
         do j_lattice = -periodic_cell_j, periodic_cell_j,1         !$6  
            do k_lattice = -periodic_cell_k, periodic_cell_k,1        !$5 
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then!$4
                do i_atom=1,n_atoms_SL,1 !$3
                 if(myid.eq.task_list_SL(i_atom)) then !$ tasks
                  ! LOOP GOES OVER UPPPER TRIANGLE OF HAMILTONIAN 
                  do j_atom=i_atom,n_atoms_SL,1 !$2

                      r_ij=0.0d0
                      TPP=0.0d0
                      dxyz=0.d0
                      coord_curr(:) = coords_SL(:,i_atom) + i_lattice*lattice_vector_SL(:,1) + &
                                                            j_lattice*lattice_vector_SL(:,2) + &
                                                            k_lattice*lattice_vector_SL(:,3)

                      dxyz(:)= coords_SL(:,j_atom)- coord_curr(:)
                      r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )

                      if(r_ij.le.mbd_cfdm_dip_cutoff/bohr)  then
                        Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
                        sigma=(r_ij/Rvdw_12)**beta
                        call MOD_CFDM_TENSOR(dxyz,r_ij,Rvdw_12,beta,TPP) 
                        CFDM_prefactor=omega_cfdm_SL(i_atom)*&
                                       omega_cfdm_SL(j_atom)*&
                                       sqrt(alpha_eff_SL(i_atom)*&
                                            alpha_eff_SL(j_atom))
                        ! Transfer each dipole matrix to CFDM hamiltonian 3x3
                        ! 
                        do i_index=1,3,1
                          do j_index=1,3,1
                          !FILL UPPER BLOCK
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)=&
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          !FILL LOWER BLOCK ALSO MAKE SURE THAT YOU DONT ADD FIELD
                          !CONTRIBUTION TWO TIMES IN DIAGONAL BLOCK below is the
                          !check
                          if(i_atom.NE.j_atom) then
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)=&
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          endif
                         enddo
                        enddo
                      endif 

                  enddo !$2
                 endif !$tasks
                enddo !$3
               endif !$4
            enddo !$5
         enddo !$6
      enddo !$7
   call sync_tensors(cfdm_hamiltonian_periodic,3*n_atoms_SL)    
   cfdm_hamiltonian = cfdm_hamiltonian + cfdm_hamiltonian_periodic
endif

    errorflag=0
    LWORK=3*n_atoms_SL*(3+(3*n_atoms_SL)/2)
    call print_dsyev_info("get_mbd_enengy_at_TS",3*n_atoms_SL)
    call DSYEV('V','U',3*n_atoms_SL,cfdm_hamiltonian,3*n_atoms_SL,cfdm_eigenvalues,WORK,LWORK,errorflag)
    call check_info(errorflag,"get_mbd_enengy_at_TS","DSYEV")

    E_int=0.0d0
    E_nonint =0.0d0
    ene_mbd_ts =   0.0
    NIMAG=0 
    do i_atom =1,n_atoms_SL
      E_nonint = E_nonint + omega_cfdm_SL(i_atom)    
    enddo

    do i_atom =1,3*n_atoms_SL
      if(cfdm_eigenvalues(i_atom).ge.0.d0) then
        E_int = E_int + dsqrt(cfdm_eigenvalues(i_atom))
      else
        NIMAG= NIMAG +1  
      endif 
    enddo

    ene_mbd_ts = ((0.5*E_int)-(1.5* E_nonint))/(SL_i*SL_j*SL_k)
    if(NIMAG.gt.0) then
      write(info_str,'(A,I4,A)')"***WARNING: found ",NIMAG," negative eigenvalues in MBD_SCS energy calculation."
      call localorb_info(info_str, use_unit,'(A)',OL_norm)  
    endif
   !
   if(allocated(cfdm_hamiltonian))          deallocate(cfdm_hamiltonian)
   if(allocated(cfdm_hamiltonian_periodic)) deallocate(cfdm_hamiltonian_periodic) 
   if(allocated(cfdm_eigenvalues))          deallocate(cfdm_eigenvalues)
   if(allocated(coords_SL))                 deallocate(coords_SL)
   if(allocated(omega_cfdm_SL))             deallocate(omega_cfdm_SL)
   if(allocated(R_vdw_SL))                  deallocate(R_vdw_SL)
   if(allocated(alpha_eff_SL))              deallocate(alpha_eff_SL)  
  
   return
  endsubroutine get_mbd_enengy_at_TS

      subroutine MOD_CFDM_TENSOR(dxyz,r_ij,Rvdw12,beta,Tij)
      ! o Tij- Erf((rij/(beta*Rvdw))^4)*Grad_i X Grad_j (1/r)
      ! tensor derived from bare coulomb potential damped with
      ! Erf((rij/(beta*Rvdw))^4)

      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: Rvdw12
      real*8,intent(in) :: beta
      real*8,dimension(3,3),intent(out) :: Tij
      ! This needs declarding for the PGI Architecture
      real*8 :: erf
      ! local vars
      integer :: i_index, j_index

      
      
             do i_index=1,3,1
               do j_index=1,3,1
                Tij(i_index,j_index)=3.d0*dxyz(i_index)*dxyz(j_index)
               enddo
             enddo

             do i_index=1,3,1
               Tij(i_index,i_index) = Tij(i_index,i_index) - r_ij**2
             enddo
             Tij=Tij/r_ij**5
             Tij=Tij*erf((r_ij/(Rvdw12*beta))**4)
             Tij=-1.d0*Tij
      return
      endsubroutine MOD_CFDM_TENSOR


      subroutine MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw12,beta,Tij)
      ! o Tij- Fermi type*Grad_i X Grad_j (1/r)
      ! tensor derived from bare coulomb potential damped with
      ! Fermi type damping

      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: Rvdw12
      real*8,intent(in) :: beta
      real*8,dimension(3,3),intent(out) :: Tij
      ! This needs declarding for the PGI Architecture
      real*8 :: erf
      real*8 :: d_param,fermi_param
      ! local vars
      integer :: i_index, j_index

      d_param = 6.0
      fermi_param =r_ij/(beta*Rvdw12)

             do i_index=1,3,1
               do j_index=1,3,1
                Tij(i_index,j_index)=3.d0*dxyz(i_index)*dxyz(j_index)
               enddo
             enddo

             do i_index=1,3,1
               Tij(i_index,i_index) = Tij(i_index,i_index) - r_ij**2
             enddo
             Tij=Tij/r_ij**5
             Tij=-1.d0*Tij

             Tij=Tij*(1.0/(1.0+exp(-d_param*(fermi_param-1.0))))
      return
      endsubroutine MBD_TENSOR_MBD_rsSCS 



   subroutine self_consistent_screening_iso_damp(ene_mbd_rsSCS)
      ! local variable 
      real*8 :: C6_free
      real*8 :: alpha_free
      real*8 :: R_free
      real*8 :: mol_c6_coeff
      real*8 :: C6_atom
      real*8 :: R_vdw_free
      integer :: i_freq
      integer :: i_myatom
      real*8 :: mol_pol_eigen(3)
      real*8 :: WORK(9)
      real*8 :: ene_mbd_rsSCS
      integer::errorflag,LWORK
      ! o INPUT
      ! o hirshfeld_volume from hirshfeld partioning 
      ! o relay_matrix is 3N x 3N (imaginary) frequency dependent matrix
      ! o R_p is effective dipole spread of quantum harmonic oscilator defined
      !   with respect to Mayer representation A. Mayer, Phys. Rev. B, 75,045407(2007)
      ! o alpha_omega is single pole polarizability of atom in molecule
      ! o coupled_atom_pol contains atom resolved screened polarizabilities at every
      !   imaginary frequency 
      ! o alpha_eff is atom resolved screened polarizabilty at ZERO frequency 
      ! o C6_eff is atom resolved screened C6 coefficients calculated by
      !   integrating coupled_atom_pol over entire frequency range 



      if(.not.allocated(relay_matrix))     allocate(relay_matrix(3*n_atoms,3*n_atoms))

      if(n_periodic .gt. 0) then
      if(.not.allocated(relay_matrix_periodic)) allocate(relay_matrix_periodic(3*n_atoms,3*n_atoms))
      endif

      if(.not.allocated(R_p))              allocate(R_p(n_atoms))
      if(.not.allocated(Rvdw_iso))         allocate(Rvdw_iso(n_atoms))
      if(.not.allocated(alpha_omega))      allocate(alpha_omega(n_atoms))
      if(.not.allocated(coupled_atom_pol)) allocate(coupled_atom_pol(20,n_atoms))
      if(.not.allocated(alpha_eff))        allocate(alpha_eff(n_atoms))
      if(.not.allocated(C6_eff))           allocate(C6_eff(n_atoms))


      casmir_omega        = 0.d0
      casmir_omega_weight = 0.d0 
!     call gauleg(0.d0,30.d0,casmir_omega(1:20),casmir_omega_weight(1:20),20)
      ! currently I am using   
      call gauss_legendre_grid(casmir_omega,casmir_omega_weight)

      coupled_atom_pol = 0.d0
          mol_c6_coeff = 0.d0
             alpha_eff = 0.d0
                C6_eff = 0.d0 
 

     write(info_str,'(2x,A)')"| Many-Body Dispersion (MBD@rsSCS) energy "
     call localorb_info(info_str, use_unit,'(A)',OL_norm)
      if(n_periodic .eq. 0) then
        write(info_str,'(2x,A)')"| Dynamic molecular polarizability alpha(iw) (bohr^3)"
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
      else
        write(info_str,'(2x,A)')"| Dynamic polarizability of unit cell alpha(iw) (bohr^3)"
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
      endif
      write(info_str,'(2x,A)')"| ----------------------------------------------------------------------------------"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,A)')"|  omega(Ha)   alpha_iso      alpha_xx       alpha_yy       alpha_zz"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)


    ! Loop over Casimir-Polder frequencies
      do i_freq=0,20,1
              !! zeroout array before getting actual frequency dependent parameters
              R_p = 0.0d0
              alpha_omega = 0.0d0
              relay_matrix= 0.0d0
             mol_pol_tensor= 0.0d0
             Rvdw_iso=0.d0 
              ! loop over atoms
              do i_myatom=1,n_atoms,1
                   if (.not.vdw_hirshfeld_data_external(species(i_myatom))) then
                   call get_vdw_param(species_element(species(i_myatom)),&
                                 species_z(species(i_myatom)),C6_free,alpha_free,R_vdw_free)
                   call get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),C6_free,&
                              alpha_free,casmir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                  ! 
                   Rvdw_iso(i_myatom)= (hirshfeld_volume(i_myatom)**(1.d0/3.d0))*R_vdw_free
                   else !! Case when parameter(s) C_6, alpha and R_vdw are specified in control.in
                   call get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),vdw_hirshfeld_C6(species(i_myatom)),&
                             vdw_hirshfeld_alpha(species(i_myatom)),casmir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                   Rvdw_iso(i_myatom)= (hirshfeld_volume(i_myatom)**(1.d0/3.d0))* &
                             vdw_hirshfeld_R0(species(i_myatom))
   
                   endif 
              enddo ! end loop over atoms


              ! compute screened polarizability using isotropically damped
              ! tensor Ambrosetti et al
              ! matrix(INPUT:-alpha_omega,R_p,Rvdw_iso) 
              call get_alpha_omega_relay_matrix_iso_damp() 

              ! Contract relay_matrix^-1 for polarizability
              call contract_matrix(mol_pol_tensor)  

              ! get coupled atomic "local" polarizability by summing over row or columns of
              ! relay tensor
              if (i_freq.eq.0) then
                 ! Store static polarizability for post SCS routines  
                 call get_coupled_atom_polarizability(alpha_eff)
              endif

              if (i_freq.ge.1) then
                 call get_coupled_atom_polarizability(coupled_atom_pol(i_freq,:))
              endif 

              ! Casmir-Polder integral for molecular C6-coefficent 
              LWORK=9 
              call DSYEV('V','U',3,mol_pol_tensor,3,mol_pol_eigen,WORK,LWORK,errorflag) 
              call check_info(errorflag,"self_consistent_screening_iso_damp","DSYEV")
              if(i_freq.ge.1)  then
                  mol_c6_coeff = mol_c6_coeff + (casmir_omega_weight(i_freq)* (sum(mol_pol_eigen)/3.0)**2)
              endif

              write(info_str,'(2x,"| ",F10.6,4(e15.6))')casmir_omega(i_freq),&
                     sum(mol_pol_eigen)/3.0,mol_pol_eigen(1),mol_pol_eigen(2),mol_pol_eigen(3)
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
      enddo ! end   over freq


              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  

              if(n_periodic .eq. 0) then   
              Write(info_str,'(2x,A,f20.8,2x,A)')&
              "| Molecular C6 coefficient            :  ",0.954929658551372d0*mol_c6_coeff,"    hartree bohr^6"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              else
              Write(info_str,'(2x,A,f20.8,2x,A)')&
              "| Unit cell C6 coefficient            :  ",0.954929658551372d0*mol_c6_coeff,"    hartree bohr^6"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              endif

              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              if(n_periodic .eq. 0) then
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha (bohr^3)  of atom in molecule"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
              else
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha (bohr^3)  of atom in unit cell"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  
           
              endif  
              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)  

              do i_myatom=1,n_atoms,1
                 C6_atom = 0.0d0 
                 do i_freq=1,20,1
                 C6_atom = C6_atom + (casmir_omega_weight(i_freq)*coupled_atom_pol(i_freq,i_myatom)**2)
                 enddo
                 write(info_str,'(2x,A,I4,2x,A,f12.6,6x,f12.6)')"| ATOM ",i_myatom,&
                                      trim(species_name(species(i_myatom))),&
                                      C6_atom*0.954929658551372d0,alpha_eff(i_myatom)
                 call localorb_info(info_str, use_unit,'(A)',OL_norm)  

                 !store effctive C6 for post SCS routines 
                  C6_eff(i_myatom)=C6_atom*0.954929658551372d0
              enddo 


              write(info_str,'(2x,A)')&
              "| ----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)

              !  
              ! Compute many-body-dispersion energy new alpha_eff and C6_eff
              call get_mbd_enengy_at_iso_damp_scs(ene_mbd_rsSCS)

   

              if(allocated(R_p))                   deallocate(R_p)
              if(allocated(Rvdw_iso))              deallocate(Rvdw_iso)
              if(allocated(alpha_omega))           deallocate(alpha_omega)
              if(allocated(relay_matrix))          deallocate(relay_matrix)
              if(allocated(relay_matrix_periodic)) deallocate(relay_matrix_periodic)
              if(allocated(coupled_atom_pol))      deallocate(coupled_atom_pol)
              if(allocated(alpha_eff))             deallocate(alpha_eff)
              if(allocated(C6_eff))                deallocate(C6_eff)



      return
      endsubroutine self_consistent_screening_iso_damp

    subroutine get_alpha_omega_relay_matrix_iso_damp()

    integer ::i_index,j_index
    real*8,dimension(3,3)::TPP
    real*8,dimension(3) :: dxyz
    real*8,dimension(3) :: coord_curr
    real*8 :: r_ij
    real*8 :: r_pp
    real*8 :: Rvdw12
    real*8 :: beta
    integer :: i_row, i_col
    integer :: i_lattice, j_lattice, k_lattice
    integer :: errorflag,periodic_cell_i,periodic_cell_j,periodic_cell_k

    !For LAPACK
    integer,dimension(3*n_atoms):: IPIV
    real*8,dimension(3*n_atoms):: WORK

    ! initio values
    relay_matrix=0.0d0

        select case (flag_xc)
          case (1) !PBE0
               beta=0.85
          case (6) ! PBE
               beta=0.83
          case (7) !HSE 
               beta=0.85
          case default
               beta=1.0
               write(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"WITHOUT short range damping"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
        end select



     ! compute relay matrix of  cluster or unit cell
       do i_row=1,n_atoms,1 !#1
         if(myid.eq.task_list(i_row)) then      
         do i_col=i_row,n_atoms,1 !#2
         TPP=0.d0
            if(i_row.eq.i_col) then  !$1
               do i_index=1,3,1
                  do j_index=1,3,1
                     if(i_index.eq.j_index) then
                        relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=1.d0/alpha_omega(i_row)
                     else 
                        relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=0.d0
                     endif   
                  enddo
               enddo

            else
               dxyz(:) = coords(:,i_col)-coords(:,i_row)
               r_ij = dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
               r_pp = dSqrt(R_p(i_row)**2 + R_p(i_col)**2)
               Rvdw12 = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
               call SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
               do i_index=1,3,1
                  do j_index=1,3,1
                   relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=TPP(i_index,j_index)
                   relay_matrix(3*i_col-3+j_index,3*i_row-3+i_index)=TPP(i_index,j_index)
                  enddo
               enddo

            endif !$1
         enddo   !#2
        endif ! task 
       enddo  !#1
  call sync_tensors(relay_matrix,3*n_atoms)
  if (n_periodic .gt. 0) then
   
      relay_matrix_periodic=0.d0 

      if(.NOT.mbd_scs_vacuum_axis(1)) then     
      periodic_cell_i = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 +& 
                                                               lattice_vector(2,1)**2 +&
                                                               lattice_vector(3,1)**2)) 
      else
      periodic_cell_i = 1
      endif 
      if(.NOT.mbd_scs_vacuum_axis(2)) then     
      periodic_cell_j = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 +& 
                                                               lattice_vector(2,2)**2 +&
                                                               lattice_vector(3,2)**2)) 
      else
      periodic_cell_j = 1
      endif 
      if(.NOT.mbd_scs_vacuum_axis(3)) then     
      periodic_cell_k = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 +& 
                                                               lattice_vector(2,3)**2 +&
                                                               lattice_vector(3,3)**2))
      else
      periodic_cell_k = 1
      endif 
            
      do i_lattice = -periodic_cell_i, periodic_cell_i,1
         do j_lattice = -periodic_cell_j, periodic_cell_j,1
            do k_lattice = -periodic_cell_k, periodic_cell_k,1
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then !#1
                  do i_row = 1, n_atoms, 1 ! atom1 loop
                    if(myid.eq.task_list(i_row)) then
                     do i_col = i_row, n_atoms, 1 ! atom2 loop
                          coord_curr = 0.d0
                          dxyz = 0.d0
                                r_ij = 0.d0
                                TPP  = 0.d0
                               Rvdw12= 0.d0  
                          ! find the coordinate of images
                          coord_curr(:) = coords(:,i_col) + i_lattice*lattice_vector(:,1) + &
                                                            j_lattice*lattice_vector(:,2) + &
                                                            k_lattice*lattice_vector(:,3)

                          dxyz(:) = coords(:,i_row)- coord_curr(:)
                          r_ij    = sqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                          r_pp    = sqrt(R_p(i_row)**2 + R_p(i_col)**2)
                          Rvdw12  = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
                          if(r_ij.le.mbd_scs_dip_cutoff/bohr) then
                            call SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
                            do i_index=1,3,1
                            do j_index=1,3,1
                            relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index)=&
                            relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index) + TPP(i_index,j_index)
                             if(i_col.NE.i_row) then
                               relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index)=&
                               relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index) + TPP(i_index,j_index)
                              endif
                            enddo
                            enddo
                          endif
                     enddo !atom2 loop
                   endif! parallel
                  enddo !atom1 loop
               endif  !#1
            enddo
         enddo
      enddo
   call sync_tensors(relay_matrix_periodic,3*n_atoms)
   relay_matrix = relay_matrix + relay_matrix_periodic  
   endif

   call DGETRF(3*n_atoms, 3*n_atoms, relay_matrix, 3*n_atoms, IPIV, errorflag)
   if(errorflag.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call localorb_info(info_str, use_unit,'(A)',OL_norm)
   endif
   call check_info(errorflag,"get_alpha_omega_relay_matrix_iso_damp","DGETRF")

   call DGETRI(3*n_atoms, relay_matrix, 3*n_atoms, IPIV, WORK,3*n_atoms,errorflag )
   if(errorflag.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call localorb_info(info_str, use_unit,'(A)',OL_norm)
   endif
   call check_info(errorflag,"get_alpha_omega_relay_matrix_iso_damp","DGETRI")
   return
   endsubroutine  get_alpha_omega_relay_matrix_iso_damp


  subroutine get_mbd_enengy_at_iso_damp_scs(ene_mbd_rsSCS)

     ! local vars  
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian 
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian_periodic
     real*8,dimension(:,:),allocatable:: coords_SL
     real*8,dimension(:),allocatable:: cfdm_eigenvalues
     real*8,dimension(:),allocatable:: omega_cfdm_SL
     real*8,dimension(:),allocatable:: R_vdw_SL
     real*8,dimension(:),allocatable::alpha_eff_SL
     real*8,dimension(:),allocatable:: WORK
     real*8,dimension(:),allocatable:: task_list_SL  
     real*8,dimension(3,3):: TPP,lattice_vector_SL
     real*8,dimension(3)::dxyz,coord_curr
     real*8:: r_ij
     real*8:: C6_free
     real*8:: alpha_free
     real*8:: R_vdw_free
     real*8:: Rvdw_12
     real*8:: beta
     real*8:: CFDM_prefactor
     real*8:: E_int
     real*8:: E_nonint
     real*8:: sigma
     real*8:: ene_mbd_rsSCS
     integer :: errorflag,i_atom,j_atom
     integer :: SL_i,SL_j,SL_k ! MBD super cell index's
     integer :: i_index, j_index ! 3x3 block index
     integer :: i_lattice,j_lattice,k_lattice ! lattice increament index
     integer :: periodic_cell_i,periodic_cell_j,periodic_cell_k
     integer :: NIMAG
     integer ::n_atoms_SL
     integer :: LWORK
!     Begin work       
        select case (flag_xc)
          case (1) !PBE0
               beta=0.85
          case (6) ! PBE
               beta=0.83
          case (7) !HSE 
               beta=0.85
          case default
               beta=1.0
               write(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"WITHOUT short range damping"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
        end select

      if ( n_periodic .eq. 0) then
      ! For Cluster calculation NO super cell needed 
      ! intiate this index so as to avoide any invalid floting point when
      ! normalizing energy in case of cluster/molecule
      SL_i = 1 
      SL_j = 1
      SL_k = 1

!     number of atoms in cluster /molecule
      n_atoms_SL = n_atoms
 
      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK((3*n_atoms_SL)*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))

      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian=0.0d0 
      cfdm_eigenvalues=0.0d0

!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN CLUSTER/MOLECULE
      do i_atom=1,n_atoms_SL
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks)  
      enddo

          do i_atom =1,n_atoms
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0
                  C6_free=0.d0
                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                C6_free,alpha_free,R_vdw_free)

                  R_vdw_SL(i_atom)= R_vdw_free*((alpha_eff(i_atom)/alpha_free)**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  coords_SL(:,i_atom)   = coords(:,i_atom)
                  alpha_eff_SL(i_atom) = alpha_eff(i_atom)
                  else
                  !! In case when parameter(s) C_6, alpha, R_vdw are specified via hirshfeld_param keyword in control.in      
                  R_vdw_SL(i_atom)= vdw_hirshfeld_R0(species(i_atom))* & 
                      ((alpha_eff(i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  coords_SL(:,i_atom)   = coords(:,i_atom)
                  alpha_eff_SL(i_atom) = alpha_eff(i_atom)

                  endif 
         enddo
      else

      if(.NOT.mbd_scs_vacuum_axis(1)) then
      SL_i = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2)) 
      else
      SL_i = 1 
      endif

      if(.NOT.mbd_scs_vacuum_axis(2)) then
      SL_j = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)) 
      else
      SL_j = 1 
      endif

      if(.NOT.mbd_scs_vacuum_axis(3)) then
      SL_k = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)) 
      else
      SL_k = 1 
      endif
       
!     number of atoms in super cell 
      n_atoms_SL = n_atoms*SL_i*SL_j*SL_k 
      write(info_str,'(2X,A,i3,A,i3,A,i3,A,f6.2,A)') &
         "| Creating super cell of dimension",  SL_i," X ", SL_j," X ", &
         SL_k, " in MBD calculation using" ,mbd_supercell_cutoff, &
         "  Angstrom radius"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,"| containing", I6,A)')n_atoms_SL, "  atom"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
        


      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_hamiltonian_periodic)) allocate(cfdm_hamiltonian_periodic(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK(3*n_atoms_SL*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))
      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian = 0.0d0 
      cfdm_eigenvalues = 0.0d0
      cfdm_hamiltonian_periodic = 0.d0 
!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN SUPERCELL/MOLECULE
      do i_atom=1,n_atoms_SL 
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks) 
      enddo

      ! Get all the parameter  required for MBD with super cell calculation
       j_atom = 0
       do i_lattice = 0, SL_i-1,1
         do j_lattice = 0, SL_j-1,1
            do k_lattice = 0, SL_k-1,1
              do i_atom =1,n_atoms !<--- atom index
                  j_atom = j_atom +1  !<--- dummy index supercell
                  coords_SL(:,j_atom) = coords(:,i_atom) + i_lattice*lattice_vector(:,1) + &
                                                           j_lattice*lattice_vector(:,2) + &
                                                           k_lattice*lattice_vector(:,3)   
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0  
                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then  
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                C6_free,alpha_free,R_vdw_free)

                  R_vdw_SL(j_atom)=R_vdw_free*((alpha_eff(i_atom)/alpha_free)**(1.d0/3.d0))
                  alpha_eff_SL(j_atom) = alpha_eff(i_atom)
                  omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  else
                   !! In case when parameter(s) C_6, alpha, R_vdw are specified via hirshfeld_param keyword in control.in    
                  R_vdw_SL(j_atom)=vdw_hirshfeld_R0(species(i_atom))* & 
                     ((alpha_eff(i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  alpha_eff_SL(j_atom) = alpha_eff(i_atom)
                  omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
  
                  endif  
              enddo 
            enddo            
         enddo            
       enddo            
      lattice_vector_SL(:,1)  = lattice_vector(:,1)*SL_i 
      lattice_vector_SL(:,2)  = lattice_vector(:,2)*SL_j 
      lattice_vector_SL(:,3)  = lattice_vector(:,3)*SL_k 
      endif 

      ! Construct cfdm hamiltonian     
      do i_atom=1,n_atoms_SL,1 !$1
         if(myid.eq.task_list_SL(i_atom)) then
         do j_atom=i_atom,n_atoms_SL,1 !$2

         if(i_atom.eq.j_atom) then  !#1

               do i_index=1,3,1
                  do j_index=1,3,1
                  if(i_index.eq.j_index) then
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = omega_cfdm_SL(i_atom)**2.0
                  endif  
                  enddo
               enddo

         else
            r_ij=0.0d0
            TPP=0.0d0
            dxyz=0.d0
            dxyz(:)= coords_SL(:,i_atom)-coords_SL(:,j_atom)
            r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
            Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
            sigma=(r_ij/Rvdw_12)**beta
            call MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw_12,beta,TPP)
            CFDM_prefactor=omega_cfdm_SL(i_atom)*omega_cfdm_SL(j_atom)*&
                           dsqrt(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))

                   ! Transfer each dipole matrix to CFDM hamiltonian i.j accordingly         
               do i_index=1,3,1
                  do j_index=1,3,1
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = TPP(i_index,j_index)*CFDM_prefactor
                   cfdm_hamiltonian(3*j_atom-3+j_index,3*i_atom-3+i_index) = TPP(i_index,j_index)*CFDM_prefactor
                  enddo
               enddo

         endif !#1
         enddo !$2
       endif ! tasks 
      enddo !$1
      call sync_tensors(cfdm_hamiltonian,3*n_atoms_SL)

! Adds dipole field due to image cells based spherical cutoff mbd_cfdm_dip_cutoff
if (n_periodic .gt. 0) then 
      if(.NOT.mbd_scs_vacuum_axis(1)) then
      periodic_cell_i = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,1)**2 +&
                                                                lattice_vector_SL(2,1)**2 +&
                                                                lattice_vector_SL(3,1)**2))
      else
      periodic_cell_i = 1
      endif   

      if(.NOT.mbd_scs_vacuum_axis(2)) then
      periodic_cell_j = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,2)**2 +&
                                                                lattice_vector_SL(2,2)**2 +&
                                                                lattice_vector_SL(3,2)**2)) 
      else
      periodic_cell_j = 1 
      endif   

      if(.NOT.mbd_scs_vacuum_axis(3)) then
      periodic_cell_k = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,3)**2 +&
                                                                lattice_vector_SL(2,3)**2 +&
                                                                lattice_vector_SL(3,3)**2)) 
      else
      periodic_cell_k = 1 
      endif  
 
      do i_lattice = -periodic_cell_i, periodic_cell_i,1          !$7
         do j_lattice = -periodic_cell_j, periodic_cell_j,1         !$6  
            do k_lattice = -periodic_cell_k, periodic_cell_k,1        !$5 
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then!$4
                do i_atom=1,n_atoms_SL,1 !$3
                 if(myid.eq.task_list_SL(i_atom)) then !$ tasks
                  ! LOOP GOES OVER UPPPER TRIANGLE OF HAMILTONIAN 
                  do j_atom=i_atom,n_atoms_SL,1 !$2

                      r_ij=0.0d0
                      TPP=0.0d0
                      dxyz=0.d0
                      coord_curr(:) = coords_SL(:,i_atom) + i_lattice*lattice_vector_SL(:,1) + &
                                                            j_lattice*lattice_vector_SL(:,2) + &
                                                            k_lattice*lattice_vector_SL(:,3)

                      dxyz(:)= coords_SL(:,j_atom)- coord_curr(:)
                      r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )

                      if(r_ij.le.mbd_cfdm_dip_cutoff/bohr)  then
                        Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
                        sigma=(r_ij/Rvdw_12)**beta
                        call MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw_12,beta,TPP) 
                        CFDM_prefactor=omega_cfdm_SL(i_atom)*&
                                       omega_cfdm_SL(j_atom)*&
                                       sqrt(alpha_eff_SL(i_atom)*&
                                            alpha_eff_SL(j_atom))
                        do i_index=1,3,1
                          do j_index=1,3,1
                          !FILL UPPER BLOCK
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)=&
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          !FILL LOWER BLOCK ALSO MAKE SURE THAT YOU DONT ADD FIELD
                          !CONTRIBUTION TWO TIMES IN DIAGONAL BLOCK below is the
                          !check
                          if(i_atom.NE.j_atom) then
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)=&
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          endif
                         enddo
                        enddo
                      endif 

                  enddo !$2
                 endif !$tasks
                enddo !$3
               endif !$4
            enddo !$5
         enddo !$6
      enddo !$7
   call sync_tensors(cfdm_hamiltonian_periodic,3*n_atoms_SL)    
   cfdm_hamiltonian = cfdm_hamiltonian + cfdm_hamiltonian_periodic
endif

    errorflag=0
    LWORK=3*n_atoms_SL*(3+(3*n_atoms_SL)/2)
    call print_dsyev_info("get_mbd_enengy_at_iso_damp_scs",3*n_atoms_SL)
    call DSYEV('V','U',3*n_atoms_SL,cfdm_hamiltonian,3*n_atoms_SL,cfdm_eigenvalues,WORK,LWORK,errorflag)
    call check_info(errorflag,"get_mbd_enengy_at_iso_damp_scs","DSYEV")

    E_int=0.0d0
    E_nonint =0.0d0
    ene_mbd_rsSCS =   0.0
    NIMAG=0 
    do i_atom =1,n_atoms_SL
      E_nonint = E_nonint + omega_cfdm_SL(i_atom)    
    enddo

    do i_atom =1,3*n_atoms_SL
      if(cfdm_eigenvalues(i_atom).ge.0.d0) then
        E_int = E_int + dsqrt(cfdm_eigenvalues(i_atom))
      else
        NIMAG= NIMAG +1  
      endif 
    enddo

    ene_mbd_rsSCS = ((0.5*E_int)-(1.5* E_nonint))/(SL_i*SL_j*SL_k)
    if(NIMAG.gt.0) then
      write(info_str,'(A,I4,A)')"***WARNING: found ",NIMAG," negative eigenvalues in MBD_rsSCS energy calculation."
      call localorb_info(info_str, use_unit,'(A)',OL_norm)  
    endif
   !
   if(allocated(cfdm_hamiltonian))          deallocate(cfdm_hamiltonian)
   if(allocated(cfdm_hamiltonian_periodic)) deallocate(cfdm_hamiltonian_periodic) 
   if(allocated(cfdm_eigenvalues))          deallocate(cfdm_eigenvalues)
   if(allocated(coords_SL))                 deallocate(coords_SL)
   if(allocated(omega_cfdm_SL))             deallocate(omega_cfdm_SL)
   if(allocated(R_vdw_SL))                  deallocate(R_vdw_SL)
   if(allocated(alpha_eff_SL))              deallocate(alpha_eff_SL)  
  
   return

  endsubroutine get_mbd_enengy_at_iso_damp_scs

  subroutine check_info(info,caller,lapackcall)

   implicit none

   !Arguments

   integer,    INTENT(IN)   :: info
   character(len=*), INTENT(IN)   :: caller
   character(len=*), INTENT(IN)   :: lapackcall

   if (info .ne. 0) then
     write(info_str,'(2X,A,A)')"***Error in LAPACK CALL ", trim(lapackcall)
     call aims_stop (info_str,caller)
   endif

  endsubroutine check_info

  subroutine print_dsyev_info(caller, rank)
   implicit none

   !Arguments

   integer,    INTENT(IN)   :: rank
   character(len=*), INTENT(IN)   :: caller

   write(info_str,'(2X,A,A,A,I10,A,I10,A)') "| ",trim(caller), &
     " solves ", rank, " x ", rank, " eigenvalue problem."
   call localorb_info(info_str, use_unit,'(A)',OL_norm)

  endsubroutine print_dsyev_info


endmodule SCS_CFDM
