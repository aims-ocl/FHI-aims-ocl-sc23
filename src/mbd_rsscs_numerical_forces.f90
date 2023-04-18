      module mbd_rsSCS_numerical_forces 
      use physics, only: forces_on
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
      use elpa1_2013
      use elpa2_2013
      use scs_cfdm , only:gauss_legendre_grid , SCS_TENSOR_MBD_rsSCS ,MBD_TENSOR_MBD_rsSCS ,get_alpha_omega_and_Rp
      implicit none
      PRIVATE ! DEFAULT PRIVATE
      PUBLIC :: mbd_rsSCS_energy_and_forces                !<-- compute finite difference forces for MBD@rsSCS

      character*200 :: info_str  
      integer :: info
      real*8,dimension(:,:),allocatable:: coupled_atom_pol                 
      real*8,dimension(:,:,:),allocatable:: coupled_molecule_pol
      real*8,dimension(:,:),allocatable:: relay_matrix
      real*8,dimension(:,:), allocatable :: relay_matrix_local
      real*8,dimension(:),allocatable:: alpha_omega
      real*8,dimension(:),allocatable:: R_p
      real*8,dimension(:),allocatable:: Rvdw_iso
      real*8,dimension(:,:,:),allocatable::fd_conf_cords
      real*8,dimension(:),allocatable::fd_conf_ene_mbd
      real*8 :: fd_step_size
      
      real*8 :: casimir_omega(0:20)
      real*8 :: casimir_omega_weight(0:20)
      ! scalapack stuff        
      integer,EXTERNAL :: numroc   ! blacs routine
      integer :: me, procs, icontxt, prow, pcol, max_procs, myrow, mycol,nblocks  ! blacs data
      integer :: myArows, myAcols   ! size of local subset of global array
      integer :: mbd_comm_rows, mbd_comm_cols ! ELPA
      integer, dimension(9)   :: ides_a, ides_mbd, ides_mbd_vec  ! scalapack array desc
      integer, dimension(9)   :: ALPHA_WORK
      integer :: iproc,jproc, myi,myj, local_tensor_i,local_tensor_j
      real*8,dimension(3,3)::inv_alpha_tensor
      integer,dimension(:),allocatable:: local_ipvi
      real*8,dimension(:),allocatable:: local_work
      integer,dimension(:),allocatable:: local_iwork
      real*8,dimension(:,:),allocatable:: cfdm_ham_local
      real*8,dimension(:),allocatable:: cfdm_ham_local_eigen
      real*8,dimension(:,:),allocatable:: cfdm_ham_local_vectors 
      logical :: use_scalapack_inversion
      
 

      contains



      subroutine allocate_mpi_scs_data()
  
      if(.not.allocated(relay_matrix)) then
        allocate(relay_matrix(3*n_atoms,3*n_atoms), stat=info)
        call check_allocation(info, 'relay_matrix: memory allocation request denied')
        relay_matrix=0.d0
      endif

      if(.not.allocated(R_p))  then
        allocate(R_p(n_atoms), stat=info)
        call check_allocation(info, 'R_p: memory allocation request denied')
      endif
      if(.not.allocated(Rvdw_iso)) then
        allocate(Rvdw_iso(n_atoms), stat=info)
        call check_allocation(info, 'Rvdw_iso: memory allocation request denied')
      endif
      if(.not.allocated(alpha_omega)) then
        allocate(alpha_omega(n_atoms), stat=info)
        call check_allocation(info, 'alpha_omega: memory allocation request denied')
      endif
      if(.not.allocated(coupled_atom_pol)) then
        allocate(coupled_atom_pol(0:20,n_atoms), stat=info)
        call check_allocation(info, 'coupled_atom_pol: memory allocation request denied')
        coupled_atom_pol=0.d0
      endif
      if(.not.allocated(coupled_molecule_pol)) then
        allocate(coupled_molecule_pol(3,3,0:20), stat=info)
        call check_allocation(info, 'coupled_molecule_pol: memory allocation request denied')
        coupled_molecule_pol=0.d0
      endif


endsubroutine allocate_mpi_scs_data


subroutine deallocate_mpi_scs_data

      if(allocated(relay_matrix)) then
        deallocate(relay_matrix, stat=info)
        call check_allocation(info, 'relay_matrix: memory deallocation request denied')
      endif

      if(allocated(R_p))  then
        deallocate(R_p, stat=info)
        call check_allocation(info, 'R_p: memory deallocation request denied')
      endif
      if(allocated(Rvdw_iso)) then
        deallocate(Rvdw_iso, stat=info)
        call check_allocation(info, 'Rvdw_iso: memory deallocation request denied')
      endif
      if(allocated(alpha_omega)) then
        deallocate(alpha_omega, stat=info)
        call check_allocation(info, 'alpha_omega: memory deallocation request denied')
      endif
      if(allocated(coupled_atom_pol)) then
        deallocate(coupled_atom_pol, stat=info)
        call check_allocation(info, 'coupled_atom_pol: memory deallocation request denied')      
      endif
      if(allocated(coupled_molecule_pol)) then
        deallocate(coupled_molecule_pol, stat=info)
        call check_allocation(info, 'coupled_molecule_pol: memory deallocation request denied')       
      endif
      if(allocated(fd_conf_ene_mbd)) then
        deallocate(fd_conf_ene_mbd, stat=info)
        call check_allocation(info, 'fd_conf_ene_mbd: memory deallocation request denied')       
      endif

endsubroutine deallocate_mpi_scs_data



subroutine allocate_scalapack_scs_data()
  implicit none
    ! allocate local block  
    ! assign associated descriptor for the local relay
    ! Tensor 3x3--->block(i,j) for I, J atoms
    nblocks = 3
    myArows = numroc(3*n_atoms,nblocks,myrow,0,prow)
    myAcols = numroc(3*n_atoms,nblocks,mycol,0,pcol)
  
    Write(info_str,'(2X,A,I6,A,I6,A,I6,A)')"| Maping polarizability tensor on",prow*pcol,&
    " blacs processes each containing",myArows/3," X ",myAcols/3," tensor blocks(atoms)"
    call localorb_info(info_str,use_unit,'(A)')
    Write(info_str,'(2X,A,F12.6,A)')"| Each processes allocates ",(8*myArows*myAcols)/2.0**20," MB memory"    
    call localorb_info(info_str,use_unit,'(A)')
        
    
    if(.not.allocated(relay_matrix_local))  then
       allocate(relay_matrix_local(myArows,myAcols),stat=info)
       call check_allocation(info, 'relay_matrix_local: memory allocation request denied')
    endif


    !Descriptor for relay_matrix_local tensor 
    call descinit(ides_a, 3*n_atoms, 3*n_atoms, nblocks, nblocks, 0, 0, &
                    icontxt, myArows, info)

    if(.not.allocated(local_ipvi))  then
       allocate(local_ipvi(3*n_atoms),stat=info)
       call check_allocation(info, 'local_ipvi: memory allocation request denied')
    endif

    if(.not.allocated(local_work))  then
       allocate(local_work(myArows*myAcols),stat=info)
       call check_allocation(info, 'local_work: memory allocation request denied')
    endif
    if(.not.allocated(local_iwork))  then
       allocate(local_iwork(myArows*myAcols),stat=info)
       call check_allocation(info, 'local_iwork: memory allocation request denied')
    endif 
    if(.not.allocated(R_p))  then
      allocate(R_p(n_atoms), stat=info)
      call check_allocation(info, 'R_p: memory allocation request denied')
    endif
    if(.not.allocated(Rvdw_iso)) then
      allocate(Rvdw_iso(n_atoms), stat=info)
      call check_allocation(info, 'Rvdw_iso: memory allocation request denied')
    endif
    if(.not.allocated(alpha_omega)) then
      allocate(alpha_omega(n_atoms), stat=info)
      call check_allocation(info, 'alpha_omega: memory allocation request denied')
    endif
    if(.not.allocated(coupled_atom_pol)) then
      allocate(coupled_atom_pol(0:20,n_atoms), stat=info)
      call check_allocation(info, 'coupled_atom_pol: memory allocation request denied')
      coupled_atom_pol=0.d0
    endif
    if(.not.allocated(coupled_molecule_pol)) then
      allocate(coupled_molecule_pol(3,3,0:20), stat=info)
      call check_allocation(info, 'coupled_molecule_pol: memory allocation request denied')
      coupled_molecule_pol=0.d0
    endif

endsubroutine allocate_scalapack_scs_data


subroutine deallocate_scalapack_scs_data()
    
    
    if(allocated(relay_matrix_local))  then
       deallocate(relay_matrix_local,stat=info)
       call check_allocation(info, 'relay_matrix_local: memory deallocation request denied')
    endif

    if(allocated(local_ipvi))  then
       deallocate(local_ipvi,stat=info)
       call check_allocation(info, 'local_ipvi: memory deallocation request denied')
    endif

    if(allocated(local_work))  then
       deallocate(local_work,stat=info)
       call check_allocation(info, 'local_work: memory deallocation request denied')
    endif
    if(allocated(local_iwork))  then
       deallocate(local_iwork,stat=info)
       call check_allocation(info, 'local_iwork: memory deallocation request denied')
    endif 
    if(allocated(R_p))  then
      deallocate(R_p, stat=info)
      call check_allocation(info, 'R_p: memory deallocation request denied')
    endif
    if(allocated(Rvdw_iso)) then
      deallocate(Rvdw_iso, stat=info)
      call check_allocation(info, 'Rvdw_iso: memory deallocation request denied')
    endif
    if(allocated(alpha_omega)) then
      deallocate(alpha_omega, stat=info)
      call check_allocation(info, 'alpha_omega: memory deallocation request denied')
    endif
    if(allocated(coupled_atom_pol)) then
      deallocate(coupled_atom_pol, stat=info)
      call check_allocation(info, 'coupled_atom_pol: memory deallocation request denied')
      
    endif
    if(allocated(coupled_molecule_pol)) then
      deallocate(coupled_molecule_pol, stat=info)
      call check_allocation(info, 'coupled_molecule_pol: memory deallocation request denied')
    endif

endsubroutine deallocate_scalapack_scs_data


subroutine setup_blacs()
  implicit none
    call blacs_pinfo (me,procs)
    !Blockcyclic 2D  
    do pcol = nint(sqrt(real(procs))),2,-1
       if(mod(procs,pcol) == 0 ) exit
    enddo
    prow = procs/pcol
    call blacs_get(0,0,icontxt)
    call blacs_gridinit(icontxt,'R',prow,pcol)
    call blacs_gridinfo(icontxt,prow,pcol,myrow,mycol)
    if (use_elpa_mbd) then
      call get_elpa_row_col_comms_2013(mpi_comm_world, myrow, mycol, &
           mbd_comm_rows, mbd_comm_cols)
    endif  

   ! Need information for inversion
   ! nblocks must be 3 (hardwired in polarizability code)
   max_procs = max(prow,pcol)
   nblocks = 3
   if (nblocks * max_procs > 3*n_atoms) then
     use_scalapack_inversion = .false.
   else
     use_scalapack_inversion = .true.
   endif
    

endsubroutine setup_blacs

subroutine exit_blacs()
  implicit none

  if (use_elpa_mbd) then
     call mpi_comm_free( mbd_comm_rows,info )
     call mpi_comm_free( mbd_comm_cols,info )
  endif  
  call BLACS_Gridexit ( icontxt )

endsubroutine exit_blacs



subroutine map_global_index_to_local_index(i_tensor,proc_i,i_nblocks,proc_local,tensor_local)
   implicit none
   integer :: i_tensor, proc_i, i_nblocks, proc_local, tensor_local, Iminus
     
   Iminus = (3*i_tensor-2)-1                      
   proc_local   = mod((Iminus/i_nblocks),proc_i)
   tensor_local  = (Iminus/(proc_i*i_nblocks))*i_nblocks + mod(Iminus,i_nblocks) + 1
   return
endsubroutine



  subroutine allocate_mbd_data(number_of_oscillators)
    implicit none
    integer,intent(in) :: number_of_oscillators

    nblocks=3
    myArows = numroc(3*number_of_oscillators,nblocks,myrow,0,prow)
    myAcols = numroc(3*number_of_oscillators,nblocks,mycol,0,pcol)
    
    call descinit(ides_mbd,3*number_of_oscillators,3*number_of_oscillators,nblocks,nblocks,0,0,icontxt,myArows,info)  
    call descinit(ides_mbd_vec,3*number_of_oscillators,3*number_of_oscillators,nblocks,nblocks,0,0,icontxt,myArows,info)       
 
    if(.not.allocated(cfdm_ham_local_vectors)) then
      allocate(cfdm_ham_local_vectors(myArows,myAcols) , stat=info)
      call check_allocation(info, 'cfdm_ham_local_vectors')
    endif
    if(.not.allocated(cfdm_ham_local_eigen)) then
      allocate(cfdm_ham_local_eigen(3*number_of_oscillators), stat=info)
      call check_allocation(info, 'cfdm_ham_local_eigen')
    endif
        
    if(.not.allocated(cfdm_ham_local)) then
      allocate(cfdm_ham_local(myArows,myAcols), stat=info)
      call check_allocation(info, 'cfdm_ham_local')
    endif
    
endsubroutine allocate_mbd_data



subroutine deallocate_mbd_data()
  implicit none
    if(allocated(cfdm_ham_local)) then
      deallocate(cfdm_ham_local, stat=info)
      call check_allocation(info, 'cfdm_ham_local: Deallocation request denied')
    endif
     if(allocated(cfdm_ham_local_vectors)) then
      deallocate(cfdm_ham_local_vectors, stat=info)
      call check_allocation(info, 'cfdm_ham_local_vectors: Deallocation request denied')
    endif
    if(allocated(cfdm_ham_local_eigen)) then
      deallocate(cfdm_ham_local_eigen, stat=info)
      call check_allocation(info, 'cfdm_ham_local_eigen: Deallocation request denied')
    endif 


end subroutine deallocate_mbd_data

real*8 function C6AA(alpha)
    implicit none
    real*8,dimension(0:20):: alpha
    integer :: i_frequency
    real*8  :: C6_atom
    C6_atom=0.d0
     do i_frequency=1,20,1
          C6_atom = C6_atom + (casimir_omega_weight(i_frequency)*alpha(i_frequency)*alpha(i_frequency))
     enddo
     C6AA=C6_atom*3.0/pi
endfunction


subroutine get_identity(matrix)
   implicit none
   real*8,dimension(3,3)::matrix
   integer :: i
   matrix=0.0d0

   do i=1,3,1
   matrix(i,i)=1.d0   
   enddo   
   return
endsubroutine get_identity





subroutine mbd_rsSCS_energy_and_forces(ene_mbd_rsSCS,mbd_rsSCS_forces)
    implicit none
      ! local variable 
      real*8,dimension(3,n_atoms)::mbd_rsSCS_forces
      real*8                     :: ene_mbd_rsSCS
      real*8                     :: C6_free
      real*8                     :: alpha_free
      real*8                     :: R_vdw_free
      integer                    :: i_freq
      integer                    :: i_myatom
      integer                    :: i_config
      integer                    :: i_config_counter 
      integer                    :: i_atom
      integer                    :: i_coord
      ! o INPUT
      ! o hirshfeld_volume from hirshfeld partioning 
      ! o relay_matrix is 3N x 3N (imaginary) frequency dependent matrix
      ! o R_p is effective dipole spread of quantum harmonic oscilator defined
      !       with respect to Mayer representation A. Mayer, Phys. Rev. B, 75,045407(2007)
      ! o alpha_omega is single pole polarizability of atom in molecule
      ! o coupled_atom_pol contains atom resolved screened polarizabilities at every
      !   imaginary frequency 
      
      ! Begin Work        

       call initialize_finite_diff_infrastructure()      
      
       if(use_scalapack_mbd) then
          call  setup_blacs()
       endif
          
       if(use_scalapack_mbd .and. use_scalapack_inversion) then
        call  allocate_scalapack_scs_data()
       else 
        call  allocate_mpi_scs_data()
       endif

!      write(info_str,'(A)')"  Computing MBD@rsSCS forces using finite differences..."
!      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      
  

      call gauss_legendre_grid(casimir_omega,casimir_omega_weight)
  
      write(info_str,'(2x,A)')"| Many-Body Dispersion (MBD@rsSCS) energy "
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      if(n_periodic .eq. 0) then
        write(info_str,'(2x,A)')"| Dynamic molecular polarizability alpha(iw)(bohr^3)"
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
      else
        write(info_str,'(2x,A)')"| Dynamic polarizability of unit cell alpha(iw)(bohr^3)"
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
      endif
      write(info_str,'(2x,A)')"|----------------------------------------------------------------------------------"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,A)')"|  omega(Ha)        XX             YY            ZZ              YZ            XZ            XY"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      
      
      i_config_counter=0 

      do i_config=0,6*n_atoms
          ! Loop over Casimir-Polder frequencies
         
          ! We only want to evaluate forces for unconstrained atoms/directions
          if (use_relaxation_constraints .and. i_config > 0) then
             i_atom = ceiling(i_config/6.d0)
             i_coord = ceiling((mod(i_config-1,6) + 1)/2d0)
             
             if (constrain_coords(i_coord,i_atom)) then
                fd_conf_ene_mbd(i_config) = 0.d0
                
                if(forces_on) then 
                  if (i_config.gt.0.AND.MOD(i_config,6).eq.0) then
                     i_config_counter = i_config_counter + 1 
                    write(info_str,'(2x,A,I6,3(2x,E13.6))')"| MBD@rsSCS force on ATOM",i_config_counter,&
                      -1.d0*(fd_conf_ene_mbd(i_config-5)-fd_conf_ene_mbd(i_config-4))/(2.0d0*fd_step_size)*(hartree / bohr),&
                      -1.d0*(fd_conf_ene_mbd(i_config-3)-fd_conf_ene_mbd(i_config-2))/(2.0d0*fd_step_size)*(hartree / bohr),&
                      -1.d0*(fd_conf_ene_mbd(i_config-1)-fd_conf_ene_mbd(i_config  ))/(2.0d0*fd_step_size)*(hartree / bohr)
                    call localorb_info(info_str, use_unit,'(A)',OL_norm)
                  endif
                endif

                cycle 
             
             endif
             
          endif
           
          do i_freq=0,20,1
              R_p = 0.0d0
              alpha_omega = 0.0d0
              Rvdw_iso=0.d0 
              ! loop over atoms
              do i_myatom=1,n_atoms,1
                if (.not.vdw_hirshfeld_data_external(species(i_myatom))) then
                  call get_vdw_param(species_element(species(i_myatom)),&
                                 species_z(species(i_myatom)),C6_free,alpha_free,R_vdw_free)
                  call get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),C6_free,&
                              alpha_free,casimir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                  Rvdw_iso(i_myatom)= (hirshfeld_volume(i_myatom)**(1d0/3d0))*R_vdw_free
                else
                  call get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),vdw_hirshfeld_C6(species(i_myatom)),&
                              vdw_hirshfeld_alpha(species(i_myatom)),casimir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                  Rvdw_iso(i_myatom)= (hirshfeld_volume(i_myatom)**(1d0/3d0)) & 
                                      *vdw_hirshfeld_R0(species(i_myatom))
                endif 
              enddo ! end loop over atoms
              
              if(use_scalapack_mbd .and. use_scalapack_inversion) then  
                 call get_non_local_polarizability_tensor_MBD_rs_SCS_p(i_config,i_freq,alpha_omega,R_p,Rvdw_iso)
              else   
                 call get_non_local_polarizability_tensor_MBD_rs_SCS(i_config,i_freq,alpha_omega,R_p,Rvdw_iso)
              endif   


              if(i_config.eq.0) then
                
                write(info_str,'(2x,"| ",F10.6,6(e15.6))')casimir_omega(i_freq),&
                   coupled_molecule_pol(1,1,i_freq),coupled_molecule_pol(2,2,i_freq),coupled_molecule_pol(3,3,i_freq),&
                   coupled_molecule_pol(2,3,i_freq),coupled_molecule_pol(1,3,i_freq),coupled_molecule_pol(1,2,i_freq)
                call localorb_info(info_str, use_unit,'(A)',OL_norm)
              endif

          enddo ! end   over freq

          if(i_config.eq.0) then
            write(info_str,'(2x,A)')&
            "|----------------------------------------------------------------------------------"
            call localorb_info(info_str, use_unit,'(A)',OL_norm)

            if(n_periodic .eq. 0) then
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha(bohr^3)  of atom in molecule"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)
            else
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha(bohr^3)  of atom in unit cell"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)

            endif
              write(info_str,'(2x,A)')&
              "|----------------------------------------------------------------------------------"
              call localorb_info(info_str, use_unit,'(A)',OL_norm)
          
           
         
          do i_myatom=1,n_atoms,1
              write(info_str,'(2x,A,I4,2x,A,f12.6,6x,f12.6)')"| ATOM",i_myatom,trim(species_name(species(i_myatom))),&
                                      C6AA(coupled_atom_pol(:,i_myatom)),coupled_atom_pol(0,i_myatom)
              call localorb_info(info_str, use_unit,'(A)',OL_norm)
          enddo

          write(info_str,'(2x,A)')&
                "|----------------------------------------------------------------------------------"
                call localorb_info(info_str, use_unit,'(A)',OL_norm)
          endif      
          
          if(use_scalapack_mbd) then
             call get_mbd_energy_MBD_rsSCS_p(i_config,fd_conf_ene_mbd(i_config))
          else
             call get_mbd_energy_MBD_rsSCS(i_config,fd_conf_ene_mbd(i_config))
          endif 
            

          if(forces_on) then 
                if (i_config.gt.0.AND.MOD(i_config,6).eq.0) then
                     i_config_counter = i_config_counter + 1 
                 write(info_str,'(2x,A,I6,3(2x,E13.6))')"| MBD@rsSCS force on ATOM",i_config_counter,&
                 ((-1.d0*(fd_conf_ene_mbd(i_config-5)-fd_conf_ene_mbd(i_config-4)))/(2.0d0*fd_step_size))*(hartree / bohr),&
                 ((-1.d0*(fd_conf_ene_mbd(i_config-3)-fd_conf_ene_mbd(i_config-2)))/(2.0d0*fd_step_size))*(hartree / bohr),&
                 ((-1.d0*(fd_conf_ene_mbd(i_config-1)-fd_conf_ene_mbd(i_config  )))/(2.0d0*fd_step_size))*(hartree / bohr)
                 call localorb_info(info_str, use_unit,'(A)',OL_norm)
               endif 
          else
                ! No force evaluation requested    
                exit 
          endif
 

      enddo ! end over configurations
        
      ! MBD@rsSCS energy
      ene_mbd_rsSCS=fd_conf_ene_mbd(0)          

      write(info_str,'(2x,A)')&
      "|----------------------------------------------------------------------------------"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
     
     mbd_rsSCS_forces=0.d0
     if(forces_on) then
       ! Calculate finite difference forces
       i_config=0
       do i_myatom=1,n_atoms

         mbd_rsSCS_forces(1,i_myatom) = -1.d0*(fd_conf_ene_mbd(i_config+1)-fd_conf_ene_mbd(i_config+2))/(2.0d0*fd_step_size)
         mbd_rsSCS_forces(2,i_myatom) = -1.d0*(fd_conf_ene_mbd(i_config+3)-fd_conf_ene_mbd(i_config+4))/(2.0d0*fd_step_size)
         mbd_rsSCS_forces(3,i_myatom) = -1.d0*(fd_conf_ene_mbd(i_config+5)-fd_conf_ene_mbd(i_config+6))/(2.0d0*fd_step_size)
         i_config = i_config + 6  
       
       enddo 
     endif      
   
     if(use_scalapack_mbd .and. use_scalapack_inversion) then
       call deallocate_scalapack_scs_data()
     else
       call deallocate_mpi_scs_data()   
     endif

     if(use_scalapack_mbd) then
       call  exit_blacs()
     endif

     return

endsubroutine mbd_rsSCS_energy_and_forces
      
    subroutine initialize_finite_diff_infrastructure()
    !local variables  
    integer :: i_row, i_col, i_index,i_config   

    real*8,dimension(6,3) :: fd_operator
!   Compute all configuration required for finite difference discritizations 
!   Parallelization over finite difference configs 
    
      if(.not.allocated(fd_conf_cords))  then
        allocate(fd_conf_cords(0:6*n_atoms,3,n_atoms), stat=info)
        call check_allocation(info, 'fd_conf_cords')
      endif
      if(.not.allocated(fd_conf_ene_mbd)) then
        allocate(fd_conf_ene_mbd(0:6*n_atoms), stat=info)
        call check_allocation(info, 'fd_conf_ene_mbd')
      endif

    fd_step_size= 0.000002d0


    ! FIXME this looks ugly i know 
    fd_conf_cords   = 0.d0
    fd_operator     = 0.d0
    fd_operator(1,1)= fd_step_size  ; fd_operator(1,2)= 0.0d0         ; fd_operator(1,3)= 0.0d0  
    fd_operator(2,1)=-fd_step_size  ; fd_operator(2,2)= 0.0d0         ; fd_operator(2,3)= 0.0d0  
    fd_operator(3,1)= 0.0d0         ; fd_operator(3,2)= fd_step_size  ; fd_operator(3,3)= 0.0d0  
    fd_operator(4,1)= 0.0d0         ; fd_operator(4,2)=-fd_step_size  ; fd_operator(4,3)= 0.0d0  
    fd_operator(5,1)= 0.0d0         ; fd_operator(5,2)= 0.0d0         ; fd_operator(5,3)= fd_step_size
    fd_operator(6,1)= 0.0d0         ; fd_operator(6,2)= 0.0d0         ; fd_operator(6,3)=-fd_step_size 

    i_config=0
    !unpertubed geom  
    fd_conf_cords(0,:,:)=coords(:,:)
    do i_row=1,n_atoms
        do i_index=1,6,1
           i_config= i_config+1  
          do  i_col=1,n_atoms
              if(i_row.eq.i_col) then
                 fd_conf_cords(i_config,1,i_col)= coords(1,i_col) + fd_operator(i_index,1)
                 fd_conf_cords(i_config,2,i_col)= coords(2,i_col) + fd_operator(i_index,2) 
                 fd_conf_cords(i_config,3,i_col)= coords(3,i_col) + fd_operator(i_index,3) 
              else   
                 fd_conf_cords(i_config,:,i_col)=coords(:,i_col)
              endif   
          enddo
        enddo    
    enddo

    return
    endsubroutine initialize_finite_diff_infrastructure


    subroutine get_non_local_polarizability_tensor_MBD_rs_SCS(i_config,i_freq,alpha_omega,R_p,Rvdw_iso)

    implicit none
    real*8,dimension(n_atoms):: alpha_omega
    real*8,dimension(n_atoms):: R_p
    real*8,dimension(n_atoms):: Rvdw_iso
    integer ::i_index,j_index
    real*8,dimension(3,3)::TPP,matrix
    real*8,dimension(3) :: dxyz
    real*8,dimension(3) :: coord_curr, eigen
    real*8 :: r_ij
    real*8 :: r_pp
    real*8 :: Rvdw12
    real*8 :: beta
    integer :: i_row, i_col
    integer :: i_lattice, j_lattice, k_lattice
    integer :: periodic_cell_i,periodic_cell_j,periodic_cell_k
    integer :: i_config, i_freq

    !For LAPACK
    integer,dimension(3*n_atoms):: IPIV
    real*8,dimension(3*n_atoms):: WORK
    real*8,dimension(9):: work_tensor

    ! initio values


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

       !!!!!!!!!!!! 
       relay_matrix=0.d0
       !!!!!!!!!!!!


     ! compute relay matrix of  cluster or unit cell
       do i_row=1,n_atoms,1 !#1
       if(myid.eq.task_list(i_row)) then 
         do i_col=i_row,n_atoms,1 !#2
         TPP=0.d0
            if(.not. empty(i_row) .and. .not. empty(i_col)) then !empty    
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
               dxyz(:) =fd_conf_cords(i_config,:,i_col)-fd_conf_cords(i_config,:,i_row)
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
           endif!empty
         enddo   !#2
       endif
       enddo  !#1
       
  if (n_periodic .gt. 0) then
   
      
      ! Setting periodic_cell_* = 1 for respective vaccum direction 
      ! just to make sure real space integration of dipole field include all atoms in first/unit cell.

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
                     if(.not. empty(i_row) .and. .not. empty(i_col)) then !empty  
                     do i_col = i_row, n_atoms, 1 ! atom2 loop
                          coord_curr = 0.d0
                          dxyz = 0.d0
                                r_ij = 0.d0
                                TPP  = 0.d0
                               Rvdw12= 0.d0  
                          ! find the coordinate of images fd_conf_cords(i_config,:,i_col)
                            coord_curr(:) = fd_conf_cords(i_config,:,i_col) + i_lattice*lattice_vector(:,1) + &
                                                            j_lattice*lattice_vector(:,2) + &
                                                            k_lattice*lattice_vector(:,3)

                            dxyz(:) = fd_conf_cords(i_config,:,i_row)- coord_curr(:)
                            r_ij    = sqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                            r_pp    = sqrt(R_p(i_row)**2 + R_p(i_col)**2)
                            Rvdw12  = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
                            if(r_ij.le.mbd_scs_dip_cutoff/bohr) then
                            call SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
                              do i_index=1,3,1
                              do j_index=1,3,1
                                 relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=&
                                 relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index) + TPP(i_index,j_index)
                                 if(i_col.NE.i_row) then
                                  relay_matrix(3*i_col-3+j_index,3*i_row-3+i_index)=&
                                  relay_matrix(3*i_col-3+j_index,3*i_row-3+i_index) + TPP(i_index,j_index)
                                 endif
                              enddo
                              enddo
                            endif
                     enddo !atom2 loop
                    endif!empty
                   endif ! task
                  enddo !atom1 loop
               endif  !#1
            enddo
         enddo
      enddo
   
   endif ! Periodic check
   call sync_tensors(relay_matrix,3*n_atoms)
   call DGETRF(3*n_atoms, 3*n_atoms, relay_matrix, 3*n_atoms, IPIV, info)
   if(info.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call localorb_info(info_str, use_unit,'(A)',OL_norm)
   endif
   call check_info(info,"get_non_local_polarizability_tensor_MBD_rs_SCS","DGETRF")


   call DGETRI(3*n_atoms, relay_matrix, 3*n_atoms, IPIV, WORK,3*n_atoms,info )
   if(info.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call localorb_info(info_str, use_unit,'(A)',OL_norm)
   endif
   call check_info(info,"get_non_local_polarizability_tensor_MBD_rs_SCS","DGETRI")

 
 ! Tensor contraction 
   do i_row=1,n_atoms,1
      matrix=0.0
      do i_col=1,n_atoms,1
        do i_index=1,3,1
          do j_index=1,3,1
           matrix(i_index,j_index) = matrix(i_index,j_index) +&
                                     relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)
          enddo
        enddo
      enddo



     coupled_molecule_pol(:,:,i_freq) = coupled_molecule_pol(:,:,i_freq) + matrix(:,:) 
     call DSYEV('N','U',3,matrix,3,eigen,work_tensor,9,info)
     call check_info(info,"get_non_local_polarizability_tensor_MBD_rs_SCS","DSYEV")
     coupled_atom_pol(i_freq,i_row) =(sum(eigen))/3.d0
    
   enddo
  
   




   return
   endsubroutine get_non_local_polarizability_tensor_MBD_rs_SCS

  subroutine get_non_local_polarizability_tensor_MBD_rs_SCS_p(i_config,i_freq,alpha_omega,R_p,Rvdw_iso)
    real*8,dimension(n_atoms):: alpha_omega
    real*8,dimension(n_atoms):: R_p
    real*8,dimension(n_atoms):: Rvdw_iso
    integer ::i_index,j_index
    real*8,dimension(3,3)::TPP
    real*8,dimension(3,3)::MACROPOL, MICROPOL
    real*8,dimension(3):: MICROPOL_EIGEN
    real*8,dimension(3) :: dxyz
    real*8,dimension(3) :: coord_curr
     
    real*8 :: r_ij
    real*8 :: r_pp
    real*8 :: Rvdw12
    real*8 :: beta
    integer :: i_row, i_col
    integer:: i_freq
    integer :: i_lattice, j_lattice, k_lattice
    integer :: i_config,periodic_cell_i,periodic_cell_j,periodic_cell_k
    integer :: mbdwork1, mbdwork2

    
       select case (flag_xc)
          case (1) !PBE0
               beta=0.85
          case (6) ! PBE
               beta=0.83
          case (7) !HSE 
               beta=0.85

         case default  !FIXME add well described error
              beta=1.0
               write(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"this fxc defaulting it to 1.0 , which will give MBD energy"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
               write(info_str,'(A)')"WITHOUT short range damping"
               call localorb_info(info_str, use_unit,'(A)',OL_norm)
       end select

      ! compute relay tensor of  cluster or unit cell    
      do i_row=1,n_atoms
         call map_global_index_to_local_index(i_row,prow,nblocks,iproc,myi)
         if (myrow==iproc) then ! process grid block
            do i_col=1,n_atoms
            call map_global_index_to_local_index(i_col,pcol,nblocks,jproc,myj)
               if (mycol==jproc) then ! process grid block
                   if(.not. empty(i_row) .and. .not. empty(i_col)) then !empty
                   if(i_row.eq.i_col) then! Check atom -eq
                      !(myi,myj)<----is first coordinate  of local(1,1)_(atom_row,atom_col) for above i_row atom and i_col
                      !atom index tensors, respectively.
                      call get_identity(inv_alpha_tensor)
                      inv_alpha_tensor(:,:) = inv_alpha_tensor(:,:)*(1.d0/alpha_omega(i_row))
                      local_tensor_i=0
                      do i_index=myi,myi+2
                      local_tensor_i= local_tensor_i +1
                         local_tensor_j=0
                         do j_index=myj,myj+2
                         local_tensor_j=local_tensor_j +1
                         relay_matrix_local(i_index,j_index)=inv_alpha_tensor(local_tensor_i,local_tensor_j)
                         enddo
                      enddo
 
                   else! check and fill off-diagongal tensor for i_row and i_col index
                     dxyz(:) =fd_conf_cords(i_config,:,i_col)-fd_conf_cords(i_config,:,i_row)
                     r_ij = dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 +(dxyz(3))**2.0d0 )
                     r_pp = dSqrt(R_p(i_row)**2 + R_p(i_col)**2)
                     Rvdw12 = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
                     call SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
                     
                      local_tensor_i=0  
                      do i_index=myi,myi+2
                      local_tensor_i= local_tensor_i +1
                         local_tensor_j=0
                         do j_index=myj,myj+2
                         local_tensor_j=local_tensor_j +1
                         relay_matrix_local(i_index,j_index)=TPP(local_tensor_i,local_tensor_j)
                         enddo
                      enddo

                  endif! Check atom block    
                  endif!empty
               endif! Check on process block col
            enddo !loop atom col
         endif !Check on process block row
      enddo !loop atom row

        

     
  if (n_periodic .gt. 0) then
   


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
                    call map_global_index_to_local_index(i_row,prow,nblocks,iproc,myi)
                     if (myrow==iproc) then
                      if(.not. empty(i_row) .and. .not. empty(i_col)) then !empty
                       do i_col = 1, n_atoms, 1 ! atom2 loop
                          call map_global_index_to_local_index(i_col,pcol,nblocks,jproc,myj)
                          if (mycol==jproc) then

                          coord_curr = 0.d0
                          dxyz = 0.d0
                                r_ij = 0.d0
                                TPP  = 0.d0
                               !Rvdw12= 0.d0  
                          ! find the coordinate of images
                          coord_curr(:) = fd_conf_cords(i_config,:,i_col) + i_lattice*lattice_vector(:,1) + &
                                                            j_lattice*lattice_vector(:,2) + &
                                                            k_lattice*lattice_vector(:,3)

                          dxyz(:) = fd_conf_cords(i_config,:,i_row)- coord_curr(:)
                          r_ij    = sqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                          r_pp    = sqrt(R_p(i_row)**2 + R_p(i_col)**2)
                          Rvdw12  = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
                          if(r_ij.le.mbd_scs_dip_cutoff/bohr) then
                          call SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)

                         local_tensor_i=0  
                         do i_index=myi,myi+2
                            local_tensor_i= local_tensor_i +1
                            local_tensor_j=0
                               do j_index=myj,myj+2
                                local_tensor_j=local_tensor_j +1
                                relay_matrix_local(i_index,j_index)=relay_matrix_local(i_index,j_index)+TPP(local_tensor_i,local_tensor_j)
                              enddo
                          enddo

  
                          endif
                        endif ! jproc
                      
                     enddo !atom2 loop
                     endif!empty
                   endif! iproc
                  enddo !atom1 loop
               endif  !#1

            enddo ! loop imag -k--- +k  real space integration
         enddo ! loop imag -j--- +j  real space integration
      enddo ! loop imag -i--- +i  real space integration
   
   
   endif

   call pdgetrf (3*n_atoms,3*n_atoms,relay_matrix_local,1,1,ides_a,local_ipvi,info)
   call check_info(info,"get_non_local_polarizability_tensor_MBD_rs_SCS_p","Matrix factorization failed in pdgetrf")

   call pdgetri (3*n_atoms,relay_matrix_local,1,1,ides_a,local_ipvi,local_work,-1,local_iwork,-1,info)

   mbdwork1=local_work(1)
   mbdwork2=local_iwork(1)

   if (allocated(local_work)) deallocate(local_work, stat=info)
   if (info /= 0) print *, "local_WORK: Deallocation request denied"

   if (allocated(local_iwork)) deallocate(local_iwork, stat=info)
   if (info /= 0) print *, "local_IWORK: Deallocation request denied"

   allocate(local_work(mbdwork1))
   allocate(local_iwork(mbdwork2))

   call pdgetri (3*n_atoms,relay_matrix_local,1,1,ides_a,local_ipvi,local_work,mbdwork1,local_iwork,mbdwork2,info)
   call check_info(info,"get_non_local_polarizability_tensor_MBD_rs_SCS_p","Matrix inversion failed in pdgetri")


   MACROPOL=0.d0    
   !Tensor contraction to get \alpha_cluster_or_unit_cell(3,3)
   do i_row=1,n_atoms
     MICROPOL(:,:)=0.d0
     call map_global_index_to_local_index(i_row,prow,nblocks,iproc,myi)
     if (myrow==iproc) then
       do i_col=1,n_atoms
         call map_global_index_to_local_index(i_col,pcol,nblocks,jproc,myj)
         if (mycol==jproc) then
           !!!!!!!!!!!!!!!!!!!! POL MOL OR UNIT CELL !!!!!!!!!!!!!!!!!!!!!!!
           local_tensor_i=0
           do i_index=myi,myi+2
             local_tensor_i= local_tensor_i +1
             local_tensor_j=0
             do j_index=myj,myj+2
               local_tensor_j=local_tensor_j +1
               !write(use_unit,*)i_row,i_col,myi,myj,local_tensor_i,local_tensor_j
               MACROPOL(local_tensor_i,local_tensor_j)=MACROPOL(local_tensor_i,local_tensor_j)+&
                                                              relay_matrix_local(i_index,j_index)
             enddo
           enddo

           !!!!!!!!!!!!!!!!!!!! Progected \alpha_isotropic MOL OR UNIT CELL !!!!!!!!!!!!!!!!!!!!!!!
           local_tensor_i=0
           do i_index=myi,myi+2
             local_tensor_i= local_tensor_i + 1
             local_tensor_j=0
             do j_index=myj,myj+2
               local_tensor_j=local_tensor_j + 1
               MICROPOL(local_tensor_i,local_tensor_j)=MICROPOL(local_tensor_i,local_tensor_j)+&
                                                                  relay_matrix_local(i_index,j_index)
             enddo
           enddo
         endif
       enddo
     endif
     call sync_tensors(MICROPOL,3)
     call DSYEV('N','U',3,MICROPOL,3,MICROPOL_EIGEN,ALPHA_WORK,9,info)
     call check_info(info,"get_non_local_polarizability_tensor_MBD_rs_SCS_p","DSYEV")
     ! only isotropic value
     coupled_atom_pol(i_freq,i_row)=sum(MICROPOL_EIGEN)/3.0
   enddo !
   call sync_tensors(MACROPOL,3)
   coupled_molecule_pol(:,:,i_freq)= MACROPOL(:,:) 
   return
   endsubroutine  get_non_local_polarizability_tensor_MBD_rs_SCS_p

 
   subroutine get_mbd_energy_MBD_rsSCS(i_config,ene_mbd_rsSCS)

     ! local vars ! optimized for finite difference forces  
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian 
     real*8,dimension(:,:),allocatable:: coords_SL
     real*8,dimension(:),allocatable:: cfdm_eigenvalues
     real*8,dimension(:),allocatable:: omega_cfdm_SL
     real*8,dimension(:),allocatable:: R_vdw_SL
     real*8,dimension(:),allocatable::alpha_eff_SL
     real*8,dimension(:),allocatable:: WORK
     real*8,dimension(:),allocatable:: task_list_SL  
     real*8,dimension(3,3):: TPP,lattice_vector_SL
     real*8,dimension(3)::dxyz,coord_curr
     real*8:: r_ij,ene_mbd_rsSCS
     real*8:: C6_free
     real*8:: alpha_free
     real*8:: R_vdw_free
     real*8:: Rvdw_12
     real*8:: beta
     real*8:: CFDM_prefactor
     real*8:: E_int
     real*8:: E_nonint
     real*8:: sigma
     integer :: info,i_atom,j_atom,i,j
     integer :: i_config
     integer :: SL_i,SL_j,SL_k ! MBD super cell index's
     integer :: i_index, j_index ! 3x3 block index
     integer :: i_lattice,j_lattice,k_lattice ! lattice increament index
     integer :: periodic_cell_i,periodic_cell_j,periodic_cell_k
     integer :: NIMAG
     integer :: n_atoms_SL
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
      if(i_config.eq.0) then
      if(forces_on) then
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy and forces..."
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      else
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy ..."          
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      endif
      endif

      if(.NOT.allocated(cfdm_hamiltonian)) then
        allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL), stat=info) 
        call check_allocation(info, 'cfdm_hamiltonian','get_mbd_energy_MBD_rsSCS',3*n_atoms_SL,3*n_atoms_SL)
      endif
      if(.NOT.allocated(cfdm_eigenvalues)) then
        allocate(cfdm_eigenvalues(3*n_atoms_SL), stat=info)
        call check_allocation(info, 'cfdm_eigenvalues')
      endif 
      if(.NOT.allocated(coords_SL)) then
        allocate(coords_SL(3,n_atoms_SL), stat=info)
        call check_allocation(info, 'coords_SL')
      endif 
      if(.NOT.allocated(omega_cfdm_SL)) then
        allocate(omega_cfdm_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'omega_cfdm_SL')
      endif
      if(.NOT.allocated(R_vdw_SL)) then
        allocate(R_vdw_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'R_vdw_SL')
      endif
      if(.NOT.allocated(alpha_eff_SL)) then
        allocate(alpha_eff_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'alpha_eff_SL')
      endif
      if(.NOT.allocated(WORK)) then
        allocate(WORK((3*n_atoms_SL)*(3+(3*n_atoms_SL)/2)), stat=info)
        call check_allocation(info, 'WORK','get_mbd_energy_MBD_rsSCS',(3*n_atoms_SL)*(3+(3*n_atoms_SL)/2))
      endif
      if(.NOT.allocated(task_list_SL)) then
        allocate(task_list_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'task_list_SL')
      endif

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

                  R_vdw_SL(i_atom)= R_vdw_free*((coupled_atom_pol(0,i_atom)/alpha_free)**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))
                  coords_SL(:,i_atom)   = fd_conf_cords(i_config,:,i_atom)
                  alpha_eff_SL(i_atom) = coupled_atom_pol(0,i_atom)
                  else
                  R_vdw_SL(i_atom)= vdw_hirshfeld_R0(species(i_atom))* & 
                     ((coupled_atom_pol(0,i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))
                  coords_SL(:,i_atom)   = fd_conf_cords(i_config,:,i_atom)
                  alpha_eff_SL(i_atom) = coupled_atom_pol(0,i_atom)
                  endif

                  if (alpha_eff_SL(i_atom) < 0d0) then
                    write(info_str,'(A,I6,A,E16.7)') &
                      "*** WARNING: Negative polarizability for atom ", &
                      i_atom, ": ", alpha_eff_SL(i_atom)
                      call localorb_info(info_str, use_unit,'(A)',OL_norm)
                      alpha_eff_SL(i_atom) = 0.001d0
                      omega_cfdm_SL(i_atom) = &
                          (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom)) / &
                          (alpha_eff_SL(i_atom)**2)
                      R_vdw_SL(i_atom) = R_vdw_free
                  endif 
         enddo
      else ! periodic case

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
      if(i_config.eq.0) then
      if(forces_on) then
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy and forces..."
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      else
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy ..."          
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      endif
  


      write(info_str,'(2X,A,i3,A,i3,A,i3,A,f6.2,A)') &
         "| Creating super cell of dimension",  SL_i," X ", SL_j," X ", &
         SL_k, " in MBD calculation using" ,mbd_supercell_cutoff, &
         "  Angstrom radius"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,"| containing", I6,A)')n_atoms_SL, "  atom"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

      endif  


      if(.NOT.allocated(cfdm_hamiltonian)) then
        allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL), stat=info)
        call check_allocation(info, 'cfdm_hamiltonian','get_mbd_energy_MBD_rsSCS',3*n_atoms_SL,3*n_atoms_SL)
      endif
 
      if(.NOT.allocated(cfdm_eigenvalues)) then
        allocate(cfdm_eigenvalues(3*n_atoms_SL), stat=info)
        call check_allocation(info, 'cfdm_eigenvalues')
      endif 
      if(.NOT.allocated(coords_SL)) then
        allocate(coords_SL(3,n_atoms_SL), stat=info)
        call check_allocation(info, 'coords_SL')
      endif 
      if(.NOT.allocated(omega_cfdm_SL)) then
        allocate(omega_cfdm_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'omega_cfdm_SL')
      endif
      if(.NOT.allocated(R_vdw_SL)) then
        allocate(R_vdw_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'R_vdw_SL')
      endif
      if(.NOT.allocated(alpha_eff_SL)) then
        allocate(alpha_eff_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'alpha_eff_SL')
      endif
      if(.NOT.allocated(WORK)) then
        allocate(WORK(3*n_atoms_SL*(3+(3*n_atoms_SL)/2)), stat=info)
        call check_allocation(info, 'WORK','get_mbd_energy_MBD_rsSCS',3*n_atoms_SL*(3+(3*n_atoms_SL)/2))
      endif
      if(.NOT.allocated(task_list_SL)) then
        allocate(task_list_SL(n_atoms_SL), stat=info)
        call check_allocation(info, 'task_list_SL')
      endif
      
      cfdm_hamiltonian = 0.0d0 
      cfdm_eigenvalues = 0.0d0
 
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
                                       
                  coords_SL(:,j_atom) =fd_conf_cords(i_config,:,i_atom) + i_lattice*lattice_vector(:,1) + &
                                                           j_lattice*lattice_vector(:,2) + &
                                                           k_lattice*lattice_vector(:,3)   
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0 
                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then   
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                C6_free,alpha_free,R_vdw_free)

                  R_vdw_SL(j_atom)=R_vdw_free*((coupled_atom_pol(0,i_atom)/alpha_free)**(1.d0/3.d0))
                  alpha_eff_SL(j_atom) = coupled_atom_pol(0,i_atom)
                  omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))
                  else
                  R_vdw_SL(j_atom)=vdw_hirshfeld_R0(species(i_atom)) & 
                     *((coupled_atom_pol(0,i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  alpha_eff_SL(j_atom) = coupled_atom_pol(0,i_atom)
                  omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))

                  endif 
                  if (alpha_eff_SL(i_atom) < 0d0) then
                    write(info_str,'(A,I6,A,E16.7)') &
                      "*** WARNING: Negative polarizability for atom ", &
                      i_atom, ": ", alpha_eff_SL(i_atom)
                      call localorb_info(info_str, use_unit,'(A)',OL_norm)
                      alpha_eff_SL(i_atom) = 0.001d0
                      omega_cfdm_SL(i_atom) = &
                          (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom)) / &
                          (alpha_eff_SL(i_atom)**2)
                      R_vdw_SL(i_atom) = R_vdw_free
                  endif   
              enddo 
            enddo            
         enddo            
       enddo            
      lattice_vector_SL(:,1)  = lattice_vector(:,1)*SL_i 
      lattice_vector_SL(:,2)  = lattice_vector(:,2)*SL_j 
      lattice_vector_SL(:,3)  = lattice_vector(:,3)*SL_k 
      endif ! periodic vs cluster 

      ! Construct cfdm hamiltonian     
      do i_atom=1,n_atoms_SL,1 !$1
        if(myid.eq.task_list_SL(i_atom)) then
         do j_atom=i_atom,n_atoms_SL,1 !$2
         if(.not. empty(i_atom) .and. .not. empty(j_atom)) then !empty      
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
        endif!empty
         enddo !$2
         endif
      enddo !$1
    
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
                  if(myid.eq.task_list_SL(i_atom)) then
                  ! LOOP GOES OVER UPPPER TRIANGLE OF HAMILTONIAN 
                  if(.not. empty(i_atom) .and. .not. empty(j_atom)) then !empty
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
                        CFDM_prefactor = omega_cfdm_SL(i_atom)*&
                                         omega_cfdm_SL(j_atom)*&
                                         sqrt(alpha_eff_SL(i_atom)*&
                                              alpha_eff_SL(j_atom))
                        do i_index=1,3,1
                          do j_index=1,3,1
                             cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index)=&
                             cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          enddo
                        enddo
                      endif 

                  enddo !$2
                 endif!empty
                endif
                enddo !$3
               endif !$4
            enddo !$5
         enddo !$6
      enddo !$7  
   endif

    call sync_tensors(cfdm_hamiltonian,3*n_atoms_SL)
    info=0
    LWORK=3*n_atoms_SL*(3+(3*n_atoms_SL)/2)

    !call print_dsyev_info("get_mbd_energy_MBD_rsSCS",3*n_atoms_SL)  
    call DSYEV('N','U',3*n_atoms_SL,cfdm_hamiltonian,3*n_atoms_SL,cfdm_eigenvalues,WORK,LWORK,info)
    call check_info(info,"get_mbd_energy_MBD_rsSCS","DSYEV")
    

    E_int=0.0d0
    E_nonint =0.0d0
    ene_mbd_rsSCS =   0.0
    NIMAG=0 
    if (info.eq.0) then
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
    else
           write(info_str,'(A)')"***Error:- Digonalization of CFDM hamiltonian failed"
           call localorb_info(info_str, use_unit,'(A)',OL_norm)
           !THIS OUTPUT IS ONLY FOR BL for info what is going wrong
           write(info_str,'(A,I6)')"***Errorcode = ", info
           call localorb_info(info_str, use_unit,'(A)',OL_norm)
           call aims_stop ()

    endif
   
      if(allocated(cfdm_hamiltonian))          deallocate(cfdm_hamiltonian)
      if(allocated(cfdm_eigenvalues))          deallocate(cfdm_eigenvalues)
      if(allocated(coords_SL))                 deallocate(coords_SL)
      if(allocated(omega_cfdm_SL))             deallocate(omega_cfdm_SL)
      if(allocated(R_vdw_SL))                  deallocate(R_vdw_SL)
      if(allocated(alpha_eff_SL))              deallocate(alpha_eff_SL)  
  return
  end subroutine get_mbd_energy_MBD_rsSCS 

subroutine get_mbd_energy_MBD_rsSCS_p(i_config,ene_mbd_rsSCS)

     ! local vars  
     real*8,dimension(:,:),allocatable:: coords_SL
     real*8,dimension(:),allocatable:: omega_cfdm_SL
     real*8,dimension(:),allocatable:: R_vdw_SL
     real*8,dimension(:),allocatable::alpha_eff_SL
     real*8,dimension(:),allocatable:: cfdm_ham_local_WORK
     integer::mbd_lwork, i_config 
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
     integer :: info,i_atom,j_atom
     integer :: SL_i,SL_j,SL_k ! MBD super cell index's
     integer :: i_index, j_index ! 3x3 block index
     integer :: i_lattice,j_lattice,k_lattice ! lattice increament index
     integer :: periodic_cell_i,periodic_cell_j,periodic_cell_k
     integer :: NIMAG
     integer :: n_atoms_SL
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
               call localorb_info(info_str, use_unit,'(A)',OL_high)
               write(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
               call localorb_info(info_str, use_unit,'(A)',OL_high)
               write(info_str,'(A)')"WITHOUT short range damping"
               call localorb_info(info_str, use_unit,'(A)',OL_high)
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
      if(i_config.eq.0) then
      if(forces_on) then
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy and forces..."
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      else
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy ..."          
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      endif
 
      Write(info_str,'(2x,A,I6,A,I6,A,I6,A)')"| Maping MBD tensor on",prow*pcol,&
      " blacs processes each containing",numroc(3*n_atoms_SL,3,myrow,0,prow)/3," X ",numroc(3*n_atoms_SL,3,mycol,0,pcol)/3," tensor blocks(atoms)"
      call localorb_info(info_str, use_unit,'(A)',OL_norm) 
      endif


      call allocate_mbd_data(n_atoms_SL)
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))

          do i_atom =1,n_atoms
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0
                  C6_free=0.d0
                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then  
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                C6_free,alpha_free,R_vdw_free)

                  R_vdw_SL(i_atom)= R_vdw_free*((coupled_atom_pol(0,i_atom)/alpha_free)**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))
                  coords_SL(:,i_atom)   = fd_conf_cords(i_config,:,i_atom)
                  alpha_eff_SL(i_atom) = coupled_atom_pol(0,i_atom)
                  else
                  R_vdw_SL(i_atom)= vdw_hirshfeld_R0(species(i_atom))* & 
                     ((coupled_atom_pol(0,i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))
                  coords_SL(:,i_atom)   = fd_conf_cords(i_config,:,i_atom)
                  alpha_eff_SL(i_atom) = coupled_atom_pol(0,i_atom)
                  endif   
                  
                  if (alpha_eff_SL(i_atom) < 0d0) then
                    write(info_str,'(A,I6,A,E16.7)') &
                      "*** WARNING: Negative polarizability for atom ", &
                      i_atom, ": ", alpha_eff_SL(i_atom)
                      call localorb_info(info_str, use_unit,'(A)',OL_norm)
                      alpha_eff_SL(i_atom) = 0.001d0
                      omega_cfdm_SL(i_atom) = &
                          (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom)) / &
                          (alpha_eff_SL(i_atom)**2)
                      R_vdw_SL(i_atom) = R_vdw_free
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
        
      if(i_config.eq.0) then
      if(forces_on) then
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy and forces..."  
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      else
      write(info_str,'(2X,A)')"| Computing MBD@rsSCS energy ..."  
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      endif
      write(info_str,'(2X,A,i3,A,i3,A,i3,A,f6.2,A)') &
         "| Creating super cell of dimension",  SL_i," X ", SL_j," X ", &
         SL_k, " in MBD calculation using" ,mbd_supercell_cutoff, &
         "  Angstrom radius"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2x,"| containing", I6,A)')n_atoms_SL, "  atom"
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

      Write(info_str,'(2x,A,I6,A,I6,A,I6,A)')"| Maping MBD tensor on",prow*pcol,&
      " blacs processes each containing",numroc(3*n_atoms_SL,3,myrow,0,prow)/3," X ",numroc(3*n_atoms_SL,3,mycol,0,pcol)/3," tensor blocks(atoms)"
      call localorb_info(info_str, use_unit,'(A)',OL_norm) 
      endif  

      call allocate_mbd_data(n_atoms_SL)
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))

      ! Get all the parameter  required for MBD with super cell calculation
       j_atom = 0
       do i_lattice = 0, SL_i-1,1
         do j_lattice = 0, SL_j-1,1
            do k_lattice = 0, SL_k-1,1
              do i_atom =1,n_atoms !<--- atom index
                  j_atom = j_atom +1  !<--- dummy index supercell
                  coords_SL(:,j_atom) =fd_conf_cords(i_config,:,i_atom) + i_lattice*lattice_vector(:,1) + &
                                                              j_lattice*lattice_vector(:,2) + &
                                                           k_lattice*lattice_vector(:,3)     
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0 
                  if (.not.vdw_hirshfeld_data_external(species(i_atom))) then   
                  call get_vdw_param(species_element(species(i_atom)),species_z(species(i_atom)),&
                                                C6_free,alpha_free,R_vdw_free)

                  R_vdw_SL(j_atom)=R_vdw_free*((coupled_atom_pol(0,i_atom)/alpha_free)**(1.d0/3.d0))
                  alpha_eff_SL(j_atom) = coupled_atom_pol(0,i_atom)
                  omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))
                  else
                  R_vdw_SL(j_atom)=vdw_hirshfeld_R0(species(i_atom)) & 
                     *((coupled_atom_pol(0,i_atom)/vdw_hirshfeld_alpha(species(i_atom)))**(1.d0/3.d0))
                  alpha_eff_SL(j_atom) = coupled_atom_pol(0,i_atom)
                  omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom))/(coupled_atom_pol(0,i_atom)*coupled_atom_pol(0,i_atom))
                  endif

                  if (alpha_eff_SL(i_atom) < 0d0) then
                    write(info_str,'(A,I6,A,E16.7)') &
                      "*** WARNING: Negative polarizability for atom ", &
                      i_atom, ": ", alpha_eff_SL(i_atom)
                      call localorb_info(info_str, use_unit,'(A)',OL_norm)
                      alpha_eff_SL(i_atom) = 0.001d0
                      omega_cfdm_SL(i_atom) = &
                          (4.d0/3.d0)*C6AA(coupled_atom_pol(:,i_atom)) / &
                          (alpha_eff_SL(i_atom)**2)
                      R_vdw_SL(i_atom) = R_vdw_free
                  endif 

              enddo 
            enddo            
         enddo            
       enddo            
      lattice_vector_SL(:,1)  = lattice_vector(:,1)*SL_i 
      lattice_vector_SL(:,2)  = lattice_vector(:,2)*SL_j 
      lattice_vector_SL(:,3)  = lattice_vector(:,3)*SL_k 
      endif 

      ! Construct CFDM hamiltonian matrix 
      do i_atom=1,n_atoms_SL,1 !#1
         call map_global_index_to_local_index(i_atom,prow,nblocks,iproc,myi)
         if (myrow==iproc) then !#2 process grid block
           ! BL: ELPA needs full matrix
           ! do j_atom=i_atom,n_atoms_SL,1 !#1
           do j_atom=1,n_atoms_SL,1 !#1
             call map_global_index_to_local_index(j_atom,pcol,nblocks,jproc,myj)
             if (mycol==jproc) then !#4 process grid block  
              if(.not. empty(i_atom) .and. .not. empty(j_atom)) then !empty
               if(i_atom.eq.j_atom) then  !#5
                 call get_identity(inv_alpha_tensor)
                 inv_alpha_tensor(:,:) = inv_alpha_tensor(:,:)*omega_cfdm_SL(i_atom)**2.0
                 local_tensor_i=0
                 do i_index=myi,myi+2
                    local_tensor_i= local_tensor_i +1
                    local_tensor_j=0
                    do j_index=myj,myj+2
                      local_tensor_j=local_tensor_j +1
                      cfdm_ham_local(i_index,j_index)=inv_alpha_tensor(local_tensor_i,local_tensor_j)
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
                   sqrt(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))
               
                 local_tensor_i=0           
                 do i_index=myi,myi+2
                   local_tensor_i= local_tensor_i +1
                   local_tensor_j=0
                     do j_index=myj,myj+2
                        local_tensor_j=local_tensor_j +1
                        cfdm_ham_local(i_index,j_index)=TPP(local_tensor_i,local_tensor_j)*CFDM_prefactor
                     enddo
                 enddo
               endif !#5
              endif!empty
             endif ! #4 process grid block
           enddo !#3
         endif ! #2 process grid block
       enddo ! #1
      

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
 
      do i_lattice = -periodic_cell_i, periodic_cell_i,1          !#7
         do j_lattice = -periodic_cell_j, periodic_cell_j,1         !#6  
            do k_lattice = -periodic_cell_k, periodic_cell_k,1        !#5 
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then!#4

                do i_atom=1,n_atoms_SL,1 !#3
                call map_global_index_to_local_index(i_atom,prow,nblocks,iproc,myi)
                 if (myrow==iproc) then !#2 process grid block
                   ! BL: ELPA needs full matrix
                   ! do j_atom=i_atom,n_atoms_SL,1 !#1
                   if(.not. empty(i_atom) .and. .not. empty(j_atom)) then !empty 
                   do j_atom=1,n_atoms_SL,1 !#1
                     call map_global_index_to_local_index(j_atom,pcol,nblocks,jproc,myj)
                       if (mycol==jproc) then !#0 

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
                        local_tensor_i=0 
                        do i_index=myi,myi+2
                           local_tensor_i= local_tensor_i +1
                           local_tensor_j=0
                           do j_index=myj,myj+2
                              local_tensor_j=local_tensor_j +1
                              cfdm_ham_local(i_index,j_index)= cfdm_ham_local(i_index,j_index) +&
                                                                        TPP(local_tensor_i,local_tensor_j)*CFDM_prefactor
                              
                           enddo
                        enddo

                      endif
                    endif !#0

                  enddo !#1
                 endif!empty
                 endif !#2process grid block
                enddo !#3

               endif !#4
            enddo !#5
         enddo !#6
      enddo !#7

   endif

   if (use_elpa_mbd) then
     !call print_dsyev_info("get_mbd_energy_MBD_rsSCS_elpa",3*n_atoms_SL) 
     call solve_evp_real_2stage_2013(3*n_atoms_SL, 3*n_atoms_SL, &
          cfdm_ham_local, myArows, cfdm_ham_local_eigen, &
          cfdm_ham_local_vectors, myArows, nblocks, mbd_comm_rows, &
          mbd_comm_cols, mpi_comm_global)
   else

     allocate (cfdm_ham_local_WORK(1))
     
     CALL PDSYEV('V','U',3*n_atoms_SL,cfdm_ham_local,1,1,ides_mbd,cfdm_ham_local_eigen,&
                 cfdm_ham_local_vectors,1,1,ides_mbd_vec,cfdm_ham_local_WORK,-1,info)
     call check_info(info,"get_mbd_energy_MBD_rsSCS_p","PDSYEV")
     mbd_lwork=cfdm_ham_local_WORK(1) 

     if (allocated(cfdm_ham_local_WORK)) deallocate(cfdm_ham_local_WORK, stat=info)
     if (info /= 0) print *, "cfdm_ham_local_WORK: Deallocation request denied"

     allocate (cfdm_ham_local_WORK(mbd_lwork))   
     !call print_dsyev_info("get_mbd_energy_MBD_rsSCS_p",3*n_atoms_SL)  
     CALL PDSYEV('V','U',3*n_atoms_SL,cfdm_ham_local,1,1,ides_mbd,cfdm_ham_local_eigen,&
                 cfdm_ham_local_vectors,1,1,ides_mbd_vec,cfdm_ham_local_WORK,mbd_lwork,info)
     call check_info(info,"get_mbd_energy_MBD_rsSCS_p","PDSYEV")
   end if
     
             E_int = 0.d0
         E_nonint  = 0.d0
    ene_mbd_rsSCS  = 0.d0
            NIMAG  = 0 
    do i_atom =1,n_atoms_SL
      E_nonint = E_nonint + omega_cfdm_SL(i_atom)    
    enddo

    do i_atom =1,3*n_atoms_SL
      if(cfdm_ham_local_eigen(i_atom).ge.0.d0) then
        E_int = E_int + dsqrt(cfdm_ham_local_eigen(i_atom))
      else
        NIMAG= NIMAG +1  
      endif 
    enddo

    ene_mbd_rsSCS = ((0.5*E_int)-(1.5* E_nonint))/(SL_i*SL_j*SL_k)
    if(NIMAG.gt.0) then
      write(info_str,'(A,I4,A)')"***WARNING: found ",NIMAG," negative eigenvalues in MBD energy calculation."
      call localorb_info(info_str, use_unit,'(A)',OL_high)  
    endif

    if(allocated(coords_SL))                 deallocate(coords_SL)
    if(allocated(omega_cfdm_SL))             deallocate(omega_cfdm_SL)
    if(allocated(R_vdw_SL))                  deallocate(R_vdw_SL)
    if(allocated(alpha_eff_SL))              deallocate(alpha_eff_SL)
    call deallocate_mbd_data  
   
    return
  
  endsubroutine get_mbd_energy_MBD_rsSCS_p

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

   write(info_str,'(2X,A,A,A,I10,A,I10,A)') "| ", trim(caller), &
     " solves ", rank, " x ", rank, " eigenvalue problem."
   call localorb_info(info_str, use_unit,'(A)',OL_norm)


  endsubroutine print_dsyev_info

   
  endmodule mbd_rsSCS_numerical_forces 


