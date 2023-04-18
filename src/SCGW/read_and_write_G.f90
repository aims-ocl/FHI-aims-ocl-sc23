      subroutine write_green_fn ()

      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use gw_para
      use constants
      use mpi_tasks
      use synchronize_mpi
      use numerical_utilities 
      use scgw_grid
      use poles_fit
      use gt
      use scgw_allocations 

      implicit none

      
      if(myid.eq.0) then 
         open(55, file="Gt.dat")
         open(56, file="Gw.dat")
         do i_spin = 1, n_spin, 1 
           i_index = 0
           do i_basis = 1, n_basis, 1 
             do j_basis = 1, i_basis, 1 
               do i_tau = -ntau, ntau, 1 
                  write (55,*) green_fn_time (i_basis, j_basis, i_tau, i_spin)
               enddo
               do i_freq = 1, nomega, 1 
                  write (56,*) green_fn_freq (i_basis, j_basis, i_freq, i_spin)
               enddo
             enddo
           enddo
         enddo
         close(55)
         close(56)
      endif

      end subroutine write_green_fn

!-------------------------------------------------
      subroutine read_green_fn ()

      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use gw_para
      use constants
      use mpi_tasks
      use synchronize_mpi
      use numerical_utilities 
      use scgw_grid
      use poles_fit
      use gt
      use scgw_allocations 

      implicit none

      open(55, file="Gt.dat", status='old', action='read')
      open(56, file="Gw.dat", status='old', action='read')
      do i_spin = 1, n_spin, 1
        do i_basis = 1, n_basis, 1
          do j_basis = 1, i_basis, 1
            do i_tau = -ntau, ntau, 1
               read (55,*) green_fn_time_ens (i_basis, j_basis, i_tau, i_spin)
               green_fn_time_ens (j_basis, i_basis, i_tau, i_spin) =&
               green_fn_time_ens (i_basis, j_basis, i_tau, i_spin)
            enddo
            do i_freq = 1, nomega, 1
               read (56,*) green_fn_freq_ens (i_basis, j_basis, i_freq, i_spin)
               green_fn_freq_ens (j_basis, i_basis, i_freq, i_spin) =&
               green_fn_freq_ens (i_basis, j_basis, i_freq, i_spin)
            enddo
          enddo
        enddo
      enddo
      close(55)
      close(56)

      end subroutine read_green_fn
