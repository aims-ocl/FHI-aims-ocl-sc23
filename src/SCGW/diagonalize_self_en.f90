      subroutine diagonalize_self_en &
          ( self_energy_freq, &
           self_energy_diag)


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
      use scgw_grid



     !ARGUMENTS
      implicit none
      complex*16 self_energy_freq(n_basis, n_basis, nomega)
!output
      complex*16 self_energy_diag (n_states, n_states, nomega ) !the diagonalized green's function

!auxiliary
      complex*16 aux_green_fn (n_basis, n_basis, nomega)
      complex*16 green_tmp (n_basis, n_basis)
      complex*16 ovlp_full (n_basis, n_basis)
!      complex*16 DUMMY(n_basis, n_basis)
!      complex*16 , dimension(:), allocatable :: WORK
!      complex*16 RWORK(2*n_basis)
!      integer LWORK
!      integer INFO
      character*2 iter
      character*14 filename
      logical output
      complex*16 aux_matr (n_basis, n_states)
      complex*16 aux_KS_eigenvector(n_basis, n_states)
 
!counters
      integer i_freq 
      integer n, i, i_task
      integer i_state, k,i_basis,j_basis
      real*8 max_noise
      real *4 ran
      integer i_seed
      
!      if(myid.eq.0)then
!        write(use_unit,*) " "
!        write(use_unit,*) "Transform G(w) in the KS basis"
!        write(use_unit,*) " "
!      endif

!      aux_matr (:,:)= 0.d0
!      aux_green_fn (:,:,:) = 0.d0
!      diagonal_green_fn(:,:,:) = 0.d0
      aux_KS_eigenvector(:,:) = KS_eigenvector(:,:,1,1)
!
!      n = 0
!      ovlp_full(:,:) = 0.d0
!      do i = 1, n_basis
!        ovlp_full(1:i,i) = overlap_matrix(n+1:n+i)
!        ovlp_full(i,1:i-1) = overlap_matrix(n+1:n+i-1)
!        n = n+i
!      enddo

!      if(myid.eq.0.and..true.)then
!        print *, n_loc_grid
!        do i_task = 1, n_tasks,1 
!          do i_freq = 1, n_loc_grid
!            write(use_unit,*) i_task, i_freq, map_index (i_task, i_freq)
!          enddo
!        enddo
!      endif

      do k = 1, n_loc_grid, 1
        i_freq = map_index (myid+1, k)
        if(i_freq.gt.0)then
!         green_tmp (:,:)= 0.d0
!         call zgemm ('N','N',n_basis,n_basis, n_basis,&
!           1.d0, green_fn_freq (1,1,i_freq), n_basis, &
!           ovlp_full, n_basis, 0.d0,        &
!           green_tmp, n_basis)
!
!        aux_green_fn (:,:,i_freq) = 0.d0
!        call zgemm ('N','N',n_basis,n_basis, n_basis,&
!           1.d0, ovlp_full, n_basis, &
!           green_tmp, n_basis, 0.d0, &
!           aux_green_fn (1,1,i_freq), n_basis)
   
       aux_matr (:,:)= (0.d0,0.d0)
       call zgemm ('N','N',n_basis,n_states, n_basis,&
           (1.d0,0.d0), self_energy_freq (1,1,i_freq), n_basis, &
           aux_KS_eigenvector, n_basis, (0.d0,0.d0),        &
           aux_matr , n_basis)
   
       self_energy_diag (:,:,i_freq) = (0.d0,0.d0)
       call zgemm ('T','N',n_states,n_basis, n_basis,&
           (1.d0,0.d0),aux_KS_eigenvector , n_basis, &
           aux_matr, n_basis, (0.d0,0.d0),        &
           self_energy_diag(1,1,i_freq)  , n_states)
        endif
      enddo

      do i_state = 1, n_states, 1
        call sync_matrix_complex (&
         self_energy_diag  (i_state, 1:n_states, 1:nomega), &
         n_states, nomega )
      enddo


      output = .false.
      if(myid.eq.0)then
      if(output) then
        do i_basis = 1, n_states, 1
          if( i_basis.lt.10 ) then
             write(iter,'(A,I1)') "0",i_basis
          else
             write(iter,'(I2)') i_basis
          endif
          filename = "gree_diag_"//iter//".dat"
          open(77, file=filename)
            do i_freq = 1, nomega, 1
              write(77,*) omega(i_freq), &
                       real(self_energy_diag(i_basis,i_basis,i_freq)),&
                      aimag(self_energy_diag(i_basis,i_basis,i_freq))
            enddo
          close(77)
        enddo
      endif
      endif


      end subroutine diagonalize_self_en
