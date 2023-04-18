      subroutine diagonalize_green_fn_dmft_GW &
          ( green_fn_freq, &
           hamiltonian_GW,&
           diagonal_green_fn, i_k_point,new_chemical_potential)


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
      use localorb_io, only: use_unit


     !ARGUMENTS
      implicit none
      complex*16 green_fn_freq(n_basis, n_basis, nomega)
      complex*16 hamiltonian_GW(n_basis, n_basis, nomega)
!output
      complex*16 diagonal_green_fn (n_states, n_states, nomega ) !the diagonalized green's function
      integer i_k_point

!auxiliary
!      complex*16 aux_green_fn (n_basis, n_basis, nomega)
!      complex*16 green_tmp (n_basis, n_basis)
!      complex*16 ovlp_full (n_basis, n_basis)
!      complex*16 aux_matr (n_basis, n_states)
      complex*16, dimension(:,:), allocatable ::  green_tmp
      complex*16, dimension(:,:), allocatable ::  ovlp_full
      complex*16, dimension(:,:), allocatable ::  aux_matr
      complex*16, dimension(:,:,:), allocatable ::  aux_green_fn
      complex*16, dimension(:,:,:), allocatable ::  diag_ham_GW
      complex*16, dimension(:), allocatable :: ipiv
      complex*16, dimension(:), allocatable :: work

!      complex*16 DUMMY(n_basis, n_basis)
!      complex*16 , dimension(:), allocatable :: WORK
!      complex*16 RWORK(2*n_basis)
!      integer LWORK
!      integer INFO
      character*2 iter
      character*14 filename
      logical output
      complex*16 aux_KS_eigenvector(n_basis, n_states)
      real*8 new_chemical_potential
 
!counters
      integer i_freq 
      integer info 
      integer n, i, i_task
      integer i_state, k,i_basis,j_basis, i_index
      
      allocate (green_tmp(n_basis, n_basis))       
      allocate (ovlp_full(n_basis, n_basis))       
      allocate (aux_matr(n_basis, n_states))       
      allocate (aux_green_fn (n_basis, n_basis, nomega))       
      allocate (diag_ham_GW (n_states, n_states, nomega))       
!      real*8 max_noise
!      real *4 ran
!      integer i_seed
      
!      if(myid.eq.0)then
!        write(use_unit,*) " "
!        write(use_unit,*) "Transform G(w) in the KS basis"
!        write(use_unit,*) " "
!      endif

!      print *, 'here 0'
      !aux_matr (:,:) = (0.d0,0.d0)
      aux_green_fn (:,:,:) =  (0.d0,0.d0)
      diagonal_green_fn(:,:,:) = (0.d0,0.d0)
      diag_ham_GW(:,:,:) = (0.d0,0.d0)

!      print *, 'green_fn_freq ', green_fn_freq 
!      print *, '#######################################################'
!      print *, ' '

      do i_basis = 1, n_basis
        do i_state = 1, n_states

if(.not.real_eigenvectors)then
           aux_KS_eigenvector(i_basis,i_state) = &
           KS_eigenvector_complex(i_basis,i_state,1,i_k_point)

else
           aux_KS_eigenvector(i_basis,i_state) = &
           DCMPLX(KS_eigenvector(i_basis,i_state,1,1))
endif
        enddo
      enddo

!      n = 0
      ovlp_full(:,:) = 0.d0
!      do i = 1, n_basis
!        ovlp_full(1:i,i) = overlap_matrix(n+1:n+i)
!        ovlp_full(i,1:i-1) = overlap_matrix(n+1:n+i-1)
!        n = n+i
!      enddo
      i_index = 0
      do j_basis = 1, n_basis, 1
        do i_basis = 1, j_basis, 1
         i_index = i_index+1
         ovlp_full (j_basis , i_basis) =  DCMPLX(overlap_matrix (i_index))
         ovlp_full (i_basis , j_basis) =  ovlp_full (j_basis , i_basis)
        enddo
      enddo
!      print*, ovlp_full

!      if(myid.eq.0.and..true.)then
!        print *, n_loc_grid
!        do i_task = 1, n_tasks,1 
!          do i_freq = 1, n_loc_grid
!            write(use_unit,*) i_task, i_freq, map_index (i_task, i_freq)
!          enddo
!        enddo
!      endif
!        do i_state = 1, n_states, 1
!          write(use_unit, *) i_state, diagonal_green_fn (i_state,i_state, 1)
!        enddo
!      print *, 'Transfrom G in the KS basis'
      !do k = 1, n_loc_grid, 1
      do i_freq = 1, nomega, 1
       ! i_freq = map_index (myid+1, k)
        !print *, i_freq
!        if(i_freq.gt.0)then

         green_tmp (:,:)= (0.d0,0.d0)
         aux_matr (:,:)=  (0.d0,0.d0)
!       print *,'green_tmp' , green_tmp
!       print *,'aux_matr' ,  aux_matr
!      print *, '#######################################################'
!      print *, ' '
 
         call zgemm ('N','N',n_basis,n_basis, n_basis,&
           (1.d0,0.d0), green_fn_freq (:,:,i_freq), n_basis, &
           ovlp_full, n_basis, (0.d0,0.d0),        &
           green_tmp, n_basis)
!         if(i_freq.eq.1)print * ,'green_tmp', green_tmp
!      print *, '#######################################################'
!      print *, ' '

         call zgemm ('N','N',n_basis,n_basis, n_basis,&
            (1.d0,0.d0), ovlp_full, n_basis, &
            green_tmp, n_basis, (0.d0,0.d0), &
            aux_green_fn (:,:,i_freq), n_basis)
    
!       if(i_freq.eq.1)print * ,'aux_green_fn', aux_green_fn (:,:,i_freq)

        call zgemm ('N','N',n_basis,n_states, n_basis,&
            !(1.d0,0.d0), aux_green_fn (:,:,i_freq), n_basis, &
            !(1.d0,0.d0), green_fn_freq (:,:,i_freq), n_basis, &
            (1.d0,0.d0), hamiltonian_GW (:,:,i_freq), n_basis, &
            aux_KS_eigenvector, n_basis, (0.d0,0.d0),        &
            aux_matr,n_basis)
    
!        print *, 'here 3' 
!        if(i_freq.eq.1)print * ,'aux_matr', aux_matr
 
        call zgemm ('C','N',n_states,n_states, n_basis,&
            (1.d0,0.d0),aux_KS_eigenvector , n_basis, &
            aux_matr, n_basis, (0.d0,0.d0), &
            diag_ham_GW (:,:,i_freq) , n_states)
!        print *, 'here 4'
!        do i_state = 1, n_states, 1
!          write(use_unit, *) i_state, diagonal_green_fn (i_state,i_state, i_freq)
!        enddo
!        endif

     diagonal_green_fn(:,:,i_freq)  = (0.d0,0.d0)

 do i_basis = 1, n_states, 1
 do j_basis = 1, n_states, 1
  diagonal_green_fn(i_basis,j_basis,i_freq) = diagonal_green_fn(i_basis,j_basis,i_freq) + &
       1./(((0.d0,1.d0)*omega(i_freq))-diag_ham_GW (i_basis,j_basis,i_freq))!+(0.d0,0.001))!-DCMPLX(chemical_potential))
!  diagonal_green_fn(i_basis,j_basis,i_freq) = diagonal_green_fn(i_basis,j_basis,i_freq) + 1./(diag_ham_GW (i_basis,j_basis,i_freq))!+(0.d0,0.001))!-DCMPLX(chemical_potential))
!  diagonal_green_fn(i_basis,j_basis,i_freq) = diagonal_green_fn(i_basis,j_basis,i_freq) + 1./(((0.d0,1.d0)*omega(i_freq))-diag_ham_GW (i_basis,j_basis,i_freq)+(0.d0,0.001))!+DCMPLX(chemical_potential))!-DCMPLX(chemical_potential))
!  diagonal_green_fn(i_basis,j_basis,i_freq) = diagonal_green_fn(i_basis,j_basis,i_freq) + 1.d0/((0.d0,1.d0)*omega(i_freq)-diag_ham_GW (i_basis,j_basis,i_freq)+(0.d0,0.0001))
!  diagonal_green_fn(i_basis,j_basis,i_freq) = diagonal_green_fn(i_basis,j_basis,i_freq) + 1.d0/((0.d0,1.d0)*omega(i_freq)-KS_eigenvalue (i_basis,1,i_k_point))
 enddo 
 enddo


if(.false.)then
       if (.not. allocated(ipiv))then
          allocate(ipiv(n_states))
       endif

       if (.not. allocated(work))then
          allocate(work (n_states))
       endif

     diagonal_green_fn(:,:,i_freq)  = diag_ham_GW(:,:,i_freq)

         call zgetrf( n_states, n_states, diagonal_green_fn(:,:,i_freq), &
                    n_states , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
  !          stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call zgetri(n_states, diagonal_green_fn(:,:,i_freq), n_states, &
                 ipiv, work, n_states, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif
       if (allocated(ipiv))then
          deallocate(ipiv)
       endif

       if (allocated(work))then
          deallocate(work)
       endif
endif

      enddo
      !print *, 'Transfrom G in the KS basis - DONE'

!      if(myid.eq.0)then
!        open (22, file = "test_diag_G.dat")
!        do i_state = 1, n_states, 1
!          write(use_unit, *) i_state, diag_ham_GW (i_state,i_state,1), KS_eigenvalue(i_state,1,i_k_point)
!        enddo
!        close(22)
!      endif

!      do i_state = 1, n_states, 1
!        call sync_matrix_complex (&
!         diagonal_green_fn  (i_state, 1:n_states, 1:nomega), &
!         n_states, nomega )
!      enddo


      deallocate (green_tmp)
      deallocate (ovlp_full)
      deallocate (aux_matr)
      deallocate (aux_green_fn )

      !add some noise to G to test the analytic continuation
      !set to FALSE for any calculation!
!      seed = 1
!      if(.false.)then
!        do i_basis = 1, n_states, 1
!          do j_basis = 1, n_states, 1
!            do i_freq = 1, nomega, 1
!               max_noise = abs(diagonal_green_fn(i_basis,i_basis,1))
!               diagonal_green_fn(i_basis,i_basis,i_freq) = &
!               diagonal_green_fn(i_basis,i_basis,i_freq) +&
!               max_noise * ( ran(i_seed)-0.5d0)/500
!            enddo
!          enddo
!        enddo
!      endif
!
!      output = .false.
!      if(myid.eq.0)then
!      if(output) then
!        do i_basis = 1, n_states, 1
!          if( i_basis.lt.10 ) then
!             write(iter,'(A,I1)') "0",i_basis
!          else
!             write(iter,'(I2)') i_basis
!          endif
!          filename = "gree_diag_"//iter//".dat"
!          open(77, file=filename)
!            do i_freq = 1, nomega, 1
!              write(77,*) omega(i_freq), &
!                       real(diagonal_green_fn(i_basis,i_basis,i_freq)),&
!                      aimag(diagonal_green_fn(i_basis,i_basis,i_freq))
!            enddo
!          close(77)
!        enddo
!      endif
!      endif


      end subroutine diagonalize_green_fn_dmft_GW
