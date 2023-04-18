      subroutine diagonalize_green_fn_dmft &
          ( green_fn_freq, diagonal_green_fn, i_spin, k_summed_overlap_matr, &
            inv_k_summed_overlap_matr,i_k_point, aux_KS_eigenvector, &
            full_hamiltonian, on_site_xc_matr, loc_self_enrg, hartree_pot_LDA, &
            full_k_xc_matr, diag_ham_k)


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
!output
      complex*16 diagonal_green_fn (n_states, n_states, nomega ) !the diagonalized green's function
      complex*16, dimension(:,:,:), allocatable ::test_diagonal_green_fn !the diagonalized green's function
      complex*16 k_summed_overlap_matr (n_basis, n_basis ) !the diagonalized green's function
      complex*16 inv_k_summed_overlap_matr (n_basis, n_basis ) !the diagonalized green's function
      complex*16 full_hamiltonian (n_basis, n_basis ) 
      complex*16 full_k_xc_matr (n_basis, n_basis ) 
      complex*16 diag_ham_k (n_states, n_states ) 
      complex*16 loc_self_enrg (n_basis, n_basis ) 
      real*8 on_site_xc_matr (n_basis, n_basis ) 
      real*8 hartree_pot_LDA(n_basis, n_basis ) 
      integer i_spin
      integer i_k_point
      integer j_k_point
      integer i_k
   integer info

!auxiliary
!      complex*16 aux_green_fn (n_basis, n_basis, nomega)
!      complex*16 green_tmp (n_basis, n_basis)
!      complex*16 ovlp_full (n_basis, n_basis)
!      complex*16 aux_matr (n_basis, n_states)
      complex*16, dimension(:,:), allocatable ::  green_tmp
      complex*16, dimension(:,:), allocatable ::  ham_tmp
      complex*16, dimension(:,:), allocatable ::  aux_hamiltonian
      complex*16, dimension(:,:), allocatable ::  ovlp_full
      complex*16, dimension(:,:), allocatable ::  aux_matr
      complex*16, dimension(:,:), allocatable ::  aux_ham_matr
      complex*16, dimension(:,:,:), allocatable ::  aux_green_fn
      complex*16, dimension(:,:), allocatable ::  aux_ham
      complex*16, dimension(:,:), allocatable ::  diagonal_ham 
      complex*16, dimension(:,:,:), allocatable :: KS_eigenvector_k_sum 
      complex*16, dimension(:,:,:), allocatable :: aux_diagonal_green_fn 
      complex*16, dimension(:,:,:), allocatable :: diagonal_green_fn_temp 
      complex*16, dimension(:,:), allocatable :: output_ham 

!      complex*16 DUMMY(n_basis, n_basis)
!      complex*16 , dimension(:), allocatable :: WORK
!      complex*16 RWORK(2*n_basis)
!      integer LWORK
!      integer INFO
      character*2 iter
      character*14 filename
      logical output
      complex*16 aux_KS_eigenvector(n_basis, n_states)
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work
 
!counters
      integer i_freq 
      integer n, i, i_task
      integer i_state, k,i_basis,j_basis, i_index
      
      allocate (green_tmp(n_basis, n_basis))       
      allocate (ham_tmp(n_basis, n_basis))       
      allocate (ovlp_full(n_basis, n_basis))       
      allocate (aux_matr(n_basis, n_states))       
      allocate (aux_ham_matr(n_basis, n_states))       
      allocate (test_diagonal_green_fn (n_states, n_states, nomega))       
      allocate (aux_green_fn (n_basis, n_basis, nomega))       
      allocate (diagonal_green_fn_temp (n_states, n_states, nomega))       
      allocate (aux_ham (n_basis, n_basis))       
      allocate (diagonal_ham (n_states, n_states))       
      allocate (aux_hamiltonian (n_basis, n_basis))       
      allocate (KS_eigenvector_k_sum(n_basis, n_states, n_spin))       
      allocate (aux_diagonal_green_fn(n_states, n_states, nomega))       
      allocate (output_ham(n_basis, n_basis))       
!      real*8 max_noise
!      real *4 ran
!      integer i_seed
       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif
      
!      if(myid.eq.0)then
!        write(use_unit,*) " "
!        write(use_unit,*) "Transform G(w) in the KS basis"
!        write(use_unit,*) " "
!      endif

!      print *, 'here 0'
      !aux_matr (:,:) = (0.d0,0.d0)
      aux_green_fn (:,:,:) =  (0.d0,0.d0)
      diagonal_green_fn(:,:,:) = (0.d0,0.d0)

!      print *, 'green_fn_freq ', green_fn_freq 
!      print *, '#######################################################'
!      print *, ' '


!      KS_eigenvector_k_sum(:,:,:) =(0.d0,0.d0)
!              i_k = 0
!              do i_k_point = 1, n_k_points
!                 if (myid == modulo(i_k_point, n_tasks) .and. myid <= n_k_points) then
!                    i_k = i_k + 1
!                     do i_basis = 1, n_basis
!                        do i_state = 1, n_states
!                        KS_eigenvector_k_sum(i_basis,i_state,1) = KS_eigenvector_k_sum(i_basis,i_state,i_spin) + &
!                        DCMPLX(KS_eigenvector_complex(i_basis,i_state,1,i_k))
!                       enddo
!                       enddo
!                  end if
!               enddo

!     KS_eigenvector_k_sum = (1./n_k_points)*KS_eigenvector_k_sum

!if(myid.eq. 0) write(use_unit,*) 'KS_eigenvector_k_sum', KS_eigenvector_k_sum

!      do i_basis = 1, n_basis
!        do i_state = 1, n_states
!           aux_KS_eigenvector(i_basis,i_state) = &
!!           ((1./(n_k_points**(0.5)))*KS_eigenvector_complex(i_basis,i_state,1,1))
!           (KS_eigenvector_complex(i_basis,i_state,1,1))!i_k_point))
!        enddo
!      enddo



!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!if(myid.eq.0) then
!write(use_unit,*) 'use_scalapack', use_scalapack
!write(use_unit,*) 'collect_eigenvectors', collect_eigenvectors
!endif

!          if (use_scalapack) then

!             In the ScaLapack case (n_tasks >= n_k_points), each node only
!             cares about scalapack_wrapper:my_k_point, and works on it in
!             scalapack_wrapper:eigenvec.  Obviously, n_k_points_task=1.
!             if (collect_eigenvectors) then

!                      do i_basis = 1, n_basis
!                        do i_state = 1, n_states
!                          aux_KS_eigenvector(i_basis,i_state) = &
!           ((1./(n_k_points**(0.5)))*KS_eigenvector_complex(i_basis,i_state,1,1))
!                          (KS_eigenvector_complex(i_basis,i_state,1,1))!i_k_point))
!                        enddo
!                       enddo




!                The array KS_eigenvectors{_complex}(:,:,:,1:1) contains the
!                corresponding eigencoefficients.
!             else
!                It doesn't!
!             end if
!
!          else
!
!              In this case, loops over the k-points can be done as follows:
!
!              i_k = 0
!              do j_k_point = 1, n_k_points
!                 if (myid == modulo(j_k_point, n_tasks) .and. myid <= n_k_points) then
!                    i_k = i_k + 1
!               do i_basis = 1, n_basis
!                 do i_state = 1, n_states
!                   aux_KS_eigenvector(i_basis,i_state) = &
!           ((1./(n_k_points**(0.5)))*KS_eigenvector_complex(i_basis,i_state,1,1))
!                   (KS_eigenvector_complex(i_basis,i_state,1,i_k))!i_k_point))
!                  enddo
!                enddo

!!                    ... occ_numbers(:,:, i_k_point) ... KS_eigenvector(:,:,:, i_k) ...
!                 end if
!              end do
!
!              This is how it is done in get_KS_orbitals_p0(), which is /the/
!              reference by definition.  Be aware that myid starts at 0 whereas
!              i_k_point starts at 1.
!
!              Together with the condition (myid <= n_k_points) [nota bene:
!              "<="] this leads to the following distribution:
!              if (n_tasks <= n_k_points) then
!                 Each task is responsible for all k-points with
!                 myid == modulo(i_k_point, n_tasks).
!              else
!                 Only tasks with (1 <= myid <= n_k_points) are responsible for any
!                 k-point, namely for the one with (myid == i_k_point).
!              end if

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------












!if(myid.eq. 0)then
!write(use_unit,*) 'n_tasks', n_tasks
!write(use_unit,*) 'n_k_poits', n_k_points
!endif
!write(use_unit,*) 'KS_eigenvector_complex', KS_eigenvector_complex(:,:,1, i_k_point)



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
!        do i_basis = 1, n_basis, 1
!          write(use_unit, *) i_basis, on_site_xc_matr(i_basis,i_basis), &
!                               hartree_pot_LDA(i_basis,i_basis), &
!                               loc_self_enrg(i_basis,i_basis)
!
!        enddo
!      print *, 'Transfrom G in the KS basis'
!      do k = 1, n_loc_grid, 1
!       i_freq = map_index (myid+1, k)
!!       !print *, i_freq
!        if(i_freq.gt.0)then
!      	 green_tmp (:,:)= (0.d0,0.d0)
!         aux_matr (:,:)=  (0.d0,0.d0)
!       print *,'green_tmp' , green_tmp
!       print *,'aux_matr' ,  aux_matr
!      print *, '#######################################################'
!      print *, ' '

           aux_hamiltonian(:,:) =(0.d0,0.d0)

     do i_basis= 1, n_basis
        do j_basis= 1, n_basis

            aux_hamiltonian(i_basis,j_basis) = aux_hamiltonian(i_basis,j_basis) &
                                            + (full_hamiltonian(i_basis,j_basis)- &
           !                                    full_k_xc_matr(i_basis,j_basis) -  &
                                               on_site_xc_matr(i_basis,j_basis) +  &
                                               !hartree_pot_LDA(i_basis,j_basis) +  &
                                               loc_self_enrg(i_basis,j_basis))
       enddo
     enddo

!test of hermiticity of the Hamiltonian
output_ham(:,:) = (0.d0,0.d0)
         call zgemm ('C','N',n_basis,n_basis, n_basis,&
           (1.d0,0.d0), aux_hamiltonian (:,:), n_basis, &
            aux_hamiltonian (:,:), n_basis, (0.d0,0.d0),&
           output_ham, n_basis)

!do i_basis = 1 , n_basis,1
!do j_basis = 1 , n_basis,1

!write(use_unit,*) "output_ham",i_basis,j_basis, output_ham(i_basis,j_basis)

!enddo
!enddo
!stop


ham_tmp(:,:) =(0.d0,0.d0)
aux_ham(:,:) =(0.d0,0.d0)
aux_ham_matr(:,:) =(0.d0,0.d0)
diagonal_ham(:,:) =(0.d0,0.d0)
         call zgemm ('N','N',n_basis,n_basis, n_basis,&
           (1.d0,0.d0), aux_hamiltonian (:,:), n_basis, &
            inv_k_summed_overlap_matr, n_basis, (0.d0,0.d0),&
            !ovlp_full, n_basis, (0.d0,0.d0),&
            ham_tmp, n_basis)

         call zgemm ('N','N',n_basis,n_basis, n_basis,&
            (1.d0,0.d0),  inv_k_summed_overlap_matr, n_basis, &
            !(1.d0,0.d0),  ovlp_full, n_basis, &
            ham_tmp, n_basis, (0.d0,0.d0), &
            aux_ham (:,:), n_basis)

        call zgemm ('N','N',n_basis,n_states, n_basis,&
!           (1.d0,0.d0), aux_ham (:,:), n_basis, &
           (1.d0, 0.d0), aux_hamiltonian (:,:), n_basis, &
            aux_KS_eigenvector(:,:), n_basis, (0.d0,0.d0),&
            aux_ham_matr,n_basis)

        call zgemm ('T','N',n_states,n_states, n_basis,&
            (1.d0,0.d0),aux_KS_eigenvector , n_basis, &
            aux_ham_matr, n_basis, (0.d0,0.d0), &
            diagonal_ham (:,:) , n_states)


      test_diagonal_green_fn(:,:,:)= (0.d0,0.d0)
!      do i_k_point = 1, n_k_point
        do i_freq = 1, nomega
           do i_state = 1, n_states
            test_diagonal_green_fn(i_state,i_state,i_freq) = &
                 test_diagonal_green_fn(i_state,i_state,i_freq) &
               + 1.d0/((0.d0,1.d0)*omega(i_freq) &
                     + chemical_potential -  KS_eigenvalue(i_state,1,i_k_point))
!diagonal_ham(i_state,i_state)+((0.d0,0.0)))
           enddo
        enddo
!if(i_k_point.eq.12) then
!do i_basis = 1 , n_states,1
!do j_basis = 1 , n_states,1

!write(use_unit,*) "diagonal ham",i_basis,j_basis, diagonal_ham(i_basis,j_basis)

!enddo
!enddo
!endif
!stop


!diag_ham_k(:,:) = diagonal_ham(:,:)




     aux_green_fn (:,:,:) =  (0.d0,0.d0)
     diagonal_green_fn_temp(:,:,:) = (0.d0,0.d0)

      do i_freq = 1, nomega, 1
         green_tmp (:,:)= (0.d0,0.d0)
         aux_matr (:,:)=  (0.d0,0.d0)



         call zgemm ('N','N',n_basis,n_basis, n_basis,&
           (1.d0,0.d0), green_fn_freq (:,:,i_freq), n_basis, &
            k_summed_overlap_matr, n_basis, (0.d0,0.d0),        &
           green_tmp, n_basis)

         call zgemm ('N','N',n_basis,n_basis, n_basis,&
            (1.d0,0.d0),  k_summed_overlap_matr, n_basis, &
            green_tmp, n_basis, (0.d0,0.d0), &
            aux_green_fn (:,:,i_freq), n_basis)

        call zgemm ('N','N',n_basis,n_states, n_basis,&
!            (1.d0,0.d0), aux_green_fn (:,:,i_freq), n_basis, &
            (1.d0,0.d0), green_fn_freq (:,:,i_freq), n_basis, &
            !(1.d0,0.d0), green_tmp (:,:), n_basis, &
           aux_KS_eigenvector(:,:), n_basis, (0.d0,0.d0),        &
            aux_matr,n_basis)


        call zgemm ('C','N',n_states,n_states, n_basis,&
            (1.d0,0.d0),aux_KS_eigenvector , n_basis, &
            aux_matr, n_basis, (0.d0,0.d0), &
            !diagonal_green_fn_temp (:,:,i_freq) , n_states)
            diagonal_green_fn (:,:,i_freq) , n_states)
      enddo
!diag_ham_k(:,:) = diagonal_green_fn(:,:,100)



!      diagonal_green_fn(:,:,:) = diagonal_green_fn_temp(:,:,:)
!
!
!! at each frequency point evaluate the LU factorization, and check for errors
!       do i_freq = 1, nomega, 1
!       !do i_freq = 1, nomega, 1
!         call zgetrf( n_states, n_states, diagonal_green_fn(:,:,i_freq), &
!                    n_states , ipiv, info )

!        if (info.ne.0) then
!          if (myid.eq.0)then
!            write(use_unit,*) " * Failure of LU decomposition at",&
!                            " frequency point " , i_freq
!            write(use_unit,*) " * Error info = ", info
!!  !          stop
!          endif
!        endif
!! this inverts the matrix in the LU factorized form
!        call zgetri(n_states, diagonal_green_fn(:,:,i_freq), n_states, &
!     !   call zgetri(n_basis, gf_non_loc(:,:,i_k_point,i_freq), n_basis, &
!                 ipiv, work, n_states, info)
!
!      if (info.ne.0) then
!         if(myid.eq.0)then
!            write(use_unit,*) " * Failure of matrix inversion at",&
!                            " frequency point " , i_freq
!            write(use_unit,*) " * Error info = ", info
!          endif
!        endif
!
!       enddo
    

diag_ham_k(:,:) = diagonal_green_fn(:,:,100)







!diag_ham_k(:,:) = diagonal_ham(:,:)

!diagonal_green_fn(:,:,:) = (0.d0,0.d0)
!        !do i_state = 1, n_states, 1
!        do i_freq = 1, nomega, 1
 
!          diagonal_green_fn(:,:,i_freq) = diagonal_green_fn(:,:,i_freq) + aux_diagonal_green_fn(1,1,i_freq)
   
!        enddo


!diagonal_green_fn (:,:,:) = diagonal_green_fn (1,1,:)


!        do i_state = 1, n_states, 1
!          write(use_unit, *)i_k_point, i_state, real(diagonal_ham (i_state,i_state)), KS_eigenvalue(i_state,1,i_k_point)
!        enddo
if(i_k_point.eq.1) then
        do i_state = 1, n_states, 1
          write(use_unit, *) i_k_point, i_state, &
             real(test_diagonal_green_fn (i_state,i_state,1)), &
             real(diagonal_green_fn (i_state,i_state,1))
        enddo
endif
!stop
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


      end subroutine diagonalize_green_fn_dmft
