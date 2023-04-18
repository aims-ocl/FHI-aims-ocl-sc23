           subroutine transform_to_KS_basis_complex_parallel(NAO_matr,&
                                                    KS_matr)






     use dimensions
!     use runtime_choices
     use species_data
     use pbc_lists
     use physics
     use prodbas
     use constants
     use mpi_tasks
!     use synchronize_mpi
     use localized_basbas
     use timing

      implicit none

      integer i_k_point,i_k,i
      complex*16, dimension(n_basis,n_basis,n_k_points) :: NAO_matr
      complex*16, dimension(n_states,n_states) :: KS_matr 
      complex*16,dimension(:,:), allocatable :: aux_matr
      complex*16,dimension(:,:,:), allocatable :: KS_eigenvector_tmp
      integer, allocatable ::  myid_tmp(:)
     !find the node where the eigenvector is stored
     allocate(myid_tmp(n_k_points))
     allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin))

    i_k = 0


!        KS_matr(:,:) = (0.d0,0.d0)

   do i_k_point =1, n_k_points,1

!           myid_tmp(i_k_point) = MOD(i_k_point, n_tasks)

!       if (myid .eq. myid_tmp(i_k_point)) then
           if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
!if(myid.eq.0)&
!write(use_unit,*) ' i_k_point,myid, MOD(i_k_point, n_tasks)',i_k_point,myid, MOD(i_k_point, n_tasks)
!          endif    
            i_k = i_k + 1

             KS_eigenvector_tmp(:,:,:) = KS_eigenvector_complex(:,:,:,i_k)

           else
             ! zero temp. KS eigenvector on all other threads, prior to allreduce
             KS_eigenvector_tmp(:,:,:) = 0.d0
           end if
           call sync_eigenvector_complex(KS_eigenvector_tmp)

          
!            KS_eigenvector_tmp(:,:,:,i_k) = KS_eigenvector_complex(:,:,:,i_k)
       !call sync_eigenvector_complex(KS_eigenvector_tmp(:,:,:,i_k))
!          endif    
!       call sync_eigenvector_complex(KS_eigenvector_tmp(:,:,:,i_k_point))
!          enddo    
        KS_matr(:,:) = (0.d0,0.d0)
!   do i_k_point =1, n_k_points,1

       if (.not. allocated(aux_matr))then
          allocate(aux_matr(n_basis,n_states))
       endif

!       call sync_eigenvector_complex(KS_eigenvector_tmp(:,:,:,i_k_point))

        aux_matr(:,:) = (0.d0,0.d0)

      if(real_eigenvectors)then

        call zgemm ('N','N',n_basis,n_states, n_basis,&
            (1.d0,0.d0),NAO_matr(:,:,i_k_point), n_basis, &
            DCMPLX(KS_eigenvector(:,:,1,i_k)), n_basis, (0.d0,0.d0),&
            aux_matr(:,:),n_basis)


        call zgemm ('T','N',n_states,n_states, n_basis,&
            (1.d0,0.d0),DCMPLX(KS_eigenvector(:,:,1,i_k)) , n_basis, &
            aux_matr(:,:), n_basis, (0.d0,0.d0), &
            KS_matr(:,:) , n_states)

       else

        call zgemm ('N','N',n_basis,n_states, n_basis,&
            (1.d0,0.d0),NAO_matr(:,:,i_k_point), n_basis, &
            !(1.d0,0.d0),NAO_matr(:,:,i_k), n_basis, &
            !KS_eigenvector_complex(:,:,1,i_k), n_basis, (0.d0,0.d0),&
            KS_eigenvector_tmp(:,:,1), n_basis, (0.d0,0.d0),&
            aux_matr(:,:),n_basis)


        call zgemm ('C','N',n_states,n_states, n_basis,&
            !(1.d0,0.d0),KS_eigenvector_complex(:,:,1,i_k) , n_basis, &
            (1.d0,0.d0),KS_eigenvector_tmp(:,:,1) , n_basis, &
            aux_matr(:,:), n_basis, (0.d0,0.d0), &
            KS_matr(:,:) , n_states)

       endif

       if (allocated(aux_matr))then
          deallocate(aux_matr)
       endif
!        KS_eigenvector_complex(1:n_basis,1:n_states,1,i_k) = KS_matr(1:n_basis,1:n_states)

      !endif
       !call sync_eigenvector_complex(KS_eigenvector_complex(:,:,:,i_k_point))

! if(i_k_point.eq.1)then
  do i = 1,n_states,1
    if(myid.eq.0)&
     write(use_unit,*) 'KS_matr',i_k_point,i, real(KS_matr(i,i)), KS_eigenvalue(i,1,i_k_point)
   enddo

! endif

   enddo

  deallocate(myid_tmp)
 end subroutine transform_to_KS_basis_complex_parallel
