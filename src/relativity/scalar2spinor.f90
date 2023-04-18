
! This file contains subroutines relating to fully-relativistic integrations.
! The integrations we generated in integrate_hamiltonian_matrix_p2 are scalar integrations, which should be assembled 
! into spinor integrations here. Spinor integrations are what finally go into the diagonalization process.
! -- Rundong Zhao, Nov. 2018 @ Durham, NC
!
 subroutine transform_scalar2spinor(iop, n_spin, n_centers_basis, ld_dirac_scalar, mat_scalar, mat_spinor)
  use pbc_lists
  use dimensions,  only: n_basis, n_k_points, n_centers_basis_I
  use rel_x2c_mod, only: dim_matrix_rel, n_basis_small, n_centers_basis_I_small, &
                         large_comp_idx, small_comp_idx, large_comp_clb, small_comp_clb 
  use localorb_io, only: localorb_info, use_unit, OL_norm
  use mpi_tasks,   only: myid,n_tasks
  use synchronize_mpi_basic, only: sync_vector_complex
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  implicit none
  integer,intent(in) :: iop ! For large comp. part (VTS matrix): iop=1,2,3; for small comp. part (W and Ss matrix): iop=4,5.
  integer,intent(in) :: n_spin, n_centers_basis, ld_dirac_scalar
  real*8,intent(in) :: mat_scalar(ld_dirac_scalar,n_spin)  ! the input scalar integration
  complex*16,intent(out) :: mat_spinor(dim_matrix_rel,dim_matrix_rel,3,n_k_points,n_spin) ! the output spinor integration

  integer :: s=1 ! Comes from BDF code, in most cases, s=1
  integer :: n_useful_cells ! Not every cell in n_cells (see pbc_lists) is useful, we need to filter the useless ones out,
                            ! or the memory will explode.
  integer :: i_spin, i_index, i_k_point, i_basis_1, i_basis_2, i_cell_1, i_cell_2, i_cell, i_full_cell_1, i_full_cell_2
  integer :: n, i,j, temp
  integer :: row, column
  real*8,allocatable :: mat_cell(:,:) ! a temporary square matrix, mat_cell(n,n), n=n_basis or n=n_basis_small
  real*8,allocatable :: mat_full(:,:) ! a temporary square matrix, mat_full(n_centers_basis_I,n_centers_basis_I)
  real*8,allocatable :: mat_sq(:,:,:) ! a temporary square matrix that will be transformed to the upper triangular matrix mat_spinor
  complex*16,allocatable :: A(:,:), B(:,:),test(:,:)

  integer,allocatable :: usefulcell_to_fullcell(:), usefulcell_to_fullcell_s(:)
  integer,allocatable :: cbasis_to_useful_cell(:), cbasis_to_useful_cell_s(:)
  integer,allocatable :: cell_index_tmp(:)
  integer,allocatable :: cellbasis_to_fullbasis(:,:), cellbasis_to_fullbasis_s(:,:)

  real*8 :: dis,coords(3), t1,t2
  real*8 :: thresh=1.d-8
  character*200 :: info_str
  complex*16 :: temp_comp

  if(iop.eq.1 .or. iop.eq.2 .or. iop.eq.3) n=n_basis
  if(iop.eq.4 .or. iop.eq.5) n=n_basis_small
  allocate( cell_index_tmp(n_cells) )

  if (iop.eq.1 .or. iop.eq.2 .or. iop.eq.3)then

     cell_index_tmp = 0
     n_useful_cells = 0
    ! firstly, I need to get the value of n_useful_cells
     do i_basis_1=1, n_centers_basis

        temp = center_to_cell( Cbasis_to_center(i_basis_1) )
        do i=1, n_useful_cells
           if( temp.eq.cell_index_tmp(i) ) goto 1001 ! this cell_index is already contained
        enddo

        n_useful_cells = n_useful_cells + 1
        cell_index_tmp(n_useful_cells) = center_to_cell( Cbasis_to_center(i_basis_1) )

 1001 continue
     enddo

     allocate( usefulcell_to_fullcell(n_useful_cells), cbasis_to_useful_cell(n_centers_basis) )
     allocate( cellbasis_to_fullbasis(n,n_useful_cells) )
    ! then, I calculate usefulcell_to_fullcell and cbasis_to_useful_cell
     usefulcell_to_fullcell = 0; cbasis_to_useful_cell = 0; cellbasis_to_fullbasis = -100
     i_cell_1 = 0
     do i_basis_1=1, n_centers_basis

        temp = center_to_cell( Cbasis_to_center(i_basis_1) )
        do i=1, i_cell_1
           if( temp.eq.usefulcell_to_fullcell(i) )then ! this cell_index is already contained
               cbasis_to_useful_cell(i_basis_1) = i
               cellbasis_to_fullbasis(Cbasis_to_basis(i_basis_1),i) = i_basis_1
               goto 1002 
           endif
        enddo

        i_cell_1 = i_cell_1 + 1
        cbasis_to_useful_cell(i_basis_1) = i_cell_1
        usefulcell_to_fullcell(i_cell_1) = center_to_cell( Cbasis_to_center(i_basis_1) )
        cellbasis_to_fullbasis(Cbasis_to_basis(i_basis_1),i_cell_1) = i_basis_1

 1002 continue
     enddo

  else if (iop.eq.4 .or. iop.eq.5)then

     cell_index_tmp = 0
     n_useful_cells = 0
     do i_basis_1=1, n_centers_basis

        temp = center_to_cell( Cbasis_to_center_s(i_basis_1) )
        do i=1, n_useful_cells
           if( temp.eq.cell_index_tmp(i) ) goto 1003 ! this cell_index is already contained
        enddo

        n_useful_cells = n_useful_cells + 1
        cell_index_tmp(n_useful_cells) = center_to_cell( Cbasis_to_center_s(i_basis_1) )

 1003 continue
     enddo

     allocate( usefulcell_to_fullcell_s(n_useful_cells), cbasis_to_useful_cell_s(n_centers_basis) )
     allocate( cellbasis_to_fullbasis_s(n,n_useful_cells) )

     usefulcell_to_fullcell_s = 0; cbasis_to_useful_cell_s = 0; cellbasis_to_fullbasis_s = -100
     i_cell_1 = 0
     do i_basis_1=1, n_centers_basis

        temp = center_to_cell( Cbasis_to_center_s(i_basis_1) )
        do i=1, i_cell_1
           if( temp.eq.usefulcell_to_fullcell_s(i) )then ! this cell_index is already contained
               cbasis_to_useful_cell_s(i_basis_1) = i
               cellbasis_to_fullbasis_s(Cbasis_to_basis_s(i_basis_1),i) = i_basis_1
               goto 1004 
           endif
        enddo

        i_cell_1 = i_cell_1 + 1
        cbasis_to_useful_cell_s(i_basis_1) = i_cell_1
        usefulcell_to_fullcell_s(i_cell_1) = center_to_cell( Cbasis_to_center_s(i_basis_1) )
        cellbasis_to_fullbasis_s(Cbasis_to_basis_s(i_basis_1),i_cell_1) = i_basis_1

 1004 continue
     enddo

  endif

  deallocate(cell_index_tmp)
  
  if(iop.eq.1) then
     write(info_str,'(2X,A,I5)') "n_useful_cells=", n_useful_cells
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
  endif

  call aims_allocate( mat_cell, n, n, "mat_cell" )
  call aims_allocate( mat_sq, 4, dim_matrix_rel, dim_matrix_rel, "mat_sq" )

  do i_spin=1, n_spin

 !!! allocate ( ham_index(n_centers_basis,n_centers_basis) )
 !!! i_index = 0
 !!! do i_basis_2 = 1, n_centers_basis
 !!! do i_basis_1 = 1, i_basis_2
 !!!    i_index = i_index + 1
 !!!    ham_index(i_basis_1, i_basis_2) = i_index 
 !!! end do
 !!! end do

    ! First, expand the packed upper triangular matrix to square form:

 !!! call aims_allocate( mat_full,           n_centers_basis,n_centers_basis, "mat_full" )
 !!! call aims_allocate( mat_full_transpose, n_centers_basis,n_centers_basis, "mat_full_transpose" )
 !!! mat_full = 0.d0

 !!! do i_basis_1 = 1, n_centers_basis
 !!!    mat_full( 1:i_basis_1, i_basis_1 ) = mat_scalar( ham_index(1:i_basis_1,i_basis_1), i_spin )
 !!! end do

 !!! deallocate ( ham_index ) ! this could be a large array, deallocate it to save the memory

 !!! do i=1, n_centers_basis
 !!! do j=1, n_centers_basis
 !!!    mat_full_transpose(j,i)=mat_full(i,j)
 !!! enddo
 !!! enddo

 !!! mat_full = mat_full + mat_full_transpose

 !!! do i_basis_2 = 1, n_centers_basis
 !!!    mat_full(i_basis_2, i_basis_2) = mat_full(i_basis_2, i_basis_2)/2
 !!! end do

    ! Then, devide the full matrix to different cells:

                ! if(i_spin.eq.1)then
                !   write(6,*)'mat_full:'
                !   do i=1, n_centers_basis
                !      write(6,"(20(f13.7))")mat_full(1:n_centers_basis,i)
                !   enddo
                ! endif

    ! Transform scalar integrations to spinor integrations:

    ! (Rundong) The serial code here used to be: 
    ! do i_cell_2 = 1, n_useful_cells
    ! do i_cell_1 = 1, n_useful_cells
    ! For MPI parallelization, we switch it to a single loop like now:
     do i_cell = 1, n_useful_cells*n_useful_cells
     if( (myid == mod(i_cell,n_tasks) .and. myid <= n_useful_cells*n_useful_cells) ) then
        i_cell_2 = int((i_cell-1.d-5)/n_useful_cells)+1
        i_cell_1 = mod(i_cell,n_useful_cells)
        if(i_cell_1.eq.0)i_cell_1=n_useful_cells

        allocate ( A(dim_matrix_rel,dim_matrix_rel), B(dim_matrix_rel,dim_matrix_rel) )
        mat_cell = 0.d0; mat_sq = 0.d0

       ! large component part
        if(iop.eq.1 .or. iop.eq.2 .or. iop.eq.3)then
          do i_basis_2=1, n
             if(cellbasis_to_fullbasis(i_basis_2,i_cell_2).eq.-100)cycle
             do i_basis_1=1, n
                if(cellbasis_to_fullbasis(i_basis_1,i_cell_1).eq.-100)cycle
                row    = cellbasis_to_fullbasis(i_basis_1,i_cell_1)
                column = cellbasis_to_fullbasis(i_basis_2,i_cell_2)
                if(row.gt.column)then
                   temp = row*(row-1)/2
                  !mat_cell(i_basis_1,i_basis_2) = mat_full(row, column)
                   mat_cell(i_basis_1,i_basis_2) = mat_scalar( temp+column, i_spin )
                else
                   temp = column*(column-1)/2
                   mat_cell(i_basis_1,i_basis_2) = mat_scalar( temp+row, i_spin )
                endif
             enddo
          enddo
          do i_basis_1=1, n
             if(cellbasis_to_fullbasis(i_basis_1,i_cell_1).eq.-100) mat_cell(i_basis_1,:) = 0.d0
          enddo

          call trans_aomat_nr2r_fockmat_calc_sq &
               (s, n_basis, dim_matrix_rel, large_comp_idx, large_comp_clb, mat_cell, mat_sq)
        endif

       ! small component part
        if(iop.eq.4 .or. iop.eq.5)then
          do i_basis_2=1, n
             if(cellbasis_to_fullbasis_s(i_basis_2,i_cell_2).eq.-100)cycle
             do i_basis_1=1, n
                if(cellbasis_to_fullbasis_s(i_basis_1,i_cell_1).eq.-100)cycle
                row    = cellbasis_to_fullbasis_s(i_basis_1,i_cell_1)
                column = cellbasis_to_fullbasis_s(i_basis_2,i_cell_2)
                if(row.gt.column)then
                   temp = row*(row-1)/2
                  !mat_cell(i_basis_1,i_basis_2) = mat_full(row, column)
                   mat_cell(i_basis_1,i_basis_2) = mat_scalar( temp+column, i_spin )
                else
                   temp = column*(column-1)/2
                   mat_cell(i_basis_1,i_basis_2) = mat_scalar( temp+row, i_spin )
                endif
             enddo
          enddo
          do i_basis_1=1, n
             if(cellbasis_to_fullbasis_s(i_basis_1,i_cell_1).eq.-100) mat_cell(i_basis_1,:) = 0.d0
          enddo

          call trans_aomat_nr2r_fockmat_calc_sq &
               (s, n_basis_small, dim_matrix_rel, small_comp_idx, small_comp_clb, mat_cell, mat_sq)
        endif

        ! now, transform the data in mat_sq into mat_spinor for output
        do i=1, dim_matrix_rel
        do j=1, dim_matrix_rel
           A(j,i)=dcmplx( mat_sq(1,j,i),mat_sq(2,j,i) ) ! matrix A
           B(j,i)=dcmplx( mat_sq(3,j,i),mat_sq(4,j,i) ) ! matrix B
        enddo
        enddo

        ! Multiply k phase factor to the matrix to obtain the final matrix for diagonalization
        if(iop.eq.1 .or. iop.eq.2 .or. iop.eq.3)then
           i_full_cell_1 = usefulcell_to_fullcell(i_cell_1);   i_full_cell_2 = usefulcell_to_fullcell(i_cell_2)
        elseif(iop.eq.4 .or. iop.eq.5)then
           i_full_cell_1 = usefulcell_to_fullcell_s(i_cell_1); i_full_cell_2 = usefulcell_to_fullcell_s(i_cell_2)
        endif

        do i_k_point = 1, n_k_points

          temp_comp = k_phase(i_full_cell_1,i_k_point) * dconjg( k_phase(i_full_cell_2,i_k_point) )

          do i=1, dim_matrix_rel
          do j=1, dim_matrix_rel

             if( dabs(dble(A(j,i))).gt.thresh .or. dabs(aimag(A(j,i))).gt.thresh )then
                mat_spinor(j,i,1,i_k_point,i_spin) = mat_spinor(j,i,1,i_k_point,i_spin) + A(j,i) * temp_comp
                mat_spinor(j,i,3,i_k_point,i_spin) = mat_spinor(j,i,3,i_k_point,i_spin) + dconjg(A(j,i)) * temp_comp
             endif
             if( dabs(dble(B(j,i))).gt.thresh .or. dabs(aimag(B(j,i))).gt.thresh )then
                mat_spinor(j,i,2,i_k_point,i_spin) = mat_spinor(j,i,2,i_k_point,i_spin) + B(j,i) * temp_comp
             endif

          enddo
          enddo
        enddo ! end of i_k_point

        deallocate ( A,B )

     endif
     enddo ! end of i_cell

     do i_k_point =1 ,n_k_points
        do j=1, 3
           call sync_vector_complex(mat_spinor(1,1,j,i_k_point,i_spin),dim_matrix_rel*dim_matrix_rel)
        enddo
     enddo
                              if(i_spin.eq.1)then
                             !if(iop.eq.1)then
                               !write(6,*)'A(k=1):'
                               !do i=1, dim_matrix_rel
                               !   write(6,"(10(f11.5,f11.5,2x))")mat_spinor(1:dim_matrix_rel,i,1,1,i_spin)
                               !enddo
                               !write(6,*)'B(k=1):'
                               !do i=1, dim_matrix_rel
                               !   write(6,"(10(f11.5,f11.5,2x))")mat_spinor(1:dim_matrix_rel,i,2,1,i_spin)
                               !enddo
                               !write(6,*)'C(k=1):'
                               !do i=1, dim_matrix_rel
                               !   write(6,"(10(f11.5,f11.5,2x))")mat_spinor(1:dim_matrix_rel,i,3,1,i_spin)
                               !enddo
                             !endif
                              endif

 !!! call aims_deallocate( mat_full, "mat_full" )
 !!! call aims_deallocate( mat_full_transpose, "mat_full_transpose" )
  enddo ! end of n_spin

  if(iop.eq.1 .or. iop.eq.2 .or. iop.eq.3)then
     deallocate(usefulcell_to_fullcell,cbasis_to_useful_cell)
     deallocate( cellbasis_to_fullbasis )
  endif
  if(iop.eq.4 .or. iop.eq.5)then
     deallocate(usefulcell_to_fullcell_s,cbasis_to_useful_cell_s)
     deallocate( cellbasis_to_fullbasis_s )
  endif
  call aims_deallocate( mat_cell, "mat_cell" )
  call aims_deallocate( mat_sq,   "mat_sq" )

 end subroutine transform_scalar2spinor


 subroutine trans_aomat_nr2r_fockmat_calc_sq(s,n1,n2,idx,cgs,fnr,ab)
 !> when transform fock matrix from 1c picture to 2c or 4c picture
 !  the result relativistic Fock matrix will be in the format:
 !   ( A   B  )
 !   (-B^* A^*)
 !   \f$A_{uv}= <\chi_v^{\dag}|F|\chi_u>\f$
 !   \f$B_{uv}= <\chi_v^{\dag}|F|\bar{\chi}_u>\f$
 !  here, we get A and B matrix
 !  here, for solid calculation, A and B not have any symmetry
  implicit none
  integer,intent(in)    :: n2,n1,s,idx(4,n2)  ! Nr. of R-AO and NR-AO basis, spins; index link them
  real*8, intent(in)    :: cgs(4,n2),fnr(s,n1,n1) ! CG coef; Matrix in NR AO picture
  real*8, intent(inout) :: ab(4,s,*)  ! A (1:2,:,:) and B(3:4,:,:) matrix, square form
  real*8  :: a1(s),a2(s),b1(s),b2(s)  ! s is a very small number (<10), can use auto array
  integer :: u,v,m,k,l
  real*8 :: thresh=1.d-8
  l = 0
  do u=1,n2
  do v=1,n2
    a1=0.d0  ! Real and Imag part of A and B
    a2=0.d0
    b1=0.d0
    b2=0.d0
    if(idx(1,u).gt.0)then
      if( idx(1,v).gt.0 )then
          if( dabs(fnr(1,idx(1,v),idx(1,u))).gt.thresh ) a1(:)=a1(:)+fnr(:,idx(1,v),idx(1,u))*cgs(1,v)*cgs(1,u)
      endif
      if( idx(2,v).gt.0 )then
          if( dabs(fnr(1,idx(2,v),idx(1,u))).gt.thresh ) b1(:)=b1(:)+fnr(:,idx(2,v),idx(1,u))*cgs(2,v)*cgs(1,u)
      endif
      if( idx(3,v).gt.0 )then
          if( dabs(fnr(1,idx(3,v),idx(1,u))).gt.thresh ) a2(:)=a2(:)-fnr(:,idx(3,v),idx(1,u))*cgs(3,v)*cgs(1,u)
      endif
      if( idx(4,v).gt.0 )then
          if( dabs(fnr(1,idx(4,v),idx(1,u))).gt.thresh ) b2(:)=b2(:)-fnr(:,idx(4,v),idx(1,u))*cgs(4,v)*cgs(1,u)
      endif
    endif
    if(idx(2,u).gt.0)then
      if( idx(1,v).gt.0 )then
          if( dabs(fnr(1,idx(1,v),idx(2,u))).gt.thresh ) b1(:)=b1(:)-fnr(:,idx(1,v),idx(2,u))*cgs(1,v)*cgs(2,u)
      endif
      if( idx(2,v).gt.0 )then
          if( dabs(fnr(1,idx(2,v),idx(2,u))).gt.thresh ) a1(:)=a1(:)+fnr(:,idx(2,v),idx(2,u))*cgs(2,v)*cgs(2,u)
      endif
      if( idx(3,v).gt.0 )then
          if( dabs(fnr(1,idx(3,v),idx(2,u))).gt.thresh ) b2(:)=b2(:)+fnr(:,idx(3,v),idx(2,u))*cgs(3,v)*cgs(2,u)
      endif
      if( idx(4,v).gt.0 )then
          if( dabs(fnr(1,idx(4,v),idx(2,u))).gt.thresh ) a2(:)=a2(:)-fnr(:,idx(4,v),idx(2,u))*cgs(4,v)*cgs(2,u)
      endif
    endif
    if(idx(3,u).gt.0)then
      if( idx(1,v).gt.0 )then 
          if( dabs(fnr(1,idx(1,v),idx(3,u))).gt.thresh ) a2(:)=a2(:)+fnr(:,idx(1,v),idx(3,u))*cgs(1,v)*cgs(3,u)
      endif
      if( idx(2,v).gt.0 )then 
          if( dabs(fnr(1,idx(2,v),idx(3,u))).gt.thresh ) b2(:)=b2(:)-fnr(:,idx(2,v),idx(3,u))*cgs(2,v)*cgs(3,u)
      endif
      if( idx(3,v).gt.0 )then 
          if( dabs(fnr(1,idx(3,v),idx(3,u))).gt.thresh ) a1(:)=a1(:)+fnr(:,idx(3,v),idx(3,u))*cgs(3,v)*cgs(3,u)
      endif
      if( idx(4,v).gt.0 )then 
          if( dabs(fnr(1,idx(4,v),idx(3,u))).gt.thresh ) b1(:)=b1(:)-fnr(:,idx(4,v),idx(3,u))*cgs(4,v)*cgs(3,u)
      endif
    endif
    if(idx(4,u).gt.0)then
      if( idx(1,v).gt.0 )then 
          if ( dabs(fnr(1,idx(1,v),idx(4,u))).gt.thresh ) b2(:)=b2(:)+fnr(:,idx(1,v),idx(4,u))*cgs(1,v)*cgs(4,u)
      endif
      if( idx(2,v).gt.0 )then 
          if ( dabs(fnr(1,idx(2,v),idx(4,u))).gt.thresh ) a2(:)=a2(:)+fnr(:,idx(2,v),idx(4,u))*cgs(2,v)*cgs(4,u)
      endif
      if( idx(3,v).gt.0 )then 
          if ( dabs(fnr(1,idx(3,v),idx(4,u))).gt.thresh ) b1(:)=b1(:)+fnr(:,idx(3,v),idx(4,u))*cgs(3,v)*cgs(4,u)
      endif
      if( idx(4,v).gt.0 )then 
          if ( dabs(fnr(1,idx(4,v),idx(4,u))).gt.thresh ) a1(:)=a1(:)+fnr(:,idx(4,v),idx(4,u))*cgs(4,v)*cgs(4,u)
      endif
    endif
    l = l+1
    ab(1,:,l) = ab(1,:,l)+a1(:)  ! real part
    ab(2,:,l) = ab(2,:,l)+a2(:)  ! imag part
    ab(3,:,l) = ab(3,:,l)+b1(:)
    ab(4,:,l) = ab(4,:,l)+b2(:)
  enddo
  enddo
 end subroutine trans_aomat_nr2r_fockmat_calc_sq

 subroutine trans_aomat_r2nr_denvu(s,n,idxu,cgsu,idxv,cgsv,ab,d)
  implicit none
  integer,  intent(in)    :: idxu(4),idxv(4),s,n
  real*8,   intent(in)    :: cgsu(4),cgsv(4),ab(s,4)
  real*8,   intent(inout) :: d(s,n,n)
  if(idxu(1).gt.0)then
    if(idxv(1).gt.0) d(:,idxv(1),idxu(1))=d(:,idxv(1),idxu(1))+ab(:,1)*cgsv(1)*cgsu(1)
    if(idxv(2).gt.0) d(:,idxv(2),idxu(1))=d(:,idxv(2),idxu(1))-ab(:,3)*cgsv(2)*cgsu(1)
    if(idxv(3).gt.0) d(:,idxv(3),idxu(1))=d(:,idxv(3),idxu(1))+ab(:,2)*cgsv(3)*cgsu(1)
    if(idxv(4).gt.0) d(:,idxv(4),idxu(1))=d(:,idxv(4),idxu(1))+ab(:,4)*cgsv(4)*cgsu(1)
  endif
  if(idxu(2).gt.0)then
    if(idxv(1).gt.0) d(:,idxv(1),idxu(2))=d(:,idxv(1),idxu(2))+ab(:,3)*cgsv(1)*cgsu(2)
    if(idxv(2).gt.0) d(:,idxv(2),idxu(2))=d(:,idxv(2),idxu(2))+ab(:,1)*cgsv(2)*cgsu(2)
    if(idxv(3).gt.0) d(:,idxv(3),idxu(2))=d(:,idxv(3),idxu(2))-ab(:,4)*cgsv(3)*cgsu(2)
    if(idxv(4).gt.0) d(:,idxv(4),idxu(2))=d(:,idxv(4),idxu(2))+ab(:,2)*cgsv(4)*cgsu(2)
  endif
  if(idxu(3).gt.0)then
    if(idxv(1).gt.0) d(:,idxv(1),idxu(3))=d(:,idxv(1),idxu(3))-ab(:,2)*cgsv(1)*cgsu(3)
    if(idxv(2).gt.0) d(:,idxv(2),idxu(3))=d(:,idxv(2),idxu(3))+ab(:,4)*cgsv(2)*cgsu(3)
    if(idxv(3).gt.0) d(:,idxv(3),idxu(3))=d(:,idxv(3),idxu(3))+ab(:,1)*cgsv(3)*cgsu(3)
    if(idxv(4).gt.0) d(:,idxv(4),idxu(3))=d(:,idxv(4),idxu(3))+ab(:,3)*cgsv(4)*cgsu(3)
  endif
  if(idxu(4).gt.0)then
    if(idxv(1).gt.0) d(:,idxv(1),idxu(4))=d(:,idxv(1),idxu(4))-ab(:,4)*cgsv(1)*cgsu(4)
    if(idxv(2).gt.0) d(:,idxv(2),idxu(4))=d(:,idxv(2),idxu(4))-ab(:,2)*cgsv(2)*cgsu(4)
    if(idxv(3).gt.0) d(:,idxv(3),idxu(4))=d(:,idxv(3),idxu(4))-ab(:,3)*cgsv(3)*cgsu(4)
    if(idxv(4).gt.0) d(:,idxv(4),idxu(4))=d(:,idxv(4),idxu(4))+ab(:,1)*cgsv(4)*cgsu(4)
  endif
 end subroutine trans_aomat_r2nr_denvu



 subroutine relativistic_get_clbidx_oneorb(do_small_comp,qq,ll,jj,mms,mva,idx,cclb)
 ! get CG coefficients for one orbital, large or small component
  implicit none
  logical,intent(in)    :: do_small_comp
  integer,intent(in)    :: ll,jj,mms,mva,qq
  integer,intent(out)   :: idx(4)
  real*8, intent(out)   :: cclb(4)
  real*8, parameter :: d2=0.70710678118654752440d0 !< \f$\sqrt(2)/2\f$
  real*8, external  :: clb
  integer :: mm1,qqq,mm,m1
  mm1 = mms+1
  mm  = iabs(mms) ! m.alpha = mj-1/2
  m1  = iabs(mm1) ! m.beta  = mj+1/2
  qqq = qq+ll+1
  idx = 0
  cclb(1)=clb(ll,jj,mva, 1)
  cclb(2)=clb(ll,jj,mva,-1)
  if(mm.ne.0) cclb(1)=cclb(1)*d2
  if(m1.ne.0) cclb(2)=cclb(2)*d2
  ! add the factor (-1)^{ (|m|+m)/2 } under Condon-Shortley convention ....
  if(mms.gt.0.and.mod(mm,2).eq.1) cclb(1)=-cclb(1)
  if(mm1.gt.0.and.mod(m1,2).eq.1) cclb(2)=-cclb(2)
  cclb(3:4) = cclb(1:2)
  if(do_small_comp)then ! small component
    if(mms.ge.0) cclb(1)=-cclb(1)
    if(mm1.ge.0) cclb(2)=-cclb(2)
    if(mm.le.ll) idx(3) = qqq+mm
    if(m1.le.ll) idx(4) = qqq+m1
    if(mm.ne.0.and.mm.le.ll) idx(1) = qqq-mm
    if(m1.ne.0.and.m1.le.ll) idx(2) = qqq-m1
  else                  ! large component
    if(mms.lt.0) cclb(3)=-cclb(3)
    if(mm1.lt.0) cclb(4)=-cclb(4)
    if(mm.le.ll) idx(1) = qqq+mm
    if(m1.le.ll) idx(2) = qqq+m1
    if(mm.ne.0.and.mm.le.ll) idx(3) = qqq-mm
    if(m1.ne.0.and.m1.le.ll) idx(4) = qqq-m1
  endif
 end subroutine relativistic_get_clbidx_oneorb
!...clb1...vector coupling coeffic....................
!... clebsh-gordon coeffic....................
      function clb(l,j,mj,is)
      use localorb_io, only: use_unit
      implicit double precision(a-h,o-z)
      data d1,d2,d5/1.0d0,2.0d0,0.5d0/
      pl=l
      u=mj*d5
      if(j.gt.2*l.and.is.gt.0) go to 1001
      if(j.gt.2*l.and.is.lt.0) go to 1002
      if(j.lt.2*l.and.is.gt.0) go to 1003
      if(j.lt.2*l.and.is.lt.0) go to 1004
      write(use_unit,2000)
      write(use_unit,"(a)")'clb'
      stop
 1001 clb= dsqrt((pl+u+d5)/(d2*pl+d1))
      return
 1002 clb= dsqrt((pl-u+d5)/(d2*pl+d1))
      return
 1003 clb= -dsqrt((pl-u+d5)/(d2*pl+1))
      return
 1004 clb= dsqrt((pl+u+d5)/(d2*pl+d1))
      return
 2000 format('  the combination of j and mj is not satisfied in clb and caused stop')
      end function clb




 subroutine transform_spinor2scalar_denmat(denmat_spinor,kdensmat_complex)
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  use dimensions, only : n_basis
  use rel_x2c_mod
  implicit none

  complex*16, intent(in)  :: denmat_spinor(2*dim_matrix_rel,2*dim_matrix_rel)
  complex*16, intent(out) :: kdensmat_complex(n_basis,n_basis)

  integer :: naos(2)
  integer :: i,j,k
  real*8,allocatable :: denmat_scalar(:,:)

  call aims_allocate( denmat_scalar, n_basis, n_basis, "denmat_scalar" )
  naos(1)=n_basis; naos(2)=n_basis_small

  call trans_aomat_r2nr_denmat(201,1, dim_matrix_rel, naos, large_comp_idx,large_comp_clb, small_comp_idx,small_comp_clb, &
       denmat_spinor, denmat_scalar)

  do i=1, n_basis
  do j=1, n_basis
     kdensmat_complex(j,i) = dcmplx(denmat_scalar(j,i),0.d0)
  enddo
  enddo

  if(allocated(denmat_scalar)) call aims_deallocate( denmat_scalar, "denmat_scalar" )

 end subroutine transform_spinor2scalar_denmat


 subroutine trans_aomat_r2nr_denmat(iop0,s,nao,naos,lidx,lcgs,sidx,scgs,dr,dnr)
 ! input is of form  ( A  B ) density matrix, CnC^{\dag}
 !                   ( E  C )
 ! note, here use the symmetry:
 ! <u|v> = <\bar{u}|\bar{v}>^*, <\bar{u}|v> = - <u|\bar{v}>^*
 ! and the notation:
 !  D_{vu} = A_{vu}, D_{\bar{v}\bar{u}} = C_{vu},
 !  D_{\bar{v}u} = E_{vu}, D_{v\bar{u}} = B_{vu}
 !  u|v = u^{\dag}v, u and v are spinor, i and j are scalar basis
 !
 ! FOR ELECTRON DENSITY CALCULATION: (note, electron density must be REAL, not complex)
 !   rho = Re (\sum_{vu}  (  D_{vu} u|v + D_{\bar{v}\bar{u}} \bar{u}|\bar{v} 
 !                         + D_{\bar{v}u} u|\bar{v} + D_{v\bar{u}} \bar{u}|v  ) )
 !  use symmetry, one can get (equation 1)
 !   rho = Re (\sum_{vu} (  (A_{vu}+C_{vu}^*) u|v + (B_{vu}-E_{vu}^*) \bar{u}|v ) )
 !       = \sum_{ij} p_{ij} i*j
 !
 ! and when the input density matrix and basis overlap have the symmetry form: 
 !  <u|v> = <v|u>^*, A_{uv} = A_{vu}^{*}, C_{uv} = C_{vu}^{*}, E_{uv} = B_{vu}^* 
 !
 ! then (equation 2)
 !   rho = \sum_{v<u} (  D_{vu} u|v + D_{uv} v|u  
 !                     + D_{\bar{u}\bar{v}} \bar{v}|\bar{u}
 !                     + D_{\bar{v}\bar{u}} \bar{u}|\bar{v}
 !                     + D_{\bar{v}u} u|\bar{v} + D_{u\bar{v}} \bar{v}|u
 !                     + D_{\bar{u}v} v|\bar{u} + D_{v\bar{u}} \bar{u}|v )
 !         + \sum_{u} D_{uu} u|u + D_{\bar{u}\bar{u}} \bar{u}|\bar{u}
 !       = Re( \sum_{u} (A_{uu}+C_{uu})u|u 
 !           + \sum_{v<u} [ (A_{vu}+A_{uv}^*+C_{uv}+C_{vu}^*)u|v
 !                        + (B_{vu}+E_{uv}^*-B_{uv}-E_{vu}^*)\bar{u}|v ]
 !       = \sum_{ij} p_{ij} i*j
 !
 !  INPUT control parameter iop0 == iopt*100+iop
 !   iop : = 1/2 only get large/small component density matrix in scalar basis
 !         else both
 !   iopt: = 0/1 use equation 2/1, and the output density matrix will be symmetrize
 !           i.e., dij(out) = (p_{ij}+p_{ji})/2
 !         else, use equation 1 and do not symmetrize
  implicit none
  integer,intent(in)  :: s,iop0       ! Nr. of spins; iop=1/2 only large/small comp., else both
  integer,intent(in)  :: nao,naos(2)  ! Dim of A and B matrix (Nr. half R-AO), NR-AOs(L/S comp.)
  integer,intent(in)  :: lidx(*),sidx(*)  ! index   link NR-AOs and R-AOs
  real*8, intent(in)  :: lcgs(*),scgs(*)  ! CG coef link NR-AOs and R-AOs
  complex*16, intent(in)  :: dr(s,nao+nao,nao+nao)
  real*8,     intent(out) :: dnr(*)
  real*8,     allocatable :: ab(:,:,:)
  complex*16  :: aa(s)
  real*8  :: ss
  integer :: i,j,k,iop,iopt

  if(nao.le.0) return
  if(naos(1).le.0.and.naos(2).le.0) return
  iop  = mod(iop0,100)
  iopt = iabs(iop0)/100
  if(iop.eq.2.and.naos(2).le.0) return
  if(iop.eq.1.and.naos(1).le.0) return

  if( iopt.eq.0 )then
    allocate( ab(s,4,nao*(nao+1)/2) )
    k = 0
    do i=1,nao
      do j=1,i-1
        k = k+1
        aa= dr(:,i,j)+dr(:,j+nao,i+nao)+dconjg(dr(:,j,i)+dr(:,i+nao,j+nao))
        ab(:,1,k) = dreal( aa )
        ab(:,2,k) = dimag( aa )
        aa= dr(:,i,j+nao)-dr(:,j,i+nao)+dconjg( dr(:,j+nao,i)-dr(:,i+nao,j) )
        ab(:,3,k) = dreal( aa )
        ab(:,4,k) = dimag( aa )
      enddo
      k = k+1
      ab(:,  1,k) = dreal( dr(:,i,i) ) + dreal( dr(:,i+nao,i+nao) )
      ab(:,2:4,k) = 0.d0
    enddo
  else
    allocate( ab(s,4,nao*nao) )
    k = 0
    do i=1,nao
    do j=1,nao
      k = k+1
      aa= dr(:,i,j)+dconjg(dr(:,i+nao,j+nao))
      ab(:,1,k) = dreal( aa )
      ab(:,2,k) = dimag( aa )
      aa= dr(:,i,j+nao)-dconjg( dr(:,i+nao,j) )
      ab(:,3,k) = dreal( aa )
      ab(:,4,k) = dimag( aa )
    enddo
    enddo
  endif

  k = 1
 ! large component part matrix
  if(iop.ne.2.and.naos(1).gt.0)then
    call trans_aomat_r2nr_denmat_calc(iopt,s,naos(1),nao,lidx,lcgs,ab,dnr)
    k = naos(1)*naos(1)*s+1
  endif

 ! small component part matrix
  if(iop.ne.1.and.naos(2).gt.0)then
    call trans_aomat_r2nr_denmat_calc(iopt,s,naos(2),nao,sidx,scgs,ab,dnr(k))
  endif

  deallocate( ab )
 end subroutine trans_aomat_r2nr_denmat

 subroutine trans_aomat_r2nr_denmat_calc(iop,s,n1,n2,idx,cgs,ab,d)
  implicit none
  integer,intent(in)  :: iop
  integer,intent(in)  :: n2,n1,s,idx(4,n2)  ! Nr. of R-AO and NR-AO basis, spins; index link them
  real*8, intent(in)  :: cgs(4,n2)          ! CG coef
  real*8, intent(in)  :: ab(s,4,*)          ! A and B matrix
  real*8, intent(out) :: d(s,n1,n1)         ! Matrix in NR AO picture
  real*8  :: a1(s),a2(s),b1(s),b2(s)
  integer :: u,v,m,k

  d = 0.d0
  k = 0
  if(iop.eq.0)then
    do u=1,n2
      do v=1,u-1
        k = k+1
        call trans_aomat_r2nr_denvu(s,n1,idx(:,u),cgs(:,u),idx(:,v),cgs(:,v),ab(:,:,k),d)
      enddo
      k = k+1
      do m=1,4
        if(idx(m,u).gt.0) d(:,idx(m,u),idx(m,u))=d(:,idx(m,u),idx(m,u)) &
          +ab(:,1,k)*cgs(m,u)*cgs(m,u)
      enddo
    enddo
  else
    do u=1,n2
    do v=1,n2
      k = k+1
      call trans_aomat_r2nr_denvu(s,n1,idx(:,u),cgs(:,u),idx(:,v),cgs(:,v),ab(:,:,k),d)
    enddo
    enddo
  endif
  if( iop.le.1 )then
    do u = 2,n1
    do v = 1,u-1
      d(:,v,u) = 0.5d0*(d(:,v,u)+d(:,u,v))
      d(:,u,v) = d(:,v,u)
    enddo
    enddo
  endif
 end subroutine trans_aomat_r2nr_denmat_calc


