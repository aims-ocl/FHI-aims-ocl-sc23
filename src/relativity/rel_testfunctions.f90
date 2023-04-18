
 subroutine test_construct_mat_sum_packed( n_spin, n_centers_basis, ld_dirac_scalar, mat_scalar, mat_spinor)
  use pbc_lists
  use dimensions,  only: n_basis, n_k_points, n_centers_basis_I
  use localorb_io, only: use_unit
  implicit none
  integer,intent(in) :: n_spin, n_centers_basis, ld_dirac_scalar
  real*8,intent(in) :: mat_scalar(ld_dirac_scalar,n_spin)  ! the input scalar integration
  complex*16,intent(out) :: mat_spinor(n_basis,n_basis,n_k_points,n_spin) ! the output spinor integration

  integer :: n_useful_cells
  integer :: i_spin, i_index, i_k_point, i_basis_1, i_basis_2, i_cell_1, i_cell_2, i_full_cell_1, i_full_cell_2, n, i,j
  integer,allocatable :: ham_index(:,:)
  real*8,allocatable :: mat_cell(:,:,:,:) ! a temporary square matrix, mat_cell(n_basis,n_basis,n_useful_cells,n_useful_cells)
  real*8,allocatable :: mat_full(:,:) ! a temporary square matrix, mat_full(n_centers_basis_I,n_centers_basis_I)
  real*8,allocatable :: mat_full_transpose(:,:)
  integer,allocatable :: full_cell_index(:) ! links n_centers_basis to n_cells
  integer,allocatable :: cbasis_to_cell(:)
  integer,allocatable :: cell_index_tmp(:)
  complex*16 :: temp

  allocate( cell_index_tmp(n_cells) )

  cell_index_tmp = 0
  n_useful_cells = 0
  do i_basis_1=1, n_centers_basis

     temp = center_to_cell( Cbasis_to_center(i_basis_1) )
     do i=1, n_useful_cells
        if( temp.eq.cell_index_tmp(i) ) goto 2001 ! this cell_index is already contained
     enddo

     n_useful_cells = n_useful_cells + 1
     cell_index_tmp(n_useful_cells) = center_to_cell( Cbasis_to_center(i_basis_1) )

 2001 continue
  enddo

  allocate( full_cell_index(n_useful_cells), cbasis_to_cell(n_centers_basis) )

  full_cell_index = 0; cbasis_to_cell = 0
  i_cell_1 = 0
  do i_basis_1=1, n_centers_basis

     temp = center_to_cell( Cbasis_to_center(i_basis_1) )
     do i=1, i_cell_1
        if( temp.eq.full_cell_index(i) )then ! this cell_index is already contained
            cbasis_to_cell(i_basis_1) = i
            goto 2002 
        endif
     enddo

     i_cell_1 = i_cell_1 + 1
     cbasis_to_cell(i_basis_1) = i_cell_1
     full_cell_index(i_cell_1) = center_to_cell( Cbasis_to_center(i_basis_1) )

 2002 continue
  enddo

  deallocate(cell_index_tmp)


  n=n_basis
  allocate ( mat_cell(n,n,n_useful_cells,n_useful_cells), ham_index(n_centers_basis,n_centers_basis) )
  mat_spinor = (0.d0,0.d0)

  mat_cell = 0.d0

  i_index = 0
  do i_basis_2 = 1, n_centers_basis
  do i_basis_1 = 1, i_basis_2
     i_index = i_index + 1
     ham_index(i_basis_1, i_basis_2) = i_index 
  end do
  end do

  do i_spin=1, n_spin

    ! First, expand the packed upper triangular matrix to square form:

     allocate ( mat_full(n_centers_basis,n_centers_basis), mat_full_transpose(n_centers_basis,n_centers_basis) )
     mat_full = 0.d0

     do i_basis_1 = 1, n_centers_basis
        mat_full( 1:i_basis_1, i_basis_1 ) = mat_scalar( ham_index(1:i_basis_1,i_basis_1), i_spin )
     end do

     do i=1, n_centers_basis
     do j=1, n_centers_basis
        mat_full_transpose(j,i)=mat_full(i,j)
     enddo
     enddo

     mat_full = mat_full + mat_full_transpose

     do i_basis_2 = 1, n_centers_basis
        mat_full(i_basis_2, i_basis_2) = mat_full(i_basis_2, i_basis_2)/2
     end do


     do i_basis_2=1, n_centers_basis

        i_cell_2 = cbasis_to_cell(i_basis_2)

        do i_basis_1=1, n_centers_basis

           i_cell_1 = cbasis_to_cell(i_basis_1)

           mat_cell( Cbasis_to_basis(i_basis_1), Cbasis_to_basis(i_basis_2), i_cell_1, i_cell_2 ) = mat_full(i_basis_1, i_basis_2)
        enddo
     enddo

     deallocate ( mat_full,mat_full_transpose)

    ! Transform scalar integrations to spinor integrations:

     do i_cell_2 = 1, n_useful_cells
     do i_cell_1 = 1, n_useful_cells

        i_full_cell_1 = full_cell_index(i_cell_1);   i_full_cell_2 = full_cell_index(i_cell_2)

        ! Multiply k phase factor to the matrix to obtain the final matrix for diagonalization

        do i_k_point = 1, n_k_points
          do i=1, n_basis
          do j=1, n_basis

             temp = mat_cell(j,i,i_cell_1,i_cell_2) * k_phase(i_full_cell_1,i_k_point) * dconjg( k_phase(i_full_cell_2,i_k_point) )
             mat_spinor(j,i,i_k_point,i_spin) = mat_spinor(j,i,i_k_point,i_spin) + temp

          enddo
          enddo
        enddo ! end of i_k_point

     enddo ! end of cell_2
     enddo ! end of cell_1

  enddo

 end subroutine test_construct_mat_sum_packed

 subroutine density_matrix_con_test01(n_compute,i_basis,dm)
  use dimensions, only: n_centers_basis_I
  implicit none
  integer :: n_compute, i_basis(n_compute)
  real*8 :: dm(n_compute,n_compute)

  integer :: i,j
  real*8 :: dm_full(n_centers_basis_I,n_centers_basis_I)

  dm_full = 0.d0
  do i=1, n_compute
  do j=1, n_compute
     dm_full(i_basis(j),i_basis(i)) = dm(j,i)
  enddo
  enddo

  do i=1, 20
    write(6,"(20f13.7)") dm(1:20,i)
  enddo

 end subroutine density_matrix_con_test01

 subroutine writecc(n,c)
  implicit none
  integer :: n
  real*8 :: c(4*n,4*n)
  integer i,j
  do i=1, 4*n
     write(6,"(24f11.5)")c(:,i)
     write(6,"(a)")'    '
  enddo
 end subroutine writecc

 subroutine test_print_modrho(np,r,rho,pot)
  implicit none
  integer :: np
  real*8 :: r(np), rho(np), pot(np)

  integer :: i,j,k, np_use
  real*8 :: r_use(np), rho_use(np), pot_use(np)
  real*8 :: temp

  do i=1, np
  do j=i+1, np
     if(r(j).lt.r(i))then
        temp = r(j)
        r(j) = r(i)
        r(i) = temp
        temp = rho(j)
        rho(j) = rho(i)
        rho(i) = temp
        temp = pot(j)
        pot(j) = pot(i)
        pot(i) = temp
     endif
  enddo
  enddo

  np_use=1
  r_use(1) = r(1)
  rho_use(1) = rho(1)
  pot_use(1) = pot(1)
  do i=2, np
     do j=1, i-1
        if(dabs(r(i)-r(j)).lt.1.d-10)then ! this point has alreay been counted
          goto 1001
        endif
     enddo
     np_use = np_use + 1
     r_use(np_use) = r(i)
     rho_use(np_use) = rho(i)
     pot_use(np_use) = pot(i)
1001 continue
  enddo

  write(6,"('np_use=',i5,3x,'np=',i5)")np_use,np
  do i=1, np_use
     write(6,"('i=',i5,3x,'r=',f13.7,3x,'rho=',f16.7,3x,'hartree_pot=',f16.7)")i,r_use(i),rho_use(i),pot_use(i)
  enddo

 end subroutine test_print_modrho




