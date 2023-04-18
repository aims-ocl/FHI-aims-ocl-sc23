!***h* FHI-aims
!  NAME
!    tetrahedron_integration
!  SYNOPSIS
module tetrahedron_integration
  !  PURPOSE
  !    Perform tetrahedron integration for DOS and spectral functions. 
  !    Designed to be as independent as possible.
  !    ref.
  !    Extensions of the tetrahedron method for evaluating spectral properties of solids
  !    A H MacDonald etc. J. Phys. C: Solid State Phys., Vol. 12, 2991 (1979)
  !    DOI: 10.1088/0022-3719/12/15/008
  !  USES
use mpi_tasks
!use mpi_utilities
!use synchronize_mpi
implicit none
  !  AUTHOR
  !    Yi Yao
  !  HISTORY
  !    April 2019 created
  !  COPYRIGHT
  !  TODO
  !  SOURCE
  public :: ltidos, ltispectral
  private :: ltidos_single, ltispectral_single, matinv3, det3, sort

contains

  subroutine ltidos &
    ( cell, kmesh_n1, kmesh_n2, kmesh_n3, nbands, &
      eigs, nenergies, energies, dos_output, comm )
    !  PURPOSE
    !    linear tetrahedron integration method for density of states
    !  USES
    implicit none
    !  ARGUMENTS
    real*8,      intent(in) :: cell(3,3)
    integer,     intent(in) :: kmesh_n1, kmesh_n2, kmesh_n3, nbands
    real*8,      intent(in) :: eigs(nbands,kmesh_n1,kmesh_n2,kmesh_n3)
    integer,     intent(in) :: nenergies
    real*8,      intent(in) :: energies(nenergies)
    real*8,      intent(out) :: dos_output(nenergies)
    integer,     intent(in) :: comm
    !  INPUTS
    !    o cell - 3d cell vector
    !    o kmesh_n1, kmesh_n2, kmesh_n3 - number of k points in 
    !      each direction
    !    o nbands - number of states in the calculation
    !    o eigs - eigenvalues ( KS or GW or any other states )
    !    o nenergies - number of points for energies for output
    !    o energies - list of energies for output
    !  OUTPUT
    !    o dos_output - list of density of states for output
    !  AUTHOR
    !    Yi Yao
    !  HISTORY
    !    Created in April 2019
    !  SOURCE
    !  local variables
    integer  :: indices(3,8)
    integer  :: kmesh_size(3)
    real*8   :: kpts_submesh(3,8)
    real*8   :: kpts(3,4)
    real*8   :: invcell(3,3), B(3,3)
    integer  :: main_diagonal(2,4)
    integer  :: remain_points(6,4)
    integer  :: list_i_r(6)
    real*8   :: main_diagonal_lengths(4)
    integer  :: index_diagonal_shortest
    integer  :: i_a, i_b
    integer  :: tetra_indecies(4,6)
    integer  :: s(4)
    real*8   :: volume, detkpts
    real*8   :: det_invcell
    real*8   :: E(nbands,4)
    integer  :: my_rank, mpi_size, mpierr


    !  counters
    integer  :: ii, jj
    integer  :: i_k, j_k, k_k
    integer  :: ii_k, jj_k, kk_k
    integer  :: n
    

    !

    call mpi_comm_rank(comm,my_rank,mpierr)
    call mpi_comm_size(comm,mpi_size,mpierr)

    !
    kmesh_size = (/kmesh_n1,kmesh_n2,kmesh_n3/)
    call matinv3(cell,invcell)
    B = invcell / spread(kmesh_size,1,3)
    indices = reshape((/0,0,0,&
               0,0,1,&
               0,1,0,&
               0,1,1,&
               1,0,0,&
               1,0,1,&
               1,1,0,&
               1,1,1/),shape(indices))
    
    do ii = 1,8
      do jj = 1,3
        kpts_submesh(jj,ii) =  dot_product(indices(:,ii),B(:,jj))
      end do
    end do
    main_diagonal = reshape((/1,8,2,7,4,5,3,6/),shape(main_diagonal))
    remain_points = reshape((/2,4,3,7,5,6,&
                     1,3,4,8,6,5,&
                     1,2,6,8,7,3,&
                     1,2,4,8,7,5/),shape(remain_points))
    do ii = 1,4
      main_diagonal_lengths(ii) = sqrt(sum((kpts_submesh(:,main_diagonal(1,ii)) &
                                     - kpts_submesh(:,main_diagonal(2,ii)))**2))
    end do
    index_diagonal_shortest = minloc(main_diagonal_lengths,1)
    i_a = main_diagonal(1,index_diagonal_shortest)
    i_b = main_diagonal(2,index_diagonal_shortest)
    list_i_r = remain_points(:,index_diagonal_shortest)
    do ii = 1,6
      tetra_indecies(1,ii) = i_a
      tetra_indecies(2,ii) = i_b
      tetra_indecies(3,ii) = list_i_r(ii)
      tetra_indecies(4,ii) = list_i_r(mod(ii,6)+1)
    end do

    dos_output = 0.0

    do ii = 1,6
      s = tetra_indecies(:,ii)
      do jj = 1,4
        kpts(:,jj) = kpts_submesh(:,s(jj))
      end do
      call det3(kpts(:,1:3) - kpts(:,2:4), detkpts)
      volume = abs(detkpts / 6.0)
      n = 0
      do i_k = 1,kmesh_n1
        do j_k = 1,kmesh_n2
          do k_k = 1,kmesh_n3
            n = n + 1
            if (mod(n, mpi_size) .ne. my_rank) then
              cycle
            end if
            !if n % world.size .ne. world.rank then
            !  continue
            !end if
            do jj = 1,4
              ii_k = mod((i_k + indices(1,s(jj)) - 1) , kmesh_n1) + 1
              jj_k = mod((j_k + indices(2,s(jj)) - 1) , kmesh_n2) + 1
              kk_k = mod((k_k + indices(3,s(jj)) - 1) , kmesh_n3) + 1
              E(:,jj) = eigs(:,ii_k,jj_k,kk_k)
            end do
            call ltidos_single(nenergies, energies, dos_output, volume, nbands, E)
          end do
        end do
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE, dos_output, nenergies, MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpierr)

    call det3(invcell, det_invcell)
    dos_output = dos_output / abs(det_invcell)
  end subroutine ltidos

  subroutine ltidos_single(nenergies, energies, dos, volume, nbands, E)
    !  PURPOSE
    !    linear tetrahedron integration method for density of states
    !    helper subroutine for single tetrahedron
    !  USES
    implicit none
    !  ARGUMENTS
    integer,     intent(in) :: nenergies
    real*8,      intent(in) :: energies(nenergies)
    real*8,     intent(inout) :: dos(nenergies)
    real*8,      intent(in) :: volume
    integer,     intent(in) :: nbands
    real*8,      intent(in) :: E(nbands,4)
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    Yi Yao
    !  HISTORY
    !    Created in April 2019
    !  SOURCE
    ! local variables
    real*8        :: zero, de
    real*8        :: ee(4)
    real*8        :: v, ni, gi, delta
    real*8        :: f12, f13, f14, f23, f24, f34
    real*8        :: f21, f31, f41, f32, f42, f43
    ! counter
    integer       :: ii, kk, m, n, j

    zero = energies(1)
    de = (energies(nenergies) - zero) / (nenergies - 1)
    do ii = 1,nbands
      ee = E(ii,:)
      call sort(ee)
      do j = 1,3
        m = max(1, int((ee(j) - zero) / de) + 2 -1) 
        n = min(nenergies , int((ee(j + 1) - zero) / de) + 1+1)
        if (n > m) then
          !do kk = 1,nenergies
          do kk = m,n
            v = energies(kk)
            if (v >= ee(j) .and. v <= ee(j+1)) then
              v = energies(kk)
              f12 = (v-ee(2)) / (ee(1)-ee(2))
              f13 = (v-ee(3)) / (ee(1)-ee(3))
              f14 = (v-ee(4)) / (ee(1)-ee(4))
              f23 = (v-ee(3)) / (ee(2)-ee(3))
              f24 = (v-ee(4)) / (ee(2)-ee(4))
              f34 = (v-ee(4)) / (ee(3)-ee(4))
              f21 = 1.0 - f12
              f31 = 1.0 - f13
              f41 = 1.0 - f14
              f32 = 1.0 - f23
              f42 = 1.0 - f24
              f43 = 1.0 - f34
              if (j == 1) then
                  ni = f21*f31*f41
                  gi = 3.0 * ni / (v - ee(1))
              else if (j == 2) then
                  delta = ee(4)-ee(1)
                  ni = f42*f32+f41*f24*f32+f41*f31*f23
                  gi = (3.0 / delta)*(f23*f31+f32*f24)
              else if (j == 3) then
                  ni = (1.0 - f14*f24*f34)
                  gi = 3.0 * (1.0 - ni) / (ee(4) - v)
              end if
              !if (IEEE_IS_FINITE(gi)) then
              dos(kk) = dos(kk) + volume * gi
              !end if
            end if
          end do
        end if
      end do
    end do
  end subroutine ltidos_single



  subroutine ltispectral &
    ( cell, kmesh_n1, kmesh_n2, kmesh_n3, nbands, &
      eigs, func, nenergies, energies, spectral_output, comm )
    !  PURPOSE
    !    linear tetrahedron integration method for spectral function
    !  USES
    implicit none
    !  ARGUMENTS
    real*8,      intent(in) :: cell(3,3)
    integer,     intent(in) :: kmesh_n1, kmesh_n2, kmesh_n3, nbands
    real*8,      intent(in) :: eigs(nbands,kmesh_n1,kmesh_n2,kmesh_n3)
    real*8,      intent(in) :: func(nbands,kmesh_n1,kmesh_n2,kmesh_n3)
    integer,     intent(in) :: nenergies
    real*8,      intent(in) :: energies(nenergies)
    real*8,      intent(out) :: spectral_output(nenergies)
    integer,     intent(in) :: comm
    !  INPUTS
    !    o cell - 3d cell vector
    !    o kmesh_n1, kmesh_n2, kmesh_n3 - number of k points in 
    !      each direction
    !    o nbands - number of states in the calculation
    !    o eigs - eigenvalues ( omega in 1.1 in the reference)
    !    o func - function ( F in 1.1 in the reference )
    !    o nenergies - number of points for energies for output
    !    o energies - list of energies for output
    !  OUTPUT
    !    o spectral_output - list of spectral function for output
    !  AUTHOR
    !    Yi Yao
    !  HISTORY
    !    Created in April 2019
    !  SOURCE
    !  local variables
    integer  :: indices(3,8)
    integer  :: kmesh_size(3)
    real*8   :: kpts_submesh(3,8)
    real*8   :: kpts(3,4)
    real*8   :: invcell(3,3), B(3,3)
    integer  :: main_diagonal(2,4)
    integer  :: remain_points(6,4)
    integer  :: list_i_r(6)
    real*8   :: main_diagonal_lengths(4)
    integer  :: index_diagonal_shortest
    integer  :: i_a, i_b
    integer  :: tetra_indecies(4,6)
    integer  :: s(4)
    real*8   :: volume, detkpts
    real*8   :: det_invcell
    real*8   :: E(nbands,4)
    real*8   :: F(nbands,4)
    integer  :: my_rank, mpi_size, mpierr


    !  counters
    integer  :: ii, jj
    integer  :: i_k, j_k, k_k
    integer  :: ii_k, jj_k, kk_k
    integer  :: n
    !

    call mpi_comm_rank(comm,my_rank,mpierr)
    call mpi_comm_size(comm,mpi_size,mpierr)
    

    !
    kmesh_size = (/kmesh_n1,kmesh_n2,kmesh_n3/)
    call matinv3(cell,invcell)
    B = invcell / spread(kmesh_size,1,3)
    indices = reshape((/0,0,0,&
               0,0,1,&
               0,1,0,&
               0,1,1,&
               1,0,0,&
               1,0,1,&
               1,1,0,&
               1,1,1/),shape(indices))
    
    do ii = 1,8
      do jj = 1,3
        kpts_submesh(jj,ii) =  dot_product(indices(:,ii),B(:,jj))
      end do
    end do
    main_diagonal = reshape((/1,8,2,7,4,5,3,6/),shape(main_diagonal))
    remain_points = reshape((/2,4,3,7,5,6,&
                     1,3,4,8,6,5,&
                     1,2,6,8,7,3,&
                     1,2,4,8,7,5/),shape(remain_points))
    do ii = 1,4
      main_diagonal_lengths(ii) = sqrt(sum((kpts_submesh(:,main_diagonal(1,ii)) &
                                     - kpts_submesh(:,main_diagonal(2,ii)))**2))
    end do
    index_diagonal_shortest = minloc(main_diagonal_lengths,1)
    i_a = main_diagonal(1,index_diagonal_shortest)
    i_b = main_diagonal(2,index_diagonal_shortest)
    list_i_r = remain_points(:,index_diagonal_shortest)
    do ii = 1,6
      tetra_indecies(1,ii) = i_a
      tetra_indecies(2,ii) = i_b
      tetra_indecies(3,ii) = list_i_r(ii)
      tetra_indecies(4,ii) = list_i_r(mod(ii,6)+1)
    end do

    spectral_output = 0.0

    do ii = 1,6
      s = tetra_indecies(:,ii)
      do jj = 1,4
        kpts(:,jj) = kpts_submesh(:,s(jj))
      end do
      call det3(kpts(:,1:3) - kpts(:,2:4), detkpts)
      volume = abs(detkpts / 6.0)
      n = 0
      do i_k = 1,kmesh_n1
        do j_k = 1,kmesh_n2
          do k_k = 1,kmesh_n3
            n = n + 1
            if (mod(n, mpi_size) .ne. my_rank) then
              cycle
            end if
            do jj = 1,4
              ii_k = mod((i_k + indices(1,s(jj)) - 1) , kmesh_n1) + 1
              jj_k = mod((j_k + indices(2,s(jj)) - 1) , kmesh_n2) + 1
              kk_k = mod((k_k + indices(3,s(jj)) - 1) , kmesh_n3) + 1
              E(:,jj) = eigs(:,ii_k,jj_k,kk_k)
              F(:,jj) = func(:,ii_k,jj_k,kk_k)
            end do
            call ltispectral_single(nenergies, energies, spectral_output, volume, nbands, E, F)
          end do
        end do
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE, spectral_output, nenergies, MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpierr)

  
    call det3(invcell, det_invcell)
    spectral_output = spectral_output / abs(det_invcell)
  end subroutine ltispectral

  subroutine ltispectral_single(nenergies, energies, spectral, volume, nbands, E, F)
    !  PURPOSE
    !    linear tetrahedron integration method for spectral function
    !    helper subroutine for single tetrahedron
    !  USES
    implicit none
    !  ARGUMENTS
    integer,     intent(in) :: nenergies
    real*8,      intent(in) :: energies(nenergies)
    real*8,     intent(inout) :: spectral(nenergies)
    real*8,      intent(in) :: volume
    integer,     intent(in) :: nbands
    real*8,      intent(in) :: E(nbands,4)
    real*8,      intent(in) :: F(nbands,4)
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    Yi Yao
    !  HISTORY
    !    Created in April 2019
    !  SOURCE
    ! local variables
    real*8        :: zero, de
    real*8        :: ee(4), ff(4)
    real*8        :: v, ni, gi, delta
    real*8        :: f12, f13, f14, f23, f24, f34
    real*8        :: f21, f31, f41, f32, f42, f43
    real*8        :: I1i, I2i, I3i, I4i
    integer       :: isort_ee(4)
    ! counter
    integer       :: ii, kk, m, n, j

    zero = energies(1)
    de = (energies(nenergies) - zero) / (nenergies - 1)
    do ii = 1,nbands
      ee = E(ii,:)
      isort_ee = rargsort(ee)
      call sort(ee)
      ff = F(ii,isort_ee)

      do j = 1,3
        m = max(1, int((ee(j) - zero) / de) + 2 -1) 
        n = min(nenergies , int((ee(j + 1) - zero) / de) + 1+1)
        if (n > m) then
          !do kk = 1,nenergies
          do kk = m,n
            v = energies(kk)
            if (v >= ee(j) .and. v <= ee(j+1)) then
              v = energies(kk)
              f12 = (v-ee(2)) / (ee(1)-ee(2))
              f13 = (v-ee(3)) / (ee(1)-ee(3))
              f14 = (v-ee(4)) / (ee(1)-ee(4))
              f23 = (v-ee(3)) / (ee(2)-ee(3))
              f24 = (v-ee(4)) / (ee(2)-ee(4))
              f34 = (v-ee(4)) / (ee(3)-ee(4))
              f21 = 1.0 - f12
              f31 = 1.0 - f13
              f41 = 1.0 - f14
              f32 = 1.0 - f23
              f42 = 1.0 - f24
              f43 = 1.0 - f34
              if (j == 1) then
                  ni = f21*f31*f41
                  gi = 3.0d0 * ni / (v - ee(1))
                  I1i = 1.0d0 / 3 * ( f12 + f13 + f14)
                  I2i = 1.0d0 / 3 * f21
                  I3i = 1.0d0 / 3 * f31
                  I4i = 1.0d0 / 3 * f41
              else if (j == 2) then
                  delta = ee(4)-ee(1)
                  ni = f42*f32+f41*f24*f32+f41*f31*f23
                  gi = (3.0d0 / delta)*(f23*f31+f32*f24)
                  I1i = f14 / 3 + f13*f31*f23/(gi*delta)
                  I2i = f23 / 3 + f24*f24*f32/(gi*delta)
                  I3i = f32 / 3 + f31*f31*f23/(gi*delta)
                  I4i = f41 / 3 + f42*f24*f32/(gi*delta)
              else if (j == 3) then
                  ni = (1.0d0 - f14*f24*f34)
                  gi = 3.0d0 * (1.0d0 - ni) / (ee(4) - v)
                  I1i = 1.0d0 / 3 * f14
                  I2i = 1.0d0 / 3 * f24
                  I3i = 1.0d0 / 3 * f34
                  I4i = 1.0d0 / 3 * (f41 + f42 + f43)
              end if
              spectral(kk) = spectral(kk) + volume * gi * &
                             ( I1i * ff(1) + I2i * ff(2) + &
                               I3i * ff(3) + I4i * ff(4) )
            end if
          end do
        end if
      end do
    end do
  end subroutine ltispectral_single

  !!!! supporting subroutines 
  subroutine matinv3(A,B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    implicit none
    real*8, intent(in) :: A(3,3)   !! Matrix
    real*8, intent(out) :: B(3,3)   !! Inverse matrix
    real*8             :: detinv
  
    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end subroutine
  
  subroutine det3(A,DET)
    !! calculate the determinant of a 3x3 matrix
  
    implicit none
  
    real*8, intent(in)  :: A(3,3)
  
    real*8 :: DET
  
  
    DET =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)
  
  end subroutine

  subroutine sort(a)
    !! selectrion sort ( I only need to sort array with 4 elements)
    !! sort array a from small to large
    real*8, intent(in out) :: a(:)
    integer :: i, minIndex
    real*8 :: temp
 
    do i = 1, SIZE(a)-1
       minIndex = MINLOC(a(i:), 1) + i - 1
       if (a(i) > a(minIndex)) then
          temp = a(i)
          a(i) = a(minIndex)
          a(minIndex) = temp
       end if
    end do
  end subroutine 

  function rargsort(a) result(b)
  ! Returns the indices that would sort an array.
  !
  ! Arguments
  ! ---------
  !
  real*8, intent(in):: a(:)   ! array of numbers
  integer :: b(size(a))         ! indices into the array 'a' that sort it
  !
  ! Example
  ! -------
  !
  ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]
  
  integer :: N                           ! number of numbers/vectors
  integer :: i,imin                      ! indices: i, i of smallest
  integer :: temp1                       ! temporary
  real*8 :: temp2
  real*8 :: a2(size(a))
  a2 = a
  N=size(a)
  do i = 1, N
      b(i) = i
  end do
  do i = 1, N-1
      ! find ith smallest in 'a'
      imin = minloc(a2(i:),1) + i - 1
      ! swap to position i in 'a' and 'b', if not already there
      if (imin /= i) then
          temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
          temp1 = b(i); b(i) = b(imin); b(imin) = temp1
      end if
  end do
  end function

  !!!!

end module tetrahedron_integration
