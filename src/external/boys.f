      module boys

      use dimensions,   only: n_states

      implicit none

c     the full transformation matrix (scnr*)
      real*8, allocatable, dimension(:,:,:) :: backstreet_boys
c     the boys-centers global matrix
      real*8, allocatable, dimension(:,:,:,:) :: boys_centers_mat
c     the boys-transformed KS_eigenvector
      real*8, dimension(:,:,:,:), allocatable    :: KS_temp
c     lowest state considered in localization
      integer, allocatable, dimension(:,:)  :: boys_sub_min_KS_state
c     highest state considered in localization
      integer, allocatable, dimension(:,:)  :: boys_sub_max_KS_state
c     flag controls when the localization is triggered
c     0 = at beginning of SCF
c     1 = throughout the entire cycle
c     2 = at the end, before writing restart info
      integer, allocatable, dimension(:)  :: boys_sub_flags
      integer                             :: boys_sub_flag = 0

      contains

c     allocations ------------------------------------------------------|
      subroutine allocate_boys()
      use dimensions,   only: n_sub_boys, n_spin

      if (.not. allocated(boys_sub_min_KS_state)) then
        allocate(boys_sub_min_KS_state(n_sub_boys,n_spin))
      end if
      if (.not. allocated(boys_sub_max_KS_state)) then
        allocate(boys_sub_max_KS_state(n_sub_boys,n_spin))
      end if
      if (.not. allocated(boys_sub_flags)) then
        allocate(boys_sub_flags(n_sub_boys))
      end if

      end subroutine

c     old subroutine to calculate boys-centers, extended to spin --------|
      subroutine boys_centers(dipx, dipy, dipz, nstates, transfmat)
      use dimensions,           only: n_spin
      implicit none
      integer :: nstates
c     CAVE: This routine is still the old one only supporting spin=none
      complex*16 ::
     &  dipx((nstates*(nstates-1)/2+nstates)*n_spin*(n_spin+1)/2)
      complex*16 ::
     &  dipy((nstates*(nstates-1)/2+nstates)*n_spin*(n_spin+1)/2)
      complex*16 ::
     &  dipz((nstates*(nstates-1)/2+nstates)*n_spin*(n_spin+1)/2)
      real*8            :: x(nstates, nstates, n_spin)
      real*8            :: y(nstates, nstates, n_spin)
      real*8            :: z(nstates, nstates, n_spin)
      integer           :: i_st, j_st, num, i_spin, offset
      real*8, optional, intent(out)  :: transfmat(nstates, nstates)


      num=0
      do i_st=1, nstates
         do j_st=i_st, nstates
            num=num+1
            do i_spin=1, n_spin

               ! We need an offset in case of spin since the dipmat matrix is
               ! stored in a list of a/a, a/b, b/b pairs
               ! Yes, it is -1
               offset = (nstates*(nstates-1)/2+nstates) *
     &                  (i_spin*(i_spin+1)/2 - 1)

               x(i_st, j_st, i_spin)=real(dipx(num+offset))
               x(j_st, i_st, i_spin)=x(i_st, j_st, i_spin)
               y(i_st, j_st, i_spin)=real(dipy(num+offset))
               y(j_st, i_st, i_spin)=y(i_st, j_st, i_spin)
               z(i_st, j_st, i_spin)=real(dipz(num+offset))
               z(j_st, i_st, i_spin)=z(i_st, j_st, i_spin)
            enddo
         enddo
      enddo

      do i_spin=1, n_spin
         call boyscalc(x(:,:,i_spin),
     &                 y(:,:,i_spin),
     &                 z(:,:,i_spin), transfmat, nstates)
      end do
      write(*,*) 'Boys calculated!' 
      call output_xyzboys(x, y, z, nstates)

      end subroutine

c     main routine called from outside doing a localization -------------|
      subroutine boys_transform_subspace(dipx, dipy, dipz, nstates)
      use dimensions,           only: n_spin, n_basis, n_states
      use physics,              only: KS_eigenvector
      use synchronize_mpi,      only: bcast_eigenvector_boys
      use mpi_tasks,            only: myid
      implicit none
      complex*16 ::
     &  dipx((nstates*(nstates-1)/2+nstates)*n_spin*(n_spin+1)/2)
      complex*16 ::
     &  dipy((nstates*(nstates-1)/2+nstates)*n_spin*(n_spin+1)/2)
      complex*16 ::
     &  dipz((nstates*(nstates-1)/2+nstates)*n_spin*(n_spin+1)/2)
      integer    :: nstates

      if (myid.eq.0) then
        call boys_on_subspaces(dipx, dipy, dipz, nstates)
        call boys_transform()
      end if

      call bcast_eigenvector_boys(KS_eigenvector,
     &                            n_states*n_basis*n_spin)

      end subroutine

c     Assemble the entire transfmat ------------------------------------|
      subroutine boys_on_subspaces(dipx, dipy, dipz, nstates)
      use dimensions,   only: n_states, n_sub_boys, n_spin

      implicit none

      complex*16, dimension(:) :: dipx
      complex*16, dimension(:) :: dipy
      complex*16, dimension(:) :: dipz

      integer    :: nstates
      integer    :: nstates_subspace(n_sub_boys, n_spin)
      integer    :: i_sub_boys, i_spin, num

      do i_sub_boys = 1, n_sub_boys
        do i_spin = 1, n_spin
          nstates_subspace(i_sub_boys,i_spin) = 1 +
     &    boys_sub_max_KS_state(i_sub_boys,i_spin) -
     &    boys_sub_min_KS_state(i_sub_boys,i_spin)
       end do
      end do

      call boys_matrix_setup()

      do i_spin = 1, n_spin
        call boys_prepare_dipelms(nstates, dipx, dipy, dipz, i_spin)
      end do

      call output_xyzboys(boys_centers_mat(:,:,:,1),
     &                    boys_centers_mat(:,:,:,2),
     &                    boys_centers_mat(:,:,:,3),
     &                    n_states)

c   ! This is purely a debug print to check the transfmat matrix elements ---------
c
c     open(unit=120, file='backstreet_boys.csv',ACTION='WRITE')
c     do i_spin=1, n_spin
c       do num=1, nstates
c         write(120,'(2X,<nstates>F12.8)')
c    &    backstreet_boys(num, 1:nstates, i_spin)
c       enddo
c     end do
c     close(unit=120)
c
c   ! --------------------------------------------------------------------------
      end subroutine

c     Prepare the <phi_a|r|phi_b> submatrices --------------------------|
      subroutine boys_prepare_dipelms(nstates, dipx, dipy, dipz, i_spin)
      use dimensions, only: n_spin, n_sub_boys
      use physics, only: occ_numbers

      implicit none

      integer,intent(in)   :: nstates, i_spin
      complex*16, dimension(:), intent(in) :: dipx
      complex*16, dimension(:), intent(in) :: dipy
      complex*16, dimension(:), intent(in) :: dipz
      logical                              :: log_col, log_row

c     shell of the full transformation matrix
      real*8    :: x(nstates, nstates), 
     &             y(nstates, nstates),
     &             z(nstates, nstates)
      integer   :: i_st, j_st, num, minb, maxb, offset, i_sub_boys,
     &             n_states_occupied

      ! We need an offset in case of spin since the dipmat matrix is
      ! stored in a list of a/a, a/b, b/b pairs
      ! Yes, it is -1
      offset = (nstates*(nstates-1)/2+nstates) * 
     &         (i_spin*(i_spin+1)/2 - 1)

      x = 0.d0
      y = 0.d0
      z = 0.d0

      do i_st = nstates, 1, -1
        n_states_occupied = i_st

        ! First not non occupied state.
        if(ANY(occ_numbers(i_st, i_spin,:) > 1.d-5)) then
          exit
        end if
      enddo


      do i_sub_boys = 1, n_sub_boys

        minb = boys_sub_min_KS_state(i_sub_boys,i_spin)
        maxb = boys_sub_max_KS_state(i_sub_boys,i_spin)

c       ! Set up the correct input matrix for all subspaces
        num=0
        do i_st=1, nstates
           do j_st=i_st, nstates
              num=num+1
              log_col = (i_st .ge. minb) .and. (i_st .le. maxb) .and.
     &             (j_st .le. n_states_occupied)
              log_row = (j_st .ge. minb) .and. (j_st .le. maxb) .and.
     &             (i_st .le. n_states_occupied)
              if (log_col .or. log_row) then
                x(i_st,j_st)=real(dipx(offset + num))
                x(j_st,i_st)=x(i_st,j_st)
                y(i_st,j_st)=real(dipy(offset + num))
                y(j_st,i_st)=y(i_st,j_st)
                z(i_st,j_st)=real(dipz(offset + num))
                z(j_st,i_st)=z(i_st,j_st)
              end if
           enddo
        enddo
      enddo ! n_sub_boys

c   ! This is purely a debug print to check the dipole matrix elements ---------
c
c     if (i_spin .eq. 1) then
c       open(unit=120, file='subspace_mat_before_1.csv',ACTION='WRITE')
c     else
c       open(unit=120, file='subspace_mat_before_2.csv',ACTION='WRITE')
c     end if
c     do num=1, nstates
c       write(120,'(2X,<nstates>(F12.8,1X))')
c    &  x(num, 1:nstates)
c     enddo
c     do num=1, nstates
c       write(120,'(2X,<nstates>(F12.8,1X))')
c    &  y(num, 1:nstates)
c     enddo
c     do num=1, nstates
c       write(120,'(2X,<nstates>(F12.8,1X))')
c    &  z(num, 1:nstates)
c     enddo
c     close(unit=120)
c
c   ! --------------------------------------------------------------------------
      call boyscalc(x, y, z, backstreet_boys(:,:,i_spin), nstates)
      write(*,*) 'Boys calculated! Spin: ', i_spin

c   ! This is purely a debug print to check the dipole matrix elements ---------
c
c      if (i_spin .eq. 1) then
c        open(unit=120, file='subspace_mat_after_1.csv',ACTION='WRITE')
c      else
c        open(unit=120, file='subspace_mat_after_2.csv',ACTION='WRITE')
c      end if
c      do num=1, nstates
c        write(120,'(2X,<nstates>(F12.8,1X))')
c     &  x(num, 1:nstates)
c      enddo
c      do num=1, nstates
c        write(120,'(2X,<nstates>(F12.8,1X))')
c     &  y(num, 1:nstates)
c      enddo
c      do num=1, nstates
c        write(120,'(2X,<nstates>(F12.8,1X))')
c     &  z(num, 1:nstates)
c      enddo
c      close(unit=120)
c
c   ! --------------------------------------------------------------------------
      boys_centers_mat(:,:,i_spin,1) = x
      boys_centers_mat(:,:,i_spin,2) = y
      boys_centers_mat(:,:,i_spin,3) = z

      end subroutine

c     Apply the Boys transformation to the actually chosen subspace ----|
      subroutine boys_transform()
      use physics,      only: KS_eigenvector
      use dimensions,   only: n_basis, n_spin, n_k_points, n_states
      implicit none
c     KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
      integer                    :: i_basis, i_spin, i_k, i

c     only for debug purposes ------------------------------------
c     real*8  :: backstreet_boys_inspect(n_states, n_states, n_spin)
c     real*8  :: KS_eigenvector_inspect(n_basis, n_states, n_spin,
c    &                                  n_k_points)
c     backstreet_boys_inspect = backstreet_boys
c     KS_eigenvector_inspect = KS_eigenvector
c     only for debug purposes ------------------------------------


      write(*,*) 'Applying Boys transformation ... '

      if (.not. allocated(KS_temp)) then
        allocate(KS_temp(n_basis, n_states, n_spin, n_k_points))
        KS_temp =0.d0
      end if

      do i_spin = 1, n_spin
        do i_k = 1, n_k_points
          call DGEMM('N', 'N', n_basis, n_states, n_states, 1.d0,
     &               KS_eigenvector(:, :, i_spin, i_k), n_basis,
     &               backstreet_boys(:, :, i_spin), n_states, 0.d0,
     &               KS_temp(:, :, i_spin, i_k), n_basis)
        end do
      end do

      KS_eigenvector(:,:,:,:) = KS_temp(:,:,:,:)

      end subroutine

c     setup a n_states x n_states rotation matrix ----------------------|
      subroutine boys_matrix_setup()
      use dimensions, only: n_states, n_spin
      implicit none
      integer :: i, i_spin

      if (.not. allocated(backstreet_boys)) then
        allocate(backstreet_boys(n_states, n_states, n_spin))
      end if
      if (.not. allocated(boys_centers_mat)) then
        allocate(boys_centers_mat(n_states, n_states, n_spin, 3))
      end if


c     Setup an identity matrix
      backstreet_boys = 0.d0
      boys_centers_mat = 0.d0

      do i_spin = 1, n_spin
        do i = 1, n_states
          backstreet_boys(i,i,i_spin) = 1.d0
        end do
      end do

      end subroutine

c    Output boys centers (diagonals of transformed matrices) -----------|
      subroutine output_xyzboys(x, y, z, nstates)
      use geometry
      use constants
      use dimensions
      use species_data 
      implicit none
      integer :: nstates
c     CAVE: This has changed to spin-polarized from initially spin=none
      real*8 :: x(nstates, nstates, n_spin), y(nstates, nstates, n_spin)
      real*8 :: z(nstates, nstates, n_spin)
      integer :: i_st, j_st, i_at, ic, i_spin
      open(unit=120, file='geometry_boys.xyz',ACTION='WRITE')
      write(120, '(2X,I4)') nstates+n_atoms
      write(120, '(2X, A)') '# Atomic geometry and Boys centers '
      do i_at=1, n_atoms
        write(120,'(2X,A,3(2X,F16.8))') 
     &                trim(species_name(species(i_at))),
     &                (coords(ic, i_at)*bohr, ic=1,3,1)
      enddo
      do i_spin=1, n_spin
         write(120, '(2X, A, I2)') '# Spin channel ', i_spin
         do i_st=1, nstates
            write(120,'(2X,A,3(2X,F16.8))') "BOYC ", 
     &            x(i_st, i_st, i_spin)*bohr,
     &            y(i_st, i_st, i_spin)*bohr,
     &            z(i_st, i_st, i_spin)*bohr
         enddo
      end do
      close(unit=120)
      end subroutine


c    The following routine is courtesy of David Manolopoulos. ----------|
      subroutine boyscalc (x,y,z,v,n)
      implicit double precision (a-h,o-z)
      integer n
c
c     ------------------------------------------------------------------ 
c     Transforms the matrices x, y, z from the occupied MO basis to a 
c     Boys localised orbital basis. 
c
c     The coordinates of the Boys orbital centres are returned in the 
c     diagonals of x, y, z, and the orthogonal transformation from the 
c     MO basis to the Boys orbital basis is returned in v. The number of
c     Jacobi sweeps required for convergence is returned in nsweep.
c
c     This subroutine could be made roughly twice as fast by exploiting
c     the symmetry of x, y, z (see the Numerical Recipes jacobi routine)
c     but I have not bothered with this.
c     ------------------------------------------------------------------ 
c
      dimension x(n,n),y(n,n),z(n,n)
      dimension v(n,n)

      integer i, j, k, nsweep
c
      v = 0.d0
      objfn = 0.d0
      do k = 1,n
         v(k,k) = 1.d0
         objfn = objfn+x(k,k)**2+y(k,k)**2+z(k,k)**2
      enddo
c
      do nsweep = 1,10000
         do j = 2,n
            do i = 1,j-1
               dx = x(i,i)-x(j,j)
               dy = y(i,i)-y(j,j)
               dz = z(i,i)-z(j,j)
               a = x(i,j)**2+y(i,j)**2+z(i,j)**2
               b = x(i,j)*dx+y(i,j)*dy+z(i,j)*dz
               c = dx**2+dy**2+dz**2
               a = a-c/4
               if (a.ne.0.d0 .or. b.ne.0.d0) then
                  gamma = datan2(b,-a)/4
                  c = dcos(gamma)
                  s = dsin(gamma)
                  do k = 1,n
                     xki = x(k,i)*c+x(k,j)*s
                     xkj = x(k,j)*c-x(k,i)*s
                     x(k,i) = xki
                     x(k,j) = xkj
                     yki = y(k,i)*c+y(k,j)*s
                     ykj = y(k,j)*c-y(k,i)*s
                     y(k,i) = yki
                     y(k,j) = ykj
                     zki = z(k,i)*c+z(k,j)*s
                     zkj = z(k,j)*c-z(k,i)*s
                     z(k,i) = zki
                     z(k,j) = zkj
                     vki = v(k,i)*c+v(k,j)*s
                     vkj = v(k,j)*c-v(k,i)*s
                     v(k,i) = vki
                     v(k,j) = vkj
                  enddo
                  do k = 1,n
                     xik = c*x(i,k)+s*x(j,k)
                     xjk = c*x(j,k)-s*x(i,k)
                     x(i,k) = xik
                     x(j,k) = xjk
                     yik = c*y(i,k)+s*y(j,k)
                     yjk = c*y(j,k)-s*y(i,k)
                     y(i,k) = yik
                     y(j,k) = yjk
                     zik = c*z(i,k)+s*z(j,k)
                     zjk = c*z(j,k)-s*z(i,k)
                     z(i,k) = zik
                     z(j,k) = zjk
                  enddo
               endif
            enddo
         enddo
         oldfn = objfn
         objfn = 0.d0
         do k = 1,n
            objfn = objfn+x(k,k)**2+y(k,k)**2+z(k,k)**2
         enddo
c         ! uncommented to hide unneccessary output
c         write(*,*) objfn
        if (objfn .le. oldfn) return
      enddo
      stop 'Boys calc did not converge in 10000 steps. Check calculation
     &      and change increase nsweeps in external/boys.f'! Not converged after 10000 Jacobi sweeps
      end subroutine

c     deallocations ------------------------------------------------------|
      subroutine deallocate_boys()
        
        if (allocated(backstreet_boys)) deallocate(backstreet_boys)
        if (allocated(boys_centers_mat)) deallocate(boys_centers_mat)
        if (allocated(KS_temp)) deallocate(KS_temp)
        if (allocated(boys_sub_min_KS_state)) then
          deallocate(boys_sub_min_KS_state)
        end if
        if (allocated(boys_sub_max_KS_state)) then
          deallocate(boys_sub_max_KS_state)
        end if
        if (allocated(boys_sub_flags)) then
          deallocate(boys_sub_flags)
        end if

      end subroutine

      end module
