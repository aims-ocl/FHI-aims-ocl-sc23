program naive_fourier_transform
! Reads a file that has THREE real values in each line, corresponding to
! time, autocorrelation function, auto correlation function already windowed

 implicit none
 integer :: i, j, k, n_max, i_code
 logical :: eof
 character*40 :: trash
 real*8 :: pi, omega_0, omega, dt, c, factor
 real*8, dimension(:), allocatable :: autocorrel
! real*8, dimension(:), allocatable :: windowed_autocorr
 real*8, dimension(:), allocatable :: fourier
 real*8, dimension(:), allocatable :: imag_fourier
 real*8, dimension(:), allocatable :: time

 pi = 2.*acos(0.0)
 n_max = 0
 eof = .false.
 open(10, file='autocorr.dat')
  read (10,*,iostat = i_code) trash, trash, trash
    if (i_code.ne.0) then
        eof = .true.
    end if
 
 do while (.not.eof)
  read (10,*,iostat = i_code) trash, trash, trash
    if (i_code.ne.0) then
        eof = .true.
    end if
   n_max = n_max+1
 end do

 close(10)

! write(6,*) n_max

 n_max = n_max - 1

 allocate(autocorrel(0:n_max))
! allocate(windowed_autocorr(0:n_max))
 allocate(time(0:n_max))
 allocate(fourier(0:2*n_max))
 allocate(imag_fourier(0:2*n_max))

 eof = .false.
 open(11, file='autocorr.dat')

! reads the file and keeps the time and the autocorrelation function already windowed
 i=0
 do while (.not.eof)
  read (11,*,iostat = i_code) time(i), trash, autocorrel(i)
    if (i_code.ne.0) then
        eof = .true.
    end if
    i = i+1
 end do

 close(11)

 fourier(:) = 0.0
 imag_fourier(:) = 0.0


! now we do a FT of a _symmetric_ function in the range [-n_max;n_max]

 omega_0 = 2.*pi/(2.*time(n_max))
 dt=time(1) ! in ps
 c = 29979245800.0 ! in cm/s
 factor = 1e-12*c !transforms from ps in cm-1

 ! window the function - actually, no, this is now done in the python script
 !call window_the_thing(autocorrel, windowed_autocorr, n_max)

 open(12, file='raw_fourier_transform.dat') 
 
 ! k is the frequency grid
 do k=0, 2*n_max, 1
    omega = omega_0 * k

   ! i = 0 treated separately
   fourier(k) = dt*autocorrel(0)
!   imag_fourier(k) = 0.d0

   do i=1, n_max, 1

      ! the factor of 2 accounts for the weight of points at -i and +i
        fourier(k) = fourier(k) + dt*2.*autocorrel(i)*cos(omega*i*dt)

!      ! this is -i
!      imag_fourier(k) = imag_fourier(k) - autocorrel(i)*sin(omega*i)
!      ! this is +i
!      imag_fourier(k) = imag_fourier(k) + autocorrel(i)*sin(omega*i)

   end do

  if(k.ge.2.and.(omega/(factor*2*pi)).le.5000) then
    write(12,*) omega/(factor*2*pi), fourier(k)
  endif

 end do

 close(12) 

 deallocate(autocorrel)
! deallocate(windowed_autocorr)
 deallocate(fourier)
 deallocate(imag_fourier)
 deallocate(time)

end program naive_fourier_transform

!subroutine window_the_thing(autocorrel, windowed_autocorr, n_max)
! implicit none
! integer :: window_function, n, n_max
! real*8, dimension(0:n_max) :: autocorrel
! real*8, dimension(0:n_max) :: windowed_autocorr
! real*8 :: pi
!
! window_function = 1
! pi = 2.*acos(0.0)
! 
! select case(window_function)
! case(0)
! ! do nothing
!  do n=0, n_max, 1
!    windowed_autocorr(n) = autocorrel(n)
!  end do 
! case (1)
! ! the Hann window, as w_0
!  do n=0, n_max, 1
!    windowed_autocorr(n) = autocorrel(n)*0.5*(1+cos(2.*pi*n/(2*n_max-1.)))
!  end do
! end select
!
!end subroutine window_the_thing
