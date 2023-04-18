program broadener

implicit none

integer, parameter      :: wp = selected_real_kind(4)   ! 4/8/10

real(wp), parameter     :: pi = 4_wp*datan(1d0)

character*64    :: fname

integer :: io_error, read_error, i, k, n, steps, spec_size, nm_lo, nm_hi

real(wp)        :: dummy, spec_inc, alpha, beta, gamma, delta, t_eps, div
real(wp), allocatable, dimension(:)     :: full_spectrum, engy, ostr

open(100, file='TDDFT_LR_Spectrum_Singlet.dat', action='read', iostat=io_error)
if(io_error==0) then
  n=1
  do
    read(100, *, iostat=read_error) dummy, dummy
  if(read_error.gt.0) then
    stop 'TDDFT Data File read Error'
  elseif(read_error.lt.0) then
    exit
  else
    n=n+1 !Count lines
  endif
  enddo
else
  write(*,*) 'TDDFT Data File Error (does TDDFT_LR_Spectrum_Singlet.dat exist?)'
  stop
endif

write(*,*) ''
write(*,*) 'File has ',n-1,' data entries.'
write(*,*) ''

allocate(engy(1:n-1))
allocate(ostr(1:n-1))

rewind(100,iostat=io_error)
if(io_error==0) then
  n=1
  do
    read(100, *, iostat=read_error) engy(n), ostr(n)
  if(read_error.gt.0) then
    stop 'File read Error'
  elseif(read_error.lt.0) then
    exit
  else
    n=n+1
  endif
  enddo
else
  write(*,*) 'TDDFT Data File Error (does TDDFT_LR_Spectrum_Singlet.dat exist?)'
  stop
endif
close(100)

open(100, file='casida_parameters.dat', action='read', iostat=io_error)
if(io_error==0) then
  do
    read(100, *, iostat=read_error) alpha, beta, gamma, delta, t_eps, div, nm_lo, nm_hi
  if(read_error.gt.0) then
    stop 'Parameter File read Error'
  elseif(read_error.lt.0) then
    exit
  else
    write(*,*) ''
    write(*,'(a)') '  --  Input Parameters  --'
    write(*,'(a,f6.2,a,f6.2)') ' -  Lorentz:  Width = ', dble(alpha), ' | Weight = ', dble(beta)
    write(*,'(a,f6.2,a,f6.2)') ' - Gaussian:  Width = ', dble(gamma), ' | Weight = ', dble(delta)
    write(*,'(a,f6.2)') ' - Temperature broadening: ', dble(t_eps)
    write(*,'(a,f6.2)') ' - Total Spectrum Divisor: ', dble(div)
    write(*,'(a,f8.2,a,f8.2,a)') ' - Wavelength plotting within ', dble(nm_lo), ' and ', dble(nm_hi),' nm.'
    write(*,*) ''
  endif
  enddo
else
  write(*,*) 'Parameter File Error (does casida_parameters.dat exist?)'
  stop
endif

engy = ( 1_wp + t_eps/100_wp ) * engy

open(110,file='Sticks_nm.dat')
do i=1, n-1
    write(110,'(f10.3,4x,f24.10)') 1236.7_wp/engy(i), ostr(i)
enddo
close(110)

steps = 1000 ! resolution (steps per eV)
spec_size = 40 * steps ! total spectrum is 0 to 40eV wide
spec_inc = 1.0_wp/dble(steps)

allocate(full_spectrum(1:spec_size))
full_spectrum = 0_wp

do k=1, n-1
  if(mod(k,n/10)==0) write(*,'(f5.2,a)') dble(k)/dble(n)*100_wp,'% done.'
  do i=1, spec_size
    full_spectrum(i) = full_spectrum(i) + &
      beta * ostr(k) / pi * alpha / ( ( engy(k) - dble(i) * spec_inc )**2 + alpha**2 )
    full_spectrum(i) = full_spectrum(i) + &
      delta * ostr(k) / (gamma * sqrt(2_wp * pi) ) * dexp(- ( dble(i) * spec_inc - engy(k) )**2 / (2.0 * gamma**2) )
  enddo
enddo

full_spectrum = full_spectrum / div

open(110,file='Broadened_Spectrum_nm.dat')
do i=1, spec_size
  if( ((1236.7_wp/(float(i)/float(steps))) < dble(nm_hi)) .and. ((1236.7_wp/(float(i)/float(steps))) > dble(nm_lo)) ) then
    write(110,'(f10.3,4x,f24.10)') 1236.7_wp/(float(i)/float(steps)), full_spectrum(i)
  endif
enddo
close(110)

open(110,file='Broadened_Spectrum_eV.dat')
do i=1, spec_size
    write(110,'(f10.3,4x,f24.10)') float(i)/float(steps), full_spectrum(i)
enddo
close(110)

end program
