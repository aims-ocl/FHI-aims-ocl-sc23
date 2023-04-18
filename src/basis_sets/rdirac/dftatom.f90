! Copyright (c) 2006-2012 Ondřej Čertík, Jiří Vackář, John Pask
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.


module types_dftatom

! This module defines the 'dp' double precision type.

implicit none
private
public dp

integer, parameter :: dp=kind(0.d0)          ! double precision

end module ! types





module utils_dftatom

! Various utilities for general use in Fortran programs.

use types_dftatom, only: dp
implicit none
private
public upcase, lowcase, white_char, blank, num_strings, get_string, &
    stop_error, savetxt, newunit, assert, str

interface str
    module procedure str_int, str_real, str_real_n
end interface

contains

function upcase(s) result(t)
! Returns string 's' in uppercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z')) then
        ! if lowercase, make uppercase
        t(i:i) = char(ichar(t(i:i)) + diff)
    end if
end do
end function

function lowcase(s) result(t)
! Returns string 's' in lowercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
        ! if uppercase, make lowercase
        t(i:i) = char(ichar(t(i:i)) - diff)
    end if
end do
end function

logical function white_char(char) ! white character
! returns .true. if char is space (32) or tab (9), .false. otherwise
character, intent(in) :: char
if (iachar(char) == 32 .or. iachar(char) == 9) then
    white_char = .true.
else
    white_char = .false.
end if
end function

logical function blank(string)
! Returns true if string contains only white characters
character(*), intent(in) :: string
integer :: i
do i = 1, len(string)
    if (.not. white_char(string(i:i))) exit
end do
blank = (i>len(string))
end function

integer function num_strings(s) result(n)
! Returns number of substrings contained in input string 's' delimited
! by white space.
character(*), intent(in) :: s    ! input string
character(len(s)+2) :: t         ! temporary string to facilitate analysis
integer :: i
t = " " // s // " "
n = 0
do i = 1, len(t)-1
    if (white_char(t(i:i)) .and. .not. white_char(t(i+1:i+1))) n = n + 1
end do
end function

!--------------------------------------------------------------------------------------------------!

subroutine get_string(s,is,ss)
! Returns first substring ss in string s, delimited by white space, starting at
! index is in s. If ss is found, is is set to (index of last character of ss in
! s) + 1; else is is set to 0. If is is out of range on input, routine
! terminates with is = -1.
character(*), intent(in) :: s   ! input string
integer, intent(inout) :: is    ! on input: starting index for search for ss in
                                ! s on output: (index of last character of ss in
                                ! s) + 1
character(*), intent(out) :: ss ! first substring in s, starting from index is
character(len(s)+1) :: t        ! temporary string to facilitate search
integer i, i1, i2
logical prevwhite, curwhite
if (is <= 0 .or. is > len(s)) then
    ss = ""; is = -1; return
end if
t = s // " "
if (is == 1) then
    prevwhite = .true.
else
    prevwhite = white_char(t(is-1:is-1))
end if
i1 = 0; i2 = 0
do i = is, len(t)
    curwhite = white_char(t(i:i))
    if (prevwhite .and. .not. curwhite) i1 = i   ! beginning of substring
    if (i1>0 .and. curwhite) then                ! end of substring
        i2 = i-1; exit
    end if
    prevwhite=curwhite
end do
if (i2 > 0) then
    ss = t(i1:i2); is = i2+1
else
    ss = ""; is = 0
end if
end subroutine

integer function newunit(unit) result(n)
! Returns lowest i/o unit number not in use (to be used in older compilers).
!
! Starting at 10 to avoid lower numbers which are sometimes reserved.
! Note: largest valid unit number may be system-dependent.
!
! Arguments
! ---------
!
! If present, the new unit will be returned into it
integer, intent(out), optional :: unit
!
! Example
! -------
!
! integer :: u
! open(newunit(u), file="log.txt", status="old")
! read(u, *) a, b
! close(u)
!
! In new compilers, just use the "newunit" keyword argument:
!
! integer :: u
! open(newunit=u, file="log.txt", status="old")
! read(u, *) a, b
! close(u)

logical inuse
integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
integer, parameter :: nmax=999  ! may be system-dependent
do n = nmin, nmax
    inquire(unit=n, opened=inuse)
    if (.not. inuse) then
        if (present(unit)) unit=n
        return
    end if
end do
call stop_error("newunit ERROR: available unit not found.")
end function

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine

subroutine savetxt(filename, d)
! Saves a 2D array into a textfile.
!
! Arguments
! ---------
!
character(len=*), intent(in) :: filename  ! File to save the array to
real(dp), intent(in) :: d(:, :)           ! The 2D array to save
!
! Example
! -------
!
! real(dp) :: data(3, 2)
! call savetxt("log.txt", data)

integer :: s, i
open(newunit(s), file=filename, status="replace")
do i = 1, size(d, 1)
    write(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine assert(condition)
! If condition == .false., it aborts the program.
!
! Arguments
! ---------
!
logical, intent(in) :: condition
!
! Example
! -------
!
! call assert(a == 5)

if (.not. condition) call stop_error("Assert failed.")
end subroutine

pure integer function str_int_len(i) result(sz)
! Returns the length of the string representation of 'i'
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, '(i0)') i
sz = len_trim(s)
end function

pure function str_int(i) result(s)
! Converts integer "i" to string
integer, intent(in) :: i
character(len=str_int_len(i)) :: s
write(s, '(i0)') i
end function

pure integer function str_real_len(r, fmt) result(sz)
! Returns the length of the string representation of 'i'
real(dp), intent(in) :: r
character(len=*), intent(in) :: fmt
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fmt) r
sz = len_trim(s)
end function

pure function str_real(r) result(s)
! Converts the real number "r" to string with 7 decimal digits.
real(dp), intent(in) :: r
character(len=*), parameter :: fmt="(f0.6)"
character(len=str_real_len(r, fmt)) :: s
write(s, fmt) r
end function

pure function str_real_n(r, n) result(s)
! Converts the real number "r" to string with 'n' decimal digits.
real(dp), intent(in) :: r
integer, intent(in) :: n
character(len=str_real_len(r, "(f0." // str_int(n) // ")")) :: s
write(s, "(f0." // str_int(n) // ")") r
end function

end module ! utils





module mesh_dftatom

! Contains mesh utilities (creating the exponential mesh and its derivatives).

use types_dftatom, only: dp
use utils_dftatom, only: stop_error

implicit none

private
public mesh_exp, mesh_exp_deriv, get_mesh_exp_params, mesh_exp_deriv2, &
    linspace, meshgrid

contains

function mesh_exp(r_min, r_max, a, N) result(mesh)
! Generates exponential mesh of N elements on [r_min, r_max]
!
! Arguments
! ---------
!
! The domain [r_min, r_max], the mesh will contain both endpoints:
real(dp), intent(in) :: r_min, r_max
!
! The fraction of the rightmost vs. leftmost elements of the mesh (for a > 1
! this means the "largest/smallest"); The only requirement is a > 0. For a == 1
! a uniform mesh will be returned:
real(dp), intent(in) :: a
!
! The number of elements in the mesh:
integer, intent(in) :: N
!
! Returns
! -------
!
! The generated mesh:
real(dp) :: mesh(N+1)
!
! Note: Every exponential mesh is fully determined by the set of parameters
! (r_min, r_max, a, N). Use the get_mesh_exp_params() subroutine to obtain them
! from the given mesh.
!
! Example
! -------
!
! real(dp) :: r(11)
! r = mesh_exp(0._dp, 50._dp, 1e9_dp, 10)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    call stop_error("mesh_exp: a > 0 required")
else if (abs(a - 1) < tiny(1._dp)) then
    alpha = (r_max - r_min) / N
    do i = 1, N+1
        mesh(i) = alpha * (i-1.0_dp) + r_min
    end do
else
    if (N > 1) then
        beta = log(a) / (N-1)
        alpha = (r_max - r_min) / (exp(beta*N) - 1)
        do i = 1, N+1
            mesh(i) = alpha * (exp(beta*(i-1)) - 1) + r_min
        end do
    else if (N == 1) then
        mesh(1) = r_min
        mesh(2) = r_max
    else
        call stop_error("mesh_exp: N >= 1 required")
    end if
end if
end function

function mesh_exp_deriv(r_min, r_max, a, N) result(Rp)
! Generates dR/dt where R(t) is the mesh returned by mesh_exp()
!
! Input parameters the same as for mesh_exp(). The variable "t" is defined by:
! t = 1, 2, ..., N+1
! So it describes a uniform mesh, with a step size 1, and the corresponding
! physical points are given by the R(t) array.
!
! Output parameters:
!     Rp(N+1) ....... dR/dt
real(dp), intent(in) :: r_min
real(dp), intent(in) :: r_max
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: Rp(N+1)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    call stop_error("mesh_exp_deriv: a > 0 required")
else if (abs(a - 1) < tiny(1._dp)) then
    call stop_error("mesh_exp_deriv: a == 1 not implemented")
else
    if (N > 1) then
        beta = log(a)/(N-1)
        alpha = (r_max - r_min) / (exp(beta*N) - 1)
        do i = 1, N+1
            Rp(i) = alpha * beta * exp(beta*(i-1))
        end do
    else
        call stop_error("mesh_exp_deriv: N > 1 required")
    end if
end if
end function

function mesh_exp_deriv2(r_min, r_max, a, N) result(Rpp)
! Generates d^R/dt^2 where R(t) is the mesh returned by mesh_exp()
!
! Input parameters the same as for mesh_exp(). The variable "t" is defined by:
! t = 1, 2, ..., N+1
! So it describes a uniform mesh, with a step size 1, and the corresponding
! physical points are given by the R(t) array.
!
! Output parameters:
!     Rp(N+1) ....... d^2R/dt^2
real(dp), intent(in) :: r_min
real(dp), intent(in) :: r_max
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: Rpp(N+1)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    call stop_error("mesh_exp_deriv2: a > 0 required")
else if (abs(a - 1) < tiny(1._dp)) then
    call stop_error("mesh_exp_deriv2: a == 1 not implemented")
else
    if (N > 1) then
        beta = log(a)/(N-1)
        alpha = (r_max - r_min) / (exp(beta*N) - 1)
        do i = 1, N+1
            Rpp(i) = alpha * beta**2 * exp(beta*(i-1))
        end do
    else
        call stop_error("mesh_exp_deriv2: N > 1 required")
    end if
end if
end function

subroutine get_mesh_exp_params(R, r_min, r_max, a, N)
! Given any exponential mesh R, it determines the get_mesh()'s parameters
!
! This only looks at the number of elements, the leftmost and the rightmost
! elements (so the middle elements are not checked/taken into account).
real(dp), intent(in) :: R(:)
real(dp), intent(out) :: r_min, r_max, a
integer, intent(out) :: N
r_min = R(1)
r_max = R(size(R))
a = (R(size(R)) - R(size(R)-1)) / (R(2) - R(1))
N = size(R) - 1
end subroutine

function linspace(a, b, n) result(s)
real(dp), intent(in) :: a, b
integer, intent(in) :: n
real(dp) :: s(n)
s = mesh_exp(a, b, 1.0_dp, n-1)
end function

subroutine meshgrid(x, y, x2, y2)
real(dp), intent(in) :: x(:), y(:)
real(dp), intent(out) :: x2(:, :), y2(:, :)
x2 = spread(x, 1, size(y))
y2 = spread(y, 2, size(x))
end subroutine

end module ! mesh





module ode1d_dftatom

! General utilities for solving 1D ODEs. the Adams and rk4 subroutines
! are used by Schroedinger, Dirac and Poisson solvers. The integrate
! function is used at other places in dftatom to calculate integrals of the
! radial density/orbitals.

use types_dftatom, only: dp
use utils_dftatom, only: stop_error

implicit none
private
public integrate, normalize, parsefunction, get_n_nodes, get_min_idx, &
        adams_interp_outward, adams_extrapolation_outward, &
        adams_interp_outward_implicit, &
        adams_interp_inward_implicit, &
        get_midpoints, rk4_integrate3, rk4_integrate4, rk4_integrate, &
        integrate_trapz_1, integrate_trapz_3, integrate_trapz_5, &
        integrate_trapz_7, integrate_simpson, integrate_adams

contains

real(dp) function integrate(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

! Choose one from the integration rules below:
!s = integrate_trapz_1(Rp, f)
!s = integrate_trapz_3(Rp, f)
!s = integrate_trapz_5(Rp, f)
s = integrate_trapz_7(Rp, f)
!s = integrate_simpson(Rp, f)
!s = integrate_adams(Rp, f)
end function

real(dp) function integrate_trapz_1(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (g(1) + g(N)) / 2
s = s + sum(g(2:N-1))
end function

real(dp) function integrate_trapz_3(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (9 * (g(1) + g(N)) + 28 * (g(2) + g(N-1)) + 23 * (g(3) + g(N-2))) / 24
s = s + sum(g(4:N-3))
end function

real(dp) function integrate_trapz_5(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (  475 * (g(1) + g(N  )) &
    + 1902 * (g(2) + g(N-1)) &
    + 1104 * (g(3) + g(N-2)) &
    + 1586 * (g(4) + g(N-3)) &
    + 1413 * (g(5) + g(N-4)) &
    ) / 1440
s = s + sum(g(6:N-5))
end function

real(dp) function integrate_trapz_7(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (  36799 * (g(1) + g(N  )) &
    + 176648 * (g(2) + g(N-1)) &
    +  54851 * (g(3) + g(N-2)) &
    + 177984 * (g(4) + g(N-3)) &
    +  89437 * (g(5) + g(N-4)) &
    + 130936 * (g(6) + g(N-5)) &
    + 119585 * (g(7) + g(N-6)) &
    ) / 120960
s = s + sum(g(8:N-7))
end function

real(dp) function integrate_simpson(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: i, N
N = size(Rp)
g = f * Rp
s = 0
do i = 2, N-1, 2
    s = s + g(i-1) + 4*g(i) + g(i+1)
end do
s = s / 3
if (modulo(N, 2) == 0) then
    ! If N is even, add the last slice separately
    s = s + (5*g(N) + 8*g(N-1) - g(N-2)) / 12
end if
end function

real(dp) function integrate_adams(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: i
s = integrate_trapz_1(Rp(:4), f(:4))
g = f * Rp
do i = 4, size(Rp)-1
    s = s + adams_interp_outward(g, i)
end do
end function

subroutine normalize(Rp, Y)
!     normalizes Y inplace on the grid R
!     I.e. int |Y|^2 dr = 1
!     input parameters:
!     Y (array, in, out) .... the unnormalized wavefunction
!     R (array, in) .... the radial grid
!     output parameters:
!     Y is normalized inplace
real(dp), intent(in) :: Rp(:)
real(dp), intent(inout) :: Y(size(Rp))

real(dp) :: S
S = integrate(Rp, Y**2)
S = sqrt(abs(S))
if (S > 0) then
    Y = Y / S
else
    call stop_error("normalize: zero function")
end if
end subroutine

subroutine parsefunction(y, nodes, minidx, positive)
! parses the function y, returns:
! nodes: the number of intersection with the x-axis, not counting
!   the beginning and infinity
! minidx: the index of the last minimum, i.e. the place, where the
!   function was last closest to the x-axis (in absolute value), so for
!   indexes > minidx, the function goes to +- infinity, or possibly to
!   zero in infinity
! positive: true or false, depending if the y is approaching the x-axis
!   in the infinity from above (positive) or below (negative).
! This information is is used in the routine checke to determine, if the
! function lies below or above a certain energy.
real(dp), intent(in) :: y(:)
integer, intent(out) :: nodes, minidx
logical, intent(out) :: positive

nodes = get_n_nodes(y)
minidx = get_min_idx(y)
positive = y(size(y)) > 0
end subroutine

integer function get_n_nodes(y) result(nodes)
! Returns the number of nodes of the function 'y'
real(dp), intent(in) :: y(:)

integer :: last_sign, last_i, i, isy
nodes = 0
last_sign = int(sign(1.0_dp, y(1)))
last_i = -1

do i = 2, size(y)
  isy = int(sign(1.0_dp, y(i)))
  if (isy == -last_sign) then
      last_sign = isy
      last_i = i - 1
      nodes = nodes + 1
  end if
end do
end function

integer function get_min_idx(y) result(k)
! Returns the index of the last minimum of the function 'y'
real(dp), intent(in) :: y(:)
k = size(y)
do while (abs(y(k-1)) < abs(y(k)))
  k = k-1
  if (k == 1) exit
end do
k = k - 1
end function

real(dp) function adams_extrapolation_outward(y, i) result(r)
! Adams extrapolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(55*y(i) - 59*y(i-1) + 37*y(i-2) - 9*y(i-3)) / 24
end function

real(dp) function adams_interp_outward(y, i) result(r)
! Adams interpolation formula
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(9*y(i+1) + 19*y(i) - 5*y(i-1) + y(i-2)) / 24
end function

real(dp) function adams_interp_outward_implicit(y, i) result(r)
! Adams interpolation formula for implicit method
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = +(19*y(i) - 5*y(i-1) + y(i-2)) / 24
end function

real(dp) function adams_interp_inward_implicit(y, i) result(r)
! Adams interpolation formula for implicit method
real(dp), intent(in) :: y(:)
integer, intent(in) :: i
r = -(19*y(i) - 5*y(i+1) + y(i+2)) / 24
end function

subroutine interp(t, X, n, Y, val)
!     Interpolates the x, y points to estimate the value at the point "t".
!     val (out) ...... the estimated function value at the point "t".
!     t (in).......... the "x" value at which to estimate the value
!     X (array, in) .. the x-values of the points
!     Y (array, in) .. the y-values of the points
!     n (in) ......... the length of the arrays X and Y
real(dp), intent(in) :: t
integer, intent(in) :: n
real(dp), intent(in) :: X(n), Y(n)
real(dp), intent(out) :: val
real(dp) :: f, denum
integer :: i, j
val = 0
do j = 1, n
    f = 1
    denum = 1
    do i=1, n
        if (i == j) cycle
        f = f * (t-X(i))
        denum = denum * (X(j)-X(i))
    end do
    val = val + Y(j) * f / denum
end do
end subroutine


subroutine get_val(f, x, R, n, V, i)
! returns the interpolated value of f(x), where f is defined on the
! grid R
! input parameters
! R ... grid
! n ... length of the grid
! f ... function defined on the grid R
! x ... point at which to interpolate
! i ... the integer number of the radial point at which "x" lies
! output parameters
! V ... the interpolated function value
integer, intent(in) :: n
real(dp), intent(in) :: f(n), x, R(n)
real(dp), intent(out) :: V
integer, intent(in) :: i

integer :: j1, j2, n1, n2
if (n < 4) call stop_error("get_val: n >= 4 required")
j1 = i-1
j2 = j1+1

n1 = j1-1
n2 = j2+1
if (n1 < 1) then
    n2 = n2-n1+1
    n1 = 1
end if
if (n2 > n) then
    n1 = n1-(n2-n)
    n2 = n
end if
call interp(x, r(n1:n2), n2-n1+1, f(n1:n2), V)
end subroutine

function get_midpoints(R, V) result(Vmid)
real(dp), intent(in) :: R(:), V(:)
real(dp) :: Vmid(size(R)-1)
integer :: i
if (.not.(size(R) == size(V))) then
    call stop_error("get_midpoints: incorrect array sizes")
end if
do i = 1, size(R) - 1
    call get_val(V, (R(i) + R(i+1))/2, R, size(R), Vmid(i), i+1)
end do
end function

subroutine rk4_integrate3(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
! Integrates the following set of equations outwards:
! dy1/dx =           y2
! dy2/dx = C1 + C2 * y2
! The above system can be equivalently written as a second order ODE:
! y1'' - C2 * y1' = C1
! The coefficients C1 and C2 are passed in as arrays, so they can have any
! dependence on R. For example for a Poisson equation of the form:
! V'' + 2 * V' / R = -4*pi*n
! we would get C1 = -4*pi*n, C2 = -2/R
real(dp), intent(in) :: R(:) ! Grid
real(dp), intent(in) :: y0(:) ! Initial condition
! Coefficients C1 and C2 at grid points and midpoints:
real(dp), intent(in) :: C1(:), C2(:), C1mid(:), C2mid(:)
! Maximum value (if y1 > max_val, the integration stops)
real(dp), intent(in) :: max_val
! Solution y1 and y2
real(dp), intent(out) :: y1(:), y2(:)
! The integration stops at R(imax)
integer, intent(out) :: imax

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(1) = y(1)
y2(1) = y(2)
do i = 2, n
    ! rk4 step size
    h = R(i)-R(i-1)

    ! Do rk4 step:
    dydx(1) =                      y(2)
    dydx(2) = C1(i-1)  + C2(i-1) * y(2)
    yt = y + h/2 * dydx
    dyt(1) =                           yt(2)
    dyt(2) = C1mid(i-1) + C2mid(i-1) * yt(2)
    yt = y + h/2 * dyt
    dym(1) =                           yt(2)
    dym(2) = C1mid(i-1) + C2mid(i-1) * yt(2)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) =                 yt(2)
    dyt(2) = C1(i) + C2(i) * yt(2)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val) then
        imax = i
        return
    end if
end do
imax = n
end subroutine

subroutine rk4_integrate4(R, y0, C, Cmid, max_val, y1, y2, imax)
! Integrates the following set of equations outwards:
! dy1/dx = C(:, 1, 1) * y1 + C(:, 1, 2) * y2
! dy2/dx = C(:, 2, 1) * y1 + C(:, 2, 2) * y2
real(dp), intent(in) :: R(:) ! Grid
real(dp), intent(in) :: y0(:) ! Initial condition
! Coefficients C1 and C2 at grid points and midpoints:
real(dp), intent(in) :: C(:, :, :), Cmid(:, :, :)
! Maximum value (if y1 > max_val, the integration stops)
real(dp), intent(in) :: max_val
! Solution y1 and y2
real(dp), intent(out) :: y1(:), y2(:)
! The integration stops at R(imax)
integer, intent(out) :: imax

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(1) = y(1)
y2(1) = y(2)
do i = 2, n
    ! rk4 step size
    h = R(i)-R(i-1)

    ! Do rk4 step:
    dydx(1) = C(i-1, 1, 1) * y(1) + C(i-1, 1, 2) * y(2)
    dydx(2) = C(i-1, 2, 1) * y(1) + C(i-1, 2, 2) * y(2)
    !dydx = matmul(C(i-1, :, :), y)
    yt = y + h/2 * dydx
    dyt(1) = Cmid(i-1, 1, 1) * yt(1) + Cmid(i-1, 1, 2) * yt(2)
    dyt(2) = Cmid(i-1, 2, 1) * yt(1) + Cmid(i-1, 2, 2) * yt(2)
    !dyt = matmul(Cmid(i-1, :, :), yt)
    yt = y + h/2 * dyt
    dym(1) = Cmid(i-1, 1, 1) * yt(1) + Cmid(i-1, 1, 2) * yt(2)
    dym(2) = Cmid(i-1, 2, 1) * yt(1) + Cmid(i-1, 2, 2) * yt(2)
    !dym = matmul(Cmid(i-1, :, :), yt)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) = C(i, 1, 1) * yt(1) + C(i, 1, 2) * yt(2)
    dyt(2) = C(i, 2, 1) * yt(1) + C(i, 2, 2) * yt(2)
    !dyt = matmul(C(i, :, :), yt)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val .or. abs(y(2)) >= max_val) then
        imax = i
        return
    end if
end do
imax = n
end subroutine

subroutine rk4_integrate(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
! Integrates the following set of equations outwards:
! dy1/dx =                y2
! dy2/dx = C1 * y1 + C2 * y2
real(dp), intent(in) :: R(:) ! Grid
real(dp), intent(in) :: y0(:) ! Initial condition
! Coefficients C1 and C2 at grid points and midpoints:
real(dp), intent(in) :: C1(:), C2(:), C1mid(:), C2mid(:)
! Maximum value (if y1 > max_val, the integration stops)
real(dp), intent(in) :: max_val
! Solution y1 and y2
real(dp), intent(out) :: y1(:), y2(:)
! The integration stops at R(imax)
integer, intent(out) :: imax

integer :: i
integer :: n
real(dp), dimension(size(y0)) :: dym, dyt, yt, dydx, y
real(dp) :: h

n = size(R)
y = y0
y1(1) = y(1)
y2(1) = y(2)
do i = 2, n
    ! rk4 step size
    h = R(i)-R(i-1)

    ! Do rk4 step:
    dydx(1) =                            y(2)
    dydx(2) = C1(i-1) * y(1) + C2(i-1) * y(2)
    yt = y + h/2 * dydx
    dyt(1) =                                   yt(2)
    dyt(2) = C1mid(i-1) * yt(1) + C2mid(i-1) * yt(2)
    yt = y + h/2 * dyt
    dym(1) =                                   yt(2)
    dym(2) = C1mid(i-1) * yt(1) + C2mid(i-1) * yt(2)
    yt = y + h * dym
    dym = dyt + dym
    dyt(1) =                         yt(2)
    dyt(2) = C1(i) * yt(1) + C2(i) * yt(2)
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    if (abs(y(1)) >= max_val) then
        imax = i
        return
    end if
end do
imax = n
end subroutine

end module ! ode1d





module rschroed_dftatom

! Routines in this module solve the radial Schroedinger equation outward and
! inward using the implicit Adams method.

use types_dftatom, only: dp
use ode1d_dftatom, only: adams_interp_outward_implicit, &
    adams_interp_inward_implicit, get_midpoints, rk4_integrate
use utils_dftatom, only: stop_error

implicit none
private
public schroed_outward_adams, schroed_inward_adams

contains

subroutine schroed_outward_adams(l, Z, E, R, Rp, V, P, Q, imax)
! Integrates the Schrodinger equation outward using Adams 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form outwards using
! predictor-corrector:
! P' = Q
! Q' = P'' = C * P
! where C = 2*(V-E) + l*(l+1)/R**2
!
! Returns P(r), Q(r).
!
! It tranforms the problem to a uniform mesh 1 <= t <= N + 1, defined by
! R(t) and Rp(t) = dR/dt and Rpp(t) = d^2R/dt^2:
!
! u(t)   = P(R(t))
! up(t)  = du/dt  = dP/dR * dR/dt             = Q * Rp
! upp(t) = dup/dt = dQ/dR * Rp^2 + Q * dRp/dt = C*P*Rp^2 + Q * Rpp
!
! So
!
! up  = u * Rp
! upp = u * C * Rp^2 + up * Rpp / Rp
!
! For example if Rp = al*R, Rpp = al^2 * R, we get:
!
! up  = u * al * R
! upp = u * C * al^2 * R^2 + up * al

integer, intent(in) :: l
real(dp), intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: P(:), Q(:)
real(dp), parameter :: max_val = 1e6_dp
integer, intent(out) :: imax ! The integration was carried to R(imax)

real(dp), dimension(size(R)) :: C, u1, u2, u1p, u2p, Vmid
integer :: i
real(dp) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp

if (size(R) < 5) call stop_error("size(R) < 5")
if (.not. (size(R) == size(Rp) .and. size(R) == size(V) .and. &
    size(R) == size(P) .and. size(R) == size(P) .and. size(R) == size(Q))) then
    call stop_error("Array sizes mismatch")
end if
C = (l*(l+1)/R**2 + 2 * (V-E))
Vmid(:3) = get_midpoints(R(:4), V(:4))
call integrate_rschroed_rk4(l, Z, E, R(:4), V(:4), Vmid(:3), &
    u1(:4), u2(:4), imax)
!u1(1:4) = R(1:4) ** (l+1)
!u2(1:4) = (l+1) * R(1:4) ** l
u1p(1:4) = Rp(1:4)          * u2(1:4)
u2p(1:4) = Rp(1:4) * C(1:4) * u1(1:4)

do i = 4, size(R)-1
    u1p(i) = Rp(i)        * u2(i)
    u2p(i) = Rp(i) * C(i) * u1(i)
    u1_tmp  = u1(i) + adams_interp_outward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_outward_implicit(u2p, i)

    lambda = 9.0_dp / 24
    Delta = 1 - lambda**2 * C(i+1) * Rp(i+1)**2
    M(1, 1) = 1 / Delta
    M(2, 1) = lambda * C(i+1) * Rp(i+1) / Delta
    M(1, 2) = lambda * Rp(i+1) / Delta
    M(2, 2) = 1 / Delta

    u1(i+1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i+1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
    if (abs(u1(i+1)) >= max_val .or. abs(u2(i+1)) >= max_val) then
        P = u1
        Q = u2
        imax = i
        return
    end if
end do
P = u1
Q = u2
imax = size(R)
end subroutine

subroutine schroed_inward_adams(l, E, R, Rp, V, P, Q, imin)
! Integrates the Schrodinger equation inward using Adams 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form inwards using
! predictor-corrector:
! P' = Q
! Q' = P'' = (2*(V-E) + l*(l+1)/R**2) * P
!
! Returns P(r), Q(r).
!
! Important note: This only works for exponential meshes, where the mesh
! marameter (as defined by mesh_exp()) is equal to "a = (rmax/rmin)**((N-1)/N)".
! Otherwise it will produce wrong answer.
integer, intent(in) :: l
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: P(:), Q(:)
real(dp), parameter :: max_val = 1e300_dp
integer, intent(out) :: imin

real(dp), dimension(size(R)) :: C, u1, u2, u1p, u2p
integer :: i, i_max
real(dp) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp
real(dp) :: R_max

C = (l*(l+1)/R**2 + 2 * (V-E))

i_max = size(R)-4
if (i_max < 4) call stop_error("imax < 4")
lambda = sqrt(-2*E)
! We require that exp(-lambda*(R-R(1)) ~ epsilon(1.0_dp),
! if we start further from
! the origin, it might sometimes blow up, if we start closer, we might not get
! as precise asymptotic.
! It follows that R ~ R(1)-log(epsilon(1.0_dp)) / lambda
R_max = R(1)-log(epsilon(1.0_dp))/lambda
do while (R(i_max) > R_max)
    if (i_max == 2) then
        call stop_error("Can't start the inward integration")
    end if
    i_max = i_max - 1
end do

u1(i_max:i_max+4) = exp(-lambda * (R(i_max:i_max+4)-R(1)))
u2(i_max:i_max+4) = - lambda * u1(i_max:i_max+4)
u1p(i_max:i_max+4) = Rp(i_max:i_max+4)                * u2(i_max:i_max+4)
u2p(i_max:i_max+4) = Rp(i_max:i_max+4) * C(i_max:i_max+4) * u1(i_max:i_max+4)

do i = i_max, 2, -1
    u1p(i) = Rp(i)        * u2(i)
    u2p(i) = Rp(i) * C(i) * u1(i)
    u1_tmp  = u1(i) + adams_interp_inward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_inward_implicit(u2p, i)

    lambda = -9.0_dp / 24
    Delta = 1 - lambda**2 * C(i-1) * Rp(i-1)**2
    M(1, 1) = 1 / Delta
    M(2, 1) = lambda * C(i-1) * Rp(i-1) / Delta
    M(1, 2) = lambda * Rp(i-1) / Delta
    M(2, 2) = 1 / Delta

    u1(i-1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i-1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
    if (abs(u1(i-1)) >= max_val .or. abs(u2(i-1)) >= max_val) then
        P = u1
        Q = u2
        P(i_max+1:) = 0
        Q(i_max+1:) = 0
        imin = i
        return
    end if
end do
P = u1
Q = u2
P(i_max+4:) = 0
Q(i_max+4:) = 0
imin = 1
end subroutine

subroutine integrate_rschroed_rk4(l, Z, E, R, V, Vmid, P, Q, imax)
! Integrates the Schrodinger equation outward using Runge-Kutta 4th order method
!
! It integrates the Schroedinger equation in the R(r) form outwards:
! R'' = -2/R * R' + (2*(V-E) + l*(l+1)/R**2) * R
!
! Returns P(r), Q(r), where these are defined as:
! P(r) = r*R(r)
! Q(r) = P'(r) = r * R'(r) + R(r)
! where R(r) is the radial solution.
integer, intent(in) :: l
real(dp), intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:)
real(dp), intent(in) :: V(:), Vmid(:)
real(dp), intent(out) :: P(:), Q(:)
integer, intent(out) :: imax ! The integration was carried to R(imax)

real(dp), parameter :: max_val = 1e6_dp
real(dp) :: y0(2), y1(size(R)), y2(size(R))
real(dp), dimension(size(R)) :: C1, C2
real(dp), dimension(size(R)-1) :: C1mid, C2mid
real(dp), dimension(size(R)-1) :: Rmid

! y(:) are values of the components (2) of the equation, in our case:
! y(1) = R, y(2) = R'
! dydx(:) are the derivatives of the components (2) of the equation

if (l == 0) then
    y0(1) = 1-Z*R(1)
    y0(2) = -Z
else
    y0(1) = R(1)**l
    y0(2) = l*R(1)**(l-1)
end if
if (size(V) /= size(Vmid) + 1) call stop_error("Vmid size is wrong")

C1 = 2*(V-E) + l*(l+1)/R**2
C2 = -2/R
Rmid = (R(:size(R)-1) + R(2:)) / 2
C1mid = 2*(Vmid-E) + l*(l+1)/Rmid**2
C2mid = -2/Rmid
call rk4_integrate(R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
P(:imax) = y1(:imax)*R(:imax) ! P(r) = r * R(r)
Q(:imax) = y2(:imax)*R(:imax) + y1(:imax) ! Q(r) = P'(r) = r * R'(r) + R(r)
end subroutine

end module ! rschroed





module rdirac_dftatom

! Routines in this module solve the radial Dirac equation outward and
! inward using the implicit Adams method.

use types_dftatom, only: dp
use ode1d_dftatom, only: adams_interp_outward_implicit, &
    adams_interp_inward_implicit, get_midpoints, rk4_integrate4
use utils_dftatom, only: stop_error

implicit none
private
public dirac_outward_adams, dirac_inward_adams


contains

subroutine dirac_outward_adams(c, kappa, Z, E, R, Rp, V, P, Q, imax)
!integrates the Dirac eq, returns r*R, where R is the radial
!solution.
! input parameters:
! kappa .... the kappa in the Dirac eq.
! Z .... the nucleus charge in Hartree atomic units
! E .... the energy at which to integration the equation
! R .... radial grid
! V .... potential on the radial grid
! c .... speed of light
! output parameters:
! Q .... f-component of the radial dirac wave function
! P .... g-component of the radial dirac wave function
real(dp), intent(in) :: c
integer, intent(in) :: kappa
real(dp), intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: Q(:)
real(dp), intent(out) :: P(:)
integer, intent(out) :: imax

real(dp), parameter :: max_val = 1e6_dp
real(dp), dimension(size(R), 2, 2) :: Ctot

real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: i
real(dp) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp
real(dp) :: Vmid(3)

if (size(R) < 4) call stop_error("size(R) <= 4")
Vmid = get_midpoints(R(:4), V(:4))
call integrate_radial_dirac_r_rk4(c, kappa, Z, E, R(:4), V(:4), Vmid, &
    u1(:4), u2(:4), imax)
if (imax /= 4) call stop_error("rk4 failed")

Ctot(:, 1, 1) = -kappa / R
Ctot(:, 2, 2) = +kappa / R
Ctot(:, 1, 2) = +(E-V)/c + 2*c
Ctot(:, 2, 1) = -(E-V)/c

u1p(:4) = Rp(:4) * (Ctot(:4, 1, 1) * u1(:4) + Ctot(:4, 1, 2) * u2(:4))
u2p(:4) = Rp(:4) * (Ctot(:4, 2, 1) * u1(:4) + Ctot(:4, 2, 2) * u2(:4))

do i = 4, size(R)-1
    u1p(i) = Rp(i) * (Ctot(i, 1, 1)*u1(i) + Ctot(i, 1, 2)*u2(i))
    u2p(i) = Rp(i) * (Ctot(i, 2, 1)*u1(i) + Ctot(i, 2, 2)*u2(i))
    u1_tmp  = u1(i) + adams_interp_outward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_outward_implicit(u2p, i)

    lambda = 9.0_dp / 24
    Delta = 1 - lambda**2 * Rp(i+1)**2 * (Ctot(i+1, 1, 2) * Ctot(i+1, 2, 1) &
        -Ctot(i+1, 1, 1) * Ctot(i+1, 2, 2))
    M(1, 1) = (1 - lambda * Rp(i+1) * Ctot(i+1, 2, 2)) / Delta
    M(2, 1) = lambda * Rp(i+1) * Ctot(i+1, 2, 1) / Delta
    M(1, 2) = lambda * Rp(i+1) * Ctot(i+1, 1, 2) / Delta
    M(2, 2) = (1 - lambda * Rp(i+1) * Ctot(i+1, 1, 1)) / Delta

    u1(i+1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i+1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
                      ! write(6,"('   i=',i4,3x,'u1=',f12.4,10x,'M1=',f12.4,3x,'M2=',f12.4,3x,'u1_tmp=',e20.8,3x,'u2_tmp=',e20.8)")i,u1(i+1),M(1,1),M(1,2),u1_tmp,u2_tmp
    if (abs(u1(i+1)) >= max_val .or. abs(u2(i+1)) >= max_val) then
        P = u1
        Q = u2
        imax = i
        return
    end if
end do
P = u1
Q = u2
imax = size(R)
end subroutine

subroutine dirac_inward_adams(c, kappa, E, R, Rp, V, P, Q, imin)
!integrates the Dirac eq. inwards, returns r*R, where R is the radial
!solution.
! input parameters:
! kappa .... the kappa in the Dirac eq.
! E .... the energy at which to integration the equation
! R .... radial grid
! V .... potential on the radial grid
! c .... speed of light
! output parameters:
! f .... f-component of the radial dirac wave function
! P .... g-component of the radial dirac wave function
real(dp), intent(in) :: c
integer, intent(in) :: kappa
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:), Rp(:)
real(dp), intent(in) :: V(:)
real(dp), intent(out) :: Q(:)
real(dp), intent(out) :: P(:)
integer, intent(out) :: imin

integer :: nr
integer :: i_max
real(dp), parameter :: max_val = 1e20_dp
real(dp) :: lambda
real(dp), dimension(size(R), 2, 2) :: Ctot
real(dp), dimension(size(R)) :: u1, u2, u1p, u2p
integer :: i
real(dp) :: Delta, M(2, 2), u1_tmp, u2_tmp
real(dp) :: R_max

nr = size(R)

if (E > 0) call stop_error("E < 0 required")

i_max = nr-4
if (i_max < 2) call stop_error("size(R) too small to start inward integraion")
lambda = sqrt(-2*E-E**2/c**2)
! We require that exp(-lambda*(R-R(1)) ~ epsilon(1.0_dp),
! if we start further from
! the origin, it might sometimes blow up, if we start closer, we might not get
! as precise asymptotic.
! It follows that R ~ R(1)-log(epsilon(1.0_dp)) / lambda
R_max = R(1)-log(epsilon(1.0_dp))/lambda
do while (R(i_max) > R_max)
    if (i_max == 2) then
        print *, "E =", E, "lambda =", lambda
        call stop_error("Can't start the inward integration")
    end if
    i_max = i_max - 1
end do
! (Rundong) The original version of dft_atom does not support our cutoff strategy very well.
! The reason is: the inward integration starts from a point where the potential could be 
! very large, and thus, finally, the value on the innermost point could be too large that,
! I believe, it exceeds the pricision of double precision.
! I added the following IF sentence and the code works well now. But, perhaps, under some
! extreme situations, bugs could still occur. Let's see.
if(V(i_max).gt.90000.d0)then
  do while (V(i_max) > 90000.d0)
    i_max = i_max - 1
  end do
  i_max=i_max+1
 !i_max=i_max-4
endif
                      ! write(6,"('E=',f12.6,4x,'lambda=',f12.6,4x,'R_max=',f12.6,10x,'i_max=',i3,4x,'R(1)=',f12.6,4x,'R(i_max)=',f12.6)")E,lambda,R_max,i_max,R(1),R(i_max)

u1(i_max+4:) = 0
u2(i_max+4:) = 0

u1(i_max:i_max+4) = exp(-lambda * (R(i_max:i_max+4)-R(1))) / sqrt(-E/(E+2*c**2))
u2(i_max:i_max+4) = -exp(-lambda * (R(i_max:i_max+4)-R(1)))

Ctot(:, 1, 1) = -kappa / R
Ctot(:, 2, 2) = +kappa / R
Ctot(:, 1, 2) = +(E-V)/c + 2*c
Ctot(:, 2, 1) = -(E-V)/c
                      ! write(6,"('V:')")
                      ! write(6,"(20f12.4)")V(:i_max+4)
                      ! write(6,"('Ctot:')")
                      ! write(6,"(20f12.4)")Ctot(:i_max+4,2,1)

u1p(i_max:i_max+4) = Rp(i_max:i_max+4) * &
    (Ctot(i_max:i_max+4, 1, 1)*u1(i_max:i_max+4) &
        + Ctot(i_max:i_max+4, 1, 2) * u2(i_max:i_max+4))
u2p(i_max:i_max+4) = Rp(i_max:i_max+4) * &
    (Ctot(i_max:i_max+4, 2, 1)*u1(i_max:i_max+4) &
        + Ctot(i_max:i_max+4, 2, 2) * u2(i_max:i_max+4))

do i = i_max, 2, -1
    u1p(i) = Rp(i) * (Ctot(i, 1, 1)*u1(i) + Ctot(i, 1, 2)*u2(i))
    u2p(i) = Rp(i) * (Ctot(i, 2, 1)*u1(i) + Ctot(i, 2, 2)*u2(i))
    u1_tmp  = u1(i) + adams_interp_inward_implicit(u1p, i)
    u2_tmp  = u2(i) + adams_interp_inward_implicit(u2p, i)

    lambda = -9.0_dp / 24
    Delta = 1 - lambda**2 * Rp(i-1)**2 * (Ctot(i-1, 1, 2) * Ctot(i-1, 2, 1) &
        -Ctot(i-1, 1, 1) * Ctot(i-1, 2, 2))
    M(1, 1) = (1 - lambda * Rp(i-1) * Ctot(i-1, 2, 2)) / Delta
    M(2, 1) = lambda * Rp(i-1) * Ctot(i-1, 2, 1) / Delta
    M(1, 2) = lambda * Rp(i-1) * Ctot(i-1, 1, 2) / Delta
    M(2, 2) = (1 - lambda * Rp(i-1) * Ctot(i-1, 1, 1)) / Delta

    u1(i-1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i-1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
                      ! write(6,"('   i=',i4,3x,'u1=',f12.4,10x,'M1=',f12.4,3x,'M2=',f12.4,3x,'u1_tmp=',f12.4,3x,'u2_tmp=',f12.4)")i,u1(i-1),M(1,1),M(1,2),u1_tmp,u2_tmp
                      ! write(6,"('   i=',i4,3x,'u1_tmp=',e20.8,3x,'u2_tmp=',e20.8)")i,u1_tmp,u2_tmp
    if (abs(u1(i-1)) >= max_val .or. abs(u2(i-1)) >= max_val) then
        P = u1
        Q = u2
        imin = i
        return
    end if
end do
P = u1
Q = u2
imin = 1
end subroutine

real(dp) recursive function a_k(n, lambda, sigma, zeta, kappa, c) result(r)
integer, intent(in) :: n, kappa
real(dp), intent(in) :: lambda, sigma, c, zeta
if (n < 1) then
    r = 0
    call stop_error("a_k: n >= 1 required")
else
    r = c/(n*lambda) * (kappa + (n-sigma)*sigma*lambda/zeta - zeta*lambda/c**2)
    r = r * b_k(n, lambda, sigma, zeta, kappa, c)
end if
end function

real(dp) recursive function b_k(k, lambda, sigma, zeta, kappa, c) result(r)
integer, intent(in) :: k, kappa
real(dp), intent(in) :: lambda, sigma, c, zeta
integer :: n
if (k < 1) then
    r = 0
    call stop_error("b_k: n >=1 required")
else if (k == 1) then
    r = 1/(2*c) * (kappa + zeta/lambda)
else
    n = k-1
    r = 1/(2*n*lambda) * (kappa**2 - (n-sigma)**2 - zeta**2/c**2) * &
        b_k(n, lambda, sigma, zeta, kappa, c)
end if
end function

function get_asymptotic(r, E, c, kappa, zeta, n_terms, norm)
! Calculates the Dirac asymptotic with 'n_terms' terms
! n_terms ... 0, 1, 2, 3, ...
! With n_terms=0, we get the simplest asymptotic, with high n_terms, we get a
! very precise value
!
! Example (calculate 20 terms):
!    zeta = -V(i_max)*R(i_max)
!    Y = get_asymptotic(R(i_max), E, c, kappa, zeta, 20, .true.)
real(dp), intent(in) :: r, E, c, zeta
integer, intent(in) :: kappa, n_terms
logical, intent(in) :: norm ! .true. ... normalize the asymptotics
real(dp) :: get_asymptotic(2)

real(dp) :: lambda, sigma, a_term, b_term
real(dp), parameter :: norm_constant = 1e-12_dp
integer :: i
lambda = sqrt(c**2 - E**2/c**2)
sigma = E*zeta/(c**2*lambda)
a_term = 1
b_term = 0
do i = 1, n_terms
    a_term = a_term + a_k(i, lambda, sigma, zeta, kappa, c) / r**i
    b_term = b_term + b_k(i, lambda, sigma, zeta, kappa, c) / r**i
end do
a_term = a_term * sqrt((c**2+E)/(2*c**2))
b_term = b_term * sqrt((c**2-E)/(2*c**2))
get_asymptotic(1) = (a_term + b_term)
get_asymptotic(2) = (a_term - b_term)
if (norm) then
    ! Normalized, so that P(r) = norm_constant:
    get_asymptotic = get_asymptotic / get_asymptotic(1) * norm_constant
else
    ! Correct asymptotic form:
    get_asymptotic = get_asymptotic * r**sigma * exp(-lambda*r)
end if
end function

subroutine integrate_radial_dirac_r_rk4(c, kappa, Z, E, R, V, Vmid, P, Q, imax)
!integrates the Dirac eq, returns r*R, where R is the radial
!solution.
! input parameters:
! kappa .... the kappa in the Dirac eq.
! Z .... the nucleus charge in Hartree atomic units
! E .... the energy at which to integration the equation
! R .... radial grid
! V .... potential on the radial grid
! c .... speed of light
! output parameters:
! Q .... f-component of the radial dirac wave function
! P .... g-component of the radial dirac wave function
real(dp), intent(in) :: c
integer, intent(in) :: kappa
real(dp), intent(in) :: Z
real(dp), intent(in) :: E
real(dp), intent(in) :: R(:)
real(dp), intent(in) :: V(:), Vmid(:)
real(dp), intent(out) :: Q(:)
real(dp), intent(out) :: P(:)
integer, intent(out) :: imax

real(dp), parameter :: max_val = 1e6_dp
real(dp) :: y(2)
real(dp), dimension(size(R), 2, 2) :: Ctot
real(dp), dimension(size(R)-1, 2, 2) :: Cmid
real(dp), dimension(size(R)-1) :: Rmid

real(dp) :: beta, Z1
integer :: l

beta = sqrt(kappa**2-(Z/c)**2)
if (dabs(Z) .gt. 1.d-10) then
    y(1) = R(1)**beta
    y(2) = R(1)**beta * c * (beta + kappa) / Z
else
    Z1 = V(1)
    if (kappa < 0) then
        l = -kappa-1
        y(1) = +R(1)**(l+1)
        y(2) = +R(1)**(l+2) * (E + Z1) / (c*(2*l+3))
    else
        l = kappa
        y(1) = -R(1)**(l+2) * (E + Z1) / (c*(2*l+1))
        y(2) = +R(1)**(l+1)
    end if
end if

Ctot(:, 1, 1) = -kappa / R
Ctot(:, 2, 2) = +kappa / R
Ctot(:, 1, 2) = +(E-V)/c + 2*c
Ctot(:, 2, 1) = -(E-V)/c

Rmid = (R(:size(R)-1) + R(2:)) / 2
Cmid(:, 1, 1) = -kappa / Rmid
Cmid(:, 2, 2) = +kappa / Rmid
Cmid(:, 1, 2) = +(E-Vmid)/c + 2*c
Cmid(:, 2, 1) = -(E-Vmid)/c
call rk4_integrate4(R, y, Ctot, Cmid, max_val, P, Q, imax)
end subroutine


end module ! rdirac_dftatom





module reigen_dftatom

! Solves the radial Schroedinger/Dirac eigenproblem

use types_dftatom, only: dp
use ode1d_dftatom, only: integrate, normalize, parsefunction, get_n_nodes, get_min_idx
use rschroed_dftatom, only: schroed_outward_adams, schroed_inward_adams
use rdirac_dftatom, only: dirac_outward_adams, dirac_inward_adams
use utils_dftatom, only: stop_error

implicit none
private
public solve_radial_eigenproblem, integrate_rproblem_outward


contains

subroutine integrate_rproblem_outward(l, E, R, Rp, V, Z, c, relat, &
    P, Q, imax)
! Integrates the radial Schroedinger/Dirac equations outward for the given E
!
! Input parameters:
! l .... the quantum number l
! E .... the energy at which to integration the equation
! R(:) .... radial grid
! V(:) .... potential on the radial grid
! Z .... the nucleus charge
! c .... speed of light
! relat ... the relativistic mode (1, 2 or 3), see solve_radial_eigenproblem()
!           for more information
!
! Output parameters:
! P(:), Q(:) ..... The P and Q components.
! For Schroedinger equation P(r) = r*R(r), Q(r) = P'(r), for Dirac equation,
! P(r) = r*g(r), Q(r) = r*f(r), where "g" and "f" are the large and small
! components of the Dirac equation
!
! The components are not normalized.

integer, intent(in) :: l, relat
real(dp), intent(in) :: Z, E, c, R(:), Rp(:), V(:)
real(dp), intent(out) :: P(:), Q(:)
integer, intent(out) :: imax

integer :: kappa

if (relat == 0) then
    call schroed_outward_adams(l, Z, E, R, Rp, V, P, Q, imax)
else if (relat == 1) then
    call stop_error("Scalar relativistic case not implemented")
else if (relat == 2 .or. relat == 3) then
    if (relat == 3) then
        if (l == 0) then
              call stop_error("for l=0 only spin up (relat=2) is allowed")
        end if
        kappa = l
    else
        kappa = -l-1
    end if
    call dirac_outward_adams(c, kappa, Z, E, R, Rp, V, P, Q, imax)
else
    call stop_error("wrong value of relat.")
end if
end subroutine

subroutine integrate_radial_problem_inward(l, E, R, Rp, V, c, &
    relat, P, Q, imin)
! Integrates the radial Schroedinger/Dirac equations inward for the given E
!
! Input parameters:
! l .... the quantum number l
! E .... the energy at which to integration the equation
! R(:) .... radial grid
! V(:) .... potential on the radial grid
! c .... speed of light
! relat ... the relativistic mode (1, 2 or 3), see solve_radial_eigenproblem()
!           for more information
!
! Output parameters:
! P(:), Q(:) ..... The P and Q components.
! For Schroedinger equation P(r) = r*R(r), Q(r) = P'(r), for Dirac equation,
! P(r) = r*g(r), Q(r) = r*f(r), where "g" and "f" are the large and small
! components of the Dirac equation
!
! The components are not normalized.

integer, intent(in) :: l, relat
real(dp), intent(in) :: E, c, R(:), Rp(:), V(:)
real(dp), intent(out) :: P(:), Q(:)
integer, intent(out) :: imin

integer :: kappa

if (relat == 0) then
    call schroed_inward_adams(l, E, R, Rp, V, P, Q, imin)
else if (relat == 1) then
    call stop_error("Scalar relativistic case not implemented")
else if (relat == 2 .or. relat == 3) then
    if (relat == 3) then
        if (l == 0) then
              call stop_error("for l=0 only spin up (relat=2) is allowed")
        end if
        kappa = l
    else
        kappa = -l-1
    end if
    call dirac_inward_adams(c, kappa, E, R, Rp, V, P, Q, imin)
else
    call stop_error("Wrong value of relat.")
end if
end subroutine

logical function is_E_above(n, l, nods_actual)
integer, intent(in) :: n, l, nods_actual
is_E_above = nods_actual > n-l-1
end function

subroutine solve_radial_eigenproblem(ngrid,n, l, Ein, eps, max_iter, &
    R, Rp, V, Z, c, relat, perturb, Emin_init, Emax_init, converged, E, P, Q)
! Solves the radial Dirac (Schroedinger) equation and returns the eigenvalue
! (E) and normalized eigenvector (P, Q) for the given "n" and "l".
!
!    Finds the wavefunction with defined "n" and "l". The potential is "V".
!    relat ... 0 nonrelat (runge-kutta)
!              2 relat (runge-kutta) spin up
!              3 relat (runge-kutta) spin down
!    Ein gives the initial energy for the perturbation method
!    eps ... the solver halves the difference between Emin and Emax until
!            |Emax-Emin|<eps
!    max_iter ... The maximum allowed number of shooting iterations to solve
!                 the eigenproblem
!    R  ... the grid on which to solve the equation
!    Rp ... the grid derivatives on which to solve the equation
!    V  ... the potential V(r) on the grid R
!    Z is the coulomb term in the potential V=-Z/r+..., it is used for
!    the asymptotic
!    c is the speed of light in atomic units
!    perturb ... If .true., use perturbation correction (faster, but assumes
!        the potential behaves as -Z/r for r->oo), otherwise use bisection
!        (slower, but works for any potential)
!    Emin_init, Emax_init ... The range in which to find the energy. If the
!        energy is not in the range, the "converged" variable will be equal to
!        either 9 (Emin_init too big) or 10 (Emin_init too small).
!
!    Returns eigenvalue E and the normalized wavefunctions P, Q. If it doesn't
!    converge, then converged /= 0, and E is undefined (P and Q then contains
!    the latest integration of the shooting solver if available --- this is
!    useful for debugging why it did not converge).
!
!    Conceptually, the algorithm is to use bisection to converge energy (using
!    the number of nodes as the criterion) until we get close enough, so that
!    the perturbation theory starts to converge, and then use perturbation to
!    finish it. If perturbation does not converge, we need to fail over to
!    bisection.
integer, intent(in) :: ngrid, n, l, relat, max_iter
real(dp), intent(in) :: Z, R(ngrid), Rp(ngrid), V(ngrid), eps, Ein, c
logical, intent(in) :: perturb
real(dp), intent(in) :: Emin_init, Emax_init
integer, intent(out) :: converged
real(dp), intent(out) :: P(ngrid), Q(ngrid), E


real(dp) :: Emin, Emax, dE, Pr(size(R)), Qr(size(R)), factor, S
integer :: minidx, ctp, iter
logical :: isbig
integer :: nnodes
logical :: last_bisect
integer :: imin, imax
E = Ein
if (.not.(n > 0)) call stop_error("n > 0 not satisfied")
if (.not.((0 <= l).and.(l < n))) call stop_error("0 <= l < n not satisfied")

Emax = Emax_init
Emin = Emin_init
if (E > Emax .or. E < Emin) E = (Emin + Emax) / 2

iter = 0
last_bisect = .true.
do while (iter < max_iter)
    iter = iter + 1

    ! See if bisection is converged
    if (abs(Emax - Emin) < eps) then
        if (.not. last_bisect) then
            ! The perturbation theory correction was used in the last
            ! iteration and in that case, the consistent stopping criterion is
            ! to converge with abs(dE), not abs(Emax - Emin).
            ! As such we fail, because abs(Emax - Emin) is converged, but
            ! abs(dE) isn't.
            converged = 6
            return
        end if
        if (abs(Emax - Emax_init) < tiny(1._dp)) then
            ! The algorithm didn't change Emax, so Emax_init was set
            ! incorrectly.
            converged = 10
            return
        end if
        if (abs(Emin - Emin_init) < tiny(1._dp)) then
            ! The algorithm didn't change Emin, so Emin_init was set
            ! incorrectly.
            converged = 9
            return
        end if
        call integrate_rproblem_outward(l, E, R, Rp, V, &
            Z, c, relat, P, Q, imax)
        minidx = get_min_idx(P(:imax))
        if (minidx <= 0) then
            ! The wavefunction doesn't have a peak
            converged = 4
            return
        end if

        ! Trim the wavefunction after the last minimum:
        P(minidx:) = 0
        Q(minidx:) = 0

        ! To make sure the zeros from above are not counted as nodes, we
        ! substract 1 from minidx here:
        nnodes = get_n_nodes(P(:minidx-1))

        if (nnodes /= n - l - 1) then
            ! Wrong number of nodes for the converged energy
            converged = 5
            return
        end if

        exit
    end if

    ctp = find_ctp(V + l*(l+1)/(2*R**2), E)
    ! If the classical turning point is too large (or cannot be found at all),
    ! we can't use inward integration to correct the energy, so we use
    ! bisection. Also use bisection if the user requests it.
    if (.not. perturb) then
        ctp = size(R)
    else if (ctp == 0) then
        ctp = size(R)
    else if (R(ctp) / R(size(R)) > 0.5_dp) then
        ctp = size(R)
    else if (size(R) - ctp <= 10) then
        ctp = size(R)
    else if (E >= 0) then
        ! Also do bisection for positive energies, as we cannot use inward
        ! integration for these
        ctp = size(R)
    end if
    call integrate_rproblem_outward(l, E, R(:ctp), Rp(:ctp), V(:ctp), &
        Z, c, relat, P(:ctp), Q(:ctp), imax)
    nnodes = get_n_nodes(P(:imax))

                      ! write(6,"('iter=',i3,10x,'l=',i3,3x,'E=',f14.5,3x,'sizeR=',i5,3x,'ctp=',i5)")iter,l,E,size(R),ctp
                      ! write(6,"('outward P:')")
                      ! write(6,"(20f12.4)")P(:ctp)
    if (nnodes /= n-l-1 .or. ctp == size(R) .or. imax < ctp) then
        ! If the number of nodes is not correct, or we didn't manage to
        ! integrate all the way to "ctp", or if "ctp" was too large, we just
        ! use bisection:
        isbig = is_E_above(n, l, nnodes)
        if (isbig) then
            Emax = E
        else
            Emin = E
        end if
        E = (Emin + Emax) / 2
        last_bisect = .true.
        cycle
    end if


    ! Perturbation theory correction
    call integrate_radial_problem_inward(l, E, R(ctp:), Rp(ctp:), V(ctp:), &
        c, relat, Pr(ctp:), Qr(ctp:), imin)
                      ! write(6,"('inward P:')")
                      ! write(6,"(20f12.4)")Pr(ctp:)
    if (imin > 1) then
        ! The inward integration didn't integrate to the ctp
        converged = 8
        return
    end if

    ! Normalize the inward solution to match the outward one:
    factor = P(ctp) / Pr(ctp)
    if (abs(factor) > 1e9) then
        ! Normalization factor for inward/outward is too large
        converged = 7
        return
    end if
    Pr = Pr * factor
    Qr = Qr * factor

    P(ctp+1:) = Pr(ctp+1:)
    Q(ctp+1:) = Qr(ctp+1:)
    if (relat == 2 .or. relat == 3) then
        S = integrate(Rp, P**2 + Q**2)
    else
        S = integrate(Rp, P**2)
    end if
    dE = P(ctp) * (Q(ctp) - Qr(ctp)) / (2 * S)
    if (relat == 2 .or. relat == 3) then
        dE = 2 * c * dE
    end if

    ! The only stopping criterion for perturbation theory correction:
    if (abs(dE) < eps) exit

    ! We always trust the sign of dE to drive bisection
    isbig = dE < 0
    if (isbig) then
        Emax = E
    else
        Emin = E
    end if

    ! If the dE prediction is out of the trust region, we don't trust the value
    ! of dE, and we do bisection
    if (E + dE > Emax .or. E + dE < Emin) then
        E = (Emin + Emax) / 2
        last_bisect = .true.
    else
        E = E + dE
        last_bisect = .false.
    end if
end do
if (iter == max_iter) then
    ! We didn't converge in 'max_iter' iterations
    converged = 2
    return
end if

! Normalize the wavefunction:
if (relat == 0) then
    S = integrate(Rp, P**2)
else
    S = integrate(Rp, P**2 + Q**2)
end if
S = sqrt(abs(S))
if (S > 0) then
    P = P / S
    Q = Q / S
else
    ! This would happen if the function is zero, but we already check this
    ! above (converged == 4), so we fail laudly here.
    call stop_error("solve_radial_eigenproblem: zero function")
end if

converged = 0
end subroutine


integer function find_ctp(V, E) result(ctp)
! Finds the classical turning point for the potential 'V' and energy 'E'
!
! Classical turning point 'ctp' is defined as E = V(ctp)
! The function returns the integer index into the array V.
real(dp), intent(in) :: V(:), E
integer :: i
do i = size(V), 1, -1
    if (V(i)-E <= 0) then
        ctp = i
        return
    end if
end do
ctp = 0
end function

end module ! reigen






