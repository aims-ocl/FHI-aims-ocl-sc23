module rlsy_helpers
!! Somewhat practical functions to use every now and then
use rlsy_constants, only: r8,rl_exitcode_blaslapack,rl_exitcode_param,rl_status

implicit none
private
public :: tochar
public :: rl_invert3x3matrix
public :: rl_cross
public :: rl_sqnorm
public :: rl_chop
public :: rl_clean_fractional_coordinates
public :: rl_determ
public :: rl_mom_real
public :: rl_mom_complex
public :: rl_mean
public :: rl_stddev
public :: norm2

! Functions that slice away small numbers
interface rl_chop
    module procedure rl_chop_real
    module procedure rl_chop_complex
end interface

! convert things to character
interface tochar
    module procedure tochar_int
    module procedure tochar_intarr
    module procedure tochar_real
    module procedure tochar_realarr
end interface

! determinant of 3x3 matrices
interface rl_determ
    module procedure rl_determ_int
    module procedure rl_determ_real
end interface rl_determ

  interface norm2
     module procedure norm2_vector, norm2_matrix, norm2_matrix_1
  end interface norm2

! type that can contain matrices of matrices. mom is matrix of matrices.
type rl_mom_real
    real(r8), dimension(:,:), allocatable :: m
end type
type rl_mom_complex
    complex(r8), dimension(:,:), allocatable :: m
end type

contains

!> Squared norm
pure function rl_sqnorm(a) result(nrm)
    !> vector
    real(r8), dimension(3), intent(in) :: a
    !> squared norm
    real(r8) :: nrm

    nrm=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
end function

!> Cross product
pure function rl_cross(b,c) result(a)
    !> first vector
    real(r8), dimension(3), intent(in) :: b
    !> second vector
    real(r8), dimension(3), intent(in) :: c
    !> cross product
    real(r8), dimension(3) :: a

    a(1)=b(2)*c(3)-b(3)*c(2)
    a(2)=b(3)*c(1)-b(1)*c(3)
    a(3)=b(1)*c(2)-b(2)*c(1)
end function

!> Invert 3x3 matrix analytically. No safety checks. Use wisely.
pure function rl_invert3x3matrix(m) result(n)
    !> 3x3 matrix to invert
    real(r8), dimension(3,3), intent(in) :: m
    !> inverse
    real(r8), dimension(3,3) :: n

    real(r8) :: det

    det =  m(1,1)*m(2,2)*m(3,3) - m(1,1)*m(2,3)*m(3,2)&
         - m(1,2)*m(2,1)*m(3,3) + m(1,2)*m(2,3)*m(3,1)&
         + m(1,3)*m(2,1)*m(3,2) - m(1,3)*m(2,2)*m(3,1)
    det=1.0_r8/det
    ! Calculate the inverse of the matrix
    n(1,1)=+det*(m(2,2)*m(3,3)-m(2,3)*m(3,2))
    n(2,1)=-det*(m(2,1)*m(3,3)-m(2,3)*m(3,1))
    n(3,1)=+det*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
    n(1,2)=-det*(m(1,2)*m(3,3)-m(1,3)*m(3,2))
    n(2,2)=+det*(m(1,1)*m(3,3)-m(1,3)*m(3,1))
    n(3,2)=-det*(m(1,1)*m(3,2)-m(1,2)*m(3,1))
    n(1,3)=+det*(m(1,2)*m(2,3)-m(1,3)*m(2,2))
    n(2,3)=-det*(m(1,1)*m(2,3)-m(1,3)*m(2,1))
    n(3,3)=+det*(m(1,1)*m(2,2)-m(1,2)*m(2,1))
end function

!> Determinant of a 3x3 matrix, doubles
pure function rl_determ_real(a) result(det)
    !> 3x3 matrix
    real(r8), dimension(3,3), intent(in) :: a
    !> the determinant
    real(r8) :: det
    !
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
        + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
end function

!> Determinant of a 3x3 matrix, integers
pure function rl_determ_int(a) result(det)
    !> 3x3 matrix
    integer, dimension(3,3), intent(in) :: a
    !> the determinant
    integer :: det
    !
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
        + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
end function

! Clean up fractional coordinates in a way that I like, so that they become 0-1, inclusive 0 not inclusive 1. Some
! quick testing suggested this was significantly faster than mod or anint if I am close to 0-1, which almost always
! is the case. Also means that you should never use this for large numbers.
elemental function rl_clean_fractional_coordinates(x,tol) result(y)
    !> unclean coordinate
    real(r8), intent(in) :: x
    !> tolerance
    real(r8), intent(in), optional :: tol
    !> the clean coordinate
    real(r8) :: y

    real(r8) :: dl
    if ( present(tol) ) then
        dl=tol
    else
        dl=1E-13_r8
    endif
    y=x
    do
        if ( abs(y) .lt. dl ) y=0.0_r8
        if ( abs(y-1.0_r8) .lt. dl ) y=0.0_r8
        if ( y .gt. 1.0_r8 ) y=y-1.0_r8
        if ( y .lt. 0.0_r8 ) y=y+1.0_r8
        if ( y .gt. -dl .and. y .lt. 1.0_r8 ) exit
    enddo
end function

!> Chop number to a short list of well-defined numbers, within a tolerance. Contains all n/m with m from 1 to 10 and n<m, and sqrt(3)/2. No problem adding more if that is needed. Only tests between -1 and 1
elemental function rl_chop_real(x,tol) result(y)
    !> number to be chopped
    real(r8), intent(in) :: x
    !> tolerance
    real(r8), intent(in) :: tol
    !> chopped number
    real(r8) :: y

    integer :: i
    real(r8), dimension(69), parameter :: welldefined=[&
          0.0000000000000000_r8,&
          0.1000000000000000_r8,&
         -0.1000000000000000_r8,&
          0.1111111111111111_r8,&
         -0.1111111111111111_r8,&
          0.1250000000000000_r8,&
         -0.1250000000000000_r8,&
          0.1428571428571428_r8,&
         -0.1428571428571428_r8,&
          0.1666666666666667_r8,&
         -0.1666666666666667_r8,&
          0.2000000000000000_r8,&
         -0.2000000000000000_r8,&
          0.2222222222222222_r8,&
         -0.2222222222222222_r8,&
          0.2500000000000000_r8,&
         -0.2500000000000000_r8,&
          0.2857142857142857_r8,&
         -0.2857142857142857_r8,&
          0.3000000000000000_r8,&
         -0.3000000000000000_r8,&
          0.3333333333333333_r8,&
         -0.3333333333333333_r8,&
          0.3750000000000000_r8,&
         -0.3750000000000000_r8,&
          0.4000000000000000_r8,&
         -0.4000000000000000_r8,&
          0.4285714285714285_r8,&
         -0.4285714285714285_r8,&
          0.4444444444444444_r8,&
         -0.4444444444444444_r8,&
          0.5000000000000000_r8,&
         -0.5000000000000000_r8,&
          0.5555555555555556_r8,&
         -0.5555555555555556_r8,&
          0.5714285714285714_r8,&
         -0.5714285714285714_r8,&
          0.6000000000000000_r8,&
         -0.6000000000000000_r8,&
          0.6250000000000000_r8,&
         -0.6250000000000000_r8,&
          0.6666666666666667_r8,&
         -0.6666666666666667_r8,&
          0.7000000000000000_r8,&
         -0.7000000000000000_r8,&
          0.7071067811865475_r8,&
         -0.7071067811865475_r8,&
          0.7142857142857143_r8,&
         -0.7142857142857143_r8,&
          0.7500000000000000_r8,&
         -0.7500000000000000_r8,&
          0.7777777777777778_r8,&
         -0.7777777777777778_r8,&
          0.8000000000000000_r8,&
         -0.8000000000000000_r8,&
          0.8333333333333334_r8,&
         -0.8333333333333334_r8,&
          0.8571428571428571_r8,&
         -0.8571428571428571_r8,&
          0.8660254037844386_r8,&
         -0.8660254037844386_r8,&
          0.8750000000000000_r8,&
         -0.8750000000000000_r8,&
          0.8888888888888888_r8,&
         -0.8888888888888888_r8,&
          0.9000000000000000_r8,&
         -0.9000000000000000_r8,&
          1.0000000000000000_r8,&
         -1.0000000000000000_r8]

    ! my list of well defined values. Could very well grow a little bit.
    ! starts with numbers that are common in symmetry operations
    y=x
    do i=1,69
        if ( abs(y-welldefined(i)) .lt. tol ) then
            y=welldefined(i)
            return
        endif
    enddo
end function

!> cutoff at small numbers for complex numbers
elemental function rl_chop_complex(x,tol) result(y)
    !> number to be chopped
    complex(r8), intent(in) :: x
    !> tolerance
    real(r8), intent(in) :: tol
    !> chopped number
    complex(r8) :: y

    real(r8) :: re,im
    re=real(x)
    im=aimag(x)
    if ( abs(re) .lt. tol ) re=0.0_r8
    if ( abs(im) .lt. tol ) im=0.0_r8
    y=cmplx(re,im,r8)
end function

!> convert one int to a character
pure function tochar_int(i,padding) result(s)
    !> integer to convert
    integer, intent(in) :: i
    !> pad the integer? Positive number means zer-padding, negative means pad with whitespace
    integer, intent(in), optional :: padding
    !> resulting string
    character(len=:), allocatable :: s

    character(len=range(i)) :: tmp
    character(len=range(i)+5) :: ttmp
    integer :: j,k

    if ( present(padding) ) then
        write(tmp,'(i0)') i ! get i to a string
        ! build a string that contains the padding
        do j=1,range(i)+5
            if ( padding > 0 ) then
                ttmp(j:j)='0'
            else
                ttmp(j:j)=' '
            endif
        enddo
        ! wrap it together to a nice string
        k=len(trim(adjustl(tmp))) ! how many digits where relevant
        s=ttmp(1:abs(padding)-k)//trim(adjustl(tmp))
    else
        ! much easier
        write(tmp,'(i0)') i
        s=trim(adjustl(tmp))
    endif
end function

!> convert array of integers to a character
pure function tochar_intarr(i,padding) result(s)
    integer, dimension(:), intent(in) :: i
    integer, intent(in), optional :: padding
    character(len=:), allocatable :: s

    character(len=size(i,1)*100) :: tmp
    integer :: j
    tmp=''
    do j=1,size(i,1)
        if ( present(padding) ) then
            tmp=trim(tmp)//' '//tochar_int(i(j),padding)//' '
        else
            tmp=trim(tmp)//' '//tochar_int(i(j))//' '
        endif
    enddo
    s=trim(tmp)
end function

!> convert a real to a character
pure function tochar_real(f,ndecimals,frmt) result(s)
    !> number to convert
    real(r8), intent(in) :: f
    !> how many decimals
    integer, intent(in), optional :: ndecimals
    !> maybe format instead
    character(len=*), intent(in), optional :: frmt
    !> resulting string
    character(len=:), allocatable :: s

    integer :: ndec
    character(len=100) :: tmp

    ! how many decimal places
    if ( present(ndecimals) ) then
        ndec=ndecimals
    else
        ndec=5
    endif

    ! convert to string
    if ( present(frmt) ) then
        write(tmp,trim(frmt)) f
        s=trim(tmp)
    else
        write(tmp,"(F30."//tochar_int(ndec)//")") f
        s=trim(adjustl(tmp))
    endif
end function

!> convert an array of reals to a character
pure function tochar_realarr(f,ndecimals,frmt) result(s)
    !> number to convert
    real(r8), dimension(:), intent(in) :: f
    !> how many decimals
    integer, intent(in), optional :: ndecimals
    !> maybe format instead
    character(len=*), intent(in), optional :: frmt
    !> resulting string
    character(len=:), allocatable :: s

    integer :: j,ndec
    character(len=size(f,1)*50) :: tmp

    ! how many decimal places
    if ( present(ndecimals) ) then
        ndec=ndecimals
    else
        ndec=5
    endif

    tmp=''
    do j=1,size(f,1)
        if ( present(frmt) ) then
            tmp=trim(tmp)//tochar_real(f(j),frmt=frmt)
        else
            tmp=trim(tmp)//' '//tochar_real(f(j),ndec)//' '
        endif
    enddo
    if ( present(frmt) ) then
        s=trim(tmp)
    else
        s=trim(adjustl(tmp))
    endif
end function

!> mean
pure function rl_mean(x) result(m)
    !> values
    real(r8), dimension(:), intent(in) :: x
    !> mean
    real(r8) :: m
    m=sum(x)/real(size(x),r8)
end function

!> Standard deviation
pure function rl_stddev(x) result(s)
    !> values
    real(r8), dimension(:), intent(in) :: x
    !> standard deviation
    real(r8) :: s

    real(r8) :: mean, sigma
    integer :: n,i

    n=size(x)
    mean=sum(x)/real(n,r8)
    sigma=0.0_r8
    do i=1,n
        sigma=sigma+(mean-x(i))**2
    enddo
    s=sqrt(sigma/real(n,r8))
end function


  !!  FUNCTION
  !!
  !!  Reimplementation of the norm2 intrinsic (F2008).
  !!
  pure real(r8) function norm2_vector(x) result(y)
    real(r8), intent(in) :: x(:)
    y = sqrt(sum(x**2))
  end function norm2_vector

  pure function norm2_matrix_1(x) result(yy)
    real(r8), intent(in) :: x(:,:)
    ! Note that here dim can only be 1, so this is not exactly the
    ! same as the intrinsic function.
    integer :: dim
    real(r8) :: y(size(x,2))
    real(r8) :: yy
    integer :: i
    dim = 1
    if (dim == 1) y = [(norm2_vector(x(:,i)), i=1,size(x,2))]
    yy = sqrt(sum(y**2))
  end function norm2_matrix_1

  pure function norm2_matrix(x, dim) result(y)
    real(r8), intent(in) :: x(:,:)
    ! Note that here dim can only be 1, so this is not exactly the
    ! same as the intrinsic function.
    integer, intent(in) :: dim
    real(r8) :: y(size(x,2))
    integer :: i
    if (dim == 1) y = [(norm2_vector(x(:,i)), i=1,size(x,2))]
  end function norm2_matrix

end module rlsy_helpers
