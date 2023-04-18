module rlsy_geometry
!!
!! A collection of routines for geometry, like in high school -- not structure and molecules and stuff.
!!
use rlsy_constants, only: r8, rl_huge, rl_hugeint, rl_tol, &
                          rl_sqtol, rl_pi
use rlsy_helpers, only: rl_cross,rl_sqnorm,norm2
use rlsy_sorting, only: rl_qsort
implicit none
private
! some geometric objects
public :: rl_plane
! the public functions
public :: rl_rotation_matrix_from_axis_and_angle
public :: rl_improper_rotation_matrix_from_axis_and_angle
public :: rl_inscribed_sphere_in_box
public :: rl_bounding_sphere_of_box

!> A plane on Hessian normal form, that is $$ \hat{\mathbf{n}} \cdot \mathbf{x} + p = 0 $$ Also contains contains an orthornmal coordinate system aligned with the plane, useful for projections.
type :: rl_plane
    !> The plane normal
    real(r8), dimension(3) :: normal=rl_huge
    !> The distance from the origin
    real(r8) :: p=rl_huge
    !> first vector that spans the plane. Orthogonal to the normal and v2.
    real(r8), dimension(3) :: v1=rl_huge
    !> second vector that spans the plane. Orthogonal to the normal and v1.
    real(r8), dimension(3) :: v2=rl_huge
    contains
        !> create the plane from a set of points.
        procedure :: generate => create_plane_from_points
        !> signed distance to a point
        procedure :: distance_to_point => rl_plane_point_distance
        !> sort a set of points clockwise on this plane. Or counterclockwise, not sure. In order at least.
        procedure :: anglesort => anglesort_with_plane
end type

contains

!> Returns a rotation matrix. Defined as rotation about an axis u with angle alpha. Returns 3x3 rotation matrix. Use as `matmul(rotationmatrix,vector)`.
pure function rl_rotation_matrix_from_axis_and_angle(u,alpha) result(m)
    !> the axis
    real(r8), dimension(3), intent(in) :: u
    !> the angle, in radians
    real(r8), intent(in) :: alpha
    !> the rotation matrix
    real(r8), dimension(3,3) :: m
    !
    real(r8) :: st,ct,a,b,c,invnrm

    invnrm=norm2(u)
    ! If the vector is really tiny, I can get strange stuff.  makes little sense if the vector is supersmall
    ! not sure what to do there.
    if ( invnrm .gt. rl_tol ) then
        invnrm=1.0_r8/invnrm
        a=u(1)*invnrm
        b=u(2)*invnrm
        c=u(3)*invnrm
    else
        a=0.0_r8
        b=0.0_r8
        c=0.0_r8
    endif
    st=sin(alpha)
    ct=cos(alpha)
    ! Apparently it's called Rodriguez formula.
    m(1,1)=ct+a*a*(1.0_r8-ct)
    m(1,2)=a*b*(1.0_r8-ct)-c*st
    m(1,3)=a*c*(1.0_r8-ct)+b*st
    !
    m(2,1)=b*a*(1.0_r8-ct)+c*st
    m(2,2)=ct+b*b*(1.0_r8-ct)
    m(2,3)=b*c*(1.0_r8-ct)-a*st
    !
    m(3,1)=a*c*(1.0_r8-ct)-b*st
    m(3,2)=c*b*(1.0_r8-ct)+a*st
    m(3,3)=ct+c*c*(1.0_r8-ct)
end function

!> Returns an improper rotation matrix. Defined as rotation about an axis u with angle alpha, followed by reflection through the origin. Use as `matmul(rotationmatrix,vector)`.
pure function rl_improper_rotation_matrix_from_axis_and_angle(u,alpha) result(m)
    !> the axis
    real(r8), dimension(3), intent(in) :: u
    !> the angle, in radians
    real(r8), intent(in) :: alpha
    !> the rotation matrix
    real(r8), dimension(3,3) :: m
    !
    real(r8) :: st,ct,a,b,c,invnrm

    invnrm=norm2(u)
    ! If the vector is really tiny, I can get strange stuff.  makes little sense if the vector is supersmall
    ! not sure what to do there. Now I return NaN perhaps?
    if ( invnrm .gt. rl_tol ) then
        invnrm=1.0_r8/invnrm
        a=u(1)*invnrm
        b=u(2)*invnrm
        c=u(3)*invnrm
    else
        a=0.0_r8
        b=0.0_r8
        c=0.0_r8
    endif
    !
    st=sin(alpha)
    ct=cos(alpha)
    ! quarternion thing from wikipedia
    m(1,1)=ct-a*a*(1.0_r8+ct)
    m(1,2)=-a*b*(1.0_r8+ct)-c*st
    m(1,3)=-a*c*(1.0_r8+ct)+b*st
    !
    m(2,1)=-b*a*(1.0_r8+ct)+c*st
    m(2,2)=ct-b*b*(1.0_r8+ct)
    m(2,3)=-b*c*(1.0_r8+ct)-a*st
    !
    m(3,1)=-a*c*(1.0_r8+ct)-b*st
    m(3,2)=-c*b*(1.0_r8+ct)+a*st
    m(3,3)=ct-c*c*(1.0_r8+ct)
end function

!> given a parallelepiped, what is the radius of the largest sphere that fits inside?
pure function rl_inscribed_sphere_in_box(box) result(r)
    !> vectors that define the box
    real(r8), dimension(3,3), intent(in) :: box
    !> largest sphere that fit inside box
    real(r8) :: r

    real(r8), dimension(3) :: na,nb,nc
    ! the normals of the faces of the box
    na=rl_cross(box(:,2),box(:,3))
    nb=rl_cross(box(:,3),box(:,1))
    nc=rl_cross(box(:,1),box(:,2))
    na=na/norm2(na)
    nb=nb/norm2(nb)
    nc=nc/norm2(nc)
    ! distances between opposing planes
    r=rl_huge
    r=min(r,abs(dot_product(na,box(:,1))))
    r=min(r,abs(dot_product(nb,box(:,2))))
    r=min(r,abs(dot_product(nc,box(:,3))))
    r=r*0.5_r8
end function

!> given a parallelepiped, what is the radius of the smallest sphere that fits the box?
pure function rl_bounding_sphere_of_box(box) result(r)
    !> vectors that define the box
    real(r8), dimension(3,3), intent(in) :: box
    !> radius of smallest sphere that contains the box
    real(r8) :: r

    real(r8), dimension(3) :: a,b,c,v
    real(r8), dimension(4) :: diags

    a=box(:,1)
    b=box(:,2)
    c=box(:,3)
    ! had to do this in an awkward way to avoid segfaults on intel compilers. Odd.
    v=a-b-c
    diags(1)=dot_product(v,v)
    v=a+b-c
    diags(2)=dot_product(v,v)
    v=a-b+c
    diags(3)=dot_product(v,v)
    v=a+b+c
    diags(4)=dot_product(v,v)
    r=maxval(diags)
    r=sqrt(r)*0.5_r8
end function

!> Signed distance from plane to a point.
pure function rl_plane_point_distance(plane,point) result(r)
    !> plane
    class(rl_plane), intent(in) :: plane
    !> point
    real(r8), dimension(3), intent(in) :: point
    !> (signed) distance
    real(r8) :: r

    r=dot_product(plane%normal,point)+plane%p
end function

!> Creates a plane, either from three points or a point and a normal.
subroutine create_plane_from_points(plane,points,point,normal)
    !> The plane
    class(rl_plane), intent(out) :: plane
    !> Three points defining a plane
    real(r8), dimension(3,3), intent(in), optional :: points
    !> Point in the plane
    real(r8), dimension(3), intent(in), optional :: point
    !> Normal of the plane
    real(r8), dimension(3), intent(in), optional :: normal
    !
    real(r8), dimension(3) :: v0,v1
    real(r8), dimension(3,6), parameter :: stupidpoints=reshape([&
        1.0_r8,0.0_r8,0.0_r8, &
        0.0_r8,1.0_r8,0.0_r8, &
        0.0_r8,0.0_r8,1.0_r8, &
        0.0_r8,0.7071067811865475_r8,0.7071067811865475_r8,&
        0.7071067811865475_r8,0.0_r8,0.7071067811865475_r8,&
        0.7071067811865475_r8,0.7071067811865475_r8,0.0_r8],[3,6])
    integer :: i

    ! If it's defined from three points
    if ( present(points) ) then
        v0=points(:,1)-points(:,3)
        v1=points(:,2)-points(:,3)
        plane%normal=rl_cross(v0,v1)
        plane%normal=plane%normal/norm2(plane%normal)
        plane%p=-dot_product(plane%normal,points(:,3))
    endif

    ! Define it by a normal and a point. I don't assume that the
    ! normal is of norm 1.
    if ( present(point) .and. present(normal) ) then
        plane%normal=normal/norm2(normal)
        plane%p=-dot_product(plane%normal,point)
    endif

    ! Create the vectors that span the plane. First get a vector orthogonal
    ! To the normal. I want to keep the routine pure, so I generate some strange numbers.
    ! It's impossible that they are all parallel to the normal.
    do i=1,6
        v1=rl_cross(plane%normal,stupidpoints(:,i))
        if ( rl_sqnorm(v1) .gt. rl_tol ) exit
    enddo
    v1=v1/norm2(v1)
    ! And then a third orthogonal vector
    v0=rl_cross(v1,plane%normal)
    v0=v0/norm2(v0)
    ! I think this becoms a right handed system.
    plane%v1=v0
    plane%v2=v1
end subroutine

!> Sort points clockwise with reference to this plane. Or counterclockwise, I dunno. In order at least.
subroutine anglesort_with_plane(plane,points,order)
    !> The plane the points will get projected to
    class(rl_plane), intent(in) :: plane
    !> Points to sort
    real(r8), dimension(:,:), intent(inout) :: points
    !> The order of the points
    integer, dimension(:), intent(out), optional :: order

    integer :: i,n
    real(r8), dimension(3) :: ctr
    real(r8) :: x,y
    real(r8), dimension(size(points,2)) :: angle
    integer, dimension(size(points,2)) :: dumind

    n=size(points,2)
    ! get the center of mass of the points
    ctr=0.0_r8
    do i=1,n
        ctr=ctr+points(:,i)
    enddo
    ctr=ctr/real(n,r8)
    ! project onto plane
    do i=1,n
        x=dot_product(plane%v1,points(:,i)-ctr)
        y=dot_product(plane%v2,points(:,i)-ctr)
        angle(i)=atan2(x,y)
    enddo
    ! Get a signed angle
    call rl_qsort(angle,dumind)
    points=points(:,dumind)
    if ( present(order) ) order=dumind
end subroutine

end module
