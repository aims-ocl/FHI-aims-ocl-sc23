!****h* FHI-aims/analytic_multipole_coefficients
!  NAME
!    analytic_multipole_coefficients
!  SYNOPSIS
 
module analytic_multipole_coefficients


!  PURPOSE
!   
!    The module includes the analytic constants and their calculation routines for the Hartree potential.
!    They are used in the cluster systems for the far distance part of the potentials (outside multipole radius).
!    In periodic systems they are used for the all three far distance parts:  Fourier part and real space part
!    inside and outside multipole radius.
!
!    The method is published in:
!
!    Bernard Delley, "Fast Calculation of Electrostatics in Crystals and Large Molecules", 
!    J. Phys. Chem. 100, 6107-6110 (1996).
!
!    Here we use same notation for the constants than in the paper ( a, d and c constants).
!
!  USES

  implicit none

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  ! --- d coefficient matrix

  ! Relates spherical harmonics to derivatives of a radially symmetric
  ! function.

  real*8,  dimension(:), allocatable:: dd
  integer, dimension(:,:), allocatable:: index_dd
  integer:: n_dd_lm_ijk

  ! Y(L, M) = 0.d0
  ! do n = 1, n_dd_lm_ijk
  !   L = index_dd(n, 1)
  !   M = index_dd(n, 2)
  !   i = index_dd(n, 3)
  !   j = index_dd(n, 4)
  !   k = index_dd(n, 5)
  !   Y(L,M) = Y(L,M) + sqrt((2*l+1)/4pi) * [d^L / dx^i dy^j dz^k (1/r)]_{r=1}
  ! end do

  ! --- c coefficient matrix

  ! Relates spherical harmonics to polynomials.

  real*8,  dimension(:), allocatable:: cc
  integer, dimension(:,:), allocatable:: index_cc
  integer, dimension(:,:),allocatable:: index_ijk_max_cc
  integer, dimension(:), allocatable:: n_cc_lm_ijk

  ! Y(L, M) = 0.d0
  ! do n = 1, n_cc_lm_ijk
  !   L = index_cc(n, 1)
  !   M = index_cc(n, 2)
  !   i = index_cc(n, 3)
  !   j = index_cc(n, 4)
  !   k = index_cc(n, 5)
  !   Y(L,M) = Y(L,M) + sqrt((2*l+1)/4pi) * [x^i y^j z^k F^p(r)]_{r=1}
  ! end do
  ! F^p(r) = [(1/r) (d/dr)]^p 1/r.

  ! this is a derived dimension and only affects the present module.
  integer, private :: l_max_analytic_multipole
  integer, private, parameter :: l_max_multipole = 13
  integer, private :: n_total_cartesian 
  integer, dimension(:), allocatable:: n_total_cartesian_il

!******
contains
!---------------------------------------------------------------------
!****s* analytic_multipole_coefficients/initialize_analytic_multipole_coefficients
!  NAME
!    initialize_analytic_multipole_coefficients
!  SYNOPSIS

  subroutine initialize_analytic_multipole_coefficients


!  PURPOSE
!
!    The subroutine allocates and initializes analytic multipole coefficients 
!    for the far distance part of the hartree potential. This routine calls other
!    initialization subroutines:
!    o  initialize_F_p_functions
!    o  evaluate_coeff_a
!    o  evaluate_coeff_d
!    o  evaluate_coeff_c
!  USES
    use dimensions,            only : use_forces, l_pot_max
    use mpi_tasks,             only : check_allocation
    use Hartree_F_p_functions, only : initialize_F_p_functions
    use runtime_choices,       only : hartree_d_matrix_method
    implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    real*8,  dimension(:,:,:), allocatable :: dd_matrix
    integer, dimension(:,:,:), allocatable :: index_ijk
    integer                                :: info
    ! WPH, December 12 2017:  a was formerly defined to be an automatic array,
    ! but I ran into issues when executing MPI-compiled executables on ORNL's
    ! summitdev that I traced down to this definition.  Moving it to a dynamic 
    ! array fixed the issue.  This smells like a ulimit bug, but I know for 
    ! certain that I'm running "ulimit -s unlimited"
    integer,                   allocatable :: a(:,:,:,:,:,:)
    integer :: i_l

    ! Begin with allocations
    ! Allocations are done in the module text, as they depend on dimensions
    ! that are determined along the way

    if (use_forces) then
      l_max_analytic_multipole = l_pot_max + 1
    else
      l_max_analytic_multipole = l_pot_max
    end if

    allocate( n_total_cartesian_il(0:l_max_analytic_multipole),stat=info)
    call check_allocation(info, ' n_total_cartesian_il                     ')

    ! YY: parameters for switch between pseudo inverse and submatrix method for 
    !     D matrix evaluation
    !if (use_d_matrix_pseudo_inverse) then
    if (hartree_d_matrix_method == 2) then
      n_total_cartesian = (l_max_analytic_multipole+1)*(l_max_analytic_multipole+2)/2
      do i_l = 0, l_max_analytic_multipole,1
        n_total_cartesian_il(i_l) = (i_l+1)*(i_l+2)/2
      end do
    else if (hartree_d_matrix_method == 1) then ! normal way submatrix method by Delley
      n_total_cartesian = l_max_analytic_multipole*2+1
      do i_l = 0, l_max_analytic_multipole,1
        n_total_cartesian_il(i_l) = (i_l*2+1)
      end do
    end if

    ! globally used index lists for the later evaluation of
    ! ylm-like functions and their derivatives
    allocate( n_cc_lm_ijk(0:l_max_analytic_multipole),stat=info)
    call check_allocation(info, ' n_cc_lm_ijk                     ')

    allocate(index_ijk_max_cc(3,0:l_max_analytic_multipole),stat=info)
    call check_allocation(info, 'index_ijk_max_cc              ')

    ! locally used arrays
    allocate(dd_matrix((l_max_analytic_multipole*2+1), n_total_cartesian,0:l_max_analytic_multipole ),stat=info)
    call check_allocation(info, 'dd_matrix                     ')

    allocate(index_ijk(3,n_total_cartesian,0:l_max_analytic_multipole),stat=info)
    call check_allocation(info, 'index_ijk                     ')

    allocate(a(0:l_max_multipole,0:l_max_multipole, 0:l_max_multipole, 0:l_max_multipole, &
         0:l_max_multipole, 0:l_max_multipole), stat=info)
    call check_allocation(info, 'a                             ')

    ! initialize the radial functions which define a multipole decomposition
    call  initialize_F_p_functions

    call evaluate_coeff_a(a)
    call evaluate_coeff_d(dd_matrix,index_ijk,a) 
    call evaluate_coeff_c(dd_matrix,index_ijk,a)

    ! deallocate locally used arrays 

    deallocate(index_ijk)
    deallocate(dd_matrix)
    deallocate(a)

  end subroutine initialize_analytic_multipole_coefficients
!******
!--------------------------------------------------------------------------------------------
!****s* analytic_multipole_coefficients/evaluate_coeff_a
!  NAME
!    evaluate_coeff_a
!  SYNOPSIS

  subroutine evaluate_coeff_a(a)

!  PURPOSE
!    Evaluates so called a-coefficients for the far distance hartree potential.
!    See: Bernard Delley, "Fast Calculation of Electrostatics in Crystals and Large Molecules", 
!    J. Phys. Chem. 100, 6107-6110 (1996).
!  USES
    implicit none
!  ARGUMENTS
    integer, intent(out) :: a(0:l_max_multipole,0:l_max_multipole, 0:l_max_multipole, 0:l_max_multipole, &
         0:l_max_multipole, 0:l_max_multipole)
!  INPUTS
!    none
!  OUTPUTS
!    o a 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer:: x,y,z, xx,yy,zz, k

!    write(use_unit,*) 'Initializing coefficients a for periodic Hartree potential'
    
    a = 0
    a(0,0,0,0,0,0) = 1
!    l = 1
    a(1,0,0,1,0,0) = -1
!    a(0,1,0,0,1,0) = -1
!    a(0,0,1,0,0,1) = -1

    y = 0
    z = 0
    yy = 0
    zz = 0


    do x = 2,l_max_multipole-1,1
       do xx = x,0,-2

          k =  (x + y + z + xx + yy + zz)+1 - 2
          k = (0,1)**(k*2)*k
         
          if( xx > 0)then
             a(x,y,z , xx, yy, zz) =  a(x-1, y, z, xx+1, yy, zz) * (xx+1) &
                  + a(x-1, y, z, xx-1, yy, zz) * k 
          else
             a(x,y,z , xx, yy, zz) = a(x-1, y, z, xx+1, yy, zz) * (xx+1)
          end if

       end do
    end do


    do x = 0,l_max_multipole-1,1
       do y = 0,l_max_multipole-1,1
          do yy = y,0,-2
             do xx = x,0,-2

                k = (x + y + z + xx + yy + zz)+1 - 2
                k = (0,1)**(k*2)*k

                if(y>0)then
                   a(x,y,z , xx, yy, zz) =  a(x,y,z , xx, yy, zz) &
                        + a(x, y-1, z, xx, yy+1, zz) * (yy+1)

                   if(yy>0)then
                      a(x,y,z , xx, yy, zz) =  a(x,y,z , xx, yy, zz) &
                           + a(x, y-1, z, xx, yy-1, zz) * k 
                   end if
                end if
             end do
          end do
       end do
    end do

    
    do x = 0,l_max_multipole-1,1
       do y = 0,l_max_multipole-1,1
          do z = 1,l_max_multipole-1,1
             do zz = z,0,-2
                do yy = y,0,-2
                   do xx = x,0,-2

                      k = (x + y + z + xx + yy + zz)+1 - 2
                      k = (0,1)**(k*2)*k

                      if(z>0)then
                         a(x,y,z , xx, yy, zz) =  a(x,y,z , xx, yy, zz) &
                              + a(x, y, z-1, xx, yy, zz+1) * (zz+1)

                         if(zz>0)then
                            a(x,y,z , xx, yy, zz) =  a(x,y,z , xx, yy, zz) &
                                 + a(x, y, z-1, xx, yy, zz-1) * k 
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do

!    write(use_unit,*) '  done'


  end subroutine evaluate_coeff_a
  
!******
!-----------------------------------------------------------------------------------------
!****s* analytic_multipole_coefficients/evaluate_coeff_d
!  NAME
!    evaluate_coeff_d
!  SYNOPSIS

  subroutine evaluate_coeff_d(dd_matrix,index_ijk,a)

!  PURPOSE
!    Evaluates so called d-coefficients for the far distance hartree potential.
!    See: Bernard Delley, "Fast Calculation of Electrostatics in Crystals and Large Molecules", 
!    J. Phys. Chem. 100, 6107-6110 (1996).
!  USES
    use grids,       only : lebedev_grid_floor
    use dimensions,  only : n_max_angular
    use localorb_io, only : use_unit
    use mpi_tasks,   only : check_allocation, aims_stop
    use numerical_utilities, only: pseudo_inverse
    implicit none
! ARGUMENTS
    real*8,  intent(out) :: dd_matrix((l_max_analytic_multipole*2+1), n_total_cartesian,0:l_max_analytic_multipole)
    !real*8,  intent(out) :: dd_matrix((l_max_analytic_multipole*2+1), (l_max_analytic_multipole*2+1),0:l_max_analytic_multipole)
    integer, intent(out) :: index_ijk(3,n_total_cartesian,0:l_max_analytic_multipole)
    integer, intent(in)  :: a(0:l_max_multipole,0:l_max_multipole, 0:l_max_multipole, 0:l_max_multipole, &
         0:l_max_multipole, 0:l_max_multipole)
!  INPUTS
!    o a
!  OUTPUTS
!    o dd_matrix
!    o index_ijk
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    real*8 :: dd_matrix_copy((l_max_analytic_multipole*2+1), n_total_cartesian,0:l_max_analytic_multipole)
    real*8 :: dd_matrix_inv( n_total_cartesian,(l_max_analytic_multipole*2+1),0:l_max_analytic_multipole)
    integer:: i,j,k,n
    real*8, dimension(:,:), allocatable::  r_angular
    real*8, dimension(:), allocatable::  w_angular
    integer :: n_lebedev, n_angular, i_angular
    real*8,dimension(3):: coord_current, dir_tab
    real*8:: dist_tab, deriv, value
    real*8,dimension(4):: trigonom_tab
    real*8,dimension(:), allocatable:: ylm_tab
    integer, dimension(:,:), allocatable :: index_lm
    integer:: i_m, i_l, nn, n_max, aa,bb,cc

    real*8::maximi
    integer,dimension(:),allocatable:: ipivot
    integer:: info
    integer:: matrix_size, matrix_size_2, position
    real*8,allocatable, dimension(:,:):: unit_matrix

    integer*8,dimension(:,:),allocatable:: ijk_order
!    logical:: converged
    character(*), parameter :: func = 'evaluate_coeff_d'


    !write(use_unit,*) 'Initializing coefficients d for periodic Hartree potential'

    dd_matrix = 0.0
    allocate(ylm_tab( (l_max_analytic_multipole+1)**2),stat=info)
    call check_allocation(info, 'ylm_tab                       ')

    allocate( index_lm( -l_max_analytic_multipole:l_max_analytic_multipole,0:l_max_analytic_multipole),stat=info)
    call check_allocation(info, 'index_lm                      ')

    allocate(ipivot(l_max_analytic_multipole*2+1),stat=info)
    call check_allocation(info, 'ipivot                        ')

    allocate(unit_matrix( l_max_analytic_multipole*2+1,l_max_analytic_multipole*2+1),stat=info) 
    call check_allocation(info, 'unit_matrix                   ')

    allocate(ijk_order((l_max_analytic_multipole*2+1)*(l_max_analytic_multipole*2+1)*(l_max_analytic_multipole*2+1),4),stat=info)
    call check_allocation(info, 'ijk_order                     ')
  
    n = 0
    do i_l = 0, l_max_analytic_multipole,1
       do i_m = -i_l, i_l
          n = n+1
          index_lm(i_m,i_l) = n
       end do
    end do

    n_max = n

    ijk_order = 1e15
    index_ijk = -1
    index_ijk_max_cc = 0

    do i_l = 0, l_max_analytic_multipole,1
       nn = 0
       ijk_order = 10000
       do i = 0,i_l
          do j = 0,i_l
             do k = 0,i_l

                if( i_l== i+j+k)then

                   nn = nn + 1
                   !ijk_order(nn,1) = i**2 + j**2 + k**2
                   ijk_order(nn,1) = i + j + k
                   ijk_order(nn,2) = i
                   ijk_order(nn,3) = j
                   ijk_order(nn,4) = k

                   index_ijk_max_cc(1,i_l) = max( index_ijk_max_cc(1,i_l),i)
                   index_ijk_max_cc(2,i_l) = max( index_ijk_max_cc(2,i_l),j)
                   index_ijk_max_cc(3,i_l) = max( index_ijk_max_cc(3,i_l),k)

                end if
             end do
          end do
       end do

       do nn = 1,n_total_cartesian_il(i_l),1

          position = minloc(ijk_order(1:((i_l*2+1)*(i_l*2+1)*(i_l*2+1)),1),1)

          index_ijk(1,nn,i_l) = ijk_order( position,2)
          index_ijk(2,nn,i_l) = ijk_order( position,3)
          index_ijk(3,nn,i_l) = ijk_order( position,4)
          
          ijk_order(position,1) = 1e15

       end do
       if (any(index_ijk_max_cc(:, i_l) /= i_l)) then
          call aims_stop('Assertion on index_ijk_max_cc failed', func)
       end if
    end do



    n_lebedev  =  lebedev_grid_floor(n_max_angular)
    allocate(r_angular(3,n_max_angular),stat=info)
    call check_allocation(info, 'r_angular                     ')

    allocate(w_angular(n_max_angular),stat=info)
    call check_allocation(info, 'w_angular                     ')


    n_angular = n_max_angular
    call get_angular_grid ( n_lebedev, n_angular, r_angular(1,1),  w_angular(1) )


    dd_matrix = 0

    do i_angular = 1, n_angular
                 

       coord_current(:) =  r_angular(:, i_angular)
       
       dist_tab = 1.0
       
       dir_tab =  coord_current


       call tab_trigonom2(dir_tab, trigonom_tab, 1)
          
       call ylm_real ( dir_tab, l_max_analytic_multipole, ylm_tab(1) )


       do i_l = 0, l_max_analytic_multipole,1
          do i_m = -i_l, i_l,1
             !do n = 1, (i_l*2+1),1
             do n = 1, n_total_cartesian_il(i_l) ,1
                

                i = index_ijk(1,n,i_l)
                j = index_ijk(2,n,i_l)
                k = index_ijk(3,n,i_l)


                deriv = 0.0
                value =  w_angular(i_angular) * ylm_tab( index_lm(i_m,i_l))
                if( abs(value) < 1e-30) value = 0

                   
                do aa = i,0,-2
                   do bb= j,0,-2
                      do cc = k,0,-2



                         deriv = deriv + a(i,j,k,aa,bb,cc)* &
                              coord_current(1)**aa * coord_current(2)**bb * coord_current(3)**cc

                         if(abs(deriv) < 1e-30)  deriv = 0
                         
                      end do
                   end do
                end do

                value = value*deriv
                dd_matrix(i_l+ i_m+1 ,n,i_l) =  dd_matrix(i_l+ i_m+1 ,n,i_l) + value

             end do
          end do
       end do
    end do
    

    dd_matrix_copy = dd_matrix
    dd_matrix_inv = 0

    do i_l = 0, l_max_analytic_multipole

       matrix_size = i_l*2+1
       matrix_size_2 = n_total_cartesian_il(i_l)
       do i=1,matrix_size
          do j=1, matrix_size_2
             if(abs(dd_matrix(i,j,i_l)) < 1e-15)  dd_matrix(i,j,i_l) = 0.0
             if(abs(dd_matrix_copy(i,j,i_l)) < 1e-15)  dd_matrix_copy(i,j,i_l) = 0.0
          end do   
       end do

       !call DGETRF(matrix_size, matrix_size, dd_matrix(1:matrix_size,1:matrix_size,i_l), matrix_size, ipivot, info)
       !
       !unit_matrix = 0.0
       !do i=1, matrix_size
       !   unit_matrix(i,i) = 1.0
       !end do

       !do i=1, matrix_size

       !   call DGETRS('N', matrix_size, 1, dd_matrix(1:matrix_size,1:matrix_size,i_l) , matrix_size, &
       !        ipivot, unit_matrix(1:matrix_size,i), matrix_size, info)

       !   if( info /= 0) then
       !      write(use_unit,*) 'Error in DGETRS',info
       !      stop
       !   end if
       !end do


       !do i=1,matrix_size
       !   maximi = maxval(abs(dd_matrix(1:matrix_size,i,i_l)))

       !   do j=1, matrix_size
       !      if(abs(unit_matrix(i,j)) < 1e-10*maximi)  unit_matrix(i,j) = 0.0
       !   end do
       !end do
       !     

       !dd_matrix(1:matrix_size,1:matrix_size,i_l) =  transpose(unit_matrix(1:matrix_size,1:matrix_size))

       call pseudo_inverse(func, matrix_size, matrix_size_2, dd_matrix_copy(1:matrix_size,1:matrix_size_2,i_l), &
       dd_matrix_inv(1:matrix_size_2,1:matrix_size,i_l), 1d-15)
       do i=1,matrix_size
          maximi = maxval(abs(dd_matrix_inv(1:matrix_size_2,i,i_l)))

          do j=1, matrix_size_2
             if(abs(dd_matrix_inv(j,i,i_l)) < 1e-10*maximi)  dd_matrix_inv(j,i,i_l) = 0.0
          end do
       end do
       dd_matrix(1:matrix_size,1:matrix_size_2,i_l) = transpose(dd_matrix_inv(1:matrix_size_2,1:matrix_size,i_l))
       


    end do



    deallocate(w_angular)
    deallocate(r_angular)
    deallocate(ijk_order)
    deallocate(unit_matrix)
    deallocate(ipivot)
    deallocate(index_lm)
    deallocate(ylm_tab)


    ! Put the data in a sparse form which can be used in actual calculations

    n_dd_lm_ijk = 0

    do i_l = 0, l_max_analytic_multipole,1
       do i_m = -i_l, i_l,1
          !do n = 1, (i_l*2+1),1
          do n = 1, n_total_cartesian_il(i_l),1
             if(dd_matrix(i_l+ i_m+1 ,n,i_l) /= 0  )then

                n_dd_lm_ijk = n_dd_lm_ijk + 1

             end if
          end do
       end do
    end do


    ! only at this point are we ready to allocate dd / index_dd to their exact length ...
    allocate(dd(n_dd_lm_ijk),stat=info)
    call check_allocation(info, 'dd                            ')

    allocate(index_dd(n_dd_lm_ijk,5),stat=info)
    call check_allocation(info, 'index_dd                      ')

  

    nn = 0
    do i_l = 0, l_max_analytic_multipole,1
       do i_m = -i_l, i_l,1
          !do n = 1, (i_l*2+1),1
          do n = 1, n_total_cartesian_il(i_l),1
             if(dd_matrix(i_l+ i_m+1 ,n,i_l) /= 0  )then
             
                nn = nn + 1

                index_dd(nn, 1) = i_l
                index_dd(nn, 2) = i_m
                index_dd(nn, 3) = index_ijk(1,n,i_l)
                index_dd(nn, 4) = index_ijk(2,n,i_l)
                index_dd(nn, 5) = index_ijk(3,n,i_l)

                dd(nn) = dd_matrix(i_l+ i_m+1 ,n,i_l)

!                write(use_unit,*) dd(nn), i_l, i_m

             end if
          end do
       end do
    end do





  end subroutine evaluate_coeff_d
!******
!----------------------------------------------------------------------
!****s* analytic_multipole_coefficients/evaluate_coeff_c
!  NAME
!    evaluate_coeff_c
!  SYNOPSIS

  subroutine evaluate_coeff_c( dd_matrix,index_ijk,a)

!  PURPOSE
!    Evaluates so called c-coefficients for the far distance hartree potential.
!    See: Bernard Delley, "Fast Calculation of Electrostatics in Crystals and Large Molecules", 
!    J. Phys. Chem. 100, 6107-6110 (1996).
!  USES
    use Hartree_F_p_functions, only : F_1_div_r_coeff
    use localorb_io,           only : use_unit
    use mpi_tasks,             only : check_allocation
    implicit none
!  ARGUMENTS
    real*8,  intent(in) :: dd_matrix((l_max_analytic_multipole*2+1), n_total_cartesian,0:l_max_analytic_multipole)
    integer, intent(in) :: index_ijk(3,n_total_cartesian,0:l_max_analytic_multipole)
    integer, intent(in) :: a(0:l_max_multipole,0:l_max_multipole, 0:l_max_multipole, 0:l_max_multipole, &
         0:l_max_multipole, 0:l_max_multipole)
!  INPUTS
!    o dd_matrix
!    o index_ijk
!    o a
!  OUTPUT
!    none (modifies module variables)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    real*8,  dimension(:,:,:), allocatable:: c
!    integer, dimension(:,:),   allocatable:: index_lm
    integer:: i_l, i_lm, n, i,j,k, ii, jj, kk, nn

    
!    write(use_unit,*) 'Initializing coefficients c for periodic Hartree potential'

    allocate( c(0:l_max_analytic_multipole , 0:l_max_analytic_multipole, 0:l_max_analytic_multipole),stat=nn)
    call check_allocation(nn, 'c                             ')

!    allocate( index_lm( -l_max_analytic_multipole:l_max_analytic_multipole,0:l_max_analytic_multipole))


    ! Calculate the number of none zero  C coefficients.
    n_cc_lm_ijk = 0 

    do i_l = 0, l_max_analytic_multipole

       if(i_l/=0) n_cc_lm_ijk(i_l) = n_cc_lm_ijk(i_l-1)

       do i_lm = 1, (i_l*2+1),1

          ! C depends on l, m, i, j, k
          ! this temp-memory is allocated only i,j,k.
          c = 0.0

          !do n = 1, (i_l*2+1),1
          do n = 1, n_total_cartesian_il(i_l),1

             i = index_ijk(1,n,i_l)
             j = index_ijk(2,n,i_l)
             k = index_ijk(3,n,i_l)

             do ii = i,0,-2
                do jj = j,0,-2
                   do kk = k,0,-2


                      c(ii,jj,kk)  =   c(ii,jj,kk) +  dd_matrix(i_lm, n, i_l) *  (a(i,j,k, ii,jj,kk) & 
                                                      /  F_1_div_r_coeff((i+j+k+ii+jj+kk)/2))

                   end do
                end do
             end do
          end do

          do ii = 0,l_max_analytic_multipole, 1
             do jj = 0,l_max_analytic_multipole, 1 
                do kk = 0,l_max_analytic_multipole, 1


                   if( c(ii,jj,kk)  /= 0 ) then

                      ! number of none zero coefficients
                      n_cc_lm_ijk(i_l) = n_cc_lm_ijk(i_l) + 1
                      !write(use_unit,*) c(ii,jj,kk), i_l, i_lm - i_l - 1, ii, jj, kk
                      !write(use_unit,*) i_l, i_lm - i_l - 1, ii, jj, kk

                   end if
                   !if ((i_l .eq. 4) .and. (i_lm -i_l -1 .eq.2)) then
                   !   write(use_unit,*) c(ii,jj,kk)
                   !end if
                end do
             end do
          end do
       end do
    end do


    ! Only at this point are we ready to dimension the cc / index_cc arrays correctly
    allocate(cc(n_cc_lm_ijk(l_max_analytic_multipole)),stat=nn)
    call check_allocation(nn, 'cc                            ')


    allocate(index_cc(n_cc_lm_ijk(l_max_analytic_multipole),6),stat=nn)
    call check_allocation(nn, 'index_cc                      ')


    ! Same again, but this time the actual data is saved.
    nn = 0

    do i_l = 0, l_max_analytic_multipole
       do i_lm = 1, (i_l*2+1),1

          c = 0.0

          !do n = 1, (i_l*2+1),1
          do n = 1, n_total_cartesian_il(i_l),1

             i = index_ijk(1,n,i_l)
             j = index_ijk(2,n,i_l)
             k = index_ijk(3,n,i_l)

             do ii = i,0,-2
                do jj = j,0,-2
                   do kk = k,0,-2


                      c(ii,jj,kk)  =   c(ii,jj,kk) +  dd_matrix(i_lm, n, i_l) *  (a(i,j,k, ii,jj,kk) & 
                                                      /  F_1_div_r_coeff((i+j+k+ii+jj+kk)/2))


                   end do
                end do
             end do
          end do

          do ii = 0,l_max_analytic_multipole, 1
             do jj = 0,l_max_analytic_multipole, 1 
                do kk = 0,l_max_analytic_multipole, 1

                   if( c(ii,jj,kk)  /= 0 ) then

                      nn = nn + 1
                      cc(nn) = c(ii,jj,kk)
                      index_cc(nn, 1) = i_l                 ! l
                      index_cc(nn, 2) = i_lm - i_l - 1      ! m
                      index_cc(nn, 3) = ii                  ! i
                      index_cc(nn, 4) = jj                  ! j
                      index_cc(nn, 5) = kk                  ! k
                      index_cc(nn, 6) = (i+j+k+ii+jj+kk)/2  ! p for Fp functions

                      !write(use_unit,*) cc(nn), i_l, i_lm - i_l - 1 , ii, jj, kk
                      !write(use_unit,*) i_l, i_lm - i_l - 1 , ii, jj, kk

                   end if
                   !if ((i_l .eq. 4) .and. (i_lm -i_l -1 .eq.2)) then
                   !   write(use_unit,*) c(ii,jj,kk)
                   !end if
                end do
             end do
          end do
       end do
    end do

!    deallocate( index_lm )
    deallocate( c )




  end subroutine evaluate_coeff_c
!******
!----------------------------------------------------------------------
!****s* analytic_multipole_coefficients/cleanup_analytic_multipole_coefficients
!  NAME
!    cleanup_analytic_multipole_coefficients
!  SYNOPSIS

  subroutine cleanup_analytic_multipole_coefficients

!  PURPOSE
!    Deallocated the variables of the module analytic_multipole_coefficients
!  USES
    implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none     
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



    if (allocated(dd)) then
       deallocate(dd)
    end if
    if (allocated(index_dd)) then
       deallocate(index_dd)
    end if
    if (allocated(cc)) then
       deallocate(cc)
    end if
    if (allocated(index_cc)) then
       deallocate(index_cc)
    end if
    if (allocated(n_cc_lm_ijk)) then
      deallocate(n_cc_lm_ijk)
    end if
    if (allocated(index_ijk_max_cc)) then
      deallocate(index_ijk_max_cc)
    end if
    if (allocated(n_total_cartesian_il)) then
      deallocate(n_total_cartesian_il)
    end if
    
  end subroutine cleanup_analytic_multipole_coefficients
!******
!----------------------------------------------------------------------
end module analytic_multipole_coefficients
