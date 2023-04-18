!****s* FHI-aims/get_angular_grid
!  NAME
!    get_angular_grid
!  SYNOPSIS

subroutine get_angular_grid &
     ( i_lebedev, n_angular, r_angular, w_angular )

!  PURPOSE
!  subroutine get_angular_grid obtains a Lebedev angular integration grid.
!
!  This subroutine is not part of grids.f as we do not use any modules / module data
!  in this version - rather, all arguments are passed along explicitly.
!
!  Note:
!  This subroutine returns either an original Lebedev-Laikov grid
!  (of which there are many but not so accurate) or a corrected
!  Delley-style grid (decidedly more accurate, but some possible grids are
!  missing).
!  The handling of these two is currently atrocious - need to look at the
!  entire angular grid infrastructure again and make consistent.
!
!  USES

  use dimensions,      only : n_max_angular
  use runtime_choices, only : force_lebedev
  use localorb_io,     only : use_unit
  implicit none
!  ARGUMENTS

  integer i_lebedev
  integer n_angular
  real*8 r_angular(3,n_max_angular)
  real*8 w_angular(n_max_angular)

!  INPUTS
!   o i_lebedev -- index of lebedev grid
!
!  OUTPUT
!   o n_angular -- number of angular grid points on this shell
!   o r_angular -- location of each "angular" grid point on unit sphere
!   o w_angular -- weights of the integration grid points
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





  !  local variables

  double precision x(n_max_angular)
  double precision y(n_max_angular)
  double precision z(n_max_angular)
  double precision w(n_max_angular)

  integer n_leb

  !  counters

  integer i_ang

  !  begin work

  !  Look for correct Lebedev grid. If the number of grid points
  !  exceeds the preset maximum, only set n_angular, but do not touch
  !  any arrays. This verification must be done in a fine-grained way,
  !  unfortunately.

  !  There are two versions of Lebedev's grids in this code. The first is
  !  due to the original paper of Lebedev-Laikov, but it is accurate to
  !  only 12 digits, and has a few flaws. The second version is due to
  !  Bernard Delley, is accurate to 17 grids, and flaws are corrected. We
  !  use the second version unless the first version has denser grids
  !  available.
  if (force_lebedev.eq.2) then
      if (i_lebedev.eq.0) then
          n_angular = 0
      else if (i_lebedev.eq.1) then
          n_angular = 50
          call estd00050(x,y,z,w)
      else if (i_lebedev.eq.2) then
          n_angular = 62
          call estd00062(x,y,z,w)
      else if (i_lebedev.eq.3) then
          n_angular = 72
          call estd00072(x,y,z,w)
      else if (i_lebedev.eq.4) then
          n_angular = 86
          call estd00086(x,y,z,w)
      else if (i_lebedev.eq.5) then
          n_angular = 98
          call estd00098(x,y,z,w)
      else if (i_lebedev.eq.6) then
          n_angular = 114
          call estd00114(x,y,z,w)
      else if (i_lebedev.eq.7) then
          n_angular = 128
          call estd00128(x,y,z,w)
      else if (i_lebedev.eq.8) then
          n_angular = 146
          call estd00146(x,y,z,w)
      else if (i_lebedev.eq.9) then
          n_angular = 163
          call estd00163(x,y,z,w)
      else if (i_lebedev.eq.10) then
          n_angular = 182
          call estd00182(x,y,z,w)
      else if (i_lebedev.eq.11) then
          n_angular = 201
          call estd00201(x,y,z,w)
      else if (i_lebedev.eq.12) then
          n_angular = 222
          call estd00222(x,y,z,w)
      else if (i_lebedev.eq.13) then
          n_angular = 243
          call estd00243(x,y,z,w)
      else if (i_lebedev.eq.14) then
          n_angular = 266
          call estd00266(x,y,z,w)
      else if (i_lebedev.eq.15) then
          n_angular = 289
          call estd00289(x,y,z,w)
      else if (i_lebedev.eq.16) then
          n_angular = 314
          call estd00314(x,y,z,w)
      else if (i_lebedev.eq.17) then
          n_angular = 339
          call estd00339(x,y,z,w)
      else if (i_lebedev.eq.18) then
          n_angular = 366
          call estd00366(x,y,z,w)
      else if (i_lebedev.eq.19) then
          n_angular = 393
          call estd00393(x,y,z,w)
      else if (i_lebedev.eq.20) then
          n_angular = 422
          call estd00422(x,y,z,w)
      else if (i_lebedev.eq.21) then
          n_angular = 451
          call estd00451(x,y,z,w)
      else if (i_lebedev.eq.22) then
          n_angular = 482
          call estd00482(x,y,z,w)
      else if (i_lebedev.eq.23) then
          n_angular = 513
          call estd00513(x,y,z,w)
      else if (i_lebedev.eq.24) then
          n_angular = 546
          call estd00546(x,y,z,w)
      else if (i_lebedev.eq.25) then
          n_angular = 579
          call estd00579(x,y,z,w)
      else if (i_lebedev.eq.26) then
          n_angular = 614
          call estd00614(x,y,z,w)
      else if (i_lebedev.eq.27) then
          n_angular = 649
          call estd00649(x,y,z,w)
      else if (i_lebedev.eq.28) then
          n_angular = 686
          call estd00686(x,y,z,w)
      else if (i_lebedev.eq.29) then
          n_angular = 723
          call estd00723(x,y,z,w)
      else if (i_lebedev.eq.30) then
          n_angular = 762
          call estd00762(x,y,z,w)
      else if (i_lebedev.eq.31) then
          n_angular = 801
          call estd00801(x,y,z,w)
      else if (i_lebedev.eq.32) then
          n_angular = 842
          call estd00842(x,y,z,w)
      else if (i_lebedev.eq.33) then
          n_angular = 883
          call estd00883(x,y,z,w)
      else if (i_lebedev.eq.34) then
          n_angular = 926
          call estd00926(x,y,z,w)
      else if (i_lebedev.eq.35) then
          n_angular = 969
          call estd00969(x,y,z,w)
      else if (i_lebedev.eq.36) then
          n_angular = 1014
          call estd01014(x,y,z,w)
      else if (i_lebedev.eq.37) then
          n_angular = 1059
          call estd01059(x,y,z,w)
      else if (i_lebedev.eq.38) then
          n_angular = 1106
          call estd01106(x,y,z,w)
      else if (i_lebedev.eq.39) then
          n_angular = 1153
          call estd01153(x,y,z,w)
      else if (i_lebedev.eq.40) then
          n_angular = 1202
          call estd01202(x,y,z,w)
      else if (i_lebedev.eq.41) then
          n_angular = 1251
          call estd01251(x,y,z,w)
      else if (i_lebedev.eq.42) then
          n_angular = 1302
          call estd01302(x,y,z,w)
      else if (i_lebedev.eq.43) then
          n_angular = 1862
          call estd01862(x,y,z,w)
      else if (i_lebedev.eq.44) then
          n_angular = 3282
          call estd03282(x,y,z,w)
      else if (i_lebedev.eq.45) then
          n_angular = 5102
          call estd05102(x,y,z,w)
      else if (i_lebedev.eq.46) then
          n_angular = 7322
          call estd07322(x,y,z,w)
      else if (i_lebedev.eq.47) then
          n_angular = 9942
          call estd09942(x,y,z,w)
      else if (i_lebedev.eq.48) then
          n_angular = 12962
          call estd12962(x,y,z,w)
      else if (i_lebedev.eq.49) then
          n_angular = 16382
          call estd16382(x,y,z,w)
      end if

     if (n_angular.le.n_max_angular) then
        do i_ang = 1, n_angular, 1
           r_angular (1, i_ang) = x(i_ang)
           r_angular (2, i_ang) = y(i_ang)
           r_angular (3, i_ang) = z(i_ang)
           w_angular (i_ang)    = w(i_ang)
        enddo
     end if

  else if ((force_lebedev.eq.1).or.(i_lebedev.gt.20)) then
     !       use the historic Lebedev-Laikov grids from subroutine
     !       Lebedev-Laikov.f

     if (i_lebedev.eq.0) then
        n_angular = 0
     else if (i_lebedev.eq.1) then
        n_angular = 6
        if (n_angular.le.n_max_angular) then
           call LD0006(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.2) then
        n_angular = 14
        if (n_angular.le.n_max_angular) then
           call LD0014(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.3) then
        n_angular = 26
        if (n_angular.le.n_max_angular) then
           call LD0026(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.4) then
        n_angular = 38
        if (n_angular.le.n_max_angular) then
           call LD0038(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.5) then
        n_angular = 50
        if (n_angular.le.n_max_angular) then
           call LD0050(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.6) then
        n_angular = 78
        if (n_angular.le.n_max_angular) then
           ! check external/HeoXu.f for explanation of this grid. YY
           call HX0078(x,y,z,w,n_leb) 
        end if
     else if (i_lebedev.eq.7) then
        n_angular = 86
        if (n_angular.le.n_max_angular) then
           call LD0086(x,y,z,w,n_leb)
        end if
     !else if ((i_lebedev.eq.6).or.(i_lebedev.eq.7)) then
     !   !VB Patch until we can figure out what is wrong with this grid
     !   !          n_angular = 74
     !   !          if (n_angular.le.n_max_angular) then
     !   !            call LD0074(x,y,z,w,n_leb)
     !   !          end if
     !   i_lebedev = 7
     !   n_angular = 86
     !   if (n_angular.le.n_max_angular) then
     !      call LD0086(x,y,z,w,n_leb)
     !   end if
     !   !        else if (i_lebedev.eq.7) then
     !   !          n_angular = 86
     !   !          if (n_angular.le.n_max_angular) then
     !   !            call LD0086(x,y,z,w,n_leb)
     !   !          end if
     else if (i_lebedev.eq.8) then
        n_angular = 110
        if (n_angular.le.n_max_angular) then
           call LD0110(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.9) then
        n_angular = 146
        if (n_angular.le.n_max_angular) then
           call LD0146(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.10) then
        n_angular = 170
        if (n_angular.le.n_max_angular) then
           call LD0170(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.11) then
        n_angular = 194
        if (n_angular.le.n_max_angular) then
           call LD0194(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.12) then
        n_angular = 248
        if (n_angular.le.n_max_angular) then
           ! check external/HeoXu.f for explanation of this grid. YY
           call HX0248(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.13) then
        n_angular = 288
        if (n_angular.le.n_max_angular) then
           ! check external/HeoXu.f for explanation of this grid. YY
           call HX0288(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.14) then
        n_angular = 302
        if (n_angular.le.n_max_angular) then
           call LD0302(x,y,z,w,n_leb)
        end if
     !else if ( (i_lebedev.eq.12) .or. &
     !     (i_lebedev.eq.13) .or. &
     !     (i_lebedev.eq.14) ) then
     !   !VB Patch until we can figure out what is wrong with these grids
     !   !          n_angular = 230
     !   !          if (n_angular.le.n_max_angular) then
     !   !            call LD0230(x,y,z,w,n_leb)
     !   !          end if
     !   !        else if (i_lebedev.eq.13) then
     !   !          n_angular = 266
     !   !          if (n_angular.le.n_max_angular) then
     !   !            call LD0266(x,y,z,w,n_leb)
     !   !          end if
     !   i_lebedev = 14
     !   n_angular = 302
     !   if (n_angular.le.n_max_angular) then
     !      call LD0302(x,y,z,w,n_leb)
     !   end if
     !   !        else if (i_lebedev.eq.14) then
     !   !          n_angular = 302
     !   !          if (n_angular.le.n_max_angular) then
     !   !            call LD0302(x,y,z,w,n_leb)
     !   !          end if
     else if (i_lebedev.eq.15) then
        n_angular = 350
        if (n_angular.le.n_max_angular) then
           call LD0350(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.16) then
        n_angular = 434
        if (n_angular.le.n_max_angular) then
           call LD0434(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.17) then
        n_angular = 590
        if (n_angular.le.n_max_angular) then
           call LD0590(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.18) then
        n_angular = 770
        if (n_angular.le.n_max_angular) then
           call LD0770(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.19) then
        n_angular = 974
        if (n_angular.le.n_max_angular) then
           call LD0974(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.20) then
        n_angular = 1202
        if (n_angular.le.n_max_angular) then
           call LD1202(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.21) then
        n_angular = 1454
        if (n_angular.le.n_max_angular) then
           call LD1454(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.22) then
        n_angular = 1730
        if (n_angular.le.n_max_angular) then
           call LD1730(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.23) then
        n_angular = 2030
        if (n_angular.le.n_max_angular) then
           call LD2030(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.24) then
        n_angular = 2354
        if (n_angular.le.n_max_angular) then
           call LD2354(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.25) then
        n_angular = 2702
        if (n_angular.le.n_max_angular) then
           call LD2702(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.26) then
        n_angular = 3074
        if (n_angular.le.n_max_angular) then
           call LD3074(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.27) then
        n_angular = 3470
        if (n_angular.le.n_max_angular) then
           call LD3470(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.28) then
        n_angular = 3890
        if (n_angular.le.n_max_angular) then
           call LD3890(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.29) then
        n_angular = 4334
        if (n_angular.le.n_max_angular) then
           call LD4334(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.30) then
        n_angular = 4802
        if (n_angular.le.n_max_angular) then
           call LD4802(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.eq.31) then
        n_angular = 5294
        if (n_angular.le.n_max_angular) then
           call LD5294(x,y,z,w,n_leb)
        end if
     else if (i_lebedev.ge.32) then
        n_angular = 5810
        if (n_angular.le.n_max_angular) then
           call LD5810(x,y,z,w,n_leb)
        end if
     else
        write(use_unit,*) "* Incorrect angular integration grid requested."
        write(use_unit,*) "* Check call to get_angular_grid."
        stop
     end if

     if (n_angular.le.n_max_angular) then
        do i_ang = 1, n_angular, 1
           r_angular (1, i_ang) = x(i_ang)
           r_angular (2, i_ang) = y(i_ang)
           r_angular (3, i_ang) = z(i_ang)
           w_angular (i_ang)    = w(i_ang)
        enddo
     end if

  else

     !       n_leb is the angular momentum up to which the chosen
     !       grid integrated a given integrand exactly.

     !       Notice that Delley does not provide all the angular grids
     !       in the original Lebedev-Laikov subroutine ... we have to
     !       patch :-(

     if (i_lebedev.eq.0) then
        n_leb = 0
     else if (i_lebedev.le.1) then
        n_leb = 3
     else if (i_lebedev.le.2) then
        n_leb = 5
     else if (i_lebedev.le.3) then
        n_leb = 7
     else if (i_lebedev.le.5) then
        i_lebedev = 5
        n_leb = 11
     else if (i_lebedev.le.8) then
        i_lebedev = 8
        n_leb = 17
     else if (i_lebedev.le.11) then
        i_lebedev = 11
        n_leb = 23
     else if (i_lebedev.le.14) then
        i_lebedev = 14
        n_leb = 29
     else if (i_lebedev.le.16) then
        i_lebedev = 16
        n_leb = 35
     else if (i_lebedev.le.17) then
        n_leb = 41
     else if (i_lebedev.eq.18) then
        n_leb = 47
     else if (i_lebedev.eq.19) then
        n_leb = 53
     else if (i_lebedev.eq.20) then
        n_leb = 59
     else
        !         we should never reach this point
        write(use_unit,'(1X,A,A)') &
             "Requested Delley-grid which does not exist. This should", &
             " not be possible - sorry."
        stop
     end if
     if (n_leb.gt.0) then
        call anginc &
             (r_angular, w_angular, n_angular, n_max_angular, n_leb)
        if (n_angular.lt.0) then
           !           dimension exceeded, this will be handled outside get_angular_grid
           n_angular = -n_angular
        end if
     end if

  end if

end subroutine get_angular_grid
!******	
