!****h* FHI-aims/hartree_potential_real_p0
!  NAME
!   hartree_potential_real_p0
!  SYNOPSIS

module hartree_potential_real_p0

!  PURPOSE
!    The module contains routines and variables for the far distance part of Hartree potential in cluster systems,
!    and the analytic real space part of Hartree potential in periodic systems.
!
!  USES

  use dimensions
  use runtime_choices
  use geometry
  use species_data
  use analytic_multipole_coefficients
  use pbc_lists
  use constants
  use synchronize_mpi

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



  implicit none
  
  real*8,allocatable,dimension(:,:) :: multipole_c
  real*8,allocatable,dimension(:,:),private:: Fp
  integer,private:: n_hartree_atoms
  real*8,allocatable,dimension(:),private:: hartree_atoms

  integer :: hartree_force_l_add ! set to 0 if no forces and set to 1 if (forces_on) ... 

  ! The following variables are used in:
  !   far_distance_hartree_Fp_periodic_single_atom()               and
  !   far_distance_real_gradient_hartree_potential_single_atom()
  ! They belong to the non-periodic Ewald method.
  integer, parameter,                private  ::  pmaxab = 30
  real*8, dimension(0:pmaxab),       private  ::  b0, b2, b4, b6
  real*8, dimension(0:pmaxab),       private  ::  a_save






contains






!****s*  hartree_potential_real_p0/hartree_potential_real_coeff
!  NAME
!   hartree_potential_real_coeff
!  SYNOPSIS

  subroutine hartree_potential_real_coeff( index_lm, multipole_moments, l_hartree_max_far_distance, n_atom_list)

!  PURPOSE
!  The subroutine allocates and calculates the coefficients for the real space part of the analytic Hartree potential
!
    implicit none
!  ARGUMENTS

    integer :: index_lm(-l_pot_max:l_pot_max, 0:l_pot_max)
    real*8  ::  multipole_moments( ( l_pot_max+1)**2, n_atoms)
    integer :: l_hartree_max_far_distance(n_atoms)
    integer :: n_atom_list

!  INPUTS
!    o index_lm -- order of l and m indexis in multipole_moments table
!    o multipole_moments -- multipole moments
!    o l_hartree_max_far_distance -- maximum l component in far distance Hartree potential
!    o n_atom_list -- number of atoms included in Hartree potential includid periodic mirror images.
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: info
    integer:: n, i_atom, i_l, i_m


    if(.not.allocated(multipole_c))then

       allocate(multipole_c( n_cc_lm_ijk(l_pot_max), n_atoms),stat=info)
       call check_allocation(info, 'multipole_c                   ') 
       allocate(Fp( 0:(l_pot_max+1), n_atom_list),stat=info)
       call check_allocation(info, 'Fp                            ') 
       Fp = 0.0d0
       allocate(hartree_atoms( n_atom_list),stat=info)
       call check_allocation(info, 'hartree_atoms                 ') 

       hartree_atoms = 0

    end if

    do i_atom = 1,n_atoms

       do n = 1,n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)),1

          i_l = index_cc(n, 1)
          i_m = index_cc(n, 2)


          if(i_l <=  l_hartree_max_far_distance(i_atom) .and. abs(multipole_moments(  index_lm(i_m, i_l), i_atom )) > 1e-10  )then
             multipole_c( n,i_atom) =  cc(n) * multipole_moments(  index_lm(i_m, i_l), i_atom)
          else
             multipole_c( n,i_atom) =  0.0
          end if

       end do
    end do

  end subroutine hartree_potential_real_coeff
  !******





!---------------------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_hartree_Fp_cluster
!  NAME
!   far_distance_hartree_Fp_cluster
!  SYNOPSIS

  subroutine far_distance_hartree_Fp_cluster(atom_list, n_atoms_list, dist_tab, &
       l_hartree_max_far_distance, forces_on )


!  PURPOSE
!  The subroutine calculates and saves Fp-functions for non-periodic systems.
!  These are needed for the analytic real space Hartree potential (the far
!  distance part).
!
    implicit none
!  ARGUMENTS

  integer:: n_atoms_list
  integer:: atom_list(n_atoms_list)
  real*8 :: dist_tab(n_atoms)
  integer:: l_hartree_max_far_distance(n_atoms)
  logical :: forces_on

!  INPUTS
!  o n_atoms_list -- number of relevant atoms in this routine
!  o atom_list --  relevant atoms in this routine
!  o dist_tab -- distance to relevant atoms
!  o l_hartree_max_far_distance -- maximum l in far distance Hartree potential
!  o forces_on -- are we calculating forces in this round or not?
!
!  OUTPUT
!    none 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  !local variables
  real*8  :: dist_sq
  integer :: one_minus_2l
  integer :: l_max

  ! counters
  integer:: i_atom_list, i_l

  ! Recursive tabulation of Fp ...
  do i_atom_list = 1,n_atoms_list
     l_max = l_hartree_max_far_distance(atom_list(i_atom_list)) &
          + hartree_force_l_add

     ! needed constant
     dist_sq = dist_tab(i_atom_list)**2

     ! initialize everything
     one_minus_2l = 1
     Fp(0,i_atom_list) = 1.d0/dist_tab(i_atom_list)

     do i_l = 1, l_max
        one_minus_2l = one_minus_2l-2
        Fp(i_l,i_atom_list) = Fp(i_l-1,i_atom_list) * dble(one_minus_2l)/dist_sq
     end do

  end do

!!$  if (forces_on) then
!!$    ! one more Fp needed per atom
!!$    do i_atom_list = 1,n_atoms_list
!!$     l_max = l_hartree_max_far_distance(atom_list(i_atom_list))
!!$     ! needed constant
!!$     dist_sq = dist_tab(i_atom_list)**2
!!$
!!$     Fp(l_max+1,i_atom_list) = Fp(l_max,i_atom_list) * dble(-1-2*l_max)/dist_sq     
!!$
!!$    enddo
!!$  end if

  n_hartree_atoms = n_atoms_list
  hartree_atoms(1:n_atoms_list) = atom_list

end subroutine far_distance_hartree_Fp_cluster
!******





!---------------------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_hartree_Fp_cluster_single_atom
!  NAME
!   far_distance_hartree_Fp_cluster_single_atom
!  SYNOPSIS

subroutine far_distance_hartree_Fp_cluster_single_atom(current_atom, &
     dist_tab, l_hartree_max_far_distance, forces_on )

!  PURPOSE
!  The subroutine calculates and saves Fp-functions for non-periodic systems.
!  These are needed for the analytic real space Hartree potential (the far
!  distance part).
!
!  This is otherwise same than far_distance_hartree_Fp_cluster, but only single
!  atom contribution is calculated.
!
    implicit none
!  ARGUMENTS

  integer:: current_atom
  real*8 :: dist_tab
  integer:: l_hartree_max_far_distance(n_atoms)
  logical :: forces_on

!  INPUTS
!  o current_atom -- atom which multipole component is now calculated.
!  o dist_tab -- distance to the atom
!  o l_hartree_max_far_distance -- maximum l in far distance Hartree potential
!  o forces_on -- are we calculating forces in this round or not?

!  OUTPUT
!    none 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  !local variables
  real*8  :: dist_sq
  integer :: one_minus_2l
  integer :: l_max

   ! counters
  integer:: i_l

  ! Recursive tabulation of Fp ...
  l_max = l_hartree_max_far_distance(current_atom)+hartree_force_l_add

  ! needed constant
  dist_sq = dist_tab**2

  ! initialize everything
  one_minus_2l = 1
  Fp(0,current_atom) = 1.d0/dist_tab

  do i_l = 1, l_max
     one_minus_2l = one_minus_2l-2
     Fp(i_l,current_atom) = Fp(i_l-1,current_atom) * dble(one_minus_2l)/dist_sq
  end do

!!$  if (forces_on) then
!!$     ! one more Fp needed per atom
!!$     l_max = l_hartree_max_far_distance(current_atom)
!!$     ! needed constant
!!$     dist_sq = dist_tab**2
!!$
!!$     Fp(l_max+1,current_atom) = Fp(l_max,current_atom) * dble(-1-2*l_max)/dist_sq     
!!$
!!$  end if

  n_hartree_atoms = 1
  hartree_atoms(1) = current_atom

end subroutine far_distance_hartree_Fp_cluster_single_atom
!******





!---------------------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_hartree_Fp_cluster_single_atom_p2
!  NAME
!   far_distance_hartree_Fp_cluster_single_atom_p2
!  SYNOPSIS

subroutine far_distance_hartree_Fp_cluster_single_atom_p2 &
     (  dist_tab, l_max, forces_on )


!  PURPOSE
!  The subroutine calculates and saves Fp-functions for non-periodic systems.
!  These are needed for the analytic real space Hartree potential (the far
!  distance part).
!
!  This is updated version of otherwise same than far_far_distance_hartree_Fp_cluster_single_atom
!
    implicit none
!  ARGUMENTS

  real*8 :: dist_tab
  logical :: forces_on
  integer :: l_max

!  INPUTS
!  o dist_tab -- distance to the atom
!  o l_max -- maximum l in far distance Hartree potential
!  o forces_on -- are we calculating forces in this round or not?

!  OUTPUT
!    none 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  !local variables
  real*8  :: dist_sq
  integer :: one_minus_2l

  ! counters
  integer:: i_l

  ! needed constant
  dist_sq = dist_tab**2

  ! initialize everything
  one_minus_2l = 1
  Fp(0,1) = 1.d0/dist_tab

  do i_l = 1, l_max+hartree_force_l_add
     one_minus_2l = one_minus_2l-2
     Fp(i_l,1) = Fp(i_l-1,1) * dble(one_minus_2l)/dist_sq
  end do

!!$  if (forces_on) then
!!$
!!$     Fp(l_max+1,1) = Fp(l_max,1) * dble(-1-2*l_max)/dist_sq     
!!$
!!$  end if

end subroutine far_distance_hartree_Fp_cluster_single_atom_p2
!******





!---------------------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_hartree_Fp_periodic
!  NAME
!   far_distance_hartree_Fp_periodic
!  SYNOPSIS

subroutine far_distance_hartree_Fp_periodic(atom_list_in, n_atoms_list_in, &
     atom_list_out, n_atoms_list_out, & 
     dist_tab_in, dist_tab_out,  l_hartree_max_far_distance, forces_on )

!  PURPOSE
!
!  The subroutine calculates and saves Fp-functions for periodic systems.
!  These are needed for the analytic real space Hartree potential (the far
!  distance part).
!
    use Hartree_F_p_functions, only: F_erf, F_erfc
    implicit none
!  ARGUMENTS

    integer:: n_atoms_list_in,  n_atoms_list_out
    integer:: atom_list_in(n_atoms_list_in)
    integer:: atom_list_out(n_atoms_list_out)
    real*8 :: dist_tab_in(n_atoms_list_in)
    real*8 :: dist_tab_out(n_atoms_list_out)
    integer:: l_hartree_max_far_distance(n_atoms)
    logical :: forces_on

!  INPUTS
!    o n_atoms_list_in -- number of atoms for potential inside multipole radius
!    o n_atoms_list_out -- number of atoms for potential  outside multipole radius
!    o atom_list_in -- atoms for potential inside multipole radius
!    o atom_list_out --  atoms for potential  outside multipole radius
!    o dist_tab_in -- distances to atoms for potential inside multipole radius
!    o dist_tab_out -- distances to atoms for potential outside multipole radius
!    o l_hartree_max_far_distance -- maximum l in for distance Hartree potential
!    o forces_on -- is the forces calculated
!   
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  integer:: i_atom_list,  i_atom, i_l 

  do i_atom_list = 1,n_atoms_list_in
     i_atom = center_to_atom(atom_list_in(i_atom_list))

!!$     do i_l = 0, l_hartree_max_far_distance(i_atom)+hartree_force_l_add
!!$        Fp(i_l,i_atom_list) = -F_erf(  dist_tab_in(i_atom_list), i_l )
!!$     end do
     i_l = l_hartree_max_far_distance(i_atom)+hartree_force_l_add
     call F_erf(Fp(0:i_l,i_atom_list),dist_tab_in(i_atom_list),i_l)
     Fp(0:i_l,i_atom_list) = -Fp(0:i_l,i_atom_list)
  end do

  do i_atom_list = 1,n_atoms_list_out
     i_atom = center_to_atom(atom_list_out(i_atom_list))

!!$     do i_l = 0, l_hartree_max_far_distance(i_atom)+hartree_force_l_add
!!$        Fp(i_l,i_atom_list + n_atoms_list_in) = F_erfc(  dist_tab_out(i_atom_list), i_l )
!!$     end do
     i_l = l_hartree_max_far_distance(i_atom)+hartree_force_l_add
     call F_erfc(Fp(0:i_l,i_atom_list+n_atoms_list_in),dist_tab_out(i_atom_list),i_l)
  end do

  n_hartree_atoms = n_atoms_list_in +  n_atoms_list_out

  if(n_atoms_list_in >0)then
     hartree_atoms(1:n_atoms_list_in) = atom_list_in
  end if

  if(n_atoms_list_out > 0)then
     hartree_atoms(n_atoms_list_in+1:n_atoms_list_in+n_atoms_list_out) = &
          atom_list_out(1:n_atoms_list_out)
  end if
     

!!$  if(forces_on)then
!!$     do i_atom_list = 1,n_atoms_list_in
!!$        i_atom = center_to_atom(atom_list_in(i_atom_list))
!!$
!!$        i_l = l_hartree_max_far_distance(i_atom)+1
!!$        Fp(i_l,i_atom_list) = -F_erf(  dist_tab_in(i_atom_list), i_l )
!!$        
!!$     end do
!!$
!!$     do i_atom_list = 1,n_atoms_list_out
!!$        i_atom = center_to_atom(atom_list_out(i_atom_list))
!!$        
!!$        i_l = l_hartree_max_far_distance(i_atom)+1
!!$        Fp(i_l,i_atom_list + n_atoms_list_in) = F_erfc(  dist_tab_out(i_atom_list), i_l )
!!$     end do
!!$  end if
     
end subroutine far_distance_hartree_Fp_periodic
!******





!--------------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_hartree_Fp_periodic_single_atom
!  NAME
!   far_distance_hartree_Fp_periodic_single_atom
!  SYNOPSIS

subroutine far_distance_hartree_Fp_periodic_single_atom(current_atom, i_center, &
     dist, l_hartree_max_far_distance, inside, forces_on,                       &
     multipole_radius_sq, adap_outer_radius, non_peri_extd )

!  PURPOSE
!
!  The subroutine calculates and saves Fp-functions for periodic systems.
!  These are needed for the analytic real space Hartree potential (the far
!  distance part).
! 
!  This is otherwise same routine than far_distance_hartree_Fp_periodic, but it calculates only
!  single atom contribution.
!
    use Hartree_F_p_functions, only: F_erf, F_erfc
    implicit none
!  ARGUMENTS

    integer:: current_atom, i_center
    real*8 :: dist
    integer:: l_hartree_max_far_distance(n_atoms)
    logical :: inside
    logical :: forces_on
    real*8,  intent(in)           :: multipole_radius_sq, adap_outer_radius
    logical, intent(in), optional :: non_peri_extd

!  INPUTS
!    o current_atom -- atom which contribution is calculated
!    o i_center -- the atom center  which contribution is calculated
!    o dist -- distance to atom center
!    o l_hartree_max_far_distance -- maximum l in for distance Hartree potential
!    o inside -- is the potential of the atom center inside or outside multipole radius.
!    o forces_on -- is the forces calculated
!    o multipolse_radius_sq -- see sum_up_whole_potential_p1.f90
!    o adap_outer_radius -- see sum_up_whole_potential_p1.f90
!    o non_peri_extd -- used for the non-periodic Ewald method
!   
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  ! General variables:
  integer :: lmax

  ! Variables used in section "B)" for the non-periodic Ewald method:
  integer        :: n, p
  logical, save  :: firstcall = .true.
  real*8,  save  :: a(0:pmaxab)
  real*8         :: doublefact, drel


  lmax  =  l_hartree_max_far_distance(current_atom) + hartree_force_l_add



  if ( .not. use_hartree_non_periodic_ewald ) then


    ! A) Standard case

    if (inside) then

!!$     if(forces_on)then
!!$
!!$        do i_l = 0, l_hartree_max_far_distance(current_atom)+ 1
!!$           Fp(i_l,i_center) = -F_erf(  dist, i_l )
!!$        end do
!!$        
!!$     else

!!$        do i_l = 0, l_hartree_max_far_distance(current_atom)+hartree_force_l_add
!!$           Fp(i_l,i_center) = -F_erf(  dist, i_l )
!!$        end do

          call F_erf(Fp(0:lmax,i_center),dist,lmax)
          Fp(0:lmax,i_center) = -Fp(0:lmax,i_center)

!     end if

   else

!     if(forces_on)then 

!!$        do i_l = 0, l_hartree_max_far_distance(current_atom)+hartree_force_l_add
!!$           Fp(i_l,i_center) = F_erfc(  dist, i_l )
!!$        end do
        
          call F_erfc(Fp(0:lmax,i_center),dist,lmax)
        
!!$     else
!!$
!!$        do i_l = 0, l_hartree_max_far_distance(current_atom)
!!$           Fp(i_l,i_center) = F_erfc(  dist, i_l )
!!$        end do
!!$
!!$     end if

    end if  ! inside




  else  ! use_hartree_non_periodic_ewald = .true.


    ! B) Using polynomials rather than the error function for the non-periodic
    !    Ewald method.


    ! 1. If the routine is called for the first time, the coefficient vectors
    !    are initialized.

    if (firstcall) then

      if ( l_hartree_far_distance + 1 > pmaxab ) then   ! "+1" for forces
        call aims_stop('pmaxab insufficient')
      end if

      do p = 0, pmaxab

        a(p)  =  (-1)**p  *  product((/  ( real(n,8),  n = 1, 2*p-1, 2 )  /))
        a_save(p) = a(p)

        doublefact  =  product((/  ( real(n,8),  n = 1, 2*p+7, 2 )  /))
        b0(p)  =    (-1)**p * doublefact / (2*p+1) / 48
        b2(p)  =  - (-1)**p * doublefact / (2*p+3) / 16
        b4(p)  =    (-1)**p * doublefact / (2*p+5) / 16
        b6(p)  =  - (-1)**p * doublefact / (2*p+7) / 48

      end do

      firstcall = .false.

    end if


    ! 2. Delete the Fp-values for the cases where section 3 is skipped but
    !    section 4 is evaluated.
    Fp( 0:lmax, i_center )  =  0


    ! 3. The interval [0, adap_outer_radius[ is considered where the bridge
    !    potentials have to be applied. According to the convention used in
    !    sum_up_whole_potential, the upper endpoint of the interval is excluded.

    if ( dist < adap_outer_radius ) then


      ! 3a) Calculate the radial functions for the bridge potentials.

      drel = dist / adap_outer_radius

      do p = 0, lmax
        Fp( p, i_center )  =  1 / adap_outer_radius**(2*p+1) *                          &
                              (  b0(p) + b2(p)*drel**2 + b4(p)*drel**4 + b6(p)*drel**6  )
      end do


      ! 3b) These values are correct for the extended part, but for the
      !     localized part the opposite sign is needed.
      if ( .not. present(non_peri_extd) ) then
        Fp( 0:lmax, i_center )  =  - Fp( 0:lmax, i_center )
      end if


    end if


    ! 4. Next, the outer interval [multipole_radius, infinity[, that is, the
    !    region outside of the density cloud, is considered. Note that this
    !    region generally intersects with the previous one. In the outer region,
    !    the radial part of the potentials is a(l)/dist^(2l+1). These functions
    !    have to be added in all cases except for the following: we are in the
    !    interspace [multipole_radius_sq, adap_outer_radius[ and the extended
    !    part is considered. Inclusion or exclusion of the endpoints of the
    !    respective intervals follows the convention used in
    !    sum_up_whole_potential.

    if ( dist**2 >= multipole_radius_sq  .and.       &
         .not.  (  dist < adap_outer_radius          &
                   .and. present(non_peri_extd)  )   &
       ) then

      do p = 0, lmax
        Fp( p, i_center )  =  Fp( p, i_center )  +  a(p) / dist**(2*p+1)
      end do

    end if


  end if  ! .not. use_hartree_non_periodic_ewald



!!$  n_hartree_atoms = n_atoms_list_in +  n_atoms_list_out
!!$
!!$  if(n_atoms_list_in >0)then
!!$     hartree_atoms(1:n_atoms_list_in) = atom_list_in
!!$  end if
!!$
!!$  if(n_atoms_list_out > 0)then
!!$     hartree_atoms(n_atoms_list_in+1:n_atoms_list_in+n_atoms_list_out) = &
!!$          atom_list_out(1:n_atoms_list_out)
!!$  end if

  ! n_hartree_atoms = 1
  ! hartree_atoms(1) = i_center


end subroutine far_distance_hartree_Fp_periodic_single_atom
!******


subroutine far_distance_hartree_Fp_periodic_single_atom_init()

   ! Variables used in section "B)" for the non-periodic Ewald method:
   integer        :: n, p
   logical, save  :: firstcall = .true.
   real*8,  save  :: a(0:pmaxab)
   real*8         :: doublefact, drel
 
     if (firstcall) then
 
       if ( l_hartree_far_distance + 1 > pmaxab ) then   ! "+1" for forces
         call aims_stop('pmaxab insufficient')
       end if
 
       do p = 0, pmaxab
 
         a(p)  =  (-1)**p  *  product((/  ( real(n,8),  n = 1, 2*p-1, 2 )  /))
         a_save(p) = a(p)
 
         doublefact  =  product((/  ( real(n,8),  n = 1, 2*p+7, 2 )  /))
         b0(p)  =    (-1)**p * doublefact / (2*p+1) / 48
         b2(p)  =  - (-1)**p * doublefact / (2*p+3) / 16
         b4(p)  =    (-1)**p * doublefact / (2*p+5) / 16
         b6(p)  =  - (-1)**p * doublefact / (2*p+7) / 48
 
       end do
 
       firstcall = .false.
 
     end if
 
end subroutine far_distance_hartree_Fp_periodic_single_atom_init



!---------------------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/update_outer_radius_l
!  NAME
!    update_outer_radius_l
!  SYNOPSIS

subroutine update_outer_radius_l( outer_radius_sq, multipole_moments, multipole_radius_sq, l_max, index_lm)

!  PURPOSE
!  The subroutine calculates the outer radius of  the real space part of the Hartree potential.
!  This is only for periodic systems.
!
  use Hartree_F_p_functions, only: F_erfc_single
  implicit none
!  ARGUMENTS

  real*8:: outer_radius_sq(0:l_pot_max, n_atoms)
  real*8, dimension( ( l_pot_max+1)**2, n_atoms) :: multipole_moments
  real*8, dimension(n_atoms) :: multipole_radius_sq
  integer :: l_max(n_atoms)
  integer :: index_lm(-l_pot_max:l_pot_max, 0:l_pot_max)

!  INPUTS
!   o multipole_moments -- multipole moments
!   o multipole_radius_sq -- (radius of multipole components)**2
!   o l_max -- maximum l component
!   o index_lm -- order of l and m components in multipole moments
!
!  OUTPUT
!   o outer_radius_sq -- outer radius of the real space part of the Hartree potential.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  integer:: l_hartree_max
  real*8::  pot, radius_max, radius_min, radius, ra, rb
  integer:: i_atom, i_l, n, i_div
  integer:: ll, mm, ii, jj, kk
  ! Number of interval division steps for [radius_min ... radius_max]
  ! for determining outer radius
  integer, parameter :: N_DIVISION_STEPS = 9


  if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then

     outer_radius_sq(:,:) = 0.

     do i_atom = 1,n_atoms

        if(mod(i_atom-1,n_tasks) /= myid) cycle ! distribute work over tasks

        radius_max = sqrt(atom_radius_hartree_sq(species(i_atom)))
        radius_min = sqrt(multipole_radius_sq(i_atom))
        outer_radius_sq(:,i_atom) = radius_min
        
        do n = 1,n_cc_lm_ijk(l_max(i_atom)),1
           
           ll =  index_cc(n, 1)
           mm =  index_cc(n, 2)
           ii =  index_cc(n, 3)
           jj =  index_cc(n, 4)
           kk =  index_cc(n, 5)

           ! Check if abs(pot) at current outer_radius_sq(ll, i_atom) is
           ! already smaller than threshold, in this case there needs nothing to be done.
           ! This also catches the case if abs(pot(radius_min)) < threshold
           ! since the initial value of outer_radius_sq(ll, i_atom) is radius_min

           radius = outer_radius_sq(ll, i_atom)
           pot    = multipole_c( n,i_atom) &
                  * F_erfc_single( radius, ll ) * radius**(ii+jj+kk)
           if(abs(pot) < far_distance_adaptive_hartree_radius_threshold) cycle 

           ! do a binary search for radius with
           ! abs(pot) < far_distance_adaptive_hartree_radius_threshold

           ra = radius_min ! we could also use outer_radius_sq(ll, i_atom) here
           rb = radius_max
           do i_div = 1, N_DIVISION_STEPS
              radius = 0.5*(ra+rb)
              pot    = multipole_c( n,i_atom) &
                     * F_erfc_single( radius, ll ) * radius**(ii+jj+kk)
              if (abs(pot) < far_distance_adaptive_hartree_radius_threshold) then
                 rb = radius
              else
                 ra = radius
              endif
           enddo

           outer_radius_sq(ll, i_atom) =   max( outer_radius_sq(ll, i_atom),  rb)

        end do
        
       
        do ll = 0,l_max(i_atom)
           outer_radius_sq(ll, i_atom) = outer_radius_sq(ll, i_atom)**2
        end do

     end do

     call sync_matrix(outer_radius_sq, l_pot_max+1, n_atoms)

  end if

end subroutine update_outer_radius_l
!******





!---------------------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_real_hartree_potential
!  NAME
!    far_distance_real_hartree_potential
!  SYNOPSIS

subroutine far_distance_real_hartree_potential( potential,  l_hartree_max_far_distance,  coord_current )

!  PURPOSE
!  The subroutine calculates the real space part of the far distance periodic Hartree potential.
!
    implicit none
!  ARGUMENTS

  real*8:: potential
  integer:: l_hartree_max_far_distance(n_atoms) 
  real*8 :: coord_current(3)

!  INPUTS
!  o l_hartree_max_far_distance - maximum l component
!  o coord_current -- coordinates
!
!  OUTPUT
!  o potential -- Hartree potential is added here
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  real*8:: dpot, c_pot
  integer i_atom, i_atom_list
  integer i_l,i_l2

  integer i_coord

  integer:: ii, jj, kk, n, l_max

  real*8:: coord_c(0:l_pot_max,3)
  real*8:: coord_mat(0:l_pot_max,0:l_pot_max)
  real*8::  rest_mat(0:l_pot_max,0:l_pot_max)
  real*8:: vector(n_cc_lm_ijk(l_pot_max))


  real*8 :: dir(3)


  real*8, external :: ddot


  ! Begin work.

  c_pot = 0.

  do i_atom_list = 1,n_hartree_atoms
  
     i_atom =  hartree_atoms(i_atom_list)
     l_max =  l_hartree_max_far_distance(center_to_atom(i_atom))


     dpot = 0.

     coord_c(0,1:3) = 1.0
     dir = coord_current(1:3)-coords_center(1:3,i_atom)
     

     do i_coord = 1,3,1
       do i_l = 1, index_ijk_max_cc(i_coord, l_max ) !l_pot_max

         coord_c(i_l,i_coord) =  dir(i_coord) *  coord_c(i_l-1,i_coord)

       end do
     end do
   
     ! x**n1 * y**n2
     do i_l = 0, index_ijk_max_cc(1,l_max )
        do i_l2 = 0,  index_ijk_max_cc(2,l_max )

           coord_mat(i_l,i_l2) =  coord_c(i_l,1)*coord_c(i_l2,2)

        end do
     end do

     ! z**n * Fp
     do i_l = 0,  index_ijk_max_cc(3,l_max )
        do i_l2 = 0, l_max

              rest_mat(i_l,i_l2) =  coord_c(i_l,3)* Fp(i_l2,i_atom_list)

        end do
     end do

     ! x**n1 * y**n2 * z**n3 * Fp
     do n = 1,n_cc_lm_ijk(l_max),1

           ii =  index_cc(n, 3)
           jj =  index_cc(n, 4)
           kk =  index_cc(n, 5)

           vector(n) = coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))


     end do

     ! Now sum up all multipole components from present atom
     !dpot = ddot(n_cc_lm_ijk(l_max), vector, 1, multipole_c(1,center_to_atom(i_atom)), 1)
     dpot = dot_product(vector(1:n_cc_lm_ijk(l_max)), multipole_c(1:n_cc_lm_ijk(l_max),center_to_atom(i_atom)) )

!     write(use_unit,*) vector


     if(abs(dpot) .gt. 1e-30)then
          c_pot = c_pot + dpot
     end if
     dpot = 0.0
        
  end do

  potential = potential + c_pot

end subroutine far_distance_real_hartree_potential
!******





!-------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_real_hartree_potential_single_atom
!  NAME
!   far_distance_real_hartree_potential_single_atom
!  SYNOPSIS

subroutine far_distance_real_hartree_potential_single_atom &
     ( current_center, i_center,  potential, l_hartree_max_far_distance,coord_current )

!  PURPOSE
!  The subroutine calculates the real space part of the far distance periodic Hartree potential.
!  This routine calculates single atom contribution.
!
    implicit none
!  ARGUMENTS

  integer :: current_center
  integer :: i_center
  integer:: l_hartree_max_far_distance(n_atoms) 
  real*8:: potential

!  INPUTS
!  o current_center -- the atom center which contribution is calculated, absolute index
!  o i_center --  the atom center which contribution is calculated, relative index
!  o l_hartree_max_far_distance -- maximum l-component
!
!  OUTPUT
!  o potential -- potential is added here.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  ! local variables
  real*8:: dpot, c_pot

  integer i_l,i_l2

  integer i_coord
  integer:: ii, jj, kk, n, l_max

  real*8:: coord_c(0:l_pot_max,3)
  real*8:: coord_mat(0:l_pot_max,0:l_pot_max)
  real*8::  rest_mat(0:l_pot_max,0:l_pot_max)
  real*8:: vector(n_cc_lm_ijk(l_pot_max))

  real*8 :: coord_current(3), dir(3)


  real*8, external :: ddot

  ! Begin work.

  c_pot = 0.0d0

  l_max =  l_hartree_max_far_distance(center_to_atom(current_center))

  dpot = 0.0d0

  coord_c(0,1:3) = 1.0d0
  dir = coord_current(1:3)-coords_center(1:3,current_center)
     

  do i_coord = 1,3,1
     do i_l = 1, index_ijk_max_cc(i_coord, l_max ) !l_pot_max

        coord_c(i_l,i_coord) =  dir(i_coord) *  coord_c(i_l-1,i_coord)

     end do
  end do
   
     ! x**n1 * y**n2
  do i_l = 0, index_ijk_max_cc(1,l_max )
     do i_l2 = 0,  index_ijk_max_cc(2,l_max )

        coord_mat(i_l,i_l2) =  coord_c(i_l,1)*coord_c(i_l2,2)
        
     end do
  end do

     ! z**n * Fp
  do i_l = 0,  index_ijk_max_cc(3,l_max )
     do i_l2 = 0, l_max

        rest_mat(i_l,i_l2) =  coord_c(i_l,3)* Fp(i_l2,i_center)
        
     end do
  end do

     ! x**n1 * y**n2 * z**n3 * Fp
  do n = 1,n_cc_lm_ijk(l_max),1
     
     ii =  index_cc(n, 3)
     jj =  index_cc(n, 4)
     kk =  index_cc(n, 5)
     
     vector(n) = coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))

  end do

  ! Now sum up all multipole components from present atom
  !dpot = ddot(n_cc_lm_ijk(l_max), vector, 1, multipole_c(1,center_to_atom(current_center)), 1)
  dpot = dot_product(vector(1:n_cc_lm_ijk(l_max)) , multipole_c(1:n_cc_lm_ijk(l_max),center_to_atom(current_center)) ) 
  
  if(abs(dpot) .gt. 1e-30)then
     c_pot = c_pot + dpot
  end if
  dpot = 0.0d0
  
  potential = potential + c_pot

end subroutine far_distance_real_hartree_potential_single_atom
!******





!-------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_real_hartree_potential_single_atom_p2
!  NAME
!   far_distance_real_hartree_potential_single_atom_p2
!  SYNOPSIS

subroutine far_distance_real_hartree_potential_single_atom_p2 &
     ( i_center, potential, l_max, coord_current )

!  PURPOSE
!  The subroutine calculates the real space part of the far distance periodic Hartree potential.
!  This routine calculates single atom contribution.
!  This is updated version of far_distance_real_hartree_potential_single_atom
!
    implicit none
!  ARGUMENTS

  integer :: i_center
  real*8  :: potential
  integer :: l_max 
  real*8  :: coord_current(3)

!  INPUTS
!  o i_center --  the atom center which contribution is calculated, relative index
!  o l_max -- maximum l-component
!  o coord_current -- coordinates
!
!  OUTPUT
!  o potential -- potential is added here.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  ! local variables
  real*8:: dpot
  integer:: i_l, ii, jj, kk, n, nn
  real*8:: coord_c(0:l_pot_max,3)
  real*8 :: dir(3)
  real*8, external :: ddot

  ! Begin work.

  dpot = 0.0d0

  coord_c(0,1:3) = 1.0d0
  dir = coord_current(1:3)-coords_center(1:3,i_center)
     

  
  do i_l = 1, MAXVAL(index_ijk_max_cc(:, l_max ))

     coord_c(i_l,1) =  dir(1) *  coord_c(i_l-1,1)
     coord_c(i_l,2) =  dir(2) *  coord_c(i_l-1,2)
     coord_c(i_l,3) =  dir(3) *  coord_c(i_l-1,3)

  end do
   
     ! x**n1 * y**n2 * z**n3 * Fp
  dpot = 0
  do n = 1,n_cc_lm_ijk(l_max),1 !sum over all i,j,k tuple, for which l=i+j+k
     
     ii =  index_cc(n, 3) !lambda in Delley
     jj =  index_cc(n, 4) !mu in Delley
     kk =  index_cc(n, 5) !nu in Delley
     nn =  index_cc(n, 6) !p=(i+j+k+lambda+mu+nu)/2 in Delley
     
     dpot = dpot + coord_c(ii,1)*coord_c(jj,2)*coord_c(kk,3)*Fp(nn,1)*multipole_c(n,center_to_atom(i_center))

  end do

  potential = potential + dpot

end subroutine far_distance_real_hartree_potential_single_atom_p2
!******





!-------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_real_gradient_hartree_potential
!  NAME
!   far_distance_real_gradient_hartree_potential
!  SYNOPSIS

subroutine far_distance_real_gradient_hartree_potential( dir_tab, gradient,l_hartree_max_far_distance )

!  PURPOSE
!  The subroutine calculates gradient of the real space part of the far distance Hartree potential
!  in periodic systems.
!
    implicit none
!  ARGUMENTS

  real*8, dimension(3, n_atoms) :: dir_tab
  real*8:: gradient(3,n_atoms)
  integer:: l_hartree_max_far_distance(n_atoms) 

!  INPUTS
!  o dir_tab -- directios to atoms
!  o l_hartree_max_far_distance -- maximum l component
!
!  OUTPUT
!  o gradient -- gradient is added here
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  integer i_atom, i_atom_list
  integer i_l,i_l2
  integer i_coord
  integer:: ii, jj, kk, n

  real*8:: coord_c(0:l_pot_max+1,3)
  real*8:: coord_mat(0:l_pot_max+1,0:l_pot_max+1)
  real*8:: rest_mat(0:l_pot_max+1,0:l_pot_max+1)
  real*8:: vector(n_cc_lm_ijk(l_pot_max+1))

  real*8, external :: ddot

! begin work

!  write(use_unit,*) 

  do i_atom_list = 1,n_hartree_atoms
  
     i_atom =  hartree_atoms(i_atom_list)

!     write(use_unit,*) "i_atom = ", i_atom
!     write(use_unit,*) "Fp = ", Fp

        ! x**n, y**n, and z**n

        coord_c(0,1:3) = 1.0
        do i_coord = 1,3,1
           do i_l = 1, index_ijk_max_cc(i_coord, l_hartree_max_far_distance(i_atom)+1) !l_pot_max

           coord_c(i_l,i_coord) =  dir_tab(i_coord,i_atom) *  coord_c(i_l-1,i_coord)

          end do
        end do

        
        ! x**n1 * y**n2

        do i_l = 0, index_ijk_max_cc(1,l_hartree_max_far_distance(i_atom)+1 )
         do i_l2 = 0,  index_ijk_max_cc(2,l_hartree_max_far_distance(i_atom)+1)

           coord_mat(i_l,i_l2) =  coord_c(i_l,1)*coord_c(i_l2,2)

         end do
        end do


        ! z**n * Fp

        do i_l = 0,  index_ijk_max_cc(3,l_hartree_max_far_distance(i_atom)+1)
           do i_l2 = 0, l_hartree_max_far_distance(i_atom)+1

              rest_mat(i_l,i_l2) =  coord_c(i_l,3)* Fp(i_l2,i_atom_list)

           end do
        end do

! X direction:

        ! x**n1 * y**n2 * z**n3 * Fp
        do n = 1,n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)),1

              ii =  index_cc(n, 3)-1
              jj =  index_cc(n, 4)
              kk =  index_cc(n, 5)

              if(ii >= 0)then
                 vector(n) = (ii+1)* coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))
              else
                 vector(n) = 0.d0
              end if

              ii =  index_cc(n, 3)+1

              vector(n) = vector(n) + coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6)+1)
              
           end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(1,i_atom) =   gradient(1,i_atom) + dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom))), &
             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)), i_atom  ))

!        write(use_unit,*) "vector: "
!        write(use_unit,*) n_cc_lm_ijk(l_hartree_max_far_distance(i_atom))
!        write(use_unit,*) vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)))
!        write(use_unit,*) "Gradient(1): "
!        write(use_unit,*)  gradient(1,i_atom)
!        write(use_unit,*) - dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom))), &
!             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)), i_atom  ))
!        write(use_unit,*)

! Y direction:

        do n = 1,n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)),1


              ii =  index_cc(n, 3)
              jj =  index_cc(n, 4)-1
              kk =  index_cc(n, 5)

              if(jj >= 0)then
                 vector(n) = (jj+1)* coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))
              else
                 vector(n) = 0.d0
              end if

              jj =  index_cc(n, 4)+1

              vector(n) = vector(n) + coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6)+1)

           end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(2,i_atom) = gradient(2,i_atom)  + dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom))), &
             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)), i_atom  ))

! Z direction:

        do n = 1,n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)),1


              ii =  index_cc(n, 3)
              jj =  index_cc(n, 4)
              kk =  index_cc(n, 5)-1

              if(kk >= 0)then
                 vector(n) = (kk+1)* coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))
              else
                 vector(n) = 0.d0
              end if

              kk =  index_cc(n, 5)+1

              vector(n) = vector(n) + coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6)+1)

           end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(3,i_atom) =  gradient(3,i_atom) + dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom))), &
             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)), i_atom  ))

     end do

   end subroutine far_distance_real_gradient_hartree_potential
!******





!-------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_real_gradient_hartree_potential_single_atom
!  NAME
!    far_distance_real_gradient_hartree_potential_single_atom
!  SYNOPSIS

   subroutine far_distance_real_gradient_hartree_potential_single_atom( current_atom, &
        i_center, dir_tab, gradient, l_hartree_max_far_distance, distsq, adap_outer_radius_sq )


!  PURPOSE
!  The subroutine calculates gradient of the real space part of the far distance Hartree potential
!  in periodic systems. This routine calculates only single atom contribution.
!
    implicit none
!  ARGUMENTS

    integer :: current_atom, i_center
    real*8, dimension(3) :: dir_tab
    real*8:: gradient(3)
    integer:: l_hartree_max_far_distance(n_atoms) 
    real*8, intent(in) :: distsq, adap_outer_radius_sq

!  INPUTS
!   o current_atom -- atom which contribution is calculated, absolute index
!   o i_center -- center which contribution is calculated, relative index
!   o dir_tab -- direction to center
!   o l_hartree_max_far_distance -- maximum l component
!   o distsq -- squared distance to atom center
!   o adap_outer_radius_sq -- see sum_up_whole_potential_p1.f90
!
!  OUTPUT
!   o  gradient -- gradient is added here.
!  
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  ! local variables

     integer i_coord
     integer :: n, ii, jj, kk, p

!  real*8:: F(0:l_pot_max)
     real*8:: coord_c(0:l_pot_max+1,3)
     real*8:: coord_mat(0:l_pot_max+1,0:l_pot_max+1)
     real*8:: rest_mat(0:l_pot_max+1,0:l_pot_max+1)
     real*8:: vector(n_cc_lm_ijk(l_pot_max+1))

     real*8 :: drelsq
     real*8 :: Fp_deriv_red(0:pmaxab)

     real*8, external :: ddot

        ! begin work

        ! x**n, y**n, and z**n

        coord_c(0,1:3) = 1.0
        do i_coord = 1,3
           do n = 1, index_ijk_max_cc &
                       ( i_coord, l_hartree_max_far_distance(current_atom)+1 )  !l_pot_max

              coord_c(n,i_coord)  =  dir_tab(i_coord) * coord_c(n-1,i_coord)

           end do
        end do


        ! x**ii * y**jj

        do ii = 0, index_ijk_max_cc(1,l_hartree_max_far_distance(current_atom)+1 )
           do jj = 0, index_ijk_max_cc(2,l_hartree_max_far_distance(current_atom)+1)

              coord_mat(ii,jj)  =  coord_c(ii,1) * coord_c(jj,2)

           end do
        end do


        ! z**kk * Fp

        do kk = 0, index_ijk_max_cc(3,l_hartree_max_far_distance(current_atom)+1)
           do p = 0, l_hartree_max_far_distance(current_atom)+1

              rest_mat(kk,p)  =  coord_c(kk,3) * Fp(p,i_center)

           end do
        end do


        if (use_hartree_non_periodic_ewald) then

           drelsq = distsq / adap_outer_radius_sq

           do p = 0, l_hartree_max_far_distance(current_atom)

              Fp_deriv_red(p)  =  - 1 / sqrt(adap_outer_radius_sq)**( 2*p + 3 )           &
                                    *  ( 2*b2(p) + 4*b4(p)*drelsq + 6*b6(p)*drelsq**2 )
           end do

        end if


! X direction:

        ! x**n1 * y**n2 * z**n3 * Fp
        do n = 1,n_cc_lm_ijk(l_hartree_max_far_distance(current_atom)),1

           ii =  index_cc(n, 3)
           jj =  index_cc(n, 4)
           kk =  index_cc(n, 5)
           p  =  index_cc(n, 6)

           if (ii >= 1) then
              vector(n)  =  ii * coord_mat(ii-1,jj) * rest_mat(kk,p)
           else
              vector(n) = 0.d0
           end if

           if ( .not. use_hartree_non_periodic_ewald ) then
              vector(n) = vector(n) + coord_mat(ii+1,jj) * rest_mat(kk,p+1)
           else
              vector(n) = vector(n) + coord_mat(ii+1,jj) * coord_c(kk,3) * Fp_deriv_red(p)
           end if

        end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(1) =   gradient(1) + &
             dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(current_atom))), &
             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(current_atom)), current_atom  ))

!        write(use_unit,*) "vector: "
!        write(use_unit,*) n_cc_lm_ijk(l_hartree_max_far_distance(i_atom))
!        write(use_unit,*) vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)))
!        write(use_unit,*) "Gradient(1): "
!        write(use_unit,*)  gradient(1,i_atom)
!        write(use_unit,*) - dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom))), &
!             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(i_atom)), i_atom  ))
!        write(use_unit,*)


! Y direction:

        do n = 1,n_cc_lm_ijk(l_hartree_max_far_distance(current_atom)),1

           ii =  index_cc(n, 3)
           jj =  index_cc(n, 4)
           kk =  index_cc(n, 5)
           p  =  index_cc(n, 6)

           if (jj >= 1) then
              vector(n)  =  jj * coord_mat(ii,jj-1) * rest_mat(kk,p)
           else
              vector(n) = 0.d0
           end if

           if ( .not. use_hartree_non_periodic_ewald ) then
              vector(n) = vector(n) + coord_mat(ii,jj+1) * rest_mat(kk,p+1)
           else
              vector(n) = vector(n) + coord_mat(ii,jj+1) * coord_c(kk,3) * Fp_deriv_red(p)
           end if

        end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(2) = gradient(2) &
             + dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(current_atom))), &
             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(current_atom)), current_atom  ))


! Z direction:

        do n = 1,n_cc_lm_ijk(l_hartree_max_far_distance(current_atom)),1

           ii =  index_cc(n, 3)
           jj =  index_cc(n, 4)
           kk =  index_cc(n, 5)
           p  =  index_cc(n, 6)

           if (kk >= 1) then
              vector(n)  =  coord_mat(ii,jj) * kk * rest_mat( kk-1, p )
           else
              vector(n) = 0.d0
           end if

           if ( .not. use_hartree_non_periodic_ewald ) then
              vector(n) = vector(n) + coord_mat(ii,jj) * rest_mat(kk+1,p+1)
           else
              vector(n) = vector(n) + coord_mat(ii,jj) * coord_c(kk+1,3) * Fp_deriv_red(p)
           end if

        end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(3) =  gradient(3) &
             + dot_product(vector(1:n_cc_lm_ijk(l_hartree_max_far_distance(current_atom))), &
             multipole_c(1:n_cc_lm_ijk(l_hartree_max_far_distance(current_atom)), current_atom  ))

      end subroutine far_distance_real_gradient_hartree_potential_single_atom
!******





!-------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/far_distance_real_gradient_hartree_potential_single_atom_p2
!  NAME
!   far_distance_real_gradient_hartree_potential_single_atom_p2
!  SYNOPSIS

   subroutine far_distance_real_gradient_hartree_potential_single_atom_p2 &
     ( current_atom,  dir_tab, gradient, l_max )

!  PURPOSE
!  The subroutine calculates gradient of the real space part of the far distance Hartree potential
!  in periodic systems. This routine calculates only single atom contribution.
!  This is updated version of far_distance_real_gradient_hartree_potential_single_atom
!
    implicit none
!  ARGUMENTS

    integer :: current_atom
    real*8, dimension(3) :: dir_tab
    integer :: l_max
    real*8:: gradient(3)

!  INPUTS
!   o current_atom -- atom which contribution is calculated, absolute index
!   o dir_tab -- direction to center
!   o l_max -- maximum l component
!
!  OUTPUT
!   o  gradient -- gradient is added here.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  ! local variables
     integer i_l,i_l2

     integer i_coord
     integer:: ii, jj, kk, n
     
!  real*8:: F(0:l_pot_max)
     real*8:: coord_c(0:l_pot_max+1,3)
     real*8:: coord_mat(0:l_pot_max+1,0:l_pot_max+1)
     real*8:: rest_mat(0:l_pot_max+1,0:l_pot_max+1)
     real*8:: vector(n_cc_lm_ijk(l_pot_max+1))

     real*8, external :: ddot

! begin work

        ! x**n, y**n, and z**n

     coord_c(0,1:3) = 1.0
	do i_coord = 1,3,1
	  do i_l = 1, index_ijk_max_cc &
             (i_coord, l_max+1) !l_pot_max
           coord_c(i_l,i_coord) =  dir_tab(i_coord) *  coord_c(i_l-1,i_coord)
	   
          end do
        end do

        
        ! x**n1 * y**n2

        do i_l = 0, index_ijk_max_cc(1,l_max+1 )
           do i_l2 = 0,  index_ijk_max_cc(2,l_max+1)

              coord_mat(i_l,i_l2) =  coord_c(i_l,1)*coord_c(i_l2,2)

           end do
        end do

        ! z**n * Fp

        do i_l = 0,  index_ijk_max_cc(3,l_max+1)
           do i_l2 = 0, l_max+1

              rest_mat(i_l,i_l2) =  coord_c(i_l,3)* Fp(i_l2,1)

           end do
        end do

! X direction:

        ! x**n1 * y**n2 * z**n3 * Fp
        do n = 1,n_cc_lm_ijk(l_max),1

           ii =  index_cc(n, 3)-1
           jj =  index_cc(n, 4)
           kk =  index_cc(n, 5)

           if(ii >= 0)then
              vector(n) = (ii+1)* coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))
           else
              vector(n) = 0.d0
           end if
           ii =  index_cc(n, 3)+1

           vector(n) = vector(n) + coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6)+1)
           
        end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(1) =   gradient(1) + &
             dot_product(vector(1:n_cc_lm_ijk(l_max)), &
             multipole_c(1:n_cc_lm_ijk(l_max), current_atom  ))

! Y direction:
        do n = 1,n_cc_lm_ijk(l_max),1


           ii =  index_cc(n, 3)
           jj =  index_cc(n, 4)-1
           kk =  index_cc(n, 5)

           if(jj >= 0)then
              vector(n) = (jj+1)* coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))
           else
              vector(n) = 0.d0
           end if
              
           jj =  index_cc(n, 4)+1

           vector(n) = vector(n) + coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6)+1)
           
        end do
        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(2) = gradient(2) &
             + dot_product(vector(1:n_cc_lm_ijk(l_max)), &
             multipole_c(1:n_cc_lm_ijk(l_max), current_atom  ))

        ! Z direction:

        do n = 1,n_cc_lm_ijk(l_max),1

           ii =  index_cc(n, 3)
           jj =  index_cc(n, 4)
           kk =  index_cc(n, 5)-1

           if(kk >= 0)then
              vector(n) = (kk+1)* coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6))
           else
              vector(n) = 0.d0
           end if

           kk =  index_cc(n, 5)+1

           vector(n) = vector(n) + coord_mat(ii,jj) * rest_mat(kk,index_cc(n, 6)+1)

        end do

        ! x**n1 * y**n2 * z**n3 * Fp * cc

        gradient(3) =  gradient(3) &
             + dot_product(vector(1:n_cc_lm_ijk(l_max)), &
             multipole_c(1:n_cc_lm_ijk(l_max), current_atom  ))

      end subroutine far_distance_real_gradient_hartree_potential_single_atom_p2
!******





!-------------------------------------------------------------------------------





!****s* hartree_potential_real_p0/cleanup_hartree_potential_real
!  NAME
!   cleanup_hartree_potential_real
!  SYNOPSIS

      subroutine cleanup_hartree_potential_real

!  PURPOSE
!  Deallocates the variables in module hartree_potential_real_p0
!
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


       if (allocated(multipole_c)) then
         deallocate(multipole_c)
       end if
       if (allocated(Fp)) then
         deallocate(Fp)
       end if
       if (allocated(hartree_atoms)) then
         deallocate(hartree_atoms)
       end if

end subroutine cleanup_hartree_potential_real
!******




end module hartree_potential_real_p0
!******
