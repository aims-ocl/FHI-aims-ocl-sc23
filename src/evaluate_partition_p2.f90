!****s* FHI-aims/evaluate_partition_p2
!  NAME
!   evaluate_partition_p2
!  SYNOPSIS

subroutine evaluate_partition_p2(partition_type_temp,i_center,i_center_L, &
     dist_current,dir_current,dist_current_sq,dens_current, dist_tab,i_r,&
     dir_tab,angular_weight,partition_tab,weight_tab,partition_deriv,n_atom_list,n_compute_atoms, &
     center_index,temp_free_rho,dist_tab_sq,atom_atom_tab,atom_atom_index,min_atom_atom_tab, &
     n_atom_atom_tab,atom_atom_dist_list,atom_idx_A,atom_idx_B, &
     write_out,partition_deriv_delley )
!  PURPOSE
!  organizes the calculation of the hartree partition tab for a given point, in the *p2 formalism
!
!  USES
!
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use spline
      use free_atoms
      use constants
      use pbc_lists
      use localorb_io, only: use_unit

      implicit none


!  ARGUMENTS
      integer :: n_atom_list
      integer :: atom_atom_index(n_atom_list)
      integer :: i_center, i_center_L
      integer :: partition_type_temp
      real*8  :: dist_tab(n_atom_list)
      real*8  :: dist_tab_sq(n_atom_list)
      real*8  :: i_r(n_atom_list)

      !Be care: in fact, dir_tab only have n_compute_atoms 
      real*8  :: dir_tab( 3, n_atom_list ) 

      real*8  :: angular_weight
      real*8, dimension(n_atom_list, n_atom_list) :: atom_atom_tab
      real*8, dimension(n_atom_list) :: min_atom_atom_tab
      real*8  :: dist_current, dens_current, dist_current_sq, dir_current(3)

      integer, intent(in)                :: n_atom_atom_tab
      real*8, dimension(n_atom_atom_tab) :: atom_atom_dist_list
      integer, dimension(n_atom_atom_tab):: atom_idx_A
      integer, dimension(n_atom_list+1):: atom_idx_B

      integer:: n_compute_atoms
      integer:: center_index(n_compute_atoms)
      real*8 :: temp_free_rho(n_compute_atoms)
      
      real*8 :: partition_tab, weight_tab
      real*8 :: partition_deriv
      real*8 :: partition_deriv_delley(3,n_atoms)

      integer :: write_out

! INPUTS 
! o partition_type_temp -- type of partition tab to be used
! o i_center -- absolute index of atom in question
! o i_center_L -- atomic index for: dist_tab, i_r, temp_free_rho, dist_tab_sq etc
! o dist_current -- atomic quantities for current atom in question; might not be in the atom list at all in periodic BC ... 
! o dir_current -- direction to current atom
! o dist_current_sq -- square distance (from current atom)
! o dens_current -- density at current point
! o dist_tab -- distances to all n_compute_atoms relevant atoms   
! o i_r --  log grid index for all n_compute_atoms relevant atoms
! o dir_tab -- directions to all relevant atoms
! o angular_weight -- angular weight of point in its integration shell
! o n_atom_list -- number of entries in all list
! o n_compute_atoms -- number of entries actually relevant here
! o center_index -- index array from all i_compute_atoms to absolute center number
! o temp_free_rho -- free atom density for all relevant atoms
! o dist_tab_sq -- dist_tab*dist_tab for all relevant atoms
! o atom_atom_tab -- interatomic distances for all n_atom_list atoms, must be indexed with atom_atom_index 
! o atom_atom_index -- translation index for all i_compute_atoms to the relevant i_atom_list atoms
! o min_atom_atom_tab -- nearest-neighbour distances between atoms
! OUTPUTS
! o partition_tab -- integration weight
! o parition_deriv -- its derivative
!  
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
!
!     local variables

      logical :: valid

!     variables for trial partition function
      real*8 aux_dens(n_atom_list)
      real*8 wt_dens
      real*8 partition_norm
      real*8 partition_change

      real*8 aux_dens_deriv

      !Be care: in fact, wt_dens_deriv only have n_compute_atoms 
      real*8 wt_dens_deriv(n_atom_list)

      real*8 wt_dens_deriv_sum

      real*8, external :: stratmann_partition_weight
      real*8, external :: stratmann_weight_restricted
      real*8, external :: stratmann_weight_restricted_LN
      real*8, external :: smooth_partition_edge

!     counters

      integer i_center_2, i_center_3, i_compute_atom_2

      integer :: i_compute_atom, i_atom_atom_tab, i_atom_atom_tab_2

      real*8 :: sigma = 0.2
      real*8 :: mu = 0.3
      integer :: i_center_2nd, i_center_3rd


!     begin work 
      ! First, create safety net: Any integration point that lies inside the innermost logarithmic
      ! grid shell of another atom must be avoided because else we would have to extrapolate to
      ! zero on the logarithmic grid (but zero is the point at -infty on a log. grid, and therefore this
      ! extrapolation would have horrible consequences ...)

      ! Example of how write_out might be used. Kept as a reminder only.
      if (write_out.gt.0) then
        write(use_unit,*) "Grid point ", write_out, ": In evaluate_partition_p2()."
      end if

      valid = .true.

      do i_compute_atom = 1, n_compute_atoms, 1 
         i_center_2 = center_index(i_compute_atom)
         if (dist_tab(i_compute_atom).lt.r_grid_min(species_center(i_center_2))) then
            valid = .false.
            partition_tab = 0.d0
         end if
      end do

      ! continue only if the grid point is safe
      if (valid) then

      select case(partition_type_temp)

          case(1)
             wt_dens =  &
                  dens_current / dist_current_sq
! 		  1d0/(sigma*sqrt(pi*2d0))*exp(-0.5*((dist_current-mu)/sigma)**2)
          case(2)
             wt_dens =  &
                  dens_current / dist_current
          case(3)
             wt_dens =  &
                 dens_current
          case(4)
             wt_dens =  &
               1.0/ ( 1+  &
               exp( (dist_current- hartree_partition_parameters(1))/ &
                    hartree_partition_parameters(2) ) )
          case(5) 
             ! Stratmann partition tab
             i_atom_atom_tab = atom_atom_index(i_center_L)
             ! first check whether the point is equal to 1 
             if (dist_current.lt.min_atom_atom_tab(i_atom_atom_tab)) then
                wt_dens = 1d0
             else
                wt_dens = 1d0
                ! then check if it is equal to 0
                do i_compute_atom = 1, n_compute_atoms
                   i_atom_atom_tab_2 = atom_atom_index(i_compute_atom)
                   if (dist_tab(i_compute_atom).lt.min_atom_atom_tab(i_atom_atom_tab_2)) wt_dens = 0d0
                end do
                if (wt_dens.gt.0d0) then
                   ! No? alright, then do the long & winded calculation of its actual weight ... 
                   wt_dens = stratmann_partition_weight ( i_center_L,  dist_tab, stratmann_a, & 
                        atom_atom_tab, n_atom_list, n_compute_atoms, i_atom_atom_tab, atom_atom_index )
                end if
             end if
          case(6)
             wt_dens =  &
                  dens_current / dist_current_sq
          case(7) 
             ! Stratmann_smooth partition tab
             ! almost VERBATIM copy of case(5) - not ideal but wanted to impact remaining code as
             !   little as possible while programming on a bus.
             i_atom_atom_tab = atom_atom_index(i_center_L)
             ! first check whether the point is equal to 1 
             if (dist_current.lt.min_atom_atom_tab(i_atom_atom_tab)) then
                ! in this case, all other atoms have zero weight. Further below,
                ! we will simply skip their computation.
                wt_dens = 1d0
             else
                wt_dens = 1d0
                ! then check if it is equal to 0
                do i_compute_atom = 1, n_compute_atoms
                   i_atom_atom_tab_2 = atom_atom_index(i_compute_atom)
                   if (dist_tab(i_compute_atom).lt.min_atom_atom_tab(i_atom_atom_tab_2)) wt_dens = 0d0
                end do
                if (wt_dens.gt.0d0) then
                   ! No? alright, then do the long & winded calculation of its actual weight ... 
                   wt_dens = stratmann_partition_weight ( i_center_L,  dist_tab, stratmann_a, & 
                        atom_atom_tab, n_atom_list, n_compute_atoms, i_atom_atom_tab, atom_atom_index )
                   ! Now let the weight fall off smoothly towards the edge of the minimum atomic
                   ! radius needed for the purpose of partitioning the Hartree potential
                   ! i_center_L indexes the atomic center to which the present point belongs in list of locally relevant centers
                   ! center_index translates into the index of that center in entire system
                   ! and finally, species_center translates that to the species index in question.
                   i_center_2 = center_index(i_center_L) 
                   wt_dens = wt_dens * smooth_partition_edge ( dist_current, species_center (i_center_2) )
                end if
             end if
          case(8) 
             ! Stratmann_smoother partition tab
             ! This version ensures that atoms can not simply vanish
             ! from the list of contributing atoms without their contribution
             ! vanishing smoothly in the weight as well.
             i_atom_atom_tab = atom_atom_index(i_center_L)
             ! first check whether the point is equal to 1 
             if (dist_current.lt.min_atom_atom_tab(i_atom_atom_tab)) then
                ! in this case, all other atoms have zero weight. Further below,
                ! we will simply skip their computation.
                wt_dens = 1d0
             else
                wt_dens = 1d0
                ! then check if it is equal to 0
                do i_compute_atom = 1, n_compute_atoms
                   i_atom_atom_tab_2 = atom_atom_index(i_compute_atom)
                   if (dist_tab(i_compute_atom).lt.min_atom_atom_tab(i_atom_atom_tab_2)) wt_dens = 0d0
                end do
                if (wt_dens.gt.0d0) then
                   ! No? alright, then do the long & winded calculation of its actual weight ... 
                   wt_dens = stratmann_weight_restricted ( i_center_L,  dist_tab, stratmann_a, & 
                        atom_atom_tab, n_atom_list, n_compute_atoms, center_index, i_atom_atom_tab, & 
                        atom_atom_index, write_out ) 
                   i_center_2 = center_index(i_center_L) 
                end if
             end if
!         LN introduced as preparation for a new atom_atom_tab
          case(9)
             ! Stratmann_smoother partition tab
             ! This version ensures that atoms can not simply vanish
             ! from the list of contributing atoms without their contribution
             ! vanishing smoothly in the weight as well.
             i_atom_atom_tab = atom_atom_index(i_center_L)
             ! first check whether the point is equal to 1 
             if (dist_current.lt.min_atom_atom_tab(i_atom_atom_tab)) then
                ! in this case, all other atoms have zero weight. Further below,
                ! we will simply skip their computation.
                wt_dens = 1d0
             else
                wt_dens = 1d0
                ! then check if it is equal to 0
                do i_compute_atom = 1, n_compute_atoms
                   i_atom_atom_tab_2 = atom_atom_index(i_compute_atom)
                   if (dist_tab(i_compute_atom).lt.min_atom_atom_tab(i_atom_atom_tab_2)) wt_dens = 0d0
                end do
                if (wt_dens.gt.0d0) then
                   ! No? alright, then do the long & winded calculation of its actual weight ... 
                   wt_dens = stratmann_weight_restricted_LN ( i_center_L,  dist_tab, stratmann_a, &
                        n_atom_atom_tab, atom_atom_dist_list, atom_idx_A, atom_idx_B, &
                        n_atom_list, n_compute_atoms, center_index, i_atom_atom_tab, &
                        atom_atom_index, write_out )
                   i_center_2 = center_index(i_center_L) 
                end if
             end if
             
      end select

      ! obtain sum over all density weights at current integration point for normalisation.
      !     obtain sum over all density weights at current integration point for normalisation.
      if (partition_type_temp.eq.5) then      ! becke partition weight is treated explicitly ... 
         if ((wt_dens.gt.0d0).and.(wt_dens.lt.1d0)) then   ! only calculate interstitial weights, rest does not matter anyway - know the outcome already
            partition_norm = 0d0
            do i_compute_atom = 1, n_compute_atoms, 1 
               i_atom_atom_tab = atom_atom_index(i_compute_atom)
               partition_norm = partition_norm +  &
                    stratmann_partition_weight ( i_compute_atom,  dist_tab, stratmann_a, &
                    atom_atom_tab, n_atom_list, n_compute_atoms, i_atom_atom_tab, atom_atom_index )
            end do
         else
            partition_norm = 1d0  ! use neutral element for division for trivial cases ... 
         end if
      else if (partition_type_temp.eq.7) then      ! becke partition weight is treated explicitly ... 
         if ((wt_dens.gt.0d0).and.(wt_dens.lt.1d0)) then   ! only calculate interstitial weights, rest does not matter anyway - know the outcome already
            partition_norm = 0d0
            do i_compute_atom = 1, n_compute_atoms, 1 
               i_atom_atom_tab = atom_atom_index(i_compute_atom)
               ! Now let the weight fall off smoothly towards the edge of the minimum atomic
               ! radius needed for the purpose of partitioning the Hartree potential
               i_center_2 = center_index(i_compute_atom) 
               partition_norm = partition_norm +  &
                    stratmann_partition_weight ( i_compute_atom,  dist_tab, stratmann_a, &
                    atom_atom_tab, n_atom_list, n_compute_atoms, i_atom_atom_tab, atom_atom_index ) &
                    * smooth_partition_edge ( dist_tab(i_compute_atom), species_center (i_center_2) )
            end do
         else
            partition_norm = 1d0  ! use neutral element for division for trivial cases ... 
         end if
      else if (partition_type_temp.eq.8) then      ! becke partition weight is treated explicitly ...
         ! Stratmann_smoother partition tab
         ! This version ensures that atoms can not simply vanish
         ! from the list of contributing atoms without their contribution
         ! vanishing smoothly in the weight as well.
         if ((wt_dens.gt.0d0).and.(wt_dens.lt.1d0)) then   ! only calculate interstitial weights, rest does not matter anyway - know the outcome already
            partition_norm = 0d0
            do i_compute_atom = 1, n_compute_atoms, 1 
               i_atom_atom_tab = atom_atom_index(i_compute_atom)
               i_center_2 = center_index(i_compute_atom) 
               partition_change = stratmann_weight_restricted ( i_compute_atom,  dist_tab, stratmann_a, &
                    atom_atom_tab, n_atom_list, n_compute_atoms, center_index, i_atom_atom_tab, & 
                    atom_atom_index, write_out )
               partition_norm = partition_norm + partition_change

            end do
         else
            partition_norm = 1d0  ! use neutral element for division for trivial cases ...
         end if

      else if (partition_type_temp.eq.9) then      ! becke partition weight is treated explicitly ...
         ! Stratmann_smoother partition tab
         ! This version ensures that atoms can not simply vanish
         ! from the list of contributing atoms without their contribution
         ! vanishing smoothly in the weight as well.
         if ((wt_dens.gt.0d0).and.(wt_dens.lt.1d0)) then   ! only calculate interstitial weights, rest does not matter anyway - know the outcome already
            partition_norm = 0d0
            do i_compute_atom = 1, n_compute_atoms, 1
               i_atom_atom_tab = atom_atom_index(i_compute_atom)
               i_center_2 = center_index(i_compute_atom)
               partition_change = stratmann_weight_restricted_LN ( i_compute_atom,  dist_tab, stratmann_a, &
                         n_atom_atom_tab, atom_atom_dist_list, atom_idx_A, atom_idx_B, &
                         n_atom_list, n_compute_atoms, center_index, i_atom_atom_tab, &
                         atom_atom_index, write_out )
               partition_norm = partition_norm + partition_change
            end do

         else
            partition_norm = 1d0  ! use neutral element for division for trivial cases ... 
         end if
      else
         partition_norm = 0d0
         do i_compute_atom = 1, n_compute_atoms, 1 
            i_center_2 = center_index(i_compute_atom)

            !  add contribution from i_atom_2 to partition_norm               
            select case(partition_type_temp)
               
            case(1)
               partition_norm = partition_norm + &
! 		    1d0/(sigma*sqrt(pi*2d0))*exp(-0.5*((abs(dist_tab(i_compute_atom))-mu)/sigma)**2)
                    temp_free_rho(i_compute_atom) / dist_tab_sq(i_compute_atom)
            case(2)
               partition_norm = partition_norm + &
                    temp_free_rho(i_compute_atom)/ abs(dist_tab(i_compute_atom))
            case(3)
               partition_norm = partition_norm + &
                    temp_free_rho(i_compute_atom)
            case(4)
               partition_norm = partition_norm + &
                    1.0/ ( 1+  &
                    exp((dist_tab(i_compute_atom)-hartree_partition_parameters(1))/ &
                    hartree_partition_parameters(2) ) )
            case(6)
               partition_norm = partition_norm + &
                    temp_free_rho(i_compute_atom) / dist_tab_sq(i_compute_atom)
            end select	
         end do
      end if

      

      if (partition_norm.gt.partition_acc) then
        ! tabulate partition function
    
        partition_tab =  &
          wt_dens / partition_norm


        ! Comment - this section is obsolete and only kept for
        ! archival purposes. It is possible that someone, sometime will
        ! have to look into derivatives of partition tables for
        ! otehr reasons than the original ones.

        ! if requested, compute radial derivative of partition_tab
        if (multipole_interpolation_style.eq.1) then

          wt_dens_deriv_sum = 0.d0
  
          do i_compute_atom = 1, n_compute_atoms, 1 
             i_center_2 = center_index(i_compute_atom)

             ! only ever needed for the hartree partition table at this point!
             ! However, must compute same derivative for normal partition tab if
             ! this ever becomes relevant.
             aux_dens_deriv =  &
                  val_spline &
                  ( i_r(i_compute_atom), hartree_partition_drho_dr_spl(1,1,species_center(i_center_2)), &
                  n_grid(species_center(i_center_2)) )

            select case(partition_type_temp)
  
            case(1)

              wt_dens_deriv(i_compute_atom) =  &
               ( aux_dens_deriv -  &
                 2 * temp_free_rho(i_compute_atom) / dist_tab(i_compute_atom) ) / &
               dist_tab_sq(i_compute_atom)

            case(2)

              wt_dens_deriv(i_compute_atom) =  &
               ( aux_dens_deriv -  &
                 temp_free_rho(i_compute_atom) / dist_tab(i_compute_atom) ) / &
               dist_tab(i_compute_atom)

            case(3)

             wt_dens_deriv(i_compute_atom) = aux_dens_deriv 

            case(4)

              wt_dens_deriv(i_compute_atom) =  &
               1.0/ ( 1 + exp(  &
                 (dist_tab(i_compute_atom) - hartree_partition_parameters(1))/ &
                 hartree_partition_parameters(2) ) )

              wt_dens_deriv(i_compute_atom) = - wt_dens_deriv(i_compute_atom) * &
               1.0/ ( hartree_partition_parameters(2) * ( 1 + exp( &
               -(dist_tab(i_compute_atom)-hartree_partition_parameters(1))/ &
                    hartree_partition_parameters(2) ) ) )

            case(6)

              wt_dens_deriv(i_compute_atom) =  &
               ( aux_dens_deriv -  &
                 2 * temp_free_rho(i_compute_atom) / dist_tab(i_compute_atom) ) / &
               dist_tab_sq(i_compute_atom)

            end select

            wt_dens_deriv_sum = wt_dens_deriv_sum +  &
            wt_dens_deriv(i_compute_atom) *  &
            dot_product(dir_tab(:,i_compute_atom), dir_current(:))
          end do

          partition_deriv =  &
          ( wt_dens_deriv(i_center_L) - partition_tab*wt_dens_deriv_sum ) /  &
          partition_norm          

        !  include angular integration weights into partition_deriv
          partition_deriv =  &
            partition_deriv * &
            angular_weight * &
            4.0d0*pi

        end if


       !----------shanghui add for gradient_partition--------------
        if(use_partition_deriv) then

          do i_compute_atom = 1, n_compute_atoms, 1
             i_center_2 = center_index(i_compute_atom)

             aux_dens_deriv =  &
             val_spline &
           ( i_r(i_compute_atom), free_drho_dr_spl(1,1,species_center(i_center_2)), &
           n_grid(species_center(i_center_2)) )

             wt_dens_deriv(i_compute_atom) =  &
             ( aux_dens_deriv -  &
             2 * temp_free_rho(i_compute_atom) / dist_tab(i_compute_atom) ) / &
             dist_tab_sq(i_compute_atom)
          enddo

          do i_compute_atom = 1, n_compute_atoms, 1
             i_center_2 = center_index(i_compute_atom)

           if(i_center_2.ne.i_center)  then 
                
           partition_deriv_delley(1:3,i_center_2)= -angular_weight*4.0d0*pi* &
                          !wt_dens_deriv(i_center_2)*                        &
                          wt_dens_deriv(i_compute_atom)*                        &
                          partition_tab*                                    &
                          !(-dir_tab(1:3,i_center_2))/partition_norm
                          (-dir_tab(1:3,i_compute_atom))/partition_norm

           else  
                
             do i_compute_atom_2 = 1, n_compute_atoms, 1
                 i_center_3 = center_index(i_compute_atom_2)

                 if(i_center_3.ne.i_center ) then 
           partition_deriv_delley(1:3,i_center_2) = &  
           partition_deriv_delley(1:3,i_center_2) & 
                         -angular_weight*4.0d0*pi* &
                          !wt_dens_deriv(i_center_3)*                        &
                          wt_dens_deriv(i_compute_atom_2)*                        &
                          partition_tab*                                    &
                          !(dir_tab(1:3,i_center_3))/partition_norm
                          (dir_tab(1:3,i_compute_atom_2))/partition_norm
 
                 endif
             enddo

           endif
          
          enddo 
       endif 

       !----------shanghui add for gradient_partition--------------

       weight_tab = angular_weight * 4d0*pi
!       and multiply angular integration weights into partition_tab
        partition_tab =  &
          partition_tab * &
          weight_tab

      else
!       this grid point is too far away - we should not be integrating here!
        partition_tab = 0.d0
        if (multipole_interpolation_style.eq.1) then
           partition_deriv = 0.d0
        end if


        if(use_partition_deriv) then
           partition_deriv_delley(1:3,1:n_atoms) = 0.0d0 
        endif

      end if


      end if


    end subroutine evaluate_partition_p2

!******

