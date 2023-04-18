 !****s* FHI-aims/evaluate_partition_tab_p2
!  NAME
!   evaluate_partition_tab_p2
!  SYNOPSIS

 subroutine evaluate_partition_tab_p2 &
           ( i_center,        &  ! absolute index of atom in question
             i_center_L,      &  ! atomic index for: dist_tab, i_r, temp_free_rho, dist_tab_sq etc
             dist_current,    &  ! quantities for the current atom in question
             dist_current_sq, &
             dens_current,    &
             dist_tab,        &  ! distances to all n_compute_atoms relevant atoms
             i_r,             &  ! log grid index for all n_compute_atoms relevant atoms
             radial_weight,   &  ! radial & angular weights
             angular_weight,  &  
             partition_tab,   &  ! THE OUTPUT result ... 
             weight_tab,      &  ! only the integration weights without the partition of unity
                                 !function
             n_atom_list,     &  ! number of entries in all the lists
             n_compute_atoms, &  ! number of relevant atoms
             center_index,    &  ! index array from all i_compute_atoms to absolute center number
             temp_free_rho,   &  ! free atom density for all relevant atoms
             dist_tab_sq,     &  ! dist_tab*dist_tab for all relevant atoms
             atom_atom_tab,   &  ! interatomic distances for all n_atom_list atoms, must be indexed with atom_atom_index 
             atom_atom_index, &  ! translation index for all i_compute_atoms to the relevant i_atom_list atoms
             min_atom_atom_tab,& ! NN distances of all atoms
             n_atom_atom_tab, &
             atom_atom_dist_list,&
             atom_idx_A,      &
             atom_idx_B)

!  PURPOSE
!     calculate partition tab for integrations, given a certain point and all the relevant centers that might affect it 
!
!  USES
!
        use dimensions
        use runtime_choices
        use grids
        use pbc_lists
        use spline
        use free_atoms
        use constants
        
        implicit none

 !  ARGUMENTS
        integer :: n_atom_list
        integer :: i_center
        integer :: i_center_L
        integer :: i_atom_atom_tab, i_atom_atom_tab_2
        integer :: atom_atom_index(n_atom_list)
        real*8  :: dir_current(3), dist_current, dist_current_sq, dens_current
        real*8  :: dist_tab(n_atom_list)
        real*8  :: dist_tab_sq(n_atom_list)
        real*8  :: i_r(n_atom_list)
        real*8  :: radial_weight
        real*8  :: angular_weight
        real*8  :: atom_atom_tab(n_atom_list, n_atom_list)

        integer, intent(in)                :: n_atom_atom_tab
        real*8, dimension(n_atom_atom_tab) :: atom_atom_dist_list
        integer, dimension(n_atom_atom_tab):: atom_idx_A
        integer, dimension(n_atom_list+1):: atom_idx_B

        real*8  :: min_atom_atom_tab(n_atom_list)
        integer :: n_compute_atoms
        integer :: center_index(n_compute_atoms)
        real*8  :: temp_free_rho(n_compute_atoms)
        real*8  :: partition_tab
        real*8  :: weight_tab

!  INPUTS
! o i_center -- absolute index of atom in question
! o i_center_L -- atomic index for: dist_tab, i_r, temp_free_rho, dist_tab_sq etc
! o dist_current -- distance to for the current atom in question
! o dist_current_sq -- distance to the atom squared
! o dens_current -- density at current point
! o dist_tab -- distances to all n_compute_atoms relevant atoms
! o i_r -- log grid index for all n_compute_atoms relevant atoms
! o radial_weight -- radial & angular weights
! o angular_weight -- Lebedev weight of point in integration shell
! o n_atom_list -- number of entries in all the lists
! o n_compute_atoms -- number of relevant atoms
! o center_index -- index array from all i_compute_atoms to absolute center number
! o temp_free_rho -- free atom density for all relevant atoms
! o dist_tab_sq -- dist_tab*dist_tab for all relevant atoms
! o atom_atom_tab -- interatomic distances for all n_atom_list atoms, must be indexed with atom_atom_index 
! o atom_atom_index -- translation index for all i_compute_atoms to the relevant i_atom_list atoms
! o min_atom_atom_tab -- NN distances of all atoms
!  OUTPUT
! o partition_tab -- integration weight of point in question 
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
        
        integer :: write_out = 0

        !     variables for trial partition function
        real*8 :: aux_dens
        real*8 :: wt_dens
        real*8 :: partition_norm
        real*8 :: partition_change
      
        !     counters
        
        integer :: i_center_2
        
        integer :: i_compute_atom
        real*8, external :: stratmann_partition_weight
        real*8, external :: stratmann_weight_restricted
        real*8, external :: stratmann_weight_restricted_LN
        real*8, external :: smooth_partition_edge


        !     begin work
        
        ! First, create safety net: Any integration point that lies inside the innermost logarithmic
        ! grid shell of another atom must be avoided because else we would have to extrapolate to
        ! zero on the logarithmic grid (but zero is the point at -infty on a log. grid, and therefore this
        ! extrapolation would have horrible consequences ...)

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
           
           ! initialize numerator / denominator of partition function
           select case(partition_type)
	   
            
           case(1)
              wt_dens = &
                   dens_current / dist_current_sq 
           case(2)
              wt_dens = &
                   dens_current / dist_current
           case(3)
              wt_dens = &
                   dens_current
           case(4)
              wt_dens = &
                   1.0/ ( 1+ &
                   exp( (dist_current- integral_partition_parameters(1))/ &
                   integral_partition_parameters(2) ) )
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
              wt_dens = &
                   dens_current / dist_current_sq
           case(7) 
             ! Stratmann_smooth partition tab
             ! almost VERBATIM copy of case(5) - not ideal but wanted to impact remaining code as
             !   little as possible while programming on a bus.
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
                   wt_dens = stratmann_weight_restricted_LN ( i_compute_atom,  dist_tab, stratmann_a, &
                         n_atom_atom_tab, atom_atom_dist_list, atom_idx_A, atom_idx_B, &
                         n_atom_list, n_compute_atoms, center_index, i_atom_atom_tab, &
                         atom_atom_index, write_out )
                   i_center_2 = center_index(i_center_L)
                end if
             end if

           end select
           
           !     obtain sum over all density weights at current integration point for normalisation.
           if (partition_type.eq.5) then      ! becke partition weight is treated explicitly ... 
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
           else if (partition_type.eq.7) then      ! becke partition weight is treated explicitly ... 
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
           else if (partition_type.eq.8) then      ! becke partition weight is treated explicitly ...
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
           else if (partition_type.eq.9) then      ! becke partition weight is treated explicitly ...
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
              partition_norm = wt_dens
              do i_compute_atom = 1, n_compute_atoms, 1 
                 i_center_2 = center_index(i_compute_atom)
                                     
                 select case(partition_type)
                    
                 case(1)
                    partition_norm = partition_norm + &
                         temp_free_rho(i_compute_atom) / dist_tab_sq(i_compute_atom) 
                 case(2)
                    partition_norm = partition_norm + &
                         temp_free_rho(i_compute_atom) / dist_tab(i_compute_atom)
                 case(3)
                    partition_norm = partition_norm + &
                         temp_free_rho(i_compute_atom)
                 case(4)
                    partition_norm = partition_norm + &
                         1.0/ ( 1+  &
                         exp((dist_tab(i_compute_atom)-integral_partition_parameters(1))/ &
                         integral_partition_parameters(2) ) )
                 case(6)
                    partition_norm = partition_norm + &
                         temp_free_rho(i_compute_atom) / dist_tab_sq(i_compute_atom)
                 end select
              enddo
           end if
           if (partition_norm.gt.partition_acc) then
              !       tabulate partition function
              partition_tab = &
                   wt_dens / partition_norm
              weight_tab = radial_weight * &
                  dist_tab_sq(i_center_L) * &
                  angular_weight * &
                  4.0d0*pi
              !       and multiply integration weights into partition_tab
              partition_tab = &
                   partition_tab * &
                   weight_tab
           else
              !       this grid point is too far away - we should not be integrating here!
              partition_tab = 0.d0
           end if
           
           ! end if (valid)
        end if
 end subroutine evaluate_partition_tab_p2
      
!******

