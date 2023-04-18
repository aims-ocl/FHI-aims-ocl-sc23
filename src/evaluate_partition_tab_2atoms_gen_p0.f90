!****s* FHI-aims/evaluate_partition_tab_2atoms_gen_p0
!  NAME
!   evaluate_partition_tab_2atoms_gen_p0
!  SYNOPSIS

      subroutine evaluate_partition_tab_2atoms_gen_p0 &
      ( i_atom, dist_tab, i_r, &
        radial_weight, angular_weight, &
        partition_tab, partition_type_temp, &
        n_atoms_list_max, n_atoms_list, atoms_list,&
        atoms_list_index)

!  PURPOSE
!     Subroutine evaluate_partition_tab_2atoms_gen_p0
!     creates partition_tab for the present integration point among
!     2 atoms. These 2 atoms can be either identical or different.

!     Hrmph. We must use module grids here explicitly, although some
!     of the explicit input to subroutine evaluate_partition is also
!     contained in module grids. Hopefully this will not be too confusing.

!     What needs to be done is, first determine the set of atoms which
!     matters at this point.

!  USES

      use dimensions
      use runtime_choices
      use grids
      use geometry
      use spline
      use free_atoms
      use constants
      use pbc_lists

      implicit none

!  ARGUMENTS

      integer i_atom
      integer partition_type_temp
      real*8  dist_tab(n_atoms_list_max)
      real*8  i_r(n_atoms_list_max)
      real*8  radial_weight
      real*8  angular_weight
      integer n_atoms_list_max
      integer n_atoms_list
      integer atoms_list(n_atoms_list_max)
      integer atoms_list_index(n_centers)

      real*8 partition_tab(n_atoms_list_max)
      

!  INPUTS
!  o  i_atom -- the current atom that is under consideration
!  o  partition_type_temp -- integer number, tells what partition type is used here
!  o  dist_tab(n_atoms) -- the distances of the current integration point to all atoms 
!  o  i_r -- the distancs of the current points to all atoms with the scale of
!          the logrithmic grid.
!  o  radial_weight -- the weight of the radial integration for the current points
!  o  augular_weight -- the weight of the angular integration for the current points
!  OUTPUT
!  o  partition_tab -- the 2-center partition functions
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

!    offset parameter for partition function g(r) = rho(r)/r^2 + par_offset
      real*8, parameter :: par_offset = 0.d-2

!     variables for trial partition function
      real*8 aux_dens
      real*8 wt_dens
      real*8 aux_partition_norm
      real*8 aux_partition_norm_2
      real*8 partition_norm_start
      real*8 partition_norm

!     counters

      integer i_atom_1, i_atom_1_ind
      integer i_atom_2, i_atom_2_ind
      integer i_index
      integer i_atom_ind

!     begin work

!     initialize numerator / denominator of partition function

      i_atom_ind = atoms_list_index(i_atom)

      aux_dens = &
        val_spline &
        ( i_r(i_atom_ind), partition_rho_spl(1,1,species_center(i_atom)), &
          n_grid(species_center(i_atom)) )


         select case(partition_type_temp)

          case(1)
             wt_dens = &
                  aux_dens / (dist_tab(i_atom_ind)**2.0d0) + par_offset
          case(2)
             wt_dens = &
                  aux_dens / abs(dist_tab(i_atom_ind))
          case(3)
             wt_dens = &
                  aux_dens
          case(4)
             wt_dens = &
               1.0/ ( 1+ &
               exp( (dist_tab(i_atom_ind)- hartree_partition_parameters(1))/ &
                    hartree_partition_parameters(2) ) )
          case(6)
             wt_dens = &
                  aux_dens / (dist_tab(i_atom_ind)**2.0d0) + par_offset
          end select


!      wt_dens = aux_dens /
!     +  (dist_tab(i_atom_ind))**2.0d0

!     obtain sum over all density weights at current integration point for normalisation.

      partition_norm_start = wt_dens

      i_index = 0
!    loop over the first atom
      do i_atom_1_ind = 1, n_atoms_list, 1
         i_atom_1 = atoms_list(i_atom_1_ind)
        if (i_atom_1 .ne. i_atom ) then
           aux_dens =   val_spline &
            ( i_r(i_atom_1_ind), partition_rho_spl(1,1,species_center(i_atom_1)), &
             n_grid(species_center(i_atom_1)) )

           select case(partition_type_temp)

            case(1)
               aux_partition_norm = &
                   aux_dens / (dist_tab(i_atom_1_ind)**2.0d0) + par_offset
            case(2)
               aux_partition_norm = &
                   aux_dens / abs(dist_tab(i_atom_1_ind))
            case(3)
               aux_partition_norm = &
                    aux_dens
            case(4)
               aux_partition_norm = &
                  1.0/ ( 1+ exp( (dist_tab(i_atom_1_ind)- &
                        hartree_partition_parameters(1))/ &
                        hartree_partition_parameters(2) ) )
            case(6)
               aux_partition_norm = &
                   aux_dens / (dist_tab(i_atom_1_ind)**2.0d0) + par_offset
           end select

!     end of if (i_atom_1 .ne. i_atom)
        endif

        i_atom_2_ind = i_atom_1_ind
        i_atom_2 = i_atom_1
        partition_norm = partition_norm_start
        
        i_index = i_atom_1_ind
           
        if (i_atom_1.ne.i_atom ) then
           
           partition_norm = partition_norm_start + aux_partition_norm
           
        endif
        
!          write(use_unit,*)i_index
!          write(use_unit,*)i_atom_1,i_atom_2
!          write(use_unit,*)partition_norm_start
!          write(use_unit,*)aux_partition_norm
!          write(use_unit,*)aux_partition_norm_2
!          write(use_unit,*)partition_norm


        if (partition_norm.gt.partition_acc) then
           !              tabulate partition function
           partition_tab(i_index) = &
                wt_dens / partition_norm
           
           !         and multiply integration weights into partition_tab
           partition_tab(i_index) = &
                partition_tab(i_index) * &
                radial_weight * &
                (dist_tab(i_atom_ind))**2.d0 * &
                angular_weight * &
                4.0d0*pi
        else
           
           !         this grid point is too far away - we should not be integrating here!
           partition_tab(i_index) = 0.d0
           
        end if
        
        
        !     end loop over other atoms, i_atom_1_ind
     enddo

   end subroutine evaluate_partition_tab_2atoms_gen_p0

!---------------------------------------------------------------------
!******
