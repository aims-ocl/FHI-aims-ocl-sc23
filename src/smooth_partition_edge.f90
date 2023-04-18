!----------------------------------------------------------------------------------------------------------------
!
!  function smooth_partition_edge gives a function that is 
!
!  - one inside a sphere of radius = factor * multipole_radius_free , 0 < factor < 1
!  - zero outside radius = multipole_radius_free
!  - smoothly decays to zero for factor * multipole_radius_free < radius < multipole_radius_free
!
!  VB: To avoid discontinuities, multipole_radius_free has been replaced by the rigorous outer bound
!      for each species, free_r_cut(i_species) + w_cut(i_species), codified in outer_partition_radius(i_species).
!      If this causes any trouble (the multipole_radius_free assumption is wired deeply into the code, especially
!      in pbc_lists.f90), we have to revisit this choice.
!
!  This enables the edge of the Stratmann partition table for a given atom to decay to zero smoothly
!  at a bounded distance multipole_radius_free, regardless of the interatomic distance.
!
!----------------------------------------------------------------------------------------------------------------

     real*8 function smooth_partition_edge & 
     ( current_atom_distance, i_species )

       use species_data, only : outer_partition_radius

       implicit none

       real*8 :: current_atom_distance
       integer :: i_species

       real*8, parameter :: factor = 0.8d0
       real*8, parameter :: pi      = 3.14159265358979323846d0 

       if ( current_atom_distance.le.(factor*outer_partition_radius(i_species)) ) then

         smooth_partition_edge = 1.d0

       else if ( current_atom_distance.ge.outer_partition_radius(i_species) ) then

         smooth_partition_edge = 0.d0

       else
         ! this is the interesting case

         smooth_partition_edge = current_atom_distance / outer_partition_radius(i_species)
         smooth_partition_edge =  pi * ( smooth_partition_edge - factor ) / (1.d0 - factor) 
         smooth_partition_edge =  ( cos(smooth_partition_edge) + 1.d0 ) / 2.d0 

       end if

     end function smooth_partition_edge
