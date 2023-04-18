subroutine embedding_potential(coord_current,pot_ion_embed,force)

! calculates potential due to a homogeneous field or 
! due to a "image-plane" (1/4z) potential or
! interaction of point charges with ionic charges 

use dimensions
use geometry
use grids
use localorb_io
use runtime_choices !OTH - for vacuum-z-level; to be shifted to geometry
use mpi_tasks, only: aims_stop

implicit none

! Imported variables

real*8,dimension(3) :: coord_current
real*8 :: pot_ion_embed 
real*8, dimension(3) :: force

! Local Variables

real*8 :: dipole_tmp 
real*8 :: distance 
real*8, dimension(3) :: coord_diff
 
!Counters

integer :: i_multipole
integer :: i_coord

!OTH: For periodic part
real*8 :: dip_length
real*8 :: gradient
real*8 :: dip_origin
real*8 :: pot_jump
real*8, dimension(3) :: coord_temp
character*100 :: info_str

!_______________________________________

  ! VB Fixme: Can this loop be vectorized explicitly?
  pot_ion_embed = 0.d0
  force(:) = 0.0

!Make sure this part is not accidently executed for PBC
if (n_periodic.eq.0) then
  do i_multipole = 1, n_multipoles, 1

     distance=0.d0

     do i_coord = 1,3,1
           coord_diff(i_coord) = coord_current(i_coord) - multipole_coords(i_coord,i_multipole) 

           distance = distance + coord_diff(i_coord)**2.0d0
     end do

     distance = sqrt(distance)

     selectcase(multipole_order(i_multipole))
       case(0)
         ! negative sign because positive charges must be attractive!
         pot_ion_embed = pot_ion_embed &
                - multipole_charge(i_multipole)/distance
         ! AT: force
         if (use_forces) then 
           do i_coord = 1,3,1
             force(i_coord) = force(i_coord) - multipole_charge(i_multipole) &
                              * coord_diff(i_coord) / distance**3.d0
           end do
         endif

       case(1)
         dipole_tmp=0.d0
         do i_coord = 1,3,1
           dipole_tmp = dipole_tmp + ((coord_current(i_coord) - multipole_coords(i_coord,i_multipole)) * &
                      multipole_data(i_coord,i_multipole))
         end do

         pot_ion_embed = pot_ion_embed &
                + multipole_charge(i_multipole)*dipole_tmp/distance**3.d0

         ! AT: force 
         if (use_forces) then
           do i_coord = 1,3,1
             force(i_coord) = force(i_coord) - (multipole_charge(i_multipole) * & 
                              multipole_data(i_coord,i_multipole)*coord_diff(i_coord)) / distance**5.d0
           end do
         endif 
        endselect
  end do
endif !n_periodic.eq.0 

!OTH: Adapting to periodic fields
if (n_periodic.eq.0) then
  
   pot_ion_embed = pot_ion_embed &
    - dot_product(homogeneous_field,coord_current) 

elseif (n_periodic.eq.3) then
! BL: In case of periodic calculations a asymmetric sawtooth potential is
!     created with a potential jump at z = vacuum_z_level
!     The field has to be solely in z direction, which is checked now. 
     coord_temp = coord_current
     call map_to_center_cell(coord_temp) !Just to be sure...
     if (allocated(homogeneous_field)) then !apply effects of field

         ! Ensure that field is in z direction
         if (dabs(homogeneous_field(1)) + dabs(homogeneous_field(2)) .gt. 1.0d-10) then
            write(info_str,'(2X,2A)') & 
            "** Homogeneous field in other directions then z not supported!",&
            " Please adjust simulation cell"
            call aims_stop(info_str,"embedding_potential.f90")
         end if

         gradient = homogeneous_field(3)
         ! BL: Compare z coordinates of all 3 lattice vectors to find cell extension in z direction
         ! length = MAX_Z(a1,a2,a3) - MIN_Z(a1,a2,a3) 
         dip_length = dabs((maxval(lattice_vector(3,:)) - minval(lattice_vector(3,:))))
         ! BL: Put dipole origin (V = 0) half length away from vacuum level (V jump)
         dip_origin = vacuum_z_level + 0.5d0 * dip_length
         ! BL: V_jump = E_z * d
         pot_jump = gradient * dip_length
 
        ! E = - grad phi 
        ! z smaller than vacuum level: phi = - 0.5 V_jump - (z - z0) * gradient
        ! "From Zero to negative jump"
        if (coord_current(3) .lt. vacuum_z_level) then
          pot_ion_embed= -0.5d0 * pot_jump - (coord_current(3) - vacuum_z_level) * gradient
        ! z larger than vacuum level: phi = + 0.5 V_jump - (z - z0) * gradient
        ! "From positive jump to Zero"
        elseif (coord_current(3) .ge. vacuum_z_level) then
            pot_ion_embed= 0.5d0 * pot_jump - (coord_current(3) - vacuum_z_level) * gradient
        endif

     endif

     !OTH apply effects of mirror image potential
     if (potential_well_requested) then !add a potential well as requested by the user. Like the dipole correction, this is in z only.
         if (coord_temp(3) .gt. potential_well_start .and. coord_temp(3) .lt. potential_well_end) then
               pot_ion_embed = pot_ion_embed + potential_well_depth
         endif
          
     endif

endif

end subroutine

subroutine output_embedding_potential_z(z_min, z_max, n_points)

!  PURPOSE
!    The subroutine plots the embedded potential along the z axis!
!  USES

use constants, only: hartree_over_bohr
use dimensions
use geometry
use grids
use localorb_io
use runtime_choices

  implicit none

! ARGUMENTS
  real*8, intent(IN)   :: z_min 
  real*8, intent(IN)   :: z_max 
  integer, intent(IN)  :: n_points 

! LOCAL VARIABLE
  integer              :: i
  real*8, dimension(3) :: coord_temp
  real*8               :: pot_ion_embed
  real*8, dimension(3) :: force
  character*100 :: info_str

! START

write(info_str,*) "Embedding Potential:"
call localorb_info(info_str,use_unit,'(A)')
write(info_str,*) "z (A)      Potential (eV)     Force_z (eV/A)"
call localorb_info(info_str,use_unit,'(A)') 

do i = 1, n_points
  coord_temp(1) = 0.0
  coord_temp(2) = 0.0
  coord_temp(3) = (i-1) * (z_max-z_min) / (n_points-1) + z_min
  call embedding_potential(coord_temp, pot_ion_embed, force)

  write(info_str,*) coord_temp(3)*bohr, pot_ion_embed*hartree, force(3)*hartree_over_bohr
  call localorb_info(info_str,use_unit,'(A)') 

end do


end subroutine
