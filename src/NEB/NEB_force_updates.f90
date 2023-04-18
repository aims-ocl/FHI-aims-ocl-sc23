!
! A collection of subroutines that do all the technical stuff in the NEB
!
!  contains:
!  subroutine obtain_variable_spring_constants
!  subroutine get_max_force_image
!  subroutine update_force_PEB
!  subroutine update_force_NEB
!  subroutine update_force_CINEB
!


!---------------------------------------------------------------------------------------------------
! calculate variable spring constants according to prescription in Henkelman,Uberuaga,Jonsson; JCP v113 p9901
subroutine obtain_variable_spring_constants(n_images, spring_constants, start_energy, end_energy, energies, kmax, kmin)
implicit none
integer :: n_images, i_image
real*8, dimension(n_images) :: spring_constants, energies
real*8 :: start_energy, end_energy, kmax, kmin, Emax, Emin
! determine reference_energies
Emin = max(start_energy,end_energy)
Emax = Emin
do i_image = 1, n_images
   Emax = max(Emax, energies(i_image))
end do
write(use_unit,*) 'Spring constants on the images: '
do i_image = 1, n_images
   if (energies(i_image).lt.Emin) then
      spring_constants(i_image) = kmin
   else
      spring_constants(i_image) = kmin + (kmax-kmin)*(energies(i_image)-Emin)/(Emax-Emin)
   end if
   write(use_unit,*) 'Image: ',i_image,'spring constant ', spring_constants(i_image)
end do
end subroutine obtain_variable_spring_constants

!---------------------------------------------------------------------------------------------------
!  calculates the maximal force on i_image^th image, helps to decide the convergence of transition
!  state!
subroutine get_max_force_image(n_atoms,n_images,i_image,forces,max_force_image)
  implicit none
  integer :: n_atoms, n_images, i_image, i_coords
  real*8 :: max_force_image
  real*8, dimension(3*n_atoms,n_images) :: forces
  max_force_image = 0d0
  do i_coords = 1, 3*n_atoms
     max_force_image = max(max_force_image,abs(forces(i_coords,i_image)))
  end do
end subroutine get_max_force_image

!---------------------------------------------------------------------------------------------------
! calculates net force on the system.
! in a cluster, net force should ideally vanish
! dft-forces are noisy, so it doesn't
! noise is removed by subroutine "remove_translation_and_rotation"
!
! R.Gehrke (2007); adapted for NEB by Felix Hanke (2007)
subroutine get_net_force(n_atoms, n_images, total_forces, net_force)
  implicit none
  ! imported variables
  real*8, dimension(3, n_atoms, n_images), intent(in) :: total_forces
  real*8, dimension(3, n_images), intent(out) :: net_force
  integer :: i_atom, i_image, n_atoms, n_images

  net_force = 0.d0
  do i_image = 1, n_images
     do i_atom = 1, n_atoms, 1
        net_force(:,i_image) = net_force(:,i_image) + total_forces(:, i_atom, i_image)
     end do
  end do
end subroutine get_net_force

!------------------------------------------------------------------------------------------------------
! Clean the translational force components from each image, causing the images to remain where they are
! Adapted from the original AIMS counterpart, based on work by R Gehrke ... 
subroutine clean_translational_forces(n_atoms, n_images, total_forces)
  implicit none
  integer :: n_atoms, n_images, i_atom, i_coord, i_image
  real*8, dimension(3,n_atoms,n_images) :: total_forces
  real*8, dimension(3,n_atoms,3)        :: translation_vectors
  real*8 :: inv_norm, dnrm2, ddot, translation_component
  translation_vectors = 0.d0
  do i_atom = 1, n_atoms, 1
     do i_coord = 1, 3, 1
        translation_vectors(i_coord, i_atom, i_coord) = 1.d0
     end do
  end do
  do i_coord =1, 3, 1
     inv_norm = 1.d0 / dnrm2(3*n_atoms, translation_vectors(:, :, i_coord), 1)
     translation_vectors(:, :, i_coord) = translation_vectors(:, :, i_coord) * inv_norm
  end do
  do i_image = 1, n_images
     do i_coord = 1, 3, 1
        translation_component = ddot(3*n_atoms, translation_vectors(1,1,i_coord), 1, total_forces(1,1,i_image), 1)
        total_forces(:,:,i_image) = total_forces(:,:,i_image) - translation_component * translation_vectors(:,:,i_coord)
     end do
  end do
end subroutine clean_translational_forces

!---------------------------------------------------------------------------------------------------
! A straight forward routine implementing the pure elastic band method, with externally specified
! spring constant
subroutine update_force_PEB(n_atoms, n_images, start_coords, coords, end_coords, &
     start_forces, total_forces, end_forces, &
     start_energy, energies, end_energy, spring_constants, object_function)
  implicit none
  integer :: n_atoms, n_images, i_images, i_coords
  real*8  :: start_energy, end_energy, object_function, buf, ddot
  real*8, dimension(3*n_atoms) :: start_coords, end_coords, start_forces, end_forces
  real*8, dimension(n_images)  :: energies, spring_constants
  real*8, dimension(3*n_atoms,n_images) :: coords, total_forces, spring_forces
  ! calculate total energy of elastic band ... first the DFT energies, 
  object_function = 0d0
  do i_images = 1, n_images
     object_function = object_function + energies(i_images)
  end do
  ! use the lower one of the end points as reference:
  buf = min(start_energy,end_energy)
  object_function = object_function - dble(n_images)*buf
  ! ... then the spring_energies.
  ! Spring 1: start <-> first image
  buf = spring_constants(1)*ddot(3*n_atoms,start_coords(:)-coords(:,1),1,start_coords(:)-coords(:,1),1)
  ! intermediate images
  do i_images = 2, n_images
     buf = buf + spring_constants(i_images)*ddot(3*n_atoms,coords(:,i_images)-coords(:,i_images-1),1,coords(:,i_images)-coords(:,i_images-1),1)
  end do
  ! last image
  buf = buf + spring_constants(n_images)*ddot(3*n_atoms,coords(:,n_images)-end_coords(:),1,coords(:,n_images)-end_coords(:),1)
  object_function = object_function + buf/2d0
  ! calculate spring forces
  if (n_images.gt.1) then
     ! first image
     spring_forces(:,1) = spring_constants(1)*(start_coords(:)+coords(:,2)-2d0*coords(:,1))
     ! intermediate images
     do i_images = 2, n_images-1
        spring_forces(:,i_images) = spring_constants(i_images)*(coords(:,i_images+1)+coords(:,i_images-1)-2d0*coords(:,i_images))
     end do
     ! last image
     spring_forces(:,n_images) = spring_constants(n_images)*(end_coords(:)+coords(:,n_images-1)-2d0*coords(:,n_images))
  else 
     ! one image only:
     spring_forces(:,1) = spring_constants(1)*(start_coords(:)+end_coords(:)-2d0*coords(:,1))
  end if
  ! ... and now, add the results
  total_forces(:,:) = total_forces(:,:) + spring_forces(:,:)
  call clean_translational_forces(n_atoms,n_images,total_forces)
end subroutine update_force_PEB

!---------------------------------------------------------------------------------------------------
! nudged elastic band routine according to the prescription by Jonsson, Mills, and Jacobsen
subroutine update_force_NEB(n_atoms, n_images, start_coords, coords, end_coords, &
     start_forces, total_forces, end_forces, &
     start_energy, energies, end_energy, spring_constants)
  implicit none
  integer n_atoms, n_images
  real*8 :: start_energy, end_energy, pi
  real*8, dimension(3*n_atoms) :: start_coords, end_coords, start_forces, end_forces, tauplus, tauminus
  real*8, dimension(n_images) :: energies, spring_constants
  real*8, dimension(3*n_atoms,n_images) :: coords, total_forces
  ! working variables
  real*8, dimension(3*n_atoms,n_images) :: spring_forces, correction_forces, tangent, force_perp, force_par, force_corr
  real*8, dimension(n_images) :: tangent_magnitude, cos_phi, f_phi
  integer :: i_coords, i_images, i_atoms
  real*8 :: V_first, V_second, ddot, buf

  ! calculate tangent vectors, according to the 'improved tangent vector' prescription 
  ! by Henkelmann and Jonsson, JCP v113 p9978, eqns (8) and eqns (10)
  tangent  = 0d0
  ! image 1
  if ((start_energy.lt.energies(1)).and.(energies(1).lt.energies(2))) then
     ! V_start < V_1 < V_2; use tau+ = R_2 - R_1
     tangent(:,1) = coords(:,2) - coords(:,1)
  else if ((start_energy.gt.energies(1)).and.(energies(1).gt.energies(2))) then
     ! V_start> V_1 > V_2; use tau- = R_1 - R_start
     tangent(:,1) = coords(:,1) - start_coords(:)
  else
     ! use energy-difference weighed average of the two vectors in the previous cases. 
     tauplus(:)  = coords(:,2) - coords(:,1)
     tauminus(:) = coords(:,1) - start_coords(:)
     if (energies(2).gt.start_energy) then
        V_first  = max(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
        V_second = min(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
     else
        V_first  = min(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
        V_second = max(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
     end if
     tangent(:,1) = V_first*tauplus(:) + V_second*tauminus(:)
  end if
  ! images 2 ... n-1
  do i_images = 2, n_images-1
     if ((energies(i_images-1).lt.energies(i_images)).and.(energies(i_images).lt.energies(i_images+1))) then
        tangent(:,i_images) = coords(:,i_images+1)-coords(:,i_images)
     else if ((energies(i_images-1).gt.energies(i_images)).and.(energies(i_images).gt.energies(i_images+1))) then
        tangent(:,i_images) = coords(:,i_images)-coords(:,i_images-1)
     else
        tauplus(:)  = coords(:,i_images+1)-coords(:,i_images  )
        tauminus(:) = coords(:,i_images  )-coords(:,i_images-1)
        if (energies(i_images+1).gt.energies(i_images-1)) then
           V_first  = max(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
           V_second = min(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
        else
           V_first  = min(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
           V_second = max(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
        end if
!        write(use_unit,*) 'DEBUG: ', i_images, V_first, V_second
        tangent(:,i_images) = V_first*tauplus(:) + V_second*tauminus(:)
!        write(use_unit,*) 'DEBUG: intermediate tangent vector!'
!       write(use_unit,'("DEBUG intermediate: ",3E18.8E4)')tangent(:,i_images)

     end if
  end do
  ! last image
  if ((energies(n_images-1).lt.energies(n_images)).and.(energies(n_images).lt.end_energy)) then
     tangent(:,n_images) = end_coords(:) - coords(:,n_images)
  else if ((energies(n_images-1).gt.energies(n_images)).and.(energies(n_images).gt.end_energy)) then
     tangent(:,n_images) = coords(:,n_images) - coords(:,n_images-1)
  else
     tauplus (:) = end_coords(:) - coords(:,n_images)
     tauminus(:) = coords(:,n_images) - coords(:,n_images-1)
     if (end_energy.gt.energies(n_images-1)) then
        V_first  = max(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))
        V_second = min(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))
     else
        V_first  = min(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))
        V_second = max(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))        
     end if
     tangent(:,i_images) = V_first*tauplus(:) + V_second*tauminus(:)
  end if
  ! normalize tangents:
  do i_images = 1, n_images
     tangent_magnitude(i_images) = sqrt(ddot(3*n_atoms,tangent(:,i_images),1,tangent(:,i_images),1))
     tangent(:,i_images) = tangent(:,i_images)/tangent_magnitude(i_images)
  end do

! write(use_unit,*) 'DEBUG: Tangent vectors:'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') tangent 

  ! calculate spring forces
  ! first image
  spring_forces(:,1) = spring_constants(1)*(start_coords(:)+coords(:,2)-2d0*coords(:,1))
  ! intermediate images
  do i_images = 2, n_images-1
     spring_forces(:,i_images) = spring_constants(i_images)*(coords(:,i_images+1)+coords(:,i_images-1)-2d0*coords(:,i_images))
  end do
  ! last image
  spring_forces(:,n_images) = spring_constants(n_images)*(end_coords(:)+coords(:,n_images-1)-2d0*coords(:,n_images))

! write(use_unit,*) 'DEBUG: spring forces:'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') spring_forces
  ! calculate (8) in book chapter
  do i_images = 1, n_images
     ! the perpendicular force, first term
     buf                    = ddot(3*n_atoms,total_forces(:,i_images),1,tangent(:,i_images),1)
     force_perp(:,i_images) = total_forces(:,i_images) - buf*tangent(:,i_images)
     ! the parallel force, second term: based on spring forces
     buf                    = ddot(3*n_atoms,spring_forces(:,i_images),1,tangent(:,i_images),1)
     force_par(:,i_images)  = buf*tangent(:,i_images)
  end do

! write(use_unit,*) 'DEBUG: perpendicular forces:'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') force_perp

! write(use_unit,*) 'DEBUG: parallel forces: '
! write(use_unit,'("DEBUG: ", 3E18.8E4)') force_par
  
  ! calculate correction force (9) in book chapter
  ! calculate cos phi for Eqn (10):
  ! first image
  tauplus(:)  = coords(:,2)-coords(:,1)
  tauminus(:) = coords(:,1)-start_coords(:) 
  buf         = sqrt(ddot(3*n_atoms,tauplus,1,tauplus,1)*ddot(3*n_atoms,tauminus,1,tauminus,1))
  cos_phi(1)  = ddot(3*n_atoms,tauplus,1,tauminus,1)/buf
  ! intermediate images
  do i_images = 2, n_images-1
     tauplus(:)        = coords(:,i_images+1)-coords(:,i_images)
     tauminus(:)       = coords(:,i_images)-coords(:,i_images-1)
     buf               = sqrt(ddot(3*n_atoms,tauplus,1,tauplus,1)*ddot(3*n_atoms,tauminus,1,tauminus,1))
     cos_phi(i_images) = ddot(3*n_atoms,tauplus,1,tauminus,1)/buf
  end do
  ! last image
  tauplus(:)        = end_coords(:) - coords(:,n_images)
  tauminus(:)       = coords(:,n_images) - coords(:,n_images-1)
  buf               = sqrt(ddot(3*n_atoms,tauplus,1,tauplus,1)*ddot(3*n_atoms,tauminus,1,tauminus,1))
  cos_phi(n_images) = ddot(3*n_atoms,tauplus,1,tauminus,1)/buf
  ! phi-dependent factor f(phi):
  pi = 2d0*dasin(1d0)
  f_phi(:) = (1d0+cos(pi*cos_phi(:)))/2d0

! write(use_unit,*) 'DEBUG: f(phi)'
! write(use_unit,'("DEBUG: ",E18.8E4)') f_phi

  ! calculate the correction force, for each image separately:
  do i_images = 1, n_images
     buf = ddot(3*n_atoms,spring_forces(:,i_images),1,tangent(:,i_images),1)
     force_corr(:,i_images) = f_phi(i_images)*(spring_forces(:,i_images)-buf*tangent(:,i_images))
  end do

! write(use_unit,*) 'DEBUG: correction forces'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') force_corr

  ! sum up total forces
  total_forces(:,:) = force_perp(:,:)! + force_par(:,:) + force_corr(:,:)
  call clean_translational_forces(n_atoms,n_images,total_forces)

! write(use_unit,*) 'DEBUG: total forces'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') total_forces

! write(use_unit,*) 'DEBUG: end of iteration, NEXT!!!'
end subroutine update_force_NEB

!---------------------------------------------------------------------------------------------------
! climbing image nudged elastic band routine according to the prescription by 
! Henkelman, Uberuaga, and Jonsson; JCP v113, p9901
subroutine update_force_CINEB(n_atoms, n_images, start_coords, coords, end_coords, &
     start_forces, total_forces, end_forces, &
     start_energy, energies, end_energy, spring_constants)
  implicit none
  integer n_atoms, n_images
  real*8 :: start_energy, end_energy, pi
  real*8, dimension(3*n_atoms) :: start_coords, end_coords, start_forces, end_forces, tauplus, tauminus
  real*8, dimension(n_images) :: energies, spring_constants
  real*8, dimension(3*n_atoms,n_images) :: coords, total_forces
  ! working variables
  real*8, dimension(3*n_atoms,n_images) :: spring_forces, correction_forces, tangent, force_perp, force_par, force_corr, f_CI
  real*8, dimension(n_images) :: tangent_magnitude, cos_phi, f_phi
  integer :: i_coords, i_images, i_atoms, i_maximum 
  real*8 :: V_first, V_second, ddot, buf

  ! regular elastic band method:
  ! calculate tangent vectors, according to the 'improved tangent vector' prescription 
  ! by Henkelmann and Jonsson, JCP v113 p9978, eqns (8) and eqns (10)
  tangent  = 0d0
  ! image 1
  if ((start_energy.lt.energies(1)).and.(energies(1).lt.energies(2))) then
     ! V_start < V_1 < V_2; use tau+ = R_2 - R_1
     tangent(:,1) = coords(:,2) - coords(:,1)
  else if ((start_energy.gt.energies(1)).and.(energies(1).gt.energies(2))) then
     ! V_start> V_1 > V_2; use tau- = R_1 - R_start
     tangent(:,1) = coords(:,1) - start_coords(:)
  else
     ! use energy-difference weighed average of the two vectors in the previous cases. 
     tauplus(:)  = coords(:,2) - coords(:,1)
     tauminus(:) = coords(:,1) - start_coords(:)
     if (energies(2).gt.start_energy) then
        V_first  = max(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
        V_second = min(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
     else
        V_first  = min(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
        V_second = max(abs(start_energy-energies(1)),abs(energies(2)-energies(1)))
     end if
     tangent(:,1) = V_first*tauplus(:) + V_second*tauminus(:)
  end if
  ! images 2 ... n-1
  do i_images = 2, n_images-1
     if ((energies(i_images-1).lt.energies(i_images)).and.(energies(i_images).lt.energies(i_images+1))) then
        tangent(:,i_images) = coords(:,i_images+1)-coords(:,i_images)
     else if ((energies(i_images-1).gt.energies(i_images)).and.(energies(i_images).gt.energies(i_images+1))) then
        tangent(:,i_images) = coords(:,i_images)-coords(:,i_images-1)
     else
        tauplus(:)  = coords(:,i_images+1)-coords(:,i_images  )
        tauminus(:) = coords(:,i_images  )-coords(:,i_images-1)
        if (energies(i_images+1).gt.energies(i_images-1)) then
           V_first  = max(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
           V_second = min(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
        else
           V_first  = min(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
           V_second = max(abs(energies(i_images+1)-energies(i_images)),abs(energies(i_images)-energies(i_images-1)))
        end if
!        write(use_unit,*) 'DEBUG: ', i_images, V_first, V_second
        tangent(:,i_images) = V_first*tauplus(:) + V_second*tauminus(:)
!        write(use_unit,*) 'DEBUG: intermediate tangent vector!'
!       write(use_unit,'("DEBUG intermediate: ",3E18.8E4)')tangent(:,i_images)

     end if
  end do
  ! last image
  if ((energies(n_images-1).lt.energies(n_images)).and.(energies(n_images).lt.end_energy)) then
     tangent(:,n_images) = end_coords(:) - coords(:,n_images)
  else if ((energies(n_images-1).gt.energies(n_images)).and.(energies(n_images).gt.end_energy)) then
     tangent(:,n_images) = coords(:,n_images) - coords(:,n_images-1)
  else
     tauplus (:) = end_coords(:) - coords(:,n_images)
     tauminus(:) = coords(:,n_images) - coords(:,n_images-1)
     if (end_energy.gt.energies(n_images-1)) then
        V_first  = max(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))
        V_second = min(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))
     else
        V_first  = min(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))
        V_second = max(abs(end_energy-energies(n_images)),abs(energies(n_images)-energies(n_images-1)))        
     end if
     tangent(:,i_images) = V_first*tauplus(:) + V_second*tauminus(:)
  end if
  ! normalize tangents:
  do i_images = 1, n_images
     tangent_magnitude(i_images) = sqrt(ddot(3*n_atoms,tangent(:,i_images),1,tangent(:,i_images),1))
     tangent(:,i_images) = tangent(:,i_images)/tangent_magnitude(i_images)
  end do

! write(use_unit,*) 'DEBUG: Tangent vectors:'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') tangent 

  ! calculate spring forces
  ! first image
  spring_forces(:,1) = spring_constants(1)*(start_coords(:)+coords(:,2)-2d0*coords(:,1))
  ! intermediate images
  do i_images = 2, n_images-1
     spring_forces(:,i_images) = spring_constants(i_images)*(coords(:,i_images+1)+coords(:,i_images-1)-2d0*coords(:,i_images))
  end do
  ! last image
  spring_forces(:,n_images) = spring_constants(n_images)*(end_coords(:)+coords(:,n_images-1)-2d0*coords(:,n_images))

! write(use_unit,*) 'DEBUG: spring forces:'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') spring_forces
  ! calculate (8) in book chapter
  do i_images = 1, n_images
     ! the perpendicular force, first term
     buf                    = ddot(3*n_atoms,total_forces(:,i_images),1,tangent(:,i_images),1)
     force_perp(:,i_images) = total_forces(:,i_images) - buf*tangent(:,i_images)
     ! the parallel force, second term: based on spring forces
     buf                    = ddot(3*n_atoms,spring_forces(:,i_images),1,tangent(:,i_images),1)
     force_par(:,i_images)  = buf*tangent(:,i_images)
  end do

! write(use_unit,*) 'DEBUG: perpendicular forces:'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') force_perp

! write(use_unit,*) 'DEBUG: parallel forces: '
! write(use_unit,'("DEBUG: ", 3E18.8E4)') force_par
  
  ! calculate correction force (9) in book chapter
  ! calculate cos phi for Eqn (10):
  ! first image
  tauplus(:)  = coords(:,2)-coords(:,1)
  tauminus(:) = coords(:,1)-start_coords(:) 
  buf         = sqrt(ddot(3*n_atoms,tauplus,1,tauplus,1)*ddot(3*n_atoms,tauminus,1,tauminus,1))
  cos_phi(1)  = ddot(3*n_atoms,tauplus,1,tauminus,1)/buf
  ! intermediate images
  do i_images = 2, n_images-1
     tauplus(:)        = coords(:,i_images+1)-coords(:,i_images)
     tauminus(:)       = coords(:,i_images)-coords(:,i_images-1)
     buf               = sqrt(ddot(3*n_atoms,tauplus,1,tauplus,1)*ddot(3*n_atoms,tauminus,1,tauminus,1))
     cos_phi(i_images) = ddot(3*n_atoms,tauplus,1,tauminus,1)/buf
  end do
  ! last image
  tauplus(:)        = end_coords(:) - coords(:,n_images)
  tauminus(:)       = coords(:,n_images) - coords(:,n_images-1)
  buf               = sqrt(ddot(3*n_atoms,tauplus,1,tauplus,1)*ddot(3*n_atoms,tauminus,1,tauminus,1))
  cos_phi(n_images) = ddot(3*n_atoms,tauplus,1,tauminus,1)/buf
  ! phi-dependent factor f(phi):
  pi = 2d0*dasin(1d0)
  f_phi(:) = (1d0+cos(pi*cos_phi(:)))/2d0

! write(use_unit,*) 'DEBUG: f(phi)'
! write(use_unit,'("DEBUG: ",E18.8E4)') f_phi

  ! calculate the correction force, for each image separately:
  do i_images = 1, n_images
     buf = ddot(3*n_atoms,spring_forces(:,i_images),1,tangent(:,i_images),1)
     force_corr(:,i_images) = f_phi(i_images)*(spring_forces(:,i_images)-buf*tangent(:,i_images))
  end do

! write(use_unit,*) 'DEBUG: correction forces'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') force_corr

  ! now, work on the climbing image force
  ! determine the single image with the maximal energy: it is supposed to be converged to the TS
  i_maximum = 1
  do i_images = 2, n_images
     if (energies(i_images).gt.energies(i_maximum)) i_maximum = i_images
  end do
  ! delete all other force components on this particular image
  f_CI(:,:)               = 0d0
  force_perp(:,i_maximum) = 0d0
  force_par (:,i_maximum) = 0d0
  force_corr(:,i_maximum) = 0d0
  ! F_total,imax . tangent_imax
  buf = ddot(3*n_atoms,total_forces(:,i_maximum),1,tangent(:,i_maximum),1)
  ! Climbing image force
  f_CI(:,i_maximum) = total_forces(:,i_maximum) - 2d0*buf*tangent(:,i_maximum)
  ! sum up remaining forces
  total_forces(:,:) = force_perp(:,:) + force_par(:,:) + force_corr(:,:) + f_CI(:,:)
  call clean_translational_forces(n_atoms,n_images,total_forces)

! write(use_unit,*) 'DEBUG: total forces'
! write(use_unit,'("DEBUG: ", 3E18.8E4)') total_forces

! write(use_unit,*) 'DEBUG: end of iteration, NEXT!!!'
end subroutine update_force_CINEB
