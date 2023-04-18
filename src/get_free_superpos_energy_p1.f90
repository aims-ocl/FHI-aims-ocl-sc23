!****s* FHI-aims/get_free_superpos_energy_p1
!  NAME
!   get_free_superpos_energy_p1
!  SYNOPSIS

subroutine  get_free_superpos_energy_p1 &
     ( partition_tab, rho, rho_gradient, kinetic_density, pot_ion_embed, &
     hartree_energy_free, en_elec_free, en_xc, en_pot_xc,  &
     en_ion_ion, en_ion_embed, en_density_embed, en_vdw, en_ll_vdw, en_ll_vdw_err, &
     en_lda_c, en_pbe_c, ionic_forces, ext_charge_forces &
     )

!  PURPOSE
!  Subroutine get_free_superpos_energy.f
!
!  calculates:
!  * free atom contribution of the hartree energy 
!  * XC energies for the superposition of free atoms
!  * ion-ion interaction energy
!
!  USES

  use dimensions
  use runtime_choices
  use geometry
  use species_data
  use grids
  use spline
  use free_atoms
  use localorb_io
  use mpi_utilities
  use synchronize_mpi
  use constants
  use vdw_correction
  use ll_vdwdf
  use pseudodata
  use energy_density
  use pbc_lists
  use physics,  ONLY: elec_free_atom 
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: partition_tab
  real*8, dimension(n_spin, n_full_points) :: rho
  real*8, dimension(3, n_spin, n_full_points) :: rho_gradient
  real*8, dimension(n_spin, n_full_points) :: kinetic_density
  real*8, dimension(n_full_points) :: pot_ion_embed

  real*8 :: hartree_energy_free
  real*8 :: en_xc
  real*8 :: en_pot_xc
  real*8 :: en_ion_ion
  real*8 :: en_elec_free
  real*8 :: en_ion_embed
  real*8 :: en_density_embed
  real*8 :: en_vdw
  real*8 :: en_ll_vdw, en_ll_vdw_err
  real*8 :: en_lda_c, en_pbe_c

  real*8, dimension(3, n_atoms) :: ionic_forces
  real*8, dimension(3, n_multipoles) :: ext_charge_forces

! INPUTS
! o partition_tab -- values of partition function
! o rho -- electron density
! o rho_gradient -- gradient of electron density  (this should only ever be referenced if use_density_gradient is true.)
! o pot_ion_embed -- embedded ionic potential (this should only ever be referenced if use_embedding_potential or 
!                    use_embedding_pp is true.)
! o kinetic_density -- kinetic-energy density (only referenced for use_meta_gga)
!
! OUTPUT
! o hartree_energy_free -- free atoms Hartree potential
! o en_elec_free -- energy of free atoms electron density
! o en_xc -- exchange correlation energy of free atoms
! o en_pot_xc -- exchange correlation potential energy of free atoms
! o en_ion_ion -- ion-ion energy
! o en_ion_embed -- energy of embedded ions
! o en_density_embed -- energy of enbedded charge density
! o en_vdw -- energy of Van der Waals interaction
! o en_ll_vdw  -- non_local correlation energy from Langreth-Lundquist van der Waals density functional
! o en_ll_vdw_err  -- error of en_ll_vdw
! o en_lda_c -- pwlda correlation energy (will be used for ll_vdwdf energy calculation)
! o en_pbe_c -- pbe correlation energy (will be used for ll_vdwdf energy calculation)
! o ionic_forces -- ionic forces of free atoms
! o ext_charge_forces -- forces of external charge density
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



  !  local variables

  real*8, dimension(3) :: coord_current
  real*8 :: distance
  real*8 :: distance_squared
  real*8, dimension(3) :: direction
  real*8 :: i_r
  integer, dimension(:), allocatable :: center_index
  real*8, dimension(:), allocatable :: temp_hartree_sum
  real*8, dimension(:), allocatable :: temp_free_pot_es
  real*8, dimension(:), allocatable :: temp_free_rho
  real*8, dimension(:), allocatable :: temp_nucl
  real*8, dimension(:), allocatable :: atomic_hartree_energy
  !  real*8 :: temp_hartree
  real*8, dimension(n_max_grid) :: free_hartree_electronic

  real*8, dimension(n_spin) :: local_rho
  real*8, dimension(3,n_spin) :: local_rho_gradient
  real*8, dimension(n_spin) :: local_kinetic_density
  real*8, dimension(n_spin) :: local_xc_derivs
  real*8, dimension(3,n_spin) :: local_xc_gradient_deriv
  real*8, dimension(n_spin) :: local_xc_tau_deriv
  real*8 :: en_density_xc
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c
  real*8 :: temp_hartree_ion_ion

  !  real*8 :: squared_rho_gradient

  real*8 :: pot_ion_embed_nuc
  real*8, dimension(3) :: embedding_force, pp_embedding_field
  real*8 :: full_rho

  character*100 :: info_str

  integer :: n_compute_atoms


  !  counters

  integer :: i_species
  integer :: i_coords
  integer :: i_atom
  integer :: i_atom_2
  !  integer :: i_atom_3
  !  integer :: i_angular
  ! integer :: i_radial
  integer :: i_coord
  integer :: i_grid
  integer :: i_full_points
  integer :: i_spin
  integer :: i_center, i_center_2, i_center_L

  integer :: i_index
  !  integer :: current_atom, current_radial, current_angular
  integer :: i_my_batch

  integer :: i_compute_atom
  !  integer :: i_compute_atom_2
  integer :: i_pp_atoms, i_multipole
  !     functions

  real*8 :: get_distance
  real*8 :: int_log_mesh, rmin
  real*8 :: ddot
  real*8 :: ion_pseudo_pot

!  real*8, dimension(:,:), allocatable  :: rho_inc_partialcore
!  real*8, dimension(:,:,:), allocatable  :: rho_gradient_inc_partialcore
  !  begin work


  write (info_str,'(2X,A,A)')   &
       "Calculating total energy contributions",  &
       " from superposition of free atom densities."

  call localorb_info(info_str,use_unit,'(A)',OL_norm)
  
  if (flag_out_locpot_atom) then
     elec_free_atom = 0.0d0
  !   elec_atomic = 0.0d0 
  endif

  ! CC: Allocate arrays for energy density and set them to zero
  if (flag_energy_density) then
    if(.not.allocated(ed_hartree_energy_free))        allocate(ed_hartree_energy_free       (n_full_points))
    if(.not.allocated(ed_hartree_energy_free_nuclei)) allocate(ed_hartree_energy_free_nuclei(n_atoms))
    ed_hartree_energy_free(:)        = 0.0d0
    ed_hartree_energy_free_nuclei(:) = 0.0d0
  end if
 
  !     Preparation: Calculate free-atom Hartree energies for each species on logarithmic grids

  allocate(atomic_hartree_energy(n_species),stat=i_atom)
  call check_allocation(i_atom, 'atomic_hartree_energy         ')
  atomic_hartree_energy(:)=0d0



 ! if(use_embedding_pp.and.use_nonlinear_core) then
 !     allocate(rho_inc_partialcore(n_spin,n_full_points))
 !!     allocate(rho_gradient_inc_partialcore(3,n_spin,n_full_points))
 !     do i_spin = 1,n_spin
 !        rho_inc_partialcore(i_spin,:) = rho(i_spin,:) + partial_core_rho(:)
 !!        rho_gradient_inc_partialcore(:,i_spin,:) = rho_gradient(:,i_spin,:)
 !!        if(use_density_gradient) then
 !!           rho_gradient_inc_partialcore(:,i_spin,:) = &
 !!              rho_gradient_inc_partialcore(:,i_spin,:) + partial_core_rho_grad(:,:)
 !!        endif
 !     enddo
 ! endif




  do i_species = 1, n_species, 1
! DB: 0912812: actually we don't have to cycle here, since free_pot_es and
!              free_rho are set zero for pseudoized species, but we can.
   if(species_pseudoized(i_species).or.no_basis(i_species)) cycle

     !        calculate atomic hartree-energies on logarithmic grid
     !        factor 4*Pi*r^2 !!!
     !        ionic part is subtracted ( - ( - z/r) -> + z/r )

     if ( (n_periodic .eq. 0) .and. (.not.(force_new_functional)) ) then

        do i_grid = 1, n_grid(i_species), 1

           free_hartree_electronic(i_grid) =  &
                (free_pot_es(i_grid, i_species) *  &
                r_grid(i_grid, i_species) +  &
                species_z(i_species) ) *  &
                r_grid(i_grid, i_species)

        end do

     else

        do i_grid = 1, n_grid(i_species), 1

           free_hartree_electronic(i_grid) =  &
                (free_pot_es(i_grid, i_species) *  &
                r_grid(i_grid, i_species)) *  &
                r_grid(i_grid, i_species)

        end do
     end if

     !        evaluate diagonal part of free atomic hartree energy on
     !        logarithmic grid

     atomic_hartree_energy(i_species) =   &
          int_log_mesh(free_rho(1, i_species),   &
          free_hartree_electronic, n_grid(i_species),  &
          r_grid(1, i_species) )

  enddo

  !     begin energy summations for actual structure

  allocate(center_index(n_centers_basis_integrals),stat=i_atom)
  call check_allocation(i_atom, 'center_index                  ')

  allocate(temp_hartree_sum(n_centers_basis_integrals),stat=i_atom)
  call check_allocation(i_atom, 'temp_hartree_sum              ')

  allocate(temp_free_pot_es(n_centers_basis_integrals),stat=i_atom)
  call check_allocation(i_atom, 'temp_free_pot_es              ')

  allocate(temp_free_rho(n_centers_basis_integrals),stat=i_atom)
  call check_allocation(i_atom, 'temp_free_rho                 ')

  allocate(temp_nucl(n_centers_basis_integrals),stat=i_atom)
  call check_allocation(i_atom, 'temp_nucl                     ')


  hartree_energy_free = 0.d0
  en_ion_ion = 0.d0
  en_ion_embed = 0.d0
  en_density_embed = 0.d0
  ! en_vdw = 0.d0 !Here is commented because with the SC vdW scheme the initialization to zero is already 
  ! done in integrate_hamiltonian. If we keep the initialization to zero here, we would get a vdW energy 
  ! set equal to zero in the first SC cycle. Probably this doesn't affect the final result at convergency, 
  ! but still it's preferable to have the correct starting point.
  en_elec_free = 0.d0

  if (use_forces) then
     ionic_forces = 0.d0
     if (use_qmmm) then
        ext_charge_forces = 0.0
     end if
  end if

  en_xc      = 0.d0
  en_pot_xc  = 0.d0


  ! WARNING : some cycles will be split for real atoms an ghost atoms for optimization reasons,
  !           the first one would run over occupied centers - AS


  do i_atom = 1, n_occ_atoms, 1

     if (myid.eq.task_list(i_atom)) then

        !test
        !       write(use_unit,'(2X,A,I4)') "Atom ", i_atom
        !test end

        !         first, increment ion-ion interaction energy

        do i_center_L =  i_atom + 1, n_occ_atoms, 1
           i_center = new_centers_basis_integrals(i_center_L)

           distance = get_distance(coords(1, i_atom),   &
                coords_center(1, i_center))

           en_ion_ion = en_ion_ion +  &
                species_z(species(i_atom)) * species_z(species_center(i_center))  &
                / distance

        end do

        ! VB: FIXME: WHEN WE GET TO CONSISTENT FORCES, THIS PART MUST TAKE
        !     force_new_functional INTO ACCOUNT.

        ! ionic force contributions if needed
        if (use_forces .and. (n_periodic == 0 .and. (.not. force_new_functional))) then
           do i_atom_2 = 1, n_occ_atoms, 1
              if (i_atom .ne. i_atom_2) then

                 distance_squared = 0.d0
                 do i_coords = 1, 3, 1
                    direction(i_coords) =   &
                         coords(i_coords, i_atom)-coords(i_coords, i_atom_2)  
                    distance_squared = distance_squared +   &
                         direction(i_coords) * direction(i_coords)
                 end do

                 do i_coords = 1, 3, 1
                    ionic_forces(i_coords, i_atom) =   &
                         ionic_forces(i_coords, i_atom) +   &
                         species_z(species(i_atom_2)) *   &
                         direction(i_coords) / (distance_squared ** 1.5d0)
                 end do
              end if
           enddo
           ionic_forces(:, i_atom) = ionic_forces(:, i_atom) * species_z(species(i_atom))
        end if
        !           end if

        ! Forces on external charges coming from the "quantum" atoms
        if (use_forces .and. (n_periodic == 0) .and. use_qmmm) then
           do i_atom_2 = 1, n_multipoles, 1
              if (multipole_order(i_atom_2) == 0) then

                 distance_squared = 0.d0
                 do i_coords = 1, 3, 1
                    direction(i_coords) =   &
                         coords(i_coords, i_atom) - multipole_coords(i_coords,i_atom_2)
                    distance_squared = distance_squared +   &
                         direction(i_coords) * direction(i_coords)
                 end do

                 do i_coords = 1, 3, 1
                    ext_charge_forces(i_coords, i_atom_2) =   &
                         ext_charge_forces(i_coords, i_atom_2) -   & !!! (AT - Minus !) 
                         multipole_charge(i_atom_2) * species_z(species(i_atom)) * &
                         direction(i_coords) / (distance_squared ** 1.5d0)

                 end do
                          !       write(use_unit,'(2X,A,I4,1X,3(E30.15,1X))') "|",i_atom, ext_charge_forces(:,i_atom_2) 
              end if
           enddo
        end if

        !         Nadia
        !         second, if (use_embedding_potential) get ion-embedded charge interaction energy        
        if (use_embedding_potential) then
           call embedding_potential(coords(1,i_atom),pot_ion_embed_nuc,embedding_force)
           ! needs negative sign because the embedding potential already
           ! contains the negative charge of the electron.
           en_ion_embed = en_ion_embed  &
                - species_z(species(i_atom)) * pot_ion_embed_nuc

           if (use_forces .and. n_periodic == 0) then
              do i_coords = 1, 3, 1
                 ionic_forces(i_coords, i_atom) = ionic_forces(i_coords, i_atom)  -   &  !!! Minus 
                      embedding_force(i_coords) * species_z(species(i_atom))
              end do
                          !    write (info_str,'(2X,A)')   &
                          !   "Embedding forces on atoms:"
                          !    call localorb_info(info_str,use_unit,'(A)')
                          !    write(use_unit,'(2X,A,I4,1X,3(E30.15,1X))') "|",i_atom,ionic_forces(:, i_atom)

           end if
        end if

        ! BL: Add homogeneous field contribution always 
        ! if a homogeneous field exists
        ! Paula: the forces of homogeneous fields.
        ! BL: The sign convention puzzles me: F = E * q
        ! why is there a minus, is this correct? 
        if (use_forces .and. allocated(homogeneous_field)) then
           do i_coords = 1, 3, 1
              ionic_forces(i_coords, i_atom) = &
                 ionic_forces(i_coords, i_atom) - & 
                 homogeneous_field(i_coords) * species_z(species(i_atom))
           end do
        end if

        !DB: adding the interaction of the nuclei with the local potential part of the pseudopots
        if (use_embedding_pp) then

           call ion_pseudopot_potential(coords(:,i_atom), ion_pseudo_pot,0, pp_embedding_field)

           en_ion_embed = en_ion_embed &
             - species_z(species(i_atom)) * ion_pseudo_pot

           if (use_forces .and. n_periodic == 0) then
              do i_coords = 1, 3, 1
                 ionic_forces(i_coords, i_atom) = ionic_forces(i_coords, i_atom)  -   &  !!! Minus 
                      pp_embedding_field(i_coords) * species_z(species(i_atom))
              end do
                      !        write (info_str,'(2X,A)')   &
                      !       "Embedding forces on atoms:"
                      !        call localorb_info(info_str,use_unit,'(A)')
                      !        write(use_unit,'(2X,A,I4,1X,3(E30.15,1X))') "|",i_atom,pp_embedding_field(:) 
           end if



        end if

        !         next, add on-site Hartree energy of free atoms
        hartree_energy_free = hartree_energy_free +  &
             atomic_hartree_energy(species(i_atom))
        
        ! CC: Save nuclear contributions to hartree_energy_free
        ! FIXME: Actually this is the energy of a free atom in its own field
        ! Construction of the associated energy density field might be more correct
        ! However, this should not matter for heat fluxes, since this field moves with
        ! the nuclear coordinate
        if (flag_energy_density) then
          ed_hartree_energy_free_nuclei(i_atom) = atomic_hartree_energy(species(i_atom))
        end if

     end if


  end do ! end loop over atoms



!FIXME: DB 16/02: in analogy to multipol embedding we leave out the interaction between different pseudocores
!  if (use_embedding_pp) then
!      do i_pp_atoms = 1, n_pp_atoms, 1 
!  
!      call ion_pseudopot_potential_v1(pp_coords(:,i_pp_atoms), ion_pseudo_pot, i_pp_atoms)
!  
!           ! we also need the minus here .. we checked that!
!
!      en_ion_embed = en_ion_embed &
!          - pp_charge(pp_species(i_pp_atoms)) * ion_pseudo_pot 
!
!      enddo
!
!      if (n_multipoles.gt.0) then
!         do i_multipole = 1,n_multipoles
!
!           call ion_pseudopot_potential(multipole_coords(:,i_multipole), ion_pseudo_pot,0)
!
!           en_ion_embed = en_ion_embed &
!             - multipole_charge(i_multipole) * ion_pseudo_pot
!
!         end do
!      endif
!
!  end if


  ! DEBUG 1
  !  write(use_unit,*) "hartree_energy_free after sum of atomic_hartree_energy:", hartree_energy_free 
  ! END DEBUG
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  Here starts the ghost atoms cycle
  if (n_occ_atoms .lt. n_atoms) then
     do i_atom = n_occ_atoms + 1, n_atoms, 1

        if (myid.eq.task_list(i_atom)) then
           if (species_pseudoized(species(i_atom))) cycle
           do i_center_L =  i_atom + 1, n_atoms, 1

              i_center = new_centers_basis_integrals(i_center_L)

              distance = get_distance(coords(1, i_atom),   &
                   coords_center(1, i_center))
           end do

           ! VB: FIXME: WHEN WE GET TO CONSISTENT FORCES, THIS PART MUST TAKE
           !     force_new_functional INTO ACCOUNT.
           ! ionic force contributions if needed
           if (use_forces .and. (n_periodic == 0 .and.(.not.force_new_functional ))) then
              ionic_forces(:, i_atom) = ionic_forces(:, i_atom) * species_z(species(i_atom))
           end if
        end if

     end do ! end loop over ghost atoms
  endif


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  temp_nucl = 0.d0

  temp_free_pot_es = 0.d0

  temp_free_rho = 0.d0

  !       now obtain "off-diagonal" Hartree terms and the XC energy by 3D integration

  i_full_points = 0

  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1

           !           execute only if partition_tab.gt.0, i.e. if the integration point
           !           makes sense
           ! Note that partition_tab is zero for instance just when we are right on another
           ! atom!
           if (partition_tab(i_full_points).gt.0.d0) then

              !             get current coordinates
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

              !             tabulate non-zero free-atom density values and, for each atom,
              !             the sum of free-atom hartree potentials of all _other_ atoms 
              !             at current integration point

              n_compute_atoms = 0

              do i_center_L = 1, n_occ_centers_basis_integr, 1
                 i_center = new_centers_basis_integrals(i_center_L)

                 distance = get_distance(coord_current, coords_center(1, i_center))
!NEC_CB Move into if-branch below
!NEC_CB          i_r = invert_log_grid_p2 ( distance, species_center(i_center))

                 ! only in non-periodic case
                 if ( (n_periodic.eq.0) .and. (.not.force_new_functional) ) then
                    temp_nucl(i_center_L) = species_z(species_center(i_center))/distance
                 end if

                 if (distance.lt.multipole_radius_free(species_center(i_center))) then
                    ! This atom is relevant in sum over free-atom densities.

                    n_compute_atoms = n_compute_atoms+1

                    center_index(n_compute_atoms) = i_center_L

                    i_r = invert_log_grid_p2 ( distance, species_center(i_center))

                    temp_free_rho(n_compute_atoms) = pi4_inv  &
                         * val_spline  &
                         (i_r, renormalized_free_rho_spl(1,1,species_center(i_center)),  &
                         n_grid(species_center(i_center)))  

                    temp_free_pot_es(n_compute_atoms) = &
                         val_spline &
                         ( i_r, free_pot_es_spl(1,1, &
                         species_center(i_center)), &
                         n_grid(species_center(i_center))) 

                 end if

              end do

              ! For each non-zero atom, sum up the _potential_ contribution from all _other_ atoms
              temp_hartree_sum(1:n_compute_atoms) = 0.d0
              do i_compute_atom = 1, n_compute_atoms, 1
                 !                do i_compute_atom_2 = 1, n_compute_atoms, 1
                 ! if (i_compute_atom_2.ne.i_compute_atom) then
                 !                    temp_hartree_sum(i_compute_atom) = & 
                 !                    temp_hartree_sum(i_compute_atom) + temp_free_pot_es(i_compute_atom_2)
                 ! end if
                 !                enddo

                 temp_hartree_sum(i_compute_atom) = & 
                      temp_hartree_sum(i_compute_atom) + sum(temp_free_pot_es(1:i_compute_atom-1))

                 temp_hartree_sum(i_compute_atom) = & 
                      temp_hartree_sum(i_compute_atom) + sum(temp_free_pot_es(i_compute_atom+1:n_compute_atoms))

              enddo

              ! For the non-periodic energy functional, must for now add the contributions from
              ! the nuclei. FIXME: THIS SHOULD GO AWAY WHEN THE FUNCTIONAL IS SWITCHED FOR THE CLUSTER CASE!!
              if ( (n_periodic.eq.0) .and. (.not.force_new_functional) ) then
                 do i_compute_atom = 1, n_compute_atoms, 1

                    temp_hartree_sum(i_compute_atom) = temp_hartree_sum(i_compute_atom) + &
                         sum(temp_nucl(1:(center_index(i_compute_atom)-1)))

                    temp_hartree_sum(i_compute_atom) = temp_hartree_sum(i_compute_atom) + &
                         sum(temp_nucl((center_index(i_compute_atom)+1):n_occ_centers_basis_integr))

                    !                  do i_center_L = 1, n_occ_centers_basis_integr, 1
                    !                    if (i_center_L.ne.center_index(i_compute_atom)) then
                    !                      temp_hartree_sum(i_compute_atom) = & 
                    !                      temp_hartree_sum(i_compute_atom) + temp_nucl(i_center_L)
                    !                    end if
                    !                  enddo

                 enddo
              end if

              temp_hartree_sum(1:n_compute_atoms) = &
                   temp_hartree_sum(1:n_compute_atoms) * partition_tab(i_full_points)

              hartree_energy_free = hartree_energy_free +  &
                   ddot(n_compute_atoms, temp_free_rho,1,temp_hartree_sum,1)

              ! CC: Save hartree_energy_free at each point of the grid
              if (flag_energy_density) then
                   ed_hartree_energy_free(i_full_points) = &
                    ddot(n_compute_atoms, temp_free_rho,1,temp_hartree_sum,1) &
                    / partition_tab(i_full_points)
              end if

              !             increment XC energies
              do i_spin = 1, n_spin,1
                 if(use_embedding_pp.and.use_nonlinear_core) then
                    local_rho(i_spin) = rho(i_spin,i_full_points) +  partial_core_rho(i_full_points)
                 else
                    local_rho(i_spin) = rho(i_spin,i_full_points)
                 endif
              end do

              if (use_gga) then
                 do i_spin = 1, n_spin,1 
                    do i_coord = 1,3,1
                       ! Put this statement here to remove repitive coding below
                       if(use_embedding_pp.and.use_nonlinear_core) then
                          local_rho_gradient(i_coord,i_spin) =   &
                                 rho_gradient(i_coord,i_spin,i_full_points)+&
                                 partial_core_rho_grad(i_coord,i_full_points)
                       else
                          local_rho_gradient(i_coord,i_spin) =   &
                               rho_gradient(i_coord,i_spin,i_full_points)
                       endif
                    enddo
                 enddo
              end if

              if (use_meta_gga) then
                 !if (.not. first_integration) first_integration = .true. 
                 !if (myid.eq.0) write(use_unit,*) kinetic_density(:,:)
                 ! kinetic_density = 0.0

                 do i_spin = 1, n_spin,1
                    local_kinetic_density(i_spin) =   &
                          kinetic_density(i_spin,i_full_points)
                 enddo
              endif

              call evaluate_xc( &
                  local_rho, &
                  local_rho_gradient, local_kinetic_density, &
                  en_density_xc, &
                  en_density_x, en_density_c, &
                  local_xc_derivs, local_xc_gradient_deriv, &
                  local_xc_tau_deriv, &
                  (use_hartree_fock .or. use_meta_gga) .and. first_integration, &
                  coord_current(:))

              !             calculate exchange correlation energy and mean value of 
              !             exchange correlation potential
              !
              ! Since we already have the pieces, add terms of XC energy here.

              ! if(use_embedding_pp.and.use_nonlinear_core) then
              !
              !     call evaluate_xc_energy_shell  &
              !          ( 1, partition_tab(i_full_points), en_density_xc, local_xc_derivs,  &
              !          ! local_xc_gradient_deriv, local_xc_tau_deriv, rho_inc_partialcore(:,i_full_points), local_rho_gradient,  &
              !          local_xc_gradient_deriv, local_xc_tau_deriv, local_rho, local_rho_gradient,  &
              !          local_kinetic_density, en_xc, en_pot_xc  )
              !
              ! else

                   call evaluate_xc_energy_shell  &
                        ( 1, partition_tab(i_full_points), en_density_xc, local_xc_derivs,  &
                        ! local_xc_gradient_deriv, local_xc_tau_deriv, rho(:,i_full_points), local_rho_gradient,  &
                        local_xc_gradient_deriv, local_xc_tau_deriv, local_rho, local_rho_gradient,  &
                        local_kinetic_density, en_xc, en_pot_xc  )

              ! endif

              !if (use_gga) then
              ! !                add partially integrated gradient pieces of v_xc to potential energy density 
              !   
              !   do i_spin = 1, n_spin, 1
              !      do i_coord = 1, 3, 1
              !
              !         en_pot_xc = en_pot_xc +  &
              !              2.d0 * local_xc_gradient_deriv(i_coord,i_spin) *  &
              !              local_rho_gradient(i_coord,i_spin) *  &
              !              partition_tab(i_full_points)
              !
              !      enddo
              !   enddo
              !end if

              ! Adding meta-gga evaluation. This should all be removed and replaced with
              ! a call to evaluate_xc_energy_shell! AJL
              !
              ! To be uncommented once I've got a superposition of the initial kinetic density
              ! if (use_meta_gga) then
              !  do i_spin = 1, n_spin, 1
              !    en_pot_xc = en_pot_xc + &
              !       local_xc_tau_deriv(i_spin) * &
              !       local_kinetic_density(i_spin) * &
              !       partition_tab(i_full_points)
              !  enddo
              ! endif

              if (use_embedding_potential) then
                 full_rho = 0.d0
                 do i_spin = 1,n_spin, 1
                    full_rho = full_rho + rho(i_spin, i_full_points)
                 enddo
                 en_density_embed = en_density_embed +   &
                      partition_tab(i_full_points) * full_rho * &
                      pot_ion_embed(i_full_points)

              end if

              if (use_embedding_pp) then
                 full_rho = 0.d0
                 do i_spin = 1,n_spin, 1
                    full_rho = full_rho + rho(i_spin, i_full_points)
                 enddo
                 en_density_embed = en_density_embed +   &
                      partition_tab(i_full_points) * full_rho * &
                      whole_local_pseudpot_on_intgrid(i_full_points)

              end if

              ! end check for zero partition tab.
           end if

           !     end loop over points in a batch
        end do

        !     end distribution of tasks
     ! end if

     !     end loop over grid batches


  end do

  ! en_elec_free is the term written in the last line of Eq. (64) of 
  ! Blum et al., Computer Physics Communications 180 (2009) 2175-2196
  ! when evaluated for the superposition density of spherical free atoms.

  en_elec_free = hartree_energy_free

  !-------------------------------------
  if ( (n_periodic .gt. 0) .or. force_new_functional ) then

     !   Ghost coords
     if (n_occ_atoms.lt.n_atoms) then
        do i_atom=n_occ_atoms+1, n_atoms

           do i_coord = 1, 3, 1
              coord_current(i_coord) = coords(i_coord, i_atom) 
           end do
           temp_free_pot_es(i_atom) = 0.d0

           temp_nucl(i_atom) = 0.d0

        enddo
     endif

     do i_atom = 1, n_occ_atoms, 1

        if (myid.eq.task_list(i_atom)) then

           rmin = r_grid_min(species(i_atom))

           !             calculate values of atomic hartree potential of
           !             all atoms at their respective rmin

           !             contribution from current atom

           i_r = &
                invert_log_grid &
                ( rmin, &
                r_grid_min(species(i_atom)), &
                r_grid_inc(species(i_atom)) &
                )  

           temp_free_pot_es(i_atom) = &
                val_spline &
                (i_r, free_pot_es_spl(1,1,species(i_atom)), &
                n_grid(species(i_atom)))

           temp_nucl(i_atom) = 0.d0
           temp_free_pot_es(i_atom) = free_pot_es_at_zero(species(i_atom))

           do i_coord = 1, 3, 1
              coord_current(i_coord) = coords(i_coord, i_atom)
           end do


           do i_atom_2 = 1,   n_occ_centers_basis_integr, 1
              if (i_atom_2 .ne. i_atom) then  

                 i_center_2 =  new_centers_basis_integrals( i_atom_2)

                 distance = &
                      get_distance(coord_current,  &
                      coords_center(1, i_center_2))

                 i_r =  &
                      invert_log_grid &
                      ( distance, &
                      r_grid_min(species_center(i_center_2)), &
                      r_grid_inc(species_center(i_center_2)) &
                      )  

                 temp_free_pot_es(i_atom_2) = &
                      val_spline &
                      (i_r, free_pot_es_spl(1,1, &
                      species_center(i_center_2)), &
                      n_grid(species_center(i_center_2))) 

                 temp_nucl(i_atom_2) = &
                      species_z(species_center(i_center_2))/ distance    

              end if
           end do


           !        calculate non-diagonal elements of free atom hartree energy

           temp_hartree_ion_ion = 0.d0

           do i_atom_2 = 1, n_occ_atoms

              ! VB: FIXME: This term, by way of adding the nuclear energy,
              !            simply recomputes the ion-ion contribution computed
              !            before. In principle, this is not necessary!

              ! orig version            temp_hartree_ion_ion = temp_hartree_ion_ion +  &
              !                 (temp_free_pot_es(i_atom_2)+ &
              !                 temp_nucl(i_atom_2)) * &
              !                 species_z(species(i_atom))

              ! new version - hartee_energy_free no longer includes the nuclear contribution


              temp_hartree_ion_ion = temp_hartree_ion_ion +  &
                   temp_free_pot_es(i_atom_2) * &
                   species_z(species(i_atom))
           end do

           do i_atom_2 = n_atoms+1, n_occ_centers_basis_integr, 1

              i_center_2 =  new_centers_basis_integrals( i_atom_2)

              ! temp_hartree_ion_ion is the term written in the last line of Eq. (66) of 
              ! Blum et al., Computer Physics Communications 180 (2009) 2175-2196
              ! where the counter i_atom is synonmous with the summation index "at" in the paper.
              temp_hartree_ion_ion = temp_hartree_ion_ion + &
                   temp_free_pot_es(i_atom_2) * &
                   species_z(species(i_atom)) 
           end do
           
           !TZ store the free atom contribution 
           !if (flag_out_locpot_atom) then 
           !elec_free_atom(i_atom) =  temp_hartree_ion_ion/species_z(species(i_atom))
           !endif 
           ! hartree_energy_free is the quantity given in Eq. (61) of the 
           ! FHI-aims CPC paper, Blum et al., Computer Physics Communications 180 (2009) 2175-2196
           ! when evaluated for the superposition density of spherical free atoms.
           ! 
           ! temp_hartree_ion_ion is the term written in the last line of Eq. (66) of 
           ! Blum et al., Computer Physics Communications 180 (2009) 2175-2196
           ! where the counter i_atom is synonmous with the summation index "at" in the paper.

	   hartree_energy_free = hartree_energy_free + temp_hartree_ion_ion
           ! CC: Save nuclear contributions to hartree_energy_free
           ! This is the electrostatic field of nucleus i_atom in the field of all other atoms
           if (flag_out_locpot_atom) then 
              elec_free_atom(i_atom) = free_pot_es_at_zero(species(i_atom))
           endif

           if (flag_energy_density) then
             ed_hartree_energy_free_nuclei(i_atom) = ed_hartree_energy_free_nuclei(i_atom) + temp_hartree_ion_ion
           end if
           

        end if
     end do

  end if

  if (use_vdw_correction) then
     ! call vdW correction
     write (info_str,'(2X,A)')   &
          "Calculating vdW correction."
     call localorb_info(info_str,use_unit,'(A)')

     call calc_vdw(en_vdw,ionic_forces)

     !  write (info_str,'(2X,A,F10.6,A,F7.3,A)')   &
     !       "vdW energy: ", en_vdw, " Ha     ", en_vdw*hartree, " eV"
     !  call localorb_info(info_str,use_unit,'(A)')
  end if

!  if (use_ll_vdwdf) then
!     ! call vdW correction
!     write (info_str,'(2X,A)')   &
!          "Calculating LL_vdwdf correction."
!     call localorb_info(info_str,use_unit,'(A)')
!
!     call calc_ll_vdw(en_ll_vdw,en_ll_vdw_err,en_lda_c,en_pbe_c)
!
!  end if

  call sync_get_free_superpos_energy  &
       (  &
       hartree_energy_free, en_xc, en_pot_xc, en_ion_ion,  &
       en_elec_free, &
       en_ion_embed, en_density_embed, en_vdw, en_ll_vdw, en_ll_vdw_err, &
       en_lda_c, en_pbe_c, ionic_forces  &
       )
  if (flag_out_locpot_atom) then
      call sync_vector(elec_free_atom, n_atoms)
      !call sync_vector(elec_atomic, n_atoms)
      !do i_atom = 1, n_atoms, 1 
      !   elec_free_atom(i_atom) = elec_free_atom(i_atom)/species_z(species(i_atom))
      !   elec_atomic(i_atom) = elec_atomic(i_atom)/species_z(species(i_atom))
      !enddo
  endif
  
  

  !     clean up for consistency

  if (allocated(temp_nucl)) then
     deallocate(temp_nucl)
  end if
  if (allocated(temp_free_rho)) then
     deallocate(temp_free_rho)
  end if
  if(allocated(temp_free_pot_es)) then
     deallocate(temp_free_pot_es)
  end if
  if(allocated(temp_hartree_sum)) then
     deallocate(temp_hartree_sum)
  end if
  if(allocated(center_index)) then
     deallocate(center_index)
  end if
  if (allocated(atomic_hartree_energy)) then
     deallocate(atomic_hartree_energy)
  end if

!  if(allocated( rho_inc_partialcore  )) then 
!     deallocate( rho_inc_partialcore  )
!  endif

!  if(allocated( rho_gradient_inc_partialcore )) then 
!     deallocate( rho_gradient_inc_partialcore )
!  endif
end subroutine get_free_superpos_energy_p1
!******	
