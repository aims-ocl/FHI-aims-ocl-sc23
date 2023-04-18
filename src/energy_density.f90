!****h* FHI-aims/energy_density
!  NAME
!    energy_density - energy_density for FHI-aims
!  SYNOPSIS

module energy_density

!  PURPOSE
!    This module takes care of all routines related to the energy density
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  AUTHOR
!    Christian Carbogno,  Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Yet unpublished and unreleased
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  SOURCE
  implicit none

save

!Kinetic energy density (DIM = n_full_points)
real*8,dimension(:),allocatable      :: ed_kinetic_energy_density
real*8,dimension(:,:,:),allocatable  :: ed_kinetic_batch

!XC energy density (DIM = n_full_points)
real*8,dimension(:),allocatable      :: ed_xc_energy_density

! Full local electrostatic Potential as it enters the KS equations (DIM = n_full_points)
real*8,dimension(:),allocatable      :: ed_local_potential

! Hartree Potential (free atoms / delta / corrections) (DIM = n_full_points)
real*8,dimension(:),allocatable      :: ed_energy_density_elec
real*8,dimension(:),allocatable      :: ed_hartree_energy_free
real*8,dimension(:),allocatable      :: ed_hartree_delta_energy
real*8,dimension(:),allocatable      :: ed_hartree_multipole_correction

! Nuclear contributions to the Hartree potential (DIM = n_atoms)
! i.e. energy of the single nuclei in the field of the free / delta atoms
real*8,dimension(:),allocatable      :: ed_energy_density_nuclei
real*8,dimension(:),allocatable      :: ed_hartree_delta_energy_nuclei
real*8,dimension(:),allocatable      :: ed_hartree_energy_free_nuclei

! EEV density (DIM = i_spin,n_full_points)
real*8,dimension(:,:),allocatable    :: ed_eev_energy_density
real*8,dimension(:,:,:),allocatable  :: ed_eev_times_kweights
real*8                               :: ed_eev_shift              
real*8,dimension(:,:),allocatable    :: ed_null1                  
real*8,dimension(:,:,:),allocatable  :: ed_null2                  
real*8,dimension(:,:,:),allocatable  :: ed_null3                  
real*8,dimension(:,:), allocatable   :: ed_null4
real*8,dimension(:,:), allocatable   :: ed_null5

! xc-potential density (DIM = n_full_points)
real*8,dimension(:),allocatable      :: ed_xc_pot_energy_density
integer,dimension(:),allocatable     :: ed_i_point_to_i_full_point_map



! Flag to make sure that update density and forces does not
! override virial while updating ed_eev_energy_density
logical                                :: flag_updating_ed

! Total values, needed for heat flux
real*8                          :: ed_energy_density_elec_sum           
real*8                          :: ed_energy_density_nuclei_sum         

! Counter
integer :: timestep = 0

contains


! Write Chetty & Martin Values to file
subroutine ed_write_chetty_martin_values(partition_tab)
!  PURPOSE
!  Subroutine ed_write_chetty_martin_values
!
!  Writes the various components of the CM energy density to ed_energy_density.dat 
!  -- This subroutine is basically a copy of output_potential_p1.f90  
!
!  USES
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use mpi_tasks
  use mpi_utilities
  use synchronize_mpi_basic, only: sync_real_number
  use localorb_io, only: localorb_info, OL_norm, use_unit
!  use heat_flux
  implicit none

!  ARGUMENTS
  real*8,dimension(n_full_points) :: partition_tab

!  INPUTS
!   o partition_tab -- integrate to check for correctness
!  OUTPUT
!    none
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
  real*8 :: grid_coord(3)
  real*8 :: ed_kinetic_energy
  real*8 :: ed_xc_energy
  real*8 :: ed_hartree_multipole_correction_sum
  real*8 :: ed_local_potential_sum
  real*8 :: ed_hartree_energy_free_sum
  real*8 :: ed_hartree_energy_free_nuclei_sum
  real*8 :: ed_hartree_delta_energy_sum
  real*8 :: ed_hartree_delta_energy_nuclei_sum   
  real*8 :: ed_eev_energy_density_sum   
  real*8 :: ed_xc_pot_energy_density_sum

!  counters
  integer :: i_coord,i_point, i_my_batch, i_index, i_atom, i_spin, current_atom
 
! info string
  character*120 :: info_str

  ! begin work

  if (.not.(allocated(ed_energy_density_elec  ))) allocate(ed_energy_density_elec  (n_full_points))
  if (.not.(allocated(ed_energy_density_nuclei))) allocate(ed_energy_density_nuclei(n_atoms      ))

  !Set everything to 0.0d0
  ed_kinetic_energy                    = 0.0d0
  ed_xc_energy                         = 0.0d0
  ed_hartree_multipole_correction_sum  = 0.0d0
  ed_local_potential_sum               = 0.0d0
  ed_hartree_energy_free_sum           = 0.0d0
  ed_hartree_energy_free_nuclei_sum    = 0.0d0
  ed_hartree_delta_energy_sum          = 0.0d0
  ed_hartree_delta_energy_nuclei_sum   = 0.0d0
  ed_eev_energy_density_sum            = 0.0d0
  ed_xc_pot_energy_density_sum         = 0.0d0 
  ed_energy_density_elec(:)            = 0.0d0               
  ed_energy_density_nuclei(:)          = 0.0d0               
  ed_energy_density_elec_sum           = 0.0d0 
  ed_energy_density_nuclei_sum         = 0.0d0


  ! Loop over batches
  i_point = 0
  do i_my_batch = 1, n_my_batches, 1

        ! Loop over points in a single batch
        do i_index = 1, batches(i_my_batch)%size, 1
           
           i_point = i_point + 1

           ! calculate grid point coordinate
           grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)
           current_atom  = batches(i_my_batch) % points(i_index) % index_atom

           ! Integrate energy contributions
           ed_xc_energy = ed_xc_energy &
              + ed_xc_energy_density(i_point) * partition_tab(i_point)     
           ed_hartree_energy_free_sum = ed_hartree_energy_free_sum &
              + ed_hartree_energy_free(i_point) * partition_tab(i_point)
           ed_hartree_delta_energy_sum = &
              ed_hartree_delta_energy_sum &
            + ed_hartree_delta_energy(i_point) * partition_tab(i_point)
           ed_hartree_multipole_correction_sum = &
              ed_hartree_multipole_correction_sum &
            + ed_hartree_multipole_correction(i_point) * partition_tab(i_point)

           if (flag_chetty_martin_energy_density) then
             ed_kinetic_energy = ed_kinetic_energy &
                + ed_kinetic_energy_density(i_point) * partition_tab(i_point)
             ed_local_potential_sum = ed_local_potential_sum &
                + ed_local_potential(i_point) * partition_tab(i_point)
             ed_energy_density_elec(i_point) = &
                ed_kinetic_energy_density(i_point) &
              + ed_local_potential(i_point) + ed_xc_energy_density(i_point) &
              - 0.5d0*( ed_hartree_energy_free(i_point) &
                      + ed_hartree_delta_energy(i_point) )
             ed_energy_density_elec_sum = &
                ed_energy_density_elec_sum &
              + ed_energy_density_elec(i_point) * partition_tab(i_point)
           end if

           ! N.B. Due to the peculiarieties of update_density_* i_spin=contains the eev for the total density and i_spin=2 the ones for the magnetization
           !      density
           if (flag_harris_foulkes_energy_density) then
             ed_eev_energy_density_sum = &
                ed_eev_energy_density_sum &
              + ed_eev_energy_density(i_point,1) * partition_tab(i_point)
             ed_xc_pot_energy_density_sum = &
                ed_xc_pot_energy_density_sum &
              + ed_xc_pot_energy_density(i_point) &
              * partition_tab(i_point)
             ed_energy_density_elec(i_point)  = -1.0d0 &
                * ed_eev_energy_density(i_point,1) &
                + ed_xc_energy_density(i_point) &
                -  ed_xc_pot_energy_density(i_point) &
                - 0.5d0 * ( ed_hartree_energy_free(i_point) &
                          + ed_hartree_delta_energy(i_point) )
             ed_energy_density_elec_sum = &
                ed_energy_density_elec_sum &
              + ed_energy_density_elec(i_point) * partition_tab(i_point)
           end if
           
        end do ! end loop over points in a single batch
        
  end do !end loop over batches

  !Sync
  if (use_mpi) then
    call sync_real_number(ed_xc_energy                         )
    call sync_real_number(ed_hartree_multipole_correction_sum  )
    call sync_real_number(ed_hartree_energy_free_sum           )
    call sync_real_number(ed_hartree_delta_energy_sum          )
    call sync_real_number(ed_energy_density_elec_sum           )

    ! Not required, we sync ed_hartree_delta_energy_nuclei in 
    ! sum_up_whole_potential_p1 already and ed_hartree_energy_free_nuclei_sum
    ! right above
    ! sync_real_number(ed_hartree_delta_energy_nuclei_sum   )
    ! sync_real_number(ed_hartree_energy_free_nuclei_sum    )

    if (flag_chetty_martin_energy_density) then
      call sync_real_number(ed_kinetic_energy                    )
      call sync_real_number(ed_local_potential_sum               )
    end if

    if (flag_harris_foulkes_energy_density) then
      call sync_real_number(ed_eev_energy_density_sum            )
      call sync_real_number(ed_xc_pot_energy_density_sum         )
    end if

  end if

  ! Loop over atoms
  do i_atom=1,n_atoms

    ! Sync ed_hartree_energy_free_nuclei among processors for 
    ! consistency: ed_hartree_delta_energy_nuclei is synced in 
    ! sum_up_whole_potential_p1 already!
    if (use_mpi) then
      call sync_real_number(ed_hartree_energy_free_nuclei(i_atom)   )
    end if

    ! Integrate nuclear contributions
    ed_hartree_delta_energy_nuclei_sum = ed_hartree_delta_energy_nuclei_sum &
       + ed_hartree_delta_energy_nuclei(i_atom)
    ed_hartree_energy_free_nuclei_sum  = ed_hartree_energy_free_nuclei_sum &
       + ed_hartree_energy_free_nuclei(i_atom)

    !Compute full nuclear energy "density"
    ed_energy_density_nuclei(i_atom) = -0.5d0 & 
       * (ed_hartree_delta_energy_nuclei(i_atom) &
        + ed_hartree_energy_free_nuclei(i_atom))
    ed_energy_density_nuclei_sum = &
       ed_energy_density_nuclei_sum + ed_energy_density_nuclei(i_atom)

  end do


  write (info_str,'(2X,2A)')  "--------------------------------------------", &
     "-------------------------------------------------------------" 
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )


      
  ! Print integrated quantities to standard out
  if (flag_harris_foulkes_energy_density) then
    write (info_str,'(2X,2A)') "* Harris-Foulkes type decomposition of ", &
       "the energy density:                                               " 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,2A)')  "|-----------------------------------------", &
       "---------------------------------------------------------------" 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated eigenvalue              energy density   : ", &
       -1.0d0 * ed_eev_energy_density_sum, " Ha", & 
       -1.0d0 * ed_eev_energy_density_sum * hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') "| Integrated xc                      energy density   : ", &
             & ed_xc_energy,                                                                            " Ha",& 
             & ed_xc_energy*hartree,                                                                    " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated xc potential correction energy density   : ", &
       -1.0d0 * ed_xc_pot_energy_density_sum, " Ha",& 
       -1.0d0 * ed_xc_pot_energy_density_sum * hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated free-atom electrostatic energy density   : ", & 
        -0.5d0*(ed_hartree_energy_free_sum &
              + ed_hartree_energy_free_nuclei_sum), " Ha", & 
        -0.5d0*(ed_hartree_energy_free_sum &
              + ed_hartree_energy_free_nuclei_sum) * hartree,       " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated Hartree correction      energy density   : ", & 
       -0.5d0*(ed_hartree_delta_energy_sum &
             + ed_hartree_delta_energy_nuclei_sum), " Ha", & 
       -0.5d0*(ed_hartree_delta_energy_sum &
             + ed_hartree_delta_energy_nuclei_sum) * hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,2A)')  "|-----------------------------------------", &
       "---------------------------------------------------------------" 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Resulting total energy:                             : ", & 
       -1.0d0*ed_eev_energy_density_sum + ed_xc_energy &
       - ed_xc_pot_energy_density_sum & 
       -0.5d0*( ed_hartree_energy_free_sum &
              + ed_hartree_energy_free_nuclei_sum &
              + ed_hartree_delta_energy_sum &
              + ed_hartree_delta_energy_nuclei_sum ), " Ha", &
       (-1.0d0*ed_eev_energy_density_sum + ed_xc_energy &
        - ed_xc_pot_energy_density_sum & 
        -0.5d0*( ed_hartree_energy_free_sum &
               + ed_hartree_energy_free_nuclei_sum &
               + ed_hartree_delta_energy_sum &
               + ed_hartree_delta_energy_nuclei_sum ))*hartree," eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,2A)')  "|-----------------------------------------", &
     "---------------------------------------------------------------" 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Resulting total energy elec:                        : ", & 
       ed_energy_density_elec_sum,         " Ha", &
       ed_energy_density_elec_sum*hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Resulting total energy nucl:                        : ", & 
       ed_energy_density_nuclei_sum,         " Ha", &
       ed_energy_density_nuclei_sum*hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Resulting total energy totl:                        : ", & 
       (ed_energy_density_elec_sum+ed_energy_density_nuclei_sum), " Ha",&
       (ed_energy_density_elec_sum+ed_energy_density_nuclei_sum)*hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,2A)')  "------------------------------------------", &
       "---------------------------------------------------------------" 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  end if

  if (flag_chetty_martin_energy_density) then
    write (info_str,'(2X,2A)') "* Chetty-Martin type decomposition of the ", &
       "energy density:                                               " 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,2A)')  "|-----------------------------------------", &
       "---------------------------------------------------------------" 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated kinetic                 energy density   : ", &
       ed_kinetic_energy, " Ha",& 
       ed_kinetic_energy*hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated electrostatic           energy density   : ", &
       ed_local_potential_sum & 
       -0.5d0*( ed_hartree_energy_free_sum &
              + ed_hartree_energy_free_nuclei_sum &
              + ed_hartree_delta_energy_sum &
              + ed_hartree_delta_energy_nuclei_sum ), " Ha", &
       (ed_local_potential_sum & 
        -0.5d0*( ed_hartree_energy_free_sum &
               + ed_hartree_energy_free_nuclei_sum &
               + ed_hartree_delta_energy_sum &
               + ed_hartree_delta_energy_nuclei_sum ))*hartree," eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated xc                      energy density   : ", &
       ed_xc_energy, " Ha", &
       ed_xc_energy * hartree, " eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,2A)')  "|-----------------------------------------", &
       "---------------------------------------------------------------" 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Resulting total energy:                             : ", & 
       ed_kinetic_energy+ed_local_potential_sum+ed_xc_energy &
       -0.5d0*( ed_hartree_energy_free_sum &
              + ed_hartree_energy_free_nuclei_sum &
              + ed_hartree_delta_energy_sum &
              + ed_hartree_delta_energy_nuclei_sum ),         " Ha", &
       (ed_kinetic_energy + ed_local_potential_sum + ed_xc_energy &
       -0.5d0*( ed_hartree_energy_free_sum &
              + ed_hartree_energy_free_nuclei_sum &
              + ed_hartree_delta_energy_sum &
              + ed_hartree_delta_energy_nuclei_sum ))*hartree," eV"
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
    write (info_str,'(2X,2A)')  "------------------------------------------", &
       "---------------------------------------------------------------" 
    call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  end if

  ! FIXME These analysis is mainly for debug purpose: It might be eliminated to unclutter the output
  write (info_str,'(2X,A)') "* Components of the electrostatic energy density:"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,2A)')  "|-------------------------------------------", &
     "-------------------------------------------------------------" 
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  if (flag_chetty_martin_energy_density) then
    write(info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
       "| Integrated local potential                          : ", &
       ed_local_potential_sum - ed_hartree_multipole_correction_sum, " Ha", & 
       (ed_local_potential_sum - ed_hartree_multipole_correction_sum) &
       * hartree, " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  end if
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
     "| Multipole correction to local potential             : ", &
     ed_hartree_multipole_correction_sum," Ha", &
     ed_hartree_multipole_correction_sum *hartree, " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,2A)')  "--------------------------------------------", &
   "-------------------------------------------------------------" 
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
     "| Electronic hartree energy:                          : ", &
     ed_hartree_energy_free_sum + ed_hartree_delta_energy_sum, " Ha", &
    (ed_hartree_energy_free_sum + ed_hartree_delta_energy_sum) * hartree, " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
     "|-> from free atoms                                   : ", &
     ed_hartree_energy_free_sum, " Ha", &
     ed_hartree_energy_free_sum * hartree, " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
     "|-> Hartree energy correction                         : ", &
     ed_hartree_delta_energy_sum, " Ha", &
     ed_hartree_delta_energy_sum * hartree, " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,2A)')  "|-------------------------------------------", &
     "-------------------------------------------------------------" 
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
     "| Nuclear    hartee  energy:                          : ", &
     ed_hartree_energy_free_nuclei_sum + ed_hartree_delta_energy_nuclei_sum, &
     " Ha", &
     (ed_hartree_energy_free_nuclei_sum + ed_hartree_delta_energy_nuclei_sum) &
     * hartree,  " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,A,A,1X,F20.8,A,1X,F20.8,A)') "|-> from free atoms     ", & 
     "                             : ", ed_hartree_energy_free_nuclei_sum,  & 
     " Ha", ed_hartree_energy_free_nuclei_sum   *hartree, " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
     "|-> Hartree energy correction                         : ", &
     ed_hartree_delta_energy_nuclei_sum, " Ha", &
     ed_hartree_delta_energy_nuclei_sum  * hartree, " eV"
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )
  write (info_str,'(2X,2A)')  "--------------------------------------------", &
     "-------------------------------------------------------------" 
  call localorb_info ( info_str, use_unit,'(A)', OL_norm  )

  call ed_deallocations()

end subroutine ed_write_chetty_martin_values

! Deallocate all arrays used for the calculation of the energy density
subroutine ed_deallocations
  implicit none

! TODO : CHECK THAT REALLY EVERYTHING GETS DEALLOCATED - IF NOT, CODE WILL
! FAIL DURING MD/RELAXATION DUE TO WRONG SIZED ALLOCS
!Kinetic energy density (DIM = n_full_points)
  if (allocated(ed_kinetic_energy_density))           deallocate(ed_kinetic_energy_density)      
  if (allocated(ed_kinetic_batch))                    deallocate(ed_kinetic_batch)               
                                                                                                  
!XC energy density (DIM = n_full_points)
  if (allocated(ed_xc_energy_density))                deallocate(ed_xc_energy_density)           
                                                                                                  
! Full local electrostatic Potential as it enters the KS equations (DIM = n_full_points)
  if (allocated(ed_local_potential))  deallocate(ed_local_potential) 

! Hartree Potential (free atoms / delta / corrections) (DIM = n_full_points)
  if (allocated(ed_energy_density_elec)) deallocate(ed_energy_density_elec) 
  if (allocated(ed_hartree_energy_free)) deallocate(ed_hartree_energy_free) 
  if (allocated(ed_hartree_delta_energy)) deallocate(ed_hartree_delta_energy)
  if (allocated(ed_hartree_multipole_correction)) &
     deallocate(ed_hartree_multipole_correction)

! Nuclear contributions to the Hartree potential (DIM = n_atoms)
! i.e. energy of the single nuclei in the field of the free / delta atoms
  if (allocated(ed_energy_density_nuclei)) &
     deallocate(ed_energy_density_nuclei)      
  if (allocated(ed_hartree_delta_energy_nuclei)) &
     deallocate(ed_hartree_delta_energy_nuclei) 
  if (allocated(ed_hartree_energy_free_nuclei)) &
     deallocate(ed_hartree_energy_free_nuclei)  

! EEV density (DIM = i_spin,n_full_points)
  if (allocated(ed_eev_energy_density)) deallocate(ed_eev_energy_density)
  if (allocated(ed_eev_times_kweights)) deallocate(ed_eev_times_kweights)
  if (allocated(ed_null1)) deallocate(ed_null1)
  if (allocated(ed_null2)) deallocate(ed_null2)
  if (allocated(ed_null3)) deallocate(ed_null3)
  if (allocated(ed_null4)) deallocate(ed_null4)
  if (allocated(ed_null5)) deallocate(ed_null5)

! xc-potential density (DIM = n_full_points)
  if (allocated(ed_xc_pot_energy_density)) deallocate(ed_xc_pot_energy_density)
  if (allocated(ed_i_point_to_i_full_point_map)) &
     deallocate(ed_i_point_to_i_full_point_map)

end subroutine ed_deallocations

subroutine ed_construct_pseudo_occupation( occ_numbers, KS_eigenvalue, pseudo_occupation, e_shift )
  use constants, only: hartree
  use dimensions
  use mpi_tasks
  use localorb_io
  implicit none
  real*8, dimension(n_states, n_spin, n_k_points),intent(INOUT)  :: occ_numbers
  real*8, dimension(n_states, n_spin, n_k_points),intent(IN)     :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points),intent(OUT)    :: pseudo_occupation
  real*8, intent(OUT)                                            :: e_shift
  
  ! local 
  integer       :: i_k_point,i_spin,i_state
  character*240 :: info_str
  real*8        :: high_occ_eev
 
  ! initialize
  pseudo_occupation(:,:,:) = 0.0d0
  
  !lowest/high_occ eev
  e_shift = 10.0d0
  high_occ_eev = KS_eigenvalue(1,1,1) - 1000.0d0

  !check for lowest eigenvalue and highest occupied eev
  do i_k_point = 1, n_k_points, 1
     do i_spin = 1, n_spin, 1
       do i_state = 1, n_states, 1
          if ( KS_eigenvalue(i_state,i_spin,i_k_point) .lt. e_shift ) &
             e_shift = KS_eigenvalue(i_state,i_spin,i_k_point)
          if ( ( KS_eigenvalue(i_state,i_spin,i_k_point) .gt. high_occ_eev ) &
               .and. &
               (occ_numbers(i_state,i_spin,i_k_point) .gt. 0.0d0) ) &
          then
             high_occ_eev = KS_eigenvalue(i_state,i_spin,i_k_point)
          end if
       end do
     end do
  end do

  !Same sign of low. and high. state?
  if ( ( e_shift * high_occ_eev ) .gt. 0.0d0 ) then
    !Do noting, e_shift becomes 0.0d0
    e_shift = 0.0d0
    ! loop over all three indexes
    do i_k_point = 1, n_k_points, 1
       do i_spin = 1, n_spin, 1
         do i_state = 1, n_states, 1
            pseudo_occupation(i_state,i_spin,i_k_point)  = -1.0d0 &
               * occ_numbers(i_state,i_spin,i_k_point) &
               * KS_eigenvalue(i_state,i_spin,i_k_point)
            
            ! Check for negative pseudo_occupation: 
            if ( pseudo_occupation(i_state,i_spin,i_k_point) .lt. 0.0d0 ) then
              write(info_str,'(A,I8,A,I2,A,I8,A)') &
                 ' *** ERROR IN ENERGY DENS.: Occupied (!) state ',i_state, &
                 ' for spin ', i_spin, ' at k-point ', i_k_point, &
                 ' still has a positive eigenvalue.'
              call aims_stop(info_str,'ed_construct_pseudo_occupation')
            end if 

         end do
       end do
    end do
  else
   ! Shift by low. eev + safety range
    e_shift = -1.0d0*(e_shift - 0.1d0)
    write(info_str,'(A,F20.8,A,F20.8,A)') &
       "  = Highest occupied state has positive EEV (", &
       high_occ_eev*hartree, " eV ) – Shifting whole spectra by", &
       e_shift*hartree," eV."
    call localorb_info(info_str)
    ! loop over all three indexes
    do i_k_point = 1, n_k_points, 1
       do i_spin = 1, n_spin, 1
         do i_state = 1, n_states, 1
            pseudo_occupation(i_state,i_spin,i_k_point) = &
               occ_numbers(i_state,i_spin,i_k_point) &
             * (KS_eigenvalue(i_state,i_spin,i_k_point) + e_shift)
            
            ! Check for negative pseudo_occupation: 
            if ( pseudo_occupation(i_state,i_spin,i_k_point) .lt. 0.0d0 ) then
              write(info_str,'(A,I8,A,I2,A,I8,A)') &
                 ' *** ERROR IN SHIFT ENERGY DENS.: Occupied (!) state ', &
                 i_state,' for spin ', i_spin, ' at k-point ', i_k_point, &
                 ' still has a positive eigenvalue.'
              call aims_stop(info_str,'ed_construct_pseudo_occupation')
            end if 

         end do
       end do
    end do
  end if

end subroutine ed_construct_pseudo_occupation

subroutine ed_eev_shift_back(shifted_eev_energy_density,e_shift,density)
  use dimensions
  use grids
  use mpi_tasks
  use mpi_utilities
  implicit none
  real*8, dimension(n_full_points,n_spin),intent(INOUT)  :: shifted_eev_energy_density
  real*8, intent(IN)                                     :: e_shift
  real*8, dimension(n_full_points,n_spin),intent(IN)     :: density
  
  ! local 
  integer       :: i_point,i_index,i_my_batch,i_spin
  character*240 :: info_str

  ! Loop over batches
  i_point = 0
  do i_my_batch = 1, n_my_batches, 1

        ! Loop over points in a single batch
        do i_index = 1, batches(i_my_batch)%size, 1
           
           i_point = i_point + 1
           
           do i_spin=1,n_spin
             shifted_eev_energy_density(i_point,i_spin) = -1.0d0 &
                * (shifted_eev_energy_density(i_point,i_spin) &
                - ( e_shift * density(i_point,i_spin) ) )
           end do 
           
        end do ! end loop over points in a single batch
        
  end do !end loop over batches

end subroutine ed_eev_shift_back


end module energy_density
