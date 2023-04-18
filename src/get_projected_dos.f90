!****s* FHI-aims/get_species_projected_dos
!  NAME
!   get_species_projected_dos
!  SYNOPSIS

      subroutine get_species_projected_dos &
                 (dostype, i_first_level, n_levels,energy_levels )

!  PURPOSE
!  The  subroutine calculates the projected density of states.
!  This means the density of states separated to contributions
!  of the different species.
!
!  USES

      use dimensions
      use basis
      use species_data
      use geometry
      use physics
      use mpi_tasks
      use runtime_choices
      use localorb_io, only : use_unit
      implicit none

!  ARGUMENTS


      character*2 dostype
      integer :: i_first_level
      integer :: n_levels
      real*8 :: energy_levels(n_levels,n_spin)


!  INPUTS
!  o dostype -- character for file name
!  o i_first_level -- the first of the states included for DOS calculations
!  o n_levels -- number of energy levels printed out
!  o energy_levels -- energy levels
! 
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

!    local variables
      real*8, dimension(:,:,:), allocatable :: KS_dos
      real*8, dimension(:,:), allocatable :: KS_levels
      real*8, dimension(:,:,:), allocatable :: projected_weight

!  parameters
      real*8  :: de
      real*8  :: en

      integer i_e
      integer i_species
      integer i_species_1
      integer i_species_2
      integer i_basis_1
      integer i_basis_2
      integer i_index
      integer i_spin
      integer i_state
      integer i_KS_state

      character*50 filename

      if (myid.eq.0) then
       write(use_unit,*) "-------------------------------------------------- "
       write(use_unit,*)
       write(use_unit,*) "  Calculating projected density of states ... "
       write(use_unit,*)
      endif

      if(.not. allocated (projected_weight)) then
        allocate(projected_weight(n_levels,n_species,n_spin))
      endif
      if(.not. allocated (KS_dos)) then
        allocate(KS_dos(n_spin,dos_n_en_points, n_species))
      endif
      if(.not. allocated (KS_levels)) then
        allocate(KS_levels(n_levels, n_spin))
      endif

      projected_weight(:,:,:) =0.d0
      do i_spin = 1, n_spin, 1
        do i_state = 1, n_levels , 1
          i_index = 0
          i_KS_state = i_state + i_first_level - 1
          do i_basis_1 = 1, n_basis, 1

           i_species_1 = species(basis_atom(i_basis_1))
           do i_basis_2 = 1, i_basis_1, 1

               i_species_2 = species(basis_atom(i_basis_2))
               i_index = i_index + 1
               projected_weight(i_state,i_species_1,i_spin) = &
                  projected_weight(i_state,i_species_1,i_spin) + &
                  KS_eigenvector(i_basis_1,i_KS_state,i_spin,1) * &
                  KS_eigenvector(i_basis_2,i_KS_state,i_spin,1) * &
                  overlap_matrix(i_index)

               if(i_basis_1.ne.i_basis_2) then
                 projected_weight(i_state,i_species_2,i_spin) = &
                    projected_weight(i_state,i_species_2,i_spin) + &
                    KS_eigenvector(i_basis_1,i_KS_state,i_spin,1) * &
                    KS_eigenvector(i_basis_2,i_KS_state,i_spin,1) * &
                    overlap_matrix(i_index)
               endif
            enddo
         enddo
       enddo
      enddo

      KS_levels(:,:) = energy_levels(:,:)
      call  evaluate_QP_dos (n_levels, dos_n_en_points, &
                         n_spin, n_species, &
                         dos_low_energy,dos_high_energy, &
                         dos_alpha, &
                         KS_levels, &
                         projected_weight, &
                         KS_dos)

      if(myid.eq.0) then
       de= (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)
       do i_species = 1, n_species

         write(filename,'(A,A1,A,A7)') &
               trim(species_name(i_species)), '_', &
               trim(dostype), 'dos.dat'

         open (110, file = filename)
         do i_e = 1, dos_n_en_points ,1

            en = dos_low_energy + dble(i_e -1) *de
            write(110,'(3f20.12)') en, &
                (KS_dos(i_spin,i_e,i_species), i_spin=1,n_spin)
         enddo
         close(110)
       enddo
      endif

      if(allocated(KS_dos)) then
        deallocate(KS_dos)
      endif
      if(allocated(KS_levels)) then
        deallocate(KS_levels)
      endif
      if(allocated(projected_weight)) then
        deallocate(projected_weight)
      endif
    end subroutine get_species_projected_dos
!******		
!-------------------------------------------------------------------------------
!****s* FHI-aims/evaluate_QP_dos
!  NAME
!   evaluate_QP_dos
!  SYNOPSIS

    subroutine  evaluate_QP_dos (n_levels,dos_n_en_points, &
         n_spin,n_species, &
         dos_low_energy,dos_high_energy, &
         dos_alpha, &
         KS_levels,projected_weight, &
         KS_dos)


!  PURPOSE
!   This subroutine is to calculate the density of states of the KS states
!
!  USES
      use constants
      implicit none
!  ARGUMENTS

      integer :: n_levels
      integer :: dos_n_en_points
      integer :: n_species
      integer :: n_spin

      real*8  :: dos_low_energy
      real*8  :: dos_high_energy
      real*8  :: dos_alpha
      real*8  :: KS_levels(n_levels,n_spin)
      real*8  :: projected_weight(n_levels,n_species,n_spin)

      real*8  :: KS_dos(n_spin,dos_n_en_points,n_species)

!  INPUTS
!  o n_levels -- number of energy levels
!  o dos_n_en_points -- number of energy points in the DOS
!  o n_species -- number of species
!  o n_spin -- number of spin
!  o dos_low_energy -- lower bound of the energy range
!  o dos_high_energy -- higher bound of the energy range
!  o dos_alpha --  Gaussian broadening factor
!  o KS_levels --  Kohn-Sham eigenvalues
!  o projected_weight -- projection weight
!  
!  OUTPUT
!  o KS_dos -- density of states
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






!    local variables
!     n_en    : number of energy points
!     en_low  : low bound of the energy spectrum
!     en_high : upper bound of the energy spectrum
!     de      : energy interval

      real*8  de
      real*8  en

!     counter
      integer i_e
      integer i_state
      integer i_species
      integer i_spin

!     begin work
!      dos_alpha=1.0

      KS_levels(:,:) = KS_levels(:,:)*hartree

      de= (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)
      KS_dos =0.d0
      do i_spin = 1, n_spin, 1
       do i_species = 1, n_species, 1
         do i_e = 1, dos_n_en_points, 1
          en= dos_low_energy + (i_e-1)*de

          do i_state = 1, n_levels, 1
            KS_dos(i_spin,i_e,i_species) = &
             KS_dos (i_spin,i_e,i_species) + &
             abs(projected_weight(i_state,i_species,i_spin)) * &
             exp(-((en-KS_levels(i_state,i_spin))/dos_alpha)**2/2.d0)/ &
             sqrt(2.d0*pi)/dos_alpha

!            write(use_unit,*) en, KS_dos(i_e,i_species)
          enddo
        enddo
       enddo
      enddo

      return
    end subroutine evaluate_QP_dos
!******		
