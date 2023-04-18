!****h* FHI-aims/vdw_correction
!  NAME
!    vdw_correction - module to handle $C_6R^{-6}$ vdW correction to DFT
!    energies and forces
!  SYNOPSIS

module vdw_correction

  !  PURPOSE
  !
  !    This module takes care of
  !    * vdW correction
  !
  !    C6/R^6 van der Waals correction. Two versions currently implemented:
  !    1) Empirical (C6 and R_vdW must be specified in the control file).
  !    Specify "use_vdw_correction" in the control file.
  !    Parameters should be provided in the control.in file, as follows:
  !    * vdw_correction
  !    * vdw_pairs N                  N is the number of atomic pairs
  !    * vdw_coeff species1_name species2_name C6 R0 d
  !    * .... (repeat N times)
  !    
  !    C6, R0 and d should be given appropriately for each functional.
  !  
  !    2) Non-empirical correction based on Hirshfeld volume partitioning,
  !    from which we extract the effective polarizability AND vdW radii.
  !    Only "use_vdw_correction_hirshfeld" must be specified in the control file.
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
  
  use mbd, only: mbd_calc_t
  
  implicit none

  ! vdw variables
  character*20, dimension(:),   allocatable :: vdw_species1
  character*20, dimension(:),   allocatable :: vdw_species2
  character*20, dimension(:,:), allocatable :: vdw_pairs_ignore
  real*8,  dimension(:), allocatable :: vdw_C6
  real*8,  dimension(:), allocatable :: vdw_R0
  real*8,  dimension(:), allocatable :: vdw_d
  real*8,  dimension(:), allocatable :: vdw_hirshfeld_C6
  real*8,  dimension(:), allocatable :: vdw_hirshfeld_alpha
  real*8,  dimension(:), allocatable :: vdw_hirshfeld_R0
  real*8,  dimension(:,:), allocatable :: vdw_C6_prefactor ! 1 by default, 0 if interaction (species1, species2) was turned off.
  logical, dimension(:), allocatable :: vdw_hirshfeld_data_external
  integer :: vdw_n1
  integer :: vdw_n2
  integer :: vdw_n3
 
  ! logical check for the hirshfeld analysis
  logical :: hirsh_analysis
  
  real*8, dimension(:), allocatable :: hirshfeld_volume
  real*8, dimension(:), allocatable :: hirshfeldw
  real*8, dimension(:), allocatable :: freeintegral

  integer :: lattice_translations
  real*8, dimension(:,:,:,:,:), allocatable :: table_f_ab
  real*8, dimension(:), allocatable :: table_f_atom
  !******

  type(mbd_calc_t) :: mbd_calc

contains
!-----------------------------------------------------

  subroutine allocate_vdw()

    use dimensions

    ! allocations for empirical vdw correction
    if (vdw_pairs .gt. 0) then
       if (.not.allocated(vdw_species1)) allocate (vdw_species1(vdw_pairs))
       if (.not.allocated(vdw_species2)) allocate (vdw_species2(vdw_pairs))
       if (.not.allocated(vdw_C6))       allocate (vdw_C6(vdw_pairs))
       if (.not.allocated(vdw_R0))       allocate (vdw_R0(vdw_pairs))
       if (.not.allocated(vdw_d))        allocate (vdw_d(vdw_pairs))
    end if

    ! these variables can also be declared as species defaults
    if (.not.allocated(vdw_hirshfeld_C6           )) allocate(vdw_hirshfeld_C6(n_species))
    if (.not.allocated(vdw_hirshfeld_alpha        )) allocate(vdw_hirshfeld_alpha(n_species))
    if (.not.allocated(vdw_hirshfeld_R0           )) allocate(vdw_hirshfeld_R0(n_species))
    if (.not.allocated(vdw_hirshfeld_data_external)) allocate(vdw_hirshfeld_data_external(n_species))
    if (.not.allocated(vdw_C6_prefactor           )) allocate(vdw_C6_prefactor(n_species, n_species))
    if ((.not.allocated(vdw_pairs_ignore)).and.(n_vdw_pairs_ignore.gt.0)) then
       allocate(vdw_pairs_ignore(2,n_vdw_pairs_ignore))
    end if
    vdw_hirshfeld_data_external(:) = .false.
    vdw_C6_prefactor(:,:) = 1d0

    ! only allocate if hirshfeld correction is actually being used
    if (use_vdw_correction_hirshfeld.or.use_mbd_old&
        .or. use_mbd_dev .or. use_libmbd) then
       if (.not.allocated(hirshfeld_volume)) allocate (hirshfeld_volume(n_atoms))
    endif
    if (use_vdw_correction_hirshfeld_sc .or. use_mbd_std) then
       if (.not.allocated(hirshfeld_volume)) allocate (hirshfeld_volume(n_atoms))
       if (.not.allocated(hirshfeldw)) allocate (hirshfeldw(n_atoms))
       if (.not.allocated(freeintegral)) allocate (freeintegral(n_atoms))
    endif

  end subroutine allocate_vdw

  subroutine deallocate_vdw()

    use dimensions

    if(allocated (vdw_species1)    ) deallocate (vdw_species1)
    if(allocated (vdw_species2)    ) deallocate (vdw_species2)
    if(allocated (vdw_C6)          ) deallocate (vdw_C6)
    if(allocated (vdw_R0)          ) deallocate (vdw_R0)
    if(allocated (vdw_d)           ) deallocate (vdw_d)
    if(allocated (hirshfeld_volume)) deallocate (hirshfeld_volume)
    if(allocated (hirshfeldw)) deallocate (hirshfeldw)
    if(allocated (freeintegral)) deallocate (freeintegral)
    if (allocated(vdw_hirshfeld_C6           )) deallocate(vdw_hirshfeld_C6)
    if (allocated(vdw_hirshfeld_alpha        )) deallocate(vdw_hirshfeld_alpha)
    if (allocated(vdw_hirshfeld_R0           )) deallocate(vdw_hirshfeld_R0)
    if (allocated(vdw_hirshfeld_data_external)) deallocate(vdw_hirshfeld_data_external)
    if (allocated(vdw_pairs_ignore           )) deallocate(vdw_pairs_ignore)
    if (allocated(vdw_C6_prefactor           )) deallocate(vdw_C6_prefactor)


end subroutine deallocate_vdw


subroutine read_vdw(vdw_pairs)

use localorb_io

implicit none

! Imported variables

integer :: vdw_pairs

! local variables
integer :: i_code
character*20 :: desc_str
character*20 :: str
character*160 :: info_str
real*8 :: C6
real*8 :: R0
real*8 :: d
character*2  :: atom1
character*2  :: atom2

!  counters
integer :: i_pair

read(7,*,iostat=i_code) desc_str
if (desc_str/="vdw_pairs") then
   write(use_unit,*)"Syntax error in control.in - expected 'vdw_pairs' line"
   stop
endif
do i_pair = 1, vdw_pairs, 1
   read(7,*) str,atom1,atom2,C6,R0,d
   vdw_species1(i_pair) = atom1
   vdw_species2(i_pair) = atom2
   vdw_C6(i_pair) = C6
   vdw_R0(i_pair) = R0
   vdw_d(i_pair) = d
   write(info_str,'(2X,A,A,A,A,A,F10.5,A,F10.5,A,F10.5)') &
      "VDW parameters for species ",atom1," and ",atom2," (C6,R0,d): "&
      ,C6," ",R0," ",d
   call localorb_info(info_str)
enddo

end subroutine read_vdw

! Free-atom C6, polarizability and vdW radius
subroutine get_vdw_param(element, dummy, C6, alpha, R0)
    use dimensions
    use species_data
    use localorb_io


    !VVG: The majority of values such as isotropic static polarizability 
    !(in bohr^3), the homo-atomic van der Waals coefficient(in hartree·bohr^6), 
    !and vdW Radii (in bohr) for neutral free atoms are taken from Ref. Chu, X. & Dalgarno, 
    !J. Chem. Phys. 121, 4083 (2004) and Mitroy, et al. 
    !J. Phys. B: At. Mol. Opt. Phys. 43, 202001 (2010) 
    !and for rest of the elements they are calculated using linear response coupled cluster
    !single double theory with accrate basis. The vdW radii for respective element are 
    !defined as discussed in Tkatchenko, A. & Scheffler, M. Phys. Rev. Lett. 102, 073005 (2009). 
    !
    !For lathanides and actinide, the C_{6} coefficients are constructed to 
    !satisfy valence electronic sum rule using the static polarizabilities from Ref
    !Dzuba, et al Phys. Rev. A 89, 042507 (2014).
    !using one-term formula which is equivalent to the two-point zeroth-order Pad\'e approximant.

    implicit none

    character(len=2), intent(in) :: element
    real(8), intent(in) :: dummy  ! keep old function signature
    real(8), intent(out) :: C6, alpha, R0

    select case (periodic_table(periodic_table(element)))  ! capitalize
    case ('H');  alpha =   4.5000; C6 =    6.5000; R0 = 3.1000
    case ('He'); alpha =   1.3800; C6 =    1.4600; R0 = 2.6500
    case ('Li'); alpha = 164.2000; C6 = 1387.0000; R0 = 4.1600
    case ('Be'); alpha =  38.0000; C6 =  214.0000; R0 = 4.1700
    case ('B');  alpha =  21.0000; C6 =   99.5000; R0 = 3.8900
    case ('C');  alpha =  12.0000; C6 =   46.6000; R0 = 3.5900
    case ('N');  alpha =   7.4000; C6 =   24.2000; R0 = 3.3400
    case ('O');  alpha =   5.4000; C6 =   15.6000; R0 = 3.1900
    case ('F');  alpha =   3.8000; C6 =    9.5200; R0 = 3.0400
    case ('Ne'); alpha =   2.6700; C6 =    6.3800; R0 = 2.9100
    case ('Na'); alpha = 162.7000; C6 = 1556.0000; R0 = 3.7300
    case ('Mg'); alpha =  71.0000; C6 =  627.0000; R0 = 4.2700
    case ('Al'); alpha =  60.0000; C6 =  528.0000; R0 = 4.3300
    case ('Si'); alpha =  37.0000; C6 =  305.0000; R0 = 4.2000
    case ('P');  alpha =  25.0000; C6 =  185.0000; R0 = 4.0100
    case ('S');  alpha =  19.6000; C6 =  134.0000; R0 = 3.8600
    case ('Cl'); alpha =  15.0000; C6 =   94.6000; R0 = 3.7100
    case ('Ar'); alpha =  11.1000; C6 =   64.3000; R0 = 3.5500
    case ('K');  alpha = 292.9000; C6 = 3897.0000; R0 = 3.7100
    case ('Ca'); alpha = 160.0000; C6 = 2221.0000; R0 = 4.6500
    case ('Sc'); alpha = 120.0000; C6 = 1383.0000; R0 = 4.5900
    case ('Ti'); alpha =  98.0000; C6 = 1044.0000; R0 = 4.5100
    case ('V');  alpha =  84.0000; C6 =  832.0000; R0 = 4.4400
    case ('Cr'); alpha =  78.0000; C6 =  602.0000; R0 = 3.9900
    case ('Mn'); alpha =  63.0000; C6 =  552.0000; R0 = 3.9700
    case ('Fe'); alpha =  56.0000; C6 =  482.0000; R0 = 4.2300
    case ('Co'); alpha =  50.0000; C6 =  408.0000; R0 = 4.1800
    case ('Ni'); alpha =  48.0000; C6 =  373.0000; R0 = 3.8200
    case ('Cu'); alpha =  42.0000; C6 =  253.0000; R0 = 3.7600
    case ('Zn'); alpha =  40.0000; C6 =  284.0000; R0 = 4.0200
    case ('Ga'); alpha =  60.0000; C6 =  498.0000; R0 = 4.1900
    case ('Ge'); alpha =  41.0000; C6 =  354.0000; R0 = 4.2000
    case ('As'); alpha =  29.0000; C6 =  246.0000; R0 = 4.1100
    case ('Se'); alpha =  25.0000; C6 =  210.0000; R0 = 4.0400
    case ('Br'); alpha =  20.0000; C6 =  162.0000; R0 = 3.9300
    case ('Kr'); alpha =  16.8000; C6 =  129.6000; R0 = 3.8200
    case ('Rb'); alpha = 319.2000; C6 = 4691.0000; R0 = 3.7200
    case ('Sr'); alpha = 199.0000; C6 = 3170.0000; R0 = 4.5400
    case ('Y');  alpha = 126.7370; C6 = 1968.580;  R0 = 4.8151
    case ('Zr'); alpha = 119.97;   C6 = 1677.91;   R0 = 4.53
    case ('Nb'); alpha = 101.603;  C6 = 1263.61;   R0 = 4.2365
    case ('Mo'); alpha =  88.4225; C6 = 1028.73;   R0 = 4.099
    case ('Tc'); alpha =  80.083;  C6 = 1390.87;   R0 = 4.076
    case ('Ru'); alpha =  65.8950; C6 =  609.754;  R0 = 3.9953
    case ('Rh'); alpha =  56.1;    C6 =  469.0;    R0 = 3.95
    case ('Pd'); alpha =  23.6800; C6 =  157.5000; R0 = 3.6600
    case ('Ag'); alpha =  50.6000; C6 =  339.0000; R0 = 3.8200
    case ('Cd'); alpha =  39.7;    C6 =  452.0;    R0 = 3.99
    case ('In'); alpha =  70.2200; C6 =  707.0460; R0 = 4.2319
    case ('Sn'); alpha =  55.9500; C6 =  587.4170; R0 = 4.3030
    case ('Sb'); alpha =  43.6719; C6 =  459.322;  R0 = 4.2760
    case ('Te'); alpha =  37.65;   C6 =  396.0;    R0 = 4.22
    case ('I');  alpha =  35.0000; C6 =  385.0000; R0 = 4.1700
    case ('Xe'); alpha =  27.3000; C6 =  285.9000; R0 = 4.0800
    case ('Cs'); alpha = 427.12;   C6 = 6582.08;   R0 = 3.78
    case ('Ba'); alpha = 275.0;    C6 = 5727.0;    R0 = 4.77
    case ('La'); alpha = 213.70;   C6 = 3884.5;    R0 = 3.14
    case ('Ce'); alpha = 204.7;    C6 = 3708.33;   R0 = 3.26
    case ('Pr'); alpha = 215.8;    C6 = 3911.84;   R0 = 3.28
    case ('Nd'); alpha = 208.4;    C6 = 3908.75;   R0 = 3.3
    case ('Pm'); alpha = 200.2;    C6 = 3847.68;   R0 = 3.27
    case ('Sm'); alpha = 192.1;    C6 = 3708.69;   R0 = 3.32
    case ('Eu'); alpha = 184.2;    C6 = 3511.71;   R0 = 3.40
    case ('Gd'); alpha = 158.3;    C6 = 2781.53;   R0 = 3.62
    case ('Tb'); alpha = 169.5;    C6 = 3124.41;   R0 = 3.42
    case ('Dy'); alpha = 164.64;   C6 = 2984.29;   R0 = 3.26
    case ('Ho'); alpha = 156.3;    C6 = 2839.95;   R0 = 3.24
    case ('Er'); alpha = 150.2;    C6 = 2724.12;   R0 = 3.30
    case ('Tm'); alpha = 144.3;    C6 = 2576.78;   R0 = 3.26
    case ('Yb'); alpha = 138.9;    C6 = 2387.53;   R0 = 3.22
    case ('Lu'); alpha = 137.2;    C6 = 2371.80;   R0 = 3.20
    case ('Hf'); alpha =  99.52;   C6 = 1274.8;    R0 = 4.21
    case ('Ta'); alpha =  82.53;   C6 = 1019.92;   R0 = 4.15
    case ('W');  alpha =  71.041;  C6 =  847.93;   R0 = 4.08
    case ('Re'); alpha =  63.04;   C6 =  710.2;    R0 = 4.02
    case ('Os'); alpha =  55.055;  C6 =  596.67;   R0 = 3.84
    case ('Ir'); alpha =  42.51;   C6 =  359.1;    R0 = 4.00
    case ('Pt'); alpha =  39.68;   C6 =  347.1;    R0 = 3.92
    case ('Au'); alpha =  36.5;    C6 =  298.0;    R0 = 3.86
    case ('Hg'); alpha =  33.9;    C6 =  392.0;    R0 = 3.98
    case ('Tl'); alpha =  69.92;   C6 =  717.44;   R0 = 3.91
    case ('Pb'); alpha =  61.8;    C6 =  697.0;    R0 = 4.31
    case ('Bi'); alpha =  49.02;   C6 =  571.0;    R0 = 4.32
    case ('Po'); alpha =  45.013;  C6 =  530.92;   R0 = 4.097
    case ('At'); alpha =  38.93;   C6 =  457.53;   R0 = 4.07
    case ('Rn'); alpha =  33.54;   C6 =  390.63;   R0 = 4.23
    case ('Fr'); alpha = 317.8;    C6 = 4224.44;   R0 = 3.90
    case ('Ra'); alpha = 246.2;    C6 = 4851.32;   R0 = 4.98
    case ('Ac'); alpha = 203.3;    C6 = 3604.41;   R0 = 2.75
    case ('Th'); alpha = 217.0;    C6 = 4047.54;   R0 = 2.85
    case ('Pa'); alpha = 154.4;    C6 = 2367.42;   R0 = 2.71
    case ('U');  alpha = 127.8;    C6 = 1877.10;   R0 = 3.00
    case ('Np'); alpha = 150.5;    C6 = 2507.88;   R0 = 3.28
    case ('Pu'); alpha = 132.2;    C6 = 2117.27;   R0 = 3.45
    case ('Am'); alpha = 131.20;   C6 = 2110.98;   R0 = 3.51
    case ('Cm'); alpha = 143.6;    C6 = 2403.22;   R0 = 3.47
    case ('Bk'); alpha = 125.3;    C6 = 1985.82;   R0 = 3.56
    case ('Cf'); alpha = 121.5;    C6 = 1891.92;   R0 = 3.55
    case ('Es'); alpha = 117.5;    C6 = 1851.1;    R0 = 3.76
    case ('Fm'); alpha = 113.4;    C6 = 1787.07;   R0 = 3.89
    case ('Md'); alpha = 109.4;    C6 = 1701.0;    R0 = 3.93
    case ('No'); alpha = 105.4;    C6 = 1578.18;   R0 = 3.78
    case default;alpha = 0.;       C6 = 0.;        R0 = 4.0
    end select
end subroutine get_vdw_param

subroutine build_C6_prefactor_matrix()
  use dimensions
  use species_data
  use localorb_io
  use mpi_tasks, only: aims_stop
  implicit none
  integer :: i_species, i_species1, i_species2, i_vdw_ignore
  logical :: found1, found2
  character*120 :: info_str
  ! build vdw_C6_prefactor matrix: coefficients are equal to unity by default and zero for all ignored pairs of species
  vdw_C6_prefactor(:,:) = 1d0
  if (n_vdw_pairs_ignore .gt. 0) then
     do i_vdw_ignore = 1, n_vdw_pairs_ignore   ! loop through all relevant pairs of species
        ! determine the species:
        found1 = .false.
        found2 = .false.
        do i_species = 1, n_species
           if (trim(species_name(i_species)).eq.vdw_pairs_ignore(1,i_vdw_ignore)) then
              found1 = .true.
              i_species1 = i_species
           end if
           if (trim(species_name(i_species)).eq.vdw_pairs_ignore(2,i_vdw_ignore)) then
              found2 = .true.
              i_species2 = i_species
           end if
        end do
        if ((.not.found1).or.(.not.found2)) then
           write(info_str,'(2X,4A)') '*** WARNING: Did not find following species for ignoring their vdw interaction: ', &
                trim(vdw_pairs_ignore(1,i_vdw_ignore)),' and ', trim(vdw_pairs_ignore(2,i_vdw_ignore))
           call localorb_info(info_str)
           write(info_str,'(2X,A)') '*** Please check input in control.in. Aborting execution. '
           call localorb_info(info_str)
           call aims_stop()
        end if
        vdw_C6_prefactor(i_species1, i_species2) = 0d0
        vdw_C6_prefactor(i_species2, i_species1) = 0d0
     end do
  end if
end subroutine build_C6_prefactor_matrix


  ! Empirical vdW routine
  subroutine calc_vdw(vdw_energy, ionic_forces)

    use dimensions
    use geometry
    use species_data
    use localorb_io
    use mpi_utilities
    use synchronize_mpi
    use constants

    implicit none

    ! Imported variables
    real*8 :: vdw_energy
    real*8, dimension(3, n_atoms) :: ionic_forces

    ! Local Variables
    real*8, dimension(3) :: coord_diff
    real*8 :: dist
    real*8 :: dist2
    real*8 :: dist6
    real*8 :: rpow
    real*8 :: rexp
    real*8 :: aexp
    real*8 :: rphi
    real*8 :: rfac
    real*8 :: ffac
    real*8 :: C6
    real*8 :: R0
    real*8 :: d
    character*160 :: info_str

    !Logical
    logical :: found

    !counters
    integer :: i_pair
    integer :: i_atom1
    integer :: i_atom2
    integer :: i_coord



    vdw_energy = 0.d0

    do i_atom1 = 1, n_atoms, 1
       if (myid.eq.task_list(i_atom1)) then
          do i_atom2 =  i_atom1 + 1, n_atoms, 1
             !          Find if the atoms pair is defined for vdW correction     
             found = .false.
             do i_pair = 1, vdw_pairs, 1
                if (.not. empty(i_atom1)) then
                   if ( (species_name(species(i_atom1)).eq.vdw_species1(i_pair)) .and. &
                        (species_name(species(i_atom2)).eq.vdw_species2(i_pair)) ) then
                      C6 = vdw_C6(i_pair)*vdw_C6_prefactor(species(i_atom1),species(i_atom2))
                      R0 = vdw_R0(i_pair)
                      d = vdw_d(i_pair)
                      found = .true.
                   endif

                   if ( (species_name(species(i_atom2)).eq.vdw_species1(i_pair)) .and. &
                        (species_name(species(i_atom1)).eq.vdw_species2(i_pair)) ) then
                      C6 = vdw_C6(i_pair)*vdw_C6_prefactor(species(i_atom1),species(i_atom2))
                      R0 = vdw_R0(i_pair)
                      d = vdw_d(i_pair)
                      found = .true.
                   endif

                endif
             end do

             !          If the atoms pair is defined, proceed with calculation
             if (found) then 
                !              write (info_str,'(2X,A,A,A,A,A,F10.5,F10.5,F10.5)')   &
                !                    "vdW param for species:" ,species_name(species(i_atom1)),"and ",&
                !                    species_name(species(i_atom2))," : ",C6,R0,d
                !             call localorb_info(info_str,use_unit,'(A)')

                do i_coord = 1, 3, 1
                   coord_diff(i_coord) = coords(i_coord,i_atom1) - coords(i_coord,i_atom2)
                enddo
                dist2 = coord_diff(1)*coord_diff(1) + &
                        coord_diff(2)*coord_diff(2) + &
                        coord_diff(3)*coord_diff(3)
             dist = dsqrt(dist2)
             dist6 = dist2*dist2*dist2 
             rpow=d*(dist/R0-1.0D0)
             rexp=dexp(-rpow)
             aexp=1.0d0+rexp
             rphi = -C6/(dist6*aexp)
             vdw_energy = vdw_energy + rphi 

!            vdW forces
             if (use_forces .and. n_periodic == 0) then
               rfac=-6.0d0/dist + d*rexp/(aexp*R0)
               ffac=rphi*rfac/dist
               do i_coord = 1, 3, 1
                 ionic_forces(i_coord, i_atom1) =   &
                 ionic_forces(i_coord, i_atom1) - coord_diff(i_coord)*ffac
                 ionic_forces(i_coord, i_atom2) =   &
                 ionic_forces(i_coord, i_atom2) + coord_diff(i_coord)*ffac
               end do
             endif
           endif ! found 
        end do
  endif
enddo

end subroutine calc_vdw


subroutine calc_vdw_hirshfeld(vdw_energy, vdw_forces, vdW_Hessian)

use dimensions
use runtime_choices
use geometry
use species_data
use localorb_io
use mpi_utilities
use synchronize_mpi
use constants
use analytical_stress
use heat_flux

implicit none

! Imported variables
real*8 :: vdw_energy
real*8, dimension(3, n_atoms) :: vdw_forces
real*8, dimension(n_atoms*3,n_atoms*3), optional, intent(out) :: vdw_Hessian


! Local Variables
real*8, dimension(3) :: coord_diff
real*8, dimension(3) :: coord_curr
real*8, dimension(3, n_atoms) :: vdw_forces_change
real*8 :: vdw_energy_change
real*8, dimension(3, 3) :: AS_vdw_change
real*8 :: dist
real*8 :: dist2
real*8 :: dist6
real*8 :: rpow
real*8 :: rexp
real*8 :: aexp
real*8 :: rphi
real*8 :: rfac
real*8 :: ffac
real*8 :: C6
real*8 :: R0
real*8 :: d
real*8 :: Sr
logical :: periodic_converged
integer :: periodic_cell
real*8 :: maxchange

real*8 :: alpha_array, R0_array, C6_array

real*8, dimension(n_atoms) :: alpha1,C61, R01 

character*160 :: info_str
character*2 :: atom
real*8 :: nucleus

real*8 :: dR0, dR0raexp, dEdr, d2Edr2
real*8 :: signum, dRdx1, dRdx2, d2Rdx1dx2

!counters
integer :: i_pair
integer :: i_atom1
integer :: i_atom2
integer :: i_coord
integer :: i_coord2
integer :: i1
integer :: i2
integer :: i3


vdw_energy = 0.d0
if (use_analytical_stress) then
   AS_vdw_stress(:,:) = 0.0d0
   if (compute_heat_flux) then
     if (.not.allocated(HF_stress_per_atom_vdw))        allocate(HF_stress_per_atom_vdw  (1:3,1:3,n_atoms))
     if (.not.allocated(HF_stress_per_atom_vdw_change)) allocate(HF_stress_per_atom_vdw_change  (1:3,1:3,n_atoms))
     HF_stress_per_atom_vdw(:,:,:) = 0.0d0
     HF_stress_per_atom_vdw_change(:,:,:) = 0.0d0
   end if
end if

if (.not.hirsh_analysis) then
       write(info_str,'(1X,A)') & 
    "Since the Hirshfeld analysis is not performed, the vdW energy is set to zero"
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
   return
endif


! Damping function parameter
d=20.d0
! Scaling factor for the onset of the vdW correction
select case (flag_xc) 
 case (1) !PBE0
  Sr=0.96
 case (6) ! PBE
  Sr=0.94
 case (7) !HSE
  Sr=0.96
 case (9) ! BLYP 
  Sr=0.62
 case (10) ! B3LYP 
  Sr=0.84
 case (12) ! revPBE
  Sr=0.60
 case (17) ! PBEsol
  Sr = 1.13  ! value from email from Alex T. to Jan H. from 2017-01-25
 case (20) ! AM05
  Sr=0.84
! Following values taken from DOI: 10.1021/ct2005616 
 case (25) ! M06-L
  Sr=1.27
 case (27) ! M06
  Sr=1.16
 case (51) ! TPSS
  Sr=0.86
 case default
  Sr=1.0

 select case (flag_post_xc)
  case (1) ! M06L
    Sr=1.27
  case (2) ! M06
    Sr=1.16
  case (6)
    Sr=0.86
  case default
 end select 
end select
 
if (n_periodic == 0) then
 vdw_n1=0
 vdw_n2=0
 vdw_n3=0
end if

! first: obtain explicit parameters for all the species BEFORE starting the vdw calculation
! that leaves them open to fudging later, within sort of one single big loop.
alpha1(:) = 0d0
R01(:)    = 0d0
C61(:)    = 0d0
do i_atom1 = 1, n_atoms
   if (.not.vdw_hirshfeld_data_external(species(i_atom1))) then
      atom=species_element(species(i_atom1))
      nucleus=species_z(species(i_atom1))
      call get_vdw_param(atom,nucleus,C6_array,alpha_array,R0_array)
   else
      C6_array    = vdw_hirshfeld_C6   (species(i_atom1))
      alpha_array = vdw_hirshfeld_alpha(species(i_atom1))
      R0_array    = vdw_hirshfeld_R0   (species(i_atom1))
   end if
   ! precompute per-atom quantities
   alpha1(i_atom1) = alpha_array*hirshfeld_volume(i_atom1)
   C61(i_atom1) = C6_array*(hirshfeld_volume(i_atom1)*hirshfeld_volume(i_atom1))
   R01(i_atom1) = R0_array*hirshfeld_volume(i_atom1)**0.333333333333333333333333
   ! safety net that should not be necessary but replaces older "if" inside double loop
   if (C61(i_atom1) .le. 0.d0) then
     C61(i_atom1) = 1.d0
   end if
   if (alpha1(i_atom1) .le. 0.d0) then
     alpha1(i_atom1) = 1.d0
   end if
   if (R01(i_atom1) .le. 0.d0) then
     R01(i_atom1) = 1.d0
   end if
!   if (isnan(R01(i_atom1))) then
!   isnan doesn't compile with IBM XLF
!   nan can be checked as follows, since nan is not equal to anything, even itself.
   if (R01(i_atom1) /= R01(i_atom1)) then
     R01(i_atom1) = 1.d0
   end if

end do

if(present(vdw_Hessian)) then
   vdW_Hessian = 0d0
endif

! calculation of all atoms once, this gives the entire vdw interaction for cluster systems:
do i_atom1 = 1, n_atoms
   if (myid.eq.task_list(i_atom1)) then
      do i_atom2 =  1, n_atoms
         C6=2.0*C61(i_atom1)*C61(i_atom2)*vdw_C6_prefactor(species(i_atom1),species(i_atom2)) & 
            / (alpha1(i_atom2)/alpha1(i_atom1)*C61(i_atom1)+alpha1(i_atom1)/alpha1(i_atom2)*C61(i_atom2))
         R0=Sr*(R01(i_atom1)+R01(i_atom2))          
         coord_diff(:) = coords(:,i_atom1) - coords(:,i_atom2)
         dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
         if (dist2 .gt. 1e-6 .and. .not. empty(i_atom1) .and. .not. empty(i_atom2)) then
            dist  = dsqrt(dist2)
            dist6 = dist2*dist2*dist2 
            rpow=d*(dist/R0-1.0D0)
            rexp=dexp(-rpow)
            aexp=1.0d0+rexp
            rphi = -C6/(dist6*aexp)
            vdw_energy = vdw_energy + rphi 
            !  vdW forces
            if (use_forces) then
               rfac=-6.0d0/dist + d*rexp/(aexp*R0)
               ffac=rphi*rfac/dist
               vdw_forces(:,i_atom1) = vdw_forces(:,i_atom1) - coord_diff(:)*ffac
               if (use_analytical_stress) then
                  do i_coord = 1, 3, 1
                     do i_coord2 = 1, 3, 1
                        AS_vdw_stress(i_coord,i_coord2) = AS_vdw_stress(i_coord,i_coord2) & 
                        + ffac * coord_diff(i_coord) * coord_diff(i_coord2)
                     end do
                  end do
                  if (compute_heat_flux) then
                    do i_coord = 1, 3, 1
                       do i_coord2 = 1, 3, 1
                          HF_stress_per_atom_vdw(i_coord,i_coord2,i_atom1) = HF_stress_per_atom_vdw(i_coord,i_coord2,i_atom1) &
                          + ffac * coord_diff(i_coord) * coord_diff(i_coord2)
                          HF_stress_per_atom_vdw(i_coord,i_coord2,i_atom2) = HF_stress_per_atom_vdw(i_coord,i_coord2,i_atom2) &
                          + ffac * coord_diff(i_coord) * coord_diff(i_coord2)
                       end do
                    end do
                  end if
               end if !use_analytical_stress
            end if !use_forces
			if (use_DFPT_phonon_reduce_memory .and. present(vdW_Hessian) .and. use_forces) then !vdW correction requested for frequencies
                                !OTH: This code is currently untested, because there
                                !are still problems with the underlying DFPT.
                                !I don't want to loose it, but I don't want
                                !peopole to use it either. Hence, I'm putting in
                                !a stop for now
                                call aims_stop('DFPT + vdW is currently untested. Please contact the developers if necessary. ')
				dR0=d/R0
				dR0raexp = dR0*rexp/(aexp)
				dEdr=rphi*rfac
				d2Edr2= rphi*dR0raexp*(-12/dist+2*dR0raexp-dR0)-42*C6/(dist6*dist2*aexp)
				if (i_atom1 .eq. i_atom2) then
				     signum = 1.0
				else 
				     signum = -1.0
				endif
				do i_coord = 1,3,1
					do i_coord2 = 1,3,1
						dRdx1 = coord_diff(i_coord)/dist
						dRdx2 = coord_diff(i_coord2)/dist
						d2Rdx1dx2 = -dRdx1*dRdx2
						if (i_coord .eq. i_coord2) then
							d2Rdx1dx2 = d2Rdx1dx2+1						
						endif
						d2Rdx1dx2 = d2Rdx1dx2 / dist
						
						vdW_Hessian(3*(i_atom1-1)+i_coord,3*(i_atom2-1)+i_coord2) = vdW_Hessian(3*(i_atom1-1)+i_coord,3*(i_atom2-1)+i_coord2) + signum*(d2Edr2*dRdx1*dRdx2+dEdr*d2Rdx1dx2)
					enddo
				enddo
				
			
			endif !vdW correction for frequencies
         end if
      end do
   end if
end do
vdw_energy = 0.5*vdw_energy ! factor of 0.5 due to complete ij summation
call sync_vdw_correction(vdw_energy,vdw_forces)

if (use_analytical_stress) then
   if (use_mpi) then
      call sync_matrix(AS_vdw_stress, 3, 3)
      if (compute_heat_flux) then
        call sync_matrix( HF_stress_per_atom_vdw(1,:,:) , 3, n_atoms )
        call sync_matrix( HF_stress_per_atom_vdw(2,:,:) , 3, n_atoms )
        call sync_matrix( HF_stress_per_atom_vdw(3,:,:) , 3, n_atoms )
      end if
   end if
   ! FK: 1/2 due to complete summation over i and j 
   AS_vdw_stress(1:3,1:3) = 0.5d0 * AS_vdw_stress(1:3,1:3)
   if (compute_heat_flux) then
     HF_stress_per_atom_vdw(1:3,1:3,1:n_atoms) = 0.25d0 * HF_stress_per_atom_vdw(1:3,1:3,1:n_atoms)
   end if
end if

! VB: The loop below can become very expensive for large systems if a small value of vdw_convergence_threshold is chosen.
!     I addressed this issue partly on May 11, 2016, by modifying the way in which the default is chosen.
!     The code path below can likely be further optimized if needed, by choosing a more sophisticated supercell
!     convergence strategy, or by simply applying precomputed pairwise lattice sums 
!     (would be a much more elegant approach).

! do we need to converge the periodic results ?
if (n_periodic .gt. 0) then
   write(info_str,'(2X,A)') "Converging periodic supercell size for vdw calculation."
   call localorb_info(info_str)
   write(info_str,'(2X,A)') "Convergence accuracy (and time!) can be adjusted using the vdw_convergence_threshold keyword."
   call localorb_info(info_str)
   periodic_converged = .false.
   periodic_cell      = 0


   do while (.not.periodic_converged)
      periodic_cell = periodic_cell + 1 
      write(info_str,'(2X,A,I5)') "| Shell number: ", periodic_cell
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      vdw_energy_change = 0d0
      vdw_forces_change = 0d0


      if (use_analytical_stress) then
         AS_vdw_change(:,:) = 0.0d0
         if (compute_heat_flux) then
           HF_stress_per_atom_vdw_change(:,:,:) = 0.0d0
         end if
      end if

      do i1 = -periodic_cell, periodic_cell
         do i2 = -periodic_cell, periodic_cell
            do i3 = -periodic_cell, periodic_cell 
               ! make sure that only those cells in the "outer shell" of currently treated cells are included
               if ((abs(i1).eq.periodic_cell).or.(abs(i2).eq.periodic_cell).or.(abs(i3).eq.periodic_cell)) then
                  do i_atom1 = 1, n_atoms, 1
                     if (myid.eq.task_list(i_atom1)) then
                        do i_atom2 =  1, n_atoms, 1
                           C6=2.0*C61(i_atom1)*C61(i_atom2)*vdw_C6_prefactor(species(i_atom1),species(i_atom2)) & 
                             / (alpha1(i_atom2)/alpha1(i_atom1)*C61(i_atom1)+alpha1(i_atom1)/alpha1(i_atom2)*C61(i_atom2))
                           R0=Sr*(R01(i_atom1)+R01(i_atom2))          
                           coord_curr(:) = coords(:,i_atom2)+i1*lattice_vector(:,1)+i2*lattice_vector(:,2)+i3*lattice_vector(:,3) 
                           coord_diff(:) = coords(:,i_atom1) - coord_curr(:)
                           dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
                           dist6 = dist2*dist2*dist2 
                           dist  = sqrt(dist2)
                           rpow=d*(dist/R0-1.0D0)
                           rexp=dexp(-rpow)
                           aexp=1.0d0+rexp
                           rphi = -C6/(dist6*aexp)
                           vdw_energy_change = vdw_energy_change + rphi 
                           if (use_forces) then
                              rfac=-6.0d0/dist + d*rexp/(aexp*R0)
                              ffac=rphi*rfac/dist
                              vdw_forces_change(:,i_atom1) = vdw_forces_change(:,i_atom1) - coord_diff(:)*ffac

                              if (use_analytical_stress) then
                                 do i_coord = 1, 3, 1
                                    do i_coord2 = 1, 3, 1
                                       AS_vdw_change(i_coord,i_coord2) = &
                                          AS_vdw_change(i_coord,i_coord2) &
                                          + ffac &
                                          * coord_diff(i_coord) &
                                          * coord_diff(i_coord2)
                                    end do
                                 end do
                                 if (compute_heat_flux) then
                                   do i_coord = 1, 3, 1
                                      do i_coord2 = 1, 3, 1
                                         HF_stress_per_atom_vdw_change(i_coord,i_coord2,i_atom1) = &
                                         HF_stress_per_atom_vdw_change(i_coord,i_coord2,i_atom1)   &
                                         + ffac * coord_diff(i_coord) * coord_diff(i_coord2)
                                         HF_stress_per_atom_vdw_change(i_coord,i_coord2,i_atom2) = & 
                                         HF_stress_per_atom_vdw_change(i_coord,i_coord2,i_atom2) &
                                         + ffac * coord_diff(i_coord) * coord_diff(i_coord2)
                                      end do
                                   end do
                                 end if
                              end if !use_analytical_stress
                           end if !use_forces
						   if (use_DFPT_phonon_reduce_memory .and. present(vdW_Hessian) .and. use_forces) then !vdW correction requested for frequencies
							    if (i_atom1 .eq. i_atom2) cycle 	
								dR0=d/R0
								dR0raexp = dR0*rexp/(aexp)
								dEdr=rphi*rfac
								d2Edr2= rphi*dR0raexp*(-12/dist+2*dR0raexp-dR0)-42*C6/(dist6*dist2*aexp)
								signum = -1.0
								do i_coord = 1,3,1
									do i_coord2 = 1,3,1
										dRdx1 = coord_diff(i_coord)/dist
										dRdx2 = coord_diff(i_coord2)/dist
										d2Rdx1dx2 = -dRdx1*dRdx2
										if (i_coord .eq. i_coord2) then
											d2Rdx1dx2 = d2Rdx1dx2+1						
										endif
										d2Rdx1dx2 = d2Rdx1dx2 / dist
										
										vdW_Hessian(3*(i_atom1-1)+i_coord,3*(i_atom2-1)+i_coord2) = vdW_Hessian(3*(i_atom1-1)+i_coord,3*(i_atom2-1)+i_coord2) + signum*(d2Edr2*dRdx1*dRdx2+dEdr*d2Rdx1dx2)
									enddo
								enddo
							endif !vdW correction for frequencies
                        end do
                     end if
                  end do
               end if
            end do
         end do
      end do
      ! synchronize only the change ... 
      vdw_energy_change = 0.5*vdw_energy_change ! factor of 0.5 due to complete ij summation
      call sync_vdw_correction(vdw_energy_change,vdw_forces_change)
      vdw_energy = vdw_energy + vdw_energy_change

      ! update forces if necessary 
      if (use_forces) then
         vdw_forces(:,:) = vdw_forces(:,:) + vdw_forces_change(:,:)
         if (use_analytical_stress) then
            if (use_mpi) then
               call sync_matrix(AS_vdw_change, 3, 3)
               if (compute_heat_flux) then
                 call sync_matrix( HF_stress_per_atom_vdw_change(1,:,:) , 3, n_atoms )
                 call sync_matrix( HF_stress_per_atom_vdw_change(2,:,:) , 3, n_atoms )
                 call sync_matrix( HF_stress_per_atom_vdw_change(3,:,:) , 3, n_atoms )
               end if
            end if
            ! FK: 1/2 due to complete summation over i and j 
            AS_vdw_stress(1:3,1:3) = AS_vdw_stress(1:3,1:3) + 0.5 * AS_vdw_change(1:3,1:3)
            if (compute_heat_flux) then
              HF_stress_per_atom_vdw(1:3,1:3,1:n_atoms) = HF_stress_per_atom_vdw(1:3,1:3,1:n_atoms) + & 
                0.25d0 * HF_stress_per_atom_vdw_change(1:3,1:3,1:n_atoms)
            end if
         end if !use_analytical_stress
      end if !use_forces

      ! check convergence
      periodic_converged = (dabs(vdw_energy_change).le.vdw_convergence_threshold)

      write(info_str,'(2X,A,E14.6,A)') "|   VDW energy change: ", vdw_energy_change*hartree, & 
           " eV."
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      
      ! Safety check and stop calculation if something went wrong
!      if (isnan(vdw_energy_change)) then
      if (vdw_energy_change /= vdw_energy_change) then
         write(use_unit,*) '*** Error: vdw_energy_change = ', vdw_energy_change
         call aims_stop("NaN","calc_vdw_hirshfeld in vdw_correction")
      end if

      ! check force convergence if necessary (note that stress is not tested)
      if (periodic_converged.and.use_forces) then
         if (sc_accuracy_forces.gt.0d0) then
            ! no need to converge forces if the accuracy has been set to a negative value & checking was turned off explicitly
            maxchange = 0d0
            do i_atom1 = 1, n_atoms
               do i_coord = 1, 3
                  maxchange = max(maxchange,dabs(vdw_forces_change(i_coord,i_atom1)))
               end do
            end do
            maxchange = maxchange*hartree/bohr
            periodic_converged = periodic_converged.and.(maxchange.lt.sc_accuracy_forces)            
            write(info_str,'(2X,A,E14.6,A)') "|   max. VDW force change: ", maxchange, " eV/AA."
            call localorb_info(info_str, use_unit, '(A)', OL_norm)
         end if ! sc_accuracy_forces
      end if !use_forces

      if (periodic_converged) then
         write(info_str,'(2X,A)') "|   Converged. "
         call localorb_info(info_str, use_unit, '(A)', OL_norm)
      else
         write(info_str,'(2X,A)') "|   Not converged. "
         call localorb_info(info_str, use_unit, '(A)', OL_norm)
      end if

   end do   ! end periodic convergence loop

   if (compute_heat_flux) then
     if (allocated(HF_stress_per_atom_vdw_change))        deallocate(HF_stress_per_atom_vdw_change)
   end if

end if  ! periodic convergence 

end subroutine calc_vdw_hirshfeld





!****************************************************************
!****************************************************************
!*************   SELF CONSISTENT ROUTINES   *********************
!****************************************************************
!****************************************************************




!****************************************************************
!*************   BEGIN ENERGY ROUTINE   *************************
!****************************************************************
subroutine en_vdw_hirshfeld(vdw_energy)
use dimensions
use runtime_choices
use geometry
use species_data
use localorb_io
use mpi_utilities
use synchronize_mpi
use constants

implicit none

! Imported variables
real*8 :: vdw_energy

! Local Variables
real*8, dimension(3) :: coord_diff
real*8, dimension(3) :: coord_curr
real*8 :: vdw_energy_change
real*8 :: dist
real*8 :: dist2
real*8 :: dist6
real*8 :: rpow
real*8 :: rexp
real*8 :: aexp
real*8 :: rphi
real*8 :: C6
real*8 :: R0
real*8 :: d
real*8 :: Sr
logical :: periodic_converged
integer :: periodic_cell
real*8 :: nucleus

real*8 :: alpha_array, R0_array, C6_array
real*8, dimension(n_atoms) :: alpha1,C61, R01 

character*160 :: info_str
character*2 :: atom

!counters
integer :: i_atom1
integer :: i_atom2
integer :: i_coord
integer :: i1
integer :: i2
integer :: i3
integer :: param

! The standard TS energy (written above) is computed in the post-process fashion with a call for the hirshfeld 
! volumes before it. Here we need to call the hirshfeld analysis each time we compute the energy and the 
! potential. Moreover, we also need the free volumes and the hirshfeld weigts. The simplest way is to insert 
! the call here.
write(info_str,'(2X,A)') &
     "Evaluating non-empirical van der Waals self-consistent correction (Tkatchenko/Scheffler 2009)."
call localorb_info ( info_str, use_unit,'(A)', OL_norm )
call hirshfeld_analysis( )


vdw_energy = 0.d0

if (.not.hirsh_analysis) then
       write(info_str,'(1X,A)') & 
    "Since the Hirshfeld analysis is not performed, the vdW energy is set to zero"
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
   return
endif

! Damping function parameter
d=20.d0
! Scaling factor for the onset of the vdW correction
select case (flag_xc) 
 case (1) !PBE0
  Sr=0.96 
 case (6) ! PBE
  Sr=0.94
 case (7) !HSE
  Sr=0.96
 case (9) ! BLYP 
  Sr=0.62
 case (10) ! B3LYP 
  Sr=0.84
 case (12) ! revPBE
  Sr=0.60
 case (17) ! PBEsol
  Sr = 1.13  ! value from email from Alex T. to Jan H. from 2017-01-25
 case (20) ! AM05
  Sr=0.84
! Following values taken from DOI: 10.1021/ct2005616 
 case (25) ! M06-L
  Sr=1.27
 case (27) ! M06
  Sr=1.16
 case (51) ! TPSS
  Sr=0.86
 case default
  Sr=1.0 

 select case (flag_post_xc)
  case (1) ! M06L
    Sr=1.27
  case (2) ! M06
    Sr=1.16
  case (6)
    Sr=0.86
  case default
 end select 
end select
 
if (n_periodic == 0) then
 vdw_n1=0
 vdw_n2=0
 vdw_n3=0
end if

! first: obtain explicit parameters for all the species BEFORE starting the vdw calculation
! that leaves them open to fudging later, within sort of one single big loop.
alpha1(:) = 0d0
R01(:)    = 0d0
C61(:)    = 0d0
do i_atom1 = 1, n_atoms
   if (.not.vdw_hirshfeld_data_external(species(i_atom1))) then
      atom=species_name(species(i_atom1))
      nucleus=species_z(species(i_atom1))
      call get_vdw_param(atom,nucleus,C6_array,alpha_array,R0_array)
   else
      C6_array    = vdw_hirshfeld_C6   (species(i_atom1))
      alpha_array = vdw_hirshfeld_alpha(species(i_atom1))
      R0_array    = vdw_hirshfeld_R0   (species(i_atom1))
   end if
   ! precompute per-atom quantities
   alpha1(i_atom1) = alpha_array*hirshfeld_volume(i_atom1)
   C61(i_atom1) = C6_array*(hirshfeld_volume(i_atom1)*hirshfeld_volume(i_atom1))
   R01(i_atom1) = R0_array*hirshfeld_volume(i_atom1)**0.333333333333333333333333
   ! safety net that should not be necessary but replaces older "if" inside double loop
   if (C61(i_atom1) .le. 0.d0) then
     C61(i_atom1) = 1.d0
   end if
   if (alpha1(i_atom1) .le. 0.d0) then
     alpha1(i_atom1) = 1.d0
   end if
   if (R01(i_atom1) .le. 0.d0) then
     R01(i_atom1) = 1.d0
   end if
!   if (isnan(R01(i_atom1))) then
   if (R01(i_atom1) /= R01(i_atom1)) then
     R01(i_atom1) = 1.d0
   end if
end do

! calculation of all atoms once, this gives the entire vdw interaction for cluster systems:
do i_atom1 = 1, n_atoms
   if (first_integration) then
      param = 0
   else
      param = task_list(i_atom1)
   end if
   if (myid.eq.param) then
      do i_atom2 =  1, n_atoms
         C6=2.0*C61(i_atom1)*C61(i_atom2)*vdw_C6_prefactor(species(i_atom1),species(i_atom2)) & 
            / (alpha1(i_atom2)/alpha1(i_atom1)*C61(i_atom1)+alpha1(i_atom1)/alpha1(i_atom2)*C61(i_atom2))
         R0=Sr*(R01(i_atom1)+R01(i_atom2))          
         coord_diff(:) = coords(:,i_atom1) - coords(:,i_atom2)
         dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
         if (dist2 .gt. 1e-6 .and. .not. empty(i_atom1) .and. .not. empty(i_atom2)) then
            dist  = dsqrt(dist2)
            dist6 = dist2*dist2*dist2 
            rpow=d*(dist/R0-1.0D0)
            rexp=dexp(-rpow)
            aexp=1.0d0+rexp
            rphi = -C6/(dist6*aexp)
            vdw_energy = vdw_energy + rphi 
         end if
      end do
   end if
end do
!default value if we use a non-periodic structure
lattice_translations=0
!calculating the total energy
vdw_energy = 0.5*vdw_energy ! factor of 0.5 due to complete ij summation
call sync_real_number(vdw_energy)


! do we need to converge the periodic results ?
if (n_periodic .gt. 0) then
   write(info_str,'(2X,A)') "Converging periodic supercell size for vdw calculation."
   call localorb_info(info_str)
   periodic_converged = .false.
   periodic_cell      = 0
   do while (.not.periodic_converged)
      periodic_cell = periodic_cell + 1 
      vdw_energy_change = 0d0
      do i1 = -periodic_cell, periodic_cell
         do i2 = -periodic_cell, periodic_cell
            do i3 = -periodic_cell, periodic_cell 
               ! make sure that only those cells in the "outer shell" of currently treated cells are included
               if ((abs(i1).eq.periodic_cell).or.(abs(i2).eq.periodic_cell).or.(abs(i3).eq.periodic_cell)) then
                  do i_atom1 = 1, n_atoms, 1
                     if (first_integration) then
                        param = 0
                     else
                        param = task_list(i_atom1)
                     end if
                     if (myid.eq.param) then
                        do i_atom2 =  1, n_atoms, 1
                           C6=2.0*C61(i_atom1)*C61(i_atom2)*vdw_C6_prefactor(species(i_atom1),species(i_atom2)) & 
                             / (alpha1(i_atom2)/alpha1(i_atom1)*C61(i_atom1)+alpha1(i_atom1)/alpha1(i_atom2)*C61(i_atom2))
                           R0=Sr*(R01(i_atom1)+R01(i_atom2))          
                           coord_curr(:) = coords(:,i_atom2)+i1*lattice_vector(:,1)+i2*lattice_vector(:,2)+i3*lattice_vector(:,3) 
                           coord_diff(:) = coords(:,i_atom1) - coord_curr(:)
                           dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
                           dist6 = dist2*dist2*dist2 
                           dist  = sqrt(dist2)
                           rpow=d*(dist/R0-1.0D0)
                           rexp=dexp(-rpow)
                           aexp=1.0d0+rexp
                           rphi = -C6/(dist6*aexp)
                           vdw_energy_change = vdw_energy_change + rphi 
                        end do
                     end if
                  end do
               end if
            end do
         end do
      end do
      ! synchronize only the change ... 
      vdw_energy_change = 0.5*vdw_energy_change ! factor of 0.5 due to complete ij summation
      call sync_real_number(vdw_energy_change)
      periodic_converged = (dabs(vdw_energy_change).le.vdw_convergence_threshold)
      vdw_energy = vdw_energy + vdw_energy_change
   end do   ! end periodic convergence loop
   !store the number of replicas needed to converge the energy, we'll use it in the potential routine
   lattice_translations = periodic_cell
end if  ! periodic convergence 

end subroutine en_vdw_hirshfeld
!*************   END ENERGY ROUTINE   ***************************


!****************************************************************
!*************   BEGIN POTENTIAL ROUTINE   **********************
!****************************************************************
subroutine calc_vdw_potential(vdw_potential)

use dimensions
use geometry
use species_data
use localorb_io
use mpi_utilities
use synchronize_mpi
use runtime_choices
use grids
use pbc_lists
use free_atoms
use spline

implicit none

real*8, dimension(:), allocatable :: vdw_potential

! Local Variables
real*8, dimension(3) :: coord_diff
real*8, dimension(3) :: coord_curr
real*8, dimension(3) :: coord_add
real*8, dimension(3) :: coord_b

real*8 :: vdw_energy_change
real*8 :: dist
real*8 :: dist2
real*8 :: dist6
real*8 :: rpow
real*8 :: rexp
real*8 :: aexp
real*8 :: rphi
real*8 :: rfac
real*8 :: ffac
real*8 :: C6
real*8 :: R0
real*8 :: d
real*8 :: Sr
real*8 :: C61
real*8 :: alpha1
real*8 :: R01
real*8 :: C62
real*8 :: alpha2
real*8 :: R02
logical :: periodic_converged
integer :: periodic_cell

integer :: i_full_points
integer :: i_my_batch
integer :: i_index
real*8, dimension(3) :: coord_current
logical :: point_on_atom
real*8 dist_tab_sq(n_centers_basis_integrals)
real*8 dist_tab(n_centers_basis_integrals)
real*8 dir_tab(3,n_centers_basis_integrals)
real*8 dir_tab_norm(3,n_centers_basis_integrals)
real*8, dimension(n_species) :: r_grid_min_sq
real*8 :: partition_norm
real*8 :: aux_dens
real*8 :: aux_dens1
real*8 :: aux_dens2
real*8 :: full_dens
integer :: i_center_L
integer :: i_center
real*8 :: radial1
real*8 :: radial2
real*8 i_r(n_centers_basis_integrals)
real*8, dimension(:,:,:), allocatable :: reference_rho_spl
real*8, dimension(:,:,:,:), allocatable :: potential_table
real*8, dimension(:), allocatable :: potential_table_atom
real*8 :: distance1
real*8 :: weight1
integer atom2
real*8 :: distance2
real*8 :: weight2
real*8 :: totalweight1
real*8 :: totalweight2
real*8 :: C6free
real*8 :: fdamp_der
real*8 :: fdamp_der_a
real*8 :: fdamp_der_b
real*8 :: r1_der
real*8 :: r2_der
real*8 :: couple_potential
real*8 :: couple_potential_a
real*8 :: couple_potential_b
real*8 :: grid_potential
real*8, dimension(n_atoms) :: alpha_array, R0_array, C6_array
real*8 :: cut_fac
real*8 i_r1
real*8 :: nucleus

character*160 :: info_str
character*2 :: atom
real*8 :: current_i_r
!counters
integer :: i_pair
integer :: i_atom1
integer :: i_atom2
integer :: i_coord
integer :: i1
integer :: i2
integer :: i3
integer :: current_atom
integer :: current_radial
integer :: current_angular 
integer :: param
integer :: i_p1, i_p2, i_p3
integer :: trasl_index1, trasl_index2, trasl_index3

!This cutoff involves only the replica of the derivative of the hirshfeld volume
integer, dimension(3) :: smart_cut
real*8 :: cutoff

! Damping function parameter
d=20.d0
! Scaling factor for the onset of the vdW correction
select case (flag_xc) 
 case (1) !PBE0
  Sr=0.96 
 case (6) ! PBE
  Sr=0.94
 case (7) !HSE
  Sr=0.96 
 case (9) ! BLYP 
  Sr=0.62
 case (10) ! B3LYP 
  Sr=0.84
 case (12) ! revPBE
  Sr=0.60
 case (17) ! PBEsol
  Sr = 1.13  ! value from email from Alex T. to Jan H. from 2017-01-25
 case (20) ! AM05
  Sr=0.84
! Following values taken from DOI: 10.1021/ct2005616 
 case (25) ! M06-L
  Sr=1.27
 case (27) ! M06
  Sr=1.16
 case (51) ! TPSS
  Sr=0.86
 case default
  Sr=1.0 

 select case (flag_post_xc)
  case (1) ! M06L
    Sr=1.27
  case (2) ! M06
    Sr=1.16
  case (6)
    Sr=0.86
  case default
 end select 
end select
 
if (n_periodic == 0) then
 vdw_n1=0
 vdw_n2=0
 vdw_n3=0
end if

vdw_potential = 0.d0

if (.not.hirsh_analysis) then
       write(info_str,'(1X,A)') & 
    "Since the Hirshfeld analysis is not performed, the vdW potential is set to zero"
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
   return
endif

! first: obtain explicit parameters for all the species BEFORE starting the vdw calculation
! that leaves them open to fudging later, within sort of one single big loop.
alpha_array(:) = 0d0
R0_array(:)    = 0d0
C6_array(:)    = 0d0

do i_atom1 = 1, n_atoms
   if (.not.vdw_hirshfeld_data_external(species(i_atom1))) then
      atom=species_name(species(i_atom1))
      nucleus=species_z(species(i_atom1))
      call get_vdw_param(atom,nucleus,C6_array(i_atom1),alpha_array(i_atom1),R0_array(i_atom1))
   else
      C6_array(i_atom1)    = vdw_hirshfeld_C6   (species(i_atom1))
      alpha_array(i_atom1) = vdw_hirshfeld_alpha(species(i_atom1))
      R0_array(i_atom1)    = vdw_hirshfeld_R0   (species(i_atom1))
   end if
end do

if(.not. use_distributed_spline_storage) then
   allocate (reference_rho_spl(n_max_spline,n_max_grid, n_atoms))
   do i_atom1 = 1, n_atoms, 1
      reference_rho_spl(:,:,i_atom1) = free_rho_spl(:,:,species(i_atom1))
   end do
end if

!In case of NON-periodic system one needs only 1,1,1 lattice dimensions, and here i'm allocating 2,2,2 cause it's the easiest way valid for both cases
!There is an extra case during the periodic cycle: when all the indexes meet 0 -> original unit cell. In this case i just store 0 in the table, but it's an extra position that's the reason for +2 instead of +1 below

!Periodic case
if (n_periodic .gt. 0) then
   !The cutoff below corresponds to 10A, the limit that we take to replicate the derivative of the hirshfeld volumes
   cutoff=18.897259886
   do  i1=1, 3, 1
      smart_cut(i1)=cutoff/sqrt(lattice_vector(1,i1)**2+lattice_vector(2,i1)**2+lattice_vector(3,i1)**2)
      smart_cut(i1)=smart_cut(i1)+1
   end do
   periodic_cell = lattice_translations
   
   
   allocate (potential_table(n_atoms,2*smart_cut(1)+1,2*smart_cut(2)+1,2*smart_cut(3)+1))
   
   potential_table=0.d0
   
   
   !Store non-periodic elements
   do i_atom1 = 1, n_atoms, 1
      totalweight1=hirshfeldw(i_atom1)/freeintegral(i_atom1)
      alpha1=alpha_array(i_atom1)*totalweight1
      C61=C6_array(i_atom1)*(totalweight1*totalweight1)
      R01=R0_array(i_atom1)*(totalweight1**0.333333333333333333333333)
      do i_atom2 = 1, n_atoms, 1
         totalweight2=hirshfeldw(i_atom2)/freeintegral(i_atom2)
         alpha2=alpha_array(i_atom2)*totalweight2
         C62=C6_array(i_atom2)*(totalweight2*totalweight2)
         R02=R0_array(i_atom2)*(totalweight2**0.333333333333333333333333)
         C6free=2.0*C6_array(i_atom1)*C6_array(i_atom2)*vdw_C6_prefactor(species(i_atom1),species(i_atom2))
         C6free=C6free/((alpha_array(i_atom2)/alpha_array(i_atom1))*C6_array(i_atom1)+(alpha_array(i_atom1)/alpha_array(i_atom2))*C6_array(i_atom2))
         if ( (C61 .gt. 0.d0) .or. (C62 .gt. 0.d0) ) then  
            C6=2.0*C61*C62*vdw_C6_prefactor(species(i_atom1),species(i_atom2))/(alpha2/alpha1*C61+alpha1/alpha2*C62)
         else 
            C6 = 0.d0
         end if
         R0=(R01+R02)
         coord_diff(:) = coords(:,i_atom1) - coords(:,i_atom2)
         dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
         
         !Store periodic elements
         trasl_index1=0
         do i1 = -periodic_cell, periodic_cell, 1
            trasl_index1 = trasl_index1 +1
            trasl_index2 = 0
            do i2 = -periodic_cell, periodic_cell, 1
               trasl_index2 = trasl_index2 +1
               trasl_index3 = 0
               do i3 = -periodic_cell, periodic_cell, 1
                  trasl_index3 = trasl_index3 +1
                  ! make sure that it treats everything but the original unit cell
                  !Displacement in case of periodic structure
                  couple_potential_a= 0.d0      
                  fdamp_der_a = 0.d0
                  fdamp_der_b= 0.d0
                  couple_potential = 0.d0
                  fdamp_der = 0.d0
                  couple_potential_b= 0.d0
                  coord_add(:) = i1*lattice_vector(:,1)+i2*lattice_vector(:,2)+i3*lattice_vector(:,3)
                  !only atom2 displaced
                  coord_curr(:) = coords(:,i_atom2) + coord_add(:)
                  coord_diff(:) = coords(:,i_atom1) - coord_curr(:)
                  dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
                  if (dist2 .gt. 1e-6 .and. .not. empty(i_atom1) .and. .not. empty(i_atom2)) then
                     dist6 = dist2*dist2*dist2
                     dist  = dsqrt(dist2)
                     rpow=d*(dist/(Sr*R0)-1.0D0)
                     rexp=dexp(-rpow)
                     aexp=1.0d0+rexp
                     couple_potential = (C6free/(dist6*aexp))
                     couple_potential_a = couple_potential*totalweight2/freeintegral(i_atom1)
                     couple_potential_b = couple_potential*totalweight1/freeintegral(i_atom2)
                     fdamp_der = ((-1.d0*rexp*d*dist)/(Sr*R0*R0*aexp*aexp))*(C6/dist6)
                     fdamp_der_a = fdamp_der*R01/(3.0d0*hirshfeldw(i_atom1))
                     fdamp_der_b = fdamp_der*R02/(3.0d0*hirshfeldw(i_atom2))
                  endif
                  potential_table(i_atom1,smart_cut(1)+1,smart_cut(2)+1,smart_cut(3)+1) = potential_table(i_atom1,smart_cut(1)+1,smart_cut(2)+1,smart_cut(3)+1)+0.5*(couple_potential_a + fdamp_der_a)
                  if (abs(i1) .le. smart_cut(1) .and. abs(i2) .le. smart_cut(2)  .and. abs(i3) .le. smart_cut(3)) then
                     potential_table(i_atom2,i1+smart_cut(1)+1,i2+smart_cut(2)+1,i3+smart_cut(3)+1) = potential_table(i_atom2,i1+smart_cut(1)+1,i2+smart_cut(2)+1,i3+smart_cut(3)+1)+0.5*(couple_potential_b + fdamp_der_b)
                  endif
               end do ! periodicity 1
            end do ! 2
         end do ! 3
         !endif
      end do!end loop on atom2
   end do!end loop atom1
   
   
   !loop over the grid - r
   i_full_points = 0
   r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)
   do i_my_batch = 1, n_my_batches, 1
      do i_index = 1, batches(i_my_batch)%size, 1
         i_full_points = i_full_points + 1
         coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
         call tab_atom_centered_coords_p0 &
              ( coord_current, &
              dist_tab_sq, &
              dir_tab, &
              n_centers_basis_integrals, centers_basis_integrals )
         !Check this control, should be done for every couple of atoms not only one atom
         point_on_atom = .false.
         do i_center = 1, n_centers_basis_integrals, 1
            !Need to do exactly the same check as in the initialization of the partition tab: a point is defined to be "on an atom" when it is inside the innermost logarithmic grid shell of that atom.
            !At this point it is impossible to spline the density of the logarithmic grid.
            if ( dist_tab_sq(i_center).lt.r_grid_min_sq(species_center(centers_basis_integrals(i_center)))) then
               point_on_atom = .true.
               exit ! exit the loop
            end if
         end do
         if (.not.point_on_atom) then
            current_atom = batches(i_my_batch) % points(i_index) % index_atom
            current_radial = batches(i_my_batch) % points(i_index) % index_radial
            current_angular = batches(i_my_batch) % points(i_index) % index_angular
            call tab_global_geometry_p0 &
                 ( dist_tab_sq, &
                 dir_tab, &
                 dist_tab, &
                 i_r, &
                 dir_tab_norm, &
                 n_centers_basis_integrals, centers_basis_integrals )
            partition_norm = 0.d0
            if(.not. use_distributed_spline_storage) then
               do i_center_L = 1, n_centers_basis_integrals, 1
                  i_center = centers_basis_integrals(i_center_L)
                  aux_dens = val_spline &
                       ( i_r(i_center_L), reference_rho_spl(1,1,center_to_atom(i_center)), &
                    n_grid(species_center(i_center)) )         
                  partition_norm = partition_norm + aux_dens
               enddo
            else
               do i_center_L = 1, n_centers_basis_integrals, 1                
                  i_center = centers_basis_integrals(i_center_L)         
                  aux_dens = val_spline &
                       ( i_r(i_center_L), free_rho_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )       
                  partition_norm = partition_norm + aux_dens
               enddo
            endif
            
            grid_potential= 0.d0
            trasl_index1 = 0
            do i1 = -smart_cut(1), smart_cut(1), 1
               trasl_index1 = trasl_index1+1
               trasl_index2 = 0
               do i2 = -smart_cut(2), smart_cut(2), 1
                  trasl_index2 = trasl_index2+1
                  trasl_index3 = 0
                  do i3 = -smart_cut(3), smart_cut(3), 1
                     trasl_index3 = trasl_index3+1
                     ! make sure that it treats everything but the original unit cell
                     do i_atom2 = 1, n_atoms, 1
                        coord_add(:) = i1*lattice_vector(:,1)+i2*lattice_vector(:,2)+i3*lattice_vector(:,3)
                        coord_curr(:) = coord_current(:) - coord_add(:)
                        coord_b(:) = coord_curr(:) - coords(:,i_atom2)
                        distance2 = dsqrt(coord_b(1)**2 + coord_b(2)**2 + coord_b(3)**2)
                        call tab_i_r_p0( distance2, i_atom2, i_r1 ) !find i_r1 log coord for distance2
                        aux_dens2 = val_spline &   !find the density
                             ( i_r1, reference_rho_spl(1,1,center_to_atom(i_atom2)), &
                             n_grid(species_center(i_atom2)) )
                        weight2=aux_dens2*distance2*distance2*distance2
                        weight2=weight2/partition_norm
                        grid_potential = grid_potential + potential_table(i_atom2,i1+smart_cut(1)+1,i2+smart_cut(2)+1,i3+smart_cut(3)+1)*weight2
                     enddo
                  enddo
               enddo
            enddo
         endif ! point on atom
         vdw_potential(i_full_points) = -1.d0*grid_potential
      end do  ! points on batch
   end do ! batches

!Non periodic case
else
   allocate (potential_table_atom(n_atoms))
   potential_table_atom=0.d0
   do i_atom1 = 1, n_atoms, 1
      totalweight1=hirshfeldw(i_atom1)/freeintegral(i_atom1)
      alpha1=alpha_array(i_atom1)*totalweight1
      C61=C6_array(i_atom1)*(totalweight1*totalweight1)
      R01=R0_array(i_atom1)*(totalweight1**0.333333333333333333333333)
      do i_atom2 = 1, n_atoms, 1
         totalweight2=hirshfeldw(i_atom2)/freeintegral(i_atom2)
         alpha2=alpha_array(i_atom2)*totalweight2
         C62=C6_array(i_atom2)*(totalweight2*totalweight2)
         R02=R0_array(i_atom2)*(totalweight2**0.333333333333333333333333)
         C6free=2.0*C6_array(i_atom1)*C6_array(i_atom2)*vdw_C6_prefactor(species(i_atom1),species(i_atom2))
         C6free=C6free/((alpha_array(i_atom2)/alpha_array(i_atom1))*C6_array(i_atom1)+(alpha_array(i_atom1)/alpha_array(i_atom2))*C6_array(i_atom2))
         if ( (C61 .gt. 0.d0) .or. (C62 .gt. 0.d0) ) then
            C6=2.0*C61*C62*vdw_C6_prefactor(species(i_atom1),species(i_atom2))/(alpha2/alpha1*C61+alpha1/alpha2*C62)
         else 
            C6 = 0.d0
         end if
         R0=(R01+R02)
         coord_diff(:) = coords(:,i_atom1) - coords(:,i_atom2)
         dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
         couple_potential_a= 0.d0
         fdamp_der_a = 0.d0
         fdamp_der_b= 0.d0
         couple_potential = 0.d0
         fdamp_der = 0.d0
         couple_potential_b= 0.d0
         dist2 = coord_diff(1)*coord_diff(1)+coord_diff(2)*coord_diff(2)+coord_diff(3)*coord_diff(3)
         if (dist2 .gt. 1e-6 .and. .not. empty(i_atom1) .and. .not. empty(i_atom2)) then
            dist6 = dist2*dist2*dist2 
            dist  = dsqrt(dist2)
            rpow=d*(dist/(Sr*R0)-1.0D0)
            rexp=dexp(-rpow)
            aexp=1.0d0+rexp
            couple_potential = (C6free/(dist6*aexp))
            couple_potential_a = couple_potential*totalweight2/freeintegral(i_atom1)
            fdamp_der = ((-1.d0*rexp*d*dist)/(Sr*R0*R0*aexp*aexp))*(C6/dist6)
            fdamp_der_a = fdamp_der*R01/(3.0d0*hirshfeldw(i_atom1))
         endif
         potential_table_atom(i_atom1) = potential_table_atom(i_atom1)+couple_potential_a+fdamp_der_a
      enddo
   enddo
   
   !loop over the grid - r
   i_full_points = 0
   r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)
   do i_my_batch = 1, n_my_batches, 1
      do i_index = 1, batches(i_my_batch)%size, 1
         i_full_points = i_full_points + 1
         coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
         call tab_atom_centered_coords_p0 &
              ( coord_current, &
              dist_tab_sq, &
              dir_tab, &
              n_centers_basis_integrals, centers_basis_integrals )
         !Check this control, should be done for every couple of atoms not only one atom
         point_on_atom = .false.
         do i_center = 1, n_centers_basis_integrals, 1
            ! need to do exactly the same check as in the initialization of the partition tab: a point is defined to be "on an atom" when it is inside the innermost logarithmic grid shell of that atom.
            !At this point it is impossible to spline the density of the logarithmic grid.
            if ( dist_tab_sq(i_center).lt.r_grid_min_sq(species_center(centers_basis_integrals(i_center)))) then
               point_on_atom = .true.
               exit ! exit the loop
            end if
         end do
         if (.not.point_on_atom) then
            current_atom = batches(i_my_batch) % points(i_index) % index_atom
            current_radial = batches(i_my_batch) % points(i_index) % index_radial
            current_angular = batches(i_my_batch) % points(i_index) % index_angular
            call tab_global_geometry_p0 &
                 ( dist_tab_sq, &
                 dir_tab, &
                 dist_tab, &
                 i_r, &
                 dir_tab_norm, &
                 n_centers_basis_integrals, centers_basis_integrals )
            partition_norm = 0.d0
            if(.not. use_distributed_spline_storage) then
               do i_center_L = 1, n_centers_basis_integrals, 1
                  i_center = centers_basis_integrals(i_center_L)
                  aux_dens = val_spline &
                       ( i_r(i_center_L), reference_rho_spl(1,1,center_to_atom(i_center)), &
                       n_grid(species_center(i_center)) )
                  partition_norm = partition_norm + aux_dens
               enddo
            else
               do i_center_L = 1, n_centers_basis_integrals, 1
                  i_center = centers_basis_integrals(i_center_L)
                  aux_dens = val_spline &
                       ( i_r(i_center_L), free_rho_spl(1,1,species_center(i_center)), &
                       n_grid(species_center(i_center)) )
                  partition_norm = partition_norm + aux_dens
               enddo
            endif
            grid_potential= 0.d0
            do i_atom1 = 1, n_atoms, 1
               distance1 = dsqrt(dist_tab_sq(i_atom1))
               aux_dens1 = val_spline &
                    ( i_r(i_atom1), reference_rho_spl(1,1,center_to_atom(i_atom1)), &
                    n_grid(species_center(i_atom1)) )
               weight1 = abs(aux_dens1)
               weight1 = distance1*distance1*distance1*weight1
               weight1=weight1/partition_norm
               grid_potential = grid_potential+potential_table_atom(i_atom1)*weight1
            enddo
         endif
         vdw_potential(i_full_points) = -1.d0*grid_potential
      enddo
   enddo
endif

if(allocated(potential_table_atom)) deallocate(potential_table_atom)
if(allocated(potential_table)) deallocate(potential_table)



end subroutine calc_vdw_potential
!*************   END POTENTIAL ROUTINE   ************************


!****************************************************************
!*************   BEGIN HIRSHFELD ROUTINE   **********************
!****************************************************************
!This routine is used only inside scf_solver (and the related initialize, reinitialize) to update the hirshfeld quantities for the next step in the SC cycle
subroutine hirsh()

  use physics
  use runtime_choices, only: use_local_index, transport_calculation, &
      transport_lead_calculation, use_density_matrix
  use mpi_utilities
  use synchronize_mpi
  use dimensions, only: n_hamiltonian_matrix_size
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  use pbc_lists
  use species_data, only: l_shell_max

  implicit none
  character*160 :: info_str

  real*8, allocatable :: my_density_matrix(:)

     call get_n_max_for_densities( partition_tab, hartree_partition_tab)
     call kweight_occs('scf_solver:missing', occ_numbers)
     if (use_density_matrix) then
        if (use_local_index .or. transport_calculation .or. transport_lead_calculation) then
           call aims_allocate( my_density_matrix, n_hamiltonian_matrix_size, "my_density_matrix" )
           call update_missing_density_densmat &
                ( KS_eigenvector, KS_eigenvector_complex, &
                occ_numbers, partition_tab, hartree_partition_tab, rho, &
                l_shell_max, my_density_matrix )
           call aims_deallocate( my_density_matrix,                          "my_density_matrix" )
        else
           call update_missing_density_densmat &
                ( KS_eigenvector, KS_eigenvector_complex, &
                occ_numbers, partition_tab,  hartree_partition_tab, rho, &
                l_shell_max, hamiltonian(1,1) )
        end if
     else
        call update_missing_density_orbital &
             ( KS_eigenvector, KS_eigenvector_complex, &
             KS_eigenvalue, occ_numbers, &
             partition_tab, hartree_partition_tab, rho, &
             l_shell_max )
     end if
     call de_kweight_occs('scf_solver:missing', occ_numbers)
     
end subroutine hirsh
!*************   END HIRSHFELD ROUTINE   ************************


end module vdw_correction


