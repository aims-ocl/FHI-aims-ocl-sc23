!****h* FHI-aims/basis
!  NAME
!     basis
!  SYNOPSIS

 module basis

!  PURPOSE
!
!  Module basis contains all arrays related to the final (shrunken) basis
!  which is used to calculate the actual Hamiltonian.
!
!  Subroutines:
!  * allocate_basis
!  * cleanup_basis
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  SOURCE


      implicit none

!  global variable declarations

!     basis_wave_spl : Basis functions, cubic splined on the logarithmic grid
!     basis_wave_s_spl : For relativistic cases (X2C and 4C-DKS), basis_wave_spl saves the large component part, while basis_wave_s_spl saves the small comp part.
!     basis_deriv_spl : First derivative of basis functions, cubic splined on the logarithmic grid
!     basis_kinetic_spl : Part of kinetic operator on basis function, cubic splined on the logarithmic grid
!                         [e - v_basis(r)] * u_basis(r)
!                         which is also:  u"(r) + [l*(l+1)/r^2]*u(r)
!     basis_atom : atom to which a basis fn corresponds
!     basis_l    : angular momentum quantum number of given basis function
!     basis_m    : magnetic quantum number of given basis function
!     basis_fn   : The radial basis function number which corresponds to
!                  basis function i_basis [There are more basis functions (at,type,n,l,m)
!                  than there are necessary radial functions (at,type, n,l).]
!     basisfn_species : i_basis_fn -> i_species
!     basisfn_l       : i_basis_fn -> i_l
!     basisfn_type : Stores the type of a basis function:
!                  confined basis fn, fixed ionic basis fn, fixed hydrogenic basis fn
!     basisfn_n    : main quantum number of given basis function
!     outer_radius(i_basis_fn) : radius at which the basis function has a
!                                magnitude of wave_threshold (~ 1d-6).
!     outer_radius_sq(i_basis_fn) == outer_radius(i_basis_fn)**2
!     atom_radius(i_species): maximum outer_radius of basis_fn
!                                of that species
!        Explanation:
!        (formerly outer_radius_of_atoms) Radius beyond which all basis functions are exactly zero.
!
!        Please note that in the code it is always assumed that
!            multipole_radius >= multipole_radius_free >= atom_radius.
!        Any violation to this will lead to unexpected results.
!
!     atom_radius_sq(i_species): atom_radius(i_species)**2

!     Any quantities containing the "ext" label refer to the basis set created
!     by the additional basis functions added using the for_aux keyword.
!     This is the OBS+ basis set, which is described in and around Fig. 1 of 
!     New Journal of Physics 17, 093020 (2015)

!
!

      real*8, dimension(:,:,:), allocatable :: basis_wave_spl
      real*8, dimension(:,:,:), allocatable :: basis_wave_s_spl ! for fully_relativistic small component
      real*8, dimension(:,:,:), allocatable :: basis_deriv_spl
      real*8, dimension(:,:,:), allocatable :: basis_deriv_s_spl ! for small comp.
      real*8, dimension(:,:,:), allocatable :: basis_kinetic_spl
      real*8, dimension(:,:,:), allocatable :: basis_kinetic_s_spl ! for small comp.
      real*8, dimension(:,:,:), allocatable :: &
                                      basis_kinetic_scaled_zora_spl
                                      
      real*8, dimension(:,:,:), allocatable :: ext_wave_spl
      real*8, dimension(:,:,:), allocatable :: ext_deriv_spl
      real*8, dimension(:,:,:), allocatable :: ext_kinetic_spl
      real*8, dimension(:,:,:), allocatable :: &
                                      ext_kinetic_scaled_zora_spl
                                      
      integer, dimension(:), allocatable :: basis_atom
      integer, dimension(:), allocatable :: basis_small_atom ! for small component
      integer, dimension(:), allocatable :: basis_l
      integer, dimension(:), allocatable :: basis_small_l ! for small component
      integer, dimension(:), allocatable :: basis_k       ! needed only for relativistic large comp.
      integer, dimension(:), allocatable :: basis_small_k ! needed only for small component
      integer, dimension(:), allocatable :: basis_m
      integer, dimension(:), allocatable :: basis_small_m ! for small component
      integer, dimension(:), allocatable :: basis_fn
      integer, dimension(:), allocatable :: basis_small_fn ! for small component
      
      integer, dimension(:), allocatable :: ext_atom
      integer, dimension(:), allocatable :: ext_l
      integer, dimension(:), allocatable :: ext_m
      integer, dimension(:), allocatable :: ext_fn

      integer, dimension(:), allocatable :: basisfn_species
      integer, dimension(:), allocatable :: basisfn_l
      character*8, dimension(:), allocatable :: basisfn_type
      integer, dimension(:), allocatable :: basisfn_n
      integer, dimension(:), allocatable :: basisfn_k

      integer, dimension(:), allocatable :: extfn_species
      integer, dimension(:), allocatable :: extfn_l
      character*8, dimension(:), allocatable :: extfn_type
      integer, dimension(:), allocatable :: extfn_n


      real*8, dimension(:), allocatable :: outer_radius
      real*8, dimension(:), allocatable :: outer_radius_sq
      real*8, dimension(:), allocatable :: atom_radius
      real*8, dimension(:), allocatable :: atom_radius_sq

      real*8, dimension(:), allocatable :: outer_ext_radius
      real*8, dimension(:), allocatable :: outer_ext_radius_sq
      real*8, dimension(:), allocatable :: atom_ext_radius
      real*8, dimension(:), allocatable :: atom_ext_radius_sq


      integer, dimension(:), allocatable :: perm_basis_fns_spl
      integer, dimension(:), allocatable :: perm_ext_fns_spl
      integer, dimension(:), allocatable :: i_radial_fn
      integer, dimension(:), allocatable :: i_ext_radial_fn

      integer, dimension(:), allocatable :: basis_fn_start_spl
      integer, dimension(:), allocatable :: ext_fn_start_spl

      real*8, dimension(:,:,:), allocatable :: basis_wave_ordered
      real*8, dimension(:,:,:), allocatable :: basis_wave_s_ordered
      real*8, dimension(:,:,:), allocatable :: basis_deriv_ordered
      real*8, dimension(:,:,:), allocatable :: basis_deriv_s_ordered
      real*8, dimension(:,:,:), allocatable :: basis_kinetic_ordered
      real*8, dimension(:,:,:), allocatable :: basis_kinetic_s_ordered
      real*8, dimension(:,:,:), allocatable :: &
                                     basis_kinetic_scaled_zora_ordered

      real*8, dimension(:,:,:), allocatable :: ext_wave_ordered
      real*8, dimension(:,:,:), allocatable :: ext_deriv_ordered
      real*8, dimension(:,:,:), allocatable :: ext_kinetic_ordered
      real*8, dimension(:,:,:), allocatable :: &
                                     ext_kinetic_scaled_zora_ordered

      integer, dimension(:), allocatable :: n_basis_fn_species
      integer, dimension(:), allocatable :: n_basis_atom
      logical, dimension(:,:), allocatable :: basis_fn_atom

      integer, dimension(:), allocatable :: n_ext_fn_species
      integer, dimension(:), allocatable :: n_ext_atom
      logical, dimension(:,:), allocatable :: ext_fn_atom

      integer, dimension(:), allocatable :: spline_offset
      integer,    allocatable, dimension(:)   :: basis_mapping      

      integer, dimension(:), allocatable :: ext_spline_offset
      integer,    allocatable, dimension(:)   :: ext_mapping  

      real*8, allocatable, dimension(:) :: outradsq_save
      integer, allocatable, dimension(:) :: invlist_save

      ! --- JW: Arrays to loop basis functions by their angular momentum

      ! The consistency check loop at the end of generate_full_bas.f90
      ! exemplifies how these arrays are to be used.

      integer :: max_basis_L
      integer :: max_n_basis_fnLsp
      
      integer :: max_ext_L
      integer :: max_n_ext_fnLsp

      ! Group different radial fns by angular momentum L and i_species:

      ! L, i_species -> Number of radial parts in this L-channel
      integer, allocatable :: Lsp2n_basis_fnLsp(:,:)
      integer, allocatable :: Lsp2n_ext_fnLsp(:,:)
      ! i_fnLsp, L, i_species -> i_basis_fn
      integer, allocatable :: Lsp2basis_fn(:,:,:)
      integer, allocatable :: Lsp2ext_fn(:,:,:)

      ! The basis functions of a given atom make up a closed block in the list
      ! of all basis functions, and always have the same sequence for each
      ! atom of a given species.  Additionally, different Ms of otherwise
      ! identical quantum numbers follow in order.  Therefore, a given radial
      ! part uniquely defines the position within such a block.

      ! i_fnLsp, L, i_species -> i_basis_sp   [see atom2basis_off]
      integer, allocatable :: Lsp2basis_sp(:,:,:)
      integer, allocatable :: Lsp2ext_sp(:,:,:)
      ! i_basis == atom2basis_off(i_atom) + i_basis_sp + M [see Lsp2basis_sp]
      integer, allocatable :: atom2basis_off(:)
      integer, allocatable :: atom2ext_off(:)
      ! i_species -> Number of basis function per atom of this kind
      integer :: max_n_basis_sp
      integer :: max_n_basis_sp2 ! In case of atom splitting in calculate_fock
      integer, allocatable :: sp2n_basis_sp(:)

      integer :: max_n_ext_sp
      integer, allocatable :: sp2n_ext_sp(:)


      contains
!******
!---------------------------------------------------------------------
!****s* basis/allocate_basis
!  NAME
!    allocate_basis
!  SYNOPSIS

        subroutine allocate_basis( )

!  PURPOSE
!    Subroutine allocate_basis allocates all basis-related arrays, after
!    the relevant dimensions n_basis, n_basis_fn are known.    

!  USES

        use dimensions,      only : n_atoms, n_basis, n_basis_fns, n_max_spline, n_species, &
                                    n_max_grid, use_basis_gradients
        use rel_x2c_mod,     only : n_basis_small
        use runtime_choices, only : flag_rel, REL_atomic_zora, REL_x2c, REL_4c_dks
        use mpi_tasks,       only : check_allocation

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



        integer:: info
        character(*), parameter :: func = 'allocate_basis'

        if (.not.allocated(basis_wave_spl)) then
           allocate ( basis_wave_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info )
           call check_allocation(info, 'basis_wave_spl                ')
        end if
        if (use_basis_gradients) then
           if (.not.allocated(basis_deriv_spl)) then
              allocate( basis_deriv_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info)
              call check_allocation(info, 'basis_deriv_spl               ')
           end if
        end if
        if (.not.allocated(basis_kinetic_spl)) then
           allocate( basis_kinetic_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info)
           call check_allocation(info, 'basis_kinetic_spl             ')
        end if
        if ((.not.allocated(basis_kinetic_scaled_zora_spl)) .and. (flag_rel==REL_atomic_zora)) then
           allocate( basis_kinetic_scaled_zora_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info)
           call check_allocation(info, 'basis_kinetic_scaled_zora_spl ')
        end if

        if (.not.allocated(basis_wave_ordered)) then
           allocate( basis_wave_ordered(n_basis_fns,n_max_spline, n_max_grid),stat=info)
           call check_allocation(info, 'basis_wave_ordered            ')
        end if
        if (use_basis_gradients) then
           if (.not.allocated(basis_deriv_ordered)) then
              allocate( basis_deriv_ordered(n_basis_fns,n_max_spline, n_max_grid ),stat=info)
              call check_allocation(info, 'basis_deriv_ordered           ')
           end if
        end if

        if (.not.allocated(basis_kinetic_ordered)) then
           allocate( basis_kinetic_ordered(n_basis_fns,n_max_spline, n_max_grid),stat=info)
           call check_allocation(info, 'basis_kinetic_ordered         ')
        end if

        if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
           if (.not.allocated(basis_wave_s_spl)) then
              allocate ( basis_wave_s_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info )
              call check_allocation(info, 'basis_wave_s_spl              ')
           end if
           if (.not.allocated(basis_wave_s_ordered)) then
              allocate( basis_wave_s_ordered(n_basis_fns,n_max_spline, n_max_grid),stat=info)
              call check_allocation(info, 'basis_wave_s_ordered          ')
           end if
           if (.not.allocated(basis_deriv_spl)) then
              allocate( basis_deriv_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info)
              call check_allocation(info, 'basis_deriv_spl               ')
           end if
           if (.not.allocated(basis_deriv_s_spl)) then
              allocate( basis_deriv_s_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info)
              call check_allocation(info, 'basis_deriv_s_spl             ')
           end if
           if (.not.allocated(basis_deriv_ordered)) then
              allocate( basis_deriv_ordered(n_basis_fns,n_max_spline, n_max_grid ),stat=info)
              call check_allocation(info, 'basis_deriv_ordered           ')
           end if
           if (.not.allocated(basis_deriv_s_ordered)) then
              allocate( basis_deriv_s_ordered(n_basis_fns,n_max_spline, n_max_grid ),stat=info)
              call check_allocation(info, 'basis_deriv_s_ordered         ')
           end if
           if (.not.allocated(basis_kinetic_s_spl)) then
              allocate( basis_kinetic_s_spl(n_max_spline,n_max_grid,n_basis_fns),stat=info)
              call check_allocation(info, 'basis_kinetic_s_spl           ')
           end if
           if (.not.allocated(basis_kinetic_s_ordered)) then
              allocate( basis_kinetic_s_ordered(n_basis_fns,n_max_spline, n_max_grid),stat=info)
              call check_allocation(info, 'basis_kinetic_s_ordered       ')
           end if
        endif

        if (.not.allocated(basis_kinetic_scaled_zora_ordered).and. flag_rel==REL_atomic_zora) then
           allocate( basis_kinetic_scaled_zora_ordered(n_basis_fns,n_max_spline,  n_max_grid),stat=info)
           call check_allocation(info, 'basis_kinetic_scaled_zora_orde')
        end if

        if (.not.allocated(outer_radius)) then
           allocate (outer_radius(n_basis_fns),stat=info)
           call check_allocation(info, 'outer_radius                  ')
        end if

        if (.not.allocated(basis_atom)) then
           allocate (basis_atom(n_basis),stat=info)
           call check_allocation(info, 'basis_atom                    ')
        end if
        if (.not.allocated(basis_small_atom)) then
           allocate (basis_small_atom(n_basis_small),stat=info)
           call check_allocation(info, 'basis_small_atom              ')
        end if
        if (.not.allocated(basis_l)) then
           allocate (basis_l(n_basis),stat=info)
           call check_allocation(info, 'basis_l                       ')
        end if
        if (.not.allocated(basis_small_l)) then
           allocate (basis_small_l(n_basis_small),stat=info)
           call check_allocation(info, 'basis_small_l                 ')
        end if
        if (.not.allocated(basis_k)) then
           allocate (basis_k(n_basis),stat=info)
           call check_allocation(info, 'basis_k                       ')
        end if
        if (.not.allocated(basis_small_k)) then
           allocate (basis_small_k(n_basis_small),stat=info)
           call check_allocation(info, 'basis_small_k                 ')
        end if
        if (.not.allocated(basis_m)) then
           allocate (basis_m(n_basis),stat=info)
           call check_allocation(info, 'basis_m                       ')
        end if
        if (.not.allocated(basis_small_m)) then
           allocate (basis_small_m(n_basis_small),stat=info)
           call check_allocation(info, 'basis_small_m                 ')
        end if
        if (.not.allocated(basis_fn)) then
           allocate (basis_fn(n_basis),stat=info)
           call check_allocation(info, 'basis_fn                      ')
        end if
        if (.not.allocated(basis_small_fn)) then
           allocate (basis_small_fn(n_basis_small),stat=info)
           call check_allocation(info, 'basis_small_fn                ')
        end if

        if (.not. allocated(basisfn_species)) then
           allocate(basisfn_species(n_basis_fns), stat=info)
           call check_allocation(info, 'basisfn_species', func)
        end if
        if (.not. allocated(basisfn_l)) then
           allocate(basisfn_l(n_basis_fns), stat=info)
           call check_allocation(info, 'basisfn_l', func)
        end if
        if (.not. allocated(basisfn_type)) then
           allocate(basisfn_type(n_basis_fns), stat=info)
           call check_allocation(info, 'basisfn_type', func)
        end if
        if (.not. allocated(basisfn_n)) then
           allocate(basisfn_n(n_basis_fns), stat=info)
           call check_allocation(info, 'basisfn_n', func)
        end if
        if (.not. allocated(basisfn_k)) then
           allocate(basisfn_k(n_basis_fns), stat=info)
           call check_allocation(info, 'basisfn_k', func)
        end if

        if (.not. allocated(Lsp2n_basis_fnLsp)) then
           allocate(Lsp2n_basis_fnLsp(0:max_basis_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2n_basis_fnLsp', func)
        end if
        if (.not. allocated(Lsp2basis_fn)) then
           allocate(Lsp2basis_fn(max_n_basis_fnLsp, 0:max_basis_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2basis_fn', func)
        end if
        if (.not. allocated(Lsp2basis_sp)) then
           allocate(Lsp2basis_sp(max_n_basis_fnLsp, 0:max_basis_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2basis_sp', func)
        end if
        if (.not. allocated(atom2basis_off)) then
           allocate(atom2basis_off(n_atoms), stat=info)
           call check_allocation(info, 'atom2basis_off', func)
        end if
        if (.not. allocated(sp2n_basis_sp)) then
           allocate(sp2n_basis_sp(n_species), stat=info)
           call check_allocation(info, 'sp2n_basis_sp', func)
        end if

        if (.not.allocated(outer_radius_sq)) then
           allocate (outer_radius_sq(n_basis_fns),stat=info)
           call check_allocation(info, 'outer_radius_sq               ')
        end if

        if (.not.allocated(atom_radius)) then
           allocate (atom_radius(n_species),stat=info)
           call check_allocation(info, 'atom_radius                ')
        end if
        if (.not.allocated(atom_radius_sq)) then
           allocate (atom_radius_sq(n_species),stat=info)
           call check_allocation(info, 'atom_radius_sq                ')
        end if
        
        if (.not.allocated(perm_basis_fns_spl)) then
           allocate (perm_basis_fns_spl(n_basis_fns),stat=info)
           call check_allocation(info, 'perm_basis_fns_spl            ')
        end if
        if (.not.allocated(i_radial_fn)) then
           allocate (i_radial_fn(n_basis_fns),stat=info)
           call check_allocation(info, 'i_radial_fn                   ')
        end if
        
        if (.not.allocated(basis_fn_start_spl)) then
           allocate (basis_fn_start_spl(n_species),stat=info)
           call check_allocation(info, 'basis_fn_start_spl            ')
        end if
        
        if (.not.allocated(n_basis_fn_species)) then
           allocate (n_basis_fn_species(n_species),stat=info)
           call check_allocation(info, 'n_basis_fn_species            ')
        end if
        if (.not.allocated(basis_fn_atom)) then
           allocate (basis_fn_atom(n_basis_fns,n_atoms),stat=info)
           call check_allocation(info, 'basis_fn_atom                 ')
        end if
        if (.not.allocated(n_basis_atom)) then
           allocate (n_basis_atom(n_atoms),stat=info)
           call check_allocation(info, 'n_basis_atom                  ')
        end if
        if (.not.allocated(spline_offset)) then
           allocate(spline_offset(n_species),stat=info)
           call check_allocation(info, 'spline_offset                 ')
        end if
        if (.not.allocated(basis_mapping)) then
           allocate (basis_mapping(n_basis),stat=info)
           call check_allocation(info, 'basis_mapping                 ')
        end if           

        end subroutine allocate_basis
!*****        
!---------------------------------------------------------------------
!***** basis/allocate_ext
!  NAME
!    allocate_ext
!  SYNOPSIS

        subroutine allocate_ext( )

!  PURPOSE
!    Subroutine allocate_ext allocates all auxiliary basis-related arrays, 
!    including the radial functions enabled by the for_aux keyword, 
!    after the standard dimensions n_basis, n_basis_fn (without for_aux) are known.   
!    These are the additional variables for the OBS+ basis set of Fig. 1 of 
!    New Journal of Physics 17, 093020 (2015)

!  USES

        use dimensions,       only : n_atoms, n_ext, n_ext_fns, n_max_grid, &
                                     n_species, n_max_spline, use_basis_gradients
        use runtime_choices,  only : flag_rel, REL_atomic_zora
        use localorb_io,      only : localorb_info, use_unit
        use mpi_tasks,        only : check_allocation

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        character*200 :: info_str
        real*8 memory

        integer:: info
        character(*), parameter :: func = 'allocate_ext'

        memory = 0.d0

        if (.not.allocated(ext_wave_spl)) then
           allocate ( ext_wave_spl(n_max_spline,n_max_grid,n_ext_fns),stat=info )
           call check_allocation(info, 'ext_wave_spl                ')
           memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
        end if
        if (use_basis_gradients) then
           if (.not.allocated(ext_deriv_spl)) then
              allocate( ext_deriv_spl(n_max_spline,n_max_grid,n_ext_fns),stat=info)
              call check_allocation(info, 'ext_deriv_spl               ')
              memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
           end if
        end if
        if (.not.allocated(ext_kinetic_spl)) then
           allocate( ext_kinetic_spl(n_max_spline,n_max_grid,n_ext_fns),stat=info)
           call check_allocation(info, 'ext_kinetic_spl             ')
           memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
        end if
        if ((.not.allocated(ext_kinetic_scaled_zora_spl)) .and. (flag_rel==REL_atomic_zora)) then
           allocate( ext_kinetic_scaled_zora_spl(n_max_spline,n_max_grid,n_ext_fns),stat=info)
           call check_allocation(info, 'ext_kinetic_scaled_zora_spl ')
           memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
        end if      

        if (.not.allocated(ext_wave_ordered)) then
           allocate( ext_wave_ordered(n_ext_fns,n_max_spline, n_max_grid),stat=info)
           call check_allocation(info, 'ext_wave_ordered            ')
           memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
        end if
        if (use_basis_gradients) then
           if (.not.allocated(ext_deriv_ordered)) then
              allocate( ext_deriv_ordered(n_ext_fns,n_max_spline, n_max_grid ),stat=info)
              call check_allocation(info, 'ext_deriv_ordered           ')
              memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
           end if
        end if
        if (.not.allocated(ext_kinetic_ordered)) then
           allocate( ext_kinetic_ordered(n_ext_fns,n_max_spline, n_max_grid),stat=info)
           call check_allocation(info, 'ext_kinetic_ordered         ')
           memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
        end if

        if (.not.allocated(ext_kinetic_scaled_zora_ordered).and. flag_rel==REL_atomic_zora) then
           allocate( ext_kinetic_scaled_zora_ordered(n_ext_fns,n_max_spline,  n_max_grid),stat=info)
           call check_allocation(info, 'ext_kinetic_scaled_zora_orde')
           memory = memory + 8.d0 * n_max_spline*n_max_grid*n_ext_fns
        end if

      if (.not.allocated(outer_ext_radius)) then
           allocate (outer_ext_radius(n_ext_fns),stat=info)
           call check_allocation(info, 'outer_ext_radius                  ')
           memory = memory + 8.d0 * n_ext_fns
        end if

        if (.not.allocated(ext_atom)) then
           allocate (ext_atom(n_ext),stat=info)
           call check_allocation(info, 'ext_atom                    ')
           memory = memory + 4.d0 * n_ext
        end if
        if (.not.allocated(ext_l)) then
           allocate (ext_l(n_ext),stat=info)
           call check_allocation(info, 'ext_l                       ')
           memory = memory + 4.d0 * n_ext
        end if
        if (.not.allocated(ext_m)) then
           allocate (ext_m(n_ext),stat=info)
           call check_allocation(info, 'ext_m                       ')
           memory = memory + 4.d0 * n_ext
        end if
        if (.not.allocated(ext_fn)) then
           allocate (ext_fn(n_ext),stat=info)
           call check_allocation(info, 'ext_fn                      ')
           memory = memory + 4.d0 * n_ext
        end if

        if (.not. allocated(extfn_species)) then
           allocate(extfn_species(n_ext_fns), stat=info)
           call check_allocation(info, 'extfn_species', func)
           memory = memory + 4.d0 * n_ext_fns
        end if
        if (.not. allocated(extfn_l)) then
           allocate(extfn_l(n_ext_fns), stat=info)
           call check_allocation(info, 'extfn_l', func)
           memory = memory + 4.d0 * n_ext_fns
        end if
        if (.not. allocated(extfn_type)) then
           allocate(extfn_type(n_ext_fns), stat=info)
           call check_allocation(info, 'extfn_type', func)
           memory = memory + 4.d0 * n_ext_fns
        end if
        if (.not. allocated(extfn_n)) then
           allocate(extfn_n(n_ext_fns), stat=info)
           call check_allocation(info, 'extfn_n', func)
           memory = memory + 4.d0 * n_ext_fns
        end if

        if (.not. allocated(Lsp2n_ext_fnLsp)) then
           allocate(Lsp2n_ext_fnLsp(0:max_ext_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2n_ext_fnLsp', func)
           memory = memory + 4.d0 * (max_ext_L+1)*n_species
        end if
        if (.not. allocated(Lsp2ext_fn)) then
           allocate(Lsp2ext_fn(max_n_ext_fnLsp, 0:max_ext_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2ext_fn', func)
           memory = memory + 4.d0 * max_n_ext_fnLsp* (max_ext_L+1)* n_species
        end if
        if (.not. allocated(Lsp2ext_sp)) then
           allocate(Lsp2ext_sp(max_n_ext_fnLsp, 0:max_ext_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2ext_sp', func)
           memory = memory + 4.d0 * max_n_ext_fnLsp* (max_ext_L+1)* n_species
        end if
        if (.not. allocated(atom2ext_off)) then
           allocate(atom2ext_off(n_atoms), stat=info)
           call check_allocation(info, 'atom2ext_off', func)
           memory = memory + 4.d0 * n_atoms
        end if
        if (.not. allocated(sp2n_ext_sp)) then
           allocate(sp2n_ext_sp(n_species), stat=info)
           call check_allocation(info, 'sp2n_ext_sp', func)
           memory = memory + 4.d0 * n_species
        end if

        if (.not.allocated(outer_ext_radius_sq)) then
           allocate (outer_ext_radius_sq(n_ext_fns),stat=info)
           call check_allocation(info, 'outer_ext_radius_sq               ')
           memory = memory + 8.d0 * n_ext_fns
        end if

        if (.not.allocated(atom_ext_radius)) then
           allocate (atom_ext_radius(n_species),stat=info)
           call check_allocation(info, 'atom_ext_radius                ')
           memory = memory + 8.d0 * n_species
        end if
        if (.not.allocated(atom_ext_radius_sq)) then
           allocate (atom_ext_radius_sq(n_species),stat=info)
           call check_allocation(info, 'atom_ext_radius_sq                ')
           memory = memory + 8.d0 * n_species
        end if
        
        if (.not.allocated(perm_ext_fns_spl)) then
           allocate (perm_ext_fns_spl(n_ext_fns),stat=info)
           call check_allocation(info, 'perm_ext_fns_spl            ')
           memory = memory + 4.d0 * n_ext_fns
        end if
        if (.not.allocated(i_ext_radial_fn)) then
           allocate (i_ext_radial_fn(n_ext_fns),stat=info)
           call check_allocation(info, 'i_ext_radial_fn                   ')
           memory = memory + 4.d0 * n_ext_fns
        end if
        
        if (.not.allocated(ext_fn_start_spl)) then
           allocate (ext_fn_start_spl(n_species),stat=info)
           call check_allocation(info, 'ext_fn_start_spl            ')
           memory = memory + 4.d0 * n_species
        end if
        
        if (.not.allocated(n_ext_fn_species)) then
           allocate (n_ext_fn_species(n_species),stat=info)
           call check_allocation(info, 'n_ext_fn_species            ')
           memory = memory + 4.d0 * n_species
        end if
        if (.not.allocated(ext_fn_atom)) then
           allocate (ext_fn_atom(n_ext_fns,n_atoms),stat=info)
           call check_allocation(info, 'ext_fn_atom                 ')
           memory = memory + 4.d0 * n_ext_fns * n_atoms
        end if
        if (.not.allocated(n_ext_atom)) then
           allocate (n_ext_atom(n_atoms),stat=info)
           call check_allocation(info, 'n_ext_atom                  ')
           memory = memory + 4.d0 * n_atoms
        end if
        if (.not.allocated(ext_spline_offset)) then
           allocate(ext_spline_offset(n_species),stat=info)
           call check_allocation(info, 'ext_spline_offset                 ')
           memory = memory + 4.d0 * n_species
        end if
        if (.not.allocated(ext_mapping)) then
           allocate (ext_mapping(n_ext),stat=info)
           call check_allocation(info, 'ext_mapping                 ')
           memory = memory + 4.d0 * n_ext
        end if 

        write (info_str,'(2X,A)')         'Per-task memory consumption for arrays in subroutine allocate_ext:'
        call localorb_info(info_str,use_unit,'(A)')        
        write (info_str,'(2X,A,F18.6,A)') '| ', memory/1000000, 'MB.'
        call localorb_info(info_str,use_unit,'(A)')        

        end subroutine allocate_ext     
        
        
!******
!---------------------------------------------------------------------
!****s* basis/cleanup_basis
!  NAME
!    cleanup_basis
!  SYNOPSIS

        subroutine cleanup_basis ( )

!  PURPOSE
!    Subroutine cleanup_basis deallocates all basis-related arrays
!  USES
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




        if (allocated(basis_wave_spl)) then
          deallocate(basis_wave_spl)
        end if
        if (allocated(basis_wave_s_spl)) then
          deallocate(basis_wave_s_spl)
        end if
        if (allocated(basis_deriv_spl)) then
          deallocate(basis_deriv_spl)
        end if
        if (allocated(basis_deriv_s_spl)) then
          deallocate(basis_deriv_s_spl)
        end if
        if (allocated(basis_kinetic_spl)) then
          deallocate(basis_kinetic_spl)
        end if
        if (allocated(basis_kinetic_s_spl)) then
          deallocate(basis_kinetic_s_spl)
        end if
        if (allocated(basis_kinetic_scaled_zora_spl)) then
          deallocate(basis_kinetic_scaled_zora_spl)
        end if

        if (allocated(basis_kinetic_scaled_zora_ordered))then
           deallocate( basis_kinetic_scaled_zora_ordered)
        end if

        if (allocated(outer_radius)) then
          deallocate(outer_radius)
        end if
        if (allocated(outer_radius_sq)) then
          deallocate(outer_radius_sq)
        end if
        if (allocated(atom_radius)) then
          deallocate(atom_radius)
        end if
        if (allocated(atom_radius_sq)) then
          deallocate(atom_radius_sq)
        end if

        if (allocated(basis_atom)) then
          deallocate(basis_atom)
        end if
        if (allocated(basis_small_atom)) then
          deallocate(basis_small_atom)
        end if
        if (allocated(basis_l)) then
          deallocate(basis_l)
        end if
        if (allocated(basis_small_l)) then
          deallocate(basis_small_l)
        end if
        if (allocated(basis_k)) then
          deallocate(basis_k)
        end if
        if (allocated(basis_small_k)) then
          deallocate(basis_small_k)
        end if
        if (allocated(basis_m)) then
          deallocate(basis_m)
        end if
        if (allocated(basis_small_m)) then
          deallocate(basis_small_m)
        end if
        if (allocated(basis_fn)) then
          deallocate(basis_fn)
        end if
        if (allocated(basis_small_fn)) then
          deallocate(basis_small_fn)
        end if

        if (allocated(basisfn_species)) then
           deallocate(basisfn_species)
        end if
        if (allocated(basisfn_l)) then
           deallocate(basisfn_l)
        end if
        if (allocated(basisfn_type)) then
           deallocate(basisfn_type)
        end if
        if (allocated(basisfn_n)) then
           deallocate(basisfn_n)
        end if
        if (allocated(basisfn_k)) then
           deallocate(basisfn_k)
        end if

        if (allocated(Lsp2n_basis_fnLsp)) then
           deallocate(Lsp2n_basis_fnLsp)
        end if
        if (allocated(Lsp2basis_fn)) then
           deallocate(Lsp2basis_fn)
        end if
        if (allocated(Lsp2basis_sp)) then
           deallocate(Lsp2basis_sp)
        end if
        if (allocated(atom2basis_off)) then
           deallocate(atom2basis_off)
        end if
        if (allocated(sp2n_basis_sp)) then
           deallocate(sp2n_basis_sp)
        end if

        if (allocated(perm_basis_fns_spl)) then
           deallocate(perm_basis_fns_spl)
        end if
        if (allocated(i_radial_fn)) then
           deallocate(i_radial_fn)
        end if

        if (allocated(basis_fn_start_spl)) then
           deallocate(basis_fn_start_spl)
        end if

        if (allocated(n_basis_fn_species)) then
           deallocate(n_basis_fn_species)
        end if
        if (allocated(basis_fn_atom)) then
           deallocate(basis_fn_atom)
        end if
        if (allocated(spline_offset)) then
           deallocate(spline_offset)
        end if
        if (allocated(n_basis_atom)) then
           deallocate(n_basis_atom)
        end if
        if (allocated(basis_deriv_ordered)) then
           deallocate(basis_deriv_ordered)
        end if
        if (allocated(basis_deriv_s_ordered)) then
           deallocate(basis_deriv_s_ordered)
        end if
        if (allocated(basis_wave_ordered)) then
           deallocate(basis_wave_ordered)
        end if
        if (allocated(basis_wave_s_ordered)) then
           deallocate(basis_wave_s_ordered)
        end if
        if (allocated(basis_kinetic_ordered)) then
           deallocate(basis_kinetic_ordered)
        end if
        if (allocated(basis_kinetic_s_ordered)) then
           deallocate(basis_kinetic_s_ordered)
        end if
        if (allocated(basis_mapping)) then
          deallocate(basis_mapping)
        end if
        if (allocated(outradsq_save)) then
          deallocate(outradsq_save)
        end if
        if (allocated(invlist_save)) then
          deallocate(invlist_save)
        end if
        end subroutine cleanup_basis
!******
!---------------------------------------------------------------------
!****s* basis/cleanup_basis
!  NAME
!    cleanup_ext_basis
!  SYNOPSIS

        subroutine cleanup_ext_basis( )

!  PURPOSE
!    Subroutine cleanup_ext_basis deallocates all array associated with the
!    OBS+ basis set of Fig. 1 of Ref. New Journal of Physics 17, 093020 (2015)
!  USES
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




        if (allocated(ext_wave_spl)) then
          deallocate(ext_wave_spl)
        end if
        if (allocated(ext_deriv_spl)) then
          deallocate(ext_deriv_spl)
        end if
        if (allocated(ext_kinetic_spl)) then
          deallocate(ext_kinetic_spl)
        end if
        if (allocated(ext_kinetic_scaled_zora_spl)) then
          deallocate(ext_kinetic_scaled_zora_spl)
        end if

        if (allocated(ext_kinetic_scaled_zora_ordered))then
           deallocate( ext_kinetic_scaled_zora_ordered)
        end if

        if (allocated(outer_ext_radius)) then
          deallocate(outer_ext_radius)
        end if
        if (allocated(outer_ext_radius_sq)) then
          deallocate(outer_ext_radius_sq)
        end if
        if (allocated(atom_ext_radius)) then
          deallocate(atom_ext_radius)
        end if
        if (allocated(atom_ext_radius_sq)) then
          deallocate(atom_ext_radius_sq)
        end if
        if (allocated(ext_atom)) then
          deallocate(ext_atom)
        end if
        if (allocated(ext_l)) then
          deallocate(ext_l)
        end if
        if (allocated(ext_m)) then
          deallocate(ext_m)
        end if
        if (allocated(ext_fn)) then
          deallocate(ext_fn)
        end if

        if (allocated(extfn_species)) then
           deallocate(extfn_species)
        end if
        if (allocated(extfn_l)) then
           deallocate(extfn_l)
        end if
        if (allocated(extfn_type)) then
           deallocate(extfn_type)
        end if
        if (allocated(extfn_n)) then
           deallocate(extfn_n)
        end if

        if (allocated(Lsp2n_ext_fnLsp)) then
           deallocate(Lsp2n_ext_fnLsp)
        end if
        if (allocated(Lsp2ext_fn)) then
           deallocate(Lsp2ext_fn)
        end if
        if (allocated(Lsp2ext_sp)) then
           deallocate(Lsp2ext_sp)
        end if
        if (allocated(atom2ext_off)) then
           deallocate(atom2ext_off)
        end if
        if (allocated(sp2n_ext_sp)) then
           deallocate(sp2n_ext_sp)
        end if

        if (allocated(perm_ext_fns_spl)) then
           deallocate(perm_ext_fns_spl)
        end if
        if (allocated(i_ext_radial_fn)) then
           deallocate(i_ext_radial_fn)
        end if

        if (allocated(ext_fn_start_spl)) then
           deallocate(ext_fn_start_spl)
        end if

        if (allocated(n_ext_fn_species)) then
           deallocate(n_ext_fn_species)
        end if
        if (allocated(ext_fn_atom)) then
           deallocate(ext_fn_atom)
        end if
        if (allocated(ext_spline_offset)) then
           deallocate(ext_spline_offset)
        end if
        if (allocated(n_ext_atom)) then
           deallocate(n_ext_atom)
        end if
        if (allocated(ext_deriv_ordered)) then
           deallocate(ext_deriv_ordered)
        end if
        if (allocated(ext_wave_ordered)) then
           deallocate(ext_wave_ordered)
        end if
        if (allocated(ext_kinetic_ordered)) then
           deallocate(ext_kinetic_ordered)
        end if
        if (allocated(ext_mapping)) then
          deallocate(ext_mapping)
        end if
        end subroutine cleanup_ext_basis
    
    
    !******
!---------------------------------------------------------------------
!****s* basis/cleanup_basis
!  NAME
!    cleanup_basis
!  SYNOPSIS

        subroutine cleanup_nonused_ext ( )

!  PURPOSE
!    Subroutine cleanup_basis deallocates all auxiliary basis-related arrays
!  USES
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




        if (allocated(ext_wave_spl)) then
          deallocate(ext_wave_spl)
        end if
        if (allocated(ext_deriv_spl)) then
          deallocate(ext_deriv_spl)
        end if
        if (allocated(ext_kinetic_spl)) then
          deallocate(ext_kinetic_spl)
        end if
        if (allocated(ext_kinetic_scaled_zora_spl)) then
          deallocate(ext_kinetic_scaled_zora_spl)
        end if

        if (allocated(ext_kinetic_scaled_zora_ordered))then
           deallocate( ext_kinetic_scaled_zora_ordered)
        end if

        if (allocated(outer_ext_radius)) then
          deallocate(outer_ext_radius)
        end if
        if (allocated(outer_ext_radius_sq)) then
          deallocate(outer_ext_radius_sq)
        end if
        if (allocated(atom_ext_radius)) then
          deallocate(atom_ext_radius)
        end if
        if (allocated(atom_ext_radius_sq)) then
          deallocate(atom_ext_radius_sq)
        end if
        if (allocated(ext_atom)) then
          deallocate(ext_atom)
        end if
        if (allocated(ext_l)) then
          deallocate(ext_l)
        end if
        if (allocated(ext_m)) then
          deallocate(ext_m)
        end if
        if (allocated(ext_fn)) then
          deallocate(ext_fn)
        end if

        if (allocated(extfn_species)) then
           deallocate(extfn_species)
        end if
        if (allocated(extfn_l)) then
           deallocate(extfn_l)
        end if
        if (allocated(extfn_type)) then
           deallocate(extfn_type)
        end if
        if (allocated(extfn_n)) then
           deallocate(extfn_n)
        end if

        if (allocated(Lsp2n_ext_fnLsp)) then
           deallocate(Lsp2n_ext_fnLsp)
        end if
        if (allocated(Lsp2ext_fn)) then
           deallocate(Lsp2ext_fn)
        end if
        if (allocated(Lsp2ext_sp)) then
           deallocate(Lsp2ext_sp)
        end if
        if (allocated(atom2ext_off)) then
           deallocate(atom2ext_off)
        end if
        if (allocated(sp2n_ext_sp)) then
           deallocate(sp2n_ext_sp)
        end if

        if (allocated(perm_ext_fns_spl)) then
           deallocate(perm_ext_fns_spl)
        end if
        if (allocated(i_ext_radial_fn)) then
           deallocate(i_ext_radial_fn)
        end if

        if (allocated(ext_fn_start_spl)) then
           deallocate(ext_fn_start_spl)
        end if

        if (allocated(n_ext_fn_species)) then
           deallocate(n_ext_fn_species)
        end if
        if (allocated(ext_fn_atom)) then
           deallocate(ext_fn_atom)
        end if
        if (allocated(ext_spline_offset)) then
           deallocate(ext_spline_offset)
        end if
        if (allocated(n_ext_atom)) then
           deallocate(n_ext_atom)
        end if
        if (allocated(ext_deriv_ordered)) then
           deallocate(ext_deriv_ordered)
        end if
        if (allocated(ext_wave_ordered)) then
           deallocate(ext_wave_ordered)
        end if
        if (allocated(ext_kinetic_ordered)) then
           deallocate(ext_kinetic_ordered)
        end if
        if (allocated(ext_mapping)) then
          deallocate(ext_mapping)
        end if
        end subroutine cleanup_nonused_ext
        
!******
!---------------------------------------------------------------------
 end module basis
