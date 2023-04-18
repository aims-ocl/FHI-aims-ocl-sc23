c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           march-may 2008
c     last revision:  march 2012
c###########################################################

      module globalvars
c      **********************************************       
c      description of global variables and data-types
c      **********************************************       
       implicit none

c      http://www.mathwithmrherte.com/pi_digits.htm
       double precision, parameter :: pi=3.14159265358979323846d0

       complex(8), parameter :: ione  = (0.0d0,1.0d0)
       complex(8), parameter :: cone  = (1.0d0,0.0d0)
       complex(8), parameter :: czero = (0.0d0,0.0d0)

c      http://en.wikipedia.org/wiki/Hartree >>>
c           Mohr, Peter J.; Taylor, Barry N.; Newell, David B.
c           "CODATA Recommended Values of the Fundamental Physical
c           Constants: 2006". Rev. Mod. Phys. 80: 633--730 (2008)
       double precision, parameter :: hartree = 27.211384d0
       double precision, parameter :: bohr_radius = 0.529177d0

c      -- global variables --

c      system specific parameters to be found 
c      after reading/importing external files
       integer  num_atoms,      ! number of atoms
     &          num_atom_types, ! number of atom types
     &          num_atom_sorts, ! number of different "sorts" of atoms
     &          num_at_gto,     ! max. number of basis functions at atom   
     &          num_val_gto,    ! max. number of basis functions due to valence electons   
c    &          num_mos,        ! = nsaos, number of molecular orbitals
     &	        nsaos,          ! nsaos: from mos file
     &          nsaos0,         ! saves value of nsaos in case of ecp 
     &          nvalence,       ! reduced dimension of space, if core states are excluded
     &          nspin,          ! amount of spin degrees of freedom
     &          nelectr,        ! number of electrons, will be set later  
     &          nelectr0        ! saves value of nelectr in case of ecp

c      names of external files
       character(8), parameter :: tcntrl_file_name = 'tcontrol'
       character(64) :: crd_file_name      = 'geometry.in',       ! 'coord' 
     &                  crd_tmpfile_name   = 'geometry.tmp',
     &                  basis_file_name    = 'basis-indices.out', ! 'basis'
     &                  tmpbasis_file_name = 'aims.pseudo.basis.tmp',
     &                  mos_file_name      = 'mos.aims',          !  'mos'
     &                  alpha_file_name    = 'alpha.aims',        !  'alpha'
     &                  beta_file_name     = 'beta.aims',         !  'beta'
     &                  output_file_name   = 'TE.dat',            !  'TEV.dat'
     &                  smat12_file_name   = 'smat12',
     &                  smat12inv_file_name= 'smat12inv',
     &                  hsrc_file_name     = 'hsource',
     &                  hion_file_name     = 'hion.data',
     &                  self_energy_file   = 'self.energy.in'
     
c      next file names a fixed to have compartibility with 
c      TURBOMOLE ridft module which reads density matrix ;
c      same remark is applied to the development version of the fhi-aims                   
       character(4), parameter  :: charge_dens_file = 'dmat'
       character(4), parameter  :: spin_dens_file   = 'smat'

c      occupation numbers file name: was needed for testing   
       character(64), parameter :: occ_nmbrs_file = 'occnmbrs'

       character(64) :: outomat_file_name, inomat_file_name
       character(64) :: outhmat_file_name, inhmat_file_name
      
c      this should be adjusted if needed ...
       integer, parameter :: contr  = 10 ! max. number of GTOs inside CGTO
c      ! max. number of diff. orbitals within given l
       integer, parameter :: orbdim = 8

c      max. length allowed to mark atoms 
       integer, parameter :: atsymbol_len = 16 
	       
c      atoms info	    
       type type_atom
        character(atsymbol_len)  symbol  ! atomic symbol ("c","h", etc.)
        integer                  charge  ! atom charge
        integer                  atype   ! atom type (1:c,2:h,3:au,etc.)
        integer                  asort   ! atom sort (1:"c",2:"c ring",3:"h",etc.)
        double precision         pos(3)  ! x,y,z coordinates
        logical                  llead   ! .true. if atom is coupled
        logical                  rlead   ! to L or R leads, respectively
        double precision         exfield ! local exchange field [Hartree]  
        double precision         imse    ! local self-energy (imag.piece) [Hartree]  
       end type
       type(type_atom), allocatable :: atom(:)
c      allocate(atom(1:num_atoms),stat=ierr)

c      different atom names are collected in auxiliary array
c      like 1 - h,  2 - c, ..., 10 - au, etc.
c      here atomic attributes (e.g. "c ring" or "c chain")
c      are not distinguished
       character(atsymbol_len), allocatable :: atom_type(:)
c      allocate(atom_type(1:num_atom_types),stat=ierr)

c      different atomic "sorts" are collected in another auxiliary array
c      where atomic attributes (e.g. "c ring" or "c chain") are distinguished
       character(atsymbol_len), allocatable :: atom_sort(:)
c      allocate(atom_type(1:num_atom_sorts),stat=ierr)

c      auxiliary array keeping similar information as the one above,
c      atomic "sorts" will refer as "c_chain", "c_string", etc.
       character(atsymbol_len), allocatable :: ref_sort(:)
     
c      MO energies,  n = MO index
       double precision, allocatable :: mo_en(:,:)  
c                                     ! (En,ispin)
c      MO expansion coefficients over CGTOs
       double precision, allocatable :: mo_coeff(:,:,:) 
c                                     !     (icoeff,n,ispin)
c      after reading mos/alpha/beta files we get
c       nspin = 1 ! mos
c       nspin = 2 ! alpha, beta spins
       
c      allocate(mo_en(num_mos,nspin),stat=ierr)
c      allocate(mo_coeff(nsoas,num_mos,nspin),stat=ierr)

c      expansion coefficients in orthogonal basis
       double precision, allocatable :: mo_orth(:,:,:)

c      spectrum info >>>
c      (this type is used for the density-matrix related routines) 
       type spectrum
        double precision  en       ! mo energy
        integer           spin     ! =1 or 2 : spin index
        integer           moindex  ! original index of mo
       end type

c      contracted GTOs info       
       type cgto 
        integer           ngto           ! degree of contraction == number of GTO  
        double precision  icoeff(contr)  ! expansion coeffs 
        double precision  xi(contr)      ! exponents 
        integer           lm             ! lm = 1 2  3  4   5 6 7 8 9
c                                        !      s px py pz  d  ...  d
        logical           valence        ! .false. --> core, .true. --> valence state 
       end type

c      >>> TURBOMOLE functions:
c      >>> s-functions:  lm = 1
c      >>> p-functions:  lm = 2   3   4
c                             px  py  pz
c      >>> d-functions   lm = 5  (-xx-yy+2zz)/sqrt(12)
c                             6  xz
c                             7  yz
c                             8  xy
c                             9  (xx-yy)/2
c      >>> f-functions   lm =10  (2zzz-3xxz-3yyz)/sqrt(60)
c                            11  (-xxx-xyy+4xzz)/sqrt(40)
c                            12  (-yyy-xxy+4yzz)/sqrt(40)
c                            13  xyz
c                            14  (xxz-yyz)/2
c                            15  (xxx-3xyy)/sqrt(24)
c                            16  (yyy-3xxy)/sqrt(24)
c      integer, parameter :: lmax0 = 3  ! up to f-functions

c      >>> for FHI-aims code, larger lmax is required 
c                             up to to h-functions :
       integer, parameter :: lmax0 = 5 


       type(cgto), allocatable :: aos(:,:)
c      allocate(aos(1:num_atom_types,1:num_at_gto),stat=ierr)
c               ...like   aos(atype,orb_number)        
c      we get 'aos' from file 'basis'       
      
c      cartesian CGTOs info       
       type cart_cgto        
        double precision icoeff(contr) ! expansion coeffs 
        integer          pows(3)       ! 1:px,2:py,3:pz ==> pows of polynoms        
       end type 
      
       integer, allocatable :: n_basis_func(:)
c              ! number of basis functions of given atom type
c      allocate(n_basis_func(1:num_atom_types),stat=ierr)

       integer, allocatable :: n_val_basis_func(:)
c      ! in case of fhi-aims input, core states may be integrated out >>>
c      ! number of basis functions at atom associated with valence electrons 

c      ! maximum n quantum number of the core state for given atom type
       integer, allocatable :: ncore(:)
       integer :: ncore_states = 0

c      overlap matrix       
       double precision, allocatable :: smat(:,:)
c      square-root of overlap matrix
       double precision, allocatable :: smat12(:,:)
c      and inv(smat12)
       double precision, allocatable :: smat12inv(:,:)

c      normalizations of CGTO : orbnorm = |<mu|mu>|
       double precision, allocatable :: orbnorm(:)

c      equilibrium density matrix (per spin); last index == ispin
c      double precision, allocatable :: densmat0(:,:,:)
c      non-equilibrium density matrix (per spin); last index == ispin
c      double precision, allocatable :: neqdmat(:,:,:)

c      import/export overlap matrix: relevant for the fhi-aims 
       logical :: save_omat  = .false.
       logical :: read_omat  = .false.
c      export ks-hamiltonian matrix
       logical :: save_hmat  = .false.

c      the '$aims_input on' flag is used 
c      (i)  to read a 'geometry.in' file in the fhi-aims format
c      (ii) to avoid the use of turbomole 'basis' file by
c           reading  the fhi-aims file 'basis-indices.out' instead
       logical :: aims_input = .false.
c      to protect FHI-aims users, is set here to .true. anyway !

c                                       hamiltonian matrices
       double precision, allocatable :: h0mat(:,:,:)
c                                        ! bare hamiltonian
       complex(8), allocatable       :: hmat(:,:,:),
c                                        ! including a self-energy sigma
     &                                  sigmann(:),
c                                        ! a diagonal part of self-energy
     &                                  bsigmann(:)
c                          ! a diagonal part of self-energy for b-channel
       double precision, allocatable :: gamma_left(:),
c                                        ! gamma_L: left electrode
     &                                  gamma_right(:)
c                                        ! gamma_R: right electrode

c      imaginary part of the self-energy at given atom: 
       double precision, allocatable :: at_self_energy(:)

       double precision :: eta = 1.0d-10
c                          ! infinitesimal eta for the Green's function

c      --   arrays needed for non-equilibrium implementation   --
c           complex eigen-energies of the "extended hamiltonian" h+sigma  
       complex(8), allocatable :: zpoles(:,:),      ! (en,ispin)
c                                 corresponding eigenvectors             
     &                            heigvec(:,:,:),   ! (n,p,ispin)
c                                 and inverse of eigenvectors matrix            
     &                            invheigvec(:,:,:) ! (p,n,ispin)

c      non-equilibrium density matrix (per spin); last index == ispin
       double precision, allocatable :: neqdmat(:,:,:)

c      testing flag
       logical :: testing   = .false. 

       logical :: zspectrum = .false. 
c      if .true., outputs a complex-valued  
c      spectrum to zmos/zalpha/zbeta files  

c      a flag '$landauer' is introduced 
c      (i)   for testing purposes;  
c      (ii)  it fixes the only choice available so far for
c            the fhi-aims user in the first common release 
c      (iii) it means:
c            take a self-energy with imaginarty piece only,
c            build up a Green's function and compute ballistic 
c            Landauer transmission function      

c      this original line(s) can be commented out  ...
c      logical :: landauer    = .false.
c      logical :: do_landauer = .false.
c      >>>  becomes '.true.' if '$landauer on' is found in <tcontrol>
 
c      instead, next lines can be used ...
c      ! CAUTION : next line is to protect a default fhi-aims user ;
c      !           in the first publicly available release, it
c      !           prevents using any advanced options except of a basic
c      !           option "$landauer on" to compute transmission function
       logical ::  landauer    = .true.
       logical ::  do_landauer = .false.

c      ! in case of fhi-aims input, core states can be integrated out
       logical ::  ecp = .false. 
        
c      flag for a self-consitency cycle for the noneq.dens.matrix
       logical :: do_dmat_cycle = .false.  
c      ! becomes .true. if keyword "$densmat_cycle on" is found

c      default value of admixing factor for the density matrix
       double precision :: defdmix = 0.01d0
c      admixing  factor ($dmix keyword)
       double precision dmix
       
c      self.-cons.-cycle iteration number
       integer :: iter = 0 

c      flag for transmission calculation       
       logical :: do_trans_calc = .false.  
c      ! becomes .true. if keyword "$transmission  on" is found

c      >>> inserted by RK & AB :..
c      flag for transmission eigenchannels
       logical :: do_cond_channels = .false.
c      ! becomes .true. if keyword "$transmission channels" is found
c      >>> done with insert by RK & AB.

c      if .plt files with scattering wave functions are required ...
       logical :: do_scatt_wave_func = .false.
c      ! becomes .true. if keyword "$transmission wave_functions" is found

c      ! scaling factor for scattering wave functions output (optional) 
c      ! will be modified by input or estimated automatically        
       double precision :: wf_scaling_fctr = -1.0d0

c      interval [in a.u.] along z-axis for computation of 
c      scattering wave functions
       double precision :: wf_zmin = -1.0d0, wf_zmax  = 1.0d0
c      ! modified by 'wf_z_range <zmin> <zmax>' in $ldos group   

c      global array contains eigenchannels as solution of the 
c      eigenvalue problem for the transmissio matrix tt+
       complex(8), allocatable :: utt(:,:,:)

c      flag for optional calculation of local density of states
       logical :: do_ldos = .false.
c      ! becomes .true. if keyword '$ldos on' is found

c      flag for optional calculation of lm-projected local density of states
       logical :: do_lmdos = .false.
c      ! becomes .true. if keyword '$ldos lm-projected' is found

c      ================== space-resolved LDOS ==============================
c      flags for optional calculation of r-resolved density of states
c      ! becomes .true. if keyword '$ldos 3d' is found
       logical :: do_ldos3d = .false.
c      ! second option is to integrate over given energy window
c      ! becomes .true. if keyword '$ldos 3d_ewindow' is found
       logical :: do_ldos3d_ewindow = .false.
     
c      interval [in a.u.] along z-axis for computation of the r-resolved DOS
       double precision :: zmin = -1.0d0, zmax  = 1.0d0
c      ! modified by 'z_range <zmin> <zmax>' in $ldos group   

c      ! scaling factor (optional) 
c      ! will be modified by input or estimated automatically        
       double precision :: scaling_fctr = -1.0d0
c      ================== space-resolved LDOS ==============================

c      optional flag for (lm-decomposed) population analysis
       logical :: do_pop = .false.
c                ! becomes .true. if keyword '$population on' is found
       logical :: do_lmpop = .false.
c      ! becomes .true. if keyword '$population lm-projected' is found

c      "pseudo-occupation numbers" were needed only for testing purposes!
c      flag for saving pseudo-occupation numbers       
c      logical :: store_occn = .true.  
c      ! .true./.false. depending on keyword "$saveoccn on/off" in <tcontrol>

c      flag for import the overlap matrices in 
c      a binary format: may be helpful for the NEGF cycle
       logical :: importsmat = .false. 
c      if becomes .true. after 1st iteration, then sqrt(overlap) and its 
c      inverse are stored on hard drive and read later on for the next 
c      interations -- in that case, however, precision is slightly lost 
c      (??? check why it is so!) so procedure is justified only for big 
c      systems where one can indeed gain in cpu-time 

       logical :: storesmat = .false. 
c      if .true. allows above flag 'importsmat' to do its job ;
c      when .false., sqrts of overlap matrices are calculated at each iteration

c      next flag controls output of loewdin charges/magn.moments 
       logical :: qoutput = .false. 
c      becomes active if "$output_charges on" is found in tcontrol-file
             
c      specification of left/right surfaces
       integer left_satoms(3), right_satoms(3)
c      first option: keyword "$dist"
c      thickness of a slab (coupling region) 
c      connected to reservoirs defined in atomic units 
       double precision :: dist = 4.0d0
       logical ::          dist_found = .false.

c      second option : keyword "$nlayers"
c      "thickness" of a slab (coupling region) connected 
c      to reservoirs defined in atomic layers
       integer :: nlayers = 1  
       logical :: nlayers_found = .false.

       integer, parameter :: deep = 3
c                           ! layers 3,4, etc ... are treated as one

c      tolerance parameter to identify group of atoms 
c      belonging to the same crystallographic plane
       double precision :: dz_tolerance = 3.0d-2

c      model self-energy    ! absorbing boundary conditions     
       double precision rsigma(deep), isigma(deep)
       complex(8)       sigma(deep)

c      import self-energy flag, becomes true if a group 
c      "$self_energy file=<...>" if found in <tcontrol>
       logical :: import_self_energy = .false.

c      >>> adjusting real piece of the self-energy for given EF
       logical :: tune_rsigma = .false.
c      ! becomes .true. if keyword '$adjust_rsigma  on' is found
       logical :: expert = .false.

c      initial guess for a real piece ( r.sigma = -rguess*|im.sigma| )
       double precision :: rguess = 0.1d0
c      el.number conv.criterium for the re.sigma search
       double precision :: nelcnv  = 1.0d-5
c      maximum number of allowed iterations
       integer :: iter_limit = 15
c      in case of expert run (adjust rsigma + self-cons. cycle for density)
c      we take into account previous guess for re(sigma), and 
c      use as the second point rfactor * re(sigma)
       double precision :: rfactor = 0.9d0
c      <<<< end of $adjust_rsigma group

c      flag to allow spin-polarized treatment of electrodes
       logical :: sp_leads = .false.
c      ! switches to .true. if a group '$hexchange on/off' is active (i.e. "on")

c      exchange field in left/right leads [Hartree]
       double precision :: hex_left = 0.0d0, hex_right = 0.0d0

c      optional flag: if keyword $hsource is found, 
c      local exchange fields on atoms are introduced;
c      requires input-file <hsrc_file_name>
       logical :: localexf = .false.

c      optional flag: helps building the "domain-wall" solution
c      in case of antiparallel magnetizations in the electrodes;
c      if switched on, matrix elements of the hamiltonian
c      are symmetrized appropriately, so that 
c      "junction's left side/spin-up" equals "junction's right-side/spin-down"   
       logical :: z_symmetrize  = .false. 
       logical :: ap_symmetrize = .false. 
c      becomes .true. if flag "$symmetrize ap" is found in <tcontrol> 
c
c      similar flag to insure symmetric solution without domain-wall
       logical :: p_symmetrize = .false. 
c      becomes .true. if flag "$symmetrize p" is found in <tcontrol> 

c      further flags to insure symmetric solutions
       logical :: symmetrize   = .false. 
       logical :: x_symmetrize = .false.
       logical :: y_symmetrize = .false.

c      optional flag to perform calculations for the 
c      isolated molecule ($no-leads on)
       logical :: no_leads = .false.

       logical :: mo_restart = .false.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      ... variables for the 'lda+u' type of calculation ...

c      "hubbard-u" related dimensions >>>
       integer ion_size,       ! dimension of "ion" block in hamiltonian
     &         res_size        ! dimension of "reservoirs" block

c      lda+u:  on/off == true/false
       logical :: ldau_old = .false. 
       logical :: ldau     = .false. 

c      >>> specification of reservoirs >>>
       integer :: num_ratoms = 0 
c      number of different reservoirs' atom types
       integer, allocatable :: ratom_types(:) 
c      specific atom sorts which belong to reservoirs

c      parameters below accept temporary values and  
c      are to be updated once <tcontrol> file is read 

       double precision ::  u_j = 4.0d0
c      screened coulomb interaction [eV]

       integer :: nlocorb = 5 ! at least, one d-ion
c      number of localized orbitals (per spin)
       double precision :: loc_cutoff = 0.9d0
c      only orbitals with localization degree above 
c      'cut-off' value are corrected

c      parameters related to $hion group
       logical :: save_hion = .false.
c      if .true., a hamiltonian of the 3d-ion and its density matrix
c      is saved to external file : may be used further to approach 
c      a required magnetic solution for molecules absorbed on surfaces 
       logical :: read_hion = .false.
c      if .true. hamiltonian of the ion and its dens.matrix
c      are read from external file

c      lm-decomposed local density of states at correlated 
c      site(s) to output into external file
       logical :: spectr_func = .false.

c      flag for importing density matrix 
       logical :: read_dmat = .false.

c      projector operators
       double precision, allocatable :: osmat(:,:)
       complex(8), allocatable       :: cmplx_osmat(:,:)
   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      ... variables for the 'mol_subsystem' type of calculation ...
c
c      i.e. if the density matrix, hamiltonian and the local DOS to be 
c      projected on the "molecule" excluding reservoirs' degrees of freedom:
c      this option is expected to be used as an addtional analysis tool

c      "mol_subsystem" related dimensions >>>

       integer molecule_size,       ! dimension of "molecular" block in hamiltonian
     &         reservoir_size       ! dimension of "reservoirs" block

c      mol_subsystem flag:  on/off == true/false
       logical :: mol_subsystem = .false.

c      >>> specification of reservoirs >>>
       integer :: num_res_atoms = 0 
c      number of different reservoirs' atom types
       integer, allocatable :: res_atom_types(:) 
c      specific atom sorts which belong to reservoirs

c      orbital projected local density of states at the "molecule"
c      to be output into external file
       logical :: mol_spectr_func = .false.
       
c      amount of occupied and empty molecular orbitals 
c      onto which the local DOS is to be projected ...       
       integer :: occ_mol_orbitals = 10, empty_mol_orbitals = 10 

c      import molecular hamiltonian from external file
       logical :: read_hmat  = .false.

c      projector operators
c      double precision, allocatable :: osmatrix(:,:)
c      complex(8), allocatable       :: cmplx_osmatrix(:,:)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      energy points [Hartree]
       double precision ener,            ! 1st point 
     &                  estep,           ! step
     &                  eend             ! last point  
       logical :: conductance = .false.  ! .true.  if energy mesh is not defined 
                                         ! .false. otherwise

c      flag to switch to eV units for energy points,
c      becomes .true. if $evunits is set to "on"
       logical :: evunits = .false.

c      bias-voltage [eV]
       double precision :: bias = 0.0d0

c      testing >>>
       integer :: occ=0, a_occ=0, b_occ=0 ! number of occupied orbitals
       logical :: calc_densmat0 = .false.  
c                 ! if .true. -- equilibrium dens.matrix D is evaluated, 
c                 ! and a number of electrons is calculated as N = tr(DS)

c      flag for the fermi energy
       logical :: efermi_guess = .false.       
c      ! becomes .true. if value of $efermi is found in tcontrol       

c      ! fixed efermi during iterative procedure >>>>
       logical :: fixed_efermi = .false.

c      efermi = (\mu_left+\mu_right)/2: 
c      here is set to zero to be on the safe side ...
       double precision :: efermi = 0.0d0

c      parameters for the efermi search
       double precision :: nelconv = 1.0d-7   
c                          ! conv.cretirium for the electrons number
       integer          :: nclstmo = 2
c                          ! number of orbitals per spin taken up/down
c                          ! from the efermi-guess to estimate density of states  

c      convergence checks
       double precision :: densconv = 1.0d-4, dnorm = 1.0d0
       logical          :: dnormfound = .false.    
       double precision, parameter :: fctr = 0.1d0
     
      end module globalvars

