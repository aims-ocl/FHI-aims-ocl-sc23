!****h* FHI-aims/ll_vdwdf
!  NAME
!    ll_vdwdf 
!  SYNOPSIS

module ll_vdwdf

!  PURPOSE
!   Module that connects to the LL van der Waals density functional developed
!   by Dion and coworkers ("Langreth-Lundqvist functional) as implemented by
!   the group of Claudia Ambrosch-Draxl and coworkers. 
!
!   !!! We are very grateful for code provided by the group of Claudia Ambrosch-Draxl
!   !!! and coworkers that has been used in the implementation below. If you use this
!   !!! functionality, please cite their work.
!
!    Subroutines included in this module:
!    * allocate_ll_vdw
!    * read_ll_vdw_parameters
!    * calc_ll_vdw
!    * local_x_c_total_energy
!    * local_lda_x_c_total_energy
!    * local_pbe_x_c_total_energy
!    * local_x_c_total_energy_even_spacing_grid 
!    * local_lda_x_c_energy_even_spacing_grid
!    * local_lda_x_c_total_energy_even_spacing_grid 
!    * local_pbe_x_c_total_energy_even_spacing_grid 
!    * obtain_charge_density
!    * charge_density
!    * charge_density_even_grid
!    * charge_density_even_grid0
!    * charge_density_even_grid0_p1
!    * charge_density_trans
!    * charge_density_grad_even_grid  
!    * density_interpolation_even_grid
!    * cube_cell_dimension
!    * integrand
!    * integrand_even_spacing_grid
!    * calc_kernel
!    * read_kernel_data
!    * cleanup_ll_vdw
!    * get_center
!    * output_ll_vdwdf
!  AUTHOR
!    Significant portions of this code have been developed by the group of
!    Claudia Ambrosch-Draxl, Peter Puschnig, Dmitrii Nabok and others in Leoben.
!    Please cite their work when using this implementation.
!    
!    Portions of this module developed by FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Portions of this version:
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement.
!
!   Copyright by Claudia Ambrosch-Draxl, Peter Puschnig, Dmitrii Nabok and others in Leoben
!   for significant portions of the code upon which this module is based. 
!  SEE ALSO
!
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009), submitted.
!
!    M. Dion, H. Rydberg, E. Schroeder, D. C. Langreth, and B. I.
!    Lundqvist, Phys. Rev. Lett. 92, 246401 2004
!
!    Claudia Ambrosch-Draxl, P. Puschnig, Dmitrii Nabok, code developed at University of Leoben. 
!    !!! A full citation will follow when their implementation is published - please check back with us !!!
!
!  HISTORY
!    Release version, FHI-aims (2009).
!  SOURCE
   use localorb_io,only:use_unit
   implicit none

    real*8,  dimension(:), allocatable   :: cell_origin

    real*8    :: lat_vec(3,3)
    real*8    :: lat(3)
    real*8    :: offset_coord(3)

    integer   :: cell_size_vdw(3)=0


! Cuba Monte-Carlo integration variables (input from control.in)

    integer          :: ndim, ndim_eg
    integer          :: ncomp
    integer          :: verbose
    integer          :: key1, key2, key3
    integer          :: maxpass

    double precision  :: epsrel, epsabs

    parameter (ndim = 2)
    parameter (ndim_eg = 6)
    parameter (ncomp = 1)

    integer          :: last
    integer          :: ngiven
    integer          :: ldxgiven, ldxgiven_eg
    integer*8        :: nextra

    double precision :: border, maxchisq, mindeviation
    integer*8        :: mineval, mineval_eg
    integer*8        :: maxeval

    parameter (last = 4)
    parameter (ldxgiven_eg = ndim_eg)
    parameter (ldxgiven = ndim)
    parameter (mineval = 16*ndim)
    parameter (mineval_eg = 16*ndim_eg)
    parameter (border = 0D0)
    parameter (maxchisq = 10D0)
    parameter (mindeviation = .25D0)
  
! read kernel data variables

    character*80     :: kernelfile 

    real*8,  dimension(:,:), allocatable :: kernel
    real*8           :: dmin,dmax,dstep
    real*8           :: deltamin,deltamax,deltastep
    real*8           :: eps
    real*8           :: dcutoff

    parameter (eps = 1.0d-10)
    parameter (dcutoff = 0.02)

    integer          :: i_count
    integer          :: time0, time1, corate
    integer          :: tot_time, minutes, seconds

!---------------------------------------------------------------------------------------------
    contains
!---------------------------------------------------------------------------------------------

!******	

subroutine allocate_ll_vdw ()
!---------------------------------------------------------------------------------------------
! Subroutine allocate_ll_vdw allocates the variables related to even grid charge density
! and their gradients. 
!---------------------------------------------------------------------------------------------
!
!    use dimensions
!    use physics

    implicit none
!
! allocate cube parameters
!
    allocate ( cell_origin(3) )
!    allocate ( cell_edge_units(3))
!    allocate ( cell_edge_steps(3))
!
end subroutine allocate_ll_vdw


subroutine read_ll_vdw_parameters ()
!---------------------------------------------------------------------------------------------
! Read parameters for ll_vdw_functional
!FIXME: if "integration_grid" option is not included in the control.in file in case of even spacing grids.
!---------------------------------------------------------------------------------------------
!
    use constants
    use runtime_choices
    use localorb_io
    use dimensions
    use physics

    implicit none

    character*20 desc_str
    integer i_code

    logical :: flag_es  = .false.
    logical :: flag_eu  = .false.

    logical :: flag_eos = .false.
    logical :: flag_eof = .false.

    logical :: flag_edge = .false.
    integer :: i_edge = 0

    integer :: i_coord, nelect
    character*140 :: info_str

    character*80     :: xsffile 

    verbose = 0
    key1 = 5000000
    key2 = 1
    key3 = 1
    maxpass = 5

    epsrel = 1.d-16
    epsabs = 1.d-02

    flag_eg = .true.

    write(info_str,'(4X,A)') "| Read ll_vdw_functional parameters: "
    call localorb_info(info_str)

! read ll_vdw_functional parameters

    read (7,*,iostat=i_code) desc_str

    if (i_code.ne.0) then
      write(use_unit,*) "End-of-file instead of ll_vdw_functional parameters?"
      stop
    end if

    do while (.not.flag_eos)

      if (desc_str(1:1).eq."#") then

          continue

      else if (desc_str.eq.'vdwdf') then

          backspace (7)

          read(7,*) desc_str, desc_str 

          if (desc_str.eq.'cell_grid') then
     
              backspace (7)
     
              read(7,*) desc_str, desc_str, desc_str
     
              if (desc_str .eq. 'even_spacing') then 
     
                 flag_eg = .true.
     
              else
     
                 write(use_unit,'(1X,A,A,A)') &
                 "* Unknown even_spacing specification ", desc_str, &
                 ". Please correct."
                 stop
     
              endif

          else if (desc_str.eq.'cell_origin') then

              backspace (7)

              read(7,*) desc_str, desc_str, &
             (cell_origin(i_coord), i_coord=1,3,1)
              do i_coord=1,3,1
              cell_origin(i_coord) = cell_origin(i_coord)/bohr
              enddo
              flag_cell = .true.

          else if (desc_str.eq.'cell_edge_steps') then

              backspace (7)

              read(7,*) desc_str, desc_str, &
             (cell_edge_steps(i_coord), i_coord=1,3,1)
              if ( (cell_edge_steps(1).le.0.d0) .or. &
                   (cell_edge_steps(2).le.0.d0) .or. &
                   (cell_edge_steps(3).le.0.d0) ) then
                   write(use_unit,'(1X,A,A)') &
                   "* Zero or less grid steps are not possible. ", &
                   "Please correct."
                   stop
              end if
              flag_es = .true.

          else if (desc_str.eq.'cell_edge_units') then

              backspace (7)

              read(7,*) desc_str, desc_str, &
             (cell_edge_units(i_coord), i_coord=1,3,1)
              do i_coord=1,3,1
              cell_edge_units(i_coord) = cell_edge_units(i_coord)/bohr
              enddo
              if (n_periodic <1) then
              if ( (cell_edge_units(1).eq.0.d0) .and. &
                   (cell_edge_units(2).eq.0.d0) .and. &
                   (cell_edge_units(3).eq.0.d0) ) &
                then
                write(use_unit,'(1X,A,A)') &
                "* Zero-length grid units are not possible. ", &
                "Please correct."
                stop
              end if
              endif
              flag_eu = .true.

          else if (desc_str.eq.'cell_size') then

              backspace(7)
      
              read(7,*) desc_str, desc_str, &
              (cell_size_vdw(i_coord), i_coord=1,3,1)
      
          else

             write(use_unit,'(1X,A,A,A)') &
             "* Unknown vdwdf specification ", desc_str, &
             ". Please correct."
             stop
    
          endif

      else if (desc_str.eq.'mc_int') then

           backspace(7)
    
           read(7,*) desc_str, desc_str
    
           if (desc_str.eq.'kernel_data') then

               flag_kernel = .true.
    
               backspace(7)
      
               read(7,*) desc_str, desc_str, kernelfile 
               call read_kernel_data()
      
           else if ((desc_str.eq.'output_flag').or.(desc_str.eq.'Output_flag')) then

               backspace(7)
      
               read(7,*) desc_str, desc_str, verbose

           else if ((desc_str.eq.'number_of_MC').or.(desc_str.eq.'Number_of_MC')) then

               backspace(7)
      
               read(7,*) desc_str, desc_str, key1
      
           else if (desc_str.eq.'key2') then

               backspace(7)
      
               read(7,*) desc_str, desc_str, key2
      
           else if (desc_str.eq.'key3') then

               backspace(7)
      
               read(7,*) desc_str, desc_str, key3
      
           else if (desc_str.eq.'Maxpass') then

               backspace(7)
      
               read(7,*) desc_str, desc_str, maxpass 
      
           else if ((desc_str.eq.'relative_accuracy').or.(desc_str.eq.'Relative_accuracy')) then

               backspace(7)
      
               read(7,*) desc_str, desc_str, epsrel 
      
           else if ((desc_str.eq.'absolute_accuracy').or.(desc_str.eq.'Absolute_accuracy')) then

               backspace(7)
      
               read(7,*) desc_str, desc_str, epsabs
      
           else
    
             write(use_unit,'(1X,A,A,A)') &
             "* Unknown mc_int specification ", desc_str, &
             ". Please correct."
             stop
    
           end if

     else

       backspace(7)
       flag_eos = .true.
!
     end if
!
!     next line

    read (7,*,iostat=i_code) desc_str
!
    if (i_code.ne.0) then
      flag_eos = .true.
      flag_eof = .true.
    end if
!
  end do
!
   if (.not.flag_eof) then
     backspace(7)
   end if

    if (n_periodic>1) then
       if ((.not.flag_eu).and.(.not.flag_es)) then
          write(use_unit,'(1X,A)') &
          "* Please offer cell_edge_steps or cell_edge_units. "
          stop
       endif
    else
       if ((.not.flag_eu).or.(.not.flag_es)) then
          write(use_unit,'(1X,A)') &
          "* Please offer cell_edge_steps and cell_edge_units. "
          stop
       endif
    endif

    if (.not.flag_kernel) then
      write(use_unit,'(1X,A)') &
      "* Please offer kernel data in order to calculate LL-vdw term. "
      stop
    end if

!
! set up basic parameters for cube_cell
!
!    call cube_cell_dimension()

end subroutine read_ll_vdw_parameters


subroutine calc_ll_vdw(ll_vdw_energy,ll_vdw_energy_err,lda_en_c,pbe_en_c)
!---------------------------------------------------------------------------------------------
! Calculate energy contributaion from the Langreth-Lundqvist van der Waals functional 
!---------------------------------------------------------------------------------------------

    use physics
    use dimensions
    use localorb_io
    use mpi_utilities
    use synchronize_mpi
    use constants
    use xc

    use runtime_choices
    use grids
    use geometry 
    use species_data
    use pbc_lists

    implicit none

! output
!    real*8, dimension(3, n_atoms) :: i_forces
    real*8 :: ll_vdw_energy, ll_vdw_energy_err, func(6), x(6)
    real*8 :: lda_en_c, pbe_en_c
    
! local variables

    integer*8        :: neval
    integer          :: nregions, fail
    double precision :: integral(ncomp), error(ncomp), prob(ncomp)

    integer          :: i,j,k
    integer          :: cell(3)

    real*8           :: en_c, en_c2, en_x, en_xc2, en_xc0
    real*8           :: dens(n_spin), dens_grad(3,n_spin)
    real*8           :: en_density_c, en_density_c2, en_density_x, en_density_xc, en_density_xc2
    real*8, dimension(n_spin) :: pot_x, pot_c, pot_xc
    real*8, dimension(n_spin) :: pot_x2, pot_c2, pot_xc2
    real*8           :: en_pot_x, en_pot_c, en_pot_xc0
    real*8           :: en_pot_xc2
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8           :: lda_en_x, pbe_en_x

    real*8 time_ll_vdwdf
    real*8 tot_time_ll_vdwdf
    real*8 rtime

    character*160    :: info_str

! count
    integer          :: i_point, i_spin, i_dir, i_count

!    external integrand


    maxeval=1
    do i_count=1,32
       maxeval=maxeval*2
    enddo
! FIXME force calculation..

!    i_forces = i_forces
 
! calculating the non_local part of the correlation energy

    ll_vdw_energy = 0.d0
    ll_vdw_energy_err = 0.d0

    ngiven=0
    nextra=0

      call cpu_time(time_ll_vdwdf)

    if(myid.eq.0) then
    i_count=0
    if (flag_eg) then

        call lldivonne(ndim_eg, ncomp, integrand_even_spacing_grid,  &
        epsrel, epsabs, verbose, mineval_eg, maxeval,	 &
        key1, key2, key3, maxpass,			 &
        border, maxchisq, mindeviation,		 &
        ngiven, ldxgiven_eg, 0, nextra, 0,		 &
        nregions, neval, fail, integral, error, prob)

    else

        call lldivonne(ndim, ncomp, integrand,  &
        epsrel, epsabs, verbose, mineval, maxeval,	 &
        key1, key2, key3, maxpass,			 &
        border, maxchisq, mindeviation,		 &
        ngiven, ldxgiven, 0, nextra, 0,		 &
        nregions, neval, fail, integral, error, prob)

    endif 
   endif

    if (myid == 0) then
        ll_vdw_energy = integral(1)
        ll_vdw_energy_err = error(1)
    endif

    call cpu_time(rtime)
    tot_time_ll_vdwdf = rtime- time_ll_vdwdf

    call local_lda_x_c_total_energy (lda_en_x,lda_en_c)
    call local_pbe_x_c_total_energy (pbe_en_x,pbe_en_c)
 
      if(myid.eq.0) then
       write(use_unit,'(10X,A)')"---------------------------------------------"
       write(use_unit,'(A)') " "
       write(use_unit,'(A)') " "
       write(use_unit,'(A)') " "

       write(use_unit,'(10X,A)')"---------------------------------------------"

       write(use_unit,'(10X,A)') "| LL_vdwdf calculation comes to stop ...  "
       write(use_unit,'(A)') " "
       write(use_unit,'(10X,A)')"---------------------------------------------"
       write(use_unit,'(10X,A,F12.3,A)') &
          "| Total time for calculating non-local correlation energy from ll_vdwdf: ", &
           tot_time_ll_vdwdf, " s"

       write(use_unit,'(A)')
       write(use_unit,'(10X,A)')"---------------------------------------------"

      endif

end subroutine calc_ll_vdw


subroutine local_x_c_total_energy &
           (lda_en_x,lda_en_c,pbe_en_x,pbe_en_c)
!---------------------------------------------------------------------------------------------
! Calculate local exchange and correlation energy and potential for the case of 
! pwlda and pbe or revpbe
!---------------------------------------------------------------------------------------------

    implicit none

! output
    real*8 :: lda_en_x,lda_en_c
    real*8 :: pbe_en_x,pbe_en_c
    
    lda_en_c =0.0
    lda_en_x =0.0

    pbe_en_c =0.0
    pbe_en_x =0.0

    call local_lda_x_c_total_energy (lda_en_x,lda_en_c)
    call local_pbe_x_c_total_energy (pbe_en_x,pbe_en_c)

end subroutine local_x_c_total_energy 

subroutine local_lda_x_c_total_energy (lda_en_x,lda_en_c)
!---------------------------------------------------------------------------------------------
! Calculate local exchange and correlation energy and potential for the case of pwlda 
!---------------------------------------------------------------------------------------------

    use physics
    use dimensions
    use localorb_io
!    use mpi_utilities
!    use synchronize_mpi
    use constants
    use xc

    use runtime_choices
    use grids
    use geometry 
    use species_data
    use pbc_lists

    implicit none

! output
    real*8 :: lda_en_x,lda_en_c
    
! local variables
    real*8           :: dens(n_spin), dens_grad(3,n_spin)
    real*8           :: lda_en_density_c, lda_en_density_x, lda_en_density_xc
    real*8           :: lda_en_xc0,lda_en_pot_c,lda_en_pot_x,lda_en_pot_xc0 
    real*8, dimension(n_spin)   :: lda_pot_x, lda_pot_c, lda_pot_xc
    real*8 :: tmp

    logical          :: xc_undefined

! count
    integer          :: i_point, i_spin, i_dir

    lda_en_c =0.0
    lda_en_x =0.0
    lda_en_xc0=0.0
    lda_en_pot_c=0.0
    lda_en_pot_x=0.0
    lda_en_pot_xc0=0.0


    do i_point =1,n_full_points

       if (partition_tab(i_point).gt.0.0d0) then

              xc_undefined = .true.

              do i_spin = 1, n_spin, 1
                 ! This is true if both densities are below zero
                 xc_undefined = xc_undefined .and. &
                      (rho(i_spin,i_point).le.0.d0)
              enddo

              do i_spin = 1, n_spin, 1
                 ! This is true if one of the densities is less than zero
                 xc_undefined = xc_undefined .or. &
                      (rho(i_spin,i_point).lt.0.d0)
              enddo

              if (xc_undefined) then

                 lda_en_density_x = 0.d0
                 lda_en_density_c = 0.d0
                 lda_en_density_xc = 0.d0

              else

                 call x_cpot_pw91_lda &
                 (rho(1:n_spin,i_point), lda_pot_x, lda_pot_c, lda_pot_xc, &
                 lda_en_density_x, lda_en_density_c, lda_en_density_xc)

              endif
 
            do i_spin=1,n_spin
     
              lda_en_c   = lda_en_c     + lda_en_density_c*rho(i_spin,i_point)*partition_tab(i_point)
              lda_en_x   = lda_en_x     + lda_en_density_x*rho(i_spin,i_point)*partition_tab(i_point)
              lda_en_xc0 = lda_en_xc0   + lda_en_density_xc*rho(i_spin,i_point)*partition_tab(i_point)
 
              lda_en_pot_c   = lda_en_pot_c   + lda_pot_c(i_spin)*dens(i_spin)*partition_tab(i_point)
              lda_en_pot_x   = lda_en_pot_x   + lda_pot_x(i_spin)*dens(i_spin)*partition_tab(i_point)
              lda_en_pot_xc0 = lda_en_pot_xc0 + lda_pot_xc(i_spin)*dens(i_spin)*partition_tab(i_point)
 
            enddo

       endif

    enddo

end subroutine local_lda_x_c_total_energy 

subroutine local_lda_x_c_energy_even_spacing_grid (density,lda_en_x,lda_en_c)
!---------------------------------------------------------------------------------------------
! Calculate local exchange and correlation energy and potential for the case of pwlda 
! from density of evan_spacing_grids
!---------------------------------------------------------------------------------------------

    use physics
    use dimensions
    use localorb_io
!    use mpi_utilities
!    use synchronize_mpi
    use constants
    use xc

    use runtime_choices
    use grids
    use geometry 
    use species_data
    use pbc_lists

    implicit none

! output
    real*8 :: lda_en_x, lda_en_c
    
! local variables
    real*8           :: density(n_spin)
    real*8           :: lda_en_density_c, lda_en_density_x, lda_en_density_xc, lda_en_xc
    real*8           :: lda_en_pot_c, lda_en_pot_x
    real*8, dimension(n_spin)   :: lda_pot_x, lda_pot_c, lda_pot_xc
    real*8 :: tmp

! count
    integer          :: i_point, i_spin, i_dir

    lda_en_c =0.0
    lda_en_x =0.0
    lda_en_pot_c=0.0
    lda_en_pot_x=0.0

     call x_cpot_pw91_lda &
    (density, lda_pot_x, lda_pot_c, lda_pot_xc, &
     lda_en_density_x, lda_en_density_c, lda_en_density_xc)

     lda_en_c=lda_en_density_c
     lda_en_x=lda_en_density_x
     lda_en_xc=lda_en_density_xc
 
end subroutine local_lda_x_c_energy_even_spacing_grid

subroutine local_pbe_x_c_total_energy (pbe_en_x,pbe_en_c)
!---------------------------------------------------------------------------------------------
! Calculate local exchange and correlation energy and potential for the case of pbe or revpbe
!---------------------------------------------------------------------------------------------

    use physics
    use dimensions
    use localorb_io
!    use mpi_utilities
!    use synchronize_mpi
    use constants
    use xc

    use runtime_choices
    use grids
    use geometry 
    use species_data
    use pbc_lists

    implicit none

! output
    real*8 :: pbe_en_x,pbe_en_c
    
! local variables
    real*8           :: dens(n_spin), dens_grad(3,n_spin)
    real*8           :: pbe_en_xc0,pbe_en_pot_c,pbe_en_pot_x,pbe_en_pot_xc0 
    real*8           :: pbe_en_density_c, pbe_en_density_x, pbe_en_density_xc
    real*8, dimension(n_spin)   :: pbe_pot_x, pbe_pot_c, pbe_pot_xc
    real*8, dimension(n_spin)   :: pbe_local_xc_derivs
    real*8, dimension(3,n_spin) :: pbe_xc_gradient_deriv
    real*8 :: tmp

    logical          :: xc_undefined


! count
    integer          :: i_point, i_spin, i_dir

    pbe_en_c =0.0
    pbe_en_x =0.0
    pbe_en_xc0=0.0
    pbe_en_pot_c=0.0
    pbe_en_pot_x=0.0
    pbe_en_pot_xc0=0.0

    do i_point =1,n_full_points

       if (partition_tab(i_point).gt.0.0d0) then

           xc_undefined = .true.

           do i_spin = 1, n_spin, 1
              ! This is true if both densities are below zero
              xc_undefined = xc_undefined .and. &
                   (rho(i_spin,i_point).le.0.d0)
           enddo

           do i_spin = 1, n_spin, 1
              ! This is true if one of the densities is less than zero
              xc_undefined = xc_undefined .or. &
                   (rho(i_spin,i_point).lt.0.d0)
           enddo

           if (xc_undefined) then

              pbe_en_density_x = 0.d0
              pbe_en_density_c = 0.d0
              pbe_en_density_xc = 0.d0

           else

              if (flag_xc.eq.6) then
                 call x_c_partials_pbe &
                 (rho(1:n_spin,i_point), rho_gradient(1:3,1:n_spin,i_point), pbe_en_density_x, pbe_en_density_c, &
                  pbe_en_density_xc, pbe_local_xc_derivs(1), pbe_xc_gradient_deriv(1,1))
              else if (flag_xc.eq.12) then
                 call x_c_partials_revpbe &
                 (rho(1:n_spin,i_point), rho_gradient(1:3,1:n_spin,i_point), pbe_en_density_x, pbe_en_density_c, &
                  pbe_en_density_xc, pbe_local_xc_derivs(1), pbe_xc_gradient_deriv(1,1))
              endif
    
           endif
    
           do i_spin=1,n_spin
    
              pbe_en_c   = pbe_en_c   + pbe_en_density_c*rho(i_spin,i_point)*partition_tab(i_point)
              pbe_en_x   = pbe_en_x   + pbe_en_density_x*rho(i_spin,i_point)*partition_tab(i_point)
              pbe_en_xc0 = pbe_en_xc0 + pbe_en_density_xc*rho(i_spin,i_point)*partition_tab(i_point)
    
              pbe_en_pot_xc0=pbe_en_pot_xc0+pbe_local_xc_derivs(i_spin)*rho(i_spin,i_point)*partition_tab(i_point)
    
              do i_dir=1,3
                 pbe_en_pot_xc0=pbe_en_pot_xc0+2.d0*pbe_xc_gradient_deriv(i_dir,i_spin)* &
                           rho_gradient(i_dir,i_spin,i_point)*partition_tab(i_point)
              enddo
    
           enddo

       endif

    enddo

end subroutine local_pbe_x_c_total_energy 

subroutine local_x_c_total_energy_even_spacing_grid &
           (lda_en_x,lda_en_c,pbe_en_x,pbe_en_c)
!---------------------------------------------------------------------------------------------
! Calculate local exchange and correlation energy and potential of even_spacing_grids 
! for the case of pwlda and pbegga
!------------------------------------------------------------------------------------------

    implicit none

! output
    real*8 :: pbe_en_x,pbe_en_c
    real*8 :: lda_en_x,lda_en_c
    
    lda_en_c =0.0
    lda_en_x =0.0

    pbe_en_c =0.0
    pbe_en_x =0.0

    call local_lda_x_c_total_energy_even_spacing_grid (lda_en_x,lda_en_c)
    call local_pbe_x_c_total_energy_even_spacing_grid (pbe_en_x,pbe_en_c)

end subroutine local_x_c_total_energy_even_spacing_grid

subroutine local_lda_x_c_total_energy_even_spacing_grid (lda_en_x,lda_en_c)
!---------------------------------------------------------------------------------------------
! Calculate local exchange and correlation energy and potential of even_spacing_grids 
! for the case of pwlda
!------------------------------------------------------------------------------------------

    use physics
    use dimensions
    use localorb_io
    use mpi_utilities
    use synchronize_mpi
    use constants
    use xc

    use runtime_choices
    use grids
    use geometry 
    use species_data
    use pbc_lists

    implicit none

! output
    real*8 :: lda_en_x,lda_en_c
    
! local variables
    real*8           :: dens(n_spin), dens_grad(3,n_spin)
    real*8           :: lda_en_density_c, lda_en_density_x, lda_en_density_xc
    real*8           :: lda_en_xc0,lda_en_pot_c,lda_en_pot_x,lda_en_pot_xc0 
    real*8, dimension(n_spin)   :: lda_pot_x, lda_pot_c, lda_pot_xc
    real*8           :: weight

! count
    integer          :: i_point, i_spin, i_dir
    integer          :: i, j, k

    lda_en_c =0.0
    lda_en_x =0.0
    lda_en_xc0=0.0
    lda_en_pot_c=0.0
    lda_en_pot_x=0.0
    lda_en_pot_xc0=0.0

!FIXME spin components..

    dens=0.0
    dens_grad=0.0
    weight=1./dble(cell_edge_steps(1)*cell_edge_steps(2)*cell_edge_steps(3))

    do i=1,cell_edge_steps(1)
    do j=1,cell_edge_steps(2)
    do k=1,cell_edge_steps(3)

       if (rho_e(i,j,k).gt.eps) then
          dens(1)=rho_e(i,j,k)
          do i_dir=1,3
             dens_grad(i_dir,1)=rho_e_gradient(i_dir,i,j,k)
          enddo

          call x_cpot_pw91_lda &
          (dens, lda_pot_x, lda_pot_c, lda_pot_xc, &
           lda_en_density_x, lda_en_density_c, lda_en_density_xc)

          do i_spin=1,n_spin
    
             lda_en_c   = lda_en_c     + lda_en_density_c*dens(i_spin)*weight
             lda_en_x   = lda_en_x     + lda_en_density_x*dens(i_spin)*weight
             lda_en_xc0 = lda_en_xc0   + lda_en_density_xc*dens(i_spin)*weight

             lda_en_pot_c   = lda_en_pot_c   + lda_pot_c(i_spin)*dens(i_spin)*weight
             lda_en_pot_x   = lda_en_pot_x   + lda_pot_x(i_spin)*dens(i_spin)*weight
             lda_en_pot_xc0 = lda_en_pot_xc0 + lda_pot_xc(i_spin)*dens(i_spin)*weight

          enddo
       endif

    enddo
    enddo
    enddo

end subroutine local_lda_x_c_total_energy_even_spacing_grid

subroutine local_pbe_x_c_total_energy_even_spacing_grid (pbe_en_x,pbe_en_c)
!---------------------------------------------------------------------------------------------
! Calculate local exchange and correlation energy and potential of even_spacing_grids 
! for the case of pbe
!------------------------------------------------------------------------------------------

    use physics
    use dimensions
    use localorb_io
    use mpi_utilities
    use synchronize_mpi
    use constants
    use xc

    use runtime_choices
    use grids
    use geometry 
    use species_data
    use pbc_lists

    implicit none

! output
    real*8 :: pbe_en_x,pbe_en_c
    
! local variables
    real*8           :: dens(n_spin), dens_grad(3,n_spin)
    real*8           :: pbe_en_density_c, pbe_en_density_x, pbe_en_density_xc
    real*8           :: pbe_en_xc0,pbe_en_pot_c,pbe_en_pot_x,pbe_en_pot_xc0 
    real*8, dimension(n_spin)   :: pbe_pot_x, pbe_pot_c, pbe_pot_xc
    real*8, dimension(n_spin)   :: pbe_local_xc_derivs
    real*8, dimension(3,n_spin) :: pbe_xc_gradient_deriv
    real*8           :: weight

! count
    integer          :: i_point, i_spin, i_dir
    integer          :: i, j, k

    pbe_en_c =0.0
    pbe_en_x =0.0
    pbe_en_xc0=0.0
    pbe_en_pot_c=0.0
    pbe_en_pot_x=0.0
    pbe_en_pot_xc0=0.0

!FIXME spin components..

    dens=0.0
    dens_grad=0.0
    weight=1./dble(cell_edge_steps(1)*cell_edge_steps(2)*cell_edge_steps(3))

    if (flag_xc.eq.6) then

       do i=1,cell_edge_steps(1)
       do j=1,cell_edge_steps(2)
       do k=1,cell_edge_steps(3)

          if (rho_e(i,j,k).gt.eps) then

             dens(1)=rho_e(i,j,k)
             do i_dir=1,3
                dens_grad(i_dir,1)=rho_e_gradient(i_dir,i,j,k)
             enddo

             call x_c_partials_pbe &
             (dens, dens_grad, pbe_en_density_x, pbe_en_density_c, pbe_en_density_xc, &
              pbe_local_xc_derivs, pbe_xc_gradient_deriv)

             do i_spin=1,n_spin
       
                pbe_en_c       = pbe_en_c       + pbe_en_density_c*dens(i_spin)*weight
                pbe_en_x       = pbe_en_x       + pbe_en_density_x*dens(i_spin)*weight
                pbe_en_xc0     = pbe_en_xc0     + pbe_en_density_xc*dens(i_spin)*weight
                pbe_en_pot_xc0 = pbe_en_pot_xc0 + pbe_local_xc_derivs(i_spin)*dens(i_spin)*weight

                do i_dir=1,3
                   pbe_en_pot_xc0=pbe_en_pot_xc0+2.d0*pbe_xc_gradient_deriv(i_dir,i_spin)* &
                             dens_grad(i_dir,i_spin)*weight
                enddo

             enddo

          endif

       enddo
       enddo
       enddo

    else if (flag_xc.eq.12) then

       do i=1,cell_edge_steps(1)
       do j=1,cell_edge_steps(2)
       do k=1,cell_edge_steps(3)

          if (rho_e(i,j,k).gt.eps) then
    
             dens(1)=rho_e(i,j,k)
             do i_dir=1,3
                dens_grad(i_dir,1)=rho_e_gradient(i_dir,i,j,k)
             enddo

             call x_c_partials_revpbe &
             (dens, dens_grad, pbe_en_density_x, pbe_en_density_c, pbe_en_density_xc, &
              pbe_local_xc_derivs, pbe_xc_gradient_deriv)
    
             do i_spin=1,n_spin
       
                pbe_en_c       = pbe_en_c       + pbe_en_density_c*dens(i_spin)*weight
                pbe_en_x       = pbe_en_x       + pbe_en_density_x*dens(i_spin)*weight
                pbe_en_xc0     = pbe_en_xc0     + pbe_en_density_xc*dens(i_spin)*weight
                pbe_en_pot_xc0 = pbe_en_pot_xc0 + pbe_local_xc_derivs(i_spin)*dens(i_spin)*weight
    
                do i_dir=1,3
                   pbe_en_pot_xc0=pbe_en_pot_xc0+2.d0*pbe_xc_gradient_deriv(i_dir,i_spin)* &
                             dens_grad(i_dir,i_spin)*weight
                enddo
    
             enddo
    
          endif

       enddo
       enddo
       enddo

    endif

end subroutine local_pbe_x_c_total_energy_even_spacing_grid

subroutine obtain_charge_density( )
!---------------------------------------------------------------------------------------------
! Obtain charge densities and their gradient projected on spacing grids
!
! output
! rho_e ! charge density projected on even grids 
! rho_e_gradient ! gradient of charge densities 
! rho_e_2gradient ! square of the gradients of charge densities 
!FIXME spin components..
!---------------------------------------------------------------------------------------------

    use runtime_choices

    implicit none

    if (flag_eg) then
    
        call charge_density_even_grid()

    end if

end subroutine obtain_charge_density


subroutine charge_density_even_grid( )
!---------------------------------------------------------------------------------------------
! Obtain charge densities and their gradient projected on even grids
!
! output
! rho_e ! charge density projected on even grids
! rho_e_gradient ! gradient of charge densities projected on even grids
! rho_e_2gradient ! square of the gradients of charge densities 
!---------------------------------------------------------------------------------------------

    use physics
    use dimensions
    use mpi_utilities
    
    implicit none

    integer          :: i, j, k

!    call get_my_task()

    if (n_periodic>1) then
       call charge_density_even_grid0_p1( )
    else
       call charge_density_even_grid0( )
    endif

!
! Translation of charge density to the unit_cell
!
   call charge_density_trans( )
!
! Checking input density on negative values and zeros  
!
    do i=1,cell_edge_steps(1)
    do j=1,cell_edge_steps(2)
    do k=1,cell_edge_steps(3)
!
       if ((rho_e(i,j,k).lt.0d0).or.(abs(rho_e(i,j,k)).lt.eps)) rho_e(i,j,k)=eps
!
    end do
    end do
    end do
!
! Obtaining gradient of charge densities
!
   call charge_density_grad_even_grid()

end subroutine charge_density_even_grid 

subroutine charge_density(i_point,dens,dens_2gradient)
!---------------------------------------------------------------------------------------------
! density and square of density gradients at a given point
!---------------------------------------------------------------------------------------------
!
    use physics
    use dimensions
    use geometry

    implicit none

! input
    integer :: i_point

! output
    real*8 :: dens       
    real*8 :: grad(3), dens_2gradient 


! local variables
    integer :: i, j, k, i_spin
    integer :: i_dir 

    real*8  :: rho_local
  
!FIXME: check with where n_angular and n_division should play a role..
!How to determine pick.. - how to obtain i_full_points from above information

    dens=0.0
    grad=0.0
    if (n_spin.gt.1) then
       do i_spin=1,n_spin
          dens =dens+rho(i_spin,i_point)
          do i_dir=1,3
             grad(i_dir)=grad(i_dir)+rho_gradient(i_dir,i_spin,i_point)
          enddo
       enddo
    else
       dens =rho(1,i_point)
       do i_dir=1,3
          grad(i_dir)=rho_gradient(i_dir,1,i_point)
       enddo
    endif

    if ((dens.lt.0d0).or.(abs(dens).lt.eps)) dens=eps
    dens_2gradient=grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3)


return

end subroutine charge_density

 
subroutine charge_density_even_grid0( )
!---------------------------------------------------------------------------------------------
! Obtain charge densities projected on even grids
!---------------------------------------------------------------------------------------------
!
! TODO: figure out automatically which region to plot
! FIXME: ATM all plots use the same grid  (grid of the first plot)
! ATM only canonical cartesian basis is supported

      use physics
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use species_data
      use mpi_utilities
      use localorb_io
      use basis
      use cartesian_ylm
      use constants

      implicit none

!  local variables

      real*8 inv_bohr_3
      real*8 coord_current(3)

      real*8 dist_tab(n_atoms, cell_edge_steps(3) )
      real*8 dist_tab_sq(n_atoms, cell_edge_steps(3))
      real*8 i_r(n_atoms, cell_edge_steps(3) )
      real*8 dir_tab(3, n_atoms, cell_edge_steps(3) )
      real*8 trigonom_tab(4, n_atoms, cell_edge_steps(3) )
      real*8 radial_wave(n_basis, cell_edge_steps(3) )
      real*8 wave(n_basis, cell_edge_steps(3) )

      integer :: n_compute
      integer :: i_basis(n_basis)
      integer :: orb_idx=0

      integer :: n_compute_fns
      integer :: i_basis_fns(n_basis_fns*n_atoms)
      integer :: i_basis_fns_inv(n_basis_fns,n_atoms)
      integer :: i_atom_fns(n_basis_fns*n_atoms)

      integer :: n_compute_atoms
      integer :: atom_index(n_atoms)
      integer :: atom_index_inv(n_atoms)

      integer :: spline_array_start(n_atoms)
      integer :: spline_array_end(n_atoms)

!     other local variables

      integer, dimension(n_spin) :: max_occ_number
      real*8, dimension(n_states,n_spin) :: occ_numbers_sqrt

      integer :: l_ylm_max
      integer :: n_points

      real*8, dimension(n_basis, n_states, n_spin) :: &
        KS_vec_times_occ_sqrt

      real*8, dimension(n_states, n_basis, n_spin) :: &
        KS_ev_compute

      real*8, dimension(:,:,:), allocatable ::  KS_orbital
      real*8, dimension(:,:), allocatable :: local_rho
      real*8, dimension(:,:), allocatable ::  local_orb_up
      real*8, dimension(:,:), allocatable ::  local_orb_down

      integer, dimension(:,:), allocatable :: index_lm
      real*8, dimension(:,:,:), allocatable :: ylm_tab

      real*8, dimension(:,:,:), allocatable :: dir_tab_global
      real*8, dimension(:,:), allocatable :: dist_tab_sq_global

      real*8 cell_units(3)
      real*8, dimension(:,:), allocatable :: coord_tmp, coord_tmp1 

!     counters
      integer :: coord_x, coord_y, coord_z
      integer :: i_l
      integer :: i_m
      integer :: i_state
      integer :: i_point

      integer :: i_full_points

      integer :: i_spin = 1
      integer :: i_index

      inv_bohr_3 = 1.0d0/(bohr**3)
!     begin work


!      if (n_spin.ne.1) then
!         write(use_unit,*) "vdwdf currently supports non-spin-polarized",&
!              "ensities only! "
!         return
!      endif

      cell_units(1)= cell_edge_units(1)
      cell_units(2)= cell_edge_units(2)
      cell_units(3)= cell_edge_units(3)

      l_ylm_max = l_wave_max

      allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
       cell_edge_steps(3) ) )
      allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )
      allocate( coord_tmp(3,n_atoms), coord_tmp1(3,n_atoms) )

!     translate center of current coordinates to offset_coord

        call trans_cell(coord_tmp,offset_coord)

        coord_tmp1=coords
        coords=coord_tmp
        coord_tmp=coord_tmp1

!     initialize index_lm

        i_index = 0
        
        do i_l = 0, l_ylm_max, 1
            do i_m = -i_l, i_l
             i_index = i_index + 1
             index_lm(i_m, i_l) = i_index
            enddo
        enddo


!       find the maximal occupation number

        do i_spin = 1, n_spin, 1
        ! initialize
            max_occ_number(i_spin) = 0
            do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,1)).gt.0.d0) then
              max_occ_number(i_spin) = i_state
              exit
             endif
            enddo
        enddo

!allocate

        allocate( local_rho(1:cell_edge_steps(3),1:2) )
        allocate(local_orb_up(max_occ_number(1), &
            1:cell_edge_steps(3)))

        local_orb_up = 0.d0
        if (n_spin .gt. 1) then
              allocate(local_orb_down(max_occ_number(2), &
                        1:cell_edge_steps(3)))
              local_orb_down = 0.d0
        endif

        allocate( KS_orbital(n_states,1:cell_edge_steps(3),1:2) )
        KS_orbital = 0.d0


!     allocate the sqrt array for occupation numbers and fill it
!     up to max_occ_number

      do i_spin = 1, n_spin, 1
        do i_state = 1, max_occ_number(i_spin), 1

          occ_numbers_sqrt(i_state, i_spin) = &
            sqrt(occ_numbers(i_state,i_spin,1))

          KS_vec_times_occ_sqrt(:,i_state,i_spin) = &
          KS_eigenvector(:,i_state,i_spin,1) * &
          occ_numbers_sqrt(i_state,i_spin)

        enddo
      enddo

      i_full_points = 0

      do coord_x =   1,cell_edge_steps(1),1

      if(myid.eq.0) then
         print *, 'evaluating charge density of grid x:',coord_x
      endif

        do coord_y = 1,cell_edge_steps(2),1

           n_compute = 0
           i_basis = 0
           i_point = 0

           local_rho = 0.d0

           do coord_z = 1,cell_edge_steps(3),1
            i_point = i_point+1

!    generate output grid
              coord_current(1)=cell_units(1)*(coord_x-1)+cell_origin(1)
              coord_current(2)=cell_units(2)*(coord_y-1)+cell_origin(2)
              coord_current(3)=cell_units(3)*(coord_z-1)+cell_origin(3)
!              coord_current(1)=cell_units(1)*(coord_x-1)-offset_coord(1)+cell_origin(1)
!              coord_current(2)=cell_units(2)*(coord_y-1)-offset_coord(2)+cell_origin(2)
!              coord_current(3)=cell_units(3)*(coord_z-1)-offset_coord(3)+cell_origin(3)

!     compute atom-centered coordinates of current integration point as viewed from all atoms
                 call tab_atom_centered_coords_v2 &
                      ( coord_current,dist_tab_sq(1,i_point), &
                      dir_tab(1,1,i_point))

!     determine which basis functions are relevant at current grid point,
!     and tabulate their indices

!     next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                 call prune_basis_v2 &
                      (dist_tab_sq(1,i_point), n_compute, i_basis)

!          end loop over the z component
           enddo

           n_points = i_point
           if (n_compute.gt.0) then

             ! Determine all radial functions, ylm functions and their derivatives that
             ! are best evaluated strictly locally at each individual grid point.
             i_point = 0
             do coord_z = 1,cell_edge_steps(3),1
                   i_point = i_point+1
                   n_compute_atoms = 0
                   n_compute_fns = 0
                   i_basis_fns_inv = 0

                ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
                ! without any copying and without doing any unnecessary operations.
                ! The price is that the interface is no longer explicit in terms of physical
                ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
               call prune_radial_basis_v2 &
                 ( dist_tab_sq(1,i_point), &
                   n_compute_atoms, atom_index, atom_index_inv, &
                   n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                   i_atom_fns, spline_array_start, spline_array_end )

               ! Tabulate distances, unit vectors, and inverse logarithmic grid units
               ! for all atoms which are actually relevant
               call tab_local_geometry &
               ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                 dir_tab(1,1,i_point), dist_tab(1,i_point), &
                 i_r(1,i_point) &
               )

               ! Now evaluate radial functions u(r) from the previously stored compressed spline arrays
               call evaluate_radial_functions_v2 &
               (   spline_array_start, spline_array_end, &
                   n_compute_atoms, n_compute_fns, &
                   dist_tab(1,i_point), i_r(1,i_point), &
                   atom_index, i_basis_fns_inv, &
                   basis_wave_ordered, radial_wave(1,i_point), &
                   .false. &
                    )

                 call tab_trigonom_v2 &
                       ( n_compute_atoms, dir_tab(1,1,i_point), &
                       trigonom_tab(1,1,i_point) &
                       )

                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm_v2 &
                       ( n_compute_atoms, atom_index, &
                       trigonom_tab(1,1,i_point), l_shell_max, &
                       l_ylm_max, &
                       ylm_tab(1,1,i_point) )


               ! tabulate total wave function value for each basis function in all cases -
               ! but only now are we sure that we have ylm_tab ...
               call evaluate_waves_v2 &
                   (l_ylm_max, ylm_tab(1,1,i_point), &
                    dist_tab(1,i_point), &
                   index_lm, n_compute, &
                   i_basis, radial_wave(1,i_point), &
                   wave(1,i_point), n_compute_atoms, &
                   atom_index_inv, &
                   i_basis_fns_inv &
                   )
!        end loop over z component
         enddo

      end if

      if (n_compute.gt.0) then

          do i_spin = 1, n_spin, 1
            if (max_occ_number(i_spin).gt.0) then
                  call evaluate_KS_density_v2 &
                      (n_points, wave(1,1), n_compute, &
                      i_basis, KS_vec_times_occ_sqrt(1,1,i_spin), &
                      KS_ev_compute(1,1,i_spin), &
                      max_occ_number(i_spin), &
                      occ_numbers_sqrt(1,i_spin), &
                      KS_orbital(1,1,i_spin), local_rho(1,i_spin) &
                       )
            else
                if (i_spin .eq. 1) then
                    call evaluate_KS_orbital_density &
                       (orb_idx,n_points, wave(1,1), n_compute, &
                       i_basis, KS_eigenvector(1,1,i_spin,1), &
                       KS_ev_compute(1,1,i_spin), &
                       max_occ_number(i_spin), &
                       local_orb_up(1,i_spin), &
                       local_rho(1,i_spin) &
                       )
                else
                    call evaluate_KS_orbital_density &
                       (orb_idx,n_points, wave(1,1), n_compute, &
                       i_basis, KS_eigenvector(1,1,i_spin,1), &
                       KS_ev_compute(1,1,i_spin), &
                       max_occ_number(i_spin), &
                       local_orb_down(1,i_spin), &
                       local_rho(1,i_spin) &
                       )
                endif
            endif
          end do

      else

            local_rho = 0.d0

      end if

! assign the density

      do coord_z = 1,cell_edge_steps(3),1
         if (n_spin.eq.1) then
             rho_e(coord_x,coord_y,coord_z)=                    &
             local_rho(coord_z,1) 
         else
             rho_e(coord_x,coord_y,coord_z)=                    &
             (local_rho(coord_z,1) + local_rho(coord_z,2))
         endif

      enddo

!     end loop over y component integration loop
      end do

!     end loop over x component
      end do

!     recover to the original coordinates
      coords=coord_tmp

!     deallocation

      if (allocated(ylm_tab)) then
        deallocate(ylm_tab)
      end if
      if (allocated(index_lm)) then
        deallocate(index_lm)
      end if
      if (allocated(dir_tab_global)) then
         deallocate(dir_tab_global)
      end if
      if (allocated(dist_tab_sq_global)) then
         deallocate(dist_tab_sq_global)
      end if
      if (allocated(local_rho)) then
         deallocate(local_rho)
      end if
      if (allocated(KS_orbital)) then
         deallocate(KS_orbital)
      end if
      if (allocated(local_orb_up)) then
         deallocate(local_orb_up)
      end if
      if (allocated(local_orb_down)) then
         deallocate(local_orb_down)
      end if

end subroutine charge_density_even_grid0

subroutine charge_density_even_grid0_p1()
!---------------------------------------------------------------------------------------------
! Obtain charge densities of periodic system projected on even grids
!---------------------------------------------------------------------------------------------
!
! FIXME: ATM all plots use the same grid  (grid of the first plot)
! FIXME: ATM only canonical cartesian basis is supported 
! FIXME: only non-spin-polarized calculations are supported 

      use mpi_tasks
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use species_data
      use mpi_utilities
      use localorb_io
      use basis
      use cartesian_ylm
      use constants
      use plot
      use physics
      use synchronize_mpi
      use density_matrix_evaluation
      use pbc_lists
      
      implicit none

! local variables
      real*8 inv_bohr_3
      real*8 sqrt_inv_bohr_3
     
      real*8 coord_current(3), coord_current_temp(3)
      real*8 dist_tab(n_centers_basis_integrals)
      real*8 dist_tab_sq(n_centers_basis_integrals, cell_edge_steps(3))
      real*8 i_r(n_centers_basis_integrals)
      real*8 dir_tab(3, n_centers_basis_integrals, cell_edge_steps(3))
      real*8 dir_tab_norm(3, n_centers_basis_integrals)
      real*8 trigonom_tab(4, n_centers_basis_integrals)
      real*8 en_upper_limit,en_lower_limit


      real*8,dimension(:),allocatable:: radial_wave
      real*8,dimension(:,:),allocatable:: wave

      integer :: n_compute, n_compute_a
      integer,dimension(:),allocatable :: i_basis
      integer :: i_cell
     
      integer :: n_compute_fns
      integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
      integer :: i_basis_fns_inv(n_basis_fns,n_centers)
      integer :: i_atom_fns(n_basis_fns*n_centers_integrals)
     
      integer :: n_compute_atoms, i_k_point, i_k
     
      integer :: atom_index(n_centers_integrals)
      integer :: atom_index_inv(n_centers)
     
      integer :: spline_array_start(n_centers_integrals)
      integer :: spline_array_end(n_centers_integrals)

! other local variables
      integer, dimension(n_spin,n_k_points) :: max_occ_number
      real*8, dimension(n_states,n_spin) :: occ_numbers_sqrt
     
      integer :: l_ylm_max
      integer :: n_points
     

      logical :: if_stm

  real*8,     dimension(:,:,:),allocatable :: KS_vec
  integer :: cell_edge_steps_max
  complex*16, dimension(:,:,:),allocatable :: KS_vec_complex
  integer :: eig_dens_count, n_compute_basis_local, n_compute_basis
  integer :: n_points_per_task
  real*8,     dimension(:,:,:),allocatable :: densmat
  real*8,     dimension(:,:),allocatable :: densmat_tmp
  real*8,     dimension(:,:),allocatable :: densmat_sparse
  real*8,     dimension(:),allocatable :: densmat_sparse_tmp
  real*8,     dimension(:,:),allocatable :: density_matrix_con
  real*8,     dimension(:,:),allocatable :: work

  real*8,     dimension(:,:,:),allocatable :: KS_ev_compute
  complex*16, dimension(:,:,:),allocatable :: KS_ev_compute_complex

  real*8, dimension(:,:,:), allocatable ::  KS_orbital
  complex*16, dimension(:,:,:), allocatable ::  KS_orbital_complex
  real*8, dimension(:,:), allocatable :: local_rho
  real*8, dimension(:,:), allocatable :: local_free_rho
  real*8, dimension(:),   allocatable :: local_rho_temp
  real*8, dimension(:,:,:), allocatable ::  local_orb
  logical, dimension(:), allocatable :: point_in_cell

  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:,:), allocatable :: dir_tab_global
  real*8, dimension(:,:), allocatable :: dist_tab_sq_global

  real*8 cell_units(3)

! counters
      integer :: coord_x, coord_y, coord_z, i_z
      integer :: i_l
      integer :: i_m
      integer :: i_state
      integer :: i_point
     
      integer :: i_full_points
     
      integer :: i_spin = 1
      integer :: i_index, i_bas
      integer:: info, n_batch, n_batch_size, i_batch_z


!     begin work



  cell_units(1)= cell_edge_units(1)
  cell_units(2)= cell_edge_units(2)
  cell_units(3)= cell_edge_units(3)

  n_cube =1
  i_cell=1
  cell_edge_steps_max=0
  do i_z=1,3,1
     cell_edge_steps_max=MAX(cell_edge_steps(i_z),cell_edge_steps_max)
  enddo
 
  inv_bohr_3 = 1.0d0/(bohr**3)
  sqrt_inv_bohr_3 = sqrt(inv_bohr_3)
  
  l_ylm_max = l_wave_max

! First allocations

  if(.not. allocated( ylm_tab))then
     allocate( ylm_tab( (l_ylm_max+1)**2,n_centers_integrals),stat=info )
     if(info/=0)then
        write(use_unit,*)'Error in allocation: ylm_tab'
        stop
     end if
  end if


 if(.not. allocated( index_lm))then
     allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info) 
     if(info/=0)then
        write(use_unit,*)'Error in allocation: index_lm'
        stop
     end if
 end if

!     initialize index_lm
  i_index = 0
  do i_l = 0, l_ylm_max, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo

 do i_k = 1, n_k_points
  !     find the maximal occupation number
   do i_spin = 1, n_spin, 1
      ! initialize
      max_occ_number(i_spin,i_k) = 0
      do i_state = n_states, 1, -1
         if (dabs(occ_numbers(i_state,i_spin,i_k)).gt.0.d0) then
            max_occ_number(i_spin,i_k) = i_state
            exit
         endif
      enddo
   enddo
 end do

 if(.not. allocated(i_basis))then
     allocate(i_basis(n_centers_basis_T),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: i_basis'
        stop
     end if
 end if

 if(.not. allocated(densmat))then
    allocate(densmat(n_centers_basis_T,n_centers_basis_T,n_cube),stat=info)
    if(info/=0)then
       write(use_unit,*)'Error in allocation: densmat'
       stop
    end if
 end if
 if(.not. allocated(densmat_tmp))then
    allocate(densmat_tmp(n_centers_basis_T,n_centers_basis_T),stat=info)
    if(info/=0)then
       write(use_unit,*)'Error in allocation: densmat_tmp'
       stop
    end if
 end if
 if(.not. allocated(densmat_sparse))then
    allocate(densmat_sparse(n_hamiltonian_matrix_size,n_cube),stat=info)
    if(info/=0)then
       write(use_unit,*)'Error in allocation: densmat_sparse'
       stop
    end if
 end if
 if(.not. allocated(densmat_sparse_tmp))then
    allocate(densmat_sparse_tmp(n_hamiltonian_matrix_size),stat=info)
    if(info/=0)then
       write(use_unit,*)'Error in allocation: densmat_sparse'
       stop
    end if
 end if


 allocate( local_rho(1:cell_edge_steps(3),1:2) ) 
 allocate( local_free_rho(1:cell_edge_steps(3),1:2) ) 
 allocate( local_rho_temp(cell_edge_steps(3)))
 allocate( point_in_cell(1:cell_edge_steps(3))) 

! Calculate density matrices if needed
 n_cube=1
 n_spin=1
 densmat(1:n_centers_basis_T,1:n_centers_basis_T,i_cell)=0.0
 densmat_sparse(1:n_hamiltonian_matrix_size,i_cell)=0.0
 en_lower_limit=-1.d10
 en_upper_limit=1.d10

 do i_spin=1,n_spin,1
    call evaluate_densmat_part(KS_eigenvector, KS_eigenvector_complex, occ_numbers, &
         densmat_tmp, densmat_sparse_tmp, i_spin, &
         KS_eigenvalue, en_lower_limit, en_upper_limit)
    if(packed_matrix_format == PM_none)then
       densmat(1:n_centers_basis_T,1:n_centers_basis_T,i_cell)=&
            densmat(1:n_centers_basis_T,1:n_centers_basis_T,i_cell)+&
            densmat_tmp(1:n_centers_basis_T,1:n_centers_basis_T)
    else
       densmat_sparse(1:n_hamiltonian_matrix_size,i_cell)=&
            densmat_sparse(1:n_hamiltonian_matrix_size,i_cell)+&
            densmat_sparse_tmp(1:n_hamiltonian_matrix_size)
    endif
 enddo
 
! First, determine the size of the arrays that have to be allocated

 n_compute_basis_local=0
 n_compute_basis=0

 n_points_per_task=int(ceiling(dble(cell_edge_steps(3))/dble(n_tasks)))

 do coord_x =   1,cell_edge_steps(1),1

   do coord_y = 1,cell_edge_steps(2),1
      i_k = 0
      local_rho = 0.d0
      local_free_rho = 0.d0
      point_in_cell = .true.
      
      
      i_full_points = 0
      n_compute = 0
      i_basis = 0
      i_point = 0
      
         if (myid.lt.cell_edge_steps(3)) then
            do coord_z = myid*n_points_per_task+1,MIN((myid+1)*n_points_per_task,cell_edge_steps(3)),1
                       i_point = i_point+1
                       !    generate output grid 
                       coord_current(1)=cell_units(1)*(coord_x-1)-offset_coord(1)+cell_origin(1)
                       coord_current(2)=cell_units(2)*(coord_y-1)-offset_coord(2)+cell_origin(2)
                       coord_current(3)=cell_units(3)*(coord_z-1)-offset_coord(3)+cell_origin(3)


                       if(n_periodic > 0)then
                          coord_current_temp = coord_current
                          call map_to_center_cell(coord_current_temp(1:3) )   
                          if( abs(coord_current(1) -coord_current_temp(1))> 1e-3 .or. & 
                              abs(coord_current(2)-coord_current_temp(2))>1e-3 .or. &
                                abs(coord_current(3) - coord_current_temp(3))>1e-3 )then
!                             point_in_cell(i_point) = .false.
                          end if
                          coord_current = coord_current_temp
                       end if

                       if(point_in_cell(i_point))then

                       !     compute atom-centered coordinates of current integration point as viewed from all atoms
                       call tab_atom_centered_coords_p0 &
                            ( coord_current,  &
                            dist_tab_sq(1,i_point),  &
                            dir_tab(1,1,i_point), &
                            n_centers_basis_integrals, centers_basis_integrals )
 
                       call prune_basis_p0 &
                            (dist_tab_sq(1,i_point), &
                            n_compute_a, n_compute, i_basis,  &
                            n_centers_basis_T, n_centers_integrals, inv_centers_basis_integrals )
                       

                    end if ! point_in_cell 
                 enddo     ! coord_z

                 n_max_compute_dens = MAX(n_compute, n_max_compute_dens)
                 n_compute_basis_local=MAX(n_compute,n_compute_basis_local)

                 n_points = i_point
                 if (n_compute.gt.0) then

                    ! Determine all radial functions, ylm functions and their derivatives that
                    ! are best evaluated strictly locally at each individual grid point.
                    i_point = 0

                    do coord_z = myid*n_points_per_task+1,MIN((myid+1)*n_points_per_task,cell_edge_steps(3)),1

                          i_point = i_point+1
                          n_compute_atoms = 0
                          n_compute_fns = 0
                          i_basis_fns_inv = 0


                          if(point_in_cell(i_point))then

                          ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                          ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                          ! without any copying and without doing any unnecessary operations. 
                          ! The price is that the interface is no longer explicit in terms of physical 
                          ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
                          call prune_radial_basis_p0 &
                               ( dist_tab_sq(1,i_point),  &
                               dist_tab(1), &
                               dir_tab(1,1,i_point), &
                               n_compute_atoms, atom_index, atom_index_inv, &
                               n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                               i_atom_fns, spline_array_start, spline_array_end, &
                               n_centers_integrals, centers_basis_integrals)


                               n_max_compute_fns_dens = MAX(n_compute_fns, n_max_compute_fns_dens)

                       end if ! point_in_cell

                 enddo ! coord_z


              end if !n_compute
           end if

        enddo ! coord_y
     enddo ! coord_x 


 call sync_find_max(n_compute_basis_local,n_compute_basis)

! now allocate necessary quantities for this subroutine

 if(.not.allocated(radial_wave))then
     allocate(radial_wave(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave'
        stop
     end if
 else
     deallocate(radial_wave)
     allocate(radial_wave(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave'
        stop
     end if
 end if

 if(.not. allocated(wave))then
     allocate(wave(n_compute_basis,  cell_edge_steps(3)),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: wave'
        stop
     end if
 else
   deallocate(wave)
   allocate(wave(n_compute_basis,  cell_edge_steps(3)),stat=info)
   if(info/=0)then
      write(use_unit,*)'Error in allocation: wave'
      stop
   end if   
 end if

 if(.not. allocated(density_matrix_con))then
    allocate(density_matrix_con(n_compute_basis,n_compute_basis),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: density_matrix_con'
        stop
     end if
  else
     deallocate(density_matrix_con)
     allocate(density_matrix_con(n_compute_basis,n_compute_basis),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: density_matrix_con'
        stop
     end if
 end if

 if(.not. allocated(work))then
    allocate(work(n_compute_basis,cell_edge_steps(3)),stat=info)
    call check_allocation(info, 'work            ')
 end if
 

 if(real_eigenvectors)then

     if(.not. allocated( KS_ev_compute))then
        allocate( KS_ev_compute(n_states, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute'
           stop
        end if
     end if
     
     if(.not. allocated( KS_vec))then
        allocate( KS_vec(n_centers_basis_T, n_states, n_cube),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_vec'
           stop
        end if
     end if
     KS_vec = 0.d0

 else

     if(.not. allocated( KS_vec_complex))then
        allocate( KS_vec_complex(n_centers_basis_T,n_states,n_cube),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation:  KS_vec_complex'
           stop
        end if
     end if
     KS_vec_complex = (0.d0, 0.d0)


     if(.not. allocated( KS_ev_compute_complex))then
        allocate( KS_ev_compute_complex(n_states, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute_complex'
           stop
        end if
     end if


  end if

  if(real_eigenvectors)then
     allocate( KS_orbital(1,cell_edge_steps(3),n_cube) ) 
  else
     allocate( KS_orbital_complex(1,cell_edge_steps(3),n_cube) ) 
  end if



  ! ************************************************
  ! ************************************************

  do coord_x =   1,cell_edge_steps(1),1
if (myid == 0) print *, 'evaluating charge density of grid x:',coord_x
     do coord_y = 1,cell_edge_steps(2),1
        
           local_rho = 0.d0
           local_rho_temp = 0.d0
           local_free_rho = 0.d0
           work = 0.0
           if(real_eigenvectors) then
              KS_orbital=0.0
           else
              KS_orbital_complex=(0.d0,0.d0)
           endif

           point_in_cell = .true.
           
           
                 i_full_points = 0
                 n_compute = 0
                 i_basis = 0
                 i_point = 0
                 

                 if (myid.lt.cell_edge_steps(3)) then
                    do coord_z = myid*n_points_per_task+1,MIN((myid+1)*n_points_per_task,cell_edge_steps(3)),1

                       
                       i_point = i_point+1
                       !    generate output grid 
                       coord_current(1)=cell_units(1)*(coord_x-1)-offset_coord(1)+cell_origin(1)
                       coord_current(2)=cell_units(2)*(coord_y-1)-offset_coord(2)+cell_origin(2)
                       coord_current(3)=cell_units(3)*(coord_z-1)-offset_coord(3)+cell_origin(3)


                       if(n_periodic > 0)then
                          coord_current_temp = coord_current
                          call map_to_center_cell(coord_current_temp(1:3) )   
                          if( abs(coord_current(1) -coord_current_temp(1))> 1e-3 &
                              .or. abs(coord_current(2)-coord_current_temp(2))>1e-3 .or. &
                                abs(coord_current(3) - coord_current_temp(3))>1e-3 )then
!                             point_in_cell(i_point) = .false.
                          end if
                          coord_current = coord_current_temp
                       end if
                       

                       if(point_in_cell(i_point))then

                          !     compute atom-centered coordinates of current integration point as viewed from all atoms
                          call tab_atom_centered_coords_p0 &
                               ( coord_current,  &
                               dist_tab_sq(1,i_point),  &
                               dir_tab(1,1,i_point), &
                               n_centers_basis_integrals, centers_basis_integrals )
                          
                       !     determine which basis functions are relevant at current grid point,
                       !     and tabulate their indices

                       !     next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
 
                       call prune_basis_p0 &
                            (dist_tab_sq(1,i_point), &
                            n_compute_a, n_compute, i_basis,  &
                            n_centers_basis_T, n_centers_integrals, inv_centers_basis_integrals )

                    end if ! point_in_cell 
                 enddo     !          end loop over the z component

                 n_points = i_point
                 if (n_compute.gt.0) then

                    ! Determine all radial functions, ylm functions and their derivatives that
                    ! are best evaluated strictly locally at each individual grid point.
                    i_point = 0

                    do coord_z = myid*n_points_per_task+1,MIN((myid+1)*n_points_per_task,cell_edge_steps(3)),1

                          i_point = i_point+1
                          n_compute_atoms = 0
                          n_compute_fns = 0
                          i_basis_fns_inv = 0


                          if(.not. point_in_cell(i_point))then
                             
                             wave(:,i_point) = 0.d0
                          else

                          ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                          ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                          ! without any copying and without doing any unnecessary operations. 
                          ! The price is that the interface is no longer explicit in terms of physical 
                          ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
                          call prune_radial_basis_p0 &
                               ( dist_tab_sq(1,i_point),  &
                               dist_tab(1), &
                               dir_tab(1,1,i_point), &
                               n_compute_atoms, atom_index, atom_index_inv, &
                               n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                               i_atom_fns, spline_array_start, spline_array_end, &
                               n_centers_integrals, centers_basis_integrals)


                          ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                          ! for all atoms which are actually relevant

                          call tab_local_geometry_p0 &
                               ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                               dir_tab(1,1,i_point), dist_tab(1),  &
                               i_r(1) &
                               )
!                          local_free_rho(myid*n_points_per_task+i_point,1) = 0.d0


                          ! Now evaluate radial functions u(r) from the previously stored compressed spline arrays  

                          call evaluate_radial_functions_p0 &
                               (   spline_array_start, spline_array_end, &
                               n_compute_atoms, n_compute_fns,  &
                               dist_tab(1), i_r(1), &
                               atom_index, i_basis_fns_inv, &
                               basis_wave_ordered, radial_wave(1), &
                               .false., n_compute, n_max_compute_fns_dens   &
                               )

                          call tab_trigonom_p0 &
                               ( n_compute_atoms, dir_tab(1,1,i_point),  &
                               trigonom_tab(1,1) &
                               )

                          ! tabulate distance and Ylm's w.r.t. other atoms            
                          call tab_wave_ylm_p0 &
                               ( n_compute_atoms, atom_index,  &
                               trigonom_tab(1,1), l_shell_max,  &
                               l_ylm_max, &
                               ylm_tab(1,1) )


                          ! tabulate total wave function value for each basis function in all cases -
                          ! but only now are we sure that we have ylm_tab ...

                          call evaluate_waves_p0 &
                               (l_ylm_max, ylm_tab(1,1),  &
                               dist_tab(1),  &
                               index_lm, n_compute,  &
                               i_basis, radial_wave(1),  &
                               wave(1,i_point), n_compute_atoms,   &
                               atom_index_inv, n_compute_fns,  &
                               i_basis_fns_inv, n_max_compute_fns_dens  &
                               )
                       end if ! .not. point_in_cell

                    enddo !        end loop over z component

                 end if ! if (n_compute.gt.0)
                ! write(use_unit,*) 'free',local_free_rho(:,1)


                    if (n_compute.gt.0) then

                                if(packed_matrix_format /= PM_none )then
                                   call  prune_density_matrix_sparse(densmat_sparse(1,i_cell), density_matrix_con, &
                                        n_compute, i_basis)
                                else
                                   call  prune_density_matrix(densmat(1,1,i_cell), density_matrix_con, &
                                        n_compute, i_basis)
                                end if

                                call evaluate_KS_density_densmat(n_points, wave(1,1), n_compute,   &
                                     local_rho(myid*n_points_per_task+1:MIN((myid+1)* &
                                     n_points_per_task,cell_edge_steps(3)),1), &
                                     n_compute_basis, n_centers_basis_T, &
                                     density_matrix_con, work)

                    else
                                local_rho(myid*n_points_per_task+1:MIN((myid+1)* &
                                     n_points_per_task,cell_edge_steps(3)),1) = 0.d0
                    end if ! n_compute
                endif ! myid

           call sync_vector(local_rho(1,1),cell_edge_steps(3)*n_spin)
           if (real_eigenvectors) then
              call sync_vector(KS_orbital(1,1,1),cell_edge_steps(3))
           else
              call sync_vector_complex(KS_orbital_complex(1,1,1),cell_edge_steps(3))
           endif
           call sync_vector(local_free_rho(1,1),cell_edge_steps(3)*2)

           if(myid ==0)then

              ! plot data
                 do i_z = 1,cell_edge_steps(3),1
                             rho_e(coord_x,coord_y,i_z)=local_rho(i_z,1)
                 enddo ! i_z
              endif

           end do            !     end loop over y component integration loop
        end do        !     end loop over x component
 !end do ! i_spin

  !     finally, deallocate stuff.
  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if
  if (allocated(dir_tab_global)) then
     deallocate(dir_tab_global)
  end if
  if (allocated(dist_tab_sq_global)) then
     deallocate(dist_tab_sq_global)
  end if
  if (allocated(local_rho)) then
     deallocate(local_rho)
  end if
  if (allocated(KS_orbital)) then
     deallocate(KS_orbital)
  end if
  if (allocated(KS_orbital_complex)) then
     deallocate(KS_orbital_complex)
  end if
  if (allocated(local_orb)) then
     deallocate(local_orb)
  end if

  if (allocated(point_in_cell)) deallocate(point_in_cell)

  if(allocated(radial_wave))deallocate(radial_wave)

  if(allocated(wave)) deallocate(wave)

  if(allocated(i_basis)) deallocate(i_basis)

  if(allocated( KS_ev_compute)) deallocate( KS_ev_compute )

  if(allocated( KS_ev_compute_complex)) deallocate( KS_ev_compute_complex)

  if(allocated(local_free_rho)) deallocate(local_free_rho) 

  if(allocated( KS_vec))deallocate( KS_vec)

  if(allocated( KS_vec_complex))deallocate( KS_vec_complex)

  if(allocated( densmat))deallocate( densmat)

  if(allocated( densmat_tmp))deallocate( densmat_tmp)

  if(allocated( densmat_sparse))deallocate( densmat_sparse)

  if(allocated( densmat_sparse_tmp))deallocate( densmat_sparse_tmp)

  if(allocated( work))deallocate( work)

  if(allocated( density_matrix_con))deallocate( density_matrix_con)

end subroutine charge_density_even_grid0_p1

subroutine charge_density_grad_even_grid( )
!---------------------------------------------------------------------------------------------
! Obtain gradients of densities and squares of the gradients of densities
!---------------------------------------------------------------------------------------------
 
    use constants
    use dimensions
    use physics

    implicit none
   
    real*8   :: h(3)
    real*8   :: grad(3)
  
    integer  :: cell(3) 
    integer  :: nx, ny, nz
    integer  :: i, j, k, i_dir

    cell=cell_edge_steps
    cell = cell-1  
!    
    h=lat/dble(cell)
  
    nx=cell(1) 
    ny=cell(2) 
    nz=cell(3)

! calculate [ grad(n(r)) ]^2
  
    do i=1, nx
    do j=1, ny
    do k=1, nz
    
       if (i.gt.1) then
         grad(1)=(rho_e(i+1,j  ,k  )-rho_e(i-1,j  ,k  ))/(2*h(1))     
       else 
         grad(1)=(rho_e(2  ,j  ,k  )-rho_e(nx ,j  ,k  ))/(2*h(1))
       end if
     
       if (j.gt.1) then
         grad(2)=(rho_e(i  ,j+1,k  )-rho_e(i  ,j-1,k  ))/(2*h(2))       
       else 
         grad(2)=(rho_e(i  ,2  ,k  )-rho_e(i  ,ny ,k  ))/(2*h(2))
       end if
     
       if (k.gt.1) then
         grad(3)=(rho_e(i  ,j  ,k+1)-rho_e(i  ,j  ,k-1))/(2*h(3))       
       else 
         grad(3)=(rho_e(i  ,j  ,2  )-rho_e(i  ,j  ,nz ))/(2*h(3))
       end if
 
       do i_dir=1,3 
          rho_e_gradient(i_dir,i,j,k)=grad(i_dir)
       enddo 
       rho_e_2gradient(i,j,k)=grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3)
        
    end do
    end do
    end do
   
    do i=1,cell(1)+1
    do j=1,cell(2)+1
       do i_dir=1,3 
          rho_e_gradient(i_dir,i,j,cell(3)+1)=rho_e_gradient(i_dir,i,j,1) 
       enddo 
       rho_e_2gradient(i,j,cell(3)+1)=rho_e_2gradient(i,j,1) 
    end do  
    end do
    
    do j=1,cell(2)+1
    do k=1,cell(3)+1
       do i_dir=1,3 
          rho_e_gradient(i_dir,cell(1)+1,j,k)=rho_e_gradient(i_dir,1,j,k) 
       enddo 
       rho_e_2gradient(cell(1)+1,j,k)=rho_e_2gradient(1,j,k) 
    end do  
    end do
    
    do i=1,cell(1)+1
    do k=1,cell(3)+1
       do i_dir=1,3 
          rho_e_gradient(i_dir,i,cell(2)+1,k)=rho_e_gradient(i_dir,i,1,k) 
       enddo 
       rho_e_2gradient(i,cell(2)+1,k)=rho_e_2gradient(i,1,k) 
    end do  
    end do

end subroutine charge_density_grad_even_grid


subroutine density_interpolation_even_grid(xx,dens,gr2dens)
!---------------------------------------------------------------------------------------------
!
! BLEND_103 extends scalar point data into a cube.
!
!  Diagram:
!
!    011--------------111 
!      |               |
!      |               | 
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               | 
!      |               |
!    000--------------100 
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the 
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!---------------------------------------------------------------------------------------------
!
    use physics
    use dimensions

    implicit none

! input
    real*8 :: xx(3)                 ! input point

! output
    real*8 :: dens       ! interpolated value of density
    real*8 :: gr2dens    ! interpolated value of grad(n)**2

    integer :: i, j, k
    real*8  :: r,s,t
    real*8  :: V000, V100, V010, V001, V110, V101, V011, V111
    integer :: nx, ny, nz

integer:: coord_x, coord_y, coord_z

    nx=cell_edge_steps(1)-1
    ny=cell_edge_steps(2)-1
    nz=cell_edge_steps(3)-1

! location of indexes
    i = 1 + int(xx(1)*nx)
    j = 1 + int(xx(2)*ny)
    k = 1 + int(xx(3)*nz)
 
    r=(xx(1)*dble(nx)-dble(i-1))
    s=(xx(2)*dble(ny)-dble(j-1))
    t=(xx(3)*dble(nz)-dble(k-1))

! Density interpolation
    V000 = rho_e(i  ,j  ,k  )
    V100 = rho_e(i+1,j  ,k  )
    V010 = rho_e(i  ,j+1,k  )
    V110 = rho_e(i+1,j+1,k  )
    V001 = rho_e(i  ,j  ,k+1)
    V101 = rho_e(i+1,j  ,k+1)
    V011 = rho_e(i  ,j+1,k+1)
    V111 = rho_e(i+1,j+1,k+1)
   
    dens  = 								      &
      1	      * ( + V000 )						      &
    + r	      * ( - V000 + V100 )					      & 
    +	s     * ( - V000	+ V010 )				      &
    +	    t * ( - V000	       + V001 ) 			      &
    + r * s     * ( + V000 - V100 - V010  		     + V110 )	      &
    + r	  * t * ( + V000 - V100        - V001	     + V101 )		      &
    +	s * t * ( + V000	- V010 - V001 + V011 )			      &
    + r * s * t * ( - V000 + V100 + V010 + V001 - V011 - V101 - V110 + V111 )

! Squared gradient density interpolation
    V000 = rho_e_2gradient(i  ,j  ,k  )
    V100 = rho_e_2gradient(i+1,j  ,k  )
    V010 = rho_e_2gradient(i  ,j+1,k  )
    V110 = rho_e_2gradient(i+1,j+1,k  )
    V001 = rho_e_2gradient(i  ,j  ,k+1)
    V101 = rho_e_2gradient(i+1,j  ,k+1)
    V011 = rho_e_2gradient(i  ,j+1,k+1)
    V111 = rho_e_2gradient(i+1,j+1,k+1)
   
    gr2dens  = 								      &
      1	      * ( + V000 )						      &
    + r	      * ( - V000 + V100 )					      & 
    +	s     * ( - V000	+ V010 )				      &
    +	    t * ( - V000	       + V001 ) 			      &
    + r * s     * ( + V000 - V100 - V010  		     + V110 )	      &
    + r	  * t * ( + V000 - V100        - V001	     + V101 )		      &
    +	s * t * ( + V000	- V010 - V001 + V011 )			      &
    + r * s * t * ( - V000 + V100 + V010 + V001 - V011 - V101 - V110 + V111 )
   
if (gr2dens.gt.1e+20) then
      write(use_unit,*) "density gradient diverge -- problem"
      stop
endif

return

end subroutine density_interpolation_even_grid 

subroutine cube_cell_dimension()
!---------------------------------------------------------------------------------------------
! Define the dimensions of the cube cells
!---------------------------------------------------------------------------------------------
!
    use constants
    use dimensions
    use geometry 
    use runtime_choices

    implicit none

    integer i_dir, i, count
    integer nelect
    
    real*8, dimension(:,:), allocatable :: coord_tmp, coord_tmp1 
    lat_vec=0.0

    if (.not. flag_cell) cell_origin=0.0

! determining cell parameters

    if (n_periodic>1) then
      do i_dir=1,3
          lat_vec(i_dir,i_dir)=lattice_vector(i_dir,i_dir)
         if (cell_edge_steps(i_dir).eq.0) then
             cell_edge_steps(i_dir)=lattice_vector(i_dir,i_dir)/cell_edge_units(i_dir)+1.0
         else
             cell_edge_units(i_dir)=lattice_vector(i_dir,i_dir)/dble(cell_edge_steps(i_dir)-1)
         endif 
       enddo
    else
       do i_dir=1,3
          if ((cell_edge_steps(i_dir).eq.0).or.(cell_edge_units(i_dir).eq.0.0)) then
              write(use_unit,'(a)') 'please offer valid cell_edge_steps and cell_edge_units' 
              stop
          else 
              lat_vec(i_dir,i_dir)= cell_edge_units(i_dir)*dble(cell_edge_steps(i_dir)-1)
          endif 
       enddo
    endif

      offset_coord = (real(cell_edge_steps(1:3)-1)/2.d0)*cell_edge_units

    do i_dir=1,3
       lat(i_dir)= &
       dsqrt(lat_vec(1,i_dir)**2d0+lat_vec(2,i_dir)**2d0+lat_vec(3,i_dir)**2d0)
    enddo

end subroutine cube_cell_dimension


subroutine integrand(ndim, x, ncomp, func)
!---------------------------------------------------------------------------------------------
! This subroutine gives the integrand function for Cuba program
!
! Integrand = n(r1)*phi(r1,r2)*n(r2)
!
! Input:   ndim    : number of dimensions of the integral
!          x       : point in ndim space
!                   ( x(i) in [0,1] )
!          ncomp   : number of components of the integral
! Output:  func    : value of the integrant at x
!---------------------------------------------------------------------------------------------
!
    use physics
    use constants
    use runtime_choices
    use dimensions
    use grids
    use geometry 
    use species_data
    use pbc_lists

    implicit none

! input
    integer :: ndim         
    real*8  :: x(ndim)      
    integer :: ncomp        

! output
    real*8  :: func(ncomp)  

    real*8  :: coord1(3), coord2(3), coord1b(3), coord2b(3), dr(3)
    real*8  :: distance, q01, q02
    real*8  :: delta, phi
    real*8  :: density1, density2
    real*8  :: grad2den1, grad2den2
    real*8  :: weight1, weight2
    
    real*8    :: volume, const
    real*8    :: tmp

    integer   :: i_point1, i_point2  
    integer   :: i_atom1, i_atom2
    integer   :: i_radial1, i_radial2
    integer   :: i_angular1, i_angular2
    integer   :: i_species1, i_species2

    func=0.0d0

    i_point1     = 1 + int(x(1)*(n_full_points-1)) 
    i_point2     = 1 + int(x(2)*(n_full_points-1)) 

    if ((partition_tab(i_point1).gt.0.0d0).and.(partition_tab(i_point2).gt.0.0d0)) then
       continue
    else
       return
    endif

    coord1(:) = batches(i_point1)%points(1)%coords(:)
    coord2(:) = batches(i_point2)%points(1)%coords(:)

    dr        = coord2-coord1
    distance  = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

    call charge_density(i_point1,density1,grad2den1)
    call charge_density(i_point2,density2,grad2den2)

    q01=wave_vector(density1,grad2den1)
    q02=wave_vector(density2,grad2den2)

    call calc_kernel(q01,q02,distance,phi)

    const=0.5d0

    weight1=partition_tab(i_point1)
    weight2=partition_tab(i_point2)
    weight1=1.0
    weight2=1.0

    func(ncomp) = const*density1*weight1*phi*density2*weight2
          
return

end subroutine integrand


subroutine integrand_even_spacing_grid(ndim, x, ncomp, func)
!---------------------------------------------------------------------------------------------
! This subroutine gives the integrand function of even spacing grids for Cuba program
!
! Integrand = n(r1)*phi(r1,r2)*n(r2)
!
! Input:   ndim    : number of dimensions of the integral
!          x       : point in ndim space
!                   ( x(i) in [0,1] )
!          ncomp   : number of components of the integral
! Output:  func    : value of the integrant at x
!---------------------------------------------------------------------------------------------
!
    use physics
    use constants
    use runtime_choices
!
    implicit none
!
! input
    integer :: ndim         
    real*8  :: x(ndim)      
    integer :: ncomp        

! output
    real*8  :: func(ncomp)  

    real*8  :: xsc(6), xpr(6), r1(3), r2(3), dr(3)
    real*8  :: distance, q01, q02
    real*8  :: d1, d2, D, delta, phi
    real*8  :: density1, density2
    real*8  :: grad2den1, grad2den2
    integer :: i, j, k 
    
    real*8    :: volume, const
    real*8    :: tmp
   
    real*8    :: small, c2, CvdW 
    real*8    :: Zab

    parameter (small = 1.d-8)
    parameter (c2  = 3.093667726280136d0)    ! c2=(3 pi^2)^(1/3)
    parameter (Zab=-0.8491d0  )

    func=0.0d0

    xsc(1)=x(1);    xsc(4)=(2.0*dble(cell_size_vdw(1))+1.0)*x(4)-dble(cell_size_vdw(1))
    xsc(2)=x(2);    xsc(5)=(2.0*dble(cell_size_vdw(2))+1.0)*x(5)-dble(cell_size_vdw(2))
    xsc(3)=x(3);    xsc(6)=(2.0*dble(cell_size_vdw(3))+1.0)*x(6)-dble(cell_size_vdw(3))

    do i=1,6
      if   ( xsc(i).lt.0.0d0 )     then
             xpr(i)=xsc(i)-dble(int(xsc(i)))+1.0d0
      else
             xpr(i)=xsc(i)-dble(int(xsc(i)))
      end if
    end do

    call density_interpolation_even_grid(xpr(1:3),density1,grad2den1)
!
    call density_interpolation_even_grid(xpr(4:6),density2,grad2den2)

! calculation real space radius vectors

    r1(:)=xsc(1)*lat_vec(:,1)+xsc(2)*lat_vec(:,2)+xsc(3)*lat_vec(:,3)
    r2(:)=xsc(4)*lat_vec(:,1)+xsc(5)*lat_vec(:,2)+xsc(6)*lat_vec(:,3)

    dr   = r2 - r1
    distance  = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

    q01=wave_vector(density1,grad2den1)
    q02=wave_vector(density2,grad2den2)

    call calc_kernel(q01,q02,distance,phi)

    volume=cal_cell_volume(lat_vec)

    const=0.5d0*(2.0*cell_size_vdw(1)+1.0)*(2.0*cell_size_vdw(2)+1.0)* &
                (2.0*cell_size_vdw(3)+1.0)*volume*volume

    func(ncomp) = const*density1*phi*density2

return

end subroutine integrand_even_spacing_grid


subroutine calc_kernel(q01,q02,dis,phi)
!---------------------------------------------------------------------------------------------
! This subroutine gives parameters for calculating the integrand kernel
!---------------------------------------------------------------------------------------------
!
    use physics
    use constants
    use runtime_choices
!
    implicit none
!
! input
    real*8  :: q01, q02, dis

! output
    real*8  :: phi

!local variables
    real*8  :: delta, D
    real*8  :: d1, d2    
    real*8    :: phi11,phi12,phi21,phi22,phi1,phi2
    real*8    :: D0,delta0
    real*8    :: tmp
    integer   :: iD,idelta

    real*8    :: CvdW 
    parameter (CvdW = 32.6650486837726d0)    ! CvdW = 12 (4*pi/9)^3

    idelta=0

    d1    = dis*q01
    d2    = dis*q02
    D     = 0.5d0*(d1 + d2)

    if ((d1.lt.eps).and.(d2.lt.eps)) then
        delta=0.d0
    else
        delta = (d1 - d2)/(d1 + d2)
    end if

    delta = dabs(delta)
  
    if (delta.lt.deltamin .or. delta.gt.deltamax) then
        write(21,*) 'delta out of range',delta
        stop
    end if 
    if (dcutoff.lt.dmin) then
        write(21,*) 'dcutoff too small! dmin =',dmin
        stop
    end if

   
! use empirical relation to determine whether the asymptotic for phi can be used
! (compare the plots for the exact phi and its asymptotic form)
    
 
      if ( ( (0.167*D - 0.1296*D*delta).gt.1.0d0 ) .or. ( D.gt.dmax ) ) then
      
      ! use asymptotic form of Phi(D,delta)     
      
           phi = - CvdW/( (d1*d2)**2d0*(d1**2d0+d2**2d0) )
           
      else if (D .le. dcutoff) then
  
           phi = kernel(idelta,1)

      else
  
           ! bilinear interpolation between tabulated values
           iD     = 1 + int((D     - Dmin)/dstep)
           idelta = 1 + int((delta - deltamin)/deltastep)
            
           D0     = Dmin + dble(iD-1)*dstep
           delta0 = deltamin + dble(idelta-1)*deltastep
           
           phi11  = kernel(idelta    ,iD)
           phi12  = kernel(idelta + 1,iD)
           phi21  = kernel(idelta    ,iD + 1)
           phi22  = kernel(idelta + 1,iD + 1)                  
            
           phi1   = phi11 + (D - D0)*(phi21 - phi11)/dstep
           phi2   = phi12 + (D - D0)*(phi22 - phi12)/dstep     
           
           phi    = phi1  + (delta - delta0)*(phi2 - phi1)/deltastep
       
      end if
return

end subroutine calc_kernel


subroutine read_kernel_data ()
!---------------------------------------------------------------------------------------------
! Read the kernal data
!---------------------------------------------------------------------------------------------

    use constants
    
    implicit none
 
    integer          :: i
    integer          :: in 
    integer          :: nd, ndelta

    real*8           :: d
  
    parameter (in=10) 
  
    open(in,file=kernelfile,status='old')
    read(in,*) dmin,dmax,dstep
    read(in,*) deltamin,deltamax,deltastep 
    read(in,*) nd, ndelta
    
    allocate(kernel(ndelta,nd))
   
    read(in,*) kernel
    
!    do i = 1, nd
!        d        = dmin + (i-1)*dstep
!        kernel(:,i) = kernel(:,i)/(4.0d0*pi*d*d) 
!    end do
    
    close(in)

end subroutine read_kernel_data


subroutine cleanup_ll_vdw()
!---------------------------------------------------------------------------------------------
!  Deallocate arrays
!---------------------------------------------------------------------------------------------
!
    implicit none

    if ( allocated (cell_origin) ) then
        deallocate (cell_origin)
    end if
    if ( allocated (kernel) ) then
        deallocate (kernel)
    end if
     

end subroutine cleanup_ll_vdw


real*8 function cal_cell_volume(vec)
!---------------------------------------------------------------------------------------------
! function calculating the volume of the cell
!---------------------------------------------------------------------------------------------

    use constants

    implicit none

    real*8  :: vec(3,3)
    real*8  :: crossab(3)
    real*8  :: vol
 
! calculate unit cell volume  
    crossab(1) = vec(2,1)*vec(3,2) - vec(3,1)*vec(2,2)
    crossab(2) = vec(3,1)*vec(1,2) - vec(1,1)*vec(3,2)
    crossab(3) = vec(1,1)*vec(2,2) - vec(2,1)*vec(1,2)

    vol        = crossab(1)*vec(1,3) + crossab(2)*vec(2,3) + crossab(3)*vec(3,3)

    cal_cell_volume=vol     

    return

end function cal_cell_volume


real*8 function wave_vector(density,grad2den)
!---------------------------------------------------------------------------------------------
! function calculating the wave vector
!---------------------------------------------------------------------------------------------
!
    use physics
    use constants
    use dimensions
    use runtime_choices
!
    implicit none

! input 
    real*8  :: density
    real*8  :: dens(n_spin)
    real*8  :: grad2den
!
    real*8  :: lda_en_x, lda_en_c
    real*8  :: kf, ec, q0
    
    real*8    :: c2, CvdW 
    real*8    :: Zab

    parameter (c2  = 3.093667726280136d0)    ! c2=(3 pi^2)^(1/3)
    parameter (CvdW = 32.6650486837726d0)    ! CvdW = 12 (4*pi/9)^3
    parameter (Zab=-0.8491d0  )

    dens=0.0
    dens(1)=density
    call local_lda_x_c_energy_even_spacing_grid (dens,lda_en_x,lda_en_c)

    kf = c2 * density**third
    q0 = kf - 4.0d0*pi/3.0d0*lda_en_c - &
          Zab/(36.0d0*kf)*grad2den/(density*density)

    wave_vector=q0

    return

end function wave_vector

subroutine trans_cell(coord_tmp,center)
!---------------------------------------------------------------------------------------------
! translating the center of the cell to the given center
!---------------------------------------------------------------------------------------------
!
    use dimensions
    use geometry 
!
    implicit none

!input
    real*8, dimension(3) :: center

!output
    real*8, dimension(3,n_atoms) :: coord_tmp

!local variable
    real*8, dimension(3) :: center0, trans 

!counter
    integer*4 :: i_atom, i_dir

    center0=0.0
    do i_atom=1, n_atoms
       center0=center0+coords(:,i_atom)
    enddo
    center0=center0/dble(n_atoms)
    trans=center-center0
!
    do i_atom=1, n_atoms
       do i_dir=1,3
          coord_tmp(i_dir,i_atom)=coords(i_dir,i_atom)+trans(i_dir)
       enddo
    enddo

    return

end subroutine trans_cell
    
subroutine get_center(center)
!---------------------------------------------------------------------------------------------
! translating the center of the cell to the given center
!---------------------------------------------------------------------------------------------
!
    use dimensions
    use geometry 
!
    implicit none

!output
    real*8, dimension(3) :: center

!local variable
    real*8, dimension(3) :: center0, trans 

!counter
    integer*4 :: i_atom, i_dir

    center=0.0
    do i_atom=1, n_atoms
       center=center+coords(:,i_atom)
    enddo
    center=center/dble(n_atoms)

    return

end subroutine get_center 

subroutine charge_density_trans( )
!---------------------------------------------------------------------------------------------
! translating the center of the cell to the given center
!---------------------------------------------------------------------------------------------
!
    use dimensions
    use physics
!
    implicit none

    real*8 origin(3)
    real*8, dimension(:,:,:), allocatable :: current_rho1, current_rho2

! counter
    integer*4 i_dir, i, j, k
    integer*4 nx, ny, nz, di, dj, dk 

!    if (n_periodic>1) then
!       origin=-(offset_coord-cell_origin)
!    else
       origin=cell_origin
!    endif
!
    nx=cell_edge_steps(1)
    ny=cell_edge_steps(2)
    nz=cell_edge_steps(3)
!
    allocate(current_rho1(nx,ny,nz),current_rho2(nx,ny,nz))
!    
    current_rho1=rho_e
!
    if (dabs(origin(1))>eps) then  
     di=int(dble(nx)*origin(1)/lat(1))
     do i=1,nx
     do j=1,ny
     do k=1,nz
        if (di>0) then
           if ( (i-di) > 0 ) then
              current_rho1(i,j,k)=rho_e(i-di,j,k)
           else
              current_rho1(i,j,k)=rho_e(i-di+nx,j,k)
           end if
        else
           if ( (i-di) > nx ) then
              current_rho1(i,j,k)=rho_e(i-di,j,k)
           else
              current_rho1(i,j,k)=rho_e(i-di-nx,j,k)
           end if
        end if
     end do
     end do
     end do
    endif
    
    current_rho2=current_rho1
    
    if (dabs(origin(2))>eps) then   
     dj=int(dble(ny)*origin(2)/lat(2))
     do i=1,nx
     do j=1,ny
     do k=1,nz
        if (dj>0) then
           if ( (j-dj) > 0 ) then
              current_rho1(i,j,k)=current_rho2(i,j-dj,k)
           else
              current_rho1(i,j,k)=current_rho2(i,j-dj+ny,k)
           end if
        else
           if ( (j-dj) > ny ) then
              current_rho1(i,j,k)=current_rho2(i,j-dj,k)
           else
              current_rho1(i,j,k)=current_rho2(i,j-dj-ny,k)
           end if
        end if
     end do
     end do
     end do
    endif
    
    current_rho2=current_rho1
    
    if (dabs(origin(3))>eps) then  
     dk=int(dble(nz)*origin(3)/lat(3))
     do i=1,nx
     do j=1,ny
     do k=1,nz
        if (dk>0) then
           if ( (k-dk) > 0 ) then
              current_rho1(i,j,k)=current_rho2(i,j,k-dk)
           else
              current_rho1(i,j,k)=current_rho2(i,j,k-dk+nz)
           end if
        else
           if ( (k-dk) > nz ) then
              current_rho1(i,j,k)=current_rho2(i,j,k-dk)
           else
              current_rho1(i,j,k)=current_rho2(i,j,k-dk-nz)
           end if
        end if
     end do
     end do
     end do
    endif
    
     rho_e=current_rho1
     deallocate(current_rho1,current_rho2)

end subroutine charge_density_trans

!---------------------------------------------------------------------------------------------
end module ll_vdwdf
!---------------------------------------------------------------------------------------------
    
