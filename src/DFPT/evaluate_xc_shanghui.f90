!****s* FHI-aims/evaluate_xc
!  NAME
!  evaluate_xc
!  SYNOPSIS
      subroutine evaluate_xc_shanghui &
      ( rho, rho_gradient, en_density_xc, &
       en_density_x, en_density_c, &
        local_xc_derivs, xc_gradient_deriv,dVxc_drho,coord_current &
      )

!  PURPOSE
!  Subroutine evaluate_xc evaluates the exchange correlation potential
!  and energy contribution for a given density at one integration point.
!
!  VB: The present version supports LDA and gradient functionals.
!
!  USES
!
      use runtime_choices
      use dimensions
      use xc
      use constants
      use localorb_io, only: use_unit

     !-------------shanghui add for libxc---------------------
     ! use xc_f90_types_m
     ! use xc_f90_lib_m 
     !-------------shanghui end add for libxc-----------------

      implicit none

!  ARGUMENTS

      real*8 rho(n_spin)
      real*8 rho_gradient(3,n_spin)
      real*8 en_density_xc
      real*8 en_density_x(n_spin)
      real*8 en_density_c
      real*8 local_xc_derivs(n_spin)
      real*8 xc_gradient_deriv(3,n_spin)
      real*8  dVxc_drho(n_spin)
      real*8, optional :: coord_current(3)  !SAG

!  INPUTS
!   o flag_xc : Determines type of XC functional, see read_control.f
!   o rho     : Density at current point, given as spin up and spin down if spin-polarized
!   o rho_gradient : Density gradient at current point
!  OUTPUT
!   o en_density_xc : Exchange-correlation energy density at current point
!   o local_xc_derivs : Local parts of the XC potential:
!               Partial derivative of the exchange-correlation
!               energy functional by the density plus spin-weighted
!               partial derivative of the XC energy functional by spin
!   o xc_gradient_deriv : Partial derivative of the exchange-correlation
!               energy functional by the modulus square of the density gradient
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

      logical :: xc_undefined

!     counters

      integer :: i_spin
   
     !-------------shanghui add for siesta_pbe test---------
     ! real*8   DECDD(n_spin), DECDGD(3,n_spin),   & 
     !          DEXDD(n_spin), DEXDGD(3,n_spin)

     !-------------shanghui end add for siesta_pbc test-----


     !-------------shanghui add for DFT-repository test-------
      integer ideriv, npt
      real*8  rhoa, rhob, sigma_aa, sigma_bb, sigma_ab 

     !-----exchange------->
      real*8   zk_x,vrhoa_x,vrhob_x,vsigma_aa_x,vsigma_bb_x,vsigma_ab_x
      real*8   v2rhoa2_x,v2rhob2_x,v2rhoab_x
      real*8   v2rhoasigma_aa_x,v2rhoasigma_ab_x,v2rhoasigma_bb_x
      real*8   v2rhobsigma_bb_x,v2rhobsigma_ab_x,v2rhobsigma_aa_x
      real*8   v2sigma_aa2_x,v2sigma_aaab_x,v2sigma_aabb_x
      real*8   v2sigma_ab2_x,v2sigma_abbb_x,v2sigma_bb2_x

     !----correlation----->
      real*8   zk_c,vrhoa_c,vrhob_c,vsigma_aa_c,vsigma_bb_c,vsigma_ab_c
      real*8   v2rhoa2_c,v2rhob2_c,v2rhoab_c
      real*8   v2rhoasigma_aa_c,v2rhoasigma_ab_c,v2rhoasigma_bb_c
      real*8   v2rhobsigma_bb_c,v2rhobsigma_ab_c,v2rhobsigma_aa_c
      real*8   v2sigma_aa2_c,v2sigma_aaab_c,v2sigma_aabb_c
      real*8   v2sigma_ab2_c,v2sigma_abbb_c,v2sigma_bb2_c
     !-------------shanghui end add for DFT-repository test-------






     !-------------shanghui add for libxc---------------------
     ! TYPE(xc_f90_pointer_t) :: xc_func
     ! TYPE(xc_f90_pointer_t) :: xc_info
     
 
     ! integer ::  func_id  ! here 1 is slater-x, 9 is pz-correlation
     ! real*8 zk(n_spin),sigma(n_spin),    &   
     !        ex(n_spin),vx(n_spin),fx(n_spin),kx(n_spin),  &
     !        ec(n_spin),vc(n_spin),fc(n_spin),kc(n_spin) 
     !
     ! real*8 vsigma_x(n_spin), v2rho2_x(n_spin), & 
     !        v2rhosigma_x(n_spin), v2sigma2_x(n_spin)
     ! real*8 vsigma_c(n_spin), v2rho2_c(n_spin), & 
     !        v2rhosigma_c(n_spin), v2sigma2_c(n_spin)

 
      !------shanghui test for libxc-interface-------PZ-LDA 
      !rho(1)=0.1d0      
    
      !func_id=1
      !call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
      !call xc_f90_lda(xc_func, 1, rho(1),ex(1),vx(1),fx(1),kx(1))
      !func_id=9
      !call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
      !call xc_f90_lda(xc_func, 1, rho(1),ec(1),vc(1),fc(1),kc(1))

      !write(use_unit,*) 'libxc:',rho,ex+ec,vx+vc,fx+fc


      !------shanghui test for libxc-interface-------GGA-PBE 
     ! rho(1)=0.1d0
     ! rho_gradient(1,1)=0.1d0 
     ! rho_gradient(2,1)=0.2d0
     ! rho_gradient(3,1)=0.3d0     
     ! sigma(1)=rho_gradient(1,1)*rho_gradient(1,1)+  & 
     !          rho_gradient(2,1)*rho_gradient(2,1)+  & 
     !          rho_gradient(3,1)*rho_gradient(3,1)      

     ! func_id=101
     ! call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
     ! call xc_f90_gga(xc_func, 1, rho(1), sigma(1),ex(1),vx(1), &  
     !      vsigma_x(1), v2rho2_x(1),v2rhosigma_x(1), v2sigma2_x(1))
     ! func_id=130
     ! call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
     ! call xc_f90_gga(xc_func, 1, rho(1), sigma(1),ec(1),vc(1), &  
     !      vsigma_c(1), v2rho2_c(1),v2rhosigma_c(1), v2sigma2_c(1))

     ! write(use_unit,*) 'libxc-pbe:',rho,ex+ec,vx+vc, & 
     !             (vsigma_x(1)+vsigma_c(1))*2.0d0*rho_gradient(1,1), & 
     !             (vsigma_x(1)+vsigma_c(1))*2.0d0*rho_gradient(2,1), & 
     !             (vsigma_x(1)+vsigma_c(1))*2.0d0*rho_gradient(3,1)
     !-------------shanghui end add for libxc---------------------

 


!  begin work

      !---------shanghui add  to initial dVxc_drho to 0 ----
        dVxc_drho(1:n_spin) = 0.0d0
      !---------shanghui end add to initial dVxc_drho to 0 ----

        xc_undefined = .true.

        do i_spin = 1, n_spin, 1
          ! This is true if both densities are below zero
          xc_undefined = xc_undefined .and. (rho(i_spin).le.0.d0)
        enddo

        do i_spin = 1, n_spin, 1
          ! This is true if one of the densities is less than zero
          xc_undefined = xc_undefined .or. (rho(i_spin).lt.0.d0)
        enddo

        if (xc_undefined) then

          en_density_xc = 0.d0
          en_density_x = 0.d0
          en_density_c = 0.d0

          local_xc_derivs = 0.d0
          xc_gradient_deriv = 0.d0

        else

!         calculate xc-potential, thus initializing new potential
          select case(flag_xc)
           case(0)
!        Hartree-Fock calculation, set the XC contribution to zero here.
            en_density_xc = 0.d0
            en_density_x = 0.d0
            en_density_c = 0.d0
            local_xc_derivs = 0.d0
            xc_gradient_deriv = 0.d0

          case(1)
!           Hybrid-PBE0 functional: here only 3/4 of the PBE exchange are included.
            call xc_partials_pbe0 &
            ( rho, rho_gradient, en_density_xc, &
              en_density_x, en_density_c, local_xc_derivs, &
              xc_gradient_deriv )

          case(3)
!           Perdew-Zunger LDA
            
            call xcpot_pz_lda_shanghui(rho, en_density_xc, en_density_x, en_density_c,  &
                         local_xc_derivs, dVxc_drho )
          !  write(use_unit,*) rho, en_density_xc,local_xc_derivs,dVxc_drho
          !  stop


          case(4)
!           PW91_gga
!           FIXME: The following routine does not use the Hesse matrix of
!           the density. Therefore, it returns only the energy density and
!           its partial derivatives. To compute the local potential explicitly,
!           use another routine (which actually requires the Hessian of the density).
!           Johan has already supplied that subroutine. Need to integrate it here.
            call xc_partials_pw91gga &
            ( rho, rho_gradient, en_density_xc,  en_density_x, &
              en_density_c, local_xc_derivs, xc_gradient_deriv )

         case(5)  
            !PBE with VDW. SAG
            
            
           call xc_partials_pbe_vdw &
            ( rho, rho_gradient, en_density_xc,  en_density_x, &
              en_density_c, local_xc_derivs, xc_gradient_deriv, &
              coord_current)
      

          case(6)
!           Perdew-Burke-Ernzerhoff PRL 77, 3865 (1996)
!           FIXME: The following routine does not use the Hesse matrix of
!           the density. Therefore, it returns only the energy density and
!           its partial derivatives. To compute the local potential explicitly,
!           use another routine (which actually requires the Hessian of the density).
!           Johan has already supplied that subroutine. Need to integrate it here.
            call xc_partials_pbe &
            ( rho, rho_gradient, en_density_xc,  en_density_x, &
              en_density_c, local_xc_derivs, xc_gradient_deriv )
           !write(use_unit,*) 'aims-PBE'
           !write(use_unit,'(7F15.10)')   & 
           !           rho, en_density_x,en_density_c, & 
           !           local_xc_derivs,xc_gradient_deriv*2.0d0  
           !            ! 2*gradient_n*(d(E)/d(grad_n.grad_n) 


           !-------------shanghui begin DFT_repository-------------------

           if(n_spin.eq.1) then ! closed shell 
           
             ideriv = 2
             npt = 1
 
             rhoa = rho(1)/2.d0
             rhob = rho(1)/2.d0
             
             sigma_aa =  (rho_gradient(1,1)*rho_gradient(1,1)+  &
                         rho_gradient(2,1)*rho_gradient(2,1)+  &
                         rho_gradient(3,1)*rho_gradient(3,1))/4.0d0
             sigma_bb =  (rho_gradient(1,1)*rho_gradient(1,1)+  &
                         rho_gradient(2,1)*rho_gradient(2,1)+  &
                         rho_gradient(3,1)*rho_gradient(3,1))/4.0d0
             sigma_ab =  (rho_gradient(1,1)*rho_gradient(1,1)+  &
                         rho_gradient(2,1)*rho_gradient(2,1)+  &
                         rho_gradient(3,1)*rho_gradient(3,1))/4.0d0

             call  uks_x_pbe                                         &  
            (ideriv,npt,rhoa,rhob,sigma_aa,sigma_bb,sigma_ab,       &  
             zk_x,vrhoa_x,vrhob_x,vsigma_aa_x,vsigma_bb_x,vsigma_ab_x, &
             v2rhoa2_x,v2rhob2_x,v2rhoab_x,                         &
             v2rhoasigma_aa_x,v2rhoasigma_ab_x,v2rhoasigma_bb_x,    & 
             v2rhobsigma_bb_x,v2rhobsigma_ab_x,v2rhobsigma_aa_x,    &
             v2sigma_aa2_x,v2sigma_aaab_x,v2sigma_aabb_x,           &
             v2sigma_ab2_x,v2sigma_abbb_x,v2sigma_bb2_x)
             
             call  uks_c_pbe                                       &
            (ideriv,npt,rhoa,rhob,sigma_aa,sigma_bb,sigma_ab,       &
             zk_c,vrhoa_c,vrhob_c,vsigma_aa_c,vsigma_bb_c,vsigma_ab_c,&
             v2rhoa2_c,v2rhob2_c,v2rhoab_c,                          &
             v2rhoasigma_aa_c,v2rhoasigma_ab_c,v2rhoasigma_bb_c,     &
             v2rhobsigma_bb_c,v2rhobsigma_ab_c,v2rhobsigma_aa_c,     &
             v2sigma_aa2_c,v2sigma_aaab_c,v2sigma_aabb_c,            &
             v2sigma_ab2_c,v2sigma_abbb_c,v2sigma_bb2_c)


             !------(0) energy density: exc----------------
             en_density_x = zk_x/(rhoa+rhob)
             en_density_c = zk_c/(rhoa+rhob)

             !------(1.1) d exc/d rho---------------------- 
             local_xc_derivs =  vrhoa_x + vrhob_c

             !------(1.2) 0.5*d exc/d gradient_rho---------
             xc_gradient_deriv(1,1) =   &     ! x-coord
             ( (vsigma_aa_x)*rho_gradient(1,1) +            &
               (vsigma_ab_x)*rho_gradient(1,1)/2.0d0 +      & 
               (vsigma_aa_c)*rho_gradient(1,1) +            &  
               (vsigma_ab_c)*rho_gradient(1,1)/2.0d0 )/2.0d0     
             xc_gradient_deriv(2,1) =   &     ! y-coord
             ( (vsigma_aa_x)*rho_gradient(2,1) +            &
               (vsigma_ab_x)*rho_gradient(2,1)/2.0d0 +      & 
               (vsigma_aa_c)*rho_gradient(2,1) +            &  
               (vsigma_ab_c)*rho_gradient(2,1)/2.0d0 )/2.0d0     
             xc_gradient_deriv(3,1) =   &     ! z-coord
             ( (vsigma_aa_x)*rho_gradient(3,1) +            &
               (vsigma_ab_x)*rho_gradient(3,1)/2.0d0 +      & 
               (vsigma_aa_c)*rho_gradient(3,1) +            &  
               (vsigma_ab_c)*rho_gradient(3,1)/2.0d0 )/2.0d0     

           !write(use_unit,*) 'aims-PBE-DFT-repository'
           !write(use_unit,'(7F15.10)')   & 
           !           rho, en_density_x,en_density_c, & 
           !           local_xc_derivs,xc_gradient_deriv*2.0d0  
      

          else  ! open shell 
            
             ideriv = 2
             npt = 1
 
             rhoa = rho(1)
             rhob = rho(2)
             
             sigma_aa =  rho_gradient(1,1)*rho_gradient(1,1)+  &
                         rho_gradient(2,1)*rho_gradient(2,1)+  &
                         rho_gradient(3,1)*rho_gradient(3,1)
             sigma_bb =  rho_gradient(1,2)*rho_gradient(1,2)+  &
                         rho_gradient(2,2)*rho_gradient(2,2)+  &
                         rho_gradient(3,2)*rho_gradient(3,2)
             sigma_ab =  rho_gradient(1,1)*rho_gradient(1,2)+  &
                         rho_gradient(2,1)*rho_gradient(2,2)+  &
                         rho_gradient(3,1)*rho_gradient(3,2)

             call  uks_x_pbe                                         &  
            (ideriv,npt,rhoa,rhob,sigma_aa,sigma_bb,sigma_ab,       &  
             zk_x,vrhoa_x,vrhob_x,vsigma_aa_x,vsigma_bb_x,vsigma_ab_x, &
             v2rhoa2_x,v2rhob2_x,v2rhoab_x,                         &
             v2rhoasigma_aa_x,v2rhoasigma_ab_x,v2rhoasigma_bb_x,    & 
             v2rhobsigma_bb_x,v2rhobsigma_ab_x,v2rhobsigma_aa_x,    &
             v2sigma_aa2_x,v2sigma_aaab_x,v2sigma_aabb_x,           &
             v2sigma_ab2_x,v2sigma_abbb_x,v2sigma_bb2_x)
             
              call  uks_c_pbe                                       &
            (ideriv,npt,rhoa,rhob,sigma_aa,sigma_bb,sigma_ab,       &
             zk_c,vrhoa_c,vrhob_c,vsigma_aa_c,vsigma_bb_c,vsigma_ab_c,&
             v2rhoa2_c,v2rhob2_c,v2rhoab_c,                          &
             v2rhoasigma_aa_c,v2rhoasigma_ab_c,v2rhoasigma_bb_c,     &
             v2rhobsigma_bb_c,v2rhobsigma_ab_c,v2rhobsigma_aa_c,     &
             v2sigma_aa2_c,v2sigma_aaab_c,v2sigma_aabb_c,            &
             v2sigma_ab2_c,v2sigma_abbb_c,v2sigma_bb2_c)

           endif


           !-------------shanghui end DFT_repository-------------------


             

 
          case(7)
            call xc_partials_hse &
            ( hse_omega_pbe, rho, rho_gradient, en_density_xc, &
              en_density_x, en_density_c, local_xc_derivs, &
              xc_gradient_deriv )

          case(8)
!           Perdew-Wang 1991 LDA
            call xcpot_pw91_lda_shanghui(rho,en_density_xc,en_density_x,en_density_c, & 
                 local_xc_derivs, dVxc_drho)

          case(9)
!         BLYP
             call xc_partials_blyp &
            ( rho, rho_gradient, en_density_xc, &
              en_density_x, en_density_c, local_xc_derivs, &
               xc_gradient_deriv )
          case(10)
!         B3LYP - rpa parametrized
             call xc_partials_b3lyp &
            ( rho, rho_gradient, en_density_xc, &
              en_density_x, en_density_c, local_xc_derivs, &
              xc_gradient_deriv )
          case(11)
!         RPBE
           call xc_partials_rpbe &
               ( rho, rho_gradient, en_density_xc, &
              en_density_x, en_density_c, local_xc_derivs, &
                xc_gradient_deriv )
          case(12)
!         revPBE
           call xc_partials_revpbe &
               ( rho, rho_gradient, en_density_xc, &
                 en_density_x, en_density_c, local_xc_derivs, &
                xc_gradient_deriv )

          case(14)
!         revPBE_vdw.   SAG
           call xc_partials_revpbe_vdw &
                ( rho, rho_gradient, en_density_xc, &
                en_density_x, en_density_c, local_xc_derivs, &
                xc_gradient_deriv, coord_current )
           
         case(15)
!         VWN5
             call xc_vwn5_lda &
            ( rho, &
              local_xc_derivs, en_density_xc, &
              en_density_x, en_density_c)
          case(16)
!         VWN5
             call xc_vwn_lda &
            ( rho, &
              local_xc_derivs, en_density_xc, &
              en_density_x, en_density_c)
          case(17)
!         PBEsol
             call xc_partials_pbesol &
            ( rho, rho_gradient, en_density_xc, &
                 en_density_x, en_density_c, local_xc_derivs, &
              xc_gradient_deriv)
!@@edu>
          case(21) ! PBEint
             call xc_partials_pbeint &
            ( rho, rho_gradient, en_density_xc, &
                 en_density_x, en_density_c, local_xc_derivs, &
              xc_gradient_deriv)
!@@edu<
          case(13)
!         PBEsol0
             call xc_partials_pbesol0 &
            ( rho, rho_gradient, en_density_xc, &
                 en_density_x, en_density_c, local_xc_derivs, &
              xc_gradient_deriv)
          case(20)
!         AM05
!             if(spin_treatment.eq.1) then
!              write(use_unit,*) "Chosen functional does not work with ",&
!                          " spin-polarization (yet)."
!              stop
!             endif
             call xc_partials_am05 &
            ( rho, rho_gradient, en_density_xc, local_xc_derivs, &
              xc_gradient_deriv)
         case default
            write(use_unit,*) "Chosen type of XC is not yet implemented."
            stop
          end select

        end if

      end subroutine evaluate_xc_shanghui
!---------------------------------------------------------------------
!******

 

  subroutine xcpot_pz_lda_shanghui &
       ( density, en_density_xc,en_density_x, en_density_c, pot_xc, dVxc_drho)
    !  PURPOSE
    !    Calculates pz_lda xc energy and its derivatives with respect to the density.
    !  USES

    use constants
    use runtime_choices
    use dimensions
    implicit none

    !  ARGUMENTS

    real*8 :: density(n_spin)    
    real*8 :: en_density_xc,en_density_x, en_density_c
    real*8, dimension(n_spin) :: pot_xc
    real*8, dimension(n_spin) :: dVxc_drho

    !  INPUTS
    !   o density --  Full density rho(r)
    !  OUTPUT
    !   o pot_xc --  Exchange-correlation potential at current grid point
    !   o en_density_xc -- Exchange-correlation energy density at current grid point
    !   
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !  local variables:
    !       en_density_x : Exchange energy density at current grid point
    !       pot_x :        Exchange potential at current grid point



    !       for spin-polarized version

    real*8 :: rho
    real*8 :: zeta
    real*8 :: vxc_av
    real*8 :: dvxc

    real*8 :: EX,EC,VX,VC,DVXDN(n_spin),DVCDN(n_spin)
    integer:: i_spin

    !  begin work

    if (spin_treatment.eq.0) then
       ! unpolarized version - use routine written by shanghui

       call pz_lda(density(1),en_density_xc,en_density_x, en_density_c, pot_xc(1),dVxc_drho(1))
       

    else if (spin_treatment.eq.1) then
       ! Spin-polarized version.

       rho = density(1)+density(2)
       zeta = (density(1)-density(2)) / rho

       call stvxc_spin &
            (rho, zeta, vxc_av, dvxc, en_density_xc)

       pot_xc(1) = vxc_av + dvxc/2.d0
       pot_xc(2) = vxc_av - dvxc/2.d0

    end if

  end subroutine xcpot_pz_lda_shanghui


  subroutine xcpot_pw91_lda_shanghui &
       ( density, en_density_xc, en_density_x, en_density_c, &
         pot_xc, dVxc_drho  & 
       )

    !  PURPOSE
    !    Calculates pw_lda xc energy and its derivatives with respect to the density. for nspin=1 
    !    ref: Phys. Rev. B 45, 13244–13249 (1992) 
    !  USES

    use constants
    use dimensions
    use runtime_choices

    implicit none
    !  ARGUMENTS

    real*8, dimension(n_spin) ::  density
    real*8 en_density_xc
    real*8, dimension(n_spin) :: pot_xc
    real*8, dimension(n_spin) :: dVxc_drho

    !  INPUTS
    !   o density -- Full density rho(r)
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density at current grid point
    !   o pot_xc -- XC potential (possibly different for spin up, spin down 
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !  local variables
    !       aux_density : Extra variable to specifiy the actual input density which enters each subroutine
    !       en_density_c : Correlation energy density at current grid point
    !       pot_c :        correlation potential at current grid point (spin-up and spin-down)
    !       inv_k_fermi :  Inverse of the Fermi wave vector of a homogeneous electron gas of given density
    !       r_seitz:       Wigner-Seitz radius
    !       zeta :         would be spin polarisation
    !       ec_rs :        derivative of en_density_c w.r.t. r_seitz
    !       ec_zeta :      derivative of en_density_c w.r.t. zeta
    !       alpha_c :      "correlation contribution to spin stiffness"

    real*8 :: aux_density
    real*8, dimension(n_spin) :: en_density_x

    real*8 :: en_density_c
    real*8, dimension(n_spin) :: pot_c
    real*8 :: inv_k_f
    real*8 :: r_seitz
    real*8 :: zeta = 0.d0
    real*8 :: ec_rs = 0.d0
    real*8 :: ec_zeta = 0.d0
    real*8 :: alpha_c = 0.d0

    !       counters

    integer :: i_spin

    !  begin work


    !------------shanghui add here for an UNPOLARIZED pw-lda version for dVxc_drho----------
     if(spin_treatment.eq.0 ) then 
       call pw_lda(density(1),en_density_xc,en_density_x, en_density_c,pot_xc(1),dVxc_drho(1))
       return
     endif
    !-----------shanghui end add for an UNPOLARIZED pw-lda version ----------

    ! Exchange energy first, for both spin channels
    do i_spin = 1, n_spin, 1
       aux_density = dble(n_spin) * density(i_spin)
       call xlda( aux_density, pot_xc(i_spin), &
            en_density_x(i_spin) )
    enddo

    if (spin_treatment.eq.1) then

       aux_density = density(1)+density(2)
       zeta = ( density(1)-density(2) ) / aux_density

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*density(i_spin)
       enddo

       en_density_xc = en_density_xc / aux_density

    else

       en_density_xc = en_density_x(1)
       aux_density = density(1)
       zeta = 0.d0

    end if

    ! Correlation energy is next

    inv_k_f = (pisq3*aux_density)**third
    r_seitz = const_rs/inv_k_f

    if (spin_treatment.eq.0) then
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c(1), &
            pot_c(1), ec_rs, ec_zeta, alpha_c )
    else
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c(1), &
            pot_c(2), ec_rs, ec_zeta, alpha_c )
    end if

    en_density_xc = &
         en_density_xc + en_density_c

    do i_spin = 1, n_spin, 1
       pot_xc(i_spin) = pot_xc(i_spin) + pot_c(i_spin)
    enddo

  end subroutine xcpot_pw91_lda_shanghui
