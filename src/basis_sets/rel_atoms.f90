 !
 ! subroutine rel_atoms performs relativistic atom SCF by adopting the open
 ! source code dftatom (Computer Physics Communications 184, 1777-1791, 2013)
 ! to solve the radial Dirac equation.
 ! Currently for four-component Dirac-Kohn-Sham (4c-DKS) method, and exact 
 ! two-component (X2C) method.
 ! -- Rundong Zhao April, 2018
 !
 subroutine rel_atoms &
!   inputs     
     ( i_species, nuclear_charge, nspin, n_atomic, nqn, kappa, lnumber, occ_number,&
       n_grid, radial_grid, cut_atom, cut_type, r_cut, &
       cut_scale, cut_width, use_grad, i_exch,   &
!   outputs       
       atom_pot, atom_pot_es, atom_density, dens_1stderiv, dens_2ndderiv, density_diff, &
       upper_comp, upper_comp_deriv, lower_comp, lower_comp_deriv, eigenvalues )
       
    use constants,       only : light_speed, pi4_inv
    use dimensions,      only : n_max_grid, n_max_ind_fns
    use grids,           only : r_grid_min, r_grid_inc, log_r_grid_inc
    use localorb_io,     only : localorb_info, use_unit, OL_norm
    use mpi_tasks,       only : aims_stop_coll
    use reigen_dftatom
    implicit none
!   Inputs
    real*8  :: nuclear_charge                        ! charge of the species
    integer :: i_species, nspin, n_atomic            ! nspin: spin channels. n_atomic: number of atomic basis functions
    integer, dimension(n_atomic) :: nqn              ! primary quantum number
    integer, dimension(n_atomic) :: kappa            ! quantum number kappa
    integer, dimension(n_atomic) :: lnumber          ! l quantum number
    real*8 , dimension(n_atomic) :: occ_number       ! occupation number for orbitals
    integer :: n_grid                                ! number of grid points in radial grid
    real*8, dimension(n_grid) :: radial_grid         ! radial grid
    integer :: i_exch                                ! type of exchange-correlation functional to be used
    logical :: cut_atom                              ! do we use cutoff potential (confinement potential) or not
    integer :: cut_type
    real*8  :: r_cut
    real*8  :: cut_scale
    real*8  :: cut_width
    integer :: use_grad
!  Outputs    
    real*8, dimension(n_max_grid,nspin) :: atom_pot                   ! effective atomic radial potential (hartree potential + XC)
    real*8, dimension(n_max_grid,nspin) :: atom_pot_es                ! electrostatic potential
    real*8, dimension(n_max_grid,nspin) :: atom_density               ! electron density
    real*8, dimension(n_max_grid,nspin) :: density_diff               ! electron density difference: rho_4c - rho_large
    real*8, dimension(n_max_grid,nspin) :: dens_1stderiv              ! 1st derivative of density
    real*8, dimension(n_max_grid,nspin) :: dens_2ndderiv              ! 2nd derivative of density
    real*8, dimension(n_max_grid,n_max_ind_fns) :: upper_comp         ! upper(large) component of dirac spinor
    real*8, dimension(n_max_grid,n_max_ind_fns) :: upper_comp_deriv   ! derivative of the upper component
    real*8, dimension(n_max_grid,n_max_ind_fns) :: lower_comp         ! lower(small) component of dirac spinor
    real*8, dimension(n_max_grid,n_max_ind_fns) :: lower_comp_deriv   ! derivative of the lower component
    real*8, dimension(n_atomic)  :: eigenvalues

!   external cutoff function for potential
    real*8, external :: cutoff_pot
!    
!   local variables     
    real*8, allocatable :: density(:)       ! electron density
    real*8, allocatable :: density_large(:) ! large component electron density
    real*8, allocatable :: density_1st(:)   ! 1st derivative of electron density
    real*8, allocatable :: density_2nd(:)   ! 2nd derivative of electron density
    real*8, allocatable :: potential(:)     ! total potential (Hartree potential + XC potential )
    real*8, allocatable :: prev_e(:)        ! eigenvalues of the last iteration
    real*8, allocatable :: nuc_pot(:)       ! -Z/r
    real*8, allocatable :: hartree(:)       ! Harthree potential
    real*8, allocatable :: xc_pot(:)        ! exchange-correlation potential
    real*8, allocatable :: v1(:),v2(:),v3(:)! temporary array
    real*8, allocatable :: xc_energy(:)
    real*8, allocatable :: x_energy(:)
    real*8, allocatable :: c_energy(:)
    real*8, allocatable :: r_prime(:)       ! grid derivatives 
    real*8, allocatable :: weight(:)        ! radial integration weight
    real*8, allocatable :: spline_coef(:,:) ! spline_coef for density_1st and density_2nd

    integer :: i,j,k,l,m,n
    integer :: iter, i_outer, relat, converged_dftatom

    real*8 :: sum_elec, temp1, conv_test
    real*8 :: r_outer
    real*8 :: mix_param = 0.5d0  ! Parameter for density mixer

    integer, parameter :: maxscf = 100
    real*8,  parameter :: eps = 1.d-6
    integer, parameter :: order = 5
    character*200 :: info_str

    allocate( density(n_grid), density_large(n_grid), density_1st(n_grid), density_2nd(n_grid), &
              potential(n_grid), prev_e(n_atomic), nuc_pot(n_grid), &
              hartree(n_grid), xc_pot(n_grid), v1(n_grid), v2(n_grid), v3(n_grid), &
              xc_energy(n_grid), x_energy(n_grid), c_energy(n_grid), &
              r_prime(n_grid), weight(n_grid), spline_coef(4,n_grid) )
 
   ! nuclear potential:
    do i=1,n_grid
      nuc_pot(i) = -nuclear_charge/radial_grid(i)
    end do

   ! Initialize the potential with nuclear potential 
    potential = nuc_pot
    eigenvalues=0.d0
    prev_e=0.d0
    v1=0.d0; v2=0.d0

  ! Rundong: dftatom requires that we provide the derivative at each grid. In aims, the radial grids are practically generated through:
  ! r_i = r_min * e^{alpha*(i-1)}, in which the e^{alpha} is saved in r_grid_inc, and alpha is in log_r_grid_inc. Thus,the grid derivative
  ! r_i^{prime} = r_min / e^{alpha} * alpha * e^{alpha*i} = alpha * r_i
    do i=1, n_grid
      r_prime(i) = log_r_grid_inc(i_species) * radial_grid(i)
      ! For logarithmic grids used by aims, the radial integration weight is:
      weight(i) = log_r_grid_inc(i_species) * radial_grid(i)
    enddo

    write(info_str,'(2X,A)') "Use dftatom to solve the Dirac equation of the free atoms"
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,l3)') "cut_atom: ", cut_atom
    call localorb_info(info_str,use_unit,'(A)',OL_norm)


  !------------------ SCF loop ------------------
    do iter=1,maxscf

       ! Add cutoff potential to the radial potential
        if (cut_atom) then
       ! Set outermost radius for radial function integration
       !  r_outer = r_cut + cut_width
          do j=1,n_grid
            potential(j) = potential(j) + cutoff_pot(radial_grid(j), cut_type, r_cut, cut_width, cut_scale )
          end do
       !  i_outer = 1 + ( log(r_outer/radial_grid(1)))/(log(radial_grid(2)/radial_grid(1)))
       !else
       !  i_outer = n_grid
        end if
 
      ! Solve radial dirac equation for each occipied orbital
        do j=1,n_atomic

          if( kappa(j) .eq. -lnumber(j)-1 )then
            relat=2  ! spin up
          elseif( kappa(j) .eq. lnumber(j) )then
            relat=3  ! spin down
          else
            call aims_stop_coll('illegal kappa','rel_atom')
          endif

          call solve_radial_eigenproblem(n_grid,nqn(j), lnumber(j), -10.d0, 1.d-10, 100, radial_grid, r_prime, potential, nuclear_charge,&
               light_speed, relat, .true., -1.d4, 10.d0, converged_dftatom, eigenvalues(j), upper_comp(1,j), lower_comp(1,j))
         !write(use_unit,"('iter:',i4,5x,'orb:',i3,5x,'n:',i3,5x,'l:',i3,5x,'k:',i3,5x,'Z:',f8.4,5x,'converge:',i3,5x,'eigen:',e12.6)")iter,j,nqn(j),lnumber(j),kappa(j),nuclear_charge,converged_dftatom,eigenvalues(j)
        enddo


      ! Calculate electron density. Note, the theoretical expression of rho should be:
      ! rho(i) = 1/4pi * \sum_{orbitals} {occ_nuber * (P**2+Q**2)/r**2} (see the cpc paper in 2013).
      ! However, the rho generated here does not involve 1/4pi, since in both subroutine vestat and vexcor,
      ! they request a rho without 1/4pi.
        density(:) = 0.0d0
        do j=1,n_atomic
        do k=1,n_grid
          density(k) = density(k) + (occ_number(j)*(upper_comp(k,j)**2+lower_comp(k,j)**2)) / (radial_grid(k)**2)
        enddo
        enddo

        sum_elec=0.d0 ! apporximately equals to nuclear_charge
        do k=1, n_grid
          sum_elec = sum_elec + weight(k) * density(k) * radial_grid(k)**2
        enddo

       ! DO NOT normalize the density here, like writing:
       ! density = nuclear_charge/sum_elec * density
       ! For neutral atoms, this can be done (actually not necessary, as it brings in very little difference for atomic systems).
       ! However, we also use rel_atoms to solve ions, such a normalization can lead to errors.

       ! Compute density gradiants if needed
        if(use_grad.eq.1) then
          call fderiv(1,n_grid,radial_grid,density,density_1st,spline_coef)
          call fderiv(2,n_grid,radial_grid,density,density_2nd,spline_coef)
        end if 


       ! Calculate Hartree potential (including the nuclear potential and the Coulomb potential from other electrons)
        call vestat(n_grid, sum_elec-nuclear_charge, temp1, radial_grid, density, hartree, .false.)


       ! Exchage-correlation potential
        call vexcor(i_exch,n_grid,radial_grid,density,density_1st,density_2nd,xc_pot,xc_energy,x_energy,c_energy,.false.)
 

       ! New electronic potential for the next iteration
        v3 =  hartree + xc_pot

        call anderson(mix_param,iter,n_grid,radial_grid,v1,v3,v2,potential)

       ! Calculate mean difference between new and old densities
        temp1=0.d0
        do k=1, n_atomic
          temp1 = temp1 + (eigenvalues(k)-prev_e(k))**2
        end do
        conv_test=sqrt(temp1/dble(n_atomic))

        prev_e = eigenvalues
        
        write(info_str,'(2X,A,E14.6)') "Convergence of radial Dirac equation:", conv_test
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
        if ((iter.gt.2).and.(conv_test.lt.eps)) exit   

    end do 
  !------------------ end SCF loop ------------------

    if (iter .eq. maxscf) then
      do i=1, 10
        write(use_unit,*) "WARNING!: Atom SCF is not converged "
      enddo
    end if 

!
! Calcualte derivatives of wave functions is required
!
    if(use_grad.eq.1) then
      do j= 1, n_atomic
        call fderiv(1,n_grid,radial_grid,upper_comp(:,j),upper_comp_deriv(:,j),spline_coef)
        call fderiv(1,n_grid,radial_grid,lower_comp(:,j),lower_comp_deriv(:,j),spline_coef)
      end do
    end if 

    density_large(:) = 0.0d0
    do j=1,n_atomic
    do k=1,n_grid
      density_large(k) = density_large(k) + occ_number(j)*upper_comp(k,j)**2 / (radial_grid(k)**2)
    enddo
    enddo

    atom_pot(1:n_grid,1) = potential(1:n_grid)
    atom_pot_es(1:n_grid,1) = hartree(1:n_grid)
    atom_density(1:n_grid,1) = density(1:n_grid)
    dens_1stderiv(1:n_grid,1) = density_1st(1:n_grid)
    dens_2ndderiv(1:n_grid,1) = density_2nd(1:n_grid)
    density_diff(1:n_grid,1) = density(1:n_grid)-density_large(1:n_grid)


! Print the wavefunctions    
    if(.false.) then
      write(use_unit,*)'--------------------------------------------'
      write(use_unit,*)'     Output of the Radial Dirac Equation'
      write(use_unit,*)'--------------------------------------------'
      do j = 1,n_atomic
        write(use_unit,"('n, l, kappa, E:',i4,i4,i4,f20.13)")nqn(j),lnumber(j),kappa(j),eigenvalues(j)
      enddo 
    endif

    deallocate( density, density_large, density_1st, density_2nd, potential, prev_e, nuc_pot, hartree, xc_pot, &
              v1, v2, v3, xc_energy, x_energy, c_energy, r_prime, weight, spline_coef )

 end subroutine rel_atoms



