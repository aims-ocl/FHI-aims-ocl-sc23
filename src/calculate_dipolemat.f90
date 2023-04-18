!****h* FHI-aims/calculate_dipolemat
!  NAME
!    calculate_dipolemat 
!  SYNOPSIS

module calculate_dipolemat 
!  PURPOSE
!  This module contains all routines related to the evaluation of the Momentum-
!  matarix-element as postprocessing
!
!  AUTHOR
!    Bjoern Bieniek
!  HISTORY
!    Development version, FHI-aims (2010).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  SOURCE
  implicit none

save
! the momentum matrix as function of space p = <psi |grad  psi>

real*8, dimension(:), allocatable :: dipole_mat_full_oned
real*8, dimension(:), allocatable :: dipole_mat_full_oned_two
real*8, dimension(:), allocatable :: die_el
real*8, dimension(:), allocatable :: seebeck
real*8, dimension(:), allocatable :: abtew
real*8, dimension(:), allocatable :: abtewseebeck
real*8, dimension(:), allocatable :: fermideriv

real*8, dimension(:,:), allocatable :: dipolemat_full_w
complex*16, dimension(:,:), allocatable :: dipolemat_full_w_complex
complex*16, dimension(:), allocatable :: dipelement_one
complex*16, dimension(:), allocatable :: dipelement_two
complex*16, dimension(:), allocatable :: dipelement_three
real*8:: omegapl, abtewintegral_cond, abtewintegral_seeb
integer :: n_state_min,n_state_max

contains
subroutine allocate_spectra
!  PURPOSE
!    allocation of arrays containing spectral values of dielectric function, optical conductivity, seebeck coeff ...
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: check_allocation

  integer :: info
    if (.not.allocated( die_el)) then
      allocate( die_el ( n_omega),stat=info)
      call check_allocation(info, 'die_el                      ')
    end if


  if (.not.allocated( seebeck)) then
     allocate( seebeck ( n_omega),stat=info)
     call check_allocation(info, 'seebeck                      ')
  end if

  if (.not.allocated( abtew)) then
     allocate( abtew (n_greenenergy),stat=info)
     call check_allocation(info, 'abtew                      ')
  end if

  if (.not.allocated( abtewseebeck)) then
     allocate( abtewseebeck (n_greenenergy),stat=info)
     call check_allocation(info, 'abtewseebeck                      ')
  end if

  if (.not.allocated( fermideriv)) then
     allocate( fermideriv (n_greenenergy),stat=info)
     call check_allocation(info, 'fermideriv                      ')
  end if
end subroutine allocate_spectra


subroutine allocate_dipole_mat
!  PURPOSE
!    allocation of momentum matrix
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: check_allocation
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
  integer :: info

!  if(out_mommat) then  !!! this is would be useful to separate spectrum and momentum-matrix tasks
!    if (.not.allocated(dipole_mat_full_oned)) then
!      allocate(dipole_mat_full_oned(n_cells_in_hamiltonian*n_basis*n_basis),stat=info)
!      call check_allocation(info, 'dipole_mat_full_oned                      ')
!    end if
!  else
    if (.not.allocated(dipole_mat_full_oned)) then
      allocate(dipole_mat_full_oned(n_cells_in_hamiltonian*n_basis*n_basis),stat=info)
      call check_allocation(info, 'dipole_mat_full_oned                      ')
    end if

end subroutine allocate_dipole_mat

subroutine allocate_dipole_mat_two
!  PURPOSE
!    allocation of momentum matrix
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: check_allocation
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
  integer :: info
  if (.not.allocated(dipole_mat_full_oned)) then
     allocate(dipole_mat_full_oned(n_cells_in_hamiltonian*n_basis*n_basis),stat=info)
     call check_allocation(info, 'dipole_mat_full_oned                      ')
  end if

  if (.not.allocated(dipole_mat_full_oned_two)) then
     allocate(dipole_mat_full_oned_two(n_cells_in_hamiltonian*n_basis*n_basis),stat=info)
     call check_allocation(info, 'dipole_mat_full_oned_two                      ')
  end if

end subroutine allocate_dipole_mat_two

subroutine allocate_dipole_mat_k
!  PURPOSE
!    allocation of dipole matrix, integrated, Fourier transformed, basis transformed
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: check_allocation
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
  integer :: info
  if (real_eigenvectors) then
      if (.not.allocated( dipolemat_full_w)) then
	allocate( dipolemat_full_w ( n_basis,n_basis),stat=info)
	call check_allocation(info, 'dipolemat_full_w                      ')
      end if
  else
      if (.not.allocated( dipolemat_full_w_complex)) then
	allocate( dipolemat_full_w_complex ( n_basis,n_basis),stat=info)
	call check_allocation(info, 'dipolemat_full_w_complex                      ')
      end if
  endif
  if (ep1==ep2) then
    if (.not.allocated( dipelement_one)) then
      allocate( dipelement_one ((n_states+1)*n_states/2),stat=info)
      call check_allocation(info, 'dipelement_one                      ')
    end if
  else
    if (.not.allocated( dipelement_one)) then
      allocate( dipelement_one ((n_states+1)*n_states/2),stat=info)
      call check_allocation(info, 'dipelement_one                      ')
    end if
    if (.not.allocated( dipelement_two)) then
      allocate( dipelement_two ((n_states+1)*n_states/2),stat=info)
      call check_allocation(info, 'dipelement_two                      ')
    end if
  endif
  if(out_mommat) then
      if (.not.allocated( dipelement_one)) then
	allocate( dipelement_one ((n_states+1)*n_states/2),stat=info)
	call check_allocation(info, 'dipelement_one                      ')
      end if
      if (.not.allocated( dipelement_two)) then
	allocate( dipelement_two ((n_states+1)*n_states/2),stat=info)
	call check_allocation(info, 'dipelement_two                      ')
      end if
      if (.not.allocated( dipelement_three)) then
	allocate( dipelement_three ((n_states+1)*n_states/2),stat=info)
	call check_allocation(info, 'dipelement_three                     ')
      end if
  endif
end subroutine allocate_dipole_mat_k


subroutine clean_dipole_mat
!  NAME
!    clean_dipole_mat - deallocation of dipole matrix
!  SYNOPSIS
!  PURPOSE
!    deallocation of local variables
!  USES
   use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  if (allocated(dipelement_one))  deallocate(dipelement_one)
  if (allocated(dipelement_two))  deallocate(dipelement_two)
  if (allocated(dipelement_three))  deallocate(dipelement_three)
  if (allocated(dipolemat_full_w_complex))  deallocate(dipolemat_full_w_complex)
  if (allocated(dipolemat_full_w))  deallocate(dipolemat_full_w)
end subroutine clean_dipole_mat

subroutine clean_dipole_mat_final
!  PURPOSE
!    deallocation of local variables
!  USES
   use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  if (allocated(dipole_mat_full_oned))  deallocate(dipole_mat_full_oned)
  if (allocated(dipole_mat_full_oned_two))  deallocate(dipole_mat_full_oned_two)
  if (allocated(die_el))  deallocate(die_el)
  if (allocated(seebeck))  deallocate(seebeck)
  if (allocated(abtew))  deallocate(abtew)
  if (allocated(abtewseebeck))  deallocate(abtewseebeck)
  if (allocated(fermideriv))  deallocate(fermideriv)


end subroutine clean_dipole_mat_final


subroutine evaluate_dipole_mat( & 
     n_points, partition, n_compute, gradient_basis_wave, &
     n_basis_list, wave, dipolemat_shell,i_coord)

!  PURPOSE
!  Evaluates mom_mat for one grid batch,
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: n_basis_list
  integer :: n_points ! number of points in this batch
  integer :: n_compute ! number of non zero basis functions in this batch
  real*8  :: partition(n_points)
  real*8  :: wave(n_basis_list, n_points) ! input: basis wave
  real*8  :: gradient_basis_wave(n_basis_list,3, n_points) ! input: the gradient of the basis (x,y or z)
  real*8  :: dipolemat_shell( n_compute, n_compute)

! INPUTS
! o gradient_basis_wave
! o wave
!  OUTPUT
! o dipolemat_shell - see descriptions above
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


  !  begin work
      real*8 wave_compute(n_compute,n_points)
      real*8 grad_wave_compute(n_compute,n_points) 
!     counters

      integer :: i_point, i_coord
!     begin work

      dipolemat_shell = 0.0
!     Condense basis functions to only those that are used at present integration points

      do i_point = 1, n_points, 1
         wave_compute(1:n_compute, i_point) =  &
              sqrt(partition(i_point))* &
                wave(1:n_compute, i_point)
         grad_wave_compute(1:n_compute, i_point) =  &
              sqrt(partition(i_point))*gradient_basis_wave(1:n_compute,i_coord, i_point)   
      enddo

        call dgemm('N', 'T', n_compute, n_compute,  &
            n_points, 1.0d0,  &
             wave_compute, n_compute, grad_wave_compute,  &
             n_compute, 0.0d0, dipolemat_shell, n_compute )
end subroutine evaluate_dipole_mat


function ergebnis(abtew, abtewseebeck, zellgr, switch)
!  use pbc_lists
!  use geometry
!  use dimensions
  use runtime_choices
!  use basis
  implicit none

  real*8:: abtew(n_greenenergy)
  real*8:: abtewseebeck(n_greenenergy)
  real*8:: summe_cond, summe_seeb
  real*8:: Emin_ha, Emax_ha, ergebnis, prefactorAbtewCond, prefactorAbtewSeeb, zellgr, energie
  integer :: i, maxcounter, switch


    Emin_ha=Emin/hartree
    Emax_ha=Emax/hartree
    summe_cond=0.0d0
    summe_seeb=0.0d0   
    maxcounter=n_greenenergy/3 ! for more detailed simpsons rule; modulo abfrage in read_control. 

select case (switch)
  case(1) ! Conducivity Integral

!  prefactorAbtew=(2*pi*((-1.602176565E-19)**2) * (1.519829877E-16/(2*pi))&
! * 9.109382E-31)/(100 * 1 * (zellgr * bohr**3)*1E-30 * (9.109382E-31)**2  ) 
! (e**2 * hbar * m) / (V * m**2) is in units (Coulumb**2 (Ha*s) * kg )/(meter**3 * kg**2 ) 
!upper kg has to be added to due to absence of electronic mass in the 
! momentum matrix-elements, hbar-Hatree is compatible with 1/Ha unit of the summation
    prefactorAbtewCond = 42827.9032811 * (1/(zellgr * bohr**3))   ! is the conductivity prefactor. 


   do i=1, maxcounter
     summe_cond= summe_cond + abtew(3*i) + 3*abtew(3*i-1) + 3*abtew(3*i-2) + abtew(3*i-3) ! gives 1/Ha unit
   enddo
   ergebnis=prefactorAbtewCond * 0.375d0*((Emax_ha-Emin_ha)/n_greenenergy) * summe_cond ! Simpson Rule. has Hartree dimension!

  case(2) ! Seebeck Integral

!  prefactorAbtew=(2*pi*((-1.602176565E-19)**2) * (1.519829877E-16/(2*pi)) &
!* 9.109382E-31)  /  (100 * 1 * (zellgr * bohr**3)*1E-30 * (9.109382E-31)**2  ) 
! (e**2 * hbar * m) / (V * m**2) is (Coulumb**2 (Ha*s) * kg )/(meter**3 * kg**2 ) upper kg has 
!to be added to due to absence of electronic mass in the momentum matrix-elements
! prefactorAbtew-seebeck= (2pi * (-1.602e-19)(A*s) * (6.62606957E-34/2pi)(V*A*s*s)  * 
!  * (9.109382E-31)kg ) / (1 * (zellgr*bohr**3)*1E-30 * (9.109382E-31)**2(kg**2)   ) 
    prefactorAbtewSeeb = -1.165406543E8 * (1/(zellgr * bohr**3))

 
   do i=1, maxcounter
     summe_seeb= summe_seeb + abtewseebeck(3*i) + 3*abtewseebeck(3*i-1) + 3*abtewseebeck(3*i-2) + abtewseebeck(3*i-3) 
   enddo
                     
    ergebnis=prefactorAbtewSeeb * 0.375d0*((Emax_ha-Emin_ha)/n_greenenergy) * summe_seeb
endselect
end function

subroutine calc_fermideriv(fermideriv, chemical_potential)

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none
  integer :: i
  real*8:: Emin_ha, Emax_ha, energie, chemical_potential, widthf
  real*8:: fermideriv(n_greenenergy)

   widthf=widthtwo/hartree
   Emin_ha=Emin/hartree
   Emax_ha=Emax/hartree

  do i=1, n_greenenergy

   energie=Emin_ha+i*((Emax_ha - Emin_ha)/n_greenenergy)
   fermideriv(i)=(exp((energie+chemical_potential)/widthf))/ &
                 (widthf*((exp(chemical_potential/widthf)+exp(energie/widthf))**2))

  enddo 

end subroutine 


subroutine calc_diel(dipelementxi,dipelementxj,die_el,seebeck ,abtew, abtewseebeck, fermideriv, omegapl,&
 chemical_potential,KS_eigen,KS_eigen_full,k_weight)

  !  PURPOSE
  !   Sum up Momentummatrix elements to dielectric function (Fermis Golden rule)
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  real*8:: omega,omegaev,fermione,fermitwo,dfermi,energie
  real*8:: omegapl
  real*8:: k_weight
  real*8:: die_el( n_omega)
  real*8:: seebeck( n_omega)
  real*8:: abtew (n_greenenergy)
  real*8:: abtewseebeck (n_greenenergy)
  real*8:: fermideriv (n_greenenergy)
  real*8:: Emin_ha, Emax_ha
  real*8:: widthtolerance
  complex*16:: dipelementxi ((n_states+1)*n_states/2)
  complex*16:: dipelementxj ((n_states+1)*n_states/2)
  real*8:: dipmultfeld ((n_states+1)*n_states/2)
  real*8:: scaling, scaling_one, scaling_two
  real*8 , dimension(n_states, n_spin), INTENT(IN) ::     KS_eigen
  real*8 , dimension(n_states, n_spin,n_k_points), INTENT(IN) ::     KS_eigen_full
  real*8 :: weight
  real*8:: chemical_potential,width,ohm,widthf, dipmult
  real*8:: fermioccupations(n_states)
  
  integer:: n_state,m_state,num,omegaind, energieind


     width=widthone/hartree
     widthf=widthtwo/hartree

    Emin_ha=Emin/hartree
    Emax_ha=Emax/hartree

num = 0
do n_state = n_state_min, n_state_max
	 do m_state = n_state, n_state_max
     num = num + 1
       dipmultfeld(num)=abs(dipelementxi(num)*dipelementxj(num))
     if (n_state==m_state) then
         omegapl=omegapl+sqrt(1.0/(2.0*pi))*(1.0/(width))*exp(-0.5*((KS_eigen(n_state,1)&
         -chemical_potential)**2/((width)**2)))*k_weight*(abs(dipelementxi(num)*dipelementxj(num)))/pi
     endif
     enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! d.c. value calculation (Abtew et al. 2007) 

if (flag_out_dclimit) then


! region-restriction wrt Fermi-deriv can be set up already here. 

   do energieind=0, n_greenenergy-1
      energie=Emin_ha+energieind*((Emax_ha - Emin_ha)/n_greenenergy)
      num=0
    if(fermideriv(energieind+1)<1E-30) then
      abtew(energieind+1)=0.0d0
    else
      do n_state = n_state_min, n_state_max
			scaling_one=sqrt(1.0/(2.0*pi))*(1.0/width)*exp(-0.5*((( KS_eigen(n_state,1)-energie  )**2)/((width)**2)))  
         do m_state = n_state, n_state_max
            num = num + 1
		  if (n_state/=m_state) then
       		scaling_two=sqrt(1.0/(2.0*pi))*(1.0/width)*exp(-0.5*((( KS_eigen(m_state,1)-energie  )**2)/((width)**2)))
            abtew(energieind+1)=abtew(energieind+1) + k_weight*dipmultfeld(num)* &
                                scaling_one*scaling_two*fermideriv(energieind+1)
            abtewseebeck(energieind+1)=abtewseebeck(energieind+1) + &
                k_weight*dipmultfeld(num)*scaling_one*scaling_two*fermideriv(energieind+1) * &
                ((KS_eigen(n_state,1)+KS_eigen(m_state,1))/2 - chemical_potential)
		  endif ! n!=m condition
	    enddo ! end m-state
	  enddo ! end n-state
    endif ! end of fermideriv>1e-30 condition
    enddo ! end of energy loop
endif ! END of Abtew case 


do n_state= n_state_min, n_state_max
fermioccupations(n_state)=(1.0/(1.0+exp((KS_eigen(n_state,1)-chemical_potential)/widthf)))
enddo


    do omegaind = 0, n_omega-1
        omegaev=omega_min+omegaind*((omega_max-omega_min)/n_omega)
        num=0
        omega = omegaev/hartree
	  do n_state = n_state_min, n_state_max
	    do m_state = n_state, n_state_max
	      num = num + 1
           if(n_state/=m_state) then 
              ohm=KS_eigen(m_state,1)-KS_eigen(n_state,1)
              scaling=sqrt(1.0/(2.0*pi))*(1.0/width)*exp(-0.5*(((ohm-omega)**2)/((width)**2)))
	          dfermi=(fermioccupations(n_state)-fermioccupations(m_state))              
              die_el(omegaind+1) = die_el(omegaind+1)+k_weight*dipmultfeld(num)*dfermi*scaling
              seebeck(omegaind+1)=seebeck(omegaind+1)+k_weight*dipmultfeld(num)*dfermi*scaling*&
              ((KS_eigen(n_state,1)+KS_eigen(m_state,1))/2-chemical_potential)
          endif
	    enddo
	  enddo
    enddo


end subroutine calc_diel


subroutine out_die_el(die_el,seebeck,abtew, abtewseebeck, fermideriv, ep1_in,ep2_in, zellgr, chempot)
  !  PURPOSE
  !   Write Dielectric function \epsilon_ep1_ep2 to file
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  !  ARGUMENTS
  real*8, INTENT(IN) :: die_el(n_omega)
  real*8, INTENT(IN) :: seebeck(n_omega)
  real*8, INTENT(IN) :: abtew(n_greenenergy)
  real*8, INTENT(IN) :: abtewseebeck(n_greenenergy)
  real*8, INTENT(IN) :: fermideriv(n_greenenergy)

  CHARACTER(len=1), INTENT(IN):: ep1_in
  CHARACTER(len=1), INTENT(IN):: ep2_in
  real*8:: zellgr
  real*8:: chempot

  ! Write dipelement to file
  !  INPUTS
  !    o die_el -- n_omega Dielectric function
  !  OUTPUT
  !    o file
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: num, ierror
  integer:: omega
  real*8:: omega_value, energie_value, sigma, abtewvalue
  real*8:: abtewseebeckvalue
  real*8:: prefactorL11
  real*8:: prefactorL12

  CHARACTER(len=55) :: fmt1
  CHARACTER(len=55) :: fmt2
  CHARACTER(len=55) :: fmt3
  CHARACTER(len=45) :: nameOmega
  CHARACTER(len=55) :: nameDClimit
  CHARACTER(len=10) :: fermi
  ! file indices can be defined here already

  write(fermi, '(f7.4)')(chempot*hartree)

  nameOmega = "dielectric_full_matrix_"//trim(ep1_in)//"_"//trim(ep2_in)//"_"//trim(fermi)//".out"
  nameDClimit= "DC_limit_transport_full_matrix_"//trim(ep1_in)//"_"//trim(ep2_in)//"_"//trim(fermi)//".out"

! equals (2pi * e*e * hbar )/(1 * CellVol * mass * 100(cm) ) - 
! EXPLICIT prefactorL11=(2*pi * ((-1.602176565E-19)**2) * (4.135667E-15/(2*pi)) * 9.109382E-31) /
! (100 * 1*(zellgr * bohr**3)*1E-30*(9.109382E-31)**2) ! is (Coulumb**2 (Ha*s) * kg )/(meter**3 * kg**2 ), 100 is for m->cm
 prefactorL11= 1.16540652466953E6 * (1/(1*(zellgr * bohr**3))) 
! factor "1" before cellsize to be changed to "3" when sum runs over the three spatial directions, xx, yy, zz. 
! EXPLICIT prefactorL12= (2*pi * -1.602176565E-19 * (4.135667E-15/(2*pi)) 
!(hbar in eV*s - as forthcoming omega-division is in eV as well) * &
! 4.35974434E-18 (Ha-in-Joule) * 9.109382E-31 )/(1 * (zellgr * bohr**3) * 1E-30*(9.109382E-31)**2))
 prefactorL12 = -3.17123256619787E9 *1E6 * 0.01 * (1/(1*(zellgr * bohr**3)))  
! equals (2pi * e * hbar * Hartree-to-Joule ) / (1 * CellVol * mass  ), factor 1e6 to give muV, 
!factor 0.01 to correct the Ohm*cm unit of sigma that enters S (not L12) , negative sign is for single elementary charge


if(flag_out_dclimit) then

  fmt1 = '(5ES17.4)'
  open(unit=8, file=nameOmega,ACTION='WRITE')
  WRITE (8,122) trim(ep1_in), trim(ep2_in)
  122 FORMAT ('# omega (eV), Im(\epsilon_{',A1,'_',A1, &
     '})   L11-Conductivity (Ohm-1*cm-1)       L12      (L12/(T*L11))-', &
     'Seebeck(muV/K)')
  num=0
     do omega=0,n_omega-1
         num=num+1
         omega_value=omega_min+omega*((omega_max-omega_min)/n_omega) ! omega in eV
         energie_value=Emin+omega*((Emax-Emin)/n_omega)
         sigma=prefactorL11*die_el(num)/omega_value
             if ( sigma < 1.0E-30)then
               sigma=0.0d0
             endif
         WRITE (8,fmt1) omega_value, die_el(num)/(pi*((omega_value/hartree)*(omega_value/hartree ))), sigma,& 
                (prefactorL12*seebeck(num)/omega_value) ,(prefactorL12*seebeck(num)/omega_value)/(11604.519*widthtwo* &
                 prefactorL11*die_el(num)/omega_value)
!        WRITE (8,fmt1) omega_value, die_el(num)/(pi*(omega_value**2)), sigma, (prefactorL12*seebeck(num)/omega_value), 
!(prefactorL12*seebeck(num)/omega_value)/(11604.519*widthtwo* prefactorL11*die_el(num)/omega_value) 
! differs in hartree or eV for dielectric function

     end do
  close(unit=8)

  fmt2 = '(4ES17.4)'
  open(unit=9, file=nameDClimit,ACTION='WRITE')
  WRITE (9,125)
  125 FORMAT ('# Energy (eV)     AbtewSigmaValue   AbtewSeebeckValue   Fermideriv')
  num=0
     do omega=0,n_greenenergy-1
         num=num+1
         energie_value=Emin+omega*((Emax-Emin)/n_greenenergy)
         abtewvalue=abtew(num)
         abtewseebeckvalue=abtewseebeck(num)
             if (abtewvalue<1.0E-30)then
               abtewvalue=0.0d0
             endif   
             if (abs(abtewseebeckvalue) < 1.0E-30)then
               abtewseebeckvalue=0.0d0
             endif   
         WRITE (9,fmt2) energie_value, abtewvalue, abtewseebeckvalue, fermideriv(num)
     end do
  close(unit=9)


else

  fmt3 = '(5ES17.4)'
  open(unit=8, file=nameOmega,ACTION='WRITE')
  WRITE (8,124) trim(ep1_in), trim(ep2_in)
  124 FORMAT ('# omega (eV), Im(\epsilon_{',A1,'_',A1, &
     '})   L11-Conductivity (Ohm-1*cm-1)       L12      (L12/(T*L11))-', &
     'Seebeck(muV/K)')
  num=0
     do omega=0,n_omega-1
         num=num+1
         omega_value=omega_min+omega*((omega_max-omega_min)/n_omega)
         sigma=prefactorL11*die_el(num)/omega_value
          if ( sigma<1.0E-30)then
            sigma=0.0d0
          endif
         WRITE (8,fmt3) omega_value, die_el(num)/(pi*((omega_value/hartree)*(omega_value/hartree ))), sigma, &
                (prefactorL12*seebeck(num)/omega_value) ,&
                 (prefactorL12*seebeck(num)/omega_value)/(11604.519*widthtwo* prefactorL11*die_el(num)/omega_value)
     end do
  close(unit=8)
endif


end subroutine out_die_el


subroutine construct_dipolemat_p1(dipole_mat_full_oned, dipolemat_full_w, &
     dipolemat_full_w_complex,k_phase_ac, work_ham )


  !  PURPOSE
  !   Contruct Momentum matrix for the lapack eigenvalue solver.
  !   This is needed in periodic systems and if the packed matrix format is in use.
  !
  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop
  implicit none

  !  ARGUMENTS

  real*8::     dipole_mat_full_oned (n_cells_in_hamiltonian*n_basis*n_basis)
  real*8::     dipolemat_full_w         ( n_basis,n_basis)
  complex*16:: dipolemat_full_w_complex ( n_basis,n_basis)
  real*8, dimension( n_centers_basis_I , n_centers_basis_I ),optional :: work_ham
  complex*16:: k_phase_ac(n_cells_in_hamiltonian)


  ! JW: In principle, optional dummy arguments are only allowed with an
  ! explicit interface (e.g. within a module subroutine).  While I do not
  ! really expect trouble, this might be a possible explanation if there is
  ! trouble with some compilers.

  !  INPUTS
  !    o dipole_mat -- Momentum matrix with centers terms separately
  !    o k_point -- k-point wanted to be calculated 
  !    o work_ham -- work space needed if for non-packed format is in use
  !  OUTPUT
  !    o dipolemat_w -- Momentum matrix ready to go lapack eigenvalue solver if real eigenvector are used
  !    o dipolemat_w_complex -- Momentum matrix ready to go lapack eigenvalue solver  if real eigenvector are used
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !
  ! SOURCE

  complex*16, dimension(:,:),allocatable :: work

  integer:: n_first_part, i_cell,      i_basis, cell_shift, i_index_real, i_size,num
  integer:: i_index_new,  i_index_old, i_index, i_basis_2,  i_basis_1
  integer:: i_spin

  select case( packed_matrix_format)

  case(PM_index) !---------------------------------------------------------------------------


     if(real_eigenvectors)then

        dipolemat_full_w = 0.0
     else
        dipolemat_full_w_complex = (0.0,0.0)
     end if

     do i_cell = 1, n_cells_in_hamiltonian-1

        do i_basis_2 = 1, n_basis

           do i_basis_1 = 1, n_basis

                 num = (i_cell-1)*(n_basis*n_basis)+(i_basis_2-1)*n_basis+i_basis_1
                 if(real_eigenvectors)then


                    dipolemat_full_w(i_basis_1,i_basis_2) =   &
                         dipolemat_full_w(i_basis_1,i_basis_2) &
                         + dble(k_phase_ac( i_cell))  &
                         * dipole_mat_full_oned( num)

                 else ! complex eigenvectors
                    dipolemat_full_w_complex(i_basis_1,i_basis_2) =  &
                         dipolemat_full_w_complex(i_basis_1,i_basis_2) &
                         + k_phase_ac( i_cell)  &
                         * dipole_mat_full_oned(  num )
                 end if  ! complex eigenvectos ends
           end do  ! i_size
        end do  ! i_basis_2
     end do  ! i_cell

  case(PM_none)  !-----------------------------------------------------------------------------

  case default
     
     call localorb_info('Invalid packing type')
     call aims_stop

  end select  ! packed matrix format
end subroutine construct_dipolemat_p1
!******
subroutine calc_dipelement_p0(dipelement, dipolemat_full_w, dipolemat_full_w_complex, KS_vec, KS_vec_complex, k_point)


  !  PURPOSE
  !   Calculete Momentum-Matrix elements by summing up of KS EV coefficents (basis transformation)
  !
  ! USES
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  use mpi_tasks, only: aims_stop
  use localorb_io, only: localorb_info
  implicit none
  !  ARGUMENTS

  complex*16, INTENT(IN)::     KS_vec_complex         (n_basis, n_states, n_spin)
  real*8, INTENT(IN)::     KS_vec         (n_basis, n_states, n_spin)
  real*8, INTENT(IN)::     dipolemat_full_w         (n_basis,n_basis)
  complex*16, INTENT(IN):: dipolemat_full_w_complex (n_basis,n_basis)
  integer, INTENT(IN):: k_point
  complex*16, dimension((n_states+1)*n_states/2), INTENT(OUT) ::     dipelement
  !  INPUTS
  !    o KS_eigenvector_complex
  !    o KS_eigenvector
  !    o dipolemat_w 
  !    o dipolemat_w_complex
  !    o k_point -- k-point wanted to be calculated 
  !  OUTPUT
  !    o dipelement
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: i_spin,i_basis_2, i_basis_1, n_state, m_state,i_basis,i_state
  complex*16:: dummy(n_basis) 
  real*8:: dummy_real(n_basis) 
  integer::  num
  select case( packed_matrix_format)

  case(PM_index)
     i_spin = 1
     num = 0
     dipelement = (0.0,0.0)
     do n_state = n_state_min, n_state_max
       do m_state = n_state, n_state_max
          num = num + 1
          dipelement(num) = (0.0,0.0)
          if(real_eigenvectors)then
	      dummy_real=0.0
	      call dgemm('N', 'N', 1, n_basis, n_basis, (1.d0,0.d0), &
			    (KS_vec(1:n_basis,n_state,i_spin)), 1 , &
			    dipolemat_full_w, n_basis, (0.d0,0.d0), &
			    dummy_real, 1)
	      call dgemm('N', 'N', 1, 1, n_basis, (1.d0,0.d0), &
			    dummy_real, 1 , &
			    KS_vec(1:n_basis,m_state,i_spin), n_basis, (0.d0,0.d0), &
			    dipelement(num), 1)
         else
	      dummy=(0.0,0.0)
              call zgemm('N', 'N', 1, n_basis, n_basis, (1.d0,0.d0), &
			conjg(KS_vec_complex(1:n_basis,n_state,i_spin)), 1 , &
			 dipolemat_full_w_complex, n_basis, (0.d0,0.d0), &
			 dummy, 1)
	      call zgemm('N', 'N', 1, 1, n_basis, (1.d0,0.d0), &
			 dummy, 1 , &
			 (KS_vec_complex(1:n_basis,m_state,i_spin)), n_basis, (0.d0,0.d0), &
			 dipelement(num), 1)
	 endif
      end do
    end do
  case default
     
     call localorb_info('Invalid packing type')
     call aims_stop

  end select
end subroutine calc_dipelement_p0

subroutine get_state_minmax(KS_eigen)

  !  PURPOSE
  !   Get the highest and lowest states at Gamma point that are smaller/bigger than E_max/E_min 
  !
  ! USES
    use dimensions, only: n_states, n_spin, n_k_points
    use pbc_lists
    use runtime_choices
    implicit none
    real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) ::     KS_eigen
    real*8 :: Emin_ha
    real*8 :: Emax_ha
    integer :: i, i_k_point, i_state, i_spin
    integer :: n_homo_k(n_k_points,n_spin)
    n_state_min=1
    n_state_max=1
    Emin_ha=Emin/hartree
    Emax_ha=Emax/hartree
    do i = 1,n_states
       if (KS_eigen(i,1,1)<=Emin_ha) then
          n_state_min=n_state_min+1
       endif
       if (KS_eigen(i,1,1)<=Emax_ha) then
          n_state_max=n_state_max+1
       endif
    enddo
    if (n_state_min>n_states) then
       n_state_min=n_states-1
    endif
    if (n_state_max>n_states) then
       n_state_max=n_states
    endif
    if (n_state_min>=n_state_max) then
       n_state_min=n_state_max-1
    endif

end subroutine get_state_minmax

subroutine out_dipelement(dipelement_one,dipelement_two,dipelement_three, KS_Eigenvalue, k_point)

  !  PURPOSE
  !   Write Momentummatrixelements (x,y,z) to file for one k-point
  !
  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  !  ARGUMENTS
  complex*16, dimension((n_states+1)*n_states/2), INTENT(IN) ::     dipelement_one
  complex*16, dimension((n_states+1)*n_states/2), INTENT(IN) ::     dipelement_two
  complex*16, dimension((n_states+1)*n_states/2), INTENT(IN) ::     dipelement_three
  integer, INTENT(IN) :: k_point
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) ::     KS_eigenvalue
  ! Write dipelement to file
  !  INPUTS
  !    o dipelement_i -- n_states+1 x n_states/2 Momentummatrix
  !    o k_point -- k-point wanted to be calculated 
  !    o  KS_Eigenvalue
  !  OUTPUT
  !    o file dipelement_k_point.dat
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: n_state, m_state, num, ierror
  CHARACTER(len=5) :: value
  CHARACTER(len=1) :: valuek1
  CHARACTER(len=2) :: valuek2
  CHARACTER(len=3) :: valuek3
  CHARACTER(len=4) :: valuek4
  CHARACTER(len=5) :: valuek5
  CHARACTER(len=55) :: fmt
  CHARACTER(len=25) :: name
  CHARACTER(len=25) :: name2
  real*8 :: PP

  write(value,'(I3)') 7
  fmt = '(I3, F10.4, I3, F10.4, ' // value // 'ES14.4)'
  if (k_point<=9) then
     write(valuek1,'(I1)') k_point
     name = 'element_k_'//trim(valuek1)//'.dat'
     name2 = 'element_k_'//trim(valuek1)
  elseif (k_point<=99) then
     write(valuek2,'(I2)') k_point
     name = 'element_k_'//trim(valuek2)//'.dat'
     name2 = 'element_k_'//trim(valuek2)
  elseif (k_point<=999) then
     write(valuek3,'(I3)') k_point
     name = 'element_k_'//trim(valuek3)//'.dat'
     name2 = 'element_k_'//trim(valuek3)
  elseif (k_point<=9999) then
     write(valuek4,'(I4)') k_point
     name = 'element_k_'//trim(valuek4)//'.dat'
     name2 = 'element_k_'//trim(valuek4)
  elseif (k_point<=99999) then
     write(valuek5,'(I5)') k_point
     name = 'element_k_'//trim(valuek5)//'.dat'
     name2 = 'element_k_'//trim(valuek5)
  endif

  open(unit=9, file=name,ACTION='WRITE')
     WRITE (9,120) k_point, k_point_list(k_point,1:3)
     120 FORMAT ('# K-point:',I5,' at', 3F10.6)
     WRITE (9,121)
     121 FORMAT ('# KS state i, KS energy i [eV], KS state j, ', &
        'KS energy j [eV], Re(p^x_ij), Im(p^x_ij), Re(p^y_ij), Im(p^y_ij), ', &
        'Re(p^z_ij), Im(p^z_ij), (p^x_ij^2+p^y_ij^2+p^z_ij^2)/3')
     num = 0
     do n_state = n_state_min,n_state_max
       do m_state = n_state, n_state_max
          num = num + 1
          PP = dble((dipelement_one(num)*conjg(dipelement_one(num))+dipelement_two(num)* &
               conjg(dipelement_two(num))+dipelement_three(num)*conjg(dipelement_three(num)))) 
          WRITE (9,fmt) n_state, KS_eigenvalue(n_state, 1, k_point)*hartree, m_state, KS_eigenvalue(m_state, 1, k_point)*hartree, &
                        dipelement_one(num),dipelement_two(num),dipelement_three(num), (PP)/3.
      end do
     end do
  close(unit=9)

end subroutine out_dipelement

!-----------------------------------------------------
!****s* FHI-aims/calculate_dipolemat_p0
!  NAME
!   calculate_dipolemat_p0
!  SYNOPSIS

subroutine calculate_dipolemat_p0( partition_tab_std, basis_l_max, dipole_mat_full_oned, direction, count1)

!  PURPOSE
!  The subroutine integrates the dipole matrix
!  using a fixed basis set (no adaptive modifications).
!  
!  Adapted from integrate_olvp_matrix and integrate_hamitonian_matrix
!
!  USES
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use plus_u
  use mpi_utilities
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
  use constants
  use species_data, only: species_name
  use load_balancing
  use pbc_lists, only: centers_basis_integrals, inv_centers_basis_integrals, &
      n_cells_in_hamiltonian
  use physics, only: deallocate_vdw_splines
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points) :: partition_tab_std
  integer basis_l_max (n_species)
  integer coord
  integer count1
  character direction
  ! when this routine is called, dipole_mat has either the dimension
  ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
  ! so we declare it here as a 1D assumed size array
  real*8 :: dipole_mat_full_oned(*)


!  INPUTS
!  o partition_tab_std -- values of partition function
!  o basis_l_max -- maximum l component of basis functions
!
!  OUTPUT
!  o overlap_matrix -- overlap matrix
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
  real*8, dimension(:,:,:),allocatable :: work_ham
  real*8, dimension(:,:),allocatable :: work_ovl

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab
  real*8, dimension(:), allocatable :: dist_tab_full
  real*8, dimension(:,:), allocatable :: dir_tab_full_norm
  real*8, dimension(:), allocatable :: i_r_full

  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:,:)  ,allocatable:: wave
  real*8, dimension(:,:,:), allocatable :: gradient_basis_wave
  real*8 coord_current(3)
  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 i_r(n_centers_integrals)
  real*8 dir_tab(3, n_centers_integrals, n_max_batch_size)
  real*8 trigonom_tab(4, n_centers_integrals)


  !    Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
  !     The hope is that such a separate treatment will allow to minimize numerical noise
  !     introduced through ZORA
  real*8, dimension(:,:), allocatable :: dipolemat_shell

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points

  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  !     for pruning of atoms, radial functions, and basis functions, to only the relevant ones ...

  integer :: n_compute_a, n_compute_c
  integer :: i_basis(n_centers_basis_I)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_max_compute_atoms)
  integer :: spline_array_end(n_max_compute_atoms)

! VB - renewed index infrastructure starts here

      real*8 one_over_dist_tab(n_max_compute_atoms)

      ! indices for basis functions that are nonzero at current point
!
! VB: Warning - n_max_compute_fns_ham is outrightly wrong for the density
!     update!!
!
!     This should MUCH rather be counted in prune_basis for each batch and then
!     passed here, that would be the appropriate maximum!
!

      integer :: rad_index(n_max_compute_atoms)
      integer :: wave_index(n_max_compute_fns_ham)
      integer :: l_index(n_max_compute_fns_ham)
      integer :: l_count(n_max_compute_fns_ham)
      integer :: fn_atom(n_max_compute_fns_ham)

      ! indices for known zero basis functions at current point
      integer :: n_zero_compute
      integer :: zero_index_point(n_max_compute_ham)

      ! active atoms in current batch
      integer :: n_batch_centers
      integer :: batch_center(n_centers_integrals)

  !     for splitting of angular shells into "octants"

  integer division_low
  integer division_high

  !  counters

  integer i_basis_1
  integer i_basis_2
  integer i_cell
  integer i_grid
  integer i_index, i_l, i_m
  integer i_coord
  integer i_center, i_center_L
  integer i_division

  integer i_species

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2

  integer :: i_my_batch

  character*200 :: info_str
  integer :: info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  integer dim_overlap_matrix ! actual dimension of overlap matrix

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i, j, i_off, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all


  !  begin work
select case (direction) 
     case ('a')
        coord=count1
     case ('x') 
       coord=1
     case ('y') 
       coord=2
     case ('z') 
       coord=3
endselect



  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A,A)')"Integrating Momentum Matrix with loadbalancing"
  else
    write(info_str,'(2X,A,A)')"Integrating momentum matrix."
  endif
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  !     begin with general allocations
  l_ylm_max = l_wave_max

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) )
  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )
  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')
  allocate(radial_wave_deriv(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_deriv             ')

  allocate(wave(n_max_compute_ham, n_max_batch_size), STAT=info )
  call check_allocation(info, 'wave                          ')

  allocate (gradient_basis_wave(n_max_compute_ham,3 ,n_max_batch_size),STAT=info)
  call check_allocation(info, 'gradient_basis_wave           ')

  allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
  call check_allocation(info, 'dylm_dtheta_tab               ')

  allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
  call check_allocation(info, 'scaled_dylm_dphi_tab          ')

  allocate(dist_tab_full(n_centers_integrals),STAT=info )
  call check_allocation(info, 'dist_tab_full                 ')

  allocate(dir_tab_full_norm(3,n_centers_integrals),STAT=info )
  call check_allocation(info, 'dir_tab_full_norm             ')

  allocate(i_r_full(n_centers_integrals),STAT=info )
  call check_allocation(info, 'i_r_full                      ')

  allocate ( dipolemat_shell(n_max_compute_ham,n_max_compute_ham), stat=info )
  call check_allocation(info,'dipolemat_shell','integrate_ovlp_matrix_p2')

  allocate(work_ham(1, 1, n_spin),stat=info)
  call check_allocation(info, 'work_ham                      ')
  allocate(work_ovl(1, 1),stat=info)
  call check_allocation(info, 'work_ovl                      ')
  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))

    dim_overlap_matrix = batch_perm(n_bp)%n_local_matrix_size

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std

    dim_overlap_matrix = n_hamiltonian_matrix_size

  endif

  if(get_batch_weights) allocate(batch_times(n_my_batches_work))

  !-----------------------------------------------------------------------------

  ! initialize

  dipole_mat_full_oned(1:n_cells_in_hamiltonian*n_basis*n_basis) = 0.d0
  i_basis_fns_inv = 0
  gradient_basis_wave=0
  ! initialize index_lm

  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo

  i_full_points = 0
  i_full_points_2 = 0

  call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  do i_my_batch = 1, n_my_batches_work, 1

        if(get_batch_weights) time_start = mpi_wtime()

        n_compute_c = 0
        n_compute_a = 0
        i_basis = 0

        i_point = 0

        ! loop over one batch
        do i_index = 1, batches_work(i_my_batch)%size, 1

           i_full_points_2 = i_full_points_2 + 1

           if (partition_tab(i_full_points_2).gt.0.d0) then

              i_point = i_point+1

              !     get current integration point coordinate
              coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if

              !     compute atom-centered coordinates of current integration point,
              !     as viewed from all atoms
              call tab_atom_centered_coords_p0( coord_current, dist_tab_sq(1,i_point), &
                   dir_tab(1,1,i_point), n_centers_integrals, centers_basis_integrals )

              !    determine which basis functions are relevant at current integration point,
              !     and tabulate their indices

              ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
              if (.not.prune_basis_once) then
                 call prune_basis_p2( dist_tab_sq(1,i_point), n_compute_c, &
                      i_basis, n_centers_basis_I, n_centers_integrals, &
                      inv_centers_basis_integrals )
              end if
           end if

        enddo ! end loop over one batch

        if (prune_basis_once) then
           n_compute_c = batches_work(i_my_batch)%batch_n_compute
           i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
        end if
        ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
        call collect_batch_centers_p2 &
        ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, &
          n_batch_centers, batch_center &
        )

        n_points = i_point

        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_c.gt.0) then

           i_point = 0

           ! loop over one batch of integration points
           do i_index = 1, batches_work(i_my_batch)%size, 1

              ! Increment the (global) counter for the grid, to access storage arrays
              i_full_points = i_full_points + 1

              if (partition_tab(i_full_points).gt.0.d0) then

                 i_point = i_point+1

                 ! for all integrations
                 partition(i_point) = partition_tab(i_full_points)

                 n_compute_atoms = 0
                 n_compute_fns = 0

                 ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                 ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
                 ! without any copying and without doing any unnecessary operations.
                 ! The price is that the interface is no longer explicit in terms of physical
                 ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                 call tab_global_geometry_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dir_tab(1,1,i_point), &
                      dist_tab_full, &
                      i_r_full, &
                      dir_tab_full_norm, &
                      n_centers_integrals,  centers_basis_integrals)

                 call prune_radial_basis_p2 &
                   ( n_max_compute_atoms, n_max_compute_fns_ham, &
                     dist_tab_sq(1,i_point), dist_tab(1,i_point), dir_tab(1,1,i_point), &
                     n_compute_atoms, atom_index, atom_index_inv, &
                     n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                     i_atom_fns, spline_array_start, spline_array_end, &
                     n_centers_integrals, centers_basis_integrals, &
                     n_compute_c, i_basis, &
                     n_batch_centers, batch_center, &
                     one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                     fn_atom, n_zero_compute, zero_index_point &
                    )

                 ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                 ! for all atoms which are actually relevant
                 call tab_local_geometry_p2 &
                      ( n_compute_atoms, atom_index, &
                        dist_tab(1,i_point),  &
                      i_r )

                 ! compute trigonometric functions of spherical coordinate angles
                 ! of current integration point, viewed from all atoms
                 call tab_trigonom_p0 &
                      ( n_compute_atoms, dir_tab(1,1,i_point),  &
                      trigonom_tab )

                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm_p0 &
                      ( n_compute_atoms, atom_index,  &
                      trigonom_tab, basis_l_max,  &
                      l_ylm_max, &
                      ylm_tab )

                 ! Now evaluate radial functions
                 ! from the previously stored compressed spline arrays
                 call evaluate_radial_functions_p0  &
                      ( spline_array_start, spline_array_end,  &
                      n_compute_atoms, n_compute_fns,   &
                      dist_tab(1,i_point), i_r,  &
                      atom_index, i_basis_fns_inv,  &
                      basis_wave_ordered, radial_wave,  &
                      .false. , n_compute_c, n_max_compute_fns_ham )

                 ! tabulate total wave function value for each basis function
                 call evaluate_waves_p2  &
                   ( n_compute_c, n_compute_atoms, n_compute_fns, &
                     l_ylm_max, ylm_tab, one_over_dist_tab,   &
                     radial_wave(1), wave(1,i_point), &
                     rad_index, wave_index, l_index, l_count, fn_atom, &
                     n_zero_compute, zero_index_point &
                   )

                 call tab_gradient_ylm_p0  &
                      ( trigonom_tab(1,1), basis_l_max,   &
                      l_ylm_max, n_compute_atoms, atom_index,  &
                      ylm_tab(1,1),   &
                      dylm_dtheta_tab(1,1),   &
                      scaled_dylm_dphi_tab(1,1)  )
                 call evaluate_radial_functions_p0  &
                      ( spline_array_start, spline_array_end,  &
                      n_compute_atoms, n_compute_fns,   &
                      dist_tab(1,i_point), i_r,  &
                      atom_index, i_basis_fns_inv,  &
                      basis_deriv_ordered,   &
                      radial_wave_deriv(1), .true.,  &
                      n_compute_c, n_max_compute_fns_ham )

                 ! and finally, assemble the actual gradients
                 call evaluate_wave_gradient_p2  &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
                   l_ylm_max, ylm_tab,  &
                   dylm_dtheta_tab,  &
                   scaled_dylm_dphi_tab,  &
                   radial_wave,  &
                   radial_wave_deriv,  &
                   gradient_basis_wave(1:n_compute_c,1:3,i_point),  &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point &
                 )

                 ! Reset i_basis_fns_inv
                 i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

              end if

              ! end loop over one batch
           enddo

           !<psi|grad_i_coord|psi>
           !do i_coord = 1, 3
	      call evaluate_dipole_mat( & 
		  n_points, partition, n_compute_c, gradient_basis_wave, &
		  n_max_compute_ham, wave, dipolemat_shell,coord)

              !keep upper and lower Triangle
	      call update_full_matrix_p0XXX &
		      ( n_compute_c, n_compute_c,  &
		      i_basis, &
		      dipolemat_shell,  &
		      dipole_mat_full_oned)
           !enddo
        else
           i_full_points = i_full_points + batches_work(i_my_batch)%size
        end if ! n_compute.gt.0

        if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  enddo ! end loop over batches

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_global,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integration: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)

  !     synchronise the hamiltonian
  if(.not. use_local_index) call sync_vector(dipole_mat_full_oned, n_cells_in_hamiltonian*n_basis*n_basis )

  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)


  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if

  if (allocated(dipolemat_shell)) then
     deallocate( dipolemat_shell )
  end if

  if (allocated(ins_idx)) then
     deallocate(ins_idx)
  endif

  if (allocated(batch_times)) then
     deallocate(batch_times)
  endif

  if (allocated(radial_wave)) then
     deallocate(radial_wave)
  endif

  if (allocated(radial_wave_deriv)) then
     deallocate(radial_wave_deriv)
  endif

  if (allocated(wave)) then
     deallocate(wave)
  endif

  if (allocated(gradient_basis_wave)) then
     deallocate(gradient_basis_wave)
  endif

  if (allocated(dylm_dtheta_tab)) then
     deallocate(dylm_dtheta_tab)
  endif

  if (allocated(scaled_dylm_dphi_tab)) then
     deallocate(scaled_dylm_dphi_tab)
  endif

  if (allocated(dist_tab_full)) then
     deallocate(dist_tab_full)
  endif

  if (allocated(dir_tab_full_norm)) then
     deallocate(dir_tab_full_norm)
  endif

  if (allocated(i_r_full)) then
     deallocate(i_r_full)
  endif

end subroutine calculate_dipolemat_p0
subroutine update_full_matrix_p0XXX &
     ( n_compute_c, n_compute_a, i_basis, &
     matrix_shell, matrix)

!  PURPOSE
!  Subroutine update_full_matrix adds a part of the integrals in a
!  matrix (overlap or Hamiltonian matrix)
!  (only for the n_compute basis functions that are nonzero at the
!  current integration shell) to the full matrix (which contains
!  all n_basis basis functions). The link between i_compute = 1 ... n_compute
!  and i_basis = 1 ... n_basis is provided by the index array 
!  i_basis(i_compute).
!
!  Modified to keep upper and lower triangle of the matrix_shell
!
!  USES

  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: aims_stop
  use localorb_io, only: localorb_info
  implicit none

!  ARGUMENTS

  integer n_compute_c
  integer n_compute_a
  integer i_basis(n_compute_c)
  real*8 matrix_shell(n_compute_c,n_compute_a)
  real*8 matrix (n_cells_in_hamiltonian*n_basis*n_basis)
!  INPUTS
!   o n_compute_c == n_compute_a -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
!   o matrix_shell -- hamiltonian / overlap_matrix of relevant basis functions
!
!  OUTPUT
!   o matrix -- date from matrix_shell is added here
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
!
  !  local variables

  !     auxiliary matrices for Level 3 Blas matrix multiplications

  !     counters

  integer :: i_compute_1
  integer :: i_compute_2
  integer :: i_offset
  integer :: i_index_real
  integer :: i_offset_first_part
  integer :: i_cell_index, i_cell_1
  integer :: i_max_basis, i_min_basis
  integer :: i_one_part, i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: offset(n_cells) !_in_hamiltonian)
  integer :: offset_end(n_cells) !_in_hamiltonian)
  integer :: help
!  integer :: direct_citing(n_compute_a,n_compute_a )
!  real*8 :: data(n_cells_in_hamiltonian, n_basis)
!  real*8,dimension(:,:), allocatable :: data
  integer:: i_cell_old, num

!NEC_CB
  integer::help1(n_compute_a)
  integer::help2(n_compute_a)

  !     begin work

  !      now add the aux. ham. matrix to the actual ham. matrix
  !      this requires translating between the actually computed matrix elements and the
  !      full ham. matrix ...
  !      Version to keep upper and lower triangle for momentummatrix calculation 
  !      <|psi_1|\nabla|\psi_2>\neq<|psi_2|\nabla|\psi_1>
  ! When basis is smaller than n_basis

  !      write(use_unit,*) n_compute_c,n_compute_a

  select case(packed_matrix_format)


  case(PM_index) !------------------------------------------------------------------------------------


        do i_compute_1 = 1, n_compute_a, 1
          help1(i_compute_1)=Cbasis_to_basis(i_basis(i_compute_1))
          help2(i_compute_1)=center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
        end do

        do i_compute_1 = 1, n_compute_a, 1

           i_basis_1 = help1(i_compute_1)!Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_old = help2(i_compute_1)!center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))



           do i_compute_2 = 1, n_compute_a
              
              i_basis_2 = help1(i_compute_2)!Cbasis_to_basis(i_basis(i_compute_2))
              i_cell    =  position_in_hamiltonian( i_cell_old, help2(i_compute_2) )!center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))
              num = (i_cell-1)*(n_basis*n_basis)+(i_basis_1-1)*n_basis+i_basis_2
              matrix(num) = matrix(num) + matrix_shell(i_compute_2, i_compute_1)
           end do
        end do
     

  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select


end subroutine update_full_matrix_p0XXX
end module calculate_dipolemat

