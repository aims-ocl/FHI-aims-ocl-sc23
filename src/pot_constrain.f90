      MODULE POTCONSTRAIN
      use dimensions
      use dimensions
      use geometry
      use physics
      use species_data
      use runtime_choices
      use mpi_tasks
      IMPLICIT NONE
      CONTAINS
  
      subroutine apply_constrain(free_energy)
      implicit none
      real*8 :: constrain_forces(3, n_atoms), constrain_pot, free_energy
      character(*), parameter :: func = 'apply_constrain'
      free_energy=0.d0
      if (potconstype.eq.'nanotubecylinder') then
         call nanotubecylinder(constrain_pot, constrain_forces)
      else if (potconstype.eq.'sphere') then
         call sphere(potconspar, constrain_pot, constrain_forces)
      else
         call aims_stop_coll('Constraining potential not implemented', func)
      endif
      total_forces=total_forces+constrain_forces
      free_energy=free_energy+constrain_pot
      end subroutine

      subroutine nanotubecylinder(pot, forces)
      ! Cylinder along x direction that mimics a 6,6 carbon nanotube
      ! Analytical form: V(r)=a2*r^2+a4*r^4+a6*r^6+a8*r^8
      ! and parameters following Dellago and Naor, CPC 169 (2005) 36--39
      implicit none
      real*8 :: a2, a4, a6, a8
      real*8  :: pot, forces(3, n_atoms), r
      integer i_atom, i_dim
      a2=-0.000101790427
      a4=0.0001362104651
      a6=8.1919580588d-06
      a8=3.188645093e-06
      pot=0.d0
      forces=0.d0
      do i_atom=1, n_atoms
         if (species_name(species(i_atom)).eq.'O') then
            call cylindrical_distance(coords(:,i_atom), r)
            pot=pot+a2*(r**2)+a4*(r**4)+a6*(r**6)+a8*(r**8)
!            forces(i_atom, 1)=0.d0
            if (r.gt.0.d0) then
               do i_dim=2, 3
                  forces(i_dim,i_atom)=-(2*a2*r+4*a4*(r**3)+6*a6*(r**5)+8*a8*(r**7))*coords(i_dim, i_atom)/r
               enddo
            endif
         endif
      enddo
      end subroutine

      subroutine sphere(radius, pot, forces)
      ! Cylinder along x direction that mimics a 6,6 carbon nanotube
      ! Analytical form: V(r)=a2*r^2+a4*r^4+a6*r^6+a8*r^8
      ! and parameters following Dellago and Naor, CPC 169 (2005) 36--39
      implicit none
      real*8 :: radius
      real*8  :: pot, forces(3, n_atoms), r
      integer i_atom, i_dim
      real*8 :: kT, a2 

      kT=0.00095 ! kT at 300K
      a2=4.0*kT
      ! This potential simply rises to 10kT over a distance of 0.5 Bohr beyond the radius that is set.
      ! It is zero between -R and R.
      ! It is just a random parametrization, feel free to make it more sensible to your needs
      ! It is also applied only to heavy atoms (no hydrogens)
      pot=0.d0
      forces=0.d0
      do i_atom=1, n_atoms
         if (species_name(species(i_atom)).ne.'H') then
            call spherical_distance(coords(:,i_atom), r)
            if (r.gt.-radius .and. r.lt.radius) then
               pot=pot+0.0
               forces(:,i_atom)=0.0
            else if (r.le.-radius) then
               pot=pot+a2*(r+radius)**2
               do i_dim=1, 3
                  forces(i_dim,i_atom)=-(2*a2*(r+radius))*coords(i_dim, i_atom)/r
               enddo
            else if (r.ge.radius) then
               pot=pot+a2*(r-radius)**2
               do i_dim=1, 3
                  forces(i_dim,i_atom)=-(2*a2*(r-radius))*coords(i_dim, i_atom)/r
               enddo
            endif
         endif
      enddo
      end subroutine

      subroutine cylindrical_distance(q, r)
      implicit none
      real*8 :: q(3), r
         r=dsqrt(q(2)**2 + q(3)**2)
      end subroutine 
      
      subroutine spherical_distance(q, r)
      implicit none
      real*8 :: q(3), r
         r=dsqrt(q(1)**2 + q(2)**2 + q(3)**2)
      end subroutine 
      
      END MODULE
