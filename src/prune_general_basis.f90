!****s* FHI-aims/prune_general_basis
!  NAME
!    prune_general_basis
!  SYNOPSIS

subroutine prune_general_basis(n_points, points, &
&                        n_centers, center2coords, center2species, &
&                        n_species, max_bas_L, max_n_bas_fnLsp, n_max_comp, &
&                        Lsp2n_bas_fnLsp, Lsp2bas_fn, Lsp2bas_sp, &
&                        n_bas_fns, basfn2radius, &
&                        n_comp, comp2center, comp2bas_sp, &
&                        comp2fn, comp2l, comp2m)

  !  PURPOSE
  !
  !     Figure out which bas functions reach at least one of the points.
  !
  !     Even if the large number of parameters makes this procedure look as if
  !     it is hard to use, it is not so bad in practice, as almost any of the
  !     input parameters is readily available in one or the other core module.
  !
  !  USES

  use localorb_io, only: use_unit
  use mpi_tasks
  use geometry !For debug purpose
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  integer, intent(IN) :: n_centers
  real*8, intent(IN) :: center2coords(3, n_centers)
  integer, intent(IN) :: center2species(n_centers)
  integer, intent(IN) :: n_species, max_bas_L, max_n_bas_fnLsp, n_max_comp
  integer, intent(IN) :: Lsp2n_bas_fnLsp(0:max_bas_L, n_species)
  integer, intent(IN) :: Lsp2bas_fn(max_n_bas_fnLsp, 0:max_bas_L, n_species)
  integer, intent(IN) :: Lsp2bas_sp(max_n_bas_fnLsp, 0:max_bas_L, n_species)
  integer, intent(IN) :: n_bas_fns
  real*8, intent(IN) :: basfn2radius(n_bas_fns)
  integer, intent(OUT) :: n_comp
  integer, intent(OUT) :: comp2center(n_max_comp)
  integer, intent(OUT) :: comp2bas_sp(n_max_comp)
  integer, intent(OUT) :: comp2fn(n_max_comp)
  integer, intent(OUT) :: comp2l(n_max_comp)
  integer, intent(OUT) :: comp2m(n_max_comp)

  !  INPUTS
  !   Grid points
  !    o n_points -- Number of grid points
  !    o points -- Grid point positions
  !   Centers (from pbc_lists.f90)
  !    o n_centers -- Number of atomic positions
  !    o center2coords -- Atomic positions  (pbc_lists:coords_center)
  !    o center2species -- Atomic species numbers (pbc_lists:species_center)
  !    o n_species -- Number of species
  !   Basis properties
  !         (replace "bas" with either "basis" [-> basis.f90]
  !                                 or "basbas" [-> prodbas.f90])
  !    o max_bas_L -- Maximum angular momentum
  !    o max_n_bas_fnLsp -- Maximum number of radial parts in L-channel
  !    o n_max_comp -- Array dimension
  !                    (very safe choice: n_centers * max_n_bas_fnLsp
  !                                       * (max_bas_L+1) * (2*max_bas_L+1))
  !    o Lsp2n_bas_fnLsp -- Number of radial parts i_fn in L-channel
  !    o Lsp2bas_fn -- (i_fn, L, i_species) -> i_bas_fn (see basis.f90)
  !    o Lsp2bas_sp -- (i_fn, L, i_species) -> i_bas_sp (see basis.f90)
  !    o n_bas_fns -- Total number of radial parts
  !    o basfn2radius -- Radii of radial parts (when are they zero)
  !                      (basis:outer_radius; prodbas:charge_radius_basbas_fn)
  !  OUTPUTS
  !    o n_comp -- Number of significant basis functions
  !    o comp2center -- Position of basis functions
  !    o comp2bas_sp -- Index of basis functions within center/atom
  !    o comp2fn, comp2l, comp2m -- Nature of basis function
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  !OTH: I'm adding now all basis functions which are within the unit cell

  real*8 :: species2radius(n_species)
  real*8 :: coord(3), dist, distvec(3), spec_rad, fn_rad
  real*8 :: point_min(3), point_max(3)
  logical :: is_in(max_n_bas_fnLsp, 0:max_bas_L)
  integer :: i_center, i_species, L, i_fnLsp, M
  integer :: i_bas_sp, i_bas_fn, i
  integer :: i_point
  integer :: n_comp_uptonow
  integer :: info
  character(*), parameter :: func = 'prune_general_basis'


  ! --- species2radius

  species2radius = 0.d0
  do i_species = 1, n_species
     do L = 0, max_bas_L
        do i_fnLsp = 1, Lsp2n_bas_fnLsp(L, i_species)
           i_bas_fn = Lsp2bas_fn(i_fnLsp, L, i_species)
           fn_rad = basfn2radius(i_bas_fn)
           species2radius(i_species) = max(species2radius(i_species), fn_rad)
        end do
     end do
  end do

  ! --- point_min, point_max

  do i = 1, 3
     point_min(i) = minval(points(i, :))
     point_max(i) = maxval(points(i, :))
     point_min(i) = min(point_min(i),sum(-lattice_vector(i,:))*0.5d0)
     point_max(i) = max(point_max(i),sum(lattice_vector(i,:))*0.5d0)
  end do

  ! --- center & point loops

  n_comp_uptonow = 0
  ! Loop over centers first to make sure that functions from the same center
  ! are subsequent.
  CENTER_LOOP: do i_center = 1, n_centers
     i_species = center2species(i_center)
     spec_rad = species2radius(i_species)
     coord = center2coords(:, i_center)
     ! Is there any point or are we completely off?
     if (any(coord + spec_rad < point_min)) cycle CENTER_LOOP
     if (any(coord - spec_rad > point_max)) cycle CENTER_LOOP

     ! Prepare is_in(:,:) to be .false. for all valid entry and
     ! .true. otherwise for the all(is_in)-check below to work.
     is_in = .true.
     do L = 0, max_bas_L
        do i_fnLsp = 1, Lsp2n_bas_fnLsp(L, i_species)
           is_in(i_fnLsp, L) = .false.
        end do
     end do

     POINT_LOOP: do i_point = 1, n_points
        distvec = coord - points(:, i_point)
        if (sum(distvec**2).lt.1e-8) then
          dist=0d0
        else
           dist = sqrt(sum(distvec**2))
        endif
        if (dist.ne.dist .or. dist.gt.1e8) then  !First part tests for NaN
                                                 !Second part thought to test for infinity
          write(use_unit,*) 'Dist:, ', dist
          call aims_stop('Dist has failed in prune_general_basis.f90')
        endif
        if (dist > spec_rad) cycle POINT_LOOP

        L_LOOP: do L = 0, max_bas_L
           do i_fnLsp = 1, Lsp2n_bas_fnLsp(L, i_species)
              if (is_in(i_fnLsp, L)) cycle !Checks whether this basis function has already been accounted for.
              i_bas_fn = Lsp2bas_fn(i_fnLsp, L, i_species)
              fn_rad = basfn2radius(i_bas_fn)
              if (dist .le. fn_rad) then
                 is_in(i_fnLsp, L) = .true.
                 i_bas_sp = Lsp2bas_sp(i_fnLsp, L, i_species)
                 ! Loop over M last to make sure that functions with same
                 ! center and i_bas_fn are subsequent.
                 do M = -L, L
                    n_comp_uptonow = n_comp_uptonow + 1
                    if (n_comp_uptonow > n_max_comp) then
                       call aims_stop('Too many basis functions', func)
                    end if
                    comp2center(n_comp_uptonow) = i_center
                    comp2bas_sp(n_comp_uptonow) = i_bas_sp + M
                    comp2fn(n_comp_uptonow) = i_bas_fn
                    comp2l(n_comp_uptonow) = L
                    comp2m(n_comp_uptonow) = M
                 end do
              end if
           end do
        end do L_LOOP
        if (all(is_in)) exit POINT_LOOP   ! Nothing to add any more.
     end do POINT_LOOP
  end do CENTER_LOOP
  n_comp = n_comp_uptonow

end subroutine prune_general_basis
!******
