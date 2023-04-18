!  Module control contains all necessary control parameters
!  
!  R.Gehrke (2005)
!

module control

  character*30 :: simulation_flag
  character*30 :: diff_norm
  real*8 :: hard_sphere_radius
  character*30 :: potential_flag
  logical :: use_pair_pot
  real*8 :: energy_interval
  logical :: spin_polarized
  logical :: charged
  logical :: verbose

end module control
