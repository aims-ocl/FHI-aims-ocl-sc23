!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/integrate_free_atom_hessian_p1
!  NAME
!   integrate_free_atom_hessian_p1
!  SYNOPSIS

subroutine integrate_free_atom_hessian_p1 &
     (hellman_feynman_hessian_free_part)

!  PURPOSE
!
! get d2V_free(|R1-R2|)/dR1x dR2x  for new_hessian_calculation

! note by shanghui:
!                                                    d V^{I}_{free}(R_J-R_{I+R}) 
! (1) free_force:         F(J) =   Z_J * \sum_{I,R}  --------------------------
!                                                          d R_J 
!

!                                            d^2  V^{I}_{free}(R_{J+R}-R_I)
! (2) Hessian: Hessian(J+R, I) = - Z_J * ----------------------------------
!                                                     d R_{J+R} d R_I

 
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use free_atoms
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use pbc_lists
  use load_balancing

  implicit none

!  ARGUMENTS


  real*8, dimension(3, n_centers_in_sc_DFPT, 3, n_atoms)          ::  hellman_feynman_hessian_free_part

!  INPUTS
!
!  OUTPUT
! o hellman_feynman_hessian_free_part : d2V_free(|R1-R2|)/dR1x dR2x
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
  integer i_center, j_center
  integer i_atom, j_atom
  integer current_center_i, current_center_j
  integer i_coord, j_coord

  real*8 coord_current(3)
  real*8 dist_tab_sq
  real*8 dist_tab_in
  real*8 dir_tab(3)
  real*8 dir_tab_in(3)
  real*8 log_weight
  real*8 radial_weight
  real*8 trigonom_tab(4)
  real*8 i_r
  real*8 i_r_log

  real*8 :: d_v_hartree_free_d_r,hartree_free_2nd_deriv


  !  external functions
  character*100 :: info_str
  integer :: mpierr
  integer :: info



  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all


  hellman_feynman_hessian_free_part = 0.0d0


  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()


    do i_center = 1, n_atoms  

      current_center_i   = center_in_sc_DFPT_to_center(i_center)
      i_atom = center_in_sc_DFPT_to_atom(i_center)

          do j_center = 1, n_centers_in_sc_DFPT, 1

             current_center_j = center_in_sc_DFPT_to_center(j_center)
             j_atom         = center_in_sc_DFPT_to_atom(j_center)

           do i_coord = 1, 3, 1
              coord_current(i_coord) = coords_center(i_coord, current_center_j)
            end do

            call tab_single_atom_centered_coords_p0 &
                 ( current_center_i, coord_current,  &
                 dist_tab_sq,  &
                 dir_tab )


            if (j_center.eq. i_center ) then
              dist_tab_sq = r_grid_min(species(j_atom))**2 + 1e-15
            end if


         !-----keep the same as free_rho_superpos-----
            if(dist_tab_sq.lt.multipole_radius_free_sq(species(i_atom))) then
             

               call tab_single_atom_centered_coords_radial_log_p0 &
                    ( current_center_i, dist_tab_sq, dir_tab,  &
                    dist_tab_in, i_r, i_r_log, dir_tab_in )

               call tab_single_radial_weights_v2 &
                    ( i_atom, dist_tab_in, i_r, &
                      log_weight, radial_weight )


                  ! splines are now derivatives df/di where i is the grid point index
                  ! must convert to df/dr = df/di * di/dr

                  ! obtain gradients of the free atom density and potential explicitly

                  ! radial derivative of hartree potential of free atom # i_atom (without core potential!!!)
                  ! subtract core potential => - ( - Z/r ) = + Z/r , d/dr => - Z/r^2

               if(j_center == i_center) then

                   d_v_hartree_free_d_r =  &
                   val_spline_deriv( i_r_log, &
                   free_pot_es_spl(1,1,species(i_atom)),  &
                   n_grid(species(i_atom))) * log_weight

              
                  hellman_feynman_hessian_free_part(1:3,i_center,1:3,j_center) = 0.0d0

               else
                 
                 !-----------shanghui add for new_hellman_feynman_hessian-----
                  d_v_hartree_free_d_r =  &
                          val_spline_deriv( i_r_log, &
                          free_pot_es_spl(1,1,species(i_atom)),  &
                          n_grid(species(i_atom))) * log_weight

                 ! note: free_pot_es use n_spine_coeff=4 ===>  cubic spline
                 ! hartree_free_2nd_deriv = d^2 V_free(|r-R|)/ d|r-R|^2
                  hartree_free_2nd_deriv= val_spline_2nd_deriv(i_r_log, &
                        free_pot_es_spl(1,1,species(i_atom)),  &
                        n_grid(species(i_atom))) * log_weight*log_weight &
                        - val_spline_deriv( i_r_log, &
                        free_pot_es_spl(1,1,species(i_atom)),  &
                        n_grid(species(i_atom))) * log_weight/dist_tab_in  

                  do i_coord = 1,3 
                  do j_coord = 1,3 

                   if(i_coord.eq.j_coord) then 
                   hellman_feynman_hessian_free_part(j_coord,j_center,i_coord,i_center)= & 
                   - species_z(species(j_atom)) * ( d_v_hartree_free_d_r* & 
                    (-1.0d0+dir_tab_in(i_coord)*dir_tab_in(j_coord))/dist_tab_in &
                     -hartree_free_2nd_deriv*dir_tab_in(i_coord)*dir_tab_in(j_coord) )
                  !----here '-' is hessian_free_part = - d^2V_free/dRIdRJ 
                   else 
                   hellman_feynman_hessian_free_part(j_coord,j_center,i_coord,i_center)= & 
                   - species_z(species(j_atom)) * ( d_v_hartree_free_d_r* & 
                    (dir_tab_in(i_coord)*dir_tab_in(j_coord))/dist_tab_in &
                     - hartree_free_2nd_deriv*dir_tab_in(i_coord)*dir_tab_in(j_coord) )    
                  !----here '-' is hessian_free_part = - d^2V_free/dRIdRJ 
                   endif
 
                  enddo
                  enddo

               end if ! j_atom = i_center


          ! for free_atom_hessian, only near part. 
            end if 


          end do ! loop over j_atom
         !end if
      end do  ! loop over i_center




  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for free_atom_hessian: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

end subroutine  integrate_free_atom_hessian_p1
!******
!------------------------------------------------------------------------------------------------------------
