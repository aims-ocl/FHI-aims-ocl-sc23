!****s* FHI-aims/get_occupation_numbers
!  NAME
!   get_occupation_numbers
!  SYNOPSIS

subroutine update_chemical_pot_cluster(inv_full_ovlp_matrix_sqrt,&
                                       full_ovlp_matrix_sqrt, &
                                       full_ovlp_matrix, &
                                       free_cluster_ham, hybrid_func,&
                                       loc_self_enrg, &
                                       hartree_pot_LDA,&
                                       KS_eigenvalue, &
                                       n_electrons, &
                                       t_out,&
                                       new_chem_pot,&
                                       embed_part_number_new,&
                                       k_summed_overlap_matr,&
                                       self_energy_freq)

! PURPOSE
! Subroutine get_occupation_numbers determines occupation numbers according to two available schemes:
! o 1) gaussian smearing according to Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
! o 2) fermi smearing
! o 3) methfessel-paxton (methfessel et al., PRB 40, 6, 40 (1989))
!
!  USES

  use dimensions
  use runtime_choices
  use localorb_io
  use constants
  use gw_para 
  use scgw_grid


  implicit none

!  ARGUMENTS

  real*8, dimension(n_states, n_spin,n_k_points), intent(in) :: KS_eigenvalue
  !real*8, dimension(n_states, n_spin), intent(out) :: occ_numbers
  real*8, intent(in)  :: n_electrons
  real*8, intent(out) :: new_chem_pot 
  real*8, intent(out) :: embed_part_number_new 
  logical :: t_out

   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: &
   inv_full_ovlp_matrix_sqrt
   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: &
   full_ovlp_matrix_sqrt
   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: &
   full_ovlp_matrix
   real*8, dimension(n_basis,n_basis), intent(in) :: free_cluster_ham
   real*8, dimension(n_basis,n_basis), intent(in) :: k_summed_overlap_matr
   real*8, dimension(n_basis,n_basis), intent(in) :: loc_self_enrg
   complex*16, dimension(n_basis,n_basis,nomega) :: hybrid_func
   complex*16, dimension(n_basis,n_basis,nomega) :: self_energy_freq
   real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA




!  INPUTS
!  o KS_eigenvalue -- Kohn-Sham eigenvalues
!  o n_electrons -- number of electrons
!  o t_out -- is the information printed out or not ?
!
!  OUTPUT
!  o occ_numbers -- occupation weights of different KS states
!  o chemical_potential -- chemical potential
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
  real*8 :: chemical_potential_l
  real*8 :: chemical_potential_r
  real*8 :: diff_l
  real*8 :: diff_r
  real*8 :: diff_electrons
  real*8 :: diff_electrons_thr
  real*8 :: lowest_eigenvalue
  real*8 :: highest_eigenvalue

  character*120 :: info_str

  !  counters
  integer :: i_state,  number_reiterations
  integer :: i_k_points 
  integer :: i_spin
  integer :: i_counter

!write(use_unit,*) 'n_electrons', n_electrons
  if (t_out) then
     write(info_str,'(2X,A)') &
          "Determining occupation numbers for Kohn-Sham eigenstates."
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if

  i_counter = 0
  !  write(use_unit,*) "get_occ_numbers.... n_spin=", n_spin, "n_states=", n_states


  lowest_eigenvalue = KS_eigenvalue(1,1,1)

  do i_k_points = 1, n_k_points,1
        do i_state = 1, n_states, 1
           if (KS_eigenvalue(i_state, 1,i_k_points) .lt. lowest_eigenvalue) then
              lowest_eigenvalue = KS_eigenvalue(i_state, 1, i_k_points)
           end if
        end do
  end do

  highest_eigenvalue = KS_eigenvalue(n_states,1,1)
  do i_k_points = 1,n_k_points,1
        do i_state = 1, n_states, 1
           if (KS_eigenvalue(i_state, 1, i_k_points) .gt. highest_eigenvalue) then
              highest_eigenvalue = KS_eigenvalue(i_state, 1, i_k_points)
           end if
        end do
  end do




  if (lowest_eigenvalue.ne.highest_eigenvalue) then
     chemical_potential_r = highest_eigenvalue
  else
     chemical_potential_r = 0.d0
  end if
  chemical_potential_l = lowest_eigenvalue

  call check_norm_dmft_p2_cluster(inv_full_ovlp_matrix_sqrt,&
                                  full_ovlp_matrix_sqrt, &
                                  full_ovlp_matrix, & 
                                  free_cluster_ham, &
                                  hybrid_func,&
                                  loc_self_enrg,&
                                  hartree_pot_LDA,& 
                                  n_electrons, &
                                  diff_l,&
                                  chemical_potential_l, &
                                  i_counter,&
                                  embed_part_number_new,&
                                  k_summed_overlap_matr,&
                                  self_energy_freq) 

  call check_norm_dmft_p2_cluster(inv_full_ovlp_matrix_sqrt,&
                                  full_ovlp_matrix_sqrt, &
                                  full_ovlp_matrix, &
                                  free_cluster_ham, &
                                  hybrid_func, &
                                  loc_self_enrg, &
                                  hartree_pot_LDA,&
                                  n_electrons, &
                                  diff_r, &
                                  chemical_potential_r,&
                                  i_counter,&
                                  embed_part_number_new,&
                                  k_summed_overlap_matr,&
                                  self_energy_freq) 

  do while (diff_l * diff_r .gt. 0.d0)
     !     interval for chemical potential still not found
     !     Must extend intervals both ways: 
     !     * For zero electrons, need Fermi level below lowest state
     !     * For all electrons in one channel, need Fermi level above highest state
     chemical_potential_l = chemical_potential_l - &
                           0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
     chemical_potential_r = chemical_potential_r + &
                           0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
     diff_l = diff_r

  call check_norm_dmft_p2_cluster(inv_full_ovlp_matrix_sqrt,&
                                  full_ovlp_matrix_sqrt, &
                                  full_ovlp_matrix, &
                                  free_cluster_ham,&
                                  hybrid_func,&
                                  loc_self_enrg,&  
                                  hartree_pot_LDA,&
                                  n_electrons, &
                                  diff_l, &
                                  chemical_potential_l,&
                                  i_counter,&
                                  embed_part_number_new,&
                                  k_summed_overlap_matr,&
                                  self_energy_freq)   

  call check_norm_dmft_p2_cluster(inv_full_ovlp_matrix_sqrt, &
                                  full_ovlp_matrix_sqrt, &
                                  full_ovlp_matrix, & 
                                  free_cluster_ham,&
                                  hybrid_func,&
                                  loc_self_enrg,&
                                  hartree_pot_LDA,& 
                                  n_electrons, &
                                  diff_r, &
                                  chemical_potential_r,&
                                  i_counter,&
                                  embed_part_number_new,&
                                  k_summed_overlap_matr,&
                                  self_energy_freq) 


!stop
     !     If the right interval cannot be found, stop potentially
     !     infinite loop
     if (i_counter .gt. max_zeroin) then
        write (info_str,'(1X,A)') "* Determination of electronic &
                          &chemical potential (Fermi level):"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* Too many iterations needed to&
                          &find good interval around E_F."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
             "* According to this check, you may not be able to find one. "
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
             "Please look carefully at your settings."
        call localorb_info(info_str,use_unit,'(A)')
        stop
     end if
  end do
     !     write(use_unit,*)'chem_l, chem_r', chemical_potential_l, chemical_potential_r

  call zeroin_dmft_cluster(inv_full_ovlp_matrix_sqrt, &
                           full_ovlp_matrix_sqrt, &
                           full_ovlp_matrix, &
                           free_cluster_ham, &
                           hybrid_func,&
                           loc_self_enrg,&
                           hartree_pot_LDA,&
                           chemical_potential_l, &
                           chemical_potential_r, &
                           diff_l, &
                           diff_r, &
                           fermi_acc, &
                           1.d-10,&
                            max_zeroin, &
                           n_basis, &
                           n_electrons,&
                           new_chem_pot,&
                           diff_electrons, &
                           i_counter, &
                           n_k_points, &
                           n_basis, &
                           nomega,&
                           embed_part_number_new,&
                           k_summed_overlap_matr,&
                           self_energy_freq)


  !     call check_norm one more time, just to make sure we really have the 
  !     occupation numbers for the final fermi level [this is critical for H]
!  call check_norm( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter)
  call check_norm_dmft_p2_cluster(inv_full_ovlp_matrix_sqrt,&
                                  full_ovlp_matrix_sqrt,&
                                  full_ovlp_matrix, &
                                  free_cluster_ham, &
                                  hybrid_func,&
                                  loc_self_enrg,&
                                  hartree_pot_LDA,&
                                  n_electrons, &
                                  diff_electrons,&
                                  new_chem_pot,&
                                  i_counter,&
                                  embed_part_number_new,&
                                  k_summed_overlap_matr,&
                                  self_energy_freq)   

  !     and decrement i_counter which was falsely incremented ...
  i_counter = i_counter - 1
  !       If zeroin did not converge to the desired accuracy, issue
  !       warning and stop, for now; we can create nicer behavior after
  !       we gain more experience.
  if (i_counter .gt. max_zeroin) then
     write (info_str,'(1X,A)') "* Determination of electronic chemical &
                                  &potential (Fermi level):"
     call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') "* Too many iterations needed to find &
                                  &correct occupation of levels."
     call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') "* Check variables occupation_acc, &
                                  &max_zeroin before continuing your calculation."
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* Status of get_occ_numbers when aborted: "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* State #         Eigenvalue       Occ. number"
     call localorb_info(info_str,use_unit,'(A)')
     do i_spin = 1, n_spin, 1
        write (info_str,*) "* spin # ", i_spin
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,*) "----------"
        call localorb_info(info_str,use_unit,'(A)')
     end do
     stop
  end if

  !     finally, occupation thresholding - if needed
!  if (occupation_thr.gt.0.d0) then
!     call threshold_occ_numbers( n_electrons, occ_numbers, diff_electrons_thr)
!  end if

  if (t_out) then
     write (info_str, '(2X, A, E10.4)') "| Chemical potential &
                               &(Fermi level) in eV                 : ", &
          new_chem_pot * hartree
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
     write (info_str, '(2X, A, E10.4)') "| Error in electron count &
                           &due to remaining E_F inaccuracy: ", diff_electrons
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
     if (occupation_thr.gt.0.d0) then
        write (info_str, '(2X, A, E10.4)') "| Error in &
                    &electron count after thresholding : ", diff_electrons_thr
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
     end if
  end if

!  call check_n_states(n_states, n_spin, 1, occ_numbers, t_out)

  if (t_out) then
     write(info_str,*)
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if

end subroutine update_chemical_pot_cluster 

!-------------------------------------------------------------------------------
!******	
