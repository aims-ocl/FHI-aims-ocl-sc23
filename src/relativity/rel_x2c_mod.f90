

 module rel_x2c_mod
!
! This module saves variables and functions for BOTH four-component Dirac-Kohn-Sham (4c-DKS) method and exact two-component (X2C) method.
! -- Rundong Zhao, June. 2018 @ Durham, NC
!
  implicit none

  integer :: scf_iteration=0
  logical :: upw = .false.  ! whether update W matrix in SCF iterations

 !----------------------- For free atoms -----------------------
  integer, dimension(:), allocatable :: dim_mat_atom ! dim_mat_atom(n_species): the dim_matrix_rel for free atoms
  real*8, dimension(:,:), allocatable :: atom_occ_shell ! electron occupation of different shells of free atoms
  real*8, dimension(:,:), allocatable :: free_den_diff ! allocated in free_atoms.f90; for x2c, this is rho_4c - rho_2c of free atoms;
                                                       ! for q4c, this is rho_4c - rho_large, viz. rho_small
  real*8, dimension(:,:,:), allocatable :: free_den_diff_spl ! cubic spline of free_den_diff
  real*8, dimension(:,:,:), allocatable :: free_drho_dr_diff_spl ! cubic spline of 1st derivative of free_den_diff
 !----------------------- For free atoms -----------------------


 !-------------------------- For basis --------------------------
  integer, dimension(:),  allocatable :: rel_shell ! rel_shell(n_species), allocated in subroutine shrink_fixed_basis_phi_thresh
  ! Be careful with the definition of the rel_shell here. Take neon as an example, there are 4 shells in total:
  ! 1s, 2s, 2p_1/2, 2p_3/2
  integer, dimension(:,:), allocatable :: l_rel     ! l.  l_rel(rel_shell,n_species)
  integer, dimension(:,:), allocatable :: j_rel     ! j*2.  j_rel(rel_shell,n_species).  Note, it is j times 2 here!

  integer :: dim_matrix_rel   ! dim_matrix_rel here is HALF of the dimension of the H/S matrix that will be used in the diagonalization.
  ! The reason is that the relativistic H/S matrices are symmetrical, in the following shape:
  ! ( A   B  )
  ! (-B^* A^*)
  ! Thus, we only save matrices A and B. dim_matrix_rel is the dimension of A and B.
  ! Take neon as an example, the nonrel single zeta basis set is: 1s, 2s, 2p_y, 2p_z, 2p_x. Totally 5 basis.
  ! For rel cases, the basis set splits into: 1s_(1/2,1/2), 1s_(1/2,-1/2), 2s_(1/2,1/2), 2s_(1/2,-1/2), 2p_(1/2,1/2), 2p_(1/2,-1/2),
  ! 2p_(3/2,3/2), 2p_(3/2,1/2), 2p_(3/2,-1/2), 2p_(3/2,-3/2). Totally 10 basis. Here, dim_matrix_rel = 5.
  ! For periodic systems, the final symmetrical form is:
  ! [ A(k)   B(k) ]
  ! [ B^+(k) C(k) ]
  integer, dimension(:,:), allocatable :: large_comp_idx
  integer, dimension(:,:), allocatable :: small_comp_idx
  real*8, dimension(:,:), allocatable  :: large_comp_clb ! large component CG coefficients
  real*8, dimension(:,:), allocatable  :: small_comp_clb
 !-------------------------- For basis --------------------------
 

 !--------------- For scalar basis used in scalar integrations---------------
  integer :: n_basis_small  ! Total number of small scalar basis. See the definition of n_basis in module dimensions for example.
                            ! For relativistic cases, the n_basis denotes the total number of large scalar basis.
  integer :: n_centers_basis_I_small ! Total number of small scalar basis (that can touch the 0th cell) in SUPERCELL. 
                                     ! See the definition of n_centers_basis_I in module dimensions for example.
 !--------------- For scalar basis used in scalar integrations---------------


 !------------------- For integrations --------------------
  ! As for dirac_v,t,w: the dimension is: (dim_matrix_rel,dim_matrix_rel,3,n_k_points). Because:
  ! 1-3 saves the matrices A(k), B(k), and C(k). For the definition of A,B,C, see the comments of dim_matrix_rel above.
  complex*16, dimension(:,:,:,:), allocatable :: dirac_v ! See varialbe "hamiltonian" and "ham_ovlp_work" in subroutine initialize_scf for comparison
  complex*16, dimension(:,:,:,:), allocatable :: dirac_t ! <chi_large|c \dot \sigma \dot p|chi_small>
  complex*16, dimension(:,:,:,:), allocatable :: dirac_w ! Ditto
  complex*16, dimension(:,:,:,:), allocatable :: dirac_ss   ! small component overlap matrix
  complex*16, dimension(:,:,:,:), allocatable :: dirac_ssc2 ! <chi_small| 2 m c^2 |chi_small>
  complex*16, dimension(:,:,:,:), allocatable :: dirac_s ! Overlap matrix. (dim_matrix_rel,dim_matrix_rel,3,n_k_points)

  real*8, dimension(:), allocatable :: dirac_v_sum ! corresponding to *hamiltonian* of none relativistic cases (see integrate_hamiltonian_matrix_p2)
  real*8, dimension(:), allocatable :: dirac_t_sum ! corresponding to *hamiltonian* of none relativistic cases (see integrate_hamiltonian_matrix_p2)
  real*8, dimension(:), allocatable :: dirac_w_sum
  real*8, dimension(:), allocatable :: dirac_ss_sum
  real*8, dimension(:), allocatable :: dirac_ssc2_sum
  real*8, dimension(:), allocatable :: dirac_s_sum

  real*8, dimension(:,:,:), allocatable :: x2c_v ! x2c_v(2n*(2n+1)/2,n_k_points,2), 1 and 2 save the real and imag part, respectively. n=dim_matrix_rel
  real*8, dimension(:,:,:), allocatable :: x2c_s
  real*8, dimension(:,:,:), allocatable :: x2c_v_old
  real*8, dimension(:,:,:), allocatable :: x2c_t
  real*8, dimension(:,:,:), allocatable :: x2c_w
  real*8, dimension(:,:,:), allocatable :: x2c_x  ! x2c_x(2n*2n,n_k_points,2)
  real*8, dimension(:,:,:), allocatable :: x2c_fw ! x2c_fw(2n*2n,n_k_points,2)
 !------------------- For integrations --------------------


  contains


  subroutine initialize_basis_index_x2c()
   use dimensions,   only : n_species, use_ionic, use_hydro, use_ext_basis
   use species_data, only : ionic_in_large_basis, hydro_in_large_basis, include_min_basis, &
                            l_shell_max, no_basis, species_pseudoized, n_atomic, atomic_l, &
                            n_conf, conf_l, n_conf, n_ionic, ionic_l, n_hydro, hydro_l, &
                            n_sto, sto_l, n_gaussian, gaussian_l, atoms_in_structure
   use aims_memory_tracking, only: aims_allocate, aims_deallocate
   use localorb_io
   implicit none
   logical need_section, spin_down
   integer i_species, i_l, i_l_last, i_atomic, i_conf, i_ionic, i_hydro, i_sto, i_gaussian, max_shell, &
           i_basis, i,j
   character*200 info_str

   if (.not.allocated(rel_shell)) call aims_allocate( rel_shell, n_species, "rel_shell" )

   rel_shell = 0; max_shell = 0
   do i_species = 1, n_species, 1
   if(.not.species_pseudoized(i_species).and.(.not.(no_basis(i_species)))) then
      do i_l = 0, l_shell_max(i_species), 1
 
         if (include_min_basis(i_species)) then
            ! treat atomic-like wave functions first
            do i_atomic = 1, n_atomic(i_species), 1
               if (atomic_l(i_species,i_atomic).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  endif
               end if
            enddo
         end if
 
         ! treat confined basis functions next
         if (n_conf(i_species).gt.0) then
            do i_conf = 1, n_conf(i_species), 1
               if (conf_l(i_species,i_conf).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  endif
               end if
            enddo
         end if
 
         ! treat ionic basis functions next
         if (use_ionic .and. n_ionic(i_species).gt.0) then
            ! JW: Avoid test for ionic_in_large_basis if .not. use_ionic.
            need_section = (.not. use_ext_basis) .or. any(.not. ionic_in_large_basis(i_species, 1:n_ionic(i_species)))
         else
            need_section = .false.
         end if
 
         if (need_section) then
            do i_ionic = 1, n_ionic(i_species), 1
               if (ionic_l(i_species, i_ionic).eq.i_l .and. (.not. ionic_in_large_basis(i_species,i_ionic)) ) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  endif
               end if
            enddo
         end if
         
         ! treat hydrogenic basis functions next
         ! check also that does all hydrogenic function belong to the large or small basis
         if (use_hydro .and. n_hydro(i_species).gt.0) then 
            ! JW: Avoid test for hydro_in_large_basis if .not. use_hydro.
            need_section = (.not. use_ext_basis) .or. any(.not. hydro_in_large_basis(i_species, 1:n_hydro(i_species)))
         else
            need_section = .false.
         end if
         if (need_section) then
            do i_hydro = 1, n_hydro(i_species), 1
               if (hydro_l(i_species, i_hydro).eq.i_l .and. .not. hydro_in_large_basis(i_species,i_hydro)) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  endif
               end if
            enddo
         end if
 
         ! treat STO basis functions next
         if (n_sto(i_species).gt.0) then
            do i_sto = 1, n_sto(i_species), 1
               if (sto_l(i_species,i_sto).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  endif
               end if
            enddo
         end if
 
         ! treat Gaussian basis functions next
         if (n_gaussian(i_species).gt.0) then
            do i_gaussian = 1, n_gaussian(i_species), 1
               if (gaussian_l(i_species,i_gaussian).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                  endif
               end if
            enddo
         end if
      enddo      ! end loop over l-shells

      if (rel_shell(i_species).gt.max_shell) max_shell = rel_shell(i_species)

   end if   ! .not.species_pseudoized
   enddo    ! first loop over species 

   if (.not.allocated(l_rel)) call aims_allocate( l_rel, max_shell, n_species, "l_rel" )
   if (.not.allocated(j_rel)) call aims_allocate( j_rel, max_shell, n_species, "j_rel" )

   ! Now, do the loop again to obtain the values of l_rel, j_rel, and dim_matrix_rel

   rel_shell = 0
   l_rel = 0
   j_rel = 0
   dim_matrix_rel = 0
   spin_down = .true.
   do i_species = 1, n_species, 1
   if(.not.species_pseudoized(i_species).and.(.not.(no_basis(i_species)))) then
      i_basis=0
      do i_l = 0, l_shell_max(i_species), 1
 
         if (include_min_basis(i_species)) then
            ! treat atomic-like wave functions first
            do i_atomic = 1, n_atomic(i_species), 1
               if (atomic_l(i_species,i_atomic).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l + 1
                     i_basis = i_basis + i_l*2 + 2
                  elseif (spin_down)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 - 1
                     i_basis = i_basis + i_l*2 - 1 + 1
                     spin_down = .false.  ! go to the lower "else" next time
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 + 1
                     i_basis = i_basis + i_l*2 + 2
                     spin_down = .true.  ! go to the upper "elseif" next time
                  endif
               end if
            enddo
         end if
 
         ! treat confined basis functions next
         if (n_conf(i_species).gt.0) then
            do i_conf = 1, n_conf(i_species), 1
               if (conf_l(i_species,i_conf).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l + 1
                     i_basis = i_basis + i_l*2 + 2
                  elseif (spin_down)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 - 1
                     i_basis = i_basis + i_l*2 - 1 + 1
                     spin_down = .false.  ! go to the lower "else" next time
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 + 1
                     i_basis = i_basis + i_l*2 + 2
                     spin_down = .true.  ! go to the upper "elseif" next time
                  endif
               end if
            enddo
         end if
 
         ! treat ionic basis functions next
         if (use_ionic .and. n_ionic(i_species).gt.0) then
            ! JW: Avoid test for ionic_in_large_basis if .not. use_ionic.
            need_section = (.not. use_ext_basis) .or. any(.not. ionic_in_large_basis(i_species, 1:n_ionic(i_species)))
         else
            need_section = .false.
         end if
 
         if (need_section) then
            do i_ionic = 1, n_ionic(i_species), 1
               if (ionic_l(i_species, i_ionic).eq.i_l .and. (.not. ionic_in_large_basis(i_species,i_ionic)) ) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l + 1
                     i_basis = i_basis + i_l*2 + 2
                  elseif (spin_down)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 - 1
                     i_basis = i_basis + i_l*2 - 1 + 1
                     spin_down = .false.  ! go to the lower "else" next time
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 + 1
                     i_basis = i_basis + i_l*2 + 2
                     spin_down = .true.  ! go to the upper "elseif" next time
                  endif
               end if
            enddo
         end if
         
         ! treat hydrogenic basis functions next
         ! check also that does all hydrogenic function belong to the large or small basis
         if (use_hydro .and. n_hydro(i_species).gt.0) then 
            ! JW: Avoid test for hydro_in_large_basis if .not. use_hydro.
            need_section = (.not. use_ext_basis) .or. any(.not. hydro_in_large_basis(i_species, 1:n_hydro(i_species)))
         else
            need_section = .false.
         end if
         if (need_section) then
            do i_hydro = 1, n_hydro(i_species), 1
               if (hydro_l(i_species, i_hydro).eq.i_l .and. .not. hydro_in_large_basis(i_species,i_hydro)) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l + 1
                     i_basis = i_basis + i_l*2 + 2
                  elseif (spin_down)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 - 1
                     i_basis = i_basis + i_l*2 - 1 + 1
                     spin_down = .false.  ! go to the lower "else" next time
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 + 1
                     i_basis = i_basis + i_l*2 + 2
                     spin_down = .true.  ! go to the upper "elseif" next time
                  endif
               end if
            enddo
         end if

         ! treat STO basis functions next
         if (n_sto(i_species).gt.0) then
            do i_sto = 1, n_sto(i_species), 1
               if (sto_l(i_species,i_sto).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l + 1
                     i_basis = i_basis + i_l*2 + 2
                  elseif (spin_down)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 - 1
                     i_basis = i_basis + i_l*2 - 1 + 1
                     spin_down = .false.  ! go to the lower "else" next time
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 + 1
                     i_basis = i_basis + i_l*2 + 2
                     spin_down = .true.  ! go to the upper "elseif" next time
                  endif
               end if
            enddo
         end if
 
         ! treat Gaussian basis functions next
         if (n_gaussian(i_species).gt.0) then
            do i_gaussian = 1, n_gaussian(i_species), 1
               if (gaussian_l(i_species,i_gaussian).eq.i_l) then
                  if (i_l.eq.0)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l + 1
                     i_basis = i_basis + i_l*2 + 2
                  elseif (spin_down)then
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 - 1
                     i_basis = i_basis + i_l*2 - 1 + 1
                     spin_down = .false.  ! go to the lower "else" next time
                  else
                     rel_shell(i_species) = rel_shell(i_species) + 1
                     l_rel(rel_shell(i_species),i_species) = i_l
                     j_rel(rel_shell(i_species),i_species) = i_l*2 + 1
                     i_basis = i_basis + i_l*2 + 2
                     spin_down = .true.  ! go to the upper "elseif" next time
                  endif
               end if
            enddo
         end if
      enddo      ! end loop over l-shells

      if (rel_shell(i_species).gt.max_shell) max_shell = rel_shell(i_species)

      dim_matrix_rel = dim_matrix_rel + i_basis*atoms_in_structure(i_species)

   end if   ! .not.species_pseudoized
   enddo    ! first loop over species 


  ! CG coef and link index between NR AO and 2c&4c AO basis
   dim_matrix_rel = dim_matrix_rel/2
   write(info_str,'(2X,A,I10)') "dim_matrix_rel=:", dim_matrix_rel
   call localorb_info(info_str,use_unit,'(A)',OL_norm)
   if (.not.allocated(large_comp_idx)) call aims_allocate( large_comp_idx, 4, dim_matrix_rel, "large_comp_idx" )
   if (.not.allocated(small_comp_idx)) call aims_allocate( small_comp_idx, 4, dim_matrix_rel, "small_comp_idx" )
   if (.not.allocated(large_comp_clb)) call aims_allocate( large_comp_clb, 4, dim_matrix_rel, "large_comp_clb" )
   if (.not.allocated(small_comp_clb)) call aims_allocate( small_comp_clb, 4, dim_matrix_rel, "small_comp_clb" )

 ! generate large_comp_idx, large_comp_clb, small_comp_idx, small_comp_clb:
   call integral_transform_nr2r_genindex(.true.)

  end subroutine initialize_basis_index_x2c


  subroutine integral_transform_nr2r_genindex(need_small_comp)
  ! generate transformation index from 1c picture to 2c or 4c picture
   use geometry,    only : species
   use dimensions,  only : n_atoms, n_species
   implicit none
   logical,intent(in)    :: need_small_comp
  ! output: large_comp_idx(4,*),small_comp_idx(4,*),large_comp_clb(4,*),small_comp_clb(4,*)
   integer :: i,j,m,k,kk,ql,qs,sig
   integer :: ll,jj,ls,mva,mva0,mm1,mms
   integer :: i_atom, i_species, i_shell

   kk=0; ql=0; qs=0
   do i_atom=1, n_atoms
     i_species = species(i_atom)
     do i_shell=1, rel_shell(i_species)
       ll=l_rel(i_shell,i_species); jj=j_rel(i_shell,i_species); ls=jj-ll  ! L; 2J; small comp. L==2j-l=l\pm 1
       mva0=1; mm1=(jj+1)/2                  ! only half |LJ-J>==>|LJJ>, others==>time reverse
       sig=-1; if(mod(mm1-ll,2).eq.1)sig=1
       do k=1,mm1
         mva = sign(mva0,sig); mms=(mva-1)/2
         mva0= mva0+2; sig=-sig
         kk  = kk+1
         call relativistic_get_clbidx_oneorb(.false.,ql,ll,jj,mms, mva,large_comp_idx(1,kk),large_comp_clb(1,kk))
         if(.not.need_small_comp)cycle
         call relativistic_get_clbidx_oneorb(.true., qs,ls,jj,mms, mva,small_comp_idx(1,kk),small_comp_clb(1,kk))
       enddo
       ql=ql+ll+ll+1
       qs=qs+ls+ls+1
     enddo       ! end of loop over all radial ao functions
   enddo         ! end of loop over atoms
  end subroutine integral_transform_nr2r_genindex



 end module rel_x2c_mod


 subroutine prepare_x2c()
   use dimensions, only : use_initial_rho, n_max_ind_fns, n_species
   use runtime_choices, only : packed_matrix_format, default_initial_moment,&
                               spin_treatment, real_eigenvectors, &
                               use_density_matrix,&
                               use_lapack_fast, use_elsi, flag_KS, &
                               use_scalapack
   use rel_x2c_mod, only : atom_occ_shell
   use aims_memory_tracking, only : aims_allocate, aims_deallocate
   implicit none

   use_initial_rho = .false.
   packed_matrix_format = 0  ! .false.
   default_initial_moment = 0
   spin_treatment = 0
   real_eigenvectors = .false.
   use_density_matrix = .true.

  ! This is equivalent to setting "KS_method serial" in the input file:
   use_lapack_fast = .true.
   use_elsi = .true.
   flag_KS = -1
   use_scalapack = .false.
  ! -------------------------------------------------------------------

   call aims_allocate( atom_occ_shell, n_max_ind_fns, n_species, 'atom_occ_shell' )

 end subroutine prepare_x2c


 subroutine rel_matrix_initialize()
  use rel_x2c_mod
  use dimensions, only: n_k_points, n_spin, n_centers_basis_I
  use runtime_choices, only: flag_rel, REL_x2c, REL_4c_dks
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  implicit none

  integer :: ld_large, ld_small

  ld_large = (n_centers_basis_I+1)*n_centers_basis_I/2
  ld_small = (n_centers_basis_I_small+1)*n_centers_basis_I_small/2

  if(.not.(allocated(dirac_v_sum)))call aims_allocate( dirac_v_sum, ld_large, "dirac_v_sum" )
 !call aims_allocate( dirac_w_sum,    ld_small, "dirac_w_sum" )
 !call aims_allocate( dirac_s_sum,    ld_large, "dirac_s_sum" )
 !call aims_allocate( dirac_ss_sum,   ld_small, "dirac_ss_sum" )
 !call aims_allocate( dirac_ssc2_sum, ld_small, "dirac_ssc2_sum" )
  if(flag_rel.eq.REL_x2c)&
  call aims_allocate( dirac_t_sum,    ld_large, "dirac_t_sum" )

  dirac_v_sum=0.d0 !; dirac_w_sum=0.d0; dirac_ss_sum=0.d0; dirac_s_sum=0.d0; dirac_ssc2_sum=0.d0
  if(flag_rel.eq.REL_x2c) dirac_t_sum=0.d0

  if( allocated(dirac_v) )   call aims_deallocate( dirac_v, "dirac_v" )
  if( allocated(dirac_t) )   call aims_deallocate( dirac_t, "dirac_t" )
  if( allocated(dirac_w) )   call aims_deallocate( dirac_w, "dirac_w")
  if( allocated(dirac_s) )   call aims_deallocate( dirac_s, "dirac_s" )
  if( allocated(dirac_ss) )  call aims_deallocate( dirac_ss, "dirac_ss" )
  if( allocated(dirac_ssc2) )call aims_deallocate( dirac_ssc2, "dirac_ssc2" )
  call aims_allocate( dirac_v,    dim_matrix_rel*dim_matrix_rel, 3, n_k_points, n_spin, "dirac_v" )
  call aims_allocate( dirac_w,    dim_matrix_rel*dim_matrix_rel, 3, n_k_points, n_spin, "dirac_w" )
  call aims_allocate( dirac_s,    dim_matrix_rel*dim_matrix_rel, 3, n_k_points, 1,      "dirac_s" )
  call aims_allocate( dirac_ss,   dim_matrix_rel*dim_matrix_rel, 3, n_k_points, n_spin, "dirac_ss")
  call aims_allocate( dirac_ssc2, dim_matrix_rel*dim_matrix_rel, 3, n_k_points, n_spin, "dirac_ssc2")
  if(flag_rel.eq.REL_x2c)&
  call aims_allocate( dirac_t,    dim_matrix_rel*dim_matrix_rel, 3, n_k_points, n_spin, "dirac_t" )

  dirac_v    = (0.d0,0.d0)
  dirac_w    = (0.d0,0.d0)
  dirac_s    = (0.d0,0.d0)
  dirac_ss   = (0.d0,0.d0)
  dirac_ssc2 = (0.d0,0.d0)
  if(flag_rel.eq.REL_x2c) dirac_t = (0.d0,0.d0)

 end subroutine rel_matrix_initialize

 subroutine cleanup_relativity()
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use rel_x2c_mod
  if( allocated(dim_mat_atom) )          call aims_deallocate( dim_mat_atom, "dim_mat_atom" )
  if( allocated(atom_occ_shell) )        call aims_deallocate( atom_occ_shell, "atom_occ_shell" )
 ! these arrays are deallocated in subroutine cleanup_free_atoms
 !if( allocated(free_den_diff) )         deallocate( free_den_diff )
 !if( allocated(free_den_diff_spl) )     deallocate( free_den_diff_spl )
 !if( allocated(free_drho_dr_diff_spl) ) deallocate( free_drho_dr_diff_spl )
  if( allocated(rel_shell) )             call aims_deallocate( rel_shell, "rel_shell" )
  if( allocated(l_rel) )                 call aims_deallocate( l_rel, "l_rel" )
  if( allocated(j_rel) )                 call aims_deallocate( j_rel, "j_rel" )
  if( allocated(large_comp_idx) )        call aims_deallocate( large_comp_idx, "large_comp_idx" )
  if( allocated(small_comp_idx) )        call aims_deallocate( small_comp_idx, "small_comp_idx" )
  if( allocated(large_comp_clb) )        call aims_deallocate( large_comp_clb, "large_comp_clb" )
  if( allocated(small_comp_clb) )        call aims_deallocate( small_comp_clb, "small_comp_clb" )
  if( allocated(dirac_v) )               call aims_deallocate( dirac_v, "dirac_v" )
  if( allocated(dirac_t) )               call aims_deallocate( dirac_t, "dirac_t" )
  if( allocated(dirac_w) )               call aims_deallocate( dirac_w, "dirac_w")
  if( allocated(dirac_ss) )              call aims_deallocate( dirac_ss, "dirac_ss" )
  if( allocated(dirac_ssc2) )            call aims_deallocate( dirac_ssc2, "dirac_ssc2" )
  if( allocated(dirac_s) )               call aims_deallocate( dirac_s, "dirac_s" )
  if( allocated(dirac_v_sum) )           call aims_deallocate( dirac_v_sum, "dirac_v_sum" )
  if( allocated(dirac_t_sum) )           call aims_deallocate( dirac_t_sum, "dirac_t_sum" )
  if( allocated(dirac_w_sum) )           call aims_deallocate( dirac_w_sum, "dirac_w_sum" )
  if( allocated(dirac_ss_sum) )          call aims_deallocate( dirac_ss_sum, "dirac_ss_sum" )
  if( allocated(dirac_ssc2_sum) )        call aims_deallocate( dirac_ssc2_sum, "dirac_ssc2_sum" )
  if( allocated(dirac_s_sum) )           call aims_deallocate( dirac_s_sum, "dirac_s_sum" )
  if( allocated(x2c_v) )                 call aims_deallocate( x2c_v, "x2c_v" )
  if( allocated(x2c_s) )                 call aims_deallocate( x2c_s, "x2c_s" )
  if( allocated(x2c_v_old) )             call aims_deallocate( x2c_v_old, "x2c_v_old" )
  if( allocated(x2c_t) )                 call aims_deallocate( x2c_t, "x2c_t" )
  if( allocated(x2c_w) )                 call aims_deallocate( x2c_w, "x2c_w" )
  if( allocated(x2c_x) )                 call aims_deallocate( x2c_x, "x2c_x" )
  if( allocated(x2c_fw) )                call aims_deallocate( x2c_fw, "x2c_fw" )
 end subroutine cleanup_relativity


