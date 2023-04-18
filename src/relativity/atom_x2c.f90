
! X2C for free atoms:
! Before starting the molecular/solid SCF, we should firstly solve the X2C for all the free atoms
! so as to generate the difference between 4C and 2C density -- the so-called picture change error
! (PCE) correction.
! -- Rundong Zhao, Nov. 2018 @ Durham, NC
!
 subroutine atom_x2c_pce()
  use constants, only: pi4, light_speed, light_speed_sq
  use dimensions, only: n_species
  use grids, only: n_grid, r_grid, log_r_grid_inc, r_grid_inc
  use species_data, only: n_atomic, atomic_n, atomic_l, l_shell_max
  use free_atoms, only: free_potential, free_wave, free_wave_small, free_wave_eigenval, &
                        free_wave_deriv, free_wave_small_deriv
  use rel_x2c_mod, only: dim_mat_atom, atom_occ_shell
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use localorb_io, only: use_unit
  implicit none

  integer :: i_species, i_atomic, i_basis, i_l, i_j, j_l, j_j, dim_mat, shell
  integer :: i,j,k,m,n
  logical :: spin_down
  real*8 :: sum_elec
  integer,allocatable :: n_atom(:), l_atom(:), j_atom(:), m_atom(:)
  real*8,allocatable :: occ_temp(:), occ(:), atom_e(:), atom_e_temp(:)
  real*8,allocatable :: weight(:)
  real*8,allocatable :: wave(:,:), wave_small(:,:), wave_kinetic(:,:), small_kinetic(:,:)
  real*8,allocatable :: wave_deriv(:,:), small_deriv(:,:)
  real*8,allocatable :: sa(:,:), ssa(:,:), ta(:,:), tsa(:,:), va(:,:), wa(:,:) ! A matrices of the STVW matices
  real*8,allocatable :: t1a(:,:), t2a(:,:)
  real*8,allocatable :: s_square(:,:), ss_square(:,:), t_square(:,:), ts_square(:,:), v_square(:,:), w_square(:,:) ! STVW matrices of the free atom, in square form
  real*8,allocatable :: s(:,:), ss(:,:), t(:,:), ts(:,:), v(:,:), w(:,:) ! STVW matrices of the free atom, in triangular form
  real*8,allocatable :: eigenvector(:,:), eigenval(:), free_eigenval(:), denmat(:,:)
  real*8,allocatable :: rho_2c(:)
  real*8,allocatable :: term1t1(:,:), term1t2(:,:), term2t1(:,:), term2t2(:,:)
  real*8 :: alpha, factor

  if (.not.allocated(dim_mat_atom)) call aims_allocate( dim_mat_atom, n_species, 'dim_mat_atom' )

  ! Obtain dim_mat_atom:
  dim_mat_atom = 0; spin_down = .true.
  do i_species = 1, n_species
     do i_l = 0, l_shell_max(i_species)
        do i_atomic = 1, n_atomic(i_species)
           if (atomic_l(i_species,i_atomic).eq.i_l) then
              if (i_l.eq.0)then
                 dim_mat_atom(i_species) = dim_mat_atom(i_species) + 1
              elseif (spin_down)then
                 dim_mat_atom(i_species) = dim_mat_atom(i_species) + i_l
                 spin_down = .false.  ! go to the lower "else" next time
              else
                 dim_mat_atom(i_species) = dim_mat_atom(i_species) + i_l + 1
                 spin_down = .true.  ! go to the upper "elseif" next time
              endif
           end if
        enddo ! i_atomic
     enddo ! i_l
  enddo ! i_species


  do i_species=1, n_species

     write(use_unit,"(' Do atomic X2C...  Species:',i4)") i_species
     alpha = log_r_grid_inc(i_species)
     write(use_unit,"(a,e20.12,5x,a,e20.12,5x,a,i6)") 'alpha:',alpha, 'r_grid_inc:',r_grid_inc(i_species), 'n_grid:',n_grid(i_species)
     write(use_unit,"(a,e20.12,5x,a,e20.12,5x)") 'r_grid_min:',r_grid(1,i_species), 'r_grid_max:',r_grid(n_grid(i_species),i_species)

     call aims_allocate( weight, n_grid(i_species), 'weight' )
     do n=1, n_grid(i_species)
        weight(n) = r_grid(n,i_species)*alpha
     enddo

     call aims_allocate( n_atom,     dim_mat_atom(i_species), 'n_atom' )
     call aims_allocate( l_atom,     dim_mat_atom(i_species), 'l_atom' )
     call aims_allocate( j_atom,     dim_mat_atom(i_species), 'j_atom' )
     call aims_allocate( m_atom,     dim_mat_atom(i_species), 'm_atom' )
     call aims_allocate( occ_temp,   dim_mat_atom(i_species), 'occ_temp' )
     call aims_allocate( occ,        2*dim_mat_atom(i_species), 'occ' )
     call aims_allocate( atom_e_temp,dim_mat_atom(i_species), 'atom_e_temp' )
     call aims_allocate( atom_e,     2*dim_mat_atom(i_species), 'atom_e' )

     call aims_allocate( wave,         n_grid(i_species), dim_mat_atom(i_species), 'wave' )
     call aims_allocate( wave_deriv,   n_grid(i_species), dim_mat_atom(i_species), 'wave_deriv' )
     call aims_allocate( small_deriv,  n_grid(i_species), dim_mat_atom(i_species), 'small_deriv' )
     call aims_allocate( wave_small,   n_grid(i_species), dim_mat_atom(i_species), 'wave_small' )
     call aims_allocate( wave_kinetic, n_grid(i_species), dim_mat_atom(i_species), 'wave_kinetic' )
     call aims_allocate( small_kinetic,n_grid(i_species), dim_mat_atom(i_species), 'small_kinetic' )

     dim_mat=0; shell=0; l_atom=-1; j_atom=-1; m_atom=0; occ_temp=0.d0; occ=0.d0
     spin_down = .true.
     do i_l = 0, l_shell_max(i_species)
        do i_atomic = 1, n_atomic(i_species)
           if (atomic_l(i_species,i_atomic).eq.i_l) then

              shell = shell + 1

              if (i_l.eq.0)then

                 dim_mat = dim_mat + 1

                 n_atom(dim_mat)      = atomic_n(i_species,i_atomic)
                 l_atom(dim_mat)      = i_l
                 j_atom(dim_mat)      = i_l + 1 ! j_atom denotes J*2
                 m_atom(dim_mat)      = -1
                 occ_temp(dim_mat)    = atom_occ_shell(shell,i_species)
                 atom_e_temp(dim_mat) = free_wave_eigenval(i_species,i_atomic)

                 do n=1, n_grid(i_species)
                   ! Note that the basis function here should contain a 1/r.
                   ! How ever, since the integration is spherical symmetrical, the final integrand is: psi(r) \dot operator \dot psi(r) r^2
                   ! the r^2 cancels the 1/r. Thus, here we do not need to generate a basis containing 1/r.
                    wave(n,         dim_mat) = free_wave(n,i_species,i_atomic)
                    wave_small(n,   dim_mat) = free_wave_small(n,i_species,i_atomic)
                    wave_kinetic(n, dim_mat) = (free_wave_eigenval(i_species,i_atomic) - free_potential(n,i_species)) * wave(n,dim_mat)
                    small_kinetic(n,dim_mat) = (free_wave_eigenval(i_species,i_atomic) - free_potential(n,i_species) + 2*light_speed_sq) &
                                               * wave_small(n,dim_mat)
                    wave_deriv(n,   dim_mat) = free_wave_deriv(n,i_species,i_atomic)
                    small_deriv(n,  dim_mat) = free_wave_small_deriv(n,i_species,i_atomic)
                 enddo

              elseif (spin_down)then ! dim_mat = dim_mat + i_l

                 do m= -i_l*2+1, 0, 2
                    dim_mat = dim_mat + 1
                    n_atom(dim_mat)      = atomic_n(i_species,i_atomic)
                    l_atom(dim_mat)      = i_l
                    j_atom(dim_mat)      = i_l*2 - 1 ! j_atom denotes J*2
                    m_atom(dim_mat)      = m
                    occ_temp(dim_mat)    = atom_occ_shell(shell,i_species) / i_l
                    atom_e_temp(dim_mat) = free_wave_eigenval(i_species,i_atomic)

                    do n=1, n_grid(i_species)
                       wave(n,         dim_mat) = free_wave(n,i_species,i_atomic)
                       wave_small(n,   dim_mat) = free_wave_small(n,i_species,i_atomic)
                       wave_kinetic(n, dim_mat) = (free_wave_eigenval(i_species,i_atomic) - free_potential(n,i_species)) * wave(n,dim_mat)
                       small_kinetic(n,dim_mat) = (free_wave_eigenval(i_species,i_atomic) - free_potential(n,i_species) + 2*light_speed_sq) &
                                                  * wave_small(n,dim_mat)
                       wave_deriv(n,   dim_mat) = free_wave_deriv(n,i_species,i_atomic)
                       small_deriv(n,  dim_mat) = free_wave_small_deriv(n,i_species,i_atomic)
                    enddo
                 enddo

                 spin_down = .false.  ! go to the lower "else" next time
              else ! dim_mat = dim_mat + i_l + 1

                 do m= -i_l*2-1, 0, 2
                    dim_mat = dim_mat + 1
                    n_atom(dim_mat)      = atomic_n(i_species,i_atomic)
                    l_atom(dim_mat)      = i_l
                    j_atom(dim_mat)      = i_l*2 + 1 ! j_atom denotes J*2
                    m_atom(dim_mat)      = m
                    occ_temp(dim_mat)    = atom_occ_shell(shell,i_species) / (i_l+1)
                    atom_e_temp(dim_mat) = free_wave_eigenval(i_species,i_atomic)

                    do n=1, n_grid(i_species)
                       wave(n,         dim_mat) = free_wave(n,i_species,i_atomic)
                       wave_small(n,   dim_mat) = free_wave_small(n,i_species,i_atomic)
                       wave_kinetic(n, dim_mat) = (free_wave_eigenval(i_species,i_atomic) - free_potential(n,i_species)) * wave(n,dim_mat)
                       small_kinetic(n,dim_mat) = (free_wave_eigenval(i_species,i_atomic) - free_potential(n,i_species) + 2*light_speed_sq) &
                                                  * wave_small(n,dim_mat)
                       wave_deriv(n,   dim_mat) = free_wave_deriv(n,i_species,i_atomic)
                       small_deriv(n,  dim_mat) = free_wave_small_deriv(n,i_species,i_atomic)
                    enddo
                 enddo

                 spin_down = .true.  ! go to the upper "elseif" next time
              endif
           end if
        enddo ! i_atomic
     enddo ! i_l

     m=0; n=0
     do i_basis=1, dim_mat_atom(i_species)
        n=m
        m=m+2
        occ(n+1:m) = occ_temp(i_basis)*0.5
        atom_e(n+1:m) = atom_e_temp(i_basis)
     enddo


     call aims_allocate( rho_2c, n_grid(i_species), 'rho_2c')

     call aims_allocate( ta, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_t_a')
     call aims_allocate(tsa, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_ts_a')
     call aims_allocate(t1a, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_t1_a')
     call aims_allocate(t2a, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_t2_a')
     allocate( term1t1 (dim_mat_atom(i_species), dim_mat_atom(i_species)) )
     allocate( term1t2 (dim_mat_atom(i_species), dim_mat_atom(i_species)) )
     allocate( term2t1 (dim_mat_atom(i_species), dim_mat_atom(i_species)) )
     allocate( term2t2 (dim_mat_atom(i_species), dim_mat_atom(i_species)) )
     term1t1=0.d0; term1t2=0.d0; term2t1=0.d0; term2t2=0.d0
     call aims_allocate( va, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_v_a')
     call aims_allocate( sa, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_s_a')
     call aims_allocate(ssa, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_ss_a')
     call aims_allocate( wa, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_w_a')
     ta=0.d0; tsa=0.d0; va=0.d0; sa=0.d0; ssa=0.d0; wa=0.d0; t1a=0.d0; t2a=0.d0
     call aims_allocate( t_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_t_square')
     call aims_allocate( v_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_v_square')
     call aims_allocate( s_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_s_square')
     call aims_allocate(ts_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_ts_square')
     call aims_allocate(ss_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_ss_square')
     call aims_allocate( w_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_w_square')
     t_square=0.d0; v_square=0.d0; s_square=0.d0; ss_square=0.d0; w_square=0.d0; ts_square=0.d0
     call aims_allocate( t, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_t')
     call aims_allocate( v, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_v')
     call aims_allocate( s, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_s')
     call aims_allocate( w, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_w')
     call aims_allocate(ts, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_ts')
     call aims_allocate(ss, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_ss')
     t=0.d0; v=0.d0; s=0.d0; ss=0.d0; w=0.d0; ts=0.d0

     call aims_allocate( eigenvector, 4*dim_mat_atom(i_species)**2, 2, 'eigenvector')
     call aims_allocate( eigenval, 2*dim_mat_atom(i_species), 'eigenval')
     call aims_allocate( denmat, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'denmat')
     eigenvector=0.d0; eigenval=0.d0; denmat=0.d0


    ! Integration: potential and overlap matrices.
     do i=1, dim_mat_atom(i_species)
        do j=1, i
           if( l_atom(j) .ne. l_atom(i) )then ! s-p, s-d, s-f, p-d, p-f, d-f
             cycle
           elseif( j_atom(j) .ne. j_atom(i) )then ! p_1/2 - p_3/2, d_3/2 - d_5/2, f_5/2 - f_7/2
             cycle
           elseif( m_atom(j) .ne. m_atom(i) )then ! e.g., p_(3/2,-3/2) - p(3/2,-1/2)
             cycle
           endif
           do n=1, n_grid(i_species)
              va(j,i)  = va(j,i)  + wave(n,j)*free_potential(n,i_species)*wave(n,i) * weight(n) 
              sa(j,i)  = sa(j,i)  + wave(n,j)*wave(n,i) * weight(n) 
              wa(j,i)  = wa(j,i)  + wave_small(n,j)*free_potential(n,i_species)*wave_small(n,i) * weight(n) 
              ssa(j,i) = ssa(j,i) + wave_small(n,j)*wave_small(n,i) * weight(n) 
           enddo
           va(i,j)=va(j,i)
           sa(i,j)=sa(j,i)
           wa(i,j)=wa(j,i)
           ssa(i,j)=ssa(j,i)
        enddo
     enddo

    ! For T integration, we should integrate a full matrix (instead of a triangular matrix) and then symmetrize it.
    ! This can avoid the numerical noise.
     do i=1, dim_mat_atom(i_species)
        do j=1, dim_mat_atom(i_species)
           if( l_atom(j) .ne. l_atom(i) )then ! s-p, s-d, s-f, p-d, p-f, d-f
             cycle
           elseif( j_atom(j) .ne. j_atom(i) )then ! p_1/2 - p_3/2, d_3/2 - d_5/2, f_5/2 - f_7/2
             cycle
           elseif( m_atom(j) .ne. m_atom(i) )then ! e.g., p_(3/2,-3/2) - p(3/2,-1/2)
             cycle
           endif
           factor   = l_atom(i)*(l_atom(i)+1) - dble(j_atom(i))*(dble(j_atom(i))/2+1)/2 - 0.25d0
          !write(use_unit,"('j=',i3,5x,'factor=',f12.6)")j_atom(i),factor
           do n=1, n_grid(i_species)
              ! write(use_unit,"('grid=',i5,5x,'r=',f12.6,5x,'wave_small=',f16.10,5x,'small_deriv=',f16.10)")n,r_grid(n,i_species),wave_small(n,i),small_deriv(n,i)
              ta(j,i)  = ta(j,i)  + wave(n,j)*wave_kinetic(n,i) * weight(n)
              tsa(j,i) = tsa(j,i) + wave_small(n,j)*small_kinetic(n,i) * weight(n)

              t1a(j,i) = t1a(j,i) + wave(n,j)      * light_speed * (-small_deriv(n,i) + wave_small(n,i)/r_grid(n,i_species) *factor) * weight(n)
              t2a(j,i) = t2a(j,i) + wave_small(n,j)* light_speed * ( wave_deriv(n,i)  + wave(n,i)/      r_grid(n,i_species) *factor) * weight(n)
              term1t1(j,i)= term1t1(j,i) - wave(n,j)*small_deriv(n,i)*light_speed*weight(n)
              term1t2(j,i)= term1t2(j,i) + wave_small(n,j)*wave_deriv(n,i)*light_speed*weight(n)
              term2t1(j,i)= term2t1(j,i) + wave(n,j)*wave_small(n,i)/r_grid(n,i_species)*factor*light_speed*weight(n)
              term2t2(j,i)= term2t2(j,i) + wave_small(n,j)*wave(n,i)/r_grid(n,i_species)*factor*light_speed*weight(n)
           enddo
        enddo
     enddo
     ! Symmetrize:
     do i=1, dim_mat_atom(i_species)
        do j=1, i
           ta(j,i)  = 0.5d0*(ta(j,i)  +  ta(i,j))
           ta(i,j)  = ta(j,i)
           tsa(j,i) = 0.5d0*(tsa(j,i) + tsa(i,j))
           tsa(i,j) = tsa(j,i)
           t1a(j,i) = 0.5d0*(t1a(j,i) + t1a(i,j))
           t1a(i,j) = t1a(j,i)
           t2a(j,i) = 0.5d0*(t2a(j,i) + t2a(i,j))
           t2a(i,j) = t2a(j,i)
        enddo
     enddo

     do i=1, dim_mat_atom(i_species)
        do j=1, dim_mat_atom(i_species)
           t_square(j,i) = ta(j,i)
           v_square(j,i) = va(j,i)
           s_square(j,i) = sa(j,i)
           w_square(j,i) = wa(j,i)
           ts_square(j,i)=tsa(j,i)
           ss_square(j,i)=ssa(j,i)
        enddo
     enddo
     do i=dim_mat_atom(i_species)+1, 2*dim_mat_atom(i_species)
        do j=dim_mat_atom(i_species)+1, 2*dim_mat_atom(i_species)
           s_square(j,i) = sa( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           t_square(j,i) = ta( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           v_square(j,i) = va( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           w_square(j,i) = wa( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           ts_square(j,i) = tsa( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           ss_square(j,i) = ssa( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
        enddo
     enddo

     m = 0
     do i=1, 2*dim_mat_atom(i_species)
     do j=1, i
       m=m+1
       s(m,1)=s_square(j,i)
       t(m,1)=t_square(j,i)
       v(m,1)=v_square(j,i)
       w(m,1)=w_square(j,i)
       ts(m,1)=ts_square(j,i)
       ss(m,1)=ss_square(j,i)
     enddo
     enddo

                write(6,*)'sa:'
                do i=1, dim_mat_atom(i_species)
                   write(use_unit,"(20(f11.5))")sa(1:dim_mat_atom(i_species),i)
                enddo

                write(6,*)'va:'
                do i=1, dim_mat_atom(i_species)
                   write(use_unit,"(20(f11.5))")va(1:dim_mat_atom(i_species),i)
                enddo

                write(6,*)'ta:'
                do i=1, dim_mat_atom(i_species)
                   write(use_unit,"(20(f11.5))")ta(1:dim_mat_atom(i_species),i)
                enddo

                write(6,*)'tsa:'
                do i=1, dim_mat_atom(i_species)
                   write(use_unit,"(20(f11.5))")tsa(1:dim_mat_atom(i_species),i)
                enddo

                write(6,*)'wa:'
                do i=1, dim_mat_atom(i_species)
                   write(use_unit,"(20(f11.5))")wa(1:dim_mat_atom(i_species),i)
                enddo

                write(6,*)'ssa:'
                do i=1, dim_mat_atom(i_species)
                   write(use_unit,"(20(f11.5))")ssa(1:dim_mat_atom(i_species),i)
                enddo

   ! call atom_test_ss(dim_mat_atom(i_species),s_square,t_square,v_square,w_square,ss_square) !!!!!

     call atom_x2c_routine(dim_mat_atom(i_species),s,ss,t,ts,v,w)

     call atom_diagonalization(dim_mat_atom(i_species),s,v,eigenvector,eigenval)

     call atom_occ(dim_mat_atom(i_species),occ,atom_e)

     call atom_denmat(dim_mat_atom(i_species),occ,eigenvector,denmat)

     call atom_density_generating(n_grid(i_species), dim_mat_atom(i_species), denmat, wave, r_grid(1,i_species), rho_2c)

     sum_elec=0.d0
     do k=1, n_grid(i_species)
       sum_elec = sum_elec + weight(k) * rho_2c(k) * r_grid(k,i_species)**2
     enddo
     write(use_unit,*)'sum of two-component elec:',sum_elec

     call atom_density_pce(n_grid(i_species), i_species, rho_2c)

     call aims_deallocate( weight, 'weight')

     call aims_deallocate( n_atom,     'n_atom')
     call aims_deallocate( l_atom,     'l_atom')
     call aims_deallocate( j_atom,     'j_atom')
     call aims_deallocate( m_atom,     'm_atom')
     call aims_deallocate( occ_temp,   'occ_temp')
     call aims_deallocate( occ,        'occ')
     call aims_deallocate( atom_e_temp,'atom_e_temp')
     call aims_deallocate( atom_e,     'atom_e')

     call aims_deallocate( wave,         'wave')
     call aims_deallocate( wave_small,   'wave_small')
     call aims_deallocate( wave_kinetic, 'wave_kinetic')
     call aims_deallocate( small_kinetic,'small_kinetic')
     call aims_deallocate( wave_deriv,   'wave_deriv')
     call aims_deallocate( small_deriv,  'small_deriv')

     call aims_deallocate( rho_2c, 'rho_2c')

     call aims_deallocate( ta, 'atom_t_a')
     call aims_deallocate(tsa, 'atom_ts_a')
     call aims_deallocate(t1a, 'atom_t1_a')
     call aims_deallocate(t2a, 'atom_t2_a')
     deallocate (term1t1,term1t2,term2t1,term2t2)
     call aims_deallocate( va, 'atom_v_a')
     call aims_deallocate( sa, 'atom_s_a')
     call aims_deallocate(ssa, 'atom_ss_a')
     call aims_deallocate( wa, 'atom_w_a')

     call aims_deallocate( t_square, 'atom_t_square')
     call aims_deallocate( v_square, 'atom_v_square')
     call aims_deallocate( s_square, 'atom_s_square')
     call aims_deallocate( w_square, 'atom_w_square')
     call aims_deallocate(ss_square, 'atom_ss_square')
     call aims_deallocate(ts_square, 'atom_ts_square')

     call aims_deallocate( t, 'atom_t')
     call aims_deallocate( v, 'atom_v')
     call aims_deallocate( s, 'atom_s')
     call aims_deallocate( w, 'atom_w')
     call aims_deallocate(ts, 'atom_ts')
     call aims_deallocate(ss, 'atom_ss')

     call aims_deallocate( eigenvector, 'eigenvector')
     call aims_deallocate( eigenval, 'eigenval')
     call aims_deallocate( denmat, 'denmat')

  enddo ! i_species

 end subroutine atom_x2c_pce


 subroutine atom_x2c_pce_sto()
  use constants, only: pi4, light_speed
  use dimensions, only: n_species
  use grids, only: n_grid, r_grid, log_r_grid_inc
  use species_data, only: n_atomic, atomic_l, l_shell_max, n_sto, sto_n, sto_l, sto_k, &
                          sto_wave, sto_wave_small, sto_kinetic
  use free_atoms, only: free_potential
  use rel_x2c_mod, only: dim_mat_atom, atom_occ_shell
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use localorb_io, only: use_unit
  implicit none

  integer :: i_species, i_atomic, i_sto, i_basis, i_l, i_j, j_l, j_j, dim_mat, shell
  integer :: i,j,k,m,n
  logical :: spin_down
  real*8 :: sum_elec, sum_large, sum_small
  integer,allocatable :: n_atom(:), l_atom(:), j_atom(:), m_atom(:)
  real*8,allocatable :: occ_temp(:), occ(:)
  real*8,allocatable :: weight(:)
  real*8,allocatable :: wave(:,:), wave_small(:,:), wave_kinetic(:,:)

  real*8,allocatable :: sa(:,:), ta(:,:), va(:,:), wa(:,:) ! A matrices of the STVW matices
  real*8,allocatable :: s_square(:,:), t_square(:,:), v_square(:,:), w_square(:,:) ! STVW matrices of the free atom, in square form
  real*8,allocatable :: s(:,:), t(:,:), v(:,:), w(:,:) ! STVW matrices of the free atom, in triangular form

  real*8,allocatable :: clarge(:,:,:), csmall(:,:,:)  ! four-component large and small coefficients
  real*8,allocatable :: c2c(:,:) ! two-component coefficients after FW transformation
  real*8,allocatable :: eigenval(:), denmat(:,:), denmat_large(:,:), denmat_small(:,:)
  real*8,allocatable :: rho_large(:), rho_small(:), rho_2c(:)
  real*8 :: alpha

  if (.not.allocated(dim_mat_atom)) call aims_allocate( dim_mat_atom, n_species, 'dim_mat_atom' )

  ! Obtain dim_mat_atom:
  dim_mat_atom = 0; spin_down = .true.
  do i_species = 1, n_species
     do i_l = 0, l_shell_max(i_species)

        do i_sto = 1, n_sto(i_species)
           if (sto_l(i_species,i_sto).eq.i_l) then
              if (i_l.eq.0)then
                 dim_mat_atom(i_species) = dim_mat_atom(i_species) + 1
              elseif (spin_down)then
                 dim_mat_atom(i_species) = dim_mat_atom(i_species) + i_l
                 spin_down = .false.  ! go to the lower "else" next time
              else
                 dim_mat_atom(i_species) = dim_mat_atom(i_species) + i_l + 1
                 spin_down = .true.  ! go to the upper "elseif" next time
              endif
           end if
        enddo ! i_sto

     enddo ! i_l
  enddo ! i_species


  do i_species=1, n_species

     alpha = log_r_grid_inc(i_species)
     write(use_unit,"(a,5x,f12.8,5x,a,i6)")'alpha:',alpha,'number of grids:',n_grid(i_species)

     call aims_allocate( weight, n_grid(i_species), 'weight' )
     do n=1, n_grid(i_species)
        weight(n) = r_grid(n,i_species)*alpha
     enddo

     call aims_allocate( n_atom,     dim_mat_atom(i_species), 'n_atom' )
     call aims_allocate( l_atom,     dim_mat_atom(i_species), 'l_atom' )
     call aims_allocate( j_atom,     dim_mat_atom(i_species), 'j_atom' )
     call aims_allocate( m_atom,     dim_mat_atom(i_species), 'm_atom' )
     call aims_allocate( occ_temp,   dim_mat_atom(i_species), 'occ_temp' )
     call aims_allocate( occ,        2*dim_mat_atom(i_species), 'occ' )

     call aims_allocate( wave,         n_grid(i_species), dim_mat_atom(i_species), 'wave' )
     call aims_allocate( wave_small,   n_grid(i_species), dim_mat_atom(i_species), 'wave_small' )
     call aims_allocate( wave_kinetic, n_grid(i_species), dim_mat_atom(i_species), 'wave_kinetic' )

     ! the preliminary occupation can be obtained from the minimal NAO basis:
     dim_mat=0; shell=0; occ_temp=0.d0
     spin_down = .true.
     do i_l = 0, l_shell_max(i_species)
        do i_atomic = 1, n_atomic(i_species)
           if (atomic_l(i_species,i_atomic).eq.i_l) then
              shell = shell + 1
              if (i_l.eq.0)then
                 dim_mat = dim_mat + 1
                 occ_temp(dim_mat)    = atom_occ_shell(shell,i_species)
              elseif (spin_down)then ! dim_mat = dim_mat + i_l
                 do m= -i_l*2+1, 0, 2
                    dim_mat = dim_mat + 1
                    occ_temp(dim_mat)    = atom_occ_shell(shell,i_species) / i_l
                 enddo
                 spin_down = .false.  ! go to the lower "else" next time
              else ! dim_mat = dim_mat + i_l + 1
                 do m= -i_l*2-1, 0, 2
                    dim_mat = dim_mat + 1
                    occ_temp(dim_mat)    = atom_occ_shell(shell,i_species) / (i_l+1)
                 enddo
                 spin_down = .true.  ! go to the upper "elseif" next time
              endif
           end if
        enddo ! i_atomic
     enddo ! i_l

     ! get wavefunctions:
     dim_mat=0; l_atom=-1; j_atom=-1; m_atom=0
     spin_down = .true.
     do i_l = 0, l_shell_max(i_species)

        do i_sto = 1, n_sto(i_species)
           if (sto_l(i_species,i_sto).eq.i_l) then

              if (i_l.eq.0)then

                 dim_mat = dim_mat + 1

                 n_atom(dim_mat) = sto_n(i_species,i_sto)
                 l_atom(dim_mat) = i_l
                 j_atom(dim_mat) = i_l + 1 ! j_atom denotes J*2
                 m_atom(dim_mat) = -1

                 do n=1, n_grid(i_species)
                    wave(n,         dim_mat) = sto_wave(n,i_species,i_sto)
                    wave_small(n,   dim_mat) = sto_wave_small(n,i_species,i_sto)
                    wave_kinetic(n, dim_mat) = sto_kinetic(n,i_species,i_sto)
                 enddo

              elseif (spin_down)then ! dim_mat = dim_mat + i_l

                 do m= -i_l*2+1, 0, 2
                    dim_mat = dim_mat + 1
                    n_atom(dim_mat) = sto_n(i_species,i_sto)
                    l_atom(dim_mat) = i_l
                    j_atom(dim_mat) = i_l*2 - 1 ! j_atom denotes J*2
                    m_atom(dim_mat) = m

                    do n=1, n_grid(i_species)
                       wave(n,         dim_mat) = sto_wave(n,i_species,i_sto)
                       wave_small(n,   dim_mat) = sto_wave_small(n,i_species,i_sto)
                       wave_kinetic(n, dim_mat) = sto_kinetic(n,i_species,i_sto)
                    enddo
                 enddo

                 spin_down = .false.  ! go to the lower "else" next time
              else ! dim_mat = dim_mat + i_l + 1

                 do m= -i_l*2-1, 0, 2
                    dim_mat = dim_mat + 1
                    n_atom(dim_mat) = sto_n(i_species,i_sto)
                    l_atom(dim_mat) = i_l
                    j_atom(dim_mat) = i_l*2 + 1 ! j_atom denotes J*2
                    m_atom(dim_mat) = m

                    do n=1, n_grid(i_species)
                       wave(n,         dim_mat) = sto_wave(n,i_species,i_sto)
                       wave_small(n,   dim_mat) = sto_wave_small(n,i_species,i_sto)
                       wave_kinetic(n, dim_mat) = sto_kinetic(n,i_species,i_sto)
                    enddo
                 enddo

                 spin_down = .true.  ! go to the upper "elseif" next time
              endif
           end if
        enddo ! i_sto

     enddo ! i_l

     m=0; n=0; occ=0.d0
     do i_basis=1, dim_mat_atom(i_species)
        n=m
        m=m+2
        occ(n+1:m) = occ_temp(i_basis)*0.5 
     enddo
     call atom_occ_sto(dim_mat_atom(i_species),occ)


     call aims_allocate( rho_large, n_grid(i_species), 'rho_large')
     call aims_allocate( rho_small, n_grid(i_species), 'rho_small')
     call aims_allocate( rho_2c, n_grid(i_species), 'rho_2c')

     call aims_allocate( ta, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_t_a')
     call aims_allocate( va, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_v_a')
     call aims_allocate( sa, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_s_a')
     call aims_allocate( wa, dim_mat_atom(i_species), dim_mat_atom(i_species), 'atom_w_a')
     ta=0.d0; va=0.d0; sa=0.d0; wa=0.d0
     call aims_allocate( t_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_t_square')
     call aims_allocate( v_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_v_square')
     call aims_allocate( s_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_s_square')
     call aims_allocate( w_square, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 'atom_w_square')
     t_square=0.d0; v_square=0.d0; s_square=0.d0; w_square=0.d0
     call aims_allocate( t, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_t')
     call aims_allocate( v, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_v')
     call aims_allocate( s, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_s')
     call aims_allocate( w, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'atom_w')
     t=0.d0; v=0.d0; s=0.d0; w=0.d0

     call aims_allocate( clarge, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 2, 'clarge')
     call aims_allocate( csmall, 2*dim_mat_atom(i_species), 2*dim_mat_atom(i_species), 2, 'csmall')
     call aims_allocate( c2c, 4*dim_mat_atom(i_species)**2, 2, 'c2c')
     call aims_allocate( eigenval, 2*dim_mat_atom(i_species), 'eigenval')
     call aims_allocate( denmat, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'denmat')
     call aims_allocate( denmat_large, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'denmat_large')
     call aims_allocate( denmat_small, dim_mat_atom(i_species)*(2*dim_mat_atom(i_species)+1), 2, 'denmat_small')
     c2c=0.d0; eigenval=0.d0; denmat=0.d0; denmat_large=0.d0; denmat_small=0.d0


     do i=1, dim_mat_atom(i_species)
        do j=1, i
           if( l_atom(j) .ne. l_atom(i) )then ! s-p, s-d, s-f, p-d, p-f, d-f
             cycle
           elseif( j_atom(j) .ne. j_atom(i) )then ! p_1/2 - p_3/2, d_3/2 - d_5/2, f_5/2 - f_7/2
             cycle
           elseif( m_atom(j) .ne. m_atom(i) )then ! e.g., p_(3/2,-3/2) - p(3/2,-1/2)
             cycle
           endif
           do n=1, n_grid(i_species)
              sa(j,i) = sa(j,i) + wave(n,j)*wave(n,i) * weight(n)
              ta(j,i) = ta(j,i) + wave(n,j)*wave_kinetic(n,i) * weight(n)
              va(j,i) = va(j,i) + wave(n,j)*free_potential(n,i_species)*wave(n,i) * weight(n)
              wa(j,i) = wa(j,i) + wave_small(n,j)*free_potential(n,i_species)*wave_small(n,i) * weight(n)
           enddo
           sa(i,j)=sa(j,i)
           ta(i,j)=ta(j,i)
           va(i,j)=va(j,i)
           wa(i,j)=wa(j,i)
        enddo
     enddo

     do i=1, dim_mat_atom(i_species)
        do j=1, dim_mat_atom(i_species)
           s_square(j,i) = sa(j,i)
           t_square(j,i) = ta(j,i)
           v_square(j,i) = va(j,i)
           w_square(j,i) = wa(j,i)
        enddo
     enddo
     do i=dim_mat_atom(i_species)+1, 2*dim_mat_atom(i_species)
        do j=dim_mat_atom(i_species)+1, 2*dim_mat_atom(i_species)
           s_square(j,i) = sa( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           t_square(j,i) = ta( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           v_square(j,i) = va( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
           w_square(j,i) = wa( j-dim_mat_atom(i_species), i-dim_mat_atom(i_species) )
        enddo
     enddo

     m = 0
     do i=1, 2*dim_mat_atom(i_species)
     do j=1, i
       m=m+1
       s(m,1)=s_square(j,i)
       t(m,1)=t_square(j,i)
       v(m,1)=v_square(j,i)
       w(m,1)=w_square(j,i)
     enddo
     enddo

              ! write(use_unit,*)'sa:'
              ! do i=1, dim_mat_atom(i_species)
              !    write(6,"(20(f11.5))")sa(1:dim_mat_atom(i_species),i)
              ! enddo

                write(use_unit,*)'va:'
                do i=1, dim_mat_atom(i_species)
                   write(6,"(20(f11.5))")va(1:dim_mat_atom(i_species),i)
                enddo

              ! write(use_unit,*)'ta:'
              ! do i=1, dim_mat_atom(i_species)
              !    write(6,"(20(f11.5))")ta(1:dim_mat_atom(i_species),i)
              ! enddo

              ! write(use_unit,*)'wa:'
              ! do i=1, dim_mat_atom(i_species)
              !    write(6,"(20(f11.5))")wa(1:dim_mat_atom(i_species),i)
              ! enddo

     call atom_x2c_routine_sto(dim_mat_atom(i_species),s,s,t,t,v,w,clarge,csmall)

     call atom_diagonalization(dim_mat_atom(i_species),s,v,c2c,eigenval)

     call atom_denmat(dim_mat_atom(i_species),occ,clarge,denmat_large)
     call atom_denmat(dim_mat_atom(i_species),occ,csmall,denmat_small)
     call atom_denmat(dim_mat_atom(i_species),occ,c2c,denmat)

     call atom_density_generating(n_grid(i_species), dim_mat_atom(i_species), denmat_large, wave,       r_grid(1,i_species), rho_large)
     call atom_density_generating(n_grid(i_species), dim_mat_atom(i_species), denmat_small, wave_small, r_grid(1,i_species), rho_small)
     call atom_density_generating(n_grid(i_species), dim_mat_atom(i_species), denmat, wave, r_grid(1,i_species), rho_2c)

     sum_elec=0.d0
     do k=1, n_grid(i_species)
       sum_elec = sum_elec + weight(k) * rho_2c(k) * r_grid(k,i_species)**2
     enddo
     write(use_unit,*)'sum of two-component electron number:',sum_elec

     sum_large=0.d0
     do k=1, n_grid(i_species)
       sum_large = sum_large + weight(k) * rho_large(k) * r_grid(k,i_species)**2
     enddo
     write(use_unit,*)'sum of large component electron number:',sum_large

     sum_small=0.d0
     do k=1, n_grid(i_species)
       sum_small = sum_small + weight(k) * rho_small(k) * r_grid(k,i_species)**2
     enddo
     write(use_unit,*)'sum of small component electron number:',sum_small

     write(use_unit,*)'sum of large and small component electron number:',sum_large+sum_small

     call atom_density_pce_sto(n_grid(i_species), i_species, rho_large, rho_small, rho_2c)

     call aims_deallocate( weight,   'weight')

     call aims_deallocate( n_atom,     'n_atom')
     call aims_deallocate( l_atom,     'l_atom')
     call aims_deallocate( j_atom,     'j_atom')
     call aims_deallocate( m_atom,     'm_atom')
     call aims_deallocate( occ_temp,   'occ_temp')
     call aims_deallocate( occ,        'occ')

     call aims_deallocate( wave,         'wave')
     call aims_deallocate( wave_small,   'wave_small')
     call aims_deallocate( wave_kinetic, 'wave_kinetic')

     call aims_deallocate( rho_large, 'rho_large')
     call aims_deallocate( rho_small, 'rho_small')
     call aims_deallocate( rho_2c, 'rho_2c')

     call aims_deallocate( ta, 'atom_t_a')
     call aims_deallocate( va, 'atom_v_a')
     call aims_deallocate( sa, 'atom_s_a')
     call aims_deallocate( wa, 'atom_w_a')

     call aims_deallocate( t_square, 'atom_t_square')
     call aims_deallocate( v_square, 'atom_v_square')
     call aims_deallocate( s_square, 'atom_s_square')
     call aims_deallocate( w_square, 'atom_w_square')

     call aims_deallocate( t, 'atom_t')
     call aims_deallocate( v, 'atom_v')
     call aims_deallocate( s, 'atom_s')
     call aims_deallocate( w, 'atom_w')

     call aims_deallocate( clarge, 'clarge')
     call aims_deallocate( csmall, 'csmall')
     call aims_deallocate( c2c, 'c2c')
     call aims_deallocate( eigenval, 'eigenval')
     call aims_deallocate( denmat, 'denmat')
     call aims_deallocate( denmat_large, 'denmat_large')
     call aims_deallocate( denmat_small, 'denmat_small')

  enddo ! i_species

 end subroutine atom_x2c_pce_sto


 subroutine atom_x2c_routine_sto(dim_mat,s,ss,t,ts,v,w,clarge,csmall)
  use constants, only : light_speed
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use localorb_io, only: use_unit
  implicit none
  integer,intent(in) :: dim_mat
  real*8,intent(inout) :: s(dim_mat*(2*dim_mat+1),2), ss(dim_mat*(2*dim_mat+1),2), &
                          t(dim_mat*(2*dim_mat+1),2), ts(dim_mat*(2*dim_mat+1),2), &
                          v(dim_mat*(2*dim_mat+1),2),  w(dim_mat*(2*dim_mat+1),2)
  real*8,intent(out) :: clarge(2*dim_mat,2*dim_mat,2), csmall(2*dim_mat,2*dim_mat,2)

  integer :: iesc=1 ! =1 NESC; =2 SESC.
  integer :: xdim, jsyml(2), norb(2)
  integer :: i,j,k,l,m,n
  real*8,allocatable :: x(:,:), fw(:,:), u(:,:)
  real*8,allocatable :: cc(:,:,:),ee(:)

  xdim=dim_mat*dim_mat*4
  jsyml=1; norb=dim_mat*2

  call aims_allocate( x, xdim, 2 ,"atom_x2c_x" )
  call aims_allocate( fw, xdim, 2 ,"atom_x2c_fw" )
  call aims_allocate( u, xdim, 2 ,"atom_sesc_u" )
  x=0.d0; fw=0.d0; u=0.d0

  allocate( cc(4*dim_mat,4*dim_mat,2), ee(4*dim_mat) )
 
 ! Solve full Dirac equation
  call x2c_calc_fdirac_solver(1,1,1, jsyml, norb, light_speed, &
       v,v(1,2), w,w(1,2), t,t(1,2), s,s(1,2), ss,ss(1,2), cc,cc(1,1,2), ee)
 !call x2c_fdirac_solver( dim_mat, v,v(1,2), w,w(1,2), t,t(1,2), s,s(1,2), ss,ss(1,2), cc, ee)
                      write(use_unit,*)'Four-component atomic eigenvalues:'
                      write(use_unit,"(10f18.8)")ee

 ! Get X matrix by equation: X = BA^{-1}
  call x2c_calc_xmat_generator(1,1,1, jsyml, norb, cc, cc(1,1,2), x, x(1,2) )

  do i=1, 2*dim_mat
    clarge(1:2*dim_mat,i,1) = cc(1:2*dim_mat,           i+2*dim_mat,1)
    clarge(1:2*dim_mat,i,2) = cc(1:2*dim_mat,           i+2*dim_mat,2)
    csmall(1:2*dim_mat,i,1) = cc(2*dim_mat+1:4*dim_mat, i+2*dim_mat,1)
    csmall(1:2*dim_mat,i,2) = cc(2*dim_mat+1:4*dim_mat, i+2*dim_mat,2)
  enddo

  deallocate( cc,ee )

 ! Generate FW transformation matrix and sesc tilde{S}S^{-1}; also NESC form kinetic matrix with new X matrix
  call x2c_gen_sescu_nesct_fw_mat( iesc, 1,1,1, jsyml, norb, light_speed, &
       s,s(1,2), t,t(1,2), x,x(1,2), fw,fw(1,2), t,t(1,2), ss, ss(1,2), u,u(1,2) )

  ! W matrix ==> XWX
  call x2c_w_to_xwx(1,1,1, jsyml,norb, x,x(1,2), w,w(1,2) )


  l = dim_mat*(2*dim_mat+1) * 2
  ! NESC:
  call apbtoc(l, v, 1, w, 1, v, 1)


 ! FW Transformation:

  ! FW Transformation for NESC T matrix:
  call x2c_w_to_xwx( 1,1,1, jsyml, norb, fw, fw(1,2), t, t(1,2) )
  ! V:
  call x2c_w_to_xwx( 1,1,1, jsyml, norb, fw, fw(1,2), v, v(1,2) )

  ! Sum V and T (both FW-Transformed) to generate NESC H matrix:
  call apbtoc(l, v, 1, t, 1, v, 1)

  call aims_deallocate( x,  "atom_x2c_x" )
  call aims_deallocate( fw, "atom_x2c_fw" )
  call aims_deallocate( u,  "atom_sesc_u" )

 end subroutine atom_x2c_routine_sto


 subroutine atom_x2c_routine(dim_mat,s,ss,t,ts,v,w)
  use constants, only : light_speed
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use localorb_io, only: use_unit
  implicit none
  integer,intent(in) :: dim_mat
  real*8,intent(inout) :: s(dim_mat*(2*dim_mat+1),2), ss(dim_mat*(2*dim_mat+1),2), &
                          t(dim_mat*(2*dim_mat+1),2), ts(dim_mat*(2*dim_mat+1),2), &
                          v(dim_mat*(2*dim_mat+1),2),  w(dim_mat*(2*dim_mat+1),2)

  integer :: iesc=1 ! =1 NESC; =2 SESC.
  integer :: xdim, jsyml(2), norb(2)
  integer :: i,j,k,l,m,n
  real*8,allocatable :: x(:,:), fw(:,:), u(:,:)
  real*8,allocatable :: cc(:,:),ee(:)

  xdim=dim_mat*dim_mat*4
  jsyml=1; norb=dim_mat*2

  call aims_allocate( x, xdim, 2 ,"atom_x2c_x" )
  call aims_allocate( fw, xdim, 2 ,"atom_x2c_fw" )
  call aims_allocate( u, xdim, 2 ,"atom_sesc_u" )
  x=0.d0; fw=0.d0; u=0.d0


  allocate( cc(16*dim_mat*dim_mat,2), ee(4*dim_mat) )
 
 ! Solve full Dirac equation
  call x2c_calc_fdirac_solver_nao(1,1,1, jsyml, norb, light_speed, &
       v,v(1,2), w,w(1,2), t,t(1,2), ts,ts(1,2), s,s(1,2), ss,ss(1,2), cc,cc(1,2), ee)
 !call x2c_fdirac_solver( dim_mat, v,v(1,2), w,w(1,2), t,t(1,2), s,s(1,2), ss,ss(1,2), cc, ee)
                      write(use_unit,*)'Four-component atomic eigenvalues:'
                      write(use_unit,"(10f18.8)")ee

 ! Get X matrix by equation: X = BA^{-1}
  call x2c_calc_xmat_generator(1, 1, 1, jsyml, norb, cc, cc(:,2), x, x(1,2) )

  deallocate( cc,ee )

 ! Generate FW transformation matrix and sesc tilde{S}S^{-1}; also NESC form kinetic matrix with new X matrix
  call x2c_gen_sescu_nesct_fw_mat_atom(iesc, 1,1,1, jsyml, norb, light_speed, &
       s,s(1,2), t,t(1,2), ts,ts(1,2), x,x(1,2), fw,fw(1,2), t,t(1,2), ss, ss(1,2), u,u(1,2) )

  ! W matrix ==> XWX
  call x2c_w_to_xwx(1,1,1, jsyml, norb, x, x(1,2), w, w(1,2) )


  l = dim_mat*(2*dim_mat+1) * 2
  ! NESC:
  if( iesc.eq.1 ) call apbtoc(l, v, 1, w, 1, v, 1)


 ! FW Transformation:

 ! FW Transformation for NESC T matrix:
  call x2c_w_to_xwx( 1,1,1, jsyml, norb, fw, fw(1,2), t, t(1,2) )
 ! V:
  call x2c_w_to_xwx( 1,1,1, jsyml, norb, fw, fw(1,2), v, v(1,2) )


  ! Sum V and T (both FW-Transformed) to generate NESC H matrix:
  call apbtoc(l, v, 1, t, 1, v, 1)

  call aims_deallocate( x,  "atom_x2c_x" )
  call aims_deallocate( fw, "atom_x2c_fw" )
  call aims_deallocate( u,  "atom_sesc_u" )

 end subroutine atom_x2c_routine


 subroutine atom_diagonalization(dim_mat,s,v,eigenvector,eigenval)
  use localorb_io, only: use_unit
  implicit none
  integer,intent(in) :: dim_mat
  real*8,intent(in) :: s(dim_mat*(2*dim_mat+1),2), v(dim_mat*(2*dim_mat+1),2)
  real*8,intent(out) :: eigenvector(2*dim_mat*2*dim_mat,2), eigenval(2*dim_mat)

  integer :: xdim
  real*8  :: eps
  real*8,   allocatable :: smihalf(:,:),e(:)
  integer :: k,l

  eps = 1.d-12
  allocate( smihalf(2*dim_mat*2*dim_mat,2), e(2*dim_mat) )

  call pbc_diagmat_SX(.false., 2*dim_mat,eps, s,s(1,2), e, smihalf,smihalf(1,2), xdim )
  write(use_unit,"(' Xdim after Lowdin orthorgonalization:',i5)")xdim
  ! solve FC = SCE
  call pbc_diagmat_F(.false., 2*dim_mat,xdim, v,v(1,2), smihalf,smihalf(1,2), eigenvector(:,1),eigenvector(:,2),eigenval)
                      write(6,*)'X2C atomic eigenvalues:'
                      write(6,"(10f18.8)")eigenval

 end subroutine atom_diagonalization


 subroutine atom_occ(n,occ,e)
  implicit none
  integer,intent(in) :: n
  real*8,intent(inout) :: occ(2*n)
  real*8,intent(inout) :: e(2*n)
  integer :: i,j,k,m,orb
  real*8 :: tmp,tmpe
 ! reorder e from lowest to highest
  do orb=1, 2*n
    do i=orb+1, 2*n
       if(e(i).le.e(orb))then
          tmpe = e(i)
          e(i) = e(orb)
          e(orb) = tmpe
          tmp = occ(i)
          occ(i) = occ(orb)
          occ(orb) = tmp
       endif
    enddo
  enddo
  write(6,*)'atom occupation:'
  write(6,"(20f9.5)")occ
 end subroutine atom_occ

 subroutine atom_occ_sto(n,occ)
  implicit none
  integer,intent(in) :: n
  real*8,intent(inout) :: occ(2*n)
  integer :: i,j,k,m,orb
  real*8 :: tmp
 ! Reorder occu from lowest orbital to highest orbital.
 ! Since I used average occupation, it is OK here to order the occupation number
 ! from the largest value (viz. 1.0) to the smallest value (fraction numbers or 1.0)
  do orb=1, 2*n
    do i=orb+1, 2*n
       if(occ(i).gt.occ(orb))then
          tmp = occ(i)
          occ(i) = occ(orb)
          occ(orb) = tmp
       endif
    enddo
  enddo
  write(6,*)'atom occupation:'
  write(6,"(20f9.5)")occ
 end subroutine atom_occ_sto

 subroutine atom_denmat(dim_mat,occ,eigenvector,denmat)
  implicit none
  integer,intent(in) :: dim_mat
  real*8,intent(in) :: occ(2*dim_mat)
  real*8,intent(in) :: eigenvector(2*dim_mat*2*dim_mat,2)
  real*8,intent(out) :: denmat(dim_mat*(2*dim_mat+1),2)

  integer :: kk,i,j,k,m,n,orb
  real*8  :: cicjr,cicji,cir,cii,cjr,cji  ! cir: ci real. cii: ci imag. 

 ! mo_coef:
 ! column: n , from 1 - 2*dim_mat
 ! row: u , from 1 - 2*dim_mat
  m=0
  do orb=1, 2*dim_mat
    if(dabs(occ(orb)).gt.1.0d-8)then
      kk=1
      do i=1, 2*dim_mat
         cir=eigenvector(m+i,1)
         cii=eigenvector(m+i,2)
         do j=1,i
            cjr=eigenvector(m+j,1)
            cji=eigenvector(m+j,2)
           ! (a-ib)*(c+id)=(ac+bd)+i(ad-bc)
            cicjr=cir*cjr+cii*cji
            cicji=cir*cji-cii*cjr
            denmat(kk,1)=denmat(kk,1)+cicjr*occ(orb)
            denmat(kk,2)=denmat(kk,2)+cicji*occ(orb)
            kk=kk+1
         enddo ! j=1, i
      enddo  ! i=1, 2*dim_mat
    endif
    m=m+2*dim_mat
  enddo

 end subroutine atom_denmat


 subroutine atom_density_generating(n_grid, dim_mat, denmat, wave, r_grid, rho_2c)
  implicit none
  integer,intent(in) :: n_grid, dim_mat
  real*8,intent(in) :: denmat(dim_mat*(2*dim_mat+1),2)
  real*8,intent(in) :: wave(n_grid, dim_mat), r_grid(n_grid)
  real*8,intent(out) :: rho_2c(n_grid)

  integer :: i,j,k,m,n
  real*8,allocatable :: wave_work(:,:), denmat_square(:,:), r_grid_inv(:)

  allocate( wave_work(n_grid,2*dim_mat), denmat_square(2*dim_mat,2*dim_mat), r_grid_inv(n_grid) )

  m=0
  do i=1, 2*dim_mat
    do j=1, i
       m=m+1
       denmat_square(j,i) = denmat(m,1)
       denmat_square(i,j) = denmat(m,1)
    enddo
  enddo

              ! write(6,*)'denmat_square:'
              ! do i=1, 2*dim_mat
              !    write(6,"(20(f12.6))")denmat_square(1:2*dim_mat,i)
              ! enddo

  wave_work( 1:n_grid, 1:dim_mat) = wave( 1:n_grid, 1:dim_mat )
  wave_work( 1:n_grid, dim_mat+1:2*dim_mat) = wave( 1:n_grid, 1:dim_mat )

  do n=1, n_grid
     r_grid_inv(n) = 1.d0/r_grid(n)**2
  enddo

  rho_2c=0.d0
  do i=1, 2*dim_mat
  do j=1, 2*dim_mat
    do n=1, n_grid
       rho_2c(n) = rho_2c(n) + wave_work(n,j)*denmat_square(j,i)*wave_work(n,i) * r_grid_inv(n)
    enddo
  enddo
  enddo

  deallocate( wave_work,denmat_square,r_grid_inv )

 end subroutine atom_density_generating


 subroutine atom_density_pce(n_grid, i_species, rho_2c)
  use spline, only: cubic_spline
  use free_atoms, only: free_rho
  use rel_x2c_mod, only: free_den_diff, free_den_diff_spl
  implicit none
  integer,intent(in) :: n_grid, i_species
  real*8,intent(in) :: rho_2c(n_grid)

  integer :: i,j,k,m,n

  do n=1, n_grid
     free_den_diff(n,i_species) = free_rho(n,i_species) - rho_2c(n)
  enddo
              ! write(6,*)'free_den_diff'
              ! write(6,"(20f12.6)")free_den_diff(:,i_species)

 ! spline rho
  call cubic_spline ( free_den_diff(1,i_species), n_grid, free_den_diff_spl(1,1,i_species) )

 end subroutine atom_density_pce


 subroutine atom_density_pce_sto(n_grid, i_species, rho_large, rho_small, rho_2c)
  use spline, only: cubic_spline
  use free_atoms, only: free_rho
  use rel_x2c_mod, only: free_den_diff, free_den_diff_spl
  implicit none
  integer,intent(in) :: n_grid, i_species
  real*8,intent(in) :: rho_large(n_grid), rho_small(n_grid), rho_2c(n_grid)

  integer :: i,j,k,m,n

  do n=1, n_grid
     free_den_diff(n,i_species) = rho_large(n) + rho_small(n) - rho_2c(n)
  enddo

 ! spline rho
  call cubic_spline ( free_den_diff(1,i_species), n_grid, free_den_diff_spl(1,1,i_species) )

 end subroutine atom_density_pce_sto



