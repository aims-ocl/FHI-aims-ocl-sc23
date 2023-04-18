Module CC_cl

!  PURPOSE
!
!  This file may at some point contain all the needed variable declarations
!  for the Coupled-Cluster (CC) calculations (hybrid algorithm)
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Tonghao Shen
!  SOURCE

  Use CC_cl_distribution
  Use CC_cl_mem
  Use CC_cl_MPI
  Use dimensions
  Use mpi_tasks
  Use runtime_choices
  Use physics
  Use localorb_io
!  Use synchronize_mpi

  Implicit None

! Parameters discribe CC equations

  Integer (kind = 8) :: CC_n_config
! Number of configurations in CC calculations

  Integer :: CC_n_s
! Number of singlet excitations in CC calculations

  Integer (kind = 8) :: CC_n_d
! Number of double-excitations in CC calculations

  Integer :: CC_n_elec
! Number of electrons (spin free)

  Integer :: CC_n_occ
! Number of occupied(valence) orbitals (spin free)

  Integer :: CC_n_vir
! Number of virtual orbitals (spin free)

  Integer :: CC_n_fc
! Number of frozen core electrons

  Integer :: CC_valence
  ! The starting index of valence state

  Integer :: CC_n_state
  ! Number of total orbitals included (=n_states if full claculation is taken
  !                                    =n_states - n_frozen if frozen core approximation is used)

  Integer :: CC_n_bas
! Numbers of auxiliary basis functions on each core

  Integer :: CC_sv_strategy
! Saving strategy
! 1 ---- Maximum memory space demand, all vectors needed are saved in memory (replicated) (default)
! 2 ---- Medium memory space demand, no replication of (ka|cd) is saved 
! 3 ---- Minimum memory space demand, only necessary vectors are saved

  Integer :: CC_i_scf
! Number of CC scf iteration

  Double precision :: E_HF, CC_E_corr
  Double precision :: CC_E_PT = 0.0D0
! E_HF : HF energy
! CC_E_corr : Correlation energy
! CC_E_PT : perturbative triples correction

! RLE parameters
  Integer :: CC_RLE_ndc, CC_RLE_m_bgn
! CC_RLE_ndc ---- the dimension of current R matrix in RLE algorithm
! CC_RLE_m_bgn ---- the start point of saved basis functions

! DIIS parameters
  Integer :: CC_DIIS_ndc, CC_DIIS_m_bgn
! CC_DIIS_ndc ---- the dimension of current B matrix in DIIS algorithm
! CC_DIIS_m_bgn ---- the start point of saved basis functions

  Double precision :: CC_Norm
! Normalization coefficient 

  Logical :: CC_flag_converge = .false.
! CCSD converge flag

  Logical :: CC_RLE_flag_singular
! the flag of sigularity of R matrix

  Logical :: CC_RLE_flag_sv_overflow
! the flag of sigularity of R matrix

  Integer :: CC_n_omp_calc = 0
! Number of omp threads used in PT calculation

  Integer :: CC_PT_start = -1
  Integer :: CC_res_pt = 0

  Integer , dimension(:) , allocatable :: CC_res_inf
Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_code(code_fg,f,g,i_step,codetype)

! This subroutine maps an index (f,g) to code_fg
! codetype : 1 ---- code_fg = (f - 1) * i_step + g
!            2 ---- code_fg = (g - 1) * g / 2 + f             (f <= g)
!            3 ---- code_fg = (g - 2) * (g - 1) / 2 + f       (f < g)

  Implicit None
  Integer , intent(out) :: code_fg
  Integer , intent(in) :: i_step
  Integer , intent(in) :: f,g
  Integer , intent(in) :: codetype

  if (codetype.eq.1) then
    code_fg = (f - 1) * i_step + g
  else if (codetype.eq.2) then
    code_fg = (g - 1) * g / 2 + f
  else if (codetype.eq.3) then
    code_fg = (g - 2) * (g - 1) / 2 + f
  end if

  End Subroutine CC_cl_code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_decode(code_fg,f,g,i_step,detype)

! This subroutine maps an index code_fg to (f,g)
! detype : 1 ---- code_fg = (f - 1) * i_step + g
!          2 ---- code_fg = (g - 1) * g / 2 + f
!          3 ---- code_fg = (g - 2) * (g - 1) / 2 + f

! Therefore we have,
! detype : 1 ---- g = Mod(code_fg,i_step) , f = Int(code_fg/i_step) + 1
!          2 ---- Note that f <= g, for c = 2 * code_fg = g**2 - g + 2 * f
!                 we have (g - 1)**2 < c < (g + 1)**2
!                 and thus int(sqrt(c)) = g or int(sqrt(c)) = g - 1
!                 if f > int(sqrt(c)), g = int(sqrt(c)) + 1
!                 else g = int(sqrt(c))
!          3 ---- Similar to 2, the only difference is that f < g, assume that h = g - 1
!                 we get the same expansion by substitute (g - 1) by h.

  Implicit None

  Integer , intent(in) :: code_fg
  Integer , intent(in) :: i_step
  Integer , intent(out) :: f,g
  Integer , intent(in) :: detype

  if (detype.eq.1) then

    g = Mod(code_fg,i_step)
    if ((g.eq.0).and.(code_fg.ne.0)) then
      g = i_step
    end if

    f = Int((code_fg - g) / i_step) + 1

  else if (detype.eq.2) then

    g = Int(sqrt(2.0D0 * dble(code_fg)))
    f = code_fg - g * (g - 1) / 2
    if (f.gt.g) then
      g = g + 1
      f = code_fg - g * (g - 1) / 2
    else if (f.eq.0) then
      g = g - 1
      f = code_fg - g * (g - 1) / 2
    end if

  else if (detype.eq.3) then

    g = Int(sqrt(2.0D0 * dble(code_fg)))
    f = code_fg - g * (g - 1) / 2
    if (f.gt.g) then
      g = g + 1
      f = code_fg - g * (g - 1) / 2
    else if (f.eq.0) then
      g = g - 1
      f = code_fg - g * (g - 1) / 2
    end if

    g = g + 1

  end if

  End Subroutine CC_cl_decode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_clean_w()

! This subroutine clean w vector

  Implicit None

  CC_w_s = 0.0D0
  CC_w_d = 0.0D0

  End Subroutine CC_cl_clean_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_ev_diff(iii,aaa,E_diff)

  ! This subroutine calculates the difference between rrr and sss eigenvalue
  ! E_diff = KS_eigenvalue(sss) - KS_eigenvalue(rrr)

  Implicit None
  Integer , intent(in) :: iii,aaa
  Double precision , intent(out) :: E_diff
  Integer :: rrr,sss

  sss = aaa + CC_n_elec
  rrr = iii + CC_n_fc
  E_diff = KS_eigenvalue(sss,1,1) - KS_eigenvalue(rrr,1,1)

  End Subroutine CC_cl_ev_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_w2t()

  Implicit None
  Integer (kind = 8) :: s_tmp
  Integer :: i_send,i_recv,req1,req2
  Integer :: i_task,i_task_run,i_id
  Integer :: iii,jjj,aaa,bbb,i_tmp
  Integer :: n_ij,ij_run,ij_start,i_ij
  Integer :: n_b,b_start,b_end,b_run,a_run,n_a,a_start,a_end
  Integer :: errnum
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2
  Double precision :: coeff
  Double precision , dimension(:,:,:,:) , allocatable :: w_recv,w_tmp,w_d_all

  ! For single excitation
  CC_t_s = CC_w_s
  CC_t_d = 0.0D0

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  ij_start = CC_index_ij_D(CC_mpi_did+1)

  n_b = CC_mem_aa_D(CC_mpi_did+1)
  b_start = CC_index_aa_D(CC_mpi_did+1)
  b_end = b_start - 1 + n_b

  Allocate(w_d_all(n_ij,CC_n_vir,CC_n_vir,1),stat=errnum)
  Call check_allocation(errnum,'w_d_all in CC1')

  w_d_all(:,a_start:a_end,:,1) = CC_w_d(:,:,:,1)

  i_send = CC_mpi_gid + 1
  i_recv = CC_mpi_gid + 1

  do i_task_run = 1, CC_mpi_group_size - 1

    i_recv = i_recv + 1
    if (i_recv.gt.CC_mpi_group_size) then
      i_recv = 1
    end if

    i_send = i_send - 1
    if (i_send.lt.1) then
      i_send = CC_mpi_group_size
    end if

    n_a = CC_mem_aa_G(i_recv)
    a_start = CC_index_aa_G(i_recv)
    a_end = a_start - 1 + n_a

    Allocate(w_recv(n_ij,CC_mem_aa_G(i_recv),CC_n_vir,1),stat=errnum)
    Call check_allocation(errnum,'w_recv in CC')

    s_tmp = Int(CC_mem_aa_G(CC_mpi_gid+1),8) * Int(n_ij*CC_n_vir,8)
    Call CC_mpi_real_isend(s_tmp,CC_w_d,i_send-1,100,req1,CC_mpi_comm_group)

    s_tmp = Int(CC_mem_aa_G(i_recv),8) * Int(n_ij*CC_n_vir,8)
    Call CC_mpi_real_irecv(s_tmp,w_recv,i_recv-1,100,req2,CC_mpi_comm_group)

    Call MPI_WAIT(req1,stat1,errnum)
    Call MPI_WAIT(req2,stat2,errnum)

    w_d_all(:,a_start:a_end,:,1) = w_recv(:,:,:,1)

    Deallocate(w_recv)

  end do

  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  ij_start = CC_index_ij_D(CC_mpi_did+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  Allocate(w_tmp(n_ij,CC_n_vir,CC_n_vir,1),stat=errnum)
  Call check_allocation(errnum,'w_tmp in CC1')

  w_tmp = w_d_all

  i_task = CC_mpi_did
  i_send = CC_mpi_did + 1
  i_recv = CC_mpi_did + 1

  do i_task_run = 1, CC_mpi_domain_size

    i_recv = i_recv + 1
    if (i_recv.gt.CC_mpi_domain_size) then
      i_recv = 1
    end if

    i_send = i_send - 1
    if (i_send.lt.1) then
      i_send = CC_mpi_domain_size
    end if

    i_task = i_task + 1
    if (i_task.gt.CC_mpi_domain_size) then
      i_task = 1
    end if

    if (i_recv.ne.CC_mpi_did+1) then

      Allocate(w_recv(CC_mem_ij_D(i_recv),CC_n_vir,CC_n_vir,1),stat=errnum)
      Call check_allocation(errnum,'w_recv in CC')

      i_tmp = CC_n_vir * CC_n_vir
      s_tmp = Int(i_tmp,8) * Int(CC_mem_ij_D(CC_mpi_did+1),8)

      Call CC_mpi_real_isend(s_tmp,w_d_all,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(CC_mem_ij_D(i_recv),8)
      Call CC_mpi_real_irecv(s_tmp,w_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    end if

    n_ij = CC_mem_ij_D(i_task)
    ij_start = CC_index_ij_D(i_task)

    do ij_run = 1, n_ij
      do a_run = 1, n_a
        do bbb = 1, CC_n_vir
          
          aaa = a_start - 1 + a_run
          i_ij = ij_start - 1 + ij_run
          Call CC_cl_decode(i_ij,iii,jjj,CC_n_occ,2)
          CC_t_d(iii,jjj,a_run,bbb) = w_tmp(ij_run,aaa,bbb,1)
          CC_t_d(jjj,iii,a_run,bbb) = w_tmp(ij_run,bbb,aaa,1)

        end do
      end do
    end do

    Deallocate(w_tmp)

    if (i_recv.ne.CC_mpi_did+1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      Allocate(w_tmp(CC_mem_ij_D(i_recv),CC_n_vir,CC_n_vir,1),stat=errnum)
      Call check_allocation(errnum,'w_tmp in CC2')

      w_tmp = w_recv

      Deallocate(w_recv)

    end if

  end do

  Deallocate(w_d_all)

  !write(90+myid,*) 't_d w2t'
  !do a_run = 1, n_a
  !  do bbb = 1, CC_n_vir
  !    do iii = 1, CC_n_occ
  !      do jjj = 1, CC_n_occ
  !        aaa = a_start - 1 + a_run
  !        write(90+myid,*) iii,jjj,aaa,bbb,CC_t_d(iii,jjj,a_run,bbb)
  !      end do
  !    end do
  !  end do
  !end do

  End Subroutine CC_cl_w2t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_pyramidal(n,n_pyramidal)

  Implicit None

  Integer , intent(in) :: n
  Integer , intent(out) :: n_pyramidal

  n_pyramidal = (n ** 3 + 3 * n ** 2 + 2 * n) / 6

  End Subroutine CC_cl_pyramidal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_code_ijk(iii,jjj,kkk,n)

! iii >= jjj >= kkk

  Implicit None

  Integer , intent(in) :: iii,jjj,kkk
  Integer , intent(out) :: n

  Integer :: f,f_pyramidal,code_jk

  if ((iii.lt.jjj).or.(jjj.lt.kkk)) then
    write(use_unit,*) 'Error in CCSD(T) calculation: i < j or j < k'
  end if

  f = iii - 1
  Call CC_cl_pyramidal(f,f_pyramidal)
  Call CC_cl_code(code_jk,kkk,jjj,CC_n_occ,2)

  n = f_pyramidal + code_jk

  End Subroutine CC_cl_code_ijk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_decode_ijk(n,iii,jjj,kkk)

! iii >= jjj >= kkk
! n = pyramidal(i-1) + code_jk, thus 6n > (i-1) ** 3
! since n = pyramidal(i-1) + code_jk <= pyramidal(i), we have,
! (i+1) ** 3 = i ** 3 + 3 * i ** 2 + 3 * i + 1 > 6 * pyramidal(i) >= 6n.
! finally, (i+1) ** 3 > 6n > (i-1) ** 3, and thus,
! Int(cbrt(6n)) = i or i - 1

  Implicit None

  Integer , intent(in) :: n
  Integer , intent(out) :: iii,jjj,kkk

  Integer :: f,f_pyramidal,code_jk
  Double precision :: r

  r = dble(n * 6) ** (1.0D0 / 3.0D0)

  f = Int(r)

  Call CC_cl_pyramidal(f,f_pyramidal)

  if (f_pyramidal.ge.n) then
    f = f - 1
  end if

  iii = f + 1

  Call CC_cl_pyramidal(f,f_pyramidal)

  code_jk = n - f_pyramidal

  Call CC_cl_decode(code_jk,kkk,jjj,CC_n_occ,2)

  End Subroutine CC_cl_decode_ijk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_decode_ijk_ne(n,iii,jjj,kkk)

! iii > jjj > kkk (iii>=3, jjj>=2, or n=0)
! n = pyramidal(i-3) + code_jk, thus 6n > (i-3) ** 3
! since n = pyramidal(i-3) + code_jk <= pyramidal(i-2), we have,
! (i-1) ** 3 = i ** 3 - 3 * i ** 2 + 3 * i - 1 > 6 * pyramidal(i-2) >= 6n.
! finally, (i-1) ** 3 > 6n > (i-3) ** 3, and thus,
! Int(cbrt(6n)) = i or i - 1

  Implicit None

  Integer , intent(in) :: n
  Integer , intent(out) :: iii,jjj,kkk

  Integer :: f,f_pyramidal,code_jk
  Double precision :: r

  r = dble(n * 6) ** (1.0D0 / 3.0D0)

  ! f = i - 3 or f = i - 2
  f = Int(r)

  ! make sure f = i - 3
  Call CC_cl_pyramidal(f,f_pyramidal)

  if (f_pyramidal.ge.n) then
    f = f - 1
  end if

  ! now f = i - 3
  iii = f + 3

  Call CC_cl_pyramidal(f,f_pyramidal)

  code_jk = n - f_pyramidal

  Call CC_cl_decode(code_jk,kkk,jjj,CC_n_occ,3)

  End Subroutine CC_cl_decode_ijk_ne



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Mat_regular(i_len,i_vec,j_len,j_vec,i_start,i_end)
  
  Implicit None

  Integer , intent(in) :: i_len, j_len, i_start, i_end
  Double precision , dimension(i_len) :: i_vec
  Double precision , dimension(j_len) :: j_vec

  i_vec(i_start:i_end) = j_vec(:)

  End Subroutine CC_Mat_regular

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_D_ijk(iii,jjj,kkk,delta_ijk,D_ijk)

  Implicit None

  Integer , intent(in) :: iii,jjj,kkk
  Double precision , intent(out) :: delta_ijk,D_ijk
  
  Integer :: i_i,i_j,i_k

  delta_ijk = 2.0D0

  if (iii.eq.jjj) then
    delta_ijk = delta_ijk - 1.0D0
  end if

  if (jjj.eq.kkk) then
    delta_ijk = delta_ijk - 1.0D0
  end if

  i_i = CC_n_fc + iii
  i_j = CC_n_fc + jjj
  i_k = CC_n_fc + kkk

  D_ijk = KS_eigenvalue(i_i,1,1) + KS_eigenvalue(i_j,1,1) + KS_eigenvalue(i_k,1,1)

  End Subroutine CC_cl_D_ijk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_D_abc(aaa,bbb,ccc,D_abc)

  Implicit None

  Integer , intent(in) :: aaa,bbb,ccc
  Double precision , intent(out) :: D_abc
  
  Integer :: i_a,i_b,i_c

  i_a = CC_n_elec + aaa
  i_b = CC_n_elec + bbb
  i_c = CC_n_elec + ccc

  D_abc = KS_eigenvalue(i_a,1,1) + KS_eigenvalue(i_b,1,1) + KS_eigenvalue(i_c,1,1)

  End Subroutine CC_cl_D_abc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_write_t(i_scf)

  Implicit None

  Integer , parameter :: ures = 100
  Integer :: utmp
  Integer , intent(in) :: i_scf
  Integer :: iii,errnum,n_s,n_d,n_rec,i_rec,n_ij
  Integer , dimension(1) :: shp1
  Character (len=9) :: fname
  Double precision , dimension(:) , allocatable :: w_tmp
  
  if (myid.eq.0) then
    open(unit=ures,file='CCrestart')
      write(ures,*) i_scf
      do iii = 1, CC_n_domain
        write(ures,*) 0
      end do
      write(ures,*) 0.0D0
    close(ures)
  end if

  n_s = CC_n_occ * CC_n_vir
  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  n_d = n_ij * CC_mem_aa_G(CC_mpi_gid+1) * CC_n_vir

  n_rec = n_s + n_d

  Allocate(w_tmp(n_rec),stat=errnum)
  Call check_allocation(errnum,'w_tmp in CC')

  i_rec = 0
  shp1(1) = n_s
  w_tmp(i_rec+1:i_rec+n_s) = reshape(CC_w_s,shp1)

  i_rec = n_s
  shp1(1) = n_d
  w_tmp(i_rec+1:i_rec+n_d) = reshape(CC_w_d,shp1)

  utmp = 200 + myid
  write(fname,"(I4)") myid
  do iii = 1, 4
    if (fname(iii:iii).eq.' ') then
      fname(iii:iii) ='0'
    end if
  end do

  fname = 'CCres'//fname
    
  open(unit=utmp,file=fname,form='unformatted',access='direct',RECL=8*n_rec)
  write(utmp,REC=1) w_tmp
  close(utmp)

  Deallocate(w_tmp)

  End Subroutine CC_cl_write_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_read_t()

  Implicit None

  Integer , parameter :: ures = 100
  Integer :: utmp
  Integer :: iii,errnum,n_s,n_d,n_rec,i_rec,n_ij
  Integer , dimension(2) :: shp0
  Integer , dimension(4) :: shp1
  Character (len=9) :: fname
  Double precision , dimension(:) , allocatable :: w_tmp
 
  n_s = CC_n_occ * CC_n_vir
  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  n_d = n_ij * CC_mem_aa_G(CC_mpi_gid+1) * CC_n_vir

  n_rec = n_s + n_d

  Allocate(w_tmp(n_rec),stat=errnum)
  Call check_allocation(errnum,'w_tmp in CC')

  utmp = 200 + myid
  write(fname,"(I4)") myid
  do iii = 1, 4
    if (fname(iii:iii).eq.' ') then
      fname(iii:iii) ='0'
    end if
  end do
  
  fname = 'CCres'//fname

  open(unit=utmp,file=fname,form='unformatted',access='direct',RECL=8*n_rec)
  read(utmp,REC=1) w_tmp
  close(utmp)

  i_rec = 0
  shp0(1) = CC_n_occ
  shp0(2) = CC_n_vir
  CC_w_s = reshape(w_tmp(i_rec+1:i_rec+n_s),shp0)

  i_rec = n_s
  shp1(1) = n_ij
  shp1(2) = CC_mem_aa_G(CC_mpi_gid+1)
  shp1(3) = CC_n_vir
  shp1(4) = 1
  CC_w_d = reshape(w_tmp(i_rec+1:i_rec+n_d),shp1)

  Deallocate(w_tmp)

  End Subroutine CC_cl_read_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module CC_cl

