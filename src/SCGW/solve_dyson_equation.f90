!-------------------------------------------
      subroutine solve_dyson_equation_v1 &
         (inv_green_fn_freq,      &
          self_energy_freq,       &
          new_green_fn_freq,      &
          exchange_self_energy,   & 
          xc_matr,                &
          hartree_pot,            &
          hartree_pot0,           &
          exchange_self_energy_0  &
          )

      use scgw_grid
      use runtime_choices, only: hybrid_coeff
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use physics

      implicit none
! INPUT
      complex*16  inv_green_fn_freq      (n_basis,n_basis,nomega,n_spin)
      complex*16  self_energy_freq       (n_basis,n_basis,nomega,n_spin)
      real*8      xc_matr                (n_basis,n_basis,n_spin)
      real*8      exchange_self_energy   (n_basis,n_basis,n_spin)
      real*8      hartree_pot            (n_basis,n_basis)
      real*8      hartree_pot0           (n_basis,n_basis)
      real*8      exchange_self_energy_0 (n_basis,n_basis,n_spin)

!OUTPUT
      complex*16  new_green_fn_freq(n_basis,n_basis,nomega,n_spin)

!AUXILIARY 
      logical      output
      complex*16   inv_new_green_fn_freq(n_basis,n_basis,nomega)
      integer      i_freq
      integer      i_spin
      integer      i_basis
      character*17 filename
      character*2  iter

      do i_spin = 1, n_spin, 1
        do i_freq = 1, nomega, 1
  
         if (.not. use_hartree_fock) then
  
           inv_new_green_fn_freq   (:,:,i_freq) =         &
              inv_green_fn_freq    (:,:,i_freq,i_spin)    &
            - self_energy_freq     (:,:,i_freq,i_spin)    &
            - exchange_self_energy (:,:,i_spin)           &
            - hartree_pot          (:,:)                  &
            + hartree_pot0         (:,:)                  &
            + xc_matr              (:,:,i_spin)           
  
         elseif(use_hartree_fock)then
  
          inv_new_green_fn_freq    (:,:,i_freq) =  &
              inv_green_fn_freq    (:,:,i_freq,i_spin)    &
            - self_energy_freq     (:,:,i_freq,i_spin)    &
            - exchange_self_energy (:,:,i_spin)           &
            - hartree_pot          (:,:)                  &
            + hartree_pot0         (:,:)                  & 
            + xc_matr              (:,:,i_spin)           &
            + hybrid_coeff *                              &
              exchange_self_energy_0 (:,:,i_spin)   
  
         endif
        enddo
  
        call invert_green (inv_new_green_fn_freq,&
             nomega, new_green_fn_freq(:,:,:,i_spin), n_basis)
      enddo

      end subroutine solve_dyson_equation_v1
!-----------------------------------------------

      subroutine solve_dyson_equation (inv_green_fn_freq, &
          self_energy_freq,omega, nomega, new_green_fn_freq,&
          exchange_self_energy, xc_matr, hartree_pot, &
!          inv_overlap_matrix, &
           exchange_self_energy_0, &
          new_chem_pot )


      use runtime_choices, only: hybrid_coeff
      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use physics


      implicit none
! INPUT
      integer nomega
      complex*16  inv_green_fn_freq(n_basis,n_basis,nomega)
      complex*16  self_energy_freq (n_basis,n_basis,nomega)
      real*8  omega(nomega)
      real*8  xc_matr(n_basis, n_basis)
      real*8  exchange_self_energy(n_basis,n_basis)
      real*8  hartree_pot(n_basis,n_basis)
      real*8  exchange_self_energy_0 (n_basis,n_basis)
!      real*8  inv_overlap_matrix (n_basis, n_basis)
      real*8  new_chem_pot
      real*8 ovlp_NAO_KS (n_states, n_basis)

!OUTPUT
      complex*16  new_green_fn_freq(n_basis,n_basis,nomega)

!AUXILIARY 
      logical output
      complex*16  inv_new_green_fn_freq(n_basis,n_basis,nomega)
      integer i_freq
      character*17 filename
      character*2 iter
      real*8 id (n_basis, n_basis)
      real*8 full_ovlp_matrix(n_basis, n_basis)
      integer i_index       
      integer i_basis
      !delta_mu = new_chem_pot - KS_eigenvalue(1,1,1)
!      inv_new_green_fn_freq(:,:,:) = dcmplx(0.d0,0.d0)
!      new_green_fn_freq    (:,:,:) = dcmplx(0.d0,0.d0)
 
!     i_index = 0

!     do i_basis = 1, n_basis, 1
!       do j_basis = 1, i_basis, 1
!         i_index = i_index + 1
!         full_ovlp_matrix(i_basis,j_basis) =&
!             overlap_matrix(i_index)
!         if(i_basis.ne.j_basis)then
!             full_ovlp_matrix(j_basis,i_basis) =&
!             full_ovlp_matrix(i_basis,j_basis)
!         endif
!       enddo
!     enddo

    
!    id(:,:) = 0.d0
!      inv_green_fn_freq(:,:,:)=0.d0
!      do i_freq = 1, nomega, 1
!       do i_basis = 1, n_basis , 1
!        do j_basis = 1, n_basis, 1
!         do i_state = 1, n_states, 1
!           inv_green_fn_freq(i_basis,j_basis,i_freq) =&
!          inv_green_fn_freq(i_basis,j_basis,i_freq)+& !id(i_basis,j_basis)+&
!            ovlp_NAO_KS (i_state,i_basis)*&
!            ovlp_NAO_KS (i_state,j_basis)*&
!            ((0.d0,1.d0)*omega(i_freq)-&
!            (KS_eigenvalue(i_state,1,1) ))
!         enddo
!        enddo
!       enddo
!      enddo

       !inv_new_green_fn_freq(:,:,:)=0.d0  
!      do i_freq = 1, nomega, 1
!       green_tmp = 0.d0
!      call dgemm ('N','N',n_basis,n_basis, n_basis,&
!      1.d0,  hartree_pot (1,1), n_basis, &
!      full_ovlp_matrix(1,1), n_basis, 0.d0,        &
!      id(1,1), n_basis)

!      hartree_pot = 0.d0
!      call dgemm ('N','N',n_basis,n_basis, n_basis,&
!      1.d0, full_ovlp_matrix(1,1), n_basis, &
!      id(1,1), n_basis, 0.d0, &
!      hartree_pot (1,1), n_basis)
!      enddo

       
      !self_energy_freq = 0.d0 !this gives hartree-fock 
!      if(myid.eq.0 )then 
!         write(use_unit,*) " HYBRID COEFF    : " , hybrid_coeff
!      endif 

      do i_freq = 1, nomega, 1

       if (.not. use_hartree_fock) then

         inv_new_green_fn_freq   (:,:,i_freq) =  &
            inv_green_fn_freq    (:,:,i_freq)    &
          - self_energy_freq     (:,:,i_freq)    &
          - exchange_self_energy (:,:)           &
          + xc_matr              (:,:)           &
          - 2.d0*hartree_pot     (:,:) 

       elseif(use_hartree_fock)then
! the xc_matr has to be replaced with the 
! starting exact exchange matrix
!        inv_new_green_fn_freq     (:,:,i_freq) = &
!           inv_green_fn_freq      (:,:,i_freq)   &
!         - self_energy_freq       (:,:,i_freq)   &
!         - exchange_self_energy   (:,:)          &
!         + exchange_self_energy_0 (:,:)          &
!         - hartree_pot(:,:)*2

        inv_new_green_fn_freq   (:,:,i_freq) =  &
            inv_green_fn_freq    (:,:,i_freq)    &
          - self_energy_freq     (:,:,i_freq)    &
          - exchange_self_energy (:,:)           &
          + xc_matr              (:,:)           &
          + (hybrid_coeff)* exchange_self_energy_0 (:,:)         &
          - 2.d0*hartree_pot     (:,:)

       endif
      enddo

      call invert_green (inv_new_green_fn_freq,&
           nomega, new_green_fn_freq, n_basis)

      output = .false.
      if(myid.eq.0)then
        if (output)then
         do i_basis = 1, n_basis, 1
          if( i_basis.lt.10 ) then
           write(iter,'(A,I1)') "0",i_basis
          else
           write(iter,'(I2)') i_basis
          endif
          filename = "new_gree_fr"//iter//".dat"
          open(77, file=filename)
           do i_freq = 1, nomega, 1
            write(77,*) omega(i_freq), &
                       real(new_green_fn_freq(i_basis,i_basis,i_freq)),&
                      aimag(new_green_fn_freq(i_basis,i_basis,i_freq))
           enddo
          close(77)
         enddo
        endif
      endif
      
      end subroutine solve_dyson_equation


