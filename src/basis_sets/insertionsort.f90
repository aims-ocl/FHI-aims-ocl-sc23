!****s* FHI-aims/insertionsort
!  NAME
!    insertionsort
!  SYNOPSIS 
      subroutine insertionsort( array, length, perm, perm_inv )
!  PURPOSE
!    Sorts a given array of double precision reals using the insertionsort
!    algorithm.
      use runtime_choices
      implicit none

!  ARGUMENTS
      integer :: length
      integer :: perm(length)
      integer :: perm_inv(length)
      real*8 :: array(length)
!  INPUTS
!    o length -- number of elements to sort
!    o array -- array of elements to sort
!    o perm -- new permutation of the entries
!    o perm_inv -- inverse permutation of perm
!  OUTPUT
!    the array is sorted on exit
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


!     locals

      integer :: i_index_1, i_index_2
      real*8 :: value
      integer :: perm_aux(length,length)

      !(Rundong) I don't reorder the basis set here.
      if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
        do i_index_1 = 1, length
           array(i_index_1) = 10.d0 + 0.1d0 * i_index_1
        enddo
      endif


         do i_index_1 = 1, length, 1
            do i_index_2 = 1, length, 1

               perm_aux(i_index_1,i_index_2) = i_index_2

            enddo
            perm(i_index_1) = i_index_1

         enddo

!      end if

      do i_index_1 = 2, length, 1

!        write(use_unit,*) i_index_1, array(i_index_1)
         value = array(i_index_1)
         i_index_2 = i_index_1 - 1

         do while ( (i_index_2.gt.1) .and. (value.lt.array(i_index_2)) &
               )

            array(i_index_2 + 1) = array(i_index_2)

!            if (present(perm)) then
               perm_aux(i_index_1,i_index_2) = i_index_2 + 1
!            endif

            i_index_2 = i_index_2 - 1
!        write(use_unit,*) i_index_2, array(i_index_2)

         enddo

         if (i_index_2.eq.1) then

           if (value.lt.array(i_index_2)) then
             ! this is the case where array(i_index_1) was the
             ! largest entry of the array

             array(i_index_2+1) = array(i_index_2)
             perm_aux(i_index_1,i_index_2) = i_index_2 + 1

             i_index_2 = i_index_2-1
           end if
        end if


         array(i_index_2 + 1) = value
         perm_aux(i_index_1, i_index_1) = i_index_2 + 1

      enddo

!      if (present(perm)) then
         do i_index_1 = 2, length, 1
            do i_index_2 = 1, length, 1

               perm(i_index_2) = &
                    perm_aux(i_index_1,perm(i_index_2))

            enddo
         enddo

         do i_index_1 = 1, length, 1
            perm_inv(perm(i_index_1)) = i_index_1
         enddo
!      end if


      end subroutine insertionsort
!******
