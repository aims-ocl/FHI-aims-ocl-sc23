! This module contains any architecture specific (p*ke!) 
!  calls for the alpha case.
!
module arch_specific
  
  implicit none
  
contains
  !
!-------------------------------------------------------------------------
  real*4 function arch_rand()
    
    real*4 :: rand
    
    arch_rand = rand()
    
    return 
    
  end function arch_rand
  
  
!-------------------------------------------------------------------------
end module arch_specific
