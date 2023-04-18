c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     last revision:  jan 2012
c###########################################################

      module periodic_table
       
       implicit none
       
       contains 

       integer function get_charge(element)
       
        character(*) element
  
        select case(trim(element))
	
	 case('h')  ; get_charge =  1
	 case('he') ; get_charge =  2
	 
	 case('li') ; get_charge =  3
	 case('be') ; get_charge =  4
	 case('b')  ; get_charge =  5
	 case('c')  ; get_charge =  6
	 case('n')  ; get_charge =  7
	 case('o')  ; get_charge =  8
	 case('f')  ; get_charge =  9
	 case('ne') ; get_charge = 10
	 
	 case('na') ; get_charge = 11
	 case('mg') ; get_charge = 12
	 case('al') ; get_charge = 13
	 case('si') ; get_charge = 14
	 case('p')  ; get_charge = 15
	 case('s')  ; get_charge = 16
	 case('cl') ; get_charge = 17
	 case('ar') ; get_charge = 18
	 
	 case('k')  ; get_charge = 19
	 case('ca') ; get_charge = 20

	 case('sc') ; get_charge = 21
	 case('ti') ; get_charge = 22
	 case('v')  ; get_charge = 23
	 case('cr') ; get_charge = 24
	 case('mn') ; get_charge = 25
	 case('fe') ; get_charge = 26
	 case('co') ; get_charge = 27
	 case('ni') ; get_charge = 28
	 case('cu') ; get_charge = 29
	 case('zn') ; get_charge = 30
	
	 case('ga') ; get_charge = 31
	 case('ge') ; get_charge = 32
	 case('as') ; get_charge = 33
	 case('se') ; get_charge = 34
	 case('br') ; get_charge = 35
	 case('kr') ; get_charge = 36
	 
	 case('rb') ; get_charge = 37
	 case('sr') ; get_charge = 38 
	
	 case('y')  ; get_charge = 39
	 case('zr') ; get_charge = 40
	 case('nb') ; get_charge = 41
	 case('mo') ; get_charge = 42
	 case('tc') ; get_charge = 43
	 case('ru') ; get_charge = 44
	 case('rh') ; get_charge = 45
	 case('pd') ; get_charge = 46
	 case('ag') ; get_charge = 47
	 case('cd') ; get_charge = 48
	
	 case('in') ; get_charge = 49
	 case('sn') ; get_charge = 50
	 case('sb') ; get_charge = 51
	 case('te') ; get_charge = 52
	 case('i')  ; get_charge = 53
	 case('xe') ; get_charge = 54

	 case('cs') ; get_charge = 55 
	 case('ba') ; get_charge = 56
	 
	 case('la') ; get_charge = 57
	 case('ce') ; get_charge = 58
	 case('pr') ; get_charge = 59
	 case('nd') ; get_charge = 60
	 case('pm') ; get_charge = 61
	 case('sm') ; get_charge = 62
	 case('eu') ; get_charge = 63
	 case('gd') ; get_charge = 64
	 case('tb') ; get_charge = 65
	 case('dy') ; get_charge = 66
	 case('ho') ; get_charge = 67
	 case('er') ; get_charge = 68
	 case('tm') ; get_charge = 69
	 case('yb') ; get_charge = 70
	
	 case('lu') ; get_charge = 71
	 case('hf') ; get_charge = 72 
	 case('ta') ; get_charge = 73
	 case('w')  ; get_charge = 74
	 case('re') ; get_charge = 75
	 case('os') ; get_charge = 76
	 case('ir') ; get_charge = 77
	 case('pt') ; get_charge = 78
	 case('au') ; get_charge = 79
	 case('hg') ; get_charge = 80
	
	 case('tl') ; get_charge = 81
	 case('pb') ; get_charge = 82
	 case('bi') ; get_charge = 83
	 case('po') ; get_charge = 84
	 case('at') ; get_charge = 85
	 case('rn') ; get_charge = 86

	 case('fr') ; get_charge = 87 
	 case('ra') ; get_charge = 88
	  
	 case('ac') ; get_charge = 89
	 case('th') ; get_charge = 90
	 case('pa') ; get_charge = 91
	 case('u')  ; get_charge = 92
	 case('np') ; get_charge = 93
	 case('pu') ; get_charge = 94
	 case('am') ; get_charge = 95
	 case('cm') ; get_charge = 96
	 case('bk') ; get_charge = 97
	 case('cf') ; get_charge = 98
	 case('es') ; get_charge = 99
	 case('fm') ; get_charge = 100 
	 case('md') ; get_charge = 101 
	 case('no') ; get_charge = 102
	 
	 case('lr') ; get_charge = 103 
	 case('rf') ; get_charge = 104
	 case('db') ; get_charge = 105
	 case('sg') ; get_charge = 106 
	 case('bh') ; get_charge = 107 
	 case('hs') ; get_charge = 108 
	 case('mt') ; get_charge = 109 
	 case('ds') ; get_charge = 110 
	 case('rg') ; get_charge = 111 
	 case('cn') ; get_charge = 112 

         case default ; get_charge = -1
 	
	end select

       end function get_charge

       integer function get_ncore(z)
         integer z

         if (z<=18) then 
          get_ncore = 0 
         else if ( (z>=19).and.(z<=36) ) then 
          get_ncore = 2
         else if ( (z>=37).and.(z<=54) ) then 
          get_ncore = 3
         else if ( (z>=55).and.(z<=86) ) then 
          get_ncore = 4
         else  
          get_ncore = 5
         end if 
      
       end function get_ncore
      
      end module periodic_table
      