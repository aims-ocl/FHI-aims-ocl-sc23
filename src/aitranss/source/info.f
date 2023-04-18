c###########################################################
c      a i t r a n s s : ab initio transport simulations
c     (c)  2003-2012   : alexej bagrets,  andreas arnold
c                        florian weigend, ferdinand evers
c     institute of nanotechnology (int) &
c     institut fuer theorie der kondensierten materie (tkm)
c     karlsruhe institute of technology (kit)
c
c     author:         alexej.bagrets <at> kit.edu
c     date:           march-may 2008
c     last revision:  jan 2012
c###########################################################

      module info
c      prints out info-message
       implicit none
       
       contains 

       subroutine info_message
	print *
        print '(1x,a)', '==============================================================='
        print *          	
	print '(1x,a)', '     a i t r a n s s :  ab initio transport simulations'      
        print *
        print '(1x,a)', '     (c)  2003-2013  :  alexei bagrets,  andreas arnold '
	print '(1x,a)', '                        florian weigend, ferdinand evers' 
        print *  
	print '(1x,a)', '     institute of nanotechnology (int) &'
	print '(1x,a)', '     institut fuer theorie der kondensierten materie (tkm)'
        print '(1x,a)', '     karlsruhe institute of technology (kit)'
        print *           
        print '(1x,a)', '==============================================================='
        print *
	print '(1x,a)', 'When using this software, please cite the following references:'
	print *
	print '(1x,a)', '   A. Arnold, F. Weigend, and F. Evers, '
        print '(1x,a)', '   "Quantum chemistry calculations for molecules coupled to '
	print '(1x,a)', '   reservoirs: Formalism, implementation, and application to '
	print '(1x,a)', '   benzenedithiol", J. Chem. Phys. 126, 174101 (2007).'
	print *
	print '(1x,a)', '   J. Wilhelm, M. Walz, M. Stendel, A. Bagrets, and F. Evers, '
        print '(1x,a)', '   "Ab initio simulations of scanning-tunneling-microscope'
	print '(1x,a)', '   images with embedding techniques and application to C58-'
	print '(1x,a)', '   dimers on Au(111)", Phys. Chem. Chem. Phys. 15, 6684 (2013).'
	print *
	print '(1x,a)', '   A. Bagrets, "Spin-polarizd electron transport across metal-'
	print '(1x,a)', '   organic molecules: a density functional theory approach",'
        print '(1x,a)', '   J. Chem. Theory Comput. 9, 2801 (2013). '
        print *
        print '(1x,a)', 'Report bugs to: alexej.bagrets <at> kit.edu'
        print *
       end subroutine info_message
      
      end module info

