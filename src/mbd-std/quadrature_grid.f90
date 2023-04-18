module quadrature_grid_module

use mbdvdw_interface_module, only: dp  ! double-precision kind

implicit none

private

save

real(dp), dimension(:), allocatable, public :: casimir_omega
real(dp), dimension(:), allocatable, public :: casimir_omega_weight

public :: generate_grid

contains

subroutine generate_grid(npts_in, npts_out)
   integer, intent(in) :: npts_in
   integer, intent(out) :: npts_out
   if(.not.allocated(casimir_omega)) allocate(casimir_omega(npts_in+1))
   if(.not.allocated(casimir_omega_weight)) allocate(casimir_omega_weight(npts_in+1))
   npts_out = npts_in+1

   select case(npts_in)
      case(0)
         npts_out = 21
         if(allocated(casimir_omega))        deallocate(casimir_omega)
         if(allocated(casimir_omega_weight)) deallocate(casimir_omega_weight)

         if(.not.allocated(casimir_omega)) allocate(casimir_omega(npts_out))
         if(.not.allocated(casimir_omega_weight)) allocate(casimir_omega_weight(npts_out))

         call fhi_ames_grid()
      case(4)
          call gauss_legendre_grid4()
      case(5)
          call gauss_legendre_grid5()
      case(6)
          call gauss_legendre_grid6()
      case(7)
          call gauss_legendre_grid7()
      case(8)
          call gauss_legendre_grid8()
      case(9)
          call gauss_legendre_grid9()
      case(10)
          call gauss_legendre_grid10()
      case(11)
          call gauss_legendre_grid11()
      case(12)
          call gauss_legendre_grid12()
      case(13)
          call gauss_legendre_grid13()
      case(14)
          call gauss_legendre_grid14()
      case(15)
          call gauss_legendre_grid15()
      case(16)
          call gauss_legendre_grid16()
      case(17)
          call gauss_legendre_grid17()
      case(18)
          call gauss_legendre_grid18()
      case(19)
          call gauss_legendre_grid19()
      case(20)
          call gauss_legendre_grid20()
      case(21)
          call gauss_legendre_grid21()
      case(22)
          call gauss_legendre_grid22()
      case(23)
          call gauss_legendre_grid23()
      case(24)
          call gauss_legendre_grid24()
      case(25)
          call gauss_legendre_grid25()
      case(26)
          call gauss_legendre_grid26()
      case(27)
          call gauss_legendre_grid27()
      case(28)
          call gauss_legendre_grid28()
      case(29)
          call gauss_legendre_grid29()
      case(30)
          call gauss_legendre_grid30()
      case(31)
          call gauss_legendre_grid31()
      case(32)
          call gauss_legendre_grid32()
      case(33)
          call gauss_legendre_grid33()
      case(34)
          call gauss_legendre_grid34()
      case(35)
          call gauss_legendre_grid35()
      case(36)
          call gauss_legendre_grid36()
      case(37)
          call gauss_legendre_grid37()
      case(38)
          call gauss_legendre_grid38()
      case(39)
          call gauss_legendre_grid39()
      case(40)
          call gauss_legendre_grid40()
      case(41)
          call gauss_legendre_grid41()
      case(42)
          call gauss_legendre_grid42()
      case(43)
          call gauss_legendre_grid43()
      case(44)
          call gauss_legendre_grid44()
      case(45)
          call gauss_legendre_grid45()
      case(46)
          call gauss_legendre_grid46()
      case(47)
          call gauss_legendre_grid47()
      case(48)
          call gauss_legendre_grid48()
      case(49)
          call gauss_legendre_grid49()
      case(50)
          call gauss_legendre_grid50()
      case(55)
          call gauss_legendre_grid55()
      case(60)
          call gauss_legendre_grid60()
      case(65)
          call gauss_legendre_grid65()
      case(70)
          call gauss_legendre_grid70()
      case(75)
          call gauss_legendre_grid75()
      case(80)
          call gauss_legendre_grid80()
      case(85)
          call gauss_legendre_grid85()
      case(90)
          call gauss_legendre_grid90()
      case(95)
          call gauss_legendre_grid95()
      case(100)
          call gauss_legendre_grid100()
      case(110)
          call gauss_legendre_grid110()
      case(120)
          call gauss_legendre_grid120()
      case(130)
          call gauss_legendre_grid130()
      case(140)
          call gauss_legendre_grid140()
      case(150)
          call gauss_legendre_grid150()
      case(200)
          call gauss_legendre_grid200()
      case(250)
          call gauss_legendre_grid250()
      case(300)
          call gauss_legendre_grid300()
      case(350)
          call gauss_legendre_grid350()
      case(400)
          call gauss_legendre_grid400()
      case(450)
          call gauss_legendre_grid450()
   end select

end subroutine generate_grid

subroutine fhi_ames_grid()
   casimir_omega(1 ) = 0.0000000_DP; casimir_omega_weight(1 ) = 0.0000000_DP
   casimir_omega(2 ) = 0.0392901_DP; casimir_omega_weight(2 ) = 0.0786611_DP
   casimir_omega(3 ) = 0.1183580_DP; casimir_omega_weight(3 ) = 0.0796400_DP
   casimir_omega(4 ) = 0.1989120_DP; casimir_omega_weight(4 ) = 0.0816475_DP
   casimir_omega(5 ) = 0.2820290_DP; casimir_omega_weight(5 ) = 0.0847872_DP
   casimir_omega(6 ) = 0.3689190_DP; casimir_omega_weight(6 ) = 0.0892294_DP
   casimir_omega(7 ) = 0.4610060_DP; casimir_omega_weight(7 ) = 0.0952317_DP
   casimir_omega(8 ) = 0.5600270_DP; casimir_omega_weight(8 ) = 0.1031720_DP
   casimir_omega(9 ) = 0.6681790_DP; casimir_omega_weight(9 ) = 0.1136050_DP
   casimir_omega(10) = 0.7883360_DP; casimir_omega_weight(10) = 0.1273500_DP
   casimir_omega(11) = 0.9243900_DP; casimir_omega_weight(11) = 0.1456520_DP
   casimir_omega(12) = 1.0817900_DP; casimir_omega_weight(12) = 0.1704530_DP
   casimir_omega(13) = 1.2684900_DP; casimir_omega_weight(13) = 0.2049170_DP
   casimir_omega(14) = 1.4966100_DP; casimir_omega_weight(14) = 0.2544560_DP
   casimir_omega(15) = 1.7856300_DP; casimir_omega_weight(15) = 0.3289620_DP
   casimir_omega(16) = 2.1691700_DP; casimir_omega_weight(16) = 0.4480920_DP
   casimir_omega(17) = 2.7106200_DP; casimir_omega_weight(17) = 0.6556060_DP
   casimir_omega(18) = 3.5457300_DP; casimir_omega_weight(18) = 1.0659600_DP
   casimir_omega(19) = 5.0273400_DP; casimir_omega_weight(19) = 2.0635700_DP
   casimir_omega(20) = 8.4489600_DP; casimir_omega_weight(20) = 5.6851000_DP
   casimir_omega(21) = 25.451700_DP; casimir_omega_weight(21) = 50.955800_DP
return
endsubroutine fhi_ames_grid

subroutine gauss_legendre_grid4()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0447673888927603_DP ;  casimir_omega_weight(2) = 0.1205099587901076_DP
         casimir_omega(3) = 0.2955350568166513_DP ;  casimir_omega_weight(3) = 0.4358411270879674_DP
         casimir_omega(4) = 1.2181295981523519_DP ;  casimir_omega_weight(4) = 1.7964399307365695_DP
         casimir_omega(5) = 8.0415679561382358_DP ;  casimir_omega_weight(5) = 21.6472089833853722_DP
return
endsubroutine gauss_legendre_grid4

subroutine gauss_legendre_grid5()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0295313648167765_DP ;  casimir_omega_weight(2) = 0.0782470267057947_DP
         casimir_omega(3) = 0.1799960597963229_DP ;  casimir_omega_weight(3) = 0.2426622842865895_DP
         casimir_omega(4) = 0.6000000000000000_DP ;  casimir_omega_weight(4) = 0.6826666666666666_DP
         casimir_omega(5) = 2.0000437809992238_DP ;  casimir_omega_weight(5) = 2.6963656488905752_DP
         casimir_omega(6) = 12.1904287943877048_DP ;  casimir_omega_weight(6) = 32.3000583734505042_DP
return
endsubroutine gauss_legendre_grid5

subroutine gauss_legendre_grid6()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0209671051368779_DP ;  casimir_omega_weight(2) = 0.0550522910140742_DP
         casimir_omega(3) = 0.1223652898763400_DP ;  casimir_omega_weight(3) = 0.1568746483220375_DP
         casimir_omega(4) = 0.3688207751687427_DP ;  casimir_omega_weight(4) = 0.3659920223248676_DP
         casimir_omega(5) = 0.9760838440711290_DP ;  casimir_omega_weight(5) = 0.9685975522576772_DP
         casimir_omega(6) = 2.9420107643581699_DP ;  casimir_omega_weight(6) = 3.7717142212856860_DP
         casimir_omega(7) = 17.1697522213887446_DP ;  casimir_omega_weight(7) = 45.0817692647956818_DP
return
endsubroutine gauss_legendre_grid6

subroutine gauss_legendre_grid7()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0156662709134676_DP ;  casimir_omega_weight(2) = 0.0409005195431814_DP
         casimir_omega(3) = 0.0890488151591658_DP ;  casimir_omega_weight(3) = 0.1106673657026609_DP
         casimir_omega(4) = 0.2535790722216331_DP ;  casimir_omega_weight(4) = 0.2318335984138593_DP
         casimir_omega(5) = 0.6000000000000000_DP ;  casimir_omega_weight(5) = 0.5015510204081632_DP
         casimir_omega(6) = 1.4196755151992706_DP ;  casimir_omega_weight(6) = 1.2979323584756699_DP
         casimir_omega(7) = 4.0427264456751733_DP ;  casimir_omega_weight(7) = 5.0241868485243044_DP
         casimir_omega(8) = 22.9793038808312176_DP ;  casimir_omega_weight(8) = 59.9929282889317932_DP
return
endsubroutine gauss_legendre_grid7

subroutine gauss_legendre_grid8()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0121543689176908_DP ;  casimir_omega_weight(2) = 0.0316113918275405_DP
         casimir_omega(3) = 0.0679035954004373_DP ;  casimir_omega_weight(3) = 0.0826692623964860_DP
         casimir_omega(4) = 0.1866106234123315_DP ;  casimir_omega_weight(4) = 0.1617566291470207_DP
         casimir_omega(5) = 0.4139976952756907_DP ;  casimir_omega_weight(5) = 0.3107569334661239_DP
         casimir_omega(6) = 0.8695700582590624_DP ;  casimir_omega_weight(6) = 0.6527208431887426_DP
         casimir_omega(7) = 1.9291506207797742_DP ;  casimir_omega_weight(7) = 1.6722140241968553_DP
         casimir_omega(8) = 5.3016338513009194_DP ;  casimir_omega_weight(8) = 6.4544764883019106_DP
         casimir_omega(9) = 29.6189791866540517_DP ;  casimir_omega_weight(9) = 77.0337944274751010_DP
return
endsubroutine gauss_legendre_grid8

subroutine gauss_legendre_grid9()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0097064537286881_DP ;  casimir_omega_weight(2) = 0.0251775836703725_DP
         casimir_omega(3) = 0.0535836975808472_DP ;  casimir_omega_weight(3) = 0.0643064773559521_DP
         casimir_omega(4) = 0.1437840882005354_DP ;  casimir_omega_weight(4) = 0.1201447386947358_DP
         casimir_omega(5) = 0.3061709630438917_DP ;  casimir_omega_weight(5) = 0.2137354124241337_DP
         casimir_omega(6) = 0.6000000000000000_DP ;  casimir_omega_weight(6) = 0.3962872260015118_DP
         casimir_omega(7) = 1.1758136579019463_DP ;  casimir_omega_weight(7) = 0.8208257720036445_DP
         casimir_omega(8) = 2.5037540975876875_DP ;  casimir_omega_weight(8) = 2.0921152373342138_DP
         casimir_omega(9) = 6.7184613278475496_DP ;  casimir_omega_weight(9) = 8.0629109365624672_DP
         casimir_omega(10) = 37.0887257141087190_DP ;  casimir_omega_weight(10) = 96.2044966159523227_DP
return
endsubroutine gauss_legendre_grid9

subroutine gauss_legendre_grid10()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0079315219153048_DP ;  casimir_omega_weight(2) = 0.0205337037129960_DP
         casimir_omega(3) = 0.0434097743983572_DP ;  casimir_omega_weight(3) = 0.0515577435908860_DP
         casimir_omega(4) = 0.1145368364284179_DP ;  casimir_omega_weight(4) = 0.0932144721167302_DP
         casimir_omega(5) = 0.2371730542143752_DP ;  casimir_omega_weight(5) = 0.1572649408016086_DP
         casimir_omega(6) = 0.4445006553664406_DP ;  casimir_omega_weight(6) = 0.2686762476258545_DP
         casimir_omega(7) = 0.8098975685496360_DP ;  casimir_omega_weight(7) = 0.4895386250889388_DP
         casimir_omega(8) = 1.5178790069237980_DP ;  casimir_omega_weight(8) = 1.0064766967672136_DP
         casimir_omega(9) = 3.1430936214567891_DP ;  casimir_omega_weight(9) = 2.5579701856062806_DP
         casimir_omega(10) = 8.2930631404899433_DP ;  casimir_omega_weight(10) = 9.8496624068285676_DP
         casimir_omega(11) = 45.3885148202563826_DP ;  casimir_omega_weight(11) = 117.5051049778581103_DP
return
endsubroutine gauss_legendre_grid10

subroutine gauss_legendre_grid11()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0066032837299041_DP ;  casimir_omega_weight(2) = 0.0170701882586499_DP
         casimir_omega(3) = 0.0359089519062432_DP ;  casimir_omega_weight(3) = 0.0423185119223073_DP
         casimir_omega(4) = 0.0935806773820723_DP ;  casimir_omega_weight(4) = 0.0746797334966168_DP
         casimir_omega(5) = 0.1899434255201306_DP ;  casimir_omega_weight(5) = 0.1212628204825432_DP
         casimir_omega(6) = 0.3452219047251076_DP ;  casimir_omega_weight(6) = 0.1956677179249149_DP
         casimir_omega(7) = 0.6000000000000000_DP ;  casimir_omega_weight(7) = 0.3275101041334809_DP
         casimir_omega(8) = 1.0428075248778312_DP ;  casimir_omega_weight(8) = 0.5910510481374279_DP
         casimir_omega(9) = 1.8953011877836559_DP ;  casimir_omega_weight(9) = 1.2099895906647344_DP
         casimir_omega(10) = 3.8469480032740915_DP ;  casimir_omega_weight(10) = 3.0699612323482568_DP
         casimir_omega(11) = 10.0253552634993408_DP ;  casimir_omega_weight(11) = 11.8148287187964467_DP
         casimir_omega(12) = 54.5183297773015667_DP ;  casimir_omega_weight(12) = 140.9356603338343632_DP
return
endsubroutine gauss_legendre_grid11

subroutine gauss_legendre_grid12()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0055832858509395_DP ;  casimir_omega_weight(2) = 0.0144172198045948_DP
         casimir_omega(3) = 0.0302132896412981_DP ;  casimir_omega_weight(3) = 0.0353941356300697_DP
         casimir_omega(4) = 0.0780033826132639_DP ;  casimir_omega_weight(4) = 0.0613218171925034_DP
         casimir_omega(5) = 0.1559922048127476_DP ;  casimir_omega_weight(5) = 0.0967625865124997_DP
         casimir_omega(6) = 0.2773010424740895_DP ;  casimir_omega_weight(6) = 0.1497576608413012_DP
         casimir_omega(7) = 0.4664454067244921_DP ;  casimir_omega_weight(7) = 0.2361303180272326_DP
         casimir_omega(8) = 0.7717945011572072_DP ;  casimir_omega_weight(8) = 0.3907082766441817_DP
         casimir_omega(9) = 1.2982280801690014_DP ;  casimir_omega_weight(9) = 0.7011138464896651_DP
         casimir_omega(10) = 2.3078076268756029_DP ;  casimir_omega_weight(10) = 1.4315422710885908_DP
         casimir_omega(11) = 4.6151844694333111_DP ;  casimir_omega_weight(11) = 3.6281951995264681_DP
         casimir_omega(12) = 11.9152864277288337_DP ;  casimir_omega_weight(12) = 13.9584689022973638_DP
         casimir_omega(13) = 64.4781602825191698_DP ;  casimir_omega_weight(13) = 166.4961877659453933_DP
return
endsubroutine gauss_legendre_grid12

subroutine gauss_legendre_grid13()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0047829090901057_DP ;  casimir_omega_weight(2) = 0.0123396045114688_DP
         casimir_omega(3) = 0.0257827501765995_DP ;  casimir_omega_weight(3) = 0.0300626270559082_DP
         casimir_omega(4) = 0.0660827005903224_DP ;  casimir_omega_weight(4) = 0.0513445646366637_DP
         casimir_omega(5) = 0.1306606281516948_DP ;  casimir_omega_weight(5) = 0.0792549131176798_DP
         casimir_omega(6) = 0.2284473630547052_DP ;  casimir_omega_weight(6) = 0.1188577964640809_DP
         casimir_omega(7) = 0.3752463650656124_DP ;  casimir_omega_weight(7) = 0.1793493092495169_DP
         casimir_omega(8) = 0.6000000000000000_DP ;  casimir_omega_weight(8) = 0.2790618638770487_DP
         casimir_omega(9) = 0.9593697195096172_DP ;  casimir_omega_weight(9) = 0.4585315476110403_DP
         casimir_omega(10) = 1.5758553532254711_DP ;  casimir_omega_weight(10) = 0.8198943175616906_DP
         casimir_omega(11) = 2.7552293685749465_DP ;  casimir_omega_weight(11) = 1.6712415003253209_DP
         casimir_omega(12) = 5.4477192485187391_DP ;  casimir_omega_weight(12) = 4.2327382292080564_DP
         casimir_omega(13) = 13.9628238855115185_DP ;  casimir_omega_weight(13) = 16.2806203466391182_DP
         casimir_omega(14) = 75.2679997085300840_DP ;  casimir_omega_weight(14) = 194.1867033797393844_DP
return
endsubroutine gauss_legendre_grid13

subroutine gauss_legendre_grid14()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0041432723490365_DP ;  casimir_omega_weight(2) = 0.0106818499929661_DP
         casimir_omega(3) = 0.0222662793364750_DP ;  casimir_omega_weight(3) = 0.0258653663020621_DP
         casimir_omega(4) = 0.0567420842482167_DP ;  casimir_omega_weight(4) = 0.0436768290008541_DP
         casimir_omega(5) = 0.1111983915643437_DP ;  casimir_omega_weight(5) = 0.0662615478989440_DP
         casimir_omega(6) = 0.1919492360568343_DP ;  casimir_omega_weight(6) = 0.0969721988960306_DP
         casimir_omega(7) = 0.3097026366110883_DP ;  casimir_omega_weight(7) = 0.1415115102191963_DP
         casimir_omega(8) = 0.4829787831369910_DP ;  casimir_omega_weight(8) = 0.2103922610957832_DP
         casimir_omega(9) = 0.7453743571545055_DP ;  casimir_omega_weight(9) = 0.3246954148709923_DP
         casimir_omega(10) = 1.1624053444920233_DP ;  casimir_omega_weight(10) = 0.5311344378139596_DP
         casimir_omega(11) = 1.8754958727390176_DP ;  casimir_omega_weight(11) = 0.9474950905565596_DP
         casimir_omega(12) = 3.2374568996502977_DP ;  casimir_omega_weight(12) = 1.9291547513330358_DP
         casimir_omega(13) = 6.3444972945510765_DP ;  casimir_omega_weight(13) = 4.8836331464013432_DP
         casimir_omega(14) = 16.1679459131851431_DP ;  casimir_omega_weight(14) = 18.7813076929926659_DP
         casimir_omega(15) = 86.8878436349269521_DP ;  casimir_omega_weight(15) = 224.0072179026357162_DP
return
endsubroutine gauss_legendre_grid14

subroutine gauss_legendre_grid15()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0036240021641946_DP ;  casimir_omega_weight(2) = 0.0093377589930291_DP
         casimir_omega(3) = 0.0194272861575501_DP ;  casimir_omega_weight(3) = 0.0224989668490614_DP
         casimir_omega(4) = 0.0492780681398096_DP ;  casimir_omega_weight(4) = 0.0376452136210375_DP
         casimir_omega(5) = 0.0958870685315184_DP ;  casimir_omega_weight(5) = 0.0563236089195038_DP
         casimir_omega(6) = 0.1638582152651668_DP ;  casimir_omega_weight(6) = 0.0808455433033145_DP
         casimir_omega(7) = 0.2607386870266663_DP ;  casimir_omega_weight(7) = 0.1149344187522071_DP
         casimir_omega(8) = 0.3990059108653618_DP ;  casimir_omega_weight(8) = 0.1650309702105489_DP
         casimir_omega(9) = 0.6000000000000000_DP ;  casimir_omega_weight(9) = 0.2430938903106735_DP
         casimir_omega(10) = 0.9022422731012530_DP ;  casimir_omega_weight(10) = 0.3731722103362876_DP
         casimir_omega(11) = 1.3806926931529038_DP ;  casimir_omega_weight(11) = 0.6086136045730677_DP
         casimir_omega(12) = 2.1970213664137801_DP ;  casimir_omega_weight(12) = 1.0839821838000396_DP
         casimir_omega(13) = 3.7544165810186048_DP ;  casimir_omega_weight(13) = 2.2053264790411675_DP
         casimir_omega(14) = 7.3054811925383074_DP ;  casimir_omega_weight(14) = 5.5809087182012016_DP
         casimir_omega(15) = 18.5306376341242682_DP ;  casimir_omega_weight(15) = 21.4605477286440589_DP
         casimir_omega(16) = 99.3376890214983490_DP ;  casimir_omega_weight(16) = 255.9577387044332113_DP
return
endsubroutine gauss_legendre_grid15

subroutine gauss_legendre_grid16()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0031966603077105_DP ;  casimir_omega_weight(2) = 0.0082327662303926_DP
         casimir_omega(3) = 0.0171014158679792_DP ;  casimir_omega_weight(3) = 0.0197558527233351_DP
         casimir_omega(4) = 0.0432139419967426_DP ;  casimir_omega_weight(4) = 0.0328078139642496_DP
         casimir_omega(5) = 0.0836031596414442_DP ;  casimir_omega_weight(5) = 0.0485339777703896_DP
         casimir_omega(6) = 0.1417130971244756_DP ;  casimir_omega_weight(6) = 0.0685820714448462_DP
         casimir_omega(7) = 0.2230357965620847_DP ;  casimir_omega_weight(7) = 0.0954871526976804_DP
         casimir_omega(8) = 0.3363269938432668_DP ;  casimir_omega_weight(8) = 0.1334082654269864_DP
         casimir_omega(9) = 0.4958778910917917_DP ;  casimir_omega_weight(9) = 0.1896003320385234_DP
         casimir_omega(10) = 0.7259851799550396_DP ;  casimir_omega_weight(10) = 0.2775825130486468_DP
         casimir_omega(11) = 1.0703868752436956_DP ;  casimir_omega_weight(11) = 0.4245822041528421_DP
         casimir_omega(12) = 1.6140906775912522_DP ;  casimir_omega_weight(12) = 0.6910322261034716_DP
         casimir_omega(13) = 2.5403438870846857_DP ;  casimir_omega_weight(13) = 1.2293997484614272_DP
         casimir_omega(14) = 4.3060573493150462_DP ;  casimir_omega_weight(14) = 2.4997870004673559_DP
         casimir_omega(15) = 8.3306447726323167_DP ;  casimir_omega_weight(15) = 6.3245848740985258_DP
         casimir_omega(16) = 21.0508885801711045_DP ;  casimir_omega_weight(16) = 24.3183522168997968_DP
         casimir_omega(17) = 112.6175337215725278_DP ;  casimir_omega_weight(17) = 290.0382709844776628_DP
return
endsubroutine gauss_legendre_grid16

subroutine gauss_legendre_grid17()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0028407437354047_DP ;  casimir_omega_weight(2) = 0.0073132523944734_DP
         casimir_omega(3) = 0.0151715067977603_DP ;  casimir_omega_weight(3) = 0.0174899012524414_DP
         casimir_omega(4) = 0.0382166851601698_DP ;  casimir_omega_weight(4) = 0.0288641413461699_DP
         casimir_omega(5) = 0.0735843767577327_DP ;  casimir_omega_weight(5) = 0.0423029015299237_DP
         casimir_omega(6) = 0.1239071472818789_DP ;  casimir_omega_weight(6) = 0.0590142275493950_DP
         casimir_omega(7) = 0.1932884952868578_DP ;  casimir_omega_weight(7) = 0.0807850165079439_DP
         casimir_omega(8) = 0.2880785905540633_DP ;  casimir_omega_weight(8) = 0.1104183977062747_DP
         casimir_omega(9) = 0.4182571975440877_DP ;  casimir_omega_weight(9) = 0.1525571987988791_DP
         casimir_omega(10) = 0.6000000000000000_DP ;  casimir_omega_weight(10) = 0.2153357644274479_DP
         casimir_omega(11) = 0.8607144171429425_DP ;  casimir_omega_weight(11) = 0.3139412333276018_DP
         casimir_omega(12) = 1.2496589882212690_DP ;  casimir_omega_weight(12) = 0.4789850675582968_DP
         casimir_omega(13) = 1.8625009184624628_DP ;  casimir_omega_weight(13) = 0.7784331251622152_DP
         casimir_omega(14) = 2.9054014065954448_DP ;  casimir_omega_weight(14) = 1.3837782847271745_DP
         casimir_omega(15) = 4.8923428567623057_DP ;  casimir_omega_weight(15) = 2.8125576003943316_DP
         casimir_omega(16) = 9.4199692749699544_DP ;  casimir_omega_weight(16) = 7.1146757885921872_DP
         casimir_omega(17) = 23.7286912103643992_DP ;  casimir_omega_weight(17) = 27.3547296027453584_DP
         casimir_omega(18) = 126.7273761843634503_DP ;  casimir_omega_weight(18) = 326.2488184959813680_DP
return
endsubroutine gauss_legendre_grid17

subroutine gauss_legendre_grid18()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0025411666299899_DP ;  casimir_omega_weight(2) = 0.0065398502718662_DP
         casimir_omega(3) = 0.0135521554805430_DP ;  casimir_omega_weight(3) = 0.0155957128145076_DP
         casimir_omega(4) = 0.0340475727164809_DP ;  casimir_omega_weight(4) = 0.0256036593098562_DP
         casimir_omega(5) = 0.0652972782663843_DP ;  casimir_omega_weight(5) = 0.0372325123754240_DP
         casimir_omega(6) = 0.1093510616651265_DP ;  casimir_omega_weight(6) = 0.0513893313482539_DP
         casimir_omega(7) = 0.1693437882626418_DP ;  casimir_omega_weight(7) = 0.0693709297605707_DP
         casimir_omega(8) = 0.2500081549474851_DP ;  casimir_omega_weight(8) = 0.0931348518718314_DP
         casimir_omega(9) = 0.3585535613167681_DP ;  casimir_omega_weight(9) = 0.1257844405534405_DP
         casimir_omega(10) = 0.5062201706095452_DP ;  casimir_omega_weight(10) = 0.1724861962059178_DP
         casimir_omega(11) = 0.7111530138487372_DP ;  casimir_omega_weight(11) = 0.2423136915533056_DP
         casimir_omega(12) = 1.0040340937569268_DP ;  casimir_omega_weight(12) = 0.3522259444753406_DP
         casimir_omega(13) = 1.4399530290346687_DP ;  casimir_omega_weight(13) = 0.5364217502813420_DP
         casimir_omega(14) = 2.1258529981723462_DP ;  casimir_omega_weight(14) = 0.8708462266640207_DP
         casimir_omega(15) = 3.2921491069053670_DP ;  casimir_omega_weight(15) = 1.5471394490957384_DP
         casimir_omega(16) = 5.5132466399496378_DP ;  casimir_omega_weight(16) = 3.1436535978309821_DP
         casimir_omega(17) = 10.5734409614973952_DP ;  casimir_omega_weight(17) = 7.9511917740912859_DP
         casimir_omega(18) = 26.5640399799763287_DP ;  casimir_omega_weight(18) = 30.5696860780267521_DP
         casimir_omega(19) = 141.6672152669609659_DP ;  casimir_omega_weight(19) = 364.5893840034558480_DP
return
endsubroutine gauss_legendre_grid18

subroutine gauss_legendre_grid19()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0022866282094578_DP ;  casimir_omega_weight(2) = 0.0058831231422575_DP
         casimir_omega(3) = 0.0121798844133466_DP ;  casimir_omega_weight(3) = 0.0139956402750500_DP
         casimir_omega(4) = 0.0305316331261895_DP ;  casimir_omega_weight(4) = 0.0228750404568153_DP
         casimir_omega(5) = 0.0583586716089737_DP ;  casimir_omega_weight(5) = 0.0330459015539304_DP
         casimir_omega(6) = 0.0972827332714340_DP ;  casimir_omega_weight(6) = 0.0452033843393872_DP
         casimir_omega(7) = 0.1497444755265162_DP ;  casimir_omega_weight(7) = 0.0603123022102896_DP
         casimir_omega(8) = 0.2193527059489135_DP ;  casimir_omega_weight(8) = 0.0797811836679323_DP
         casimir_omega(9) = 0.3114634069341174_DP ;  casimir_omega_weight(9) = 0.1057606364696072_DP
         casimir_omega(10) = 0.4341630189155032_DP ;  casimir_omega_weight(10) = 0.1416800741869687_DP
         casimir_omega(11) = 0.6000000000000000_DP ;  casimir_omega_weight(11) = 0.1932653398185404_DP
         casimir_omega(12) = 0.8291816306677727_DP ;  casimir_omega_weight(12) = 0.2705861849793922_DP
         casimir_omega(13) = 1.1558340144791046_DP ;  casimir_omega_weight(13) = 0.3924754507369421_DP
         casimir_omega(14) = 1.6411924277051899_DP ;  casimir_omega_weight(14) = 0.5969211728788163_DP
         casimir_omega(15) = 2.4040953680207893_DP ;  casimir_omega_weight(15) = 0.9682929929041140_DP
         casimir_omega(16) = 3.7005539204531184_DP ;  casimir_omega_weight(16) = 1.7194989851707589_DP
         casimir_omega(17) = 6.1687490491926065_DP ;  casimir_omega_weight(17) = 3.4930862572816821_DP
         casimir_omega(18) = 11.7910495816615200_DP ;  casimir_omega_weight(18) = 8.8341404828902164_DP
         casimir_omega(19) = 29.5569307378249810_DP ;  casimir_omega_weight(19) = 33.9632262673914411_DP
         casimir_omega(20) = 157.4370501120370420_DP ;  casimir_omega_weight(20) = 405.0599695796285005_DP
return
endsubroutine gauss_legendre_grid19

subroutine gauss_legendre_grid20()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0020685270838163_DP ;  casimir_omega_weight(2) = 0.0053206999983087_DP
         casimir_omega(3) = 0.0110066968539682_DP ;  casimir_omega_weight(3) = 0.0126314155336561_DP
         casimir_omega(4) = 0.0275381209914517_DP ;  casimir_omega_weight(4) = 0.0205670909744248_DP
         casimir_omega(5) = 0.0524870458951971_DP ;  casimir_omega_weight(5) = 0.0295451545147365_DP
         casimir_omega(6) = 0.0871545984820398_DP ;  casimir_omega_weight(6) = 0.0401079258170448_DP
         casimir_omega(7) = 0.1334722656943147_DP ;  casimir_omega_weight(7) = 0.0529887326875466_DP
         casimir_omega(8) = 0.1942459517949385_DP ;  casimir_omega_weight(8) = 0.0692272502150060_DP
         casimir_omega(9) = 0.2735493056758191_DP ;  casimir_omega_weight(9) = 0.0903599093386727_DP
         casimir_omega(10) = 0.3773691388316548_DP ;  casimir_omega_weight(10) = 0.1187479666670754_DP
         casimir_omega(11) = 0.5146961792789785_DP ;  casimir_omega_weight(11) = 0.1581694585910045_DP
         casimir_omega(12) = 0.6994417570853403_DP ;  casimir_omega_weight(12) = 0.2149429673037550_DP
         casimir_omega(13) = 0.9539730808792946_DP ;  casimir_omega_weight(13) = 0.3001897928385637_DP
         casimir_omega(14) = 1.3160333165920470_DP ;  casimir_omega_weight(14) = 0.4347174301179076_DP
         casimir_omega(15) = 1.8533204768151084_DP ;  casimir_omega_weight(15) = 0.6605042688998619_DP
         casimir_omega(16) = 2.6971895481604484_DP ;  casimir_omega_weight(16) = 1.0707891653120136_DP
         casimir_omega(17) = 4.1305909988695122_DP ;  casimir_omega_weight(17) = 1.9008685743340461_DP
         casimir_omega(18) = 6.8588352394384327_DP ;  casimir_omega_weight(18) = 3.8608640186181988_DP
         casimir_omega(19) = 13.0727873594480659_DP ;  casimir_omega_weight(19) = 9.7635276929221710_DP
         casimir_omega(20) = 32.7073603258374490_DP ;  casimir_omega_weight(20) = 37.5353536820377656_DP
         casimir_omega(21) = 174.0368800662883473_DP ;  casimir_omega_weight(21) = 447.6605768032585502_DP
return
endsubroutine gauss_legendre_grid20

subroutine gauss_legendre_grid21()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0018802224684725_DP ;  casimir_omega_weight(2) = 0.0048353316170273_DP
         casimir_omega(3) = 0.0099957445042520_DP ;  casimir_omega_weight(3) = 0.0114585944343627_DP
         casimir_omega(4) = 0.0249676663374147_DP ;  casimir_omega_weight(4) = 0.0185965214879633_DP
         casimir_omega(5) = 0.0474715228170005_DP ;  casimir_omega_weight(5) = 0.0265855349248210_DP
         casimir_omega(6) = 0.0785641722558279_DP ;  casimir_omega_weight(6) = 0.0358553522675121_DP
         casimir_omega(7) = 0.1197960943496240_DP ;  casimir_omega_weight(7) = 0.0469738157489542_DP
         casimir_omega(8) = 0.1733858163134742_DP ;  casimir_omega_weight(8) = 0.0607254094743330_DP
         casimir_omega(9) = 0.2424942174884458_DP ;  casimir_omega_weight(9) = 0.0782366921338333_DP
         casimir_omega(10) = 0.3316615993429985_DP ;  casimir_omega_weight(10) = 0.1011844386881142_DP
         casimir_omega(11) = 0.4475209135511761_DP ;  casimir_omega_weight(11) = 0.1321555314958571_DP
         casimir_omega(12) = 0.6000000000000000_DP ;  casimir_omega_weight(12) = 0.1752973603796285_DP
         casimir_omega(13) = 0.8044316792780062_DP ;  casimir_omega_weight(13) = 0.2375533587547808_DP
         casimir_omega(14) = 1.0854437194813571_DP ;  casimir_omega_weight(14) = 0.3311508287387704_DP
         casimir_omega(15) = 1.4845714826876355_DP ;  casimir_omega_weight(15) = 0.4789720894983205_DP
         casimir_omega(16) = 2.0762944031658006_DP ;  casimir_omega_weight(16) = 0.7271865167653329_DP
         casimir_omega(17) = 3.0051063179851472_DP ;  casimir_omega_weight(17) = 1.1783465166658444_DP
         casimir_omega(18) = 4.5822413660483141_DP ;  casimir_omega_weight(18) = 2.0912570404156523_DP
         casimir_omega(19) = 7.5834938219230050_DP ;  casimir_omega_weight(19) = 4.2469933107496045_DP
         casimir_omega(20) = 14.4186483083735499_DP ;  casimir_omega_weight(20) = 10.7393578346665333_DP
         casimir_omega(21) = 36.0153263067962257_DP ;  casimir_omega_weight(21) = 41.2860710270519462_DP
         casimir_omega(22) = 191.4667046248343922_DP ;  casimir_omega_weight(22) = 492.3912068940510949_DP
return
endsubroutine gauss_legendre_grid21

subroutine gauss_legendre_grid22()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0017165210874463_DP ;  casimir_omega_weight(2) = 0.0044135437690269_DP
         casimir_omega(3) = 0.0091183500803553_DP ;  casimir_omega_weight(3) = 0.0104427820110461_DP
         casimir_omega(4) = 0.0227436013719973_DP ;  casimir_omega_weight(4) = 0.0168998808539809_DP
         casimir_omega(5) = 0.0431514152756930_DP ;  casimir_omega_weight(5) = 0.0240590603710706_DP
         casimir_omega(6) = 0.0712097100702361_DP ;  casimir_omega_weight(6) = 0.0322655209409634_DP
         casimir_omega(7) = 0.1081788255777496_DP ;  casimir_omega_weight(7) = 0.0419661879011242_DP
         casimir_omega(8) = 0.1558386630557780_DP ;  casimir_omega_weight(8) = 0.0537644390359346_DP
         casimir_omega(9) = 0.2166844061709273_DP ;  casimir_omega_weight(9) = 0.0685050493300289_DP
         casimir_omega(10) = 0.2942305446487530_DP ;  casimir_omega_weight(10) = 0.0874105548141960_DP
         casimir_omega(11) = 0.3934922723565226_DP ;  casimir_omega_weight(11) = 0.1123084426523416_DP
         casimir_omega(12) = 0.5217686682438419_DP ;  casimir_omega_weight(12) = 0.1460247294454859_DP
         casimir_omega(13) = 0.6899609384589543_DP ;  casimir_omega_weight(13) = 0.1930958401652003_DP
         casimir_omega(14) = 0.9148845486699239_DP ;  casimir_omega_weight(14) = 0.2611214147929032_DP
         casimir_omega(15) = 1.2235303456674129_DP ;  casimir_omega_weight(15) = 0.3634886598006601_DP
         casimir_omega(16) = 1.6614024348204408_DP ;  casimir_omega_weight(16) = 0.5252544830781409_DP
         casimir_omega(17) = 2.3100814197254005_DP ;  casimir_omega_weight(17) = 0.7969795763354163_DP
         casimir_omega(18) = 3.3278231490992036_DP ;  casimir_omega_weight(18) = 1.2909740037474824_DP
         casimir_omega(19) = 5.0554903207009545_DP ;  casimir_omega_weight(19) = 2.2906711549383716_DP
         casimir_omega(20) = 8.3427159387466574_DP ;  casimir_omega_weight(20) = 4.6514791032140801_DP
         casimir_omega(21) = 15.8286277582777650_DP ;  casimir_omega_weight(21) = 11.7616343525202254_DP
         casimir_omega(22) = 39.4808267754039193_DP ;  casimir_omega_weight(22) = 45.2153804140133104_DP
         casimir_omega(23) = 209.7265233924885308_DP ;  casimir_omega_weight(23) = 539.2518608062600833_DP
return
endsubroutine gauss_legendre_grid22

subroutine gauss_legendre_grid23()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0015733142406025_DP ;  casimir_omega_weight(2) = 0.0040446865811776_DP
         casimir_omega(3) = 0.0083519201789137_DP ;  casimir_omega_weight(3) = 0.0095570124043856_DP
         casimir_omega(4) = 0.0208059748945486_DP ;  casimir_omega_weight(4) = 0.0154281012447263_DP
         casimir_omega(5) = 0.0394024202031022_DP ;  casimir_omega_weight(5) = 0.0218837426915324_DP
         casimir_omega(6) = 0.0648610512005597_DP ;  casimir_omega_weight(6) = 0.0292046437481986_DP
         casimir_omega(7) = 0.0982178256531098_DP ;  casimir_omega_weight(7) = 0.0377476596256909_DP
         casimir_omega(8) = 0.1409191669902604_DP ;  casimir_omega_weight(8) = 0.0479847413841749_DP
         casimir_omega(9) = 0.1949644419774582_DP ;  casimir_omega_weight(9) = 0.0605618750729159_DP
         casimir_omega(10) = 0.2631224225367048_DP ;  casimir_omega_weight(10) = 0.0763909547879038_DP
         casimir_omega(11) = 0.3492651920709059_DP ;  casimir_omega_weight(11) = 0.0967979283979875_DP
         casimir_omega(12) = 0.4588950132666094_DP ;  casimir_omega_weight(12) = 0.1237701728022572_DP
         casimir_omega(13) = 0.6000000000000000_DP ;  casimir_omega_weight(13) = 0.1603854866233274_DP
         casimir_omega(14) = 0.7844931620358374_DP ;  casimir_omega_weight(14) = 0.2115883838793281_DP
         casimir_omega(15) = 1.0307354072859189_DP ;  casimir_omega_weight(15) = 0.2856656042938212_DP
         casimir_omega(16) = 1.3681844235444478_DP ;  casimir_omega_weight(16) = 0.3972178175956033_DP
         casimir_omega(17) = 1.8464905515520780_DP ;  casimir_omega_weight(17) = 0.5735760273627017_DP
         casimir_omega(18) = 2.5546560321697140_DP ;  casimir_omega_weight(18) = 0.8698923762277032_DP
         casimir_omega(19) = 3.6653224361885637_DP ;  casimir_omega_weight(19) = 1.4086785450566859_DP
         casimir_omega(20) = 5.5503263258381095_DP ;  casimir_omega_weight(20) = 2.4991161880976613_DP
         casimir_omega(21) = 9.1364946149083597_DP ;  casimir_omega_weight(21) = 5.0743252882594305_DP
         casimir_omega(22) = 17.3027220221401095_DP ;  casimir_omega_weight(22) = 12.8303599576521545_DP
         casimir_omega(23) = 43.1038602247302265_DP ;  casimir_omega_weight(23) = 49.3232835108618062_DP
         casimir_omega(24) = 228.8163360563904121_DP ;  casimir_omega_weight(24) = 588.2425392953313121_DP
return
endsubroutine gauss_legendre_grid23

subroutine gauss_legendre_grid24()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0014473168095932_DP ;  casimir_omega_weight(2) = 0.0037202521522486_DP
         casimir_omega(3) = 0.0076784560446872_DP ;  casimir_omega_weight(3) = 0.0087798954131766_DP
         casimir_omega(4) = 0.0191073389268272_DP ;  casimir_omega_weight(4) = 0.0141427267320948_DP
         casimir_omega(5) = 0.0361270795441245_DP ;  casimir_omega_weight(5) = 0.0199963555766831_DP
         casimir_omega(6) = 0.0593399398725749_DP ;  casimir_omega_weight(6) = 0.0265715449851933_DP
         casimir_omega(7) = 0.0896059521541503_DP ;  casimir_omega_weight(7) = 0.0341569002357493_DP
         casimir_omega(8) = 0.1281139628137400_DP ;  casimir_omega_weight(8) = 0.0431271023612803_DP
         casimir_omega(9) = 0.1764872057339697_DP ;  casimir_omega_weight(9) = 0.0539846913085114_DP
         casimir_omega(10) = 0.2369406010122470_DP ;  casimir_omega_weight(10) = 0.0674235046588170_DP
         casimir_omega(11) = 0.3125179117967915_DP ;  casimir_omega_weight(11) = 0.0844280475437189_DP
         casimir_omega(12) = 0.4074561261423225_DP ;  casimir_omega_weight(12) = 0.1064341433560470_DP
         casimir_omega(13) = 0.5277592467557539_DP ;  casimir_omega_weight(13) = 0.1355975265827953_DP
         casimir_omega(14) = 0.6821292136764917_DP ;  casimir_omega_weight(14) = 0.1752599026032904_DP
         casimir_omega(15) = 0.8835307089584750_DP ;  casimir_omega_weight(15) = 0.2307925396215768_DP
         casimir_omega(16) = 1.1519339737367846_DP ;  casimir_omega_weight(16) = 0.3111998788892227_DP
         casimir_omega(17) = 1.5193681389429419_DP ;  casimir_omega_weight(17) = 0.4323493920283529_DP
         casimir_omega(18) = 2.0398079198027004_DP ;  casimir_omega_weight(18) = 0.6239455173039228_DP
         casimir_omega(19) = 2.8099981617412797_DP ;  casimir_omega_weight(19) = 0.9459318539120899_DP
         casimir_omega(20) = 4.0175902531640659_DP ;  casimir_omega_weight(20) = 1.5314655574371665_DP
         casimir_omega(21) = 6.0667402220672155_DP ;  casimir_omega_weight(21) = 2.7165962936649906_DP
         casimir_omega(22) = 9.9648242964202893_DP ;  casimir_omega_weight(22) = 5.5155349506461260_DP
         casimir_omega(23) = 18.8409281574291363_DP ;  casimir_omega_weight(23) = 13.9455368081335429_DP
         casimir_omega(24) = 46.8844254502295996_DP ;  casimir_omega_weight(24) = 53.6097816493605421_DP
         casimir_omega(25) = 248.7361423662160576_DP ;  casimir_omega_weight(25) = 639.3632429654505813_DP
return
endsubroutine gauss_legendre_grid24

subroutine gauss_legendre_grid25()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0013358767331915_DP ;  casimir_omega_weight(2) = 0.0034333772047887_DP
         casimir_omega(3) = 0.0070834738127617_DP ;  casimir_omega_weight(3) = 0.0080942828231158_DP
         casimir_omega(4) = 0.0176097298283083_DP ;  casimir_omega_weight(4) = 0.0130132539478886_DP
         casimir_omega(5) = 0.0332480566374834_DP ;  casimir_omega_weight(5) = 0.0183474611322358_DP
         casimir_omega(6) = 0.0545064356941685_DP ;  casimir_omega_weight(6) = 0.0242884760371221_DP
         casimir_omega(7) = 0.0821052616930392_DP ;  casimir_omega_weight(7) = 0.0310723919565393_DP
         casimir_omega(8) = 0.1170316173923664_DP ;  casimir_omega_weight(8) = 0.0390006295224468_DP
         casimir_omega(9) = 0.1606187462467165_DP ;  casimir_omega_weight(9) = 0.0484701301312897_DP
         casimir_omega(10) = 0.2146624404574120_DP ;  casimir_omega_weight(10) = 0.0600181250342669_DP
         casimir_omega(11) = 0.2815930171944317_DP ;  casimir_omega_weight(11) = 0.0743904638654524_DP
         casimir_omega(12) = 0.3647334579808395_DP ;  casimir_omega_weight(12) = 0.0926489589359220_DP
         casimir_omega(13) = 0.4686951044920172_DP ;  casimir_omega_weight(13) = 0.1163451849938796_DP
         casimir_omega(14) = 0.6000000000000000_DP ;  casimir_omega_weight(14) = 0.1478112644720585_DP
         casimir_omega(15) = 0.7680899513345172_DP ;  casimir_omega_weight(15) = 0.1906646061021028_DP
         casimir_omega(16) = 0.9870221448642417_DP ;  casimir_omega_weight(16) = 0.2507216493782062_DP
         casimir_omega(17) = 1.2784407922708911_DP ;  casimir_omega_weight(17) = 0.3377349499255577_DP
         casimir_omega(18) = 1.6770516501764181_DP ;  casimir_omega_weight(18) = 0.4688919748360962_DP
         casimir_omega(19) = 2.2413323999368431_DP ;  casimir_omega_weight(19) = 0.6763698237660434_DP
         casimir_omega(20) = 3.0760918119506546_DP ;  casimir_omega_weight(20) = 1.0251034703955564_DP
         casimir_omega(21) = 4.3846154628421408_DP ;  casimir_omega_weight(21) = 1.6593393338113496_DP
         casimir_omega(22) = 6.6047246607709349_DP ;  casimir_omega_weight(22) = 2.9431147829042006_DP
         casimir_omega(23) = 10.8277005157089530_DP ;  casimir_omega_weight(23) = 5.9751105614844180_DP
         casimir_omega(24) = 20.4432437924905450_DP ;  casimir_omega_weight(24) = 15.1071666393551070_DP
         casimir_omega(25) = 50.8225214797034042_DP ;  casimir_omega_weight(25) = 58.0748759033261450_DP
         casimir_omega(26) = 269.4859421197828055_DP ;  casimir_omega_weight(26) = 692.6139723046327390_DP
return
endsubroutine gauss_legendre_grid25

subroutine gauss_legendre_grid26()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0012368340086875_DP ;  casimir_omega_weight(2) = 0.0031784755325699_DP
         casimir_omega(3) = 0.0065552101854690_DP ;  casimir_omega_weight(3) = 0.0074862938524137_DP
         casimir_omega(4) = 0.0162824704698250_DP ;  casimir_omega_weight(4) = 0.0120152234925946_DP
         casimir_omega(5) = 0.0307033114270580_DP ;  casimir_omega_weight(5) = 0.0168979197040197_DP
         casimir_omega(6) = 0.0502493365840785_DP ;  casimir_omega_weight(6) = 0.0222948316690473_DP
         casimir_omega(7) = 0.0755288741146963_DP ;  casimir_omega_weight(7) = 0.0284010913500546_DP
         casimir_omega(8) = 0.1073689677109616_DP ;  casimir_omega_weight(8) = 0.0354620724814886_DP
         casimir_omega(9) = 0.1468760541237260_DP ;  casimir_omega_weight(9) = 0.0437955865971495_DP
         casimir_omega(10) = 0.1955234754715632_DP ;  casimir_omega_weight(10) = 0.0538240964591008_DP
         casimir_omega(11) = 0.2552784996120011_DP ;  casimir_omega_weight(11) = 0.0661227447089209_DP
         casimir_omega(12) = 0.3287890906642715_DP ;  casimir_omega_weight(12) = 0.0814928979275316_DP
         casimir_omega(13) = 0.4196635137901181_DP ;  casimir_omega_weight(13) = 0.1010778826371524_DP
         casimir_omega(14) = 0.5328983262880461_DP ;  casimir_omega_weight(14) = 0.1265505334161409_DP
         casimir_omega(15) = 0.6755510052125963_DP ;  casimir_omega_weight(15) = 0.1604271130948436_DP
         casimir_omega(16) = 0.8578301142949564_DP ;  casimir_omega_weight(16) = 0.2066123185983827_DP
         casimir_omega(17) = 1.0949268397946874_DP ;  casimir_omega_weight(17) = 0.2713860153122147_DP
         casimir_omega(18) = 1.4102245216387814_DP ;  casimir_omega_weight(18) = 0.3652791604788864_DP
         casimir_omega(19) = 1.8412111340172961_DP ;  casimir_omega_weight(19) = 0.5068523124392346_DP
         casimir_omega(20) = 2.4510462385974905_DP ;  casimir_omega_weight(20) = 0.7308543821968978_DP
         casimir_omega(21) = 3.3529241053068879_DP ;  casimir_omega_weight(21) = 1.1074115750782592_DP
         casimir_omega(22) = 4.7663890693420417_DP ;  casimir_omega_weight(22) = 1.7923033138652953_DP
         casimir_omega(23) = 7.1642736894175600_DP ;  casimir_omega_weight(23) = 3.1786743227802452_DP
         casimir_omega(24) = 11.7251196456530575_DP ;  casimir_omega_weight(24) = 6.4530541196824407_DP
         casimir_omega(25) = 22.1096669984545962_DP ;  casimir_omega_weight(25) = 16.3152508598491117_DP
         casimir_omega(26) = 54.9181475214959320_DP ;  casimir_omega_weight(26) = 62.7185671464050500_DP
         casimir_omega(27) = 291.0657351523462353_DP ;  casimir_omega_weight(27) = 747.9947277105159174_DP
return
endsubroutine gauss_legendre_grid26

subroutine gauss_legendre_grid27()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0011484150292641_DP ;  casimir_omega_weight(2) = 0.0029509628992882_DP
         casimir_omega(3) = 0.0060840303018527_DP ;  casimir_omega_weight(3) = 0.0069445929830681_DP
         casimir_omega(4) = 0.0151005468695832_DP ;  casimir_omega_weight(4) = 0.0111288281620021_DP
         casimir_omega(5) = 0.0284425820576520_DP ;  casimir_omega_weight(5) = 0.0156163956055870_DP
         casimir_omega(6) = 0.0464793059098274_DP ;  casimir_omega_weight(6) = 0.0205427612955070_DP
         casimir_omega(7) = 0.0697282029882486_DP ;  casimir_omega_weight(7) = 0.0260707075516299_DP
         casimir_omega(8) = 0.0988879344007488_DP ;  casimir_omega_weight(8) = 0.0324021229052978_DP
         casimir_omega(9) = 0.1348852602939708_DP ;  casimir_omega_weight(9) = 0.0397945732391664_DP
         casimir_omega(10) = 0.1789417943485575_DP ;  casimir_omega_weight(10) = 0.0485848883186924_DP
         casimir_omega(11) = 0.2326693652390074_DP ;  casimir_omega_weight(11) = 0.0592235903516645_DP
         casimir_omega(12) = 0.2982076914661539_DP ;  casimir_omega_weight(12) = 0.0723264106721176_DP
         casimir_omega(13) = 0.3784262192575336_DP ;  casimir_omega_weight(13) = 0.0887533289999953_DP
         casimir_omega(14) = 0.4772257912822858_DP ;  casimir_omega_weight(14) = 0.1097330834781429_DP
         casimir_omega(15) = 0.5999999999999996_DP ;  casimir_omega_weight(15) = 0.1370650408547479_DP
         casimir_omega(16) = 0.7543598996036966_DP ;  casimir_omega_weight(16) = 0.1734571755087970_DP
         casimir_omega(17) = 0.9513082912339277_DP ;  casimir_omega_weight(17) = 0.2231129172760864_DP
         casimir_omega(18) = 1.2072123231632343_DP ;  casimir_omega_weight(18) = 0.2927937030204167_DP
         casimir_omega(19) = 1.5472599911474980_DP ;  casimir_omega_weight(19) = 0.3938390934668549_DP
         casimir_omega(20) = 2.0118273727531877_DP ;  casimir_omega_weight(20) = 0.5462357666499511_DP
         casimir_omega(21) = 2.6689350579552626_DP ;  casimir_omega_weight(21) = 0.7874035413721758_DP
         casimir_omega(22) = 3.6404845766226668_DP ;  casimir_omega_weight(22) = 1.1928596688907473_DP
         casimir_omega(23) = 5.1629037401218003_DP ;  casimir_omega_weight(23) = 1.9303602811708309_DP
         casimir_omega(24) = 7.7453824439293477_DP ;  casimir_omega_weight(24) = 3.4232770815626288_DP
         casimir_omega(25) = 12.6570787163518954_DP ;  casimir_omega_weight(25) = 6.9493672566352576_DP
         casimir_omega(26) = 23.8401961934995654_DP ;  casimir_omega_weight(26) = 17.5697906226352067_DP
         casimir_omega(27) = 59.1713029256882308_DP ;  casimir_omega_weight(27) = 67.5408560952761690_DP
         casimir_omega(28) = 313.4755213284692559_DP ;  casimir_omega_weight(28) = 805.5055095091277053_DP
return
endsubroutine gauss_legendre_grid27

subroutine gauss_legendre_grid28()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0010691524840916_DP ;  casimir_omega_weight(2) = 0.0027470487188693_DP
         casimir_omega(3) = 0.0056619809494810_DP ;  casimir_omega_weight(3) = 0.0064598479416396_DP
         casimir_omega(4) = 0.0140433932665294_DP ;  casimir_omega_weight(4) = 0.0103378832885768_DP
         casimir_omega(5) = 0.0264247815042413_DP ;  casimir_omega_weight(5) = 0.0144775450730216_DP
         casimir_omega(6) = 0.0431238592158273_DP ;  casimir_omega_weight(6) = 0.0189940462236321_DP
         casimir_omega(7) = 0.0645837979390434_DP ;  casimir_omega_weight(7) = 0.0240243319693881_DP
         casimir_omega(8) = 0.0913992158042783_DP ;  casimir_omega_weight(8) = 0.0297361203702559_DP
         casimir_omega(9) = 0.1243528853588817_DP ;  casimir_omega_weight(9) = 0.0363404483891758_DP
         casimir_omega(10) = 0.1644673185763083_DP ;  casimir_omega_weight(10) = 0.0441092155057525_DP
         casimir_omega(11) = 0.2130774152089362_DP ;  casimir_omega_weight(11) = 0.0534003016834685_DP
         casimir_omega(12) = 0.2719336443748769_DP ;  casimir_omega_weight(12) = 0.0646943752049677_DP
         casimir_omega(13) = 0.3433505216495335_DP ;  casimir_omega_weight(13) = 0.0786500904781250_DP
         casimir_omega(14) = 0.4304238815244323_DP ;  casimir_omega_weight(14) = 0.0961888746427630_DP
         casimir_omega(15) = 0.5373552788927303_DP ;  casimir_omega_weight(15) = 0.1186285735875089_DP
         casimir_omega(16) = 0.6699478243552620_DP ;  casimir_omega_weight(16) = 0.1479002029068823_DP
         casimir_omega(17) = 0.8363848184375556_DP ;  casimir_omega_weight(17) = 0.1869108985516019_DP
         casimir_omega(18) = 1.0484911986458583_DP ;  casimir_omega_weight(18) = 0.2401741731535455_DP
         casimir_omega(19) = 1.3238523715135386_DP ;  casimir_omega_weight(19) = 0.3149511059419137_DP
         casimir_omega(20) = 1.6895267837137813_DP ;  casimir_omega_weight(20) = 0.4234200037772592_DP
         casimir_omega(21) = 2.1888847165278569_DP ;  casimir_omega_weight(21) = 0.5870466455849690_DP
         casimir_omega(22) = 2.8949871083492895_DP ;  casimir_omega_weight(22) = 0.8460208164424621_DP
         casimir_omega(23) = 3.9387646472908604_DP ;  casimir_omega_weight(23) = 1.2814505970462298_DP
         casimir_omega(24) = 5.5741534485132256_DP ;  casimir_omega_weight(24) = 2.0735125088460089_DP
         casimir_omega(25) = 8.3480469175605201_DP ;  casimir_omega_weight(25) = 3.6769248372603895_DP
         casimir_omega(26) = 13.6235752769504881_DP ;  casimir_omega_weight(26) = 7.4640513147136271_DP
         casimir_omega(27) = 25.6348300704512653_DP ;  casimir_omega_weight(27) = 18.8707868790113515_DP
         casimir_omega(28) = 63.5819871546901396_DP ;  casimir_omega_weight(28) = 72.5417433423614995_DP
         casimir_omega(29) = 336.7153005362702629_DP ;  casimir_omega_weight(29) = 865.1463179714381795_DP
return
endsubroutine gauss_legendre_grid28

subroutine gauss_legendre_grid29()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0009978239879040_DP ;  casimir_omega_weight(2) = 0.0025635766011823_DP
         casimir_omega(3) = 0.0052824496665318_DP ;  casimir_omega_weight(3) = 0.0060243181042260_DP
         casimir_omega(4) = 0.0130939718702417_DP ;  casimir_omega_weight(4) = 0.0096290549024019_DP
         casimir_omega(5) = 0.0246160456041505_DP ;  casimir_omega_weight(5) = 0.0134606801208471_DP
         casimir_omega(6) = 0.0401236531401030_DP ;  casimir_omega_weight(6) = 0.0176178400933112_DP
         casimir_omega(7) = 0.0599986679983517_DP ;  casimir_omega_weight(7) = 0.0222166322337308_DP
         casimir_omega(8) = 0.0847506061611039_DP ;  casimir_omega_weight(8) = 0.0273976120244943_DP
         casimir_omega(9) = 0.1150456470924083_DP ;  casimir_omega_weight(9) = 0.0333354133221457_DP
         casimir_omega(10) = 0.1517469673282803_DP ;  casimir_omega_weight(10) = 0.0402519845687114_DP
         casimir_omega(11) = 0.1959708195247690_DP ;  casimir_omega_weight(11) = 0.0484352075825952_DP
         casimir_omega(12) = 0.2491650196112665_DP ;  casimir_omega_weight(12) = 0.0582656673381993_DP
         casimir_omega(13) = 0.3132200210546787_DP ;  casimir_omega_weight(13) = 0.0702559759529726_DP
         casimir_omega(14) = 0.3906284189202038_DP ;  casimir_omega_weight(14) = 0.0851098229569565_DP
         casimir_omega(15) = 0.4847180820471179_DP ;  casimir_omega_weight(15) = 0.1038127450794328_DP
         casimir_omega(16) = 0.5999999999999994_DP ;  casimir_omega_weight(16) = 0.1277752580619767_DP
         casimir_omega(17) = 0.7426997533898587_DP ;  casimir_omega_weight(17) = 0.1590650380600494_DP
         casimir_omega(18) = 0.9215919338258365_DP ;  casimir_omega_weight(18) = 0.2007957499438864_DP
         casimir_omega(19) = 1.1493518159784393_DP ;  casimir_omega_weight(19) = 0.2578022735359866_DP
         casimir_omega(20) = 1.4448256041785152_DP ;  casimir_omega_weight(20) = 0.3378633491415308_DP
         casimir_omega(21) = 1.8370081876118236_DP ;  casimir_omega_weight(21) = 0.4540261305926687_DP
         casimir_omega(22) = 2.3723703105130127_DP ;  casimir_omega_weight(22) = 0.6292884451750238_DP
         casimir_omega(23) = 3.1291927082720190_DP ;  casimir_omega_weight(23) = 0.9067090753212533_DP
         casimir_omega(24) = 4.2477572292010475_DP ;  casimir_omega_weight(24) = 1.3731866922422622_DP
         casimir_omega(25) = 6.0001332031219503_DP ;  casimir_omega_weight(25) = 2.2217618686271834_DP
         casimir_omega(26) = 8.9722637852281224_DP ;  casimir_omega_weight(26) = 3.9396190593912572_DP
         casimir_omega(27) = 14.6246072902669102_DP ;  casimir_omega_weight(27) = 7.9971074068088894_DP
         casimir_omega(28) = 27.4935675414246887_DP ;  casimir_omega_weight(28) = 20.2182404195425960_DP
         casimir_omega(29) = 68.1501997607038561_DP ;  casimir_omega_weight(29) = 77.7212293807953358_DP
         casimir_omega(30) = 360.7850726822511547_DP ;  casimir_omega_weight(30) = 926.9171533217428305_DP
return
endsubroutine gauss_legendre_grid29

subroutine gauss_legendre_grid30()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0009334045957258_DP ;  casimir_omega_weight(2) = 0.0023979010815464_DP
         casimir_omega(3) = 0.0049399019509120_DP ;  casimir_omega_weight(3) = 0.0056315385615115_DP
         casimir_omega(4) = 0.0122380680785943_DP ;  casimir_omega_weight(4) = 0.0089912741585551_DP
         casimir_omega(5) = 0.0229882517143179_DP ;  casimir_omega_weight(5) = 0.0125487698845914_DP
         casimir_omega(6) = 0.0374297012951682_DP ;  casimir_omega_weight(6) = 0.0163890088513474_DP
         casimir_omega(7) = 0.0558933407439577_DP ;  casimir_omega_weight(7) = 0.0206111084262841_DP
         casimir_omega(8) = 0.0788184859737716_DP ;  casimir_omega_weight(8) = 0.0253338039215627_DP
         casimir_omega(9) = 0.1067760110091904_DP ;  casimir_omega_weight(9) = 0.0307029115773563_DP
         casimir_omega(10) = 0.1405002230735112_DP ;  casimir_omega_weight(10) = 0.0369014453821285_DP
         casimir_omega(11) = 0.1809326707841716_DP ;  casimir_omega_weight(11) = 0.0441636186647028_DP
         casimir_omega(12) = 0.2292826506620191_DP ;  casimir_omega_weight(12) = 0.0527946279038177_DP
         casimir_omega(13) = 0.2871115605285195_DP ;  casimir_omega_weight(13) = 0.0631991736316367_DP
         casimir_omega(14) = 0.3564520021463336_DP ;  casimir_omega_weight(14) = 0.0759234201816370_DP
         casimir_omega(15) = 0.4399785849751038_DP ;  casimir_omega_weight(15) = 0.0917180566423294_DP
         casimir_omega(16) = 0.5412573798302813_DP ;  casimir_omega_weight(16) = 0.1116352758097951_DP
         casimir_omega(17) = 0.6651179520413798_DP ;  casimir_omega_weight(17) = 0.1371817342157401_DP
         casimir_omega(18) = 0.8182216414473230_DP ;  casimir_omega_weight(18) = 0.1705667080603292_DP
         casimir_omega(19) = 1.0099536482676561_DP ;  casimir_omega_weight(19) = 0.2151177009518480_DP
         casimir_omega(20) = 1.2538680063502361_DP ;  casimir_omega_weight(20) = 0.2760021982347564_DP
         casimir_omega(21) = 1.5701144371829014_DP ;  casimir_omega_weight(21) = 0.3615345829182477_DP
         casimir_omega(22) = 1.9896904104700461_DP ;  casimir_omega_weight(22) = 0.4856609266196846_DP
         casimir_omega(23) = 2.5622735119192241_DP ;  casimir_omega_weight(23) = 0.6729640279979591_DP
         casimir_omega(24) = 3.3715438196039571_DP ;  casimir_omega_weight(24) = 0.9694706778620220_DP
         casimir_omega(25) = 4.5674564228472629_DP ;  casimir_omega_weight(25) = 1.4680698824283427_DP
         casimir_omega(26) = 6.4408388406970856_DP ;  casimir_omega_weight(26) = 2.3751099135397586_DP
         casimir_omega(27) = 9.6180302685576766_DP ;  casimir_omega_weight(27) = 4.2113609713542610_DP
         casimir_omega(28) = 15.6601730516018058_DP ;  casimir_omega_weight(28) = 8.5485364619978643_DP
         casimir_omega(29) = 29.4164076950739783_DP ;  casimir_omega_weight(29) = 21.6121519056479769_DP
         casimir_omega(30) = 72.8759403683179272_DP ;  casimir_omega_weight(30) = 83.0793146238104612_DP
         casimir_omega(31) = 385.6848376882703633_DP ;  casimir_omega_weight(31) = 990.8180157497322398_DP
return
endsubroutine gauss_legendre_grid30

subroutine gauss_legendre_grid31()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0008750297241470_DP ;  casimir_omega_weight(2) = 0.0022477914403405_DP
         casimir_omega(3) = 0.0046296767443968_DP ;  casimir_omega_weight(3) = 0.0052760751910971_DP
         casimir_omega(4) = 0.0114637451964151_DP ;  casimir_omega_weight(4) = 0.0084152881498807_DP
         casimir_omega(5) = 0.0215178820054885_DP ;  casimir_omega_weight(5) = 0.0117276849149515_DP
         casimir_omega(6) = 0.0350012599392569_DP ;  casimir_omega_weight(6) = 0.0152868941777078_DP
         casimir_omega(7) = 0.0522021561806568_DP ;  casimir_omega_weight(7) = 0.0191780843046211_DP
         casimir_omega(8) = 0.0735015291523344_DP ;  casimir_omega_weight(8) = 0.0235022926192878_DP
         casimir_omega(9) = 0.0993916754089614_DP ;  casimir_omega_weight(9) = 0.0283822767173601_DP
         casimir_omega(10) = 0.1305016673115107_DP ;  casimir_omega_weight(10) = 0.0339703393175618_DP
         casimir_omega(11) = 0.1676319486522720_DP ;  casimir_omega_weight(11) = 0.0404589968014489_DP
         casimir_omega(12) = 0.2118015476563194_DP ;  casimir_omega_weight(12) = 0.0480958142485939_DP
         casimir_omega(13) = 0.2643130116849958_DP ;  casimir_omega_weight(13) = 0.0572044278118154_DP
         casimir_omega(14) = 0.3268427063947287_DP ;  casimir_omega_weight(14) = 0.0682149022876572_DP
         casimir_omega(15) = 0.4015681223781624_DP ;  casimir_omega_weight(15) = 0.0817084356341859_DP
         casimir_omega(16) = 0.4913502820072247_DP ;  casimir_omega_weight(16) = 0.0984845808006035_DP
         casimir_omega(17) = 0.6000000000000004_DP ;  casimir_omega_weight(17) = 0.1196646537521120_DP
         casimir_omega(18) = 0.7326748618711622_DP ;  casimir_omega_weight(18) = 0.1468548595103048_DP
         casimir_omega(19) = 0.8964855025543649_DP ;  casimir_omega_weight(19) = 0.1824109631726761_DP
         casimir_omega(20) = 1.1014472495685041_DP ;  casimir_omega_weight(20) = 0.2298815761658237_DP
         casimir_omega(21) = 1.3620214824272161_DP ;  casimir_omega_weight(21) = 0.2947779947455104_DP
         casimir_omega(22) = 1.6997042938711424_DP ;  casimir_omega_weight(22) = 0.3859681994780089_DP
         casimir_omega(23) = 2.1475619826311725_DP ;  casimir_omega_weight(23) = 0.5183272286980609_DP
         casimir_omega(24) = 2.7585854450477694_DP ;  casimir_omega_weight(24) = 0.7180757574619490_DP
         casimir_omega(25) = 3.6220337218255727_DP ;  casimir_omega_weight(25) = 1.0343075810873448_DP
         casimir_omega(26) = 4.8978572847632504_DP ;  casimir_omega_weight(26) = 1.5661017728684810_DP
         casimir_omega(27) = 6.8962668659536313_DP ;  casimir_omega_weight(27) = 2.5335579412604927_DP
         casimir_omega(28) = 10.2853440311795623_DP ;  casimir_omega_weight(28) = 4.4921515985088840_DP
         casimir_omega(29) = 16.7302711255771683_DP ;  casimir_omega_weight(29) = 9.1183392609198730_DP
         casimir_omega(30) = 31.4033497632671619_DP ;  casimir_omega_weight(30) = 23.0525218941552801_DP
         casimir_omega(31) = 77.7592086608892714_DP ;  casimir_omega_weight(31) = 88.6159994197483911_DP
         casimir_omega(32) = 411.4145954881013267_DP ;  casimir_omega_weight(32) = 1056.8489054138153733_DP
return
endsubroutine gauss_legendre_grid31

subroutine gauss_legendre_grid32()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0008219659513546_DP ;  casimir_omega_weight(2) = 0.0021113560129362_DP
         casimir_omega(3) = 0.0043478258574964_DP ;  casimir_omega_weight(3) = 0.0049533330238981_DP
         casimir_omega(4) = 0.0107609185729061_DP ;  casimir_omega_weight(4) = 0.0078933118185897_DP
         casimir_omega(5) = 0.0201851422779929_DP ;  casimir_omega_weight(5) = 0.0109856187912516_DP
         casimir_omega(6) = 0.0328042042621407_DP ;  casimir_omega_weight(6) = 0.0142943806113153_DP
         casimir_omega(7) = 0.0488704523045288_DP ;  casimir_omega_weight(7) = 0.0178932157860400_DP
         casimir_omega(8) = 0.0687159854529204_DP ;  casimir_omega_weight(8) = 0.0218686797223171_DP
         casimir_omega(9) = 0.0927678025406974_DP ;  casimir_omega_weight(9) = 0.0263248974445517_DP
         casimir_omega(10) = 0.1215682852598141_DP ;  casimir_omega_weight(10) = 0.0313896819458937_DP
         casimir_omega(11) = 0.1558027905986396_DP ;  casimir_omega_weight(11) = 0.0372227612997642_DP
         casimir_omega(12) = 0.1963368986134621_DP ;  casimir_omega_weight(12) = 0.0440270511419815_DP
         casimir_omega(13) = 0.2442670155872915_DP ;  casimir_omega_weight(13) = 0.0520643786047483_DP
         casimir_omega(14) = 0.3009897807815471_DP ;  casimir_omega_weight(14) = 0.0616778057822990_DP
         casimir_omega(15) = 0.3682984242002268_DP ;  casimir_omega_weight(15) = 0.0733239009109432_DP
         casimir_omega(16) = 0.4485184786356920_DP ;  casimir_omega_weight(16) = 0.0876202900921237_DP
         casimir_omega(17) = 0.5447021130125425_DP ;  casimir_omega_weight(17) = 0.1054171851740935_DP
         casimir_omega(18) = 0.6609117009092826_DP ;  casimir_omega_weight(18) = 0.1279074369165793_DP
         casimir_omega(19) = 0.8026425156329157_DP ;  casimir_omega_weight(19) = 0.1568001618884276_DP
         casimir_omega(20) = 0.9774682060661882_DP ;  casimir_omega_weight(20) = 0.1946024668469130_DP
         casimir_omega(21) = 1.1960538961330420_DP ;  casimir_omega_weight(21) = 0.2450913108056526_DP
         casimir_omega(22) = 1.4737970214048397_DP ;  casimir_omega_weight(22) = 0.3141329824024104_DP
         casimir_omega(23) = 1.8335830021882378_DP ;  casimir_omega_weight(23) = 0.4111669949994562_DP
         casimir_omega(24) = 2.3106132991378132_DP ;  casimir_omega_weight(24) = 0.5520273864120199_DP
         casimir_omega(25) = 2.9612986580390817_DP ;  casimir_omega_weight(25) = 0.7646255997112130_DP
         casimir_omega(26) = 3.8806567595698636_DP ;  casimir_omega_weight(26) = 1.1012214196662327_DP
         casimir_omega(27) = 5.2389556465962093_DP ;  casimir_omega_weight(27) = 1.6672837093099935_DP
         casimir_omega(28) = 7.3664143265283153_DP ;  casimir_omega_weight(28) = 2.6971070431802300_DP
         casimir_omega(29) = 10.9742030967498128_DP ;  casimir_omega_weight(29) = 4.7819918056010033_DP
         casimir_omega(30) = 17.8349002965660084_DP ;  casimir_omega_weight(30) = 9.7065164634321004_DP
         casimir_omega(31) = 33.4543930948808850_DP ;  casimir_omega_weight(31) = 24.5393508565754921_DP
         casimir_omega(32) = 82.8000043698393426_DP ;  casimir_omega_weight(32) = 94.3312840639365078_DP
         casimir_omega(33) = 437.9743460258503092_DP ;  casimir_omega_weight(33) = 1125.0098224501393815_DP
return
endsubroutine gauss_legendre_grid32

subroutine gauss_legendre_grid33()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0007735878386695_DP ;  casimir_omega_weight(2) = 0.0019869821463144_DP
         casimir_omega(3) = 0.0040909868537744_DP ;  casimir_omega_weight(3) = 0.0046594050222253_DP
         casimir_omega(4) = 0.0101210200918098_DP ;  casimir_omega_weight(4) = 0.0074187556652948_DP
         casimir_omega(5) = 0.0189732724007222_DP ;  casimir_omega_weight(5) = 0.0103126407982716_DP
         casimir_omega(6) = 0.0308097685662858_DP ;  casimir_omega_weight(6) = 0.0133971835654812_DP
         casimir_omega(7) = 0.0458524035498812_DP ;  casimir_omega_weight(7) = 0.0167363691465067_DP
         casimir_omega(8) = 0.0643920995427091_DP ;  casimir_omega_weight(8) = 0.0204048058741780_DP
         casimir_omega(9) = 0.0868011990624277_DP ;  casimir_omega_weight(9) = 0.0244914262875677_DP
         casimir_omega(10) = 0.1135500951183624_DP ;  casimir_omega_weight(10) = 0.0291043193115246_DP
         casimir_omega(11) = 0.1452294389741103_DP ;  casimir_omega_weight(11) = 0.0343771449655216_DP
         casimir_omega(12) = 0.1825798228189066_DP ;  casimir_omega_weight(12) = 0.0404778042393456_DP
         casimir_omega(13) = 0.2265316566613443_DP ;  casimir_omega_weight(13) = 0.0476203579979834_DP
         casimir_omega(14) = 0.2782591828458295_DP ;  casimir_omega_weight(14) = 0.0560816864144829_DP
         casimir_omega(15) = 0.3392544266274640_DP ;  casimir_omega_weight(15) = 0.0662251690441571_DP
         casimir_omega(16) = 0.4114297471311042_DP ;  casimir_omega_weight(16) = 0.0785349411532497_DP
         casimir_omega(17) = 0.4972621732010999_DP ;  casimir_omega_weight(17) = 0.0936663909295258_DP
         casimir_omega(18) = 0.6000000000000000_DP ;  casimir_omega_weight(18) = 0.1125221353922521_DP
         casimir_omega(19) = 0.7239641770507454_DP ;  casimir_omega_weight(19) = 0.1363689322879254_DP
         casimir_omega(20) = 0.8749974995981120_DP ;  casimir_omega_weight(20) = 0.1670221407648512_DP
         casimir_omega(21) = 1.0611504869037920_DP ;  casimir_omega_weight(21) = 0.2071450358808810_DP
         casimir_omega(22) = 1.2937578423043783_DP ;  casimir_omega_weight(22) = 0.2607501426056826_DP
         casimir_omega(23) = 1.5891818622868481_DP ;  casimir_omega_weight(23) = 0.3340699058195492_DP
         casimir_omega(24) = 1.9717403294725995_DP ;  casimir_omega_weight(24) = 0.4371332923593482_DP
         casimir_omega(25) = 2.4788362644861306_DP ;  casimir_omega_weight(25) = 0.5867633601836164_DP
         casimir_omega(26) = 3.1704068554477449_DP ;  casimir_omega_weight(26) = 0.8126152018826227_DP
         casimir_omega(27) = 4.1474081451465565_DP ;  casimir_omega_weight(27) = 1.1702135681128156_DP
         casimir_omega(28) = 5.5907479730681002_DP ;  casimir_omega_weight(28) = 1.7716168270961128_DP
         casimir_omega(29) = 7.8512787145033460_DP ;  casimir_omega_weight(29) = 2.8657581427566710_DP
         casimir_omega(30) = 11.6846057842167550_DP ;  casimir_omega_weight(30) = 5.0808823261572407_DP
         casimir_omega(31) = 18.9740595294619041_DP ;  casimir_omega_weight(31) = 10.3130686304234551_DP
         casimir_omega(32) = 35.5695371350285612_DP ;  casimir_omega_weight(32) = 26.0726391943385210_DP
         casimir_omega(33) = 87.9983272661618940_DP ;  casimir_omega_weight(33) = 100.2251688081349812_DP
         casimir_omega(34) = 465.3640892535753437_DP ;  casimir_omega_weight(34) = 1195.3007669731889564_DP
return
endsubroutine gauss_legendre_grid33

subroutine gauss_legendre_grid34()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0007293593948402_DP ;  casimir_omega_weight(2) = 0.0018732882103634_DP
         casimir_omega(3) = 0.0038562816466974_DP ;  casimir_omega_weight(3) = 0.0043909517945136_DP
         casimir_omega(4) = 0.0095367316921978_DP ;  casimir_omega_weight(4) = 0.0069860108862205_DP
         casimir_omega(5) = 0.0178680019772297_DP ;  casimir_omega_weight(5) = 0.0097003466125052_DP
         casimir_omega(6) = 0.0289935594808174_DP ;  casimir_omega_weight(6) = 0.0125833000878572_DP
         casimir_omega(7) = 0.0431093433535578_DP ;  casimir_omega_weight(7) = 0.0156907672600391_DP
         casimir_omega(8) = 0.0604713616708387_DP ;  casimir_omega_weight(8) = 0.0190874261926407_DP
         casimir_omega(9) = 0.0814059015795173_DP ;  casimir_omega_weight(9) = 0.0228497184854199_DP
         casimir_omega(10) = 0.1063231346344596_DP ;  casimir_omega_weight(10) = 0.0270697011135564_DP
         casimir_omega(11) = 0.1357351407868475_DP ;  casimir_omega_weight(11) = 0.0318600990382899_DP
         casimir_omega(12) = 0.1702797813330445_DP ;  casimir_omega_weight(12) = 0.0373610479741160_DP
         casimir_omega(13) = 0.2107524447476123_DP ;  casimir_omega_weight(13) = 0.0437492400240230_DP
         casimir_omega(14) = 0.2581485595915055_DP ;  casimir_omega_weight(14) = 0.0512505243497326_DP
         casimir_omega(15) = 0.3137210650625930_DP ;  casimir_omega_weight(15) = 0.0601575429733735_DP
         casimir_omega(16) = 0.3790589955719529_DP ;  casimir_omega_weight(16) = 0.0708548185198764_DP
         casimir_omega(17) = 0.4561963739571342_DP ;  casimir_omega_weight(17) = 0.0838550637821114_DP
         casimir_omega(18) = 0.5477654009584496_DP ;  casimir_omega_weight(18) = 0.0998527200251564_DP
         casimir_omega(19) = 0.6572156608834572_DP ;  casimir_omega_weight(19) = 0.1198045208177016_DP
         casimir_omega(20) = 0.7891338479464246_DP ;  casimir_omega_weight(20) = 0.1450534746214089_DP
         casimir_omega(21) = 0.9497202393437578_DP ;  casimir_omega_weight(21) = 0.1775244908825332_DP
         casimir_omega(22) = 1.1475161858454535_DP ;  casimir_omega_weight(22) = 0.2200418204268994_DP
         casimir_omega(23) = 1.3945458404635855_DP ;  casimir_omega_weight(23) = 0.2768607567154135_DP
         casimir_omega(24) = 1.7081652382781136_DP ;  casimir_omega_weight(24) = 0.3545910515990300_DP
         casimir_omega(25) = 2.1141676197944377_DP ;  casimir_omega_weight(25) = 0.4638690351262151_DP
         casimir_omega(26) = 2.6522240144527394_DP ;  casimir_omega_weight(26) = 0.6225367968998596_DP
         casimir_omega(27) = 3.3859046879842754_DP ;  casimir_omega_weight(27) = 0.8620459528195539_DP
         casimir_omega(28) = 4.4222838027087281_DP ;  casimir_omega_weight(28) = 1.2412851893277719_DP
         casimir_omega(29) = 5.9532312495222737_DP ;  casimir_omega_weight(29) = 1.8791020897049644_DP
         casimir_omega(30) = 8.3508578882189965_DP ;  casimir_omega_weight(30) = 3.0395120257588260_DP
         casimir_omega(31) = 12.4165506562995400_DP ;  casimir_omega_weight(31) = 5.3888237857676309_DP
         casimir_omega(32) = 20.1477479383967939_DP ;  casimir_omega_weight(32) = 10.9379962411522005_DP
         casimir_omega(33) = 37.7487814084698030_DP ;  casimir_omega_weight(33) = 27.6523872509572186_DP
         casimir_omega(34) = 93.3541771536087452_DP ;  casimir_omega_weight(34) = 106.2976538679523344_DP
         casimir_omega(35) = 493.5838251303969741_DP ;  casimir_omega_weight(35) = 1267.7217390824391714_DP
return
endsubroutine gauss_legendre_grid34

subroutine gauss_legendre_grid35()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0006888191496813_DP ;  casimir_omega_weight(2) = 0.0017690849740083_DP
         casimir_omega(3) = 0.0036412350253752_DP ;  casimir_omega_weight(3) = 0.0041451052044580_DP
         casimir_omega(4) = 0.0090017721099520_DP ;  casimir_omega_weight(4) = 0.0065902784473131_DP
         casimir_omega(5) = 0.0168571171575977_DP ;  casimir_omega_weight(5) = 0.0091415830743158_DP
         casimir_omega(6) = 0.0273347762245373_DP ;  casimir_omega_weight(6) = 0.0118425809599055_DP
         casimir_omega(7) = 0.0406084498970949_DP ;  casimir_omega_weight(7) = 0.0147423327228230_DP
         casimir_omega(8) = 0.0569043745924503_DP ;  casimir_omega_weight(8) = 0.0178972048410172_DP
         casimir_omega(9) = 0.0765097892815093_DP ;  casimir_omega_weight(9) = 0.0213732899097384_DP
         casimir_omega(10) = 0.0997841451134778_DP ;  casimir_omega_weight(10) = 0.0252495034317188_DP
         casimir_omega(11) = 0.1271738472016536_DP ;  casimir_omega_weight(11) = 0.0296216010643533_DP
         casimir_omega(12) = 0.1592316183690299_DP ;  casimir_omega_weight(12) = 0.0346074766610985_DP
         casimir_omega(13) = 0.1966420088157427_DP ;  casimir_omega_weight(13) = 0.0403542590391098_DP
         casimir_omega(14) = 0.2402552030980493_DP ;  casimir_omega_weight(14) = 0.0470479605998858_DP
         casimir_omega(15) = 0.2911321957497814_DP ;  casimir_omega_weight(15) = 0.0549267904041556_DP
         casimir_omega(16) = 0.3506057781672661_DP ;  casimir_omega_weight(16) = 0.0642998035201926_DP
         casimir_omega(17) = 0.4203638597207389_DP ;  casimir_omega_weight(17) = 0.0755734446162303_DP
         casimir_omega(18) = 0.5025648620886812_DP ;  casimir_omega_weight(18) = 0.0892899766046728_DP
         casimir_omega(19) = 0.6000000000000000_DP ;  casimir_omega_weight(19) = 0.1061841538885251_DP
         casimir_omega(20) = 0.7163254480303786_DP ;  casimir_omega_weight(20) = 0.1272685126256954_DP
         casimir_omega(21) = 0.8564009290407586_DP ;  casimir_omega_weight(21) = 0.1539646348835655_DP
         casimir_omega(22) = 1.0267942584455982_DP ;  casimir_omega_weight(22) = 0.1883102709226199_DP
         casimir_omega(23) = 1.2365516602272604_DP ;  casimir_omega_weight(23) = 0.2332954405482115_DP
         casimir_omega(24) = 1.4984066748934557_DP ;  casimir_omega_weight(24) = 0.2934253963866209_DP
         casimir_omega(25) = 1.8307380104997133_DP ;  casimir_omega_weight(25) = 0.3756983380782848_DP
         casimir_omega(26) = 2.2608575086241727_DP ;  casimir_omega_weight(26) = 0.4913758603046289_DP
         casimir_omega(27) = 2.8307706963458066_DP ;  casimir_omega_weight(27) = 0.6593490888017024_DP
         casimir_omega(28) = 3.6077875857992900_DP ;  casimir_omega_weight(28) = 0.9129190306231062_DP
         casimir_omega(29) = 4.7052802442759267_DP ;  casimir_omega_weight(29) = 1.3144372728233182_DP
         casimir_omega(30) = 6.3264028921910436_DP ;  casimir_omega_weight(30) = 1.9897403192507632_DP
         casimir_omega(31) = 8.8651500097213596_DP ;  casimir_omega_weight(31) = 3.2183693643130082_DP
         casimir_omega(32) = 13.1700364781784227_DP ;  casimir_omega_weight(32) = 5.7058167206699606_DP
         casimir_omega(33) = 21.3559647616106822_DP ;  casimir_omega_weight(33) = 11.5812997071348818_DP
         casimir_omega(34) = 39.9921255062654524_DP ;  casimir_omega_weight(34) = 29.2785953217814487_DP
         casimir_omega(35) = 98.8675538632391664_DP ;  casimir_omega_weight(35) = 112.5487394289544909_DP
         casimir_omega(36) = 522.6335536208279109_DP ;  casimir_omega_weight(36) = 1342.2727388618263831_DP
return
endsubroutine gauss_legendre_grid35

subroutine gauss_legendre_grid36()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0006515680554444_DP ;  casimir_omega_weight(2) = 0.0016733443134754_DP
         casimir_omega(3) = 0.0034437087505760_DP ;  casimir_omega_weight(3) = 0.0039193905891516_DP
         casimir_omega(4) = 0.0085107249996309_DP ;  casimir_omega_weight(4) = 0.0062274320797591_DP
         casimir_omega(5) = 0.0159301132894687_DP ;  casimir_omega_weight(5) = 0.0086302295247614_DP
         casimir_omega(6) = 0.0258155894356164_DP ;  casimir_omega_weight(6) = 0.0111663942658118_DP
         casimir_omega(7) = 0.0383217072422813_DP ;  casimir_omega_weight(7) = 0.0138791773512148_DP
         casimir_omega(8) = 0.0536491825547036_DP ;  casimir_omega_weight(8) = 0.0168179433813631_DP
         casimir_omega(9) = 0.0720519570935094_DP ;  casimir_omega_weight(9) = 0.0200401492973437_DP
         casimir_omega(10) = 0.0938464950863591_DP ;  casimir_omega_weight(10) = 0.0236138539955507_DP
         casimir_omega(11) = 0.1194239270808082_DP ;  casimir_omega_weight(11) = 0.0276209402559650_DP
         casimir_omega(12) = 0.1492658795299639_DP ;  casimir_omega_weight(12) = 0.0321613172354544_DP
         casimir_omega(13) = 0.1839651500273005_DP ;  casimir_omega_weight(13) = 0.0373584845558574_DP
         casimir_omega(14) = 0.2242528462829711_DP ;  casimir_omega_weight(14) = 0.0433670043392562_DP
         casimir_omega(15) = 0.2710342689987175_DP ;  casimir_omega_weight(15) = 0.0503826761661462_DP
         casimir_omega(16) = 0.3254367902064166_DP ;  casimir_omega_weight(16) = 0.0586565901262108_DP
         casimir_omega(17) = 0.3888744273676722_DP ;  casimir_omega_weight(17) = 0.0685148244834558_DP
         casimir_omega(18) = 0.4631360119961658_DP ;  casimir_omega_weight(18) = 0.0803864914108685_DP
         casimir_omega(19) = 0.5505072507421337_DP ;  casimir_omega_weight(19) = 0.0948443490689901_DP
         casimir_omega(20) = 0.6539423404772365_DP ;  casimir_omega_weight(20) = 0.1126647024677741_DP
         casimir_omega(21) = 0.7773094526775429_DP ;  casimir_omega_weight(21) = 0.1349175577427744_DP
         casimir_omega(22) = 0.9257487113176190_DP ;  casimir_omega_weight(22) = 0.1631053780035364_DP
         casimir_omega(23) = 1.1062056006994803_DP ;  casimir_omega_weight(23) = 0.1993820319896580_DP
         casimir_omega(24) = 1.3282453223717738_DP ;  casimir_omega_weight(24) = 0.2469080909712281_DP
         casimir_omega(25) = 1.6053307949801348_DP ;  casimir_omega_weight(25) = 0.3104459484273310_DP
         casimir_omega(26) = 1.9568923785106906_DP ;  casimir_omega_weight(26) = 0.3973933850472089_DP
         casimir_omega(27) = 2.4118036964216798_DP ;  casimir_omega_weight(27) = 0.5196551551802508_DP
         casimir_omega(28) = 3.0144712939845415_DP ;  casimir_omega_weight(28) = 0.6972014197634652_DP
         casimir_omega(29) = 3.8360516252495285_DP ;  casimir_omega_weight(29) = 0.9652354402226746_DP
         casimir_omega(30) = 4.9963944703512997_DP ;  casimir_omega_weight(30) = 1.3896706650773698_DP
         casimir_omega(31) = 6.7102606760676284_DP ;  casimir_omega_weight(31) = 2.1035322208166551_DP
         casimir_omega(32) = 9.3941534943622091_DP ;  casimir_omega_weight(32) = 3.4023307361664163_DP
         casimir_omega(33) = 13.9450621841439126_DP ;  casimir_omega_weight(33) = 6.0318615926924881_DP
         casimir_omega(34) = 22.5987093411317304_DP ;  casimir_omega_weight(34) = 12.2429793833459897_DP
         casimir_omega(35) = 42.2995690749745421_DP ;  casimir_omega_weight(35) = 30.9512636618986008_DP
         casimir_omega(36) = 104.5384572489736712_DP ;  casimir_omega_weight(36) = 118.9784256515684007_DP
         casimir_omega(37) = 552.5132746946147790_DP ;  casimir_omega_weight(37) = 1418.9537663864907699_DP
return
endsubroutine gauss_legendre_grid36

subroutine gauss_legendre_grid37()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0006172596197449_DP ;  casimir_omega_weight(2) = 0.0015851737034176_DP
         casimir_omega(3) = 0.0032618479073338_DP ;  casimir_omega_weight(3) = 0.0037116635829829_DP
         casimir_omega(4) = 0.0080588994862846_DP ;  casimir_omega_weight(4) = 0.0058939076870840_DP
         casimir_omega(5) = 0.0150779144278168_DP ;  casimir_omega_weight(5) = 0.0081610227330707_DP
         casimir_omega(6) = 0.0244206425553320_DP ;  casimir_omega_weight(6) = 0.0105473586212181_DP
         casimir_omega(7) = 0.0362250773822672_DP ;  casimir_omega_weight(7) = 0.0130912017233207_DP
         casimir_omega(8) = 0.0506699505095668_DP ;  casimir_omega_weight(8) = 0.0158359824839579_DP
         casimir_omega(9) = 0.0679806586530855_DP ;  casimir_omega_weight(9) = 0.0188319040344515_DP
         casimir_omega(10) = 0.0884370210948336_DP ;  casimir_omega_weight(10) = 0.0221379910602339_DP
         casimir_omega(11) = 0.1123833500428328_DP ;  casimir_omega_weight(11) = 0.0258246948127490_DP
         casimir_omega(12) = 0.1402414851411728_DP ;  casimir_omega_weight(12) = 0.0299772559968938_DP
         casimir_omega(13) = 0.1725276837142622_DP ;  casimir_omega_weight(13) = 0.0347001090485880_DP
         casimir_omega(14) = 0.2098745969788165_DP ;  casimir_omega_weight(14) = 0.0401227290898848_DP
         casimir_omega(15) = 0.2530600451559576_DP ;  casimir_omega_weight(15) = 0.0464074972983160_DP
         casimir_omega(16) = 0.3030450023559451_DP ;  casimir_omega_weight(16) = 0.0537604229967392_DP
         casimir_omega(17) = 0.3610242267220091_DP ;  casimir_omega_weight(17) = 0.0624459622304919_DP
         casimir_omega(18) = 0.4284944999099455_DP ;  casimir_omega_weight(18) = 0.0728077969244789_DP
         casimir_omega(19) = 0.5073477600605824_DP ;  casimir_omega_weight(19) = 0.0852984277902057_DP
         casimir_omega(20) = 0.5999999999999996_DP ;  casimir_omega_weight(20) = 0.1005220331917680_DP
         casimir_omega(21) = 0.7095724635839766_DP ;  casimir_omega_weight(21) = 0.1192976894974544_DP
         casimir_omega(22) = 0.8401508072464391_DP ;  casimir_omega_weight(22) = 0.1427545263073184_DP
         casimir_omega(23) = 0.9971629972555893_DP ;  casimir_omega_weight(23) = 0.1724781836101386_DP
         casimir_omega(24) = 1.1879423755589882_DP ;  casimir_omega_weight(24) = 0.2107419165777545_DP
         casimir_omega(25) = 1.4225872748032440_DP ;  casimir_omega_weight(25) = 0.2608816222701748_DP
         casimir_omega(26) = 1.7153100240918413_DP ;  casimir_omega_weight(26) = 0.3279240098254806_DP
         casimir_omega(27) = 2.0866216496375527_DP ;  casimir_omega_weight(27) = 0.4196775684155405_DP
         casimir_omega(28) = 2.5670007675518365_DP ;  casimir_omega_weight(28) = 0.5487081021401182_DP
         casimir_omega(29) = 3.2033214872380351_DP ;  casimir_omega_weight(29) = 0.7360948019748002_DP
         casimir_omega(30) = 4.0706934216380111_DP ;  casimir_omega_weight(30) = 1.0189960433034313_DP
         casimir_omega(31) = 5.2956238896879286_DP ;  casimir_omega_weight(31) = 1.4669860938251769_DP
         casimir_omega(32) = 7.1048026765297685_DP ;  casimir_omega_weight(32) = 2.2204784020118398_DP
         casimir_omega(33) = 9.9378669699192130_DP ;  casimir_omega_weight(33) = 3.5913966402298376_DP
         casimir_omega(34) = 14.7416268504940664_DP ;  casimir_omega_weight(34) = 6.3669588013516893_DP
         casimir_omega(35) = 23.8759811062361607_DP ;  casimir_omega_weight(35) = 12.9230355773131844_DP
         casimir_omega(36) = 44.6711118078446887_DP ;  casimir_omega_weight(36) = 32.6703924925399676_DP
         casimir_omega(37) = 110.3668871839750238_DP ;  casimir_omega_weight(37) = 125.5867126750134020_DP
         casimir_omega(38) = 583.2229883249824525_DP ;  casimir_omega_weight(38) = 1497.7648217187911541_DP
return
endsubroutine gauss_legendre_grid37

subroutine gauss_legendre_grid38()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0005855918116490_DP ;  casimir_omega_weight(2) = 0.0015037953005407_DP
         casimir_omega(3) = 0.0030940369740213_DP ;  casimir_omega_weight(3) = 0.0035200584891716_DP
         casimir_omega(4) = 0.0076422163220581_DP ;  casimir_omega_weight(4) = 0.0055866134800316_DP
         casimir_omega(5) = 0.0142926453339688_DP ;  casimir_omega_weight(5) = 0.0077294157090508_DP
         casimir_omega(6) = 0.0231366487437085_DP ;  casimir_omega_weight(6) = 0.0099791299583351_DP
         casimir_omega(7) = 0.0342978353261543_DP ;  casimir_omega_weight(7) = 0.0123697783128797_DP
         casimir_omega(8) = 0.0479359114645430_DP ;  casimir_omega_weight(8) = 0.0149397336446002_DP
         casimir_omega(9) = 0.0642516809266127_DP ;  casimir_omega_weight(9) = 0.0177330683264206_DP
         casimir_omega(10) = 0.0834935551024162_DP ;  casimir_omega_weight(10) = 0.0208012385169405_DP
         casimir_omega(11) = 0.1059659557638371_DP ;  casimir_omega_weight(11) = 0.0242052062759160_DP
         casimir_omega(12) = 0.1320401204425854_DP ;  casimir_omega_weight(12) = 0.0280181526391137_DP
         casimir_omega(13) = 0.1621680021277410_DP ;  casimir_omega_weight(13) = 0.0323289948422376_DP
         casimir_omega(14) = 0.1969002080462272_DP ;  casimir_omega_weight(14) = 0.0372470058344817_DP
         casimir_omega(15) = 0.2369092784658836_DP ;  casimir_omega_weight(15) = 0.0429079582683401_DP
         casimir_omega(16) = 0.2830201145040508_DP ;  casimir_omega_weight(16) = 0.0494823991299848_DP
         casimir_omega(17) = 0.3362500989569211_DP ;  casimir_omega_weight(17) = 0.0571869380537173_DP
         casimir_omega(18) = 0.3978625336388423_DP ;  casimir_omega_weight(18) = 0.0662998556256738_DP
         casimir_omega(19) = 0.4694386275026698_DP ;  casimir_omega_weight(19) = 0.0771829961340843_DP
         casimir_omega(20) = 0.5529757149359503_DP ;  casimir_omega_weight(20) = 0.0903129518219681_DP
         casimir_omega(21) = 0.6510231648087795_DP ;  casimir_omega_weight(21) = 0.1063262312074074_DP
         casimir_omega(22) = 0.7668734077447699_DP ;  casimir_omega_weight(22) = 0.1260858902476246_DP
         casimir_omega(23) = 0.9048351366675509_DP ;  casimir_omega_weight(23) = 0.1507818250123313_DP
         casimir_omega(24) = 1.0706316551779556_DP ;  casimir_omega_weight(24) = 0.1820851394034960_DP
         casimir_omega(25) = 1.2719943973976695_DP ;  casimir_omega_weight(25) = 0.2223917355606698_DP
         casimir_omega(26) = 1.5195690195470404_DP ;  casimir_omega_weight(26) = 0.2752176043876479_DP
         casimir_omega(27) = 1.8283373266699687_DP ;  casimir_omega_weight(27) = 0.3458609401666540_DP
         casimir_omega(28) = 2.2199200537503376_DP ;  casimir_omega_weight(28) = 0.4425520634541867_DP
         casimir_omega(29) = 2.7264440443807221_DP ;  casimir_omega_weight(29) = 0.5785357143072216_DP
         casimir_omega(30) = 3.3973175384962357_DP ;  casimir_omega_weight(30) = 0.7760301052477355_DP
         casimir_omega(31) = 4.3117100422710628_DP ;  casimir_omega_weight(31) = 1.0742015823275386_DP
         casimir_omega(32) = 5.6029662540842571_DP ;  casimir_omega_weight(32) = 1.5463841876415865_DP
         casimir_omega(33) = 7.5100272217893256_DP ;  casimir_omega_weight(33) = 2.3405793887995285_DP
         casimir_omega(34) = 10.4962892432303860_DP ;  casimir_omega_weight(34) = 3.7855675092012513_DP
         casimir_omega(35) = 15.5597296733778130_DP ;  casimir_omega_weight(35) = 6.7111086937092805_DP
         casimir_omega(36) = 25.1877795599112595_DP ;  casimir_omega_weight(36) = 13.6214685565436273_DP
         casimir_omega(37) = 47.1067534376011778_DP ;  casimir_omega_weight(37) = 34.4359820063537541_DP
         casimir_omega(38) = 116.3528435576925091_DP ;  casimir_omega_weight(38) = 132.3736006206111426_DP
         casimir_omega(39) = 614.7626944889938159_DP ;  casimir_omega_weight(39) = 1578.7059049152367152_DP
return
endsubroutine gauss_legendre_grid38

subroutine gauss_legendre_grid39()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0005563003858277_DP ;  casimir_omega_weight(2) = 0.0014285286982019_DP
         casimir_omega(3) = 0.0029388636455631_DP ;  casimir_omega_weight(3) = 0.0033429458439523_DP
         casimir_omega(4) = 0.0072571144002997_DP ;  casimir_omega_weight(4) = 0.0053028564997243_DP
         casimir_omega(5) = 0.0135674449895851_DP ;  casimir_omega_weight(5) = 0.0073314630700401_DP
         casimir_omega(6) = 0.0219520628588168_DP ;  casimir_omega_weight(6) = 0.0094562298527354_DP
         casimir_omega(7) = 0.0325220312943215_DP ;  casimir_omega_weight(7) = 0.0117074987374872_DP
         casimir_omega(8) = 0.0454205210458894_DP ;  casimir_omega_weight(8) = 0.0141193094276226_DP
         casimir_omega(9) = 0.0608270491149411_DP ;  casimir_omega_weight(9) = 0.0167305228138269_DP
         casimir_omega(10) = 0.0789629716509175_DP ;  casimir_omega_weight(10) = 0.0195862145021256_DP
         casimir_omega(11) = 0.1000985358782356_DP ;  casimir_omega_weight(11) = 0.0227394157233369_DP
         casimir_omega(12) = 0.1245618938089790_DP ;  casimir_omega_weight(12) = 0.0262533188635115_DP
         casimir_omega(13) = 0.1527506189999498_DP ;  casimir_omega_weight(13) = 0.0302041095016507_DP
         casimir_omega(14) = 0.1851464586696928_DP ;  casimir_omega_weight(14) = 0.0346846488350133_DP
         casimir_omega(15) = 0.2223343192512750_DP ;  casimir_omega_weight(15) = 0.0398093196757439_DP
         casimir_omega(16) = 0.2650268577608741_DP ;  casimir_omega_weight(16) = 0.0457204798481909_DP
         casimir_omega(17) = 0.3140965856140503_DP ;  casimir_omega_weight(17) = 0.0525971605707930_DP
         casimir_omega(18) = 0.3706181647990887_DP ;  casimir_omega_weight(18) = 0.0606669389417012_DP
         casimir_omega(19) = 0.4359247122397455_DP ;  casimir_omega_weight(19) = 0.0702223592530072_DP
         casimir_omega(20) = 0.5116836234332185_DP ;  casimir_omega_weight(20) = 0.0816439707010447_DP
         casimir_omega(21) = 0.6000000000000000_DP ;  casimir_omega_weight(21) = 0.0954331465673316_DP
         casimir_omega(22) = 0.7035597457360978_DP ;  casimir_omega_weight(22) = 0.1122596241832028_DP
         casimir_omega(23) = 0.8258306764724339_DP ;  casimir_omega_weight(23) = 0.1330316378427925_DP
         casimir_omega(24) = 0.9713501230981358_DP ;  casimir_omega_weight(24) = 0.1590014851024749_DP
         casimir_omega(25) = 1.1461442641797892_DP ;  casimir_omega_weight(25) = 0.1919280140613576_DP
         casimir_omega(26) = 1.3583528969158920_DP ;  casimir_omega_weight(26) = 0.2343330286404758_DP
         casimir_omega(27) = 1.6191832246696019_DP ;  casimir_omega_weight(27) = 0.2899173767754035_DP
         casimir_omega(28) = 1.9444066205028081_DP ;  casimir_omega_weight(28) = 0.3642579032253173_DP
         casimir_omega(29) = 2.3567825934644340_DP ;  casimir_omega_weight(29) = 0.4660178792768647_DP
         casimir_omega(30) = 2.8901294689054402_DP ;  casimir_omega_weight(30) = 0.6091388640924448_DP
         casimir_omega(31) = 3.5964562002976752_DP ;  casimir_omega_weight(31) = 0.8170080806059300_DP
         casimir_omega(32) = 4.5590989355302032_DP ;  casimir_omega_weight(32) = 1.1308526999524500_DP
         casimir_omega(33) = 5.9184196050630344_DP ;  casimir_omega_weight(33) = 1.6278654918340612_DP
         casimir_omega(34) = 7.9259328539248877_DP ;  casimir_omega_weight(34) = 2.4638356383904800_DP
         casimir_omega(35) = 11.0694192728009959_DP ;  casimir_omega_weight(35) = 3.9848437198835214_DP
         casimir_omega(36) = 16.3993699505744139_DP ;  casimir_omega_weight(36) = 7.0643115724493191_DP
         casimir_omega(37) = 26.5341042677050361_DP ;  casimir_omega_weight(37) = 14.3382785546284595_DP
         casimir_omega(38) = 49.6064937305019740_DP ;  casimir_omega_weight(38) = 36.2480323717203277_DP
         casimir_omega(39) = 122.4963262734222127_DP ;  casimir_omega_weight(39) = 139.3390895945037471_DP
         casimir_omega(40) = 647.1323931665131113_DP ;  casimir_omega_weight(40) = 1661.7770160253642189_DP
return
endsubroutine gauss_legendre_grid39

subroutine gauss_legendre_grid40()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0005291533477354_DP ;  casimir_omega_weight(2) = 0.0013587766334496_DP
         casimir_omega(3) = 0.0027950888837236_DP ;  casimir_omega_weight(3) = 0.0031788973461349_DP
         casimir_omega(4) = 0.0069004735618335_DP ;  casimir_omega_weight(4) = 0.0050402821887169_DP
         casimir_omega(5) = 0.0128963131752778_DP ;  casimir_omega_weight(5) = 0.0069637273763413_DP
         casimir_omega(6) = 0.0208568128585072_DP ;  casimir_omega_weight(6) = 0.0089739063400921_DP
         casimir_omega(7) = 0.0308820528163218_DP ;  casimir_omega_weight(7) = 0.0110979706259992_DP
         casimir_omega(8) = 0.0431007735881692_DP ;  casimir_omega_weight(8) = 0.0133662291077987_DP
         casimir_omega(9) = 0.0576739866815167_DP ;  casimir_omega_weight(9) = 0.0158130887200253_DP
         casimir_omega(10) = 0.0747996324689095_DP ;  casimir_omega_weight(10) = 0.0184782143995122_DP
         casimir_omega(11) = 0.0947185305281534_DP ;  casimir_omega_weight(11) = 0.0214079666724717_DP
         casimir_omega(12) = 0.1177219428246802_DP ;  casimir_omega_weight(12) = 0.0246572073623603_DP
         casimir_omega(13) = 0.1441611766829549_DP ;  casimir_omega_weight(13) = 0.0282915974739560_DP
         casimir_omega(14) = 0.1744598001240762_DP ;  casimir_omega_weight(14) = 0.0323905570476404_DP
         casimir_omega(15) = 0.2091292425078095_DP ;  casimir_omega_weight(15) = 0.0370511218139768_DP
         casimir_omega(16) = 0.2487888322939187_DP ;  casimir_omega_weight(16) = 0.0423930254103856_DP
         casimir_omega(17) = 0.2941917167918480_DP ;  casimir_omega_weight(17) = 0.0485654733423537_DP
         casimir_omega(18) = 0.3462586700452960_DP ;  casimir_omega_weight(18) = 0.0557562786314628_DP
         casimir_omega(19) = 0.4061226076222396_DP ;  casimir_omega_weight(19) = 0.0642043356351066_DP
         casimir_omega(20) = 0.4751878209980847_DP ;  casimir_omega_weight(20) = 0.0742168770243187_DP
         casimir_omega(21) = 0.5552097261891438_DP ;  casimir_omega_weight(21) = 0.0861936873113640_DP
         casimir_omega(22) = 0.6484036266276763_DP ;  casimir_omega_weight(22) = 0.1006616001284906_DP
         casimir_omega(23) = 0.7575951741436805_DP ;  casimir_omega_weight(23) = 0.1183244716910904_DP
         casimir_omega(24) = 0.8864318145392656_DP ;  casimir_omega_weight(24) = 0.1401369061218422_DP
         casimir_omega(25) = 1.0396851577836490_DP ;  casimir_omega_weight(25) = 0.1674152313321087_DP
         casimir_omega(26) = 1.2236918290079322_DP ;  casimir_omega_weight(26) = 0.2020083146766133_DP
         casimir_omega(27) = 1.4470102885273268_DP ;  casimir_omega_weight(27) = 0.2465671122173110_DP
         casimir_omega(28) = 1.7214235354319549_DP ;  casimir_omega_weight(28) = 0.3049820883005119_DP
         casimir_omega(29) = 2.0635126243637050_DP ;  casimir_omega_weight(29) = 0.3831159002271197_DP
         casimir_omega(30) = 2.4972049221804467_DP ;  casimir_omega_weight(30) = 0.4900758865452866_DP
         casimir_omega(31) = 3.0580535061007059_DP ;  casimir_omega_weight(31) = 0.6405183062380625_DP
         casimir_omega(32) = 3.8007346397017452_DP ;  casimir_omega_weight(32) = 0.8590293794038457_DP
         casimir_omega(33) = 4.8128578726591202_DP ;  casimir_omega_weight(33) = 1.1889499548322040_DP
         casimir_omega(34) = 6.2419822300124110_DP ;  casimir_omega_weight(34) = 1.7114304814243180_DP
         casimir_omega(35) = 8.3525182967671157_DP ;  casimir_omega_weight(35) = 2.5902475498099737_DP
         casimir_omega(36) = 11.6572561461889936_DP ;  casimir_omega_weight(36) = 4.1892256016665801_DP
         casimir_omega(37) = 17.2605470664306537_DP ;  casimir_omega_weight(37) = 7.4265677025394767_DP
         casimir_omega(38) = 27.9149548485002441_DP ;  casimir_omega_weight(38) = 15.0734657762871951_DP
         casimir_omega(39) = 52.1703324814048486_DP ;  casimir_omega_weight(39) = 38.1065437363390913_DP
         casimir_omega(40) = 128.7973352462429375_DP ;  casimir_omega_weight(40) = 146.4831796898290577_DP
         casimir_omega(41) = 680.3320843394149051_DP ;  casimir_omega_weight(41) = 1746.9781550900265756_DP
return
endsubroutine gauss_legendre_grid40

subroutine gauss_legendre_grid41()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0005039463422178_DP ;  casimir_omega_weight(2) = 0.0012940130824617_DP
         casimir_omega(3) = 0.0026616219978600_DP ;  casimir_omega_weight(3) = 0.0030266567240921_DP
         casimir_omega(4) = 0.0065695505208698_DP ;  casimir_omega_weight(4) = 0.0047968244187812_DP
         casimir_omega(5) = 0.0122739835567158_DP ;  casimir_omega_weight(5) = 0.0066232021425162_DP
         casimir_omega(6) = 0.0198420785729151_DP ;  casimir_omega_weight(6) = 0.0085280203410540_DP
         casimir_omega(7) = 0.0293642659395816_DP ;  casimir_omega_weight(7) = 0.0105356532111916_DP
         casimir_omega(8) = 0.0409566451611388_DP ;  casimir_omega_weight(8) = 0.0126731825402479_DP
         casimir_omega(9) = 0.0547640741740768_DP ;  casimir_omega_weight(9) = 0.0149711894692575_DP
         casimir_omega(10) = 0.0709641379023702_DP ;  casimir_omega_weight(10) = 0.0174647254815844_DP
         casimir_omega(11) = 0.0897721946174145_DP ;  casimir_omega_weight(11) = 0.0201945068417632_DP
         casimir_omega(12) = 0.1114477567917254_DP ;  casimir_omega_weight(12) = 0.0232084028151600_DP
         casimir_omega(13) = 0.1363025457606166_DP ;  casimir_omega_weight(13) = 0.0265633135462962_DP
         casimir_omega(14) = 0.1647106715059835_DP ;  casimir_omega_weight(14) = 0.0303275675689447_DP
         casimir_omega(15) = 0.1971215413980226_DP ;  casimir_omega_weight(15) = 0.0345840167963453_DP
         casimir_omega(16) = 0.2340763116846360_DP ;  casimir_omega_weight(16) = 0.0394340751588515_DP
         casimir_omega(17) = 0.2762289879640847_DP ;  casimir_omega_weight(17) = 0.0450030457447398_DP
         casimir_omega(18) = 0.3243736932684883_DP ;  casimir_omega_weight(18) = 0.0514472256488111_DP
         casimir_omega(19) = 0.3794802114650230_DP ;  casimir_omega_weight(19) = 0.0589634917292714_DP
         casimir_omega(20) = 0.4427407667178324_DP ;  casimir_omega_weight(20) = 0.0678023923811701_DP
         casimir_omega(21) = 0.5156322531810280_DP ;  casimir_omega_weight(21) = 0.0782862623902921_DP
         casimir_omega(22) = 0.6000000000000008_DP ;  casimir_omega_weight(22) = 0.0908346427767593_DP
         casimir_omega(23) = 0.6981719971531953_DP ;  casimir_omega_weight(23) = 0.1060004990484964_DP
         casimir_omega(24) = 0.8131169005935145_DP ;  casimir_omega_weight(24) = 0.1245226897773763_DP
         casimir_omega(25) = 0.9486660677514176_DP ;  casimir_omega_weight(25) = 0.1474033747998376_DP
         casimir_omega(26) = 1.1098310605047226_DP ;  casimir_omega_weight(26) = 0.1760245364736965_DP
         casimir_omega(27) = 1.3032665494426920_DP ;  casimir_omega_weight(27) = 0.2123273323862159_DP
         casimir_omega(28) = 1.5379599815508787_DP ;  casimir_omega_weight(28) = 0.2590951176020441_DP
         casimir_omega(29) = 1.8262844204992146_DP ;  casimir_omega_weight(29) = 0.3204127292507146_DP
         casimir_omega(30) = 2.1856507335465629_DP ;  casimir_omega_weight(30) = 0.4024357966468406_DP
         casimir_omega(31) = 2.6411832441651888_DP ;  casimir_omega_weight(31) = 0.5147268398875027_DP
         casimir_omega(32) = 3.2302130645192917_DP ;  casimir_omega_weight(32) = 0.6726746965419456_DP
         casimir_omega(33) = 4.0101503760070178_DP ;  casimir_omega_weight(33) = 0.9020945689241842_DP
         casimir_omega(34) = 5.0729848997147533_DP ;  casimir_omega_weight(34) = 1.2484938345566927_DP
         casimir_omega(35) = 6.5736526259109151_DP ;  casimir_omega_weight(35) = 1.7970795718149597_DP
         casimir_omega(36) = 8.7897824292889446_DP ;  casimir_omega_weight(36) = 2.7198154726094366_DP
         casimir_omega(37) = 12.2597990612371248_DP ;  casimir_omega_weight(37) = 4.3987134435388020_DP
         casimir_omega(38) = 18.1432604793438053_DP ;  casimir_omega_weight(38) = 7.7978773167490942_DP
         casimir_omega(39) = 29.3303309668379377_DP ;  casimir_omega_weight(39) = 15.8270304015510721_DP
         casimir_omega(40) = 54.7982695096686143_DP ;  casimir_omega_weight(40) = 40.0115162302073486_DP
         casimir_omega(41) = 135.2558704013768249_DP ;  casimir_omega_weight(41) = 153.8058709885951885_DP
         casimir_omega(42) = 714.3617679924358299_DP ;  casimir_omega_weight(42) = 1834.3093221486356015_DP
return
endsubroutine gauss_legendre_grid41

subroutine gauss_legendre_grid42()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0004804987935535_DP ;  casimir_omega_weight(2) = 0.0012337732987616_DP
         casimir_omega(3) = 0.0025374998123691_DP ;  casimir_omega_weight(3) = 0.0028851154167672_DP
         casimir_omega(4) = 0.0062619254177873_DP ;  casimir_omega_weight(4) = 0.0045706639508027_DP
         casimir_omega(5) = 0.0116958181530699_DP ;  casimir_omega_weight(5) = 0.0063072481998333_DP
         casimir_omega(6) = 0.0189001084895776_DP ;  casimir_omega_weight(6) = 0.0081149524173616_DP
         casimir_omega(7) = 0.0279567195275157_DP ;  casimir_omega_weight(7) = 0.0100157233844229_DP
         casimir_omega(8) = 0.0389706371107618_DP ;  casimir_omega_weight(8) = 0.0120338393820173_DP
         casimir_omega(9) = 0.0520725642180128_DP ;  casimir_omega_weight(9) = 0.0141965797241183_DP
         casimir_omega(10) = 0.0674223172993084_DP ;  casimir_omega_weight(10) = 0.0165350419067731_DP
         casimir_omega(11) = 0.0852131264368989_DP ;  casimir_omega_weight(11) = 0.0190851397730751_DP
         casimir_omega(12) = 0.1056770463512751_DP ;  casimir_omega_weight(12) = 0.0218888377309039_DP
         casimir_omega(13) = 0.1290917498258622_DP ;  casimir_omega_weight(13) = 0.0249956957156133_DP
         casimir_omega(14) = 0.1557890619818371_DP ;  casimir_omega_weight(14) = 0.0284648252550933_DP
         casimir_omega(15) = 0.1861657109921162_DP ;  casimir_omega_weight(15) = 0.0323673925939264_DP
         casimir_omega(16) = 0.2206969304319428_DP ;  casimir_omega_weight(16) = 0.0367898550611298_DP
         casimir_omega(17) = 0.2599537683429559_DP ;  casimir_omega_weight(17) = 0.0418381885653249_DP
         casimir_omega(18) = 0.3046252645316335_DP ;  casimir_omega_weight(18) = 0.0476434676510655_DP
         casimir_omega(19) = 0.3555470899083413_DP ;  casimir_omega_weight(19) = 0.0543693109808279_DP
         casimir_omega(20) = 0.4137388593012837_DP ;  casimir_omega_weight(20) = 0.0622219295745220_DP
         casimir_omega(21) = 0.4804532236695494_DP ;  casimir_omega_weight(21) = 0.0714638527745239_DP
         casimir_omega(22) = 0.5572411620739249_DP ;  casimir_omega_weight(22) = 0.0824329228762321_DP
         casimir_omega(23) = 0.6460398558142436_DP ;  casimir_omega_weight(23) = 0.0955689515309775_DP
         casimir_omega(24) = 0.7492925060434267_DP ;  casimir_omega_weight(24) = 0.1114517016411391_DP
         casimir_omega(25) = 0.8701140632716076_DP ;  casimir_omega_weight(25) = 0.1308559124905941_DP
         casimir_omega(26) = 1.0125241078271978_DP ;  casimir_omega_weight(26) = 0.1548324811440139_DP
         casimir_omega(27) = 1.1817798518895211_DP ;  casimir_omega_weight(27) = 0.1848306647537865_DP
         casimir_omega(28) = 1.3848616324924874_DP ;  casimir_omega_weight(28) = 0.2228861789018722_DP
         casimir_omega(29) = 1.6311962259530117_DP ;  casimir_omega_weight(29) = 0.2719180217487407_DP
         casimir_omega(30) = 1.9337610459062766_DP ;  casimir_omega_weight(30) = 0.3362101571880811_DP
         casimir_omega(31) = 2.3108169175700635_DP ;  casimir_omega_weight(31) = 0.4222183439477494_DP
         casimir_omega(32) = 2.7887142322078673_DP ;  casimir_omega_weight(32) = 0.5399713961589183_DP
         casimir_omega(33) = 3.4066054306944196_DP ;  casimir_omega_weight(33) = 0.7056086071692779_DP
         casimir_omega(34) = 4.2247012291772190_DP ;  casimir_omega_weight(34) = 0.9462041451798708_DP
         casimir_omega(35) = 5.3394782976955444_DP ;  casimir_omega_weight(35) = 1.3094847663090841_DP
         casimir_omega(36) = 6.9134294691689018_DP ;  casimir_omega_weight(36) = 1.8848131276042064_DP
         casimir_omega(37) = 9.2377242634451289_DP ;  casimir_omega_weight(37) = 2.8525397140854802_DP
         casimir_omega(38) = 12.8770473104213661_DP ;  casimir_omega_weight(38) = 4.6133074999150132_DP
         casimir_omega(39) = 19.0475097113078249_DP ;  casimir_omega_weight(39) = 8.1782406202451199_DP
         casimir_omega(40) = 30.7802323265013555_DP ;  casimir_omega_weight(40) = 16.5989725892597306_DP
         casimir_omega(41) = 57.4903046557236195_DP ;  casimir_omega_weight(41) = 41.9629499681018103_DP
         casimir_omega(42) = 141.8719316727444948_DP ;  casimir_omega_weight(42) = 161.3071635632805396_DP
         casimir_omega(43) = 749.2214441114584815_DP ;  casimir_omega_weight(43) = 1923.7705172330970527_DP
return
endsubroutine gauss_legendre_grid42

subroutine gauss_legendre_grid43()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0004586506601604_DP ;  casimir_omega_weight(2) = 0.0011776454400337_DP
         casimir_omega(3) = 0.0024218691718231_DP ;  casimir_omega_weight(3) = 0.0027532921792167_DP
         casimir_omega(4) = 0.0059754570274421_DP ;  casimir_omega_weight(4) = 0.0043601937340014_DP
         casimir_omega(5) = 0.0111577191544596_DP ;  casimir_omega_weight(5) = 0.0060135408156727_DP
         casimir_omega(6) = 0.0180240672313338_DP ;  casimir_omega_weight(6) = 0.0077315257803709_DP
         casimir_omega(7) = 0.0266489002046902_DP ;  casimir_omega_weight(7) = 0.0095339658899985_DP
         casimir_omega(8) = 0.0371273997589811_DP ;  casimir_omega_weight(8) = 0.0114426939191506_DP
         casimir_omega(9) = 0.0495778201384977_DP ;  casimir_omega_weight(9) = 0.0134821268355892_DP
         casimir_omega(10) = 0.0641444070287795_DP ;  casimir_omega_weight(10) = 0.0156799569330503_DP
         casimir_omega(11) = 0.0810010783426413_DP ;  casimir_omega_weight(11) = 0.0180679905190586_DP
         casimir_omega(12) = 0.1003560349267405_DP ;  casimir_omega_weight(12) = 0.0206831774894625_DP
         casimir_omega(13) = 0.1224575199882184_DP ;  casimir_omega_weight(13) = 0.0235688904120879_DP
         casimir_omega(14) = 0.1476010139297846_DP ;  casimir_omega_weight(14) = 0.0267765312416472_DP
         casimir_omega(15) = 0.1761382420387593_DP ;  casimir_omega_weight(15) = 0.0303675705106352_DP
         casimir_omega(16) = 0.2084884949181558_DP ;  casimir_omega_weight(16) = 0.0344161611416075_DP
         casimir_omega(17) = 0.2451529284893114_DP ;  casimir_omega_weight(17) = 0.0390125216834883_DP
         casimir_omega(18) = 0.2867327405176131_DP ;  casimir_omega_weight(18) = 0.0442673589273721_DP
         casimir_omega(19) = 0.3339524414584327_DP ;  casimir_omega_weight(19) = 0.0503177083722313_DP
         casimir_omega(20) = 0.3876898901307878_DP ;  casimir_omega_weight(20) = 0.0573347296806606_DP
         casimir_omega(21) = 0.4490154116302550_DP ;  casimir_omega_weight(21) = 0.0655342294433268_DP
         casimir_omega(22) = 0.5192432518619617_DP ;  casimir_omega_weight(22) = 0.0751910373017078_DP
         casimir_omega(23) = 0.6000000000000004_DP ;  casimir_omega_weight(23) = 0.0866589020325587_DP
         casimir_omega(24) = 0.6933166655687315_DP ;  casimir_omega_weight(24) = 0.1003984145691559_DP
         casimir_omega(25) = 0.8017542175065576_DP ;  casimir_omega_weight(25) = 0.1170167960526409_DP
         casimir_omega(26) = 0.9285772189688852_DP ;  casimir_omega_weight(26) = 0.1373255408317209_DP
         casimir_omega(27) = 1.0779978083939512_DP ;  casimir_omega_weight(27) = 0.1624254612776128_DP
         casimir_omega(28) = 1.2555245674076978_DP ;  casimir_omega_weight(28) = 0.1938347067280835_DP
         casimir_omega(29) = 1.4684711384783495_DP ;  casimir_omega_weight(29) = 0.2336858159700130_DP
         casimir_omega(30) = 1.7267139855430484_DP ;  casimir_omega_weight(30) = 0.2850366721446430_DP
         casimir_omega(31) = 2.0438491711571714_DP ;  casimir_omega_weight(31) = 0.3523751179744405_DP
         casimir_omega(32) = 2.4390076356200110_DP ;  casimir_omega_weight(32) = 0.4424641973317461_DP
         casimir_omega(33) = 2.9397949593837578_DP ;  casimir_omega_weight(33) = 0.5658101294097092_DP
         casimir_omega(34) = 3.5872282146539471_DP ;  casimir_omega_weight(34) = 0.7393205392486525_DP
         casimir_omega(35) = 4.4443852769116221_DP ;  casimir_omega_weight(35) = 0.9913585434826135_DP
         casimir_omega(36) = 5.6123365492876180_DP ;  casimir_omega_weight(36) = 1.3719231256924638_DP
         casimir_omega(37) = 7.2613115904314585_DP ;  casimir_omega_weight(37) = 1.9746314699082623_DP
         casimir_omega(38) = 9.6963429256290983_DP ;  casimir_omega_weight(38) = 2.9884205452944621_DP
         casimir_omega(39) = 13.5090002677348178_DP ;  casimir_omega_weight(39) = 4.8330079955004628_DP
         casimir_omega(40) = 19.9732943391466691_DP ;  casimir_omega_weight(40) = 8.5676577944400893_DP
         casimir_omega(41) = 32.2646586651283442_DP ;  casimir_omega_weight(41) = 17.3892924799921502_DP
         casimir_omega(42) = 60.2464377781841307_DP ;  casimir_omega_weight(42) = 43.9608450516775378_DP
         casimir_omega(43) = 148.6455190017598511_DP ;  casimir_omega_weight(43) = 168.9870574780304082_DP
         casimir_omega(44) = 784.9111126845202762_DP ;  casimir_omega_weight(44) = 2015.3617403742064198_DP
return
endsubroutine gauss_legendre_grid43

subroutine gauss_legendre_grid44()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0004382596945793_DP ;  casimir_omega_weight(2) = 0.0011252635003871_DP
         casimir_omega(3) = 0.0023139721858847_DP ;  casimir_omega_weight(3) = 0.0026303159041813_DP
         casimir_omega(4) = 0.0057082450544349_DP ;  casimir_omega_weight(4) = 0.0041639897832930_DP
         casimir_omega(5) = 0.0106560548937059_DP ;  casimir_omega_weight(5) = 0.0057400255315396_DP
         casimir_omega(6) = 0.0172079079620483_DP ;  casimir_omega_weight(6) = 0.0073749423753385_DP
         casimir_omega(7) = 0.0254315282148877_DP ;  casimir_omega_weight(7) = 0.0090866827824966_DP
         casimir_omega(8) = 0.0354134204572702_DP ;  casimir_omega_weight(8) = 0.0108949380575062_DP
         casimir_omega(9) = 0.0472608531530478_DP ;  casimir_omega_weight(9) = 0.0128216333683443_DP
         casimir_omega(10) = 0.0611043769740797_DP ;  casimir_omega_weight(10) = 0.0148915150596542_DP
         casimir_omega(11) = 0.0771009887935204_DP ;  casimir_omega_weight(11) = 0.0171328589595107_DP
         casimir_omega(12) = 0.0954380782972279_DP ;  casimir_omega_weight(12) = 0.0195783339846621_DP
         casimir_omega(13) = 0.1163383345963656_DP ;  casimir_omega_weight(13) = 0.0222660673758898_DP
         casimir_omega(14) = 0.1400658436774485_DP ;  casimir_omega_weight(14) = 0.0252409728424756_DP
         casimir_omega(15) = 0.1669336783918929_DP ;  casimir_omega_weight(15) = 0.0285564231335158_DP
         casimir_omega(16) = 0.1973133774209765_DP ;  casimir_omega_weight(16) = 0.0322763765168344_DP
         casimir_omega(17) = 0.2316468375630617_DP ;  casimir_omega_weight(17) = 0.0364781057112498_DP
         casimir_omega(18) = 0.2704613182014333_DP ;  casimir_omega_weight(18) = 0.0412557329568376_DP
         casimir_omega(19) = 0.3143884974726432_DP ;  casimir_omega_weight(19) = 0.0467248535848666_DP
         casimir_omega(20) = 0.3641888552863846_DP ;  casimir_omega_weight(20) = 0.0530286440351278_DP
         casimir_omega(21) = 0.4207831319986607_DP ;  casimir_omega_weight(21) = 0.0603460163374061_DP
         casimir_omega(22) = 0.4852932884337938_DP ;  casimir_omega_weight(22) = 0.0689026272059224_DP
         casimir_omega(23) = 0.5590963734143120_DP ;  casimir_omega_weight(23) = 0.0789859200901815_DP
         casimir_omega(24) = 0.6438961458496638_DP ;  casimir_omega_weight(24) = 0.0909659442286693_DP
         casimir_omega(25) = 0.7418194493516316_DP ;  casimir_omega_weight(25) = 0.1053245742131283_DP
         casimir_omega(26) = 0.8555476030848732_DP ;  casimir_omega_weight(26) = 0.1226971465989105_DP
         casimir_omega(27) = 0.9884981233621474_DP ;  casimir_omega_weight(27) = 0.1439327819956035_DP
         casimir_omega(28) = 1.1450800614336263_DP ;  casimir_omega_weight(28) = 0.1701833834365791_DP
         casimir_omega(29) = 1.3310591044737869_DP ;  casimir_omega_weight(29) = 0.2030376074815984_DP
         casimir_omega(30) = 1.5540898541383952_DP ;  casimir_omega_weight(30) = 0.2447270792920200_DP
         casimir_omega(31) = 1.8245088331336217_DP ;  casimir_omega_weight(31) = 0.2984518071010944_DP
         casimir_omega(32) = 2.1565450630930525_DP ;  casimir_omega_weight(32) = 0.3689082629785979_DP
         casimir_omega(33) = 2.5702197662766988_DP ;  casimir_omega_weight(33) = 0.4631739303208044_DP
         casimir_omega(34) = 3.0944228422128894_DP ;  casimir_omega_weight(34) = 0.5922435432246497_DP
         casimir_omega(35) = 3.7720793044347793_DP ;  casimir_omega_weight(35) = 0.7738109332928677_DP
         casimir_omega(36) = 4.6692008187352076_DP ;  casimir_omega_weight(36) = 1.0375581472146458_DP
         casimir_omega(37) = 5.8915583109980885_DP ;  casimir_omega_weight(37) = 1.4358092440787460_DP
         casimir_omega(38) = 7.6172979534285847_DP ;  casimir_omega_weight(38) = 2.0665348824749841_DP
         casimir_omega(39) = 10.1656376410851461_DP ;  casimir_omega_weight(39) = 3.1274582060864820_DP
         casimir_omega(40) = 14.1556573776504688_DP ;  casimir_omega_weight(40) = 5.0578151293761380_DP
         casimir_omega(41) = 20.9206139871260568_DP ;  casimir_omega_weight(41) = 8.9661290002272214_DP
         casimir_omega(42) = 33.7836097496681589_DP ;  casimir_omega_weight(42) = 18.1979901985305155_DP
         casimir_omega(43) = 63.0666687514245154_DP ;  casimir_omega_weight(43) = 46.0052015712308773_DP
         casimir_omega(44) = 155.5766323363848471_DP ;  casimir_omega_weight(44) = 176.8455527899603510_DP
         casimir_omega(45) = 821.4307737003614420_DP ;  casimir_omega_weight(45) = 2109.0829915967838133_DP
return
endsubroutine gauss_legendre_grid44

subroutine gauss_legendre_grid45()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0004191991207750_DP ;  casimir_omega_weight(2) = 0.0010763013204747_DP
         casimir_omega(3) = 0.0022131337340298_DP ;  casimir_omega_weight(3) = 0.0025154110920849_DP
         casimir_omega(4) = 0.0054585982600069_DP ;  casimir_omega_weight(4) = 0.0039807866300479_DP
         casimir_omega(5) = 0.0101875974264646_DP ;  casimir_omega_weight(5) = 0.0054848811074563_DP
         casimir_omega(6) = 0.0164462651502385_DP ;  casimir_omega_weight(6) = 0.0070427295498089_DP
         casimir_omega(7) = 0.0242963865237822_DP ;  casimir_omega_weight(7) = 0.0086706183564513_DP
         casimir_omega(8) = 0.0338167636309127_DP ;  casimir_omega_weight(8) = 0.0103863567483083_DP
         casimir_omega(9) = 0.0451049386856762_DP ;  casimir_omega_weight(9) = 0.0122096920605363_DP
         casimir_omega(10) = 0.0582793753630270_DP ;  casimir_omega_weight(10) = 0.0141628110577925_DP
         casimir_omega(11) = 0.0734821894649248_DP ;  casimir_omega_weight(11) = 0.0162709410250888_DP
         casimir_omega(12) = 0.0908825415499163_DP ;  casimir_omega_weight(12) = 0.0185630779280569_DP
         casimir_omega(13) = 0.1106808361808671_DP ;  casimir_omega_weight(13) = 0.0210728785002401_DP
         casimir_omega(14) = 0.1331139148236604_DP ;  casimir_omega_weight(14) = 0.0238397646706747_DP
         casimir_omega(15) = 0.1584614851779747_DP ;  casimir_omega_weight(15) = 0.0269103042025897_DP
         casimir_omega(16) = 0.1870541036047721_DP ;  casimir_omega_weight(16) = 0.0303399525634237_DP
         casimir_omega(17) = 0.2192831261626030_DP ;  casimir_omega_weight(17) = 0.0341952703031518_DP
         casimir_omega(18) = 0.2556131773232365_DP ;  casimir_omega_weight(18) = 0.0385567710854084_DP
         casimir_omega(19) = 0.2965978677400232_DP ;  casimir_omega_weight(19) = 0.0435226131853143_DP
         casimir_omega(20) = 0.3428997439245241_DP ;  casimir_omega_weight(20) = 0.0492134295416164_DP
         casimir_omega(21) = 0.3953158034739666_DP ;  casimir_omega_weight(21) = 0.0557787102164709_DP
         casimir_omega(22) = 0.4548104045873023_DP ;  casimir_omega_weight(22) = 0.0634053247470795_DP
         casimir_omega(23) = 0.5225581061958876_DP ;  casimir_omega_weight(23) = 0.0723290292034214_DP
         casimir_omega(24) = 0.6000000000000000_DP ;  casimir_omega_weight(24) = 0.0828501897950784_DP
         casimir_omega(25) = 0.6889186020301623_DP ;  casimir_omega_weight(25) = 0.0953555462908481_DP
         casimir_omega(26) = 0.7915386199809261_DP ;  casimir_omega_weight(26) = 0.1103487579517588_DP
         casimir_omega(27) = 0.9106643266886437_DP ;  casimir_omega_weight(27) = 0.1284939310203636_DP
         casimir_omega(28) = 1.0498695504399087_DP ;  casimir_omega_weight(28) = 0.1506786810544704_DP
         casimir_omega(29) = 1.2137646259667336_DP ;  casimir_omega_weight(29) = 0.1781071749317887_DP
         casimir_omega(30) = 1.4083780960351693_DP ;  casimir_omega_weight(30) = 0.2124401895832738_DP
         casimir_omega(31) = 1.6417131874207795_DP ;  casimir_omega_weight(31) = 0.2560106980711059_DP
         casimir_omega(32) = 1.9245768633905329_DP ;  casimir_omega_weight(32) = 0.3121640723975104_DP
         casimir_omega(33) = 2.2718454241146930_DP ;  casimir_omega_weight(33) = 0.3858101632426535_DP
         casimir_omega(34) = 2.7044505488167894_DP ;  casimir_omega_weight(34) = 0.4843480468039645_DP
         casimir_omega(35) = 3.2525955930773089_DP ;  casimir_omega_weight(35) = 0.6192720809529166_DP
         casimir_omega(36) = 3.9611568279290870_DP ;  casimir_omega_weight(36) = 0.8090801778657992_DP
         casimir_omega(37) = 4.8991463458208271_DP ;  casimir_omega_weight(37) = 1.0848032951465028_DP
         casimir_omega(38) = 6.1771423896967024_DP ;  casimir_omega_weight(38) = 1.5011434147569267_DP
         casimir_omega(39) = 7.9813876371441195_DP ;  casimir_omega_weight(39) = 2.1605236168128950_DP
         casimir_omega(40) = 10.6456077207493340_DP ;  casimir_omega_weight(40) = 3.2696529093391460_DP
         casimir_omega(41) = 14.8170181457896248_DP ;  casimir_omega_weight(41) = 5.2877290784373008_DP
         casimir_omega(42) = 21.8894683207013507_DP ;  casimir_omega_weight(42) = 9.3736543807074799_DP
         casimir_omega(43) = 35.3370853725354124_DP ;  casimir_omega_weight(43) = 19.0250658559496344_DP
         casimir_omega(44) = 65.9509974635029721_DP ;  casimir_omega_weight(44) = 48.0960196071840400_DP
         casimir_omega(45) = 162.6652716302360773_DP ;  casimir_omega_weight(45) = 184.8826495499031637_DP
         casimir_omega(46) = 858.7804271497717536_DP ;  casimir_omega_weight(46) = 2204.9342709264428777_DP
return
endsubroutine gauss_legendre_grid45

subroutine gauss_legendre_grid46()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0004013556576563_DP ;  casimir_omega_weight(2) = 0.0010304674915718_DP
         casimir_omega(3) = 0.0021187508427405_DP ;  casimir_omega_weight(3) = 0.0024078855122674_DP
         casimir_omega(4) = 0.0052250074103613_DP ;  casimir_omega_weight(4) = 0.0038094565411414_DP
         casimir_omega(5) = 0.0097494696790513_DP ;  casimir_omega_weight(5) = 0.0052464882896188_DP
         casimir_omega(6) = 0.0157343640457003_DP ;  casimir_omega_weight(6) = 0.0067326953385416_DP
         casimir_omega(7) = 0.0232361770857424_DP ;  casimir_omega_weight(7) = 0.0082828965806994_DP
         casimir_omega(8) = 0.0323268530756292_DP ;  casimir_omega_weight(8) = 0.0099132414037131_DP
         casimir_omega(9) = 0.0430952965945313_DP ;  casimir_omega_weight(9) = 0.0116415665776729_DP
         casimir_omega(10) = 0.0556492685563832_DP ;  casimir_omega_weight(10) = 0.0134878259685100_DP
         casimir_omega(11) = 0.0701177518402246_DP ;  casimir_omega_weight(11) = 0.0154746029757422_DP
         casimir_omega(12) = 0.0866538794913656_DP ;  casimir_omega_weight(12) = 0.0176277275065003_DP
         casimir_omega(13) = 0.1054385441161237_DP ;  casimir_omega_weight(13) = 0.0199770269866896_DP
         casimir_omega(14) = 0.1266848409091292_DP ;  casimir_omega_weight(14) = 0.0225572499172769_DP
         casimir_omega(15) = 0.1506435409169750_DP ;  casimir_omega_weight(15) = 0.0254092123839814_DP
         casimir_omega(16) = 0.1776098491988476_DP ;  casimir_omega_weight(16) = 0.0285812340585184_DP
         casimir_omega(17) = 0.2079317795482823_DP ;  casimir_omega_weight(17) = 0.0321309523361868_DP
         casimir_omega(18) = 0.2420205805378205_DP ;  casimir_omega_weight(18) = 0.0361276338337311_DP
         casimir_omega(19) = 0.2803637870211371_DP ;  casimir_omega_weight(19) = 0.0406551451757638_DP
         casimir_omega(20) = 0.3235416615298846_DP ;  casimir_omega_weight(20) = 0.0458158052479164_DP
         casimir_omega(21) = 0.3722480525749608_DP ;  casimir_omega_weight(21) = 0.0517354270396537_DP
         casimir_omega(22) = 0.4273170631657828_DP ;  casimir_omega_weight(22) = 0.0585699812554852_DP
         casimir_omega(23) = 0.4897574398973116_DP ;  casimir_omega_weight(23) = 0.0665144952302132_DP
         casimir_omega(24) = 0.5607973319357400_DP ;  casimir_omega_weight(24) = 0.0758150694575677_DP
         casimir_omega(25) = 0.6419431397031878_DP ;  casimir_omega_weight(25) = 0.0867852982759610_DP
         casimir_omega(26) = 0.7350577462906576_DP ;  casimir_omega_weight(26) = 0.0998289989628985_DP
         casimir_omega(27) = 0.8424657731496522_DP ;  casimir_omega_weight(27) = 0.1154721137887720_DP
         casimir_omega(28) = 0.9670970674252375_DP ;  casimir_omega_weight(28) = 0.1344081706430562_DP
         casimir_omega(29) = 1.1126851432292217_DP ;  casimir_omega_weight(29) = 0.1575641467110769_DP
         casimir_omega(30) = 1.2840460026061047_DP ;  casimir_omega_weight(30) = 0.1861976441500090_DP
         casimir_omega(31) = 1.4874768054849072_DP ;  casimir_omega_weight(31) = 0.2220431718877267_DP
         casimir_omega(32) = 1.7313370797964391_DP ;  casimir_omega_weight(32) = 0.2675373110818479_DP
         casimir_omega(33) = 2.0269146200161066_DP ;  casimir_omega_weight(33) = 0.3261740350134186_DP
         casimir_omega(34) = 2.3897473320705420_DP ;  casimir_omega_weight(34) = 0.4030813212104344_DP
         casimir_omega(35) = 2.8416975339474693_DP ;  casimir_omega_weight(35) = 0.5059869910445172_DP
         casimir_omega(36) = 3.4143111802029211_DP ;  casimir_omega_weight(36) = 0.6468961342329487_DP
         casimir_omega(37) = 4.1544591207352939_DP ;  casimir_omega_weight(37) = 0.8451286168268527_DP
         casimir_omega(38) = 5.1342205155168239_DP ;  casimir_omega_weight(38) = 1.1330942875717707_DP
         casimir_omega(39) = 6.4690877227838559_DP ;  casimir_omega_weight(39) = 1.5679258981010817_DP
         casimir_omega(40) = 8.3535798207194922_DP ;  casimir_omega_weight(40) = 2.2565978965127211_DP
         casimir_omega(41) = 11.1362525500942535_DP ;  casimir_omega_weight(41) = 3.4150048445335828_DP
         casimir_omega(42) = 15.4930821310056643_DP ;  casimir_omega_weight(42) = 5.5227500003062175_DP
         casimir_omega(43) = 22.8798570412112845_DP ;  casimir_omega_weight(43) = 9.7902340635084233_DP
         casimir_omega(44) = 36.9250853483375749_DP ;  casimir_omega_weight(44) = 19.8705195513848736_DP
         casimir_omega(45) = 68.8994238144248499_DP ;  casimir_omega_weight(45) = 50.2332992313544580_DP
         casimir_omega(46) = 169.9114368418757408_DP ;  casimir_omega_weight(46) = 193.0983478033149936_DP
         casimir_omega(47) = 896.9600730237266362_DP ;  casimir_omega_weight(47) = 2302.9155783824999162_DP
return
endsubroutine gauss_legendre_grid46

subroutine gauss_legendre_grid47()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0003846278310591_DP ;  casimir_omega_weight(2) = 0.0009875010042838_DP
         casimir_omega(3) = 0.0020302836209986_DP ;  casimir_omega_weight(3) = 0.0023071196852828_DP
         casimir_omega(4) = 0.0050061222292099_DP ;  casimir_omega_weight(4) = 0.0036489918576262_DP
         casimir_omega(5) = 0.0093391005192281_DP ;  casimir_omega_weight(5) = 0.0050234033740230_DP
         casimir_omega(6) = 0.0150679439450653_DP ;  casimir_omega_weight(6) = 0.0064428908013967_DP
         casimir_omega(7) = 0.0222443994234592_DP ;  casimir_omega_weight(7) = 0.0079209686980645_DP
         casimir_omega(8) = 0.0309342887849509_DP ;  casimir_omega_weight(8) = 0.0094723178287980_DP
         casimir_omega(9) = 0.0412188233371235_DP ;  casimir_omega_weight(9) = 0.0111130929175960_DP
         casimir_omega(10) = 0.0531962575215053_DP ;  casimir_omega_weight(10) = 0.0128612924548743_DP
         casimir_omega(11) = 0.0669839456798175_DP ;  casimir_omega_weight(11) = 0.0147371974482650_DP
         casimir_omega(12) = 0.0827208790707953_DP ;  casimir_omega_weight(12) = 0.0167638966156647_DP
         casimir_omega(13) = 0.1005708009268746_DP ;  casimir_omega_weight(13) = 0.0189679217688973_DP
         casimir_omega(14) = 0.1207260244892822_DP ;  casimir_omega_weight(14) = 0.0213800242130626_DP
         casimir_omega(15) = 0.1434121141806867_DP ;  casimir_omega_weight(15) = 0.0240361322178050_DP
         casimir_omega(16) = 0.1688936360191266_DP ;  casimir_omega_weight(16) = 0.0269785420155575_DP
         casimir_omega(17) = 0.1974812438185193_DP ;  casimir_omega_weight(17) = 0.0302574116258260_DP
         casimir_omega(18) = 0.2295404479494901_DP ;  casimir_omega_weight(18) = 0.0339326498905008_DP
         casimir_omega(19) = 0.2655025209093222_DP ;  casimir_omega_weight(19) = 0.0380763250312049_DP
         casimir_omega(20) = 0.3058781392933670_DP ;  casimir_omega_weight(20) = 0.0427757616155834_DP
         casimir_omega(21) = 0.3512745602616617_DP ;  casimir_omega_weight(21) = 0.0481375577013735_DP
         casimir_omega(22) = 0.4024174045191481_DP ;  casimir_omega_weight(22) = 0.0542928436187289_DP
         casimir_omega(23) = 0.4601785000031226_DP ;  casimir_omega_weight(23) = 0.0614042333088574_DP
         casimir_omega(24) = 0.5256117799397673_DP ;  casimir_omega_weight(24) = 0.0696751083853498_DP
         casimir_omega(25) = 0.6000000000000001_DP ;  casimir_omega_weight(25) = 0.0793621555483869_DP
         casimir_omega(26) = 0.6849161562574849_DP ;  casimir_omega_weight(26) = 0.0907924997944053_DP
         casimir_omega(27) = 0.7823051272442263_DP ;  casimir_omega_weight(27) = 0.1043874204285811_DP
         casimir_omega(28) = 0.8945935139911939_DP ;  casimir_omega_weight(28) = 0.1206956389361178_DP
         casimir_omega(29) = 1.0248393727454643_DP ;  casimir_omega_weight(29) = 0.1404407549565411_DP
         casimir_omega(30) = 1.1769392897173496_DP ;  casimir_omega_weight(30) = 0.1645899723637259_DP
         casimir_omega(31) = 1.3559193290031000_DP ;  casimir_omega_weight(31) = 0.1944554986159561_DP
         casimir_omega(32) = 1.5683510388514055_DP ;  casimir_omega_weight(32) = 0.2318471850262246_DP
         casimir_omega(33) = 1.8229579328092078_DP ;  casimir_omega_weight(33) = 0.2793074799562364_DP
         casimir_omega(34) = 2.1315190346143749_DP ;  casimir_omega_weight(34) = 0.3404821945202946_DP
         casimir_omega(35) = 2.5102481896782551_DP ;  casimir_omega_weight(35) = 0.4207221804888328_DP
         casimir_omega(36) = 2.9819585422690671_DP ;  casimir_omega_weight(36) = 0.5280911560350566_DP
         casimir_omega(37) = 3.5795677938545736_DP ;  casimir_omega_weight(37) = 0.6751160501313501_DP
         casimir_omega(38) = 4.3519846989524797_DP ;  casimir_omega_weight(38) = 0.8819565554151765_DP
         casimir_omega(39) = 5.3744221297562325_DP ;  casimir_omega_weight(39) = 1.1824313914730560_DP
         casimir_omega(40) = 6.7673933613556496_DP ;  casimir_omega_weight(40) = 1.6361569259338873_DP
         casimir_omega(41) = 8.7338737706218854_DP ;  casimir_omega_weight(41) = 2.3547579209049267_DP
         casimir_omega(42) = 11.6375715796360719_DP ;  casimir_omega_weight(42) = 3.5635141807858361_DP
         casimir_omega(43) = 16.1838489386385049_DP ;  casimir_omega_weight(43) = 5.7628780358065557_DP
         casimir_omega(44) = 23.8917798813485156_DP ;  casimir_omega_weight(44) = 10.2158681627527628_DP
         casimir_omega(45) = 38.5476095110883463_DP ;  casimir_omega_weight(45) = 20.7343513735430847_DP
         casimir_omega(46) = 71.9119477146307986_DP ;  casimir_omega_weight(46) = 52.4170405080467248_DP
         casimir_omega(47) = 177.3151279341557540_DP ;  casimir_omega_weight(47) = 201.4926475907967358_DP
         casimir_omega(48) = 935.9697113144808327_DP ;  casimir_omega_weight(48) = 2403.0269139838537740_DP
return
endsubroutine gauss_legendre_grid47

subroutine gauss_legendre_grid48()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0003689245270563_DP ;  casimir_omega_weight(2) = 0.0009471675200488_DP
         casimir_omega(3) = 0.0019472474980282_DP ;  casimir_omega_weight(3) = 0.0022125578850667_DP
         casimir_omega(4) = 0.0048007316901601_DP ;  casimir_omega_weight(4) = 0.0034984899276215_DP
         casimir_omega(5) = 0.0089541864174159_DP ;  casimir_omega_weight(5) = 0.0048143357390188_DP
         casimir_omega(6) = 0.0144431928875407_DP ;  casimir_omega_weight(6) = 0.0061715781641159_DP
         casimir_omega(7) = 0.0213152476270855_DP ;  casimir_omega_weight(7) = 0.0075825691344524_DP
         casimir_omega(8) = 0.0296306921476553_DP ;  casimir_omega_weight(8) = 0.0090606859365837_DP
         casimir_omega(9) = 0.0394638665797801_DP ;  casimir_omega_weight(9) = 0.0106205974538734_DP
         casimir_omega(10) = 0.0509045566078159_DP ;  casimir_omega_weight(10) = 0.0122785836212843_DP
         casimir_omega(11) = 0.0640597878064507_DP ;  casimir_omega_weight(11) = 0.0140529126257697_DP
         casimir_omega(12) = 0.0790560316983307_DP ;  casimir_omega_weight(12) = 0.0159642899371684_DP
         casimir_omega(13) = 0.0960419045514988_DP ;  casimir_omega_weight(13) = 0.0180363983902296_DP
         casimir_omega(14) = 0.1151914618576148_DP ;  casimir_omega_weight(14) = 0.0202965541416653_DP
         casimir_omega(15) = 0.1367082196959489_DP ;  casimir_omega_weight(15) = 0.0227765105384920_DP
         casimir_omega(16) = 0.1608300707911366_DP ;  casimir_omega_weight(16) = 0.0255134515359063_DP
         casimir_omega(17) = 0.1878353108533617_DP ;  casimir_omega_weight(17) = 0.0285512292384972_DP
         casimir_omega(18) = 0.2180500536971317_DP ;  casimir_omega_weight(18) = 0.0319419177142580_DP
         casimir_omega(19) = 0.2518573971862026_DP ;  casimir_omega_weight(19) = 0.0357477793138129_DP
         casimir_omega(20) = 0.2897088140233422_DP ;  casimir_omega_weight(20) = 0.0400437730245732_DP
         casimir_omega(21) = 0.3321383928663959_DP ;  casimir_omega_weight(21) = 0.0449207808730042_DP
         casimir_omega(22) = 0.3797807621445858_DP ;  casimir_omega_weight(22) = 0.0504897939548313_DP
         casimir_omega(23) = 0.4333938144821661_DP ;  casimir_omega_weight(23) = 0.0568873931884049_DP
         casimir_omega(24) = 0.4938877480331296_DP ;  casimir_omega_weight(24) = 0.0642829948586299_DP
         casimir_omega(25) = 0.5623625034184802_DP ;  casimir_omega_weight(25) = 0.0728885283243630_DP
         casimir_omega(26) = 0.6401564788043962_DP ;  casimir_omega_weight(26) = 0.0829715056635573_DP
         casimir_omega(27) = 0.7289105701319233_DP ;  casimir_omega_weight(27) = 0.0948728827932959_DP
         casimir_omega(28) = 0.8306532949256337_DP ;  casimir_omega_weight(28) = 0.1090317835941886_DP
         casimir_omega(29) = 0.9479153129482233_DP ;  casimir_omega_weight(29) = 0.1260202032539090_DP
         casimir_omega(30) = 1.0838855360657196_DP ;  casimir_omega_weight(30) = 0.1465924617652154_DP
         casimir_omega(31) = 1.2426270191800057_DP ;  casimir_omega_weight(31) = 0.1717568534391769_DP
         casimir_omega(32) = 1.4293802922685115_DP ;  casimir_omega_weight(32) = 0.2028813599060202_DP
         casimir_omega(33) = 1.6509970710671540_DP ;  casimir_omega_weight(33) = 0.2418527842408025_DP
         casimir_omega(34) = 1.9165725462612440_DP ;  casimir_omega_weight(34) = 0.2913217002272561_DP
         casimir_omega(35) = 2.2383873751291028_DP ;  casimir_omega_weight(35) = 0.3550889925809060_DP
         casimir_omega(36) = 2.6333456817788430_DP ;  casimir_omega_weight(36) = 0.4387331340128868_DP
         casimir_omega(37) = 3.1252316291027387_DP ;  casimir_omega_weight(37) = 0.5506608905070951_DP
         casimir_omega(38) = 3.7483638176600662_DP ;  casimir_omega_weight(38) = 0.7039321371494273_DP
         casimir_omega(39) = 4.5537322360641728_DP ;  casimir_omega_weight(39) = 0.9195642653828473_DP
         casimir_omega(40) = 5.6197501166831572_DP ;  casimir_omega_weight(40) = 1.2328148448917922_DP
         casimir_omega(41) = 7.0720584558578707_DP ;  casimir_omega_weight(41) = 1.7058367052258379_DP
         casimir_omega(42) = 9.1222688296957735_DP ;  casimir_omega_weight(42) = 2.4550038681676369_DP
         casimir_omega(43) = 12.1495643168257086_DP ;  casimir_omega_weight(43) = 3.7151810694301570_DP
         casimir_omega(44) = 16.8893182147480623_DP ;  casimir_omega_weight(44) = 6.0081133110724139_DP
         casimir_omega(45) = 24.9252366012882547_DP ;  casimir_omega_weight(45) = 10.6505567807398069_DP
         casimir_omega(46) = 40.2046577118162105_DP ;  casimir_omega_weight(46) = 21.6165614019983749_DP
         casimir_omega(47) = 74.9885690837269010_DP ;  casimir_omega_weight(47) = 54.6472434949705601_DP
         casimir_omega(48) = 184.8763448737379349_DP ;  casimir_omega_weight(48) = 210.0655489489396643_DP
         casimir_omega(49) = 975.8093420151944883_DP ;  casimir_omega_weight(49) = 2505.2682777467657615_DP
return
endsubroutine gauss_legendre_grid48

subroutine gauss_legendre_grid49()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0003541637479512_DP ;  casimir_omega_weight(2) = 0.0009092561655811_DP
         casimir_omega(3) = 0.0018692065536502_DP ;  casimir_omega_weight(3) = 0.0021237004147582_DP
         casimir_omega(4) = 0.0046077471062333_DP ;  casimir_omega_weight(4) = 0.0033571402057641_DP
         casimir_omega(5) = 0.0085926586133562_DP ;  casimir_omega_weight(5) = 0.0046181286773661_DP
         casimir_omega(6) = 0.0138566918681216_DP ;  casimir_omega_weight(6) = 0.0059172037568903_DP
         casimir_omega(7) = 0.0204435226310256_DP ;  casimir_omega_weight(7) = 0.0072656782359208_DP
         casimir_omega(8) = 0.0284085745709226_DP ;  casimir_omega_weight(8) = 0.0086757690810704_DP
         casimir_omega(9) = 0.0378200346786067_DP ;  casimir_omega_weight(9) = 0.0101608284650075_DP
         casimir_omega(10) = 0.0487601232265246_DP ;  casimir_omega_weight(10) = 0.0117356207133882_DP
         casimir_omega(11) = 0.0613266642430135_DP ;  casimir_omega_weight(11) = 0.0134166478510560_DP
         casimir_omega(12) = 0.0756350103816467_DP ;  casimir_omega_weight(12) = 0.0152225351156407_DP
         casimir_omega(13) = 0.0918203896294846_DP ;  casimir_omega_weight(13) = 0.0171744920732320_DP
         casimir_omega(14) = 0.1100407590958060_DP ;  casimir_omega_weight(14) = 0.0192968694411339_DP
         casimir_omega(15) = 0.1304802739388656_DP ;  casimir_omega_weight(15) = 0.0216178373966925_DP
         casimir_omega(16) = 0.1533535088337794_DP ;  casimir_omega_weight(16) = 0.0241702186396888_DP
         casimir_omega(17) = 0.1789106074111569_DP ;  casimir_omega_weight(17) = 0.0269925194839584_DP
         casimir_omega(18) = 0.2074435847760620_DP ;  casimir_omega_weight(18) = 0.0301302157406088_DP
         casimir_omega(19) = 0.2392940736665017_DP ;  casimir_omega_weight(19) = 0.0336373684755950_DP
         casimir_omega(20) = 0.2748628917759416_DP ;  casimir_omega_weight(20) = 0.0375786698216722_DP
         casimir_omega(21) = 0.3146219243436720_DP ;  casimir_omega_weight(21) = 0.0420320537187229_DP
         casimir_omega(22) = 0.3591289738417530_DP ;  casimir_omega_weight(22) = 0.0470920548844522_DP
         casimir_omega(23) = 0.4090464440584484_DP ;  casimir_omega_weight(23) = 0.0528741676223244_DP
         casimir_omega(24) = 0.4651650232723021_DP ;  casimir_omega_weight(24) = 0.0595205534905730_DP
         casimir_omega(25) = 0.5284339461758133_DP ;  casimir_omega_weight(25) = 0.0672075874560007_DP
         casimir_omega(26) = 0.5999999999999994_DP ;  casimir_omega_weight(26) = 0.0761559376857480_DP
         casimir_omega(27) = 0.6812582776054769_DP ;  casimir_omega_weight(27) = 0.0866441787164474_DP
         casimir_omega(28) = 0.7739188932725494_DP ;  casimir_omega_weight(28) = 0.0990273958268541_DP
         casimir_omega(29) = 0.8800956596228452_DP ;  casimir_omega_weight(29) = 0.1137629384303593_DP
         casimir_omega(30) = 1.0024253853676266_DP ;  casimir_omega_weight(30) = 0.1314465685135775_DP
         casimir_omega(31) = 1.1442304942701955_DP ;  casimir_omega_weight(31) = 0.1528639738063222_DP
         casimir_omega(32) = 1.3097439151351831_DP ;  casimir_omega_weight(32) = 0.1790654017346609_DP
         casimir_omega(33) = 1.5044250552636886_DP ;  casimir_omega_weight(33) = 0.2114757760292394_DP
         casimir_omega(34) = 1.7354115837740876_DP ;  casimir_omega_weight(34) = 0.2520604600730857_DP
         casimir_omega(35) = 2.0121780659582642_DP ;  casimir_omega_weight(35) = 0.3035804105552688_DP
         casimir_omega(36) = 2.3475172021672193_DP ;  casimir_omega_weight(36) = 0.3699948209096011_DP
         casimir_omega(37) = 2.7590377390583378_DP ;  casimir_omega_weight(37) = 0.4571145309078970_DP
         casimir_omega(38) = 3.2715150545859921_DP ;  casimir_omega_weight(38) = 0.5736965048386101_DP
         casimir_omega(39) = 3.9206978041879261_DP ;  casimir_omega_weight(39) = 0.7333446702990176_DP
         casimir_omega(40) = 4.7597005432203519_DP ;  casimir_omega_weight(40) = 0.9579519893434009_DP
         casimir_omega(41) = 5.8702035149582121_DP ;  casimir_omega_weight(41) = 1.2842448606390200_DP
         casimir_omega(42) = 7.3830822438153278_DP ;  casimir_omega_weight(42) = 1.7769654212447696_DP
         casimir_omega(43) = 9.5187644077872218_DP ;  casimir_omega_weight(43) = 2.5573358979772851_DP
         casimir_omega(44) = 12.6722303190978991_DP ;  casimir_omega_weight(44) = 3.8700056462235248_DP
         casimir_omega(45) = 17.6094896411664514_DP ;  casimir_omega_weight(45) = 6.2584559393559465_DP
         casimir_omega(46) = 25.9802269853608756_DP ;  casimir_omega_weight(46) = 11.0943000093914677_DP
         casimir_omega(47) = 41.8962298165124452_DP ;  casimir_omega_weight(47) = 22.5171497082907344_DP
         casimir_omega(48) = 78.1292878493724032_DP ;  casimir_omega_weight(48) = 56.9239082440141857_DP
         casimir_omega(49) = 192.5950876306286546_DP ;  casimir_omega_weight(49) = 218.8170519105287894_DP
         casimir_omega(50) = 1016.4789651194063254_DP ;  casimir_omega_weight(50) = 2609.6396696877532122_DP
return
endsubroutine gauss_legendre_grid49

subroutine gauss_legendre_grid50()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0003402715391351_DP ;  casimir_omega_weight(2) = 0.0008735767680645_DP
         casimir_omega(3) = 0.0017957677688414_DP ;  casimir_omega_weight(3) = 0.0020400969539748_DP
         casimir_omega(4) = 0.0044261875711759_DP ;  casimir_omega_weight(4) = 0.0032242131694547_DP
         casimir_omega(5) = 0.0082526549006076_DP ;  casimir_omega_weight(5) = 0.0044337429832287_DP
         casimir_omega(6) = 0.0133053670095006_DP ;  casimir_omega_weight(6) = 0.0056783749381546_DP
         casimir_omega(7) = 0.0196245572195340_DP ;  casimir_omega_weight(7) = 0.0069684906443430_DP
         casimir_omega(8) = 0.0272612255388467_DP ;  casimir_omega_weight(8) = 0.0083152712828608_DP
         casimir_omega(9) = 0.0362780349553520_DP ;  casimir_omega_weight(9) = 0.0097308986562372_DP
         casimir_omega(10) = 0.0467504293459098_DP ;  casimir_omega_weight(10) = 0.0112287961003410_DP
         casimir_omega(11) = 0.0587680122663356_DP ;  casimir_omega_weight(11) = 0.0128239104881290_DP
         casimir_omega(12) = 0.0724362319612418_DP ;  casimir_omega_weight(12) = 0.0145330445200177_DP
         casimir_omega(13) = 0.0878784289915625_DP ;  casimir_omega_weight(13) = 0.0163752520787442_DP
         casimir_omega(14) = 0.1052383174008137_DP ;  casimir_omega_weight(14) = 0.0183723130196705_DP
         casimir_omega(15) = 0.1246829888626363_DP ;  casimir_omega_weight(15) = 0.0205493082606819_DP
         casimir_omega(16) = 0.1464065529203032_DP ;  casimir_omega_weight(16) = 0.0229353219152610_DP
         casimir_omega(17) = 0.1706345568853563_DP ;  casimir_omega_weight(17) = 0.0255643050101035_DP
         casimir_omega(18) = 0.1976293684619649_DP ;  casimir_omega_weight(18) = 0.0284761457552244_DP
         casimir_omega(19) = 0.2276967558025521_DP ;  casimir_omega_weight(19) = 0.0317180053821682_DP
         casimir_omega(20) = 0.2611939677624155_DP ;  casimir_omega_weight(20) = 0.0353459976433954_DP
         casimir_omega(21) = 0.2985397075858409_DP ;  casimir_omega_weight(21) = 0.0394273161956008_DP
         casimir_omega(22) = 0.3402265145510621_DP ;  casimir_omega_weight(22) = 0.0440429502064737_DP
         casimir_omega(23) = 0.3868362322273963_DP ;  casimir_omega_weight(23) = 0.0492911789346037_DP
         casimir_omega(24) = 0.4390594662321446_DP ;  casimir_omega_weight(24) = 0.0552921071295881_DP
         casimir_omega(25) = 0.4977202438800146_DP ;  casimir_omega_weight(25) = 0.0621936044961286_DP
         casimir_omega(26) = 0.5638075199954548_DP ;  casimir_omega_weight(26) = 0.0701791588072620_DP
         casimir_omega(27) = 0.6385157828382675_DP ;  casimir_omega_weight(27) = 0.0794783661720427_DP
         casimir_omega(28) = 0.7232978855623665_DP ;  casimir_omega_weight(28) = 0.0903810989821751_DP
         casimir_omega(29) = 0.8199344910824827_DP ;  casimir_omega_weight(29) = 0.1032568688456575_DP
         casimir_omega(30) = 0.9306263736649648_DP ;  casimir_omega_weight(30) = 0.1185816303748303_DP
         casimir_omega(31) = 1.0581185904191786_DP ;  casimir_omega_weight(31) = 0.1369754043181186_DP
         casimir_omega(32) = 1.2058697414530262_DP ;  casimir_omega_weight(32) = 0.1592558925291545_DP
         casimir_omega(33) = 1.3782860419175507_DP ;  casimir_omega_weight(33) = 0.1865161573480742_DP
         casimir_omega(34) = 1.5810501942863651_DP ;  casimir_omega_weight(34) = 0.2202392317584835_DP
         casimir_omega(35) = 1.8215916126316236_DP ;  casimir_omega_weight(35) = 0.2624706473105766_DP
         casimir_omega(36) = 2.1097719393491445_DP ;  casimir_omega_weight(36) = 0.3160840004731060_DP
         casimir_omega(37) = 2.4589063318495508_DP ;  casimir_omega_weight(37) = 0.3852000279737628_DP
         casimir_omega(38) = 2.8873225071353823_DP ;  casimir_omega_weight(38) = 0.4758666822825103_DP
         casimir_omega(39) = 3.4208072581481281_DP ;  casimir_omega_weight(39) = 0.5971982760546529_DP
         casimir_omega(40) = 4.0965684540692537_DP ;  casimir_omega_weight(40) = 0.7633538954099459_DP
         casimir_omega(41) = 4.9698885523562835_DP ;  casimir_omega_weight(41) = 0.9971199444715251_DP
         casimir_omega(42) = 6.1257814603033713_DP ;  casimir_omega_weight(42) = 1.3367216294598185_DP
         casimir_omega(43) = 7.7004640392996961_DP ;  casimir_omega_weight(43) = 1.8495432402455416_DP
         casimir_omega(44) = 9.9233599736881484_DP ;  casimir_omega_weight(44) = 2.6617541537782352_DP
         casimir_omega(45) = 13.2055691878931647_DP ;  casimir_omega_weight(45) = 4.0279880332395512_DP
         casimir_omega(46) = 18.3443629312390506_DP ;  casimir_omega_weight(46) = 6.5139060225791594_DP
         casimir_omega(47) = 27.0567508391872735_DP ;  casimir_omega_weight(47) = 11.5470979314911251_DP
         casimir_omega(48) = 43.6223257043616073_DP ;  casimir_omega_weight(48) = 23.4361163568820459_DP
         casimir_omega(49) = 81.3341039463354321_DP ;  casimir_omega_weight(49) = 59.2470348019846398_DP
         casimir_omega(50) = 200.4713561777892323_DP ;  casimir_omega_weight(50) = 227.7471565053051847_DP
         casimir_omega(51) = 1057.9785806214763397_DP ;  casimir_omega_weight(51) = 2716.1410898192470995_DP
return
endsubroutine gauss_legendre_grid50

subroutine gauss_legendre_grid55()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0002817063325806_DP ;  casimir_omega_weight(2) = 0.0007231756852760_DP
         casimir_omega(3) = 0.0014862792582862_DP ;  casimir_omega_weight(3) = 0.0016879211172511_DP
         casimir_omega(4) = 0.0036615267384160_DP ;  casimir_omega_weight(4) = 0.0026649545065384_DP
         casimir_omega(5) = 0.0068219728700001_DP ;  casimir_omega_weight(5) = 0.0036593638295202_DP
         casimir_omega(6) = 0.0109882172091306_DP ;  casimir_omega_weight(6) = 0.0046776743000900_DP
         casimir_omega(7) = 0.0161875969205868_DP ;  casimir_omega_weight(7) = 0.0057268259472585_DP
         casimir_omega(8) = 0.0224546065181317_DP ;  casimir_omega_weight(8) = 0.0068142082161809_DP
         casimir_omega(9) = 0.0298313970381577_DP ;  casimir_omega_weight(9) = 0.0079477641050404_DP
         casimir_omega(10) = 0.0383683908524614_DP ;  casimir_omega_weight(10) = 0.0091361157343046_DP
         casimir_omega(11) = 0.0481250313726308_DP ;  casimir_omega_weight(11) = 0.0103887088646427_DP
         casimir_omega(12) = 0.0591706879921860_DP ;  casimir_omega_weight(12) = 0.0117159795637604_DP
         casimir_omega(13) = 0.0715857407766485_DP ;  casimir_omega_weight(13) = 0.0131295479708551_DP
         casimir_omega(14) = 0.0854628750299021_DP ;  casimir_omega_weight(14) = 0.0146424454518194_DP
         casimir_omega(15) = 0.1009086228987581_DP ;  casimir_omega_weight(15) = 0.0162693829704242_DP
         casimir_omega(16) = 0.1180451979236539_DP ;  casimir_omega_weight(16) = 0.0180270704192689_DP
         casimir_omega(17) = 0.1370126793645646_DP ;  casimir_omega_weight(17) = 0.0199345991126018_DP
         casimir_omega(18) = 0.1579716168535163_DP ;  casimir_omega_weight(18) = 0.0220139028210735_DP
         casimir_omega(19) = 0.1811061432806937_DP ;  casimir_omega_weight(19) = 0.0242903168625776_DP
         casimir_omega(20) = 0.2066277059180269_DP ;  casimir_omega_weight(20) = 0.0267932601738932_DP
         casimir_omega(21) = 0.2347795541067915_DP ;  casimir_omega_weight(21) = 0.0295570724142885_DP
         casimir_omega(22) = 0.2658421583875168_DP ;  casimir_omega_weight(22) = 0.0326220476028624_DP
         casimir_omega(23) = 0.3001397834563233_DP ;  casimir_omega_weight(23) = 0.0360357184154973_DP
         casimir_omega(24) = 0.3380484995288755_DP ;  casimir_omega_weight(24) = 0.0398544622613407_DP
         casimir_omega(25) = 0.3800059987515602_DP ;  casimir_omega_weight(25) = 0.0441455233239985_DP
         casimir_omega(26) = 0.4265236924375540_DP ;  casimir_omega_weight(26) = 0.0489895763322191_DP
         casimir_omega(27) = 0.4782017112987825_DP ;  casimir_omega_weight(27) = 0.0544840014639708_DP
         casimir_omega(28) = 0.5357476289856359_DP ;  casimir_omega_weight(28) = 0.0607471006882464_DP
         casimir_omega(29) = 0.6000000000000000_DP ;  casimir_omega_weight(29) = 0.0679235717334725_DP
         casimir_omega(30) = 0.6719581768035264_DP ;  casimir_omega_weight(30) = 0.0761916783502344_DP
         casimir_omega(31) = 0.7528203925122937_DP ;  casimir_omega_weight(31) = 0.0857727322981484_DP
         casimir_omega(32) = 0.8440328318987965_DP ;  casimir_omega_weight(32) = 0.0969437608703504_DP
         casimir_omega(33) = 0.9473534659524159_DP ;  casimir_omega_weight(33) = 0.1100546166762360_DP
         casimir_omega(34) = 1.0649359500240865_DP ;  casimir_omega_weight(34) = 0.1255513622753253_DP
         casimir_omega(35) = 1.1994411265789022_DP ;  casimir_omega_weight(35) = 0.1440086422253782_DP
         casimir_omega(36) = 1.3541870190326610_DP ;  casimir_omega_weight(36) = 0.1661751231107073_DP
         casimir_omega(37) = 1.5333532826979082_DP ;  casimir_omega_weight(37) = 0.1930382489472412_DP
         casimir_omega(38) = 1.7422639350350173_DP ;  casimir_omega_weight(38) = 0.2259180621281415_DP
         casimir_omega(39) = 1.9877845857610770_DP ;  casimir_omega_weight(39) = 0.2666056300909108_DP
         casimir_omega(40) = 2.2788903929103941_DP ;  casimir_omega_weight(40) = 0.3175714261121123_DP
         casimir_omega(41) = 2.6274940514235809_DP ;  casimir_omega_weight(41) = 0.3822860835128047_DP
         casimir_omega(42) = 3.0496793290382809_DP ;  casimir_omega_weight(42) = 0.4657265605697762_DP
         casimir_omega(43) = 3.5675841138094717_DP ;  casimir_omega_weight(43) = 0.5751975456547638_DP
         casimir_omega(44) = 4.2123553633556288_DP ;  casimir_omega_weight(44) = 0.7217073332722904_DP
         casimir_omega(45) = 5.0289344790496751_DP ;  casimir_omega_weight(45) = 0.9223573824706177_DP
         casimir_omega(46) = 6.0840935303564629_DP ;  casimir_omega_weight(46) = 1.2046693706700946_DP
         casimir_omega(47) = 7.4805146039808204_DP ;  casimir_omega_weight(47) = 1.6148122123128803_DP
         casimir_omega(48) = 9.3827234346187129_DP ;  casimir_omega_weight(48) = 2.2341736334806370_DP
         casimir_omega(49) = 12.0678223530570641_DP ;  casimir_omega_weight(49) = 3.2151429314875606_DP
         casimir_omega(50) = 16.0323450651120076_DP ;  casimir_omega_weight(50) = 4.8652706240529069_DP
         casimir_omega(51) = 22.2392490847214575_DP ;  casimir_omega_weight(51) = 7.8677711911614603_DP
         casimir_omega(52) = 32.7623665557739301_DP ;  casimir_omega_weight(52) = 13.9469103250642341_DP
         casimir_omega(53) = 52.7706584092579618_DP ;  casimir_omega_weight(53) = 28.3066265906746146_DP
         casimir_omega(54) = 98.3196425204150160_DP ;  casimir_omega_weight(54) = 71.5595960742252402_DP
         casimir_omega(55) = 242.2155849871174951_DP ;  casimir_omega_weight(55) = 275.0767048303795832_DP
         casimir_omega(56) = 1277.9265439374362359_DP ;  casimir_omega_weight(56) = 3280.5986137351833349_DP
return
endsubroutine gauss_legendre_grid55

subroutine gauss_legendre_grid60()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0002370566545418_DP ;  casimir_omega_weight(2) = 0.0006085241548422_DP
         casimir_omega(3) = 0.0012504430575928_DP ;  casimir_omega_weight(3) = 0.0014197181192234_DP
         casimir_omega(4) = 0.0030793547801650_DP ;  casimir_omega_weight(4) = 0.0022397939642519_DP
         casimir_omega(5) = 0.0057341209517736_DP ;  casimir_omega_weight(5) = 0.0030721563671044_DP
         casimir_omega(6) = 0.0092292786881034_DP ;  casimir_omega_weight(6) = 0.0039213602088471_DP
         casimir_omega(7) = 0.0135840589923602_DP ;  casimir_omega_weight(7) = 0.0047922257110839_DP
         casimir_omega(8) = 0.0188226428954021_DP ;  casimir_omega_weight(8) = 0.0056898403216213_DP
         casimir_omega(9) = 0.0249744554133460_DP ;  casimir_omega_weight(9) = 0.0066196153491969_DP
         casimir_omega(10) = 0.0320745234657117_DP ;  casimir_omega_weight(10) = 0.0075873557179919_DP
         casimir_omega(11) = 0.0401639083704970_DP ;  casimir_omega_weight(11) = 0.0085993393690195_DP
         casimir_omega(12) = 0.0492902229220804_DP ;  casimir_omega_weight(12) = 0.0096624073085802_DP
         casimir_omega(13) = 0.0595082446752810_DP ;  casimir_omega_weight(13) = 0.0107840663697600_DP
         casimir_omega(14) = 0.0708806394354815_DP ;  casimir_omega_weight(14) = 0.0119726073553746_DP
         casimir_omega(15) = 0.0834788119164447_DP ;  casimir_omega_weight(15) = 0.0132372418298672_DP
         casimir_omega(16) = 0.0973839041387285_DP ;  casimir_omega_weight(16) = 0.0145882615393668_DP
         casimir_omega(17) = 0.1126879665459734_DP ;  casimir_omega_weight(17) = 0.0160372253225004_DP
         casimir_omega(18) = 0.1294953322128241_DP ;  casimir_omega_weight(18) = 0.0175971794853605_DP
         casimir_omega(19) = 0.1479242311666076_DP ;  casimir_omega_weight(19) = 0.0192829190196078_DP
         casimir_omega(20) = 0.1681086900814500_DP ;  casimir_omega_weight(20) = 0.0211112988300882_DP
         casimir_omega(21) = 0.1902007728648993_DP ;  casimir_omega_weight(21) = 0.0231016064225662_DP
         casimir_omega(22) = 0.2143732305130889_DP ;  casimir_omega_weight(22) = 0.0252760104365700_DP
         casimir_omega(23) = 0.2408226448091120_DP ;  casimir_omega_weight(23) = 0.0276601031989754_DP
         casimir_omega(24) = 0.2697731709691586_DP ;  casimir_omega_weight(24) = 0.0302835603997013_DP
         casimir_omega(25) = 0.3014810105173619_DP ;  casimir_omega_weight(25) = 0.0331809474318927_DP
         casimir_omega(26) = 0.3362397792556048_DP ;  casimir_omega_weight(26) = 0.0363927104174012_DP
         casimir_omega(27) = 0.3743869785678427_DP ;  casimir_omega_weight(27) = 0.0399664011771694_DP
         casimir_omega(28) = 0.4163118346965068_DP ;  casimir_omega_weight(28) = 0.0439582004147339_DP
         casimir_omega(29) = 0.4624648444909915_DP ;  casimir_omega_weight(29) = 0.0484348235813609_DP
         casimir_omega(30) = 0.5133694635977094_DP ;  casimir_omega_weight(30) = 0.0534759213051174_DP
         casimir_omega(31) = 0.5696365027143087_DP ;  casimir_omega_weight(31) = 0.0591771237966982_DP
         casimir_omega(32) = 0.6319819714582998_DP ;  casimir_omega_weight(32) = 0.0656539305049170_DP
         casimir_omega(33) = 0.7012493448229440_DP ;  casimir_omega_weight(33) = 0.0730467186657657_DP
         casimir_omega(34) = 0.7784375489043526_DP ;  casimir_omega_weight(34) = 0.0815272464478622_DP
         casimir_omega(35) = 0.8647364066948555_DP ;  casimir_omega_weight(35) = 0.0913071719402846_DP
         casimir_omega(36) = 0.9615719044960436_DP ;  casimir_omega_weight(36) = 0.1026493192759898_DP
         casimir_omega(37) = 1.0706645144634506_DP ;  casimir_omega_weight(37) = 0.1158827302210296_DP
         casimir_omega(38) = 1.1941050594935170_DP ;  casimir_omega_weight(38) = 0.1314229945667831_DP
         casimir_omega(39) = 1.3344544185276173_DP ;  casimir_omega_weight(39) = 0.1498000369679065_DP
         casimir_omega(40) = 1.4948760332955977_DP ;  casimir_omega_weight(40) = 0.1716965835310302_DP
         casimir_omega(41) = 1.6793141528835596_DP ;  casimir_omega_weight(41) = 0.1980021570462496_DP
         casimir_omega(42) = 1.8927367884867110_DP ;  casimir_omega_weight(42) = 0.2298900245804481_DP
         casimir_omega(43) = 2.1414716861191234_DP ;  casimir_omega_weight(43) = 0.2689286834602631_DP
         casimir_omega(44) = 2.4336783579056145_DP ;  casimir_omega_weight(44) = 0.3172463519003216_DP
         casimir_omega(45) = 2.7800229849856213_DP ;  casimir_omega_weight(45) = 0.3777785855618269_DP
         casimir_omega(46) = 3.1946623143042552_DP ;  casimir_omega_weight(46) = 0.4546494265019519_DP
         casimir_omega(47) = 3.6967094632718851_DP ;  casimir_omega_weight(47) = 0.5537728740926077_DP
         casimir_omega(48) = 4.3124715330200170_DP ;  casimir_omega_weight(48) = 0.6838289531976264_DP
         casimir_omega(49) = 5.0789609527674626_DP ;  casimir_omega_weight(49) = 0.8578986553318921_DP
         casimir_omega(50) = 6.0495819018761861_DP ;  casimir_omega_weight(50) = 1.0963034298040943_DP
         casimir_omega(51) = 7.3036796885479722_DP ;  casimir_omega_weight(51) = 1.4317469838534909_DP
         casimir_omega(52) = 8.9632711209062403_DP ;  casimir_omega_weight(52) = 1.9190913771186584_DP
         casimir_omega(53) = 11.2238612176061512_DP ;  casimir_omega_weight(53) = 2.6550488794755971_DP
         casimir_omega(54) = 14.4147287314870312_DP ;  casimir_omega_weight(54) = 3.8207023130710702_DP
         casimir_omega(55) = 19.1259007568984032_DP ;  casimir_omega_weight(55) = 5.7815112319062454_DP
         casimir_omega(56) = 26.5016516935377560_DP ;  casimir_omega_weight(56) = 9.3493334137748025_DP
         casimir_omega(57) = 39.0062985598261562_DP ;  casimir_omega_weight(57) = 16.5730987475844138_DP
         casimir_omega(58) = 62.7820729677230673_DP ;  casimir_omega_weight(58) = 33.6366021627332259_DP
         casimir_omega(59) = 116.9076075023437653_DP ;  casimir_omega_weight(59) = 85.0337074979211991_DP
         casimir_omega(60) = 287.8979556998091880_DP ;  casimir_omega_weight(60) = 326.8712971074673987_DP
         casimir_omega(61) = 1518.6243166040380856_DP ;  casimir_omega_weight(61) = 3898.3068438659411186_DP
return
endsubroutine gauss_legendre_grid60

subroutine gauss_legendre_grid65()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0002022390221103_DP ;  casimir_omega_weight(2) = 0.0005191272578275_DP
         casimir_omega(3) = 0.0010666083055823_DP ;  casimir_omega_weight(3) = 0.0012107504552676_DP
         casimir_omega(4) = 0.0026258589647980_DP ;  casimir_omega_weight(4) = 0.0019089827863601_DP
         casimir_omega(5) = 0.0048875481922555_DP ;  casimir_omega_weight(5) = 0.0026161511712775_DP
         casimir_omega(6) = 0.0078622299409582_DP ;  casimir_omega_weight(6) = 0.0033355306636934_DP
         casimir_omega(7) = 0.0115638296762154_DP ;  casimir_omega_weight(7) = 0.0040705760888555_DP
         casimir_omega(8) = 0.0160098090569687_DP ;  casimir_omega_weight(8) = 0.0048249089175743_DP
         casimir_omega(9) = 0.0212213473109338_DP ;  casimir_omega_weight(9) = 0.0056023490998817_DP
         casimir_omega(10) = 0.0272235596101113_DP ;  casimir_omega_weight(10) = 0.0064069558338184_DP
         casimir_omega(11) = 0.0340457589556123_DP ;  casimir_omega_weight(11) = 0.0072430737085089_DP
         casimir_omega(12) = 0.0417217669059246_DP ;  casimir_omega_weight(12) = 0.0081153843598864_DP
         casimir_omega(13) = 0.0502902790757065_DP ;  casimir_omega_weight(13) = 0.0090289645364971_DP
         casimir_omega(14) = 0.0597952924025048_DP ;  casimir_omega_weight(14) = 0.0099893517947838_DP
         casimir_omega(15) = 0.0702866025309247_DP ;  casimir_omega_weight(15) = 0.0110026193093598_DP
         casimir_omega(16) = 0.0818203812924776_DP ;  casimir_omega_weight(16) = 0.0120754615779915_DP
         casimir_omega(17) = 0.0944598462097070_DP ;  casimir_omega_weight(17) = 0.0132152931539327_DP
         casimir_omega(18) = 0.1082760362947880_DP ;  casimir_omega_weight(18) = 0.0144303629712587_DP
         casimir_omega(19) = 0.1233487112366313_DP ;  casimir_omega_weight(19) = 0.0157298873644200_DP
         casimir_omega(20) = 0.1397673944923191_DP ;  casimir_omega_weight(20) = 0.0171242055486683_DP
         casimir_omega(21) = 0.1576325849650539_DP ;  casimir_omega_weight(21) = 0.0186249621582705_DP
         casimir_omega(22) = 0.1770571670480799_DP ;  casimir_omega_weight(22) = 0.0202453224795504_DP
         casimir_omega(23) = 0.1981680550801065_DP ;  casimir_omega_weight(23) = 0.0220002273243829_DP
         casimir_omega(24) = 0.2211081159979182_DP ;  casimir_omega_weight(24) = 0.0239066961439801_DP
         casimir_omega(25) = 0.2460384235799407_DP ;  casimir_omega_weight(25) = 0.0259841890839007_DP
         casimir_omega(26) = 0.2731409096617157_DP ;  casimir_omega_weight(26) = 0.0282550413639739_DP
         casimir_omega(27) = 0.3026214927375890_DP ;  casimir_omega_weight(27) = 0.0307449868110422_DP
         casimir_omega(28) = 0.3347137833187116_DP ;  casimir_omega_weight(28) = 0.0334837918195329_DP
         casimir_omega(29) = 0.3696834894537505_DP ;  casimir_omega_weight(29) = 0.0365060267920146_DP
         casimir_omega(30) = 0.4078336764767240_DP ;  casimir_omega_weight(30) = 0.0398520096645846_DP
         casimir_omega(31) = 0.4495110743913769_DP ;  casimir_omega_weight(31) = 0.0435689660623011_DP
         casimir_omega(32) = 0.4951136771216686_DP ;  casimir_omega_weight(32) = 0.0477124638054141_DP
         casimir_omega(33) = 0.5450999439467686_DP ;  casimir_omega_weight(33) = 0.0523481970809839_DP
         casimir_omega(34) = 0.6000000000000000_DP ;  casimir_omega_weight(34) = 0.0575542192733591_DP
         casimir_omega(35) = 0.6604293469440448_DP ;  casimir_omega_weight(35) = 0.0634237555806261_DP
         casimir_omega(36) = 0.7271057468920087_DP ;  casimir_omega_weight(36) = 0.0700687705356357_DP
         casimir_omega(37) = 0.8008701464973420_DP ;  casimir_omega_weight(37) = 0.0776245263374157_DP
         casimir_omega(38) = 0.8827127840693315_DP ;  casimir_omega_weight(38) = 0.0862554527269211_DP
         casimir_omega(39) = 0.9738059996456450_DP ;  casimir_omega_weight(39) = 0.0961627687669173_DP
         casimir_omega(40) = 1.0755457884959911_DP ;  casimir_omega_weight(40) = 0.1075944674799448_DP
         casimir_omega(41) = 1.1896048649530830_DP ;  casimir_omega_weight(41) = 0.1208585204985713_DP
         casimir_omega(42) = 1.3180010290141408_DP ;  casimir_omega_weight(42) = 0.1363405197657032_DP
         casimir_omega(43) = 1.4631860941144084_DP ;  casimir_omega_weight(43) = 0.1545275066438958_DP
         casimir_omega(44) = 1.6281627581838296_DP ;  casimir_omega_weight(44) = 0.1760405408782599_DP
         casimir_omega(45) = 1.8166399213762063_DP ;  casimir_omega_weight(45) = 0.2016797874948592_DP
         casimir_omega(46) = 2.0332416134402629_DP ;  casimir_omega_weight(46) = 0.2324878050927016_DP
         casimir_omega(47) = 2.2837917685598370_DP ;  casimir_omega_weight(47) = 0.2698397369822229_DP
         casimir_omega(48) = 2.5757080276672366_DP ;  casimir_omega_weight(48) = 0.3155739853299728_DP
         casimir_omega(49) = 2.9185550168366090_DP ;  casimir_omega_weight(49) = 0.3721850128910777_DP
         casimir_omega(50) = 3.3248354143651739_DP ;  casimir_omega_weight(50) = 0.4431135779514409_DP
         casimir_omega(51) = 3.8111431941226797_DP ;  casimir_omega_weight(51) = 0.5331934846699958_DP
         casimir_omega(52) = 4.3998817203397396_DP ;  casimir_omega_weight(52) = 0.6493565762270053_DP
         casimir_omega(53) = 5.1218864909796649_DP ;  casimir_omega_weight(53) = 0.8017767992300343_DP
         casimir_omega(54) = 6.0205408408525622_DP ;  casimir_omega_weight(54) = 1.0057865433503750_DP
         casimir_omega(55) = 7.1584410867567554_DP ;  casimir_omega_weight(55) = 1.2852048526442317_DP
         casimir_omega(56) = 8.6285895036932772_DP ;  casimir_omega_weight(56) = 1.6783642088803308_DP
         casimir_omega(57) = 10.5740042531980460_DP ;  casimir_omega_weight(57) = 2.2495692429665377_DP
         casimir_omega(58) = 13.2238401280297371_DP ;  casimir_omega_weight(58) = 3.1121778660529218_DP
         casimir_omega(59) = 16.9640501484333370_DP ;  casimir_omega_weight(59) = 4.4784400201798809_DP
         casimir_omega(60) = 22.4862144650809661_DP ;  casimir_omega_weight(60) = 6.7767164685723973_DP
         casimir_omega(61) = 31.1315550366891713_DP ;  casimir_omega_weight(61) = 10.9585982403289606_DP
         casimir_omega(62) = 45.7885361663849295_DP ;  casimir_omega_weight(62) = 19.4256677272917138_DP
         casimir_omega(63) = 73.6565627261627611_DP ;  casimir_omega_weight(63) = 39.4260466124224322_DP
         casimir_omega(64) = 137.0979952945416187_DP ;  casimir_omega_weight(64) = 99.6693716495570499_DP
         casimir_omega(65) = 337.5184668222261166_DP ;  casimir_omega_weight(65) = 383.1309349715792223_DP
         casimir_omega(66) = 1780.0718983078797919_DP ;  casimir_omega_weight(66) = 4569.2657809636566526_DP
return
endsubroutine gauss_legendre_grid65

subroutine gauss_legendre_grid70()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0001745649813764_DP ;  casimir_omega_weight(2) = 0.0004480770014238_DP
         casimir_omega(3) = 0.0009205345544637_DP ;  casimir_omega_weight(3) = 0.0010447669188318_DP
         casimir_omega(4) = 0.0022657069635432_DP ;  casimir_omega_weight(4) = 0.0016464984235335_DP
         casimir_omega(5) = 0.0042157467756660_DP ;  casimir_omega_weight(5) = 0.0022548867232735_DP
         casimir_omega(6) = 0.0067785028817737_DP ;  casimir_omega_weight(6) = 0.0028723460295096_DP
         casimir_omega(7) = 0.0099643087336436_DP ;  casimir_omega_weight(7) = 0.0035014185663255_DP
         casimir_omega(8) = 0.0137860931956307_DP ;  casimir_omega_weight(8) = 0.0041447548443641_DP
         casimir_omega(9) = 0.0182594969970058_DP ;  casimir_omega_weight(9) = 0.0048051318787598_DP
         casimir_omega(10) = 0.0234030113304936_DP ;  casimir_omega_weight(10) = 0.0054854780135976_DP
         casimir_omega(11) = 0.0292381429735033_DP ;  casimir_omega_weight(11) = 0.0061889010005954_DP
         casimir_omega(12) = 0.0357896089928089_DP ;  casimir_omega_weight(12) = 0.0069187191248614_DP
         casimir_omega(13) = 0.0430855642494803_DP ;  casimir_omega_weight(13) = 0.0076784957650361_DP
         casimir_omega(14) = 0.0511578654207043_DP ;  casimir_omega_weight(14) = 0.0084720779767617_DP
         casimir_omega(15) = 0.0600423759125749_DP ;  casimir_omega_weight(15) = 0.0093036398223839_DP
         casimir_omega(16) = 0.0697793168267342_DP ;  casimir_omega_weight(16) = 0.0101777313038259_DP
         casimir_omega(17) = 0.0804136700756015_DP ;  casimir_omega_weight(17) = 0.0110993339097249_DP
         casimir_omega(18) = 0.0919956408414075_DP ;  casimir_omega_weight(18) = 0.0120739239724980_DP
         casimir_omega(19) = 0.1045811878781766_DP ;  casimir_omega_weight(19) = 0.0131075452548236_DP
         casimir_omega(20) = 0.1182326317071362_DP ;  casimir_omega_weight(20) = 0.0142068924577933_DP
         casimir_omega(21) = 0.1330193526094884_DP ;  casimir_omega_weight(21) = 0.0153794076767136_DP
         casimir_omega(22) = 0.1490185925442667_DP ;  casimir_omega_weight(22) = 0.0166333922402223_DP
         casimir_omega(23) = 0.1663163777980980_DP ;  casimir_omega_weight(23) = 0.0179781368729441_DP
         casimir_omega(24) = 0.1850085824147225_DP ;  casimir_omega_weight(24) = 0.0194240737457245_DP
         casimir_omega(25) = 0.2052021563892091_DP ;  casimir_omega_weight(25) = 0.0209829547516961_DP
         casimir_omega(26) = 0.2270165474146049_DP ;  casimir_omega_weight(26) = 0.0226680613112033_DP
         casimir_omega(27) = 0.2505853508527850_DP ;  casimir_omega_weight(27) = 0.0244944522162161_DP
         casimir_omega(28) = 0.2760582298418578_DP ;  casimir_omega_weight(28) = 0.0264792575432452_DP
         casimir_omega(29) = 0.3036031564029310_DP ;  casimir_omega_weight(29) = 0.0286420285821125_DP
         casimir_omega(30) = 0.3334090355254185_DP ;  casimir_omega_weight(30) = 0.0310051561639027_DP
         casimir_omega(31) = 0.3656887880842788_DP ;  casimir_omega_weight(31) = 0.0335943728809332_DP
         casimir_omega(32) = 0.4006829858472447_DP ;  casimir_omega_weight(32) = 0.0364393586830251_DP
         casimir_omega(33) = 0.4386641537801632_DP ;  casimir_omega_weight(33) = 0.0395744744875736_DP
         casimir_omega(34) = 0.4799418826955366_DP ;  casimir_omega_weight(34) = 0.0430396551347375_DP
         casimir_omega(35) = 0.5248689307985325_DP ;  casimir_omega_weight(35) = 0.0468815017697647_DP
         casimir_omega(36) = 0.5738485382586160_DP ;  casimir_omega_weight(36) = 0.0511546252511900_DP
         casimir_omega(37) = 0.6273432378035589_DP ;  casimir_omega_weight(37) = 0.0559233074481541_DP
         casimir_omega(38) = 0.6858855208905166_DP ;  casimir_omega_weight(38) = 0.0612635676731017_DP
         casimir_omega(39) = 0.7500908192844151_DP ;  casimir_omega_weight(39) = 0.0672657489286345_DP
         casimir_omega(40) = 0.8206733942076653_DP ;  casimir_omega_weight(40) = 0.0740377758743830_DP
         casimir_omega(41) = 0.8984659012630140_DP ;  casimir_omega_weight(41) = 0.0817092873842961_DP
         casimir_omega(42) = 0.9844436354910395_DP ;  casimir_omega_weight(42) = 0.0904369169867099_DP
         casimir_omega(43) = 1.0797547805885841_DP ;  casimir_omega_weight(43) = 0.1004110927531150_DP
         casimir_omega(44) = 1.1857584231509812_DP ;  casimir_omega_weight(44) = 0.1118648667878014_DP
         casimir_omega(45) = 1.3040726958447473_DP ;  casimir_omega_weight(45) = 0.1250854821034270_DP
         casimir_omega(46) = 1.4366362549720426_DP ;  casimir_omega_weight(46) = 0.1404296698898841_DP
         casimir_omega(47) = 1.5857874859779466_DP ;  casimir_omega_weight(47) = 0.1583440871076123_DP
         casimir_omega(48) = 1.7543675287562956_DP ;  casimir_omega_weight(48) = 0.1793929221870204_DP
         casimir_omega(49) = 1.9458556749167992_DP ;  casimir_omega_weight(49) = 0.2042956258288241_DP
         casimir_omega(50) = 2.1645493051623963_DP ;  casimir_omega_weight(50) = 0.2339791437959637_DP
         casimir_omega(51) = 2.4158059330285262_DP ;  casimir_omega_weight(51) = 0.2696512359582444_DP
         casimir_omega(52) = 2.7063731174280306_DP ;  casimir_omega_weight(52) = 0.3129049621855994_DP
         casimir_omega(53) = 3.0448446829105937_DP ;  casimir_omega_weight(53) = 0.3658700676471867_DP
         casimir_omega(54) = 3.4423016921490022_DP ;  casimir_omega_weight(54) = 0.4314363426734113_DP
         casimir_omega(55) = 3.9132288954930896_DP ;  casimir_omega_weight(55) = 0.5135898585957677_DP
         casimir_omega(56) = 4.4768507600951715_DP ;  casimir_omega_weight(56) = 0.6179305260359963_DP
         casimir_omega(57) = 5.1591218769581664_DP ;  casimir_omega_weight(57) = 0.7524888264204401_DP
         casimir_omega(58) = 5.9957653994935125_DP ;  casimir_omega_weight(58) = 0.9290512057288070_DP
         casimir_omega(59) = 7.0370410696279926_DP ;  casimir_omega_weight(59) = 1.1653801458931810_DP
         casimir_omega(60) = 8.3554667618016136_DP ;  casimir_omega_weight(60) = 1.4890698836831036_DP
         casimir_omega(61) = 10.0587855003482822_DP ;  casimir_omega_weight(61) = 1.9445284140467274_DP
         casimir_omega(62) = 12.3126834808299890_DP ;  casimir_omega_weight(62) = 2.6062523595832796_DP
         casimir_omega(63) = 15.3826358034073731_DP ;  casimir_omega_weight(63) = 3.6055663649071148_DP
         casimir_omega(64) = 19.7157676391104602_DP ;  casimir_omega_weight(64) = 5.1883610820410926_DP
         casimir_omega(65) = 26.1132718959202634_DP ;  casimir_omega_weight(65) = 7.8508906516831649_DP
         casimir_omega(66) = 36.1289487934571412_DP ;  casimir_omega_weight(66) = 12.6955693032784875_DP
         casimir_omega(67) = 53.1090723540126390_DP ;  casimir_omega_weight(67) = 22.5046202336423491_DP
         casimir_omega(68) = 85.3941233088265079_DP ;  casimir_omega_weight(68) = 45.6749622643592019_DP
         casimir_omega(69) = 158.8908035296126968_DP ;  casimir_omega_weight(69) = 115.4665902233183203_DP
         casimir_omega(70) = 391.0771173709248956_DP ;  casimir_omega_weight(70) = 443.8556194984579406_DP
         casimir_omega(71) = 2062.2692888437086367_DP ;  casimir_omega_weight(71) = 5293.4754255257776094_DP
return
endsubroutine gauss_legendre_grid70

subroutine gauss_legendre_grid75()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0001522059404776_DP ;  casimir_omega_weight(2) = 0.0003906756503321_DP
         casimir_omega(3) = 0.0008025433370984_DP ;  casimir_omega_weight(3) = 0.0009107327666231_DP
         casimir_omega(4) = 0.0019749182031068_DP ;  casimir_omega_weight(4) = 0.0014347191706923_DP
         casimir_omega(5) = 0.0036736661259948_DP ;  casimir_omega_weight(5) = 0.0019637678304156_DP
         casimir_omega(6) = 0.0059047463681733_DP ;  casimir_omega_weight(6) = 0.0024996965867747_DP
         casimir_omega(7) = 0.0086759895399686_DP ;  casimir_omega_weight(7) = 0.0030444184318680_DP
         casimir_omega(8) = 0.0119971753219177_DP ;  casimir_omega_weight(8) = 0.0035999193122382_DP
         casimir_omega(9) = 0.0158801098581114_DP ;  casimir_omega_weight(9) = 0.0041682685667020_DP
         casimir_omega(10) = 0.0203387167372781_DP ;  casimir_omega_weight(10) = 0.0047516345525367_DP
         casimir_omega(11) = 0.0253891447379611_DP ;  casimir_omega_weight(11) = 0.0053523023988197_DP
         casimir_omega(12) = 0.0310498942185800_DP ;  casimir_omega_weight(12) = 0.0059726935470068_DP
         casimir_omega(13) = 0.0373419639971980_DP ;  casimir_omega_weight(13) = 0.0066153872307762_DP
         casimir_omega(14) = 0.0442890208001906_DP ;  casimir_omega_weight(14) = 0.0072831441902323_DP
         casimir_omega(15) = 0.0519175936943116_DP ;  casimir_omega_weight(15) = 0.0079789329919348_DP
         casimir_omega(16) = 0.0602572963214735_DP ;  casimir_omega_weight(16) = 0.0087059593931432_DP
         casimir_omega(17) = 0.0693410802295746_DP ;  casimir_omega_weight(17) = 0.0094676992615007_DP
         casimir_omega(18) = 0.0792055231451624_DP ;  casimir_omega_weight(18) = 0.0102679356462333_DP
         casimir_omega(19) = 0.0898911566787843_DP ;  casimir_omega_weight(19) = 0.0111108006979157_DP
         casimir_omega(20) = 0.1014428387096094_DP ;  casimir_omega_weight(20) = 0.0120008232548953_DP
         casimir_omega(21) = 0.1139101765845520_DP ;  casimir_omega_weight(21) = 0.0129429830601446_DP
         casimir_omega(22) = 0.1273480083160445_DP ;  casimir_omega_weight(22) = 0.0139427727481208_DP
         casimir_omega(23) = 0.1418169502051866_DP ;  casimir_omega_weight(23) = 0.0150062689540490_DP
         casimir_omega(24) = 0.1573840207943822_DP ;  casimir_omega_weight(24) = 0.0161402141563908_DP
         casimir_omega(25) = 0.1741233528164491_DP ;  casimir_omega_weight(25) = 0.0173521111778806_DP
         casimir_omega(26) = 0.1921170069183034_DP ;  casimir_omega_weight(26) = 0.0186503326549301_DP
         casimir_omega(27) = 0.2114559034747180_DP ;  casimir_omega_weight(27) = 0.0200442482565257_DP
         casimir_omega(28) = 0.2322408918686883_DP ;  casimir_omega_weight(28) = 0.0215443730137983_DP
         casimir_omega(29) = 0.2545839803216985_DP ;  casimir_omega_weight(29) = 0.0231625408380734_DP
         casimir_omega(30) = 0.2786097538633682_DP ;  casimir_omega_weight(30) = 0.0249121081941684_DP
         casimir_omega(31) = 0.3044570135299776_DP ;  casimir_omega_weight(31) = 0.0268081940031432_DP
         casimir_omega(32) = 0.3322806766223510_DP ;  casimir_omega_weight(32) = 0.0288679632344920_DP
         casimir_omega(33) = 0.3622539861510823_DP ;  casimir_omega_weight(33) = 0.0311109633898823_DP
         casimir_omega(34) = 0.3945710878563023_DP ;  casimir_omega_weight(34) = 0.0335595252812858_DP
         casimir_omega(35) = 0.4294500459330673_DP ;  casimir_omega_weight(35) = 0.0362392423006938_DP
         casimir_omega(36) = 0.4671363845007030_DP ;  casimir_omega_weight(36) = 0.0391795459452747_DP
         casimir_omega(37) = 0.5079072618112066_DP ;  casimir_omega_weight(37) = 0.0424143999400172_DP
         casimir_omega(38) = 0.5520764093619511_DP ;  casimir_omega_weight(38) = 0.0459831412103612_DP
         casimir_omega(39) = 0.5999999999999996_DP ;  casimir_omega_weight(39) = 0.0499315036345703_DP
         casimir_omega(40) = 0.6520836498267703_DP ;  casimir_omega_weight(40) = 0.0543128705419701_DP
         casimir_omega(41) = 0.7087908109764637_DP ;  casimir_omega_weight(41) = 0.0591898151315277_DP
         casimir_omega(42) = 0.7706528798538882_DP ;  casimir_omega_weight(42) = 0.0646360054919823_DP
         casimir_omega(43) = 0.8382814332172830_DP ;  casimir_omega_weight(43) = 0.0707385742817431_DP
         casimir_omega(44) = 0.9123831194927980_DP ;  casimir_omega_weight(44) = 0.0776010845883060_DP
         casimir_omega(45) = 0.9937778844753906_DP ;  casimir_omega_weight(45) = 0.0853472661821702_DP
         casimir_omega(46) = 1.0834214124619497_DP ;  casimir_omega_weight(46) = 0.0941257548297336_DP
         casimir_omega(47) = 1.1824329347057512_DP ;  casimir_omega_weight(47) = 0.1041161480951609_DP
         casimir_omega(48) = 1.2921299236944384_DP ;  casimir_omega_weight(48) = 0.1155368037681252_DP
         casimir_omega(49) = 1.4140716927478916_DP ;  casimir_omega_weight(49) = 0.1286549659953023_DP
         casimir_omega(50) = 1.5501146120449321_DP ;  casimir_omega_weight(50) = 0.1438000308529504_DP
         casimir_omega(51) = 1.7024826173417402_DP ;  casimir_omega_weight(51) = 0.1613810902115488_DP
         casimir_omega(52) = 1.8738580502302324_DP ;  casimir_omega_weight(52) = 0.1819103709010714_DP
         casimir_omega(53) = 2.0674998165208267_DP ;  casimir_omega_weight(53) = 0.2060348947813739_DP
         casimir_omega(54) = 2.2873986709891558_DP ;  casimir_omega_weight(54) = 0.2345797510094282_DP
         casimir_omega(55) = 2.5384835837968391_DP ;  casimir_omega_weight(55) = 0.2686080002339551_DP
         casimir_omega(56) = 2.8268993348256726_DP ;  casimir_omega_weight(56) = 0.3095047620177251_DP
         casimir_omega(57) = 3.1603848821424969_DP ;  casimir_omega_weight(57) = 0.3590970466343306_DP
         casimir_omega(58) = 3.5487965890873521_DP ;  casimir_omega_weight(58) = 0.4198273744598983_DP
         casimir_omega(59) = 4.0048433383321509_DP ;  casimir_omega_weight(59) = 0.4950099409398652_DP
         casimir_omega(60) = 4.5451375826432994_DP ;  casimir_omega_weight(60) = 0.5892162357961729_DP
         casimir_omega(61) = 5.1917275993986660_DP ;  casimir_omega_weight(61) = 0.7088686157758305_DP
         casimir_omega(62) = 5.9743802323853856_DP ;  casimir_omega_weight(62) = 0.8631769906312597_DP
         casimir_omega(63) = 6.9340655909375037_DP ;  casimir_omega_weight(63) = 1.0656588792910995_DP
         casimir_omega(64) = 8.1284253635711252_DP ;  casimir_omega_weight(64) = 1.3366855462782166_DP
         casimir_omega(65) = 9.6406284368710562_DP ;  casimir_omega_weight(65) = 1.7079040155123919_DP
         casimir_omega(66) = 11.5942423979201976_DP ;  casimir_omega_weight(66) = 2.2302445304645575_DP
         casimir_omega(67) = 14.1792881845971213_DP ;  casimir_omega_weight(67) = 2.9891451227383428_DP
         casimir_omega(68) = 17.7002317624184045_DP ;  casimir_omega_weight(68) = 4.1352182596685516_DP
         casimir_omega(69) = 22.6698683583797624_DP ;  casimir_omega_weight(69) = 5.9504688905683558_DP
         casimir_omega(70) = 30.0070633578479935_DP ;  casimir_omega_weight(70) = 9.0040366992157246_DP
         casimir_omega(71) = 41.4938259597419261_DP ;  casimir_omega_weight(71) = 14.5602490619222866_DP
         casimir_omega(72) = 60.9679023540134182_DP ;  casimir_omega_weight(72) = 25.8099582800999805_DP
         casimir_omega(73) = 97.9947517420359162_DP ;  casimir_omega_weight(73) = 52.3833506966992246_DP
         casimir_omega(74) = 182.2860305979529301_DP ;  casimir_omega_weight(74) = 132.4253643704815602_DP
         casimir_omega(75) = 448.5739066772051160_DP ;  casimir_omega_weight(75) = 509.0453514199984966_DP
         casimir_omega(76) = 2365.2164880703394374_DP ;  casimir_omega_weight(76) = 6070.9357778907797183_DP
return
endsubroutine gauss_legendre_grid75

subroutine gauss_legendre_grid80()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0001338830723078_DP ;  casimir_omega_weight(2) = 0.0003436383074826_DP
         casimir_omega(3) = 0.0007058701257852_DP ;  casimir_omega_weight(3) = 0.0008009412915707_DP
         casimir_omega(4) = 0.0017367503704767_DP ;  casimir_omega_weight(4) = 0.0012613646053483_DP
         casimir_omega(5) = 0.0032299032786565_DP ;  casimir_omega_weight(5) = 0.0017257073388961_DP
         casimir_omega(6) = 0.0051899355797046_DP ;  casimir_omega_weight(6) = 0.0021953635961398_DP
         casimir_omega(7) = 0.0076228902581809_DP ;  casimir_omega_weight(7) = 0.0026718008514492_DP
         casimir_omega(8) = 0.0105363030679594_DP ;  casimir_omega_weight(8) = 0.0031565373382687_DP
         casimir_omega(9) = 0.0139392555531132_DP ;  casimir_omega_weight(9) = 0.0036511478989218_DP
         casimir_omega(10) = 0.0178424365400675_DP ;  casimir_omega_weight(10) = 0.0041572740883146_DP
         casimir_omega(11) = 0.0222582145496690_DP ;  casimir_omega_weight(11) = 0.0046766357665017_DP
         casimir_omega(12) = 0.0272007223611805_DP ;  casimir_omega_weight(12) = 0.0052110437989930_DP
         casimir_omega(13) = 0.0326859548405648_DP ;  casimir_omega_weight(13) = 0.0057624139042236_DP
         casimir_omega(14) = 0.0387318812500342_DP ;  casimir_omega_weight(14) = 0.0063327817984668_DP
         casimir_omega(15) = 0.0453585734326655_DP ;  casimir_omega_weight(15) = 0.0069243198373922_DP
         casimir_omega(16) = 0.0525883514836556_DP ;  casimir_omega_weight(16) = 0.0075393553900736_DP
         casimir_omega(17) = 0.0604459487738210_DP ;  casimir_omega_weight(17) = 0.0081803912182204_DP
         casimir_omega(18) = 0.0689586984840571_DP ;  casimir_omega_weight(18) = 0.0088501281750806_DP
         casimir_omega(19) = 0.0781567441475556_DP ;  casimir_omega_weight(19) = 0.0095514905870567_DP
         casimir_omega(20) = 0.0880732770875454_DP ;  casimir_omega_weight(20) = 0.0102876547385217_DP
         casimir_omega(21) = 0.0987448040918477_DP ;  casimir_omega_weight(21) = 0.0110620809483955_DP
         casimir_omega(22) = 0.1102114491934746_DP ;  casimir_omega_weight(22) = 0.0118785498081130_DP
         casimir_omega(23) = 0.1225172940430110_DP ;  casimir_omega_weight(23) = 0.0127412032472771_DP
         casimir_omega(24) = 0.1357107620808955_DP ;  casimir_omega_weight(24) = 0.0136545912088921_DP
         casimir_omega(25) = 0.1498450525667780_DP ;  casimir_omega_weight(25) = 0.0146237248546327_DP
         casimir_omega(26) = 0.1649786315243007_DP ;  casimir_omega_weight(26) = 0.0156541373871653_DP
         casimir_omega(27) = 0.1811757878438931_DP ;  casimir_omega_weight(27) = 0.0167519537773031_DP
         casimir_omega(28) = 0.1985072641914661_DP ;  casimir_omega_weight(28) = 0.0179239709265172_DP
         casimir_omega(29) = 0.2170509740438506_DP ;  casimir_omega_weight(29) = 0.0191777500897732_DP
         casimir_omega(30) = 0.2368928181700402_DP ;  casimir_omega_weight(30) = 0.0205217237419592_DP
         casimir_omega(31) = 0.2581276162720160_DP ;  casimir_omega_weight(31) = 0.0219653195087639_DP
         casimir_omega(32) = 0.2808601723787875_DP ;  casimir_omega_weight(32) = 0.0235191043192153_DP
         casimir_omega(33) = 0.3052064960629326_DP ;  casimir_omega_weight(33) = 0.0251949525969427_DP
         casimir_omega(34) = 0.3312952057590653_DP ;  casimir_omega_weight(34) = 0.0270062431223398_DP
         casimir_omega(35) = 0.3592691455830052_DP ;  casimir_omega_weight(35) = 0.0289680902087612_DP
         casimir_omega(36) = 0.3892872532999511_DP ;  casimir_omega_weight(36) = 0.0310976160951698_DP
         casimir_omega(37) = 0.4215267247505370_DP ;  casimir_omega_weight(37) = 0.0334142730333785_DP
         casimir_omega(38) = 0.4561855294740127_DP ;  casimir_omega_weight(38) = 0.0359402255287216_DP
         casimir_omega(39) = 0.4934853439285725_DP ;  casimir_omega_weight(39) = 0.0387008056952053_DP
         casimir_omega(40) = 0.5336749831944303_DP ;  casimir_omega_weight(40) = 0.0417250578629519_DP
         casimir_omega(41) = 0.5770344301273437_DP ;  casimir_omega_weight(41) = 0.0450463926308094_DP
         casimir_omega(42) = 0.6238795836160996_DP ;  casimir_omega_weight(42) = 0.0487033757616758_DP
         casimir_omega(43) = 0.6745678762102363_DP ;  casimir_omega_weight(43) = 0.0527406840374694_DP
         casimir_omega(44) = 0.7295049476729885_DP ;  casimir_omega_weight(44) = 0.0572102689186832_DP
         casimir_omega(45) = 0.7891526073066905_DP ;  casimir_omega_weight(45) = 0.0621727802630738_DP
         casimir_omega(46) = 0.8540383773129728_DP ;  casimir_omega_weight(46) = 0.0676993173740227_DP
         casimir_omega(47) = 0.9247669861993026_DP ;  casimir_omega_weight(47) = 0.0738735945514127_DP
         casimir_omega(48) = 1.0020342810563627_DP ;  casimir_omega_weight(48) = 0.0807946348935923_DP
         casimir_omega(49) = 1.0866441582671438_DP ;  casimir_omega_weight(49) = 0.0885801418659057_DP
         casimir_omega(50) = 1.1795292847429066_DP ;  casimir_omega_weight(50) = 0.0973707466884176_DP
         casimir_omega(51) = 1.2817766112970939_DP ;  casimir_omega_weight(51) = 0.1073353960431572_DP
         casimir_omega(52) = 1.3946589876715507_DP ;  casimir_omega_weight(52) = 0.1186782364181158_DP
         casimir_omega(53) = 1.5196746055070101_DP ;  casimir_omega_weight(53) = 0.1316474795343942_DP
         casimir_omega(54) = 1.6585965650965913_DP ;  casimir_omega_weight(54) = 0.1465469139924357_DP
         casimir_omega(55) = 1.8135356480091789_DP ;  casimir_omega_weight(55) = 0.1637509859476252_DP
         casimir_omega(56) = 1.9870204748892117_DP ;  casimir_omega_weight(56) = 0.1837247434992768_DP
         casimir_omega(57) = 2.1821007767722507_DP ;  casimir_omega_weight(57) = 0.2070504830633304_DP
         casimir_omega(58) = 2.4024817225084343_DP ;  casimir_omega_weight(58) = 0.2344637415545699_DP
         casimir_omega(59) = 2.6527004526391820_DP ;  casimir_omega_weight(59) = 0.2669024897144135_DP
         casimir_omega(60) = 2.9383606846035732_DP ;  casimir_omega_weight(60) = 0.3055752331846397_DP
         casimir_omega(61) = 3.2664482922098679_DP ;  casimir_omega_weight(61) = 0.3520566059024088_DP
         casimir_omega(62) = 3.6457614485228560_DP ;  casimir_omega_weight(62) = 0.4084235989226063_DP
         casimir_omega(63) = 4.0875054489247375_DP ;  casimir_omega_weight(63) = 0.4774529368149332_DP
         casimir_omega(64) = 4.6061284144633614_DP ;  casimir_omega_weight(64) = 0.5629122946889830_DP
         casimir_omega(65) = 5.2205161627757484_DP ;  casimir_omega_weight(65) = 0.6699986832165482_DP
         casimir_omega(66) = 5.9557341278083280_DP ;  casimir_omega_weight(66) = 0.8060132423346604_DP
         casimir_omega(67) = 6.8456224590323425_DP ;  casimir_omega_weight(67) = 0.9814260977728830_DP
         casimir_omega(68) = 7.9367575467164304_DP ;  casimir_omega_weight(68) = 1.2116044127112116_DP
         casimir_omega(69) = 9.2946685877717847_DP ;  casimir_omega_weight(69) = 1.5197069224561899_DP
         casimir_omega(70) = 11.0139049556913218_DP ;  casimir_omega_weight(70) = 1.9417110305037160_DP
         casimir_omega(71) = 13.2349426320299504_DP ;  casimir_omega_weight(71) = 2.5355159622925849_DP
         casimir_omega(72) = 16.1738040217315344_DP ;  casimir_omega_weight(72) = 3.3982505739457074_DP
         casimir_omega(73) = 20.1766165283300474_DP ;  casimir_omega_weight(73) = 4.7011362430647079_DP
         casimir_omega(74) = 25.8263433530061342_DP ;  casimir_omega_weight(74) = 6.7647658019368109_DP
         casimir_omega(75) = 34.1675820900357792_DP ;  casimir_omega_weight(75) = 10.2361566414619460_DP
         casimir_omega(76) = 47.2261816459501915_DP ;  casimir_omega_weight(76) = 16.5526392298427183_DP
         casimir_omega(77) = 69.3650228353103131_DP ;  casimir_omega_weight(77) = 29.3416832712797273_DP
         casimir_omega(78) = 111.4584459475640301_DP ;  casimir_omega_weight(78) = 59.5512130114534557_DP
         casimir_omega(79) = 207.2836753742452345_DP ;  casimir_omega_weight(79) = 150.5456948955465180_DP
         casimir_omega(80) = 510.0088342732555589_DP ;  casimir_omega_weight(80) = 578.7001312471962819_DP
         casimir_omega(81) = 2688.9134958889662812_DP ;  casimir_omega_weight(81) = 6901.6468382858174664_DP
return
endsubroutine gauss_legendre_grid80

subroutine gauss_legendre_grid85()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0001186802631222_DP ;  casimir_omega_weight(2) = 0.0003046120417370_DP
         casimir_omega(3) = 0.0006256714020899_DP ;  casimir_omega_weight(3) = 0.0007098776993188_DP
         casimir_omega(4) = 0.0015392263033628_DP ;  casimir_omega_weight(4) = 0.0011176624981358_DP
         casimir_omega(5) = 0.0028620218221214_DP ;  casimir_omega_weight(5) = 0.0015285301393762_DP
         casimir_omega(6) = 0.0045976758173529_DP ;  casimir_omega_weight(6) = 0.0019435669746948_DP
         casimir_omega(7) = 0.0067509266886026_DP ;  casimir_omega_weight(7) = 0.0023639177243828_DP
         casimir_omega(8) = 0.0093276758280131_DP ;  casimir_omega_weight(8) = 0.0027907634647028_DP
         casimir_omega(9) = 0.0123350249405869_DP ;  casimir_omega_weight(9) = 0.0032253247060268_DP
         casimir_omega(10) = 0.0157813186821240_DP ;  casimir_omega_weight(10) = 0.0036688680670661_DP
         casimir_omega(11) = 0.0196761945912157_DP ;  casimir_omega_weight(11) = 0.0041227140545237_DP
         casimir_omega(12) = 0.0240306411712847_DP ;  casimir_omega_weight(12) = 0.0045882455659461_DP
         casimir_omega(13) = 0.0288570648253332_DP ;  casimir_omega_weight(13) = 0.0050669171003626_DP
         casimir_omega(14) = 0.0341693663849713_DP ;  casimir_omega_weight(14) = 0.0055602647523848_DP
         casimir_omega(15) = 0.0399830280702329_DP ;  casimir_omega_weight(15) = 0.0060699171000413_DP
         casimir_omega(16) = 0.0463152118383574_DP ;  casimir_omega_weight(16) = 0.0065976071185642_DP
         casimir_omega(17) = 0.0531848702219249_DP ;  casimir_omega_weight(17) = 0.0071451852724256_DP
         casimir_omega(18) = 0.0606128709196909_DP ;  casimir_omega_weight(18) = 0.0077146339596604_DP
         casimir_omega(19) = 0.0686221365895728_DP ;  casimir_omega_weight(19) = 0.0083080835072395_DP
         casimir_omega(20) = 0.0772378015060605_DP ;  casimir_omega_weight(20) = 0.0089278299451325_DP
         casimir_omega(21) = 0.0864873869883755_DP ;  casimir_omega_weight(21) = 0.0095763548204101_DP
         casimir_omega(22) = 0.0964009977863587_DP ;  casimir_omega_weight(22) = 0.0102563473524377_DP
         casimir_omega(23) = 0.1070115419348644_DP ;  casimir_omega_weight(23) = 0.0109707292769701_DP
         casimir_omega(24) = 0.1183549769620369_DP ;  casimir_omega_weight(24) = 0.0117226827820936_DP
         casimir_omega(25) = 0.1304705857715467_DP ;  casimir_omega_weight(25) = 0.0125156820042847_DP
         casimir_omega(26) = 0.1434012860246857_DP ;  casimir_omega_weight(26) = 0.0133535286301825_DP
         casimir_omega(27) = 0.1571939774384991_DP ;  casimir_omega_weight(27) = 0.0142403922416914_DP
         casimir_omega(28) = 0.1718999321068717_DP ;  casimir_omega_weight(28) = 0.0151808561515035_DP
         casimir_omega(29) = 0.1875752337620938_DP ;  casimir_omega_weight(29) = 0.0161799696070436_DP
         casimir_omega(30) = 0.2042812728483716_DP ;  casimir_omega_weight(30) = 0.0172433073974984_DP
         casimir_omega(31) = 0.2220853054046780_DP ;  casimir_omega_weight(31) = 0.0183770380868813_DP
         casimir_omega(32) = 0.2410610850871119_DP ;  casimir_omega_weight(32) = 0.0195880023228620_DP
         casimir_omega(33) = 0.2612895792433867_DP ;  casimir_omega_weight(33) = 0.0208838029452302_DP
         casimir_omega(34) = 0.2828597818367436_DP ;  casimir_omega_weight(34) = 0.0222729089502540_DP
         casimir_omega(35) = 0.3058696382684068_DP ;  casimir_omega_weight(35) = 0.0237647757716499_DP
         casimir_omega(36) = 0.3304270998470526_DP ;  casimir_omega_weight(36) = 0.0253699848327011_DP
         casimir_omega(37) = 0.3566513289005530_DP ;  casimir_omega_weight(37) = 0.0271004059292915_DP
         casimir_omega(38) = 0.3846740794443632_DP ;  casimir_omega_weight(38) = 0.0289693867482015_DP
         casimir_omega(39) = 0.4146412830687338_DP ;  casimir_omega_weight(39) = 0.0309919747446123_DP
         casimir_omega(40) = 0.4467148754805397_DP ;  casimir_omega_weight(40) = 0.0331851777433853_DP
         casimir_omega(41) = 0.4810749061842085_DP ;  casimir_omega_weight(41) = 0.0355682710493892_DP
         casimir_omega(42) = 0.5179219824267806_DP ;  casimir_omega_weight(42) = 0.0381631606297253_DP
         casimir_omega(43) = 0.5574801091692964_DP ;  casimir_omega_weight(43) = 0.0409948141653089_DP
         casimir_omega(44) = 0.6000000000000000_DP ;  casimir_omega_weight(44) = 0.0440917745919031_DP
         casimir_omega(45) = 0.6457629502448748_DP ;  casimir_omega_weight(45) = 0.0474867743345631_DP
         casimir_omega(46) = 0.6950853839282515_DP ;  casimir_omega_weight(46) = 0.0512174730138590_DP
         casimir_omega(47) = 0.7483242118269048_DP ;  casimir_omega_weight(47) = 0.0553273472736239_DP
         casimir_omega(48) = 0.8058831701378680_DP ;  casimir_omega_weight(48) = 0.0598667689600886_DP
         casimir_omega(49) = 0.8682203502161264_DP ;  casimir_omega_weight(49) = 0.0648943177281170_DP
         casimir_omega(50) = 0.9358571820591519_DP ;  casimir_omega_weight(50) = 0.0704783870213295_DP
         casimir_omega(51) = 1.0093892012396821_DP ;  casimir_omega_weight(51) = 0.0766991593121631_DP
         casimir_omega(52) = 1.0894990155669320_DP ;  casimir_omega_weight(52) = 0.0836510489392973_DP
         casimir_omega(53) = 1.1769720003529518_DP ;  casimir_omega_weight(53) = 0.0914457408595557_DP
         casimir_omega(54) = 1.2727153986415052_DP ;  casimir_omega_weight(54) = 0.1002159939792695_DP
         casimir_omega(55) = 1.3777816973889567_DP ;  casimir_omega_weight(55) = 0.1101204324838893_DP
         casimir_omega(56) = 1.4933974094985405_DP ;  casimir_omega_weight(56) = 0.1213496235431059_DP
         casimir_omega(57) = 1.6209987389486102_DP ;  casimir_omega_weight(57) = 0.1341338433453047_DP
         casimir_omega(58) = 1.7622760764135765_DP ;  casimir_omega_weight(58) = 0.1487530779554744_DP
         casimir_omega(59) = 1.9192299152701391_DP ;  casimir_omega_weight(59) = 0.1655500093225807_DP
         casimir_omega(60) = 2.0942416648319853_DP ;  casimir_omega_weight(60) = 0.1849470274399770_DP
         casimir_omega(61) = 2.2901640754070707_DP ;  casimir_omega_weight(61) = 0.2074687291654493_DP
         casimir_omega(62) = 2.5104377372043034_DP ;  casimir_omega_weight(62) = 0.2337719774164176_DP
         casimir_omega(63) = 2.7592426129699259_DP ;  casimir_omega_weight(63) = 0.2646865031867768_DP
         casimir_omega(64) = 3.0416971828355961_DP ;  casimir_omega_weight(64) = 0.3012703995118655_DP
         casimir_omega(65) = 3.3641230982273318_DP ;  casimir_omega_weight(65) = 0.3448869448822227_DP
         casimir_omega(66) = 3.7344011811767963_DP ;  casimir_omega_weight(66) = 0.3973124401926522_DP
         casimir_omega(67) = 4.1624566602802524_DP ;  casimir_omega_weight(67) = 0.4608898856983748_DP
         casimir_omega(68) = 4.6609301790102275_DP ;  casimir_omega_weight(68) = 0.5387516373194012_DP
         casimir_omega(69) = 5.2461205361930219_DP ;  casimir_omega_weight(69) = 0.6351479226654700_DP
         casimir_omega(70) = 5.9393325961573753_DP ;  casimir_omega_weight(70) = 0.7559413743120944_DP
         casimir_omega(71) = 6.7688423135719979_DP ;  casimir_omega_weight(71) = 0.9093682509423305_DP
         casimir_omega(72) = 7.7728242128400398_DP ;  casimir_omega_weight(72) = 1.1072396804954221_DP
         casimir_omega(73) = 9.0038203051464052_DP ;  casimir_omega_weight(73) = 1.3668910403661223_DP
         casimir_omega(74) = 10.5357528712718835_DP ;  casimir_omega_weight(74) = 1.7144472235732333_DP
         casimir_omega(75) = 12.4752812588189315_DP ;  casimir_omega_weight(75) = 2.1904936044172927_DP
         casimir_omega(76) = 14.9808736868070493_DP ;  casimir_omega_weight(76) = 2.8603451226105929_DP
         casimir_omega(77) = 18.2962207621546042_DP ;  casimir_omega_weight(77) = 3.8335708732257161_DP
         casimir_omega(78) = 22.8117819081736499_DP ;  casimir_omega_weight(78) = 5.3033222306428769_DP
         casimir_omega(79) = 29.1851862265366400_DP ;  casimir_omega_weight(79) = 7.6312534948111912_DP
         casimir_omega(80) = 38.5948232590628990_DP ;  casimir_omega_weight(80) = 11.5472519268476734_DP
         casimir_omega(81) = 53.3260123543878066_DP ;  casimir_omega_weight(81) = 18.6727410309475488_DP
         casimir_omega(82) = 78.3004314139001423_DP ;  casimir_omega_weight(82) = 33.0997962113980009_DP
         casimir_omega(83) = 125.7852044374618572_DP ;  casimir_omega_weight(83) = 67.1785499971491475_DP
         casimir_omega(84) = 233.8837370525204165_DP ;  casimir_omega_weight(84) = 169.8275823745695448_DP
         casimir_omega(85) = 575.3818998240872133_DP ;  casimir_omega_weight(85) = 652.8199593470919808_DP
         casimir_omega(86) = 3033.3603122300878567_DP ;  casimir_omega_weight(86) = 7785.6086068807398988_DP
return
endsubroutine gauss_legendre_grid85

subroutine gauss_legendre_grid90()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0001059273116997_DP ;  casimir_omega_weight(2) = 0.0002718756828159_DP
         casimir_omega(3) = 0.0005584052292128_DP ;  casimir_omega_weight(3) = 0.0006335111747478_DP
         casimir_omega(4) = 0.0013735937998130_DP ;  casimir_omega_weight(4) = 0.0009972104150624_DP
         casimir_omega(5) = 0.0025536437899585_DP ;  casimir_omega_weight(5) = 0.0013633686060748_DP
         casimir_omega(6) = 0.0041014364367063_DP ;  casimir_omega_weight(6) = 0.0017328442525776_DP
         casimir_omega(7) = 0.0060207398906271_DP ;  casimir_omega_weight(7) = 0.0021065435335956_DP
         casimir_omega(8) = 0.0083162420439341_DP ;  casimir_omega_weight(8) = 0.0024853993636364_DP
         casimir_omega(9) = 0.0109935774479411_DP ;  casimir_omega_weight(9) = 0.0028703727674292_DP
         casimir_omega(10) = 0.0140593575548921_DP ;  casimir_omega_weight(10) = 0.0032624573622436_DP
         casimir_omega(11) = 0.0175212059349933_DP ;  casimir_omega_weight(11) = 0.0036626847022197_DP
         casimir_omega(12) = 0.0213877990957353_DP ;  casimir_omega_weight(12) = 0.0040721301185178_DP
         casimir_omega(13) = 0.0256689133672669_DP ;  casimir_omega_weight(13) = 0.0044919190135233_DP
         casimir_omega(14) = 0.0303754783225068_DP ;  casimir_omega_weight(14) = 0.0049232336449357_DP
         casimir_omega(15) = 0.0355196372515928_DP ;  casimir_omega_weight(15) = 0.0053673204621017_DP
         casimir_omega(16) = 0.0411148152802416_DP ;  casimir_omega_weight(16) = 0.0058254980712744_DP
         casimir_omega(17) = 0.0471757958041762_DP ;  casimir_omega_weight(17) = 0.0062991659181894_DP
         casimir_omega(18) = 0.0537188060060699_DP ;  casimir_omega_weight(18) = 0.0067898137882976_DP
         casimir_omega(19) = 0.0607616123282846_DP ;  casimir_omega_weight(19) = 0.0072990322382634_DP
         casimir_omega(20) = 0.0683236268957281_DP ;  casimir_omega_weight(20) = 0.0078285240875514_DP
         casimir_omega(21) = 0.0764260260205385_DP ;  casimir_omega_weight(21) = 0.0083801171165148_DP
         casimir_omega(22) = 0.0850918820767024_DP ;  casimir_omega_weight(22) = 0.0089557781378547_DP
         casimir_omega(23) = 0.0943463102111833_DP ;  casimir_omega_weight(23) = 0.0095576286321543_DP
         casimir_omega(24) = 0.1042166315624002_DP ;  casimir_omega_weight(24) = 0.0101879621660263_DP
         casimir_omega(25) = 0.1147325548912866_DP ;  casimir_omega_weight(25) = 0.0108492638439614_DP
         casimir_omega(26) = 0.1259263787998002_DP ;  casimir_omega_weight(26) = 0.0115442320830990_DP
         casimir_omega(27) = 0.1378332170227658_DP ;  casimir_omega_weight(27) = 0.0122758030449196_DP
         casimir_omega(28) = 0.1504912496385587_DP ;  casimir_omega_weight(28) = 0.0130471781105293_DP
         casimir_omega(29) = 0.1639420034609811_DP ;  casimir_omega_weight(29) = 0.0138618548483055_DP
         casimir_omega(30) = 0.1782306653591407_DP ;  casimir_omega_weight(30) = 0.0147236619960677_DP
         casimir_omega(31) = 0.1934064328165061_DP ;  casimir_omega_weight(31) = 0.0156367990668290_DP
         casimir_omega(32) = 0.2095229066995230_DP ;  casimir_omega_weight(32) = 0.0166058812904738_DP
         casimir_omega(33) = 0.2266385319781602_DP ;  casimir_omega_weight(33) = 0.0176359907265733_DP
         casimir_omega(34) = 0.2448170930471710_DP ;  casimir_omega_weight(34) = 0.0187327345303343_DP
         casimir_omega(35) = 0.2641282713640610_DP ;  casimir_omega_weight(35) = 0.0199023115294231_DP
         casimir_omega(36) = 0.2846482743797740_DP ;  casimir_omega_weight(36) = 0.0211515884804308_DP
         casimir_omega(37) = 0.3064605462300053_DP ;  casimir_omega_weight(37) = 0.0224881876279802_DP
         casimir_omega(38) = 0.3296565724267995_DP ;  casimir_omega_weight(38) = 0.0239205874967294_DP
         casimir_omega(39) = 0.3543367929005572_DP ;  casimir_omega_weight(39) = 0.0254582392190124_DP
         casimir_omega(40) = 0.3806116402644715_DP ;  casimir_omega_weight(40) = 0.0271117011541041_DP
         casimir_omega(41) = 0.4086027231968017_DP ;  casimir_omega_weight(41) = 0.0288927951084475_DP
         casimir_omega(42) = 0.4384441784732055_DP ;  casimir_omega_weight(42) = 0.0308147881443405_DP
         casimir_omega(43) = 0.4702842195714704_DP ;  casimir_omega_weight(43) = 0.0328926047988332_DP
         casimir_omega(44) = 0.5042869150896352_DP ;  casimir_omega_weight(44) = 0.0351430755650035_DP
         casimir_omega(45) = 0.5406342366864962_DP ;  casimir_omega_weight(45) = 0.0375852287656830_DP
         casimir_omega(46) = 0.5795284241500895_DP ;  casimir_omega_weight(46) = 0.0402406345415575_DP
         casimir_omega(47) = 0.6211947248799057_DP ;  casimir_omega_weight(47) = 0.0431338116671242_DP
         casimir_omega(48) = 0.6658845769857474_DP ;  casimir_omega_weight(48) = 0.0462927104116474_DP
         casimir_omega(49) = 0.7138793199423218_DP ;  casimir_omega_weight(49) = 0.0497492878247042_DP
         casimir_omega(50) = 0.7654945350452905_DP ;  casimir_omega_weight(50) = 0.0535401958412620_DP
         casimir_omega(51) = 0.8210851407666723_DP ;  casimir_omega_weight(51) = 0.0577076077262511_DP
         casimir_omega(52) = 0.8810513967783992_DP ;  casimir_omega_weight(52) = 0.0623002149568834_DP
         casimir_omega(53) = 0.9458460065747070_DP ;  casimir_omega_weight(53) = 0.0673744351335078_DP
         casimir_omega(54) = 1.0159825544874532_DP ;  casimir_omega_weight(54) = 0.0729958825408905_DP
         casimir_omega(55) = 1.0920455713951776_DP ;  casimir_omega_weight(55) = 0.0792411674024029_DP
         casimir_omega(56) = 1.1747025985191981_DP ;  casimir_omega_weight(56) = 0.0862001088477747_DP
         casimir_omega(57) = 1.2647187157006725_DP ;  casimir_omega_weight(57) = 0.0939784717693682_DP
         casimir_omega(58) = 1.3629741267029842_DP ;  casimir_omega_weight(58) = 0.1027013713302854_DP
         casimir_omega(59) = 1.4704855593176911_DP ;  casimir_omega_weight(59) = 0.1125175340926124_DP
         casimir_omega(60) = 1.5884324561133807_DP ;  casimir_omega_weight(60) = 0.1236046660790392_DP
         casimir_omega(61) = 1.7181892217459398_DP ;  casimir_omega_weight(61) = 0.1361762620628533_DP
         casimir_omega(62) = 1.8613651818993495_DP ;  casimir_omega_weight(62) = 0.1504903064261806_DP
         casimir_omega(63) = 2.0198544356807977_DP ;  casimir_omega_weight(63) = 0.1668604778661157_DP
         casimir_omega(64) = 2.1958985031293854_DP ;  casimir_omega_weight(64) = 0.1856706985969934_DP
         casimir_omega(65) = 2.3921656632171500_DP ;  casimir_omega_weight(65) = 0.2073941943657680_DP
         casimir_omega(66) = 2.6118522644693041_DP ;  casimir_omega_weight(66) = 0.2326187015990291_DP
         casimir_omega(67) = 2.8588132481148683_DP ;  casimir_omega_weight(67) = 0.2620801450262066_DP
         casimir_omega(68) = 3.1377319222178315_DP ;  casimir_omega_weight(68) = 0.2967081272444215_DP
         casimir_omega(69) = 3.4543430794387815_DP ;  casimir_omega_weight(69) = 0.3376881028890870_DP
         casimir_omega(70) = 3.8157295096562915_DP ;  casimir_omega_weight(70) = 0.3865474498410557_DP
         casimir_omega(71) = 4.2307208539058134_DP ;  casimir_omega_weight(71) = 0.4452762873034413_DP
         casimir_omega(72) = 4.7104372521378410_DP ;  casimir_omega_weight(72) = 0.5164996519942509_DP
         casimir_omega(73) = 5.2690411261306940_DP ;  casimir_omega_weight(73) = 0.6037269572525041_DP
         casimir_omega(74) = 5.9247934050034985_DP ;  casimir_omega_weight(74) = 0.7117200549998832_DP
         casimir_omega(75) = 6.7015636937150607_DP ;  casimir_omega_weight(75) = 0.8470472997035710_DP
         casimir_omega(76) = 7.6310318429886728_DP ;  casimir_omega_weight(76) = 1.0189364034366928_DP
         casimir_omega(77) = 8.7559678316980030_DP ;  casimir_omega_weight(77) = 1.2406202817164964_DP
         casimir_omega(78) = 10.1352386413759525_DP ;  casimir_omega_weight(78) = 1.5315210952978482_DP
         casimir_omega(79) = 11.8516652207993651_DP ;  casimir_omega_weight(79) = 1.9209085810616460_DP
         casimir_omega(80) = 14.0247463867743907_DP ;  casimir_omega_weight(80) = 2.4542536745996495_DP
         casimir_omega(81) = 16.8320264459461768_DP ;  casimir_omega_weight(81) = 3.2047337614974092_DP
         casimir_omega(82) = 20.5465309485921424_DP ;  casimir_omega_weight(82) = 4.2951075895289845_DP
         casimir_omega(83) = 25.6057219253760486_DP ;  casimir_omega_weight(83) = 5.9417776157164006_DP
         casimir_omega(84) = 32.7463923099411005_DP ;  casimir_omega_weight(84) = 8.5499331917304655_DP
         casimir_omega(85) = 43.2887833348463644_DP ;  casimir_omega_weight(85) = 12.9373236113910437_DP
         casimir_omega(86) = 59.7933155292830918_DP ;  casimir_omega_weight(86) = 20.9205553584079986_DP
         casimir_omega(87) = 87.7741263470832678_DP ;  casimir_omega_weight(87) = 37.0842978338743237_DP
         casimir_omega(88) = 140.9750261236847564_DP ;  casimir_omega_weight(88) = 75.2653622300040581_DP
         casimir_omega(89) = 262.0862150433237048_DP ;  casimir_omega_weight(89) = 190.2710272287934572_DP
         casimir_omega(90) = 644.6931030846802742_DP ;  casimir_omega_weight(90) = 731.4048359871799221_DP
         casimir_omega(91) = 3398.5569370471571347_DP ;  casimir_omega_weight(91) = 8722.8210838324012002_DP
return
endsubroutine gauss_legendre_grid90

subroutine gauss_legendre_grid95()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000951247666537_DP ;  casimir_omega_weight(2) = 0.0002441466710803_DP
         casimir_omega(3) = 0.0005014329533970_DP ;  casimir_omega_weight(3) = 0.0005688401067335_DP
         casimir_omega(4) = 0.0012333366416934_DP ;  casimir_omega_weight(4) = 0.0008952464225173_DP
         casimir_omega(5) = 0.0022925860730832_DP ;  casimir_omega_weight(5) = 0.0012236386616340_DP
         casimir_omega(6) = 0.0036815045969757_DP ;  casimir_omega_weight(6) = 0.0015547037577225_DP
         casimir_omega(7) = 0.0054031265785296_DP ;  casimir_omega_weight(7) = 0.0018891684565565_DP
         casimir_omega(8) = 0.0074612234481075_DP ;  casimir_omega_weight(8) = 0.0022277796443438_DP
         casimir_omega(9) = 0.0098603235510794_DP ;  casimir_omega_weight(9) = 0.0025713046649797_DP
         casimir_omega(10) = 0.0126057340360093_DP ;  casimir_omega_weight(10) = 0.0029205343664741_DP
         casimir_omega(11) = 0.0157035661956335_DP ;  casimir_omega_weight(11) = 0.0032762868466639_DP
         casimir_omega(12) = 0.0191607647417807_DP ;  casimir_omega_weight(12) = 0.0036394115547772_DP
         casimir_omega(13) = 0.0229851413333405_DP ;  casimir_omega_weight(13) = 0.0040107936951241_DP
         casimir_omega(14) = 0.0271854126636444_DP ;  casimir_omega_weight(14) = 0.0043913589469623_DP
         casimir_omega(15) = 0.0317712434401158_DP ;  casimir_omega_weight(15) = 0.0047820785360881_DP
         casimir_omega(16) = 0.0367532946300997_DP ;  casimir_omega_weight(16) = 0.0051839747038656_DP
         casimir_omega(17) = 0.0421432773962414_DP ;  casimir_omega_weight(17) = 0.0055981266267012_DP
         casimir_omega(18) = 0.0479540132012309_DP ;  casimir_omega_weight(18) = 0.0060256768458846_DP
         casimir_omega(19) = 0.0541995006253205_DP ;  casimir_omega_weight(19) = 0.0064678382751941_DP
         casimir_omega(20) = 0.0608949895115279_DP ;  casimir_omega_weight(20) = 0.0069259018620444_DP
         casimir_omega(21) = 0.0680570631338714_DP ;  casimir_omega_weight(21) = 0.0074012449875059_DP
         casimir_omega(22) = 0.0757037291747295_DP ;  casimir_omega_weight(22) = 0.0078953407015959_DP
         casimir_omega(23) = 0.0838545204000065_DP ;  casimir_omega_weight(23) = 0.0084097679029463_DP
         casimir_omega(24) = 0.0925306060370701_DP ;  casimir_omega_weight(24) = 0.0089462225866660_DP
         casimir_omega(25) = 0.1017549149925716_DP ;  casimir_omega_weight(25) = 0.0095065303012729_DP
         casimir_omega(26) = 0.1115522721977646_DP ;  casimir_omega_weight(26) = 0.0100926599753465_DP
         casimir_omega(27) = 0.1219495495407839_DP ;  casimir_omega_weight(27) = 0.0107067392974995_DP
         casimir_omega(28) = 0.1329758330419836_DP ;  casimir_omega_weight(28) = 0.0113510718600299_DP
         casimir_omega(29) = 0.1446626081539594_DP ;  casimir_omega_weight(29) = 0.0120281563077544_DP
         casimir_omega(30) = 0.1570439653271263_DP ;  casimir_omega_weight(30) = 0.0127407077699504_DP
         casimir_omega(31) = 0.1701568282803996_DP ;  casimir_omega_weight(31) = 0.0134916818959602_DP
         casimir_omega(32) = 0.1840412077613748_DP ;  casimir_omega_weight(32) = 0.0142843018649568_DP
         casimir_omega(33) = 0.1987404839795562_DP ;  casimir_omega_weight(33) = 0.0151220887992346_DP
         casimir_omega(34) = 0.2143017213591434_DP ;  casimir_omega_weight(34) = 0.0160088960795475_DP
         casimir_omega(35) = 0.2307760197962522_DP ;  casimir_omega_weight(35) = 0.0169489481429228_DP
         casimir_omega(36) = 0.2482189072328796_DP ;  casimir_omega_weight(36) = 0.0179468844402097_DP
         casimir_omega(37) = 0.2666907790930859_DP ;  casimir_omega_weight(37) = 0.0190078093456914_DP
         casimir_omega(38) = 0.2862573909856460_DP ;  casimir_omega_weight(38) = 0.0201373489480214_DP
         casimir_omega(39) = 0.3069904120860091_DP ;  casimir_omega_weight(39) = 0.0213417158152629_DP
         casimir_omega(40) = 0.3289680477980898_DP ;  casimir_omega_weight(40) = 0.0226277830225258_DP
         casimir_omega(41) = 0.3522757416988438_DP ;  casimir_omega_weight(41) = 0.0240031689658113_DP
         casimir_omega(42) = 0.3770069684293224_DP ;  casimir_omega_weight(42) = 0.0254763347687636_DP
         casimir_omega(43) = 0.4032641311682586_DP ;  casimir_omega_weight(43) = 0.0270566964312459_DP
         casimir_omega(44) = 0.4311595796738704_DP ;  casimir_omega_weight(44) = 0.0287547542835008_DP
         casimir_omega(45) = 0.4608167676873489_DP ;  casimir_omega_weight(45) = 0.0305822428143471_DP
         casimir_omega(46) = 0.4923715718577722_DP ;  casimir_omega_weight(46) = 0.0325523045581573_DP
         casimir_omega(47) = 0.5259737983975844_DP ;  casimir_omega_weight(47) = 0.0346796924804851_DP
         casimir_omega(48) = 0.5617889085660018_DP ;  casimir_omega_weight(48) = 0.0369810062312859_DP
         casimir_omega(49) = 0.5999999999999999_DP ;  casimir_omega_weight(49) = 0.0394749687820568_DP
         casimir_omega(50) = 0.6408100881146273_DP ;  casimir_omega_weight(50) = 0.0421827513863302_DP
         casimir_omega(51) = 0.6844447405873927_DP ;  casimir_omega_weight(51) = 0.0451283565754996_DP
         casimir_omega(52) = 0.7311551287205316_DP ;  casimir_omega_weight(52) = 0.0483390711197353_DP
         casimir_omega(53) = 0.7812215727450483_DP ;  casimir_omega_weight(53) = 0.0518460036716899_DP
         casimir_omega(54) = 0.8349576745396788_DP ;  casimir_omega_weight(54) = 0.0556847253322589_DP
         casimir_omega(55) = 0.8927151516230259_DP ;  casimir_omega_weight(55) = 0.0598960358489208_DP
         casimir_omega(56) = 0.9548895117239439_DP ;  casimir_omega_weight(56) = 0.0645268838642711_DP
         casimir_omega(57) = 1.0219267391615039_DP ;  casimir_omega_weight(57) = 0.0696314769574587_DP
         casimir_omega(58) = 1.0943312045337497_DP ;  casimir_omega_weight(58) = 0.0752726266782212_DP
         casimir_omega(59) = 1.1726750602853984_DP ;  casimir_omega_weight(59) = 0.0815233860569088_DP
         casimir_omega(60) = 1.2576094498746127_DP ;  casimir_omega_weight(60) = 0.0884690531317123_DP
         casimir_omega(61) = 1.3498779418779441_DP ;  casimir_omega_weight(61) = 0.0962096351678329_DP
         casimir_omega(62) = 1.4503327083873077_DP ;  casimir_omega_weight(62) = 0.1048628962533586_DP
         casimir_omega(63) = 1.5599541075274510_DP ;  casimir_omega_weight(63) = 0.1145681483594582_DP
         casimir_omega(64) = 1.6798745139180886_DP ;  casimir_omega_weight(64) = 0.1254909962898805_DP
         casimir_omega(65) = 1.8114074837265313_DP ;  casimir_omega_weight(65) = 0.1378293152558061_DP
         casimir_omega(66) = 1.9560836639736197_DP ;  casimir_omega_weight(66) = 0.1518208333295631_DP
         casimir_omega(67) = 2.1156952890938969_DP ;  casimir_omega_weight(67) = 0.1677528202523745_DP
         casimir_omega(68) = 2.2923516943176483_DP ;  casimir_omega_weight(68) = 0.1859745643993192_DP
         casimir_omega(69) = 2.4885490770141883_DP ;  casimir_omega_weight(69) = 0.2069135740037679_DP
         casimir_omega(70) = 2.7072588436903420_DP ;  casimir_omega_weight(70) = 0.2310968013919238_DP
         casimir_omega(71) = 2.9520404245495371_DP ;  casimir_omega_weight(71) = 0.2591787123474543_DP
         casimir_omega(72) = 3.2271866176045005_DP ;  casimir_omega_weight(72) = 0.2919787877626256_DP
         casimir_omega(73) = 3.5379126406452319_DP ;  casimir_omega_weight(73) = 0.3305321784604342_DP
         casimir_omega(74) = 3.8906045839122161_DP ;  casimir_omega_weight(74) = 0.3761589391345393_DP
         casimir_omega(75) = 4.2931495914914697_DP ;  casimir_omega_weight(75) = 0.4305598727993014_DP
         casimir_omega(76) = 4.7553800047167512_DP ;  casimir_omega_weight(76) = 0.4959510675641396_DP
         casimir_omega(77) = 5.2896787404984504_DP ;  casimir_omega_weight(77) = 0.5752556231616820_DP
         casimir_omega(78) = 5.9118164382284366_DP ;  casimir_omega_weight(78) = 0.6723814357475129_DP
         casimir_omega(79) = 6.6421276182721565_DP ;  casimir_omega_weight(79) = 0.7926310527317999_DP
         casimir_omega(80) = 7.5071923279772301_DP ;  casimir_omega_weight(80) = 0.9433186498587514_DP
         casimir_omega(81) = 8.5422876966874384_DP ;  casimir_omega_weight(81) = 1.1347197266636273_DP
         casimir_omega(82) = 9.7950402439614894_DP ;  casimir_omega_weight(82) = 1.3815697710664976_DP
         casimir_omega(83) = 11.3310012772571227_DP ;  casimir_omega_weight(83) = 1.7054962958086830_DP
         casimir_omega(84) = 13.2423960031121357_DP ;  casimir_omega_weight(84) = 2.1390925672888050_DP
         casimir_omega(85) = 15.6622922077843825_DP ;  casimir_omega_weight(85) = 2.7329926724034670_DP
         casimir_omega(86) = 18.7883941403970454_DP ;  casimir_omega_weight(86) = 3.5686831737544962_DP
         casimir_omega(87) = 22.9247290402165582_DP ;  casimir_omega_weight(87) = 4.7828618851351603_DP
         casimir_omega(88) = 28.5584321366473262_DP ;  casimir_omega_weight(88) = 6.6165034316480140_DP
         casimir_omega(89) = 36.5099581301860994_DP ;  casimir_omega_weight(89) = 9.5208058003414529_DP
         casimir_omega(90) = 48.2494596903286990_DP ;  casimir_omega_weight(90) = 14.4063724798322834_DP
         casimir_omega(91) = 66.6280892678908145_DP ;  casimir_omega_weight(91) = 23.2960828764768344_DP
         casimir_omega(92) = 97.7861063369944787_DP ;  casimir_omega_weight(92) = 41.2951886845587595_DP
         casimir_omega(93) = 157.0279101956954548_DP ;  casimir_omega_weight(93) = 83.8116501391387771_DP
         casimir_omega(94) = 291.8911089073782819_DP ;  casimir_omega_weight(94) = 211.8760297717404910_DP
         casimir_omega(95) = 717.9424438723993944_DP ;  casimir_omega_weight(95) = 814.4547613674344575_DP
         casimir_omega(96) = 3784.5033702979881127_DP ;  casimir_omega_weight(96) = 9713.2842692106696632_DP
return
endsubroutine gauss_legendre_grid95

subroutine gauss_legendre_grid100()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000858942625814_DP ;  casimir_omega_weight(2) = 0.0002204534525562_DP
         casimir_omega(3) = 0.0004527561974681_DP ;  casimir_omega_weight(3) = 0.0005135920261776_DP
         casimir_omega(4) = 0.0011135223051184_DP ;  casimir_omega_weight(4) = 0.0008081689954694_DP
         casimir_omega(5) = 0.0020696333184536_DP ;  casimir_omega_weight(5) = 0.0011043679094102_DP
         casimir_omega(6) = 0.0033229837649868_DP ;  casimir_omega_weight(6) = 0.0014027445450094_DP
         casimir_omega(7) = 0.0048760447043133_DP ;  casimir_omega_weight(7) = 0.0017038885504445_DP
         casimir_omega(8) = 0.0067318848759188_DP ;  casimir_omega_weight(8) = 0.0020084050855163_DP
         casimir_omega(9) = 0.0088941856376975_DP ;  casimir_omega_weight(9) = 0.0023169144810549_DP
         casimir_omega(10) = 0.0113672570995894_DP ;  casimir_omega_weight(10) = 0.0026300543245760_DP
         casimir_omega(11) = 0.0141560566874287_DP ;  casimir_omega_weight(11) = 0.0029484821311963_DP
         casimir_omega(12) = 0.0172662105209422_DP ;  casimir_omega_weight(12) = 0.0032728782808808_DP
         casimir_omega(13) = 0.0207040378343009_DP ;  casimir_omega_weight(13) = 0.0036039491638462_DP
         casimir_omega(14) = 0.0244765786459272_DP ;  casimir_omega_weight(14) = 0.0039424305360497_DP
         casimir_omega(15) = 0.0285916248968266_DP ;  casimir_omega_weight(15) = 0.0042890911048851_DP
         casimir_omega(16) = 0.0330577553010846_DP ;  casimir_omega_weight(16) = 0.0046447363729020_DP
         casimir_omega(17) = 0.0378843741825312_DP ;  casimir_omega_weight(17) = 0.0050102127722131_DP
         casimir_omega(18) = 0.0430817546063617_DP ;  casimir_omega_weight(18) = 0.0053864121265173_DP
         casimir_omega(19) = 0.0486610861535960_DP ;  casimir_omega_weight(19) = 0.0057742764820112_DP
         casimir_omega(20) = 0.0546345277298805_DP ;  casimir_omega_weight(20) = 0.0061748033532789_DP
         casimir_omega(21) = 0.0610152658488987_DP ;  casimir_omega_weight(21) = 0.0065890514356949_DP
         casimir_omega(22) = 0.0678175788851898_DP ;  casimir_omega_weight(22) = 0.0070181468420149_DP
         casimir_omega(23) = 0.0750569078523359_DP ;  casimir_omega_weight(23) = 0.0074632899279584_DP
         casimir_omega(24) = 0.0827499343312168_DP ;  casimir_omega_weight(24) = 0.0079257627796486_DP
         casimir_omega(25) = 0.0909146662504676_DP ;  casimir_omega_weight(25) = 0.0084069374451244_DP
         casimir_omega(26) = 0.0995705323087081_DP ;  casimir_omega_weight(26) = 0.0089082850028277_DP
         casimir_omega(27) = 0.1087384859270335_DP ;  casimir_omega_weight(27) = 0.0094313855723050_DP
         casimir_omega(28) = 0.1184411197324416_DP ;  casimir_omega_weight(28) = 0.0099779393865793_DP
         casimir_omega(29) = 0.1287027917003242_DP ;  casimir_omega_weight(29) = 0.0105497790620522_DP
         casimir_omega(30) = 0.1395497642292843_DP ;  casimir_omega_weight(30) = 0.0111488832207926_DP
         casimir_omega(31) = 0.1510103575871145_DP ;  casimir_omega_weight(31) = 0.0117773916420540_DP
         casimir_omega(32) = 0.1631151193560597_DP ;  casimir_omega_weight(32) = 0.0124376221453950_DP
         casimir_omega(33) = 0.1758970117223359_DP ;  casimir_omega_weight(33) = 0.0131320894374632_DP
         casimir_omega(34) = 0.1893916187037813_DP ;  casimir_omega_weight(34) = 0.0138635261890972_DP
         casimir_omega(35) = 0.2036373756958502_DP ;  casimir_omega_weight(35) = 0.0146349066497887_DP
         casimir_omega(36) = 0.2186758240461980_DP ;  casimir_omega_weight(36) = 0.0154494731538041_DP
         casimir_omega(37) = 0.2345518937493921_DP ;  casimir_omega_weight(37) = 0.0163107659277027_DP
         casimir_omega(38) = 0.2513142177946808_DP ;  casimir_omega_weight(38) = 0.0172226566740985_DP
         casimir_omega(39) = 0.2690154822119637_DP ;  casimir_omega_weight(39) = 0.0181893864833153_DP
         casimir_omega(40) = 0.2877128164567630_DP ;  casimir_omega_weight(40) = 0.0192156087151732_DP
         casimir_omega(41) = 0.3074682294694304_DP ;  casimir_omega_weight(41) = 0.0203064376005550_DP
         casimir_omega(42) = 0.3283490975553116_DP ;  casimir_omega_weight(42) = 0.0214675034397648_DP
         casimir_omega(43) = 0.3504287111832737_DP ;  casimir_omega_weight(43) = 0.0227050154264389_DP
         casimir_omega(44) = 0.3737868889165822_DP ;  casimir_omega_weight(44) = 0.0240258333068610_DP
         casimir_omega(45) = 0.3985106680051070_DP ;  casimir_omega_weight(45) = 0.0254375493013016_DP
         casimir_omega(46) = 0.4246950827206421_DP ;  casimir_omega_weight(46) = 0.0269485819744063_DP
         casimir_omega(47) = 0.4524440433560248_DP ;  casimir_omega_weight(47) = 0.0285682840552213_DP
         casimir_omega(48) = 0.4818713309927963_DP ;  casimir_omega_weight(48) = 0.0303070665864009_DP
         casimir_omega(49) = 0.5131017257439772_DP ;  casimir_omega_weight(49) = 0.0321765422416135_DP
         casimir_omega(50) = 0.5462722892878303_DP ;  casimir_omega_weight(50) = 0.0341896912091083_DP
         casimir_omega(51) = 0.5815338262362277_DP ;  casimir_omega_weight(51) = 0.0363610537219063_DP
         casimir_omega(52) = 0.6190525533655930_DP ;  casimir_omega_weight(52) = 0.0387069541513929_DP
         casimir_omega(53) = 0.6590120111516703_DP ;  casimir_omega_weight(53) = 0.0412457626099662_DP
         casimir_omega(54) = 0.7016152586078608_DP ;  casimir_omega_weight(54) = 0.0439982012791380_DP
         casimir_omega(55) = 0.7470874004027065_DP ;  casimir_omega_weight(55) = 0.0469877042554454_DP
         casimir_omega(56) = 0.7956785049697709_DP ;  casimir_omega_weight(56) = 0.0502408416696153_DP
         casimir_omega(57) = 0.8476669842603353_DP ;  casimir_omega_weight(57) = 0.0537878212905126_DP
         casimir_omega(58) = 0.9033635204851952_DP ;  casimir_omega_weight(58) = 0.0576630839128375_DP
         casimir_omega(59) = 0.9631156433642096_DP ;  casimir_omega_weight(59) = 0.0619060127276507_DP
         casimir_omega(60) = 1.0273130839776436_DP ;  casimir_omega_weight(60) = 0.0665617818264228_DP
         casimir_omega(61) = 1.0963940594944279_DP ;  casimir_omega_weight(61) = 0.0716823753096186_DP
         casimir_omega(62) = 1.1708526784091444_DP ;  casimir_omega_weight(62) = 0.0773278165831514_DP
         casimir_omega(63) = 1.2512477005142391_DP ;  casimir_omega_weight(63) = 0.0835676578990864_DP
         casimir_omega(64) = 1.3382129423924660_DP ;  casimir_omega_weight(64) = 0.0904827938005837_DP
         casimir_omega(65) = 1.4324696913650692_DP ;  casimir_omega_weight(65) = 0.0981676799145058_DP
         casimir_omega(66) = 1.5348415834350233_DP ;  casimir_omega_weight(66) = 0.1067330619392123_DP
         casimir_omega(67) = 1.6462725203859063_DP ;  casimir_omega_weight(67) = 0.1163093506951838_DP
         casimir_omega(68) = 1.7678483567657561_DP ;  casimir_omega_weight(68) = 0.1270508205276188_DP
         casimir_omega(69) = 1.9008232912516521_DP ;  casimir_omega_weight(69) = 0.1391408640966840_DP
         casimir_omega(70) = 2.0466521658041676_DP ;  casimir_omega_weight(70) = 0.1527986122421741_DP
         casimir_omega(71) = 2.2070302337465417_DP ;  casimir_omega_weight(71) = 0.1682873311754863_DP
         casimir_omega(72) = 2.3839424378047989_DP ;  casimir_omega_weight(72) = 0.1859251523587903_DP
         casimir_omega(73) = 2.5797248887393960_DP ;  casimir_omega_weight(73) = 0.2060988901355118_DP
         casimir_omega(74) = 2.7971421229015410_DP ;  casimir_omega_weight(74) = 0.2292819838009381_DP
         casimir_omega(75) = 3.0394849425034103_DP ;  casimir_omega_weight(75) = 0.2560580024169791_DP
         casimir_omega(76) = 3.3106953525320382_DP ;  casimir_omega_weight(76) = 0.2871517302817795_DP
         casimir_omega(77) = 3.6155275225792525_DP ;  casimir_omega_weight(77) = 0.3234706982066235_DP
         casimir_omega(78) = 3.9597571530231424_DP ;  casimir_omega_weight(78) = 0.3661612813012552_DP
         casimir_omega(79) = 4.3504566246428311_DP ;  casimir_omega_weight(79) = 0.4166853722452051_DP
         casimir_omega(80) = 4.7963606588783296_DP ;  casimir_omega_weight(80) = 0.4769265244271483_DP
         casimir_omega(81) = 5.3083581855591566_DP ;  casimir_omega_weight(81) = 0.5493389449855610_DP
         casimir_omega(82) = 5.9001627706011126_DP ;  casimir_omega_weight(82) = 0.6371598227686019_DP
         casimir_omega(83) = 6.5892397163179126_DP ;  casimir_omega_weight(83) = 0.7447169617177284_DP
         casimir_omega(84) = 7.3981086008577837_DP ;  casimir_omega_weight(84) = 0.8778826755008849_DP
         casimir_omega(85) = 8.3562056209020135_DP ;  casimir_omega_weight(85) = 1.0447570601373133_DP
         casimir_omega(86) = 9.5025985717879529_DP ;  casimir_omega_weight(86) = 1.2567197363270397_DP
         casimir_omega(87) = 10.8900316044201855_DP ;  casimir_omega_weight(87) = 1.5300895488642481_DP
         casimir_omega(88) = 12.5910997118585168_DP ;  casimir_omega_weight(88) = 1.8888179307658384_DP
         casimir_omega(89) = 14.7079379519367883_DP ;  casimir_omega_weight(89) = 2.3690003632794614_DP
         casimir_omega(90) = 17.3879125840648925_DP ;  casimir_omega_weight(90) = 3.0267116743070770_DP
         casimir_omega(91) = 20.8499716578428931_DP ;  casimir_omega_weight(91) = 3.9521943343135248_DP
         casimir_omega(92) = 25.4308108500082888_DP ;  casimir_omega_weight(92) = 5.2968346361362428_DP
         casimir_omega(93) = 31.6699091826646395_DP ;  casimir_omega_weight(93) = 7.3275004581189007_DP
         casimir_omega(94) = 40.4758810603369454_DP ;  casimir_omega_weight(94) = 10.5438720060731335_DP
         casimir_omega(95) = 53.4768503376790179_DP ;  casimir_omega_weight(95) = 15.9543991252416877_DP
         casimir_omega(96) = 73.8303321299630255_DP ;  casimir_omega_weight(96) = 25.7993240875103034_DP
         casimir_omega(97) = 108.3363704009597086_DP ;  casimir_omega_weight(97) = 45.7324691764427911_DP
         casimir_omega(98) = 173.9438560396436344_DP ;  casimir_omega_weight(98) = 92.8174140493234745_DP
         casimir_omega(99) = 323.2984183120259445_DP ;  casimir_omega_weight(99) = 234.6425902409612831_DP
         casimir_omega(100) = 795.1299220489955815_DP ;  casimir_omega_weight(100) = 901.9697356398478405_DP
         casimir_omega(101) = 4191.1996119465720767_DP ;  casimir_omega_weight(101) = 10756.9981630536349257_DP
return
endsubroutine gauss_legendre_grid100

subroutine gauss_legendre_grid110()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000710501795945_DP ;  casimir_omega_weight(2) = 0.0001823520992504_DP
         casimir_omega(3) = 0.0003744853595292_DP ;  casimir_omega_weight(3) = 0.0004247671971259_DP
         casimir_omega(4) = 0.0009209038667499_DP ;  casimir_omega_weight(4) = 0.0006682281720777_DP
         casimir_omega(5) = 0.0017113118514691_DP ;  casimir_omega_weight(5) = 0.0009128033001933_DP
         casimir_omega(6) = 0.0027470064747296_DP ;  casimir_omega_weight(6) = 0.0011588670321061_DP
         casimir_omega(7) = 0.0040296753782299_DP ;  casimir_omega_weight(7) = 0.0014068194485151_DP
         casimir_omega(8) = 0.0055614114523420_DP ;  casimir_omega_weight(8) = 0.0016570703941240_DP
         casimir_omega(9) = 0.0073447217741400_DP ;  casimir_omega_weight(9) = 0.0019100384864284_DP
         casimir_omega(10) = 0.0093825368063088_DP ;  casimir_omega_weight(10) = 0.0021661520937332_DP
         casimir_omega(11) = 0.0116782208369624_DP ;  casimir_omega_weight(11) = 0.0024258507517030_DP
         casimir_omega(12) = 0.0142355839282276_DP ;  casimir_omega_weight(12) = 0.0026895867462202_DP
         casimir_omega(13) = 0.0170588955039835_DP ;  casimir_omega_weight(13) = 0.0029578268052071_DP
         casimir_omega(14) = 0.0201528996795751_DP ;  casimir_omega_weight(14) = 0.0032310538909560_DP
         casimir_omega(15) = 0.0235228324361253_DP ;  casimir_omega_weight(15) = 0.0035097690984706_DP
         casimir_omega(16) = 0.0271744407505080_DP ;  casimir_omega_weight(16) = 0.0037944936704746_DP
         casimir_omega(17) = 0.0311140038041038_DP ;  casimir_omega_weight(17) = 0.0040857711423456_DP
         casimir_omega(18) = 0.0353483564076747_DP ;  casimir_omega_weight(18) = 0.0043841696320750_DP
         casimir_omega(19) = 0.0398849147956439_DP ;  casimir_omega_weight(19) = 0.0046902842920673_DP
         casimir_omega(20) = 0.0447317049607634_DP ;  casimir_omega_weight(20) = 0.0050047399413625_DP
         casimir_omega(21) = 0.0498973937196595_DP ;  casimir_omega_weight(21) = 0.0053281938987608_DP
         casimir_omega(22) = 0.0553913227213082_DP ;  casimir_omega_weight(22) = 0.0056613390395268_DP
         casimir_omega(23) = 0.0612235456343372_DP ;  casimir_omega_weight(23) = 0.0060049071007336_DP
         casimir_omega(24) = 0.0674048687754590_DP ;  casimir_omega_weight(24) = 0.0063596722630626_DP
         casimir_omega(25) = 0.0739468954706692_DP ;  casimir_omega_weight(25) = 0.0067264550399630_DP
         casimir_omega(26) = 0.0808620744734665_DP ;  casimir_omega_weight(26) = 0.0071061265085911_DP
         casimir_omega(27) = 0.0881637528007146_DP ;  casimir_omega_weight(27) = 0.0074996129209094_DP
         casimir_omega(28) = 0.0958662333873659_DP ;  casimir_omega_weight(28) = 0.0079079007378668_DP
         casimir_omega(29) = 0.1039848380066826_DP ;  casimir_omega_weight(29) = 0.0083320421346937_DP
         casimir_omega(30) = 0.1125359759535036_DP ;  casimir_omega_weight(30) = 0.0087731610311973_DP
         casimir_omega(31) = 0.1215372190452311_DP ;  casimir_omega_weight(31) = 0.0092324597076005_DP
         casimir_omega(32) = 0.1310073835594657_DP ;  casimir_omega_weight(32) = 0.0097112260740481_DP
         casimir_omega(33) = 0.1409666197996002_DP ;  casimir_omega_weight(33) = 0.0102108416705941_DP
         casimir_omega(34) = 0.1514365100613184_DP ;  casimir_omega_weight(34) = 0.0107327904843890_DP
         casimir_omega(35) = 0.1624401758652455_DP ;  casimir_omega_weight(35) = 0.0112786686821714_DP
         casimir_omega(36) = 0.1740023954254540_DP ;  casimir_omega_weight(36) = 0.0118501953692014_DP
         casimir_omega(37) = 0.1861497324420062_DP ;  casimir_omega_weight(37) = 0.0124492245007818_DP
         casimir_omega(38) = 0.1989106774402810_DP ;  casimir_omega_weight(38) = 0.0130777580897870_DP
         casimir_omega(39) = 0.2123158030329738_DP ;  casimir_omega_weight(39) = 0.0137379608735165_DP
         casimir_omega(40) = 0.2263979346552223_DP ;  casimir_omega_weight(40) = 0.0144321766262559_DP
         casimir_omega(41) = 0.2411923385226699_DP ;  casimir_omega_weight(41) = 0.0151629463305251_DP
         casimir_omega(42) = 0.2567369287903668_DP ;  casimir_omega_weight(42) = 0.0159330284509945_DP
         casimir_omega(43) = 0.2730724961518676_DP ;  casimir_omega_weight(43) = 0.0167454215908986_DP
         casimir_omega(44) = 0.2902429604181541_DP ;  casimir_omega_weight(44) = 0.0176033898527403_DP
         casimir_omega(45) = 0.3082956499615759_DP ;  casimir_omega_weight(45) = 0.0185104912738857_DP
         casimir_omega(46) = 0.3272816113084583_DP ;  casimir_omega_weight(46) = 0.0194706097648479_DP
         casimir_omega(47) = 0.3472559526244382_DP ;  casimir_omega_weight(47) = 0.0204879910451249_DP
         casimir_omega(48) = 0.3682782253697777_DP ;  casimir_omega_weight(48) = 0.0215672831501796_DP
         casimir_omega(49) = 0.3904128490206458_DP ;  casimir_omega_weight(49) = 0.0227135821760436_DP
         casimir_omega(50) = 0.4137295844720898_DP ;  casimir_omega_weight(50) = 0.0239324840375854_DP
         casimir_omega(51) = 0.4383040625775484_DP ;  casimir_omega_weight(51) = 0.0252301431463203_DP
         casimir_omega(52) = 0.4642183752604204_DP ;  casimir_omega_weight(52) = 0.0266133390676791_DP
         casimir_omega(53) = 0.4915617377821764_DP ;  casimir_omega_weight(53) = 0.0280895524010986_DP
         casimir_omega(54) = 0.5204312321010726_DP ;  casimir_omega_weight(54) = 0.0296670513452401_DP
         casimir_omega(55) = 0.5509326428449136_DP ;  casimir_omega_weight(55) = 0.0313549906727426_DP
         casimir_omega(56) = 0.5831813992982907_DP ;  casimir_omega_weight(56) = 0.0331635251536972_DP
         casimir_omega(57) = 0.6173036390275262_DP ;  casimir_omega_weight(57) = 0.0351039398461468_DP
         casimir_omega(58) = 0.6534374114066420_DP ;  casimir_omega_weight(58) = 0.0371888001300439_DP
         casimir_omega(59) = 0.6917340424528650_DP ;  casimir_omega_weight(59) = 0.0394321249165818_DP
         casimir_omega(60) = 0.7323596861388052_DP ;  casimir_omega_weight(60) = 0.0418495871404951_DP
         casimir_omega(61) = 0.7754970918547659_DP ;  casimir_omega_weight(61) = 0.0444587464680876_DP
         casimir_omega(62) = 0.8213476231156446_DP ;  casimir_omega_weight(62) = 0.0472793201647123_DP
         casimir_omega(63) = 0.8701335691508574_DP ;  casimir_omega_weight(63) = 0.0503334993093178_DP
         casimir_omega(64) = 0.9221007989441526_DP ;  casimir_omega_weight(64) = 0.0536463190797935_DP
         casimir_omega(65) = 0.9775218169321143_DP ;  casimir_omega_weight(65) = 0.0572460937381929_DP
         casimir_omega(66) = 1.0366992913418662_DP ;  casimir_omega_weight(66) = 0.0611649293179137_DP
         casimir_omega(67) = 1.0999701405793461_DP ;  casimir_omega_weight(67) = 0.0654393299842937_DP
         casimir_omega(68) = 1.1677102808452469_DP ;  casimir_omega_weight(68) = 0.0701109177722925_DP
         casimir_omega(69) = 1.2403401601242878_DP ;  casimir_omega_weight(69) = 0.0752272901200514_DP
         casimir_omega(70) = 1.3183312309848618_DP ;  casimir_omega_weight(70) = 0.0808430456028518_DP
         casimir_omega(71) = 1.4022135486942389_DP ;  casimir_omega_weight(71) = 0.0870210159129761_DP
         casimir_omega(72) = 1.4925847238972860_DP ;  casimir_omega_weight(72) = 0.0938337519377279_DP
         casimir_omega(73) = 1.5901205130172136_DP ;  casimir_omega_weight(73) = 0.1013653244489362_DP
         casimir_omega(74) = 1.6955873979106921_DP ;  casimir_omega_weight(74) = 0.1097135163627389_DP
         casimir_omega(75) = 1.8098575935325636_DP ;  casimir_omega_weight(75) = 0.1189925050267271_DP
         casimir_omega(76) = 1.9339270343144612_DP ;  casimir_omega_weight(76) = 0.1293361612851810_DP
         casimir_omega(77) = 2.0689370345722091_DP ;  casimir_omega_weight(77) = 0.1409021295730334_DP
         casimir_omega(78) = 2.2162005063245136_DP ;  casimir_omega_weight(78) = 0.1538769033642948_DP
         casimir_omega(79) = 2.3772338642394168_DP ;  casimir_omega_weight(79) = 0.1684821776924604_DP
         casimir_omega(80) = 2.5537960725154663_DP ;  casimir_omega_weight(80) = 0.1849828519156441_DP
         casimir_omega(81) = 2.7479367209604040_DP ;  casimir_omega_weight(81) = 0.2036971811005720_DP
         casimir_omega(82) = 2.9620555976850445_DP ;  casimir_omega_weight(82) = 0.2250097474019084_DP
         casimir_omega(83) = 3.1989770111270022_DP ;  casimir_omega_weight(83) = 0.2493881642374702_DP
         casimir_omega(84) = 3.4620431872660604_DP ;  casimir_omega_weight(84) = 0.2774047665158252_DP
         casimir_omega(85) = 3.7552325493518763_DP ;  casimir_omega_weight(85) = 0.3097650256883333_DP
         casimir_omega(86) = 4.0833107548602641_DP ;  casimir_omega_weight(86) = 0.3473451290856228_DP
         casimir_omega(87) = 4.4520252831026266_DP ;  casimir_omega_weight(87) = 0.3912421872326978_DP
         casimir_omega(88) = 4.8683585390382271_DP ;  casimir_omega_weight(88) = 0.4428420506747199_DP
         casimir_omega(89) = 5.3408604829309967_DP ;  casimir_omega_weight(89) = 0.5039120005905011_DP
         casimir_omega(90) = 5.8800906786766385_DP ;  casimir_omega_weight(90) = 0.5767290656479162_DP
         casimir_omega(91) = 6.4992129148328432_DP ;  casimir_omega_weight(91) = 0.6642601402761894_DP
         casimir_omega(92) = 7.2148056874994717_DP ;  casimir_omega_weight(92) = 0.7704186687757737_DP
         casimir_omega(93) = 8.0479829757388881_DP ;  casimir_omega_weight(93) = 0.9004365445362105_DP
         casimir_omega(94) = 9.0259688868464050_DP ;  casimir_omega_weight(94) = 1.0614128250637422_DP
         casimir_omega(95) = 10.1843490500123810_DP ;  casimir_omega_weight(95) = 1.2631397429787921_DP
         casimir_omega(96) = 11.5703527667666588_DP ;  casimir_omega_weight(96) = 1.5193741615144229_DP
         casimir_omega(97) = 13.2477427338876925_DP ;  casimir_omega_weight(97) = 1.8498439917619303_DP
         casimir_omega(98) = 15.3042794050231343_DP ;  casimir_omega_weight(98) = 2.2835042113219499_DP
         casimir_omega(99) = 17.8634343307361014_DP ;  casimir_omega_weight(99) = 2.8639907863313558_DP
         casimir_omega(100) = 21.1033592365893661_DP ;  casimir_omega_weight(100) = 3.6590927950362868_DP
         casimir_omega(101) = 25.2887413551166631_DP ;  casimir_omega_weight(101) = 4.7779047154108429_DP
         casimir_omega(102) = 30.8266134906932052_DP ;  casimir_omega_weight(102) = 6.4034380367399146_DP
         casimir_omega(103) = 38.3691540392292438_DP ;  casimir_omega_weight(103) = 8.8583103986297793_DP
         casimir_omega(104) = 49.0147906306706886_DP ;  casimir_omega_weight(104) = 12.7465871938743991_DP
         casimir_omega(105) = 64.7317687398236927_DP ;  casimir_omega_weight(105) = 19.2873874657983642_DP
         casimir_omega(106) = 89.3372210438784862_DP ;  casimir_omega_weight(106) = 31.1889490453237492_DP
         casimir_omega(107) = 131.0517478978387942_DP ;  casimir_omega_weight(107) = 55.2862002822853000_DP
         casimir_omega(108) = 210.3649312607361708_DP ;  casimir_omega_weight(108) = 112.2073708160632464_DP
         casimir_omega(109) = 390.9202827766663404_DP ;  casimir_omega_weight(109) = 283.6603856489143709_DP
         casimir_omega(110) = 961.3192901654275602_DP ;  casimir_omega_weight(110) = 1090.3948312955790243_DP
         casimir_omega(111) = 5066.8415203746617408_DP ;  casimir_omega_weight(111) = 13004.1780764116156206_DP
return
endsubroutine gauss_legendre_grid110

subroutine gauss_legendre_grid120()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000597462510774_DP ;  casimir_omega_weight(2) = 0.0001533383525200_DP
         casimir_omega(3) = 0.0003148886692287_DP ;  casimir_omega_weight(3) = 0.0003571448730181_DP
         casimir_omega(4) = 0.0007742737948095_DP ;  casimir_omega_weight(4) = 0.0005617387663876_DP
         casimir_omega(5) = 0.0014386289916722_DP ;  casimir_omega_weight(5) = 0.0007671241176854_DP
         casimir_omega(6) = 0.0023088725773103_DP ;  casimir_omega_weight(6) = 0.0009735615566897_DP
         casimir_omega(7) = 0.0033861963735281_DP ;  casimir_omega_weight(7) = 0.0011813319965853_DP
         casimir_omega(8) = 0.0046720766781627_DP ;  casimir_omega_weight(8) = 0.0013907229347456_DP
         casimir_omega(9) = 0.0061682799706518_DP ;  casimir_omega_weight(9) = 0.0016020272486196_DP
         casimir_omega(10) = 0.0078768684583749_DP ;  casimir_omega_weight(10) = 0.0018155436310052_DP
         casimir_omega(11) = 0.0098002062696703_DP ;  casimir_omega_weight(11) = 0.0020315773741769_DP
         casimir_omega(12) = 0.0119409664970850_DP ;  casimir_omega_weight(12) = 0.0022504412703723_DP
         casimir_omega(13) = 0.0143021391751351_DP ;  casimir_omega_weight(13) = 0.0024724565766071_DP
         casimir_omega(14) = 0.0168870402496755_DP ;  casimir_omega_weight(14) = 0.0026979540325322_DP
         casimir_omega(15) = 0.0196993215914816_DP ;  casimir_omega_weight(15) = 0.0029272749314494_DP
         casimir_omega(16) = 0.0227429821090122_DP ;  casimir_omega_weight(16) = 0.0031607722485016_DP
         casimir_omega(17) = 0.0260223800202704_DP ;  casimir_omega_weight(17) = 0.0033988118317227_DP
         casimir_omega(18) = 0.0295422463499067_DP ;  casimir_omega_weight(18) = 0.0036417736626582_DP
         casimir_omega(19) = 0.0333076997248084_DP ;  casimir_omega_weight(19) = 0.0038900531940200_DP
         casimir_omega(20) = 0.0373242625492504_DP ;  casimir_omega_weight(20) = 0.0041440627725891_DP
         casimir_omega(21) = 0.0415978786492945_DP ;  casimir_omega_weight(21) = 0.0044042331563759_DP
         casimir_omega(22) = 0.0461349324855225_DP ;  casimir_omega_weight(22) = 0.0046710151358382_DP
         casimir_omega(23) = 0.0509422700434658_DP ;  casimir_omega_weight(23) = 0.0049448812699368_DP
         casimir_omega(24) = 0.0560272215223752_DP ;  casimir_omega_weight(24) = 0.0052263277488119_DP
         casimir_omega(25) = 0.0613976259553274_DP ;  casimir_omega_weight(25) = 0.0055158763960465_DP
         casimir_omega(26) = 0.0670618579072558_DP ;  casimir_omega_weight(26) = 0.0058140768247411_DP
         casimir_omega(27) = 0.0730288564124510_DP ;  casimir_omega_weight(27) = 0.0061215087631281_DP
         casimir_omega(28) = 0.0793081563295823_DP ;  casimir_omega_weight(28) = 0.0064387845670285_DP
         casimir_omega(29) = 0.0859099223105111_DP ;  casimir_omega_weight(29) = 0.0067665519383376_DP
         casimir_omega(30) = 0.0928449855993366_DP ;  casimir_omega_weight(30) = 0.0071054968707285_DP
         casimir_omega(31) = 0.1001248839004638_DP ;  casimir_omega_weight(31) = 0.0074563468461148_DP
         casimir_omega(32) = 0.1077619045792775_DP ;  casimir_omega_weight(32) = 0.0078198743079909_DP
         casimir_omega(33) = 0.1157691314865709_DP ;  casimir_omega_weight(33) = 0.0081969004407073_DP
         casimir_omega(34) = 0.1241604957285527_DP ;  casimir_omega_weight(34) = 0.0085882992870415_DP
         casimir_omega(35) = 0.1329508307384428_DP ;  casimir_omega_weight(35) = 0.0089950022401555_DP
         casimir_omega(36) = 0.1421559320438287_DP ;  casimir_omega_weight(36) = 0.0094180029502587_DP
         casimir_omega(37) = 0.1517926221666034_DP ;  casimir_omega_weight(37) = 0.0098583626910535_DP
         casimir_omega(38) = 0.1618788211400433_DP ;  casimir_omega_weight(38) = 0.0103172162364714_DP
         casimir_omega(39) = 0.1724336231810907_DP ;  casimir_omega_weight(39) = 0.0107957783043472_DP
         casimir_omega(40) = 0.1834773801159744_DP ;  casimir_omega_weight(40) = 0.0112953506306489_DP
         casimir_omega(41) = 0.1950317922247996_DP ;  casimir_omega_weight(41) = 0.0118173297458447_DP
         casimir_omega(42) = 0.2071200072467693_DP ;  casimir_omega_weight(42) = 0.0123632155340429_DP
         casimir_omega(43) = 0.2197667283733600_DP ;  casimir_omega_weight(43) = 0.0129346206658693_DP
         casimir_omega(44) = 0.2329983321535463_DP ;  casimir_omega_weight(44) = 0.0135332810079021_DP
         casimir_omega(45) = 0.2468429973446036_DP ;  casimir_omega_weight(45) = 0.0141610671249900_DP
         casimir_omega(46) = 0.2613308458659130_DP ;  casimir_omega_weight(46) = 0.0148199970073541_DP
         casimir_omega(47) = 0.2764940971538407_DP ;  casimir_omega_weight(47) = 0.0155122501722014_DP
         casimir_omega(48) = 0.2923672373755098_DP ;  casimir_omega_weight(48) = 0.0162401833101590_DP
         casimir_omega(49) = 0.3089872051411659_DP ;  casimir_omega_weight(49) = 0.0170063476705158_DP
         casimir_omega(50) = 0.3263935955621960_DP ;  casimir_omega_weight(50) = 0.0178135084066936_DP
         casimir_omega(51) = 0.3446288847387360_DP ;  casimir_omega_weight(51) = 0.0186646661350314_DP
         casimir_omega(52) = 0.3637386770318435_DP ;  casimir_omega_weight(52) = 0.0195630809968153_DP
         casimir_omega(53) = 0.3837719777859489_DP ;  casimir_omega_weight(53) = 0.0205122995561698_DP
         casimir_omega(54) = 0.4047814945242264_DP ;  casimir_omega_weight(54) = 0.0215161849163059_DP
         casimir_omega(55) = 0.4268239700502838_DP ;  casimir_omega_weight(55) = 0.0225789504946933_DP
         casimir_omega(56) = 0.4499605513632802_DP ;  casimir_omega_weight(56) = 0.0237051979657518_DP
         casimir_omega(57) = 0.4742571988409528_DP ;  casimir_omega_weight(57) = 0.0248999599593897_DP
         casimir_omega(58) = 0.4997851407789684_DP ;  casimir_omega_weight(58) = 0.0261687481974153_DP
         casimir_omega(59) = 0.5266213791106874_DP ;  casimir_omega_weight(59) = 0.0275176078602185_DP
         casimir_omega(60) = 0.5548492529872152_DP ;  casimir_omega_weight(60) = 0.0289531791065324_DP
         casimir_omega(61) = 0.5845590678953967_DP ;  casimir_omega_weight(61) = 0.0304827668233681_DP
         casimir_omega(62) = 0.6158487991574865_DP ;  casimir_omega_weight(62) = 0.0321144198665103_DP
         casimir_omega(63) = 0.6488248800224935_DP ;  casimir_omega_weight(63) = 0.0338570212700612_DP
         casimir_omega(64) = 0.6836030861639854_DP ;  casimir_omega_weight(64) = 0.0357203911638797_DP
         casimir_omega(65) = 0.7203095302891588_DP ;  casimir_omega_weight(65) = 0.0377154044495122_DP
         casimir_omega(66) = 0.7590817827959423_DP ;  casimir_omega_weight(66) = 0.0398541256594823_DP
         casimir_omega(67) = 0.8000701370581940_DP ;  casimir_omega_weight(67) = 0.0421499638756976_DP
         casimir_omega(68) = 0.8434390410585160_DP ;  casimir_omega_weight(68) = 0.0446178511274997_DP
         casimir_omega(69) = 0.8893687208283527_DP ;  casimir_omega_weight(69) = 0.0472744483504962_DP
         casimir_omega(70) = 0.9380570256247107_DP ;  casimir_omega_weight(70) = 0.0501383837907946_DP
         casimir_omega(71) = 0.9897215301315988_DP ;  casimir_omega_weight(71) = 0.0532305297205476_DP
         casimir_omega(72) = 1.0446019354208118_DP ;  casimir_omega_weight(72) = 0.0565743245329494_DP
         casimir_omega(73) = 1.1029628181886293_DP ;  casimir_omega_weight(73) = 0.0601961487639783_DP
         casimir_omega(74) = 1.1650967872133347_DP ;  casimir_omega_weight(74) = 0.0641257654151001_DP
         casimir_omega(75) = 1.2313281174443755_DP ;  casimir_omega_weight(75) = 0.0683968372166332_DP
         casimir_omega(76) = 1.3020169461328364_DP ;  casimir_omega_weight(76) = 0.0730475362937721_DP
         casimir_omega(77) = 1.3775641325735912_DP ;  casimir_omega_weight(77) = 0.0781212652281150_DP
         casimir_omega(78) = 1.4584169041563844_DP ;  casimir_omega_weight(78) = 0.0836675129460798_DP
         casimir_omega(79) = 1.5450754375475890_DP ;  casimir_omega_weight(79) = 0.0897428744724214_DP
         casimir_omega(80) = 1.6381005562789217_DP ;  casimir_omega_weight(80) = 0.0964122707056027_DP
         casimir_omega(81) = 1.7381227665325594_DP ;  casimir_omega_weight(81) = 0.1037504134579651_DP
         casimir_omega(82) = 1.8458529037412157_DP ;  casimir_omega_weight(82) = 0.1118435726657969_DP
         casimir_omega(83) = 1.9620947267311530_DP ;  casimir_omega_weight(83) = 0.1207917177309077_DP
         casimir_omega(84) = 2.0877598774452784_DP ;  casimir_omega_weight(84) = 0.1307111245116017_DP
         casimir_omega(85) = 2.2238857280073727_DP ;  casimir_omega_weight(85) = 0.1417375650469082_DP
         casimir_omega(86) = 2.3716567700166209_DP ;  casimir_omega_weight(86) = 0.1540302307437127_DP
         casimir_omega(87) = 2.5324303729302531_DP ;  casimir_omega_weight(87) = 0.1677765843512468_DP
         casimir_omega(88) = 2.7077679620388113_DP ;  casimir_omega_weight(88) = 0.1831983955954268_DP
         casimir_omega(89) = 2.8994729594754070_DP ;  casimir_omega_weight(89) = 0.2005592954871770_DP
         casimir_omega(90) = 3.1096372182921574_DP ;  casimir_omega_weight(90) = 0.2201742930758347_DP
         casimir_omega(91) = 3.3406981939072660_DP ;  casimir_omega_weight(91) = 0.2424218473056835_DP
         casimir_omega(92) = 3.5955097871362653_DP ;  casimir_omega_weight(92) = 0.2677592923667230_DP
         casimir_omega(93) = 3.8774307268843149_DP ;  casimir_omega_weight(93) = 0.2967427020263470_DP
         casimir_omega(94) = 4.1904356367454678_DP ;  casimir_omega_weight(94) = 0.3300526832955611_DP
         casimir_omega(95) = 4.5392556914819888_DP ;  casimir_omega_weight(95) = 0.3685281671490563_DP
         casimir_omega(96) = 4.9295582278708858_DP ;  casimir_omega_weight(96) = 0.4132110972658836_DP
         casimir_omega(97) = 5.3681781452858965_DP ;  casimir_omega_weight(97) = 0.4654061357613040_DP
         casimir_omega(98) = 5.8634188928075401_DP ;  casimir_omega_weight(98) = 0.5267613098673054_DP
         casimir_omega(99) = 6.4254480271920782_DP ;  casimir_omega_weight(99) = 0.5993782381239071_DP
         casimir_omega(100) = 7.0668228897698055_DP ;  casimir_omega_weight(100) = 0.6859647227295468_DP
         casimir_omega(101) = 7.8031977203601421_DP ;  casimir_omega_weight(101) = 0.7900489433073745_DP
         casimir_omega(102) = 8.6542874706449666_DP ;  casimir_omega_weight(102) = 0.9162847015438560_DP
         casimir_omega(103) = 9.6452006124694449_DP ;  casimir_omega_weight(103) = 1.0708936777932563_DP
         casimir_omega(104) = 10.8083116809133060_DP ;  casimir_omega_weight(104) = 1.2623179542170659_DP
         casimir_omega(105) = 12.1859385957337398_DP ;  casimir_omega_weight(105) = 1.5022022938638520_DP
         casimir_omega(106) = 13.8342457423025049_DP ;  casimir_omega_weight(106) = 1.8069061352295328_DP
         casimir_omega(107) = 15.8290587520334487_DP ;  casimir_omega_weight(107) = 2.1998895915896521_DP
         casimir_omega(108) = 18.2747410020290637_DP ;  casimir_omega_weight(108) = 2.7155854563590824_DP
         casimir_omega(109) = 21.3181229319876415_DP ;  casimir_omega_weight(109) = 3.4058849200336678_DP
         casimir_omega(110) = 25.1710597688688047_DP ;  casimir_omega_weight(110) = 4.3514016682137564_DP
         casimir_omega(111) = 30.1483133788109647_DP ;  casimir_omega_weight(111) = 5.6818691080280521_DP
         casimir_omega(112) = 36.7339207047230900_DP ;  casimir_omega_weight(112) = 7.6149215756281796_DP
         casimir_omega(113) = 45.7034419074545113_DP ;  casimir_omega_weight(113) = 10.5342108108810404_DP
         casimir_omega(114) = 58.3631096047622506_DP ;  casimir_omega_weight(114) = 15.1580817255117832_DP
         casimir_omega(115) = 77.0535298965956628_DP ;  casimir_omega_weight(115) = 22.9362912066841069_DP
         casimir_omega(116) = 106.3139760039738064_DP ;  casimir_omega_weight(116) = 37.0894324143586473_DP
         casimir_omega(117) = 155.9202545596408811_DP ;  casimir_omega_weight(117) = 65.7454929476414378_DP
         casimir_omega(118) = 250.2382491135264786_DP ;  casimir_omega_weight(118) = 133.4352339440673916_DP
         casimir_omega(119) = 464.9518069878388360_DP ;  casimir_omega_weight(119) = 337.3244144872141419_DP
         casimir_omega(120) = 1143.2612068325900054_DP ;  casimir_omega_weight(120) = 1296.6801236170986158_DP
         casimir_omega(121) = 6025.4826622272921668_DP ;  casimir_omega_weight(121) = 15464.3608243612707156_DP
return
endsubroutine gauss_legendre_grid120

subroutine gauss_legendre_grid130()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000509401064755_DP ;  casimir_omega_weight(2) = 0.0001307361624078_DP
         casimir_omega(3) = 0.0002684652502737_DP ;  casimir_omega_weight(3) = 0.0003044760072594_DP
         casimir_omega(4) = 0.0006600744904710_DP ;  casimir_omega_weight(4) = 0.0004788260261463_DP
         casimir_omega(5) = 0.0012263092053646_DP ;  casimir_omega_weight(5) = 0.0006537543698638_DP
         casimir_omega(6) = 0.0019678380088815_DP ;  casimir_omega_weight(6) = 0.0008294474064794_DP
         casimir_omega(7) = 0.0028855265207818_DP ;  casimir_omega_weight(7) = 0.0010061080868241_DP
         casimir_omega(8) = 0.0039804459112934_DP ;  casimir_omega_weight(8) = 0.0011839440657306_DP
         casimir_omega(9) = 0.0052538767541422_DP ;  casimir_omega_weight(9) = 0.0013631664705258_DP
         casimir_omega(10) = 0.0067073125351913_DP ;  casimir_omega_weight(10) = 0.0015439900594873_DP
         casimir_omega(11) = 0.0083424634952359_DP ;  casimir_omega_weight(11) = 0.0017266336683959_DP
         casimir_omega(12) = 0.0101612609703706_DP ;  casimir_omega_weight(12) = 0.0019113207453236_DP
         casimir_omega(13) = 0.0121658622899480_DP ;  casimir_omega_weight(13) = 0.0020982799275607_DP
         casimir_omega(14) = 0.0143586562671857_DP ;  casimir_omega_weight(14) = 0.0022877456491990_DP
         casimir_omega(15) = 0.0167422693117178_DP ;  casimir_omega_weight(15) = 0.0024799587774739_DP
         casimir_omega(16) = 0.0193195721933215_DP ;  casimir_omega_weight(16) = 0.0026751672791158_DP
         casimir_omega(17) = 0.0220936874880077_DP ;  casimir_omega_weight(17) = 0.0028736269192063_DP
         casimir_omega(18) = 0.0250679977405147_DP ;  casimir_omega_weight(18) = 0.0030756019956877_DP
         casimir_omega(19) = 0.0282461543805989_DP ;  casimir_omega_weight(19) = 0.0032813661130745_DP
         casimir_omega(20) = 0.0316320874342566_DP ;  casimir_omega_weight(20) = 0.0034912029992874_DP
         casimir_omega(21) = 0.0352300160751045_DP ;  casimir_omega_weight(21) = 0.0037054073698629_DP
         casimir_omega(22) = 0.0390444600655652_DP ;  casimir_omega_weight(22) = 0.0039242858441743_DP
         casimir_omega(23) = 0.0430802521423479_DP ;  casimir_omega_weight(23) = 0.0041481579186569_DP
         casimir_omega(24) = 0.0473425514059133_DP ;  casimir_omega_weight(24) = 0.0043773570025063_DP
         casimir_omega(25) = 0.0518368577793168_DP ;  casimir_omega_weight(25) = 0.0046122315217593_DP
         casimir_omega(26) = 0.0565690276079879_DP ;  casimir_omega_weight(26) = 0.0048531460982109_DP
         casimir_omega(27) = 0.0615452904787504_DP ;  casimir_omega_weight(27) = 0.0051004828102056_DP
         casimir_omega(28) = 0.0667722673437448_DP ;  casimir_omega_weight(28) = 0.0053546425429854_DP
         casimir_omega(29) = 0.0722569900429482_DP ;  casimir_omega_weight(29) = 0.0056160464370057_DP
         casimir_omega(30) = 0.0780069223278019_DP ;  casimir_omega_weight(30) = 0.0058851374434384_DP
         casimir_omega(31) = 0.0840299824980989_DP ;  casimir_omega_weight(31) = 0.0061623819969591_DP
         casimir_omega(32) = 0.0903345677748912_DP ;  casimir_omega_weight(32) = 0.0064482718169319_DP
         casimir_omega(33) = 0.0969295805438144_DP ;  casimir_omega_weight(33) = 0.0067433258492015_DP
         casimir_omega(34) = 0.1038244566160639_DP ;  casimir_omega_weight(34) = 0.0070480923619488_DP
         casimir_omega(35) = 0.1110291956683713_DP ;  casimir_omega_weight(35) = 0.0073631512104376_DP
         casimir_omega(36) = 0.1185543940389318_DP ;  casimir_omega_weight(36) = 0.0076891162870362_DP
         casimir_omega(37) = 0.1264112800734594_DP ;  casimir_omega_weight(37) = 0.0080266381746194_DP
         casimir_omega(38) = 0.1346117522345877_DP ;  casimir_omega_weight(38) = 0.0083764070233748_DP
         casimir_omega(39) = 0.1431684202089427_DP ;  casimir_omega_weight(39) = 0.0087391556732325_DP
         casimir_omega(40) = 0.1520946492696037_DP ;  casimir_omega_weight(40) = 0.0091156630465326_DP
         casimir_omega(41) = 0.1614046081776279_DP ;  casimir_omega_weight(41) = 0.0095067578382930_DP
         casimir_omega(42) = 0.1711133209351899_DP ;  casimir_omega_weight(42) = 0.0099133225345105_DP
         casimir_omega(43) = 0.1812367227349856_DP ;  casimir_omega_weight(43) = 0.0103362977923533_DP
         casimir_omega(44) = 0.1917917204863404_DP ;  casimir_omega_weight(44) = 0.0107766872200146_DP
         casimir_omega(45) = 0.2027962583383783_DP ;  casimir_omega_weight(45) = 0.0112355625983994_DP
         casimir_omega(46) = 0.2142693886651877_DP ;  casimir_omega_weight(46) = 0.0117140695917364_DP
         casimir_omega(47) = 0.2262313490277855_DP ;  casimir_omega_weight(47) = 0.0122134339998550_DP
         casimir_omega(48) = 0.2387036456835263_DP ;  casimir_omega_weight(48) = 0.0127349686112127_DP
         casimir_omega(49) = 0.2517091442762089_DP ;  casimir_omega_weight(49) = 0.0132800807229735_DP
         casimir_omega(50) = 0.2652721684104634_DP ;  casimir_omega_weight(50) = 0.0138502804026399_DP
         casimir_omega(51) = 0.2794186068930491_DP ;  casimir_omega_weight(51) = 0.0144471895750629_DP
         casimir_omega(52) = 0.2941760305127345_DP ;  casimir_omega_weight(52) = 0.0150725520292967_DP
         casimir_omega(53) = 0.3095738193308177_DP ;  casimir_omega_weight(53) = 0.0157282444518988_DP
         casimir_omega(54) = 0.3256433015677402_DP ;  casimir_omega_weight(54) = 0.0164162886071410_DP
         casimir_omega(55) = 0.3424179052994672_DP ;  casimir_omega_weight(55) = 0.0171388648005193_DP
         casimir_omega(56) = 0.3599333243226112_DP ;  casimir_omega_weight(56) = 0.0178983267801316_DP
         casimir_omega(57) = 0.3782276997120604_DP ;  casimir_omega_weight(57) = 0.0186972182514964_DP
         casimir_omega(58) = 0.3973418187822377_DP ;  casimir_omega_weight(58) = 0.0195382912054479_DP
         casimir_omega(59) = 0.4173193333763116_DP ;  casimir_omega_weight(59) = 0.0204245262865987_DP
         casimir_omega(60) = 0.4382069996508460_DP ;  casimir_omega_weight(60) = 0.0213591554619585_DP
         casimir_omega(61) = 0.4600549418011266_DP ;  casimir_omega_weight(61) = 0.0223456872865272_DP
         casimir_omega(62) = 0.4829169424902232_DP ;  casimir_omega_weight(62) = 0.0233879351057914_DP
         casimir_omega(63) = 0.5068507631092718_DP ;  casimir_omega_weight(63) = 0.0244900485852628_DP
         casimir_omega(64) = 0.5319184974149963_DP ;  casimir_omega_weight(64) = 0.0256565490155477_DP
         casimir_omega(65) = 0.5581869625722418_DP ;  casimir_omega_weight(65) = 0.0268923689096659_DP
         casimir_omega(66) = 0.5857281321848004_DP ;  casimir_omega_weight(66) = 0.0282028964890731_DP
         casimir_omega(67) = 0.6146196165398081_DP ;  casimir_omega_weight(67) = 0.0295940257483793_DP
         casimir_omega(68) = 0.6449451960343994_DP ;  casimir_omega_weight(68) = 0.0310722128986101_DP
         casimir_omega(69) = 0.6767954146161849_DP ;  casimir_omega_weight(69) = 0.0326445401184283_DP
         casimir_omega(70) = 0.7102682410727437_DP ;  casimir_omega_weight(70) = 0.0343187876954837_DP
         casimir_omega(71) = 0.7454698071755649_DP ;  casimir_omega_weight(71) = 0.0361035158212576_DP
         casimir_omega(72) = 0.7825152330514941_DP ;  casimir_omega_weight(72) = 0.0380081575175672_DP
         casimir_omega(73) = 0.8215295517571395_DP ;  casimir_omega_weight(73) = 0.0400431244287631_DP
         casimir_omega(74) = 0.8626487469138547_DP ;  casimir_omega_weight(74) = 0.0422199275190422_DP
         casimir_omega(75) = 0.9060209194776379_DP ;  casimir_omega_weight(75) = 0.0445513150798830_DP
         casimir_omega(76) = 0.9518076023360091_DP ;  casimir_omega_weight(76) = 0.0470514308916528_DP
         casimir_omega(77) = 1.0001852445241470_DP ;  casimir_omega_weight(77) = 0.0497359959121603_DP
         casimir_omega(78) = 1.0513468905346983_DP ;  casimir_omega_weight(78) = 0.0526225175040504_DP
         casimir_omega(79) = 1.1055040845822932_DP ;  casimir_omega_weight(79) = 0.0557305309874475_DP
         casimir_omega(80) = 1.1628890349260945_DP ;  casimir_omega_weight(80) = 0.0590818792470453_DP
         casimir_omega(81) = 1.2237570796388051_DP ;  casimir_omega_weight(81) = 0.0627010372732505_DP
         casimir_omega(82) = 1.2883895027713541_DP ;  casimir_omega_weight(82) = 0.0666154899275712_DP
         casimir_omega(83) = 1.3570967589896648_DP ;  casimir_omega_weight(83) = 0.0708561729568141_DP
         casimir_omega(84) = 1.4302221758179767_DP ;  casimir_omega_weight(84) = 0.0754579894237264_DP
         casimir_omega(85) = 1.5081462160711576_DP ;  casimir_omega_weight(85) = 0.0804604163785956_DP
         casimir_omega(86) = 1.5912913994770239_DP ;  casimir_omega_weight(86) = 0.0859082199066155_DP
         casimir_omega(87) = 1.6801280026169658_DP ;  casimir_omega_weight(87) = 0.0918523008269468_DP
         casimir_omega(88) = 1.7751806810918416_DP ;  casimir_omega_weight(88) = 0.0983506985251994_DP
         casimir_omega(89) = 1.8770361884606979_DP ;  casimir_omega_weight(89) = 0.1054697869772156_DP
         casimir_omega(90) = 1.9863524045643450_DP ;  casimir_omega_weight(90) = 0.1132857053708507_DP
         casimir_omega(91) = 2.1038689333623086_DP ;  casimir_omega_weight(91) = 0.1218860763894478_DP
         casimir_omega(92) = 2.2304195900269246_DP ;  casimir_omega_weight(92) = 0.1313720788989791_DP
         casimir_omega(93) = 2.3669471722299940_DP ;  casimir_omega_weight(93) = 0.1418609594394369_DP
         casimir_omega(94) = 2.5145210059216216_DP ;  casimir_omega_weight(94) = 0.1534890898585870_DP
         casimir_omega(95) = 2.6743578775546184_DP ;  casimir_omega_weight(95) = 0.1664157084110141_DP
         casimir_omega(96) = 2.8478471208487011_DP ;  casimir_omega_weight(96) = 0.1808275211072889_DP
         casimir_omega(97) = 3.0365808278837934_DP ;  casimir_omega_weight(97) = 0.1969443924019973_DP
         casimir_omega(98) = 3.2423904166186213_DP ;  casimir_omega_weight(98) = 0.2150264241501282_DP
         casimir_omega(99) = 3.4673911305046015_DP ;  casimir_omega_weight(99) = 0.2353828157576759_DP
         casimir_omega(100) = 3.7140364992838490_DP ;  casimir_omega_weight(100) = 0.2583830260069861_DP
         casimir_omega(101) = 3.9851853932273285_DP ;  casimir_omega_weight(101) = 0.2844709316640776_DP
         casimir_omega(102) = 4.2841851122382888_DP ;  casimir_omega_weight(102) = 0.3141829192680430_DP
         casimir_omega(103) = 4.6149750465375634_DP ;  casimir_omega_weight(103) = 0.3481711832288603_DP
         casimir_omega(104) = 4.9822169424165406_DP ;  casimir_omega_weight(104) = 0.3872339782104995_DP
         casimir_omega(105) = 5.3914598728048837_DP ;  casimir_omega_weight(105) = 0.4323552509472157_DP
         casimir_omega(106) = 5.8493508958950571_DP ;  casimir_omega_weight(106) = 0.4847570539239716_DP
         casimir_omega(107) = 6.3639064559272311_DP ;  casimir_omega_weight(107) = 0.5459695719005782_DP
         casimir_omega(108) = 6.9448653992997427_DP ;  casimir_omega_weight(108) = 0.6179257092586684_DP
         casimir_omega(109) = 7.6041529091529609_DP ;  casimir_omega_weight(109) = 0.7030903700059488_DP
         casimir_omega(110) = 8.3564970513744488_DP ;  casimir_omega_weight(110) = 0.8046394273959269_DP
         casimir_omega(111) = 9.2202581210104810_DP ;  casimir_omega_weight(111) = 0.9267109434515546_DP
         casimir_omega(112) = 10.2185590614701400_DP ;  casimir_omega_weight(112) = 1.0747631785076868_DP
         casimir_omega(113) = 11.3808486635038850_DP ;  casimir_omega_weight(113) = 1.2560932967523533_DP
         casimir_omega(114) = 12.7450977980658742_DP ;  casimir_omega_weight(114) = 1.4806026851966603_DP
         casimir_omega(115) = 14.3609395423778494_DP ;  casimir_omega_weight(115) = 1.7619490305403176_DP
         casimir_omega(116) = 16.2942469515514716_DP ;  casimir_omega_weight(116) = 2.1193196786906254_DP
         casimir_omega(117) = 18.6339529880712895_DP ;  casimir_omega_weight(117) = 2.5802300804310265_DP
         casimir_omega(118) = 21.5024614224808062_DP ;  casimir_omega_weight(118) = 3.1850651156740892_DP
         casimir_omega(119) = 25.0719839866015946_DP ;  casimir_omega_weight(119) = 3.9946859382103654_DP
         casimir_omega(120) = 29.5909974500896418_DP ;  casimir_omega_weight(120) = 5.1036411974948193_DP
         casimir_omega(121) = 35.4286737689081406_DP ;  casimir_omega_weight(121) = 6.6640901509444461_DP
         casimir_omega(122) = 43.1527210404443409_DP ;  casimir_omega_weight(122) = 8.9312876315105925_DP
         casimir_omega(123) = 53.6727635861861927_DP ;  casimir_omega_weight(123) = 12.3552038178441741_DP
         casimir_omega(124) = 68.5208307781784498_DP ;  casimir_omega_weight(124) = 17.7783574720811615_DP
         casimir_omega(125) = 90.4421283501396545_DP ;  casimir_omega_weight(125) = 26.9011119704959505_DP
         casimir_omega(126) = 124.7605930519875130_DP ;  casimir_omega_weight(126) = 43.5007755716523192_DP
         casimir_omega(127) = 182.9418876834371304_DP ;  casimir_omega_weight(127) = 77.1103483063990751_DP
         casimir_omega(128) = 293.5638079084578180_DP ;  casimir_omega_weight(128) = 156.5010043262410022_DP
         casimir_omega(129) = 545.3929900294398294_DP ;  casimir_omega_weight(129) = 395.6346774098458923_DP
         casimir_omega(130) = 1340.9556716668651006_DP ;  casimir_omega_weight(130) = 1520.8256130162380941_DP
         casimir_omega(131) = 7067.1230373831776888_DP ;  casimir_omega_weight(131) = 18137.5464068574874545_DP
return
endsubroutine gauss_legendre_grid130

subroutine gauss_legendre_grid140()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000439465858040_DP ;  casimir_omega_weight(2) = 0.0001127866377870_DP
         casimir_omega(3) = 0.0002316002176323_DP ;  casimir_omega_weight(3) = 0.0002626553059691_DP
         casimir_omega(4) = 0.0005694005732888_DP ;  casimir_omega_weight(4) = 0.0004130085023146_DP
         casimir_omega(5) = 0.0010577607706115_DP ;  casimir_omega_weight(5) = 0.0005637945908504_DP
         casimir_omega(6) = 0.0016971793946794_DP ;  casimir_omega_weight(6) = 0.0007151498924354_DP
         casimir_omega(7) = 0.0024883003300048_DP ;  casimir_omega_weight(7) = 0.0008672246174906_DP
         casimir_omega(8) = 0.0034319196578740_DP ;  casimir_omega_weight(8) = 0.0010201725023574_DP
         casimir_omega(9) = 0.0045289883890402_DP ;  casimir_omega_weight(9) = 0.0011741496281056_DP
         casimir_omega(10) = 0.0057806147786660_DP ;  casimir_omega_weight(10) = 0.0013293144345932_DP
         casimir_omega(11) = 0.0071880668046461_DP ;  casimir_omega_weight(11) = 0.0014858279784129_DP
         casimir_omega(12) = 0.0087527749452324_DP ;  casimir_omega_weight(12) = 0.0016438542615170_DP
         casimir_omega(13) = 0.0104763353018858_DP ;  casimir_omega_weight(13) = 0.0018035605899057_DP
         casimir_omega(14) = 0.0123605130908658_DP ;  casimir_omega_weight(14) = 0.0019651179515824_DP
         casimir_omega(15) = 0.0144072465211868_DP ;  casimir_omega_weight(15) = 0.0021287014111972_DP
         casimir_omega(16) = 0.0166186510755181_DP ;  casimir_omega_weight(16) = 0.0022944905214159_DP
         casimir_omega(17) = 0.0189970242112501_DP ;  casimir_omega_weight(17) = 0.0024626697521010_DP
         casimir_omega(18) = 0.0215448505002501_DP ;  casimir_omega_weight(18) = 0.0026334289387947_DP
         casimir_omega(19) = 0.0242648072275069_DP ;  casimir_omega_weight(19) = 0.0028069637523223_DP
         casimir_omega(20) = 0.0271597704707348_DP ;  casimir_omega_weight(20) = 0.0029834761914565_DP
         casimir_omega(21) = 0.0302328216850810_DP ;  casimir_omega_weight(21) = 0.0031631751008300_DP
         casimir_omega(22) = 0.0334872548193096_DP ;  casimir_omega_weight(22) = 0.0033462767163790_DP
         casimir_omega(23) = 0.0369265839922514_DP ;  casimir_omega_weight(23) = 0.0035330052408599_DP
         casimir_omega(24) = 0.0405545517609020_DP ;  casimir_omega_weight(24) = 0.0037235934521046_DP
         casimir_omega(25) = 0.0443751380143585_DP ;  casimir_omega_weight(25) = 0.0039182833469342_DP
         casimir_omega(26) = 0.0483925695307975_DP ;  casimir_omega_weight(26) = 0.0041173268238842_DP
         casimir_omega(27) = 0.0526113302379824_DP ;  casimir_omega_weight(27) = 0.0043209864081181_DP
         casimir_omega(28) = 0.0570361722213063_DP ;  casimir_omega_weight(28) = 0.0045295360222253_DP
         casimir_omega(29) = 0.0616721275272214_DP ;  casimir_omega_weight(29) = 0.0047432618068952_DP
         casimir_omega(30) = 0.0665245208140690_DP ;  casimir_omega_weight(30) = 0.0049624629957860_DP
         casimir_omega(31) = 0.0715989829068268_DP ;  casimir_omega_weight(31) = 0.0051874528493212_DP
         casimir_omega(32) = 0.0769014653172361_DP ;  casimir_omega_weight(32) = 0.0054185596525334_DP
         casimir_omega(33) = 0.0824382557960986_DP ;  casimir_omega_weight(33) = 0.0056561277825441_DP
         casimir_omega(34) = 0.0882159949904001_DP ;  casimir_omega_weight(34) = 0.0059005188517912_DP
         casimir_omega(35) = 0.0942416942842790_DP ;  casimir_omega_weight(35) = 0.0061521129336558_DP
         casimir_omega(36) = 0.1005227549098414_DP ;  casimir_omega_weight(36) = 0.0064113098777966_DP
         casimir_omega(37) = 0.1070669884214420_DP ;  casimir_omega_weight(37) = 0.0066785307231434_DP
         casimir_omega(38) = 0.1138826386354160_DP ;  casimir_omega_weight(38) = 0.0069542192173226_DP
         casimir_omega(39) = 0.1209784051463965_DP ;  casimir_omega_weight(39) = 0.0072388434520855_DP
         casimir_omega(40) = 0.1283634685414191_DP ;  casimir_omega_weight(40) = 0.0075328976252870_DP
         casimir_omega(41) = 0.1360475174440520_DP ;  casimir_omega_weight(41) = 0.0078369039409862_DP
         casimir_omega(42) = 0.1440407775329536_DP ;  casimir_omega_weight(42) = 0.0081514146604250_DP
         casimir_omega(43) = 0.1523540426926295_DP ;  casimir_omega_weight(43) = 0.0084770143179081_DP
         casimir_omega(44) = 0.1609987084689156_DP ;  casimir_omega_weight(44) = 0.0088143221170763_DP
         casimir_omega(45) = 0.1699868080179786_DP ;  casimir_omega_weight(45) = 0.0091639945246675_DP
         casimir_omega(46) = 0.1793310507556035_DP ;  casimir_omega_weight(46) = 0.0095267280806367_DP
         casimir_omega(47) = 0.1890448639334011_DP ;  casimir_omega_weight(47) = 0.0099032624455332_DP
         casimir_omega(48) = 0.1991424373905697_DP ;  casimir_omega_weight(48) = 0.0102943837082912_DP
         casimir_omega(49) = 0.2096387717542136_DP ;  casimir_omega_weight(49) = 0.0107009279800546_DP
         casimir_omega(50) = 0.2205497303882560_DP ;  casimir_omega_weight(50) = 0.0111237853025402_DP
         casimir_omega(51) = 0.2318920954210147_DP ;  casimir_omega_weight(51) = 0.0115639039025553_DP
         casimir_omega(52) = 0.2436836282148911_DP ;  casimir_omega_weight(52) = 0.0120222948278862_DP
         casimir_omega(53) = 0.2559431346787718_DP ;  casimir_omega_weight(53) = 0.0125000370037428_DP
         casimir_omega(54) = 0.2686905358651805_DP ;  casimir_omega_weight(54) = 0.0129982827534945_DP
         casimir_omega(55) = 0.2819469443404069_DP ;  casimir_omega_weight(55) = 0.0135182638324780_DP
         casimir_omega(56) = 0.2957347468674686_DP ;  casimir_omega_weight(56) = 0.0140612980294812_DP
         casimir_omega(57) = 0.3100776939995209_DP ;  casimir_omega_weight(57) = 0.0146287963969080_DP
         casimir_omega(58) = 0.3250009972459952_DP ;  casimir_omega_weight(58) = 0.0152222711781385_DP
         casimir_omega(59) = 0.3405314345462871_DP ;  casimir_omega_weight(59) = 0.0158433445088412_DP
         casimir_omega(60) = 0.3566974648672847_DP ;  casimir_omega_weight(60) = 0.0164937579785685_DP
         casimir_omega(61) = 0.3735293528326142_DP ;  casimir_omega_weight(61) = 0.0171753831497755_DP
         casimir_omega(62) = 0.3910593043946636_DP ;  casimir_omega_weight(62) = 0.0178902331437008_DP
         casimir_omega(63) = 0.4093216146767527_DP ;  casimir_omega_weight(63) = 0.0186404754166559_DP
         casimir_omega(64) = 0.4283528292442397_DP ;  casimir_omega_weight(64) = 0.0194284458663249_DP
         casimir_omega(65) = 0.4481919202119408_DP ;  casimir_omega_weight(65) = 0.0202566644261125_DP
         casimir_omega(66) = 0.4688804787636278_DP ;  casimir_omega_weight(66) = 0.0211278523267191_DP
         casimir_omega(67) = 0.4904629258503658_DP ;  casimir_omega_weight(67) = 0.0220449512283819_DP
         casimir_omega(68) = 0.5129867430515824_DP ;  casimir_omega_weight(68) = 0.0230111444551947_DP
         casimir_omega(69) = 0.5365027258298781_DP ;  casimir_omega_weight(69) = 0.0240298805951491_DP
         casimir_omega(70) = 0.5610652616924049_DP ;  casimir_omega_weight(70) = 0.0251048997667956_DP
         casimir_omega(71) = 0.5867326360935822_DP ;  casimir_omega_weight(71) = 0.0262402628965107_DP
         casimir_omega(72) = 0.6135673692822848_DP ;  casimir_omega_weight(72) = 0.0274403844004336_DP
         casimir_omega(73) = 0.6416365877190305_DP ;  casimir_omega_weight(73) = 0.0287100687232111_DP
         casimir_omega(74) = 0.6710124341738268_DP ;  casimir_omega_weight(74) = 0.0300545512534257_DP
         casimir_omega(75) = 0.7017725211737108_DP ;  casimir_omega_weight(75) = 0.0314795442146348_DP
         casimir_omega(76) = 0.7340004331129226_DP ;  casimir_omega_weight(76) = 0.0329912882233247_DP
         casimir_omega(77) = 0.7677862830827785_DP ;  casimir_omega_weight(77) = 0.0345966103136292_DP
         casimir_omega(78) = 0.8032273313400269_DP ;  casimir_omega_weight(78) = 0.0363029893558615_DP
         casimir_omega(79) = 0.8404286733325939_DP ;  casimir_omega_weight(79) = 0.0381186299461539_DP
         casimir_omega(80) = 0.8795040063650111_DP ;  casimir_omega_weight(80) = 0.0400525460216514_DP
         casimir_omega(81) = 0.9205764853421869_DP ;  casimir_omega_weight(81) = 0.0421146556655258_DP
         casimir_omega(82) = 0.9637796796154948_DP ;  casimir_omega_weight(82) = 0.0443158888152541_DP
         casimir_omega(83) = 1.0092586448124725_DP ;  casimir_omega_weight(83) = 0.0466683098841480_DP
         casimir_omega(84) = 1.0571711257131120_DP ;  casimir_omega_weight(84) = 0.0491852576599521_DP
         casimir_omega(85) = 1.1076889088051431_DP ;  casimir_omega_weight(85) = 0.0518815052683858_DP
         casimir_omega(86) = 1.1609993461850099_DP ;  casimir_omega_weight(86) = 0.0547734434980323_DP
         casimir_omega(87) = 1.2173070760647924_DP ;  casimir_omega_weight(87) = 0.0578792913962594_DP
         casimir_omega(88) = 1.2768359694132949_DP ;  casimir_omega_weight(88) = 0.0612193387862640_DP
         casimir_omega(89) = 1.3398313373442932_DP ;  casimir_omega_weight(89) = 0.0648162262534309_DP
         casimir_omega(90) = 1.4065624399413077_DP ;  casimir_omega_weight(90) = 0.0686952692417746_DP
         casimir_omega(91) = 1.4773253444935395_DP ;  casimir_omega_weight(91) = 0.0728848342349334_DP
         casimir_omega(92) = 1.5524461898815345_DP ;  casimir_omega_weight(92) = 0.0774167766308957_DP
         casimir_omega(93) = 1.6322849244306734_DP ;  casimir_omega_weight(93) = 0.0823269519304169_DP
         casimir_omega(94) = 1.7172395973683432_DP ;  casimir_omega_weight(94) = 0.0876558143427810_DP
         casimir_omega(95) = 1.8077512996084653_DP ;  casimir_omega_weight(95) = 0.0934491199926063_DP
         casimir_omega(96) = 1.9043098686184072_DP ;  casimir_omega_weight(96) = 0.0997587557480009_DP
         casimir_omega(97) = 2.0074604954532731_DP ;  casimir_omega_weight(97) = 0.1066437194909758_DP
         casimir_omega(98) = 2.1178114007642561_DP ;  casimir_omega_weight(98) = 0.1141712836847336_DP
         casimir_omega(99) = 2.2360427821040987_DP ;  casimir_omega_weight(99) = 0.1224183817153734_DP
         casimir_omega(100) = 2.3629172789742872_DP ;  casimir_omega_weight(100) = 0.1314732661627331_DP
         casimir_omega(101) = 2.4992922571362759_DP ;  casimir_omega_weight(101) = 0.1414375005081176_DP
         casimir_omega(102) = 2.6461342828107526_DP ;  casimir_omega_weight(102) = 0.1524283616411179_DP
         casimir_omega(103) = 2.8045362445456070_DP ;  casimir_omega_weight(103) = 0.1645817509967912_DP
         casimir_omega(104) = 2.9757376910727387_DP ;  casimir_omega_weight(104) = 0.1780557387418156_DP
         casimir_omega(105) = 3.1611490944858165_DP ;  casimir_omega_weight(105) = 0.1930349001841676_DP
         casimir_omega(106) = 3.3623809290586486_DP ;  casimir_omega_weight(106) = 0.2097356493230040_DP
         casimir_omega(107) = 3.5812786898138946_DP ;  casimir_omega_weight(107) = 0.2284128350813748_DP
         casimir_omega(108) = 3.8199652789991565_DP ;  casimir_omega_weight(108) = 0.2493679467195907_DP
         casimir_omega(109) = 4.0808925868735644_DP ;  casimir_omega_weight(109) = 0.2729593838804458_DP
         casimir_omega(110) = 4.3669046187781690_DP ;  casimir_omega_weight(110) = 0.2996153945697881_DP
         casimir_omega(111) = 4.6813152196114558_DP ;  casimir_omega_weight(111) = 0.3298504867903062_DP
         casimir_omega(112) = 5.0280043847616414_DP ;  casimir_omega_weight(112) = 0.3642863992366057_DP
         casimir_omega(113) = 5.4115384161303881_DP ;  casimir_omega_weight(113) = 0.4036791067668786_DP
         casimir_omega(114) = 5.8373209168290501_DP ;  casimir_omega_weight(114) = 0.4489538867808397_DP
         casimir_omega(115) = 6.3117840131901302_DP ;  casimir_omega_weight(115) = 0.5012512575584892_DP
         casimir_omega(116) = 6.8426325350751185_DP ;  casimir_omega_weight(116) = 0.5619877324154424_DP
         casimir_omega(117) = 7.4391586041095179_DP ;  casimir_omega_weight(117) = 0.6329369893933197_DP
         casimir_omega(118) = 8.1126508245115510_DP ;  casimir_omega_weight(118) = 0.7163395100853730_DP
         casimir_omega(119) = 8.8769320426090221_DP ;  casimir_omega_weight(119) = 0.8150524316851830_DP
         casimir_omega(120) = 9.7490740025002864_DP ;  casimir_omega_weight(120) = 0.9327569956537095_DP
         casimir_omega(121) = 10.7503586645871057_DP ;  casimir_omega_weight(121) = 1.0742497432572202_DP
         casimir_omega(122) = 11.9075885059596818_DP ;  casimir_omega_weight(122) = 1.2458574943920449_DP
         casimir_omega(123) = 13.2548984678610466_DP ;  casimir_omega_weight(123) = 1.4560385936124298_DP
         casimir_omega(124) = 14.8363016703425572_DP ;  casimir_omega_weight(124) = 1.7162700126445174_DP
         casimir_omega(125) = 16.7093292197977483_DP ;  casimir_omega_weight(125) = 2.0423827547450850_DP
         casimir_omega(126) = 18.9503364314716620_DP ;  casimir_omega_weight(126) = 2.4566174050715439_DP
         casimir_omega(127) = 21.6624079995481935_DP ;  casimir_omega_weight(127) = 2.9908678869395802_DP
         casimir_omega(128) = 24.9874255618937440_DP ;  casimir_omega_weight(128) = 3.6919454371494096_DP
         casimir_omega(129) = 29.1250045490455669_DP ;  casimir_omega_weight(129) = 4.6303959114414228_DP
         casimir_omega(130) = 34.3631613179849964_DP ;  casimir_omega_weight(130) = 5.9158132793364864_DP
         casimir_omega(131) = 41.1298133737674547_DP ;  casimir_omega_weight(131) = 7.7245695694155927_DP
         casimir_omega(132) = 50.0830069869838681_DP ;  casimir_omega_weight(132) = 10.3525377610860900_DP
         casimir_omega(133) = 62.2771130379795750_DP ;  casimir_omega_weight(133) = 14.3212908100512930_DP
         casimir_omega(134) = 79.4879494218164098_DP ;  casimir_omega_weight(134) = 20.6074156600535154_DP
         casimir_omega(135) = 104.8975605166207572_DP ;  casimir_omega_weight(135) = 31.1818508215582497_DP
         casimir_omega(136) = 144.6770695880233859_DP ;  casimir_omega_weight(136) = 50.4229794210273283_DP
         casimir_omega(137) = 212.1166454934597709_DP ;  casimir_omega_weight(137) = 89.3807671033167281_DP
         casimir_omega(138) = 340.3416065354004445_DP ;  casimir_omega_weight(138) = 181.4046825493679194_DP
         casimir_omega(139) = 632.2438312992250076_DP ;  casimir_omega_weight(139) = 458.5911748463315689_DP
         casimir_omega(140) = 1554.4026844207357954_DP ;  casimir_omega_weight(140) = 1762.8312997723442095_DP
         casimir_omega(141) = 8191.7626458267150156_DP ;  casimir_omega_weight(141) = 21023.7348242388870858_DP
return
endsubroutine gauss_legendre_grid140

subroutine gauss_legendre_grid150()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000383003109803_DP ;  casimir_omega_weight(2) = 0.0000982951490805_DP
         casimir_omega(3) = 0.0002018387233057_DP ;  casimir_omega_weight(3) = 0.0002288955502429_DP
         casimir_omega(4) = 0.0004962065489085_DP ;  casimir_omega_weight(4) = 0.0003598887722050_DP
         casimir_omega(5) = 0.0009217258992271_DP ;  casimir_omega_weight(5) = 0.0004912128409855_DP
         casimir_omega(6) = 0.0014787763140779_DP ;  casimir_omega_weight(6) = 0.0006229694071221_DP
         casimir_omega(7) = 0.0021678466878313_DP ;  casimir_omega_weight(7) = 0.0007552719726242_DP
         casimir_omega(8) = 0.0029895409920918_DP ;  casimir_omega_weight(8) = 0.0008882367887690_DP
         casimir_omega(9) = 0.0039445802979730_DP ;  casimir_omega_weight(9) = 0.0010219817550088_DP
         casimir_omega(10) = 0.0050338043616427_DP ;  casimir_omega_weight(10) = 0.0011566263575065_DP
         casimir_omega(11) = 0.0062581732773587_DP ;  casimir_omega_weight(11) = 0.0012922918178020_DP
         casimir_omega(12) = 0.0076187693138236_DP ;  casimir_omega_weight(12) = 0.0014291013002180_DP
         casimir_omega(13) = 0.0091167989708115_DP ;  casimir_omega_weight(13) = 0.0015671801421939_DP
         casimir_omega(14) = 0.0107535952730575_DP ;  casimir_omega_weight(14) = 0.0017066560976771_DP
         casimir_omega(15) = 0.0125306203128102_DP ;  casimir_omega_weight(15) = 0.0018476595908483_DP
         casimir_omega(16) = 0.0144494680510252_DP ;  casimir_omega_weight(16) = 0.0019903239797110_DP
         casimir_omega(17) = 0.0165118673872259_DP ;  casimir_omega_weight(17) = 0.0021347858299232_DP
         casimir_omega(18) = 0.0187196855086144_DP ;  casimir_omega_weight(18) = 0.0022811851996199_DP
         casimir_omega(19) = 0.0210749315298857_DP ;  casimir_omega_weight(19) = 0.0024296659361508_DP
         casimir_omega(20) = 0.0235797604361600_DP ;  casimir_omega_weight(20) = 0.0025803759857755_DP
         casimir_omega(21) = 0.0262364773425698_DP ;  casimir_omega_weight(21) = 0.0027334677174617_DP
         casimir_omega(22) = 0.0290475420852017_DP ;  casimir_omega_weight(22) = 0.0028890982620197_DP
         casimir_omega(23) = 0.0320155741593907_DP ;  casimir_omega_weight(23) = 0.0030474298678939_DP
         casimir_omega(24) = 0.0351433580227134_DP ;  casimir_omega_weight(24) = 0.0032086302750199_DP
         casimir_omega(25) = 0.0384338487815094_DP ;  casimir_omega_weight(25) = 0.0033728731082773_DP
         casimir_omega(26) = 0.0418901782813124_DP ;  casimir_omega_weight(26) = 0.0035403382921439_DP
         casimir_omega(27) = 0.0455156616232749_DP ;  casimir_omega_weight(27) = 0.0037112124883270_DP
         casimir_omega(28) = 0.0493138041304742_DP ;  casimir_omega_weight(28) = 0.0038856895582189_DP
         casimir_omega(29) = 0.0532883087899314_DP ;  casimir_omega_weight(29) = 0.0040639710522210_DP
         casimir_omega(30) = 0.0574430841982865_DP ;  casimir_omega_weight(30) = 0.0042462667281040_DP
         casimir_omega(31) = 0.0617822530413250_DP ;  casimir_omega_weight(31) = 0.0044327951007592_DP
         casimir_omega(32) = 0.0663101611399991_DP ;  casimir_omega_weight(32) = 0.0046237840258785_DP
         casimir_omega(33) = 0.0710313870982296_DP ;  casimir_omega_weight(33) = 0.0048194713203091_DP
         casimir_omega(34) = 0.0759507525906301_DP ;  casimir_omega_weight(34) = 0.0050201054220616_DP
         casimir_omega(35) = 0.0810733333313951_DP ;  casimir_omega_weight(35) = 0.0052259460931840_DP
         casimir_omega(36) = 0.0864044707689398_DP ;  casimir_omega_weight(36) = 0.0054372651690059_DP
         casimir_omega(37) = 0.0919497845545418_DP ;  casimir_omega_weight(37) = 0.0056543473575457_DP
         casimir_omega(38) = 0.0977151858371663_DP ;  casimir_omega_weight(38) = 0.0058774910932008_DP
         casimir_omega(39) = 0.1037068914409889_DP ;  casimir_omega_weight(39) = 0.0061070094492185_DP
         casimir_omega(40) = 0.1099314389867865_DP ;  casimir_omega_weight(40) = 0.0063432311138251_DP
         casimir_omega(41) = 0.1163957030234947_DP ;  casimir_omega_weight(41) = 0.0065865014353443_DP
         casimir_omega(42) = 0.1231069122417734_DP ;  casimir_omega_weight(42) = 0.0068371835421131_DP
         casimir_omega(43) = 0.1300726678475106_DP ;  casimir_omega_weight(43) = 0.0070956595435451_DP
         casimir_omega(44) = 0.1373009631798177_DP ;  casimir_omega_weight(44) = 0.0073623318192600_DP
         casimir_omega(45) = 0.1448002046653307_DP ;  casimir_omega_weight(45) = 0.0076376244038846_DP
         casimir_omega(46) = 0.1525792342085624_DP ;  casimir_omega_weight(46) = 0.0079219844758026_DP
         casimir_omega(47) = 0.1606473531267510_DP ;  casimir_omega_weight(47) = 0.0082158839589835_DP
         casimir_omega(48) = 0.1690143477471766_DP ;  casimir_omega_weight(48) = 0.0085198212478379_DP
         casimir_omega(49) = 0.1776905167953868_DP ;  casimir_omega_weight(49) = 0.0088343230660808_DP
         casimir_omega(50) = 0.1866867007142440_DP ;  casimir_omega_weight(50) = 0.0091599464716312_DP
         casimir_omega(51) = 0.1960143130663664_DP ;  casimir_omega_weight(51) = 0.0094972810207951_DP
         casimir_omega(52) = 0.2056853741864009_DP ;  casimir_omega_weight(52) = 0.0098469511063160_DP
         casimir_omega(53) = 0.2157125472649143_DP ;  casimir_omega_weight(53) = 0.0102096184853528_DP
         casimir_omega(54) = 0.2261091770625358_DP ;  casimir_omega_weight(54) = 0.0105859850151231_DP
         casimir_omega(55) = 0.2368893314716532_DP ;  casimir_omega_weight(55) = 0.0109767956157801_DP
         casimir_omega(56) = 0.2480678461635159_DP ;  casimir_omega_weight(56) = 0.0113828414821532_DP
         casimir_omega(57) = 0.2596603725814025_DP ;  casimir_omega_weight(57) = 0.0118049635683221_DP
         casimir_omega(58) = 0.2716834295656858_DP ;  casimir_omega_weight(58) = 0.0122440563714953_DP
         casimir_omega(59) = 0.2841544589245824_DP ;  casimir_omega_weight(59) = 0.0127010720446541_DP
         casimir_omega(60) = 0.2970918852953650_DP ;  casimir_omega_weight(60) = 0.0131770248705712_DP
         casimir_omega(61) = 0.3105151806752388_DP ;  casimir_omega_weight(61) = 0.0136729961334884_DP
         casimir_omega(62) = 0.3244449340393677_DP ;  casimir_omega_weight(62) = 0.0141901394288327_DP
         casimir_omega(63) = 0.3389029265061507_DP ;  casimir_omega_weight(63) = 0.0147296864559076_DP
         casimir_omega(64) = 0.3539122125573445_DP ;  casimir_omega_weight(64) = 0.0152929533437118_DP
         casimir_omega(65) = 0.3694972078736333_DP ;  casimir_omega_weight(65) = 0.0158813475658772_DP
         casimir_omega(66) = 0.3856837844054627_DP ;  casimir_omega_weight(66) = 0.0164963755072842_DP
         casimir_omega(67) = 0.4024993733651962_DP ;  casimir_omega_weight(67) = 0.0171396507524300_DP
         casimir_omega(68) = 0.4199730769008416_DP ;  casimir_omega_weight(68) = 0.0178129031740626_DP
         casimir_omega(69) = 0.4381357892948283_DP ;  casimir_omega_weight(69) = 0.0185179889101875_DP
         casimir_omega(70) = 0.4570203286247451_DP ;  casimir_omega_weight(70) = 0.0192569013284785_DP
         casimir_omega(71) = 0.4766615799280563_DP ;  casimir_omega_weight(71) = 0.0200317830895239_DP
         casimir_omega(72) = 0.4970966510311717_DP ;  casimir_omega_weight(72) = 0.0208449394344789_DP
         casimir_omega(73) = 0.5183650423367001_DP ;  casimir_omega_weight(73) = 0.0216988528388798_DP
         casimir_omega(74) = 0.5405088320134774_DP ;  casimir_omega_weight(74) = 0.0225961991927514_DP
         casimir_omega(75) = 0.5635728782044771_DP ;  casimir_omega_weight(75) = 0.0235398656883883_DP
         casimir_omega(76) = 0.5876050400608589_DP ;  casimir_omega_weight(76) = 0.0245329706213245_DP
         casimir_omega(77) = 0.6126564196296106_DP ;  casimir_omega_weight(77) = 0.0255788853379872_DP
         casimir_omega(78) = 0.6387816268713057_DP ;  casimir_omega_weight(78) = 0.0266812585954585_DP
         casimir_omega(79) = 0.6660390703680917_DP ;  casimir_omega_weight(79) = 0.0278440436359377_DP
         casimir_omega(80) = 0.6944912766053469_DP ;  casimir_omega_weight(80) = 0.0290715283210738_DP
         casimir_omega(81) = 0.7242052410798175_DP ;  casimir_omega_weight(81) = 0.0303683687209030_DP
         casimir_omega(82) = 0.7552528149097636_DP ;  casimir_omega_weight(82) = 0.0317396266095284_DP
         casimir_omega(83) = 0.7877111311072388_DP ;  casimir_omega_weight(83) = 0.0331908113862727_DP
         casimir_omega(84) = 0.8216630752292881_DP ;  casimir_omega_weight(84) = 0.0347279270189165_DP
         casimir_omega(85) = 0.8571978057655316_DP ;  casimir_omega_weight(85) = 0.0363575246961018_DP
         casimir_omega(86) = 0.8944113303584319_DP ;  casimir_omega_weight(86) = 0.0380867619822385_DP
         casimir_omega(87) = 0.9334071448063219_DP ;  casimir_omega_weight(87) = 0.0399234693925317_DP
         casimir_omega(88) = 0.9742969427880467_DP ;  casimir_omega_weight(88) = 0.0418762254519691_DP
         casimir_omega(89) = 1.0172014053956084_DP ;  casimir_omega_weight(89) = 0.0439544414742446_DP
         casimir_omega(90) = 1.0622510808960703_DP ;  casimir_omega_weight(90) = 0.0461684575000648_DP
         casimir_omega(91) = 1.1095873667003175_DP ;  casimir_omega_weight(91) = 0.0485296510748972_DP
         casimir_omega(92) = 1.1593636073352460_DP ;  casimir_omega_weight(92) = 0.0510505608322617_DP
         casimir_omega(93) = 1.2117463243470714_DP ;  casimir_omega_weight(93) = 0.0537450271887098_DP
         casimir_omega(94) = 1.2669165965667539_DP ;  casimir_omega_weight(94) = 0.0566283528629529_DP
         casimir_omega(95) = 1.3250716121167085_DP ;  casimir_omega_weight(95) = 0.0597174864177814_DP
         casimir_omega(96) = 1.3864264170195661_DP ;  casimir_omega_weight(96) = 0.0630312326072954_DP
         casimir_omega(97) = 1.4512158893930305_DP ;  casimir_omega_weight(97) = 0.0665904940153168_DP
         casimir_omega(98) = 1.5196969731120129_DP ;  casimir_omega_weight(98) = 0.0704185493206377_DP
         casimir_omega(99) = 1.5921512106535727_DP ;  casimir_omega_weight(99) = 0.0745413745552091_DP
         casimir_omega(100) = 1.6688876218122246_DP ;  casimir_omega_weight(100) = 0.0789880149748762_DP
         casimir_omega(101) = 1.7502459833325426_DP ;  casimir_omega_weight(101) = 0.0837910166927215_DP
         casimir_omega(102) = 1.8366005745616667_DP ;  casimir_omega_weight(102) = 0.0889869291007334_DP
         casimir_omega(103) = 1.9283644663635711_DP ;  casimir_omega_weight(103) = 0.0946168914127534_DP
         casimir_omega(104) = 2.0259944452440553_DP ;  casimir_omega_weight(104) = 0.1007273195112619_DP
         casimir_omega(105) = 2.1299966825213676_DP ;  casimir_omega_weight(105) = 0.1073707128149690_DP
         casimir_omega(106) = 2.2409332802138358_DP ;  casimir_omega_weight(106) = 0.1146066052861403_DP
         casimir_omega(107) = 2.3594298520853165_DP ;  casimir_omega_weight(107) = 0.1225026902049779_DP
         casimir_omega(108) = 2.4861843312448979_DP ;  casimir_omega_weight(108) = 0.1311361552613740_DP
         casimir_omega(109) = 2.6219772364489669_DP ;  casimir_omega_weight(109) = 0.1405952732611492_DP
         casimir_omega(110) = 2.7676836798799394_DP ;  casimir_omega_weight(110) = 0.1509813048478161_DP
         casimir_omega(111) = 2.9242874623724262_DP ;  casimir_omega_weight(111) = 0.1624107838142718_DP
         casimir_omega(112) = 3.0928976813459657_DP ;  casimir_omega_weight(112) = 0.1750182737712133_DP
         casimir_omega(113) = 3.2747683767086091_DP ;  casimir_omega_weight(113) = 0.1889597084252256_DP
         casimir_omega(114) = 3.4713218668293306_DP ;  casimir_omega_weight(114) = 0.2044164582260951_DP
         casimir_omega(115) = 3.6841765884773290_DP ;  casimir_omega_weight(115) = 0.2216003060224164_DP
         casimir_omega(116) = 3.9151804622930793_DP ;  casimir_omega_weight(116) = 0.2407595668497630_DP
         casimir_omega(117) = 4.1664510736105305_DP ;  casimir_omega_weight(117) = 0.2621866565387656_DP
         casimir_omega(118) = 4.4404243073177394_DP ;  casimir_omega_weight(118) = 0.2862275067197685_DP
         casimir_omega(119) = 4.7399135323961437_DP ;  casimir_omega_weight(119) = 0.3132933488143245_DP
         casimir_omega(120) = 5.0681820348257496_DP ;  casimir_omega_weight(120) = 0.3438755592533969_DP
         casimir_omega(121) = 5.4290321997550333_DP ;  casimir_omega_weight(121) = 0.3785644904135991_DP
         casimir_omega(122) = 5.8269160200293024_DP ;  casimir_omega_weight(122) = 0.4180735326832860_DP
         casimir_omega(123) = 6.2670729649077650_DP ;  casimir_omega_weight(123) = 0.4632701009163658_DP
         casimir_omega(124) = 6.7557032335021772_DP ;  casimir_omega_weight(124) = 0.5152158700809475_DP
         casimir_omega(125) = 7.3001871655959674_DP ;  casimir_omega_weight(125) = 0.5752194855491725_DP
         casimir_omega(126) = 7.9093654175491368_DP ;  casimir_omega_weight(126) = 0.6449062732582831_DP
         casimir_omega(127) = 8.5938999252385813_DP ;  casimir_omega_weight(127) = 0.7263113749446073_DP
         casimir_omega(128) = 9.3667434153302054_DP ;  casimir_omega_weight(128) = 0.8220055492568099_DP
         casimir_omega(129) = 10.2437564380537562_DP ;  casimir_omega_weight(129) = 0.9352671140826579_DP
         casimir_omega(130) = 11.2445273730755844_DP ;  casimir_omega_weight(130) = 1.0703199760359534_DP
         casimir_omega(131) = 12.3934754597842396_DP ;  casimir_omega_weight(131) = 1.2326677522739946_DP
         casimir_omega(132) = 13.7213542542116862_DP ;  casimir_omega_weight(132) = 1.4295699229746082_DP
         casimir_omega(133) = 15.2673306827973718_DP ;  casimir_omega_weight(133) = 1.6707317094012288_DP
         casimir_omega(134) = 17.0819060308451434_DP ;  casimir_omega_weight(134) = 1.9693219476810964_DP
         casimir_omega(135) = 19.2310923083795089_DP ;  casimir_omega_weight(135) = 2.3435053503548025_DP
         casimir_omega(136) = 21.8025006837508926_DP ;  casimir_omega_weight(136) = 2.8188010734975331_DP
         casimir_omega(137) = 24.9144119858758408_DP ;  casimir_omega_weight(137) = 3.4318046478094240_DP
         casimir_omega(138) = 28.7296231960660116_DP ;  casimir_omega_weight(138) = 4.2362279372084783_DP
         casimir_omega(139) = 33.4771758522438176_DP ;  casimir_omega_weight(139) = 5.3130162378688821_DP
         casimir_omega(140) = 39.4875439452589987_DP ;  casimir_omega_weight(140) = 6.7879191954490627_DP
         casimir_omega(141) = 47.2517259902873903_DP ;  casimir_omega_weight(141) = 8.8633085303890748_DP
         casimir_omega(142) = 57.5247734514540099_DP ;  casimir_omega_weight(142) = 11.8786730180775493_DP
         casimir_omega(143) = 71.5164861676359891_DP ;  casimir_omega_weight(143) = 16.4324727293867312_DP
         casimir_omega(144) = 91.2644623269567035_DP ;  casimir_omega_weight(144) = 23.6452571206989646_DP
         casimir_omega(145) = 120.4198239637242125_DP ;  casimir_omega_weight(145) = 35.7785084816333665_DP
         casimir_omega(146) = 166.0634038471445422_DP ;  casimir_omega_weight(146) = 57.8560445756170978_DP
         casimir_omega(147) = 243.4445267839735152_DP ;  casimir_omega_weight(147) = 102.5567498437746394_DP
         casimir_omega(148) = 390.5716442403618203_DP ;  casimir_omega_weight(148) = 208.1462690116796637_DP
         casimir_omega(149) = 725.5043303880748908_DP ;  casimir_omega_weight(149) = 526.1939070885597403_DP
         casimir_omega(150) = 1783.6022449228280493_DP ;  casimir_omega_weight(150) = 2022.6971840680328114_DP
         casimir_omega(151) = 9399.4014875003595080_DP ;  casimir_omega_weight(151) = 24122.9260763834754471_DP
return
endsubroutine gauss_legendre_grid150

subroutine gauss_legendre_grid200()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000215793905679_DP ;  casimir_omega_weight(2) = 0.0000553810126550_DP
         casimir_omega(3) = 0.0001137121623953_DP ;  casimir_omega_weight(3) = 0.0001289428498122_DP
         casimir_omega(4) = 0.0002795135905008_DP ;  casimir_omega_weight(4) = 0.0002026768949584_DP
         casimir_omega(5) = 0.0005191019384759_DP ;  casimir_omega_weight(5) = 0.0002765199911217_DP
         casimir_omega(6) = 0.0008325994837745_DP ;  casimir_omega_weight(6) = 0.0003505009497414_DP
         casimir_omega(7) = 0.0012201615094267_DP ;  casimir_omega_weight(7) = 0.0004246549527141_DP
         casimir_omega(8) = 0.0016819791793673_DP ;  casimir_omega_weight(8) = 0.0004990183397904_DP
         casimir_omega(9) = 0.0022182802440426_DP ;  casimir_omega_weight(9) = 0.0005736279026367_DP
         casimir_omega(10) = 0.0028293294135225_DP ;  casimir_omega_weight(10) = 0.0006485207628374_DP
         casimir_omega(11) = 0.0035154286802077_DP ;  casimir_omega_weight(11) = 0.0007237343661051_DP
         casimir_omega(12) = 0.0042769176540341_DP ;  casimir_omega_weight(12) = 0.0007993065070516_DP
         casimir_omega(13) = 0.0051141739282755_DP ;  casimir_omega_weight(13) = 0.0008752753639685_DP
         casimir_omega(14) = 0.0060276134824193_DP ;  casimir_omega_weight(14) = 0.0009516795376748_DP
         casimir_omega(15) = 0.0070176911250413_DP ;  casimir_omega_weight(15) = 0.0010285580924878_DP
         casimir_omega(16) = 0.0080849009784018_DP ;  casimir_omega_weight(16) = 0.0011059505986329_DP
         casimir_omega(17) = 0.0092297770060647_DP ;  casimir_omega_weight(17) = 0.0011838971758566_DP
         casimir_omega(18) = 0.0104528935847092_DP ;  casimir_omega_weight(18) = 0.0012624385382000_DP
         casimir_omega(19) = 0.0117548661212935_DP ;  casimir_omega_weight(19) = 0.0013416160399381_DP
         casimir_omega(20) = 0.0131363517167697_DP ;  casimir_omega_weight(20) = 0.0014214717227531_DP
         casimir_omega(21) = 0.0145980498776201_DP ;  casimir_omega_weight(21) = 0.0015020483642166_DP
         casimir_omega(22) = 0.0161407032765679_DP ;  casimir_omega_weight(22) = 0.0015833895276731_DP
         casimir_omega(23) = 0.0177650985639044_DP ;  casimir_omega_weight(23) = 0.0016655396136077_DP
         casimir_omega(24) = 0.0194720672309808_DP ;  casimir_omega_weight(24) = 0.0017485439126212_DP
         casimir_omega(25) = 0.0212624865275136_DP ;  casimir_omega_weight(25) = 0.0018324486601046_DP
         casimir_omega(26) = 0.0231372804344687_DP ;  casimir_omega_weight(26) = 0.0019173010927399_DP
         casimir_omega(27) = 0.0250974206944032_DP ;  casimir_omega_weight(27) = 0.0020031495069369_DP
         casimir_omega(28) = 0.0271439279012723_DP ;  casimir_omega_weight(28) = 0.0020900433193455_DP
         casimir_omega(29) = 0.0292778726518345_DP ;  casimir_omega_weight(29) = 0.0021780331295655_DP
         casimir_omega(30) = 0.0315003767609272_DP ;  casimir_omega_weight(30) = 0.0022671707851983_DP
         casimir_omega(31) = 0.0338126145430259_DP ;  casimir_omega_weight(31) = 0.0023575094493940_DP
         casimir_omega(32) = 0.0362158141626589_DP ;  casimir_omega_weight(32) = 0.0024491036710337_DP
         casimir_omega(33) = 0.0387112590563966_DP ;  casimir_omega_weight(33) = 0.0025420094577281_DP
         casimir_omega(34) = 0.0413002894293119_DP ;  casimir_omega_weight(34) = 0.0026362843517937_DP
         casimir_omega(35) = 0.0439843038289836_DP ;  casimir_omega_weight(35) = 0.0027319875093885_DP
         casimir_omega(36) = 0.0467647608002974_DP ;  casimir_omega_weight(36) = 0.0028291797830068_DP
         casimir_omega(37) = 0.0496431806245008_DP ;  casimir_omega_weight(37) = 0.0029279238075271_DP
         casimir_omega(38) = 0.0526211471461753_DP ;  casimir_omega_weight(38) = 0.0030282840900357_DP
         casimir_omega(39) = 0.0557003096920093_DP ;  casimir_omega_weight(39) = 0.0031303271036421_DP
         casimir_omega(40) = 0.0588823850854918_DP ;  casimir_omega_weight(40) = 0.0032341213855474_DP
         casimir_omega(41) = 0.0621691597618909_DP ;  casimir_omega_weight(41) = 0.0033397376395911_DP
         casimir_omega(42) = 0.0655624919881460_DP ;  casimir_omega_weight(42) = 0.0034472488435742_DP
         casimir_omega(43) = 0.0690643141925810_DP ;  casimir_omega_weight(43) = 0.0035567303616293_DP
         casimir_omega(44) = 0.0726766354096405_DP ;  casimir_omega_weight(44) = 0.0036682600619494_DP
         casimir_omega(45) = 0.0764015438451641_DP ;  casimir_omega_weight(45) = 0.0037819184401926_DP
         casimir_omega(46) = 0.0802412095680528_DP ;  casimir_omega_weight(46) = 0.0038977887489136_DP
         casimir_omega(47) = 0.0841978873345242_DP ;  casimir_omega_weight(47) = 0.0040159571333816_DP
         casimir_omega(48) = 0.0882739195515462_DP ;  casimir_omega_weight(48) = 0.0041365127741797_DP
         casimir_omega(49) = 0.0924717393864264_DP ;  casimir_omega_weight(49) = 0.0042595480369912_DP
         casimir_omega(50) = 0.0967938740299694_DP ;  casimir_omega_weight(50) = 0.0043851586300252_DP
         casimir_omega(51) = 0.1012429481210691_DP ;  casimir_omega_weight(51) = 0.0045134437695390_DP
         casimir_omega(52) = 0.1058216873410832_DP ;  casimir_omega_weight(52) = 0.0046445063539698_DP
         casimir_omega(53) = 0.1105329221868629_DP ;  casimir_omega_weight(53) = 0.0047784531472041_DP
         casimir_omega(54) = 0.1153795919318593_DP ;  casimir_omega_weight(54) = 0.0049153949715562_DP
         casimir_omega(55) = 0.1203647487853102_DP ;  casimir_omega_weight(55) = 0.0050554469110666_DP
         casimir_omega(56) = 0.1254915622601624_DP ;  casimir_omega_weight(56) = 0.0051987285257707_DP
         casimir_omega(57) = 0.1307633237610310_DP ;  casimir_omega_weight(57) = 0.0053453640776362_DP
         casimir_omega(58) = 0.1361834514042352_DP ;  casimir_omega_weight(58) = 0.0054954827689074_DP
         casimir_omega(59) = 0.1417554950827202_DP ;  casimir_omega_weight(59) = 0.0056492189936595_DP
         casimir_omega(60) = 0.1474831417894883_DP ;  casimir_omega_weight(60) = 0.0058067126034183_DP
         casimir_omega(61) = 0.1533702212140504_DP ;  casimir_omega_weight(61) = 0.0059681091877433_DP
         casimir_omega(62) = 0.1594207116273647_DP ;  casimir_omega_weight(62) = 0.0061335603707724_DP
         casimir_omega(63) = 0.1656387460717258_DP ;  casimir_omega_weight(63) = 0.0063032241247748_DP
         casimir_omega(64) = 0.1720286188731683_DP ;  casimir_omega_weight(64) = 0.0064772651018166_DP
         casimir_omega(65) = 0.1785947924951124_DP ;  casimir_omega_weight(65) = 0.0066558549847812_DP
         casimir_omega(66) = 0.1853419047532305_DP ;  casimir_omega_weight(66) = 0.0068391728590238_DP
         casimir_omega(67) = 0.1922747764128607_DP ;  casimir_omega_weight(67) = 0.0070274056060493_DP
         casimir_omega(68) = 0.1993984191917443_DP ;  casimir_omega_weight(68) = 0.0072207483207380_DP
         casimir_omega(69) = 0.2067180441924122_DP ;  casimir_omega_weight(69) = 0.0074194047537043_DP
         casimir_omega(70) = 0.2142390707902285_DP ;  casimir_omega_weight(70) = 0.0076235877805456_DP
         casimir_omega(71) = 0.2219671360049028_DP ;  casimir_omega_weight(71) = 0.0078335198998465_DP
         casimir_omega(72) = 0.2299081043851975_DP ;  casimir_omega_weight(72) = 0.0080494337619490_DP
         casimir_omega(73) = 0.2380680784387168_DP ;  casimir_omega_weight(73) = 0.0082715727306666_DP
         casimir_omega(74) = 0.2464534096408280_DP ;  casimir_omega_weight(74) = 0.0085001914802903_DP
         casimir_omega(75) = 0.2550707100592874_DP ;  casimir_omega_weight(75) = 0.0087355566304025_DP
         casimir_omega(76) = 0.2639268646337241_DP ;  casimir_omega_weight(76) = 0.0089779474212546_DP
         casimir_omega(77) = 0.2730290441519851_DP ;  casimir_omega_weight(77) = 0.0092276564326411_DP
         casimir_omega(78) = 0.2823847189684374_DP ;  casimir_omega_weight(78) = 0.0094849903494711_DP
         casimir_omega(79) = 0.2920016735126231_DP ;  casimir_omega_weight(79) = 0.0097502707775082_DP
         casimir_omega(80) = 0.3018880216402651_DP ;  casimir_omega_weight(80) = 0.0100238351129975_DP
         casimir_omega(81) = 0.3120522228825371_DP ;  casimir_omega_weight(81) = 0.0103060374702444_DP
         casimir_omega(82) = 0.3225030996537027_DP ;  casimir_omega_weight(82) = 0.0105972496715528_DP
         casimir_omega(83) = 0.3332498554818349_DP ;  casimir_omega_weight(83) = 0.0108978623042763_DP
         casimir_omega(84) = 0.3443020943322911_DP ;  casimir_omega_weight(84) = 0.0112082858501653_DP
         casimir_omega(85) = 0.3556698410989996_DP ;  casimir_omega_weight(85) = 0.0115289518926319_DP
         casimir_omega(86) = 0.3673635633445041_DP ;  casimir_omega_weight(86) = 0.0118603144080590_DP
         casimir_omega(87) = 0.3793941943760816_DP ;  casimir_omega_weight(87) = 0.0122028511477895_DP
         casimir_omega(88) = 0.3917731577521802_DP ;  casimir_omega_weight(88) = 0.0125570651180438_DP
         casimir_omega(89) = 0.4045123933210029_DP ;  casimir_omega_weight(89) = 0.0129234861656630_DP
         casimir_omega(90) = 0.4176243849012919_DP ;  casimir_omega_weight(90) = 0.0133026726782628_DP
         casimir_omega(91) = 0.4311221897243572_DP ;  casimir_omega_weight(91) = 0.0136952134081969_DP
         casimir_omega(92) = 0.4450194697661978_DP ;  casimir_omega_weight(92) = 0.0141017294305775_DP
         casimir_omega(93) = 0.4593305251092991_DP ;  casimir_omega_weight(93) = 0.0145228762465379_DP
         casimir_omega(94) = 0.4740703294853710_DP ;  casimir_omega_weight(94) = 0.0149593460440043_DP
         casimir_omega(95) = 0.4892545681631422_DP ;  casimir_omega_weight(95) = 0.0154118701293651_DP
         casimir_omega(96) = 0.5048996783593380_DP ;  casimir_omega_weight(96) = 0.0158812215447231_DP
         casimir_omega(97) = 0.5210228923663530_DP ;  casimir_omega_weight(97) = 0.0163682178868584_DP
         casimir_omega(98) = 0.5376422836070183_DP ;  casimir_omega_weight(98) = 0.0168737243455188_DP
         casimir_omega(99) = 0.5547768158453438_DP ;  casimir_omega_weight(99) = 0.0173986569804841_DP
         casimir_omega(100) = 0.5724463958025090_DP ;  casimir_omega_weight(100) = 0.0179439862587155_DP
         casimir_omega(101) = 0.5906719294497187_DP ;  casimir_omega_weight(101) = 0.0185107408750380_DP
         casimir_omega(102) = 0.6094753822741883_DP ;  casimir_omega_weight(102) = 0.0191000118822349_DP
         casimir_omega(103) = 0.6288798438416542_DP ;  casimir_omega_weight(103) = 0.0197129571589987_DP
         casimir_omega(104) = 0.6489095970087505_DP ;  casimir_omega_weight(104) = 0.0203508062471879_DP
         casimir_omega(105) = 0.6695901921715968_DP ;  casimir_omega_weight(105) = 0.0210148655930957_DP
         casimir_omega(106) = 0.6909485269734906_DP ;  casimir_omega_weight(106) = 0.0217065242310988_DP
         casimir_omega(107) = 0.7130129319349399_DP ;  casimir_omega_weight(107) = 0.0224272599521303_DP
         casimir_omega(108) = 0.7358132625139983_DP ;  casimir_omega_weight(108) = 0.0231786460040752_DP
         casimir_omega(109) = 0.7593809981544282_DP ;  casimir_omega_weight(109) = 0.0239623583761638_DP
         casimir_omega(110) = 0.7837493489341616_DP ;  casimir_omega_weight(110) = 0.0247801837253617_DP
         casimir_omega(111) = 0.8089533704876665_DP ;  casimir_omega_weight(111) = 0.0256340280090749_DP
         casimir_omega(112) = 0.8350300879436762_DP ;  casimir_omega_weight(112) = 0.0265259258957789_DP
         casimir_omega(113) = 0.8620186296954081_DP ;  casimir_omega_weight(113) = 0.0274580510333772_DP
         casimir_omega(114) = 0.8899603719046489_DP ;  casimir_omega_weight(114) = 0.0284327272642312_DP
         casimir_omega(115) = 0.9188990947351244_DP ;  casimir_omega_weight(115) = 0.0294524408862141_DP
         casimir_omega(116) = 0.9488811514156762_DP ;  casimir_omega_weight(116) = 0.0305198540708038_DP
         casimir_omega(117) = 0.9799556513513037_DP ;  casimir_omega_weight(117) = 0.0316378195626372_DP
         casimir_omega(118) = 1.0121746586317817_DP ;  casimir_omega_weight(118) = 0.0328093967997079_DP
         casimir_omega(119) = 1.0455934074352691_DP ;  casimir_omega_weight(119) = 0.0340378696107215_DP
         casimir_omega(120) = 1.0802705359901448_DP ;  casimir_omega_weight(120) = 0.0353267656652584_DP
         casimir_omega(121) = 1.1162683409448202_DP ;  casimir_omega_weight(121) = 0.0366798778744871_DP
         casimir_omega(122) = 1.1536530542053263_DP ;  casimir_omega_weight(122) = 0.0381012879654349_DP
         casimir_omega(123) = 1.1924951445373411_DP ;  casimir_omega_weight(123) = 0.0395953924801177_DP
         casimir_omega(124) = 1.2328696464968631_DP ;  casimir_omega_weight(124) = 0.0411669314840302_DP
         casimir_omega(125) = 1.2748565195563479_DP ;  casimir_omega_weight(125) = 0.0428210203056486_DP
         casimir_omega(126) = 1.3185410406359617_DP ;  casimir_omega_weight(126) = 0.0445631846718528_DP
         casimir_omega(127) = 1.3640142336385699_DP ;  casimir_omega_weight(127) = 0.0463993996535593_DP
         casimir_omega(128) = 1.4113733400292112_DP ;  casimir_omega_weight(128) = 0.0483361328927191_DP
         casimir_omega(129) = 1.4607223350029952_DP ;  casimir_omega_weight(129) = 0.0503803926476705_DP
         casimir_omega(130) = 1.5121724943593002_DP ;  casimir_omega_weight(130) = 0.0525397812694427_DP
         casimir_omega(131) = 1.5658430178556948_DP ;  casimir_omega_weight(131) = 0.0548225548096488_DP
         casimir_omega(132) = 1.6218617155651756_DP ;  casimir_omega_weight(132) = 0.0572376895622915_DP
         casimir_omega(133) = 1.6803657646204571_DP ;  casimir_omega_weight(133) = 0.0597949564603525_DP
         casimir_omega(134) = 1.7415025447169661_DP ;  casimir_omega_weight(134) = 0.0625050043857536_DP
         casimir_omega(135) = 1.8054305618833391_DP ;  casimir_omega_weight(135) = 0.0653794536123769_DP
         casimir_omega(136) = 1.8723204713393760_DP ;  casimir_omega_weight(136) = 0.0684310007900314_DP
         casimir_omega(137) = 1.9423562117769007_DP ;  casimir_omega_weight(137) = 0.0716735370979783_DP
         casimir_omega(138) = 2.0157362651537118_DP ;  casimir_omega_weight(138) = 0.0751222814561906_DP
         casimir_omega(139) = 2.0926750581274964_DP ;  casimir_omega_weight(139) = 0.0787939309879833_DP
         casimir_omega(140) = 2.1734045236258357_DP ;  casimir_omega_weight(140) = 0.0827068312886235_DP
         casimir_omega(141) = 2.2581758438105313_DP ;  casimir_omega_weight(141) = 0.0868811694819625_DP
         casimir_omega(142) = 2.3472613989228557_DP ;  casimir_omega_weight(142) = 0.0913391935543566_DP
         casimir_omega(143) = 2.4409569502787667_DP ;  casimir_omega_weight(143) = 0.0961054620589529_DP
         casimir_omega(144) = 2.5395840901259215_DP ;  casimir_omega_weight(144) = 0.1012071290045160_DP
         casimir_omega(145) = 2.6434929963069287_DP ;  casimir_omega_weight(145) = 0.1066742696057121_DP
         casimir_omega(146) = 2.7530655358523801_DP ;  casimir_omega_weight(146) = 0.1125402536082497_DP
         casimir_omega(147) = 2.8687187689453331_DP ;  casimir_omega_weight(147) = 0.1188421741504127_DP
         casimir_omega(148) = 2.9909089133905615_DP ;  casimir_omega_weight(148) = 0.1256213416309323_DP
         casimir_omega(149) = 3.1201358400765451_DP ;  casimir_omega_weight(149) = 0.1329238538817419_DP
         casimir_omega(150) = 3.2569481822926587_DP ;  casimir_omega_weight(150) = 0.1408012561691468_DP
         casimir_omega(151) = 3.4019491566001281_DP ;  casimir_omega_weight(151) = 0.1493113072633462_DP
         casimir_omega(152) = 3.5558032108024213_DP ;  casimir_omega_weight(152) = 0.1585188711445983_DP
         casimir_omega(153) = 3.7192436361059031_DP ;  casimir_omega_weight(153) = 0.1684969580098193_DP
         casimir_omega(154) = 3.8930813066639924_DP ;  casimir_omega_weight(154) = 0.1793279433011548_DP
         casimir_omega(155) = 4.0782147414422241_DP ;  casimir_omega_weight(155) = 0.1911049997499448_DP
         casimir_omega(156) = 4.2756417220980101_DP ;  casimir_omega_weight(156) = 0.2039337852436175_DP
         casimir_omega(157) = 4.4864727480794313_DP ;  casimir_omega_weight(157) = 0.2179344390981553_DP
         casimir_omega(158) = 4.7119466686377134_DP ;  casimir_omega_weight(158) = 0.2332439516070209_DP
         casimir_omega(159) = 4.9534489037758469_DP ;  casimir_omega_weight(159) = 0.2500189872606435_DP
         casimir_omega(160) = 5.2125327560071852_DP ;  casimir_omega_weight(160) = 0.2684392617377285_DP
         casimir_omega(161) = 5.4909444269612218_DP ;  casimir_omega_weight(161) = 0.2887115979269866_DP
         casimir_omega(162) = 5.7906524935966344_DP ;  casimir_omega_weight(162) = 0.3110748185229690_DP
         casimir_omega(163) = 6.1138827762719510_DP ;  casimir_omega_weight(163) = 0.3358056744264720_DP
         casimir_omega(164) = 6.4631597560335754_DP ;  casimir_omega_weight(164) = 0.3632260623208857_DP
         casimir_omega(165) = 6.8413559856451736_DP ;  casimir_omega_weight(165) = 0.3937118555786046_DP
         casimir_omega(166) = 7.2517513074560611_DP ;  casimir_omega_weight(166) = 0.4277037657995447_DP
         casimir_omega(167) = 7.6981041673095065_DP ;  casimir_omega_weight(167) = 0.4657207757489939_DP
         casimir_omega(168) = 8.1847379328709309_DP ;  casimir_omega_weight(168) = 0.5083768493224324_DP
         casimir_omega(169) = 8.7166459357667350_DP ;  casimir_omega_weight(169) = 0.5564018460431684_DP
         casimir_omega(170) = 9.2996200272260179_DP ;  casimir_omega_weight(170) = 0.6106678686954429_DP
         casimir_omega(171) = 9.9404088607010817_DP ;  casimir_omega_weight(171) = 0.6722226849013757_DP
         casimir_omega(172) = 10.6469140249981908_DP ;  casimir_omega_weight(172) = 0.7423324330297056_DP
         casimir_omega(173) = 11.4284347368994670_DP ;  casimir_omega_weight(173) = 0.8225366176632308_DP
         casimir_omega(174) = 12.2959753354020744_DP ;  casimir_omega_weight(174) = 0.9147195207553257_DP
         casimir_omega(175) = 13.2626346971372104_DP ;  casimir_omega_weight(175) = 1.0212037530638296_DP
         casimir_omega(176) = 14.3441034990614060_DP ;  casimir_omega_weight(176) = 1.1448739773488681_DP
         casimir_omega(177) = 15.5593048638375624_DP ;  casimir_omega_weight(177) = 1.2893422069288696_DP
         casimir_omega(178) = 16.9312276592939099_DP ;  casimir_omega_weight(178) = 1.4591710803930356_DP
         casimir_omega(179) = 18.4880216224413729_DP ;  casimir_omega_weight(179) = 1.6601790288038092_DP
         casimir_omega(180) = 20.2644527248196837_DP ;  casimir_omega_weight(180) = 1.8998627359064393_DP
         casimir_omega(181) = 22.3038608560898091_DP ;  casimir_omega_weight(181) = 2.1879901452298687_DP
         casimir_omega(182) = 24.6608281940391940_DP ;  casimir_omega_weight(182) = 2.5374455464681964_DP
         casimir_omega(183) = 27.4048691571218122_DP ;  casimir_omega_weight(183) = 2.9654539869589254_DP
         casimir_omega(184) = 30.6256146420824535_DP ;  casimir_omega_weight(184) = 3.4953878175055619_DP
         casimir_omega(185) = 34.4402243343051708_DP ;  casimir_omega_weight(185) = 4.1594861854798966_DP
         casimir_omega(186) = 39.0041926000446324_DP ;  casimir_omega_weight(186) = 5.0030410740607607_DP
         casimir_omega(187) = 44.5274470227564194_DP ;  casimir_omega_weight(187) = 6.0910030712743808_DP
         casimir_omega(188) = 51.2989234757585990_DP ;  casimir_omega_weight(188) = 7.5187012276230529_DP
         casimir_omega(189) = 59.7251301945632207_DP ;  casimir_omega_weight(189) = 9.4297991165008614_DP
         casimir_omega(190) = 70.3926000657926778_DP ;  casimir_omega_weight(190) = 12.0474801028245508_DP
         casimir_omega(191) = 84.1727685966634169_DP ;  casimir_omega_weight(191) = 15.7309181747759066_DP
         casimir_omega(192) = 102.4057185477379761_DP ;  casimir_omega_weight(192) = 21.0826458280836881_DP
         casimir_omega(193) = 127.2386305671648188_DP ;  casimir_omega_weight(193) = 29.1648237788774765_DP
         casimir_omega(194) = 162.2878808783572140_DP ;  casimir_omega_weight(194) = 41.9662290108009088_DP
         casimir_omega(195) = 214.0335649906427591_DP ;  casimir_omega_weight(195) = 63.5005923802508647_DP
         casimir_omega(196) = 295.0429080238023403_DP ;  casimir_omega_weight(196) = 102.6843013712320385_DP
         casimir_omega(197) = 432.3807629184830148_DP ;  casimir_omega_weight(197) = 182.0201321357579047_DP
         casimir_omega(198) = 693.5054048478134519_DP ;  casimir_omega_weight(198) = 369.4228323518088928_DP
         casimir_omega(199) = 1287.9516854795951986_DP ;  casimir_omega_weight(199) = 933.9010958349734892_DP
         casimir_omega(200) = 3165.8882604691789311_DP ;  casimir_omega_weight(200) = 3589.9295721146781943_DP
         casimir_omega(201) = 16682.5841938950543408_DP ;  casimir_omega_weight(201) = 42813.9248628811910748_DP
return
endsubroutine gauss_legendre_grid200

subroutine gauss_legendre_grid250()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000138244878274_DP ;  casimir_omega_weight(2) = 0.0000354786450616_DP
         casimir_omega(3) = 0.0000728451794833_DP ;  casimir_omega_weight(3) = 0.0000825983642545_DP
         casimir_omega(4) = 0.0001790474261237_DP ;  casimir_omega_weight(4) = 0.0001298138216163_DP
         casimir_omega(5) = 0.0003324882275691_DP ;  casimir_omega_weight(5) = 0.0001770761921467_DP
         casimir_omega(6) = 0.0005332187726048_DP ;  casimir_omega_weight(6) = 0.0002243955303249_DP
         casimir_omega(7) = 0.0007813029657696_DP ;  casimir_omega_weight(7) = 0.0002717859180660_DP
         casimir_omega(8) = 0.0010768192084039_DP ;  casimir_omega_weight(8) = 0.0003192621107607_DP
         casimir_omega(9) = 0.0014198607750529_DP ;  casimir_omega_weight(9) = 0.0003668390703113_DP
         casimir_omega(10) = 0.0018105359620699_DP ;  casimir_omega_weight(10) = 0.0004145318720880_DP
         casimir_omega(11) = 0.0022489681889255_DP ;  casimir_omega_weight(11) = 0.0004623556861064_DP
         casimir_omega(12) = 0.0027352960922916_DP ;  casimir_omega_weight(12) = 0.0005103257775354_DP
         casimir_omega(13) = 0.0032696736242189_DP ;  casimir_omega_weight(13) = 0.0005584575133180_DP
         casimir_omega(14) = 0.0038522701582477_DP ;  casimir_omega_weight(14) = 0.0006067663710788_DP
         casimir_omega(15) = 0.0044832706049828_DP ;  casimir_omega_weight(15) = 0.0006552679490240_DP
         casimir_omega(16) = 0.0051628755378508_DP ;  casimir_omega_weight(16) = 0.0007039779763584_DP
         casimir_omega(17) = 0.0058913013294513_DP ;  casimir_omega_weight(17) = 0.0007529123240380_DP
         casimir_omega(18) = 0.0066687802987899_DP ;  casimir_omega_weight(18) = 0.0008020870157729_DP
         casimir_omega(19) = 0.0074955608696308_DP ;  casimir_omega_weight(19) = 0.0008515182392577_DP
         casimir_omega(20) = 0.0083719077401862_DP ;  casimir_omega_weight(20) = 0.0009012223576212_DP
         casimir_omega(21) = 0.0092981020643671_DP ;  casimir_omega_weight(21) = 0.0009512159210974_DP
         casimir_omega(22) = 0.0102744416448185_DP ;  casimir_omega_weight(22) = 0.0010015156789221_DP
         casimir_omega(23) = 0.0113012411379720_DP ;  casimir_omega_weight(23) = 0.0010521385914789_DP
         casimir_omega(24) = 0.0123788322713717_DP ;  casimir_omega_weight(24) = 0.0011031018426946_DP
         casimir_omega(25) = 0.0135075640735259_DP ;  casimir_omega_weight(25) = 0.0011544228527071_DP
         casimir_omega(26) = 0.0146878031165697_DP ;  casimir_omega_weight(26) = 0.0012061192908165_DP
         casimir_omega(27) = 0.0159199337720299_DP ;  casimir_omega_weight(27) = 0.0012582090887385_DP
         casimir_omega(28) = 0.0172043584799996_DP ;  casimir_omega_weight(28) = 0.0013107104541810_DP
         casimir_omega(29) = 0.0185414980320547_DP ;  casimir_omega_weight(29) = 0.0013636418847494_DP
         casimir_omega(30) = 0.0199317918682553_DP ;  casimir_omega_weight(30) = 0.0014170221822143_DP
         casimir_omega(31) = 0.0213756983885976_DP ;  casimir_omega_weight(31) = 0.0014708704671498_DP
         casimir_omega(32) = 0.0228736952793025_DP ;  casimir_omega_weight(32) = 0.0015252061939692_DP
         casimir_omega(33) = 0.0244262798543401_DP ;  casimir_omega_weight(33) = 0.0015800491663712_DP
         casimir_omega(34) = 0.0260339694126260_DP ;  casimir_omega_weight(34) = 0.0016354195532309_DP
         casimir_omega(35) = 0.0276973016113262_DP ;  casimir_omega_weight(35) = 0.0016913379049395_DP
         casimir_omega(36) = 0.0294168348557502_DP ;  casimir_omega_weight(36) = 0.0017478251702376_DP
         casimir_omega(37) = 0.0311931487063217_DP ;  casimir_omega_weight(37) = 0.0018049027135458_DP
         casimir_omega(38) = 0.0330268443031479_DP ;  casimir_omega_weight(38) = 0.0018625923328356_DP
         casimir_omega(39) = 0.0349185448087356_DP ;  casimir_omega_weight(39) = 0.0019209162780513_DP
         casimir_omega(40) = 0.0368688958694195_DP ;  casimir_omega_weight(40) = 0.0019798972701252_DP
         casimir_omega(41) = 0.0388785660961110_DP ;  casimir_omega_weight(41) = 0.0020395585206010_DP
         casimir_omega(42) = 0.0409482475649899_DP ;  casimir_omega_weight(42) = 0.0020999237519064_DP
         casimir_omega(43) = 0.0430786563388039_DP ;  casimir_omega_weight(43) = 0.0021610172182985_DP
         casimir_omega(44) = 0.0452705330094654_DP ;  casimir_omega_weight(44) = 0.0022228637275191_DP
         casimir_omega(45) = 0.0475246432626723_DP ;  casimir_omega_weight(45) = 0.0022854886631945_DP
         casimir_omega(46) = 0.0498417784653118_DP ;  casimir_omega_weight(46) = 0.0023489180080067_DP
         casimir_omega(47) = 0.0522227562764442_DP ;  casimir_omega_weight(47) = 0.0024131783676882_DP
         casimir_omega(48) = 0.0546684212827044_DP ;  casimir_omega_weight(48) = 0.0024782969958688_DP
         casimir_omega(49) = 0.0571796456589899_DP ;  casimir_omega_weight(49) = 0.0025443018198130_DP
         casimir_omega(50) = 0.0597573298553582_DP ;  casimir_omega_weight(50) = 0.0026112214671026_DP
         casimir_omega(51) = 0.0624024033110918_DP ;  casimir_omega_weight(51) = 0.0026790852933002_DP
         casimir_omega(52) = 0.0651158251969354_DP ;  casimir_omega_weight(52) = 0.0027479234106395_DP
         casimir_omega(53) = 0.0678985851865647_DP ;  casimir_omega_weight(53) = 0.0028177667178014_DP
         casimir_omega(54) = 0.0707517042583861_DP ;  casimir_omega_weight(54) = 0.0028886469308165_DP
         casimir_omega(55) = 0.0736762355288320_DP ;  casimir_omega_weight(55) = 0.0029605966151555_DP
         casimir_omega(56) = 0.0766732651183620_DP ;  casimir_omega_weight(56) = 0.0030336492190634_DP
         casimir_omega(57) = 0.0797439130514418_DP ;  casimir_omega_weight(57) = 0.0031078391081946_DP
         casimir_omega(58) = 0.0828893341918392_DP ;  casimir_omega_weight(58) = 0.0031832016016131_DP
         casimir_omega(59) = 0.0861107192146282_DP ;  casimir_omega_weight(59) = 0.0032597730092310_DP
         casimir_omega(60) = 0.0894092956163741_DP ;  casimir_omega_weight(60) = 0.0033375906707379_DP
         casimir_omega(61) = 0.0927863287650352_DP ;  casimir_omega_weight(61) = 0.0034166929961150_DP
         casimir_omega(62) = 0.0962431229911939_DP ;  casimir_omega_weight(62) = 0.0034971195077982_DP
         casimir_omega(63) = 0.0997810227223074_DP ;  casimir_omega_weight(63) = 0.0035789108845690_DP
         casimir_omega(64) = 0.1034014136617591_DP ;  casimir_omega_weight(64) = 0.0036621090072702_DP
         casimir_omega(65) = 0.1071057240145657_DP ;  casimir_omega_weight(65) = 0.0037467570064274_DP
         casimir_omega(66) = 0.1108954257617020_DP ;  casimir_omega_weight(66) = 0.0038328993118735_DP
         casimir_omega(67) = 0.1147720359850892_DP ;  casimir_omega_weight(67) = 0.0039205817044726_DP
         casimir_omega(68) = 0.1187371182454096_DP ;  casimir_omega_weight(68) = 0.0040098513700573_DP
         casimir_omega(69) = 0.1227922840149999_DP ;  casimir_omega_weight(69) = 0.0041007569556826_DP
         casimir_omega(70) = 0.1269391941682102_DP ;  casimir_omega_weight(70) = 0.0041933486283094_DP
         casimir_omega(71) = 0.1311795605317177_DP ;  casimir_omega_weight(71) = 0.0042876781360585_DP
         casimir_omega(72) = 0.1355151474974231_DP ;  casimir_omega_weight(72) = 0.0043837988721446_DP
         casimir_omega(73) = 0.1399477737006857_DP ;  casimir_omega_weight(73) = 0.0044817659416413_DP
         casimir_omega(74) = 0.1444793137667983_DP ;  casimir_omega_weight(74) = 0.0045816362312186_DP
         casimir_omega(75) = 0.1491117001287453_DP ;  casimir_omega_weight(75) = 0.0046834684820042_DP
         casimir_omega(76) = 0.1538469249194580_DP ;  casimir_omega_weight(76) = 0.0047873233657368_DP
         casimir_omega(77) = 0.1586870419419320_DP ;  casimir_omega_weight(77) = 0.0048932635643750_DP
         casimir_omega(78) = 0.1636341687207646_DP ;  casimir_omega_weight(78) = 0.0050013538533520_DP
         casimir_omega(79) = 0.1686904886388435_DP ;  casimir_omega_weight(79) = 0.0051116611886579_DP
         casimir_omega(80) = 0.1738582531631218_DP ;  casimir_omega_weight(80) = 0.0052242547979698_DP
         casimir_omega(81) = 0.1791397841636235_DP ;  casimir_omega_weight(81) = 0.0053392062760162_DP
         casimir_omega(82) = 0.1845374763300464_DP ;  casimir_omega_weight(82) = 0.0054565896844392_DP
         casimir_omega(83) = 0.1900537996905553_DP ;  casimir_omega_weight(83) = 0.0055764816563694_DP
         casimir_omega(84) = 0.1956913022376211_DP ;  casimir_omega_weight(84) = 0.0056989615059779_DP
         casimir_omega(85) = 0.2014526126660132_DP ;  casimir_omega_weight(85) = 0.0058241113432751_DP
         casimir_omega(86) = 0.2073404432283353_DP ;  casimir_omega_weight(86) = 0.0059520161944520_DP
         casimir_omega(87) = 0.2133575927137955_DP ;  casimir_omega_weight(87) = 0.0060827641280609_DP
         casimir_omega(88) = 0.2195069495562098_DP ;  casimir_omega_weight(88) = 0.0062164463873548_DP
         casimir_omega(89) = 0.2257914950775780_DP ;  casimir_omega_weight(89) = 0.0063531575291480_DP
         casimir_omega(90) = 0.2322143068739121_DP ;  casimir_omega_weight(90) = 0.0064929955695383_DP
         casimir_omega(91) = 0.2387785623504040_DP ;  casimir_omega_weight(91) = 0.0066360621369090_DP
         casimir_omega(92) = 0.2454875424133799_DP ;  casimir_omega_weight(92) = 0.0067824626325878_DP
         casimir_omega(93) = 0.2523446353269457_DP ;  casimir_omega_weight(93) = 0.0069323063996341_DP
         casimir_omega(94) = 0.2593533407426636_DP ;  casimir_omega_weight(94) = 0.0070857069002055_DP
         casimir_omega(95) = 0.2665172739110902_DP ;  casimir_omega_weight(95) = 0.0072427819019916_DP
         casimir_omega(96) = 0.2738401700845178_DP ;  casimir_omega_weight(96) = 0.0074036536742741_DP
         casimir_omega(97) = 0.2813258891207950_DP ;  casimir_omega_weight(97) = 0.0075684491941377_DP
         casimir_omega(98) = 0.2889784202987027_DP ;  casimir_omega_weight(98) = 0.0077373003634572_DP
         casimir_omega(99) = 0.2968018873559702_DP ;  casimir_omega_weight(99) = 0.0079103442373007_DP
         casimir_omega(100) = 0.3048005537616696_DP ;  casimir_omega_weight(100) = 0.0080877232644188_DP
         casimir_omega(101) = 0.3129788282354531_DP ;  casimir_omega_weight(101) = 0.0082695855405443_DP
         casimir_omega(102) = 0.3213412705268208_DP ;  casimir_omega_weight(102) = 0.0084560850753019_DP
         casimir_omega(103) = 0.3298925974684363_DP ;  casimir_omega_weight(103) = 0.0086473820735283_DP
         casimir_omega(104) = 0.3386376893183436_DP ;  casimir_omega_weight(104) = 0.0088436432319098_DP
         casimir_omega(105) = 0.3475815964068630_DP ;  casimir_omega_weight(105) = 0.0090450420518683_DP
         casimir_omega(106) = 0.3567295461049211_DP ;  casimir_omega_weight(106) = 0.0092517591697168_DP
         casimir_omega(107) = 0.3660869501316169_DP ;  casimir_omega_weight(107) = 0.0094639827051593_DP
         casimir_omega(108) = 0.3756594122199428_DP ;  casimir_omega_weight(108) = 0.0096819086292973_DP
         casimir_omega(109) = 0.3854527361607700_DP ;  casimir_omega_weight(109) = 0.0099057411533818_DP
         casimir_omega(110) = 0.3954729342465121_DP ;  casimir_omega_weight(110) = 0.0101356931396285_DP
         casimir_omega(111) = 0.4057262361372300_DP ;  casimir_omega_weight(111) = 0.0103719865355330_DP
         casimir_omega(112) = 0.4162190981734186_DP ;  casimir_omega_weight(112) = 0.0106148528331972_DP
         casimir_omega(113) = 0.4269582131613159_DP ;  casimir_omega_weight(113) = 0.0108645335553115_DP
         casimir_omega(114) = 0.4379505206582387_DP ;  casimir_omega_weight(114) = 0.0111212807695386_DP
         casimir_omega(115) = 0.4492032177872837_DP ;  casimir_omega_weight(115) = 0.0113853576331861_DP
         casimir_omega(116) = 0.4607237706126892_DP ;  casimir_omega_weight(116) = 0.0116570389701955_DP
         casimir_omega(117) = 0.4725199261092370_DP ;  casimir_omega_weight(117) = 0.0119366118826199_DP
         casimir_omega(118) = 0.4845997247613428_DP ;  casimir_omega_weight(118) = 0.0122243763989218_DP
         casimir_omega(119) = 0.4969715138298940_DP ;  casimir_omega_weight(119) = 0.0125206461616121_DP
         casimir_omega(120) = 0.5096439613275323_DP ;  casimir_omega_weight(120) = 0.0128257491569526_DP
         casimir_omega(121) = 0.5226260707458463_DP ;  casimir_omega_weight(121) = 0.0131400284896016_DP
         casimir_omega(122) = 0.5359271965810181_DP ;  casimir_omega_weight(122) = 0.0134638432053921_DP
         casimir_omega(123) = 0.5495570607076846_DP ;  casimir_omega_weight(123) = 0.0137975691656008_DP
         casimir_omega(124) = 0.5635257696543422_DP ;  casimir_omega_weight(124) = 0.0141415999763691_DP
         casimir_omega(125) = 0.5778438328373756_DP ;  casimir_omega_weight(125) = 0.0144963479772362_DP
         casimir_omega(126) = 0.5925221818149423_DP ;  casimir_omega_weight(126) = 0.0148622452930375_DP
         casimir_omega(127) = 0.6075721906263346_DP ;  casimir_omega_weight(127) = 0.0152397449537799_DP
         casimir_omega(128) = 0.6230056972872728_DP ;  casimir_omega_weight(128) = 0.0156293220874775_DP
         casimir_omega(129) = 0.6388350265167440_DP ;  casimir_omega_weight(129) = 0.0160314751913374_DP
         casimir_omega(130) = 0.6550730137766135_DP ;  casimir_omega_weight(130) = 0.0164467274871561_DP
         casimir_omega(131) = 0.6717330307113426_DP ;  casimir_omega_weight(131) = 0.0168756283672079_DP
         casimir_omega(132) = 0.6888290120816961_DP ;  casimir_omega_weight(132) = 0.0173187549375418_DP
         casimir_omega(133) = 0.7063754842935127_DP ;  casimir_omega_weight(133) = 0.0177767136660853_DP
         casimir_omega(134) = 0.7243875956303253_DP ;  casimir_omega_weight(134) = 0.0182501421436651_DP
         casimir_omega(135) = 0.7428811483070773_DP ;  casimir_omega_weight(135) = 0.0187397109666958_DP
         casimir_omega(136) = 0.7618726324713235_DP ;  casimir_omega_weight(136) = 0.0192461257511015_DP
         casimir_omega(137) = 0.7813792622882415_DP ;  casimir_omega_weight(137) = 0.0197701292878455_DP
         casimir_omega(138) = 0.8014190142566502_DP ;  casimir_omega_weight(138) = 0.0203125038513575_DP
         casimir_omega(139) = 0.8220106679149993_DP ;  casimir_omega_weight(139) = 0.0208740736732065_DP
         casimir_omega(140) = 0.8431738491091703_DP ;  casimir_omega_weight(140) = 0.0214557075943827_DP
         casimir_omega(141) = 0.8649290760079581_DP ;  casimir_omega_weight(141) = 0.0220583219109092_DP
         casimir_omega(142) = 0.8872978080674008_DP ;  casimir_omega_weight(142) = 0.0226828834287444_DP
         casimir_omega(143) = 0.9103024981618578_DP ;  casimir_omega_weight(143) = 0.0233304127453999_DP
         casimir_omega(144) = 0.9339666481180371_DP ;  casimir_omega_weight(144) = 0.0240019877775380_DP
         casimir_omega(145) = 0.9583148679081315_DP ;  casimir_omega_weight(145) = 0.0246987475552789_DP
         casimir_omega(146) = 0.9833729387801754_DP ;  casimir_omega_weight(146) = 0.0254218963063022_DP
         casimir_omega(147) = 1.0091678806277433_DP ;  casimir_omega_weight(147) = 0.0261727078548030_DP
         casimir_omega(148) = 1.0357280239273683_DP ;  casimir_omega_weight(148) = 0.0269525303628424_DP
         casimir_omega(149) = 1.0630830866010732_DP ;  casimir_omega_weight(149) = 0.0277627914444561_DP
         casimir_omega(150) = 1.0912642561930908_DP ;  casimir_omega_weight(150) = 0.0286050036857476_DP
         casimir_omega(151) = 1.1203042777847994_DP ;  casimir_omega_weight(151) = 0.0294807706076511_DP
         casimir_omega(152) = 1.1502375481103555_DP ;  casimir_omega_weight(152) = 0.0303917931116049_DP
         casimir_omega(153) = 1.1811002163778599_DP ;  casimir_omega_weight(153) = 0.0313398764527133_DP
         casimir_omega(154) = 1.2129302923476124_DP ;  casimir_omega_weight(154) = 0.0323269377893561_DP
         casimir_omega(155) = 1.2457677622705727_DP ;  casimir_omega_weight(155) = 0.0333550143634781_DP
         casimir_omega(156) = 1.2796547133471405_DP ;  casimir_omega_weight(156) = 0.0344262723714241_DP
         casimir_omega(157) = 1.3146354674293752_DP ;  casimir_omega_weight(157) = 0.0355430165916202_DP
         casimir_omega(158) = 1.3507567247595937_DP ;  casimir_omega_weight(158) = 0.0367077008424871_DP
         casimir_omega(159) = 1.3880677186155861_DP ;  casimir_omega_weight(159) = 0.0379229393520905_DP
         casimir_omega(160) = 1.4266203818185867_DP ;  casimir_omega_weight(160) = 0.0391915191298458_DP
         casimir_omega(161) = 1.4664695261553837_DP ;  casimir_omega_weight(161) = 0.0405164134407589_DP
         casimir_omega(162) = 1.5076730358720600_DP ;  casimir_omega_weight(162) = 0.0419007964940628_DP
         casimir_omega(163) = 1.5502920765147916_DP ;  casimir_omega_weight(163) = 0.0433480594706273_DP
         casimir_omega(164) = 1.5943913205247620_DP ;  casimir_omega_weight(164) = 0.0448618280281934_DP
         casimir_omega(165) = 1.6400391911410177_DP ;  casimir_omega_weight(165) = 0.0464459814393224_DP
         casimir_omega(166) = 1.6873081263290921_DP ;  casimir_omega_weight(166) = 0.0481046735355163_DP
         casimir_omega(167) = 1.7362748646367439_DP ;  casimir_omega_weight(167) = 0.0498423556515614_DP
         casimir_omega(168) = 1.7870207550836847_DP ;  casimir_omega_weight(168) = 0.0516638017874994_DP
         casimir_omega(169) = 1.8396320934226520_DP ;  casimir_omega_weight(169) = 0.0535741362324208_DP
         casimir_omega(170) = 1.8942004873680522_DP ;  casimir_omega_weight(170) = 0.0555788639242835_DP
         casimir_omega(171) = 1.9508232536795855_DP ;  casimir_omega_weight(171) = 0.0576839038545969_DP
         casimir_omega(172) = 2.0096038503160294_DP ;  casimir_omega_weight(172) = 0.0598956258656211_DP
         casimir_omega(173) = 2.0706523472442364_DP ;  casimir_omega_weight(173) = 0.0622208912329836_DP
         casimir_omega(174) = 2.1340859399058312_DP ;  casimir_omega_weight(174) = 0.0646670974771711_DP
         casimir_omega(175) = 2.2000295098166580_DP ;  casimir_omega_weight(175) = 0.0672422279064767_DP
         casimir_omega(176) = 2.2686162373089926_DP ;  casimir_omega_weight(176) = 0.0699549064607142_DP
         casimir_omega(177) = 2.3399882720338265_DP ;  casimir_omega_weight(177) = 0.0728144585023109_DP
         casimir_omega(178) = 2.4142974675305187_DP ;  casimir_omega_weight(178) = 0.0758309782907687_DP
         casimir_omega(179) = 2.4917061869567743_DP ;  casimir_omega_weight(179) = 0.0790154039777577_DP
         casimir_omega(180) = 2.5723881879675541_DP ;  casimir_omega_weight(180) = 0.0823796010801161_DP
         casimir_omega(181) = 2.6565295957549462_DP ;  casimir_omega_weight(181) = 0.0859364555236165_DP
         casimir_omega(182) = 2.7443299744319289_DP ;  casimir_omega_weight(182) = 0.0896999775102685_DP
         casimir_omega(183) = 2.8360035082856689_DP ;  casimir_omega_weight(183) = 0.0936854176464283_DP
         casimir_omega(184) = 2.9317803059679464_DP ;  casimir_omega_weight(184) = 0.0979093969842642_DP
         casimir_omega(185) = 3.0319078424654142_DP ;  casimir_omega_weight(185) = 0.1023900528802587_DP
         casimir_omega(186) = 3.1366525557390141_DP ;  casimir_omega_weight(186) = 0.1071472028684354_DP
         casimir_omega(187) = 3.2463016172875130_DP ;  casimir_omega_weight(187) = 0.1122025290905266_DP
         casimir_omega(188) = 3.3611648986289664_DP ;  casimir_omega_weight(188) = 0.1175797862305049_DP
         casimir_omega(189) = 3.4815771588733981_DP ;  casimir_omega_weight(189) = 0.1233050363771917_DP
         casimir_omega(190) = 3.6079004822579011_DP ;  casimir_omega_weight(190) = 0.1294069148031240_DP
         casimir_omega(191) = 3.7405269988271312_DP ;  casimir_omega_weight(191) = 0.1359169313140567_DP
         casimir_omega(192) = 3.8798819264811684_DP ;  casimir_omega_weight(192) = 0.1428698126157524_DP
         casimir_omega(193) = 4.0264269785173292_DP ;  casimir_omega_weight(193) = 0.1503038920870972_DP
         casimir_omega(194) = 4.1806641877268600_DP ;  casimir_omega_weight(194) = 0.1582615544743195_DP
         casimir_omega(195) = 4.3431402062756845_DP ;  casimir_omega_weight(195) = 0.1667897443674665_DP
         casimir_omega(196) = 4.5144511502435085_DP ;  casimir_omega_weight(196) = 0.1759405489382908_DP
         casimir_omega(197) = 4.6952480691184979_DP ;  casimir_omega_weight(197) = 0.1857718673673657_DP
         casimir_omega(198) = 4.8862431341123527_DP ;  casimir_omega_weight(198) = 0.1963481817420015_DP
         casimir_omega(199) = 5.0882166553229000_DP ;  casimir_omega_weight(199) = 0.2077414470618061_DP
         casimir_omega(200) = 5.3020250570881444_DP ;  casimir_omega_weight(200) = 0.2200321214611898_DP
         casimir_omega(201) = 5.5286099640328725_DP ;  casimir_omega_weight(201) = 0.2333103619974893_DP
         casimir_omega(202) = 5.7690085781682212_DP ;  casimir_omega_weight(202) = 0.2476774165514434_DP
         casimir_omega(203) = 6.0243655610345304_DP ;  casimir_omega_weight(203) = 0.2632472487764412_DP
         casimir_omega(204) = 6.2959466756226883_DP ;  casimir_omega_weight(204) = 0.2801484409288560_DP
         casimir_omega(205) = 6.5851544923594556_DP ;  casimir_omega_weight(205) = 0.2985264292040168_DP
         casimir_omega(206) = 6.8935465239390874_DP ;  casimir_omega_weight(206) = 0.3185461383953368_DP
         casimir_omega(207) = 7.2228562279443604_DP ;  casimir_omega_weight(207) = 0.3403950979572412_DP
         casimir_omega(208) = 7.5750174075006171_DP ;  casimir_omega_weight(208) = 0.3642871407293407_DP
         casimir_omega(209) = 7.9521926531046203_DP ;  casimir_omega_weight(209) = 0.3904668098150075_DP
         casimir_omega(210) = 8.3568066090242343_DP ;  casimir_omega_weight(210) = 0.4192146298635678_DP
         casimir_omega(211) = 8.7915850227445311_DP ;  casimir_omega_weight(211) = 0.4508534382787650_DP
         casimir_omega(212) = 9.2596007555949829_DP ;  casimir_omega_weight(212) = 0.4857560222707728_DP
         casimir_omega(213) = 9.7643282097470721_DP ;  casimir_omega_weight(213) = 0.5243543727361660_DP
         casimir_omega(214) = 10.3097079781496905_DP ;  casimir_omega_weight(214) = 0.5671509504665067_DP
         casimir_omega(215) = 10.9002239722214949_DP ;  casimir_omega_weight(215) = 0.6147324706681710_DP
         casimir_omega(216) = 11.5409958574346305_DP ;  casimir_omega_weight(216) = 0.6677868571787802_DP
         casimir_omega(217) = 12.2378903700996435_DP ;  casimir_omega_weight(217) = 0.7271242104854779_DP
         casimir_omega(218) = 12.9976560551583606_DP ;  casimir_omega_weight(218) = 0.7937028909872454_DP
         casimir_omega(219) = 13.8280872307319598_DP ;  casimir_omega_weight(219) = 0.8686621652845476_DP
         casimir_omega(220) = 14.7382246558530348_DP ;  casimir_omega_weight(220) = 0.9533633332682754_DP
         casimir_omega(221) = 15.7386026002432224_DP ;  casimir_omega_weight(221) = 1.0494418972182642_DP
         casimir_omega(222) = 16.8415549964922313_DP ;  casimir_omega_weight(222) = 1.1588742231894875_DP
         casimir_omega(223) = 18.0615973907172851_DP ;  casimir_omega_weight(223) = 1.2840633856721608_DP
         casimir_omega(224) = 19.4159069228186816_DP ;  casimir_omega_weight(224) = 1.4279506361666596_DP
         casimir_omega(225) = 20.9249301808322130_DP ;  casimir_omega_weight(225) = 1.5941614314135204_DP
         casimir_omega(226) = 22.6131593984701560_DP ;  casimir_omega_weight(226) = 1.7871985579638039_DP
         casimir_omega(227) = 24.5101324645262650_DP ;  casimir_omega_weight(227) = 2.0127001534068993_DP
         casimir_omega(228) = 26.6517336538557466_DP ;  casimir_omega_weight(228) = 2.2777882249374284_DP
         casimir_omega(229) = 29.0819030509497694_DP ;  casimir_omega_weight(229) = 2.5915449972420217_DP
         casimir_omega(230) = 31.8549082888249728_DP ;  casimir_omega_weight(230) = 2.9656723477996794_DP
         casimir_omega(231) = 35.0384003768757637_DP ;  casimir_omega_weight(231) = 3.4154174557494734_DP
         casimir_omega(232) = 38.7175788680148401_DP ;  casimir_omega_weight(232) = 3.9608919315621911_DP
         casimir_omega(233) = 43.0009516554921802_DP ;  casimir_omega_weight(233) = 4.6289830506501746_DP
         casimir_omega(234) = 48.0284272600033901_DP ;  casimir_omega_weight(234) = 5.4561736641262630_DP
         casimir_omega(235) = 53.9828850060204886_DP ;  casimir_omega_weight(235) = 6.4927871660648577_DP
         casimir_omega(236) = 61.1070423779394361_DP ;  casimir_omega_weight(236) = 7.8095216521801296_DP
         casimir_omega(237) = 69.7285838794104791_DP ;  casimir_omega_weight(237) = 9.5077611330875609_DP
         casimir_omega(238) = 80.2985212625510343_DP ;  casimir_omega_weight(238) = 11.7363085955365865_DP
         casimir_omega(239) = 93.4513897550053088_DP ;  casimir_omega_weight(239) = 14.7194143465058858_DP
         casimir_omega(240) = 110.1027323746941278_DP ;  casimir_omega_weight(240) = 18.8054543658508457_DP
         casimir_omega(241) = 131.6128082127932544_DP ;  casimir_omega_weight(241) = 24.5550779215843278_DP
         casimir_omega(242) = 160.0734068950819164_DP ;  casimir_omega_weight(242) = 32.9088024618753678_DP
         casimir_omega(243) = 198.8361499256824629_DP ;  casimir_omega_weight(243) = 45.5245977954449188_DP
         casimir_omega(244) = 253.5459858636792205_DP ;  casimir_omega_weight(244) = 65.5068267041511092_DP
         casimir_omega(245) = 334.3179590320465877_DP ;  casimir_omega_weight(245) = 99.1206847285650383_DP
         casimir_omega(246) = 460.7687616357744673_DP ;  casimir_omega_weight(246) = 160.2841232966439691_DP
         casimir_omega(247) = 675.1450220731150011_DP ;  casimir_omega_weight(247) = 284.1226398203654071_DP
         casimir_omega(248) = 1082.7451023809144317_DP ;  casimir_omega_weight(248) = 576.6471228083744336_DP
         casimir_omega(249) = 2010.6404643375817614_DP ;  casimir_omega_weight(249) = 1457.7641702106682260_DP
         casimir_omega(250) = 4941.9879606803797287_DP ;  casimir_omega_weight(250) = 5603.6669085300300139_DP
         casimir_omega(251) = 26040.7477291924733436_DP ;  casimir_omega_weight(251) = 66829.9945252838078886_DP
return
endsubroutine gauss_legendre_grid250

subroutine gauss_legendre_grid300()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000096066882805_DP ;  casimir_omega_weight(2) = 0.0000246541275687_DP
         casimir_omega(3) = 0.0000506193748488_DP ;  casimir_omega_weight(3) = 0.0000573953465602_DP
         casimir_omega(4) = 0.0001244137579577_DP ;  casimir_omega_weight(4) = 0.0000901975840592_DP
         casimir_omega(5) = 0.0002310223589981_DP ;  casimir_omega_weight(5) = 0.0001230237547932_DP
         casimir_omega(6) = 0.0003704705002375_DP ;  casimir_omega_weight(6) = 0.0001558776795078_DP
         casimir_omega(7) = 0.0005427891676228_DP ;  casimir_omega_weight(7) = 0.0001887659642331_DP
         casimir_omega(8) = 0.0007480162332249_DP ;  casimir_omega_weight(8) = 0.0002216956659691_DP
         casimir_omega(9) = 0.0009861966979529_DP ;  casimir_omega_weight(9) = 0.0002546739645496_DP
         casimir_omega(10) = 0.0012573827723148_DP ;  casimir_omega_weight(10) = 0.0002877080943173_DP
         casimir_omega(11) = 0.0015616339206273_DP ;  casimir_omega_weight(11) = 0.0003208053273348_DP
         casimir_omega(12) = 0.0018990168964757_DP ;  casimir_omega_weight(12) = 0.0003539729699922_DP
         casimir_omega(13) = 0.0022696057772505_DP ;  casimir_omega_weight(13) = 0.0003872183638141_DP
         casimir_omega(14) = 0.0026734820003757_DP ;  casimir_omega_weight(14) = 0.0004205488877994_DP
         casimir_omega(15) = 0.0031107344022342_DP ;  casimir_omega_weight(15) = 0.0004539719613958_DP
         casimir_omega(16) = 0.0035814592602355_DP ;  casimir_omega_weight(16) = 0.0004874950477688_DP
         casimir_omega(17) = 0.0040857603382384_DP ;  casimir_omega_weight(17) = 0.0005211256572307_DP
         casimir_omega(18) = 0.0046237489354519_DP ;  casimir_omega_weight(18) = 0.0005548713507665_DP
         casimir_omega(19) = 0.0051955439389070_DP ;  casimir_omega_weight(19) = 0.0005887397436394_DP
         casimir_omega(20) = 0.0058012718795587_DP ;  casimir_omega_weight(20) = 0.0006227385090514_DP
         casimir_omega(21) = 0.0064410669920827_DP ;  casimir_omega_weight(21) = 0.0006568753818642_DP
         casimir_omega(22) = 0.0071150712784187_DP ;  casimir_omega_weight(22) = 0.0006911581623795_DP
         casimir_omega(23) = 0.0078234345751200_DP ;  casimir_omega_weight(23) = 0.0007255947201681_DP
         casimir_omega(24) = 0.0085663146245678_DP ;  casimir_omega_weight(24) = 0.0007601929979661_DP
         casimir_omega(25) = 0.0093438771501064_DP ;  casimir_omega_weight(25) = 0.0007949610156301_DP
         casimir_omega(26) = 0.0101562959351677_DP ;  casimir_omega_weight(26) = 0.0008299068741533_DP
         casimir_omega(27) = 0.0110037529064461_DP ;  casimir_omega_weight(27) = 0.0008650387597546_DP
         casimir_omega(28) = 0.0118864382211972_DP ;  casimir_omega_weight(28) = 0.0009003649480370_DP
         casimir_omega(29) = 0.0128045503587309_DP ;  casimir_omega_weight(29) = 0.0009358938082166_DP
         casimir_omega(30) = 0.0137582962161740_DP ;  casimir_omega_weight(30) = 0.0009716338074374_DP
         casimir_omega(31) = 0.0147478912085872_DP ;  casimir_omega_weight(31) = 0.0010075935151616_DP
         casimir_omega(32) = 0.0157735593735147_DP ;  casimir_omega_weight(32) = 0.0010437816076501_DP
         casimir_omega(33) = 0.0168355334800617_DP ;  casimir_omega_weight(33) = 0.0010802068725300_DP
         casimir_omega(34) = 0.0179340551425827_DP ;  casimir_omega_weight(34) = 0.0011168782134604_DP
         casimir_omega(35) = 0.0190693749390867_DP ;  casimir_omega_weight(35) = 0.0011538046548923_DP
         casimir_omega(36) = 0.0202417525344539_DP ;  casimir_omega_weight(36) = 0.0011909953469372_DP
         casimir_omega(37) = 0.0214514568085721_DP ;  casimir_omega_weight(37) = 0.0012284595703378_DP
         casimir_omega(38) = 0.0226987659895017_DP ;  casimir_omega_weight(38) = 0.0012662067415527_DP
         casimir_omega(39) = 0.0239839677917843_DP ;  casimir_omega_weight(39) = 0.0013042464179631_DP
         casimir_omega(40) = 0.0253073595600178_DP ;  casimir_omega_weight(40) = 0.0013425883031981_DP
         casimir_omega(41) = 0.0266692484178167_DP ;  casimir_omega_weight(41) = 0.0013812422525837_DP
         casimir_omega(42) = 0.0280699514222963_DP ;  casimir_omega_weight(42) = 0.0014202182787353_DP
         casimir_omega(43) = 0.0295097957242081_DP ;  casimir_omega_weight(43) = 0.0014595265572781_DP
         casimir_omega(44) = 0.0309891187338750_DP ;  casimir_omega_weight(44) = 0.0014991774327194_DP
         casimir_omega(45) = 0.0325082682930712_DP ;  casimir_omega_weight(45) = 0.0015391814244691_DP
         casimir_omega(46) = 0.0340676028529966_DP ;  casimir_omega_weight(46) = 0.0015795492330144_DP
         casimir_omega(47) = 0.0356674916585124_DP ;  casimir_omega_weight(47) = 0.0016202917462599_DP
         casimir_omega(48) = 0.0373083149387992_DP ;  casimir_omega_weight(48) = 0.0016614200460395_DP
         casimir_omega(49) = 0.0389904641046126_DP ;  casimir_omega_weight(49) = 0.0017029454147981_DP
         casimir_omega(50) = 0.0407143419523147_DP ;  casimir_omega_weight(50) = 0.0017448793424655_DP
         casimir_omega(51) = 0.0424803628748725_DP ;  casimir_omega_weight(51) = 0.0017872335335128_DP
         casimir_omega(52) = 0.0442889530800116_DP ;  casimir_omega_weight(52) = 0.0018300199142153_DP
         casimir_omega(53) = 0.0461405508157357_DP ;  casimir_omega_weight(53) = 0.0018732506401143_DP
         casimir_omega(54) = 0.0480356066034164_DP ;  casimir_omega_weight(54) = 0.0019169381037022_DP
         casimir_omega(55) = 0.0499745834786732_DP ;  casimir_omega_weight(55) = 0.0019610949423212_DP
         casimir_omega(56) = 0.0519579572402770_DP ;  casimir_omega_weight(56) = 0.0020057340463054_DP
         casimir_omega(57) = 0.0539862167073055_DP ;  casimir_omega_weight(57) = 0.0020508685673530_DP
         casimir_omega(58) = 0.0560598639848039_DP ;  casimir_omega_weight(58) = 0.0020965119271587_DP
         casimir_omega(59) = 0.0581794147382024_DP ;  casimir_omega_weight(59) = 0.0021426778263028_DP
         casimir_omega(60) = 0.0603453984767602_DP ;  casimir_omega_weight(60) = 0.0021893802534131_DP
         casimir_omega(61) = 0.0625583588463098_DP ;  casimir_omega_weight(61) = 0.0022366334946062_DP
         casimir_omega(62) = 0.0648188539315913_DP ;  casimir_omega_weight(62) = 0.0022844521432286_DP
         casimir_omega(63) = 0.0671274565684754_DP ;  casimir_omega_weight(63) = 0.0023328511098960_DP
         casimir_omega(64) = 0.0694847546663874_DP ;  casimir_omega_weight(64) = 0.0023818456328532_DP
         casimir_omega(65) = 0.0718913515412527_DP ;  casimir_omega_weight(65) = 0.0024314512886665_DP
         casimir_omega(66) = 0.0743478662593039_DP ;  casimir_omega_weight(66) = 0.0024816840032528_DP
         casimir_omega(67) = 0.0768549339920978_DP ;  casimir_omega_weight(67) = 0.0025325600632725_DP
         casimir_omega(68) = 0.0794132063831065_DP ;  casimir_omega_weight(68) = 0.0025840961278863_DP
         casimir_omega(69) = 0.0820233519262601_DP ;  casimir_omega_weight(69) = 0.0026363092409073_DP
         casimir_omega(70) = 0.0846860563568375_DP ;  casimir_omega_weight(70) = 0.0026892168433473_DP
         casimir_omega(71) = 0.0874020230551133_DP ;  casimir_omega_weight(71) = 0.0027428367863805_DP
         casimir_omega(72) = 0.0901719734631844_DP ;  casimir_omega_weight(72) = 0.0027971873447500_DP
         casimir_omega(73) = 0.0929966475154263_DP ;  casimir_omega_weight(73) = 0.0028522872306165_DP
         casimir_omega(74) = 0.0958768040830360_DP ;  casimir_omega_weight(74) = 0.0029081556078881_DP
         casimir_omega(75) = 0.0988132214331387_DP ;  casimir_omega_weight(75) = 0.0029648121070290_DP
         casimir_omega(76) = 0.1018066977029697_DP ;  casimir_omega_weight(76) = 0.0030222768403869_DP
         casimir_omega(77) = 0.1048580513896365_DP ;  casimir_omega_weight(77) = 0.0030805704180512_DP
         casimir_omega(78) = 0.1079681218560177_DP ;  casimir_omega_weight(78) = 0.0031397139642541_DP
         casimir_omega(79) = 0.1111377698533499_DP ;  casimir_omega_weight(79) = 0.0031997291343608_DP
         casimir_omega(80) = 0.1143678780611022_DP ;  casimir_omega_weight(80) = 0.0032606381324512_DP
         casimir_omega(81) = 0.1176593516447394_DP ;  casimir_omega_weight(81) = 0.0033224637295258_DP
         casimir_omega(82) = 0.1210131188320177_DP ;  casimir_omega_weight(82) = 0.0033852292823665_DP
         casimir_omega(83) = 0.1244301315084775_DP ;  casimir_omega_weight(83) = 0.0034489587530719_DP
         casimir_omega(84) = 0.1279113658328194_DP ;  casimir_omega_weight(84) = 0.0035136767292983_DP
         casimir_omega(85) = 0.1314578228728890_DP ;  casimir_omega_weight(85) = 0.0035794084452406_DP
         casimir_omega(86) = 0.1350705292630226_DP ;  casimir_omega_weight(86) = 0.0036461798033741_DP
         casimir_omega(87) = 0.1387505378835350_DP ;  casimir_omega_weight(87) = 0.0037140173970052_DP
         casimir_omega(88) = 0.1424989285631654_DP ;  casimir_omega_weight(88) = 0.0037829485336457_DP
         casimir_omega(89) = 0.1463168088053361_DP ;  casimir_omega_weight(89) = 0.0038530012592656_DP
         casimir_omega(90) = 0.1502053145391104_DP ;  casimir_omega_weight(90) = 0.0039242043834476_DP
         casimir_omega(91) = 0.1541656108957775_DP ;  casimir_omega_weight(91) = 0.0039965875054910_DP
         casimir_omega(92) = 0.1581988930120284_DP ;  casimir_omega_weight(92) = 0.0040701810415020_DP
         casimir_omega(93) = 0.1623063868607373_DP ;  casimir_omega_weight(93) = 0.0041450162525136_DP
         casimir_omega(94) = 0.1664893501103958_DP ;  casimir_omega_weight(94) = 0.0042211252736880_DP
         casimir_omega(95) = 0.1707490730143040_DP ;  casimir_omega_weight(95) = 0.0042985411446277_DP
         casimir_omega(96) = 0.1750868793306603_DP ;  casimir_omega_weight(96) = 0.0043772978408770_DP
         casimir_omega(97) = 0.1795041272747547_DP ;  casimir_omega_weight(97) = 0.0044574303066338_DP
         casimir_omega(98) = 0.1840022105045157_DP ;  casimir_omega_weight(98) = 0.0045389744887499_DP
         casimir_omega(99) = 0.1885825591407151_DP ;  casimir_omega_weight(99) = 0.0046219673720656_DP
         casimir_omega(100) = 0.1932466408232027_DP ;  casimir_omega_weight(100) = 0.0047064470161389_DP
         casimir_omega(101) = 0.1979959618045957_DP ;  casimir_omega_weight(101) = 0.0047924525934362_DP
         casimir_omega(102) = 0.2028320680829177_DP ;  casimir_omega_weight(102) = 0.0048800244290535_DP
         casimir_omega(103) = 0.2077565465747463_DP ;  casimir_omega_weight(103) = 0.0049692040420256_DP
         casimir_omega(104) = 0.2127710263305038_DP ;  casimir_omega_weight(104) = 0.0050600341883152_DP
         casimir_omega(105) = 0.2178771797935988_DP ;  casimir_omega_weight(105) = 0.0051525589055374_DP
         casimir_omega(106) = 0.2230767241052056_DP ;  casimir_omega_weight(106) = 0.0052468235595222_DP
         casimir_omega(107) = 0.2283714224565474_DP ;  casimir_omega_weight(107) = 0.0053428748927809_DP
         casimir_omega(108) = 0.2337630854906423_DP ;  casimir_omega_weight(108) = 0.0054407610749792_DP
         casimir_omega(109) = 0.2392535727555686_DP ;  casimir_omega_weight(109) = 0.0055405317555068_DP
         casimir_omega(110) = 0.2448447942113764_DP ;  casimir_omega_weight(110) = 0.0056422381182383_DP
         casimir_omega(111) = 0.2505387117929135_DP ;  casimir_omega_weight(111) = 0.0057459329385995_DP
         casimir_omega(112) = 0.2563373410309044_DP ;  casimir_omega_weight(112) = 0.0058516706430380_DP
         casimir_omega(113) = 0.2622427527337565_DP ;  casimir_omega_weight(113) = 0.0059595073710217_DP
         casimir_omega(114) = 0.2682570747326766_DP ;  casimir_omega_weight(114) = 0.0060695010396832_DP
         casimir_omega(115) = 0.2743824936928100_DP ;  casimir_omega_weight(115) = 0.0061817114112413_DP
         casimir_omega(116) = 0.2806212569932438_DP ;  casimir_omega_weight(116) = 0.0062962001633332_DP
         casimir_omega(117) = 0.2869756746788553_DP ;  casimir_omega_weight(117) = 0.0064130309623945_DP
         casimir_omega(118) = 0.2934481214871355_DP ;  casimir_omega_weight(118) = 0.0065322695402574_DP
         casimir_omega(119) = 0.3000410389532684_DP ;  casimir_omega_weight(119) = 0.0066539837740933_DP
         casimir_omega(120) = 0.3067569375969161_DP ;  casimir_omega_weight(120) = 0.0067782437698995_DP
         casimir_omega(121) = 0.3135983991943223_DP ;  casimir_omega_weight(121) = 0.0069051219496935_DP
         casimir_omega(122) = 0.3205680791395394_DP ;  casimir_omega_weight(122) = 0.0070346931425906_DP
         casimir_omega(123) = 0.3276687088987680_DP ;  casimir_omega_weight(123) = 0.0071670346799882_DP
         casimir_omega(124) = 0.3349030985620081_DP ;  casimir_omega_weight(124) = 0.0073022264950359_DP
         casimir_omega(125) = 0.3422741394964242_DP ;  casimir_omega_weight(125) = 0.0074403512266338_DP
         casimir_omega(126) = 0.3497848071060679_DP ;  casimir_omega_weight(126) = 0.0075814943281918_DP
         casimir_omega(127) = 0.3574381637028284_DP ;  casimir_omega_weight(127) = 0.0077257441813661_DP
         casimir_omega(128) = 0.3652373614937455_DP ;  casimir_omega_weight(128) = 0.0078731922150832_DP
         casimir_omega(129) = 0.3731856456900727_DP ;  casimir_omega_weight(129) = 0.0080239330300817_DP
         casimir_omega(130) = 0.3812863577437809_DP ;  casimir_omega_weight(130) = 0.0081780645292913_DP
         casimir_omega(131) = 0.3895429387174810_DP ;  casimir_omega_weight(131) = 0.0083356880543412_DP
         casimir_omega(132) = 0.3979589327940624_DP ;  casimir_omega_weight(132) = 0.0084969085285437_DP
         casimir_omega(133) = 0.4065379909326853_DP ;  casimir_omega_weight(133) = 0.0086618346066813_DP
         casimir_omega(134) = 0.4152838746781241_DP ;  casimir_omega_weight(134) = 0.0088305788319811_DP
         casimir_omega(135) = 0.4242004601308308_DP ;  casimir_omega_weight(135) = 0.0090032578006456_DP
         casimir_omega(136) = 0.4332917420854910_DP ;  casimir_omega_weight(136) = 0.0091799923343812_DP
         casimir_omega(137) = 0.4425618383462786_DP ;  casimir_omega_weight(137) = 0.0093609076613427_DP
         casimir_omega(138) = 0.4520149942274594_DP ;  casimir_omega_weight(138) = 0.0095461336059629_DP
         casimir_omega(139) = 0.4616555872484722_DP ;  casimir_omega_weight(139) = 0.0097358047881668_DP
         casimir_omega(140) = 0.4714881320331369_DP ;  casimir_omega_weight(140) = 0.0099300608325005_DP
         casimir_omega(141) = 0.4815172854231680_DP ;  casimir_omega_weight(141) = 0.0101290465877066_DP
         casimir_omega(142) = 0.4917478518167560_DP ;  casimir_omega_weight(142) = 0.0103329123573895_DP
         casimir_omega(143) = 0.5021847887435887_DP ;  casimir_omega_weight(143) = 0.0105418141423411_DP
         casimir_omega(144) = 0.5128332126883329_DP ;  casimir_omega_weight(144) = 0.0107559138952446_DP
         casimir_omega(145) = 0.5236984051752971_DP ;  casimir_omega_weight(145) = 0.0109753797884550_DP
         casimir_omega(146) = 0.5347858191277176_DP ;  casimir_omega_weight(146) = 0.0112003864956098_DP
         casimir_omega(147) = 0.5461010855159215_DP ;  casimir_omega_weight(147) = 0.0114311154878903_DP
         casimir_omega(148) = 0.5576500203094199_DP ;  casimir_omega_weight(148) = 0.0116677553457859_DP
         casimir_omega(149) = 0.5694386317489150_DP ;  casimir_omega_weight(149) = 0.0119105020873115_DP
         casimir_omega(150) = 0.5814731279551305_DP ;  casimir_omega_weight(150) = 0.0121595595135954_DP
         casimir_omega(151) = 0.5937599248923786_DP ;  casimir_omega_weight(151) = 0.0124151395729708_DP
         casimir_omega(152) = 0.6063056547059010_DP ;  casimir_omega_weight(152) = 0.0126774627445915_DP
         casimir_omega(153) = 0.6191171744531236_DP ;  casimir_omega_weight(153) = 0.0129467584428657_DP
         casimir_omega(154) = 0.6322015752502292_DP ;  casimir_omega_weight(154) = 0.0132232654438870_DP
         casimir_omega(155) = 0.6455661918567648_DP ;  casimir_omega_weight(155) = 0.0135072323352830_DP
         casimir_omega(156) = 0.6592186127223948_DP ;  casimir_omega_weight(156) = 0.0137989179909355_DP
         casimir_omega(157) = 0.6731666905214327_DP ;  casimir_omega_weight(157) = 0.0140985920720735_DP
         casimir_omega(158) = 0.6874185532023862_DP ;  casimir_omega_weight(158) = 0.0144065355564755_DP
         casimir_omega(159) = 0.7019826155814614_DP ;  casimir_omega_weight(159) = 0.0147230412975251_DP
         casimir_omega(160) = 0.7168675915108466_DP ;  casimir_omega_weight(160) = 0.0150484146150283_DP
         casimir_omega(161) = 0.7320825066545467_DP ;  casimir_omega_weight(161) = 0.0153829739198498_DP
         casimir_omega(162) = 0.7476367119066630_DP ;  casimir_omega_weight(162) = 0.0157270513745473_DP
         casimir_omega(163) = 0.7635398974892927_DP ;  casimir_omega_weight(163) = 0.0160809935923817_DP
         casimir_omega(164) = 0.7798021077696622_DP ;  casimir_omega_weight(164) = 0.0164451623772069_DP
         casimir_omega(165) = 0.7964337568387033_DP ;  casimir_omega_weight(165) = 0.0168199355069515_DP
         casimir_omega(166) = 0.8134456448961168_DP ;  casimir_omega_weight(166) = 0.0172057075635985_DP
         casimir_omega(167) = 0.8308489754899823_DP ;  casimir_omega_weight(167) = 0.0176028908128176_DP
         casimir_omega(168) = 0.8486553736621820_DP ;  casimir_omega_weight(168) = 0.0180119161366007_DP
         casimir_omega(169) = 0.8668769050544667_DP ;  casimir_omega_weight(169) = 0.0184332340224888_DP
         casimir_omega(170) = 0.8855260960336883_DP ;  casimir_omega_weight(170) = 0.0188673156133523_DP
         casimir_omega(171) = 0.9046159548987791_DP ;  casimir_omega_weight(171) = 0.0193146538218628_DP
         casimir_omega(172) = 0.9241599942364573_DP ;  casimir_omega_weight(172) = 0.0197757645142294_DP
         casimir_omega(173) = 0.9441722544972734_DP ;  casimir_omega_weight(173) = 0.0202511877679968_DP
         casimir_omega(174) = 0.9646673288687435_DP ;  casimir_omega_weight(174) = 0.0207414892093144_DP
         casimir_omega(175) = 0.9856603895277136_DP ;  casimir_omega_weight(175) = 0.0212472614351560_DP
         casimir_omega(176) = 1.0071672153600846_DP ;  casimir_omega_weight(176) = 0.0217691255268434_DP
         casimir_omega(177) = 1.0292042212423327_DP ;  casimir_omega_weight(177) = 0.0223077326612802_DP
         casimir_omega(178) = 1.0517884889862119_DP ;  casimir_omega_weight(178) = 0.0228637658273072_DP
         casimir_omega(179) = 1.0749378000554544_DP ;  casimir_omega_weight(179) = 0.0234379416547140_DP
         casimir_omega(180) = 1.0986706701713791_DP ;  casimir_omega_weight(180) = 0.0240310123644965_DP
         casimir_omega(181) = 1.1230063859330697_DP ;  casimir_omega_weight(181) = 0.0246437678492958_DP
         casimir_omega(182) = 1.1479650435872428_DP ;  casimir_omega_weight(182) = 0.0252770378940704_DP
         casimir_omega(183) = 1.1735675900932556_DP ;  casimir_omega_weight(183) = 0.0259316945475496_DP
         casimir_omega(184) = 1.1998358666397964_DP ;  casimir_omega_weight(184) = 0.0266086546562045_DP
         casimir_omega(185) = 1.2267926547820220_DP ;  casimir_omega_weight(185) = 0.0273088825732854_DP
         casimir_omega(186) = 1.2544617253809536_DP ;  casimir_omega_weight(186) = 0.0280333930567859_DP
         casimir_omega(187) = 1.2828678905414037_DP ;  casimir_omega_weight(187) = 0.0287832543710531_DP
         casimir_omega(188) = 1.3120370587601871_DP ;  casimir_omega_weight(188) = 0.0295595916086038_DP
         casimir_omega(189) = 1.3419962935134040_DP ;  casimir_omega_weight(189) = 0.0303635902495634_DP
         casimir_omega(190) = 1.3727738755300984_DP ;  casimir_omega_weight(190) = 0.0311964999782987_DP
         casimir_omega(191) = 1.4043993690197389_DP ;  casimir_omega_weight(191) = 0.0320596387781180_DP
         casimir_omega(192) = 1.4369036921430460_DP ;  casimir_omega_weight(192) = 0.0329543973272472_DP
         casimir_omega(193) = 1.4703191920397123_DP ;  casimir_omega_weight(193) = 0.0338822437210656_DP
         casimir_omega(194) = 1.5046797247529147_DP ;  casimir_omega_weight(194) = 0.0348447285482235_DP
         casimir_omega(195) = 1.5400207404193043_DP ;  casimir_omega_weight(195) = 0.0358434903507014_DP
         casimir_omega(196) = 1.5763793741246193_DP ;  casimir_omega_weight(196) = 0.0368802615008007_DP
         casimir_omega(197) = 1.6137945428597016_DP ;  casimir_omega_weight(197) = 0.0379568745312595_DP
         casimir_omega(198) = 1.6523070490495544_DP ;  casimir_omega_weight(198) = 0.0390752689580763_DP
         casimir_omega(199) = 1.6919596911696131_DP ;  casimir_omega_weight(199) = 0.0402374986398331_DP
         casimir_omega(200) = 1.7327973820092331_DP ;  casimir_omega_weight(200) = 0.0414457397211024_DP
         casimir_omega(201) = 1.7748672751925589_DP ;  casimir_omega_weight(201) = 0.0427022992129938_DP
         casimir_omega(202) = 1.8182189006222647_DP ;  casimir_omega_weight(202) = 0.0440096242686084_DP
         casimir_omega(203) = 1.8629043095727436_DP ;  casimir_omega_weight(203) = 0.0453703122175485_DP
         casimir_omega(204) = 1.9089782302263609_DP ;  casimir_omega_weight(204) = 0.0467871214299694_DP
         casimir_omega(205) = 1.9564982345207489_DP ;  casimir_omega_weight(205) = 0.0482629830882173_DP
         casimir_omega(206) = 2.0055249172570431_DP ;  casimir_omega_weight(206) = 0.0498010139522054_DP
         casimir_omega(207) = 2.0561220885096860_DP ;  casimir_omega_weight(207) = 0.0514045302139117_DP
         casimir_omega(208) = 2.1083569804788436_DP ;  casimir_omega_weight(208) = 0.0530770625465906_DP
         casimir_omega(209) = 2.1623004700378190_DP ;  casimir_omega_weight(209) = 0.0548223724660571_DP
         casimir_omega(210) = 2.2180273183512429_DP ;  casimir_omega_weight(210) = 0.0566444701339755_DP
         casimir_omega(211) = 2.2756164290772114_DP ;  casimir_omega_weight(211) = 0.0585476337477046_DP
         casimir_omega(212) = 2.3351511268189022_DP ;  casimir_omega_weight(212) = 0.0605364306777064_DP
         casimir_omega(213) = 2.3967194576611504_DP ;  casimir_omega_weight(213) = 0.0626157405315941_DP
         casimir_omega(214) = 2.4604145138167479_DP ;  casimir_omega_weight(214) = 0.0647907803447514_DP
         casimir_omega(215) = 2.5263347846185580_DP ;  casimir_omega_weight(215) = 0.0670671321204629_DP
         casimir_omega(216) = 2.5945845363293554_DP ;  casimir_omega_weight(216) = 0.0694507729693808_DP
         casimir_omega(217) = 2.6652742235056519_DP ;  casimir_omega_weight(217) = 0.0719481081270982_DP
         casimir_omega(218) = 2.7385209349472990_DP ;  casimir_omega_weight(218) = 0.0745660071633506_DP
         casimir_omega(219) = 2.8144488775963916_DP ;  casimir_omega_weight(219) = 0.0773118437335391_DP
         casimir_omega(220) = 2.8931899021216823_DP ;  casimir_omega_weight(220) = 0.0801935392677960_DP
         casimir_omega(221) = 2.9748840743434366_DP ;  casimir_omega_weight(221) = 0.0832196110414573_DP
         casimir_omega(222) = 3.0596802971257646_DP ;  casimir_omega_weight(222) = 0.0863992251278063_DP
         casimir_omega(223) = 3.1477369878950299_DP ;  casimir_omega_weight(223) = 0.0897422547979274_DP
         casimir_omega(224) = 3.2392228175446802_DP ;  casimir_omega_weight(224) = 0.0932593450063053_DP
         casimir_omega(225) = 3.3343175171656947_DP ;  casimir_omega_weight(225) = 0.0969619836849937_DP
         casimir_omega(226) = 3.4332127598127520_DP ;  casimir_omega_weight(226) = 0.1008625806658761_DP
         casimir_omega(227) = 3.5361131253891944_DP ;  casimir_omega_weight(227) = 0.1049745551617353_DP
         casimir_omega(228) = 3.6432371577278353_DP ;  casimir_omega_weight(228) = 0.1093124328642407_DP
         casimir_omega(229) = 3.7548185240740404_DP ;  casimir_omega_weight(229) = 0.1138919538654041_DP
         casimir_omega(230) = 3.8711072884673947_DP ;  casimir_omega_weight(230) = 0.1187301927782931_DP
         casimir_omega(231) = 3.9923713119906643_DP ;  casimir_omega_weight(231) = 0.1238456926308868_DP
         casimir_omega(232) = 4.1188977945395360_DP ;  casimir_omega_weight(232) = 0.1292586143353009_DP
         casimir_omega(233) = 4.2509949746990783_DP ;  casimir_omega_weight(233) = 0.1349909038008969_DP
         casimir_omega(234) = 4.3889940065317585_DP ;  casimir_omega_weight(234) = 0.1410664790693784_DP
         casimir_omega(235) = 4.5332510346362591_DP ;  casimir_omega_weight(235) = 0.1475114402109110_DP
         casimir_omega(236) = 4.6841494917816888_DP ;  casimir_omega_weight(236) = 0.1543543051445503_DP
         casimir_omega(237) = 4.8421026468254631_DP ;  casimir_omega_weight(237) = 0.1616262750409731_DP
         casimir_omega(238) = 5.0075564345653651_DP ;  casimir_omega_weight(238) = 0.1693615335484227_DP
         casimir_omega(239) = 5.1809926037509229_DP ;  casimir_omega_weight(239) = 0.1775975847700533_DP
         casimir_omega(240) = 5.3629322248009270_DP ;  casimir_omega_weight(240) = 0.1863756357305094_DP
         casimir_omega(241) = 5.5539396049788001_DP ;  casimir_omega_weight(241) = 0.1957410300303815_DP
         casimir_omega(242) = 5.7546266660292638_DP ;  casimir_omega_weight(242) = 0.2057437405258155_DP
         casimir_omega(243) = 5.9656578477750921_DP ;  casimir_omega_weight(243) = 0.2164389302287359_DP
         casimir_omega(244) = 6.1877556111545440_DP ;  casimir_omega_weight(244) = 0.2278875922397172_DP
         casimir_omega(245) = 6.4217066259309066_DP ;  casimir_omega_weight(245) = 0.2401572814665410_DP
         casimir_omega(246) = 6.6683687421883340_DP ;  casimir_omega_weight(246) = 0.2533229532089792_DP
         casimir_omega(247) = 6.9286788611645607_DP ;  casimir_omega_weight(247) = 0.2674679264908348_DP
         casimir_omega(248) = 7.2036618404960029_DP ;  casimir_omega_weight(248) = 0.2826849934150835_DP
         casimir_omega(249) = 7.4944405922084467_DP ;  casimir_omega_weight(249) = 0.2990776999183689_DP
         casimir_omega(250) = 7.8022475595853615_DP ;  casimir_omega_weight(250) = 0.3167618283034769_DP
         casimir_omega(251) = 8.1284377923684623_DP ;  casimir_omega_weight(251) = 0.3358671180286986_DP
         casimir_omega(252) = 8.4745038798372097_DP ;  casimir_omega_weight(252) = 0.3565392687097505_DP
         casimir_omega(253) = 8.8420930497080139_DP ;  casimir_omega_weight(253) = 0.3789422784890022_DP
         casimir_omega(254) = 9.2330267994274031_DP ;  casimir_omega_weight(254) = 0.4032611822882943_DP
         casimir_omega(255) = 9.6493234977389832_DP ;  casimir_omega_weight(255) = 0.4297052685482723_DP
         casimir_omega(256) = 10.0932244814609327_DP ;  casimir_omega_weight(256) = 0.4585118706142450_DP
         casimir_omega(257) = 10.5672242791316471_DP ;  casimir_omega_weight(257) = 0.4899508508777682_DP
         casimir_omega(258) = 11.0741057245651593_DP ;  casimir_omega_weight(258) = 0.5243299233966083_DP
         casimir_omega(259) = 11.6169808858254200_DP ;  casimir_omega_weight(259) = 0.5620009955728669_DP
         casimir_omega(260) = 12.1993389369578313_DP ;  casimir_omega_weight(260) = 0.6033677537496209_DP
         casimir_omega(261) = 12.8251023517643485_DP ;  casimir_omega_weight(261) = 0.6488947740805597_DP
         casimir_omega(262) = 13.4986931150072973_DP ;  casimir_omega_weight(262) = 0.6991185125657128_DP
         casimir_omega(263) = 14.2251110451189362_DP ;  casimir_omega_weight(263) = 0.7546606217679097_DP
         casimir_omega(264) = 15.0100268281430917_DP ;  casimir_omega_weight(264) = 0.8162441633551006_DP
         casimir_omega(265) = 15.8598930076862885_DP ;  casimir_omega_weight(265) = 0.8847134445952627_DP
         casimir_omega(266) = 16.7820770035603637_DP ;  casimir_omega_weight(266) = 0.9610584161781313_DP
         casimir_omega(267) = 17.7850213012553517_DP ;  casimir_omega_weight(267) = 1.0464448460641500_DP
         casimir_omega(268) = 18.8784373452172147_DP ;  casimir_omega_weight(268) = 1.1422518543784033_DP
         casimir_omega(269) = 20.0735414906371155_DP ;  casimir_omega_weight(269) = 1.2501188927783942_DP
         casimir_omega(270) = 21.3833437726432152_DP ;  casimir_omega_weight(270) = 1.3720049280429731_DP
         casimir_omega(271) = 22.8230034499678816_DP ;  casimir_omega_weight(271) = 1.5102635155643453_DP
         casimir_omega(272) = 24.4102695706346431_DP ;  casimir_omega_weight(272) = 1.6677387278530493_DP
         casimir_omega(273) = 26.1660306148076280_DP ;  casimir_omega_weight(273) = 1.8478886885641745_DP
         casimir_omega(274) = 28.1150052062961642_DP ;  casimir_omega_weight(274) = 2.0549459804043742_DP
         casimir_omega(275) = 30.2866168401912894_DP ;  casimir_omega_weight(275) = 2.2941277858073756_DP
         casimir_omega(276) = 32.7161108633330500_DP ;  casimir_omega_weight(276) = 2.5719138012130025_DP
         casimir_omega(277) = 35.4459935293397663_DP ;  casimir_omega_weight(277) = 2.8964175403088595_DP
         casimir_omega(278) = 38.5279038044619995_DP ;  casimir_omega_weight(278) = 3.2778878667250639_DP
         casimir_omega(279) = 42.0250732990289464_DP ;  casimir_omega_weight(279) = 3.7293944783812454_DP
         casimir_omega(280) = 46.0155953939807958_DP ;  casimir_omega_weight(280) = 4.2677768622796535_DP
         casimir_omega(281) = 50.5968227039343930_DP ;  casimir_omega_weight(281) = 4.9149763416096457_DP
         casimir_omega(282) = 55.8913609255278629_DP ;  casimir_omega_weight(282) = 5.6999343580835058_DP
         casimir_omega(283) = 62.0553574240321950_DP ;  casimir_omega_weight(283) = 6.6613428164006994_DP
         casimir_omega(284) = 69.2901463702646510_DP ;  casimir_omega_weight(284) = 7.8517020528514623_DP
         casimir_omega(285) = 77.8588987044135195_DP ;  casimir_omega_weight(285) = 9.3434295192966239_DP
         casimir_omega(286) = 88.1108949613990404_DP ;  casimir_omega_weight(286) = 11.2382627087052782_DP
         casimir_omega(287) = 100.5176867421101861_DP ;  casimir_omega_weight(287) = 13.6820974187852329_DP
         casimir_omega(288) = 115.7282986748884213_DP ;  casimir_omega_weight(288) = 16.8890673214191338_DP
         casimir_omega(289) = 134.6558532839881650_DP ;  casimir_omega_weight(289) = 21.1818779128854047_DP
         casimir_omega(290) = 158.6178549633651755_DP ;  casimir_omega_weight(290) = 27.0618566828827198_DP
         casimir_omega(291) = 189.5717729884722758_DP ;  casimir_omega_weight(291) = 35.3358011905796658_DP
         casimir_omega(292) = 230.5277794269387925_DP ;  casimir_omega_weight(292) = 47.3571550681692273_DP
         casimir_omega(293) = 286.3089966925517729_DP ;  casimir_omega_weight(293) = 65.5118056633009047_DP
         casimir_omega(294) = 365.0387399868098441_DP ;  casimir_omega_weight(294) = 94.2670598265547568_DP
         casimir_omega(295) = 481.2729777906967001_DP ;  casimir_omega_weight(295) = 142.6387938991264832_DP
         casimir_omega(296) = 663.2409441338610350_DP ;  casimir_omega_weight(296) = 230.6555174759750173_DP
         casimir_omega(297) = 971.7372901998334100_DP ;  casimir_omega_weight(297) = 408.8642787757368637_DP
         casimir_omega(298) = 1558.2907280520776112_DP ;  casimir_omega_weight(298) = 829.8191450188559202_DP
         casimir_omega(299) = 2893.5706621919057397_DP ;  casimir_omega_weight(299) = 2097.7831336152130461_DP
         casimir_omega(300) = 7111.9013435584374747_DP ;  casimir_omega_weight(300) = 8063.9091954789473675_DP
         casimir_omega(301) = 37473.8920927397848573_DP ;  casimir_omega_weight(301) = 96171.1350635193957714_DP
return
endsubroutine gauss_legendre_grid300

subroutine gauss_legendre_grid350()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000070613133644_DP ;  casimir_omega_weight(2) = 0.0000181217525802_DP
         casimir_omega(3) = 0.0000372068856723_DP ;  casimir_omega_weight(3) = 0.0000421868164474_DP
         casimir_omega(4) = 0.0000914461690992_DP ;  casimir_omega_weight(4) = 0.0000662942875886_DP
         casimir_omega(5) = 0.0001697999250297_DP ;  casimir_omega_weight(5) = 0.0000904155077678_DP
         casimir_omega(6) = 0.0002722822227977_DP ;  casimir_omega_weight(6) = 0.0001145518833845_DP
         casimir_omega(7) = 0.0003989098907577_DP ;  casimir_omega_weight(7) = 0.0001387068639199_DP
         casimir_omega(8) = 0.0005497034100371_DP ;  casimir_omega_weight(8) = 0.0001628842249711_DP
         casimir_omega(9) = 0.0007246870870078_DP ;  casimir_omega_weight(9) = 0.0001870878259172_DP
         casimir_omega(10) = 0.0009238891055657_DP ;  casimir_omega_weight(10) = 0.0002113215585461_DP
         casimir_omega(11) = 0.0011473415513960_DP ;  casimir_omega_weight(11) = 0.0002355893335617_DP
         casimir_omega(12) = 0.0013950804286417_DP ;  casimir_omega_weight(12) = 0.0002598950769205_DP
         casimir_omega(13) = 0.0016671456747252_DP ;  casimir_omega_weight(13) = 0.0002842427292500_DP
         casimir_omega(14) = 0.0019635811752297_DP ;  casimir_omega_weight(14) = 0.0003086362463824_DP
         casimir_omega(15) = 0.0022844347795713_DP ;  casimir_omega_weight(15) = 0.0003330796003375_DP
         casimir_omega(16) = 0.0026297583177712_DP ;  casimir_omega_weight(16) = 0.0003575767805132_DP
         casimir_omega(17) = 0.0029996076184690_DP ;  casimir_omega_weight(17) = 0.0003821317949704_DP
         casimir_omega(18) = 0.0033940425282596_DP ;  casimir_omega_weight(18) = 0.0004067486717842_DP
         casimir_omega(19) = 0.0038131269323920_DP ;  casimir_omega_weight(19) = 0.0004314314604210_DP
         casimir_omega(20) = 0.0042569287768630_DP ;  casimir_omega_weight(20) = 0.0004561842331485_DP
         casimir_omega(21) = 0.0047255200919323_DP ;  casimir_omega_weight(21) = 0.0004810110864621_DP
         casimir_omega(22) = 0.0052189770170713_DP ;  casimir_omega_weight(22) = 0.0005059161425334_DP
         casimir_omega(23) = 0.0057373798273712_DP ;  casimir_omega_weight(23) = 0.0005309035506746_DP
         casimir_omega(24) = 0.0062808129614208_DP ;  casimir_omega_weight(24) = 0.0005559774888209_DP
         casimir_omega(25) = 0.0068493650506830_DP ;  casimir_omega_weight(25) = 0.0005811421650326_DP
         casimir_omega(26) = 0.0074431289503747_DP ;  casimir_omega_weight(26) = 0.0006064018190157_DP
         casimir_omega(27) = 0.0080622017718811_DP ;  casimir_omega_weight(27) = 0.0006317607236573_DP
         casimir_omega(28) = 0.0087066849167171_DP ;  casimir_omega_weight(28) = 0.0006572231865923_DP
         casimir_omega(29) = 0.0093766841120610_DP ;  casimir_omega_weight(29) = 0.0006827935517769_DP
         casimir_omega(30) = 0.0100723094478774_DP ;  casimir_omega_weight(30) = 0.0007084762010979_DP
         casimir_omega(31) = 0.0107936754156575_DP ;  casimir_omega_weight(31) = 0.0007342755559946_DP
         casimir_omega(32) = 0.0115409009487963_DP ;  casimir_omega_weight(32) = 0.0007601960791106_DP
         casimir_omega(33) = 0.0123141094646322_DP ;  casimir_omega_weight(33) = 0.0007862422759723_DP
         casimir_omega(34) = 0.0131134289081765_DP ;  casimir_omega_weight(34) = 0.0008124186966856_DP
         casimir_omega(35) = 0.0139389917975595_DP ;  casimir_omega_weight(35) = 0.0008387299376701_DP
         casimir_omega(36) = 0.0147909352712164_DP ;  casimir_omega_weight(36) = 0.0008651806434149_DP
         casimir_omega(37) = 0.0156694011368510_DP ;  casimir_omega_weight(37) = 0.0008917755082671_DP
         casimir_omega(38) = 0.0165745359221987_DP ;  casimir_omega_weight(38) = 0.0009185192782488_DP
         casimir_omega(39) = 0.0175064909276263_DP ;  casimir_omega_weight(39) = 0.0009454167529109_DP
         casimir_omega(40) = 0.0184654222805985_DP ;  casimir_omega_weight(40) = 0.0009724727872164_DP
         casimir_omega(41) = 0.0194514909920453_DP ;  casimir_omega_weight(41) = 0.0009996922934564_DP
         casimir_omega(42) = 0.0204648630146687_DP ;  casimir_omega_weight(42) = 0.0010270802432110_DP
         casimir_omega(43) = 0.0215057093032217_DP ;  casimir_omega_weight(43) = 0.0010546416693383_DP
         casimir_omega(44) = 0.0225742058767987_DP ;  casimir_omega_weight(44) = 0.0010823816680054_DP
         casimir_omega(45) = 0.0236705338831792_DP ;  casimir_omega_weight(45) = 0.0011103054007624_DP
         casimir_omega(46) = 0.0247948796652622_DP ;  casimir_omega_weight(46) = 0.0011384180966554_DP
         casimir_omega(47) = 0.0259474348296362_DP ;  casimir_omega_weight(47) = 0.0011667250543815_DP
         casimir_omega(48) = 0.0271283963173301_DP ;  casimir_omega_weight(48) = 0.0011952316444932_DP
         casimir_omega(49) = 0.0283379664767871_DP ;  casimir_omega_weight(49) = 0.0012239433116451_DP
         casimir_omega(50) = 0.0295763531391095_DP ;  casimir_omega_weight(50) = 0.0012528655768897_DP
         casimir_omega(51) = 0.0308437696956295_DP ;  casimir_omega_weight(51) = 0.0012820040400271_DP
         casimir_omega(52) = 0.0321404351778453_DP ;  casimir_omega_weight(52) = 0.0013113643819960_DP
         casimir_omega(53) = 0.0334665743397881_DP ;  casimir_omega_weight(53) = 0.0013409523673341_DP
         casimir_omega(54) = 0.0348224177428631_DP ;  casimir_omega_weight(54) = 0.0013707738466757_DP
         casimir_omega(55) = 0.0362082018432263_DP ;  casimir_omega_weight(55) = 0.0014008347593212_DP
         casimir_omega(56) = 0.0376241690817575_DP ;  casimir_omega_weight(56) = 0.0014311411358567_DP
         casimir_omega(57) = 0.0390705679766822_DP ;  casimir_omega_weight(57) = 0.0014616991008426_DP
         casimir_omega(58) = 0.0405476532189127_DP ;  casimir_omega_weight(58) = 0.0014925148755550_DP
         casimir_omega(59) = 0.0420556857701663_DP ;  casimir_omega_weight(59) = 0.0015235947808029_DP
         casimir_omega(60) = 0.0435949329639338_DP ;  casimir_omega_weight(60) = 0.0015549452398088_DP
         casimir_omega(61) = 0.0451656686093601_DP ;  casimir_omega_weight(61) = 0.0015865727811563_DP
         casimir_omega(62) = 0.0467681730981127_DP ;  casimir_omega_weight(62) = 0.0016184840418178_DP
         casimir_omega(63) = 0.0484027335143114_DP ;  casimir_omega_weight(63) = 0.0016506857702486_DP
         casimir_omega(64) = 0.0500696437475955_DP ;  casimir_omega_weight(64) = 0.0016831848295666_DP
         casimir_omega(65) = 0.0517692046094038_DP ;  casimir_omega_weight(65) = 0.0017159882008040_DP
         casimir_omega(66) = 0.0535017239525540_DP ;  casimir_omega_weight(66) = 0.0017491029862489_DP
         casimir_omega(67) = 0.0552675167942023_DP ;  casimir_omega_weight(67) = 0.0017825364128682_DP
         casimir_omega(68) = 0.0570669054422717_DP ;  casimir_omega_weight(68) = 0.0018162958358225_DP
         casimir_omega(69) = 0.0589002196254383_DP ;  casimir_omega_weight(69) = 0.0018503887420662_DP
         casimir_omega(70) = 0.0607677966267700_DP ;  casimir_omega_weight(70) = 0.0018848227540513_DP
         casimir_omega(71) = 0.0626699814211119_DP ;  casimir_omega_weight(71) = 0.0019196056335220_DP
         casimir_omega(72) = 0.0646071268163191_DP ;  casimir_omega_weight(72) = 0.0019547452854136_DP
         casimir_omega(73) = 0.0665795935984418_DP ;  casimir_omega_weight(73) = 0.0019902497618557_DP
         casimir_omega(74) = 0.0685877506809636_DP ;  casimir_omega_weight(74) = 0.0020261272662814_DP
         casimir_omega(75) = 0.0706319752582099_DP ;  casimir_omega_weight(75) = 0.0020623861576532_DP
         casimir_omega(76) = 0.0727126529630362_DP ;  casimir_omega_weight(76) = 0.0020990349548010_DP
         casimir_omega(77) = 0.0748301780289125_DP ;  casimir_omega_weight(77) = 0.0021360823408793_DP
         casimir_omega(78) = 0.0769849534565304_DP ;  casimir_omega_weight(78) = 0.0021735371679512_DP
         casimir_omega(79) = 0.0791773911850569_DP ;  casimir_omega_weight(79) = 0.0022114084616975_DP
         casimir_omega(80) = 0.0814079122681600_DP ;  casimir_omega_weight(80) = 0.0022497054262612_DP
         casimir_omega(81) = 0.0836769470549511_DP ;  casimir_omega_weight(81) = 0.0022884374492257_DP
         casimir_omega(82) = 0.0859849353759736_DP ;  casimir_omega_weight(82) = 0.0023276141067377_DP
         casimir_omega(83) = 0.0883323267343881_DP ;  casimir_omega_weight(83) = 0.0023672451687756_DP
         casimir_omega(84) = 0.0907195805025004_DP ;  casimir_omega_weight(84) = 0.0024073406045676_DP
         casimir_omega(85) = 0.0931471661237881_DP ;  casimir_omega_weight(85) = 0.0024479105881719_DP
         casimir_omega(86) = 0.0956155633205842_DP ;  casimir_omega_weight(86) = 0.0024889655042125_DP
         casimir_omega(87) = 0.0981252623075858_DP ;  casimir_omega_weight(87) = 0.0025305159537901_DP
         casimir_omega(88) = 0.1006767640113541_DP ;  casimir_omega_weight(88) = 0.0025725727605622_DP
         casimir_omega(89) = 0.1032705802959867_DP ;  casimir_omega_weight(89) = 0.0026151469770047_DP
         casimir_omega(90) = 0.1059072341951460_DP ;  casimir_omega_weight(90) = 0.0026582498908588_DP
         casimir_omega(91) = 0.1085872601506294_DP ;  casimir_omega_weight(91) = 0.0027018930317777_DP
         casimir_omega(92) = 0.1113112042576801_DP ;  casimir_omega_weight(92) = 0.0027460881781585_DP
         casimir_omega(93) = 0.1140796245172436_DP ;  casimir_omega_weight(93) = 0.0027908473641991_DP
         casimir_omega(94) = 0.1168930910953778_DP ;  casimir_omega_weight(94) = 0.0028361828871558_DP
         casimir_omega(95) = 0.1197521865900349_DP ;  casimir_omega_weight(95) = 0.0028821073148262_DP
         casimir_omega(96) = 0.1226575063054413_DP ;  casimir_omega_weight(96) = 0.0029286334932684_DP
         casimir_omega(97) = 0.1256096585343135_DP ;  casimir_omega_weight(97) = 0.0029757745547489_DP
         casimir_omega(98) = 0.1286092648481454_DP ;  casimir_omega_weight(98) = 0.0030235439259428_DP
         casimir_omega(99) = 0.1316569603958220_DP ;  casimir_omega_weight(99) = 0.0030719553363899_DP
         casimir_omega(100) = 0.1347533942108216_DP ;  casimir_omega_weight(100) = 0.0031210228272116_DP
         casimir_omega(101) = 0.1378992295272713_DP ;  casimir_omega_weight(101) = 0.0031707607601065_DP
         casimir_omega(102) = 0.1410951441051402_DP ;  casimir_omega_weight(102) = 0.0032211838266310_DP
         casimir_omega(103) = 0.1443418305648571_DP ;  casimir_omega_weight(103) = 0.0032723070577692_DP
         casimir_omega(104) = 0.1476399967316535_DP ;  casimir_omega_weight(104) = 0.0033241458338147_DP
         casimir_omega(105) = 0.1509903659899452_DP ;  casimir_omega_weight(105) = 0.0033767158945691_DP
         casimir_omega(106) = 0.1543936776480739_DP ;  casimir_omega_weight(106) = 0.0034300333498646_DP
         casimir_omega(107) = 0.1578506873137412_DP ;  casimir_omega_weight(107) = 0.0034841146904297_DP
         casimir_omega(108) = 0.1613621672804903_DP ;  casimir_omega_weight(108) = 0.0035389767991184_DP
         casimir_omega(109) = 0.1649289069255878_DP ;  casimir_omega_weight(109) = 0.0035946369624885_DP
         casimir_omega(110) = 0.1685517131196823_DP ;  casimir_omega_weight(110) = 0.0036511128827715_DP
         casimir_omega(111) = 0.1722314106486306_DP ;  casimir_omega_weight(111) = 0.0037084226902471_DP
         casimir_omega(112) = 0.1759688426478909_DP ;  casimir_omega_weight(112) = 0.0037665849560007_DP
         casimir_omega(113) = 0.1797648710498989_DP ;  casimir_omega_weight(113) = 0.0038256187051414_DP
         casimir_omega(114) = 0.1836203770448719_DP ;  casimir_omega_weight(114) = 0.0038855434304333_DP
         casimir_omega(115) = 0.1875362615554755_DP ;  casimir_omega_weight(115) = 0.0039463791063985_DP
         casimir_omega(116) = 0.1915134457258325_DP ;  casimir_omega_weight(116) = 0.0040081462038998_DP
         casimir_omega(117) = 0.1955528714253578_DP ;  casimir_omega_weight(117) = 0.0040708657052135_DP
         casimir_omega(118) = 0.1996555017679198_DP ;  casimir_omega_weight(118) = 0.0041345591196126_DP
         casimir_omega(119) = 0.2038223216468582_DP ;  casimir_omega_weight(119) = 0.0041992484994967_DP
         casimir_omega(120) = 0.2080543382863973_DP ;  casimir_omega_weight(120) = 0.0042649564570660_DP
         casimir_omega(121) = 0.2123525818100283_DP ;  casimir_omega_weight(121) = 0.0043317061815780_DP
         casimir_omega(122) = 0.2167181058264428_DP ;  casimir_omega_weight(122) = 0.0043995214572034_DP
         casimir_omega(123) = 0.2211519880336328_DP ;  casimir_omega_weight(123) = 0.0044684266815020_DP
         casimir_omega(124) = 0.2256553308417868_DP ;  casimir_omega_weight(124) = 0.0045384468845568_DP
         casimir_omega(125) = 0.2302292620156584_DP ;  casimir_omega_weight(125) = 0.0046096077487743_DP
         casimir_omega(126) = 0.2348749353370856_DP ;  casimir_omega_weight(126) = 0.0046819356293791_DP
         casimir_omega(127) = 0.2395935312883585_DP ;  casimir_omega_weight(127) = 0.0047554575756698_DP
         casimir_omega(128) = 0.2443862577572265_DP ;  casimir_omega_weight(128) = 0.0048302013529892_DP
         casimir_omega(129) = 0.2492543507642687_DP ;  casimir_omega_weight(129) = 0.0049061954655176_DP
         casimir_omega(130) = 0.2541990752134671_DP ;  casimir_omega_weight(130) = 0.0049834691798865_DP
         casimir_omega(131) = 0.2592217256668110_DP ;  casimir_omega_weight(131) = 0.0050620525496312_DP
         casimir_omega(132) = 0.2643236271438044_DP ;  casimir_omega_weight(132) = 0.0051419764405473_DP
         casimir_omega(133) = 0.2695061359467907_DP ;  casimir_omega_weight(133) = 0.0052232725569754_DP
         casimir_omega(134) = 0.2747706405130398_DP ;  casimir_omega_weight(134) = 0.0053059734690541_DP
         casimir_omega(135) = 0.2801185622945875_DP ;  casimir_omega_weight(135) = 0.0053901126409725_DP
         casimir_omega(136) = 0.2855513566668506_DP ;  casimir_omega_weight(136) = 0.0054757244602786_DP
         casimir_omega(137) = 0.2910705138670975_DP ;  casimir_omega_weight(137) = 0.0055628442682845_DP
         casimir_omega(138) = 0.2966775599638825_DP ;  casimir_omega_weight(138) = 0.0056515083916069_DP
         casimir_omega(139) = 0.3023740578586144_DP ;  casimir_omega_weight(139) = 0.0057417541748979_DP
         casimir_omega(140) = 0.3081616083204732_DP ;  casimir_omega_weight(140) = 0.0058336200148127_DP
         casimir_omega(141) = 0.3140418510559377_DP ;  casimir_omega_weight(141) = 0.0059271453952798_DP
         casimir_omega(142) = 0.3200164658142449_DP ;  casimir_omega_weight(142) = 0.0060223709241084_DP
         casimir_omega(143) = 0.3260871735301661_DP ;  casimir_omega_weight(143) = 0.0061193383710119_DP
         casimir_omega(144) = 0.3322557375055268_DP ;  casimir_omega_weight(144) = 0.0062180907070956_DP
         casimir_omega(145) = 0.3385239646309793_DP ;  casimir_omega_weight(145) = 0.0063186721458735_DP
         casimir_omega(146) = 0.3448937066495880_DP ;  casimir_omega_weight(146) = 0.0064211281858986_DP
         casimir_omega(147) = 0.3513668614638709_DP ;  casimir_omega_weight(147) = 0.0065255056550468_DP
         casimir_omega(148) = 0.3579453744879942_DP ;  casimir_omega_weight(148) = 0.0066318527565634_DP
         casimir_omega(149) = 0.3646312400469140_DP ;  casimir_omega_weight(149) = 0.0067402191169174_DP
         casimir_omega(150) = 0.3714265028243223_DP ;  casimir_omega_weight(150) = 0.0068506558355689_DP
         casimir_omega(151) = 0.3783332593613498_DP ;  casimir_omega_weight(151) = 0.0069632155367279_DP
         casimir_omega(152) = 0.3853536596080546_DP ;  casimir_omega_weight(152) = 0.0070779524231934_DP
         casimir_omega(153) = 0.3924899085298298_DP ;  casimir_omega_weight(153) = 0.0071949223323590_DP
         casimir_omega(154) = 0.3997442677709550_DP ;  casimir_omega_weight(154) = 0.0073141827945131_DP
         casimir_omega(155) = 0.4071190573776168_DP ;  casimir_omega_weight(155) = 0.0074357930934969_DP
         casimir_omega(156) = 0.4146166575828314_DP ;  casimir_omega_weight(156) = 0.0075598143298690_DP
         casimir_omega(157) = 0.4222395106558230_DP ;  casimir_omega_weight(157) = 0.0076863094866678_DP
         casimir_omega(158) = 0.4299901228185127_DP ;  casimir_omega_weight(158) = 0.0078153434979023_DP
         casimir_omega(159) = 0.4378710662319235_DP ;  casimir_omega_weight(159) = 0.0079469833199004_DP
         casimir_omega(160) = 0.4458849810554083_DP ;  casimir_omega_weight(160) = 0.0080812980056419_DP
         casimir_omega(161) = 0.4540345775817753_DP ;  casimir_omega_weight(161) = 0.0082183587822353_DP
         casimir_omega(162) = 0.4623226384515043_DP ;  casimir_omega_weight(162) = 0.0083582391316579_DP
         casimir_omega(163) = 0.4707520209494217_DP ;  casimir_omega_weight(163) = 0.0085010148749517_DP
         casimir_omega(164) = 0.4793256593873524_DP ;  casimir_omega_weight(164) = 0.0086467642600123_DP
         casimir_omega(165) = 0.4880465675764384_DP ;  casimir_omega_weight(165) = 0.0087955680531646_DP
         casimir_omega(166) = 0.4969178413929957_DP ;  casimir_omega_weight(166) = 0.0089475096346954_DP
         casimir_omega(167) = 0.5059426614419690_DP ;  casimir_omega_weight(167) = 0.0091026750985491_DP
         casimir_omega(168) = 0.5151242958222415_DP ;  casimir_omega_weight(168) = 0.0092611533563743_DP
         casimir_omega(169) = 0.5244661029982713_DP ;  casimir_omega_weight(169) = 0.0094230362461679_DP
         casimir_omega(170) = 0.5339715347827410_DP ;  casimir_omega_weight(170) = 0.0095884186457080_DP
         casimir_omega(171) = 0.5436441394351557_DP ;  casimir_omega_weight(171) = 0.0097573985910343_DP
         casimir_omega(172) = 0.5534875648815430_DP ;  casimir_omega_weight(172) = 0.0099300774002474_DP
         casimir_omega(173) = 0.5635055620607162_DP ;  casimir_omega_weight(173) = 0.0101065598028510_DP
         casimir_omega(174) = 0.5737019884027811_DP ;  casimir_omega_weight(174) = 0.0102869540749684_DP
         casimir_omega(175) = 0.5840808114459128_DP ;  casimir_omega_weight(175) = 0.0104713721806999_DP
         casimir_omega(176) = 0.5946461125976923_DP ;  casimir_omega_weight(176) = 0.0106599299199539_DP
         casimir_omega(177) = 0.6054020910476524_DP ;  casimir_omega_weight(177) = 0.0108527470830825_DP
         casimir_omega(178) = 0.6163530678380068_DP ;  casimir_omega_weight(178) = 0.0110499476126783_DP
         casimir_omega(179) = 0.6275034900999058_DP ;  casimir_omega_weight(179) = 0.0112516597728932_DP
         casimir_omega(180) = 0.6388579354629543_DP ;  casimir_omega_weight(180) = 0.0114580163267082_DP
         casimir_omega(181) = 0.6504211166461307_DP ;  casimir_omega_weight(181) = 0.0116691547215411_DP
         casimir_omega(182) = 0.6621978862386698_DP ;  casimir_omega_weight(182) = 0.0118852172836533_DP
         casimir_omega(183) = 0.6741932416799613_DP ;  casimir_omega_weight(183) = 0.0121063514218291_DP
         casimir_omega(184) = 0.6864123304479550_DP ;  casimir_omega_weight(184) = 0.0123327098408292_DP
         casimir_omega(185) = 0.6988604554661285_DP ;  casimir_omega_weight(185) = 0.0125644507651234_DP
         casimir_omega(186) = 0.7115430807395792_DP ;  casimir_omega_weight(186) = 0.0128017381735173_DP
         casimir_omega(187) = 0.7244658372314071_DP ;  casimir_omega_weight(187) = 0.0130447420452129_DP
         casimir_omega(188) = 0.7376345289911638_DP ;  casimir_omega_weight(188) = 0.0132936386179791_DP
         casimir_omega(189) = 0.7510551395477812_DP ;  casimir_omega_weight(189) = 0.0135486106590687_DP
         casimir_omega(190) = 0.7647338385801196_DP ;  casimir_omega_weight(190) = 0.0138098477496422_DP
         casimir_omega(191) = 0.7786769888789744_DP ;  casimir_omega_weight(191) = 0.0140775465834171_DP
         casimir_omega(192) = 0.7928911536151911_DP ;  casimir_omega_weight(192) = 0.0143519112803609_DP
         casimir_omega(193) = 0.8073831039293619_DP ;  casimir_omega_weight(193) = 0.0146331537163009_DP
         casimir_omega(194) = 0.8221598268594480_DP ;  casimir_omega_weight(194) = 0.0149214938693470_DP
         casimir_omega(195) = 0.8372285336236571_DP ;  casimir_omega_weight(195) = 0.0152171601841095_DP
         casimir_omega(196) = 0.8525966682768449_DP ;  casimir_omega_weight(196) = 0.0155203899547422_DP
         casimir_omega(197) = 0.8682719167598315_DP ;  casimir_omega_weight(197) = 0.0158314297279111_DP
         casimir_omega(198) = 0.8842622163621490_DP ;  casimir_omega_weight(198) = 0.0161505357268678_DP
         casimir_omega(199) = 0.9005757656199157_DP ;  casimir_omega_weight(199) = 0.0164779742978750_DP
         casimir_omega(200) = 0.9172210346718754_DP ;  casimir_omega_weight(200) = 0.0168140223803197_DP
         casimir_omega(201) = 0.9342067760979826_DP ;  casimir_omega_weight(201) = 0.0171589680019423_DP
         casimir_omega(202) = 0.9515420362663916_DP ;  casimir_omega_weight(202) = 0.0175131108006868_DP
         casimir_omega(203) = 0.9692361672162983_DP ;  casimir_omega_weight(203) = 0.0178767625748141_DP
         casimir_omega(204) = 0.9872988391057267_DP ;  casimir_omega_weight(204) = 0.0182502478630040_DP
         casimir_omega(205) = 1.0057400532551781_DP ;  casimir_omega_weight(205) = 0.0186339045562653_DP
         casimir_omega(206) = 1.0245701558199349_DP ;  casimir_omega_weight(206) = 0.0190280845437193_DP
         casimir_omega(207) = 1.0437998521259171_DP ;  casimir_omega_weight(207) = 0.0194331543942378_DP
         casimir_omega(208) = 1.0634402217061092_DP ;  casimir_omega_weight(208) = 0.0198494960763571_DP
         casimir_omega(209) = 1.0835027340769725_DP ;  casimir_omega_weight(209) = 0.0202775077187790_DP
         casimir_omega(210) = 1.1039992652967578_DP ;  casimir_omega_weight(210) = 0.0207176044140676_DP
         casimir_omega(211) = 1.1249421153502881_DP ;  casimir_omega_weight(211) = 0.0211702190684229_DP
         casimir_omega(212) = 1.1463440264077289_DP ;  casimir_omega_weight(212) = 0.0216358033003645_DP
         casimir_omega(213) = 1.1682182020078791_DP ;  casimir_omega_weight(213) = 0.0221148283916460_DP
         casimir_omega(214) = 1.1905783272198909_DP ;  casimir_omega_weight(214) = 0.0226077862938014_DP
         casimir_omega(215) = 1.2134385898408577_DP ;  casimir_omega_weight(215) = 0.0231151906939690_DP
         casimir_omega(216) = 1.2368137026905288_DP ;  casimir_omega_weight(216) = 0.0236375781440006_DP
         casimir_omega(217) = 1.2607189270685475_DP ;  casimir_omega_weight(217) = 0.0241755092571299_DP
         casimir_omega(218) = 1.2851700974439717_DP ;  casimir_omega_weight(218) = 0.0247295699766826_DP
         casimir_omega(219) = 1.3101836474516482_DP ;  casimir_omega_weight(219) = 0.0253003729218946_DP
         casimir_omega(220) = 1.3357766372750648_DP ;  casimir_omega_weight(220) = 0.0258885588159896_DP
         casimir_omega(221) = 1.3619667825008477_DP ;  casimir_omega_weight(221) = 0.0264947980023609_DP
         casimir_omega(222) = 1.3887724845359777_DP ;  casimir_omega_weight(222) = 0.0271197920549256_DP
         casimir_omega(223) = 1.4162128626852239_DP ;  casimir_omega_weight(223) = 0.0277642754892971_DP
         casimir_omega(224) = 1.4443077879931132_DP ;  casimir_omega_weight(224) = 0.0284290175819873_DP
         casimir_omega(225) = 1.4730779189623020_DP ;  casimir_omega_weight(225) = 0.0291148243052877_DP
         casimir_omega(226) = 1.5025447392681406_DP ;  casimir_omega_weight(226) = 0.0298225403862683_DP
         casimir_omega(227) = 1.5327305975979884_DP ;  casimir_omega_weight(227) = 0.0305530514988201_DP
         casimir_omega(228) = 1.5636587497531667_DP ;  casimir_omega_weight(228) = 0.0313072865985656_DP
         casimir_omega(229) = 1.5953534031615981_DP ;  casimir_omega_weight(229) = 0.0320862204111779_DP
         casimir_omega(230) = 1.6278397639602133_DP ;  casimir_omega_weight(230) = 0.0328908760855582_DP
         casimir_omega(231) = 1.6611440868180236_DP ;  casimir_omega_weight(231) = 0.0337223280242012_DP
         casimir_omega(232) = 1.6952937276837849_DP ;  casimir_omega_weight(232) = 0.0345817049042018_DP
         casimir_omega(233) = 1.7303171996559936_DP ;  casimir_omega_weight(233) = 0.0354701929035808_DP
         casimir_omega(234) = 1.7662442321883394_DP ;  casimir_omega_weight(234) = 0.0363890391485776_DP
         casimir_omega(235) = 1.8031058338600903_DP ;  casimir_omega_weight(235) = 0.0373395553991747_DP
         casimir_omega(236) = 1.8409343589588338_DP ;  casimir_omega_weight(236) = 0.0383231219915569_DP
         casimir_omega(237) = 1.8797635781425537_DP ;  casimir_omega_weight(237) = 0.0393411920578561_DP
         casimir_omega(238) = 1.9196287534691425_DP ;  casimir_omega_weight(238) = 0.0403952960451476_DP
         casimir_omega(239) = 1.9605667181046331_DP ;  casimir_omega_weight(239) = 0.0414870465579960_DP
         casimir_omega(240) = 2.0026159610465393_DP ;  casimir_omega_weight(240) = 0.0426181435507933_DP
         casimir_omega(241) = 2.0458167172261903_DP ;  casimir_omega_weight(241) = 0.0437903798984394_DP
         casimir_omega(242) = 2.0902110633839981_DP ;  casimir_omega_weight(242) = 0.0450056473767723_DP
         casimir_omega(243) = 2.1358430201440819_DP ;  casimir_omega_weight(243) = 0.0462659430870895_DP
         casimir_omega(244) = 2.1827586607508738_DP ;  casimir_omega_weight(244) = 0.0475733763618885_DP
         casimir_omega(245) = 2.2310062269690865_DP ;  casimir_omega_weight(245) = 0.0489301761930932_DP
         casimir_omega(246) = 2.2806362526915835_DP ;  casimir_omega_weight(246) = 0.0503386992274388_DP
         casimir_omega(247) = 2.3317016958465615_DP ;  casimir_omega_weight(247) = 0.0518014383783194_DP
         casimir_omega(248) = 2.3842580792470742_DP ;  casimir_omega_weight(248) = 0.0533210321080033_DP
         casimir_omega(249) = 2.4383636410824825_DP ;  casimir_omega_weight(249) = 0.0549002744396054_DP
         casimir_omega(250) = 2.4940794958135277_DP ;  casimir_omega_weight(250) = 0.0565421257638824_DP
         casimir_omega(251) = 2.5514698063013284_DP ;  casimir_omega_weight(251) = 0.0582497245126388_DP
         casimir_omega(252) = 2.6106019680755774_DP ;  casimir_omega_weight(252) = 0.0600263997776282_DP
         casimir_omega(253) = 2.6715468067303756_DP ;  casimir_omega_weight(253) = 0.0618756849621516_DP
         casimir_omega(254) = 2.7343787895275158_DP ;  casimir_omega_weight(254) = 0.0638013325611125_DP
         casimir_omega(255) = 2.7991762523879427_DP ;  casimir_omega_weight(255) = 0.0658073301759836_DP
         casimir_omega(256) = 2.8660216435637955_DP ;  casimir_omega_weight(256) = 0.0678979178814384_DP
         casimir_omega(257) = 2.9350017854066732_DP ;  casimir_omega_weight(257) = 0.0700776070739542_DP
         casimir_omega(258) = 3.0062081557845830_DP ;  casimir_omega_weight(258) = 0.0723512009458147_DP
         casimir_omega(259) = 3.0797371908512585_DP ;  casimir_omega_weight(259) = 0.0747238167437983_DP
         casimir_omega(260) = 3.1556906110396921_DP ;  casimir_omega_weight(260) = 0.0772009099899889_DP
         casimir_omega(261) = 3.2341757723383933_DP ;  casimir_omega_weight(261) = 0.0797883008609497_DP
         casimir_omega(262) = 3.3153060451163201_DP ;  casimir_omega_weight(262) = 0.0824922029442865_DP
         casimir_omega(263) = 3.3992012229935002_DP ;  casimir_omega_weight(263) = 0.0853192546165616_DP
         casimir_omega(264) = 3.4859879645121938_DP ;  casimir_omega_weight(264) = 0.0882765533139914_DP
         casimir_omega(265) = 3.5758002706503356_DP ;  casimir_omega_weight(265) = 0.0913716930000738_DP
         casimir_omega(266) = 3.6687800015406271_DP ;  casimir_omega_weight(266) = 0.0946128051688019_DP
         casimir_omega(267) = 3.7650774361175352_DP ;  casimir_omega_weight(267) = 0.0980086037642671_DP
         casimir_omega(268) = 3.8648518788169874_DP ;  casimir_omega_weight(268) = 0.1015684344416592_DP
         casimir_omega(269) = 3.9682723179047130_DP ;  casimir_omega_weight(269) = 0.1053023286478987_DP
         casimir_omega(270) = 4.0755181405161567_DP ;  casimir_omega_weight(270) = 0.1092210630588720_DP
         casimir_omega(271) = 4.1867799100607774_DP ;  casimir_omega_weight(271) = 0.1133362249777034_DP
         casimir_omega(272) = 4.3022602122850753_DP ;  casimir_omega_weight(272) = 0.1176602843748558_DP
         casimir_omega(273) = 4.4221745770134708_DP ;  casimir_omega_weight(273) = 0.1222066733391986_DP
         casimir_omega(274) = 4.5467524834026571_DP ;  casimir_omega_weight(274) = 0.1269898738080190_DP
         casimir_omega(275) = 4.6762384574704399_DP ;  casimir_omega_weight(275) = 0.1320255145604900_DP
         casimir_omega(276) = 4.8108932717079140_DP ;  casimir_omega_weight(276) = 0.1373304785881928_DP
         casimir_omega(277) = 4.9509952577718241_DP ;  casimir_omega_weight(277) = 0.1429230221100386_DP
         casimir_omega(278) = 5.0968417446056771_DP ;  casimir_omega_weight(278) = 0.1488229066707422_DP
         casimir_omega(279) = 5.2487506358758429_DP ;  casimir_omega_weight(279) = 0.1550515459637518_DP
         casimir_omega(280) = 5.4070621423622862_DP ;  casimir_omega_weight(280) = 0.1616321692511106_DP
         casimir_omega(281) = 5.5721406869476819_DP ;  casimir_omega_weight(281) = 0.1685900035214046_DP
         casimir_omega(282) = 5.7443770021403981_DP ;  casimir_omega_weight(282) = 0.1759524768371455_DP
         casimir_omega(283) = 5.9241904426958998_DP ;  casimir_omega_weight(283) = 0.1837494456859623_DP
         casimir_omega(284) = 6.1120315389201094_DP ;  casimir_omega_weight(284) = 0.1920134495709994_DP
         casimir_omega(285) = 6.3083848197125727_DP ;  casimir_omega_weight(285) = 0.2007799965673394_DP
         casimir_omega(286) = 6.5137719384158004_DP ;  casimir_omega_weight(286) = 0.2100878841469646_DP
         casimir_omega(287) = 6.7287551391662097_DP ;  casimir_omega_weight(287) = 0.2199795602491640_DP
         casimir_omega(288) = 6.9539411068063020_DP ;  casimir_omega_weight(288) = 0.2305015303673108_DP
         casimir_omega(289) = 7.1899852496411762_DP ;  casimir_omega_weight(289) = 0.2417048173542562_DP
         casimir_omega(290) = 7.4375964715620491_DP ;  casimir_omega_weight(290) = 0.2536454817542954_DP
         casimir_omega(291) = 7.6975424985015950_DP ;  casimir_omega_weight(291) = 0.2663852117742796_DP
         casimir_omega(292) = 7.9706558340507989_DP ;  casimir_omega_weight(292) = 0.2799919935569667_DP
         casimir_omega(293) = 8.2578404306259046_DP ;  casimir_omega_weight(293) = 0.2945408742645161_DP
         casimir_omega(294) = 8.5600791761521950_DP ;  casimir_omega_weight(294) = 0.3101148326844788_DP
         casimir_omega(295) = 8.8784423122195975_DP ;  casimir_omega_weight(295) = 0.3268057747066686_DP
         casimir_omega(296) = 9.2140969185513519_DP ;  casimir_omega_weight(296) = 0.3447156741863000_DP
         casimir_omega(297) = 9.5683176209876759_DP ;  casimir_omega_weight(297) = 0.3639578835238507_DP
         casimir_omega(298) = 9.9424987067493706_DP ;  casimir_omega_weight(298) = 0.3846586429014390_DP
         casimir_omega(299) = 10.3381678623903710_DP ;  casimir_omega_weight(299) = 0.4069588227030156_DP
         casimir_omega(300) = 10.7570017876611743_DP ;  casimir_omega_weight(300) = 0.4310159404464045_DP
         casimir_omega(301) = 11.2008439838471112_DP ;  casimir_omega_weight(301) = 0.4570065018548643_DP
         casimir_omega(302) = 11.6717250696827382_DP ;  casimir_omega_weight(302) = 0.4851287258681803_DP
         casimir_omega(303) = 12.1718860437853991_DP ;  casimir_omega_weight(303) = 0.5156057259108426_DP
         casimir_omega(304) = 12.7038049923198795_DP ;  casimir_omega_weight(304) = 0.5486892351831218_DP
         casimir_omega(305) = 13.2702278376116265_DP ;  casimir_omega_weight(305) = 0.5846639829207773_DP
         casimir_omega(306) = 13.8742038418695621_DP ;  casimir_omega_weight(306) = 0.6238528524377512_DP
         casimir_omega(307) = 14.5191267253601950_DP ;  casimir_omega_weight(307) = 0.6666229816367034_DP
         casimir_omega(308) = 15.2087824371303846_DP ;  casimir_omega_weight(308) = 0.7133930042435774_DP
         casimir_omega(309) = 15.9474048373944033_DP ;  casimir_omega_weight(309) = 0.7646416774285438_DP
         casimir_omega(310) = 16.7397408252919817_DP ;  casimir_omega_weight(310) = 0.8209182017369189_DP
         casimir_omega(311) = 17.5911267884841678_DP ;  casimir_omega_weight(311) = 0.8828546160961225_DP
         casimir_omega(312) = 18.5075786811006715_DP ;  casimir_omega_weight(312) = 0.9511807493629741_DP
         casimir_omega(313) = 19.4958985789479584_DP ;  casimir_omega_weight(313) = 1.0267423372314695_DP
         casimir_omega(314) = 20.5638012488204716_DP ;  casimir_omega_weight(314) = 1.1105230788136740_DP
         casimir_omega(315) = 21.7200651463087162_DP ;  casimir_omega_weight(315) = 1.2036716234682894_DP
         casimir_omega(316) = 22.9747133828461330_DP ;  casimir_omega_weight(316) = 1.3075347631562619_DP
         casimir_omega(317) = 24.3392316576871472_DP ;  casimir_omega_weight(317) = 1.4236984828680044_DP
         casimir_omega(318) = 25.8268320426897766_DP ;  casimir_omega_weight(318) = 1.5540390254898628_DP
         casimir_omega(319) = 27.4527739861790856_DP ;  casimir_omega_weight(319) = 1.7007868055275850_DP
         casimir_omega(320) = 29.2347571729786111_DP ;  casimir_omega_weight(320) = 1.8666069262410270_DP
         casimir_omega(321) = 31.1934052286917449_DP ;  casimir_omega_weight(321) = 2.0547013144098756_DP
         casimir_omega(322) = 33.3528650933646631_DP ;  casimir_omega_weight(322) = 2.2689392275875644_DP
         casimir_omega(323) = 35.7415547906808015_DP ;  casimir_omega_weight(323) = 2.5140253176765976_DP
         casimir_omega(324) = 38.3931031159444913_DP ;  casimir_omega_weight(324) = 2.7957178600671888_DP
         casimir_omega(325) = 41.3475396713607140_DP ;  casimir_omega_weight(325) = 3.1211146424257032_DP
         casimir_omega(326) = 44.6528144775032771_DP ;  casimir_omega_weight(326) = 3.4990310570050474_DP
         casimir_omega(327) = 48.3667557555708427_DP ;  casimir_omega_weight(327) = 3.9405052452548368_DP
         casimir_omega(328) = 52.5596164514691750_DP ;  casimir_omega_weight(328) = 4.4594804148839087_DP
         casimir_omega(329) = 57.3174208834503460_DP ;  casimir_omega_weight(329) = 5.0737374165099682_DP
         casimir_omega(330) = 62.7464122703799987_DP ;  casimir_omega_weight(330) = 5.8061857622733921_DP
         casimir_omega(331) = 68.9790353209870943_DP ;  casimir_omega_weight(331) = 6.6866758276810527_DP
         casimir_omega(332) = 76.1820906474642783_DP ;  casimir_omega_weight(332) = 7.7545813959917043_DP
         casimir_omega(333) = 84.5680110873966697_DP ;  casimir_omega_weight(333) = 9.0625414022568869_DP
         casimir_omega(334) = 94.4107044908105451_DP ;  casimir_omega_weight(334) = 10.6819806526328769_DP
         casimir_omega(335) = 106.0682053929857034_DP ;  casimir_omega_weight(335) = 12.7114204677481233_DP
         casimir_omega(336) = 120.0156973143437398_DP ;  casimir_omega_weight(336) = 15.2892710223059769_DP
         casimir_omega(337) = 136.8947091324740200_DP ;  casimir_omega_weight(337) = 18.6140182655159592_DP
         casimir_omega(338) = 157.5882153516826349_DP ;  casimir_omega_weight(338) = 22.9769833031473745_DP
         casimir_omega(339) = 183.3384860994507903_DP ;  casimir_omega_weight(339) = 28.8171952761367116_DP
         casimir_omega(340) = 215.9379383924104445_DP ;  casimir_omega_weight(340) = 36.8166920790265024_DP
         casimir_omega(341) = 258.0496382925509238_DP ;  casimir_omega_weight(341) = 48.0730925733503724_DP
         casimir_omega(342) = 313.7688158874943838_DP ;  casimir_omega_weight(342) = 64.4277078062844453_DP
         casimir_omega(343) = 389.6571545559793890_DP ;  casimir_omega_weight(343) = 89.1264511113906082_DP
         casimir_omega(344) = 496.7661304501882569_DP ;  casimir_omega_weight(344) = 128.2469316773791661_DP
         casimir_omega(345) = 654.8986115543225424_DP ;  casimir_omega_weight(345) = 194.0549227636805369_DP
         casimir_omega(346) = 902.4594484635399567_DP ;  casimir_omega_weight(346) = 313.7984863538891318_DP
         casimir_omega(347) = 1322.1575624763538599_DP ;  casimir_omega_weight(347) = 556.2450510219563284_DP
         casimir_omega(348) = 2120.1422788421778023_DP ;  casimir_omega_weight(348) = 1128.9389005772484325_DP
         casimir_omega(349) = 3936.7422774179403859_DP ;  casimir_omega_weight(349) = 2853.9579872466042616_DP
         casimir_omega(350) = 9675.6284084523686033_DP ;  casimir_omega_weight(350) = 10970.6564337799190980_DP
         casimir_omega(351) = 50982.0172853888652753_DP ;  casimir_omega_weight(351) = 130837.3464822293171892_DP
return
endsubroutine gauss_legendre_grid350

subroutine gauss_legendre_grid400()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000054082375593_DP ;  casimir_omega_weight(2) = 0.0000138793674598_DP
         casimir_omega(3) = 0.0000284964125459_DP ;  casimir_omega_weight(3) = 0.0000323101831690_DP
         casimir_omega(4) = 0.0000700367942871_DP ;  casimir_omega_weight(4) = 0.0000507722630229_DP
         casimir_omega(5) = 0.0001300437207159_DP ;  casimir_omega_weight(5) = 0.0000692429622521_DP
         casimir_omega(6) = 0.0002085257083295_DP ;  casimir_omega_weight(6) = 0.0000877226614055_DP
         casimir_omega(7) = 0.0003054926893935_DP ;  casimir_omega_weight(7) = 0.0001062133042517_DP
         casimir_omega(8) = 0.0004209566946811_DP ;  casimir_omega_weight(8) = 0.0001247170823377_DP
         casimir_omega(9) = 0.0005549319833652_DP ;  casimir_omega_weight(9) = 0.0001432362489622_DP
         casimir_omega(10) = 0.0007074350804268_DP ;  casimir_omega_weight(10) = 0.0001617730794006_DP
         casimir_omega(11) = 0.0008784847921827_DP ;  casimir_omega_weight(11) = 0.0001803298601466_DP
         casimir_omega(12) = 0.0010681022155700_DP ;  casimir_omega_weight(12) = 0.0001989088856761_DP
         casimir_omega(13) = 0.0012763107455774_DP ;  casimir_omega_weight(13) = 0.0002175124575674_DP
         casimir_omega(14) = 0.0015031360822950_DP ;  casimir_omega_weight(14) = 0.0002361428844755_DP
         casimir_omega(15) = 0.0017486062381279_DP ;  casimir_omega_weight(15) = 0.0002548024824439_DP
         casimir_omega(16) = 0.0020127515454111_DP ;  casimir_omega_weight(16) = 0.0002734935753753_DP
         casimir_omega(17) = 0.0022956046645309_DP ;  casimir_omega_weight(17) = 0.0002922184955637_DP
         casimir_omega(18) = 0.0025972005926081_DP ;  casimir_omega_weight(18) = 0.0003109795842804_DP
         casimir_omega(19) = 0.0029175766727700_DP ;  casimir_omega_weight(19) = 0.0003297791923711_DP
         casimir_omega(20) = 0.0032567726040326_DP ;  casimir_omega_weight(20) = 0.0003486196808690_DP
         casimir_omega(21) = 0.0036148304518003_DP ;  casimir_omega_weight(21) = 0.0003675034216244_DP
         casimir_omega(22) = 0.0039917946589982_DP ;  casimir_omega_weight(22) = 0.0003864327979334_DP
         casimir_omega(23) = 0.0043877120578397_DP ;  casimir_omega_weight(23) = 0.0004054102051841_DP
         casimir_omega(24) = 0.0048026318822385_DP ;  casimir_omega_weight(24) = 0.0004244380514977_DP
         casimir_omega(25) = 0.0052366057808717_DP ;  casimir_omega_weight(25) = 0.0004435187583892_DP
         casimir_omega(26) = 0.0056896878309007_DP ;  casimir_omega_weight(26) = 0.0004626547614230_DP
         casimir_omega(27) = 0.0061619345523552_DP ;  casimir_omega_weight(27) = 0.0004818485108835_DP
         casimir_omega(28) = 0.0066534049231904_DP ;  casimir_omega_weight(28) = 0.0005011024724488_DP
         casimir_omega(29) = 0.0071641603950196_DP ;  casimir_omega_weight(29) = 0.0005204191278731_DP
         casimir_omega(30) = 0.0076942649095350_DP ;  casimir_omega_weight(30) = 0.0005398009756745_DP
         casimir_omega(31) = 0.0082437849156198_DP ;  casimir_omega_weight(31) = 0.0005592505318344_DP
         casimir_omega(32) = 0.0088127893871627_DP ;  casimir_omega_weight(32) = 0.0005787703305008_DP
         casimir_omega(33) = 0.0094013498415818_DP ;  casimir_omega_weight(33) = 0.0005983629247023_DP
         casimir_omega(34) = 0.0100095403590651_DP ;  casimir_omega_weight(34) = 0.0006180308870743_DP
         casimir_omega(35) = 0.0106374376025414_DP ;  casimir_omega_weight(35) = 0.0006377768105833_DP
         casimir_omega(36) = 0.0112851208383849_DP ;  casimir_omega_weight(36) = 0.0006576033092779_DP
         casimir_omega(37) = 0.0119526719578655_DP ;  casimir_omega_weight(37) = 0.0006775130190294_DP
         casimir_omega(38) = 0.0126401754993577_DP ;  casimir_omega_weight(38) = 0.0006975085983017_DP
         casimir_omega(39) = 0.0133477186713158_DP ;  casimir_omega_weight(39) = 0.0007175927289196_DP
         casimir_omega(40) = 0.0140753913760237_DP ;  casimir_omega_weight(40) = 0.0007377681168501_DP
         casimir_omega(41) = 0.0148232862341372_DP ;  casimir_omega_weight(41) = 0.0007580374930021_DP
         casimir_omega(42) = 0.0155914986100244_DP ;  casimir_omega_weight(42) = 0.0007784036140282_DP
         casimir_omega(43) = 0.0163801266379181_DP ;  casimir_omega_weight(43) = 0.0007988692631419_DP
         casimir_omega(44) = 0.0171892712488948_DP ;  casimir_omega_weight(44) = 0.0008194372509566_DP
         casimir_omega(45) = 0.0180190361986915_DP ;  casimir_omega_weight(45) = 0.0008401104163215_DP
         casimir_omega(46) = 0.0188695280963698_DP ;  casimir_omega_weight(46) = 0.0008608916271819_DP
         casimir_omega(47) = 0.0197408564338500_DP ;  casimir_omega_weight(47) = 0.0008817837814540_DP
         casimir_omega(48) = 0.0206331336163190_DP ;  casimir_omega_weight(48) = 0.0009027898079078_DP
         casimir_omega(49) = 0.0215464749935391_DP ;  casimir_omega_weight(49) = 0.0009239126670708_DP
         casimir_omega(50) = 0.0224809988920567_DP ;  casimir_omega_weight(50) = 0.0009451553521418_DP
         casimir_omega(51) = 0.0234368266483427_DP ;  casimir_omega_weight(51) = 0.0009665208899260_DP
         casimir_omega(52) = 0.0244140826428721_DP ;  casimir_omega_weight(52) = 0.0009880123417844_DP
         casimir_omega(53) = 0.0254128943351581_DP ;  casimir_omega_weight(53) = 0.0010096328045980_DP
         casimir_omega(54) = 0.0264333922997643_DP ;  casimir_omega_weight(54) = 0.0010313854117510_DP
         casimir_omega(55) = 0.0274757102633091_DP ;  casimir_omega_weight(55) = 0.0010532733341359_DP
         casimir_omega(56) = 0.0285399851424800_DP ;  casimir_omega_weight(56) = 0.0010752997811708_DP
         casimir_omega(57) = 0.0296263570830768_DP ;  casimir_omega_weight(57) = 0.0010974680018383_DP
         casimir_omega(58) = 0.0307349695001068_DP ;  casimir_omega_weight(58) = 0.0011197812857471_DP
         casimir_omega(59) = 0.0318659691189473_DP ;  casimir_omega_weight(59) = 0.0011422429642114_DP
         casimir_omega(60) = 0.0330195060175991_DP ;  casimir_omega_weight(60) = 0.0011648564113487_DP
         casimir_omega(61) = 0.0341957336700505_DP ;  casimir_omega_weight(61) = 0.0011876250452057_DP
         casimir_omega(62) = 0.0353948089907768_DP ;  casimir_omega_weight(62) = 0.0012105523289017_DP
         casimir_omega(63) = 0.0366168923803931_DP ;  casimir_omega_weight(63) = 0.0012336417717957_DP
         casimir_omega(64) = 0.0378621477724896_DP ;  casimir_omega_weight(64) = 0.0012568969306766_DP
         casimir_omega(65) = 0.0391307426816674_DP ;  casimir_omega_weight(65) = 0.0012803214109834_DP
         casimir_omega(66) = 0.0404228482528061_DP ;  casimir_omega_weight(66) = 0.0013039188680434_DP
         casimir_omega(67) = 0.0417386393115823_DP ;  casimir_omega_weight(67) = 0.0013276930083359_DP
         casimir_omega(68) = 0.0430782944162678_DP ;  casimir_omega_weight(68) = 0.0013516475907909_DP
         casimir_omega(69) = 0.0444419959108396_DP ;  casimir_omega_weight(69) = 0.0013757864281097_DP
         casimir_omega(70) = 0.0458299299794201_DP ;  casimir_omega_weight(70) = 0.0014001133881099_DP
         casimir_omega(71) = 0.0472422867020820_DP ;  casimir_omega_weight(71) = 0.0014246323951078_DP
         casimir_omega(72) = 0.0486792601120497_DP ;  casimir_omega_weight(72) = 0.0014493474313257_DP
         casimir_omega(73) = 0.0501410482543200_DP ;  casimir_omega_weight(73) = 0.0014742625383307_DP
         casimir_omega(74) = 0.0516278532457425_DP ;  casimir_omega_weight(74) = 0.0014993818185069_DP
         casimir_omega(75) = 0.0531398813365825_DP ;  casimir_omega_weight(75) = 0.0015247094365575_DP
         casimir_omega(76) = 0.0546773429736096_DP ;  casimir_omega_weight(76) = 0.0015502496210454_DP
         casimir_omega(77) = 0.0562404528647376_DP ;  casimir_omega_weight(77) = 0.0015760066659615_DP
         casimir_omega(78) = 0.0578294300452554_DP ;  casimir_omega_weight(78) = 0.0016019849323365_DP
         casimir_omega(79) = 0.0594444979456864_DP ;  casimir_omega_weight(79) = 0.0016281888498839_DP
         casimir_omega(80) = 0.0610858844613069_DP ;  casimir_omega_weight(80) = 0.0016546229186838_DP
         casimir_omega(81) = 0.0627538220233728_DP ;  casimir_omega_weight(81) = 0.0016812917109021_DP
         casimir_omega(82) = 0.0644485476720827_DP ;  casimir_omega_weight(82) = 0.0017081998725576_DP
         casimir_omega(83) = 0.0661703031313290_DP ;  casimir_omega_weight(83) = 0.0017353521253198_DP
         casimir_omega(84) = 0.0679193348852661_DP ;  casimir_omega_weight(84) = 0.0017627532683590_DP
         casimir_omega(85) = 0.0696958942567557_DP ;  casimir_omega_weight(85) = 0.0017904081802352_DP
         casimir_omega(86) = 0.0715002374877172_DP ;  casimir_omega_weight(86) = 0.0018183218208306_DP
         casimir_omega(87) = 0.0733326258214392_DP ;  casimir_omega_weight(87) = 0.0018464992333362_DP
         casimir_omega(88) = 0.0751933255868968_DP ;  casimir_omega_weight(88) = 0.0018749455462770_DP
         casimir_omega(89) = 0.0770826082851235_DP ;  casimir_omega_weight(89) = 0.0019036659755935_DP
         casimir_omega(90) = 0.0790007506776846_DP ;  casimir_omega_weight(90) = 0.0019326658267700_DP
         casimir_omega(91) = 0.0809480348773099_DP ;  casimir_omega_weight(91) = 0.0019619504970177_DP
         casimir_omega(92) = 0.0829247484407331_DP ;  casimir_omega_weight(92) = 0.0019915254775115_DP
         casimir_omega(93) = 0.0849311844637965_DP ;  casimir_omega_weight(93) = 0.0020213963556811_DP
         casimir_omega(94) = 0.0869676416788746_DP ;  casimir_omega_weight(94) = 0.0020515688175570_DP
         casimir_omega(95) = 0.0890344245546779_DP ;  casimir_omega_weight(95) = 0.0020820486501853_DP
         casimir_omega(96) = 0.0911318433984924_DP ;  casimir_omega_weight(96) = 0.0021128417440888_DP
         casimir_omega(97) = 0.0932602144609200_DP ;  casimir_omega_weight(97) = 0.0021439540958051_DP
         casimir_omega(98) = 0.0954198600431853_DP ;  casimir_omega_weight(98) = 0.0021753918104775_DP
         casimir_omega(99) = 0.0976111086070657_DP ;  casimir_omega_weight(99) = 0.0022071611045203_DP
         casimir_omega(100) = 0.0998342948875229_DP ;  casimir_omega_weight(100) = 0.0022392683083477_DP
         casimir_omega(101) = 0.1020897600080979_DP ;  casimir_omega_weight(101) = 0.0022717198691756_DP
         casimir_omega(102) = 0.1043778515991442_DP ;  casimir_omega_weight(102) = 0.0023045223538998_DP
         casimir_omega(103) = 0.1066989239189732_DP ;  casimir_omega_weight(103) = 0.0023376824520392_DP
         casimir_omega(104) = 0.1090533379779873_DP ;  casimir_omega_weight(104) = 0.0023712069787653_DP
         casimir_omega(105) = 0.1114414616658778_DP ;  casimir_omega_weight(105) = 0.0024051028780111_DP
         casimir_omega(106) = 0.1138636698819731_DP ;  casimir_omega_weight(106) = 0.0024393772256547_DP
         casimir_omega(107) = 0.1163203446688146_DP ;  casimir_omega_weight(107) = 0.0024740372327980_DP
         casimir_omega(108) = 0.1188118753490515_DP ;  casimir_omega_weight(108) = 0.0025090902491227_DP
         casimir_omega(109) = 0.1213386586657407_DP ;  casimir_omega_weight(109) = 0.0025445437663424_DP
         casimir_omega(110) = 0.1239010989261414_DP ;  casimir_omega_weight(110) = 0.0025804054217503_DP
         casimir_omega(111) = 0.1264996081491073_DP ;  casimir_omega_weight(111) = 0.0026166830018548_DP
         casimir_omega(112) = 0.1291346062161613_DP ;  casimir_omega_weight(112) = 0.0026533844461183_DP
         casimir_omega(113) = 0.1318065210263644_DP ;  casimir_omega_weight(113) = 0.0026905178508013_DP
         casimir_omega(114) = 0.1345157886550730_DP ;  casimir_omega_weight(114) = 0.0027280914729038_DP
         casimir_omega(115) = 0.1372628535167016_DP ;  casimir_omega_weight(115) = 0.0027661137342251_DP
         casimir_omega(116) = 0.1400481685315912_DP ;  casimir_omega_weight(116) = 0.0028045932255240_DP
         casimir_omega(117) = 0.1428721952971047_DP ;  casimir_omega_weight(117) = 0.0028435387108087_DP
         casimir_omega(118) = 0.1457354042630638_DP ;  casimir_omega_weight(118) = 0.0028829591317315_DP
         casimir_omega(119) = 0.1486382749116473_DP ;  casimir_omega_weight(119) = 0.0029228636121201_DP
         casimir_omega(120) = 0.1515812959418796_DP ;  casimir_omega_weight(120) = 0.0029632614626234_DP
         casimir_omega(121) = 0.1545649654588351_DP ;  casimir_omega_weight(121) = 0.0030041621855009_DP
         casimir_omega(122) = 0.1575897911676958_DP ;  casimir_omega_weight(122) = 0.0030455754795397_DP
         casimir_omega(123) = 0.1606562905727938_DP ;  casimir_omega_weight(123) = 0.0030875112451096_DP
         casimir_omega(124) = 0.1637649911817899_DP ;  casimir_omega_weight(124) = 0.0031299795893746_DP
         casimir_omega(125) = 0.1669164307151257_DP ;  casimir_omega_weight(125) = 0.0031729908316411_DP
         casimir_omega(126) = 0.1701111573209100_DP ;  casimir_omega_weight(126) = 0.0032165555088639_DP
         casimir_omega(127) = 0.1733497297953853_DP ;  casimir_omega_weight(127) = 0.0032606843813187_DP
         casimir_omega(128) = 0.1766327178091522_DP ;  casimir_omega_weight(128) = 0.0033053884384263_DP
         casimir_omega(129) = 0.1799607021393002_DP ;  casimir_omega_weight(129) = 0.0033506789047632_DP
         casimir_omega(130) = 0.1833342749076353_DP ;  casimir_omega_weight(130) = 0.0033965672462273_DP
         casimir_omega(131) = 0.1867540398251697_DP ;  casimir_omega_weight(131) = 0.0034430651764102_DP
         casimir_omega(132) = 0.1902206124430686_DP ;  casimir_omega_weight(132) = 0.0034901846631323_DP
         casimir_omega(133) = 0.1937346204102377_DP ;  casimir_omega_weight(133) = 0.0035379379351935_DP
         casimir_omega(134) = 0.1972967037377503_DP ;  casimir_omega_weight(134) = 0.0035863374893065_DP
         casimir_omega(135) = 0.2009075150703244_DP ;  casimir_omega_weight(135) = 0.0036353960972530_DP
         casimir_omega(136) = 0.2045677199650496_DP ;  casimir_omega_weight(136) = 0.0036851268132453_DP
         casimir_omega(137) = 0.2082779971775968_DP ;  casimir_omega_weight(137) = 0.0037355429815114_DP
         casimir_omega(138) = 0.2120390389561183_DP ;  casimir_omega_weight(138) = 0.0037866582441105_DP
         casimir_omega(139) = 0.2158515513430867_DP ;  casimir_omega_weight(139) = 0.0038384865489892_DP
         casimir_omega(140) = 0.2197162544853136_DP ;  casimir_omega_weight(140) = 0.0038910421582761_DP
         casimir_omega(141) = 0.2236338829523817_DP ;  casimir_omega_weight(141) = 0.0039443396568394_DP
         casimir_omega(142) = 0.2276051860637744_DP ;  casimir_omega_weight(142) = 0.0039983939611005_DP
         casimir_omega(143) = 0.2316309282249470_DP ;  casimir_omega_weight(143) = 0.0040532203281219_DP
         casimir_omega(144) = 0.2357118892726343_DP ;  casimir_omega_weight(144) = 0.0041088343649807_DP
         casimir_omega(145) = 0.2398488648296658_DP ;  casimir_omega_weight(145) = 0.0041652520384306_DP
         casimir_omega(146) = 0.2440426666696013_DP ;  casimir_omega_weight(146) = 0.0042224896848622_DP
         casimir_omega(147) = 0.2482941230914837_DP ;  casimir_omega_weight(147) = 0.0042805640205895_DP
         casimir_omega(148) = 0.2526040793050292_DP ;  casimir_omega_weight(148) = 0.0043394921524429_DP
         casimir_omega(149) = 0.2569733978265878_DP ;  casimir_omega_weight(149) = 0.0043992915887103_DP
         casimir_omega(150) = 0.2614029588862157_DP ;  casimir_omega_weight(150) = 0.0044599802504231_DP
         casimir_omega(151) = 0.2658936608462080_DP ;  casimir_omega_weight(151) = 0.0045215764830013_DP
         casimir_omega(152) = 0.2704464206314652_DP ;  casimir_omega_weight(152) = 0.0045840990682746_DP
         casimir_omega(153) = 0.2750621741720727_DP ;  casimir_omega_weight(153) = 0.0046475672368802_DP
         casimir_omega(154) = 0.2797418768584816_DP ;  casimir_omega_weight(154) = 0.0047120006810801_DP
         casimir_omega(155) = 0.2844865040097055_DP ;  casimir_omega_weight(155) = 0.0047774195679799_DP
         casimir_omega(156) = 0.2892970513549569_DP ;  casimir_omega_weight(156) = 0.0048438445531793_DP
         casimir_omega(157) = 0.2941745355291522_DP ;  casimir_omega_weight(157) = 0.0049112967948793_DP
         casimir_omega(158) = 0.2991199945827570_DP ;  casimir_omega_weight(158) = 0.0049797979684433_DP
         casimir_omega(159) = 0.3041344885064309_DP ;  casimir_omega_weight(159) = 0.0050493702814413_DP
         casimir_omega(160) = 0.3092190997709630_DP ;  casimir_omega_weight(160) = 0.0051200364892003_DP
         casimir_omega(161) = 0.3143749338830158_DP ;  casimir_omega_weight(161) = 0.0051918199108585_DP
         casimir_omega(162) = 0.3196031199571932_DP ;  casimir_omega_weight(162) = 0.0052647444459727_DP
         casimir_omega(163) = 0.3249048113049859_DP ;  casimir_omega_weight(163) = 0.0053388345916679_DP
         casimir_omega(164) = 0.3302811860411619_DP ;  casimir_omega_weight(164) = 0.0054141154603794_DP
         casimir_omega(165) = 0.3357334477081920_DP ;  casimir_omega_weight(165) = 0.0054906127981845_DP
         casimir_omega(166) = 0.3412628259193199_DP ;  casimir_omega_weight(166) = 0.0055683530037641_DP
         casimir_omega(167) = 0.3468705770209167_DP ;  casimir_omega_weight(167) = 0.0056473631480145_DP
         casimir_omega(168) = 0.3525579847747750_DP ;  casimir_omega_weight(168) = 0.0057276709943228_DP
         casimir_omega(169) = 0.3583263610610359_DP ;  casimir_omega_weight(169) = 0.0058093050195538_DP
         casimir_omega(170) = 0.3641770466024546_DP ;  casimir_omega_weight(170) = 0.0058922944357516_DP
         casimir_omega(171) = 0.3701114117107491_DP ;  casimir_omega_weight(171) = 0.0059766692126047_DP
         casimir_omega(172) = 0.3761308570557985_DP ;  casimir_omega_weight(172) = 0.0060624601006904_DP
         casimir_omega(173) = 0.3822368144584988_DP ;  casimir_omega_weight(173) = 0.0061496986555412_DP
         casimir_omega(174) = 0.3884307477080935_DP ;  casimir_omega_weight(174) = 0.0062384172625494_DP
         casimir_omega(175) = 0.3947141534048617_DP ;  casimir_omega_weight(175) = 0.0063286491627608_DP
         casimir_omega(176) = 0.4010885618290448_DP ;  casimir_omega_weight(176) = 0.0064204284795907_DP
         casimir_omega(177) = 0.4075555378369654_DP ;  casimir_omega_weight(177) = 0.0065137902464749_DP
         casimir_omega(178) = 0.4141166817852959_DP ;  casimir_omega_weight(178) = 0.0066087704355375_DP
         casimir_omega(179) = 0.4207736304844983_DP ;  casimir_omega_weight(179) = 0.0067054059872730_DP
         casimir_omega(180) = 0.4275280581824882_DP ;  casimir_omega_weight(180) = 0.0068037348413109_DP
         casimir_omega(181) = 0.4343816775796169_DP ;  casimir_omega_weight(181) = 0.0069037959683042_DP
         casimir_omega(182) = 0.4413362408761182_DP ;  casimir_omega_weight(182) = 0.0070056294029701_DP
         casimir_omega(183) = 0.4483935408532002_DP ;  casimir_omega_weight(183) = 0.0071092762783565_DP
         casimir_omega(184) = 0.4555554119890342_DP ;  casimir_omega_weight(184) = 0.0072147788613694_DP
         casimir_omega(185) = 0.4628237316109162_DP ;  casimir_omega_weight(185) = 0.0073221805896055_DP
         casimir_omega(186) = 0.4702004210849546_DP ;  casimir_omega_weight(186) = 0.0074315261095804_DP
         casimir_omega(187) = 0.4776874470446816_DP ;  casimir_omega_weight(187) = 0.0075428613163553_DP
         casimir_omega(188) = 0.4852868226600471_DP ;  casimir_omega_weight(188) = 0.0076562333946817_DP
         casimir_omega(189) = 0.4930006089483120_DP ;  casimir_omega_weight(189) = 0.0077716908616867_DP
         casimir_omega(190) = 0.5008309161284398_DP ;  casimir_omega_weight(190) = 0.0078892836111662_DP
         casimir_omega(191) = 0.5087799050206278_DP ;  casimir_omega_weight(191) = 0.0080090629596010_DP
         casimir_omega(192) = 0.5168497884927058_DP ;  casimir_omega_weight(192) = 0.0081310816938891_DP
         casimir_omega(193) = 0.5250428329552088_DP ;  casimir_omega_weight(193) = 0.0082553941209621_DP
         casimir_omega(194) = 0.5333613599069924_DP ;  casimir_omega_weight(194) = 0.0083820561192824_DP
         casimir_omega(195) = 0.5418077475333593_DP ;  casimir_omega_weight(195) = 0.0085111251923702_DP
         casimir_omega(196) = 0.5503844323587329_DP ;  casimir_omega_weight(196) = 0.0086426605243974_DP
         casimir_omega(197) = 0.5590939109560277_DP ;  casimir_omega_weight(197) = 0.0087767230379865_DP
         casimir_omega(198) = 0.5679387417149295_DP ;  casimir_omega_weight(198) = 0.0089133754542690_DP
         casimir_omega(199) = 0.5769215466714336_DP ;  casimir_omega_weight(199) = 0.0090526823553410_DP
         casimir_omega(200) = 0.5860450134010606_DP ;  casimir_omega_weight(200) = 0.0091947102491941_DP
         casimir_omega(201) = 0.5953118969782984_DP ;  casimir_omega_weight(201) = 0.0093395276372536_DP
         casimir_omega(202) = 0.6047250220049346_DP ;  casimir_omega_weight(202) = 0.0094872050846297_DP
         casimir_omega(203) = 0.6142872847100472_DP ;  casimir_omega_weight(203) = 0.0096378152932223_DP
         casimir_omega(204) = 0.6240016551245678_DP ;  casimir_omega_weight(204) = 0.0097914331777712_DP
         casimir_omega(205) = 0.6338711793334540_DP ;  casimir_omega_weight(205) = 0.0099481359450482_DP
         casimir_omega(206) = 0.6438989818086474_DP ;  casimir_omega_weight(206) = 0.0101080031762688_DP
         casimir_omega(207) = 0.6540882678261467_DP ;  casimir_omega_weight(207) = 0.0102711169129299_DP
         casimir_omega(208) = 0.6644423259706811_DP ;  casimir_omega_weight(208) = 0.0104375617461949_DP
         casimir_omega(209) = 0.6749645307316174_DP ;  casimir_omega_weight(209) = 0.0106074249100165_DP
         casimir_omega(210) = 0.6856583451939287_DP ;  casimir_omega_weight(210) = 0.0107807963781598_DP
         casimir_omega(211) = 0.6965273238282100_DP ;  casimir_omega_weight(211) = 0.0109577689653113_DP
         casimir_omega(212) = 0.7075751153839397_DP ;  casimir_omega_weight(212) = 0.0111384384324820_DP
         casimir_omega(213) = 0.7188054658903634_DP ;  casimir_omega_weight(213) = 0.0113229035968931_DP
         casimir_omega(214) = 0.7302222217695975_DP ;  casimir_omega_weight(214) = 0.0115112664465728_DP
         casimir_omega(215) = 0.7418293330667814_DP ;  casimir_omega_weight(215) = 0.0117036322598849_DP
         casimir_omega(216) = 0.7536308568023278_DP ;  casimir_omega_weight(216) = 0.0119001097302304_DP
         casimir_omega(217) = 0.7656309604515563_DP ;  casimir_omega_weight(217) = 0.0121008110961913_DP
         casimir_omega(218) = 0.7778339255573071_DP ;  casimir_omega_weight(218) = 0.0123058522773438_DP
         casimir_omega(219) = 0.7902441514813261_DP ;  casimir_omega_weight(219) = 0.0125153530160790_DP
         casimir_omega(220) = 0.8028661593005888_DP ;  casimir_omega_weight(220) = 0.0127294370256781_DP
         casimir_omega(221) = 0.8157045958549579_DP ;  casimir_omega_weight(221) = 0.0129482321449843_DP
         casimir_omega(222) = 0.8287642379529623_DP ;  casimir_omega_weight(222) = 0.0131718705000078_DP
         casimir_omega(223) = 0.8420499967427544_DP ;  casimir_omega_weight(223) = 0.0134004886727674_DP
         casimir_omega(224) = 0.8555669222557496_DP ;  casimir_omega_weight(224) = 0.0136342278778250_DP
         casimir_omega(225) = 0.8693202081307284_DP ;  casimir_omega_weight(225) = 0.0138732341468160_DP
         casimir_omega(226) = 0.8833151965266887_DP ;  casimir_omega_weight(226) = 0.0141176585214292_DP
         casimir_omega(227) = 0.8975573832330863_DP ;  casimir_omega_weight(227) = 0.0143676572552878_DP
         casimir_omega(228) = 0.9120524229865773_DP ;  casimir_omega_weight(228) = 0.0146233920251862_DP
         casimir_omega(229) = 0.9268061350038654_DP ;  casimir_omega_weight(229) = 0.0148850301521155_DP
         casimir_omega(230) = 0.9418245087407373_DP ;  casimir_omega_weight(230) = 0.0151527448327110_DP
         casimir_omega(231) = 0.9571137098879237_DP ;  casimir_omega_weight(231) = 0.0154267153815543_DP
         casimir_omega(232) = 0.9726800866149689_DP ;  casimir_omega_weight(232) = 0.0157071274849774_DP
         casimir_omega(233) = 0.9885301760739067_DP ;  casimir_omega_weight(233) = 0.0159941734669829_DP
         casimir_omega(234) = 1.0046707111751640_DP ;  casimir_omega_weight(234) = 0.0162880525679066_DP
         casimir_omega(235) = 1.0211086276488077_DP ;  casimir_omega_weight(235) = 0.0165889712365289_DP
         casimir_omega(236) = 1.0378510714049163_DP ;  casimir_omega_weight(236) = 0.0168971434363695_DP
         casimir_omega(237) = 1.0549054062076764_DP ;  casimir_omega_weight(237) = 0.0172127909669588_DP
         casimir_omega(238) = 1.0722792216785617_DP ;  casimir_omega_weight(238) = 0.0175361438008755_DP
         casimir_omega(239) = 1.0899803416448139_DP ;  casimir_omega_weight(239) = 0.0178674404374718_DP
         casimir_omega(240) = 1.1080168328503788_DP ;  casimir_omega_weight(240) = 0.0182069282741979_DP
         casimir_omega(241) = 1.1263970140473520_DP ;  casimir_omega_weight(241) = 0.0185548639965102_DP
         casimir_omega(242) = 1.1451294654870989_DP ;  casimir_omega_weight(242) = 0.0189115139874329_DP
         casimir_omega(243) = 1.1642230388312049_DP ;  casimir_omega_weight(243) = 0.0192771547578969_DP
         casimir_omega(244) = 1.1836868675036438_DP ;  casimir_omega_weight(244) = 0.0196520733990311_DP
         casimir_omega(245) = 1.2035303775067407_DP ;  casimir_omega_weight(245) = 0.0200365680576697_DP
         casimir_omega(246) = 1.2237632987248304_DP ;  casimir_omega_weight(246) = 0.0204309484364698_DP
         casimir_omega(247) = 1.2443956767409057_DP ;  casimir_omega_weight(247) = 0.0208355363200194_DP
         casimir_omega(248) = 1.2654378851930272_DP ;  casimir_omega_weight(248) = 0.0212506661285319_DP
         casimir_omega(249) = 1.2869006386988686_DP ;  casimir_omega_weight(249) = 0.0216766855006805_DP
         casimir_omega(250) = 1.3087950063784219_DP ;  casimir_omega_weight(250) = 0.0221139559074085_DP
         casimir_omega(251) = 1.3311324260067356_DP ;  casimir_omega_weight(251) = 0.0225628532984813_DP
         casimir_omega(252) = 1.3539247188304484_DP ;  casimir_omega_weight(252) = 0.0230237687838673_DP
         casimir_omega(253) = 1.3771841050839133_DP ;  casimir_omega_weight(253) = 0.0234971093519468_DP
         casimir_omega(254) = 1.4009232202429664_DP ;  casimir_omega_weight(254) = 0.0239832986268981_DP
         casimir_omega(255) = 1.4251551320566243_DP ;  casimir_omega_weight(255) = 0.0244827776676733_DP
         casimir_omega(256) = 1.4498933583996223_DP ;  casimir_omega_weight(256) = 0.0249960058111003_DP
         casimir_omega(257) = 1.4751518859912669_DP ;  casimir_omega_weight(257) = 0.0255234615618956_DP
         casimir_omega(258) = 1.5009451900289918_DP ;  casimir_omega_weight(258) = 0.0260656435325631_DP
         casimir_omega(259) = 1.5272882547880704_DP ;  casimir_omega_weight(259) = 0.0266230714363573_DP
         casimir_omega(260) = 1.5541965952421848_DP ;  casimir_omega_weight(260) = 0.0271962871366454_DP
         casimir_omega(261) = 1.5816862797631044_DP ;  casimir_omega_weight(261) = 0.0277858557563301_DP
         casimir_omega(262) = 1.6097739539614166_DP ;  casimir_omega_weight(262) = 0.0283923668512722_DP
         casimir_omega(263) = 1.6384768657344082_DP ;  casimir_omega_weight(263) = 0.0290164356518265_DP
         casimir_omega(264) = 1.6678128915913812_DP ;  casimir_omega_weight(264) = 0.0296587043770102_DP
         casimir_omega(265) = 1.6978005643314691_DP ;  casimir_omega_weight(265) = 0.0303198436261195_DP
         casimir_omega(266) = 1.7284591021538895_DP ;  casimir_omega_weight(266) = 0.0310005538529085_DP
         casimir_omega(267) = 1.7598084392860496_DP ;  casimir_omega_weight(267) = 0.0317015669280388_DP
         casimir_omega(268) = 1.7918692582205737_DP ;  casimir_omega_weight(268) = 0.0324236477955690_DP
         casimir_omega(269) = 1.8246630236586112_DP ;  casimir_omega_weight(269) = 0.0331675962300713_DP
         casimir_omega(270) = 1.8582120182633943_DP ;  casimir_omega_weight(270) = 0.0339342487012650_DP
         casimir_omega(271) = 1.8925393803352639_DP ;  casimir_omega_weight(271) = 0.0347244803535625_DP
         casimir_omega(272) = 1.9276691435270426_DP ;  casimir_omega_weight(272) = 0.0355392071086158_DP
         casimir_omega(273) = 1.9636262787270433_DP ;  casimir_omega_weight(273) = 0.0363793878995943_DP
         casimir_omega(274) = 2.0004367382459942_DP ;  casimir_omega_weight(274) = 0.0372460270463118_DP
         casimir_omega(275) = 2.0381275024538357_DP ;  casimir_omega_weight(275) = 0.0381401767815746_DP
         casimir_omega(276) = 2.0767266290228910_DP ;  casimir_omega_weight(276) = 0.0390629399394835_DP
         casimir_omega(277) = 2.1162633049451878_DP ;  casimir_omega_weight(277) = 0.0400154728175015_DP
         casimir_omega(278) = 2.1567679015039958_DP ;  casimir_omega_weight(278) = 0.0409989882250127_DP
         casimir_omega(279) = 2.1982720323929028_DP ;  casimir_omega_weight(279) = 0.0420147587321960_DP
         casimir_omega(280) = 2.2408086151900957_DP ;  casimir_omega_weight(280) = 0.0430641201341791_DP
         casimir_omega(281) = 2.2844119364109874_DP ;  casimir_omega_weight(281) = 0.0441484751464451_DP
         casimir_omega(282) = 2.3291177203793731_DP ;  casimir_omega_weight(282) = 0.0452692973493293_DP
         casimir_omega(283) = 2.3749632021752429_DP ;  casimir_omega_weight(283) = 0.0464281354003823_DP
         casimir_omega(284) = 2.4219872049375550_DP ;  casimir_omega_weight(284) = 0.0476266175353529_DP
         casimir_omega(285) = 2.4702302218215406_DP ;  casimir_omega_weight(285) = 0.0488664563802553_DP
         casimir_omega(286) = 2.5197345029337241_DP ;  casimir_omega_weight(286) = 0.0501494540988385_DP
         casimir_omega(287) = 2.5705441475930093_DP ;  casimir_omega_weight(287) = 0.0514775079020269_DP
         casimir_omega(288) = 2.6227052022941995_DP ;  casimir_omega_weight(288) = 0.0528526159483203_DP
         casimir_omega(289) = 2.6762657647803407_DP ;  casimir_omega_weight(289) = 0.0542768836663796_DP
         casimir_omega(290) = 2.7312760946629644_DP ;  casimir_omega_weight(290) = 0.0557525305344060_DP
         casimir_omega(291) = 2.7877887310655312_DP ;  casimir_omega_weight(291) = 0.0572818973536120_DP
         casimir_omega(292) = 2.8458586178042764_DP ;  casimir_omega_weight(292) = 0.0588674540565496_DP
         casimir_omega(293) = 2.9055432366633052_DP ;  casimir_omega_weight(293) = 0.0605118080953108_DP
         casimir_omega(294) = 2.9669027493679132_DP ;  casimir_omega_weight(294) = 0.0622177134580479_DP
         casimir_omega(295) = 3.0300001489108062_DP ;  casimir_omega_weight(295) = 0.0639880803676890_DP
         casimir_omega(296) = 3.0949014209421883_DP ;  casimir_omega_weight(296) = 0.0658259857211671_DP
         casimir_omega(297) = 3.1616757159958229_DP ;  casimir_omega_weight(297) = 0.0677346843334869_DP
         casimir_omega(298) = 3.2303955333908614_DP ;  casimir_omega_weight(298) = 0.0697176210571101_DP
         casimir_omega(299) = 3.3011369177224799_DP ;  casimir_omega_weight(299) = 0.0717784438541946_DP
         casimir_omega(300) = 3.3739796689363311_DP ;  casimir_omega_weight(300) = 0.0739210179064062_DP
         casimir_omega(301) = 3.4490075670704115_DP ;  casimir_omega_weight(301) = 0.0761494408565468_DP
         casimir_omega(302) = 3.5263086128466177_DP ;  casimir_omega_weight(302) = 0.0784680592844421_DP
         casimir_omega(303) = 3.6059752854025713_DP ;  casimir_omega_weight(303) = 0.0808814865310849_DP
         casimir_omega(304) = 3.6881048185732892_DP ;  casimir_omega_weight(304) = 0.0833946219965227_DP
         casimir_omega(305) = 3.7727994972647245_DP ;  casimir_omega_weight(305) = 0.0860126720497235_DP
         casimir_omega(306) = 3.8601669756062500_DP ;  casimir_omega_weight(306) = 0.0887411727035335_DP
         casimir_omega(307) = 3.9503206187306854_DP ;  casimir_omega_weight(307) = 0.0915860142244031_DP
         casimir_omega(308) = 4.0433798702086925_DP ;  casimir_omega_weight(308) = 0.0945534678643644_DP
         casimir_omega(309) = 4.1394706473620255_DP ;  casimir_omega_weight(309) = 0.0976502149233557_DP
         casimir_omega(310) = 4.2387257668996252_DP ;  casimir_omega_weight(310) = 0.1008833783731321_DP
         casimir_omega(311) = 4.3412854035643473_DP ;  casimir_omega_weight(311) = 0.1042605572994631_DP
         casimir_omega(312) = 4.4472975847485152_DP ;  casimir_omega_weight(312) = 0.1077898644483133_DP
         casimir_omega(313) = 4.5569187243392673_DP ;  casimir_omega_weight(313) = 0.1114799671946000_DP
         casimir_omega(314) = 4.6703141993896091_DP ;  casimir_omega_weight(314) = 0.1153401322879679_DP
         casimir_omega(315) = 4.7876589735875710_DP ;  casimir_omega_weight(315) = 0.1193802747724687_DP
         casimir_omega(316) = 4.9091382719143297_DP ;  casimir_omega_weight(316) = 0.1236110115230931_DP
         casimir_omega(317) = 5.0349483113513278_DP ;  casimir_omega_weight(317) = 0.1280437198947250_DP
         casimir_omega(318) = 5.1652970930221445_DP ;  casimir_omega_weight(318) = 0.1326906020406465_DP
         casimir_omega(319) = 5.3004052617437418_DP ;  casimir_omega_weight(319) = 0.1375647555228773_DP
         casimir_omega(320) = 5.4405070396232658_DP ;  casimir_omega_weight(320) = 0.1426802509168860_DP
         casimir_omega(321) = 5.5858512410814116_DP ;  casimir_omega_weight(321) = 0.1480522171995280_DP
         casimir_omega(322) = 5.7367023775208255_DP ;  casimir_omega_weight(322) = 0.1536969358080894_DP
         casimir_omega(323) = 5.8933418608030506_DP ;  casimir_omega_weight(323) = 0.1596319443765841_DP
         casimir_omega(324) = 6.0560693157662095_DP ;  casimir_omega_weight(324) = 0.1658761512808831_DP
         casimir_omega(325) = 6.2252040132208180_DP ;  casimir_omega_weight(325) = 0.1724499622786498_DP
         casimir_omega(326) = 6.4010864362316715_DP ;  casimir_omega_weight(326) = 0.1793754206986603_DP
         casimir_omega(327) = 6.5840799940435240_DP ;  casimir_omega_weight(327) = 0.1866763628332005_DP
         casimir_omega(328) = 6.7745728997736530_DP ;  casimir_omega_weight(328) = 0.1943785904132996_DP
         casimir_omega(329) = 6.9729802300018768_DP ;  casimir_omega_weight(329) = 0.2025100623089622_DP
         casimir_omega(330) = 7.1797461866781553_DP ;  casimir_omega_weight(330) = 0.2111011078997880_DP
         casimir_omega(331) = 7.3953465843842610_DP ;  casimir_omega_weight(331) = 0.2201846649101197_DP
         casimir_omega(332) = 7.6202915889787981_DP ;  casimir_omega_weight(332) = 0.2297965449108316_DP
         casimir_omega(333) = 7.8551287370863916_DP ;  casimir_omega_weight(333) = 0.2399757301628021_DP
         casimir_omega(334) = 8.1004462698353628_DP ;  casimir_omega_weight(334) = 0.2507647060232721_DP
         casimir_omega(335) = 8.3568768187825544_DP ;  casimir_omega_weight(335) = 0.2622098337853790_DP
         casimir_omega(336) = 8.6251014871992329_DP ;  casimir_omega_weight(336) = 0.2743617695644203_DP
         casimir_omega(337) = 8.9058543759347017_DP ;  casimir_omega_weight(337) = 0.2872759357332177_DP
         casimir_omega(338) = 9.1999276100797758_DP ;  casimir_omega_weight(338) = 0.3010130524359107_DP
         casimir_omega(339) = 9.5081769307755248_DP ;  casimir_omega_weight(339) = 0.3156397379366003_DP
         casimir_omega(340) = 9.8315279259679915_DP ;  casimir_omega_weight(340) = 0.3312291879945023_DP
         casimir_omega(341) = 10.1709829849288624_DP ;  casimir_omega_weight(341) = 0.3478619461637186_DP
         casimir_omega(342) = 10.5276290742460610_DP ;  casimir_omega_weight(342) = 0.3656267789381978_DP
         casimir_omega(343) = 10.9026464480758971_DP ;  casimir_omega_weight(343) = 0.3846216720789509_DP
         casimir_omega(344) = 11.2973184231809274_DP ;  casimir_omega_weight(344) = 0.4049549673245511_DP
         casimir_omega(345) = 11.7130423701493527_DP ;  casimir_omega_weight(345) = 0.4267466621437209_DP
         casimir_omega(346) = 12.1513420968532948_DP ;  casimir_omega_weight(346) = 0.4501298993085213_DP
         casimir_omega(347) = 12.6138818293974122_DP ;  casimir_omega_weight(347) = 0.4752526780637042_DP
         casimir_omega(348) = 13.1024820304914655_DP ;  casimir_omega_weight(348) = 0.5022798246691069_DP
         casimir_omega(349) = 13.6191373364972215_DP ;  casimir_omega_weight(349) = 0.5313952674027976_DP
         casimir_omega(350) = 14.1660369437710596_DP ;  casimir_omega_weight(350) = 0.5628046699817101_DP
         casimir_omega(351) = 14.7455878341225031_DP ;  casimir_omega_weight(351) = 0.5967384881951940_DP
         casimir_omega(352) = 15.3604413004202893_DP ;  casimir_omega_weight(352) = 0.6334555278363457_DP
         casimir_omega(353) = 16.0135233193396829_DP ;  casimir_omega_weight(353) = 0.6732470983424110_DP
         casimir_omega(354) = 16.7080694223972976_DP ;  casimir_omega_weight(354) = 0.7164418767470071_DP
         casimir_omega(355) = 17.4476648430786412_DP ;  casimir_omega_weight(355) = 0.7634116215712567_DP
         casimir_omega(356) = 18.2362908725025754_DP ;  casimir_omega_weight(356) = 0.8145779074550850_DP
         casimir_omega(357) = 19.0783785456327060_DP ;  casimir_omega_weight(357) = 0.8704200903310001_DP
         casimir_omega(358) = 19.9788710134306733_DP ;  casimir_omega_weight(358) = 0.9314847619840412_DP
         casimir_omega(359) = 20.9432962449260707_DP ;  casimir_omega_weight(359) = 0.9983970147671912_DP
         casimir_omega(360) = 21.9778520616958311_DP ;  casimir_omega_weight(360) = 1.0718739158791881_DP
         casimir_omega(361) = 23.0895059547737382_DP ;  casimir_omega_weight(361) = 1.1527406909915912_DP
         casimir_omega(362) = 24.2861126955058424_DP ;  casimir_omega_weight(362) = 1.2419502458283669_DP
         casimir_omega(363) = 25.5765534600501745_DP ;  casimir_omega_weight(363) = 1.3406068206306190_DP
         casimir_omega(364) = 26.9709010854142726_DP ;  casimir_omega_weight(364) = 1.4499947884643525_DP
         casimir_omega(365) = 28.4806172207251898_DP ;  casimir_omega_weight(365) = 1.5716138907566590_DP
         casimir_omega(366) = 30.1187886080231593_DP ;  casimir_omega_weight(366) = 1.7072225751085726_DP
         casimir_omega(367) = 31.9004116265637876_DP ;  casimir_omega_weight(367) = 1.8588915930433756_DP
         casimir_omega(368) = 33.8427367051229311_DP ;  casimir_omega_weight(368) = 2.0290706731913444_DP
         casimir_omega(369) = 35.9656874427775861_DP ;  casimir_omega_weight(369) = 2.2206719706519906_DP
         casimir_omega(370) = 38.2923735491411108_DP ;  casimir_omega_weight(370) = 2.4371751947064588_DP
         casimir_omega(371) = 40.8497223959989455_DP ;  casimir_omega_weight(371) = 2.6827609617476762_DP
         casimir_omega(372) = 43.6692615934087129_DP ;  casimir_omega_weight(372) = 2.9624811929116515_DP
         casimir_omega(373) = 46.7880953193943441_DP ;  casimir_omega_weight(373) = 3.2824785473745628_DP
         casimir_omega(374) = 50.2501312296495399_DP ;  casimir_omega_weight(374) = 3.6502713546476486_DP
         casimir_omega(375) = 54.1076342347988444_DP ;  casimir_omega_weight(375) = 4.0751268871230755_DP
         casimir_omega(376) = 58.4232105909666615_DP ;  casimir_omega_weight(376) = 4.5685550187347772_DP
         casimir_omega(377) = 63.2723640908453220_DP ;  casimir_omega_weight(377) = 5.1449677703126460_DP
         casimir_omega(378) = 68.7468209493643627_DP ;  casimir_omega_weight(378) = 5.8225701812480928_DP
         casimir_omega(379) = 74.9588993758562907_DP ;  casimir_omega_weight(379) = 6.6245779342690323_DP
         casimir_omega(380) = 82.0473165181309980_DP ;  casimir_omega_weight(380) = 7.5809029822217147_DP
         casimir_omega(381) = 90.1849996688832789_DP ;  casimir_omega_weight(381) = 8.7305196611616473_DP
         casimir_omega(382) = 99.5897331286242036_DP ;  casimir_omega_weight(382) = 10.1248366061488326_DP
         casimir_omega(383) = 110.5388812084258774_DP ;  casimir_omega_weight(383) = 11.8325821836009215_DP
         casimir_omega(384) = 123.3900734674389668_DP ;  casimir_omega_weight(384) = 13.9470126541790620_DP
         casimir_omega(385) = 138.6107800162292847_DP ;  casimir_omega_weight(385) = 16.5967630181972119_DP
         casimir_omega(386) = 156.8214272964054317_DP ;  casimir_omega_weight(386) = 19.9625494166238866_DP
         casimir_omega(387) = 178.8596316424528823_DP ;  casimir_omega_weight(387) = 24.3035263143612994_DP
         casimir_omega(388) = 205.8782544350598869_DP ;  casimir_omega_weight(388) = 30.0000589998744438_DP
         casimir_omega(389) = 239.4992737120479092_DP ;  casimir_omega_weight(389) = 37.6253687142619455_DP
         casimir_omega(390) = 282.0629703600874336_DP ;  casimir_omega_weight(390) = 48.0699626516609300_DP
         casimir_omega(391) = 337.0463938304238241_DP ;  casimir_omega_weight(391) = 62.7669539867800168_DP
         casimir_omega(392) = 409.7965078093066609_DP ;  casimir_omega_weight(392) = 84.1204624138253081_DP
         casimir_omega(393) = 508.8806166962729662_DP ;  casimir_omega_weight(393) = 116.3685356975828427_DP
         casimir_omega(394) = 648.7281519025392527_DP ;  casimir_omega_weight(394) = 167.4464436367215399_DP
         casimir_omega(395) = 855.1948562613241620_DP ;  casimir_omega_weight(395) = 253.3690725226675511_DP
         casimir_omega(396) = 1178.4242716744699919_DP ;  casimir_omega_weight(396) = 409.7130309512374424_DP
         casimir_omega(397) = 1726.4058368836410864_DP ;  casimir_omega_weight(397) = 726.2649574010173410_DP
         casimir_omega(398) = 2768.2997534873206860_DP ;  casimir_omega_weight(398) = 1474.0063901455218911_DP
         casimir_omega(399) = 5140.1553093006523341_DP ;  casimir_omega_weight(399) = 3726.2887315043080889_DP
         casimir_omega(400) = 12633.1691550323648698_DP ;  casimir_omega_weight(400) = 14323.9086236090588500_DP
         casimir_omega(401) = 66565.1233062650135253_DP ;  casimir_omega_weight(401) = 170828.6287803665036336_DP
return
endsubroutine gauss_legendre_grid400

subroutine gauss_legendre_grid450()
implicit none
         casimir_omega(1) = 0.0000000000000000_DP ;  casimir_omega_weight(1) = 0.0000000000000000_DP
         casimir_omega(2) = 0.0000042743561939_DP ;  casimir_omega_weight(2) = 0.0000109694304064_DP
         casimir_omega(3) = 0.0000225217847829_DP ;  casimir_omega_weight(3) = 0.0000255357813897_DP
         casimir_omega(4) = 0.0000553521725603_DP ;  casimir_omega_weight(4) = 0.0000401261828701_DP
         casimir_omega(5) = 0.0001027760040973_DP ;  casimir_omega_weight(5) = 0.0000547223608228_DP
         casimir_omega(6) = 0.0001647987859212_DP ;  casimir_omega_weight(6) = 0.0000693242386442_DP
         casimir_omega(7) = 0.0002414267658168_DP ;  casimir_omega_weight(7) = 0.0000839329745471_DP
         casimir_omega(8) = 0.0003326674716759_DP ;  casimir_omega_weight(8) = 0.0000985499217370_DP
         casimir_omega(9) = 0.0004385298132268_DP ;  casimir_omega_weight(9) = 0.0001131764812075_DP
         casimir_omega(10) = 0.0005590241104959_DP ;  casimir_omega_weight(10) = 0.0001278140701249_DP
         casimir_omega(11) = 0.0006941621047948_DP ;  casimir_omega_weight(11) = 0.0001424641131492_DP
         casimir_omega(12) = 0.0008439569645922_DP ;  casimir_omega_weight(12) = 0.0001571280396935_DP
         casimir_omega(13) = 0.0010084232897458_DP ;  casimir_omega_weight(13) = 0.0001718072830508_DP
         casimir_omega(14) = 0.0011875771152452_DP ;  casimir_omega_weight(14) = 0.0001865032801879_DP
         casimir_omega(15) = 0.0013814359149077_DP ;  casimir_omega_weight(15) = 0.0002012174718143_DP
         casimir_omega(16) = 0.0015900186052016_DP ;  casimir_omega_weight(16) = 0.0002159513025628_DP
         casimir_omega(17) = 0.0018133455492911_DP ;  casimir_omega_weight(17) = 0.0002307062212320_DP
         casimir_omega(18) = 0.0020514385613321_DP ;  casimir_omega_weight(18) = 0.0002454836810567_DP
         casimir_omega(19) = 0.0023043209110517_DP ;  casimir_omega_weight(19) = 0.0002602851399909_DP
         casimir_omega(20) = 0.0025720173286153_DP ;  casimir_omega_weight(20) = 0.0002751120610067_DP
         casimir_omega(21) = 0.0028545540097932_DP ;  casimir_omega_weight(21) = 0.0002899659123941_DP
         casimir_omega(22) = 0.0031519586214322_DP ;  casimir_omega_weight(22) = 0.0003048481680660_DP
         casimir_omega(23) = 0.0034642603072343_DP ;  casimir_omega_weight(23) = 0.0003197603078721_DP
         casimir_omega(24) = 0.0037914896938467_DP ;  casimir_omega_weight(24) = 0.0003347038179076_DP
         casimir_omega(25) = 0.0041336788972686_DP ;  casimir_omega_weight(25) = 0.0003496801908305_DP
         casimir_omega(26) = 0.0044908615295707_DP ;  casimir_omega_weight(26) = 0.0003646909261826_DP
         casimir_omega(27) = 0.0048630727059410_DP ;  casimir_omega_weight(27) = 0.0003797375307071_DP
         casimir_omega(28) = 0.0052503490520462_DP ;  casimir_omega_weight(28) = 0.0003948215186769_DP
         casimir_omega(29) = 0.0056527287117269_DP ;  casimir_omega_weight(29) = 0.0004099444122172_DP
         casimir_omega(30) = 0.0060702513550134_DP ;  casimir_omega_weight(30) = 0.0004251077416398_DP
         casimir_omega(31) = 0.0065029581864785_DP ;  casimir_omega_weight(31) = 0.0004403130457732_DP
         casimir_omega(32) = 0.0069508919539221_DP ;  casimir_omega_weight(32) = 0.0004555618723007_DP
         casimir_omega(33) = 0.0074140969573970_DP ;  casimir_omega_weight(33) = 0.0004708557780981_DP
         casimir_omega(34) = 0.0078926190585698_DP ;  casimir_omega_weight(34) = 0.0004861963295785_DP
         casimir_omega(35) = 0.0083865056904328_DP ;  casimir_omega_weight(35) = 0.0005015851030363_DP
         casimir_omega(36) = 0.0088958058673606_DP ;  casimir_omega_weight(36) = 0.0005170236850006_DP
         casimir_omega(37) = 0.0094205701955183_DP ;  casimir_omega_weight(37) = 0.0005325136725834_DP
         casimir_omega(38) = 0.0099608508836277_DP ;  casimir_omega_weight(38) = 0.0005480566738439_DP
         casimir_omega(39) = 0.0105167017540903_DP ;  casimir_omega_weight(39) = 0.0005636543081451_DP
         casimir_omega(40) = 0.0110881782544762_DP ;  casimir_omega_weight(40) = 0.0005793082065219_DP
         casimir_omega(41) = 0.0116753374693801_DP ;  casimir_omega_weight(41) = 0.0005950200120505_DP
         casimir_omega(42) = 0.0122782381326500_DP ;  casimir_omega_weight(42) = 0.0006107913802237_DP
         casimir_omega(43) = 0.0128969406399915_DP ;  casimir_omega_weight(43) = 0.0006266239793287_DP
         casimir_omega(44) = 0.0135315070619558_DP ;  casimir_omega_weight(44) = 0.0006425194908323_DP
         casimir_omega(45) = 0.0141820011573129_DP ;  casimir_omega_weight(45) = 0.0006584796097713_DP
         casimir_omega(46) = 0.0148484883868164_DP ;  casimir_omega_weight(46) = 0.0006745060451427_DP
         casimir_omega(47) = 0.0155310359273644_DP ;  casimir_omega_weight(47) = 0.0006906005203049_DP
         casimir_omega(48) = 0.0162297126865640_DP ;  casimir_omega_weight(48) = 0.0007067647733850_DP
         casimir_omega(49) = 0.0169445893177000_DP ;  casimir_omega_weight(49) = 0.0007230005576825_DP
         casimir_omega(50) = 0.0176757382351221_DP ;  casimir_omega_weight(50) = 0.0007393096420928_DP
         casimir_omega(51) = 0.0184232336300429_DP ;  casimir_omega_weight(51) = 0.0007556938115228_DP
         casimir_omega(52) = 0.0191871514867700_DP ;  casimir_omega_weight(52) = 0.0007721548673229_DP
         casimir_omega(53) = 0.0199675695993621_DP ;  casimir_omega_weight(53) = 0.0007886946277193_DP
         casimir_omega(54) = 0.0207645675887237_DP ;  casimir_omega_weight(54) = 0.0008053149282539_DP
         casimir_omega(55) = 0.0215782269201496_DP ;  casimir_omega_weight(55) = 0.0008220176222353_DP
         casimir_omega(56) = 0.0224086309213103_DP ;  casimir_omega_weight(56) = 0.0008388045811873_DP
         casimir_omega(57) = 0.0232558648007034_DP ;  casimir_omega_weight(57) = 0.0008556776953156_DP
         casimir_omega(58) = 0.0241200156665659_DP ;  casimir_omega_weight(58) = 0.0008726388739716_DP
         casimir_omega(59) = 0.0250011725462579_DP ;  casimir_omega_weight(59) = 0.0008896900461317_DP
         casimir_omega(60) = 0.0258994264061294_DP ;  casimir_omega_weight(60) = 0.0009068331608792_DP
         casimir_omega(61) = 0.0268148701718704_DP ;  casimir_omega_weight(61) = 0.0009240701878940_DP
         casimir_omega(62) = 0.0277475987493596_DP ;  casimir_omega_weight(62) = 0.0009414031179558_DP
         casimir_omega(63) = 0.0286977090460132_DP ;  casimir_omega_weight(63) = 0.0009588339634479_DP
         casimir_omega(64) = 0.0296652999926499_DP ;  casimir_omega_weight(64) = 0.0009763647588767_DP
         casimir_omega(65) = 0.0306504725658720_DP ;  casimir_omega_weight(65) = 0.0009939975613933_DP
         casimir_omega(66) = 0.0316533298109792_DP ;  casimir_omega_weight(66) = 0.0010117344513293_DP
         casimir_omega(67) = 0.0326739768654195_DP ;  casimir_omega_weight(67) = 0.0010295775327418_DP
         casimir_omega(68) = 0.0337125209827887_DP ;  casimir_omega_weight(68) = 0.0010475289339609_DP
         casimir_omega(69) = 0.0347690715573847_DP ;  casimir_omega_weight(69) = 0.0010655908081569_DP
         casimir_omega(70) = 0.0358437401493338_DP ;  casimir_omega_weight(70) = 0.0010837653339121_DP
         casimir_omega(71) = 0.0369366405102910_DP ;  casimir_omega_weight(71) = 0.0011020547158004_DP
         casimir_omega(72) = 0.0380478886097308_DP ;  casimir_omega_weight(72) = 0.0011204611849836_DP
         casimir_omega(73) = 0.0391776026618321_DP ;  casimir_omega_weight(73) = 0.0011389869998142_DP
         casimir_omega(74) = 0.0403259031529789_DP ;  casimir_omega_weight(74) = 0.0011576344464512_DP
         casimir_omega(75) = 0.0414929128698768_DP ;  casimir_omega_weight(75) = 0.0011764058394849_DP
         casimir_omega(76) = 0.0426787569283052_DP ;  casimir_omega_weight(76) = 0.0011953035225786_DP
         casimir_omega(77) = 0.0438835628025111_DP ;  casimir_omega_weight(77) = 0.0012143298691159_DP
         casimir_omega(78) = 0.0451074603552626_DP ;  casimir_omega_weight(78) = 0.0012334872828650_DP
         casimir_omega(79) = 0.0463505818685682_DP ;  casimir_omega_weight(79) = 0.0012527781986530_DP
         casimir_omega(80) = 0.0476130620750795_DP ;  casimir_omega_weight(80) = 0.0012722050830551_DP
         casimir_omega(81) = 0.0488950381901898_DP ;  casimir_omega_weight(81) = 0.0012917704350963_DP
         casimir_omega(82) = 0.0501966499448386_DP ;  casimir_omega_weight(82) = 0.0013114767869667_DP
         casimir_omega(83) = 0.0515180396190424_DP ;  casimir_omega_weight(83) = 0.0013313267047504_DP
         casimir_omega(84) = 0.0528593520761585_DP ;  casimir_omega_weight(84) = 0.0013513227891712_DP
         casimir_omega(85) = 0.0542207347979019_DP ;  casimir_omega_weight(85) = 0.0013714676763491_DP
         casimir_omega(86) = 0.0556023379201288_DP ;  casimir_omega_weight(86) = 0.0013917640385755_DP
         casimir_omega(87) = 0.0570043142694015_DP ;  casimir_omega_weight(87) = 0.0014122145851043_DP
         casimir_omega(88) = 0.0584268194003490_DP ;  casimir_omega_weight(88) = 0.0014328220629534_DP
         casimir_omega(89) = 0.0598700116338459_DP ;  casimir_omega_weight(89) = 0.0014535892577306_DP
         casimir_omega(90) = 0.0613340520960169_DP ;  casimir_omega_weight(90) = 0.0014745189944712_DP
         casimir_omega(91) = 0.0628191047580914_DP ;  casimir_omega_weight(91) = 0.0014956141384931_DP
         casimir_omega(92) = 0.0643253364771222_DP ;  casimir_omega_weight(92) = 0.0015168775962732_DP
         casimir_omega(93) = 0.0658529170375892_DP ;  casimir_omega_weight(93) = 0.0015383123163386_DP
         casimir_omega(94) = 0.0674020191939017_DP ;  casimir_omega_weight(94) = 0.0015599212901744_DP
         casimir_omega(95) = 0.0689728187138235_DP ;  casimir_omega_weight(95) = 0.0015817075531610_DP
         casimir_omega(96) = 0.0705654944228372_DP ;  casimir_omega_weight(96) = 0.0016036741855159_DP
         casimir_omega(97) = 0.0721802282494685_DP ;  casimir_omega_weight(97) = 0.0016258243132705_DP
         casimir_omega(98) = 0.0738172052715907_DP ;  casimir_omega_weight(98) = 0.0016481611092568_DP
         casimir_omega(99) = 0.0754766137637325_DP ;  casimir_omega_weight(99) = 0.0016706877941222_DP
         casimir_omega(100) = 0.0771586452454060_DP ;  casimir_omega_weight(100) = 0.0016934076373627_DP
         casimir_omega(101) = 0.0788634945304826_DP ;  casimir_omega_weight(101) = 0.0017163239583778_DP
         casimir_omega(102) = 0.0805913597776368_DP ;  casimir_omega_weight(102) = 0.0017394401275551_DP
         casimir_omega(103) = 0.0823424425418809_DP ;  casimir_omega_weight(103) = 0.0017627595673692_DP
         casimir_omega(104) = 0.0841169478272139_DP ;  casimir_omega_weight(104) = 0.0017862857535106_DP
         casimir_omega(105) = 0.0859150841404153_DP ;  casimir_omega_weight(105) = 0.0018100222160418_DP
         casimir_omega(106) = 0.0877370635459983_DP ;  casimir_omega_weight(106) = 0.0018339725405682_DP
         casimir_omega(107) = 0.0895831017223605_DP ;  casimir_omega_weight(107) = 0.0018581403694535_DP
         casimir_omega(108) = 0.0914534180191477_DP ;  casimir_omega_weight(108) = 0.0018825294030410_DP
         casimir_omega(109) = 0.0933482355158680_DP ;  casimir_omega_weight(109) = 0.0019071434009223_DP
         casimir_omega(110) = 0.0952677810817768_DP ;  casimir_omega_weight(110) = 0.0019319861832188_DP
         casimir_omega(111) = 0.0972122854370652_DP ;  casimir_omega_weight(111) = 0.0019570616319020_DP
         casimir_omega(112) = 0.0991819832153824_DP ;  casimir_omega_weight(112) = 0.0019823736921425_DP
         casimir_omega(113) = 0.1011771130277188_DP ;  casimir_omega_weight(113) = 0.0020079263736871_DP
         casimir_omega(114) = 0.1031979175276858_DP ;  casimir_omega_weight(114) = 0.0020337237522669_DP
         casimir_omega(115) = 0.1052446434782196_DP ;  casimir_omega_weight(115) = 0.0020597699710446_DP
         casimir_omega(116) = 0.1073175418197478_DP ;  casimir_omega_weight(116) = 0.0020860692420883_DP
         casimir_omega(117) = 0.1094168677398472_DP ;  casimir_omega_weight(117) = 0.0021126258478834_DP
         casimir_omega(118) = 0.1115428807444309_DP ;  casimir_omega_weight(118) = 0.0021394441428790_DP
         casimir_omega(119) = 0.1136958447305016_DP ;  casimir_omega_weight(119) = 0.0021665285550693_DP
         casimir_omega(120) = 0.1158760280605013_DP ;  casimir_omega_weight(120) = 0.0021938835876151_DP
         casimir_omega(121) = 0.1180837036383069_DP ;  casimir_omega_weight(121) = 0.0022215138205031_DP
         casimir_omega(122) = 0.1203191489868990_DP ;  casimir_omega_weight(122) = 0.0022494239122404_DP
         casimir_omega(123) = 0.1225826463277500_DP ;  casimir_omega_weight(123) = 0.0022776186015984_DP
         casimir_omega(124) = 0.1248744826619739_DP ;  casimir_omega_weight(124) = 0.0023061027093869_DP
         casimir_omega(125) = 0.1271949498532739_DP ;  casimir_omega_weight(125) = 0.0023348811402831_DP
         casimir_omega(126) = 0.1295443447127380_DP ;  casimir_omega_weight(126) = 0.0023639588846944_DP
         casimir_omega(127) = 0.1319229690855227_DP ;  casimir_omega_weight(127) = 0.0023933410206724_DP
         casimir_omega(128) = 0.1343311299394730_DP ;  casimir_omega_weight(128) = 0.0024230327158742_DP
         casimir_omega(129) = 0.1367691394557275_DP ;  casimir_omega_weight(129) = 0.0024530392295666_DP
         casimir_omega(130) = 0.1392373151213493_DP ;  casimir_omega_weight(130) = 0.0024833659146812_DP
         casimir_omega(131) = 0.1417359798240448_DP ;  casimir_omega_weight(131) = 0.0025140182199257_DP
         casimir_omega(132) = 0.1442654619490109_DP ;  casimir_omega_weight(132) = 0.0025450016919394_DP
         casimir_omega(133) = 0.1468260954779707_DP ;  casimir_omega_weight(133) = 0.0025763219775073_DP
         casimir_omega(134) = 0.1494182200904458_DP ;  casimir_omega_weight(134) = 0.0026079848258264_DP
         casimir_omega(135) = 0.1520421812673307_DP ;  casimir_omega_weight(135) = 0.0026399960908355_DP
         casimir_omega(136) = 0.1546983303968128_DP ;  casimir_omega_weight(136) = 0.0026723617335859_DP
         casimir_omega(137) = 0.1573870248827109_DP ;  casimir_omega_weight(137) = 0.0027050878246979_DP
         casimir_omega(138) = 0.1601086282552840_DP ;  casimir_omega_weight(138) = 0.0027381805468574_DP
         casimir_omega(139) = 0.1628635102845802_DP ;  casimir_omega_weight(139) = 0.0027716461973878_DP
         casimir_omega(140) = 0.1656520470963798_DP ;  casimir_omega_weight(140) = 0.0028054911908778_DP
         casimir_omega(141) = 0.1684746212908128_DP ;  casimir_omega_weight(141) = 0.0028397220618895_DP
         casimir_omega(142) = 0.1713316220637061_DP ;  casimir_omega_weight(142) = 0.0028743454677256_DP
         casimir_omega(143) = 0.1742234453307382_DP ;  casimir_omega_weight(143) = 0.0029093681912667_DP
         casimir_omega(144) = 0.1771504938544725_DP ;  casimir_omega_weight(144) = 0.0029447971438908_DP
         casimir_omega(145) = 0.1801131773743418_DP ;  casimir_omega_weight(145) = 0.0029806393684663_DP
         casimir_omega(146) = 0.1831119127396658_DP ;  casimir_omega_weight(146) = 0.0030169020424137_DP
         casimir_omega(147) = 0.1861471240457687_DP ;  casimir_omega_weight(147) = 0.0030535924808575_DP
         casimir_omega(148) = 0.1892192427732969_DP ;  casimir_omega_weight(148) = 0.0030907181398577_DP
         casimir_omega(149) = 0.1923287079307990_DP ;  casimir_omega_weight(149) = 0.0031282866197238_DP
         casimir_omega(150) = 0.1954759662006754_DP ;  casimir_omega_weight(150) = 0.0031663056684211_DP
         casimir_omega(151) = 0.1986614720885699_DP ;  casimir_omega_weight(151) = 0.0032047831850623_DP
         casimir_omega(152) = 0.2018856880763053_DP ;  casimir_omega_weight(152) = 0.0032437272234968_DP
         casimir_omega(153) = 0.2051490847784468_DP ;  casimir_omega_weight(153) = 0.0032831459959878_DP
         casimir_omega(154) = 0.2084521411026094_DP ;  casimir_omega_weight(154) = 0.0033230478770017_DP
         casimir_omega(155) = 0.2117953444135834_DP ;  casimir_omega_weight(155) = 0.0033634414070845_DP
         casimir_omega(156) = 0.2151791907014038_DP ;  casimir_omega_weight(156) = 0.0034043352968492_DP
         casimir_omega(157) = 0.2186041847534569_DP ;  casimir_omega_weight(157) = 0.0034457384310819_DP
         casimir_omega(158) = 0.2220708403307426_DP ;  casimir_omega_weight(158) = 0.0034876598729307_DP
         casimir_omega(159) = 0.2255796803483948_DP ;  casimir_omega_weight(159) = 0.0035301088682466_DP
         casimir_omega(160) = 0.2291312370605879_DP ;  casimir_omega_weight(160) = 0.0035730948500088_DP
         casimir_omega(161) = 0.2327260522499401_DP ;  casimir_omega_weight(161) = 0.0036166274428982_DP
         casimir_omega(162) = 0.2363646774215460_DP ;  casimir_omega_weight(162) = 0.0036607164679785_DP
         casimir_omega(163) = 0.2400476740017571_DP ;  casimir_omega_weight(163) = 0.0037053719475132_DP
         casimir_omega(164) = 0.2437756135418513_DP ;  casimir_omega_weight(164) = 0.0037506041099263_DP
         casimir_omega(165) = 0.2475490779267189_DP ;  casimir_omega_weight(165) = 0.0037964233948822_DP
         casimir_omega(166) = 0.2513686595887125_DP ;  casimir_omega_weight(166) = 0.0038428404585236_DP
         casimir_omega(167) = 0.2552349617268003_DP ;  casimir_omega_weight(167) = 0.0038898661788483_DP
         casimir_omega(168) = 0.2591485985311770_DP ;  casimir_omega_weight(168) = 0.0039375116612448_DP
         casimir_omega(169) = 0.2631101954134792_DP ;  casimir_omega_weight(169) = 0.0039857882441728_DP
         casimir_omega(170) = 0.2671203892427734_DP ;  casimir_omega_weight(170) = 0.0040347075050181_DP
         casimir_omega(171) = 0.2711798285874730_DP ;  casimir_omega_weight(171) = 0.0040842812661065_DP
         casimir_omega(172) = 0.2752891739633589_DP ;  casimir_omega_weight(172) = 0.0041345216008931_DP
         casimir_omega(173) = 0.2794490980878735_DP ;  casimir_omega_weight(173) = 0.0041854408403241_DP
         casimir_omega(174) = 0.2836602861408780_DP ;  casimir_omega_weight(174) = 0.0042370515793931_DP
         casimir_omega(175) = 0.2879234360320445_DP ;  casimir_omega_weight(175) = 0.0042893666838742_DP
         casimir_omega(176) = 0.2922392586750949_DP ;  casimir_omega_weight(176) = 0.0043423992972575_DP
         casimir_omega(177) = 0.2966084782690636_DP ;  casimir_omega_weight(177) = 0.0043961628478909_DP
         casimir_omega(178) = 0.3010318325868112_DP ;  casimir_omega_weight(178) = 0.0044506710563184_DP
         casimir_omega(179) = 0.3055100732709835_DP ;  casimir_omega_weight(179) = 0.0045059379428473_DP
         casimir_omega(180) = 0.3100439661376505_DP ;  casimir_omega_weight(180) = 0.0045619778353298_DP
         casimir_omega(181) = 0.3146342914878347_DP ;  casimir_omega_weight(181) = 0.0046188053771803_DP
         casimir_omega(182) = 0.3192818444271837_DP ;  casimir_omega_weight(182) = 0.0046764355356254_DP
         casimir_omega(183) = 0.3239874351940089_DP ;  casimir_omega_weight(183) = 0.0047348836102066_DP
         casimir_omega(184) = 0.3287518894959561_DP ;  casimir_omega_weight(184) = 0.0047941652415309_DP
         casimir_omega(185) = 0.3335760488555559_DP ;  casimir_omega_weight(185) = 0.0048542964202879_DP
         casimir_omega(186) = 0.3384607709649328_DP ;  casimir_omega_weight(186) = 0.0049152934965426_DP
         casimir_omega(187) = 0.3434069300499333_DP ;  casimir_omega_weight(187) = 0.0049771731892983_DP
         casimir_omega(188) = 0.3484154172439781_DP ;  casimir_omega_weight(188) = 0.0050399525963649_DP
         casimir_omega(189) = 0.3534871409719215_DP ;  casimir_omega_weight(189) = 0.0051036492045182_DP
         casimir_omega(190) = 0.3586230273442231_DP ;  casimir_omega_weight(190) = 0.0051682808999704_DP
         casimir_omega(191) = 0.3638240205617619_DP ;  casimir_omega_weight(191) = 0.0052338659791688_DP
         casimir_omega(192) = 0.3690910833316041_DP ;  casimir_omega_weight(192) = 0.0053004231599257_DP
         casimir_omega(193) = 0.3744251972940760_DP ;  casimir_omega_weight(193) = 0.0053679715928907_DP
         casimir_omega(194) = 0.3798273634614850_DP ;  casimir_omega_weight(194) = 0.0054365308733808_DP
         casimir_omega(195) = 0.3852986026688597_DP ;  casimir_omega_weight(195) = 0.0055061210535825_DP
         casimir_omega(196) = 0.3908399560370759_DP ;  casimir_omega_weight(196) = 0.0055767626551372_DP
         casimir_omega(197) = 0.3964524854487677_DP ;  casimir_omega_weight(197) = 0.0056484766821154_DP
         casimir_omega(198) = 0.4021372740374183_DP ;  casimir_omega_weight(198) = 0.0057212846344140_DP
         casimir_omega(199) = 0.4078954266900595_DP ;  casimir_omega_weight(199) = 0.0057952085215666_DP
         casimir_omega(200) = 0.4137280705639951_DP ;  casimir_omega_weight(200) = 0.0058702708770015_DP
         casimir_omega(201) = 0.4196363556180159_DP ;  casimir_omega_weight(201) = 0.0059464947727545_DP
         casimir_omega(202) = 0.4256214551585579_DP ;  casimir_omega_weight(202) = 0.0060239038346504_DP
         casimir_omega(203) = 0.4316845664012890_DP ;  casimir_omega_weight(203) = 0.0061025222579831_DP
         casimir_omega(204) = 0.4378269110486241_DP ;  casimir_omega_weight(204) = 0.0061823748236945_DP
         casimir_omega(205) = 0.4440497358836847_DP ;  casimir_omega_weight(205) = 0.0062634869150813_DP
         casimir_omega(206) = 0.4503543133812359_DP ;  casimir_omega_weight(206) = 0.0063458845350603_DP
         casimir_omega(207) = 0.4567419423361653_DP ;  casimir_omega_weight(207) = 0.0064295943239765_DP
         casimir_omega(208) = 0.4632139485100663_DP ;  casimir_omega_weight(208) = 0.0065146435780176_DP
         casimir_omega(209) = 0.4697716852965417_DP ;  casimir_omega_weight(209) = 0.0066010602682298_DP
         casimir_omega(210) = 0.4764165344058282_DP ;  casimir_omega_weight(210) = 0.0066888730601699_DP
         casimir_omega(211) = 0.4831499065694038_DP ;  casimir_omega_weight(211) = 0.0067781113341988_DP
         casimir_omega(212) = 0.4899732422652269_DP ;  casimir_omega_weight(212) = 0.0068688052064750_DP
         casimir_omega(213) = 0.4968880124643193_DP ;  casimir_omega_weight(213) = 0.0069609855506395_DP
         casimir_omega(214) = 0.5038957193993926_DP ;  casimir_omega_weight(214) = 0.0070546840202334_DP
         casimir_omega(215) = 0.5109978973562858_DP ;  casimir_omega_weight(215) = 0.0071499330718921_DP
         casimir_omega(216) = 0.5181961134889655_DP ;  casimir_omega_weight(216) = 0.0072467659893058_DP
         casimir_omega(217) = 0.5254919686589146_DP ;  casimir_omega_weight(217) = 0.0073452169080150_DP
         casimir_omega(218) = 0.5328870982997301_DP ;  casimir_omega_weight(218) = 0.0074453208410569_DP
         casimir_omega(219) = 0.5403831733078089_DP ;  casimir_omega_weight(219) = 0.0075471137054891_DP
         casimir_omega(220) = 0.5479819009600186_DP ;  casimir_omega_weight(220) = 0.0076506323498453_DP
         casimir_omega(221) = 0.5556850258592859_DP ;  casimir_omega_weight(221) = 0.0077559145825320_DP
         casimir_omega(222) = 0.5634943309090907_DP ;  casimir_omega_weight(222) = 0.0078629992012281_DP
         casimir_omega(223) = 0.5714116383178602_DP ;  casimir_omega_weight(223) = 0.0079719260233141_DP
         casimir_omega(224) = 0.5794388106343351_DP ;  casimir_omega_weight(224) = 0.0080827359173731_DP
         casimir_omega(225) = 0.5875777518149813_DP ;  casimir_omega_weight(225) = 0.0081954708358069_DP
         casimir_omega(226) = 0.5958304083246058_DP ;  casimir_omega_weight(226) = 0.0083101738486250_DP
         casimir_omega(227) = 0.6041987702713444_DP ;  casimir_omega_weight(227) = 0.0084268891784122_DP
         casimir_omega(228) = 0.6126848725772682_DP ;  casimir_omega_weight(228) = 0.0085456622365935_DP
         casimir_omega(229) = 0.6212907961858709_DP ;  casimir_omega_weight(229) = 0.0086665396609654_DP
         casimir_omega(230) = 0.6300186693077849_DP ;  casimir_omega_weight(230) = 0.0087895693546137_DP
         casimir_omega(231) = 0.6388706687061233_DP ;  casimir_omega_weight(231) = 0.0089148005262444_DP
         casimir_omega(232) = 0.6478490210228586_DP ;  casimir_omega_weight(232) = 0.0090422837319762_DP
         casimir_omega(233) = 0.6569560041477835_DP ;  casimir_omega_weight(233) = 0.0091720709186795_DP
         casimir_omega(234) = 0.6661939486315935_DP ;  casimir_omega_weight(234) = 0.0093042154689140_DP
         casimir_omega(235) = 0.6755652391447333_DP ;  casimir_omega_weight(235) = 0.0094387722475278_DP
         casimir_omega(236) = 0.6850723159837068_DP ;  casimir_omega_weight(236) = 0.0095757976499983_DP
         casimir_omega(237) = 0.6947176766266231_DP ;  casimir_omega_weight(237) = 0.0097153496525690_DP
         casimir_omega(238) = 0.7045038773398222_DP ;  casimir_omega_weight(238) = 0.0098574878642919_DP
         casimir_omega(239) = 0.7144335348375130_DP ;  casimir_omega_weight(239) = 0.0100022735810199_DP
         casimir_omega(240) = 0.7245093279964200_DP ;  casimir_omega_weight(240) = 0.0101497698414437_DP
         casimir_omega(241) = 0.7347339996275327_DP ;  casimir_omega_weight(241) = 0.0103000414852953_DP
         casimir_omega(242) = 0.7451103583071615_DP ;  casimir_omega_weight(242) = 0.0104531552137341_DP
         casimir_omega(243) = 0.7556412802695449_DP ;  casimir_omega_weight(243) = 0.0106091796521101_DP
         casimir_omega(244) = 0.7663297113634072_DP ;  casimir_omega_weight(244) = 0.0107681854151190_DP
         casimir_omega(245) = 0.7771786690749379_DP ;  casimir_omega_weight(245) = 0.0109302451744965_DP
         casimir_omega(246) = 0.7881912446197853_DP ;  casimir_omega_weight(246) = 0.0110954337293688_DP
         casimir_omega(247) = 0.7993706051067645_DP ;  casimir_omega_weight(247) = 0.0112638280793713_DP
         casimir_omega(248) = 0.8107199957761032_DP ;  casimir_omega_weight(248) = 0.0114355075006034_DP
         casimir_omega(249) = 0.8222427423151699_DP ;  casimir_omega_weight(249) = 0.0116105536246716_DP
         casimir_omega(250) = 0.8339422532547728_DP ;  casimir_omega_weight(250) = 0.0117890505208126_DP
         casimir_omega(251) = 0.8458220224492404_DP ;  casimir_omega_weight(251) = 0.0119710847813476_DP
         casimir_omega(252) = 0.8578856316436481_DP ;  casimir_omega_weight(252) = 0.0121567456105584_DP
         casimir_omega(253) = 0.8701367531317070_DP ;  casimir_omega_weight(253) = 0.0123461249171566_DP
         casimir_omega(254) = 0.8825791525079961_DP ;  casimir_omega_weight(254) = 0.0125393174105331_DP
         casimir_omega(255) = 0.8952166915183858_DP ;  casimir_omega_weight(255) = 0.0127364207009037_DP
         casimir_omega(256) = 0.9080533310126579_DP ;  casimir_omega_weight(256) = 0.0129375354036089_DP
         casimir_omega(257) = 0.9210931340035508_DP ;  casimir_omega_weight(257) = 0.0131427652476938_DP
         casimir_omega(258) = 0.9343402688366300_DP ;  casimir_omega_weight(258) = 0.0133522171889981_DP
         casimir_omega(259) = 0.9477990124755841_DP ;  casimir_omega_weight(259) = 0.0135660015279736_DP
         casimir_omega(260) = 0.9614737539077894_DP ;  casimir_omega_weight(260) = 0.0137842320324237_DP
         casimir_omega(261) = 0.9753689976751989_DP ;  casimir_omega_weight(261) = 0.0140070260654504_DP
         casimir_omega(262) = 0.9894893675358294_DP ;  casimir_omega_weight(262) = 0.0142345047187891_DP
         casimir_omega(263) = 1.0038396102614335_DP ;  casimir_omega_weight(263) = 0.0144667929518305_DP
         casimir_omega(264) = 1.0184245995771493_DP ;  casimir_omega_weight(264) = 0.0147040197366226_DP
         casimir_omega(265) = 1.0332493402492293_DP ;  casimir_omega_weight(265) = 0.0149463182090902_DP
         casimir_omega(266) = 1.0483189723272459_DP ;  casimir_omega_weight(266) = 0.0151938258268155_DP
         casimir_omega(267) = 1.0636387755474894_DP ;  casimir_omega_weight(267) = 0.0154466845336713_DP
         casimir_omega(268) = 1.0792141739045715_DP ;  casimir_omega_weight(268) = 0.0157050409316934_DP
         casimir_omega(269) = 1.0950507403986454_DP ;  casimir_omega_weight(269) = 0.0159690464604814_DP
         casimir_omega(270) = 1.1111542019659688_DP ;  casimir_omega_weight(270) = 0.0162388575845584_DP
         casimir_omega(271) = 1.1275304446009697_DP ;  casimir_omega_weight(271) = 0.0165146359890634_DP
         casimir_omega(272) = 1.1441855186783387_DP ;  casimir_omega_weight(272) = 0.0167965487842183_DP
         casimir_omega(273) = 1.1611256444841467_DP ;  casimir_omega_weight(273) = 0.0170847687189569_DP
         casimir_omega(274) = 1.1783572179653945_DP ;  casimir_omega_weight(274) = 0.0173794744042627_DP
         casimir_omega(275) = 1.1958868167079451_DP ;  casimir_omega_weight(275) = 0.0176808505466580_DP
         casimir_omega(276) = 1.2137212061532232_DP ;  casimir_omega_weight(276) = 0.0179890881923750_DP
         casimir_omega(277) = 1.2318673460646843_DP ;  casimir_omega_weight(277) = 0.0183043849827604_DP
         casimir_omega(278) = 1.2503323972555442_DP ;  casimir_omega_weight(278) = 0.0186269454215564_DP
         casimir_omega(279) = 1.2691237285899406_DP ;  casimir_omega_weight(279) = 0.0189569811545535_DP
         casimir_omega(280) = 1.2882489242702695_DP ;  casimir_omega_weight(280) = 0.0192947112624053_DP
         casimir_omega(281) = 1.3077157914241726_DP ;  casimir_omega_weight(281) = 0.0196403625672238_DP
         casimir_omega(282) = 1.3275323680053024_DP ;  casimir_omega_weight(282) = 0.0199941699537034_DP
         casimir_omega(283) = 1.3477069310228222_DP ;  casimir_omega_weight(283) = 0.0203563767055644_DP
         casimir_omega(284) = 1.3682480051153385_DP ;  casimir_omega_weight(284) = 0.0207272348581198_DP
         casimir_omega(285) = 1.3891643714858437_DP ;  casimir_omega_weight(285) = 0.0211070055678982_DP
         casimir_omega(286) = 1.4104650772151579_DP ;  casimir_omega_weight(286) = 0.0214959595001682_DP
         casimir_omega(287) = 1.4321594449722932_DP ;  casimir_omega_weight(287) = 0.0218943772354121_DP
         casimir_omega(288) = 1.4542570831411832_DP ;  casimir_omega_weight(288) = 0.0223025496957969_DP
         casimir_omega(289) = 1.4767678963843325_DP ;  casimir_omega_weight(289) = 0.0227207785927100_DP
         casimir_omega(290) = 1.4997020966650321_DP ;  casimir_omega_weight(290) = 0.0231493768965655_DP
         casimir_omega(291) = 1.5230702147510633_DP ;  casimir_omega_weight(291) = 0.0235886693301570_DP
         casimir_omega(292) = 1.5468831122240299_DP ;  casimir_omega_weight(292) = 0.0240389928868661_DP
         casimir_omega(293) = 1.5711519940199501_DP ;  casimir_omega_weight(293) = 0.0245006973751479_DP
         casimir_omega(294) = 1.5958884215280424_DP ;  casimir_omega_weight(294) = 0.0249741459907523_DP
         casimir_omega(295) = 1.6211043262763905_DP ;  casimir_omega_weight(295) = 0.0254597159184351_DP
         casimir_omega(296) = 1.6468120242346240_DP ;  casimir_omega_weight(296) = 0.0259577989646950_DP
         casimir_omega(297) = 1.6730242307656935_DP ;  casimir_omega_weight(297) = 0.0264688022234649_DP
         casimir_omega(298) = 1.6997540762605698_DP ;  casimir_omega_weight(298) = 0.0269931487766379_DP
         casimir_omega(299) = 1.7270151224917940_DP ;  casimir_omega_weight(299) = 0.0275312784315345_DP
         casimir_omega(300) = 1.7548213797238545_DP ;  casimir_omega_weight(300) = 0.0280836484975632_DP
         casimir_omega(301) = 1.7831873246207246_DP ;  casimir_omega_weight(301) = 0.0286507346042322_DP
         casimir_omega(302) = 1.8121279189932697_DP ;  casimir_omega_weight(302) = 0.0292330315632747_DP
         casimir_omega(303) = 1.8416586294318362_DP ;  casimir_omega_weight(303) = 0.0298310542774390_DP
         casimir_omega(304) = 1.8717954478721430_DP ;  casimir_omega_weight(304) = 0.0304453386987080_DP
         casimir_omega(305) = 1.9025549131454635_DP ;  casimir_omega_weight(305) = 0.0310764428392670_DP
         casimir_omega(306) = 1.9339541335674044_DP ;  casimir_omega_weight(306) = 0.0317249478382098_DP
         casimir_omega(307) = 1.9660108106227909_DP ;  casimir_omega_weight(307) = 0.0323914590876869_DP
         casimir_omega(308) = 1.9987432638078821_DP ;  casimir_omega_weight(308) = 0.0330766074221293_DP
         casimir_omega(309) = 2.0321704566950705_DP ;  casimir_omega_weight(309) = 0.0337810503745592_DP
         casimir_omega(310) = 2.0663120242892226_DP ;  casimir_omega_weight(310) = 0.0345054735043614_DP
         casimir_omega(311) = 2.1011883017493469_DP ;  casimir_omega_weight(311) = 0.0352505918009904_DP
         casimir_omega(312) = 2.1368203545540854_DP ;  casimir_omega_weight(312) = 0.0360171511687121_DP
         casimir_omega(313) = 2.1732300101944686_DP ;  casimir_omega_weight(313) = 0.0368059299973809_DP
         casimir_omega(314) = 2.2104398914830812_DP ;  casimir_omega_weight(314) = 0.0376177408253030_DP
         casimir_omega(315) = 2.2484734515743958_DP ;  casimir_omega_weight(315) = 0.0384534320999312_DP
         casimir_omega(316) = 2.2873550107976355_DP ;  casimir_omega_weight(316) = 0.0393138900432321_DP
         casimir_omega(317) = 2.3271097954100322_DP ;  casimir_omega_weight(317) = 0.0402000406284584_DP
         casimir_omega(318) = 2.3677639783858622_DP ;  casimir_omega_weight(318) = 0.0411128516761256_DP
         casimir_omega(319) = 2.4093447223644135_DP ;  casimir_omega_weight(319) = 0.0420533350772816_DP
         casimir_omega(320) = 2.4518802248883169_DP ;  casimir_omega_weight(320) = 0.0430225491526560_DP
         casimir_omega(321) = 2.4953997660731684_DP ;  casimir_omega_weight(321) = 0.0440216011574862_DP
         casimir_omega(322) = 2.5399337588586519_DP ;  casimir_omega_weight(322) = 0.0450516499416826_DP
         casimir_omega(323) = 2.5855138020023500_DP ;  casimir_omega_weight(323) = 0.0461139087767946_DP
         casimir_omega(324) = 2.6321727359886795_DP ;  casimir_omega_weight(324) = 0.0472096483612504_DP
         casimir_omega(325) = 2.6799447020374836_DP ;  casimir_omega_weight(325) = 0.0483402000169007_DP
         casimir_omega(326) = 2.7288652044104702_DP ;  casimir_omega_weight(326) = 0.0495069590904133_DP
         casimir_omega(327) = 2.7789711762276759_DP ;  casimir_omega_weight(327) = 0.0507113885744640_DP
         casimir_omega(328) = 2.8303010490218341_DP ;  casimir_omega_weight(328) = 0.0519550229651951_DP
         casimir_omega(329) = 2.8828948262752285_DP ;  casimir_omega_weight(329) = 0.0532394723728185_DP
         casimir_omega(330) = 2.9367941612017896_DP ;  casimir_omega_weight(330) = 0.0545664269046126_DP
         casimir_omega(331) = 2.9920424390568070_DP ;  casimir_omega_weight(331) = 0.0559376613408838_DP
         casimir_omega(332) = 3.0486848642780404_DP ;  casimir_omega_weight(332) = 0.0573550401255898_DP
         casimir_omega(333) = 3.1067685527850162_DP ;  casimir_omega_weight(333) = 0.0588205226961599_DP
         casimir_omega(334) = 3.1663426297884798_DP ;  casimir_omega_weight(334) = 0.0603361691786479_DP
         casimir_omega(335) = 3.2274583334891367_DP ;  casimir_omega_weight(335) = 0.0619041464760857_DP
         casimir_omega(336) = 3.2901691250744496_DP ;  casimir_omega_weight(336) = 0.0635267347815901_DP
         casimir_omega(337) = 3.3545308054545364_DP ;  casimir_omega_weight(337) = 0.0652063345491959_DP
         casimir_omega(338) = 3.4206016392131269_DP ;  casimir_omega_weight(338) = 0.0669454739595680_DP
         casimir_omega(339) = 3.4884424862877581_DP ;  casimir_omega_weight(339) = 0.0687468169197903_DP
         casimir_omega(340) = 3.5581169419350021_DP ;  casimir_omega_weight(340) = 0.0706131716410667_DP
         casimir_omega(341) = 3.6296914855818887_DP ;  casimir_omega_weight(341) = 0.0725474998416722_DP
         casimir_omega(342) = 3.7032356392141565_DP ;  casimir_omega_weight(342) = 0.0745529266266362_DP
         casimir_omega(343) = 3.7788221360060890_DP ;  casimir_omega_weight(343) = 0.0766327511012142_DP
         casimir_omega(344) = 3.8565270999557861_DP ;  casimir_omega_weight(344) = 0.0787904577790096_DP
         casimir_omega(345) = 3.9364302373545712_DP ;  casimir_omega_weight(345) = 0.0810297288537672_DP
         casimir_omega(346) = 4.0186150409898334_DP ;  casimir_omega_weight(346) = 0.0833544574075836_DP
         casimir_omega(347) = 4.1031690080584964_DP ;  casimir_omega_weight(347) = 0.0857687616379527_DP
         casimir_omega(348) = 4.1901838728532654_DP ;  casimir_omega_weight(348) = 0.0882770001920495_DP
         casimir_omega(349) = 4.2797558553774664_DP ;  casimir_omega_weight(349) = 0.0908837887065266_DP
         casimir_omega(350) = 4.3719859271468371_DP ;  casimir_omega_weight(350) = 0.0935940176606338_DP
         casimir_omega(351) = 4.4669800955498413_DP ;  casimir_omega_weight(351) = 0.0964128716605498_DP
         casimir_omega(352) = 4.5648497082620763_DP ;  casimir_omega_weight(352) = 0.0993458502863521_DP
         casimir_omega(353) = 4.6657117793476948_DP ;  casimir_omega_weight(353) = 0.1023987906442831_DP
         casimir_omega(354) = 4.7696893388318005_DP ;  casimir_omega_weight(354) = 0.1055778917836082_DP
         casimir_omega(355) = 4.8769118076940963_DP ;  casimir_omega_weight(355) = 0.1088897411537538_DP
         casimir_omega(356) = 4.9875154004192526_DP ;  casimir_omega_weight(356) = 0.1123413432939822_DP
         casimir_omega(357) = 5.1016435574423280_DP ;  casimir_omega_weight(357) = 0.1159401509716625_DP
         casimir_omega(358) = 5.2194474100541637_DP ;  casimir_omega_weight(358) = 0.1196940990053071_DP
         casimir_omega(359) = 5.3410862805809236_DP ;  casimir_omega_weight(359) = 0.1236116410366123_DP
         casimir_omega(360) = 5.4667282209307588_DP ;  casimir_omega_weight(360) = 0.1277017895430910_DP
         casimir_omega(361) = 5.5965505929072927_DP ;  casimir_omega_weight(361) = 0.1319741594171508_DP
         casimir_omega(362) = 5.7307406940343393_DP ;  casimir_omega_weight(362) = 0.1364390154721324_DP
         casimir_omega(363) = 5.8694964330161730_DP ;  casimir_omega_weight(363) = 0.1411073242790254_DP
         casimir_omega(364) = 6.0130270593848172_DP ;  casimir_omega_weight(364) = 0.1459908107821974_DP
         casimir_omega(365) = 6.1615539523592258_DP ;  casimir_omega_weight(365) = 0.1511020201959771_DP
         casimir_omega(366) = 6.3153114744727104_DP ;  casimir_omega_weight(366) = 0.1564543857431390_DP
         casimir_omega(367) = 6.4745478961177598_DP ;  casimir_omega_weight(367) = 0.1620623028620961_DP
         casimir_omega(368) = 6.6395263978224435_DP ;  casimir_omega_weight(368) = 0.1679412105870755_DP
         casimir_omega(369) = 6.8105261578182086_DP ;  casimir_omega_weight(369) = 0.1741076808896208_DP
         casimir_omega(370) = 6.9878435332957221_DP ;  casimir_omega_weight(370) = 0.1805795168699250_DP
         casimir_omega(371) = 7.1717933446874609_DP ;  casimir_omega_weight(371) = 0.1873758607957537_DP
         casimir_omega(372) = 7.3627102733755496_DP ;  casimir_omega_weight(372) = 0.1945173131131842_DP
         casimir_omega(373) = 7.5609503844202521_DP ;  casimir_omega_weight(373) = 0.2020260637011178_DP
         casimir_omega(374) = 7.7668927872537710_DP ;  casimir_omega_weight(374) = 0.2099260368022720_DP
         casimir_omega(375) = 7.9809414488128576_DP ;  casimir_omega_weight(375) = 0.2182430512572827_DP
         casimir_omega(376) = 8.2035271753140204_DP ;  casimir_omega_weight(376) = 0.2270049978831682_DP
         casimir_omega(377) = 8.4351097808390829_DP ;  casimir_omega_weight(377) = 0.2362420360863420_DP
         casimir_omega(378) = 8.6761804631305051_DP ;  casimir_omega_weight(378) = 0.2459868120915991_DP
         casimir_omega(379) = 8.9272644095364093_DP ;  casimir_omega_weight(379) = 0.2562747014951714_DP
         casimir_omega(380) = 9.1889236589435495_DP ;  casimir_omega_weight(380) = 0.2671440792373651_DP
         casimir_omega(381) = 9.4617602488440795_DP ;  casimir_omega_weight(381) = 0.2786366205283061_DP
         casimir_omega(382) = 9.7464196804714192_DP ;  casimir_omega_weight(382) = 0.2907976367812569_DP
         casimir_omega(383) = 10.0435947392809641_DP ;  casimir_omega_weight(383) = 0.3036764512002945_DP
         casimir_omega(384) = 10.3540297130407453_DP ;  casimir_omega_weight(384) = 0.3173268193656384_DP
         casimir_omega(385) = 10.6785250555362143_DP ;  casimir_omega_weight(385) = 0.3318074009776099_DP
         casimir_omega(386) = 11.0179425505134105_DP ;  casimir_omega_weight(386) = 0.3471822898621595_DP
         casimir_omega(387) = 11.3732110381363576_DP ;  casimir_omega_weight(387) = 0.3635216104668109_DP
         casimir_omega(388) = 11.7453327750922369_DP ;  casimir_omega_weight(388) = 0.3809021903693214_DP
         casimir_omega(389) = 12.1353905097605903_DP ;  casimir_omega_weight(389) = 0.3994083198845556_DP
         casimir_omega(390) = 12.5445553658232800_DP ;  casimir_omega_weight(390) = 0.4191326116595547_DP
         casimir_omega(391) = 12.9740956416385984_DP ;  casimir_omega_weight(391) = 0.4401769753130819_DP
         casimir_omega(392) = 13.4253866489964935_DP ;  casimir_omega_weight(392) = 0.4626537247346959_DP
         casimir_omega(393) = 13.8999217339733239_DP ;  casimir_omega_weight(393) = 0.4866868387096482_DP
         casimir_omega(394) = 14.3993246450310206_DP ;  casimir_omega_weight(394) = 0.5124133991715614_DP
         casimir_omega(395) = 14.9253634399177280_DP ;  casimir_omega_weight(395) = 0.5399852357427026_DP
         casimir_omega(396) = 15.4799661541337414_DP ;  casimir_omega_weight(396) = 0.5695708104534924_DP
         casimir_omega(397) = 16.0652384906587287_DP ;  casimir_omega_weight(397) = 0.6013573828383335_DP
         casimir_omega(398) = 16.6834838345237486_DP ;  casimir_omega_weight(398) = 0.6355535032143199_DP
         casimir_omega(399) = 17.3372259480854645_DP ;  casimir_omega_weight(399) = 0.6723918911793137_DP
         casimir_omega(400) = 18.0292347653320135_DP ;  casimir_omega_weight(400) = 0.7121327676130629_DP
         casimir_omega(401) = 18.7625557784450194_DP ;  casimir_omega_weight(401) = 0.7550677221538220_DP
         casimir_omega(402) = 19.5405435999543116_DP ;  casimir_omega_weight(402) = 0.8015242149560359_DP
         casimir_omega(403) = 20.3669003925772358_DP ;  casimir_omega_weight(403) = 0.8518708321813386_DP
         casimir_omega(404) = 21.2457199906252399_DP ;  casimir_omega_weight(404) = 0.9065234402317066_DP
         casimir_omega(405) = 22.1815386971102164_DP ;  casimir_omega_weight(405) = 0.9659524153844690_DP
         casimir_omega(406) = 23.1793939363506922_DP ;  casimir_omega_weight(406) = 1.0306911649467847_DP
         casimir_omega(407) = 24.2448921817275611_DP ;  casimir_omega_weight(407) = 1.1013462053770973_DP
         casimir_omega(408) = 25.3842878735319495_DP ;  casimir_omega_weight(408) = 1.1786091248958208_DP
         casimir_omega(409) = 26.6045754069879159_DP ;  casimir_omega_weight(409) = 1.2632708364300018_DP
         casimir_omega(410) = 27.9135967241476592_DP ;  casimir_omega_weight(410) = 1.3562386262681383_DP
         casimir_omega(411) = 29.3201676096094239_DP ;  casimir_omega_weight(411) = 1.4585566307777269_DP
         casimir_omega(412) = 30.8342265004445615_DP ;  casimir_omega_weight(412) = 1.5714305365458352_DP
         casimir_omega(413) = 32.4670105167788989_DP ;  casimir_omega_weight(413) = 1.6962575097511718_DP
         casimir_omega(414) = 34.2312645559225501_DP ;  casimir_omega_weight(414) = 1.8346626339102245_DP
         casimir_omega(415) = 36.1414907426949412_DP ;  casimir_omega_weight(415) = 1.9885434924814325_DP
         casimir_omega(416) = 38.2142473893200219_DP ;  casimir_omega_weight(416) = 2.1601250030477273_DP
         casimir_omega(417) = 40.4685090218598020_DP ;  casimir_omega_weight(417) = 2.3520272331618925_DP
         casimir_omega(418) = 42.9261021560708613_DP ;  casimir_omega_weight(418) = 2.5673497601590309_DP
         casimir_omega(419) = 45.6122355999315801_DP ;  casimir_omega_weight(419) = 2.8097772574598165_DP
         casimir_omega(420) = 48.5561494634656441_DP ;  casimir_omega_weight(420) = 3.0837125098926563_DP
         casimir_omega(421) = 51.7919142444527409_DP ;  casimir_omega_weight(421) = 3.3944451416651189_DP
         casimir_omega(422) = 55.3594210014377666_DP ;  casimir_omega_weight(422) = 3.7483672160278290_DP
         casimir_omega(423) = 59.3056166780783087_DP ;  casimir_omega_weight(423) = 4.1532508784443740_DP
         casimir_omega(424) = 63.6860564797935425_DP ;  casimir_omega_weight(424) = 4.6186088739540514_DP
         casimir_omega(425) = 68.5668698273859576_DP ;  casimir_omega_weight(425) = 5.1561668391565485_DP
         casimir_omega(426) = 74.0272707747543137_DP ;  casimir_omega_weight(426) = 5.7804879155234676_DP
         casimir_omega(427) = 80.1627922904136909_DP ;  casimir_omega_weight(427) = 6.5098072548538148_DP
         casimir_omega(428) = 87.0894931480726910_DP ;  casimir_omega_weight(428) = 7.3671592158429755_DP
         casimir_omega(429) = 94.9494866316628929_DP ;  casimir_omega_weight(429) = 8.3819179926952430_DP
         casimir_omega(430) = 103.9182879093186216_DP ;  casimir_omega_weight(430) = 9.5919303944970569_DP
         casimir_omega(431) = 114.2146973478982375_DP ;  casimir_omega_weight(431) = 11.0465096261161975_DP
         casimir_omega(432) = 126.1142717093221250_DP ;  casimir_omega_weight(432) = 12.8107016846218915_DP
         casimir_omega(433) = 139.9679527796282628_DP ;  casimir_omega_weight(433) = 14.9714667687998730_DP
         casimir_omega(434) = 156.2282398573049136_DP ;  casimir_omega_weight(434) = 17.6467995784739315_DP
         casimir_omega(435) = 175.4866106086338107_DP ;  casimir_omega_weight(435) = 20.9994586045209282_DP
         casimir_omega(436) = 198.5280743324144623_DP ;  casimir_omega_weight(436) = 25.2580992385450465_DP
         casimir_omega(437) = 226.4124450004896971_DP ;  casimir_omega_weight(437) = 30.7506228255868734_DP
         casimir_omega(438) = 260.5984078704438502_DP ;  casimir_omega_weight(438) = 37.9582955855186412_DP
         casimir_omega(439) = 303.1382091980196378_DP ;  casimir_omega_weight(439) = 47.6063993150359579_DP
         casimir_omega(440) = 356.9929449872076930_DP ;  casimir_omega_weight(440) = 60.8216694022915121_DP
         casimir_omega(441) = 426.5620346813580568_DP ;  casimir_omega_weight(441) = 79.4173863468914192_DP
         casimir_omega(442) = 518.6108511446556122_DP ;  casimir_omega_weight(442) = 106.4354197205271504_DP
         casimir_omega(443) = 643.9793798530039339_DP ;  casimir_omega_weight(443) = 147.2380601665370250_DP
         casimir_omega(444) = 820.9248017849427015_DP ;  casimir_omega_weight(444) = 211.8655963622535126_DP
         casimir_omega(445) = 1082.1617099692987267_DP ;  casimir_omega_weight(445) = 320.5812437482459245_DP
         casimir_omega(446) = 1491.1354123557941875_DP ;  casimir_omega_weight(446) = 518.3991517563925981_DP
         casimir_omega(447) = 2184.4821124596705886_DP ;  casimir_omega_weight(447) = 918.9239983172332131_DP
         casimir_omega(448) = 3502.7631513924575302_DP ;  casimir_omega_weight(448) = 1865.0216140612294566_DP
         casimir_omega(449) = 6503.8097575365391094_DP ;  casimir_omega_weight(449) = 4714.7753667403558211_DP
         casimir_omega(450) = 15984.5235832011239836_DP ;  casimir_omega_weight(450) = 18123.6657651843961503_DP
         casimir_omega(451) = 84223.2101565876073437_DP ;  casimir_omega_weight(451) = 216144.9819604934600648_DP
return
endsubroutine gauss_legendre_grid450

end module quadrature_grid_module
