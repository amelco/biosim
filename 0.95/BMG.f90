program BMG
use variables
implicit none

!-- FORMATS
100 FORMAT (2A7,11(A8,2x)) 
101 FORMAT (A5,10(A8,2x)) 
200 FORMAT (2A5,A8,  2x,A12,  2x,A6,  2x,2A12,  2x,A8)
201 FORMAT (2I5,F8.5,2x,F12.1,2x,F6.1,2x,2F12.1,2x,F8.2)

! Variable initialization
call initializeVariables()

! read soil parameters data from input file
!call readInputs()
call readSoilData()
print*, Os1

! opens the output file
Open(1,file='swb.out', status='replace')
!write(*,101) "Year", "ET", "Eta", "Tr", "D", "Rain", "BL"
write(*,101) "Year", "SQD", "ET", "Eta", "Tr", "Wi", "Wx", "D", "Rain", "BL", "omega"
write(1,101) "Year", "SQD", "ET", "Eta", "Tr", "Wi", "Wx", "D", "Rain", "BL", "omega"
Open(2,file='bio.out', status='replace')
write(2,200) "Year", "SQD", "A", "B", "N", "Broot", "Bshoot", "rd"


! Calculates maximum soil-water storage for each depth

Open(21,file='data.prn',status='old')
!Open(21,file='Peter2.prn',status='old')
Read(21,*)
Read(21,*)

do while (eof .ne. -1)

  ! calc. root growth 
	if (i .eq. 1) then
	  root_depth = rdg * Broot
		Wx_ant_t = root_depth * (Os1-Or1)
	else
	  root_depth = root_depth + rootGrowth(reduction_fac(Wi_ant_t, Wx_ant_t, beta), dBr)
  endif


  call calcWx(Or1, Os1, z1, Or2, Os2, z2, Wx, Wxt)
	

  ! calculates Leaf Area / shoot
  !A  =  LAperShoot(rho, Am, Aant, tau)
	! Temporary using Peter's equation for A
  A = LAperShoot2(rho, Am, Aant, tau, Wi_ant_t, Wxt, beta)

  ! calc. num. of shoots
  N = numShoots(dBs, Nant, Bshoot_ant, Nmax)

  ! calc area of tree
  Lt = A*N

  ! calculates interception of rad.
  ft = 1.0 - exp(-kt*Lt)
	
	! calculates crop coefficient
	kc = ft/ft_r
	! Everton's thesis: has to vary albedo to work
	!LAI = Lt *rho	
	!kc = (1.0-av-(as-av)*exp(-kt*LAI))/(1.0-0.23)
	!print*, LAI, kc
	!read*


  ! reads weather input file
  read(21,*, iostat=eof) Year, SQD, water_in, Rad, ET

  call waterBalance()

  ! calculates total, shoot and root biomass rates
  dB = biomass(radExt, ft, Rad, rho, bioloss_fac, Bant)
	dBs = (1.0 + reduction_fac(Wi_ant_t, Wx_ant_t, beta) * root_frac / (1.0 - root_frac)) * (1.0 - root_frac) * dB
	dBr = dB * root_frac
  ! calculates total, shoot and root biomass
  B = db + Bant
  Bshoot = Bshoot_ant + dBs
  Broot = Broot_ant + dBr

	write(2,201) Year, SQD, A, B, N, Broot, Bshoot, root_depth
  
  Aant = A
  Bant = B
	Bshoot_ant = Bshoot
	Broot_ant = Broot
  Nant = N
	i = i+1
Enddo

print*
print*, "Lost water: ", lost_water
print*

contains

! read soil parameters data from input file
subroutine readSoilData()
  implicit none
  Open(31,file='soil.in',status='old')
  !Open(31,file='soilPeter.in',status='old')
  read(31,*)
  read(31,*) Or1, Os1, alpha1, n1, Ks1
  read(31,*) Or2, Os2, alpha2, n2, Ks2
  read(31,*)
  read(31,*) z1, z2
  close(31)
end subroutine

! calculates maximum soil water storage
subroutine calcWx(Or1, Os1, z1, Or2, Os2, z2, Wx, Wxt)
  implicit none
  real, intent(in)  		  :: Or1, Os1, z1, Or2, Os2, z2
  real, intent(out), dimension(2) :: Wx
  real, intent(out) :: Wxt

  Wx(1) = (Os1 - Or1) * root_depth
	Wx(1) = min(Wx(1), z1*(Os1-Or1))
	if (root_depth > z1) then
    Wx(2) = (Os2 - Or2) * (root_depth - z1)
		Wx(2) = min(Wx(2), z1*(Os1-Or1) + z2*(Os2-Or2))
	else
	  Wx(2) = 0.0
	endif
	Wxt = Wx(1) + Wx(2)
end subroutine

subroutine calcWz(Or1, Os1, z1, Or2, Os2, root_depth, Wz)
  implicit none
  real, intent(in)  :: Or1, Os1, z1, Or2, Os2, root_depth
  real, intent(out) :: Wz

	if (root_depth > z1) then
	  Wz = ((Os1 - Or1) * z1) + ((Os2 - Or2) * (root_depth - z1))
	else
    Wz = (Os1 - Or1) * root_depth
	endif
end subroutine


subroutine tempSWB(water_in)
  implicit none
  real :: water_in

  !ETa = reduction_fac(Wi_ant_t, Wx_ant_t, beta) * (ET_ant)
  ETa = reduction_fac(Wi_ant_t, Wx_ant_t, beta) * ET * kc
  g = water_in * (1.0 - Lt*inter) - ETa
  
  if (g .ge. (Wxt - Wi_ant_t)) then
    dW = Wxt - Wi_ant_t
  else if (g .lt. -Wi_ant_t) then
    dW = -Wi_ant_t 
  else
    dW = g
  endif
  
  Wit = dW + Wi_ant_t 
  
  if (g .ge. (Wxt - Wi_ant_t)) then
    D = g + Wi_ant_t - Wxt
  Else
    D = 0.0
  Endif
end subroutine

! Equation from paper Werf et al. (2007)
real function LAperShoot(rho, Am, Ant, tau)
  implicit none
  real, intent(in) :: rho, Am, Ant, tau
  real :: dA

  dA = rho * (Am - Ant) / tau
  LAperShoot = Ant + dA
  return 
end function

real function reduction_fac(Wit, Wxt, beta)
  implicit none
	integer :: i
	real, intent(in) :: Wxt
	real, intent(inout) :: Wit
	real, intent(in) :: beta

	  !if (Wit > Wxt) then
		!  lost_water = lost_water + (Wit-Wxt)
		!	Wit = Wxt
		!endif

  reduction_fac = (Wit / Wxt)**beta
	return
end function

! Equation from Peter's spreadsheet
real function LAperShoot2(rho, Am, Ant, tau, Wit, Wxt, beta)
  implicit none
	real, intent(in) :: Wxt
	real, intent(inout) :: Wit
  real, intent(in) :: rho, Am, Ant, tau, beta
  real :: dA
	real, parameter :: afac = 1.0

  dA = rho * (Am - Ant) / tau - afac * (1.0 - reduction_fac(Wi_ant_t, Wx_ant_t, beta)) * Ant
  LAperShoot2 = Ant + dA
  return 
end function

real function numShoots(dBshoot, Nant, Bshoot, Nmax)
  implicit none
  real, intent(in) :: dBshoot, Nant, Bshoot, Nmax
  real :: dN

  dN = dBshoot * Nant / Bshoot * (1.0 - Nant/Nmax)
  numShoots = dN + Nant
  return
end function

real function biomass(radExt, ft, rad, rho, bioloss_fac, B)
  real, intent(in) :: radExt, ft, rad, rho, bioloss_fac, B
	real::alpha,bb,gama
	alpha = 0
	bb = 0.9
	gama = 0.000046
  
	if (B > (ft*radExt*reduction_fac(Wi_ant_t, Wx_ant_t, beta)*17.0/(bioloss_fac*rho))) then
	  biomass = 0.0
	else
    biomass = (reduction_fac(Wi_ant_t, Wx_ant_t, beta) * radExt*ft*rad ) / rho - (bioloss_fac * B)
  !biomass = (reduction_fac(Wi_ant_t, Wx_ant_t, beta) * radExt*ft*rad ) / rho * ( alpha+((2.0*bb)/(1.0+exp(gama*B))))
	endif

  return
end function

real function rootGrowth(omega, dBr)
implicit none
  real, intent(in) :: omega, dBr

	rootGrowth = (1.0 - omega) * rdg * dBr
	return
end function

subroutine initializeVariables()
  implicit none

  ! Variable initializationa
	i						= 1
  Wii       	= 0.0
	Wi_ant_t		= wii
	Wx_ant_t		= 0.00001
  Eta       	= 0.0
  ET_ant    	= 5.65
  beta 				= .8
  eof       	= 0
  Dsum      	= 0
  rainsum   	= 0
  BLsum     	= 0
  ETavg     	= 0
  ETaavg    	= 0
  Am 					= 0.05
  Aant 				= 0.0
  A 					= 0.0   !1/day
  rho 				= 0.0156
  tau 				= 10.0  ! day
  Nant 				= 0.6
  Nmax 				= 10000.0 
  dBs 				= 0.0
  Bant 				= 500.0
  Bshoot 			= Bant*0.4
	Broot				= Bant*0.6
	Bshoot_ant	= 200.0
	Broot_ant		= 300.0
  kt 					= 0.8
  radExt 			= 1.5
  bioloss_fac = 0.005
	root_frac   = 0.6
	root_max    = 800.0  ! mm
  Wz					= 0.0
	rdg 				= 0.2		 ! mm/g
  inter				= 0.13
	ft_r				= 0.846184277466571
	av					= 0.13
	as					= 0.25
end subroutine

subroutine waterBalance()
implicit none

110 FORMAT (2I7,9(F8.2,2x)) 

  ! calculates soil water balance
  if (eof .ne. -1) then
    call tempSWB(water_in)
    BL = D + dW + ETa - water_in * (1.0 - Lt*inter)

    ! show on screen at yearly timestep
    if (int(Year/4.0)*4 == Year) then
      if (SQD == 366) then
        write(*,110) Year, SQD, ET, Eta, ETa/ET_ant, Wit, Wxt, D, water_in, BL, &
				& reduction_fac(Wi_ant_t, Wx_ant_t, beta)
        Dsum      = 0
        rainsum   = 0
        BLsum     = 0
        ETavg     = 0
        ETaavg    = 0
     endif
    else
      if (SQD == 365) then
        write(*,110) Year, SQD, ET, Eta, ETa/ET_ant, Wit, Wxt, D, water_in, BL, &
				& reduction_fac(Wi_ant_t, Wx_ant_t, beta)
        Dsum      = 0
        rainsum   = 0
        BLsum     = 0
        ETavg     = 0
        ETaavg    = 0
      endif
    endif
    ! write in output file in daily timestep
    write(1,110) Year, SQD, ET, Eta, ETa/ET_ant, Wit, Wxt, D, water_in, BL, &
		& reduction_fac(Wi_ant_t, Wx_ant_t, beta)

    ! calculates reduction function value
    !wt = ETa/ET_ant

    ! sum
    Dsum = Dsum + D
    rainsum = rainsum + water_in 
    BLsum = BLsum + BL
    ETavg = ETavg + ET
    ETaavg =  ETaavg + ETa

    !Wi_ant(1) = Wi(1)
    !Wi_ant(2) = Wi(2)
		Wi_ant_t = Wit
		Wx_ant_t = Wxt
    ET_ant = ET

  endif
end subroutine


end program BMG
