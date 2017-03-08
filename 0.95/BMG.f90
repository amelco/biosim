program BMG
use variables
implicit none

!-- FORMATS
100 FORMAT (2A7,11(A8,2x)) 
101 FORMAT (A5,10(A8,2x)) 
102 FORMAT (A7,2x,5(A8,2x)) 
200 FORMAT (2A5,A8,  2x,A12,  2x,A6,  2x,2A12,  2x,A8,  2x,A20)
201 FORMAT (2I5,F8.5,2x,F12.1,2x,F6.1,2x,2F12.1,2x,F8.2,2x,F20.4)

real :: dbsum
real :: Bx, Bxsum
real, dimension(6) :: B_last  ! stores the last 6 years of acumulated biomass
real :: avg										! average biomass growth rate from the last 5 years
logical :: rate_0							! flag to represent zero growth rate. True when avg < 500 g

rate_0 = .false.

! read soil parameters data from input file
call readInputs()
call readSoilData()
! Variable initialization
call initializeVariables()

dbsum = 0.0
Bx = 0.0
Bxsum = 0.0

! opens the output file
Open(1,file='swb.out', status='replace')
!write(*,101) "Year", "ET", "Eta", "Tr", "D", "Rain", "BL"
write(*,102) "Year", "ET", "ETa", "Drainage", "Rain", "WB"
write(1,101) "Year", "SQD", "ET", "ETa", "Tr", "Wi", "Wx", "D", "Rain", "BL", "omega"
Open(2,file='bio.out', status='replace')
write(2,200) "Year", "SQD", "A", "B", "N", "Broot", "Bshoot", "rd", "Bx"


! Calculates maximum soil-water storage for each depth

!Open(21,file='data.prn',status='old')
Open(21,file='data2.prn',status='old')
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
	!kc = ft/ft_r
	kc = kc_min + (kc_full - kc_min) * ft
	
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
	!dBs = (1.0 + reduction_fac(Wi_ant_t, Wx_ant_t, beta) * part_root / (1.0 - part_root)) * (1.0 - part_root) * dB
	dBs = dB * part_shoot
	dBr = dB * part_root
	!dBr = dB - dBs
  ! calculates total, shoot and root biomass
  B = db + Bant
  Bshoot = Bshoot_ant + dBs
  Broot = Broot_ant + dBr
	! store total biomass in array with the last 5 yearly values of biomass value
	if (SQD == 365) then
	  avg = storeBiomass(B)
		if (avg < 500.0 .and. i > 365*5) then   ! i equals cumulative days. It must be superior to 5 years
		  rate_0 = .true.
		endif
  endif

	dbsum = dbsum + dB

	Bx = (log((reduction_fac(Wi_ant_t, Wx_ant_t, beta) * radExt*ft*rad ) / rho ) + 6.9) / gama
	Bxsum = Bxsum + Bx

	write(2,201) Year, SQD, A, B, N, Broot, Bshoot, root_depth, Bx
  
  Aant = A
  Bant = B
	Bshoot_ant = Bshoot
	Broot_ant = Broot
  Nant = N
	i = i+1
Enddo

print*
print*, "Maximum water balance error: ", BLmax
print*

contains

real function storeBiomass(BioSum)
implicit none
real, intent(in) :: BioSum
real, dimension(SIZE(B_last)-1) :: diff
integer :: i

  i = 1
  do while (i < SIZE(B_last))
	  B_last(i) = B_last(i+1)
	  i = i+1
	enddo
	B_last(SIZE(B_last)) = BioSum

  i = 1
	do while (i <= SIZE(diff))
	  diff(i) = B_last(i+1) - B_last(i)
		i = i+1
	enddo
  storeBiomass = SUM(diff)/SIZE(diff)
	return 
end function

! read soil parameters data from input file
subroutine readSoilData()
  implicit none
  Open(31,file='soil.in',status='old')
  !Open(31,file='soilPeter.in',status='old')
  read(31,*)
  read(31,*) Or1, Os1
  read(31,*) Or2, Os2
  read(31,*)
  read(31,*) z1, z2
  read(31,*)
  read(31,*) Wii
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
	real::bb
	bb = 0.9
  
	! our approach of bioloss
	if (rate_0 .eqv. .false.) then
	  biomass = (reduction_fac(Wi_ant_t, Wx_ant_t, beta) * radExt*ft*rad ) / rho * ( 1.0/(exp(gama*B)))
	else
	  biomass = 0
	endif
	! Klaas
	!biomass = (reduction_fac(Wi_ant_t, Wx_ant_t, beta) * radExt*ft*rad ) / rho *  (1.0 - B/100000.0)**2
	
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
	Wi_ant_t		= wii
	Wx_ant_t		= 0.00001
  Eta       	= 0.0
  ET_ant    	= 5.65			! from weather input file (1st line)
  eof       	= 0
  Dsum      	= 0
  rainsum   	= 0
  BLsum     	= 0
  BLmax				= 0.0
  ETavg     	= 0
  ETaavg    	= 0
  A 					= 0.0   !1/day
  dBs 				= 0.0				
  Bshoot 			= Bant*part_shoot
	Broot				= Bant*part_root
	Bshoot_ant	= 200.0
	Broot_ant		= 300.0
  radExt 			= 1.5
  bioloss_fac = 0.005
	!root_frac   = 0.6
	root_max    = 800.0  ! mm
  Wz					= 0.0
  inter				= 0.13


!  ! Soil
!	Wii       	= 0.0				! input
!	! red. function
!  beta 				= .8				! input
!	ft_r				= 0.846184277466571		! input (from Peter)
!  ! plant
!  Am 					= 0.05			! input
!  Aant 				= 0.0				! input
!  Nant 				= 0.6				! input
!  Nmax 				= 10000.0 	! input
!  Bant 				= 500.0			! input
!  tau 				= 10.0  ! day (input)
!	rdg 				= 0.2		 ! mm/g (input)
!  rho 				= 0.0156		! input 
!  ! radiation
!	kt 					= 0.8				! input

end subroutine

subroutine readInputs()
implicit none
 open(7,file='params.in', status='old')
 read(7,*)
 read(7,*) beta, ft_r
 read(7,*)
 read(7,*)
 read(7,*) kt
 read(7,*)
 read(7,*)
 read(7,*) Am
 read(7,*) Aant
 read(7,*) Nant
 read(7,*) Nmax
 read(7,*) Bant
 read(7,*) tau 
 read(7,*) rdg 
 read(7,*) rho 
 read(7,*) part_root, part_shoot
 read(7,*) gama
 read(7,*)
 read(7,*)
 read(7,*) kc_min, kc_full
 close(7)
 !print*, beta, ft_r, kt, Am, Aant, Nant, Nmax, Bant, tau, rdg, rho
end subroutine

subroutine waterBalance()
implicit none

110 FORMAT (2I7,9(F8.2,2x)) 
120 FORMAT (I7,2x,5(F8.2,2x), F6.2, 2(2x,F12.2)) 

  ! calculates soil water balance
  if (eof .ne. -1) then
    call tempSWB(water_in)
    BL = D + dW + ETa - water_in * (1.0 - Lt*inter)

    ! show on screen at yearly timestep
    if (int(Year/4.0)*4 == Year) then
      if (SQD == 366) then
        !write(*,120) Year, ETavg, ETaavg, Dsum, rainsum, BLsum, dbsum/SQD, B, Bxsum/SQD
        write(*,120) Year, ETavg, ETaavg, Dsum, rainsum, BLsum, dbsum/SQD, B, avg
        Dsum      = 0
        rainsum   = 0
        BLsum     = 0
        ETavg     = 0
        ETaavg    = 0
				dbsum = 0.0
				Bxsum = 0.0
     endif
    else
      if (SQD == 365) then
        !write(*,120) Year,ETavg, ETaavg, Dsum, rainsum, BLsum, dbsum/SQD, B, Bxsum/SQD
        write(*,120) Year, ETavg, ETaavg, Dsum, rainsum, BLsum, dbsum/SQD, B, avg
        Dsum      = 0
        rainsum   = 0
        BLsum     = 0
        ETavg     = 0
        ETaavg    = 0
				dbsum = 0.0
				Bxsum = 0.0
      endif
    endif
    ! write in output file in daily timestep
    write(1,110) Year, SQD, ET, Eta, ETa/ET_ant, Wit, Wxt, D, water_in, BL, &
		& reduction_fac(Wi_ant_t, Wx_ant_t, beta)

		if (BL > BLmax) then
		  BLmax = BL
		endif

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
