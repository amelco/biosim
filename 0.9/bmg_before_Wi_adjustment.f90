program BMG
use variables
implicit none

!-- FORMATS
100 FORMAT (2A7,11(A8,2x)) 
101 FORMAT (A5,12(A8,2x)) 
200 FORMAT (2A5,A8,  2x,A12,  2x,A6,  2x,2A12,  2x,A8)
201 FORMAT (2I5,F8.5,2x,F12.1,2x,F6.1,2x,2F12.1,2x,F8.2)

! Variable initialization
call initializeVariables()

! read soil parameters data from input file
!call readInputs()
call readSoilData()

! opens the output file
Open(1,file='swb.out', status='replace')
!write(*,101) "Year", "ET", "Eta", "Tr", "D", "Rain", "BL"
write(*,101) "Year", "SQD", "ET", "Eta", "Tr", "Wi(1)", "Wi(2)", "Wx(1)", "Wx(2)", "D", "Rain", "BL", "omega"
write(1,100) "Year", "SQD", "ET", "Eta", "Tr", "Wi(1)", "Wi(2)", "Wx(1)", "Wx(2)", "D", "Rain", "BL", "omega"
Open(2,file='bio.out', status='replace')
write(2,200) "Year", "SQD", "A", "B", "N", "Broot", "Bshoot", "rd"


! Calculates maximum soil-water storage for each depth
call calcWx(Or1, Os1, z1, Or2, Os2, z2, Wx)

Open(21,file='data.prn',status='old')
Read(21,*)
Read(21,*)

do while (eof .ne. -1)

  ! calc. root growth 
	if (i .eq. 1) then
	  root_depth = rdg * Broot
	  !root_depth = rootGrowth(reduction_fac(Wi_ant(1), Wx(1), beta), Broot)
	else
	  root_depth = root_depth + rootGrowth(reduction_fac(Wi_ant, Wz, beta, root_depth), dBr)
  endif
	
	call calcWz(Or1, Os1, z1, Or2, Os2, root_depth, Wz)

  ! calculates Leaf Area / shoot
  !A  =  LAperShoot(rho, Am, Aant, tau)
	! Temporary using Peter's equation for A
  A = LAperShoot2(rho, Am, Aant, tau, Wi_ant(1), Wz, beta)

  ! calc. num. of shoots
  N = numShoots(dBs, Nant, Bshoot, Nmax)

	!print*, root_depth, Wz, i
	!read*

  ! calc area of tree
  Lt = A*N

  ! calculates interception of rad.
  ft = 1.0 - exp(-kt*Lt)

  ! reads weather input file
  read(21,*, iostat=eof) Year, SQD, water_in, Rad, ET

  call waterBalance()

  ! calculates total, shoot and root biomass rates
  dB = biomass(wt, radExt, ft, Rad, rho, bioloss_fac, Bant)
	dBs = dB * 0.4
	dBr = dB * root_frac
  !dBs = (1.0 + wt * root_frac / (1.0-root_frac)) * (1.0 - root_frac) * dB
  !dBr = 0.6 * dB
  ! calculates total, shoot and root biomass
  B = db + Bant
  Bshoot = B * 0.4
  Broot = B * root_frac

	write(2,201) Year, SQD, A, B, N, Broot, Bshoot, root_depth
  
  Aant = A
  Bant = B
  Nant = N
	i = i+1
Enddo

contains

!subroutine readInputs()
!
!end subroutine

! read soil parameters data from input file
subroutine readSoilData()
  implicit none
  Open(31,file='soil.in',status='old')
  read(31,*)
  read(31,*) Or1, Os1, alpha1, n1, Ks1
  read(31,*) Or2, Os2, alpha2, n2, Ks2
  read(31,*)
  read(31,*) z1, z2
  close(31)
end subroutine

! calculates maximum soil water storage
subroutine calcWx(Or1, Os1, z1, Or2, Os2, z2, Wx)
  implicit none
  real, intent(in)  		  :: Or1, Os1, z1, Or2, Os2, z2
  real, intent(out), dimension(2) :: Wx
  Wx(1) = (Os1 - Or1) * z1
  Wx(2) = (Os2 - Or2) * z2
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


subroutine tempSWB(layer, water_in)
  implicit none
  integer :: layer
  real :: water_in
  
  if (layer .eq. 1) then
		ETa1 = (Wi_ant(1)/Wz)**beta 
	  ETa1 = min(ETa1, 1.0) * ET_ant
    g(layer)=water_in-ETa1
  else
	  !if (root_depth > z1) then
		  ETa2 = (Wi_ant(2)/((root_depth-z1)*(Os2-Or2)))**beta
			ETa2 = min(ETa2, 1.0) * ET_ant
		!else
		!  ETa2 = 0.0
		!endif
    g(layer)=water_in-ETa2
  endif
  
  if (g(layer) .ge. (Wx(layer) - Wi_ant(layer))) then
    dW(layer) = Wx(layer) - Wi_ant(layer)
  else if (g(layer) .lt. -Wi_ant(layer)) then
    dW(layer) = -Wi_ant(layer) 
  else
    dW(layer) = g(layer)
  endif
  
  Wi(layer)= dW(layer) + Wi_ant(layer) 
  
  if (g(layer) .ge. (Wx(layer)-Wi_ant(layer))) then
    D = g(layer) + Wi_ant(layer) - Wx(layer)
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

real function reduction_fac(Wi, Wz, beta, root_depth)
  implicit none
	real, intent(in), dimension(2) :: Wi
	real, intent(in) :: Wz, beta, root_depth

  if (root_depth > z1) then
	  reduction_fac = ( (Wi(1) + Wi(2) ) / Wz )**beta
	else
	  reduction_fac = (Wi(1) / Wz )**beta
	endif
	reduction_fac = min(reduction_fac, 1.0)
	return
end function

! Equation from Peter's spreadsheet
real function LAperShoot2(rho, Am, Ant, tau, Wi, Wz, beta)
  implicit none
	real, intent(in), dimension(2) :: Wi
  real, intent(in) :: rho, Am, Ant, tau, Wz, beta
  real :: dA
	real, parameter :: afac = 1.0

  dA = rho * (Am - Ant) / tau - afac * (1.0 - reduction_fac(Wi, Wz, beta, root_depth)) * Ant
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

real function biomass(wt, radExt, ft, rad, rho, bioloss_fac, B)
  real, intent(in) :: wt, radExt, ft, rad, rho, bioloss_fac, B

  biomass = (wt * (radExt*ft*rad / rho)) - (bioloss_fac * B)
  !biomass = (wt * (radExt*ft*rad / rho)) 
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
  Wii       	= 0.5
  Wi_ant(1) 	= Wii
  Wi_ant(2) 	= Wii
  Eta       	= 0.0
  ET_ant    	= 5.65
  beta 				= 0.8
  eof       	= 0
  Dsum      	= 0
  rainsum   	= 0
  BLsum     	= 0
  ETavg     	= 0
  ETaavg    	= 0
  Am 					= 0.05
  Aant 				= 0.0
  A 					= 0.0   !1/day
  rho 				= 0.01
  tau 				= 10.0  ! day
  Nant 				= 0.6
  Nmax 				= 10000.0 
  dBs 				= 0.0
  Bant 				= 300.0
  Bshoot 			= Bant*0.4
	Broot				= Bant*0.6
  kt 					= 0.8
  radExt 			= 1.5
  bioloss_fac = 0.0003
	root_frac   = 0.6
	root_max    = 800.0  ! mm
  Wz					= 0.0
	rdg 				= 0.2		 ! mm/g
end subroutine

subroutine waterBalance()
implicit none

110 FORMAT (2I7,11(F8.2,2x)) 

  ! calculates soil water balance
  if (eof .ne. -1) then
    call tempSWB(1, water_in)
    BL = D + dW(1) + ETa1 - water_in
    if (D .gt. 0.00) then
      excess = D
      call tempSWB(2, excess)
      BL = D + dW(2) + ETa2 - excess
    endif
	  ETa = ETa1 + ETa2
	!if(D>0.0.and.root_depth>z1) then
	!	print*, ETa1, ETa2
	!read*
	!endif

    ! show on screen at yearly timestep
    if (int(Year/4.0)*4 == Year) then
      if (SQD == 366) then
        write(*,110) Year, SQD, ET, Eta, ETa/ET_ant, Wi(1), Wi(2), Wx(1), Wx(2), D, water_in, BL, &
				& reduction_fac(Wi, Wz, beta, root_depth)
        Dsum      = 0
        rainsum   = 0
        BLsum     = 0
        ETavg     = 0
        ETaavg    = 0
     endif
    else
      if (SQD == 365) then
        write(*,110) Year, SQD, ET, Eta, ETa/ET_ant, Wi(1), Wi(2), Wx(1), Wx(2), D, water_in, BL, &
				& reduction_fac(Wi, Wz, beta, root_depth)
        Dsum      = 0
        rainsum   = 0
        BLsum     = 0
        ETavg     = 0
        ETaavg    = 0
      endif
    endif
    ! write in output file in daily timestep
    write(1,110) Year, SQD, ET, Eta, ETa/ET_ant, Wi(1), Wi(2), Wx(1), Wx(2), D, water_in, BL, &
		& reduction_fac(Wi, Wz, beta, root_depth)

    ! calculates reduction function value
    wt = ETa/ET_ant

    ! sum
    Dsum = Dsum + D
    rainsum = rainsum + water_in 
    BLsum = BLsum + BL
    ETavg = ETavg + ET
    ETaavg =  ETaavg + ETa

    Wi_ant(1) = Wi(1)
    Wi_ant(2) = Wi(2)
    ET_ant = ET

  endif
end subroutine


end program BMG

