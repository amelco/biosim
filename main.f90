program biosim
use variables
use screen
use wsb
implicit none

integer :: tnd    ! total number of days (number of data lines of wheater input file)
integer :: day    ! current day
integer :: io
!integer :: sqd    ! sequential day (1 to 365)
!integer :: year   ! current year
! Global soil variables
!real :: Or1, Os1, a1, n1, Ks1    ! parameters from 1st soil layer
!real :: Or2, Os2, a2, n2, Ks2    ! parameters from 2nd soil layer
!real :: z1, z2                   ! depths of soil layer


!real :: water_in, Rad, ET

! initialization of variables
linux = .true.    ! set to .false. when compiling under Windows machine
tnd = 0
day = 1


call iniScreen()
call readInput()
call printInitial()

! opening wheater data for reading
open(21,file='data.prn',status='old')       ! wheater input values
read(21,*)
read(21,*)

! begin of loop (dt=1 day)
do while (day .lt. tnd)
  
  read(21,*) year, SQD, water_in, Rad, ET
  if (sqd == 365) then
    !print*, year
    year = year + 1
  endif
  day = day + 1
enddo
close(21)
call waterBalance()

contains

subroutine readInput()
  ! reads all input files
  open(10, file="input.txt", status="old")    ! parameters values
  read(10,*)
  read(10,*) rho, Am, tau
  read(10,*)
  read(10,*) Nm
  read(10,*)
  read(10,*)A0, N0, Bt0
  close(10)
  open(21,file='data.prn',status='old')       ! wheater input values
  ! the wheater file has 2 lines of head
  read(21,*)
  read(21,*)
  ! read initial values of year and day
  read(21,*) year, day
  tnd = 1
  do
    read(21,*,iostat=io)
    if (io .ne. 0) exit
    tnd = tnd + 1
  enddo
  close(21)
  Open(31,file='soil.in',status='old')
  read(31,*)
  read(31,*) Or1, Os1, a1, n1, Ks1
  read(31,*) Or2, Os2, a2, n2, Ks2
  read(31,*)
  read(31,*) z1, z2
  close(31)
end subroutine

subroutine printInitial()
  100 FORMAT (A30,F12.4,A10)
  101 FORMAT (A30,I12,A10)
  write(*,'(A)') "* Parameters and initial values"
  write(*,*) "---------------- forest ------------------------------------"
  write(*,100) "Tree density", rho, "tree/m2"
  write(*,100) "Max. leaf area per shoot", Am, "m2"
  write(*,100) "Time leaf unfold process", tau, "day"
  write(*,100) "Max. number of shoots", Nm, "1/tree"
  write(*,100) "Initial LAI per shoot", A0, "m2"
  write(*,100) "Initial number of shoots", N0, "1/tree"
  write(*,100) "Initial biomass of the tree", Bt0, "g/tree"
  write(*,*) "---------------- soil layer 1-------------------------------"
  write(*,100) "Residual water content", Or1, "cm3/cm3"
  write(*,100) "Saturated water content", Os1, "cm3/cm3"
  write(*,100) "alpha", a1, "1/cm"
  write(*,100) "n", n1, "-"
  write(*,101) "Number of days to simulate", tnd, ""
  write(*,*) "---------------- soil layer 2-------------------------------"
  write(*,100) "Residual water content", Or2, "cm3/cm3"
  write(*,100) "Saturated water content", Os2, "cm3/cm3"
  write(*,100) "alpha", a2, "1/cm"
  write(*,100) "n", n2, "-"
  write(*,*) "------------------------------------------------------------"
  write(*,101) "Number of days to simulate", tnd, ""
  print*
end subroutine

end program
