program biosim
implicit none

!--- program global variables
logical :: linux

!--- model globbal variables
! plant variables
real(8) :: rho    ! tree density
real(8) :: Am     ! maximum leaf area
real(8) :: tau    ! time constant of leaf unfolding precess
real(8) :: A0     ! initial leaf area index
real(8) :: Nm     ! maximum nuumber of shoots
real(8) :: Bt     ! biomass of the tree
real(8) :: N0     ! initial number of shoots
real(8) :: Bt0    ! initial bbiomass of the tree


! initialization of variables
linux = .true.    ! set to .false. when compiling under Windows machine

call iniScreen()
call readInput()
call printInitial()

contains

subroutine iniScreen()
  if (linux) then
    call system("clear")
  else
    call system("cls")
  endif
  print*, "-== BIOmass SIMulator (Caattinga forest scenario) ==-"
  print*, "version 0.1"
  print*, "by A.H.F. Bezerra and E.A.R. Pinheiro, (2016)"
  print*
end subroutine

subroutine readInput()
  open(10, file="input.txt", status="old")
  read(10,*)
  read(10,*) rho, Am, tau
  read(10,*)
  read(10,*) Nm
  read(10,*)
  read(10,*)A0, N0, Bt0
  close(10)
end subroutine

subroutine printInitial()
  100 FORMAT (A30,F12.4,A10)
  write(*,'(A)') "* Parameters and initial values"
  write(*,100) "Tree density", rho, "tree/m2"
  write(*,100) "Max. leaf area per shoot", Am, "m2"
  write(*,100) "Time leaf unfold process", tau, "day"
  write(*,100) "Max. number of shoots", Nm, "1/tree"
  write(*,100) "Initial LAI per shoot", A0, "m2"
  write(*,100) "Initial number of shoots", N0, "1/tree"
  write(*,100) "Initial biomass of the tree", Bt0, "g/tree"
  print*
end subroutine

end program
