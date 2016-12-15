program biosim
use variables
use screen
use wsb
implicit none

! initialization of variables
linux = .true.    ! set to .false. when compiling under Windows machine

call iniScreen()
call readInput()
call printInitial()
call waterBalance()

contains

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
