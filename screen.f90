module screen
use variables, only: linux
implicit none

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


end module
