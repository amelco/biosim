module wsb

implicit none

Real water_in,Rad

real :: ET
real :: ET_ant
Real,dimension(2) :: Wi                   !Actual soil-water storage (mm)
Real,dimension(2) :: Wi_ant               !Actual soil-water storage (mm)
Real,dimension(2) :: dW                   !Change in soil-water storage (mm)
Real,dimension(2) :: Wx      !Maximum soil-water storage (mm) for each depth
Real,dimension(2) :: g                    !Net water added to the soil (mm)
Real :: D                    !Drainage (mm)
real :: excess               ! amount of water from 1 to 2 in case of D from 1 is greater than 0
Real :: Wii                         !Initial soil-water storage (mm)
Real :: ETa                  !Acutal evapotranspiration_ET (mm)
Real :: BL                          !Balanço hídrico (mm)
Integer :: i,eof,Year,SQD    !

! Global soil variables
real Or1, Os1, a1, n1, Ks1    ! parameters from 1st soil layer
real Or2, Os2, a2, n2, Ks2    ! parameters from 2nd soil layer
real z1, z2                   ! depths of soil layer

contains


subroutine init()
  100 FORMAT (A5,10(A8,2x)) 
  110 FORMAT (I5,10(F8.2,2x)) 
  ! read soil parameters data from input file
  !call readInputs()
  call readSoilData()
  call calcWx()
  Open(21,file='data.prn',status='old')
  Read(21,*)
  Read(21,*)
  Open(1,file='swb.out', status='replace')
  write(*,100) "SQD", "ET", "Eta", "Tr", "Wi(1)", "Wi(2)", "Wx(1)", "Wx(2)", "D", "Rain", "BL"
  write(1,100) "SQD", "ET", "Eta", "Tr", "Wi(1)", "Wi(2)", "Wx(1)", "Wx(2)", "D", "Rain", "BL"
  
  
  ! inicializacao das variaveis
  Wii = 0.5
  Wi_ant(1) = Wii
  Wi_ant(2) = Wii
  Eta = 0.0
  ET_ant = 5.65
  eof = 0
  
  do while (eof .ne. -1)
    read(21,*, iostat=eof) Year, SQD, water_in, Rad, ET
    if (eof .ne. -1) then
      call tempSWB(1, water_in)
      BL = D + dW(1) + ETa - water_in
      if (D .gt. 0.00) then
        excess = D
        call tempSWB(2, excess)
        BL = D + dW(2) - excess
      endif
      write(*,110) SQD, ET, Eta, ETa/ET_ant, Wi(1), Wi(2), Wx(1), Wx(2), D, water_in, BL
      write(1,110) SQD, ET, Eta, ETa/ET_ant, Wi(1), Wi(2), Wx(1), Wx(2), D, water_in, BL
      Wi_ant(1) = Wi(1)
      Wi_ant(2) = Wi(2)
      ET_ant = ET
    endif
  Enddo
end subroutine

!subroutine readInputs()
!
!end subroutine

! read soil parameters data from input file
subroutine readSoilData()
  Open(31,file='soil.in',status='old')
  read(31,*)
  read(31,*) Or1, Os1, a1, n1, Ks1
  read(31,*) Or2, Os2, a2, n2, Ks2
  read(31,*)
  read(31,*) z1, z2
  close(31)
end subroutine

! calculates maximum soil water storage
subroutine calcWx()
  Wx(1) = (Os1 - Or1) * z1
  Wx(2) = (Os2 - Or2) * z2
end subroutine

subroutine tempSWB(layer, water_in)
  integer :: layer
  real :: water_in
  
  if (layer .eq. 1) then
    ETa = ((Wi_ant(layer) / Wx(layer))**0.8) * ET_ant
    g(layer)=water_in-ETa
  else
    g(layer)=water_in
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

end module
