module variables
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


end module
