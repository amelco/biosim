module variables
implicit none
! 

Real water_in,Rad	     		! water input for every soil layer

real :: ET		     		! Daily evapotranspiration (mm/day)
real :: ET_ant		     		! Previous daily evapotranspiration (mm/day)
Real :: ETa                  		! Total Actual evapotranspiration (mm/day)
Real :: ETa1                  	! Acutal evapotranspiration from upper layer (mm/day)
Real :: ETa2                  	! Acutal evapotranspiration from lower layer (mm/day)
!Real,dimension(2) :: Wi      		! Actual soil-water storage (mm/day)
!Real,dimension(2) :: Wi_ant  		! Actual soil-water storage (mm/day)
Real :: dW      		! Change in soil-water storage (mm/day)
Real,dimension(2) :: Wx      		! Maximum available water for each depth, function of root depth (mm/day)
real :: Wit
real :: Wi_ant_t
Real :: g       		! Net water added to the soil (mm/day)
real :: Wxt											! Total available water
real :: Wx_ant_t											! Total available water
real :: Wz											! available water as a function of root depth (mm)
real :: beta 										! parameter for the evapotranspiration reduction function
Real :: D                    		! Total Drainage (mm/day)
Real :: D1                    	! Drainage (mm/day)
Real :: D2                    	! Drainage (mm/day)
real :: excess               		! Amount of water from layer 1 to 2 in case of drainage of 1 is greater than 0 (mm/day)
real :: lost_water							! Adds water that was lost due to condition of Wx as function of root depth (check with drainage or something...)
Real :: Wii                  		! Initial soil-water storage (mm)
Real :: Wt                   		! water stress reduction function
Real :: BL                   		! Water balance (mm/day)
real :: BLmax
integer :: eof				! end of file
integer :: Year				! 4 digits year
integer :: SQD    			! sequential day (1 to 365)

real :: Dsum				! Yearly cumulative Drainage (mm)
real :: rainsum				! Yearly cumulative rainfall (mm)
real :: BLsum				! Yearly cumulative water balance (mm)
real :: ETavg				! Yearly average potential ET (mm)
real :: ETaavg				! Yearly average actual ET (mm)

! Global soil variables
real :: Or1, Os1    			! parameters from 1st soil layer (-)
real :: Or2, Os2    			! parameters from 2nd soil layer (-)
real :: alpha1, n1, Ks1			! parameters necessary in case of using VG equation
real :: alpha2, n2, Ks2			! parameters necessary in case of using VG equation 
real :: z1, z2                   	! depths of soil layer (mm)

! Global tree parameters
real :: rho   				! tree density (trees/m2)
real :: Am    				! maximum leaf area per shoot (m2)
real :: tau   				! time of unfolding leafs (day)
real :: kt    				! light extinction coefficient (-)
real :: ft    				! intercepted radiation fraction (-)
real :: A 						! leaf area per shoot - Werf et al. (2007)  (m2)
real :: A2						! leaf area per shoot - Peter's spreadsheet (m2)
real :: Aant     			! leaf area from previous day
real :: dA				! leaf area daily rate (m2/day)
real :: root_frac			! fraction of biomass to root growth
real :: root_depth		! root depth (mm)
real :: root_max			! maximum root depth (mm)
real :: inter					! interception factor of rain by the tree
real :: kc						! crop coefficient
real :: kc_min, kc_full  ! params for kc
real ft_r 
real :: LAI						! leaf area index
real :: gama   				! bioloss parameter
real :: alpha   				! bioloss parameter

real :: N				! number of shoot (-)
real :: Nant				! number of shoot from the previous day (-)
real :: Nmax            		! maximum number of shoots (-)
real :: B				! Total Biomass (g)
real :: Bant				! Total Biomass from the previous day (g)
real :: Bshoot				! Biomass of shoots (g)
real :: Broot     			! biomass of roots (g)
real :: Bshoot_ant				! Biomass of shoots (g)
real :: Broot_ant     			! biomass of roots (g)
real :: dB				! Daily rate of total biomass (g/day)
real :: dBs				! Daily rate of total biomass of shoots (g/day)
real :: dBr				! Daily rate of total biomass of roots (g/day)
real :: rdg    ! parameter to convert produced root biomass (g) into root growth (mm) [g/mm]
real :: part_root, part_shoot    ! biomass partitioning: root and shoot

real :: Lt    				! Leaf area index per tree (-)

real :: radExt    			! radiation needed to generate 1 gram of biomass (-)
real :: bioloss_fac	  		! factor to compute biomass loss due to plant maintenance (-)

integer :: i					! counter

end module
