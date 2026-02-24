module soilpropertiesmod

use parametersmod, only : i1,sp

implicit none

type layerinfo
  real(sp) :: zpos    ! depth of layer midpoint from soil surface (cm)
  real(sp) :: dz      ! soil layer thickness (cm)
  real(sp) :: sand    ! mass fraction
  real(sp) :: silt    ! mass fraction
  real(sp) :: clay    ! mass fraction
  real(sp) :: cfvo    ! coarse fragment content (volume fraction)
  real(sp) :: orgm    ! organic matter content (mass fraction)
  real(sp) :: bulk    ! bulk density (g m-3)
  real(sp) :: Tsat    ! porosity (fraction)
  real(sp) :: T33     ! water content at -33 KPa tension (fraction)
  real(sp) :: T1500   ! water content at -1500 KPa tension (fraction)
  real(sp) :: whc     ! water holding capacity defined as -33 - -1500 KPa tension (fraction)
  real(sp) :: lambda  ! pore size distribution index (unitless) also called 1/B
  real(sp) :: psi_e   ! tension at air entry (bubbling pressure) (kPa)
  real(sp) :: psi_f   ! capillary head at the wetting front (mm)
  real(sp) :: Ksat    ! saturated hydraulic conductivity (mm h-1)
  real(sp) :: ki      ! intrinsic permeability (m2)
end type layerinfo

type soildata
  integer(i1) :: WRB     ! WRB 2006 subgroup (code)
  integer(i1) :: USDA    ! WRB 2006 subgroup (code)
  type(layerinfo), allocatable, dimension(:) :: layer
end type soildata

real(sp), parameter :: omcf = 1.724  ! conversion factor from organic carbon to organic matter (Nelson & Sommers 1996)

contains

! ---------------------------------------------------------------------------

subroutine soilproperties(soil)

use parametersmod,   only : sp
use pedotransfermod, only : fDp,fDb,fTsat,fT33,fT1500,calcKsat,fPsi_e

implicit none

! argument

type(soildata), intent(inout) :: soil

! local variables

integer  :: usda

real(sp) :: zpos

real(sp) :: sand
real(sp) :: silt
real(sp) :: clay
real(sp) :: orgm
real(sp) :: cfvo

real(sp) :: Db
real(sp) :: Dp

real(sp) :: Tsat
real(sp) :: T33
real(sp) :: T1500

real(sp) :: lambda
real(sp) :: psi_e
real(sp) :: psi_f
real(sp) :: Ksat
real(sp) :: ki

integer :: nl
integer :: l

real(sp) :: particlesum

! ----------
! initalize output variables

soil%layer%bulk   = -9999.
soil%layer%Tsat   = -9999.
soil%layer%T33    = -9999.
soil%layer%T1500  = -9999.
soil%layer%whc    = -9999.
soil%layer%lambda = -9999.
soil%layer%psi_e  = -9999.
soil%layer%psi_f  = -9999.
soil%layer%Ksat   = -9999.
soil%layer%ki     = -9999.


nl = size(soil%layer)

usda = soil%usda

do l = 1,nl

  ! assign input variables

  zpos = soil%layer(l)%zpos
  sand = soil%layer(l)%sand
  silt = soil%layer(l)%silt
  clay = soil%layer(l)%clay
  orgm = soil%layer(l)%orgm 
  cfvo = soil%layer(l)%cfvo
    
  ! adjust particle size distribution so the mass fractions sum to 1

  particlesum = sand + silt + clay

  sand = sand / particlesum
  silt = silt / particlesum
  clay = clay / particlesum

  ! further adjust for organic matter
  
  sand = sand - sand * orgm
  silt = silt - silt * orgm
  clay = clay - clay * orgm

  if (orgm > 0.08) cycle  ! skip organic soils for the moment
  
  ! ---

  Dp = fDp(orgm,cfvo)
    
  Db = fDb(usda,clay,cfvo,zpos,orgm,Dp)

  Tsat = fTsat(Dp,Db)
    
  T33 = fT33(Tsat,clay,sand,orgm)
    
  T1500 = fT1500(T33,clay)

  psi_e = fPsi_e(sand,clay,orgm,T33,lambda)

  psi_f = (2. + 3. * lambda) / (1. + 3. * lambda) * psi_e / 2.  ! Sandoval et al. (2024) eqn 29, NB this value will be negative

  call calcKsat(sand,clay,orgm,Db,Tsat,T33,T1500,lambda,Ksat,ki)

  ! ---
  ! assign output
  
  soil%layer(l)%sand   = sand
  soil%layer(l)%silt   = silt
  soil%layer(l)%clay   = clay

  soil%layer(l)%bulk   = Db
  soil%layer(l)%Tsat   = Tsat
  soil%layer(l)%T33    = T33
  soil%layer(l)%T1500  = T1500
  soil%layer(l)%whc    = T33 - T1500
  soil%layer(l)%lambda = lambda
  soil%layer(l)%psi_e  = psi_e
  soil%layer(l)%psi_f  = psi_f
  soil%layer(l)%Ksat   = Ksat
  soil%layer(l)%ki     = ki

end do

end subroutine soilproperties

! ---------------------------------------------------------------------------

end module soilpropertiesmod