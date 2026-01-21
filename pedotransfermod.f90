module pedotransfermod

! module containing the pedotransfer functions of Balland et al. (2008) as adjusted and modified by Sandoval et al. (2024)
! references:
! Balland, V., Pollacco, J. A. P., & Arp, P. A. (2008). Modeling soil hydraulic properties for a wide range of soil conditions.
!   Ecological Modelling, 219(3-4), 300-316. doi:10.1016/j.ecolmodel.2008.07.009
! Sandoval, D., Prentice, I. C., & Nóbrega, R. L. B. (2024). 
!   Simple process-led algorithms for simulating habitats (SPLASH v.2.0): robust calculations of water and energy fluxes. 
!   Geoscientific Model Development, 17(10), 4229-4309. doi:10.5194/gmd-17-4229-2024

use parametersmod, only : sp

implicit none

public :: fDp
public :: fDb
public :: fTsat
public :: fT33
public :: fT1500
public :: calcKsat

! SoilGrids map legend for USDA soil orders
! Histosols: 5,10-13
! Spodosols: 15-19
! Inceptisols: 85,86,98-94
! Entisols: 95-99
! Aridisols: 50-56 
! Mollisols: 69,70,72,73,74,77
! Aquoll, argiaquoll, argiudoll, cryaquoll, duraquoll, natrustoll and paleustoll 71,75,76
! Alfisols: 80-84
! Oxisols: 30-34
! Ultisols: 60-64
! Vertisols: 40-45
! Andisols: 20-27

! Typical bulk densities of coarse fragments
!  Sedimentary rocks   Siliceous  Mafic Peridotite
! 1.8 2.0 2.2 2.4 2.6   2.7 2.8    3.0     3.2

real(sp), parameter :: Dcf = 2.7  ! bulk density of coarse fragments (g cm-3)

contains

! ----------------------------

real(sp) function fDp(orgm,cfvo)

! function to estimate soil particle density (g cm-3) (Sandoval et al., 2024, eqn A26)

use parametersmod, only : sp

implicit none

! argument

real(sp), intent(in) :: orgm  ! soil organic matter (mass fraction)
real(sp), intent(in) :: cfvo  ! coarse fragment content (volume fraction)

! local variable

real(sp) :: Dp

! ---

Dp = 1. / (orgm / 1.3 + (1. - orgm) / 2.65)

fDp = Dp * (1. - cfvo) + cfvo * Dcf

end function fDp

! ----------------------------

real(sp) function fDb(usda,clay,cfvo,zpos,orgm,Dp)

! function to estimate soil bulk density (Balland et al., 2008, eqn 18) in g cm-3
! This function requires information about the USDA soil type to determine the parameter set to use,
! e.g., the field SoilGrids250m 2017-03 - Predicted USDA 2014 suborder classes 
! NB not all of the soil suborders listed in Table 6 of Balland et al. (2008)
! are present in the SoilGrids250m dataset.

use parametersmod, only : sp

implicit none

! arguments

integer,  intent(in) :: usda  ! USDA 2014 suborder class (code, complete legend in SoilGrids TAXOUSDA file)
real(sp), intent(in) :: clay  ! clay content (mass fraction)
real(sp), intent(in) :: cfvo  ! coarse fragment content (volume fraction)
real(sp), intent(in) :: zpos  ! below surface (cm)
real(sp), intent(in) :: orgm  ! organic matter content (mass fraction)
real(sp), intent(in) :: Dp    ! particle density (g cm-3)

!  Parameters specific to USDA soil order and suborder

real(sp), dimension(4,19), parameter :: parstable = reshape( & 
  [ 1.18, 0.89, 0.022, 6.27,  &  !  1 Histosols
    1.18, 0.70, 0.022, 6.27,  &  !  2 Spodosols
    1.18, 0.89, 0.022, 6.27,  &  !  3 Inceptisols General
    1.50, 1.10, 0.010, 6.27,  &  !  4 Inceptisols Exceptions: dystropept, eutropept, halaquept, tropept, ustropept and xerumbrept
    1.50, 1.10, 0.022, 6.27,  &  !  5 Entisols    General
    1.50, 1.10, 0.002, 6.27,  &  !  6 Entisols    Exceptions: cryorthent
    1.50, 1.10, 0.022, 6.27,  &  !  7 Aridisols
    1.50, 1.10, 0.022, 6.27,  &  !  8 Mollisols   General
    1.50, 1.40, 0.022, 6.27,  &  !  9 Mollisols   Exceptions: aquoll, argiaquoll, argiudoll, cryaquoll, duraquoll, natrustoll and paleustoll
    1.65, 1.10, 0.022, 6.27,  &  ! 10 Alfisols    General
    1.18, 0.89, 0.022, 6.27,  &  ! 11 Alfisols    Exceptions: cryoboralf, fragiboralf, kandiudalf, kanhaplustalf and paleudalf
    1.40, 1.20, 0.002, 6.27,  &  ! 12 Oxisols     General
    1.40, 1.20, 0.001, 6.27,  &  ! 13 Oxisols     Exceptions: kandiudox
    1.55, 1.00, 0.005, 6.27,  &  ! 14 Ultisols    General
    1.55, 1.00, 0.001, 6.27,  &  ! 15 Ultisols    Exceptions: haploxerult, kandihumult, kanhaplustult, palehumult, rhodudult and udult
    1.90, 1.20, 0.010, 6.27,  &  ! 16 Vertisols
    1.40, 1.10, 0.050, 6.27,  &  ! 17 Andisols    General
    1.00, 0.00, 0.001, 2.,    &  ! 18 Andisols    Exceptions: haplaquands, hapludands, haplustands, melanudands, udivitrands, vitraquands, vitricryands, vitritorrands, vitrixerands
    1.00, 0.00, 0.001, 2. ],  &  ! 19 Andisols    Exceptions: melanudands, udivitrands, vitraquands, vitricryands, vitritorrands, vitrixerands
  shape(parstable))

! local variables

real(sp), dimension(4) :: pars
real(sp) :: a
real(sp) :: b
real(sp) :: c
real(sp) :: d

real(sp) :: Dbs

! ----

select case(usda)
case(5,10:13)             ! Histosols, including Gelisols: Histels
  pars = parstable(:,1)
case(15:19)               ! Spodosols
  pars = parstable(:,2)
case(85:86,89:94)         ! Inceptisols general
  pars = parstable(:,3)
case(95:99)               ! Entisols general
  pars = parstable(:,5)
case(6:7)                 ! Gelisols (Orthels and Turbels), using table definition for Cryorthent
  pars = parstable(:,6)
case(50:56)               ! Aridisols
  pars = parstable(:,7)
case(69:70,72:74,77)      ! Mollisols general
  pars = parstable(:,8)
case(71,75:76)            ! Mollisols exceptions
  pars = parstable(:,9)
case(80:84)               ! Alfisols general
  pars = parstable(:,10)
case(30:34)               ! Oxisols general
  pars = parstable(:,12)
case(60:64)               ! Ultisols general
  pars = parstable(:,14)
case(40:45)               ! Vertisols
  pars = parstable(:,16)
case(20:27)               ! Andisols general
  pars = parstable(:,17)
case default
  write(0,*)'USDA type',usda,' was not classified, using average parameter set'
  pars = [ 1.40, 0.95, 0.015, 5.82 ]
end select

a = pars(1)
b = pars(2)
c = pars(3)
d = pars(4)

Dbs = (a + (Dp - a - b * (1. - clay)) * (1. - exp(-c * zpos))) / (1. + d * orgm)

fDb = Dbs * (1. - cfvo) + cfvo * Dcf  ! Balland et al. (2008) eqn A.13

end function fDb

! ----------------------------

real(sp) function fTsat(Dp,Db)

! function to estimate soil porosity Theta-sat (Sandoval et al., 2024, eqn A25c) (units)

use parametersmod, only : sp

implicit none

real(sp), intent(in) :: Dp  ! particle density (g cm-3)
real(sp), intent(in) :: Db  ! bulk density (g cm-3)

fTsat = 1. - Db / Dp

end function fTsat

! ----------------------------

real(sp) function fT33(Tsat,clay,sand,orgm)

! function to estimate soil water content at field capacity Theta-33 (Sandoval et al., 2024, eqn A25b) (units)

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: Tsat  ! soil porosity (units)
real(sp), intent(in) :: clay  ! clay content (mass fraction)
real(sp), intent(in) :: sand  ! sand content (mass fraction)
real(sp), intent(in) :: orgm  ! soil organic matter content (mass fraction)

! parameters

real(sp), parameter :: a = -0.0547
real(sp), parameter :: b = -0.0010
real(sp), parameter :: c =  0.4760
real(sp), parameter :: d =  0.9402

! ----

fT33 = Tsat * (c + (d - c) * clay**0.5) * exp((a * sand - b * orgm) / Tsat)

end function fT33

! ----------------------------

real(sp) function fT1500(T33,clay)

! function to estimate soil water content at wilting point Theta-1500 (Sandoval et al., 2024, eqn A25a) (units)

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: T33   ! field capacity water content (volume fraction)
real(sp), intent(in) :: clay  ! clay content (mass fraction)

! parameters

real(sp), parameter :: c =  0.2018
real(sp), parameter :: d =  0.7809

! ----

fT1500 = T33 * (c + (d - c) * clay**0.5)

end function fT1500

! ----------------------------

subroutine calcKsat(sand,clay,orgm,Db,Tsat,T33,T1500,lambda,Ksat)

! function to estimate saturated hydraulic conductivity (Sandoval et al., 2024) (mm h-1)
! NB this equation comes from the code on github in the file splash.point.R, lines 351-363
! because eqn A25d in the GMD paper does not appear to produce valid results.

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in)  :: sand   ! sand content (mass fraction)
real(sp), intent(in)  :: clay   ! clay content (mass fraction)
real(sp), intent(in)  :: orgm   ! organic matter content (mass fraction)
real(sp), intent(in)  :: Db     ! bulk density (g cm-3) 
real(sp), intent(in)  :: Tsat   ! porosity
real(sp), intent(in)  :: T33    ! fractional water contentent at field capacity
real(sp), intent(in)  :: T1500  ! sand content (mass fraction)
real(sp), intent(out) :: lambda ! pore size distribution index (unitless)
real(sp), intent(out) :: Ksat   ! saturated hydraulic conductivity (mm h-1)

! parameters

real(sp), parameter :: Ksmax = 857.48454

real(sp), parameter :: k2 = -2.70927
real(sp), parameter :: k3 =  3.62264
real(sp), parameter :: k4 =  7.33398
real(sp), parameter :: k5 = -8.11795
real(sp), parameter :: k6 = 18.75552
real(sp), parameter :: k7 =  1.03319

real(sp), parameter :: l1500 = log(1500.)
real(sp), parameter :: l33   = log(33.)
real(sp), parameter :: num   = l1500 - l33

! local variables

real(sp) :: Tdrain  ! water content 
real(sp) :: B       ! shape coefficient of the water retention curve = 1 / lambda (Sandoval eqn 56c)

! ----

Tdrain = Tsat - T33

B = num / (log(T33) - log(T1500))  ! Saxton and Rawls eqn 15

lambda = 1. / B

Ksat = Ksmax / (1. + exp(k2 * sand + k3 * Db + k4 * clay + k5 * Tdrain + k6 * orgm + k7 * lambda))

end subroutine calcKsat

! ----------------------------

real(sp) function fPsi_e(sand,clay,orgm,T33,lambda)

! Equations to estimate the "tension at air entry (bubbling pressure)" of a soil (mm) from
! Saxton, K. E., & Rawls, W. J. (2006). Soil Water Characteristic Estimates by Texture and Organic Matter
! for Hydrologic Solutions. Soil Science Society of America Journal, 70(5), 1569. doi:10.2136/sssaj2005.0117
! used in Sandoval et al. (2024) eqn 29 (called psi_b by Sandoval)
! NB The Saxton and Rawls paper calls for inputs in mass percent, but it appears the equations require fraction 
! to return realistic outputs.

! As Sandoval notes in the rsplash code (lines 365-383) the Saxton and Rawls regression sometimes leads to
! unphysical (negative) values for air entry tension. This is also apparent in Saxton and Rawls Fig. 1D
! Negative values for air-entry pressure would imply hydrophobic soils and have been observed experimentally
! but for this model, Sandoval sets the minimum air entry pressure equal to the Saxton and Rawls "A" value,
! which represents the air entry pressure at the intercept of the log-log retention curve.

implicit none

! arguments

real(sp), intent(in) :: sand   ! sand content (mass fraction)
real(sp), intent(in) :: clay   ! clay content (mass fraction)
real(sp), intent(in) :: orgm   ! organic matter content (mass fraction)
real(sp), intent(in) :: T33    ! soil water content at field capacity (volume fraction)
real(sp), intent(in) :: lambda ! pore size distribution index (unitless)

! parameters

real(sp), parameter :: kPa2mm = -1000. / 9.80665 ! source https://en.wikipedia.org/wiki/Centimetre_or_millimetre_of_water
real(sp), parameter :: ln33 = log(33.)

! local variables

real(sp) :: A
real(sp) :: B
real(sp) :: TS33t
real(sp) :: TS33
real(sp) :: psi_et

real(sp) :: psi_e

! ----

! Saxton and Rawls eqn 3

TS33t = 0.278 * sand + 0.034 * clay + 0.022 * orgm - 0.018 * sand * orgm - 0.027 * clay * orgm - 0.584 * sand * clay + 0.078 

TS33 = TS33t + (0.636 * TS33t - 0.107)

! Saxton and Rawls eqn 4

psi_et = -21.67 * sand - 27.93 * clay - 81.97 * TS33 + 71.12 * sand * TS33 + 8.29 * clay * TS33 + 14.05 * sand * clay + 27.16

psi_e = psi_et + (0.02 * psi_et**2 - 0.113 * psi_et - 0.70)

if (psi_e < 0.) then

  B = 1. / lambda
  A = exp(ln33 + B * log(T33))  ! Saxton and Rawls eqn. 14

  psi_e = A

end if

fPsi_e = psi_e * kPa2mm

end function fPsi_e

! ----------------------------

end module pedotransfermod
