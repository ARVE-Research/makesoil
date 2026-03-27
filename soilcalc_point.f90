program soilcalc_point

! for use with ascii text input data files with four columns: sand, organic matter, observed Ksat, observed whc
! writes out the input table and calculated Ksat and whc. Writes to standard output

! compile with Makefile

use netcdf
use parametersmod,     only : i1,i2,sp,dp
use soilpropertiesmod, only : omcf,soildata,soilproperties

character(200) :: infile

integer :: io

! data structure for single-pixel values

type(soildata) :: soil

real(sp) :: sand
real(sp) :: orgm
real(sp) :: Ksat
real(sp) :: whc
real(sp) :: Ksat_calc
real(sp) :: whc_calc

! -------------------------------------------------------
! open input file

call getarg(1,infile)

open(100,file=infile,status='old')

allocate(soil%layer(1))

do

  read(100,*,iostat=io)orgm,sand,Ksat,whc
  
  if (io > 0) exit
  
  if (sand <= 0.) cycle

  soil%layer%sand = sand / 100.  ! convert from percent to fraction
  soil%layer%orgm = orgm / 100.

  call soilproperties(soil)
  
  Ksat_calc = soil%layer(1)%Ksat
  whc_calc  = soil%layer(1)%whc
  
  write(*,*)sand,orgm,Ksat,Ksat_calc,whc,whc_calc

end do

close(100)

end program soilcalc_point
