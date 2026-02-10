program soildepth

! gfortran -o makethickness makethickness.f90 -I/usr/local/include -L/usr/local/lib -lnetcdff
! calculate soil depth using files from the Pelletier et al. (2016) files

use iso_fortran_env
use netcdf

implicit none

integer, parameter :: i1 = int8
integer, parameter :: i2 = int16
integer, parameter :: sp = real32

character(100) :: lowlandfile ! = 'upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.nc'  (byte nd: -1)
character(100) :: uplandfile  ! = 'upland_hill-slope_soil_thickness.nc' (float nd: -1)
character(100) :: slpfracfile ! = 'hill-slope_valley-bottom.nc' (float nd: not defined, but is -1)

character(100) :: outfile

integer(i1), allocatable, dimension(:,:) :: lowland
real(sp),    allocatable, dimension(:,:) :: upland
real(sp),    allocatable, dimension(:,:) :: landformf

real(sp),    allocatable, dimension(:,:) :: thickness

integer :: status
integer :: ncid
integer :: dimid
integer :: varid
integer :: xlen
integer :: ylen

integer :: x
integer :: y

real(sp), dimension(2) :: actual_range

real(sp), parameter :: missing = -9999.

! ---------

call getarg(1,lowlandfile)
call getarg(2,uplandfile)
call getarg(3,slpfracfile)

status = nf90_open(lowlandfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lon',dimid)
if (status == nf90_ebaddim) status = nf90_inq_dimid(ncid,'x',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lat',dimid)
if (status == nf90_ebaddim) status = nf90_inq_dimid(ncid,'y',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

! ---

allocate(lowland(xlen,ylen))

status = nf90_inq_varid(ncid,'Band1',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lowland)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! ---

allocate(upland(xlen,ylen))

status = nf90_open(uplandfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'Band1',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,upland)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! ---

allocate(landformf(xlen,ylen))

status = nf90_open(slpfracfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'Band1',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,landformf)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! ---

allocate(thickness(xlen,ylen))

thickness = missing

do y = 1,ylen
  do x = 1,xlen

  if (lowland(x,y) < 0 .and. upland(x,y) < 0.) cycle
  
    if (upland(x,y) < 0.) then
    
      thickness(x,y) = real(lowland(x,y))
    
    else if (lowland(x,y) < 0) then
    
      thickness(x,y) = upland(x,y)
    
    else
    
      thickness(x,y) = upland(x,y) * landformf(x,y) + real(lowland(x,y)) * (1. - landformf(x,y))
      
    end if

  end do
end do

actual_range(1) = minval(thickness,mask=thickness /= missing)
actual_range(2) = maxval(thickness,mask=thickness /= missing)

write(0,*)'soil depth actual range',actual_range

! ---
  
call getarg(4,outfile)

status = nf90_open(outfile,nf90_write,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'Band1',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,thickness)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'_FillValue',missing)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! ---------

contains

subroutine handle_err(status)

!   Internal subroutine - checks error status after each netcdf call,
!   prints out text message each time an error code is returned. 

integer, intent (in) :: status

if(status /= nf90_noerr) then 
  write(0,*)'NetCDF error: ',trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

end program soildepth
