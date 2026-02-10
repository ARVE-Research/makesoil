program maskandpack

! gfortran -o mask-and-pack mask-and-pack.f90 -I/home/public/easybuild/software/netCDF-Fortran/4.6.1-gompi-2023a/include -lnetcdff

! just grab a preformatted variable and insert into the target output file

use iso_fortran_env
use netcdf

implicit none

integer, parameter :: sp = real32
integer, parameter :: dp = real64
integer, parameter :: i2 = int16

character(100) :: infile
character(100) :: outfile
character(100) :: landfile = '/home/terraces/datasets/topography/WNA1km/WNAtopo_v2026.nc'

character(50) :: ovarname

integer :: status
integer :: ncid
integer :: dimid
integer :: varid
integer :: xlen
integer :: ylen

real(sp) :: rmissing

integer(i2), parameter :: imissing = -32768

real(sp) :: varmax

real(sp) :: scale_factor
real(sp) :: add_offset
 
real(sp),    allocatable, dimension(:,:) :: var_in
integer(i2), allocatable, dimension(:,:) :: landf
integer(i2), allocatable, dimension(:,:) :: var_out

real(sp), dimension(2) :: actual_range

! ------------------------------------------------------------------------------------------------------------

call getarg(1,infile)
call getarg(2,outfile)
call getarg(3,ovarname)

! -------------------------------------

status = nf90_open(infile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'x',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'y',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

allocate(var_in(xlen,ylen))

status = nf90_inq_varid(ncid,'Band1',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,var_in)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'_FillValue',rmissing)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! -------------

allocate(landf(xlen,ylen))
allocate(var_out(xlen,ylen))

! -------------

status = nf90_open(landfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'landf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,landf)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! --

status = nf90_open(outfile,nf90_write,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,ovarname,varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

where (landf == 0_i2) var_in = rmissing

varmax = maxval(var_in)

if ((varmax - add_offset) / scale_factor > 32767.) then
  add_offset = 32767. * scale_factor
else
  add_offset = 0.
end if

write(0,*)'sf, ao: ',scale_factor,add_offset

where (var_in /= rmissing)
  var_out = nint((var_in - add_offset) / scale_factor,i2)
elsewhere
  var_out = imissing
end where

status = nf90_put_var(ncid,varid,var_out)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

actual_range(1) = real(minval(var_out,mask=var_out /= imissing)) * scale_factor + add_offset
actual_range(2) = real(maxval(var_out,mask=var_out /= imissing)) * scale_factor + add_offset

write(0,*)'in:  ',minval(var_in,mask=var_in/=rmissing),maxval(var_in)

write(0,*)'out: ',actual_range

status = nf90_put_att(ncid,varid,'actual_range',actual_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

!-------------------------------------------------------

contains

subroutine handle_err(status)

!   Internal subroutine - checks error status after each netcdf call,
!   prints out text message each time an error code is returned. 

integer, intent (in) :: status

if(status /= nf90_noerr) then 
  print *, trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

!-------------------------------------------------------

end program maskandpack
