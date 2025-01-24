program masklandmask

! gfortran -o masklandmask masklandmask.f90 -I/usr/local/include -L/usr/local/lib -lnetcdff

use iso_fortran_env
use netcdf

implicit none

integer, parameter :: sp = real32
integer, parameter :: i2 = int16

character(200) :: classfile   ! the file containing the land cover classes
character(200) :: targetfile  ! the file to be modified

real(sp),    allocatable, dimension(:,:,:) :: classfrac
integer(i2), allocatable, dimension(:,:)   :: var

integer :: status
integer :: ncid
integer :: dimid
integer :: varid
integer :: xlen
integer :: ylen
integer :: nclass

integer :: x
integer :: y

integer(i2) :: imissing

real(sp) :: landfrac

! ---------
! G3WBM (surface cover classes)
! 0: Land
! 1: Land (No Landsat observation)
! 10: Snow
! 20: Wet Soil / Wet Vegetation / Lava
! 30: Salt Marsh
! 40: Temporal Flooded Area
! 50: Permanent Water
! 51: Permanent Water (Added by SWBD)
! 99: Ocean (Given by external land/sea mask)

! ---------

call getarg(1,classfile)

status = nf90_open(classfile,nf90_nowrite,ncid)
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

status = nf90_inq_dimid(ncid,'class',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=nclass)
if (status /= nf90_noerr) call handle_err(status)

allocate(classfrac(xlen,ylen,nclass))

status = nf90_inq_varid(ncid,'classfrac',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,classfrac)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! ---------

call getarg(2,targetfile)

allocate(var(xlen,ylen))

status = nf90_open(targetfile,nf90_write,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'Band1',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,var)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'missing_value',imissing)
if (status /= nf90_noerr) call handle_err(status)

do y = 1,ylen
  do x = 1,xlen
  
    landfrac = sum(classfrac(x,y,1:2))
    
    if (landfrac <= 0.) var(x,y) = imissing

  end do
end do

status = nf90_put_var(ncid,varid,var)
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

end program masklandmask