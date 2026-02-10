program soilcalc

! generate an input file with the makesoil_minimal.sh
! and an empty output file with soildata_template.cdl

! use Makefile to compile

use netcdf
use parametersmod,     only : i1,i2,sp,dp
use soilpropertiesmod, only : omcf,soildata,soilproperties

implicit none

! parameter

real(sp) :: rmissing = -9999.

! local variables

character(200) :: infile
character(200) :: outfile

integer :: status
integer :: ifid
integer :: ofid
integer :: dimid
integer :: varid
integer :: xlen
integer :: ylen
integer :: nl

integer :: i
integer :: j
integer :: l

real(sp), allocatable, dimension(:,:) :: datacheck

real(sp), dimension(2) :: actual_range

! input variables

integer(i2) :: missing
real(sp)    :: scale_factor
real(sp)    :: add_offset

real(dp), allocatable, dimension(:) :: x   ! projection x coordinate
real(dp), allocatable, dimension(:) :: y   ! projection y coordinate

real(dp), allocatable, dimension(:,:) :: lon   ! geodetic longitude
real(dp), allocatable, dimension(:,:) :: lat   ! geodetic latitude

real(sp), allocatable, dimension(:) :: zpos   ! vertical position of layer midpoint (below soil surface, vertical down) (cm)
real(sp), allocatable, dimension(:) :: dz     ! soil layer thickness (cm)
real(sp), allocatable, dimension(:,:) :: layer_bnds     ! depth coordinate of the layer boundaries

integer(i1), allocatable, dimension(:,:)   :: usda  ! USDA soil classification (code)
real(sp), allocatable, dimension(:,:,:) :: sand  ! sand content by mass (fraction)
real(sp), allocatable, dimension(:,:,:) :: clay  ! clay content by mass (fraction)
real(sp), allocatable, dimension(:,:,:) :: cfvo  ! coarse fragment content by volume (fraction)
real(sp), allocatable, dimension(:,:,:) :: soc   ! soil organic carbon content by mass (fraction)

! output variables

real(sp), allocatable, dimension(:,:,:) :: silt   ! silt content by mass implied as residual sand and clay (fraction)
real(sp), allocatable, dimension(:,:,:) :: bulk   ! bulk density (kg m-3)
real(sp), allocatable, dimension(:,:,:) :: Tsat   ! soil porosity (Theta-sat) (fraction)
real(sp), allocatable, dimension(:,:,:) :: T33    ! soil water content at field capacity (-33 kPa) (fraction)
real(sp), allocatable, dimension(:,:,:) :: T1500  ! soil water content at permanent wilting point (psi -1500 kPa) (fraction)
real(sp), allocatable, dimension(:,:,:) :: whc    ! water holding capacity defined as T33-T1500, reduced for coarse fragment volume (fraction)
real(sp), allocatable, dimension(:,:,:) :: Ksat   ! saturated hydraulic conductivity (mm h-1)
real(sp), allocatable, dimension(:,:,:) :: ki     ! intrinsic permeability (m2)
real(sp), allocatable, dimension(:,:,:) :: lambda ! pore size distribution index (unitless)
real(sp), allocatable, dimension(:,:,:) :: psi_e  ! tension at air entry (mm)
real(sp), allocatable, dimension(:,:,:) :: psi_f  ! capillary head at the wetting front (mm)

! data structure for single-pixel values

type(soildata) :: soil

! -------------------------------------------------------
! open input and output files

call getarg(1,infile)
call getarg(2,outfile)

status = nf90_open(infile,nf90_nowrite,ifid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_open(outfile,nf90_write,ofid)
if (status /= nf90_noerr) call handle_err(status)

! ---
! get dimensions and coordinate variables from input file
! allocate memory for input variables

status = nf90_inq_dimid(ifid,'x',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'y',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'depth',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=nl)
if (status /= nf90_noerr) call handle_err(status)

allocate(zpos(nl))
allocate(dz(nl))
allocate(layer_bnds(2,nl))
allocate(x(xlen))
allocate(y(ylen))
allocate(lon(xlen,ylen))
allocate(lat(xlen,ylen))

! ---
! get input coordinate variables

status = nf90_inq_varid(ifid,'depth',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,zpos)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'dz',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,dz)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'layer_bnds',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,layer_bnds)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'x',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,x)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'y',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,y)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

! ---
! write out coordinate variables

status = nf90_inq_varid(ofid,'depth',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,zpos)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'dz',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,dz)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'layer_bnds',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,layer_bnds)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'x',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,x)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'y',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,y)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

! ---
! free up memory for variables no longer needed

deallocate(x)
deallocate(y)
deallocate(lon)
deallocate(lat)

! ---
! allocate soil state variables for one pixel

allocate(soil%layer(nl))

soil%layer%zpos = zpos
soil%layer%dz   = dz

do l = 1,nl
  write(0,*)l,zpos(l),dz(l)
end do

! ---
! allocate input and output arrays

!input

allocate(usda(xlen,ylen))
allocate(sand(xlen,ylen,nl))
allocate(clay(xlen,ylen,nl))
allocate(cfvo(xlen,ylen,nl))
allocate(soc(xlen,ylen,nl))

! output

allocate(silt(xlen,ylen,nl)) 
allocate(bulk(xlen,ylen,nl)) 
allocate(Tsat(xlen,ylen,nl)) 
allocate(T33(xlen,ylen,nl)) 
allocate(T1500(xlen,ylen,nl)) 
allocate(whc(xlen,ylen,nl)) 
allocate(Ksat(xlen,ylen,nl)) 
allocate(ki(xlen,ylen,nl)) 
allocate(lambda(xlen,ylen,nl)) 
allocate(psi_e(xlen,ylen,nl)) 
allocate(psi_f(xlen,ylen,nl)) 

! -------------------------------------------------------
! read input soil spatial data

status = nf90_inq_varid(ifid,'USDA',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,usda)
if (status /= nf90_noerr) call handle_err(status)

call getvar(ifid,'sand',sand)
call getvar(ifid,'clay',clay)
call getvar(ifid,'cfvo',cfvo)
call getvar(ifid,'soc',soc)

status = nf90_close(ifid)
if (status /= nf90_noerr) call handle_err(status)

! ---------

bulk   = rmissing
Tsat   = rmissing
T33    = rmissing
T1500  = rmissing
whc    = rmissing
Ksat   = rmissing
lambda = rmissing
psi_e  = rmissing

write(0,*)'calculating'

do j = 1,ylen
  do i = 1,xlen
  
    soil%usda = usda(i,j)
    
    if (soil%usda <= 3) cycle  ! skip all cells with no soil classification
      
    soil%layer%sand = sand(i,j,:)
    soil%layer%clay = clay(i,j,:)
    soil%layer%cfvo = cfvo(i,j,:)
    soil%layer%orgm =  soc(i,j,:) * omcf
    
    call soilproperties(soil)      
    
    bulk(i,j,:)   = soil%layer%bulk
    Tsat(i,j,:)   = soil%layer%Tsat
    T33(i,j,:)    = soil%layer%T33
    T1500(i,j,:)  = soil%layer%T1500
    whc(i,j,:)    = soil%layer%whc
    Ksat(i,j,:)   = soil%layer%Ksat
    lambda(i,j,:) = soil%layer%lambda
    psi_e(i,j,:)  = soil%layer%psi_e
    
    ! check for bad data
    if (any(soil%layer%whc <= 0)) then
      write(0,*)i,j
      do l = 1,nl
        write(0,*)l,soil%layer(l)%sand,soil%layer(l)%clay,soil%layer(l)%whc
      end do
    end if

  end do
end do

! ---------
! write output variables

write(0,*)'writing'

call putvar(ofid,'bulk',bulk)
call putvar(ofid,'Tsat',Tsat)
call putvar(ofid,'T33',T33)
call putvar(ofid,'T1500',T1500)
call putvar(ofid,'whc',whc)
call putvar(ofid,'Ksat',Ksat)
call putvar(ofid,'lambda',lambda)
call putvar(ofid,'psi_e',psi_e)

! ---------
! close output file

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

! -------------------------------------------------------

contains

! ---------

subroutine getvar(ncid,varname,var)

use parametersmod, only : i2,sp

implicit none

integer,                    intent(in)  :: ncid     ! netCDF file ID
character(*),               intent(in)  :: varname  ! netCDF variable name
real(sp), dimension(:,:,:), intent(out) :: var      ! variable values to read

real(sp)    :: scale_factor
real(sp)    :: add_offset
integer(i2) :: missing_value

integer(i2), allocatable, dimension(:,:,:) :: ivar

integer :: xlen,ylen,nl

! ----

xlen = size(var,dim=1)
ylen = size(var,dim=2)
nl   = size(var,dim=3)

allocate(ivar(xlen,ylen,nl))

status = nf90_inq_varid(ncid,varname,varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,ivar)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'missing_value',missing_value)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'add_offset',add_offset)
if (status == nf90_enotatt) then
  add_offset = 0.
else
  if (status /= nf90_noerr) call handle_err(status)
end if

var = rmissing

where (ivar /= missing_value) var = real(ivar) * scale_factor + add_offset

write(0,*)varname,': ',minval(var,mask=var/=rmissing),maxval(var,mask=var/=rmissing)

end subroutine getvar

! ---------

subroutine putvar(ncid,varname,var)

use parametersmod, only : i2,sp

implicit none

! arguments

integer,                    intent(in) :: ncid     ! netCDF file ID
character(*),               intent(in) :: varname  ! netCDF variable name
real(sp), dimension(:,:,:), intent(in) :: var      ! variable values to write

! parameter

real(sp) :: rmissing = -9999.

! ----

status = nf90_inq_varid(ncid,varname,varid)
if (status /= nf90_noerr) call handle_err(status)

actual_range = [minval(var,mask=var /= rmissing),maxval(var,mask=var /= rmissing)]

write(0,*)varname,' range: ',actual_range

status = nf90_put_att(ncid,varid,'actual_range',actual_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,var)
if (status /= nf90_noerr) call handle_err(status)

end subroutine putvar

! ---------

subroutine handle_err(status)

!   Internal subroutine - checks error status after each netcdf call,
!   prints out text message each time an error code is returned. 

integer, intent (in) :: status

if(status /= nf90_noerr) then 
  write(0,*)'NetCDF error: ',trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

end program soilcalc
