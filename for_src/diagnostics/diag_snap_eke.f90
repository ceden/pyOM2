

  

subroutine init_snap_cdf_eke
!=======================================================================
!     initialize NetCDF snapshot file for EKE variables
!=======================================================================
 use main_module   
 use eke_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_eke.cdf'

 call def_grid_cdf('pyOM_eke.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_eke.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)


       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'EKE', NCFLOAT,4,dims,iret)
       name = 'meso-scale energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id  = ncvdef (ncid,'L_Rossby', NCFLOAT,3,dims,iret)
       name = 'Rossby Radius'; unit = 'm'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'L_Rhines', NCFLOAT,4,dims,iret)
       name = 'Rhines scale'; unit = 'm'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'eke_diss_iw', NCFLOAT,4,dims,iret)
       name = 'Dissipation of EKE to IW'; unit = 'm^2/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'eke_diss_tke', NCFLOAT,4,dims,iret)
       name = 'Dissipation of EKE to TKE'; unit = 'm^2/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id  = ncvdef (ncid,'eke_bot_flux', NCFLOAT,3,dims,iret)
       name = 'flux by bottom friction'; unit = 'm^3/s^3'
       call dvcdf(ncid,id,name,len_trim(name),unit,16,spval)

       if (enable_eke_leewave_dissipation ) then

        dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
        id  = ncvdef (ncid,'c_lee', NCFLOAT,3,dims,iret)
        name = 'Lee wave dissipation coefficient'; unit = '1/s'
        call dvcdf(ncid,id,name,len_trim(name),unit,16,spval)

        dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
        id  = ncvdef (ncid,'eke_lee_flux', NCFLOAT,3,dims,iret)
        name = 'lee wave flux'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,len_trim(name),unit,16,spval)

        dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
        id  = ncvdef (ncid,'c_Ri_diss', NCFLOAT,4,dims,iret)
        name = 'Interior dissipation coefficient'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)
       endif
       call ncclos (ncid, iret)
       print'(a)',' done '
 endif
 call fortran_barrier
end subroutine init_snap_cdf_eke




subroutine diag_snap_eke
!=======================================================================
!     write eke variables to NetCDF snapshot file
!=======================================================================
 use main_module   
 use eke_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_eke.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_eke.cdf'
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   ilen=ilen+1
   fxa = itt*dt_tracer/86400.0
   iret=nf_inq_varid(ncid,'Time',itimeid)
   iret= nf_put_vara_double(ncid,itimeid,ilen,1,fxa)
   call ncclos (ncid, iret)
 endif
 call fortran_barrier

 if (my_pe == 0 ) then
    iret=nf_open('pyOM_eke.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 ! Rossby radius
 bloc(is_pe:ie_pe,js_pe:je_pe) = L_rossby(is_pe:ie_pe,js_pe:je_pe)
 where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'L_Rossby',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 bloc(is_pe:ie_pe,js_pe:je_pe) = eke_bot_flux(is_pe:ie_pe,js_pe:je_pe)
 where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then  
      iret=nf_inq_varid(ncid,'eke_bot_flux',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)  
 endif

 if (enable_eke_leewave_dissipation ) then

    bloc(is_pe:ie_pe,js_pe:je_pe) = c_lee(is_pe:ie_pe,js_pe:je_pe)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe==0) then  
      iret=nf_inq_varid(ncid,'c_lee',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)  
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = eke_lee_flux(is_pe:ie_pe,js_pe:je_pe)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe==0) then  
      iret=nf_inq_varid(ncid,'eke_lee_flux',id)
      iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)  
    endif

 endif

 ! 3D fields

 do k=1,nz
    ! EKE
    bloc(is_pe:ie_pe,js_pe:je_pe) = eke(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'EKE',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    ! Rhines scale
    bloc(is_pe:ie_pe,js_pe:je_pe) = L_Rhines(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'L_Rhines',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif


    ! Eke dissipation
    bloc(is_pe:ie_pe,js_pe:je_pe) = eke_diss_iw(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'eke_diss_iw',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    ! Eke dissipation
    bloc(is_pe:ie_pe,js_pe:je_pe) = eke_diss_tke(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'eke_diss_tke',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    if (enable_eke_leewave_dissipation ) then
     bloc(is_pe:ie_pe,js_pe:je_pe) = c_Ri_diss(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe == 0 ) then 
        iret=nf_inq_varid(ncid,'c_Ri_diss',id)
        iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
    endif
 enddo

 if (my_pe==0)   iret = nf_close (ncid)
end subroutine diag_snap_eke





