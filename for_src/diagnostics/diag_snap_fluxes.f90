
module module_diag_snap_fluxes
 implicit none
 real*8, allocatable :: adv_flux_east_temp(:,:,:)
 real*8, allocatable :: adv_flux_east_salt(:,:,:)
 real*8, allocatable :: adv_flux_north_temp(:,:,:)
 real*8, allocatable :: adv_flux_north_salt(:,:,:)
 real*8, allocatable :: iso_flux_east_temp(:,:,:)
 real*8, allocatable :: iso_flux_east_salt(:,:,:)
 real*8, allocatable :: iso_flux_north_temp(:,:,:)
 real*8, allocatable :: iso_flux_north_salt(:,:,:)
 real*8, allocatable :: skew_flux_east_temp(:,:,:)
 real*8, allocatable :: skew_flux_east_salt(:,:,:)
 real*8, allocatable :: skew_flux_north_temp(:,:,:)
 real*8, allocatable :: skew_flux_north_salt(:,:,:)
end module module_diag_snap_fluxes


subroutine diag_snap_fluxes
 use main_module
 use isoneutral_module
 use module_diag_snap_fluxes
 implicit none
 real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz,3)
 include "netcdf.inc"
 integer :: ncid,iret,itdimid,ilen,itimeid,id,k
 real*8 :: cloc(nx,ny),fxa
 real*8, parameter :: spval = -1.0d33

 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,temp(:,:,:,tau),bloc)
 adv_flux_east_temp = flux_east
 adv_flux_north_temp = flux_north

 call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,salt(:,:,:,tau),bloc)
 adv_flux_east_salt = flux_east
 adv_flux_north_salt = flux_north

 

 if (enable_neutral_diffusion) then
  P_diss_iso = 0.0; dtemp_iso = 0.0; dsalt_iso = 0.0
  !call isoneutral_diffusion_pre
  bloc = temp
  call isoneutral_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,bloc,.true.)
  iso_flux_east_temp = flux_east
  iso_flux_north_temp = flux_north
  bloc = salt
  call isoneutral_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,bloc,.false.)
  iso_flux_east_salt = flux_east
  iso_flux_north_salt = flux_north
  if (enable_skew_diffusion) then
   P_diss_skew = 0.0;
   bloc = temp
   call isoneutral_skew_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,bloc,.true.)
   skew_flux_east_temp = flux_east
   skew_flux_north_temp = flux_north
   bloc = salt
   call isoneutral_skew_diffusion(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,bloc,.false.)
   skew_flux_east_salt = flux_east
   skew_flux_north_salt = flux_north
  endif
 endif



 if (my_pe == 0 ) then
   iret=nf_open('pyOM_fluxes.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_fluxes.cdf'
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
    iret=nf_open('pyOM_fluxes.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 do k=1,nz  
    cloc(is_pe:ie_pe,js_pe:je_pe) = adv_flux_east_temp(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'adv_flux_east_temp',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = adv_flux_east_salt(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'adv_flux_east_salt',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = adv_flux_north_temp(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'adv_flux_north_temp',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = adv_flux_north_salt(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'adv_flux_north_salt',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif

!!!!!!!!!!!!!

    cloc(is_pe:ie_pe,js_pe:je_pe) = iso_flux_east_temp(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'iso_flux_east_temp',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = iso_flux_east_salt(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'iso_flux_east_salt',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = iso_flux_north_temp(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'iso_flux_north_temp',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = iso_flux_north_salt(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'iso_flux_north_salt',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
!!!!!!!!!!!!!

    cloc(is_pe:ie_pe,js_pe:je_pe) = skew_flux_east_temp(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'skew_flux_east_temp',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = skew_flux_east_salt(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'skew_flux_east_salt',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = skew_flux_north_temp(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'skew_flux_north_temp',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    
    cloc(is_pe:ie_pe,js_pe:je_pe) = skew_flux_north_salt(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) cloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,cloc)
    if (my_pe == 0 ) then      
         iret=nf_inq_varid(ncid,'skew_flux_north_salt',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),cloc)
    endif
    

 enddo

 if (my_pe==0)   iret = nf_close (ncid)

end subroutine diag_snap_fluxes



subroutine init_diag_snap_fluxes
 use main_module
 use module_diag_snap_fluxes
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim, dims(4)
 character :: name*60, unit*16
 real*8, parameter :: spval = -1.0d33
 
 allocate(adv_flux_east_temp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 adv_flux_east_temp = 0.
 allocate(adv_flux_east_salt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 adv_flux_east_salt = 0.
 allocate(adv_flux_north_temp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 adv_flux_north_temp = 0.
 allocate(adv_flux_north_salt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 adv_flux_north_salt = 0.

 allocate(iso_flux_east_temp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 iso_flux_east_temp = 0.
 allocate(iso_flux_east_salt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 iso_flux_east_salt = 0.
 allocate(iso_flux_north_temp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 iso_flux_north_temp = 0.
 allocate(iso_flux_north_salt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 iso_flux_north_salt = 0.

 allocate(skew_flux_east_temp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 skew_flux_east_temp = 0.
 allocate(skew_flux_east_salt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 skew_flux_east_salt = 0.
 allocate(skew_flux_north_temp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 skew_flux_north_temp = 0.
 allocate(skew_flux_north_salt(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) )
 skew_flux_north_salt = 0.


 if (my_pe==0) print'(a)',' preparing file pyOM_fluxes.cdf'

 call def_grid_cdf('pyOM_fluxes.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_fluxes.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'adv_flux_east_temp', NCFLOAT,4,dims,iret)
      name = 'Eastward advective temperature flux'; unit = 'm K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'adv_flux_east_salt', NCFLOAT,4,dims,iret)
      name = 'Eastward advective salinity flux'; unit = 'm /s'
      call dvcdf(ncid,id,name,60,unit,16,spval)
      
      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'adv_flux_north_temp', NCFLOAT,4,dims,iret)
      name = 'Northward advective temperature flux'; unit = 'm K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'adv_flux_north_salt', NCFLOAT,4,dims,iret)
      name = 'Northward advective salinity flux'; unit = 'm /s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'iso_flux_east_temp', NCFLOAT,4,dims,iret)
      name = 'Eastward isopycnal temperature flux'; unit = 'm K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'iso_flux_east_salt', NCFLOAT,4,dims,iret)
      name = 'Eastward isopycnal salinity flux'; unit = 'm /s'
      call dvcdf(ncid,id,name,60,unit,16,spval)
      
      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'iso_flux_north_temp', NCFLOAT,4,dims,iret)
      name = 'Northward isopycnal temperature flux'; unit = 'm K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'iso_flux_north_salt', NCFLOAT,4,dims,iret)
      name = 'Northward isopycnal salinity flux'; unit = 'm /s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'skew_flux_east_temp', NCFLOAT,4,dims,iret)
      name = 'Eastward skew temperature flux'; unit = 'm K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'skew_flux_east_salt', NCFLOAT,4,dims,iret)
      name = 'Eastward skew salinity flux'; unit = 'm /s'
      call dvcdf(ncid,id,name,60,unit,16,spval)
      
      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'skew_flux_north_temp', NCFLOAT,4,dims,iret)
      name = 'Northward skew temperature flux'; unit = 'm K/s'
      call dvcdf(ncid,id,name,60,unit,16,spval)

      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'skew_flux_north_salt', NCFLOAT,4,dims,iret)
      name = 'Northward skew salinity flux'; unit = 'm /s'
      call dvcdf(ncid,id,name,60,unit,16,spval)



      call ncclos (ncid, iret)
      print'(a)',' done '
 endif
 call fortran_barrier
 
 
end subroutine init_diag_snap_fluxes
