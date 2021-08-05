


  
subroutine init_snap_cdf_tke
!=======================================================================
!     initialize NetCDF snapshot file for tke variables
!=======================================================================
 use main_module   
 use tke_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_tke.cdf'

 call def_grid_cdf('pyOM_tke.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_tke.cdf',NF_WRITE, ncid)
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
      id  = ncvdef (ncid,'kappaH', NCFLOAT,4,dims,iret)
      name = 'Vertical diffusivity'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)
      
      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'tke', NCFLOAT,4,dims,iret)
      name = 'Turbulent kinetic energy'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'Prandtl', NCFLOAT,4,dims,iret)
      name = 'Prandtl number'; unit = ' '
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'mxl', NCFLOAT,4,dims,iret)
      name = 'Mixing length'; unit = 'm'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'tke_diss', NCFLOAT,4,dims,iret)
      name = 'TKE dissipation'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
      id  = ncvdef (ncid,'forc_tke', NCFLOAT,3,dims,iret)
      name = 'TKE surface flux'; unit = 'm^3/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
      id  = ncvdef (ncid,'tke_surf_corr', NCFLOAT,3,dims,iret)
      name = 'Correction of TKE surface flux'; unit = 'm^3/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      call ncclos (ncid, iret)
 endif

 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_cdf_tke




subroutine diag_snap_tke
!=======================================================================
!     write tke variables to NetCDF snapshot file
!=======================================================================
 use main_module   
 use tke_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_tke.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_tke.cdf'
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   ilen=ilen+1
   fxa = itt*dt_tracer/86400.0
   iret=nf_inq_varid(ncid,'Time',itimeid)
   iret= nf_put_vara_double(ncid,itimeid,(/ilen/),(/1/),(/fxa/))
   call ncclos (ncid, iret)
 endif
 call fortran_barrier

 if (my_pe == 0 ) then
    iret=nf_open('pyOM_tke.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 ! TKE forcing
 bloc(is_pe:ie_pe,js_pe:je_pe) = forc_tke_surface(is_pe:ie_pe,js_pe:je_pe)
 where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'forc_tke',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 ! surface correction of TKE forcing
 bloc(is_pe:ie_pe,js_pe:je_pe) = tke_surf_corr(is_pe:ie_pe,js_pe:je_pe)
 where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'tke_surf_corr',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 ! 3D fields

 do k=1,nz

   ! vertical diffusivity
   bloc(is_pe:ie_pe,js_pe:je_pe) = kappaH(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'kappaH',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

    ! TKE
    bloc(is_pe:ie_pe,js_pe:je_pe) = tke(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'tke',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
   
    ! Prandtl zahl
    bloc(is_pe:ie_pe,js_pe:je_pe) = Prandtlnumber(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Prandtl',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    ! mixing length
    bloc(is_pe:ie_pe,js_pe:je_pe) = mxl(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'mxl',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    ! TKE dissipation
    bloc(is_pe:ie_pe,js_pe:je_pe) = tke_diss(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'tke_diss',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
 enddo

 if (my_pe==0)   iret = nf_close (ncid)

end subroutine diag_snap_tke


