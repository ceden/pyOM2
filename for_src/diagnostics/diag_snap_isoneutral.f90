



subroutine init_snap_cdf_isoneutral
!=======================================================================
!     initialize NetCDF snapshot file for isoneutral variables
!=======================================================================
 use main_module   
 use isoneutral_module
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_isoneutral.cdf'

 call def_grid_cdf('pyOM_isoneutral.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_isoneutral.cdf',NF_WRITE, ncid)
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
      id  = ncvdef (ncid,'K_gm', NCFLOAT,4,dims,iret)
      name = 'skewness diffusivity'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)
       
      if (enable_neutral_diffusion .and. enable_skew_diffusion) then
       dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'B1_gm', NCFLOAT,4,dims,iret)
       name = 'Zonal comp. GM streamfct.'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'B2_gm', NCFLOAT,4,dims,iret)
       name = 'Meridional comp. GM streamfct.'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
      endif

      if (enable_TEM_friction) then
       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'kappa_gm', NCFLOAT,4,dims,iret)
       name = 'Vertical diffusivity'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
      endif

      if (enable_ml_para ) then
       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id  = ncvdef (ncid,'mld', NCFLOAT,3,dims,iret)
       name = 'Mixed layer depth'; unit = 'm'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id  = ncvdef (ncid,'Ri_ml', NCFLOAT,3,dims,iret)
       name = 'Richardson number'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'B1_ml', NCFLOAT,4,dims,iret)
       name = 'Zonal ML streamfct.'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'B2_ml', NCFLOAT,4,dims,iret)
       name = 'Meridional ML streamfct.'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'K_ml', NCFLOAT,4,dims,iret)
       name = 'Diffusivity'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'u_ml', NCFLOAT,4,dims,iret)
       name = 'ML eddy velocity'; unit = 'm/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'v_ml', NCFLOAT,4,dims,iret)
       name = 'ML eddy velocity'; unit = 'm/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'w_ml', NCFLOAT,4,dims,iret)
       name = 'ML eddy velocity'; unit = 'm/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'Nsqr', NCFLOAT,4,dims,iret)
       name = 'test'; unit = '1/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

      endif

      call ncclos (ncid, iret)
 endif
   
 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_cdf_isoneutral


subroutine diag_snap_isoneutral
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use isoneutral_module
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_isoneutral.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_isoneutral.cdf'
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
    iret=nf_open('pyOM_isoneutral.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 if (enable_ml_para ) then
   bloc(is_pe:ie_pe,js_pe:je_pe) = mld(is_pe:ie_pe,js_pe:je_pe)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then  
        iret=nf_inq_varid(ncid,'mld',id)
        iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)  
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = Ri_ml(is_pe:ie_pe,js_pe:je_pe)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then  
        iret=nf_inq_varid(ncid,'Ri_ml',id)
        iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)  
   endif


endif

 do k=1,nz

    ! GM diffusivity
    bloc(is_pe:ie_pe,js_pe:je_pe) = K_gm(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'K_gm',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
   if (enable_neutral_diffusion .and. enable_skew_diffusion) then
    bloc(is_pe:ie_pe,js_pe:je_pe) = B1_gm(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B1_gm',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = B2_gm(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B2_gm',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
   endif

   if (enable_TEM_friction) then
    ! vertical GM diffusivity
    bloc(is_pe:ie_pe,js_pe:je_pe) = kappa_gm(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'kappa_gm',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
   endif

   if (enable_ml_para ) then

    bloc(is_pe:ie_pe,js_pe:je_pe) = B1_ml(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B1_ml',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = B2_ml(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B2_ml',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = K_ml(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'K_ml',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = u_ml(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'u_ml',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
      
    bloc(is_pe:ie_pe,js_pe:je_pe) = v_ml(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'v_ml',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
      
    bloc(is_pe:ie_pe,js_pe:je_pe) = w_ml(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'w_ml',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
      
    bloc(is_pe:ie_pe,js_pe:je_pe) = Nsqr(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Nsqr',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
      
   endif
 enddo

 if (my_pe==0)   iret = nf_close (ncid)

end subroutine diag_snap_isoneutral








