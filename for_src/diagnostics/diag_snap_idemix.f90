



subroutine init_snap_cdf_idemix
!=======================================================================
!     initialize NetCDF snapshot file for idemix variables
!=======================================================================
 use main_module   
 use idemix_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33
 real*8 :: bloc(nx,ny)
 

 if (my_pe==0) print'(a)',' preparing file pyOM_idemix.cdf'

 call def_grid_cdf('pyOM_idemix.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_idemix.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)
             



      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'Nsqr', NCFLOAT,4,dims,iret)
      name = 'Square of stability frequency'; unit = '1/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'E_iw', NCFLOAT,4,dims,iret)
      name = 'Internal wave energy'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'c0', NCFLOAT,4,dims,iret)
      name = 'vertical IW group velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
      id  = ncvdef (ncid,'cstar', NCFLOAT,3,dims,iret)
      name = 'modal gravity wave speed'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'v0', NCFLOAT,4,dims,iret)
      name = 'hor. IW group velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'iw_diss', NCFLOAT,4,dims,iret)
      name = 'IW dissipation'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      if (enable_idemix3) then
        dims = (/Lon_udim,lat_tdim,1,1/)
        id  = ncvdef (ncid,'forc_iw_surface_u', NCFLOAT,2,dims,iret)
        name = 'IW surface forcing'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        dims = (/Lon_udim,lat_tdim,1,1/)
        id  = ncvdef (ncid,'forc_iw_bottom_u', NCFLOAT,2,dims,iret)
        name = 'IW bottom forcing'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        dims = (/Lon_tdim,lat_udim,1,1/)
        id  = ncvdef (ncid,'forc_iw_surface_v', NCFLOAT,2,dims,iret)
        name = 'IW surface forcing'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        dims = (/Lon_tdim,lat_udim,1,1/)
        id  = ncvdef (ncid,'forc_iw_bottom_v', NCFLOAT,2,dims,iret)
        name = 'IW bottom forcing'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)
      else
        dims = (/Lon_tdim,lat_tdim,1,1/)
        id  = ncvdef (ncid,'forc_iw_surface', NCFLOAT,2,dims,iret)
        name = 'IW surface forcing'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        dims = (/Lon_tdim,lat_tdim,1,1/)
        id  = ncvdef (ncid,'forc_iw_bottom', NCFLOAT,2,dims,iret)
        name = 'IW bottom forcing'; unit = 'm^3/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)
      endif
       
      if (enable_idemix3) then

         ! IDEMIX 3.0
         
       dims = (/Lon_udim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'Tc_U', NCFLOAT,4,dims,iret)
       name = 'inverse time scale'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_udim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'Tc_V', NCFLOAT,4,dims,iret)
       name = 'inverse time scale'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'Omega', NCFLOAT,4,dims,iret)
       name = 'factor'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'F_s', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_udim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'F_d', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'G_s', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_udim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'G_d', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'O_s', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_udim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'O_d', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'P_s', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_udim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'P_d', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_udim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'A_u', NCFLOAT,4,dims,iret)
       name = 'vertical viscosity for u'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_udim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'A_v', NCFLOAT,4,dims,iret)
       name = 'vertical viscosity for v'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
      endif
     


      if (enable_leewaves) then
      
       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id = ncvdef(ncid,'Lee_flux',NCFLOAT,3,dims,iret)
       name = 'Leewave energy flux' ; unit = 'm^3/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id = ncvdef(ncid,'inv_fr',NCFLOAT,3,dims,iret)
       name = 'inverse Froude number' ; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,1,1/)
       id = ncvdef(ncid,'h_rms',NCFLOAT,2,dims,iret)
       name = 'RMS height' ; unit = 'm'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_tdim,1,1/)
       id = ncvdef(ncid,'k_s',NCFLOAT,2,dims,iret)
       name = 'Topographic wavenumber' ; unit = 'm^-1'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,1,1/)
       id = ncvdef(ncid,'k_n',NCFLOAT,2,dims,iret)
       name = 'Topographic wavenumber' ; unit = 'm^-1'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_tdim,1,1/)
       id = ncvdef(ncid,'phi_s',NCFLOAT,2,dims,iret)
       name = 'Angle' ; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id = ncvdef(ncid,'u_bottom',NCFLOAT,3,dims,iret)
       name = 'Zonal bottom velocity' ; unit = 'm/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id = ncvdef(ncid,'v_bottom',NCFLOAT,3,dims,iret)
       name = 'Meridional bottom velocity' ; unit = 'm/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id = ncvdef(ncid,'stress_lee_u',NCFLOAT,3,dims,iret)
       name = 'Bottom stress' ; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id = ncvdef(ncid,'stress_lee_v',NCFLOAT,3,dims,iret)
       name = 'Bottom stress' ; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'E_lee_p', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'E_lee_m', NCFLOAT,4,dims,iret)
       name = 'Internal wave energy'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'Lambda_0', NCFLOAT,4,dims,iret)
       name = 'Parameter'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'Tc_lee', NCFLOAT,4,dims,iret)
       name = 'Parameter'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
       id  = ncvdef (ncid,'iw_diss_lee', NCFLOAT,4,dims,iret)
       name = 'Dissipation of lee waves'; unit = 'm^2/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)

      endif


 
      if (enable_idemix_M2) then

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'E_M2', NCFLOAT,3,dims,iret)
       name = 'M2 Energy'; unit = 'm^3/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_tdim,itimedim/)
       id  = ncvdef (ncid,'E_struct_M2', NCFLOAT,4,dims,iret)
       name = 'M2 structure function'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'cg_M2', NCFLOAT,3,dims,iret)
       name = 'M2 group velocity'; unit = 'm/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_udim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'kdot_x_M2', NCFLOAT,3,dims,iret)
       name = 'M2 refraction'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_udim,itimedim,1/)
       id  = ncvdef (ncid,'kdot_y_M2', NCFLOAT,3,dims,iret)
       name = 'M2 refraction'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'tau_M2', NCFLOAT,3,dims,iret)
       name = 'M2 decay time scale'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'alpha_M2_cont', NCFLOAT,3,dims,iret)
       name = 'M2-continuum coupling coef'; unit = 's/m^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'forc_M2', NCFLOAT,3,dims,iret)
       name = 'M2 forcing'; unit = 'm^3/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)

      endif

      if (enable_idemix_niw) then
       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'E_NIW', NCFLOAT,3,dims,iret)
       name = 'NIW Energy'; unit = 'm^3/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_tdim,itimedim/)
       id  = ncvdef (ncid,'E_struct_NIW', NCFLOAT,4,dims,iret)
       name = 'NIW structure function'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'cg_NIW', NCFLOAT,3,dims,iret)
       name = 'NIW group velocity'; unit = 'm/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_udim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'kdot_x_NIW', NCFLOAT,3,dims,iret)
       name = 'NIW refraction'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_udim,itimedim,1/)
       id  = ncvdef (ncid,'kdot_y_NIW', NCFLOAT,3,dims,iret)
       name = 'NIW refraction'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'tau_NIW', NCFLOAT,3,dims,iret)
       name = 'NIW decay time scale'; unit = '1/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,itimedim,1/)
       id  = ncvdef (ncid,'forc_NIW', NCFLOAT,3,dims,iret)
       name = 'NIW forcing'; unit = 'm^3/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)

      endif

      call ncclos (ncid, iret)
 endif
   
 if (my_pe==0) iret=nf_open('pyOM_idemix.cdf',NF_WRITE,ncid)

 if (enable_idemix3) then

   bloc(is_pe:ie_pe,js_pe:je_pe) = forc_iw_surface_u(is_pe:ie_pe,js_pe:je_pe)
   where( maskU(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
      iret=nf_inq_varid(ncid,'forc_iw_surface_u',id)
      iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = forc_iw_bottom_u(is_pe:ie_pe,js_pe:je_pe)
   where( maskU(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
      iret=nf_inq_varid(ncid,'forc_iw_bottom_u',id)
      iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = forc_iw_surface_v(is_pe:ie_pe,js_pe:je_pe)
   where( maskV(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
      iret=nf_inq_varid(ncid,'forc_iw_surface_v',id)
      iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = forc_iw_bottom_v(is_pe:ie_pe,js_pe:je_pe)
   where( maskV(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
      iret=nf_inq_varid(ncid,'forc_iw_bottom_v',id)
      iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
   endif

 else
   bloc(is_pe:ie_pe,js_pe:je_pe) = forc_iw_surface(is_pe:ie_pe,js_pe:je_pe)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
      iret=nf_inq_varid(ncid,'forc_iw_surface',id)
      iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = forc_iw_bottom(is_pe:ie_pe,js_pe:je_pe)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
      iret=nf_inq_varid(ncid,'forc_iw_bottom',id)
      iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
   endif
 endif


 if (enable_leewaves) then
 
    bloc(is_pe:ie_pe,js_pe:je_pe) = h_rms(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe==0) then
       iret=nf_inq_varid(ncid,'h_rms',id)
       iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = k_s(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe==0) then
       iret=nf_inq_varid(ncid,'k_s',id)
       iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
    endif
    
    bloc(is_pe:ie_pe,js_pe:je_pe) = k_n(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe==0) then
       iret=nf_inq_varid(ncid,'k_n',id)
       iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
    endif
 
    bloc(is_pe:ie_pe,js_pe:je_pe) = ph_s(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe==0) then
       iret=nf_inq_varid(ncid,'phi_s',id)
       iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
    endif
    
 endif


 if (my_pe==0) iret=nf_close(ncid)
 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_cdf_idemix


subroutine diag_snap_idemix
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use idemix_module   
 use diagnostics_module

 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33
 
 
 if (my_pe == 0 ) then
   iret=nf_open('pyOM_idemix.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_idemix.cdf'
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
    iret=nf_open('pyOM_idemix.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 ! cstar
 bloc(is_pe:ie_pe,js_pe:je_pe) = cstar(is_pe:ie_pe,js_pe:je_pe)
 where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'cstar',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 if (enable_idemix_M2) then

  bloc = 0d0
  do k=2,np-1
    bloc(is_pe:ie_pe,js_pe:je_pe) = bloc(is_pe:ie_pe,js_pe:je_pe) + &
           E_M2(is_pe:ie_pe,js_pe:je_pe,k,tau)*dphit(k)*maskTp(is_pe:ie_pe,js_pe:je_pe,k)
  enddo
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'E_M2',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = cg_M2(is_pe:ie_pe,js_pe:je_pe)
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'cg_M2',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = kdot_x_M2(is_pe:ie_pe,js_pe:je_pe)
  where( maskU(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'kdot_x_M2',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = kdot_y_M2(is_pe:ie_pe,js_pe:je_pe)
  where( maskV(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'kdot_y_M2',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = tau_M2(is_pe:ie_pe,js_pe:je_pe)
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'tau_M2',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = alpha_M2_cont(is_pe:ie_pe,js_pe:je_pe)
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'alpha_M2_cont',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc = 0d0
  do k=2,np-1
    bloc(is_pe:ie_pe,js_pe:je_pe) = bloc(is_pe:ie_pe,js_pe:je_pe) + &
           forc_M2(is_pe:ie_pe,js_pe:je_pe,k)*dphit(k)*maskTp(is_pe:ie_pe,js_pe:je_pe,k)
  enddo
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'forc_M2',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

 endif

 if (enable_idemix_niw) then

  bloc = 0d0
  do k=2,np-1
    bloc(is_pe:ie_pe,js_pe:je_pe) = bloc(is_pe:ie_pe,js_pe:je_pe) + &
           E_niw(is_pe:ie_pe,js_pe:je_pe,k,tau)*dphit(k)*maskTp(is_pe:ie_pe,js_pe:je_pe,k)
  enddo
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'E_NIW',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = cg_niw(is_pe:ie_pe,js_pe:je_pe)
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'cg_NIW',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = kdot_x_niw(is_pe:ie_pe,js_pe:je_pe)
  where( maskU(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'kdot_x_NIW',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = kdot_y_niw(is_pe:ie_pe,js_pe:je_pe)
  where( maskV(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'kdot_y_NIW',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc(is_pe:ie_pe,js_pe:je_pe) = tau_niw(is_pe:ie_pe,js_pe:je_pe)
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'tau_NIW',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif

  bloc = 0d0
  do k=2,np-1
    bloc(is_pe:ie_pe,js_pe:je_pe) = bloc(is_pe:ie_pe,js_pe:je_pe) + &
           forc_niw(is_pe:ie_pe,js_pe:je_pe,k)*dphit(k)*maskTp(is_pe:ie_pe,js_pe:je_pe,k)
  enddo
  where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
  call pe0_recv_2D(nx,ny,bloc)
  if (my_pe==0) then
    iret=nf_inq_varid(ncid,'forc_NIW',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
  endif
 endif


 ! 3D fields

 do k=1,nz


   ! stability frequency
   bloc(is_pe:ie_pe,js_pe:je_pe) = Nsqr(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Nsqr',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! Internal wave energy
   bloc(is_pe:ie_pe,js_pe:je_pe) = E_iw(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'E_iw',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! Internal wave dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = iw_diss(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'iw_diss',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! Internal wave vertical group velocity
   bloc(is_pe:ie_pe,js_pe:je_pe) = c0(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'c0',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! Internal wave horizontal group velocity
   bloc(is_pe:ie_pe,js_pe:je_pe) = v0(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'v0',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   if (enable_idemix3) then

      ! IDEMIX 3.0

    bloc(is_pe:ie_pe,js_pe:je_pe) = Tc_U(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Tc_U',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = Tc_V(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Tc_V',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
      
    bloc(is_pe:ie_pe,js_pe:je_pe) = 0.5*(Om_id3_u(is_pe:ie_pe,js_pe:je_pe,k)& 
                                        +Om_id3_u(is_pe-1:ie_pe-1,js_pe:je_pe,k))
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Omega',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
      
    bloc(is_pe:ie_pe,js_pe:je_pe) = F_s(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'F_s',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
    bloc(is_pe:ie_pe,js_pe:je_pe) = F_d(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'F_d',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = G_s(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'G_s',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
    bloc(is_pe:ie_pe,js_pe:je_pe) = G_d(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'G_d',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = O_s(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'O_s',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
    bloc(is_pe:ie_pe,js_pe:je_pe) = O_d(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'O_d',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = P_s(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'P_s',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
    bloc(is_pe:ie_pe,js_pe:je_pe) = P_d(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'P_d',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
    bloc(is_pe:ie_pe,js_pe:je_pe) = viscU_idemix(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'A_u',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = viscV_idemix(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'A_v',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
   endif ! enable_idemix3
   
   if (enable_idemix_M2) then
    bloc(is_pe:ie_pe,js_pe:je_pe) = E_struct_M2(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'E_struct_M2',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
   endif

   if (enable_idemix_niw) then
    bloc(is_pe:ie_pe,js_pe:je_pe) = E_struct_niw(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'E_struct_NIW',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
   endif
   
   if (enable_leewaves) then

    bloc(is_pe:ie_pe,js_pe:je_pe) = E_lee_p(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'E_lee_p',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = E_lee_m(is_pe:ie_pe,js_pe:je_pe,k,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'E_lee_m',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = Lambda_0(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Lambda_0',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = Tc_lee(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'Tc_lee',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif


    bloc(is_pe:ie_pe,js_pe:je_pe) = iw_diss_lee(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'iw_diss_lee',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

   endif
     
 enddo ! k=1,nz



   if (enable_leewaves) then
   
    bloc(is_pe:ie_pe,js_pe:je_pe) = flux_lee_tot(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0) then
       iret=nf_inq_varid(ncid,'Lee_flux',id)
       iret=nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),bloc)
    endif
 
    bloc(is_pe:ie_pe,js_pe:je_pe) = inv_fr(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0) then
       iret=nf_inq_varid(ncid,'inv_fr',id)
       iret=nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),bloc)
    endif
    
    bloc(is_pe:ie_pe,js_pe:je_pe) = u_bot(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0) then
       iret=nf_inq_varid(ncid,'u_bottom',id)
       iret=nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = v_bot(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0) then
       iret=nf_inq_varid(ncid,'v_bottom',id)
       iret=nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = stress_lee_u(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0) then
       iret=nf_inq_varid(ncid,'stress_lee_u',id)
       iret=nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = stress_lee_v(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0) then
       iret=nf_inq_varid(ncid,'stress_lee_v',id)
       iret=nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),bloc)
    endif


   endif ! enable_leewaves


if (my_pe ==0) iret=nf_close(ncid)

end subroutine diag_snap_idemix






