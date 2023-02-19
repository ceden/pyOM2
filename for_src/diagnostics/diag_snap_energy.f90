



subroutine init_snap_cdf_energy
!=======================================================================
!     initialize NetCDF snapshot file of energy transfer fields
!=======================================================================
 use main_module   
 use isoneutral_module
 use diagnostics_module
 use idemix_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_energy.cdf'

 call def_grid_cdf('pyOM_energy.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_energy.cdf',NF_WRITE, ncid)
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
      id  = ncvdef (ncid,'K_diss_v', NCFLOAT,4,dims,iret)
      name = 'Dissipation of kin. En.'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'K_diss_bot', NCFLOAT,4,dims,iret)
      name = 'Dissipation of kin. En.'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'K_diss_h', NCFLOAT,4,dims,iret)
      name = 'Dissipation of kin. En.'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'P_diss_v', NCFLOAT,4,dims,iret)
      name = 'Dissipation of pot. En.'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'P_diss_nonlin', NCFLOAT,4,dims,iret)
      name = 'Dissipation of pot. En.'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      if (enable_neutral_diffusion) then
        dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
        id  = ncvdef (ncid,'P_diss_iso', NCFLOAT,4,dims,iret)
        name = 'Dissipation of pot. En.'; unit = 'm^2/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)
      endif

      if (enable_skew_diffusion) then
        dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
        id  = ncvdef (ncid,'P_diss_skew', NCFLOAT,4,dims,iret)
        name = 'Dissipation of pot. En.'; unit = 'm^2/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)
      endif
     
      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'P_diss_hmix', NCFLOAT,4,dims,iret)
      name = 'Dissipation of pot. En.'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'P_diss_adv', NCFLOAT,4,dims,iret)
      name = 'Dissipation of pot. En.'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      if (enable_TEM_friction.or.enable_biharmonic_thickness_mixing) then
       dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
       id  = ncvdef (ncid,'K_diss_gm', NCFLOAT,4,dims,iret)
       name = 'Dissipation of mean en.'; unit = 'm^2/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)
      endif

      if (enable_idemix .and. enable_idemix3) then
       dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
       id  = ncvdef (ncid,'K_diss_idemix', NCFLOAT,4,dims,iret)
       name = 'Dissipation of mean en.'; unit = 'm^2/s^3'
       call dvcdf(ncid,id,name,32,unit,16,spval)
       
      endif
    

      
      call ncclos (ncid, iret)
 endif
   
 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_cdf_energy


subroutine diag_snap_energy
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use isoneutral_module
 use diagnostics_module
 use idemix_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_energy.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_energy.cdf'
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
    iret=nf_open('pyOM_energy.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 do k=1,nz

   ! dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = K_diss_v(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'K_diss_v',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = K_diss_bot(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'K_diss_bot',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = K_diss_h(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'K_diss_h',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   if (enable_TEM_friction.or.enable_biharmonic_thickness_mixing) then
     ! dissipation
     bloc(is_pe:ie_pe,js_pe:je_pe) = K_diss_gm(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe == 0 ) then   
       iret=nf_inq_varid(ncid,'K_diss_gm',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
   endif

   ! dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = P_diss_v(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'P_diss_v',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = P_diss_hmix(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'P_diss_hmix',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = P_diss_nonlin(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'P_diss_nonlin',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   if (enable_neutral_diffusion) then
     ! dissipation
     bloc(is_pe:ie_pe,js_pe:je_pe) = P_diss_iso(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'P_diss_iso',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
   endif

   if (enable_skew_diffusion) then
     ! dissipation
     bloc(is_pe:ie_pe,js_pe:je_pe) = P_diss_skew(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe == 0 ) then 
         iret=nf_inq_varid(ncid,'P_diss_skew',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif
   endif

   ! dissipation
   bloc(is_pe:ie_pe,js_pe:je_pe) = P_diss_adv(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'P_diss_adv',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   if (enable_idemix .and. enable_idemix3) then

     ! dissipation
     bloc(is_pe:ie_pe,js_pe:je_pe) = K_diss_idemix(is_pe:ie_pe,js_pe:je_pe,k)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
     call pe0_recv_2D(nx,ny,bloc)
     if (my_pe == 0 ) then 
         iret=nf_inq_varid(ncid,'K_diss_idemix',id)
         iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
     endif

   endif
    
 enddo

 if (my_pe==0)   iret = nf_close (ncid)

end subroutine diag_snap_energy






