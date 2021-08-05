



subroutine init_snap_cdf
!=======================================================================
!     initialize NetCDF snapshot file for basic variables
!=======================================================================
 use main_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,n
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33
 real*8 :: bloc(nx,ny)

 if (my_pe==0) print'(a)',' preparing file pyOM.cdf'

 call def_grid_cdf('pyOM.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)

      dims = (/Lon_tdim,lat_tdim,1,1/)
      id  = ncvdef (ncid,'ht', NCFLOAT,2,dims,iret)
      name = 'depth'; unit = 'm'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,itimedim/)
      id  = ncvdef (ncid,'temp', NCFLOAT,4,dims,iret)
      name = 'Temperature'; unit = 'deg C'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,itimedim/)
      id  = ncvdef (ncid,'salt', NCFLOAT,4,dims,iret)
      name = 'Salinity'; unit = 'g/kg'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'u', NCFLOAT,4,dims,iret)
      name = 'Zonal velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'v', NCFLOAT,4,dims,iret)
      name = 'Meridional velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'w', NCFLOAT,4,dims,iret)
      name = 'Vertical velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      if (enable_streamfunction) then
       dims = (/Lon_udim,lat_udim,iTimedim,1/)
       id  = ncvdef (ncid,'psi', NCFLOAT,3,dims,iret)
       name = 'Streamfunction'; unit = 'm^3/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_udim,lat_udim,1,1/)
       do n=1,nisle
        write(name, '("psi",i3)') n
        call replace_space_zero(name)
        id  = ncvdef (ncid,name, NCFLOAT,2,dims,iret)
        write(name, '("Boundary streamfunction ",i3)') n
        unit = 'm^3/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)   
       enddo
      else
       dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
       id  = ncvdef (ncid,'surf_press', NCFLOAT,3,dims,iret)
       name = 'Surface pressure'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'press', NCFLOAT,4,dims,iret)
       name = 'total pressure'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
      endif

      if (.not. enable_hydrostatic) then
       dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'p_hydro', NCFLOAT,4,dims,iret)
       name = 'Hydrostatic pressure'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
       id  = ncvdef (ncid,'p_non_hydro', NCFLOAT,4,dims,iret)
       name = 'Non hydrostatic pressure'; unit = 'm^2/s^2'
       call dvcdf(ncid,id,name,32,unit,16,spval)
      endif

      dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
      id  = ncvdef (ncid,'forc_temp_surface', NCFLOAT,3,dims,iret)
      name = 'Surface temperature flux'; unit = 'deg C m/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
      id  = ncvdef (ncid,'forc_salt_surface', NCFLOAT,3,dims,iret)
      name = 'Surface salinity flux'; unit = 'g/Kg m/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,itimedim,1/)
      id  = ncvdef (ncid,'taux', NCFLOAT,3,dims,iret)
      name = 'Surface wind stress'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_udim,itimedim,1/)
      id  = ncvdef (ncid,'tauy', NCFLOAT,3,dims,iret)
      name = 'Surface wind stress'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      call ncclos (ncid, iret)
 endif
   
 if (my_pe==0) iret=nf_open('pyOM.cdf',NF_WRITE,ncid)

 bloc(is_pe:ie_pe,js_pe:je_pe) = ht(is_pe:ie_pe,js_pe:je_pe)
 where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'ht',id)
    iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
 endif

 if (enable_streamfunction) then
  do n=1,nisle
   write(name, '("psi",i3)') n
   call replace_space_zero(name)
   bloc(is_pe:ie_pe,js_pe:je_pe) = psin(is_pe:ie_pe,js_pe:je_pe,n)
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe==0) then
    iret=nf_inq_varid(ncid,name,id)
    iret= nf_put_vara_double(ncid,id,(/1,1/), (/nx,ny/),bloc)
   endif
  enddo
 endif

 if (my_pe==0) iret=nf_close(ncid)
 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_cdf


subroutine diag_snap
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33

 if (my_pe == 0 ) then
   iret=nf_open('pyOM.cdf',NF_WRITE,ncid)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   ilen=ilen+1
   fxa = itt*dt_tracer/86400.0
   if (fxa <1.0) then
     print'(a,f12.2,a,i4)',' writing snapshot at ',fxa*86400,' s, time steps in file : ',ilen
   else
     print'(a,f12.2,a,i4)',' writing snapshot at ',fxa,' d, time steps in file : ',ilen
   endif
   iret=nf_inq_varid(ncid,'Time',itimeid)
   iret= nf_put_vara_double(ncid,itimeid,(/ilen/),(/1/),(/fxa/))
   call ncclos (ncid, iret)
 endif
 call fortran_barrier

 if (my_pe == 0 ) then
    iret=nf_open('pyOM.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 ! streamfunction or surface pressure
 bloc(is_pe:ie_pe,js_pe:je_pe) = psi(is_pe:ie_pe,js_pe:je_pe,tau)
 call pe0_recv_2D(nx,ny,bloc)
 if (enable_streamfunction.and.my_pe==0) then
    iret=nf_inq_varid(ncid,'psi',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif
 if (.not.enable_streamfunction.and.my_pe==0) then
    iret=nf_inq_varid(ncid,'surf_press',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif


 ! surface temperature flux
 bloc(is_pe:ie_pe,js_pe:je_pe) = forc_temp_surface(is_pe:ie_pe,js_pe:je_pe)
 where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'forc_temp_surface',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 ! surface salinity flux
 bloc(is_pe:ie_pe,js_pe:je_pe) = forc_salt_surface(is_pe:ie_pe,js_pe:je_pe)
 where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'forc_salt_surface',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 ! zonal wind stress
 bloc(is_pe:ie_pe,js_pe:je_pe) = surface_taux(is_pe:ie_pe,js_pe:je_pe)
 where( maskU(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'taux',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 ! meridional wind stress
 bloc(is_pe:ie_pe,js_pe:je_pe) = surface_tauy(is_pe:ie_pe,js_pe:je_pe)
 where( maskV(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe==0) then
    iret=nf_inq_varid(ncid,'tauy',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)
 endif

 ! 3D fields

 do k=1,nz

   if (.not. enable_hydrostatic) then
    ! hydrostatic pressure 
    bloc(is_pe:ie_pe,js_pe:je_pe) = p_hydro(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'p_hydro',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    ! non hydrostatic pressure 
    bloc(is_pe:ie_pe,js_pe:je_pe) = p_non_hydro(is_pe:ie_pe,js_pe:je_pe,k,taup1)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'p_non_hydro',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
   endif


   if (.not.enable_streamfunction) then
    ! total pressure 
    bloc(is_pe:ie_pe,js_pe:je_pe) = p_hydro(is_pe:ie_pe,js_pe:je_pe,k) + psi(is_pe:ie_pe,js_pe:je_pe,tau) 
    where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'press',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
   endif

   ! zonal velocity
   bloc(is_pe:ie_pe,js_pe:je_pe) = u(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskU(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'u',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! meridional velocity
   bloc(is_pe:ie_pe,js_pe:je_pe) = v(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskV(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'v',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! vertical velocity
   bloc(is_pe:ie_pe,js_pe:je_pe) = w(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'w',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! temperature
   bloc(is_pe:ie_pe,js_pe:je_pe) = temp(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'temp',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif

   ! salinity
   bloc(is_pe:ie_pe,js_pe:je_pe) = salt(is_pe:ie_pe,js_pe:je_pe,k,tau)
   where( maskT(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'salt',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
   endif
 enddo
 if (my_pe==0)   iret = nf_close (ncid)
end subroutine diag_snap



subroutine def_grid_cdf(filename)
!=======================================================================
!      Define standard grid in netcdf file
!=======================================================================
 use main_module   
 implicit none
 include "netcdf.inc"
 character*(*) filename
 integer :: ncid,iret,n
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lon_tid,lon_uid,itimeid
 integer :: lat_tdim,lat_udim,lat_uid,lat_tid
 integer :: z_tdim,z_udim,z_tid,z_uid
 character :: name*24, unit*16

 if (my_pe==0) then
      !iret = nf_create (filename, nf_clobber, ncid)
      iret = nf_create (filename, IOR(NF_CLOBBER,NF_64BIT_OFFSET) , ncid)
      if (iret /= 0) call halt_stop('NETCDF:'//nf_strerror(iret))
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
!     dimensions
      lon_tdim  = ncddef(ncid, 'xt', nx , iret)
      Lon_udim  = ncddef(ncid, 'xu', nx , iret)
      lat_tdim  = ncddef(ncid, 'yt', ny , iret)
      Lat_udim  = ncddef(ncid, 'yu', ny , iret)
      iTimedim  = ncddef(ncid, 'Time', nf_unlimited, iret)
!     grid variables
      Lon_tid  = ncvdef (ncid,'xt',NCFLOAT,1,(/lon_tdim/),iret)
      Lon_uid  = ncvdef (ncid,'xu',NCFLOAT,1,(/lon_udim/),iret)
      Lat_tid  = ncvdef (ncid,'yt',NCFLOAT,1,(/lat_tdim/),iret)
      Lat_uid  = ncvdef (ncid,'yu',NCFLOAT,1,(/lat_udim/),iret)
      itimeid  = ncvdef (ncid,'Time', NCFLOAT,1,(/itimedim/),iret)
!     attributes of the grid
      if (coord_degree) then
       name = 'Longitude on T grid     '; unit = 'degrees E'
       call ncaptc(ncid, Lon_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'Longitude on U grid     '; unit = 'degrees E'
       call ncaptc(ncid, Lon_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_uid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'Latitude on T grid     '; unit = 'degrees N'
       call ncaptc(ncid, Lat_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'Latitude on U grid     '; unit = 'degrees N'
       call ncaptc(ncid, Lat_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_uid, 'units',     NCCHAR, 16, unit, iret) 
      else
       name = 'zonal axis T grid     '; unit = 'km'
       call ncaptc(ncid, Lon_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'zonal axis U grid     '; unit = 'km'
       call ncaptc(ncid, Lon_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lon_uid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'meridional axis T grid'; unit = 'km'
       call ncaptc(ncid, Lat_tid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_tid, 'units',     NCCHAR, 16, unit, iret) 
       name = 'meridional axis U grid'; unit = 'km'
       call ncaptc(ncid, Lat_uid, 'long_name', NCCHAR, 24, name, iret) 
       call ncaptc(ncid, Lat_uid, 'units',     NCCHAR, 16, unit, iret) 
      endif

      name = 'Time '; unit = 'days'
      call ncaptc(ncid, itimeid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, itimeid, 'units',     NCCHAR, 16, unit, iret) 
      call ncaptc(ncid, iTimeid,'time_origin',NCCHAR, 20,'01-JAN-1900 00:00:00', iret)

      z_tdim    = ncddef(ncid, 'zt',  nz, iret)
      z_udim    = ncddef(ncid, 'zu',  nz, iret)
      z_tid  = ncvdef (ncid,'zt', NCFLOAT,1,(/z_tdim/),iret)
      z_uid  = ncvdef (ncid,'zu', NCFLOAT,1,(/z_udim/),iret)
      name = 'Height on T grid     '; unit = 'm'
      call ncaptc(ncid, z_tid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, z_tid, 'units',     NCCHAR, 16, unit, iret) 
      name = 'Height on U grid     '; unit = 'm'
      call ncaptc(ncid, z_uid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, z_uid, 'units',     NCCHAR, 16, unit, iret) 

      call ncendf(ncid, iret)
      iret= nf_put_vara_double(ncid,z_tid,(/1/),(/nz/) ,zt)
      iret= nf_put_vara_double(ncid,z_uid,(/1/),(/nz/) ,zw)
      call ncclos (ncid, iret)
 endif


 do n=0,n_pes-1
   call fortran_barrier
   if (my_pe==n) then
      iret=nf_open(filename,NF_WRITE,ncid)
      iret=nf_inq_varid(ncid,'xt',lon_tid)
      iret=nf_inq_varid(ncid,'xu',lon_uid)
      iret=nf_inq_varid(ncid,'yt',lat_tid)
      iret=nf_inq_varid(ncid,'yu',lat_uid)
      if (coord_degree) then
       iret= nf_put_vara_double(ncid,lon_Tid,(/is_pe/),(/ie_pe-is_pe+1/) ,xt(is_pe:ie_pe))
       iret= nf_put_vara_double(ncid,lon_uid,(/is_pe/),(/ie_pe-is_pe+1/) ,xu(is_pe:ie_pe))
       iret= nf_put_vara_double(ncid,lat_Tid,(/js_pe/),(/je_pe-js_pe+1/) ,yt(js_pe:je_pe))
       iret= nf_put_vara_double(ncid,lat_uid,(/js_pe/),(/je_pe-js_pe+1/) ,yu(js_pe:je_pe))
      else
       iret= nf_put_vara_double(ncid,lon_Tid,(/is_pe/),(/ie_pe-is_pe+1/) ,xt(is_pe:ie_pe)/1e3)
       iret= nf_put_vara_double(ncid,lon_uid,(/is_pe/),(/ie_pe-is_pe+1/) ,xu(is_pe:ie_pe)/1e3)
       iret= nf_put_vara_double(ncid,lat_Tid,(/js_pe/),(/je_pe-js_pe+1/) ,yt(js_pe:je_pe)/1e3)
       iret= nf_put_vara_double(ncid,lat_uid,(/js_pe/),(/je_pe-js_pe+1/) ,yu(js_pe:je_pe)/1e3)
      endif
      call ncclos (ncid, iret)
   endif
   call fortran_barrier
 enddo
end subroutine def_grid_cdf




subroutine dvcdf(ncid,ivarid,name,iname,unit,iunit,spval)
!=======================================================================
!     define some standard attributes of variable ivarid in NetCDF file ncid 
!=======================================================================
 implicit none
 integer ncid,ivarid,iname,iunit,iret
 character (len=*) :: name, unit
 real*8 :: spval, vv
 real :: v4=-1e33
 include "netcdf.inc"
 vv=spval
 call ncaptc(ncid,ivarid, 'long_name', NCCHAR,iname , name, iret) 
 if (iret.ne.0) print*,nf_strerror(iret)
 call ncaptc(ncid,ivarid, 'units',     NCCHAR,iunit, unit, iret) 
 if (iret.ne.0) print*,nf_strerror(iret)
 !call ncapt (ncid,ivarid, 'missing_value',NCDOUBLE,1,vv,iret)
 call ncapt (ncid,ivarid, 'missing_value',NCFLOAT,1,v4,iret)
 if (iret.ne.0) print*,nf_strerror(iret)
 !call ncapt (ncid,ivarid, '_FillValue', NCDOUBLE, 1,vv, iret)
 call ncapt (ncid,ivarid, '_FillValue', NCFLOAT, 1,v4, iret)
 if (iret.ne.0) print*,nf_strerror(iret)
end subroutine dvcdf




