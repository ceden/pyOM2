



subroutine init_snap_biharmonic_thickness
!=======================================================================
!     initialize NetCDF snapshot file for basic variables
!=======================================================================
 use main_module   
 use biharmonic_thickness_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33


 if ( .not. biharmonic_thickness_active) return
 
 if (my_pe==0) print'(a)',' preparing file pyOM_biha_thk.cdf'
 
 call def_grid_cdf('pyOM_biha_thk.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_biha_thk.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)
   
      id  = ncvdef (ncid,'K_thkbi', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
      name = 'Viscosity'; unit = 'm^4/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)
     
      id  = ncvdef (ncid,'B_x', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
      name = 'Streamfunction'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)
     
      id  = ncvdef (ncid,'B_y', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
      name = 'Streamfunction'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)
     
      id  = ncvdef (ncid,'mld', NCFLOAT,3,(/Lon_tdim,lat_tdim,iTimedim/),iret)
      name = 'Mixed layer depth'; unit = 'm'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      if (enable_biharmonic_thickness_backscatter_skew .or. &
          enable_biharmonic_thickness_backscatter_harmonic ) then
       id  = ncvdef (ncid,'A_thk', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
       name = 'Viscosity'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       id  = ncvdef (ncid,'B_x_back', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
       name = 'Streamfunction'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       id  = ncvdef (ncid,'B_y_back', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
       name = 'Streamfunction'; unit = 'm^2/s'
       call dvcdf(ncid,id,name,32,unit,16,spval)

       if (enable_conserve_energy) then 
        id  = ncvdef (ncid,'diss_back', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
        name = 'Dissipation'; unit = 'm^2/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval)      
       endif
             
       if (enable_biharmonic_thickness_backscatter_integrate_energy) then
         id  = ncvdef (ncid,'E_back', NCFLOAT,3,(/Lon_tdim,lat_tdim,iTimedim/),iret)
         name = 'Energy'; unit = 'm^2/s^2'
         call dvcdf(ncid,id,name,32,unit,16,spval)      
       endif       
      endif
  
      if (enable_biharmonic_thickness_bolus_skew .and. enable_conserve_energy) then       
        id  = ncvdef (ncid,'diss_biha_skew', NCFLOAT,4,(/Lon_tdim,lat_tdim,z_udim,iTimedim/),iret)
        name = 'Dissipation'; unit = 'm^2/s^3'
        call dvcdf(ncid,id,name,32,unit,16,spval) 
      endif 

      call ncclos (ncid, iret)
 endif
  
 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_biharmonic_thickness



subroutine diag_snap_biharmonic_thickness
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use biharmonic_thickness_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret
 real*8 :: bloc(nx,ny),fxa
 integer :: itdimid,ilen,itimeid,id,k
 real*8, parameter :: spval = -1.0d33

 if ( .not. biharmonic_thickness_active) return

 if (my_pe == 0 ) then
   print'(a)',' writing to file pyOM_biha_thk.cdf'
   iret=nf_open('pyOM_biha_thk.cdf',NF_WRITE,ncid)
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
    iret=nf_open('pyOM_biha_thk.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)
 endif

 bloc(is_pe:ie_pe,js_pe:je_pe) = mld(is_pe:ie_pe,js_pe:je_pe)
 where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
 call pe0_recv_2D(nx,ny,bloc)
 if (my_pe == 0 ) then 
       iret = nf_inq_varid(ncid,'mld',id)      
       iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)     
 endif

 if (enable_biharmonic_thickness_backscatter_skew .or. &
     enable_biharmonic_thickness_backscatter_harmonic ) then
  if (enable_biharmonic_thickness_backscatter_integrate_energy) then
    bloc(is_pe:ie_pe,js_pe:je_pe) = e_back(is_pe:ie_pe,js_pe:je_pe,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret = nf_inq_varid(ncid,'E_back',id)      
       iret= nf_put_vara_double(ncid,id,(/1,1,ilen/), (/nx,ny,1/),bloc)     
    endif 
  endif
 endif

 do k=1,nz

   bloc(is_pe:ie_pe,js_pe:je_pe) = K_thkbi(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'K_thkbi',id)      
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)       
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = diag_Bx(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B_x',id)      
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)       
   endif

   bloc(is_pe:ie_pe,js_pe:je_pe) = diag_By(is_pe:ie_pe,js_pe:je_pe,k)
   where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
   call pe0_recv_2D(nx,ny,bloc)
   if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B_y',id)      
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)       
   endif

   if (enable_biharmonic_thickness_backscatter_skew  .or. &
       enable_biharmonic_thickness_backscatter_harmonic) then
   
    bloc(is_pe:ie_pe,js_pe:je_pe) = A_thk(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'A_thk',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = diag_Bx_back(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B_x_back',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif

    bloc(is_pe:ie_pe,js_pe:je_pe) = diag_By_back(is_pe:ie_pe,js_pe:je_pe,k)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
    call pe0_recv_2D(nx,ny,bloc)
    if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'B_y_back',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
    endif
    
    if (enable_conserve_energy) then
      bloc(is_pe:ie_pe,js_pe:je_pe) = diss_back(is_pe:ie_pe,js_pe:je_pe,k)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'diss_back',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif  
    endif
    
   endif

   if (enable_biharmonic_thickness_bolus_skew .and. enable_conserve_energy) then
      bloc(is_pe:ie_pe,js_pe:je_pe) = diss_biha_skew(is_pe:ie_pe,js_pe:je_pe,k)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,k) == 0.) bloc(is_pe:ie_pe,js_pe:je_pe) = spval
      call pe0_recv_2D(nx,ny,bloc)
      if (my_pe == 0 ) then 
       iret=nf_inq_varid(ncid,'diss_biha_skew',id)
       iret= nf_put_vara_double(ncid,id,(/1,1,k,ilen/), (/nx,ny,1,1/),bloc)
      endif       
   endif 

 enddo
 if (my_pe==0)   iret = nf_close (ncid)
end subroutine diag_snap_biharmonic_thickness




subroutine biharmonic_thickness_diag_store_tensor(is_,ie_,js_,je_,nz_,Bx_loc,By_loc)
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use biharmonic_thickness_module
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: Bx_loc(is_:ie_,js_:je_,nz_), By_loc(is_:ie_,js_:je_,nz_)
 
 diag_Bx= Bx_loc
 diag_By= By_loc
end subroutine biharmonic_thickness_diag_store_tensor


subroutine biharmonic_thickness_diag_store_tensor_back(is_,ie_,js_,je_,nz_,Bx_loc,By_loc)
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use biharmonic_thickness_module
 implicit none
 integer :: is_,ie_,js_,je_,nz_
 real*8 :: Bx_loc(is_:ie_,js_:je_,nz_), By_loc(is_:ie_,js_:je_,nz_)

 if (enable_biharmonic_thickness_backscatter_skew) then
  diag_Bx_back= Bx_loc
  diag_By_back= By_loc
 endif
end subroutine biharmonic_thickness_diag_store_tensor_back



