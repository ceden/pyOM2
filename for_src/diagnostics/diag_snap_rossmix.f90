



subroutine init_snap_cdf_rossmix
!=======================================================================
!     initialize NetCDF snapshot file for rossmix variables
!=======================================================================
 use main_module   
 use rossmix_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,n
 integer :: lon_tdim,lon_udim,itimedim
 integer :: lat_tdim,lat_udim,id,z_tdim,z_udim,p_tdim,p_udim,p_tid,p_uid
 integer :: dims(4)
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33

 if (my_pe==0) print'(a)',' preparing file pyOM_rossmix.cdf'

 call def_grid_cdf('pyOM_rossmix.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_rossmix.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      if (iret /= 0) call halt_stop('NETCDF:'//nf_strerror(iret))
!     dimensions
      p_tdim  = ncddef(ncid, 'phit', nphi , iret)
      p_udim  = ncddef(ncid, 'phiu', nphi , iret)
!     grid variables
      p_tid  = ncvdef (ncid,'phit',NCFLOAT,1,p_tdim,iret)
      name = 'wave angle t grid'; unit = ' '
      call ncaptc(ncid, p_tid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, p_tid, 'units',     NCCHAR, 16, unit, iret) 
      p_uid  = ncvdef (ncid,'phiu',NCFLOAT,1,p_udim,iret)
      name = 'wave angle u grid'; unit = ' '
      call ncaptc(ncid, p_uid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, p_uid, 'units',     NCCHAR, 16, unit, iret) 
!     attributes of the grid
      call ncendf(ncid, iret)
      iret= nf_put_vara_double(ncid,p_Tid,1,nphi ,phit)
      iret= nf_put_vara_double(ncid,p_uid,1,nphi ,phiu)
      call ncclos (ncid, iret)
  endif

 if (my_pe==0) then
      iret=nf_open('pyOM_rossmix.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',lon_tdim)
      iret=nf_inq_dimid(ncid,'xu',lon_udim)
      iret=nf_inq_dimid(ncid,'yt',lat_tdim)
      iret=nf_inq_dimid(ncid,'yu',lat_udim)
      iret=nf_inq_dimid(ncid,'phit',p_tdim)
      iret=nf_inq_dimid(ncid,'phiu',p_udim)
      iret=nf_inq_dimid(ncid,'zt',z_tdim)
      iret=nf_inq_dimid(ncid,'zu',z_udim)
      iret=nf_inq_dimid(ncid,'Time',itimedim)

      dims = (/Lon_tdim,lat_tdim,z_udim,itimedim/)
      id  = ncvdef (ncid,'Nsqr', NCFLOAT,4,dims,iret)
      name = 'Stability frequency'; unit = '1/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      do n=1,nmodes
        write(name, '("E_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim, p_tdim, iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'Sum of Energy'; unit = 'm^3/s^2'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("E_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim, p_tdim, iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'Sum of Energy'; unit = 'm^3/s^2'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("cgu_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_udim,lat_tdim,p_tdim, iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'group velocity in x direction'; unit = 'm/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("cgu_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_udim,lat_tdim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'group velocity in x direction'; unit = 'm/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("cgv_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_udim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'group velocity in y direction'; unit = 'm/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("cgv_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_udim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'group velocity in y direction'; unit = 'm/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("phidot_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_udim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'group velocity in phi direction'; unit = 'm/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("phidot_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_udim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'group velocity in phi direction'; unit = 'm/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("gamma_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'frequency'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("gamma_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'frequency'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("nu_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'frequency'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("nu_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'frequency'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)


        write(name, '("lambda_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'frequency'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("lambda_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,p_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'frequency'; unit = '1/s'
        call dvcdf(ncid,id,name,32,unit,16,spval)


        write(name, '("phin_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'structure function'; unit = '1/m'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("R_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,iTimedim,1/)
        id  = ncvdef (ncid,name, NCFLOAT,3,dims,iret)
        name = 'Rossby radius'; unit = 'm'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("phi_ga_l_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'structure function'; unit = '1/m'
        call dvcdf(ncid,id,name,32,unit,16,spval)

        write(name, '("phi_ga_s_",i1)') n-1; !call replace_space_zero(name)
        dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
        id  = ncvdef (ncid,name, NCFLOAT,4,dims,iret)
        name = 'structure function'; unit = '1/m'
        call dvcdf(ncid,id,name,32,unit,16,spval)

      enddo

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'Bx', NCFLOAT,4,dims,iret)
      name = 'Streamfunction'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'By', NCFLOAT,4,dims,iret)
      name = 'Streamfunction'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'stress_ux', NCFLOAT,4,dims,iret)
      name = 'Zonal momentum flux'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'stress_vx', NCFLOAT,4,dims,iret)
      name = 'Meridional momentum flux'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'stress_uy', NCFLOAT,4,dims,iret)
      name = 'Zonal momentum flux'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'stress_vy', NCFLOAT,4,dims,iret)
      name = 'Meridional momentum flux'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'K_diss_ross', NCFLOAT,4,dims,iret)
      name = 'Energy transfer'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'P_diss_ross', NCFLOAT,4,dims,iret)
      name = 'Energy transfer'; unit = 'm^2/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_tdim,z_udim,iTimedim/)
      id  = ncvdef (ncid,'we', NCFLOAT,4,dims,iret)
      name = 'Eddy velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_tdim,lat_udim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'ue', NCFLOAT,4,dims,iret)
      name = 'Eddy velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      dims = (/Lon_udim,lat_tdim,z_tdim,iTimedim/)
      id  = ncvdef (ncid,'ve', NCFLOAT,4,dims,iret)
      name = 'Eddy velocity'; unit = 'm/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      call ncclos (ncid, iret)
 endif
   
 if (my_pe==0) iret=nf_open('pyOM_rossmix.cdf',NF_WRITE,ncid)

 if (my_pe==0) iret=nf_close(ncid)
 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine init_snap_cdf_rossmix


subroutine diag_snap_rossmix
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use rossmix_module   
 use diagnostics_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,npe,n
 real*8 :: a(1:ie_pe-is_pe+1,1:je_pe-js_pe+1,nphi),az(1:ie_pe-is_pe+1,1:je_pe-js_pe+1,nz),fxa
 integer :: itdimid,ilen,itimeid,id
 real*8, parameter :: spval = -1.0d33
 character :: name*32

 if (my_pe == 0 ) then
   iret=nf_open('pyOM_rossmix.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_rossmix.cdf'
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

 do npe = 0,n_pes
   if (my_pe == npe ) then
    print*,'PE ',my_pe,' writes to disk ...'
    iret=nf_open('pyOM_rossmix.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',itdimid)
    iret=nf_inq_dimlen(ncid, itdimid,ilen)

    do n=1,nmodes
      write(name, '("E_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=E_s(is_pe:ie_pe,js_pe:je_pe,:,tau,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("E_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=E_l(is_pe:ie_pe,js_pe:je_pe,:,tau,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("cgv_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=cgv_s(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskVp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("cgv_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=cgv_l(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskVp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("cgu_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=cgu_s(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskUp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("cgu_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=cgu_l(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskUp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("phidot_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=phidot_s(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskWp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("phidot_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=phidot_l(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskWp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("gamma_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=gamma_l(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("gamma_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=gamma_s(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("nu_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=nu_s(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("nu_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=nu_l(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)


      write(name, '("lambda_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=lambda_l(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("lambda_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      a=lambda_s(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) a = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),a)

      write(name, '("phin_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      az=phinW(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

      write(name, '("R_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      az(:,:,nz)=Rn(is_pe:ie_pe,js_pe:je_pe,n)
      where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,1/),az(:,:,nz))

      write(name, '("phi_ga_l_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      az=phi_ga_l(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

      write(name, '("phi_ga_s_",i1)') n-1; !call replace_space_zero(name)
      iret=nf_inq_varid(ncid,name,id)
      az=phi_ga_s(is_pe:ie_pe,js_pe:je_pe,:,n)
      where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
      iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

     enddo

    iret=nf_inq_varid(ncid,'Nsqr',id)
    az=Nsqr(is_pe:ie_pe,js_pe:je_pe,:,tau)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'Bx',id)
    az=Bx(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'By',id)
    az=By(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'stress_ux',id)
    az=stress_ux(is_pe:ie_pe,js_pe:je_pe,:)
    !where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'stress_vx',id)
    az=stress_vx(is_pe:ie_pe,js_pe:je_pe,:)
    !where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'stress_uy',id)
    az=stress_uy(is_pe:ie_pe,js_pe:je_pe,:)
    !where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'stress_vy',id)
    az=stress_vy(is_pe:ie_pe,js_pe:je_pe,:)
    !where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'K_diss_ross',id)
    az=K_diss_ross(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'P_diss_ross',id)
    az=P_diss_ross(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'ue',id)
    az=rossmix_ue(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskU(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'ve',id)
    az=rossmix_ve(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskV(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret=nf_inq_varid(ncid,'we',id)
    az=rossmix_we(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) az = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),az)

    iret = nf_close (ncid)
   endif
  call fortran_barrier
 enddo

end subroutine diag_snap_rossmix





