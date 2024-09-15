




subroutine rossmix2_init_snap_cdf
!=======================================================================
!     initialize NetCDF snapshot file for rossmix2 variables
!=======================================================================
 use main_module   
 use rossmix2_module   
 implicit none
 include "netcdf.inc"
 integer :: iret, ncid, id
 integer :: p_tdim,p_udim,p_tid,p_uid,xtdim,xudim,ytdim,yudim
 integer :: ztdim,zudim,tdim
 character :: name*32, unit*16
 real*8, parameter :: spval = -1.0d33
  
 if (my_pe==0) print'(a)',' preparing file pyOM_rossmix2.cdf'

 call def_grid_cdf('pyOM_rossmix2.cdf')

 if (my_pe==0) then
      iret=nf_open('pyOM_rossmix2.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      if (iret /= 0) call halt_stop('NETCDF:'//nf_strerror(iret))
!     dimensions
      p_tdim  = ncddef(ncid, 'phit', nphi , iret)
      p_udim  = ncddef(ncid, 'phiu', nphi , iret)
!     grid variables
      p_tid  = ncvdef (ncid,'phit',NCFLOAT,1,(/p_tdim/),iret)
      name = 'wave angle t grid'; unit = ' '
      call ncaptc(ncid, p_tid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, p_tid, 'units',     NCCHAR, 16, unit, iret) 
      p_uid  = ncvdef (ncid,'phiu',NCFLOAT,1,(/p_udim/),iret)
      name = 'wave angle u grid'; unit = ' '
      call ncaptc(ncid, p_uid, 'long_name', NCCHAR, 24, name, iret) 
      call ncaptc(ncid, p_uid, 'units',     NCCHAR, 16, unit, iret) 
!     attributes of the grid
      call ncendf(ncid, iret)
      iret= nf_put_vara_double(ncid,p_Tid,(/1/),(/nphi/) ,phit)
      iret= nf_put_vara_double(ncid,p_uid,(/1/),(/nphi/) ,phiu)
      call ncclos (ncid, iret)
  endif

 if (my_pe==0) then
      iret=nf_open('pyOM_rossmix2.cdf',NF_WRITE, ncid)
      iret=nf_set_fill(ncid, NF_NOFILL, iret)
      call ncredf(ncid, iret)
      iret=nf_inq_dimid(ncid,'xt',xtdim)
      iret=nf_inq_dimid(ncid,'xu',xudim)
      iret=nf_inq_dimid(ncid,'yt',ytdim)
      iret=nf_inq_dimid(ncid,'yu',yudim)
      iret=nf_inq_dimid(ncid,'phit',p_tdim)
      iret=nf_inq_dimid(ncid,'phiu',p_udim)
      iret=nf_inq_dimid(ncid,'zt',ztdim)
      iret=nf_inq_dimid(ncid,'zu',zudim)
      iret=nf_inq_dimid(ncid,'Time',tdim)
 
      id  = ncvdef (ncid,'f', NCFLOAT,3,(/xtdim,ytdim,tdim/),iret)
      name = 'Coriolis frequency'; unit = '1/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'R1', NCFLOAT,3,(/xtdim,ytdim,tdim/),iret)
      name = 'Rossby radius'; unit = 'm'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'beta', NCFLOAT,3,(/xtdim,ytdim,tdim/),iret)
      name = 'Change of Coriolis frequency'; unit = '1/ms'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'Nsqr', NCFLOAT,4,(/xtdim,ytdim,zudim,tdim/),iret)
      name = 'Square of stability frequency'; unit = '1/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'Bx', NCFLOAT,4,(/xudim,ytdim,zudim,tdim/),iret)
      name = 'Streamfunction'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'By', NCFLOAT,4,(/xtdim,yudim,zudim,tdim/),iret)
      name = 'Streamfunction'; unit = 'm^2/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

!      id  = ncvdef (ncid,'ue', NCFLOAT,4,(/xudim,ytdim,ztdim,tdim/),iret)
!      name = 'Eddy velocity'; unit = 'm/s'
!      call dvcdf(ncid,id,name,32,unit,16,spval)

!      id  = ncvdef (ncid,'ve', NCFLOAT,4,(/xtdim,yudim,ztdim,tdim/),iret)
!      name = 'Eddy velocity'; unit = 'm/s'
!      call dvcdf(ncid,id,name,32,unit,16,spval)

!      id  = ncvdef (ncid,'we', NCFLOAT,4,(/xtdim,ytdim,zudim,tdim/),iret)
!      name = 'Eddy velocity'; unit = 'm/s'
!      call dvcdf(ncid,id,name,32,unit,16,spval)

!      id  = ncvdef (ncid,'ug', NCFLOAT,4,(/xudim,ytdim,ztdim,tdim/),iret)
!      name = 'Background velocity'; unit = 'm/s'
!      call dvcdf(ncid,id,name,32,unit,16,spval)

!      id  = ncvdef (ncid,'vg', NCFLOAT,4,(/xtdim,ytdim,ztdim,tdim/),iret)
!      name = 'Background velocity'; unit = 'm/s'
!      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'ub_lim', NCFLOAT,4,(/xtdim,ytdim,ztdim,tdim/),iret)
      name = 'ub limiter'; unit = ' '
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'uu_lim', NCFLOAT,4,(/xtdim,ytdim,ztdim,tdim/),iret)
      name = 'uu limiter'; unit = ' '
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'E', NCFLOAT,4,(/xtdim,ytdim,p_tdim,tdim/),iret)
      name = 'Energy'; unit = 'm^3/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)
 
      id  = ncvdef (ncid,'forc', NCFLOAT,4,(/xtdim,ytdim,p_tdim,tdim/),iret)
      name = 'Energy transfer'; unit = 'm^3/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'forc_uu', NCFLOAT,4,(/xtdim,ytdim,p_tdim,tdim/),iret)
      name = 'Energy transfer'; unit = 'm^3/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'forc_diss', NCFLOAT,3,(/xtdim,ytdim,tdim/),iret)
      name = 'Energy transfer'; unit = 'm^3/s^3'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      if (enable_rossmix2_write_ub) then
       id  = ncvdef (ncid,'ub', NCFLOAT,5,(/xtdim,ytdim,ztdim,p_tdim,tdim/),iret)
       name = 'Structure function'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)

       id  = ncvdef (ncid,'uu', NCFLOAT,5,(/xtdim,ytdim,ztdim,p_tdim,tdim/),iret)
       name = 'Structure function'; unit = ' '
       call dvcdf(ncid,id,name,32,unit,16,spval)
      
      endif

      id  = ncvdef (ncid,'om_r', NCFLOAT,4,(/xtdim,ytdim,p_tdim,tdim/),iret)
      name = 'Frequency'; unit = '1/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'om_i', NCFLOAT,4,(/xtdim,ytdim,p_tdim,tdim/),iret)
      name = 'Frequency'; unit = '1/s'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'kh_appr', NCFLOAT,3,(/xtdim,ytdim,tdim/),iret)
      name = 'Wavenumber'; unit = '1/m'
      call dvcdf(ncid,id,name,32,unit,16,spval)


      id  = ncvdef (ncid,'stress_ux', NCFLOAT,4,(/xtdim,ytdim,ztdim,tdim/),iret)
      name = 'lateral stress'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'stress_uy', NCFLOAT,4,(/xtdim,ytdim,ztdim,tdim/),iret)
      name = 'lateral stress'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'stress_vx', NCFLOAT,4,(/xtdim,ytdim,ztdim,tdim/),iret)
      name = 'lateral stress'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)

      id  = ncvdef (ncid,'stress_vy', NCFLOAT,4,(/xtdim,ytdim,ztdim,tdim/),iret)
      name = 'lateral stress'; unit = 'm^2/s^2'
      call dvcdf(ncid,id,name,32,unit,16,spval)



      call ncclos (ncid, iret)
 endif
   
 if (my_pe==0) iret=nf_open('pyOM_rossmix2.cdf',NF_WRITE,ncid)

 if (my_pe==0) iret=nf_close(ncid)
 if (my_pe==0) print'(a)',' done '
 call fortran_barrier

end subroutine rossmix2_init_snap_cdf






subroutine rossmix2_diag_snap
!=======================================================================
!     write to NetCDF snapshot file
!=======================================================================
 use main_module   
 use rossmix2_module   
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,npe,n
 real*8 :: fxa,aloc(1:ie_pe-is_pe+1,1:je_pe-js_pe+1,nz)
 real*8 :: bloc(1:ie_pe-is_pe+1,1:je_pe-js_pe+1,nphi)
 integer :: tdim,ilen,id
 real*8, parameter :: spval = -1.0d33


 if (my_pe == 0 ) then
   iret=nf_open('pyOM_rossmix2.cdf',NF_WRITE,ncid)
   print'(a)',' writing to file pyOM_rossmix2.cdf'
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',tdim)
   iret=nf_inq_dimlen(ncid,tdim,ilen)
   ilen=ilen+1
   fxa = itt*dt_tracer/86400.0
   iret=nf_inq_varid(ncid,'Time',id)
   iret= nf_put_vara_double(ncid,id,(/ilen/),(/1/),(/fxa/))
   call ncclos (ncid, iret)
 endif
 call fortran_barrier

 do npe = 0,n_pes
   if (my_pe == npe ) then
    iret=nf_open('pyOM_rossmix2.cdf',NF_WRITE,ncid)
    iret=nf_inq_dimid(ncid,'Time',tdim)
    iret=nf_inq_dimlen(ncid, tdim,ilen)

    iret=nf_inq_varid(ncid,'f',id)
    aloc(:,:,1) = coriolis_t(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) aloc(:,:,1) = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,1/),aloc(:,:,1))

    iret=nf_inq_varid(ncid,'R1',id)
    aloc(:,:,1) = R1(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) aloc(:,:,1) = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,1/),aloc(:,:,1))

    iret=nf_inq_varid(ncid,'beta',id)
    aloc(:,:,1)=beta(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) aloc(:,:,1) = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,1/),aloc(:,:,1))
    
    iret=nf_inq_varid(ncid,'Nsqr',id)
    aloc=Nsqr_lim(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

    iret=nf_inq_varid(ncid,'Bx',id)
    aloc=Bx(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)
    
    iret=nf_inq_varid(ncid,'By',id)
    aloc=By(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

!    iret=nf_inq_varid(ncid,'ue',id)
!    aloc=ue(is_pe:ie_pe,js_pe:je_pe,:)
!    where( maskU(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
!    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

!    iret=nf_inq_varid(ncid,'ve',id)
!    aloc=ve(is_pe:ie_pe,js_pe:je_pe,:)
!    where( maskV(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
!    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

!    iret=nf_inq_varid(ncid,'we',id)
!    aloc=we(is_pe:ie_pe,js_pe:je_pe,:)
!    where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
!    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

!    iret=nf_inq_varid(ncid,'ug',id)
!    aloc=ug(is_pe:ie_pe,js_pe:je_pe,:)
!    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
!    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

!    iret=nf_inq_varid(ncid,'vg',id)
!    aloc=vg(is_pe:ie_pe,js_pe:je_pe,:)
!    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
!    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

    iret=nf_inq_varid(ncid,'ub_lim',id)
    aloc=ub_lim(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

    iret=nf_inq_varid(ncid,'uu_lim',id)
    aloc=uu_lim(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

    iret=nf_inq_varid(ncid,'E',id)   
    bloc=E_r(is_pe:ie_pe,js_pe:je_pe,:,tau)
    where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) bloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),bloc)

    iret=nf_inq_varid(ncid,'forc',id)   
    bloc= forcx(is_pe:ie_pe,js_pe:je_pe,:) + forcy(is_pe:ie_pe,js_pe:je_pe,:)  
    where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) bloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),bloc)

    iret=nf_inq_varid(ncid,'forc_uu',id)   
    bloc= forcuu(is_pe:ie_pe,js_pe:je_pe,:) 
    where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) bloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),bloc)

    iret=nf_inq_varid(ncid,'forc_diss',id)
    aloc(:,:,1) = forc_diss(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) aloc(:,:,1) = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,1/),aloc(:,:,1))


  if (enable_rossmix2_write_ub) then
    do n=1,nphi
     iret=nf_inq_varid(ncid,'ub',id)   
     aloc=ub(is_pe:ie_pe,js_pe:je_pe,:,n)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
     iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,n,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1,1/),aloc)

     iret=nf_inq_varid(ncid,'uu',id)   
     aloc=uu(is_pe:ie_pe,js_pe:je_pe,:,n)
     where( maskW(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
     iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,n,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1,1/),aloc)
    enddo 
  endif

    iret=nf_inq_varid(ncid,'om_r',id)   
    bloc=real(om_max(is_pe:ie_pe,js_pe:je_pe,:))
    where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) bloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),bloc)
 
    iret=nf_inq_varid(ncid,'om_i',id)   
    bloc=imag(om_max(is_pe:ie_pe,js_pe:je_pe,:))
    where( maskTp(is_pe:ie_pe,js_pe:je_pe,:) == 0.) bloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nphi,1/),bloc)

    iret=nf_inq_varid(ncid,'kh_appr',id)
    aloc(:,:,1)=kh_appr(is_pe:ie_pe,js_pe:je_pe)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,nz) == 0.) aloc(:,:,1) = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,1/),aloc(:,:,1))


    iret=nf_inq_varid(ncid,'stress_ux',id)   
    aloc = stress_ux(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)
 
    iret=nf_inq_varid(ncid,'stress_uy',id)   
    aloc = stress_uy(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

    iret=nf_inq_varid(ncid,'stress_vx',id)   
    aloc = stress_vx(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)
 
    iret=nf_inq_varid(ncid,'stress_vy',id)   
    aloc = stress_vy(is_pe:ie_pe,js_pe:je_pe,:)
    where( maskT(is_pe:ie_pe,js_pe:je_pe,:) == 0.) aloc = spval
    iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe,1,ilen/), (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/),aloc)

    iret = nf_close (ncid)
   endif
  call fortran_barrier
 enddo

end subroutine rossmix2_diag_snap


