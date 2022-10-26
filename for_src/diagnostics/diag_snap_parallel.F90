
#ifdef with_netcdf_parallel 



subroutine init_diag_snap_parallel
 use main_module
 implicit none
 
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,id,xtdim,ytdim,ztdim,tdim,nc_mode,start(4),count(4),n
 integer :: xudim,yudim,zudim
 character :: fname*80,name*32
 real*8, parameter :: spval = -1.0d33
  
 fname = 'pyOM.cdf' 
 if (my_pe==0)  print*,'preparing file ',fname(1:len_trim(fname))

 nc_mode = IOR(nf_clobber,nf_mpiio)
 nc_mode = IOR(nc_mode,nf_netcdf4 )
 nc_mode = IOR(nc_mode,nf_classic_model )
 
 iret = NF_CREATE_par(fname,nc_mode,MPI_COMM_WORLD,MPI_INFO_NULL, ncid) 
 if (my_pe==0.and.iret/=0) print*,'at opening:',nf_strerror(iret)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)
  
 iret = nf_def_dim(ncid,'xt', nx,xtdim)
 iret = nf_def_dim(ncid,'yt', ny,ytdim)
 iret = nf_def_dim(ncid,'zt', nz,ztdim) 
 
 iret = nf_def_dim(ncid,'xu', nx,xudim)
 iret = nf_def_dim(ncid,'yu', ny,yudim)
 iret = nf_def_dim(ncid,'zu', nz,zudim) 
 
 iret = nf_def_dim(ncid,'time', nf_unlimited,tdim)
 if (my_pe==0.and.iret/=0) print*,nf_strerror(iret)
  
 iret = nf_def_var(ncid, 'xt',NF_DOUBLE,1,xtdim,id) 
 if (coord_degree) then
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('Longitude'),'Longitude')
  iret = nf_put_att_text(ncid,id,'units',len_trim('degrees E'),'degrees E')
 else
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('xt'),'xt')
  iret = nf_put_att_text(ncid,id,'units',len_trim('km'),'km')
 endif
 
 iret = nf_def_var(ncid, 'xu',NF_DOUBLE,1,xudim,id) 
 if (coord_degree) then
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('Longitude'),'Longitude')
  iret = nf_put_att_text(ncid,id,'units',len_trim('degrees E'),'degrees E') 
 else
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('xu'),'xu')
  iret = nf_put_att_text(ncid,id,'units',len_trim('km'),'km')
 endif
 
 iret = nf_def_var(ncid, 'yt',NF_DOUBLE,1,ytdim,id) 
 if (coord_degree) then
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('Latitude'),'Latitude')
  iret = nf_put_att_text(ncid,id,'units',len_trim('degrees N'),'degrees N')
 else
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('yt'),'yt')
  iret = nf_put_att_text(ncid,id,'units',len_trim('km'),'km') 
 endif 
 
 iret = nf_def_var(ncid, 'yu',NF_DOUBLE,1,yudim,id) 
 if (coord_degree) then
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('Latitude'),'Latitude')
  iret = nf_put_att_text(ncid,id,'units',len_trim('degrees N'),'degrees N')
 else
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('yu'),'yu')
  iret = nf_put_att_text(ncid,id,'units',len_trim('km'),'km')
 endif 
 
 iret = nf_def_var(ncid, 'zt',NF_DOUBLE,1,ztdim,id) 
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('zt'),'zt')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m') 

 iret = nf_def_var(ncid, 'zu',NF_DOUBLE,1,zudim,id) 
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('zu'),'zu')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m') 
  
 iret = nf_def_var(ncid, 'time',NF_DOUBLE,1,tdim,id)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('Time'),'Time')
 iret = nf_put_att_text(ncid,id,'units',len_trim('days'),'days') 
 iret = nf_put_att_text(ncid,id,'time_origin',len_trim('01-JAN-1900 00:00:00'),'01-JAN-1900 00:00:00')   
   
 id  = ncvdef (ncid,'ht',NF_DOUBLE,2,(/xtdim, ytdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('depth'),'depth') 
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m')
  
 id  = ncvdef (ncid,'u',NF_DOUBLE,4,(/xudim, ytdim,ztdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity') 
 iret = nf_put_att_text(ncid,id,'units',len_trim('m/s'),'m/s')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 id  = ncvdef (ncid,'v',NF_DOUBLE,4,(/xtdim, yudim,ztdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m/s'),'m/s')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
  
 id  = ncvdef (ncid,'w',NF_DOUBLE,4,(/xtdim, ytdim,zudim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m/s'),'m/s')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 id  = ncvdef (ncid,'temp',NF_DOUBLE,4,(/xtdim, ytdim,ztdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('Temperature'),'Temperature')
 iret = nf_put_att_text(ncid,id,'units',len_trim('deg C'),'deg C')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 id  = ncvdef (ncid,'salt',NF_DOUBLE,4,(/xtdim, ytdim,ztdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('Salinity'),'Salinity')
 iret = nf_put_att_text(ncid,id,'units',len_trim('g/kg'),'g/kg')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 if (enable_streamfunction) then 
   id  = ncvdef (ncid,'psi',NF_DOUBLE,3,(/xudim, yudim,tdim/),iret)
   iret = nf_put_att_text(ncid,id,'long_name',len_trim('Streamfunction'),'Streamfunction')
   iret = nf_put_att_text(ncid,id,'units',len_trim('m^3/s'),'m^3/s') 
   iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
   do n=1,nisle
        write(name, '("psi",i3)') n
        call replace_space_zero(name)
        id  = ncvdef (ncid,name, NF_DOUBLE,2,(/xudim, yudim/),iret)
        write(name, '("Boundary streamfunction ",i3)') n       
        iret = nf_put_att_text(ncid,id,'long_name',len_trim(name),name)        
        iret = nf_put_att_text(ncid,id,'units',len_trim('m^3/s'),'m^3/s')
        iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
   enddo
 else
   id  = ncvdef (ncid,'surf_press',NF_DOUBLE,3,(/xtdim, ytdim,tdim/),iret)
   iret = nf_put_att_text(ncid,id,'long_name',len_trim('Surface pressure'),'Surface pressure')
   iret = nf_put_att_text(ncid,id,'units',len_trim('m^2/s^2'),'m^2/s^2')
   iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
   id  = ncvdef (ncid,'press',NF_DOUBLE,4,(/xtdim, ytdim,ztdim,tdim/),iret)
   iret = nf_put_att_text(ncid,id,'long_name',len_trim('total pressure'),'total pressure')
   iret = nf_put_att_text(ncid,id,'units',len_trim('m^2/s^2'),'m^2/s^2')
   iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 endif

 if (.not. enable_hydrostatic) then
   id  = ncvdef (ncid,'p_hydro',NF_DOUBLE,4,(/xtdim, ytdim,ztdim,tdim/),iret)
   iret = nf_put_att_text(ncid,id,'long_name',len_trim('Hydrostatic pressure'),'Hydrostatic pressure')
   iret = nf_put_att_text(ncid,id,'units',len_trim('m^2/s^2'),'m^2/s^2')
   iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
   id  = ncvdef (ncid,'p_non_hydro',NF_DOUBLE,4,(/xtdim, ytdim,ztdim,tdim/),iret)
   iret = nf_put_att_text(ncid,id,'long_name',len_trim('Non hydrostatic pressure'),'Non hydrostatic pressure')
   iret = nf_put_att_text(ncid,id,'units',len_trim('m^2/s^2'),'m^2/s^2')
   iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 endif

 id  = ncvdef (ncid,'forc_temp_surface',NF_DOUBLE,3,(/xtdim, ytdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('Surface temperature flux'),'Surface temperature flux')
 iret = nf_put_att_text(ncid,id,'units',len_trim('deg C m/s'),'deg C m/s')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 id  = ncvdef (ncid,'forc_salt_surface',NF_DOUBLE,3,(/xtdim, ytdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('Surface salinity flux'),'Surface salinity flux')
 iret = nf_put_att_text(ncid,id,'units',len_trim('g/Kg m/s'),'g/Kg m/s')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 id  = ncvdef (ncid,'taux',NF_DOUBLE,3,(/xudim, ytdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('Surface wind stress'),'Surface wind stress')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m^2/s^2'),'m^2/s^2')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 id  = ncvdef (ncid,'tauy',NF_DOUBLE,3,(/xtdim, yudim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('Surface wind stress'),'Surface wind stress')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m^2/s^2'),'m^2/s^2')
 iret = nf_put_att_double(ncid,id,'missing_value',nf_double,1,spval)
 
 iret= nf_enddef(ncid)
                
 iret=nf_inq_varid(ncid,'xt',id)
 if (coord_degree) then
  iret= nf_put_vara_double(ncid,id,is_pe,ie_pe-is_pe+1,xt(is_pe:ie_pe))  
 else
  iret= nf_put_vara_double(ncid,id,is_pe,ie_pe-is_pe+1,xt(is_pe:ie_pe)/1e3) 
 endif
 iret=nf_inq_varid(ncid,'xu',id)
 if (coord_degree) then
  iret= nf_put_vara_double(ncid,id,is_pe,ie_pe-is_pe+1,xu(is_pe:ie_pe))  
 else
  iret= nf_put_vara_double(ncid,id,is_pe,ie_pe-is_pe+1,xu(is_pe:ie_pe)/1e3) 
 endif
 
 iret=nf_inq_varid(ncid,'yt',id)
 if (coord_degree) then
  iret= nf_put_vara_double(ncid,id,js_pe,je_pe-js_pe+1,yt(js_pe:je_pe)) 
 else
  iret= nf_put_vara_double(ncid,id,js_pe,je_pe-js_pe+1,yt(js_pe:je_pe)/1e3) 
 endif
 
 iret=nf_inq_varid(ncid,'yu',id)
 if (coord_degree) then
  iret= nf_put_vara_double(ncid,id,js_pe,je_pe-js_pe+1,yu(js_pe:je_pe))   
 else
  iret= nf_put_vara_double(ncid,id,js_pe,je_pe-js_pe+1,yu(js_pe:je_pe)/1e3)  
 endif 

 iret=nf_inq_varid(ncid,'zt',id)
 iret= nf_put_vara_double(ncid,id,1,nz,zt(1:nz))   
 
 iret=nf_inq_varid(ncid,'zu',id)
 iret= nf_put_vara_double(ncid,id,1,nz,zw(1:nz))   
 
 iret=nf_inq_varid(ncid,'ht',id)
 iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe/),(/ie_pe-is_pe+1,je_pe-js_pe+1/),ht(is_pe:ie_pe,js_pe:je_pe))

 if (enable_streamfunction) then
  do n=1,nisle
   write(name, '("psi",i3)') n
   call replace_space_zero(name)
   iret=nf_inq_varid(ncid,name,id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,js_pe/),(/ie_pe-is_pe+1,je_pe-js_pe+1/),psin(is_pe:ie_pe,js_pe:je_pe,n))
  enddo
 endif

 iret= nf_close(ncid) 
 
 end subroutine init_diag_snap_parallel
 
 
 

subroutine diag_snap_parallel
 use main_module
 implicit none 
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,id,xtdim,ytdim,ztdim,tdim,nc_mode,start(4),count(4),ilen,k
 integer :: xudim,yudim,zudim
 real*8 :: aloc(is_pe:ie_pe,js_pe:je_pe,1:nz)
 real*8, parameter :: spval = -1.0d33
 character :: fname*80

 fname = 'pyOM.cdf' 
 if (my_pe==0)  print*,'writing to file ',fname(1:len_trim(fname)),' at itt=',itt

 nc_mode = IOR(nf_write,nf_mpiio)
 nc_mode = IOR(nc_mode,nf_netcdf4 )
 nc_mode = IOR(nc_mode,nf_classic_model )
 
 iret = NF_open_par(fname,nc_mode,MPI_COMM_WORLD,MPI_INFO_NULL, ncid) 
 if (my_pe==0.and.iret/=0) print*,'at opening:',nf_strerror(iret)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)

 iret=nf_inq_dimid(ncid,'time',id)
 iret=nf_inq_dimlen(ncid,id,ilen)
 ilen=ilen+1
 
 iret=nf_inq_varid(ncid,'time',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 iret= nf_put_vara_double(ncid,id,ilen,1,itt*dt_tracer/86400.0) 
 
 start = (/is_pe,        js_pe,        1,ilen/)
 count = (/ie_pe-is_pe+1,je_pe-js_pe+1,nz,1/)
 
 iret=nf_inq_varid(ncid,'u',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = u(is_pe:ie_pe,js_pe:je_pe,1:nz,tau)
 where (maskU(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
   
 iret=nf_inq_varid(ncid,'v',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = v(is_pe:ie_pe,js_pe:je_pe,1:nz,tau)
 where (maskV(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
                                  
 iret=nf_inq_varid(ncid,'w',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = w(is_pe:ie_pe,js_pe:je_pe,1:nz,tau)
 where (maskW(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
 iret= nf_put_vara_double(ncid,id,start,count,aloc)                                
 
 iret=nf_inq_varid(ncid,'temp',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = temp(is_pe:ie_pe,js_pe:je_pe,1:nz,tau)
 where (maskT(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
                                  
 iret=nf_inq_varid(ncid,'salt',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = salt(is_pe:ie_pe,js_pe:je_pe,1:nz,tau)
 where (maskT(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 if (.not. enable_hydrostatic) then
   iret=nf_inq_varid(ncid,'p_hydro',id)
   iret = nf_var_par_access(ncid, id, nf_collective)
   aloc = p_hydro(is_pe:ie_pe,js_pe:je_pe,1:nz)
   where (maskT(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
   iret= nf_put_vara_double(ncid,id,start,count,aloc)

   iret=nf_inq_varid(ncid,'p_non_hydro',id)
   iret = nf_var_par_access(ncid, id, nf_collective)
   aloc = p_non_hydro(is_pe:ie_pe,js_pe:je_pe,1:nz,taup1)
   where (maskT(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
   iret= nf_put_vara_double(ncid,id,start,count,aloc)
 endif

 if (.not.enable_streamfunction) then
   iret=nf_inq_varid(ncid,'press',id)
   iret = nf_var_par_access(ncid, id, nf_collective)
   aloc = p_hydro(is_pe:ie_pe,js_pe:je_pe,1:nz)
   do k=1,nz
    aloc(:,:,k) = aloc(:,:,k) + psi(is_pe:ie_pe,js_pe:je_pe,tau)
   enddo 
   where (maskT(is_pe:ie_pe,js_pe:je_pe,1:nz) == 0d0) aloc = spval
   iret= nf_put_vara_double(ncid,id,start,count,aloc)
 endif 
 
 start = (/is_pe,        js_pe,        ilen,1/)
 count = (/ie_pe-is_pe+1,je_pe-js_pe+1,1   ,1/)

 if (enable_streamfunction) then
    iret=nf_inq_varid(ncid,'psi',id)
 else    
    iret=nf_inq_varid(ncid,'surf_press',id)
 endif 
 aloc(:,:,1) = psi(is_pe:ie_pe,js_pe:je_pe,tau)
 where (maskT(is_pe:ie_pe,js_pe:je_pe,1) == 0d0) aloc(:,:,1) = spval
 iret = nf_var_par_access(ncid, id, nf_collective)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'forc_temp_surface',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc(:,:,1) = forc_temp_surface(is_pe:ie_pe,js_pe:je_pe)
 where (maskT(is_pe:ie_pe,js_pe:je_pe,1) == 0d0) aloc(:,:,1) = spval
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'forc_salt_surface',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc(:,:,1) = forc_salt_surface(is_pe:ie_pe,js_pe:je_pe)
 where (maskT(is_pe:ie_pe,js_pe:je_pe,1) == 0d0) aloc(:,:,1) = spval
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'taux',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc(:,:,1) = surface_taux(is_pe:ie_pe,js_pe:je_pe)
 where (maskU(is_pe:ie_pe,js_pe:je_pe,1) == 0d0) aloc(:,:,1) = spval
 iret= nf_put_vara_double(ncid,id,start,count,surface_taux(is_pe:ie_pe,js_pe:je_pe))

 iret=nf_inq_varid(ncid,'tauy',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc(:,:,1) = surface_tauy(is_pe:ie_pe,js_pe:je_pe)
 where (maskV(is_pe:ie_pe,js_pe:je_pe,1) == 0d0) aloc(:,:,1) = spval
 iret= nf_put_vara_double(ncid,id,start,count,surface_tauy(is_pe:ie_pe,js_pe:je_pe))
       
 iret= nf_close(ncid) 

 end subroutine diag_snap_parallel
 
 

 
#else

 subroutine init_diag_snap_parallel
 use main_module
 implicit none 
 if (my_pe==0) print*,'WARNING: no parallel netcdf available, no output possible'
 end subroutine init_diag_snap_parallel


 subroutine diag_snap_parallel
 use main_module
 implicit none 
 if (my_pe==0) print*,'WARNING: no parallel netcdf available, no output possible'
 end subroutine diag_snap_parallel
 
 
#endif



