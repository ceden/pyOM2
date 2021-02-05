


integer function closest_angle(phi)
 use main_module   
 use rossmix_module   
 implicit none
 real*8, intent(in) :: phi
 real*8 :: pp
 pp = phi
 if (pp < 0.) pp = pp + 2*pi 
 if (pp > 2*pi) pp = pp - 2*pi
 !closest_angle = minloc( (phit(2:nphi-1)  - pp)**2,1 )+1 
 closest_angle = minloc(dabs(phit(2:nphi-1)-pp),1)+1 
end function

integer function zonal_bc_long(phi,Rloc,kloc)
 implicit none
 real*8, intent(in) :: phi,Rloc,kloc
 real*8 ::  sinp,cosp,rRn
 integer :: closest_angle
 rRn = 1./max(1d-12,Rloc)
 sinp = sin(phi)
 cosp = cos(phi)
 zonal_bc_long = closest_angle( datan( kloc**2*sinp*cosp / (kloc**2*sinp*sinp+rRn**2) )  )
end function

integer function zonal_bc_short(phi)
 implicit none
 real*8, intent(in) :: phi
 real*8 ::  sinp
 integer :: closest_angle
 sinp = sin(phi)
 zonal_bc_short = closest_angle(datan( cos(phi)/dsign(max(1d-22,dabs(sinp)),sinp)  ))
end function

integer function meridional_bc(phi)
 implicit none
 real*8, intent(in) :: phi
 integer :: closest_angle
 meridional_bc = closest_angle(datan(-dtan(phi)))
end function



subroutine reflect_rossmix_initialize
 use main_module   
 use rossmix_module   
 implicit none
 integer :: i,j,n

 ! western boundary
 n=0
 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
    if (maskT(i,j,nz) == 0.0 .and. maskT(i+1,j,nz)== 1.0) n=n+1
  enddo
 enddo
 allocate( boundary_west_i(n), boundary_west_j(n) )
 n=1
 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
    if (maskT(i,j,nz) == 0.0 .and. maskT(i+1,j,nz)== 1.0) then
      boundary_west_i(n)=i; boundary_west_j(n)=j; n=n+1
    endif
  enddo
 enddo

 
 ! eastern boundary
 n=0
 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
    if (maskT(i,j,nz) == 1.0 .and. maskT(i+1,j,nz)== 0.0) n=n+1
  enddo
 enddo
 allocate( boundary_east_i(n), boundary_east_j(n) )
 n=1
 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
    if (maskT(i,j,nz) == 1.0 .and. maskT(i+1,j,nz)== 0.0) then
      boundary_east_i(n)=i; boundary_east_j(n)=j; n=n+1
    endif
  enddo
 enddo

 ! southern boundary
 n=0
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
     if (maskT(i,j,nz) == 0.0 .and. maskT(i,j+1,nz)== 1.0) n=n+1
  enddo
 enddo
 allocate( boundary_south_i(n), boundary_south_j(n) )
 n=1
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
    if (maskT(i,j,nz) == 0.0 .and. maskT(i,j+1,nz)== 1.0) then
      boundary_south_i(n)=i; boundary_south_j(n)=j; n=n+1
    endif
  enddo
 enddo  

 ! northern boundary
 n=0
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
    if (maskT(i,j,nz) == 1.0 .and. maskT(i,j+1,nz)== 0.0) n=n+1
  enddo
 enddo
 allocate( boundary_north_i(n), boundary_north_j(n) )
 n=1
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
    if (maskT(i,j,nz) == 1.0 .and. maskT(i,j+1,nz)== 0.0) then
      boundary_north_i(n)=i; boundary_north_j(n)=j; n=n+1
    endif
  enddo
 enddo  
end subroutine reflect_rossmix_initialize




subroutine reflect_rossmix(is_,ie_,js_,je_,nphi_,m,advfe_s,advfn_s,advfe_l,advfn_l,El_loc,Es_loc)
 ! refection boundary condition for advective flux
 ! always upwind at boundary
 use main_module   
 use rossmix_module   
 implicit none
 integer, intent(in) :: is_,ie_,js_,je_,nphi_,m
 real*8, intent(inout) :: advfe_s(is_:ie_,js_:je_,nphi_),advfn_s(is_:ie_,js_:je_,nphi_)
 real*8, intent(inout) :: advfe_l(is_:ie_,js_:je_,nphi_),advfn_l(is_:ie_,js_:je_,nphi_)
 real*8, intent(in) :: El_loc(is_:ie_,js_:je_,nphi_)
 real*8, intent(in) :: Es_loc(is_:ie_,js_:je_,nphi_)
 integer :: i,j,k,kk,n,zonal_bc_long,zonal_bc_short,meridional_bc
 real*8 :: flux



 ! reflexion at western boundary
 do k=2,nphi-1
  do n=1,size(boundary_west_i,1)
     i=boundary_west_i(n); j=boundary_west_j(n)
     
     if (cgu_l(i+1,j,k,m) < 0) then
      flux = cgu_l(i+1,j,k,m)*el_loc(i+1,j,k)
      kk = zonal_bc_long(phit(k),Rn(i+1,j,m),km_l(i+1,j,m) )
      advfe_l(i,j,k)  = advfe_l(i,j,k) + flux
      advfe_s(i,j,kk) = advfe_s(i,j,kk) - flux
     endif

     !if (cgu_s(i+1,j,k,m) < 0) then
      ! kk = zonal_bc_short(phit(k))
     ! flux = cgu_s(i+1,j,k,m)*es_loc(i+1,j,k)
     ! advfe_s(i,j,k)  = advfe_s(i,j,k) + flux
     ! advfe_l(i,j,kk) = advfe_l(i,j,kk) - flux
     !endif
  enddo
 enddo

 ! reflexion at eastern boundary
 do k=2,nphi-1
  do n=1,size(boundary_east_i,1)
     i=boundary_east_i(n); j=boundary_east_j(n)

     if ( cgu_s(i-1,j,k,m) > 0) then
      flux = cgu_s(i-1,j,k,m)*es_loc(i,j,k)
      kk =zonal_bc_short(phit(k))
      advfe_s(i,j,k)  = advfe_s(i,j,k) + flux
      advfe_l(i,j,kk) = advfe_l(i,j,kk) - flux
     endif

     !if (cgu_l(i-1,j,k,m) > 0) then
     ! flux = cgu_l(i-1,j,k,m)*el_loc(i,j,k)
     ! kk = zonal_bc_long(phit(k),Rn(i,j,m),km_l(i,j,m) )
     ! advfe_l(i,j,k) = advfe_l(i,j,k) + flux
     ! advfe_s(i,j,kk) = advfe_s(i,j,kk) - flux
     !endif
  enddo
 enddo




 ! reflexion at southern boundary
 do k=2,nphi-1
  do n=1,size(boundary_south_i,1)
     i=boundary_south_i(n); j=boundary_south_j(n)
     kk = meridional_bc(phit(k))

     if (cgv_l(i,j+1,k,m) < 0) then
        flux = cosu(j)*cgv_l(i,j+1,k,m)*el_loc(i,j+1,k)
        advfn_l(i,j,k ) = advfn_l(i,j,k ) + flux
        advfn_l(i,j,kk) = advfn_l(i,j,kk) - flux
     endif

     if (cgv_s(i,j+1,k,m) < 0) then
        flux = cosu(j)*cgv_s(i,j+1,k,m)*es_loc(i,j+1,k)
        advfn_s(i,j,k ) = advfn_s(i,j,k ) + flux
        advfn_s(i,j,kk) = advfn_s(i,j,kk) - flux
     endif

  enddo
 enddo


 ! reflexion at northern boundary
 do k=2,nphi-1
  do n=1,size(boundary_north_i,1)
     i=boundary_north_i(n); j=boundary_north_j(n)
     kk = meridional_bc(phit(k))

     if (cgv_l(i,j-1,k,m) > 0) then
       flux = cosu(j)*cgv_l(i,j-1,k,m)*el_loc(i,j,k)
       advfn_l(i,j,k)  = advfn_l(i,j,k) + flux
       advfn_l(i,j,kk) = advfn_l(i,j,kk) - flux
     endif

     if (cgv_s(i,j-1,k,m) > 0) then
       flux = cosu(j)*cgv_s(i,j-1,k,m)*es_loc(i,j,k)
       advfn_s(i,j,k)  = advfn_s(i,j,k) + flux
       advfn_s(i,j,kk) = advfn_s(i,j,kk) - flux
     endif

  enddo
 enddo
end subroutine reflect_rossmix



