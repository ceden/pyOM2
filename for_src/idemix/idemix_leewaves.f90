






subroutine idemix_friction_lee
!=======================================================================
!  wave-induced gravity wave drag by lee waves
!=======================================================================
 use main_module   
 use idemix_module   
 implicit none
 integer :: i,j,k,ks
 real*8 :: diss(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: flux(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,0:nz),yN,yS
 !---------------------------------------------------------------------------------
 ! explicit vertical friction of zonal momentum
 !---------------------------------------------------------------------------------
 
 do k=1,nz-1
     flux(:,:,k)= -(2/pi)*Lambda_0(:,:,k)*u_bot(:,:)/max(u0_lee_min,u0_bot(:,:))* &
                       (E_lee_p(:,:,k,tau)-E_lee_m(:,:,k,tau))*maskW(:,:,k)
 enddo
 flux(:,:,0)=0.
 flux(:,:,nz) = 0.0 
 do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks = max(1,kbot(i,j))
     flux(i,j,ks-1) = -stress_lee_u(i,j)
   enddo
 enddo 
 do i=is_pe-onx,ie_pe+onx-1
   flux(i,:,:) = 0.5*(flux(i,:,:)+flux(i+1,:,:))
 enddo 

 if (enable_leewaves_obc_taper_north) then
    yN = 0.0; if (my_blk_j == n_pes_j) yN = yu(ny)
    call global_max(yN)
    do j=js_pe-onx,je_pe+onx
      flux(:,j,:)=  flux(:,j,:)*(1- exp( -(yN-yt(j))**2/leewaves_north_taper_len**2  ) )
    enddo
 endif
 
 if (enable_leewaves_obc_taper_south) then
    yS = 1e20; if (my_blk_j == 1) yS = yu(0)
    call global_min(yS)
    do j=js_pe-onx,je_pe+onx
      flux(:,j,:)=  flux(:,j,:)*(1- exp( -(yt(j)-yS)**2/leewaves_south_taper_len**2  ) )
    enddo
 endif

 do k=1,nz
   du_mix(:,:,k) = du_mix(:,:,k) -(flux(:,:,k)-flux(:,:,k-1))/dzt(k)*maskU(:,:,k)
   diss(:,:,k) = u(:,:,k,tau)*(flux(:,:,k)-flux(:,:,k-1))/dzt(k)*maskU(:,:,k)
 enddo
 call ugrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
 K_diss_idemix = diss

 !---------------------------------------------------------------------------------
 ! explicit vertical friction of meridional momentum
 !---------------------------------------------------------------------------------
 do k=1,nz-1
     flux(:,:,k)= -(2/pi)*Lambda_0(:,:,k)*v_bot(:,:)/max(u0_lee_min,u0_bot(:,:))* &
                            (E_lee_p(:,:,k,tau)-E_lee_m(:,:,k,tau))*maskW(:,:,k)
 enddo
 flux(:,:,0)=0.
 flux(:,:,nz) = 0.0 
 do j=js_pe-onx,je_pe+onx
   do i=is_pe-onx,ie_pe+onx
     ks = max(1,kbot(i,j))
     flux(i,j,ks-1) = -stress_lee_v(i,j)
   enddo
 enddo
 do j=js_pe-onx,je_pe+onx-1
   flux(:,j,:) = 0.5*(flux(:,j,:)+flux(:,j+1,:))
 enddo

 if (enable_leewaves_obc_taper_north) then
    yN = 0.0; if (my_blk_j == n_pes_j) yN = yu(ny)
    call global_max(yN)
    do j=js_pe-onx,je_pe+onx
      flux(:,j,:)=  flux(:,j,:)*(1- exp( -(yN-yt(j))**2/leewaves_north_taper_len**2  ) )
    enddo
 endif
 
 if (enable_leewaves_obc_taper_south) then
    yS = 1e20; if (my_blk_j == 1) yS = yu(0)
    call global_min(yS)
    do j=js_pe-onx,je_pe+onx
      flux(:,j,:)=  flux(:,j,:)*(1- exp( -(yt(j)-yS)**2/leewaves_south_taper_len**2  ) )
    enddo
 endif

 
 do k=1,nz
   dv_mix(:,:,k) = dv_mix(:,:,k) - (flux(:,:,k)-flux(:,:,k-1))/dzt(k)*maskV(:,:,k)
   diss(:,:,k) = v(:,:,k,tau)*(flux(:,:,k)-flux(:,:,k-1))/dzt(k)*maskV(:,:,k)
 enddo
 call vgrid_to_tgrid(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,diss,diss)
 K_diss_idemix = K_diss_idemix + diss

end subroutine idemix_friction_lee





subroutine integrate_leewaves
 use main_module
 use idemix_module
 implicit none
 integer :: i,j,k,ks,n
 real*8 :: flux_p(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,0:nz)
 real*8 :: flux_m(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,0:nz)
 real*8 :: aloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 !real*8 :: bloc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: floc,Nloc,fxa,fxb,yN,yS,dt_loc
 real*8 :: uz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: vz(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 
 real*8 :: maskUtr,maskVtr
 maskUtr(i,j,k) = maskW(i+1,j,k)*maskW(i,j,k)
 maskVtr(i,j,k) = maskW(i,j+1,k)*maskW(i,j,k)
 logical, save :: first = .true.


 if (first) then
   !  Ah_lee *dt/dx^2 <= 0.5  
   
   if (enable_idemix_hor_diffusion_iter) then
    fxa = 0.
    do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx
       fxa = max(fxa, Ah_lee* (dt_tracer/idemix_hor_diffusion_iter)/dyt(j)**2 )
      enddo
    enddo
    call global_max(fxa)
    if (fxa>0.5) then 
     if (my_pe==0) print*,'ERROR: Ah_lee too large for grid'
     call halt_stop(' in integrate_leewaves (1)')
    endif
   endif
   
   if (enable_idemix_hor_diffusion) then
    fxa = 0.
    do j=js_pe-onx,je_pe+onx
      do i=is_pe-onx,ie_pe+onx
       fxa = max(fxa, Ah_lee* dt_tracer/dyt(j)**2 )
      enddo
    enddo
    call global_max(fxa)
    if (fxa>0.5) then 
     if (my_pe==0) print*,'ERROR: Ah_lee too large for grid'
     call halt_stop(' in integrate_leewaves (2)')
    endif
   endif

   first = .false.
  endif  

 
 do k=1,nz
   do j=js_pe-onx,je_pe+onx
     do i=is_pe-onx,ie_pe+onx
        floc = max(abs(coriolis_t(i,j)),1D-6)
        Nloc = max(floc*1.2, sqrt( max(0d0, Nsqr(i,j,k,tau))  ) )
        Lambda_0(i,j,k) = 0.65*floc/Nloc/log(10d0)
        
        !fxa = 0.5*floc/Nloc*(Nloc**2+0.5*floc**2)/(Nloc**2-floc**2 ) 
        !fxa = fxa * log( (Nloc+sqrt(Nloc**2-floc**2) )/(Nloc-sqrt(Nloc**2-floc**2)) )
        !fxa = fxa - 3./2.*floc/sqrt(Nloc**2-floc**2 ) 
        !Lambda_0(i,j,k) = fxa*(2./pi)/(1-2./pi*asin(floc/Nloc) )*maskW(i,j,k)
        
        fxa = 0.5*(1.+tanh( (Nloc/floc-Noverf_min0)/Noverf_min1)  )
        fxb = 0.5*(1.+tanh(  (abs(coriolis_t(i,j) ) -1d-6)/1d-8 )     )
        Lambda_0(i,j,k) = fxb*fxa*Lambda_0(i,j,k)*maskW(i,j,k)
    enddo
   enddo
  enddo

 
 do k=1,nz-1
  uz = (u(:,:,k+1,tau)-u(:,:,k,tau))/dzw(k)*maskU(:,:,k)*maskU(:,:,k+1)
  vz = (v(:,:,k+1,tau)-v(:,:,k,tau))/dzw(k)*maskV(:,:,k)*maskV(:,:,k+1)
  do j = js_pe,je_pe
   do i = is_pe,ie_pe
     fxa = (u_bot(i,j)*0.5*(uz(i,j)+uz(i-1,j)) &
           +v_bot(i,j)*0.5*(vz(i,j)+vz(i,j-1)) )/max(u0_bot(i,j),u0_lee_min)
     fxb =  (2/pi)*Lambda_0(i,j,k)*fxa  
     if (abs(fxb) > Tc_max) then
       Lambda_0(i,j,k) =  Tc_max/abs(fxa)*pi/2. 
       fxb =  (2/pi)*Lambda_0(i,j,k)*fxa 
     endif
     Tc_lee(i,j,k) = fxb 
   enddo
  enddo 
 enddo


 if (enable_leewaves_obc_taper_north) then
    yN = 0.0; if (my_blk_j == n_pes_j) yN = yu(ny)
    call global_max(yN)
    do j=js_pe-onx,je_pe+onx
      Tc_lee(:,j,:)=  Tc_lee(:,j,:)*(1- exp( -(yN-yt(j))**2/leewaves_north_taper_len**2  ) )
    enddo
 endif
 
 if (enable_leewaves_obc_taper_south) then
    yS = 1e20; if (my_blk_j == 1) yS = yu(0)
    call global_min(yS)
    do j=js_pe-onx,je_pe+onx
      Tc_lee(:,j,:)=  Tc_lee(:,j,:)*(1- exp( -(yt(j)-yS)**2/leewaves_south_taper_len**2  ) )
    enddo
 endif


 
 do k=1,nz-1
  c_lee(:,:,k)  = (2/pi)*0.5*(Lambda_0(:,:,k)+Lambda_0(:,:,k+1))*u0_bot(:,:)
 enddo

 flux_m=0;flux_p=0;mean_flow_to_E_lee=0.;iw_diss_lee=0.
 
 E_lee_p(:,:,:,taup1) = E_lee_p(:,:,:,tau) 
 E_lee_m(:,:,:,taup1) = E_lee_m(:,:,:,tau) 
 
 
 dt_loc = dt_tracer/idemix_lee_ts_fac
 do n=1,idemix_lee_ts_fac
 
  call leewave_superbee(flux_p, c_lee,E_lee_p(:,:,:,taup1))
  call leewave_superbee(flux_m,-c_lee,E_lee_m(:,:,:,taup1))

  do j = js_pe,je_pe
    do i = is_pe,ie_pe   
      ks = max(1,kbot(i,j))
      flux_m(i,j,ks-1) = - c_lee(i,j,ks)*E_lee_m(i,j,ks,taup1)
      flux_p(i,j,ks-1) = flux_lee_tot(i,j) - flux_m(i,j,ks-1)
     enddo
  enddo
  k=nz-1; flux_p(:,:,k) =  c_lee(:,:,k)*E_lee_p(:,:,k,taup1)
  k=nz-1; flux_m(:,:,k) = -flux_p(:,:,k)

  do k=1,nz-1
   E_lee_p(:,:,k,taup1) = E_lee_p(:,:,k,taup1) -  &
              maskW(:,:,k)*dt_loc*( flux_p(:,:,k)-flux_p(:,:,k-1) )/dzw(k) 
   E_lee_m(:,:,k,taup1) = E_lee_m(:,:,k,taup1) -  &
              maskW(:,:,k)*dt_loc*( flux_m(:,:,k)-flux_m(:,:,k-1) )/dzw(k) 
  enddo                                        
 enddo 
 
 mean_flow_to_E_lee = Tc_lee*(E_lee_p(:,:,:,taup1) - E_lee_m(:,:,:,taup1))
 aloc = maskW*dt_tracer*Tc_lee
 E_lee_p(:,:,:,taup1) = E_lee_p(:,:,:,taup1) + aloc*E_lee_p(:,:,:,taup1)
 E_lee_m(:,:,:,taup1) = E_lee_m(:,:,:,taup1) - aloc*E_lee_m(:,:,:,taup1)
  
 aloc = maskW*dt_tracer*(E_lee_p(:,:,:,taup1) - E_lee_m(:,:,:,taup1))/(2*tau_v_lee)
 E_lee_p(:,:,:,taup1) = E_lee_p(:,:,:,taup1) -  aloc
 E_lee_m(:,:,:,taup1) = E_lee_m(:,:,:,taup1) +  aloc
     
 !aloc = min(1d0/dt_tracer, alpha_c*E_iw(:,:,:,tau) )
 !bloc = aloc*E_lee_p(:,:,:,taup1)
 !iw_diss_lee = bloc
 !E_lee_p(:,:,:,taup1) = E_lee_p(:,:,:,taup1) -  maskW*dt_tracer*bloc 
 !bloc = aloc*E_lee_m(:,:,:,taup1)
 !iw_diss_lee = iw_diss_lee + bloc
 !E_lee_m(:,:,:,taup1) = E_lee_m(:,:,:,taup1) -  maskW*dt_tracer*bloc 

 aloc =  maskW*alpha_c*E_iw(:,:,:,tau)*(E_lee_p(:,:,:,taup1) + E_lee_m(:,:,:,taup1))/2.
 aloc = min( E_lee_p(:,:,:,taup1)/dt_tracer+1d-16, aloc )
 aloc = min( E_lee_m(:,:,:,taup1)/dt_tracer+1D-16, aloc )
 iw_diss_lee = aloc*2
 E_lee_p(:,:,:,taup1) = E_lee_p(:,:,:,taup1) -  dt_tracer*aloc
 E_lee_m(:,:,:,taup1) = E_lee_m(:,:,:,taup1) -  dt_tracer*aloc

 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_p(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_p(:,:,:,taup1))     
 call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_m(:,:,:,taup1)) 
 call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_m(:,:,:,taup1)) 
   
 if (enable_idemix_hor_diffusion) then
 !---------------------------------------------------------------------------------
 ! add tendency due to lateral diffusion
 !---------------------------------------------------------------------------------

  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
      flux_east(i,j,k)=Ah_lee * (E_lee_p(i+1,j,k,tau)-E_lee_p(i,j,k,tau))/(cost(j)*dxu(i))*maskUtr(i,j,k)
    enddo
   enddo
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     flux_north(i,j,k)= Ah_lee *  (E_lee_p(i,j+1,k,tau)-E_lee_p(i,j,k,tau))/dyu(j)*maskVtr(i,j,k)*cosu(j)
    enddo
   enddo
  enddo

  do j=js_pe,je_pe
    do i=is_pe,ie_pe
     E_lee_p(i,j,:,taup1)= E_lee_p(i,j,:,taup1) + dt_tracer*maskW(i,j,:)* &
                                  (( flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                  +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
  enddo
  
  do k=1,nz
   do j=js_pe,je_pe
    do i=is_pe-1,ie_pe
      flux_east(i,j,k)=Ah_lee *(E_lee_m(i+1,j,k,tau)-E_lee_m(i,j,k,tau))/(cost(j)*dxu(i))*maskUtr(i,j,k)
    enddo
   enddo
   do j=js_pe-1,je_pe
    do i=is_pe,ie_pe
     flux_north(i,j,k)= Ah_lee *(E_lee_m(i,j+1,k,tau)-E_lee_m(i,j,k,tau))/dyu(j)*maskVtr(i,j,k)*cosu(j)
    enddo
   enddo
  enddo

  do j=js_pe,je_pe
    do i=is_pe,ie_pe
     E_lee_m(i,j,:,taup1)= E_lee_m(i,j,:,taup1) + dt_tracer*maskW(i,j,:)* &
                                  (( flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                  +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
   enddo
  enddo

 endif

 if (enable_idemix_hor_diffusion_iter) then
 !---------------------------------------------------------------------------------
 ! add tendency due to lateral diffusion with iterative method in case of high resolution
 !---------------------------------------------------------------------------------

  fxa = dt_tracer/idemix_hor_diffusion_iter
  do n=1,idemix_hor_diffusion_iter

    call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_p(:,:,:,taup1)) 
    call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_p(:,:,:,taup1))
    call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_p(:,:,:,taup1)) 
    do k=1,nz
     do j=js_pe,je_pe
      do i=is_pe-1,ie_pe
        flux_east(i,j,k)=Ah_lee* &
            (E_lee_p(i+1,j,k,taup1)-E_lee_p(i,j,k,taup1))/(cost(j)*dxu(i))*maskUtr(i,j,k)
      enddo
     enddo
     do j=js_pe-1,je_pe
      do i=is_pe,ie_pe
       flux_north(i,j,k)= Ah_lee * &
            (E_lee_p(i,j+1,k,taup1)-E_lee_p(i,j,k,taup1))/dyu(j)*maskVtr(i,j,k)*cosu(j)
      enddo
     enddo
    enddo

    do j=js_pe,je_pe
      do i=is_pe,ie_pe
       E_lee_p(i,j,:,taup1)= E_lee_p(i,j,:,taup1) + fxa*maskW(i,j,:)* &
                                    (( flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                    +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
     enddo
    enddo
    
    call border_exchg_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_m(:,:,:,taup1)) 
    call setcyclic_xyz   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_m(:,:,:,taup1))
    call set_obc_boundary_xyz(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,E_lee_m(:,:,:,taup1)) 
    do k=1,nz
     do j=js_pe,je_pe
      do i=is_pe-1,ie_pe
        flux_east(i,j,k)=Ah_lee * &
            (E_lee_m(i+1,j,k,taup1)-E_lee_m(i,j,k,taup1))/(cost(j)*dxu(i))*maskUtr(i,j,k)
      enddo
     enddo
     do j=js_pe-1,je_pe
      do i=is_pe,ie_pe
       flux_north(i,j,k)= Ah_lee * &
            (E_lee_m(i,j+1,k,taup1)-E_lee_m(i,j,k,taup1))/dyu(j)*maskVtr(i,j,k)*cosu(j)
      enddo
     enddo
    enddo

    do j=js_pe,je_pe
      do i=is_pe,ie_pe
       E_lee_m(i,j,:,taup1)= E_lee_m(i,j,:,taup1) + fxa*maskW(i,j,:)* &
                                    (( flux_east(i,j,:) - flux_east(i-1,j,:))/(cost(j)*dxt(i)) &
                                    +(flux_north(i,j,:) -flux_north(i,j-1,:))/(cost(j)*dyt(j)) )
     enddo
    enddo    
 
  enddo

 endif
                                                                                                   
end subroutine 




subroutine leewave_flux
  use main_module
  use idemix_module
  implicit none
  integer   :: i,j,n,m,ks
  real*8 :: N_loc, f_loc   ! stability freq and abs of coriolis parameter
  real*8 :: u_0       ! phase velocity
  real*8 :: dk,S_fac
  real*8 :: wavnum(nk) ! wavenumbers
  real*8 :: uk(nk)     ! u_0*k 
  real*8 :: F_top(nk),freq_fac(nk)    ! topogr. spectrum, frequency factor
  real*8 :: flux_lee(nph),ustress_lee(nph),vstress_lee(nph),yN,yS,fxa,mu
  logical, save :: first = .true.


  nu = 0.8
  mu = 2*(nu+1)
  
  if (first) then
   topo_fac = 8*1.3*pi*h_rms**2*nu*k_s**(mu-2)
   
   ph(1) = -pi + 0.5*(2.0*pi/nph)
   do i = 2,nph
        ph(i) = ph(i-1) + 2.0*pi/nph
   enddo
   
   first = .false.
  endif
   
  
  do j = js_pe,je_pe
   do i = is_pe,ie_pe
    ks = max(1,kbot(i,j))
    u_bot(i,j) = 0.5*( u(i,j,ks,tau)*maskU(i,j,ks) + u(i-1,j,ks,tau)*maskU(i-1,j,ks) )
    v_bot(i,j) = 0.5*( v(i,j,ks,tau)*maskV(i,j,ks) + v(i,j-1,ks,tau)*maskV(i,j-1,ks) )
   enddo
  enddo
 
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,u_bot) 
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,u_bot)
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,v_bot) 
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,v_bot)
  u0_bot = sqrt(u_bot**2+v_bot**2)


if (enable_leewaves_slowflux) then

  ! slow version with numerical integration
    
  do j = js_pe,je_pe
   do i  = is_pe,ie_pe
    ks = kbot(i,j)
    if (ks>0) then
     if ((Nsqr(i,j,ks,tau) .le. 0) .or. (Nsqr(i,j,ks,tau) .le. coriolis_t(i,j)**2)) then
      flux_lee(:) = 0
      inv_fr(i,j) = 0. 
      ustress_lee(:) = 0
      vstress_lee(:) = 0 
     else 
      N_loc = sqrt(Nsqr(i,j,ks,tau))
      f_loc = abs(coriolis_t(i,j))
      do n = 1,nph
        u_0 = abs(u_bot(i,j)*cos(ph(n)) + v_bot(i,j)*sin(ph(n)))
        if (u_0>0) then
         !calculate wavenumber grid
         dk = (N_loc - f_loc)/(u_0*nk)
         wavnum(1) = f_loc/u_0 + .5*dk
         do m = 2,nk
           wavnum(m) = wavnum(m-1) + dk
         enddo
         !phase velocity times wavenumber
         uk(:) = u_0*wavnum(:)
         ! calculate the integrand
         where (N_loc**2>uk(:)**2 .and. uk(:)**2> f_loc**2)  
           freq_fac(:) = sqrt(N_loc**2 - uk(:)**2)*sqrt(uk(:)**2 - f_loc**2)
         else where
           freq_fac(:) = 0.
         end where       
         ! topography spectrum
         F_top(:) = h_rms(i,j)**2*nu/(pi*k_n(i,j)*k_s(i,j))*wavnum(:) &
                    /(1.+ (cos(ph(n)-ph_s(i,j))*wavnum(:)/k_s(i,j))**2  &
                        + (sin(ph(n)-ph_s(i,j))*wavnum(:)/k_n(i,j))**2)**(nu+1.0)       
         !integrate over wavenumber to get energy flux as a function of angle 
         flux_lee(n) = sum(4.0*pi**2*u_0*freq_fac(:)*F_top(:))*dk
         S_fac=1; if (u_bot(i,j)*cos(ph(n))+ v_bot(i,j)*sin(ph(n)) <0 ) S_fac=-1 
         ustress_lee(n) = sum(4.0*pi**2*cos(ph(n))*freq_fac(:)*F_top(:)*S_fac )*dk
         vstress_lee(n) = sum(4.0*pi**2*sin(ph(n))*freq_fac(:)*F_top(:)*S_fac )*dk   
         
        else ! u_0>0
         flux_lee(n) = 0. 
         ustress_lee(n) = 0
         vstress_lee(n) = 0 
        endif 
       enddo ! n=1,nph

       ! acount for critical Froud number
       fxa = k_s(i,j)*u0_bot(i,j)
       fxa = fxa**(mu-2)*h_rms(i,j)**2*nu/(mu-2)*(f_loc**(2-mu)-N_loc**(2-mu))
       inv_fr(i,j) = N_loc*sqrt(max(0d0,fxa))/max(u0_bot(i,j),1d-16)
      
      
       !inv_fr(i,j) = h_rms(i,j)*N_loc/max(u0_bot(i,j),1d-16)
       if (inv_fr(i,j) > inv_fr_c ) then
           flux_lee(:) = flux_lee(:)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
           ustress_lee(:) = ustress_lee(:)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
           vstress_lee(:) = vstress_lee(:)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
       endif    
      endif  ! N>0 and N>f
    else !ks>0
     flux_lee(:) = 0
     inv_fr(i,j) = 0.
     ustress_lee(:) = 0
     vstress_lee(:) = 0    
    endif  
    !Integrate over angle 
    flux_lee_tot(i,j) = sum(flux_lee(:))*(2*pi/nph)  
    stress_lee_u(i,j) = sum(ustress_lee(:))*(2*pi/nph)    
    stress_lee_v(i,j) = sum(vstress_lee(:))*(2*pi/nph)    
   enddo
  enddo

else

 ! faster version with approximated isotropic spectrum
 
 flux_lee_tot(:,:) = 0
 inv_fr(:,:) = 0. 
 stress_lee_u(:,:) = 0
 stress_lee_v(:,:) = 0 
  
  
  do j = js_pe,je_pe
   do i  = is_pe,ie_pe
    ks = kbot(i,j)
    if (ks>0) then
     if ((Nsqr(i,j,ks,tau) .gt. 0) .and. (Nsqr(i,j,ks,tau) .gt. coriolis_t(i,j)**2)) then
      N_loc = sqrt(Nsqr(i,j,ks,tau))
      f_loc = max(0.2e-4, abs(coriolis_t(i,j)) )
      fxa = (max(1.,N_loc/f_loc)-1.)**(nu-0.15)
      fxa = topo_fac(i,j)*N_loc**(-mu+4)*fxa
      flux_lee_tot(i,j) = fxa*u0_bot(i,j)**(mu-1)
      stress_lee_u(i,j) = fxa*u_bot(i,j)*u0_bot(i,j)**(mu-3) ! sign is wrong here, but corrected above
      stress_lee_v(i,j) = fxa*v_bot(i,j)*u0_bot(i,j)**(mu-3)
      ! account for critical Froud number
      
      
      fxa = k_s(i,j)*u0_bot(i,j)
      fxa = fxa**(mu-2)*h_rms(i,j)**2*nu/(mu-2)*(f_loc**(2-mu)-N_loc**(2-mu))
      inv_fr(i,j) = N_loc*sqrt(max(0d0,fxa))/max(u0_bot(i,j),1d-16)
      
      !inv_fr(i,j) = h_rms(i,j)*N_loc/max(u0_bot(i,j),1d-16)
      if (inv_fr(i,j) > inv_fr_c ) then
           flux_lee_tot(i,j) = flux_lee_tot(i,j)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
           stress_lee_u(i,j) = stress_lee_u(i,j)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
           stress_lee_v(i,j) = stress_lee_v(i,j)*(inv_fr_c/(inv_fr(i,j)+1d-16))**2
      endif    
     endif  ! N>0 and N>f
    endif ! ks>0
   enddo
  enddo

endif


 if (enable_leewaves_obc_taper_north) then
    yN = 0.0; if (my_blk_j == n_pes_j) yN = yu(ny)
    call global_max(yN)
    do j=js_pe-onx,je_pe+onx
      flux_lee_tot(:,j)=  flux_lee_tot(:,j)*(1- exp( -(yN-yt(j))**2/leewaves_north_taper_len**2  ) )
      stress_lee_u(:,j)=  stress_lee_u(:,j)*(1- exp( -(yN-yt(j))**2/leewaves_north_taper_len**2  ) )
      stress_lee_v(:,j)=  stress_lee_v(:,j)*(1- exp( -(yN-yt(j))**2/leewaves_north_taper_len**2  ) )
    enddo
 endif
 
 if (enable_leewaves_obc_taper_south) then
    yS = 1e20; if (my_blk_j == 1) yS = yu(0)
    call global_min(yS)
    do j=js_pe-onx,je_pe+onx
      flux_lee_tot(:,j)=  flux_lee_tot(:,j)*(1- exp( -(yt(j)-yS)**2/leewaves_south_taper_len**2  ) )
      stress_lee_u(:,j)=  stress_lee_u(:,j)*(1- exp( -(yt(j)-yS)**2/leewaves_south_taper_len**2  ) )
      stress_lee_v(:,j)=  stress_lee_v(:,j)*(1- exp( -(yt(j)-yS)**2/leewaves_south_taper_len**2  ) )
    enddo
 endif

  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,flux_lee_tot) 
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,flux_lee_tot) 
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,stress_lee_u) 
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,stress_lee_u) 
  call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,stress_lee_v) 
  call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,stress_lee_v)     
  
end subroutine 






subroutine leewave_superbee(flux,vel,var)
 use main_module
 use idemix_module
 implicit none
 integer :: k,i,j,kp2,km1
 real*8 :: flux(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,0:nz)
 real*8 :: var(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: vel(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: Rjp,Rj,Rjm,uCFL,Cr
 real*8 :: Limiter
 
 Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 
 
 do k=1,nz-2
       kp2=min(nz-1,k+2); !if (kp2>np) kp2=3
       km1=max(1,k-1) !if (km1<1) km1=np-2
       do j=js_pe,je_pe
        do i=is_pe,ie_pe
         Rjp=(var(i,j,kp2)-var(i,j,k+1))!*maskW(i,j,k+1)
         Rj =(var(i,j,k+1)-var(i,j,k  ))!*maskW(i,j,k  )
         Rjm=(var(i,j,k  )-var(i,j,km1))!*maskW(i,j,km1)
         uCFL = ABS( vel(i,j,k)*dt_tracer/dzw(k) )
         IF (Rj.NE.0.) THEN
          IF (vel(i,j,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (vel(i,j,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux(i,j,k) = vel(i,j,k)*(var(i,j,k+1)+var(i,j,k))*0.5d0   &
                                -ABS(vel(i,j,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
        enddo
       enddo
 enddo
end subroutine  leewave_superbee








