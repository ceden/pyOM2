


 
 subroutine vertical_velocity
!=======================================================================
!       vertical velocity from continuity : 
!       \int_0^z w_z dz =w(z)-w(0) = - \int dz (u_x +v_y)  
!        w(z)=-int dz u_x + v_y
!=======================================================================
    use main_module   
    implicit none
    integer :: i,j,k
    ! integrate from bottom to surface to see error in w
    k=1
    do j=js_pe-onx+1,je_pe+onx
     do i=is_pe-onx+1,ie_pe+onx
         w(i,j,k,taup1) =-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,taup1)-          u(i-1,j,k,taup1))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,taup1)-cosu(j-1)*v(i,j-1,k,taup1))/(cost(j)*dyt(j)) )
     enddo
    enddo
    do k=2,nz
     do j=js_pe-onx+1,je_pe+onx
      do i=is_pe-onx+1,ie_pe+onx
          w(i,j,k,taup1) = w(i,j,k-1,taup1)-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,taup1)          -u(i-1,j,k,taup1))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,taup1)-cosu(j-1)*v(i,j-1,k,taup1))/(cost(j)*dyt(j)) )
      enddo
     enddo
    enddo
end subroutine vertical_velocity



 
 subroutine vertical_velocity_tau(ntau)
!=======================================================================
!       vertical velocity from continuity : 
!       \int_0^z w_z dz =w(z)-w(0) = - \int dz (u_x +v_y)  
!        w(z)=-int dz u_x + v_y
!=======================================================================
    use main_module   
    implicit none
    integer, intent(in) :: ntau
    integer :: i,j,k
    ! integrate from bottom to surface to see error in w
    k=1
    do j=js_pe-onx+1,je_pe+onx
     do i=is_pe-onx+1,ie_pe+onx
         w(i,j,k,ntau) =-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,ntau)-          u(i-1,j,k,ntau))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,ntau)-cosu(j-1)*v(i,j-1,k,ntau))/(cost(j)*dyt(j)) )
     enddo
    enddo
    do k=2,nz
     do j=js_pe-onx+1,je_pe+onx
      do i=is_pe-onx+1,ie_pe+onx
          w(i,j,k,ntau) = w(i,j,k-1,ntau)-maskW(i,j,k)*dzt(k)* &
               ((        u(i,j,k,ntau)          -u(i-1,j,k,ntau))/(cost(j)*dxt(i)) &
               +(cosu(j)*v(i,j,k,ntau)-cosu(j-1)*v(i,j-1,k,ntau))/(cost(j)*dyt(j)) )
      enddo
     enddo
    enddo
end subroutine vertical_velocity_tau





