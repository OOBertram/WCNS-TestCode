subroutine  space_derivative_WCNS(rhs,u,ax,ay,az,bx,by,bz,cx,cy,cz,im,jm,km,le)  
use common1;use scheme_WCNS
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real rhs(im,jm,km,5)
real fn(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
dimension ww(-mo1:mo1+1,5)
real uu(5),ex,ey,ez ,a,b,c
real,parameter::error=1d-5,kapa=0.d0
select case(le)
case(1)
  do  k=1,km ;do  j=1,jm; do  i=-3,im+2
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call scheme_semi_WCNS_VF(fn(i,j,k,:),ww,ex,ey,ez)
    else
      call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)
    endif        
  enddo; enddo; enddo
case(2)
  do  k=1,km ;do  j=-3,jm+2; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call scheme_semi_WCNS_VF(fn(i,j,k,:),ww,ex,ey,ez)
    else
      call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)
    endif        
  enddo; enddo; enddo
case(3)
  do  k=-3,km+2 ;do  j=1,jm; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call scheme_semi_WCNS_VF(fn(i,j,k,:),ww,ex,ey,ez)
    else
      call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)
    endif        
  enddo; enddo; enddo
end select
  select case(id_scheme_WCNS_E)
  case(4) !WCNS-E4
    a=3.d0/8.d0*(3.d0-2.d0*kapa);b=1.d0/24.d0*(22.d0*kapa-1.d0)
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:)))   
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:)))   
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:)))   
      end select                  
    enddo; enddo; enddo
  case(5) !WCNS-E5
    a=75.d0/64.d0;b=-25.d0/384.d0;c=3.d0/640.d0
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:))+c*(fn(i+2,j,k,:)-fn(i-3,j,k,:))) 
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:))+c*(fn(i,j+2,k,:)-fn(i,j-3,k,:))) 
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:))+c*(fn(i,j,k+2,:)-fn(i,j,k-3,:))) 
      end select
    enddo; enddo; enddo
  end select

end subroutine

subroutine  space_derivative_WCNS_FT(rhs,u,ax,ay,az,bx,by,bz,cx,cy,cz,im,jm,km,le)  
use common1;use scheme_WCNS
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real rhs(im,jm,km,5)
real fn(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
dimension ww(-mo1:mo1+1,5)
real uu(5),ex,ey,ez ,a,b,c
real,parameter::error=1d-5,kapa=0.d0
select case(le)
case(1)
  do  k=1,km ;do  j=1,jm; do  i=-3,im+2
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call Linear(fn(i,j,k,:),ww,ex,ey,ez)
    else
      call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)
    endif        
  enddo; enddo; enddo
case(2)
  do  k=1,km ;do  j=-3,jm+2; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call Linear(fn(i,j,k,:),ww,ex,ey,ez)
    else
      call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)
    endif        
  enddo; enddo; enddo
case(3)
  do  k=-3,km+2 ;do  j=1,jm; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call Linear(fn(i,j,k,:),ww,ex,ey,ez)
    else
      call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)
    endif        
  enddo; enddo; enddo
end select
  select case(id_scheme_WCNS_E)
  case(4) !WCNS-E4
    a=3.d0/8.d0*(3.d0-2.d0*kapa);b=1.d0/24.d0*(22.d0*kapa-1.d0)
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:)))   
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:)))   
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:)))   
      end select                  
    enddo; enddo; enddo
  case(5) !WCNS-E5
    a=75.d0/64.d0;b=-25.d0/384.d0;c=3.d0/640.d0
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:))+c*(fn(i+2,j,k,:)-fn(i-3,j,k,:))) 
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:))+c*(fn(i,j+2,k,:)-fn(i,j-3,k,:))) 
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:))+c*(fn(i,j,k+2,:)-fn(i,j,k-3,:))) 
      end select
    enddo; enddo; enddo
  end select

end subroutine

subroutine  space_derivative_WCNS_CF(rhs,u,ax,ay,az,bx,by,bz,cx,cy,cz,im,jm,km,le)  
use common1;use scheme_WCNS
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real rhs(im,jm,km,5)
real fn(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
dimension ww(-mo1:mo1+1,5)
real uu(5),ex,ey,ez ,a,b,c
real,parameter::error=1d-5,kapa=0.d0
select case(le)
case(1)
  do  k=1,km ;do  j=1,jm; do  i=-3,im+2
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)       
  enddo; enddo; enddo
case(2)
  do  k=1,km ;do  j=-3,jm+2; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)       
  enddo; enddo; enddo
case(3)
  do  k=-3,km+2 ;do  j=1,jm; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez)    
  enddo; enddo; enddo
end select
  select case(id_scheme_WCNS_E)
  case(4) !WCNS-E4
    a=3.d0/8.d0*(3.d0-2.d0*kapa);b=1.d0/24.d0*(22.d0*kapa-1.d0)
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:)))   
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:)))   
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:)))   
      end select                  
    enddo; enddo; enddo
  case(5) !WCNS-E5
    a=75.d0/64.d0;b=-25.d0/384.d0;c=3.d0/640.d0
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:))+c*(fn(i+2,j,k,:)-fn(i-3,j,k,:))) 
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:))+c*(fn(i,j+2,k,:)-fn(i,j-3,k,:))) 
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:))+c*(fn(i,j,k+2,:)-fn(i,j,k-3,:))) 
      end select
    enddo; enddo; enddo
  end select

end subroutine

subroutine  space_derivative_WCNS_Zhao(rhs,u,ax,ay,az,bx,by,bz,cx,cy,cz,im,jm,km,le)!!!need to be amend
use common1;use scheme_WCNS
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real rhs(im,jm,km,5)
real fn(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
dimension ww(-mo1:mo1+1,5)
real uu(5),ex,ey,ez ,a,b,c
real,parameter::error=1d-5,kapa=0.d0
select case(le)
case(1)
  do  k=1,km ;do  j=1,jm; do  i=-3,im+2
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    call scheme_semi_WCNS_Zhao(fn(i,j,k,:),ww,ex,ey,ez)       
  enddo; enddo; enddo
case(2)
  do  k=1,km ;do  j=-3,jm+2; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    call scheme_semi_WCNS_Zhao(fn(i,j,k,:),ww,ex,ey,ez)       
  enddo; enddo; enddo
case(3)
  do  k=-3,km+2 ;do  j=1,jm; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    call scheme_semi_WCNS_Zhao(fn(i,j,k,:),ww,ex,ey,ez)    
  enddo; enddo; enddo
end select
  select case(id_scheme_WCNS_E)
  case(4) !WCNS-E4
    a=3.d0/8.d0*(3.d0-2.d0*kapa);b=1.d0/24.d0*(22.d0*kapa-1.d0)
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:)))   
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:)))   
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:)))   
      end select                  
    enddo; enddo; enddo
  case(5) !WCNS-E5
    a=75.d0/64.d0;b=-25.d0/384.d0;c=3.d0/640.d0
    do k=1,km ;do j=1,jm; do i=1,im
      select case(le)
      case(1)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:))+c*(fn(i+2,j,k,:)-fn(i-3,j,k,:))) 
      case(2)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:))+c*(fn(i,j+2,k,:)-fn(i,j-3,k,:))) 
      case(3)
        rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:))+c*(fn(i,j,k+2,:)-fn(i,j,k-3,:))) 
      end select
    enddo; enddo; enddo
  end select

end subroutine

subroutine  space_derivative_WCNS_update(rhs,u,ax,ay,az,bx,by,bz,cx,cy,cz,im,jm,km,le,th,vel)  
use common1;use scheme_WCNS
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real rhs(im,jm,km,5)
integer :: w =2
real distance(3),vel(3),vel_shock,vels
real fn(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
dimension ww(-mo1:mo1+1,5)
real uu(5),ex,ey,ez ,a,b,c,th,eqs
real,parameter::error=1d-5,kapa=0.d0
distance(:) = 0.0
eps = 1.0d-3
vel_shock = 0.0
select case(le)
case(1)
vels = 0.0
  do  k=1,km ;do  j=1,jm; do  i=-3,im+2
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call Linear(fn(i,j,k,:),ww,ex,ey,ez)
      vel_s1(i,j,k) = 0.0
    else
      call scheme_semi_WCNS_update(fn(i,j,k,:),ww,ex,ey,ez,vel_s1(i,j,k),le)
      ! call scheme_semi_WCNS(fn(i,j,k,:),ww,ex,ey,ez,le)
    endif 
  enddo; enddo; enddo

  do  k=1,km ;do  j=1,jm; do  i=-3,im+2
  if(signal(i,j,k,le) .eqv. .false.)then
    ! vel = (abs(u(i-w,j,k,2)-u(i+w,j,k,2))+eps) / (abs(u(i-w,j,k,1)-u(i+w,j,k,1))+eps)
    if(vel_s1(i,j,k) > vels) vels = vel_s1(i,j,k)
    !write(*,*) vel_s(i,j,k)
  endif
  enddo;enddo;enddo
  vel(1) = vels
  !write(*,*) vel(1)
case(2)
  do  k=1,km ;do  j=-3,jm+2; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call Linear(fn(i,j,k,:),ww,ex,ey,ez)
      vel_s2(i,j,k) = 0.0
    else
      call scheme_semi_WCNS_update(fn(i,j,k,:),ww,ex,ey,ez,vel_s2(i,j,k),le)
    endif 
  enddo; enddo; enddo
  do  k=1,km ;do  j=-3,jm+2; do  i=1,im
  if(signal(i,j,k,le) .eqv. .false.)then
    ! vel = (abs(u(i,j-w,k,3)-u(i,j+w,k,3))+eps) / (abs(u(i,j-w,k,1)-u(i,j+w,k,1))+eps)
    if(vels < vel_s2(i,j,k)) vels = vel_s2(i,j,k)
    !write(*,*) vel_s(i,j,k)
  endif
  enddo;enddo;enddo
vel(2) = vels
!  write(*,*) vel(2)
case(3)
vels = 0.0
  do  k=-3,km+2 ;do  j=1,jm; do  i=1,im
    call neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,ww,ex,ey,ez,im,jm,km,i,j,k,le)	  
    if(signal(i,j,k,le) .eqv. .true.)then
      call Linear(fn(i,j,k,:),ww,ex,ey,ez)
      vel_s3(i,j,k) = 0.0
    else
      call scheme_semi_WCNS_update(fn(i,j,k,:),ww,ex,ey,ez,vel_s3(i,j,k),le)
    endif 
  enddo; enddo; enddo
  do  k=-3,km+2 ;do  j=1,jm; do  i=1,im
  if(signal(i,j,k,le) .eqv. .false.)then
    !vel = (abs(u(i,j,k-w,4)-u(i,j,k+w,4))+eps) / (abs(u(i,j,k-w,1)-u(i,j,k+w,1))+eps)
    if(vel_s3(i,j,k) > vels) vels = vel_s3(i,j,k)
    !write(*,*) vel_s(i,j,k)
  endif
  enddo;enddo;enddo
vel(3) = vels
end select
select case(id_scheme_WCNS_E)
case(4) !WCNS-E4
  a=3.d0/8.d0*(3.d0-2.d0*kapa);b=1.d0/24.d0*(22.d0*kapa-1.d0)
  do k=1,km ;do j=1,jm; do i=1,im
    select case(le)
    case(1)
      rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:)))   
    case(2)
      rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:)))   
    case(3)
      rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:)))   
    end select                  
  enddo; enddo; enddo
case(5) !WCNS-E5
  a=75.d0/64.d0;b=-25.d0/384.d0;c=3.d0/640.d0
  do k=1,km ;do j=1,jm; do i=1,im
    select case(le)
    case(1)
      rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i-1,j,k,:))+b*(fn(i+1,j,k,:)-fn(i-2,j,k,:))+c*(fn(i+2,j,k,:)-fn(i-3,j,k,:))) 
    case(2)
      rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j-1,k,:))+b*(fn(i,j+1,k,:)-fn(i,j-2,k,:))+c*(fn(i,j+2,k,:)-fn(i,j-3,k,:))) 
    case(3)
      rhs(i,j,k,:)=-(a*(fn(i,j,k,:)-fn(i,j,k-1,:))+b*(fn(i,j,k+1,:)-fn(i,j,k-2,:))+c*(fn(i,j,k+2,:)-fn(i,j,k-3,:))) 
    end select
  enddo; enddo; enddo
end select

end subroutine


subroutine neighbor_point_WCNS(u,ax,ay,az,bx,by,bz,cx,cy,cz,w,ex,ey,ez,im,jm,km,i,j,k,le)
use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real w(-mo1:mo1+1,5)
integer i,j,k,im,jm,km,i1,j1,k1,le
i1=max(1,min(i,im));j1=max(1,min(j,jm));k1=max(1,min(k,km))
select case(le)
case(1)
    ex=ax(i1,j1,k1)
    ey=ay(i1,j1,k1)
    ez=az(i1,j1,k1)
do l=-mo1,mo1+1
  w(l,:)=u(i+l,j,k,:) 
enddo
case(2)
    ex=bx(i1,j1,k1)
    ey=by(i1,j1,k1)
    ez=bz(i1,j1,k1)
do l=-mo1,mo1+1
    w(l,:)=u(i,j+l,k,:)
enddo
case(3)
    ex=cx(i1,j1,k1)
    ey=cy(i1,j1,k1)
    ez=cz(i1,j1,k1)
do l=-mo1,mo1+1
  w(l,:)=u(i,j,k+l,:)
enddo
end select
end subroutine


subroutine Linear(fn,u,ex,ey,ez)
use input;use ghost;use scheme_WCNS;
real fn(5),u(-mo1:mo1+1,5)
real uL(5),uR(5),fL(5),fR(5),LL(5,5),RR(5,5),aLR(5),fLR(5)
real LuL(5),LuR(5),Lfn(5),LfLR(5),Ldu(5),a
integer,parameter::fluxtype=1
integer :: k
do k=1,5
  call Linear_inter(LuL(k),u(-2:3,k))
enddo
  LuR(:) = LuL(:)
  uL = LuL;uR = LuR
selectcase(fluxtype)
  case(1)!lax
alpha=max(armm(uL,ex,ey,ez),armm(uR,ex,ey,ez))
call flux(fL,uL,ex,ey,ez)
call flux(fR,uR,ex,ey,ez)
do k=1,5
  fn(k)=0.5*(fL(k)+fR(k)-alpha*(uR(k)-uL(k)))
enddo
case(2)!roe
call eigno(fLR,aLR,ll,rr,uL,uR,ex,ey,ez) 
  Ldu=matmul(ll,uR-uL)
do k=1,5
  fn(k)=fLR(k)-0.5* (rr(k,1)*(abs(aLR(1))*Ldu(1))&
                    +rr(k,2)*(abs(aLR(2))*Ldu(2))&
                    +rr(k,3)*(abs(aLR(3))*Ldu(3))&
                    +rr(k,4)*(abs(aLR(4))*Ldu(4))&
                    +rr(k,5)*(abs(aLR(5))*Ldu(5)))      
enddo
end select
endsubroutine

subroutine Linear_inter(LuL,u)
real LuL(5),u(-2:3),con
integer i
con = 1.0/256.0
  LuL = 3.0*con*u(-2)-25.0*con*u(-1)+150.0*con*u(0)  &
  +150.0*con*u(1)-25.0*con*u(2)+3.0*con*u(3)
endsubroutine
