subroutine spacederivativek_WCNS(rhs,ukm,th,le)
use common1
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
real rhs(i1,j1,k1,5)
call space_derivative_WCNS(rhs,ukm,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,i1,j1,k1,le)
end  subroutine  

subroutine spacederivativek_WCNS_FT(rhs,ukm,th,le)
use common1
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
real rhs(i1,j1,k1,5)
call space_derivative_WCNS_FT(rhs,ukm,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,i1,j1,k1,le)
end  subroutine  

subroutine spacederivativek_WCNS_CF(rhs,ukm,th,le)
use common1
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
real rhs(i1,j1,k1,5)
call space_derivative_WCNS_CF(rhs,ukm,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,i1,j1,k1,le)
end  subroutine  

subroutine spacederivativek_WCNS_update(rhs,ukm,th,le,vel)
use common1
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5),vel(3)
real rhs(i1,j1,k1,5)
call space_derivative_WCNS_update(rhs,ukm,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,i1,j1,k1,le,th,vel(:))
end  subroutine 

real function eigenmax(le)
use common1
implicit none
integer le
real,external::eigen_max
eigenmax=eigen_max(u1,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,i1,j1,k1,le)               	 
end function

function eigen_max(u,ax,ay,az,bx,by,bz,cx,cy,cz,im,jm,km,le)
use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
dimension uijk(5)
select case(le)
case(1)
  eigen_max=0.
do k=1,km;do j=1,jm;do i=1,im;do l=1,5
 	uijk(l)=u(i,j,k,l)
  eigen_max=amax1(eigen_max,armm(uijk,ax(i,j,k),ay(i,j,k),az(i,j,k)))
enddo;enddo;enddo;enddo
return
case(2)
  eigen_max=0.
do k=1,km;do j=1,jm;do i=1,im;do l=1,5
  uijk(l)=u(i,j,k,l)
	eigen_max= amax1(eigen_max,armm(uijk,bx(i,j,k),by(i,j,k),bz(i,j,k)))
enddo;enddo;enddo;enddo
case(3)
  eigen_max=0.
do k=1,km;do j=1,jm;do i=1,im;do l=1,5
  uijk(l)=u(i,j,k,l)
  eigen_max=amax1(eigen_max,armm(uijk,cx(i,j,k),cy(i,j,k),cz(i,j,k)))
enddo;enddo;enddo;enddo
end select
end function



real function armm(u,ex,ey,ez)	!ex,ey,ez������������ֵ
use input
real u(5)
! write(*,*) ex,ey,ez
ue=(ex*u(2)+ey*u(3)+ez*u(4))/u(1)
ce=sqrt((ex**2+ey**2+ez**2)*abs(gin*p(u)/u(1)))
! ce=abs(gin*p(u)/u(1))
armm=abs(ue)+ce
endfunction