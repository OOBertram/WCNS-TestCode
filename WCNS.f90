subroutine ubar_WCNS_E5(ubar,u)
implicit none
real u(-2:3),d1u(3),d2u(3),dx,ubar
real u3(3),u5
real a(3),w(3),wu,bet(3),sum_bet,is(3),wu5
real,parameter::error=1D-6
real,external::armmax,minmod
integer num_dab,fluxtype,i,j,k

d1u(1)=1.d0/2.d0*(u(-2)-4.d0*u(-1)+3.d0*u(0))
d1u(2)=1.d0/2.d0*(u(1)-u(-1))
d1u(3)=1.d0/2.d0*(-3.d0*u(0)+4.d0*u(1)-u(2))

d2u(1)=1.d0/2.d0*(u(-2)-2.d0*u(-1)+u(0))
d2u(2)=1.d0/2.d0*(u(-1)-2.d0*u(0)+u(1))
d2u(3)=1.d0/2.d0*(u(0)-2.d0*u(1)+u(2))

do k=1,3
u3(k)=u(0)+1.d0/2.d0*d1u(k)+1.d0/4.d0*d2u(k)
enddo

a(1)=1.d0/16.d0; a(2)=10.d0/16.d0; a(3)=5.d0/16.d0

u5=a(1)*u3(1)+a(2)*u3(2)+a(3)*u3(3)

do k=1,3
    is(k)=d1u(k)**2+4.d0*d2u(k)**2.d0
    bet(k)=a(k)/(error+is(k))**2.d0
enddo

sum_bet=1.d0/(bet(1)+bet(2)+bet(3))

do k=1,3; w(k)=bet(k)*sum_bet;enddo      

wu5=0.
do k=1,3; wu5=wu5+w(k)*u3(k);enddo
ubar=wu5

end subroutine
!************************WCNS-WENOZ***************************
subroutine ubar_WCNS_E5_Z(ubar,u) 
implicit none
real u(-2:3),d1u(3),d2u(3),dx,ubar,v
real u3(3),u5
real a(3),w(3),bet(3),betZ(3),sum_bet,is(3),wu5,tao5
real,parameter::error=1D-40,C=1.d0,q=1.d0
real,external::armmax,minmod
integer fluxtype,i,j,k

d1u(1)=1.d0/2.d0*(u(-2)-4.d0*u(-1)+3.d0*u(0))
d1u(2)=1.d0/2.d0*(u(1)-u(-1))
d1u(3)=1.d0/2.d0*(-3.d0*u(0)+4.d0*u(1)-u(2))

d2u(1)=1.d0/2.d0*(u(-2)-2.d0*u(-1)+u(0))
d2u(2)=1.d0/2.d0*(u(-1)-2.d0*u(0)+u(1))
d2u(3)=1.d0/2.d0*(u(0)-2.d0*u(1)+u(2))

do k=1,3; u3(k)=u(0)+1.d0/2.d0*d1u(k)+1.d0/4.d0*d2u(k); enddo
 a(1)=1.d0/16.d0; a(2)=10.d0/16.d0; a(3)=5.d0/16.d0
do k=1,3; is(k)=d1u(k)**2.0+4.d0*d2u(k)**2.0; enddo
tao5=abs(is(1)-is(3))
do k=1,3; bet(k)=a(k)*(C+(tao5/(is(k)+error)))**1.0;enddo   

sum_bet=1.d0/(bet(1)+bet(2)+bet(3))

do k=1,3; w(k)=bet(k)*sum_bet;enddo     

  wu5=0.d0
do k=1,3; wu5=wu5+w(k)*u3(k);enddo
ubar=wu5
end subroutine

! subroutine ubar_WCNS_E5_Z(ubar,u) 
! implicit none
! real u(-2:3),d1u(3),d2u(3),dx,ubar,v
! real Qp(-2:3),beta0,beta1,beta2,beta3,c1,c2,c3,c4,c5,tao
! real alfa0,alfa1,alfa2,alfa3,epsilon,alfasum
! real Qhalf0,Qhalf1,Qhalf2,Qhalf3
!   epsilon = 1.0d-40
!   Qp(:)=u(:);
!   beta0=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2;
!   beta1=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2+0.25D0*(Qp(-1)-Qp(1))**2;
!   beta2=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2;
!   tao=abs(beta0-beta2);
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   alfa0= 1.0D0*(1.0D0+tao/(beta0+epsilon));
!   alfa1=10.0D0*(1.0D0+tao/(beta1+epsilon));
!   alfa2= 5.0D0*(1.0D0+tao/(beta2+epsilon));

!   alfasum = alfa0 + alfa1 + alfa2
      
!   alfa0 = alfa0 / alfasum
!   alfa1 = alfa1 / alfasum
!   alfa2 = alfa2 / alfasum

!   Qhalf0=3.0D0*Qp(-2)-10.0D0*Qp(-1)+15.0D0*Qp(0);
!   Qhalf1=-Qp(-1)+6.0D0*Qp(0)+3.0D0*Qp(1);
!   Qhalf2=3.0D0*Qp(0)+6.0D0*Qp(1)-Qp(2);
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ubar=alfa0*Qhalf0+alfa1*Qhalf1+alfa2*Qhalf2;
!   ubar=0.125D0*ubar/(alfa0+alfa1+alfa2);
! end subroutine
!***********************WCNS-TENO**********************
subroutine ubar_WCNS_E5_T(ubar,u)
implicit none
real u(-2:3),d1u(3),d2u(3),dx,ubar,v
real u3(3),u5
real a(3),w(3),bet(3),sum_bet,is(3),wu5
real O(3),d(3),x(3),r(3),ad(3),tao5,sum_R,sum_ad
real,parameter::error=1D-8,C=1.d0,q=6.d0,Ct=1E-5
real,external::armmax,minmod
integer fluxtype,i,j,k

d1u(1)=1.d0/2.d0*(u(-2)-4.d0*u(-1)+3.d0*u(0))
d1u(2)=1.d0/2.d0*(u(1)-u(-1))
d1u(3)=1.d0/2.d0*(-3.d0*u(0)+4.d0*u(1)-u(2))

d2u(1)=1.d0/2.d0*(u(-2)-2.d0*u(-1)+u(0))
d2u(2)=1.d0/2.d0*(u(-1)-2.d0*u(0)+u(1))
d2u(3)=1.d0/2.d0*(u(0)-2.d0*u(1)+u(2))

do k=1,3; u3(k)=u(0)+1.d0/2.d0*d1u(k)+1.d0/4.d0*d2u(k); enddo
 a(1)=1.d0/16.d0; a(2)=10.d0/16.d0; a(3)=5.d0/16.d0
do k=1,3; is(k)=d1u(k)**2.d0+4.d0*d2u(k)**2.d0; enddo
!****scale separation****
tao5=abs(is(1)-is(3))
do k=1,3; r(k)=(C+(tao5/(is(k)+error)))**q;enddo   
sum_r=1.d0/(r(1)+r(2)+r(3))
do k=1,3
    x(k)=r(k)*sum_r
    if(x(k)<Ct) then; d(k)=1.0d-20; else; d(k)=1.d0; endif
enddo 
!****weight process****    
do k=1,3; ad(k)=a(k)*d(k); enddo
sum_ad=1.d0/(ad(1)+ad(2)+ad(3))
do k=1,3; w(k)=ad(k)*sum_ad; enddo    
wu5=0.d0
do k=1,3; wu5=wu5+w(k)*u3(k); enddo
ubar=wu5
! ubar=a(1)*u3(1)+a(2)*u3(2)+a(3)*u3(3)
end subroutine 


!****************WCNS-CU6***************
subroutine ubar_WCNS_E5_CU6(ubar,u) 
implicit none
real u(-2:3),d1u(3),d2u(3),dx,ubar,v
real Qp(-2:3),beta0,beta1,beta2,beta3,c1,c2,c3,c4,c5,tao
real alfa0,alfa1,alfa2,alfa3,epsilon
real Qhalf0,Qhalf1,Qhalf2,Qhalf3
  epsilon = 1.0d-40
  Qp(:)=u(:);
  beta0=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2.0+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2.0;
  beta1=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2.0+0.25D0*(Qp(-1)-Qp(1))**2.0;
  beta2=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2.0+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2.0;
  
  c1=(-3.0D0*Qp(-2)+30.0D0*Qp(-1)+20.0D0*Qp(0)-60.0D0*Qp(1)+15.0D0*Qp(2)-2.0D0*Qp(3))/60.0D0;
  c2=(Qp(-2)-16.0D0*Qp(-1)+30.0D0*Qp(0)-16.0D0*Qp(1)+Qp(2))/12.0D0;
  c3=(Qp(-2)+Qp(-1)-10.0D0*Qp(0)+14.0D0*Qp(1)-7.0D0*Qp(2)+Qp(3))/4.0D0;
  c4=-Qp(-2)+4.0D0*Qp(-1)-6.0D0*Qp(0)+4.0D0*Qp(1)-Qp(2);
  c5=Qp(-2)-5.0D0*Qp(-1)+10.0D0*Qp(0)-10.0D0*Qp(1)+5.0D0*Qp(2)-Qp(3);
  beta3=c1**2.0+c2**2.0+c3**2.0+c4**2.0+c5**2.0;
  tao=abs(beta3-(beta0+4.0D0*beta1+beta2)/6.0D0);
  alfa0=  1.0D0 / 32.0D0*(1.0D0+tao/(beta0+epsilon))**1.0;
  alfa1= 15.0D0 / 32.0D0*(1.0D0+tao/(beta1+epsilon))**1.0;
  alfa2= 15.0D0 / 32.0D0*(1.0D0+tao/(beta2+epsilon))**1.0;
  alfa3=  1.0D0 / 32.0D0*(1.0D0+tao/(beta3+epsilon))**1.0;
  Qhalf0=(3.0D0*Qp(-2)-10.0D0*Qp(-1)+15.0D0*Qp(0));
  Qhalf1=(-Qp(-1)+6.0D0*Qp(0)+3.0D0*Qp(1));
  Qhalf2=(3.0D0*Qp(0)+6.0D0*Qp(1)-Qp(2));
  Qhalf3=(15.0D0*Qp(1)-10.0D0*Qp(2)+3.0D0*Qp(3));
  ubar = alfa0*Qhalf0+alfa1*Qhalf1+alfa2*Qhalf2+alfa3*Qhalf3;
  ubar = 0.125*ubar/(alfa0+alfa1+alfa2+alfa3);
end subroutine