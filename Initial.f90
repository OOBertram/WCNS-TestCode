subroutine readinput   
use input
implicit none
namelist /inputdata/qin,pin,gin,sm,cfl,cflvis, &
id_scheme,timebegin,timepart,timestop,interval     
open(1,file='./inputdata.txt',status='old')                
read(1,inputdata) 
close(1)
timeend=timebegin + timepart
uin=sm*sqrt(gin*pin/qin)
write(*,inputdata) 
ncfl=ceiling(CFL-1.e-6)
end subroutine

subroutine initialvalue(index)  
use common1
integer ::  index
select case(index)
case(1)
    call initial_value_DoubleMach(u1,x1,y1,i1,j1,k1) 
case(2)
    call initial_value_2DRiemann
case(3)
    call initial_value_2DRT(u1,x1,y1,i1,j1,k1) 
case(4)
    call initialvalue_Richtmyer_Meshkov
case(5)
    call initialvalue_Shock_Vortex
case(6)
    call initialvalue_IRV
case(7)
    call initial_value_2DRiemann2
case(8)
    call initial_value_2DRiemann3
end select
end

subroutine initial_value_DoubleMach(u,x,y,im,jm,km)	
use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real x(im,jm,km),y(im,jm,km)
do  k=1,km	;  do  j=1,jm	; do  i=1,im	 
if(y(i,j,k)<(x(i,j,k)-1./6.)*sqrt(3.))then
    call conserve(u(i,j,k,:),1.4,0.,0.,0.,1.0)
else
    call conserve(u(i,j,k,:),8.0,7.14471,-4.125,0.,116.5)
endif
end do ;end do ; end do 
end

subroutine initial_value_2DRiemann 
use common1 
do k=1,k1; do j=1,j1; do i=1,i1 
if(x1(i,j,k)<=0.8.and.y1(i,j,k)<=0.8)  then 
    call  conserve(u1(i,j,k,:),0.138, 1.206, 1.206,0.,0.029)    
else if (x1(i,j,k)>=0.8.and.y1(i,j,k)<0.8)  then
    call  conserve(u1(i,j,k,:),0.5323,0., 1.206, 0.,0.3)   
else if (x1(i,j,k)>=0.8.and.y1(i,j,k)>=0.8)  then           
    call  conserve(u1(i,j,k,:),1.5,0.,0., 0.,1.5)          
else 
    call  conserve(u1(i,j,k,:),0.5323,1.206 ,0.,0.,0.3)
endif
end do; end do; end do
end subroutine  

subroutine initial_value_2DRiemann2   
use common1 
do k=1,k1; do j=1,j1; do i=1,i1 
if(x1(i,j,k)<=0.5.and.y1(i,j,k)<=0.5)  then 
    call  conserve(u1(i,j,k,:),1.0, -0.75, 0.5,0.,1.0)    
else if (x1(i,j,k)>=0.5.and.y1(i,j,k)<0.5)  then
    call  conserve(u1(i,j,k,:),3.0,-0.75, -0.5, 0.,1.0)
else if (x1(i,j,k)>=0.5.and.y1(i,j,k)>=0.5)  then           
    call  conserve(u1(i,j,k,:),1.0,0.75,-0.5, 0.,1.0)    
else 
    call  conserve(u1(i,j,k,:),2.0,0.75 ,0.5,0.,1.0)
endif
end do; end do; end do
end subroutine  

subroutine initial_value_2DRiemann3   
use common1 
do k=1,k1; do j=1,j1; do i=1,i1 
if(x1(i,j,k)<=0.5.and.y1(i,j,k)<=0.5)  then 
    call  conserve(u1(i,j,k,:),1.0, 0.0, 0.7276,0.,1.0)    
else if (x1(i,j,k)>=0.5.and.y1(i,j,k)<0.5)  then
    call  conserve(u1(i,j,k,:),0.8,0.0,0.0, 0.,1.0)    
else if (x1(i,j,k)>=0.5.and.y1(i,j,k)>=0.5)  then           
    call  conserve(u1(i,j,k,:),0.5313,0.0,0.0, 0.,0.4)          
else 
    call  conserve(u1(i,j,k,:),1.0,0.7276 ,0.0,0.,1.0)
endif
end do; end do; end do
end subroutine  

subroutine initial_value_2DRT(u,x,y,im,jm,km)	 
use ghost
use input
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real x(im,jm,km),y(im,jm,km)
real,parameter::pi=3.1415926
do  k=1,km	 ; do  j=1,jm ;do  i=1,im	 
if(y(i,j,k)<=0.5)   then 
  qs=2.0;us=0.0; ps=2.0*y(i,j,k)+1.0 ; vs=-0.025*sqrt(abs(gin*ps/qs))*cos(8.0*pi*x(i,j,k))
  call  conserve(u(i,j,k,:),qs ,us ,vs,0.,ps)
else
  qs=1.0;us=0.0; ps= y(i,j,k)+1.5;  vs=-0.025*sqrt(abs(gin*ps/qs))*cos(8.0*pi*x(i,j,k))
  call  conserve(u(i,j,k,:),qs ,us ,vs,0.,ps)
endif
end do; end do; end do
end subroutine

subroutine initialvalue_Richtmyer_Meshkov 
use common1 
real,parameter::pi=3.1415926535897932384626433832795
real Ma,qsc,gs,xx0,yy0,r0,xx1,ut,deltr,Ts,qsL,usL,psL,qsR,usR,psR,Ma2,qs,us,ps
gs=1.4;xx0=4.;yy0=2.5;r0=1.;xx1=5.5;ut=0.8;deltr=0.02;TsL=1.0
!Ma=1.093;qsc=0.138;
Ma=1.23;qsc=0.166;

 qsL=1.;usL=1.;psL=1.
 usR=(2.+(gs-1.)*Ma*Ma)/((gs+1.)*Ma*Ma)
 qsR=((gs+1.)*Ma*Ma)/(2.+(gs-1.)*Ma*Ma)
 psR=(2*gs*Ma*Ma-(gs-1))/(gs+1)
 TsR=(1.+0.5*(gs-1.)*Ma*Ma)*((2.*gs/(gs-1.)*Ma*Ma)-1.)/( (gs+1)**2/(2.*(gs-1.))*Ma*Ma)
 Ma2=Ma*usR/sqrt(TsR)
 
do k=1,k1; do j=1,j1; do i=1,i1 
    rr=sqrt((x1(i,j,k)-xx0)**2+(y1(i,j,k)-yy0)**2)
if(x1(i,j,k)>=xx1)  then
    qs=qsR;us=usR-ut;ps=TsR*qs
elseif(rr>=r0+0.5*deltr)  then
    qs=1.;us=1.-ut;ps=TsL*qs
elseif(rr<=r0-0.5*deltr)  then
    qs=qsc;us=1.-ut;ps=1.
else !   r0-0.5*deltr< rr < r0 + 0.5*deltr 
    qs=qsc+(1-qsc)*(rr-(r0-0.5*deltr))/deltr;us=1.-ut;ps=1.
endif
call  conserve(u1(i,j,k,:),qs,us ,0.,0.,ps)
end do; end do; end do
end subroutine 

subroutine initialvalue_Shock_Vortex	
use common1 ;use input 
!implicit none
gs = 1.4
do k=1,k1; do j=1,j1; do i=1,i1 
xc=0.25; yc=0.5 ; r=sqrt((x1(i,j,k)-xc)**2 + (y1(i,j,k)-yc)**2)
eps=0.3;rc=0.05;tau=r/rc;alpha=0.204
uu=eps*tau*exp(alpha*(1-tau**2))*(y1(i,j,k)-yc)/r
vv=-eps*tau*exp(alpha*(1-tau**2))*(x1(i,j,k)-xc)/r
TT=-(gs-1)*eps**2*exp(2*alpha*(1-tau**2))/(4*alpha*gs)
if(x1(i,j,k)<=0.5)   then
S=0.;T=1.
qs=((T+TT)/exp(S))**(1./(gs-1))
ps=(T+TT)*qs
us=1.1*sqrt(1.4)+uu
vs=vv
  call  conserve(u1(i,j,k,:),qs,us,vs,0.,ps)  
else
  qs=1.*(((gs+1.)*1.1**2)/(2+(gs-1.)*1.1**2))
  ps=1.*(2*gs*1.1**2-(gs-1))/(gs+1)
  us=1.*1.1*sqrt(1.4)/qs+uu
  vs=0.
 call  conserve(u1(i,j,k,:),qs,us,vs,0.,ps) 
endif 
enddo;enddo;enddo
end subroutine 

subroutine initialvalue_IRV    
use common1 
do k=1,k1; do j=1,j1; do i=1,i1 
if(x1(i,j,k)>=0.5.and.y1(i,j,k)>=0.5)  then 
   call  conserve(u1(i,j,k,:),1., 0.1, 0.1,0.,1.)    
else if (x1(i,j,k)<0.5.and.y1(i,j,k)>=0.5)  then
   call  conserve(u1(i,j,k,:),0.5197,-0.6259, 0.1, 0.,0.4)    
else if (x1(i,j,k)<0.5.and.y1(i,j,k)<0.5)  then           
   call  conserve(u1(i,j,k,:),0.8,0.1,0.1, 0.,0.4)          
else 
   call  conserve(u1(i,j,k,:),0.5197,0.1 ,-0.6259,0.,0.4)
endif
end do; end do; end do
end subroutine

subroutine readsolution 
use common1 	
implicit none  
call read_solution(u1,i1,j1,k1,'Consequence.plt') 
end subroutine

subroutine read_solution(u,im,jm,km,solutionfile) 
use input;use common1; use input
character cut*1,solutionfile*13
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
OPEN(1,FILE=solutionfile,STATUS='UNKNOWN')
read(1,*) cut
read(1,*) cut
do k=1,km;do j=1,jm ;do i=1,im
read(1,*) xxx,yyy,zzz,qs,ps,us,vs,ws
u(i,j,k,1)=qs
u(i,j,k,2)=qs*us
u(i,j,k,3)=qs*vs
u(i,j,k,4)=qs*ws
u(i,j,k,5)=ps/(gin-1)+0.5*qs*(us**2+vs**2+ws**2)
enddo;enddo;enddo
close(1)
do l=1,5	       
do  k=1,km;do  j=1,jm
do  i=1-mo,0
    u(i,j,k,l)=uo(l)
    enddo
    do  i=im+1,im+mo
    u(i,j,k,l)=uo(l)
    enddo
enddo;enddo
do k=1,km;do  i=1,im
    do  j=1-mo,0
    u(i,j,k,l)=uo(l)
    enddo 
    do  j=jm+1,jm+mo
    u(i,j,k,l)=uo(l)
    enddo
enddo;enddo
do i=1,im;do j=1,jm
    do k=1-mo,0
    u(i,j,k,l)=uo(l)
    enddo
    do k=km+1,km+mo
    u(i,j,k,l)=uo(l)
    enddo
enddo ;enddo 
enddo 
end  subroutine
    
subroutine  conserve(u,qs,us,vs,ws,ps)
use input
implicit none
real u(5),qs,us,vs,ws,ps
u(1)=qs
u(2)=qs*us
u(3)=qs*vs
u(4)=qs*ws
u(5)=ps/(gin-1.)+0.5*(u(2)*us+u(3)*vs+u(4)*ws)
end  subroutine

subroutine qpuvw
use common1  
implicit none
call get_qpuvw(u1,qs1,ps1,us1,vs1,ws1,i1,j1,k1) 
end subroutine

subroutine get_qpuvw(u,qs,ps,us,vs,ws,im,jm,km)
use ghost
implicit none
integer i,j,k,im,jm,km,l
real,external::p
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real qs(im,jm,km),ps(im,jm,km)
real us(im,jm,km),vs(im,jm,km),ws(im,jm,km)	
real ui(5)
do k=1,km;do j=1,jm;do i=1,im
do l=1,5
    ui(l)=u(i,j,k,l)
enddo 
qs(i,j,k)=ui(1)
us(i,j,k)=ui(2)/ui(1)
vs(i,j,k)=ui(3)/ui(1)
ws(i,j,k)=ui(4)/ui(1)
ps(i,j,k)=p(ui)
if(qs(i,j,k)<0) then 
write(*,*) '***********negative rho!!!*********'
stop
elseif(isnan(qs(i,j,k))) then
write(*,*) '***********NAN rho!!!!!*********'
stop
endif
enddo;enddo;enddo
end subroutine



subroutine  writesolution
use common1 
implicit none
call write_solution(x1,y1,z1,qs1,ps1,us1,vs1,ws1,i1,j1,k1,'Consequence.plt') 
end subroutine

subroutine write_solution(xm,ym,zm,qs,ps,us,vs,ws,im,jm,km,solutionfile)
use common1
implicit none
integer i,j,k,im,jm,km,l,num_bpx,num_bpy
character solutionfile*15
real xm(im,jm,km),ym(im,jm,km),zm(im,jm,km)
real qs(im,jm,km),ps(im,jm,km)
real us(im,jm,km),vs(im,jm,km),ws(im,jm,km)
integer sx(im,jm,km),sy(im,jm,km),sxy(im,jm,km)
num_bpx=0
num_bpy=0
open(1,FILE=solutionfile,STATUS='UNKNOWN')
write(1,*) ' variables = "x", "y", "z", "q", "p", "u", "v", "w","sx","sy","sxy","vel_shockx","vel_shocky","vel_shockz","vel_xyz"'
write(1,*) 'zone,i=',im,',j=',jm,',k=',km,',f=point'
do  k=1,km; do  j=1,jm ; do  i=1,im
    if(signal(i,j,k,1) .eqv. .false.)then 
        sx(i,j,k) = 0
        num_bpx = num_bpx + 1
    else
        sx(i,j,k) = 1
    endif
    if(signal(i,j,k,2) .eqv. .false.)then 
        num_bpy = num_bpy + 1
        sy(i,j,k) = 0
    else
        sy(i,j,k) = 1
    endif
    sxy(i,j,k) = min(sx(i,j,k),sy(i,j,k))
    write(1,9) xm(i,j,k),ym(i,j,k),zm(i,j,k),qs(i,j,k),ps(i,j,k),us(i,j,k),vs(i,j,k),ws(i,j,k),&
    sx(i,j,k),sy(i,j,k),sxy(i,j,k),vel_s1(i,j,k),vel_s2(i,j,k),vel_s3(i,j,k), &
    sqrt(vel_s1(i,j,k)**2+vel_s2(i,j,k)**2+vel_s3(i,j,k)**2)
enddo;enddo;enddo
close(1)
write(*,*) real(num_bpx)/real(imax1*jmax1*kmax1),real(num_bpy)/real(imax1*jmax1*kmax1)
9 format(15(g20.8E3,1x))
end subroutine
    