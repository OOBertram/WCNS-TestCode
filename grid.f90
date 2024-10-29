subroutine  gridassemble(index)
use common1; use input
implicit none 
real xh,yh,zh
integer i,j,k
integer index
select case(index)
case(1)
   !Double-Mach
   ! imaxo1=961;jmaxo1=241;kmaxo1=2;
   ! imaxo1=481;jmaxo1=121;kmaxo1=2;
   ! imaxo1=241;jmaxo1=61;kmaxo1=2;
   imaxo1=201;jmaxo1=51;kmaxo1=2;
   xmi=0.;xma=4.;ymi=0.;yma=1.;zmi=-0.5;zma=0.5;
   timestop = 0.2
case(2)
   !2DRiemann
   !imaxo1=451;jmaxo1=451;kmaxo1=2;
   ! imaxo1=401;jmaxo1=401;kmaxo1=2;
   ! imaxo1=201;jmaxo1=201;kmaxo1=2;
   imaxo1=101;jmaxo1=101;kmaxo1=2;
   xmi=0.;xma=1.;ymi=0.;yma=1.;zmi=-0.5;zma=0.5;
   timestop = 0.8
case(3)
   !2DRT
   ! imaxo1=241;jmaxo1=961;kmaxo1=2;
   ! imaxo1=121;jmaxo1=481;kmaxo1=2;
   imaxo1=51;jmaxo1=201;kmaxo1=2;
   xmi=0.;xma=0.25;ymi=0.;yma=1.;zmi=-0.5;zma=0.5;
   timestop = 1.95
case(4)
   !2DRichmyer
   ! imaxo1=301;jmaxo1=251;kmaxo1=2;
   imaxo1=601;jmaxo1=501;kmaxo1=2;
   ! imaxo1=151;jmaxo1=126;kmaxo1=2;
   xmi=0.;xma=6.;ymi=0.;yma=5.;zmi=-0.5;zma=0.5;
   timestop = 8.0
case(5)
   !Shock_vortex
   imaxo1=501;jmaxo1=201;kmaxo1=2
   ! imaxo1=251;jmaxo1=101;kmaxo1=2;
   ! imaxo1=201;jmaxo1=81;kmaxo1=2;
   xmi=0.;xma=2.;ymi=0.;yma=1.;zmi=-0.5;zma=0.5;
   timestop = 0.6
case(6)
   !2DIRV
   imaxo1=451;jmaxo1=451;kmaxo1=2;
   ! imaxo1=401;jmaxo1=401;kmaxo1=2;
   ! imaxo1=201;jmaxo1=201;kmaxo1=2;
   xmi=0.;xma=1.;ymi=0.;yma=1.;zmi=-0.5;zma=0.5;
   timestop = 0.3
case(7)
   !2DRiemann2
   imaxo1=501;jmaxo1=501;kmaxo1=2;
   ! imaxo1=401;jmaxo1=401;kmaxo1=2;
   ! imaxo1=201;jmaxo1=201;kmaxo1=2;
   ! imaxo1=126;jmaxo1=126;kmaxo1=2;
   xmi=0.;xma=1.;ymi=0.;yma=1.;zmi=-0.5;zma=0.5;
   timestop = 0.3
case(8)
   !2DRiemann3
   !imaxo1=451;jmaxo1=451;kmaxo1=2;
   ! imaxo1=401;jmaxo1=401;kmaxo1=2;
   imaxo1=201;jmaxo1=201;kmaxo1=2;
   xmi=0.;xma=1.;ymi=0.;yma=1.;zmi=-0.5;zma=0.5;
   timestop = 0.25
end select
imax1=imaxo1-1;jmax1=jmaxo1-1;kmax1=kmaxo1-1;
allocate(xo1(imaxo1,jmaxo1,kmaxo1),yo1(imaxo1,jmaxo1,kmaxo1),zo1(imaxo1,jmaxo1,kmaxo1))
allocate(x1(imax1,jmax1,kmax1),y1(imax1,jmax1,kmax1),z1(imax1,jmax1,kmax1))
allocate(ax1(imax1,jmax1,kmax1),ay1(imax1,jmax1,kmax1),az1(imax1,jmax1,kmax1))
allocate(bx1(imax1,jmax1,kmax1),by1(imax1,jmax1,kmax1),bz1(imax1,jmax1,kmax1))
allocate(cx1(imax1,jmax1,kmax1),cy1(imax1,jmax1,kmax1),cz1(imax1,jmax1,kmax1))
allocate(jac1(imax1,jmax1,kmax1)) 
allocate(qs1(imax1,jmax1,kmax1),ps1(imax1,jmax1,kmax1),us1(imax1,jmax1,kmax1),vs1(imax1,jmax1,kmax1),ws1(imax1,jmax1,kmax1))
allocate(u1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,5),uk1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,5)) 
allocate(uk2(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,5),uk3(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,5)) 
allocate(uk(11,1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,5)) 
allocate(sm1(1-mo:imax1+mo-1,1-mo:jmax1+mo-1,1-mo:kmax1+mo-1))
allocate(sm2(1-mo:imax1+mo-1,1-mo:jmax1+mo-1,1-mo:kmax1+mo-1))
allocate(sm3(1-mo:imax1+mo-1,1-mo:jmax1+mo-1,1-mo:kmax1+mo-1))
allocate(signal(2-mo:imax1+mo-3,2-mo:jmax1+mo-3,2-mo:kmax1+mo-3,3)) 
allocate(signal2(2-mo:imax1+mo-3,2-mo:jmax1+mo-3,2-mo:kmax1+mo-3,3)) 
allocate(vel_s1(3-mo:imax1+mo-3,3-mo:jmax1+mo-3,3-mo:kmax1+mo-3)) 
allocate(vel_s2(3-mo:imax1+mo-3,3-mo:jmax1+mo-3,3-mo:kmax1+mo-3)) 
allocate(vel_s3(3-mo:imax1+mo-3,3-mo:jmax1+mo-3,3-mo:kmax1+mo-3)) 
allocate(Averagerho(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo))
allocate(Averagerho_p(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo))
allocate(AverageE(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo))
allocate(smLiu1(3-mo:imax1+mo-3,3-mo:jmax1+mo-3,-mo:kmax1+mo-3))
allocate(smLiu2(3-mo:imax1+mo-3,3-mo:jmax1+mo-3,3-mo:kmax1+mo-3))
allocate(smLiu3(3-mo:imax1+mo-3,3-mo:jmax1+mo-3,3-mo:kmax1+mo-3))
xh=(xma-xmi)/(imaxo1-1);yh=(yma-ymi)/(jmaxo1-1);zh=(zma-zmi)/(kmaxo1-1)
do k=1,kmaxo1; do j=1,jmaxo1; do i=1,imaxo1
    xo1(i,j,k)=xmi+(i-1)*xh 
    yo1(i,j,k)=ymi+(j-1)*yh 
    zo1(i,j,k)=zmi+(k-1)*zh 
end do ;end do;end do
open(1,file='Initial.plt') 
write(1,*) ' variables = "x", "y", "z" '
write(1,*) 'zone,i=',imaxo1,',j=',jmaxo1,',k=',kmaxo1,',f=point'
do k=1,kmaxo1; do j=1,jmaxo1; do i=1,imaxo1
    write(1,*) xo1(i,j,k),yo1(i,j,k),zo1(i,j,k)
end do ;end do;end do
i1=imaxo1-1;j1=jmaxo1-1;k1=kmaxo1-1
end subroutine

subroutine midgrid
use common1
implicit none 
call middle_grid(xo1,yo1,zo1,x1,y1,z1,              &
ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,jac1,i1,j1,k1)                
end subroutine

subroutine middle_grid(x,y,z,xm,ym,zm,ax,ay,az,bx,by,bz,cx,cy,cz,jac,im,jm,km)                  	
    real x(im+1,jm+1,km+1),y(im+1,jm+1,km+1),z(im+1,jm+1,km+1)
    real xm(im,jm,km),ym(im,jm,km),zm(im,jm,km)
    real ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
    real bx(im,jm,km),by(im,jm,km),bz(im,jm,km)
    real cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
    real jac(im,jm,km)
    do  i=1,im; do  j=1,jm ; do  k=1,km
    xm(i,j,k)=(x(i,j,k)+x(i+1,j,k)+x(i,j+1,k)+x(i+1,j+1,k)+		&	
               x(i,j,k+1)+x(i+1,j,k+1)+x(i,j+1,k+1)+x(i+1,j+1,k+1))/8.
    ym(i,j,k)=(y(i,j,k)+y(i+1,j,k)+y(i,j+1,k)+y(i+1,j+1,k)+     &
               y(i,j,k+1)+y(i+1,j,k+1)+y(i,j+1,k+1)+y(i+1,j+1,k+1))/8.
    zm(i,j,k)=(z(i,j,k)+z(i+1,j,k)+z(i,j+1,k)+z(i+1,j+1,k)+     &
               z(i,j,k+1)+z(i+1,j,k+1)+z(i,j+1,k+1)+z(i+1,j+1,k+1))/8.
    xa=(x(i+1,j,k)+x(i+1,j+1,k)+x(i+1,j,k+1)+x(i+1,j+1,k+1)-    &	
       (x(i,j,k)+x(i,j+1,k)+x(i,j,k+1)+x(i,j+1,k+1)))/4.
    xb=(x(i,j+1,k)+x(i+1,j+1,k)+x(i,j+1,k+1)+x(i+1,j+1,k+1)-    & 
       (x(i,j,k)+x(i+1,j,k)+x(i,j,k+1)+x(i+1,j,k+1)))/4.
    xc=(x(i,j,k+1)+x(i,j+1,k+1)+x(i+1,j,k+1)+x(i+1,j+1,k+1)-    &   
       (x(i,j,k)+x(i,j+1,k)+x(i+1,j,k)+x(i+1,j+1,k)))/4.
    ya=(y(i+1,j,k)+y(i+1,j+1,k)+y(i+1,j,k+1)+y(i+1,j+1,k+1)-    &
       (y(i,j,k)+y(i,j+1,k)+y(i,j,k+1)+y(i,j+1,k+1)))/4.
    yb=(y(i,j+1,k)+y(i+1,j+1,k)+y(i,j+1,k+1)+y(i+1,j+1,k+1)-    &
       (y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1)))/4.
    yc=(y(i,j,k+1)+y(i,j+1,k+1)+y(i+1,j,k+1)+y(i+1,j+1,k+1)-    &
       (y(i,j,k)+y(i,j+1, k)+y(i+1,j,k)+y(i+1,j+1,k)))/4.
    za=(z(i+1,j,k)+z(i+1,j+1,k)+z(i+1,j,k+1)+z(i+1,j+1,k+1)-    &
       (z(i,j,k)+z(i,j+1,k)+z(i,j,k+1)+z(i,j+1,k+1)))/4.
    zb=(z(i,j+1,k)+z(i+1,j+1,k)+z(i,j+1,k+1)+z(i+1,j+1,k+1)-    &
       (z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1)))/4.
    zc=(z(i,j,k+1)+z(i,j+1,k+1)+z(i+1,j,k+1)+z(i+1,j+1,k+1)-    &
       (z(i,j,k)+z(i,j+1,k)+z(i+1,j,k)+z(i+1,j+1,k)))/4.
    jac(i,j,k)=(xa*yb-ya*xb)*zc+(xc*ya-yc*xa)*zb+(xb*yc-yb*xc)*za		
    ax(i,j,k)=(yb*zc-zb*yc)/jac(i,j,k)   
    ay(i,j,k)=(zb*xc-xb*zc)/jac(i,j,k)	
    az(i,j,k)=(xb*yc-yb*xc)/jac(i,j,k)	
    bx(i,j,k)=(yc*za-zc*ya)/jac(i,j,k)
    by(i,j,k)=(zc*xa-xc*za)/jac(i,j,k)
    bz(i,j,k)=(xc*ya-yc*xa)/jac(i,j,k)
    cx(i,j,k)=(ya*zb-za*yb)/jac(i,j,k)
    cy(i,j,k)=(za*xb-xa*zb)/jac(i,j,k)
    cz(i,j,k)=(xa*yb-ya*xb)/jac(i,j,k)
    enddo;enddo;enddo
    end subroutine