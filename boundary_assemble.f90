subroutine jacobi(xa,ya,za,xb,yb,zb,xc,yc,zc,    &
        	       ax,ay,az,bx,by,bz,cx,cy,cz)
xa=by*cz-bz*cy; ya=bz*cx-bx*cz; za=bx*cy-by*cx
xb=cy*az-cz*ay; yb=cz*ax-cx*az; zb=cx*ay-cy*ax
xc=ay*bz-az*by; yc=az*bx-ax*bz; zc=ax*by-ay*bx
endsubroutine

subroutine outerproduct(ab1,ab2,ab3,a1,a2,a3,b1,b2,b3)
ab1=a2*b3-a3*b2
ab2=a3*b1-a1*b3
ab3=a1*b2-a2*b1
endsubroutine

subroutine mirr(w,ax,ay,az,bx,by,bz,cx,cy,cz)
use input   
real rough
dimension w(5)               
data refuze/1./
rough = 1.0
pw=p(w)	
wa=ax*w(2)+ay*w(3)+az*w(4)		      
wc=cx*w(2)+cy*w(3)+cz*w(4)		    
wb=bx*w(2)+by*w(3)+bz*w(4)        
wb=-refuze*wb	                   
wa=rough*wa						              
wc=rough*wc
ajac=(ax*by-ay*bx)*cz+(az*bx-ax*bz)*cy+(ay*bz-az*by)*cx
w(2)=((by*cz-bz*cy)*wa+(cy*az-cz*ay)*wb+(ay*bz-az*by)*wc)/ajac	   !!qu
w(3)=((bz*cx-bx*cz)*wa+(cz*ax-cx*az)*wb+(az*bx-ax*bz)*wc)/ajac	   !!qv
w(4)=((bx*cy-by*cx)*wa+(cx*ay-cy*ax)*wb+(ax*by-ay*bx)*wc)/ajac	   !!qw
w(5)=pw/(gin-1.)+0.5*(w(2)**2+w(3)**2+w(4)**2)/w(1)
endsubroutine

subroutine boundaryk(ukm,le,index)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
integer index
select case(index)
case(1)
  select case(le)
    case(1)
    call boundk_ileft_DM(ukm)
    call boundk_iright_DM(ukm)
    case(2)
    call boundk_jleft_DM(ukm)
    call boundk_jright_DM(ukm)
    case(3)
    call boundk_kleft_DM(ukm)
    call boundk_kright_DM(ukm)
  end select
case(2)
  select case(le)
    case(1)
    call boundk_ileft_Riemann(ukm)
    call boundk_iright_Riemann(ukm)
    case(2)
    call boundk_jleft_Riemann(ukm)
    call boundk_jright_Riemann(ukm)
    case(3)
    call boundk_kleft_Riemann(ukm)
    call boundk_kright_Riemann(ukm)
  end select 
case(3)
  select case(le)
    case(1)
    call boundk_ileft_2DRT(ukm)
    call boundk_iright_2DRT(ukm)
    case(2)
    call boundk_jleft_2DRT(ukm)
    call boundk_jright_2DRT(ukm)
    case(3)
    call boundk_kleft_2DRT(ukm)
    call boundk_kright_2DRT(ukm)
  end select 
case(4)
  select case(le)
    case(1)
    call boundk_ileft_Richmyer(ukm)
    call boundk_iright_Richmyer(ukm)
    case(2)
    call boundk_jleft_Richmyer(ukm)
    call boundk_jright_Richmyer(ukm)
    case(3)
    call boundk_kleft_Richmyer(ukm)
    call boundk_kright_Richmyer(ukm)
  end select
case(5)
  select case(le)
    case(1)
    call boundk_ileft_ShockVortex(ukm)
    call boundk_iright_ShockVortex(ukm)
    case(2)
    call boundk_jleft_ShockVortex(ukm)
    call boundk_jright_ShockVortex(ukm)
    case(3)
    call boundk_kleft_ShockVortex(ukm)
    call boundk_kright_ShockVortex(ukm)
  end select
case(6)
  select case(le)
    case(1)
    call boundk_ileft_Riemann(ukm)
    call boundk_iright_Riemann(ukm)
    case(2)
    call boundk_jleft_Riemann(ukm)
    call boundk_jright_Riemann(ukm)
    case(3)
    call boundk_kleft_Riemann(ukm)
    call boundk_kright_Riemann(ukm)
  end select
case(7)
  select case(le)
    case(1)
    call boundk_ileft_Riemann(ukm)
    call boundk_iright_Riemann(ukm)
    case(2)
    call boundk_jleft_Riemann(ukm)
    call boundk_jright_Riemann(ukm)
    case(3)
    call boundk_kleft_Riemann(ukm)
    call boundk_kright_Riemann(ukm)
  end select 
case(8)
  select case(le)
    case(1)
    call boundk_ileft_Riemann(ukm)
    call boundk_iright_Riemann(ukm)
    case(2)
    call boundk_jleft_Riemann(ukm)
    call boundk_jright_Riemann(ukm)
    case(3)
    call boundk_kleft_Riemann(ukm)
    call boundk_kright_Riemann(ukm)
  end select 
end select 
endsubroutine
    

subroutine boundk_ileft_DM(ukm)   
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
do k=1,k1 ; do j=1,j1 ; do i=1-mo,0
  call  conserve(ukm(i,j,k,:),8.0,7.14471,-4.125,0.,116.5)        !i����߽�
enddo; enddo; enddo 
endsubroutine

subroutine boundk_iright_DM(ukm)   
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)  
  call bound_out1_i_right(ukm,i1,j1,k1) !i����߽�
endsubroutine

subroutine boundk_jleft_DM(ukm)                  
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_wall_j_left(ukm,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,i1,j1,k1,1,i1,1,k1)
do k=1,k1; do j=1-mo,0; do i=1,i1
  if(x1(i,1,1)<1./6.)   call conserve(ukm(i,j,k,:),8.0,7.14471,-4.125,0.,116.5)   
enddo ; enddo ; enddo 
endsubroutine

subroutine boundk_jright_DM(ukm)
use common1;use ghost;use input
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
ss=10.
do k=1,k1; do j=j1+1,j1+mo ; do i=1,i1 !j����߽�
if(x1(i,j1,k1)<=(1./6. + sqrt(3.)/3*(1.+2.*ss*timenow(2))))  then  
  call conserve(ukm(i,j,k,:),8.0,7.14471 ,-4.125 ,0. ,116.5 )
else
  call conserve(ukm(i,j,k,:),1.4  ,0. ,0. ,0.,1.0)
endif
end do ; end do; end do  
endsubroutine

subroutine boundk_kleft_DM(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_left(ukm,i1,j1,k1)                   
endsubroutine

subroutine boundk_kright_DM(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_right(ukm,i1,j1,k1)                 
endsubroutine

subroutine boundk_ileft_Riemann(ukm)     
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_i_left(ukm,i1,j1,k1) 
endsubroutine
 
subroutine boundk_iright_Riemann(ukm) 
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_i_right(ukm,i1,j1,k1) !i����߽�
endsubroutine
 
subroutine boundk_jleft_Riemann(ukm)                  
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_j_left(ukm,i1,j1,k1)
!call bound_1_period_j_left(uk1,i1,j1,k1)
endsubroutine

subroutine boundk_jright_Riemann(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_j_right(ukm,i1,j1,k1) !j����߽�
endsubroutine

subroutine boundk_kleft_Riemann(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_left(ukm,i1,j1,k1)   !k����������                    
endsubroutine

subroutine boundk_kright_Riemann(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_right(ukm,i1,j1,k1)  !k����������                    
endsubroutine

subroutine boundk_ileft_2DRT(ukm)     !i������߽�,��������
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_wall_i_left(ukm,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,    &
                                i1,j1,k1,1,j1,1,k1)   !i����������
endsubroutine

subroutine boundk_iright_2DRT(ukm)    !i�����ұ߽�,��������
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_wall_i_right(ukm,ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,    &
                                i1,j1,k1,1,j1,1,k1) !i����������
endsubroutine

subroutine boundk_jleft_2DRT(ukm)                  
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
do k=1,k1; do j=1-mo,0; do i=1,i1     
  call conserve(ukm(i,j,k,:),2.0,0.,0.,0.,1.0)
enddo;enddo; enddo
endsubroutine

subroutine boundk_jright_2DRT(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
do k=1,k1; do j=j1+1,j1+mo; do i=1,i1       
  call conserve(ukm(i,j,k,:),1.,0.,0.,0.,2.5)
enddo;enddo; enddo
endsubroutine

subroutine boundk_kleft_2DRT(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_left(ukm,i1,j1,k1)   !k����������                    
endsubroutine

subroutine boundk_kright_2DRT(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_right(ukm,i1,j1,k1)  !k����������                    
endsubroutine

subroutine boundk_ileft_Richmyer(ukm)     !i������߽�,��������
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_i_left(ukm,i1,j1,k1) !i����߽�
endsubroutine
 
subroutine boundk_iright_Richmyer(ukm)    !i�����ұ߽�,��������
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_i_right(ukm,i1,j1,k1) !i����߽�
endsubroutine
 


subroutine boundk_jleft_Richmyer(ukm)                  
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_j_left(ukm,i1,j1,k1)
!call bound_1_period_j_left(uk1,i1,j1,k1)
endsubroutine

subroutine boundk_jright_Richmyer(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_j_right(ukm,i1,j1,k1) !j����߽�
endsubroutine

subroutine boundk_kleft_Richmyer(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_left(ukm,i1,j1,k1)   !k����������                    
endsubroutine

subroutine boundk_kright_Richmyer(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_right(ukm,i1,j1,k1)  !k����������                    
endsubroutine

subroutine boundk_ileft_ShockVortex(ukm)     !i������߽�,��������
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_i_left(ukm,i1,j1,k1) !i����߽�
endsubroutine
 
subroutine boundk_iright_ShockVortex(ukm)    !i�����ұ߽�,��������
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_i_right(ukm,i1,j1,k1) !i����߽�
endsubroutine
 
subroutine boundk_jleft_ShockVortex(ukm)                  
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_j_left(ukm,i1,j1,k1)
!call bound_1_period_j_left(uk1,i1,j1,k1)
endsubroutine

subroutine boundk_jright_ShockVortex(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_out1_j_right(ukm,i1,j1,k1) !j����߽�
endsubroutine

subroutine boundk_kleft_ShockVortex(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_left(ukm,i1,j1,k1)   !k����������                    
endsubroutine

subroutine boundk_kright_ShockVortex(ukm)
use common1;use ghost
real ukm(1-mo:i1+mo,1-mo:j1+mo,1-mo:k1+mo,5)
call bound_1_period_k_right(ukm,i1,j1,k1)  !k����������                    
endsubroutine

  subroutine bound_out1_i_left(u,im,jm,km)  
  use ghost
  real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
  do k=1,km ; do j=1,jm ;do i=0,mo-1
  do l=1,5
   u(-i,j,k,l)=u(-i+1,j,k,l)
  enddo
  enddo;enddo;enddo
  endsubroutine
  
  subroutine bound_out1_i_right(u,im,jm,km) 
  use ghost
  real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
  do k=1,km ;do j=1,jm ; do i=im+1,im+mo
  do l=1,5
   u(i,j,k,l)=u(i-1,j,k,l)
  enddo
  enddo;enddo;enddo
  endsubroutine
  
  subroutine  bound_out1_j_left(u,im,jm,km)   ! j��������
  use ghost
  real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
  do k=1,km;do i=1,im;do j=0,mo-1
  do l=1,5
    u(i,-j,k,l)=u(i,-j+1,k,l)
    ! u(i,-j,k,l)=u(i,1,k,l)
  enddo
  enddo;enddo;enddo
  endsubroutine
  
  subroutine  bound_out1_j_right(u,im,jm,km)   ! j�ҷ�������
  use ghost
  real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
  do k=1,km;do i=1,im;do j=1,mo
  do l=1,5
     u(i,jm+j,k,l)=u(i,jm+j-1,k,l)
    ! u(i,jm+j,k,l)=u(i,jm,k,l)
  enddo
  enddo;enddo;enddo
  endsubroutine

subroutine bound_wall_i_left(u,ax,ay,az,bx,by,bz,cx,cy,cz,     &
                            	im,jm,km,jl,jr,kl,kr)  !i����������
use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real  ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real  bx(im,jm,km),by(im,jm,km),bz(im,jm,km) 
real  cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real ww(5)	  
do k=kl,kr
do j=jl,jr
axs=ax(1,j,k);ays=ay(1,j,k);azs=az(1,j,k)
cxs=cx(1,j,k);cys=cy(1,j,k);czs=cz(1,j,k)
bxs=bx(1,j,k);bys=by(1,j,k);bzs=bz(1,j,k)
do i=1-mo,0
do l=1,5
 ww(l)=u(1-i,j,k,l)
enddo
call jacobi(xa,ya,za,xb,yb,zb,xc,yc,zc,		  &       !!����������ϵ�i��x��ƫ�����x��i��ƫ��
            axs,ays,azs,bxs,bys,bzs,cxs,cys,czs)
call outerproduct(xa,ya,za,xb,yb,zb,xc,yc,zc)
call mirr(ww,xb,yb,zb,xa,ya,za,xc,yc,zc)	!����ֱ
! call mirr(ww,bxs,bys,bzs,axs,ays,azs,cxs,cys,czs)	!��ֱ
do l=1,5			
u(i,j,k,l)=ww(l)	    !!���ú���Ҫ��(i=-1,i=0�Ĵ���ֵ)
enddo
enddo;enddo;enddo
endsubroutine

subroutine bound_wall_i_right(u,ax,ay,az,bx,by,bz,cx,cy,cz,    &
                                im,jm,km,jl,jr,kl,kr) !i����������

use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real  ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real  bx(im,jm,km),by(im,jm,km),bz(im,jm,km) 
real  cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real ww(5)	  
do k=kl,kr
do j=jl,jr
axs=ax(im,j,k);ays=ay(im,j,k);azs=az(im,j,k)
cxs=cx(im,j,k);cys=cy(im,j,k);czs=cz(im,j,k)
bxs=bx(im,j,k);bys=by(im,j,k);bzs=bz(im,j,k)
do i=im+1,im+mo
do l=1,5
 ww(l)=u(2*im+1-i,j,k,l)
enddo
call jacobi(xa,ya,za,xb,yb,zb,xc,yc,zc,   &
             axs,ays,azs,bxs,bys,bzs,cxs,cys,czs)
call outerproduct(xa,ya,za,xb,yb,zb,xc,yc,zc)
call mirr(ww,xb,yb,zb,xa,ya,za,xc,yc,zc)	!����ֱ
!	call mirr(ww,bxs,bys,bzs,axs,ays,azs,cxs,cys,czs)	!��ֱ
do l=1,5			
u(i,j,k,l)=ww(l)
enddo
enddo;enddo;enddo
endsubroutine

subroutine bound_wall_j_left(u,ax,ay,az,bx,by,bz,cx,cy,cz,    &
                           	im,jm,km,il,ir,kl,kr) !j����������
use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
real  ax(im,jm,km),ay(im,jm,km),az(im,jm,km)
real  bx(im,jm,km),by(im,jm,km),bz(im,jm,km) 
real  cx(im,jm,km),cy(im,jm,km),cz(im,jm,km)
real ww(5)	  
do k=kl,kr
do i=il,ir
axs=ax(i,1,k);ays=ay(i,1,k);azs=az(i,1,k)
cxs=cx(i,1,k);cys=cy(i,1,k);czs=cz(i,1,k)
bxs=bx(i,1,k);bys=by(i,1,k);bzs=bz(i,1,k)
do j=1-mo,0
do l=1,5
    ww(l)=u(i,1-j,k,l)
enddo
call jacobi(xa,ya,za,xb,yb,zb,xc,yc,zc,    &
            axs,ays,azs,bxs,bys,bzs,cxs,cys,czs)
call outerproduct(xb,yb,zb,xa,ya,za,xc,yc,zc)
call mirr(ww,xa,ya,za,xb,yb,zb,xc,yc,zc)	!����ֱ
!	call mirr(ww,axs,ays,azs,bxs,bys,bzs,cxs,cys,czs)	!��ֱ
do l=1,5			
 u(i,j,k,l)=ww(l)
enddo
enddo;enddo;enddo
endsubroutine

subroutine bound_1_period_k_left(u,im,jm,km)	!k��һ�ܻ�����������
use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
do i=1,im
do j=1,jm
do k=0,mo-1
do l=1,5
  u(i,j,-k,l)=u(i,j,km-k,l)
enddo
enddo;enddo;enddo
endsubroutine

subroutine bound_1_period_k_right(u,im,jm,km) !k��һ�ܻ�������
use ghost
real u(1-mo:im+mo,1-mo:jm+mo,1-mo:km+mo,5)
do i=1,im
do j=1,jm
do k=km+1,km+mo
do l=1,5
 u(i,j,k,l)=u(i,j,k-km,l)
enddo
enddo;enddo;enddo
endsubroutine
