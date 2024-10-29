subroutine eigno(ff,aa,ll,rr,u1,u2,ex,ey,ez)
use input
real ff(5),aa(5),ll(5,5),rr(5,5),u1(5),u2(5),u12(5),f1(5),f2(5)
real g1
real cseeeeee,beta
beta=1E-10
g1=gin-1.0
us1=u1(2)/u1(1);vs1=u1(3)/u1(1);ws1=u1(4)/u1(1)
us2=u2(2)/u2(1);vs2=u2(3)/u2(1);ws2=u2(4)/u2(1)
uevewe1=ex*us1+ey*vs1+ez*ws1 
uevewe2=ex*us2+ey*vs2+ez*ws2 
if(uevewe1.gt.uevewe2) then	
   sq1=sqrt(max(u1(1),beta));   sq2=sqrt(max(u2(1),beta))
   sq12=sq1+sq2
   us=(sq1*us1+sq2*us2)/sq12
   vs=(sq1*vs1+sq2*vs2)/sq12
   ws=(sq1*ws1+sq2*ws2)/sq12
   Hs=(sq1*H(u1)+sq2*H(u2))/sq12
   call flux(f1,u1,ex,ey,ez)
   call flux(f2,u2,ex,ey,ez)
   ff=(f1+f2)/2.
else 
   u12=(u1+u2)/2.
   us=u12(2)/u12(1); vs=u12(3)/u12(1); ws=u12(4)/u12(1)
   Hs=H(u12)
   call flux(ff,u12,ex,ey,ez)
endif
eeeeee=sqrt(ex**2+ey**2+ez**2)
uevewe=ex*us+ey*vs+ez*ws
uuvvww=0.5*(us**2+vs**2+ws**2)
!cs=sqrt(abs(g1*(Hs-uuvvww)))
cs=sqrt(max(g1*(Hs-uuvvww),beta))
aa(1)=uevewe-eeeeee*cs
aa(2)=uevewe
aa(3)=uevewe+eeeeee*cs
aa(4)=uevewe
aa(5)=uevewe
csgeeeeee=cs/(g1*eeeeee)
ll(1,1)=uuvvww+uevewe*csgeeeeee;ll(1,2)=-us-ex*csgeeeeee;ll(1,3)=-vs-ey*csgeeeeee;ll(1,4)=-ws-ez*csgeeeeee;ll(1,5)=1.
ll(2,1)=uuvvww-cs*cs/g1;        ll(2,2)=-us;             ll(2,3)=-vs;             ll(2,4)=-ws;             ll(2,5)=1.
ll(3,1)=uuvvww-uevewe*csgeeeeee;ll(3,2)=-us+ex*csgeeeeee;ll(3,3)=-vs+ey*csgeeeeee;ll(3,4)=-ws+ez*csgeeeeee;ll(3,5)=1.
ll(4,5)=0.;ll(5,5)=0.
cseeeeee=cs/eeeeee;ee1=1./(eeeeee*eeeeee)
rr(1,1)=1.;                                  rr(1,2)=-2.;      rr(1,3)=1.
rr(2,1)=us-ex*cseeeeee;                      rr(2,2)=-2*us;    rr(2,3)=us+ex*cseeeeee 
rr(3,1)=vs-ey*cseeeeee;                      rr(3,2)=-2*vs;    rr(3,3)=vs+ey*cseeeeee
rr(4,1)=ws-ez*cseeeeee;                      rr(4,2)=-2*ws;    rr(4,3)=ws+ez*cseeeeee
rr(5,1)=uuvvww+cs*cs/(gin-1.)-uevewe*cseeeeee;rr(5,2)=-2*uuvvww;rr(5,3)=uuvvww+cs*cs/(gin-1.)+uevewe*cseeeeee
rr(1,4)=0.;rr(1,5)=0.
if(abs(ex).ge.abs(ey).and.abs(ex).ge.abs(ez)) then
  ll(4,2)=-ey;ll(4,3)=ex;ll(4,4)=0.
  ll(5,2)=-ez;ll(5,3)=0.;ll(5,4)=ex
  rr(2,4)=-ey*ee1;             rr(2,5)=-ez*ee1
  rr(3,4)=(ex*ex+ez*ez)*ee1/ex;rr(3,5)=-ey*ez*ee1/ex
  rr(4,4)=-ey*ez*ee1/ex;       rr(4,5)=(ex*ex+ey*ey)*ee1/ex 
else 
if(abs(ey).ge.abs(ez).and.abs(ey).ge.abs(ex)) then
  ll(4,2)=0.; ll(4,3)=-ez;ll(4,4)=ey
  ll(5,2)=ey; ll(5,3)=-ex;ll(5,4)=0.
  rr(2,4)=-ex*ez*ee1/ey;       rr(2,5)=(ey*ey+ez*ez)*ee1/ey
  rr(3,4)=-ez*ee1;             rr(3,5)=-ex*ee1
  rr(4,4)=(ex*ex+ey*ey)*ee1/ey;rr(4,5)=-ex*ez*ee1/ey 
else
  ll(4,2)=ez; ll(4,3)=0.; ll(4,4)=-ex
  ll(5,2)=0.; ll(5,3)=ez; ll(5,4)=-ey
  rr(2,4)=(ey*ey+ez*ez)*ee1/ez;rr(2,5)=-ex*ey*ee1/ez
  rr(3,4)=-ex*ey*ee1/ez;       rr(3,5)=(ex*ex+ez*ez)*ee1/ez
  rr(4,4)=-ex*ee1;             rr(4,5)=-ey*ee1
endif
endif
ll(4,1)=-ll(4,2)*us-ll(4,3)*vs-ll(4,4)*ws
ll(5,1)=-ll(5,2)*us-ll(5,3)*vs-ll(5,4)*ws
rr(5,4)=us*rr(2,4)+vs*rr(3,4)+ws*rr(4,4)
rr(5,5)=us*rr(2,5)+vs*rr(3,5)+ws*rr(4,5)
rr(1:5,1:3)=rr(1:5,1:3)*(gin-1)/(2*cs*cs)
  end subroutine 