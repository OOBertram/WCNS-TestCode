real function p(u)
use input
implicit none
real u(5),g1
g1=gin-1.
p=g1*(u(5)-0.5*(u(2)**2+u(3)**2+u(4)**2)/u(1))
endfunction

real function uo(k)
use input
implicit none 
real u(5) 
integer k
uin=sm*sqrt(gin*pin/qin)
u(1)=qin
u(2)=qin*uin	           
u(3)=0.
u(4)=0.
u(5)=pin/(gin-1)+0.5*u(2)**2/qin
uo=u(k)
end function

subroutine flux(f,u,ex,ey,ez)	      !ex,ey,ez�����ͨ��
implicit none
real ex,ey,ez
real f(5),u(5),ps
real,external::p
ps=p(u)
f(1)=ex*u(2)+ey*u(3)+ez*u(4)
f(2)=f(1)*u(2)/u(1)+ex*ps
f(3)=f(1)*u(3)/u(1)+ey*ps
f(4)=f(1)*u(4)/u(1)+ez*ps
f(5)=f(1)*(u(5)+ps)/u(1)
end subroutine


real function H(u)	 !��
implicit none
real u(5)
real,external::p
H=(u(5)+p(u))/u(1)
endfunction

function timenowminimum(timenow,lem)	
real timenow(3)
timenowminimum=timenow(lem)
do le=1,3
  if(timenow(le).le.timenowminimum) then
    lem=le
    timenowminimum=timenow(lem)
  endif
enddo
  end function

function timenowminimum1(timenow,lem)	
real timenow(4)
timenowminimum=timenow(lem)
do le=1,4
  if(timenow(le).le.timenowminimum) then
    lem=le
    timenowminimum=timenow(lem)
  endif
enddo
  end function