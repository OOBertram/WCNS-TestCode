subroutine scheme_semi_WCNS(fn,u,ex,ey,ez)
use input;use  ghost;use scheme_WCNS
real fn(5),u(-mo1:mo1+1,5),Lu(-mo1:mo1+1)
real uL(5),uR(5),fL(5),fR(5),LL(5,5),RR(5,5),aLR(5),fLR(5)
real LuL(5),LuR(5),Lfn(5),LfLR(5),Ldu(5),a
integer i,j,k
integer,parameter::fluxtype=1
! real,external :: p

uL=u(0,:);uR=u(1,:)
call eigno(fLR,aLR,ll,rr,uL,uR,ex,ey,ez) 

do k=1,5
  do i=-mo1,mo1+1
    Lu(i)=dot_product(LL(k,:),u(i,:))
  enddo
  selectcase(id_scheme_WCNS)
  case(1)
      !write(*,*) "js"
    call ubar_WCNS_E5(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5(LuR(k),Lu(3:-2:-1))
  case(2)
      !write(*,*) "z"
    call ubar_WCNS_E5_Z(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5_Z(LuR(k),Lu(3:-2:-1))
  case(3)
      !write(*,*) "t"
    call ubar_WCNS_E5_T(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5_T(LuR(k),Lu(3:-2:-1))
  case(4)
    call ubar_WCNS_E5_CU6(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5_CU6(LuR(k),Lu(3:-2:-1))
  endselect
enddo
  uL=matmul(RR,LuL);uR=matmul(RR,LuR)
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

subroutine scheme_semi_WCNS_VF(fn,u,ex,ey,ez)
use input;use  ghost;use scheme_WCNS
real fn(5),u(-mo1:mo1+1,5),Lu(-mo1:mo1+1)
real uL(5),uR(5),fL(5),fR(5),LL(5,5),RR(5,5),aLR(5),fLR(5)
real LuL(5),LuR(5),Lfn(5),LfLR(5),Ldu(5),a
integer i,j,k
integer,parameter::fluxtype=1
! real,external :: p

do k=1,5
  do i=-mo1,mo1+1
    Lu(i)=u(i,K)
  enddo
  selectcase(id_scheme_WCNS)
  case(1)
      !write(*,*) "js"
    call ubar_WCNS_E5(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5(LuR(k),Lu(3:-2:-1))
  case(2)
      !write(*,*) "z"
    call ubar_WCNS_E5_Z(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5_Z(LuR(k),Lu(3:-2:-1))
  case(3)
      !write(*,*) "t"
    call ubar_WCNS_E5_T(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5_T(LuR(k),Lu(3:-2:-1))
  case(4)
    call ubar_WCNS_E5_CU6(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
    call ubar_WCNS_E5_CU6(LuR(k),Lu(3:-2:-1))
  endselect
enddo
  uL=LuL;uR=LuR
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

subroutine scheme_semi_WCNS_Zhao(fn,u,ex,ey,ez)
  use input;use  ghost;use scheme_WCNS
  real fn(5),u(-mo1:mo1+1,5),Lu(-mo1:mo1+1)
  real uL(5),uR(5),fL(5),fR(5),LL(5,5),RR(5,5),aLR(5),fLR(5)
  real LuL(5),LuR(5),Lfn(5),LfLR(5),Ldu(5),a
  integer i,j,k
  integer,parameter::fluxtype=1
  ! real,external :: p
  
  uL=u(0,:);uR=u(1,:)
  call eigno(fLR,aLR,ll,rr,uL,uR,ex,ey,ez) 
  
  do k=1,5
    do i=-mo1,mo1+1
      Lu(i)=dot_product(LL(k,:),u(i,:))
    enddo
    selectcase(id_scheme_WCNS)
    case(1)
        !write(*,*) "js"
      call ubar_WCNS_E5(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
      call ubar_WCNS_E5(LuR(k),Lu(3:-2:-1))
    case(2)
        !write(*,*) "z"
      call ubar_WCNS_E5_Z(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
      call ubar_WCNS_E5_Z(LuR(k),Lu(3:-2:-1))
    case(3)
        !write(*,*) "t"
      call ubar_WCNS_E5_T(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
      call ubar_WCNS_E5_T(LuR(k),Lu(3:-2:-1))
    endselect
  enddo
    uL=matmul(RR,LuL);uR=matmul(RR,LuR)
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


subroutine scheme_semi_WCNS_update(fn,u,ex,ey,ez,vel_shock,le)
    use input;use  ghost;use scheme_WCNS
    real fn(5),u(-mo1:mo1+1,5),Lu(-mo1:mo1+1),SL,SR,betaL,betaR,velL,velR,u_,c_,e_,eL,eR
    real uL(5),uR(5),fL(5),fR(5),LL(5,5),RR(5,5),aLR(5),fLR(5)
    real LuL(5),LuR(5),Lfn(5),LfLR(5),Ldu(5),ak,eps,pL,pR,div
    integer i,j,k
    integer,parameter::fluxtype=1
    eps = 0.1
    
    uL=u(0,:);uR=u(1,:)
    call eigno(fLR,aLR,ll,rr,uL,uR,ex,ey,ez) 
    
    do k=1,5
      do i=-mo1,mo1+1
        Lu(i)=dot_product(LL(k,:),u(i,:))
      enddo
      selectcase(id_scheme_WCNS)
      case(1)
        call ubar_WCNS_E5(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
        call ubar_WCNS_E5(LuR(k),Lu(3:-2:-1))
      case(2)
        call ubar_WCNS_E5_Z(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
        call ubar_WCNS_E5_Z(LuR(k),Lu(3:-2:-1))
      case(3)
        call ubar_WCNS_E5_T(LuL(k),Lu(-2:3))  !f(j-1/2)+ 
        call ubar_WCNS_E5_T(LuR(k),Lu(3:-2:-1))
      endselect
    enddo
    
    uL=matmul(RR,LuL);uR=matmul(RR,LuR)
      
    pL = abs(p(uL));pR = abs(p(uR))
    velL = uL(1+le)/uL(1)
    velR = uR(1+le)/uR(1)

    e_L = abs(pL/(uL(1)*0.4))
    e_R = abs(pR/(uR(1)*0.4))

    e_ = (sqrt(uL(1))*e_L+sqrt(uR(1))*e_R) / (SQRT(uL(1))+SQRT(uR(1)))
    u_ = (sqrt(uL(1))*velL+sqrt(uR(1))*velR) / (SQRT(uL(1))+SQRT(uR(1)))
    c_ = sqrt(0.4 * (e_-0.5*u_**2) / sqrt(uL(1)*uR(1)))

    SL = min(velL - sqrt(1.4*pL/uL(1)),u_-c_)
    SR = max(velR + sqrt(1.4*pR/uR(1)),u_+c_)
    betaL = uL(1)*(SL - velL)
    betaR = uR(1)*(SR - velR)
    vel_shock = abs((betaR*velR-betaL*velL+pL-pR)/(betaR-betaL))
    selectcase(fluxtype)
    case(1)
    alpha=max(armm(uL,ex,ey,ez),armm(uR,ex,ey,ez))
    call flux(fL,uL,ex,ey,ez)
    call flux(fR,uR,ex,ey,ez)
    do k=1,5
      fn(k)=0.5*(fL(k)+fR(k)-alpha*(uR(k)-uL(k)))
    enddo
    case(2)
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