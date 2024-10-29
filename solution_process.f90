subroutine Global_Supervise(ukm,number,le,timeuse)
 use common1;use scheme_WCNS
implicit none
real :: ukm(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5),sumrho,sump,sumE
real :: timebg,timend,timeuse
integer :: number,le
  sumrho = 0.0
  sump = 0.0
  sumE = 0.0
  call CPU_TIME(timebg)
  select case(mix_or_not)
  case(1)
    call Average(ukm(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5),sumrho,sump,sumE)
    call Cal_Smoothness(ukm(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5),le,sumrho,sump,sumE)
    call Judge(number,le)
  case(2)
    call Cal_Smoothness_Liu(ukm(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5),le)
    call Judge_Liu(number,le)
  case(3)
    call Cal_Smoothness_TOT(ukm(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5),le)
  end select
  call CPU_TIME(timend)
  timeuse = timeuse + timend - timebg
endsubroutine

subroutine Average(u,sumrho,sump,sumE)
use common1;use ghost
real u(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5),sumrho,sump,pa,sumE
integer :: i,j,k
real,external ::p
do k = 1,kmax1;do j = 1,jmax1;do i = 1,imax1
  sumrho = sumrho + u(i,j,k,1)
  pa = p(u(i,j,k,:))
  sump = sump + pa
  sumE = sumE + u(i,j,k,5)
enddo; enddo; enddo
  sumrho = sumrho / (imax1 * jmax1 * kmax1)
  sump = sump / (imax1 * jmax1 * kmax1)
  sumE = sumE / (imax1 * jmax1 * kmax1)
endsubroutine

subroutine Cal_Smoothness(u,le,sumrho,sump,sumE)
use common1
implicit none
real u(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5)
real Qp(-1:1),sump,sumrho,sumprho,sumE,sumprho_inv,sumE_inv,sumrho_inv
real,external:: p
integer :: i,j,k,q,le
  sumrho_inv = 1.0/sumrho
  sumprho_inv = sumrho / sump
  sumE_inv = 1.0 / sumE
  do k=1-mo,kmax1+mo ;do  j=1-mo,jmax1+mo; do  i=1-mo,imax1+mo
    Averagerho(i,j,k) = u(i,j,k,1) * sumrho_inv
    Averagerho_p(i,j,k) = p(u(i,j,k,:)) * sumrho_inv * sumprho_inv
    Averagerho_p(i,j,k) = Averagerho(i,j,k) / Averagerho_p(i,j,k)
    AverageE(i,j,k) = u(i,j,k,5) * sumE_inv
  enddo;enddo;enddo
  select case(le)
    case(1)
    do k=1-mo,kmax1+mo ;do  j=1-mo,jmax1+mo; do  i=2-mo,imax1+mo-1
      Qp(-1) = Averagerho_p(i-1,j,k) + Averagerho(i-1,j,k) + AverageE(i-1,j,k)
      Qp(0)  = Averagerho_p(i,j,k) + Averagerho(i,j,k) + AverageE(i,j,k)
      Qp(1)  = Averagerho_p(i+1,j,k) + Averagerho(i+1,j,k) + AverageE(i+1,j,k)
      sm1(i,j,k) = (Qp(-1)-2.0D0*Qp(0)+Qp(1))**2+0.25D0*(Qp(-1)-Qp(1))**2;
      ! write(*,*) "FT"
    enddo;enddo;enddo

    case(2)
    do k=1-mo,kmax1+mo ;do  j=2-mo,jmax1+mo-1; do  i=1-mo,imax1+mo
      Qp(-1)  = Averagerho_p(i,j-1,k) + Averagerho(i,j-1,k) + AverageE(i,j-1,k)
      Qp(0 )  = Averagerho_p(i,j,k) + Averagerho(i,j,k) + AverageE(i,j,k)
      Qp(1 )  = Averagerho_p(i,j+1,k) + Averagerho(i,j+1,k) + AverageE(i,j+1,k)
      sm2(i,j,k) = (Qp(-1)-2.0D0*Qp(0)+Qp(1))**2+0.25D0*(Qp(-1)-Qp(1))**2;
    enddo;enddo;enddo
    
    case(3)
    do k=2-mo,kmax1+mo-1 ;do  j=1-mo,jmax1+mo; do  i=1-mo,imax1+mo
      Qp(-1)  = Averagerho_p(i,j,k-1) + Averagerho(i,j,k-1) + AverageE(i,j,k-1)
      Qp(0 )  = Averagerho_p(i,j,k) + Averagerho(i,j,k) + AverageE(i,j,k)
      Qp(1 )  = Averagerho_p(i,j,k+1) + Averagerho(i,j,k+1) + AverageE(i,j,k+1)
      sm3(i,j,k) = (Qp(-1)-2.0D0*Qp(0)+Qp(1))**2+0.25D0*(Qp(-1)-Qp(1))**2;
    enddo;enddo;enddo
  end select
endsubroutine


subroutine Judge(number,le)
use common1
real tao1,tao2,ct1,ct2,cm1,cm2,beta0,beta1,beta2,beta3
real eps
integer :: i,j,k,number,le
eps = 1.0d-6
signal(:,:,:,le) = .false.
signal2(:,:,:,le) = .true.
select case(le)
case(1)
do k=1,kmax1 ;do  j=2-mo,jmax1+mo-3; do  i=1,imax1;!3,imax1-2;
  beta0 = sm1(i-1,j,k);beta1 = sm1(i,j,k);beta2 = sm1(i+1,j,k);beta3 = sm1(i+2,j,k)
  tao1 = abs(beta0-beta2)
  tao2 = abs(beta1-beta3)

  ct1 = tao1*(beta0*beta1+beta1*beta2+beta2*beta0) 
  ct2 = tao2*(beta1*beta2+beta2*beta3+beta3*beta1)

  cm1 = min(2*beta0**2*min(beta1,beta2),beta0*beta1*beta2+2*beta1**2*min(beta0,beta2),2*beta2**2*min(beta0,beta1))
  cm2 = min(2*beta1**2*min(beta2,beta3),beta1*beta2*beta3+2*beta2**2*min(beta1,beta3),2*beta3**2*min(beta1,beta2))

  if(ct1 <= cm1+eps .and. ct2 <= cm2+eps)then
      signal(i,j,k,le) = .true.;!!!!!!!!!variable-based
      signal2(i,j,k,le) = .true.
   else
      signal(i,j,k,le) = .false.;!!!!!!!!!characteristic-based
      signal2(i,j,k,le) = .false.
   endif
enddo;enddo;enddo;
case(2)
do k=1,kmax1 ;do  j=1,jmax1; do  i=2-mo,imax1+mo-3!do  j=3,jmax1-2; do  i=2-mo,imax1+mo-3
  tao1 = abs(sm2(i,j-1,k)-sm2(i,j+1,k))
  tao2 = abs(sm2(i,j,k)-sm2(i,j+2,k))
  beta0 = sm2(i,j-1,k);beta1 = sm2(i,j,k);beta2 = sm2(i,j+1,k);beta3 = sm2(i,j+2,k)

  ct1 = tao1*(beta0*beta1+beta1*beta2+beta2*beta0) 
  ct2 = tao2*(beta1*beta2+beta2*beta3+beta3*beta1)

  cm1 = min(2*beta0**2*min(beta1,beta2),beta0*beta1*beta2+2*beta1**2*min(beta0,beta2),2*beta2**2*min(beta0,beta1))
  cm2 = min(2*beta1**2*min(beta2,beta3),beta1*beta2*beta3+2*beta2**2*min(beta1,beta3),2*beta3**2*min(beta1,beta2))

  if(ct1 <= cm1+eps .and. ct2 <= cm2+eps)then
      signal(i,j,k,le) = .true.;!!!!!!!!!variable-based
      signal2(i,j,k,le) = .true.
   else
      signal(i,j,k,le) = .false.;!!!!!!!!!characteristic-based
      signal2(i,j,k,le) = .false.
      number = number + 1
   endif
enddo;enddo;enddo;
case(3)
do k=2-mo,kmax1+mo-3 ;
  signal(:,:,k,le) = signal(:,:,1,le)
  signal2(:,:,k,le) = signal2(:,:,1,le)
enddo
end select
endsubroutine


!need to be fixed at 9/22
subroutine Cal_Smoothness_Liu(u,le)  
!Efficient implementation of high-order WENO schemes with sharing function for solving Euler equations
use common1
implicit none
real u(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5)
real Qp(-2:2),d1u(3),d2u(3),a(3),bet(3),is(3),tao5,C,error,sumbet,eps
real,external:: p
integer :: i,j,k,q,le,m
error = 1.0d-40
! eps = 1.0d-6
smLiu1(:,:,:) = 0
smLiu2(:,:,:) = 0
smLiu3(:,:,:) = 0
C=1.0
select case (le)
case(1)
  do k=1,kmax1; ;do  j=1,jmax1; do  i=1,imax1!do  i=3,imax1-2
    Qp(-2) = u(i-2,j,k,1) * p(u(i-2,j,k,:)) * u(i-2,j,k,5)
    Qp(-1) = u(i-1,j,k,1) * p(u(i-1,j,k,:)) * u(i-1,j,k,5)
    Qp(0 ) = u(i,j,k,1) * p(u(i,j,k,:)) * u(i,j,k,5)
    Qp(1 ) = u(i+1,j,k,1) * p(u(i+1,j,k,:)) * u(i+1,j,k,5)
    Qp(2 ) = u(i+2,j,k,1) * p(u(i+2,j,k,:)) * u(i+2,j,k,5)
    d1u(1)=1.d0/2.d0*(Qp(-2)-4.d0*Qp(-1)+3.d0*Qp(0))
    d1u(2)=1.d0/2.d0*(Qp(1)-Qp(-1))
    d1u(3)=1.d0/2.d0*(-3.d0*Qp(0)+4.d0*Qp(1)-Qp(2))

    d2u(1)=1.d0/2.d0*(Qp(-2)-2.d0*Qp(-1)+1.d0*Qp(0))
    d2u(2)=1.d0/2.d0*(Qp(1)-2.d0*Qp(0)+Qp(-1))
    d2u(3)=1.d0/2.d0*(Qp(0)-2.d0*Qp(1)+Qp(2))

    do q=1,3; is(q)=d1u(q)**2.d0+4.d0*d2u(q)**2.d0; enddo
      tao5=abs(is(1)-is(3))
    if(tao5 > min(is(1),is(3))) then
      smLiu1(i,j,k) = 1
    else
      smLiu1(i,j,k) = 0
    endif
  enddo;enddo;enddo
case(2)
  do k=1,kmax1 ;do  j=1,jmax1; do  i=1,imax1
    Qp(-2) = u(i,j-2,k,1) * p(u(i,j-2,k,:)) * u(i,j-2,k,5)
    Qp(-1) = u(i,j-1,k,1) * p(u(i,j-1,k,:)) * u(i,j-1,k,5)
    Qp(0 ) = u(i,j,k,1) * p(u(i,j,k,:)) * u(i,j,k,5)
    Qp(1 ) = u(i,j+1,k,1) * p(u(i,j+1,k,:)) * u(i,j+1,k,5)
    Qp(2 ) = u(i,j+2,k,1) * p(u(i,j+2,k,:)) * u(i,j+2,k,5)
    d1u(1)=1.d0/2.d0*(Qp(-2)-4.d0*Qp(-1)+3.d0*Qp(0))
    d1u(2)=1.d0/2.d0*(Qp(1)-Qp(-1))
    d1u(3)=1.d0/2.d0*(-3.d0*Qp(0)+4.d0*Qp(1)-Qp(2))
  
    d2u(1)=1.d0/2.d0*(Qp(-2)-2.d0*Qp(-1)+1.d0*Qp(0))
    d2u(2)=1.d0/2.d0*(Qp(1)-2.d0*Qp(0)+Qp(-1))
    d2u(3)=1.d0/2.d0*(Qp(0)-2.d0*Qp(1)+Qp(2))
  
    do q=1,3; is(q)=d1u(q)**2.d0+4.d0*d2u(q)**2.d0; enddo
      tao5=abs(is(1)-is(3))
    if(tao5 > min(is(1),is(3))) then
      smLiu2(i,j,k) = 1
    else
      smLiu2(i,j,k) = 0
    endif
  enddo;enddo;enddo
case(3)
  do k=1,kmax1 ;do  j=1,jmax1; do  i=1,imax1
    Qp(-2) = u(i,j,k-2,1) * p(u(i,j,k-2,:)) * u(i,j,k-2,5)
    Qp(-1) = u(i,j,k-1,1) * p(u(i,j,k-1,:)) * u(i,j,k-1,5)
    Qp(0 ) = u(i,j,k,1) * p(u(i,j,k,:)) * u(i,j,k,5)
    Qp(1 ) = u(i,j,k+1,1) * p(u(i,j,k+1,:)) * u(i,j,k+1,5)
    Qp(2 ) = u(i,j,k+2,1) * p(u(i,j,k+2,:)) * u(i,j,k+2,5)
    d1u(1)=1.d0/2.d0*(Qp(-2)-4.d0*Qp(-1)+3.d0*Qp(0))
    d1u(2)=1.d0/2.d0*(Qp(1)-Qp(-1))
    d1u(3)=1.d0/2.d0*(-3.d0*Qp(0)+4.d0*Qp(1)-Qp(2))
  
    d2u(1)=1.d0/2.d0*(Qp(-2)-2.d0*Qp(-1)+1.d0*Qp(0))
    d2u(2)=1.d0/2.d0*(Qp(1)-2.d0*Qp(0)+Qp(-1))
    d2u(3)=1.d0/2.d0*(Qp(0)-2.d0*Qp(1)+Qp(2))
  
    do q=1,3; is(q)=d1u(q)**2.d0+4.d0*d2u(q)**2.d0; enddo
      tao5=abs(is(1)-is(3))
    if(tao5 > min(is(1),is(3))) then
      smLiu3(i,j,k) = 1
    else
      smLiu3(i,j,k) = 0
    endif
  enddo;enddo;enddo
end select
endsubroutine


subroutine Judge_Liu(number,le)  !On shock sensors for hybrid compact/WENO schemes
use common1
real eps
integer :: i,j,k,number,le,m
eps = 0.5
signal(:,:,:,le) = .false.
select case(le)
case(1)
  do k=3-mo,kmax1+mo-3 ;do  j=3-mo,jmax1+mo-3; do  i=3-mo,imax1+mo-4
    if(smLiu1(i,j,k) == 0 .and. smLiu1(i+1,j,k) == 0) then
      signal(i,j,k,le) = .true.
    endif
  enddo;enddo;enddo;
case(2)
  do k=3-mo,kmax1+mo-3 ;do  j=3-mo,jmax1+mo-4; do  i=3-mo,imax1+mo-3
      if(smLiu2(i,j,k) == 0 .and. smLiu2(i,j+1,k) == 0) then
        signal(i,j,k,le) = .true.
      endif
  enddo;enddo;enddo;
case(3)
do k=3-mo,kmax1+mo-4 ;
  signal(:,:,k,le) = signal(:,:,1,le)
enddo
end select
endsubroutine

subroutine Cal_Smoothness_TOT(u,le)
use common1
implicit none
real u(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,1:5)
real Qp(-2:3),fr(5),beta10(5),beta11(5),beta12(5),beta_sum(3),betasum2
real fr1(5),beta110(5),beta111(5),beta112(5),beta_sum1(3),betasum21,h,Midd,Midd2
integer :: i,j,k,q,le
  h = (xma-xmi)/(imaxo1-1)
  select case(le)
    case(1)
    do k=1,kmax1 ;do  j=1,jmax1; do  i=3,imax1-2!3,imax1-2
      do q=1,5
        Qp(-2) = u(i-2,j,k,q)
        Qp(-1) = u(i-1,j,k,q)
        Qp(0 ) = u(i  ,j,k,q)
        Qp(1 ) = u(i+1,j,k,q)
        Qp(2 ) = u(i+2,j,k,q)
        Midd = sqrt((Qp(-2)**2.0+Qp(-1)**2.0+Qp(0)**2.0+Qp(1)**2.0+Qp(2)**2.0)*h)
        
        beta10(q)=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2.0+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2.0;
        beta11(q)=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2.0+0.25D0*(Qp(-1)-Qp(1))**2.0;
        beta12(q)=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2.0+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2.0;
        
        Qp(-2) = u(i+3,j,k,q)
        Qp(-1) = u(i+2,j,k,q)
        Qp(0 ) = u(i+1,j,k,q)
        Qp(1 ) = u(i  ,j,k,q)
        Qp(2 ) = u(i-1,j,k,q)
        Midd2 = sqrt((Qp(-2)**2.0+Qp(-1)**2.0+Qp(0)**2.0+Qp(1)**2+Qp(2)**2.0)*h)
        
        beta110(q)=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2.0+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2.0;
        beta111(q)=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2.0+0.25D0*(Qp(-1)-Qp(1))**2.0;
        beta112(q)=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2.0+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2.0;
        
        if(Midd<1.0d-10)then
          Midd = 1
        endif
        if(Midd2<1.0d-10)then
          Midd2 = 1
        endif
        fr(q) = 1.0/Midd
        fr1(q) = 1.0/Midd2
      enddo

    beta_sum(1) = 1.0/5.0*(fr(1)*beta10(1)+fr(2)*beta10(2)+fr(3)*beta10(3)+fr(4)*beta10(4)+fr(5)*beta10(5))
    beta_sum(2) = 1.0/5.0*(fr(1)*beta11(1)+fr(2)*beta11(2)+fr(3)*beta11(3)+fr(4)*beta11(4)+fr(5)*beta11(5))
    beta_sum(3) = 1.0/5.0*(fr(1)*beta12(1)+fr(2)*beta12(2)+fr(3)*beta12(3)+fr(4)*beta12(4)+fr(5)*beta12(5))
    betasum2 = beta_sum(1)+beta_sum(2)+beta_sum(3)

    beta_sum1(1) = 1.0/5.0*(fr1(1)*beta110(1)+fr1(2)*beta110(2)+fr1(3)*beta110(3)+fr1(4)*beta110(4)+fr1(5)*beta110(5))
    beta_sum1(2) = 1.0/5.0*(fr1(1)*beta111(1)+fr1(2)*beta111(2)+fr1(3)*beta111(3)+fr1(4)*beta111(4)+fr1(5)*beta111(5))
    beta_sum1(3) = 1.0/5.0*(fr1(1)*beta112(1)+fr1(2)*beta112(2)+fr1(3)*beta112(3)+fr1(4)*beta112(4)+fr1(5)*beta112(5))
    betasum21 = beta_sum1(1)+beta_sum1(2)+beta_sum1(3)

    if(betasum2<1.0 .and. betasum21<1.0)then
      signal(i,j,k,le) = .true.
      ! write(*,*) betasum2
    else
      signal(i,j,k,le) = .false.
      ! write(*,*) betasum2
    endif
    enddo;enddo;enddo
    case(2)
      do k=1,kmax1 ;do  j=3,jmax1-2; do  i=1,imax1
        do q=1,5
          Qp(-2) = u(i,j-2,k,q)
          Qp(-1) = u(i,j-1,k,q)
          Qp(0 ) = u(i  ,j,k,q)
          Qp(1 ) = u(i,j+1,k,q)
          Qp(2 ) = u(i,j+2,k,q)
          Midd = sqrt((Qp(-2)**2.0+Qp(-1)**2.0+Qp(0)**2.0+Qp(1)**2+Qp(2)**2.0)*h)
          
          beta10(q)=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2.0+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2.0;
          beta11(q)=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2.0+0.25D0*(Qp(-1)-Qp(1))**2.0;
          beta12(q)=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2.0+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2.0;
          
          Qp(-2) = u(i,j+3,k,q)
          Qp(-1) = u(i,j+2,k,q)
          Qp(0 ) = u(i,j+1,k,q)
          Qp(1 ) = u(i  ,j,k,q)
          Qp(2 ) = u(i,j-1,k,q)
          Midd2 = sqrt((Qp(-2)**2.0+Qp(-1)**2.0+Qp(0)**2.0+Qp(1)**2.0+Qp(2)**2.0)*h)
          
          beta110(q)=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2.0+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2.0;
          beta111(q)=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2.0+0.25D0*(Qp(-1)-Qp(1))**2.0;
          beta112(q)=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2.0+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2.0;
          
          if(Midd<1.0d-10)then
            Midd = 1
          endif
          if(Midd2<1.0d-10)then
            Midd2 = 1
          endif
          fr(q) = 1.0/Midd
          fr1(q) = 1.0/Midd2
        enddo
  
      beta_sum(1) = 1.0/5.0*(fr(1)*beta10(1)+fr(2)*beta10(2)+fr(3)*beta10(3)+fr(4)*beta10(4)+fr(5)*beta10(5))
      beta_sum(2) = 1.0/5.0*(fr(1)*beta11(1)+fr(2)*beta11(2)+fr(3)*beta11(3)+fr(4)*beta11(4)+fr(5)*beta11(5))
      beta_sum(3) = 1.0/5.0*(fr(1)*beta12(1)+fr(2)*beta12(2)+fr(3)*beta12(3)+fr(4)*beta12(4)+fr(5)*beta12(5))
      betasum2 = beta_sum(1)+beta_sum(2)+beta_sum(3)
  
      beta_sum1(1) = 1.0/5.0*(fr1(1)*beta110(1)+fr1(2)*beta110(2)+fr1(3)*beta110(3)+fr1(4)*beta110(4)+fr1(5)*beta110(5))
      beta_sum1(2) = 1.0/5.0*(fr1(1)*beta111(1)+fr1(2)*beta111(2)+fr1(3)*beta111(3)+fr1(4)*beta111(4)+fr1(5)*beta111(5))
      beta_sum1(3) = 1.0/5.0*(fr1(1)*beta112(1)+fr1(2)*beta112(2)+fr1(3)*beta112(3)+fr1(4)*beta112(4)+fr1(5)*beta112(5))
      betasum21 = beta_sum1(1)+beta_sum1(2)+beta_sum1(3)
  
      if(betasum2<1.0 .and. betasum21<1.0)then
        signal(i,j,k,le) = .true.
        ! write(*,*) betasum2
      else
        signal(i,j,k,le) = .false.
        ! write(*,*) betasum2
      endif
      enddo;enddo;enddo
    case(3)
      do k=1,kmax1 ;do  j=1,jmax1 ;do  i=1,imax1
        do q=1,5
          Qp(-2) = u(i,j,k-2,q)
          Qp(-1) = u(i,j,k-1,q)
          Qp(0 ) = u(i  ,j,k,q)
          Qp(1 ) = u(i,j,k+1,q)
          Qp(2 ) = u(i,j,k+2,q)
          Midd = sqrt((Qp(-2)**2.0+Qp(-1)**2.0+Qp(0)**2.0+Qp(1)**2.0+Qp(2)**2.0)*h)
          
          beta10(q)=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2.0+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2.0;
          beta11(q)=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2.0+0.25D0*(Qp(-1)-Qp(1))**2.0;
          beta12(q)=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2.0+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2.0;
          
          Qp(-2) = u(i,j,k+3,q)
          Qp(-1) = u(i,j,k+2,q)
          Qp(0 ) = u(i,j,k+1,q)
          Qp(1 ) = u(i,j  ,k,q)
          Qp(2 ) = u(i,j,k-1,q)
          Midd2 = sqrt((Qp(-2)**2.0+Qp(-1)**2.0+Qp(0)**2.0+Qp(1)**2.0+Qp(2)**2.0)*h)
          
          beta110(q)=(Qp(-2)-2.0D0*Qp(-1)+Qp(0))**2.0+0.25D0*(Qp(-2)-4.0D0*Qp(-1)+3.0D0*Qp(0))**2.0;
          beta111(q)=(Qp(-1)-2.0D0*Qp(0)+Qp(1))**2.0+0.25D0*(Qp(-1)-Qp(1))**2.0;
          beta112(q)=(Qp(0)-2.0D0*Qp(1)+Qp(2))**2.0+0.25D0*(3.0D0*Qp(0)-4.0D0*Qp(1)+Qp(2))**2.0;
          
          if(Midd<1.0d-10)then
            Midd = 1
          endif
          if(Midd2<1.0d-10)then
            Midd2 = 1
          endif
          fr(q) = 1.0/Midd
          fr1(q) = 1.0/Midd2
        enddo
  
      beta_sum(1) = 1.0/5.0*(fr(1)*beta10(1)+fr(2)*beta10(2)+fr(3)*beta10(3)+fr(4)*beta10(4)+fr(5)*beta10(5))
      beta_sum(2) = 1.0/5.0*(fr(1)*beta11(1)+fr(2)*beta11(2)+fr(3)*beta11(3)+fr(4)*beta11(4)+fr(5)*beta11(5))
      beta_sum(3) = 1.0/5.0*(fr(1)*beta12(1)+fr(2)*beta12(2)+fr(3)*beta12(3)+fr(4)*beta12(4)+fr(5)*beta12(5))
      betasum2 = beta_sum(1)+beta_sum(2)+beta_sum(3)
  
      beta_sum1(1) = 1.0/5.0*(fr1(1)*beta110(1)+fr1(2)*beta110(2)+fr1(3)*beta110(3)+fr1(4)*beta110(4)+fr1(5)*beta110(5))
      beta_sum1(2) = 1.0/5.0*(fr1(1)*beta111(1)+fr1(2)*beta111(2)+fr1(3)*beta111(3)+fr1(4)*beta111(4)+fr1(5)*beta111(5))
      beta_sum1(3) = 1.0/5.0*(fr1(1)*beta112(1)+fr1(2)*beta112(2)+fr1(3)*beta112(3)+fr1(4)*beta112(4)+fr1(5)*beta112(5))
      betasum21 = beta_sum1(1)+beta_sum1(2)+beta_sum1(3)
  
      if(betasum2<1.0 .and. betasum21<1.0)then
        signal(i,j,k,le) = .true.
        ! write(*,*) betasum2
      else
        signal(i,j,k,le) = .false.
        ! write(*,*) betasum2
      endif
      enddo;enddo;enddo
  end select
endsubroutine


subroutine Extend(le)
use common1
real tao1,tao2,ct1,ct2,cm1,cm2,beta0,beta1,beta2,beta3
real eps
integer :: i,j,k,number,le
select case(le)
case(1)
do k=1,kmax1 ;do  j=2-mo,jmax1+mo-3; do  i=1,imax1
  if((signal(i,j,k,le) .eqv. .false.) .and. (signal2(i,j,k,le) .eqv. .false.)) then
    signal(i+1,j,k,le)= .false.
    signal(i-1,j,k,le)= .false.
    ! signal(i+2,j,k,le)= .false.
    ! signal(i-2,j,k,le)= .false.
  endif
enddo;enddo;enddo;
case(2)
do k=1,kmax1 ;do  j=1,jmax1; do  i=2-mo,imax1+mo-3
  if((signal(i,j,k,le) .eqv. .false.) .and. (signal2(i,j,k,le) .eqv. .false.)) then
    signal(i,j+1,k,le)= .false.
    signal(i,j-1,k,le)= .false.
    ! signal(i,j+2,k,le)= .false.
    ! signal(i,j-2,k,le)= .false.
  endif
enddo;enddo;enddo;
end select
endsubroutine
