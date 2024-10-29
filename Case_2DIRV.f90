subroutine solutiontime_RK3_WCNS_2DIRV(index)   
    use input; use common1; use scheme_WCNS
    integer nt(3),num_uns,num,num_check,index
    real rhs(i1,j1,k1,5),timenowmin,distance(3),vel(3),timeuse
    real average_th,num_iter
    average_th = 0.0;
    num_iter = 0.0
    vel(:) = 0.0;
    distance(:) = 0.0
    num = 0
    num_check = 0
    timeuse = 0.0
    check_sign = 1
    signal(:,:,:,:) = .false.
    do i=1,3
        nt(i)=0
        timenow(i)=timebegin
    enddo
    timenowmin=timebegin
    le=1
    do while (timenowmin<timeend) 
        th=min(cfl/eigenmax(le),timeend,abs(timeend-timenow(le)))
        average_th = average_th + th
        num_iter = num_iter + 1.0
        if(check_sign == 1)then
            if(num ==0)then
                call Global_Supervise(u1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,1,timeuse)
                call Global_Supervise(u1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,2,timeuse)
                check_sign = 0
            else
                call Global_Supervise(u1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,le,timeuse)
                check_sign = 0
            endif
            num = num +1
        endif
        select case (mix_or_not)
        case(0)  !CF Reconstruction
            call boundaryk(u1,le,index)  
            call spacederivativek_WCNS_CF(rhs,u1,th,le)
            uk1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:) + th*rhs(1:i1,1:j1,1:k1,:)
        !***************!
            call boundaryk(uk1,le,index) 
            call spacederivativek_WCNS_CF(rhs,uk1,th,le)
            uk2(1:i1,1:j1,1:k1,:)=0.75*u1(1:i1,1:j1,1:k1,:)+0.25*uk1(1:i1,1:j1,1:k1,:)+0.25*th*rhs(1:i1,1:j1,1:k1,:)
        !***************!
            call boundaryk(uk2,le,index) 
            call spacederivativek_WCNS_CF(rhs,uk2,th,le,vel(:))
            u1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:)/3.+2*uk2(1:i1,1:j1,1:k1,:)/3.+2*th*rhs(1:i1,1:j1,1:k1,:)/3.
        case(1)  !FT Reconstruction
            call boundaryk(u1,le,index)  
            call spacederivativek_WCNS(rhs,u1,th,le)
            uk1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:) + th*rhs(1:i1,1:j1,1:k1,:)
        !***************!
            call boundaryk(uk1,le,index) 
            call spacederivativek_WCNS(rhs,uk1,th,le)
            uk2(1:i1,1:j1,1:k1,:)=0.75*u1(1:i1,1:j1,1:k1,:)+0.25*uk1(1:i1,1:j1,1:k1,:)+0.25*th*rhs(1:i1,1:j1,1:k1,:)
        !***************!
            call boundaryk(uk2,le,index)
            call spacederivativek_WCNS_update(rhs,uk2,th,le,vel(:))
            u1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:)/3.+2*uk2(1:i1,1:j1,1:k1,:)/3.+2*th*rhs(1:i1,1:j1,1:k1,:)/3.
        case(2)  !Liu Reconstruction & TOT
            call boundaryk(u1,le,index)  
            call spacederivativek_WCNS(rhs,u1,th,le)
            uk1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:) + th*rhs(1:i1,1:j1,1:k1,:)
            call Global_Supervise(uk1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,le,timeuse)
        !***************!
            call boundaryk(uk1,le,index) 
            call spacederivativek_WCNS(rhs,uk1,th,le)
            uk2(1:i1,1:j1,1:k1,:)=0.75*u1(1:i1,1:j1,1:k1,:)+0.25*uk1(1:i1,1:j1,1:k1,:)+0.25*th*rhs(1:i1,1:j1,1:k1,:)
            call Global_Supervise(uk2(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,le,timeuse)
        !***************!
            call boundaryk(uk2,le,index)
            call spacederivativek_WCNS(rhs,uk2,th,le,vel(:))
            u1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:)/3.+2*uk2(1:i1,1:j1,1:k1,:)/3.+2*th*rhs(1:i1,1:j1,1:k1,:)/3.
            call Global_Supervise(u1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,le,timeuse)
        case(3)   !TOT
            call boundaryk(u1,le,index)  
            call spacederivativek_WCNS(rhs,u1,th,le)
            uk1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:) + th*rhs(1:i1,1:j1,1:k1,:)
            call Global_Supervise(uk1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,le,timeuse)
        !***************!
            call boundaryk(uk1,le,index) 
            call spacederivativek_WCNS(rhs,uk1,th,le)
            uk2(1:i1,1:j1,1:k1,:)=0.75*u1(1:i1,1:j1,1:k1,:)+0.25*uk1(1:i1,1:j1,1:k1,:)+0.25*th*rhs(1:i1,1:j1,1:k1,:)
            call Global_Supervise(uk2(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,le,timeuse)
        !***************!
            call boundaryk(uk2,le,index)
            call spacederivativek_WCNS(rhs,uk2,th,le,vel(:))
            u1(1:i1,1:j1,1:k1,:)=u1(1:i1,1:j1,1:k1,:)/3.+2*uk2(1:i1,1:j1,1:k1,:)/3.+2*th*rhs(1:i1,1:j1,1:k1,:)/3.
            call Global_Supervise(u1(1-mo:imax1+mo,1-mo:jmax1+mo,1-mo:kmax1+mo,:),num_uns,le,timeuse)
        end select
    select case(le)
    case(1)
        if(2*(xma-xmi)/Dble(imax1) < abs(distance(1)+abs(vel(1))*th))then
        distance(1) = 0.0
        check_sign = 1
        else
        distance(1) = distance(1) + abs(vel(1))*th
        endif
    case(2)
        if(2*(yma-ymi)/Dble(jmax1) < abs(distance(2)+abs(vel(2))*th))then
        distance(2) = 0.0
        check_sign = 1
        !write(*,*) "y",vel_shock,le
        else
        distance(2) = distance(2) + abs(vel(2))*th
        endif
    case(3)
        if(2*(zma-zmi)/Dble(kmax1) < abs(distance(3)+abs(vel(3)))*th)then
        distance(3) = 0.0
        check_sign = 1
        else
        distance(3) = distance(3) + abs(vel(3))*th
        endif
    end select
        timenow(le)=timenow(le)+th
        nt(le)=nt(le)+1
        write(*,*) le,nt(le),timenow(le),th,vel(1),vel(2),distance(1)/((xma-xmi)/Dble(imax1)), &
        distance(2)/((yma-ymi)/Dble(jmax1))
        timenowmin=timenowminimum(timenow,le)
        vel(:) = 0.0
        enddo 
        write(*,*) "Superviser timeuse :",timeuse
        write(*,*) "Average    timestep:",average_th / num_iter
    end subroutine