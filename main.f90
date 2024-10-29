program WCNS_FT_TEST
use input
real timecost,timepartcost,t0,timecost_CPU,timepartcost_CPU
integer problem_index
problem_index = 7
write (*, *) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
write (*, *) "$#                                                     #$"
write (*, *) "$#       High-order 2D WCNS Numerical Test Solver      #$"
write (*, *) "$#                School of Aeronautics                #$"
write (*, *) "$#          Northwestern Polytechnical University      #$"
write (*, *) "$#              (C) Copyright@Xuan Liu                 #$"
write (*, *) "$#             email:bertram0507@163.com               #$"
write (*, *) "$#                                                     #$"
write (*, *) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
call readinput
t0 = timebegin
if(t0 == 0.00000)then
    call gridassemble(problem_index)
    call midgrid
    call initialvalue(problem_index)
    call qpuvw
    call writesolution
else
    call gridassemble(problem_index)
    call midgrid
    call readsolution
endif

!!!!!!!!!Calculation!!!!!!!!!!!
open(11,file='Tcost.plt',status='unknown')
do while(timebegin<timestop)
    call CPU_time(start_CPU) 
    select case(problem_index)
    case(1)
        call solutiontime_RK3_WCNS_2DDoubleMach(problem_index)
    case(2)
        call solutiontime_RK3_WCNS_2DRiemann(problem_index)
    case(3)
        call solutiontime_RK3_WCNS_2DRT(problem_index)
    case(4)
        call solutiontime_RK3_WCNS_2DRichmyer(problem_index)
    case(5)
        call solutiontime_RK3_WCNS_2DShockVortex(problem_index)
    case(6)
        call solutiontime_RK3_WCNS_2DIRV(problem_index)
    case(7)
        call solutiontime_RK3_WCNS_2DRiemann2(problem_index)
    case(8)
        call solutiontime_RK3_WCNS_2DRiemann3(problem_index)
    end select
    call CPU_time(finish_CPU)
    timebegin=timeend				    
    timeend=timebegin+min(timepart,abs(timestop-timebegin))	       
    n=n+1			

        timecost=finish-start
        totaltimecost=timecost+totaltimecost 
        timecost_CPU=finish_CPU-start_CPU 
        totaltimecost_CPU=timecost_CPU+totaltimecost_CPU
        
        write(11,9) n,n*timepart,timecost,totaltimecost,timecost_CPU,totaltimecost_CPU !д�в�ͺ�ʱ
    9 format(I5,5(g16.4,1x))   
    call qpuvw
    if(mod(n,interval)==0) then
        call writesolution
        !call renewinput
    endif			
enddo
call qpuvw
call writesolution
end program
