module input
common /inputdata/qin,pin,gin,sm,cfl,cflvis, &
id_scheme,timebegin,timepart,timestop,interval
real uin
real timeend
integer ncfl
real timenow(3)
endmodule

module ghost
integer,parameter:: mo=8
integer,parameter:: mo1=2
endmodule 

module common1
use ghost
integer::i1,j1,k1,check_sign
integer::imax1,jmax1,kmax1,imaxo1,jmaxo1,kmaxo1
real:: xmi,xma,ymi,yma,zmi,zma
real,allocatable,save:: xo1(:,:,:),yo1(:,:,:),zo1(:,:,:)                        
real,allocatable,save:: x1(:,:,:),y1(:,:,:),z1(:,:,:)    
real,allocatable,save:: ax1(:,:,:),ay1(:,:,:),az1(:,:,:), &
                        bx1(:,:,:),by1(:,:,:),bz1(:,:,:), &	   
                        cx1(:,:,:),cy1(:,:,:),cz1(:,:,:)    	   
real,allocatable,save:: jac1(:,:,:)                                
real,allocatable,save:: qs1(:,:,:),ps1(:,:,:),us1(:,:,:),vs1(:,:,:),ws1(:,:,:)          
real,allocatable,save:: qs1_sum(:,:,:),ps1_sum(:,:,:),us1_sum(:,:,:),vs1_sum(:,:,:),ws1_sum(:,:,:)
real,allocatable,save:: u1(:,:,:,:),uk1(:,:,:,:),uk2(:,:,:,:),uk3(:,:,:,:),uk(:,:,:,:,:),ukresidual(:,:,:,:)
real,allocatable,save:: xb1(:,:),yb1(:,:),zb1(:,:)
real,allocatable,save:: sm1(:,:,:),sm2(:,:,:),sm3(:,:,:)
integer,allocatable,save:: smLiu1(:,:,:),smLiu2(:,:,:),smLiu3(:,:,:)
real,allocatable,save:: Averagerho_p(:,:,:)
real,allocatable,save:: Averagerho(:,:,:)
real,allocatable,save:: AverageE(:,:,:)
logical,allocatable,save:: signal(:,:,:,:)
logical,allocatable,save:: signal2(:,:,:,:)
real,allocatable,save:: vel_s1(:,:,:)
real,allocatable,save:: vel_s2(:,:,:)
real,allocatable,save:: vel_s3(:,:,:)
end module

module scheme_WCNS
integer,parameter::npre=1001
real::apre(1001,3)
integer,parameter::id_scheme_WCNS=2  !1 WENO-JS; 2 WENO-Z; 3 TENO 4 CU6
integer,parameter::mix_or_not=3      !0:NOMix 1:FT 2:Liu Shengping 2020 CAF 3:TOT
integer,parameter::id_scheme_WCNS_E=5!4 E4; 5 E5
integer,parameter::eigntype=1  !1 detector; 2 u
integer,parameter::fretype=1 !1 wre; 2 3th smothness indictor
integer,parameter::iteration=1
endmodule

