subroutine solutionsource(th)
use common1
real,parameter::gg=1.
do  k=1,k1 ; do  j=1,j1; do  i=1,i1
  v=u1(i,j,k,3)/(u1(i,j,k,1)) ; v=v+gg*th 
  u1(i,j,k,3)=u1(i,j,k,1)*v 
  u1(i,j,k,5)=u1(i,j,k,5)+u1(i,j,k,3)*gg*th 
enddo; enddo; enddo
end subroutine