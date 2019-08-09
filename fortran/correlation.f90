module correlation_mod

implicit none

real*8 :: h

private :: h

contains

function correlation(y1,y2,x1,x2,hc_optional)

real*8 :: correlation
real*8, intent(in) :: y1(:),y2(:),x1(:),x2(:)
real*8, optional, intent(in) :: hc_optional
integer :: i,j,k,n,m,max_samples
integer :: start_1,start_2,end_1,end_2,n_xbase
integer :: idx1,idx2
real*8 :: hc
real*8 :: mean1,mean2,lx,ly1,ly2,xmin,xmax
real*8 :: num,den1,den2,w1,w2,den
real*8 :: dx(max(size(x1,1),size(x2,1))-1)
real*8 :: median1,median2,iqr1,iqr2,temp1,temp2
real*8 :: xbase(3*(size(x1,1)+size(x2,1)))
logical :: valid1,valid2
real*8 :: a1,a2,b1,b2
real*8 :: xmid,delta_x
real*8 :: data_for_integration(3*(size(x1,1)+size(x2,1)),6)

hc=0.4d0
if(present(hc_optional)) hc=hc_optional

n=size(y1,1)
m=size(y2,1)

! find valid domain for both series, discard any data more than one point
! outside range of other data series
xmin=max(x1(1),x2(1))
xmax=min(x1(n),x2(m))
if (xmin.eq.x1(1)) then
  start_1=1
  do start_2=2,m
    if (x2(start_2).gt.xmin) exit
  enddo
  start_2=start_2-1
else
  start_2=1
  do start_1=2,n
    if (x1(start_1).gt.xmin) exit
  enddo
  start_1=start_1-1
endif
if (xmax.eq.x1(n)) then
  end_1=n
  do end_2=m-1,1,-1
    if (x2(end_2).lt.xmax) exit
  enddo
  end_2=end_2+1
else
  end_2=m
  do end_1=n-1,1,-1
    if (x1(end_1).lt.xmax) exit
  enddo
  end_1=end_1+1
endif

! find median and interquartile range of data spacing
dx(1:end_1-start_1)=x1(start_1+1:end_1)-x1(start_1:end_1-1)
call merge_sort(dx(1:end_1-start_1))
if (mod(end_1-start_1,2).eq.1) then ! odd number of points
  median1=dx((end_1-start_1)/2+1)
  i=floor((end_1-start_1)/4d0)
  if (mod(end_1-start_1,4).eq.1) then
    temp1=(dx(i)+3*dx(i+1))/4
    temp2=(3*dx(3*i+1)+dx(3*i+2))/4
  else
    temp1=(3*dx(i+1)+dx(i+2))/4
    temp2=(dx(3*i+2)+3*dx(3*i+3))/4
  endif
else
  median1=(dx((end_1-start_1)/2)+dx((end_1-start_1)/2+1))/2
  temp1=(dx((end_1-start_1)/4)+dx((end_1-start_1)/4+1))/2
  temp2=(dx(3*(end_1-start_1)/4+1)+dx(3*(end_1-start_1)/4+2))/2
endif
iqr1=temp2-temp1

dx(1:end_2-start_2)=x2(start_2+1:end_2)-x2(start_2:end_2-1)
call merge_sort(dx(1:end_2-start_2))
if (mod(end_2-start_2,2).eq.1) then ! odd number of points
  median2=dx((end_2-start_2)/2+1)
  i=floor((end_2-start_2)/4d0)
  if (mod(end_2-start_2,4).eq.1) then
    temp1=(dx(i)+3*dx(i+1))/4
    temp2=(3*dx(3*i+1)+dx(3*i+2))/4
  else
    temp1=(3*dx(i+1)+dx(i+2))/4
    temp2=(dx(3*i+2)+3*dx(3*i+3))/4
  endif
else
  median2=(dx((end_2-start_2)/2)+dx((end_2-start_2)/2+1))/2
  temp1=(dx((end_2-start_2)/4)+dx((end_2-start_2)/4+1))/2
  temp2=(dx(3*(end_2-start_2)/4+1)+dx(3*(end_2-start_2)/4+2))/2
endif
iqr2=temp2-temp1

! set interval around points
h=hc*max(median1,median2,iqr1,iqr2)

! construct points where kernel starts, midpoint, stops
j=1
do i=start_1,end_1
  if (x1(i)-h.gt.x1(start_1)) then
    xbase(j)=x1(i)-h
    j=j+1
  endif
  xbase(j)=x1(i)
  j=j+1
  if (x1(i)+h.lt.x1(end_1)) then
    xbase(j)=x1(i)+h
    j=j+1
  endif
enddo
do i=start_2,end_2
  if (x2(i)-h.gt.x2(start_2)) then
    xbase(j)=x2(i)-h
    j=j+1
  endif
  xbase(j)=x2(i)
  j=j+1
  if (x2(i)+h.lt.x2(end_2)) then
    xbase(j)=x2(i)+h
    j=j+1
  endif
enddo
n_xbase=j-1
call merge_sort(xbase(1:n_xbase))
! remove duplicates 
i=2
do
  if (xbase(i).eq.xbase(i-1)) then
    xbase(i:n_xbase-1)=xbase(i+1:n_xbase)
    n_xbase=n_xbase-1
  else
    i=i+1
  endif
  if (i.gt.n_xbase) exit
enddo

! generate points of valid data for integration
idx1=2
idx2=2
j=1
do i=2,n_xbase
  if (xbase(i).le.max(x1(start_1),x2(start_2))) cycle
  if (xbase(i).gt.min(x1(end_1),x2(end_2))) cycle
  ! find both series have a data point within valid interval
  valid1=.false.
  valid2=.false.
  idx1=max(1,idx1-1)
  idx2=max(1,idx2-1)
  do k=idx1,n
    if (abs(x1(k)-(xbase(i-1)+xbase(i))/2d0).le.h) then
      valid1=.true.
      idx1=max(1,k-1)
      exit
    endif
  enddo
  do k=idx2,m
    if (abs(x2(k)-(xbase(i-1)+xbase(i))/2d0).le.h) then
      valid2=.true.
      idx2=max(1,k-1)
      exit
    endif
  enddo
  if (valid1.and.valid2) then
    data_for_integration(j,1:2)=xbase(i-1:i)
    ! linear interpolate
    do k=max(1,idx1-1),n-1
      if (x1(k).ge.xbase(i-1)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,3)=y1(k)+(y1(k+1)-y1(k))*(xbase(i-1)-x1(k))/(x1(k+1)-x1(k))
    do k=max(1,idx1-1),n-1
      if (x1(k).ge.xbase(i)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,4)=(y1(k)+(y1(k+1)-y1(k))*(xbase(i)-x1(k))/(x1(k+1)-x1(k))-data_for_integration(j,3))/(data_for_integration(j,2)-data_for_integration(j,1))
    do k=max(1,idx2-1),m-1
      if (x2(k).ge.xbase(i-1)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,5)=y2(k)+(y2(k+1)-y2(k))*(xbase(i-1)-x2(k))/(x2(k+1)-x2(k))
    do k=max(1,idx2-1),m-1
      if (x2(k).ge.xbase(i)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,6)=(y2(k)+(y2(k+1)-y2(k))*(xbase(i)-x2(k))/(x2(k+1)-x2(k))-data_for_integration(j,5))/(data_for_integration(j,2)-data_for_integration(j,1))
    j=j+1
  endif
enddo

j=j-1

! calculate means by integration
mean1=0d0
mean2=0d0
num=0d0
do i=1,j
  delta_x=data_for_integration(i,2)-data_for_integration(i,1)
  a1=data_for_integration(i,3)
  b1=data_for_integration(i,4)
  mean1=mean1+delta_x*(a1+delta_x*b1/2d0)
  a2=data_for_integration(i,5)
  b2=data_for_integration(i,6)
  mean2=mean2+delta_x*(a2+delta_x*b2/2d0)
  num=num+delta_x
enddo
mean1=mean1/num
mean2=mean2/num

! calculate correlation by integration
data_for_integration(1:j,3)=data_for_integration(1:j,3)-mean1
data_for_integration(1:j,5)=data_for_integration(1:j,5)-mean2
num=0d0
den1=0d0
den2=0d0
do i=1,j
  delta_x=data_for_integration(i,2)-data_for_integration(i,1)
  a1=data_for_integration(i,3)
  b1=data_for_integration(i,4)
  a2=data_for_integration(i,5)
  b2=data_for_integration(i,6)
  num=num+delta_x*(a1*a2+delta_x*(b1*a2/2d0+b2*a1/2d0+delta_x*b1*b2/3d0))
  den1=den1+delta_x*(a1**2+delta_x*(b1*a1+delta_x*b1**2/3d0))
  den2=den2+delta_x*(a2**2+delta_x*(b2*a2+delta_x*b2**2/3d0))
enddo

correlation=num/sqrt(den1*den2)

end function

function merge_arrays(b,c)

real*8, intent(in) :: b(:),c(:)
real*8 :: merge_arrays(size(b,1)+size(c,1))
integer :: i,j,k,n,m

n=size(b,1)
m=size(c,1)

i=1
j=1
do k=1,n+m
  if (i.gt.n) then
    merge_arrays(k:n+m)=c(j:m)
    exit
  elseif (j.gt.m) then
    merge_arrays(k:n+m)=b(i:n)
    exit
  endif
  if (b(i).lt.c(j)) then
    merge_arrays(k)=b(i)
    i=i+1
  else
    merge_arrays(k)=c(j)
    j=j+1
  endif
enddo

end function

subroutine merge_sort(a)

real*8, intent(inout) :: a(:)
integer :: i,m,n

n=size(a,1)

m=1
do
  i=1
  do
    a(i:min(i+2*m-1,n))=merge_arrays(a(i:min(i+m-1,n)),a(i+m:min(i+2*m-1,n)))
    if (i.ge.n-m) exit
    i=i+2*m
  enddo
  if (m.ge.n) exit
  m=m*2
enddo

end subroutine 

end module
