program example1

use correlation_mod

implicit none

integer, parameter :: n=20
real*8 :: y1(n),y2(n),t1(n),t2(n)
integer :: i

open(11,file="test_series1.txt",form="formatted")
do i=1,n
  read(11,*)t1(i),y1(i)
enddo
close(11)

open(12,file="test_series2.txt",form="formatted")
do i=1,n
  read(12,*)t2(i),y2(i)
enddo
close(12)

write(*,'(f8.4)')correlation(y1,y2,t1,t2)

end program
