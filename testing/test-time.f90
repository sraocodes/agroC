Program test
  Character(len=10) :: s(3)
  Character(8)  :: date
  Integer :: datetime(8),i
  Character(10) :: time
  Character(5)  :: zone
  Real :: x(5,4,3)
  Real, Allocatable :: a(:)
  Read(*,*) s
  Write(*,'(a,1x,a,1x,a)') s
  Call date_and_Time(date=date,time=time,zone=zone,VALUES=datetime)
  Print *,datetime
  Print *,'date=',date,' time=',time,' zone=',zone
  Print *,Shape(x)
  Print *,Ubound(x,3)
  Print *,Size(x,1)
  Print *,'ubound=',Ubound(x)
  Print *,'ubound=',Ubound(datetime)
  Allocate(a(0))
  Print *,'a=',(a(i),i=1,Ubound(a,1))
End Program test
