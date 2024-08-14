Program test

  Integer :: i,n
  Real(8) :: dt=0.1
  Real(8) :: f1,f2,rate

  n=20
  f1=0.5
  f2=0.5
  rate=-0.1
  Do i=1,n
     f1=f1+rate*dt
     f2=f2*Exp(-rate*dt)
     Print *,f1,f1
  Enddo
End Program test
