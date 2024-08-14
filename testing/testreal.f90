Program test
  Use ISO_FORTRAN_ENV 
  Real(real32) :: x32
  Real(real64) :: x64, x,y
  Real(real128) :: x128

  Print *,real32,Tiny(x32),Huge(x32),Epsilon(x32),Precision(x32)
  Print *,real64,Tiny(x64),Huge(x64),Epsilon(x64),Precision(x64)
  Print *,real128,Tiny(x128),Huge(x128),Epsilon(x128),Precision(x128)

  x=Tiny(x64)
  y=10.0_real64
  Do i=1,20
     x=x/y
     Print *,'x=',x,x==x
     y=y*10.0_real64
  Enddo
End Program test
