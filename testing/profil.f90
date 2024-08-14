Program day
  Implicit None
  Integer, Parameter :: dp=8
  Integer, Parameter :: dim=5
  Integer :: n
  Real(dp) :: z(dim) = (/ 0,1,3,6,10 /)
  Real(dp) :: rrd(dim) = (/ 0.2,0.3,0.4,0.1,0.0 /)
  Real(dp) :: dz(dim)
  Real(dp) :: mass=5,mass2
  Real(dp) :: conc(dim)
  Real(dp) :: conc2(dim)

  dz(1)=0.5*(z(2)-z(1))
  dz(dim)=0.5*(z(dim)-z(dim-1))
  Do n=2,dim-1
     dz(n)=0.5*(z(n+1)-z(n-1))
  Enddo
  conc=mass/(z(dim)-z(1))
  conc2=rrd*mass/dz
  Print *,'Mass=',mass
  Print *,'sum(rrd)=',sum(rrd)
  Print *,' n       z      dz       c     rrd      c2'
  Do n=1,dim
     Write(*,1) n,z(n),dz(n),conc(n),rrd(n),conc2(n)
  Enddo
  mass2=Sum(conc*dz)
  Print *,'Mass=',mass2
  mass2=Sum(conc2*dz)
  Print *,'Mass=',mass2
1 Format(I3,10F8.2)
End Program day
