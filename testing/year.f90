Program year
  Implicit None
  Integer, Parameter :: dp=8
  Integer :: stunde,tag
  Real(dp) :: decl,latitude,sind,time,dt,w,elevation,dsinbe
  Real(dp), Parameter :: pi=4*atan(1.0_dp)
  Real(dp), Parameter :: rad=pi/180.0_dp
  latitude=50.91_dp ! 50°54'36''
  
  ! sommer=172
  ! 27. März:     31+28+27=86  Winkel: 41.1 Grad
  ! 17. April: 31+28+31+17=107 Winkel: 49.2 Grad 5:37-19:33
  ! 18. April: 31+28+31+18=108 Winkel: 49.6 Grad 5:35-19:35

  ! yearly declination
  time=108.0
  dt=1.0_dp/24.0_dp
  Do tag=1,1
     Do stunde=1,24
        time=time+dt
        decl=Asin(0.39795_dp*Cos(0.98563_dp*(time-173)*rad))
        w=15*(stunde-12)*rad
        elevation=Asin(Sin(decl)*Sin(latitude*rad)+Cos(decl)*Cos(w)*Cos(latitude*rad))
        dsinbe = Sin(Asin(elevation)+0.4_dp*elevation)
        Print *,time,decl/rad,w/rad,elevation/rad,dsinbe
     End Do
  End Do

  

End Program year
