Program day
  Implicit None
  Integer, Parameter :: dp=8
  Integer :: stunde
  Real(dp) :: dec,tag,latitude
  latitude=50.91_dp ! 50°54'36''
  ! sommer=172
  ! 27. März:     31+28+27=86  Winkel: 41.1 Grad
  ! 17. April: 31+28+31+17=107 Winkel: 49.2 Grad 5:37-19:33 
  tag=107
  Do stunde=0,23
     Call winkel(tag, stunde, latitude)
  End Do

Contains

  Subroutine winkel(tag, stunde, lat)
    Implicit None
    Real(dp), Intent(in) :: tag, lat
    Integer, Intent(in) :: stunde
    Real(dp), Parameter :: pi=4.0_dp*Atan(1.0_dp)
    Real(dp), Parameter :: rad=pi/180.0_dp
    Real(dp) dec,sinld,cosld,sinb,sunrise,sunset
    Real(dp), Save :: sinbold=0 
    
    dec = -23.4_dp*Cos(2.0_dp*PI*(tag+10.0_dp)/365.0_dp)
    sinld= Sin(dec*rad)*Sin(lat*rad)
    cosld= Cos(dec*rad)*Cos(lat*rad)
    sinb=Max(0.0_dp,sinld+cosld*Cos(2.0_dp*PI*(stunde+12.0_dp)/24.0_dp))
    Write(*,'(I3,10F10.4)') stunde,dec*180/PI,sinb,Asin(sinb)*180/PI
  End Subroutine winkel

End Program day
