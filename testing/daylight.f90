Program daylight
  Implicit None
  Integer, Parameter :: dp=8
  Integer :: i, sommer
  Real(dp) dec,rad,pi,sinld,cosld,sunrise,sunset,dl,dayOfYear,latitude
  PI=4.0_dp*Atan(1.0_dp)
  latitude=50.866725_dp
  sommer=31+28+31+30+31+21 ! 172.5
  Do i=1,365
     dayOfYear=i
     dl = daylength(dayOfYear, latitude)
     Print *,i,dl
  End Do

Contains

  Real(dp) Function daylength(day, lat)
    Implicit None
    Real(dp), Intent(in) :: day, lat
    Real(dp), Parameter :: pi=4.0_dp*Atan(1.0_dp)
    Real(dp), Parameter :: rad=pi/180.0_dp
    Real(dp) dec,sinld,cosld,sinb,sunrise,sunset
    
    dec = -23.4_dp*Cos(2.0_dp*PI*(day+10.0_dp)/365.0_dp)
    sinld= Sin(dec*rad)*Sin(lat*rad)
    cosld= Cos(dec*rad)*Cos(lat*rad)
    sinb=Max(null_dp,sinld+cosld*Cos(2.0_dp*PI*(hour+12.0_dp)/24.0_dp))
    sunrise = 12-Acos(-Tan(lat*rad)*Tan(dec*rad))/15/rad
    sunset = 12+Acos(-Tan(lat*rad)*Tan(dec*rad))/15/rad
    daylength = sunset-sunrise
  End Function daylength

End Program daylight
