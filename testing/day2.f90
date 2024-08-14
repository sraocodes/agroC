Program day
  Implicit None
  Integer, Parameter :: dp=8
  Integer :: stunde,tag
  Real(dp) :: dec,latitude,t
  latitude=50.91_dp ! 50°54'36''
  ! sommer=172
  ! 27. März:     31+28+27=86  Winkel: 41.1 Grad
  ! 17. April: 31+28+31+17=107 Winkel: 49.2 Grad 5:37-19:33 
  tag=107
  Do tag=107,108
     Do stunde=1,24
        t=stunde
        Call winkel(tag, t, latitude)
     End Do
  End Do

Contains

  Subroutine winkel(tag, t, lat)
    Implicit None
    Real(dp), Intent(in) :: t, lat
    Integer, Intent(in) :: tag
    Real(dp), Parameter :: pi=4.0_dp*Atan(1.0_dp)
    Real(dp), Parameter :: rad=pi/180.0_dp
    Real(dp) dec,sinld,cosld,sinb,sunrise,sunset,dl
    Real(dp), Save :: sinbold=0 
    
    dec = -23.4_dp*Cos(2.0_dp*PI*(tag+10.0_dp)/365.0_dp)
    sinld= Sin(dec*rad)*Sin(lat*rad)
    cosld= Cos(dec*rad)*Cos(lat*rad)
    ! daylength (dl=h)
    sunrise = 12-Acos(-Tan(lat*rad)*Tan(dec*rad))/15/rad
    sunset = 12+Acos(-Tan(lat*rad)*Tan(dec*rad))/15/rad
    dl = sunset-sunrise
    if(.not. between(t, sunrise, sunset+1)) then
       sinbOld = 0.0_dp
    end if
!    sinb=Max(0.0_dp,sinld+cosld*Cos(2.0_dp*PI*(t+12.0_dp)/24.0_dp))
    sinb = 0.5_dp*(sinbOld+(sin(dec*rad)*sin(lat*rad)+ &
         Cos(dec*rad)*Cos(15._dp*(t-12._dp)*rad)*Cos(lat*rad)))
!    sinb = sin(dec*rad)*sin(lat*rad)+Cos(dec*rad)*Cos(15.0_dp*(t-12._dp)*rad)*Cos(lat*rad)
    sinb=Max(0.0_dp,sinb)
    sinbold=sinb
    Write(*,'(I3,10F10.4)') Nint(t),sunrise,sunset,sinb,Asin(sinb)*180/PI
  End Subroutine winkel

  Logical Function between(x, a, b)
    ! test if x is between a and b
    Real(dp), Intent(in) :: x, a, b
    between=( a<=x .and. x<=b )
  End Function between

End Program day
