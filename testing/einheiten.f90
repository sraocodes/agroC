Program einheiten
  Implicit None
  Integer, Parameter :: dp=8
  Real(dp), Parameter :: null_dp=0.0_dp
  Real(dp), Parameter :: one_dp=1.0_dp
  Real(dp) :: dec,tag,latitude,dpar,global_rad
  Integer :: stunde=12
  Real(dp), Parameter :: pi=4.0_dp*Atan(1.0_dp)
  Real(dp), Parameter :: rad=pi/180.0_dp
  Real(dp) :: sinld,cosld,sinb,dl,sunrise,sunset,dsinb,dsinbe
  Real(dp) :: frdf,apar,pardf,pardr
  Real(dp) :: atmtr,detr
  Real(dp) :: scp,parlpp,amax,amdvs,amtmp,asssl,eff
  Real(dp) :: UnitFactor, AreaFactor

  UnitFactor=10.0
  AreaFactor = 0.0001/(1000.0/UnitFactor)**2
  ! sommer=172
  ! 27. März: 31+28+27=86 Winkel: 41.1 Grad
  tag=86
  latitude=50.91_dp ! 50°54'36''
  dec = -23.4_dp*Cos(2.0_dp*PI*(tag+10.0_dp)/365.0_dp)
  sinld= Sin(dec*rad)*Sin(latitude*rad)
  cosld= Cos(dec*rad)*Cos(latitude*rad)
  sinb=Max(0.0_dp,sinld+cosld*Cos(2.0_dp*PI*(stunde+12.0_dp)/24.0_dp))
  Write(*,'(I3,10F10.4)') stunde,dec*180/PI,sinb,Asin(sinb)*180/PI
  sunrise = 12-Acos(-Tan(latitude*rad)*Tan(dec*rad))/15/rad
  sunset = 12+Acos(-Tan(latitude*rad)*Tan(dec*rad))/15/rad
  dl = sunset-sunrise
  Print *, 'daylength=',dl
  dsinb=3600.0_dp*(dl*sinld+24.0_dp*cosld*Sqrt(one_dp-(sinld/cosld)**2)/PI)
  dsinbe=3600.0_dp*(dl*(sinld+0.4_dp*(sinld*sinld+0.5_dp*cosld*cosld))+ &
       12.0_dp*cosld*(2.0_dp+3.0_dp*0.4_dp*sinld)*Sqrt(one_dp-(sinld/cosld)**2)/PI)
  detr=1370.0_dp*(one_dp+0.033_dp*Cos(2.0_dp*PI*tag/365.0_dp)) ! (J/m2/s)
  global_rad=1000.1
  dpar=0.5_dp*global_rad*(1000.0_dp/UnitFactor)**2  ! J/L2/T -> J/m2/T
  If(dsinb>null_dp) Then
     atmtr=global_rad*(1000.0_dp/UnitFactor)**2/(detr*dsinb) ! dimensionless
     If(atmtr.Le.0.07_dp) Then
        frdf=one_dp
     Else If(atmtr.Le.0.35_dp) Then
        frdf=one_dp-2.3_dp*(atmtr-0.07_dp)**2
     Else If(atmtr.Le.0.75_dp) Then
        frdf=1.333_dp -1.46_dp*atmtr
     Else  
        frdf=0.23_dp
     Endif
  Else
     frdf=null_dp
  Endif

  apar=dpar*sinb*(1.0_dp+0.4_dp*sinb)/dsinbe
  pardf=Min(apar,frdf*apar*sinb/dsinb)
  pardr=apar-pardf
  Print *,'apar=',apar,' pardf=',pardf,' pardr=',pardr

  scp=0.2 ! SucrosPlant(ip)%scp
  parlpp=pardr*(one_dp-scp)/sinb
  eff=0.45*AreaFactor ! SucrosPlant(ip)%eff*AreaFactor
  amdvs=0.9
  amtmp=0.9
  amax=40.0*amdvs*amtmp*AreaFactor ! SucrosPlant(ip)%amx*amdvs*amtmp*AreaFactor
  
  asssl=amax*(one_dp-Exp(-eff*parlpp/amax))
  Print *,'parlpp=',parlpp,' amax=',amax,' eff=',eff,' asssl=',asssl
  

End Program einheiten
