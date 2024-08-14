Program day
  Implicit None
  Integer, Parameter :: dp=8
  Integer :: stunde
  Real(dp) :: simtime,t
  Real(dp), Parameter :: timefactor=1.0_dp/24.0_dp
  Real(dp), Parameter :: eps=1e-8_dp

  simtime=0.5
  Do stunde=1,48
     t = Nint(((SimTime-1)*TimeFactor-Int((SimTime-1)*TimeFactor+eps))/TimeFactor)   !calculate daytime (from 0 to 23)
     
     Print *,simtime,t
     simtime=simtime+0.5
  End Do

End Program day
