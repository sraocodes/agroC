Program sunriseset
  Implicit None
  Real(8) :: time, t, sunrise=8.25, sunset=16.75, secperh
  Integer :: i

  Do i=0,23
     t=i
     time = t
     If(t<sunrise) Then
        secperh=0
     else If(t-sunrise < 1) Then
        secPerH = Anint(3600*(t-sunrise))
     Else If(t<sunset) Then
        secPerH = 3600
     Else If(t-sunset<1) Then
        time = sunset
        secPerH = Anint(3600*(1-(t-sunset)))
     Else
        secPerH = 0
     End If
     Print *, secperh
  end do
End Program sunriseset
