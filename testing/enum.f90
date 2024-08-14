Program bsp
  implicit none 

  integer :: wert

  enum, bind( C )
    enumerator :: MON, DIE, MIT = 3, DON, FRE = 5, SAM = 66, SON = 77
  end enum  
 
  write( *, * ) MON
  write( *, * ) MIT
  write( *, * ) SAM
  
  wert = 4
  
  If( wert == MIT + 1 ) Then
    write( *, * ) "Donnerstag"
 Endif
End Program bsp
