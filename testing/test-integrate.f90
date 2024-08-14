PROGRAM integriere
  IMPLICIT NONE
  Integer, Parameter :: NumNP=10
  Integer, Parameter :: ntab=5
  Real(8), Dimension(NumNP) :: z=(/0,-1,-2,-3,-4,-5,-6,-7,-8,-9/)
  Real(8), Dimension(NumNP) :: rrd
  Real(8), Dimension(ntab,2) :: tab
  Real(8) :: rna=-0.5
  Real(8) :: rootdepth=-8.5
  Real(8) :: x1,x2,rl,dz,dz2
  Integer :: i,i1,j
  Logical :: last
!  (0,-0.5) (-0.5,-1.5) (-1.5,-2.5) 
!  tab(:,1)=(/ 0.125,0.25,0.375,0.5,1.0 /)
!  tab(:,2)=(/ 0.8817,0.08,0.257,0.126,0.01 /)
  tab(:,1)=(/ 0.0,0.1,0.2,0.3,0.4 /)
  tab(:,2)=(/ 5,4,3,2,1 /)
  Write( *, '(2F7.4)') ((tab(i,j),j=1,2),i=1,ntab)
  Print *,'Integrate: ',IntegrateTab(tab, 0.0_8, 1.0_8)
  call NormalizeRootDistribution()
  Print *,'NormalizeRootDistribution'
  Write( *, '(2F7.4)') ((tab(i,j),j=1,2),i=1,ntab)
  Print *,'Integrate: ',IntegrateTab(tab, 0.0_8, 1.0_8)
  Print *,'CalculateRelativeRootdist'
  Call  CalculateRelativeRootdist()
  Write( *, '(10F7.4)') rrd
  Print *,Sum(rrd)
  
Contains

  Subroutine CalculateRelativeRootdist
    Integer :: i,i1
    rrd=0
    Do i=1,NumNP
       If(rna>z(i)) Exit
    End Do
    If(i>NumNP) Return
    Print *,'i=',i
    i1=Min(Max(i,2),NumNP-1)
    rl=rna-rootdepth
    x1=0
    dz=0.5*(z(i1-1)-z(i1))
    x2=(rna-z(i1-1)+dz)/rl
    Write(*,'(4F7.3)') x1,x2,dz
    ! first node
    If(rna>0.5*(z(i1)+z(i1-1))) Then
       rrd(i1-1)=IntegrateTab(tab,x1,x2) 
       Write(*,'(5F8.4)') x1,x2,x2-x1,dz,rrd(i1-1)
    Endif
    Do i=i1,NumNP-1
       x1=x2
       dz=0.5*(z(i-1)-z(i+1))
       dz2=0.5*(z(i)-z(i+1))
       x2=( rna-z(i)+dz2 )/rl
       rrd(i)=IntegrateTab(tab,x1,x2) 
       Write(*,'(5F8.4)') x1,x2,x2-x1,dz,rrd(i)
       If(rootdepth>=z(i+1)) Exit
    End Do
    ! last node
    If(i<NumNP) Then
       If(rootdepth<0.5*(z(i)+z(i+1))) Then
          x1=x2
          dz=0.5*(z(i)-z(i+1))
          x2=(rna-z(i+1))/rl
          rrd(i+1)=IntegrateTab(tab,x1,x2) 
       Write(*,'(5F8.4)') x1,x2,x2-x1,dz,rrd(i+1)
       Endif
    Endif
!    rrd=rrd/Sum(rrd)
  End Subroutine CalculateRelativeRootdist

  Subroutine NormalizeRootDistribution
    Integer :: j, n1
    Real(8) :: z, xmin
    n1=ntab
    If(n1<=1) Stop 'Invalid distribution of roots'
    xmin=tab(1,1)
    z=tab(n1,1)-xmin
    If(z==0) Stop 'Invalid distribution of roots'
    Do j=1,n1
       tab(j,1)=(tab(j,1)-xmin)/z
    End Do
    ! Integrate the root distribution curve
    ! integral = 0.5 * sum( (y(i)+y(i-1)) * (x(i)-x(i-1)) )
    z=0
    Do j=2,n1
       z=z+(tab(j,2)+tab(j-1,2)) * (tab(j,1)-tab(j-1,1))  
    End Do
    If(z==0) Stop 'Invalid integral of the roots distribution'
    z=0.5*z
    ! set integral to 1
    Do j=1,n1
       tab(j,2)=tab(j,2)/z
    End Do
  End Subroutine NormalizeRootDistribution

  Real(8) Function IntegrateTab(tab, a, b)
    Real(8), Dimension(:,:) :: tab
    Real(8) :: a,b
    Integer :: i
    IntegrateTab=0
    Do i=1,ntab
       If(a<=tab(i,1)) Exit
    End Do
    If(i>ntab .or. b==a) Then
       IntegrateTab=0
    Else If(b<=tab(i,1)) Then
       If(i<=1) Then
          IntegrateTab=0
       Else
          IntegrateTab=0.5*(InterpolateTab(tab,a)+InterpolateTab(tab,b))*(b-a)
       End If
    Else
       If(i>1) Then
          IntegrateTab=0.5*(tab(i,2)+InterpolateTab(tab,a))*(tab(i,1)-a)
       Else
          IntegrateTab=0
       End If
       Do i=i+1,ntab
          If(b<=tab(i,1)) Exit
          IntegrateTab=IntegrateTab+0.5*(tab(i,2)+tab(i-1,2)) *&
               (tab(i,1)-tab(i-1,1))
       End Do
       If(i<=ntab) IntegrateTab=IntegrateTab+0.5*&
            (InterpolateTab(tab,b)+tab(i-1,2)) * (b-tab(i-1,1))
    End If
  End Function IntegrateTab
  
  Real(8) Function InterpolateTab(tab, x)
    Implicit None
    Real(8), Intent(in) :: tab(:,:)
    Real(8), Intent(in) :: x
    Integer :: i
    If(x<=tab(1,1)) Then
       InterpolateTab=tab(1,2)
    Else If(x>=tab(ntab,1)) Then
       InterpolateTab=tab(ntab,2)
    Else
       Do i=2,ntab
          If(x<=tab(i,1)) Exit
       End Do
       InterpolateTab=(x-tab(i-1,1))/(tab(i,1)-tab(i-1,1))* &
            (tab(i,2)-tab(i-1,2))+tab(i-1,2)
    End If
  End Function InterpolateTab

END PROGRAM integriere

