Module Nitrogen
  Use datatypes
  Use TimeData, only: TimeFactor
  Use Geometry, Only: NumNP, Urea, NH4, NO3, NSnit, NitSink, NitSinkRoot
  Use Variables, Only: lNitrogen, LNPools
  Implicit None
  Logical :: lNRootUpt, lNreduc, lNfixat
  Real(dp) :: Rkfixopt    ! [M L-2 t-1] optimum F fixation rate
  Real(dp) :: denThresh   ! [-]  threshold value of effective saturation (WFP) for denitrification 
  Real(dp) :: Ro    ! [-] C/N ratio of the biomass
  !Real(dp) :: Fe_in ! [-] Synthesis efficiency constant; default WAVE = 0.3,
  !                  ! negative Value indicates the Use of the RothC Value computed internally
  Real(dp) :: Fh    ! [-] Humification coefficient; default WAVE = 0.4,
                    ! negative Value indicates the Use of the RothC Value= 0.54
  Real(dp) :: Rknitri    ! [1/T] Nitrification constant [h]
  Real(dp) :: Rkdenit    ! [1/T] Denitrification constant [h] 
  Real(dp) :: Rkhyd      ! [1/T] Ureum hydrolyse constant [h]
  Real(dp) :: Rkvol      ! [1/T] Volatilization constant [h]
  Real(dp) :: rdo        ! [1/L] D0: travel distance resistence between bulk soil and root
    
  Real(dp) :: C2Nfact    ! factor to create initial N pool from initial C pool
  Real(dp) :: UnitFactorC ! factor to create initial N pool from initial C pool
  Real(dp) :: UnitFactorWave ! factor to convert from input to wave
  Real(dp), Allocatable :: PoolNitMan(:),PoolNitHum(:),PoolNitLit(:)
  Real(dp), Allocatable :: rmin(:,:)
  Real(dp), Allocatable :: reddenit(:), redwfps(:)
  Real(dp), Allocatable :: rhurea(:), rnitri(:), rdenit(:), rvol(:), rConc(:), NFact(:)
  Real(dp), Allocatable :: rn2o(:), rn2(:)
  Real(dp), Allocatable :: unc(:,:), und(:,:)
  ! boundary for solute transport
  Real(dp), Allocatable :: N_on_top(:)
  Integer :: NTop_dim, NTop_ind=1
  Real(dp), Allocatable :: NTopTime(:), NTopMass(:,:)
  ! factors for converting concentrations
  Real(dp) :: urea_to_nh4, nh4_to_no3, nh4_to_n2o, nh4_to_nh3, no3_to_n2o, no3_to_n2 
  ! in Wave
  ! parasol: 1 bulk density,                            - chpar(1) (M/L3)
  !          2 distribution coefficient (l/kg)          - chpar(7)
  !          3 chem diffusion in pure water (mm^2/day)  - chpar(5)
  !          4 a coefficient                            - 
  !          5 b coefficient                            - 
  !          6 hydrodynamic dispersivity (mm)           - chpar(2)
  !          7 mobile total coefficient (-)             - 1
  !          8 mass transfer coefficient (1/day)        - 0
  !          9 adsorbed fraction in the mobile zone (-) - 1
  Type NitSinkBalanceType
     Real(dp) :: rndemlv
     Real(dp) :: rndemst
     Real(dp) :: rndemrt
     Real(dp) :: rndemso
     Real(dp) :: rndemcrn
     Real(dp) :: tunc
     Real(dp) :: totdem
     Real(dp) :: tund
     Real(dp) :: an_tot
     Real(dp) :: anlv
     Real(dp) :: anst
     Real(dp) :: anrt
     Real(dp) :: anso
     Real(dp) :: ancrn
     Real(dp) :: NitReduct
     Real(dp) :: anclv
     Real(dp) :: ancst
     Real(dp) :: ancrt
     Real(dp) :: ancso
     Real(dp) :: anccrn
     Real(dp) :: sinksum(NSnit)
     Real(dp) :: sinkrootsum(NSnit)
     Real(dp) :: rhurea, rnitri, rn2o, rn2, rdenit, rvol
     Real(dp) :: Rkfixact
     Real(dp) :: tunc_tot, tund_tot
  End Type NitSinkBalanceType
  Type(NitSinkBalanceType) :: NitSinkBalance
Contains

  Subroutine NitrogenIn(output)
    Use Carbon, Only: co2alf
    Use Geometry, Only: MatNum, coord
    Use Material, Only: NMat
    Use Solute, Only: NS
    Use Variables, Only: rorad, UnitFactorPlants, w0_dens


    Implicit None
    Logical, Intent(in) :: output(:)
    Integer :: i, m, nfrom,nto,nby,ndim
    Real(dp), Allocatable :: n_init(:,:)
    Real(dp) :: exp_alt, exp_neu, w

    If(debug) Print *,'NitrogenIn'
    If(.Not.lNitrogen) Then
       Print *,'Skipping block I'
       If(headline(1:12) == '*** BLOCK I:') Call SkipBlock('I')
       Return
    Endif
    Write(*,*) 'reading nitrogen input'
    read(30,*)
    Read(30,*) lNpools, lNRootUpt, lNfixat, Rkfixopt, denThresh
!    Read(30,*) lNpools, lNRootUpt, lNreduc
    lNreduc=.True. ! not used at the moment
    read(30,*)
    Read(30,*) Ro
    !Read(30,*) Fe_in
    Read(30,*) Fh
    if(Fh<null_dp) Fh=0.54_dp
    read(30,*)
    read(30,*) Rknitri
    read(30,*) Rkdenit
    read(30,*) Rkhyd
    read(30,*) Rkvol
    read(30,*) w0_dens
    read(30,*) rorad
    Read(30,*) rdo
    Read(30,*)
    Read(30,*) C2Nfact
    ! nitrogen mass boundary
    Read(30,*)
    ! read fertilizer input
    Read(30,*) NTop_dim
    Allocate(NTopTime(NTop_dim))
    Allocate(NTopMass(NSNit,NTop_dim))
    Read(30,*)
    Do i=1,NTop_dim
       Read(30,*) NTopTime(i),NTopMass(:,i)
    Enddo
    Read(30,*)
    ! read initial pools
    if(co2alf.gt.null_dp) then
       nfrom=1
       nto=1
       nby=1
       ndim=1
    Else If(equal(co2alf, -1.0_dp)) Then
       nfrom=NumNP
       nto=1
       nby=-1
       ndim=NumNP
    else 
       nfrom=1
       nto=NMat
       nby=1
       ndim=NMat
    endif
    Allocate(n_init(NSnit,ndim))
    Do i=nfrom,nto,nby
       Read(30,*) n_init(:,i)
    Enddo
    Read(30,1) headline
    If(output(1)) Then
       Write(50,'(/A)') 'Nitrogen input'
       Write(50,'(A/)') '=============='
       Write(50,*) 'lNpools: ',lNpools
       Write(50,1) 'Ro: ',Ro
       !Write(50,1) 'Fe: ',Fe_in
       Write(50,1) 'Fh: ',Fh
       Write(50,1) 'Rknitri: ',Rknitri
       Write(50,1) 'Rkdenit: ',Rkdenit
       Write(50,1) 'Rkhyd: ',Rkhyd
       Write(50,1) 'Rkvol: ',Rkvol
       Write(50,1) 'Factor C to N: ',C2Nfact
       Write(50,'(A,/,1P,(3e13.5))') 'n_init:',n_init
    End If
1   Format(A,1P,E13.5)
    ! Allocate fields
    Allocate(PoolNitMan(NumNP))
    Allocate(PoolNitHum(NumNP))
    Allocate(PoolNitLit(NumNP))
    Allocate(rmin(NumNP, NSnit))
    Allocate(unc(NumNP, NSnit))
    Allocate(und(NumNP, NSnit))
    Allocate(rhurea(NumNP))
    Allocate(rnitri(NumNP))
    Allocate(rn2o(NumNP))
    Allocate(rn2(NumNP))
    Allocate(rdenit(NumNP))
    Allocate(rvol(NumNP))
    Allocate(rConc(NSnit))
    Allocate(NFact(NS))
    Allocate(reddenit(NumNP))
    Allocate(redwfps(NumNP))
    Allocate(N_on_top(NSnit))
    N_on_top=0
    rhurea=0
    rnitri=0
    rn2o=0
    rn2=0
    NitSink=0
    NitSinkRoot=0
    ! init pools
    If(co2alf.Le.null_dp) Then
       ! equally distribute init pool material on the nodes
       Do i=1,NumNP
          If(equal(co2alf, -1.0_dp)) Then
             m=i
          Else 
             m=MatNum(i)
          Endif
          PoolNitLit(i)=n_init(1,m)
          PoolNitMan(i)=n_init(2,m)
          PoolNitHum(i)=n_init(3,m)
       End Do
    Else
       ! exponentially distribute init pool material on the nodes
       exp_alt=1.0
       exp_neu=1.0
       Do i=NumNP,2,-1
          exp_neu=Exp((coord(i-1)-coord(NumNP))*co2alf)
          w=exp_alt-exp_neu
          !write(*,*) w,coord(NumNP)-coord(i-1)
          exp_alt=exp_neu
          PoolNitLit(i)=n_init(1,1)*w
          PoolNitMan(i)=n_init(2,1)*w
          PoolNitHum(i)=n_init(3,1)*w
       End Do
       w=exp_neu
       !write(*,*) w,coord(NumNP),coord(1)
       PoolNitLit(1)=n_init(1,1)*w
       PoolNitMan(1)=n_init(2,1)*w
       PoolNitHum(1)=n_init(3,1)*w
    Endif
    Deallocate(n_init)
    !PoolNitLit=C2Nfact*PoolDPM
    !PoolNitMan=C2Nfact*PoolRPM
    !PoolNitHum=C2Nfact*PoolHum
    ! convert to Wave (mg/l)
    UnitFactorWave=Units_M(Unit_M_Input)/(Units_L(Unit_L_Input)/Units_L(3))**3
    ! convert mass unit from plants (kg) to agroc (mg,g,kg,t)
    UnitFactorPlants=Units_M(3)/Units_M(Unit_M_Input)
    rConc(Urea)=0
    rConc(NH4)=0 !1.8/UnitFactorWave ! M/(L3 water)/T
    rConc(NO3)=0 !6.2/UnitFactorWave
    NFact=1      ! for other solutes than N
    NFact(Urea)=14.0/60.0 ! factor N/Urea
    NFact(NH4)=14.0/18.0  ! factor N/NH4
    NFact(NO3)=14.0/62.0  ! factor N/NO3
    urea_to_nh4=14.0/60.0
    nh4_to_no3=62.0/18.0
    nh4_to_n2o=44.0/18.0
    nh4_to_nh3=17.0/18.0
    no3_to_n2o=44.0/62.0
    no3_to_n2=28.0/62.0
    Call reset_NitSinkBalance()
  End Subroutine NitrogenIn


  Subroutine nit_sink(dt)
    !     in   : csolio, csolmo, decsoli, decsolm, dt, imobsw, ncs, nla,
    !            nr_of_sol, parasol, wcio, wcmo
    !     out  : decsoli, decsolm, sinki, NitSink
    !     calls: nit_upt, seqtrans, sol_sink
    !###################################################################################
    Implicit None
    Real(dp), Intent(in) :: dt

    If(debug) Print *,'nit_sink'
    If(lNitrogen) Then
       Call nit_upt(dt)
       Call seqtrans(dt)
       Call nit_correct_sink(dt)
    End If
  End Subroutine nit_sink
  
  Subroutine nit_upt(dt)
    Use Geometry, Only: NumNP, ThNew, MatNum, cRootMax, delta_z, Sink
    Use Solute, Only: ChPar, Conc
    Use Plants, Only: SucrosPlant, InterpolateTab, Plant, PlantData,&
         PlantGeometry, PlantGeometryArray, Potatoe, SugarBeet, NitReduct, col_fact
    Use Material, Only: ThS
    Use Timedata, Only: Simtime
    Use Variables, Only: rorad, UnitFactorPlants, w0_dens,anlv, anst, anrt, anso, ancrn,&
         Nit_ancrt, an_tot, tunc_tot, tund_tot
    Use Carbon, Only: sfact
    
    Implicit None
    Real(dp), Intent(in) :: dt
    Real(dp) :: xncle, xncst, xncrt, xncso, xnccrn, rmncl, fndef
    Real(dp) :: rndemlv, rndemst, rndemrt, rndemso, rndemcrn, totdem
    Real(dp) :: fixfraclv, fixfracst, fixfracrt, fixfraccrn, Rkfixact, Rkfixdt
    Real(dp) :: fixactlv, fixactst, fixactrt, fixactcrn
    Real(dp) :: anclv=0, ancst=0, ancrt=0, anccrn=0, ancso=0
    Real(dp) :: dvs, temp_sum, tunc, hcsolo, tund, tpdnup, diffus_rm, Dw, TauW, rdens
    Real(dp) :: totvr, rnuplv, rnupst, rnuprt, rnupcrn
    Real(dp) :: rootdepfac, meansfac, nstressfac
    Integer :: ip,ptype, sp, m, i
    Integer, Save :: oldseason=0
    Real(dp), Parameter :: rlncl=0.005
    Real(dp) :: sinksum
    Character(128) :: line
    Logical :: neq

    If(debug) Print *,'nit_upt'
  ! uptake from plants
    If(.Not.lNRootUpt) Return
    ip=PlantGeometry%index
    ptype=Plant(ip)%type ! plant type
    If(Plant(ip)%season/=oldseason) Then
       If(oldseason==0 .Or. ip==0 .Or. Plant(ip)%season==0) Then
          neq=.False.
       Else
          neq=(SucrosPlant(ip)%emergence(oldseason) /= SucrosPlant(ip)%emergence(Plant(ip)%season))
       Endif
       If(oldseason==0 .Or. neq) Then
          ! s_ancl = 0.05 in wave
          If(ip>0) Then
             anlv = 0.05*SucrosPlant(ip)%ssl*SucrosPlant(ip)%nsl / &
                  (SucrosPlant(ip)%sla/PlantData%Area_ha_to_L2)
          Else
             anlv=0
          Endif
          anst=0
          anrt=0
          anso=0
          ancrn=0
          an_tot=0
          Call reset_NitSinkBalance()
          NitReduct=1
          NitSink=0
       Endif
       oldseason=Plant(ip)%season
    Endif
    NitSinkRoot=0
    If(ip==0 .Or. Plant(ip)%season==0) Return
    dvs=PlantGeometry%dvs
    temp_sum=PlantGeometry%temp_sum
    If(ptype==Potatoe .Or. ptype==SugarBeet) Then
       xncle=InterpolateTab(SucrosPlant(ip)%tab(13),temp_sum)
       xncst=InterpolateTab(SucrosPlant(ip)%tab(14),temp_sum)
       xncrt=InterpolateTab(SucrosPlant(ip)%tab(15),temp_sum)
       xnccrn=InterpolateTab(SucrosPlant(ip)%tab(17),temp_sum)
    Else
       xncle=InterpolateTab(SucrosPlant(ip)%tab(13),dvs)
       xncst=InterpolateTab(SucrosPlant(ip)%tab(14),dvs)
       xncrt=InterpolateTab(SucrosPlant(ip)%tab(15),dvs)
       xnccrn=InterpolateTab(SucrosPlant(ip)%tab(17),dvs)
    End If
    rmncl=0.5d0*xncle
    if (rmncl.eq.0.d0) then
       call WriteError('rmncl is zero')
    endif
    If(PlantGeometry%wlvg>0) anclv=anlv/PlantGeometry%wlvg
    If(PlantGeometry%wst>0) ancst=anst/PlantGeometry%wst 
    If(PlantGeometry%wrt>0) ancrt=anrt/PlantGeometry%wrt 
    If(PlantGeometry%wso>0) ancso=anso/PlantGeometry%wso 
    If(PlantGeometry%wcrn>0) anccrn=ancrn/PlantGeometry%wcrn 
    NitReduct = nit_reduct(anclv, rlncl, rmncl)
    !Print *,'NitReduct=',simtime%time,NitReduct,anclv,rlncl,rmncl
    fndef=1-Sqrt(1-NitReduct*NitReduct)
    !fndef=1-Sqrt(1-NitReduct)
    rndemlv=Max(0.0,(PlantGeometry%wlvg+PlantGeometry%wlvd)*xncle-anlv) ! kg N / L2
    rndemst=Max(0.0, PlantGeometry%wst*xncst-anst)
    rndemrt=Max(0.0, PlantGeometry%wrt*xncrt-anrt)
    rndemcrn=Max(0.0, PlantGeometry%wcrn*xnccrn-ancrn)
    If(ptype==Potatoe .Or. ptype==SugarBeet) Then
       xncso=InterpolateTab(SucrosPlant(ip)%tab(16),temp_sum)
    Else
       xncso=InterpolateTab(SucrosPlant(ip)%tab(16),dvs)
    Endif
    rndemso= Min( Max(0.d0, PlantGeometry%wso*xncso-anso)*fndef, rndemlv+rndemst)
    !Write(*,'(a,1p,9e13.5)') 'rndemso=', rndemso,PlantGeometry%wso,anso,xncso,(PlantGeometry%wso*xncso-anso)*fndef, rndemlv+rndemst
    totdem=Max(0.d0,rndemlv+rndemst+rndemrt+rndemcrn) ! kg N / L2 / T bei uns
    !   N fixation    M.Herbst 09/09/2022
    if (lNfixat .And.totdem>0.d0) then
       !   computation of the root density and element thickness weighted average reduction factor
       !   for soil temperature, water content and aeration (too much CO2)
       meansfac=0.d0       
       do i=1,NumNP     
          meansfac=meansfac+(sfact(1,i)*sfact(2,i)*sfact(3,i)*((PlantGeometryArray%rrd(NumNP-i+1)&
               +abs(delta_z(i)/(PlantGeometry%rootdepth/PlantData%UnitFactor)))/2.d0))
       enddo
       ! reducion according to actual root depth/ maximum root depth
       rootdepfac=PlantGeometry%rootdepth/SucrosPlant(ip)%root_max
       ! increase of N fixation according to N deficit in the crop
       nstressfac=1.d0+1.d0-NitReduct   
       ! actual N fixation rate
       Rkfixdt=Rkfixopt*rootdepfac*meansfac*col_fact(ptype)*nstressfac*dt 
       fixfraclv=rndemlv/totdem     !fraction of leaves in total N demand
       fixfracst=rndemst/totdem
       fixfracrt=rndemrt/totdem
       fixfraccrn=rndemcrn/totdem
       fixactlv=Min(rndemlv, Rkfixdt*fixfraclv)     ! actual nitrogen fixation ending up in the laeves
       fixactst=Min(rndemst, Rkfixdt*fixfracst)     ! fixation could not be larger than demand
       fixactrt=Min(rndemrt, Rkfixdt*fixfracrt)
       fixactcrn=Min(rndemcrn, Rkfixdt*fixfraccrn)
       rndemlv=rndemlv-fixactlv    ! reduce the root uptake demand for leaves by N fixation going to the leaves
       rndemst=rndemst-fixactst
       rndemrt=rndemrt-fixactrt
       rndemcrn=rndemcrn-fixactcrn
       totdem=Max(0.d0, rndemlv+rndemst+rndemrt+rndemcrn)  ! update total demand
       Rkfixact=(fixactlv+fixactst+fixactrt+fixactcrn)/dt  ! total N fixation: Mass N/L2/T
       !  add the symbiotically fixed N to the organs (Masss N/L2)
       anlv=anlv+fixactlv
       anst=anst+fixactst
       anrt=anrt+fixactrt
       ancrn=ancrn+fixactcrn
    Else
       Rkfixact=0.d0
    Endif
    !Write(*,'(a,1p,99e13.5)') 'time,wlvg,wlvd,anlv,rndemlv,rndemst,rndemrt,totdem,tunc=',simtime%time,&
    !     PlantGeometry%wlvg,PlantGeometry%wlvd,anlv,&
    !     rndemlv,rndemst,rndemrt,totdem,tunc
    tunc=0.0
    Do sp=NH4,NO3
       Do i=1,NumNP
          hcsolo = Max(0.0,Min(Conc(sp,i),cRootMax(sp))*NFact(sp)) ! M / (L3 water) 
          unc(i,sp) = Max(0.0,Sink(i))*hcsolo*dt*delta_z(i) ! dt*dx*(unit of NitSinkRoot)
          tunc=tunc+unc(i,sp)
       Enddo
    Enddo
    If (tunc>totdem) Then
       unc(:,NH4)=unc(:,NH4)*totdem/tunc
       unc(:,NO3)=unc(:,NO3)*totdem/tunc
       und(:,NH4) =0.0
       und(:,NO3) =0.0
       tunc=totdem
       tund=0.0
    else
       tpdnup=totdem-tunc
       tund=0.d0
       Do sp=NH4,NO3
          Do i=1,NumNP
             ! from solute
             m=MatNum(i)
             Dw=ChPar(5,m,sp) !TODO: noch nicht Temperatur abhaengig
             TauW=ThNew(i)**(7./3.)/ThS(m)**2
             hcsolo = Max(0.0,Conc(sp,i)*NFact(sp))*ThNew(i) ! volume soil cRootMax?
             diffus_rm = Dw*TauW
             rdens=PlantGeometryArray%rrd(NumNP-i+1)*w0_dens
             und(i,sp)=rdens*rorad*2*PI*diffus_rm*hcsolo*dt*delta_z(i)/rdo ! 
             tund=tund+und(i,sp)
          enddo
       enddo
       if (tund.gt.tpdnup) then
          Do sp = NH4,NO3
             und(:,sp)=und(:,sp)*tpdnup/tund
          enddo
          tund=tpdnup
       endif
    Endif
    !     calculate the total uptake (mg/m**2) (wave)
    totvr=tunc+tund
    an_tot=an_tot+totvr     ! total seasonal uptake
    tunc_tot=tunc_tot+tunc  ! cumlulative convective uptake for output
    tund_tot=tund_tot+tund  ! cumlulative diffusive uptake for output
    If(totdem==0.0) Then
       rnuplv=0
       rnupst=0
       rnuprt=0
       rnupcrn=0
       !rnupso=0
    Else
       rnuplv=rndemlv*totvr/totdem
       rnupst=rndemst*totvr/totdem
       rnuprt=rndemrt*totvr/totdem
       rnupcrn=rndemcrn*totvr/totdem
       !rnupso=rndemso*totvr/totdem
    Endif
    ! calculate new masses in the different parts
    If(Abs(anlv+anst)<1.0e-20) Then
       anlv=anlv+rnuplv
       anst=anst+rnupst
    Else
!       anlv=anlv+rnuplv-rnupso*anlv/(anlv+anst)
!       anst=anst+rnupst-rnupso*anst/(anlv+anst)
       anlv=anlv+rnuplv-rndemso*anlv/(anlv+anst)
       If(anlv<0.0) Then
          Write(line,'(a,1p,e13.5,a)') 'time=',SimTime%time,&
               ' anlv<0, decrease xncso or increase xncle!'
          Call WriteError(Trim(line))
       Endif
       anst=anst+rnupst-rndemso*anst/(anlv+anst)
       If(anst<0.0) Then
          Write(line,'(a,1p,e13.5,a)') 'time=',SimTime%time,&
               ' anst<0, decrease xncso or increase xncst!'
          Call WriteError(Trim(line))
       Endif
    Endif
    anrt=Max(null_dp, anrt+rnuprt-PlantGeometry%ratwrtd*ancrt)
    Nit_ancrt=ancrt
    anso=anso+rndemso
    ancrn=ancrn+rndemcrn
    Call set_NitSinkBalance(rndemlv,rndemst,rndemrt,rndemso,rndemcrn,&
         tunc,totdem,tund,an_tot,anlv,anst,anrt,anso,ancrn,NitReduct,&
         anclv,ancst,ancrt,ancso,anccrn,Rkfixact,tunc_tot,tund_tot)
    ! calculate the uptake rates
    sinksum=0
    Do sp=NH4,NO3
       Do i=1,NumNP
          NitSinkRoot(i,sp)=(unc(i,sp)+und(i,sp))*UnitFactorPlants/(NFact(sp)*delta_z(i)*dt)
          if(sp==nh4) sinksum=sinksum+NitSinkRoot(i,sp)
       Enddo
    Enddo
  End Subroutine nit_upt


  Real(dp) Function nit_reduct(anclv, rlncl, rmncl)
    Use Plants, Only: PlantGeometry, Plant, Potatoe, SugarBeet
    Implicit None
    Real(dp), Intent(in) :: anclv, rlncl,rmncl
    Integer :: ip, ptype
    
    If(debug) Print *,'nit_reduct'
    ip=PlantGeometry%index
    If(Plant(ip)%season==0 .Or. ip==0 .Or. .Not. lNreduc) Then
       nit_reduct=1.0
    Elseif (lNitrogen) Then
       ptype=Plant(ip)%type   ! plant type
       If(ptype==Potatoe .Or. ptype==SugarBeet) Then
          If(PlantGeometry%temp_sum.Le.450.0 .And. PlantGeometry%plai .Lt. 0.75) Then
             nit_reduct=1.0
          else
             nit_reduct=Min(1.0, Max(0.0,((anclv-rlncl)/(rmncl-rlncl))))
          endif
       !        winter wheat, spring wheat and maize
       else
          If(PlantGeometry%dvs.Le.0.3 .And. PlantGeometry%plai.Lt.0.75) Then
             nit_reduct=1.0
          else
             nit_reduct=Min(1.0, Max(0.0,((anclv-rlncl)/(rmncl-rlncl))))
          endif
       endif
    Else
       nit_reduct=1.0
    Endif
  End Function nit_reduct


  Subroutine nit_correct_sink(dt)
    !###################################################################################
    !     in   : csolio, csolmo, decsoli, decsolm, dt, imobsw, ncs, nla,
    !            nr_of_sol, parasol, wcio, wcmo, rhuream,rnitrim,sunm,
    !            rhureai,rnitrii,suni
    !     out  : sinki, NitSink,rhuream,rnitrim,sunm,rhureai,rnitrii,suni
    !     calls:
    !###################################################################################
    Use Geometry, Only: MatNum, ThNew, delta_z
    Use Solute, Only: ChPar, Conc
    Implicit None
    Real(dp), Intent(in) :: dt
    Integer :: i,sp
    Real(dp) :: supply, demand, deficit
    Real(dp), Save :: sinksum(NSnit), sinkrootsum(NSnit)

    If(debug) Print *,'nit_correct_sink'
    Do sp = 1, NSnit
       ! the sinkterm
       Do i=1,NumNP
          supply=Conc(sp,i)*(ThNew(i)+ChPar(7,MatNum(i),sp)*ChPar(1,MatNum(i),sp))/dt
          If(NitSink(i,sp)<0.0d0) Then
             demand=-NitSink(i,sp)
             If(demand>supply) Then
                deficit=supply-demand
                supply=supply-deficit
                NitSink(i,sp)= Min(NitSink(i,sp)-deficit,0.d0)
                if(sp.eq.1) rhurea(i)=max(rhurea(i)+deficit,0.0d0)
                if(sp.eq.2) rnitri(i)=max(rnitri(i)+deficit,0.0d0)
                if(sp.eq.3) rdenit(i)=max(rdenit(i)+deficit,0.0d0)
             Endif
          !   If(sp.Eq.2) Write(*,'(1p,e13.5)') demand,supply,deficit,NitSink(i,2),rnitri(i)
          Endif
          If(lNRootUpt .And. NitSinkRoot(i,sp)>0.0d0) Then
             demand=NitSinkRoot(i,sp)
             If(demand>supply) demand=supply
             NitSinkRoot(i,sp)=-Max(demand,0.d0) ! negative means removal
          Endif
       Enddo
       ! accumulated over 1 timestep (day or hour), set to 0 in output
       sinksum(sp)=sinksum(sp)+Sum(NitSink(:,sp)*delta_z)
       sinkrootsum(sp)=sinkrootsum(sp)+Sum(NitSinkRoot(:,sp)*delta_z)
    Enddo
    Call  set_NitSinkBalance2(sinksum, sinkrootsum)
  End Subroutine nit_correct_sink

  !###################################################################################
  subroutine seqtrans(dt)
    !     in   : ncs, rdeniti, rdenitm, rhureai, rhuream, rmini,
    !            rminm, rnitrii, rnitrim, rvoli, rvolm, suni, sunm 
    !     out  : decnorg, decsoli, decsolm,
    !     calls: denitrification, hydrolysis, miner_immob, nitrification, rate_reduction,
    !            sol_sink, volatilisation
    !###################################################################################
    Use Geometry, Only: MatNum,ThNew, delta_z
    Use Material, Only: thR, thS
    Use Solute, Only: Conc
    Use Carbon, Only: sfact, massCO2conc, relgasD
    Implicit None
    Real(dp), Intent(in) :: dt
    Integer :: i,m
    Real(dp) :: rknitri_act, rkhyd_act, rkdenit_act, rkvol_act
    Real(dp) :: rsat, thetas, thetar, k_one, pn2n2o, pno3co2, no3conc, co2conc

    If(debug) Print *,'seqtrans'
    do i= 1,NumNP
       !        calculation of the different moisture threshold values
       m=MatNum(i)
       thetas=thS(m)
       thetar=thR(m)
       if (thetas .eq.0.d0) then
          Call WriteOutput('saturated moisture content is zero')
          Call WriteError('programme stopped : check err_file')
       endif
       rsat=min((ThNew(i)-thetar)/(thetas-thetar), 0.999999999999999d0)  !rsat=ThNew(i)/thetas
       ! original WAVE approach:
       if (denThresh .gt. 0.0) then
          if ((rsat-denThresh).lt.0.0001) then
             reddenit(i)=0.d0
          else
             reddenit(i)=((rsat-denThresh)/(1.0d0-denThresh))**2
          endif
       else
          reddenit(i)=0.50d0+(atan(0.60d0*3.14159*(.1d0*(rsat*100.0d0)+denThresh)))/3.14159  ! alternative function of del Grosso et al. 2000:
          !xinflect=9.0d0-(10.0d0*(0.145d0+(-1.175*min(0.113, relgasD(i)))))                 ! Das ist sehr unklar im Paper beschrieben
          !reddenit(i)=0.50d0+(atan(0.60d0*3.14159*(.1d0*(rsat*100.0d0)-xinflect)))/3.14159  ! Darum bleibt es erst mal raus
       endif   
       redwfps(i)=max(0.10d0, 1.50d0*rsat-0.320d0)  !N20 production reduction term, del Grosso et al. 2000
    enddo
    do i=1,NumNP
       rkhyd_act=sfact(1,i)*sfact(2,i)*rkhyd
       rhurea(i)=rkhyd_act*Conc(Urea,i)*ThNew(i)
    enddo 
    do i=1,NumNP
       rknitri_act=sfact(1,i)*sfact(2,i)*rknitri
       rnitri(i)=rknitri_act*Conc(NH4,i)*ThNew(i)
       ! N20 production from nitrification;  hole-in-the-pipe” approach (Firestone and Davidson, 1989).
       !rn2o(i)=rnitri(i)*0.0006d0*nh4_to_n2o   ! Mass N2O/(L3 soil)/T (Lin et al., 2000)
       !rnitri(i)=rnitri(i)*0.9994d0            ! reduce nitrification accordingly
       rn2o(i)=rnitri(i)*0.00105d0*nh4_to_n2o    ! Mass N2O/(L3 soil)/T (Florian Stange Diss)
       rnitri(i)=rnitri(i)*0.99895d0             ! reduce nitrification accordingly
    enddo
    do i=1,NumNP
       rkdenit_act=sfact(2,i)*reddenit(i)*rkdenit  ! sfact(2) accounts for temperature
       rdenit(i)=rkdenit_act*Conc(NO3,i)*ThNew(i)
       ! N20 production from denitrification, del Grosso et al., 2000
       co2conc=max(1.0e-15, massCO2conc(i))/dt/TimeFactor  ! CO2 production in M CO2/L-3 soil/day
       no3conc=Conc(NO3,i)*ThNew(i)          ! M NO3/L-3 soil
       pno3co2=no3conc/co2conc               ! no3 to co2 production ratio 
       k_one=max(1.70d0,38.40d0-350.0d0*relgasD(i))  ! relative gas Diffusion coeff is consistent with CO2 transport routine
       pn2n2o=1.0d0/(1.0d0+max(0.160d0*k_one, k_one*exp(-0.8*pno3co2))*redwfps(i))  !   Ratio N2 to N20 converted to multiplier (or fraction)
       rn2o(i)=rn2o(i)+rdenit(i)*min(0.999999999999999d0, pn2n2o)*no3_to_n2o  ! Mass N2O/(L3 soil)/T
       rn2(i)=rdenit(i)*(1.0d0-min(0.999999999999999d0, pn2n2o))*no3_to_n2    !Mass N2/(L3 soil)/T
       !write (*,*) i, rn2o(i), rn2(i), rdenit(i), k_one, pn2n2o, relgasD(i), redwfps(i), co2conc, no3conc, no3conc/co2conc
    enddo
    Do  i=1,NumNP
       rkvol_act=sfact(1,i)*sfact(2,i)*rkvol
       rvol(i)=rkvol_act*Conc(NH4,i)*ThNew(i)
    Enddo
    Call miner_immob(dt)
    !      urea
    NitSink(:,Urea)=-rhurea ! M/(L3 soil)/T
    !      ammonia
    NitSink(:,NH4)=rhurea*urea_to_nh4 - rnitri + rmin(:,NH4) - rvol ! M/(L3 soil)/T
    !     nitrate
    NitSink(:,NO3)=rnitri*nh4_to_no3 - rdenit + rmin(:,NO3) ! M/(L3 soil)/T
    !Write(*,'(1p,99e13.5)') Sum(nitsink(:,NH4)),Sum(nitsink(:,NO3)),Sum(rhurea),&
    !     Sum(rnitri),Sum(rdenit),Sum(rvol),Sum(rmin(:,NH4)),Sum(rmin(:,NO3))
    NitSinkBalance%rhurea=sum(rhurea*delta_z)
    NitSinkBalance%rnitri=sum(rnitri*delta_z)
    NitSinkBalance%rn2o=sum(rn2o*delta_z)
    NitSinkBalance%rn2=sum(rn2*delta_z)
    NitSinkBalance%rdenit=Sum(rdenit*delta_z)
    rvol=rvol*nh4_to_nh3 ! for output
    NitSinkBalance%rvol=Sum(rvol*delta_z)
    !Write(*,'(a,1p,999e13.5)') 'urea,nitsink,rmin,rnitri,rvol=',&
    !     Sum(delta_z*NitSink(:,Urea)),Sum(delta_z*NitSink(:,NH4)),Sum(delta_z*rmin(:,NH4)),&
    !     Sum(delta_z*rnitri),Sum(delta_z*rvol)
  end subroutine seqtrans

  Subroutine miner_immob(dt)
    Use datatypes
    Use Geometry, Only: NumNP, MatNum, ThNew
    Use Carbon, Only: PoolDPM, PoolRPM, sfact, co2parm, co2rdpm,co2rrpm,co2rhum,&
                      eff_rate_lit_from_nit,eff_rate_man_from_nit
    Use Variables, Only: lNitrogen
    Use Solute, Only: Conc

    Implicit None
    Real(dp), Intent(in) :: dt
    Integer :: i
    Real(dp) :: eff_rate_hum, eff_rate_lit, eff_rate_man
    Real(dp) :: check_mim, chno3, fe
    Real(dp) :: pxl, pxm, red_fact, shortage_min_n, vimma, vimmb
    Real(dp) :: change_nh4, change_no3 !, plus

    If(debug) Print *,'miner_immob'
    If(.Not.lNitrogen) Return
    !     pxl= ratio defining mineralisation or immobilisation for litter
    !     pxm= ratio defining mineralisation or immobilisation for manure
    !          > 0 = mineralisation
    !          < 0 = immobilisation
    Do i=1,NumNP
       red_fact = sfact(1,i)*sfact(2,i)*sfact(3,i)
       fe=co2parm(2,MatNum(i)) ! dimensionless
       !       If(PoolDPM(i)<1e-10) Then
       !          Print *,i,PoolNitLit(i),PoolDPM(i),PoolNitLit(i)/PoolDPM(i),PoolNitMan(i)/PoolRPM(i),fe/ro
       !       Endif
       pxl=PoolNitLit(i)/PoolDPM(i)-(fe/ro) ! dimensionless
       pxm=PoolNitMan(i)/PoolRPM(i)-(fe/ro) ! dimensionless
       !Print *,i,pxl,pxm
       eff_rate_lit = co2rdpm*red_fact ! 1/T
       eff_rate_man = co2rrpm*red_fact
       eff_rate_hum = co2rhum*red_fact
       ! unit = unit of Pool = M/(L3 soil)
       check_mim = ( eff_rate_hum*PoolNitHum(i)+ &
            pxm*eff_rate_man*PoolRPM(i)+ &
            pxl*eff_rate_lit*PoolDPM(i) )*dt ! M/(L3 soil)
       !If(i==1) Write(*,'(a,1p,99e13.5)') 'pxl,pxm,chek_mim=',pxl,pxm,check_mim,eff_rate_lit, eff_rate_man,eff_rate_hum,&
       !     PoolNitHum(i),PoolRPM(i),PoolDPM(i),red_fact, sfact(1,i),sfact(2,i),sfact(3,i)
       If(check_mim .Ge. 0.d0) Then  ! Mineralization, no immobilization
          change_nh4 = check_mim ! unit of Pool = M/(L3 soil)
          change_no3 = 0.d0
          eff_rate_lit = co2rdpm*red_fact ! 1/T
          eff_rate_man = co2rrpm*red_fact
       Else  ! Immobilization occurs
          !Print *,'miner_immob else1, node=',i
          If(Conc(NH4,i)*NFact(NH4)*ThNew(i) .Ge. (-check_mim)) Then  ! The present NH4 concentration is higher than mineraliztation N demand ! convert from NH4 to N
             change_nh4 = check_mim  ! The N mieralization demand is taken from liquid phase NH4 
             change_no3 = 0.d0       ! No NO3 from liquid phase required 
             eff_rate_lit = co2rdpm*red_fact  ! Lit pool decomposition rate not affected
             eff_rate_man = co2rrpm*red_fact  ! Man pool decomposition rate not affected
          else  ! The present NH4 concentration is not enough to fulfill mineralization N demand
          !Print *,'miner_immob else2, node=',i
             change_nh4 = -Conc(NH4,i)*NFact(NH4)*ThNew(i)  ! The present NH4 is entirely used for N mineralization
             chno3  = Conc(NH4,i)*NFact(NH4)*ThNew(i)+check_mim ! has unit of Pool; the NO3 required to fullfil demand is computed   
             If(Conc(NO3,i)*NFact(NO3)*ThNew(i) .Ge. (-chno3)) Then  ! enough NO3 available  !convert from NO3 to N
                change_no3 = chno3    ! The amount removed from liquid phase NO3 is set to the required NO3 
                eff_rate_lit = co2rdpm*red_fact  ! Lit pool decomposition rate not affected
                eff_rate_man = co2rrpm*red_fact  ! Man pool decomposition rate not affected
             else   ! NOT enough NO3 available
          !Print *,'miner_immob else3, node=',i
                change_no3 = -Conc(NO3,i)*NFact(NO3)*ThNew(i)   ! The present NO3 is entirely used for N mineralization 
                shortage_min_n = Conc(NO3,i)*NFact(NO3)*ThNew(i)+chno3   ! The NO3 gap is estimated  
                if(pxl .ge. 0.d0  .and. pxm .lt. 0.d0) then   ! Lit pool has enough N, Man pool is missing N for mineralization
                   vimmb    = pxm*eff_rate_man*PoolRPM(i)*dt  ! this is the amount of N potentially decomposed from Man(RPM) pool ! unit of Pool = M/(L3 soil)
                   vimma    = vimmb-shortage_min_n      ! this is the amount of N that could ACTUALLY be decomposed from Man(RPM) pool accounting for the NO3 gap
                   eff_rate_lit = co2rdpm*red_fact      ! Lit pool decomposition rate not affected
                   eff_rate_man = vimma/(pxm*PoolRPM(i)*dt) ! Man(RPM) decomp. rate is computed as ratio between N that could actually be decomposed and amount of N potentially decomposed    
                elseif(pxl .lt. 0.d0  .and.pxm .ge. 0.d0) then ! Lit(DPM) pool is missing N for mineralization, Man pool has enough N
          !Print *,'miner_immob else4, node=',i
                   vimmb    = pxl*eff_rate_lit*PoolDPM(i)*dt  ! this is the amount of N potentially decomposed from Lit(DPM) pool
                   vimma    = vimmb-shortage_min_n   ! this is the amount of N that could ACTUALLY be decomposed from Lit(DPM) pool accounting for the NO3 gap
                   eff_rate_lit = vimma/(pxl*PoolDPM(i)*dt) ! Lit(DPM) decomp. rate is computed as ratio between N that could actually be decomposed and amount of N potentially decomposed
                   eff_rate_man = co2rrpm*red_fact          ! Man pool decomposition rate not affected  ! 1/T  
                elseif(pxl .lt. 0.d0 .and. pxm .lt. 0.d0) then  ! Both, Lit(DPM) and Man(RPM) pools are missing N for mineralization
          !Print *,'miner_immob else5, node=',i
                   vimmb    = pxl*eff_rate_lit*PoolDPM(i)*dt  ! this is the amount of N potentially decomposed from Lit(DPM) pool
                   vimma    = vimmb-shortage_min_n   ! this is the amount of N that could ACTUALLY be decomposed from Lit(DPM) pool accounting for the NO3 gap
                   if(vimma .gt. 0.d0) then   ! Not?? Enough NO3 to fulfill Lit(DPM) decomposition N demand, may be there is enough NO3 left for Man decomposition?
                      vimmb    = pxm*eff_rate_man*PoolRPM(i)*dt ! this is the amount of N potentially decomposed from Man(RPM) pool ! unit of Pool
                      vimma    = vimmb+vimma ! ??
                      eff_rate_lit = 0.d0  ! No Lit(DPM) decomposition
                      eff_rate_man = vimma/(pxm*PoolRPM(i)*dt) ! Man decomp. rate is reduced ??
                   else ! NOT enough NO3 to fulfill Lit(DPM) decomposition N demand
                      Print *,'troubles in rminimm1'
                      eff_rate_lit = 0.0_dp
                      eff_rate_man = 0.0_dp
                   endif
                else
                   Print *,'troubles in rminimm2'
                   eff_rate_lit = 0.0_dp
                   eff_rate_man = 0.0_dp
                endif
             endif
          endif
       Endif
       If(.Not.LNPools) Then
          eff_rate_lit = 0
          eff_rate_man = 0
          eff_rate_hum = 0
       Endif
       !If(eff_rate_lit>co2rdpm*red_fact) Print *,'eff_rate_lit>co2rdpm',eff_rate_lit*red_fact,co2rdpm
       ! transfer to rates to carbon
       eff_rate_lit_from_nit(i)=eff_rate_lit
       eff_rate_man_from_nit(i)=eff_rate_man
       If(eff_rate_lit<null_dp .Or. eff_rate_man<null_dp) Print *,'eff_rate_lit<null_dp'
       !       Print *,pxl,pxm,fe,fh,ro,eff_rate_lit,eff_rate_man
       !        the gain/losses of the nitrogen litter pool
       !        through mineralisation/immobilisation and nitrogen humification
       PoolNitLit(i)=PoolNitLit(i)+(-pxl-fe*fh/ro)*eff_rate_lit*PoolDPM(i)*dt  ! unit of Pool = M/(L3 soil)
       !plus=(-pxl-fe*fh/ro)*eff_rate_lit*PoolDPM(i)*dt  ! unit of Pool = M/(L3 soil)
       !        the gain/losses of the nitrogen litter pool
       !        through mineralisation/immobilisation and nitrogen humification
       PoolNitMan(i)=PoolNitMan(i)+(-pxm-fe*fh/ro)*eff_rate_man*PoolRPM(i)*dt
       !plus=plus+(-pxm-fe*fh/ro)*eff_rate_man*PoolRPM(i)*dt
       !        the changes in the n-humus pool through mineralisation and n-litter humification
       PoolNitHum(i)=PoolNitHum(i)+(fe*fh/ro*(eff_rate_lit*PoolDPM(i)+ &
            eff_rate_man*PoolRPM(i))-eff_rate_hum*PoolNitHum(i))*dt
       !plus=plus+(fe*fh/ro*(eff_rate_lit*PoolDPM(i)+ &
       !     eff_rate_man*PoolRPM(i))-eff_rate_hum*PoolNitHum(i))*dt
       !        retain the rate of change of the ammonia and nitrate species
       !        to be incorporated in the sequential transformation process  (mg/(day*liter))

       ! convert from N to NH4/NO3
       If(LNPools) Then
          rmin(i,NH4)=change_nh4/NFact(NH4)/dt ! M/(L3 soil)/T
          rmin(i,NO3)=change_no3/NFact(NO3)/dt
       Else
          rmin(i,NH4)=0
          rmin(i,NO3)=0
       Endif
       !Write(*,'(1p,99e13.5)') sfact(1,i),sfact(2,i),sfact(3,i),change_nh4,-plus
    Enddo
  End Subroutine miner_immob

  Subroutine NITinpFromOC()
    Use Geometry, Only: OCinpLit, OCinpMan, rnodeResN, rnodedeadwN
    Implicit None
    
    If(debug) Print *,'NITinpFromOC'
    If(lNitrogen) Then
       ! add organic C input and harvest residues
       PoolNitLit=PoolNitLit+OCinpLit*C2NFact+(rnodeResN + rnodedeadwN)*0.59 
       PoolNitMan=PoolNitMan+OCinpMan*C2NFact+(rnodeResN + rnodedeadwN)*0.41
    Endif
  End Subroutine NITinpFromOC

  Subroutine WriteNitrogenOutput(SimTime)
    Use TimeData, Only: iday_to_date, simtime_to_iday
    Implicit None
    Real(dp), Intent(in) :: SimTime
    Logical, Save :: firstCall=.True.
    Integer :: ios,yy,mm,dd
    Character(len=*), Parameter :: fnplant = "plantupt.out"
    Character(len=*), Parameter :: fnnit = "nitrogen.out"

    If(debug) Print *,'WriteNitrogenOutput'
    If(.Not.lNitrogen) Return
    If(firstCall) Then
       firstCall=.false.
       Open(64,file=Trim(fnplant),Status='REPLACE',ACTION='WRITE',&
            FORM='FORMATTED',iostat=ios)
       If(ios/=0) Call WriteError( 'Error while opening nitrogen output file: '&
            // Trim(fnplant))
       Write(64,2) '#nitrogen'
       Write(64,2) '# 1=year','# 2=month','# 3=day','# 4=simulation time',&
            '# 5=rndemlv ( (mass Nitrogen)/L2/T )','# 6=rndemst ( (mass Nitrogen)/L2/T )',&
            '# 7=rndemrt ( (mass Nitrogen)/L2/T )','# 8=rndemso ( (mass Nitrogen)/L2/T )',&
            '# 9=rndemcrn ( (mass Nitrogen)/L2/T )',&
            '#10=totdem ( (mass Nitrogen)/L2/T )','#11=ctunc ( (cum mass Nitrogen)/L2 )',&
            '#12=ctund ( (cum mass Nitrogen)/L2 )',&
            '#13=an_tot ( (mass Nitrogen)/L2 )','#14=anlv ( (mass Nitrogen)/L2 )',&
            '#15=anst ( (mass Nitrogen)/L2 )','#16=anrt ( (mass Nitrogen)/L2 )',&
            '#17=anso ( (mass Nitrogen)/L2 )','#18=ancrn ( (mass Nitrogen)/L2 )',&
            '#19=NitReduct (-)',&
            '#20=anclv ( (mass Nitrogen)/(mass dry matter) )',&
            '#21=ancst ( (mass Nitrogen)/(mass dry matter) )',&
            '#22=ancrt ( (mass Nitrogen)/(mass dry matter) )',&
            '#23=ancso ( (mass Nitrogen)/(mass dry matter) )',&
            '#24=anccrn ( (mass Nitrogen)/(mass dry matter) )',&
            '#25=rkfixact ( (mass Nitrogen)/L2/T )'
       Write(64,'(A,A,99A13)') '#     date','    time',&
            'rndemlv','rndemst','rndemrt','rndemso','rndemcrn',&
            'totdem','ctunc','ctund','tottup','anlv','anst','anrt','anso','ancrn',&
            'NitReduct','anclv','ancst','ancrt','ancso','anccrn','rkfixact'
       ! nitrogen
       Open(65,file=Trim(fnnit),Status='REPLACE',ACTION='WRITE',&
            FORM='FORMATTED',iostat=ios)
       If(ios/=0) Call WriteError( 'Error while opening nitrogen output file: '&
            // Trim(fnnit))
       Write(65,2) '#plantupt'
       Write(65,2) '# 1=year','# 2=month','# 3=day','# 4=simulation time',&
            '# 5=NitSink(NH4)  (M/(L2 soil)/T)',&
            '# 6=NitSink(NO3)  (M/(L2 soil)/T)',&
            '# 7=NitSinkRoot(NH4)  (M/(L2 soil)/T)',&
            '# 8=NitSinkRoot(NO3)  (M/(L2 soil)/T)',&
            '# 9=rhurea (M/L2/T)',&
            '#10=rnitri (M/L2/T)',&
            '#11=rdenit (M/L2/T)',&
            '#12=rvol (M/L2/T)',&
            '#13=rN2O (M/L2/T)',&
            '#14=rN2 (M/L2/T)'
       Write(65,'(A,A,99A13)') '#     date','    time',&
            'NitSink(NH4)','NitSink(NO3)','NSkRoot(NH4)','NSkRoot(NO3)',&
            'rhurea','rnitri','rdenit','rvol', 'rN2O', 'rN2'
    Else
       Open(64,file=fnplant,Status='OLD',ACTION='WRITE',&
            FORM='FORMATTED',POSITION='APPEND')
       Open(65,file=fnnit,Status='OLD',ACTION='WRITE',&
            FORM='FORMATTED',POSITION='APPEND')
    Endif
    Call iday_to_date(simtime_to_iday(SimTime),yy,mm,dd)
    Write(64,'(I4,2I3,I8,1P,99E13.5)') yy,mm,dd, Nint(SimTime),&
         NitSinkBalance%rndemlv,NitSinkBalance%rndemst,NitSinkBalance%rndemrt,&
         NitSinkBalance%rndemso,NitSinkBalance%rndemcrn,&
         NitSinkBalance%totdem,NitSinkBalance%tunc_tot,NitSinkBalance%tund_tot,NitSinkBalance%an_tot,&
         NitSinkBalance%anlv,NitSinkBalance%anst,NitSinkBalance%anrt,&
         NitSinkBalance%anso,NitSinkBalance%ancrn,&
         NitSinkBalance%NitReduct,NitSinkBalance%anclv,NitSinkBalance%ancst,&
         NitSinkBalance%ancrt,NitSinkBalance%ancso,NitSinkBalance%anccrn,NitSinkBalance%Rkfixact
    Write(65,'(I4,2I3,I8,1P,99E13.5)') yy,mm,dd, Nint(SimTime),&
         NitSinkBalance%sinksum(NH4:NO3),NitSinkBalance%sinkrootsum(NH4:NO3),&
         NitSinkBalance%rhurea,NitSinkBalance%rnitri,&
         NitSinkBalance%rdenit,NitSinkBalance%rvol,&
         NitSinkBalance%rn2o, NitSinkBalance%rn2
    NitSinkBalance%sinksum(NH4:NO3)=null_dp
    NitSinkBalance%sinkrootsum(NH4:NO3)=null_dp
2   Format(A)
  End Subroutine WriteNitrogenOutput
 
  Subroutine NitrogenTop(tnew, Prec, rSoil, cTop)
    Implicit None
    Integer, Intent(in) :: tnew
    Real(dp), Intent(in) :: Prec, rSoil
    Real(dp), Intent(out) :: cTop(:)
    Integer :: i
    Real(dp) :: nConc, vol
    
    If(debug) Print *,'NitrogenTop'
    ! boundaries for nitrogen
    If(.Not.lNitrogen) Return
    cTop(1:NSnit)=0
    !cBot(1:NSnit)=0
    !Print *,'ntop',tnew,NTop_ind,NTop_dim,NTopTime
    If(NTop_ind<=NTop_dim) Then
       If(tnew==NTopTime(NTop_ind)) Then
          N_on_top=N_on_top+NTopMass(:,NTop_ind)
          NTop_ind=NTop_ind+1
       Endif
    Endif
    ! waiting for a flow into the soil
    vol=Prec-rSoil
    Do i=1,NSnit
       If(N_on_top(i)>0 .And. vol>0.05/Units_L(Unit_L_Input)*Units_T(Unit_T_Input)) Then
          ! more than 0.05 mm/h
          nConc=N_on_top(i)/vol
          cTop(i)=nConc/NFact(i) ! split N to Urea, NH4 and NO3 and convert to conc
          Write(*,'(a,i8,a,1p,e13.5,a,e13.5,a,e13.5)')&
               'nitrogen on top, time=',tnew,' mass=',N_on_top(i),' flow=',vol,' conc=',nConc
          N_on_top(i)=0
       Endif
    Enddo
  End Subroutine NitrogenTop

  Subroutine set_NitSinkBalance(rndemlv,rndemst,rndemrt,rndemso,rndemcrn,tunc,totdem,tund,an_tot,&
       anlv,anst,anrt,anso,ancrn,NitReduct,anclv,ancst,ancrt,ancso,anccrn,Rkfixact,&
       tunc_tot,tund_tot)
    Real(dp), Intent(in) :: rndemlv,rndemst,rndemrt,rndemso,rndemcrn,&
         tunc,totdem,tund,an_tot,anlv,anst,anrt,anso,ancrn,NitReduct,anclv,&
         ancrt,ancst,ancso,anccrn,Rkfixact,tunc_tot,tund_tot

    If(debug) Print *,'set_NitSinkBalance'
    NitSinkBalance%rndemlv=rndemlv
    NitSinkBalance%rndemst=rndemst
    NitSinkBalance%rndemrt=rndemrt
    NitSinkBalance%rndemso=rndemso
    NitSinkBalance%rndemcrn=rndemcrn
    NitSinkBalance%tunc=tunc
    NitSinkBalance%totdem=totdem
    NitSinkBalance%tund=tund
    NitSinkBalance%an_tot=an_tot
    NitSinkBalance%anlv=anlv
    NitSinkBalance%anst=anst
    NitSinkBalance%anrt=anrt
    NitSinkBalance%anso=anso
    NitSinkBalance%ancrn=ancrn
    NitSinkBalance%NitReduct=NitReduct
    NitSinkBalance%anclv=anclv
    NitSinkBalance%ancrt=ancrt
    NitSinkBalance%ancst=ancst
    NitSinkBalance%ancso=ancso
    NitSinkBalance%anccrn=anccrn
    NitSinkBalance%Rkfixact=Rkfixact
    NitSinkBalance%tunc_tot=tunc_tot
    NitSinkBalance%tund_tot=tund_tot
  End Subroutine set_NitSinkBalance

  Subroutine reset_NitSinkBalance()
    If(debug) Print *,'reset_NitSinkBalance'
    NitSinkBalance%rndemlv=0
    NitSinkBalance%rndemst=0
    NitSinkBalance%rndemrt=0
    NitSinkBalance%rndemso=0
    NitSinkBalance%rndemcrn=0
    NitSinkBalance%tunc=0
    NitSinkBalance%totdem=0
    NitSinkBalance%tund=0
    NitSinkBalance%an_tot=0
    NitSinkBalance%anlv=0
    NitSinkBalance%anst=0
    NitSinkBalance%anrt=0
    NitSinkBalance%anso=0
    NitSinkBalance%ancrn=0
    NitSinkBalance%NitReduct=1
    NitSinkBalance%anclv=0     
    NitSinkBalance%ancrt=0
    NitSinkBalance%ancst=0
    NitSinkBalance%ancso=0
    NitSinkBalance%anccrn=0
    NitSinkBalance%sinksum=0
    NitSinkBalance%sinkrootsum=0
    NitSinkBalance%tunc_tot=0
    NitSinkBalance%tund_tot=0
  End Subroutine reset_NitSinkBalance
  
  Subroutine set_NitSinkBalance2(sum1, sum2)
    Real(dp), Intent(in), Dimension(:) :: sum1, sum2
    If(debug) Print *,'set_NitSinkBalance2'
    NitSinkBalance%sinksum=sum1
    NitSinkBalance%sinkrootsum=sum2
  End Subroutine set_NitSinkBalance2

End Module Nitrogen
