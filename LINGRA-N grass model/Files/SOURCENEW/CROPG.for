*------------------------------------------------------------------------- *
*  SUBROUTINE CROPG                                                        *
*  Author: Joost Wolf                                                      *
*  Date: Grass rop model developed on the basis of LINGRA.fst in March 2012* 
*  Purpose: This subroutine simulates the dry matter increase of a         *
*           grass crop as function of intercepted radiation, temperature,  *
*           radiation use efficiency, water and nitrogen availability.     *
*           Growth rate is lowest of both source and of sink limited growth*
*                                                                          *
*                                                                          *
*  FORMAL PARAMETERS:(I= input, O= output, C= control, IN= init., T-time)  *
*  name     meaning                                  units       class     *
*  ----     -------                                  -----       -----     *
*  ICROP    number of grass data set                  -           I        *
*  ISOIL    number of soil type                       -           I        *
*  INITI    indicates initialization of run           -           I,C      *
*  IOPT     indicates optimal (=1), water limited (=2)                     *
*            or water and N limited run (=3)          -           I,C      *
*  IDAY     julian day number                         -           I        *
*  IDEM     date of emergence (1 for permament grass) -           I        *
*  IDEMERG  date of emergence (1 for permanent grass) -           O        *
*  IDPL     date of planting/sowing                   -           I        *
*  IDFLOW   date of flowering                                     O        *
*  IDHALT   end date of crop growth                               O        *
*  PL       indicates planting as start of simulation -           I,C      *
*  TERMIN   indicates terminal section                -           I,C      *
*  EMERG    indicates crop emergence                  -           I,C      *
*  TMIN     minimum air temperature                   C           I        *
*  TMAX     maximum air temperature                   C           I        *
*  AVRAD    daily total irradiation                   J m-2 d-1   I        *
*  CO       atmospheric CO2 concentration             ppmv        I        *
*  TRANRF   reduction factor due to drought/wetness   -           I        *  
*  RDMSO    soil related maxiumum rooting depth       cm          I        *
*  DAYLP    photoperiodically active daylength        h           I        *
*  WLVG     weight of living leaves                   kg DM ha-1  O        *
*  WLVDD    weight of dead leaves                     kg DM ha-1  O        *
*  WST      weight of stems                           kg DM ha-1  O        *
*  WRT      weight of roots                           kg DM ha-1  O        *
*  WRE      weight of reserves (storage carbohydrates) kg DM ha-1 O        *
*  GRASS    weight of cutted grass                     kg DM ha-1 O        *
*  YIELD    weight of cutted grass and harvestable grass kg DM ha-1 O      *
*  TADRW    weight of cutted grass and living/dead leaves kg DM ha-1 O     *
*  TILLER   number of tillers                         m-2         O        *
*  RD       actual rooting depth                      cm          O        *
*  RDMCR    crop specific maximum rooting depth       cm          O        *
*  RR       root growth rate                          cm d-1      O        *
*  RDM      soil/crop related maximal rooting depth   cm          I        *
*  LAI      leaf area index                           m2 m-2      O        *
*  DEPNR   crop group number for soil water depletion -           O        *
*  CFET    crop specific correction for transpiration -           O        * 
*  IAIRDU  air ducts in roots present (=1) or not(=0)  -          O,C      *
*  TSUM     temperature sum from emergence            C d         O        *
*  DVS      development stage                         -           O        *
*  DVSEND   development stage at end of growth period -           O        *
*  TSULP    temperature sum from sowing/planting      C d         O        *
*  FINT     fractional light interception             -           O        *
*  TPARINT  total intercepted radiation (PAR)         MJ m-2      O        *
*  TPAR     total photosynthetically active radiation MJ m-2      O        *
*  TSUML    temperature sum from emergence incl. dayl.effect C.d  O        *
*  NNI      nitrogen nutrition index                  -           O        *
*  NMINT    total mineral N from soil and fertiliser  kg N ha-1   O        *
*  NMIN     mineral N available from soil for crop    kg N ha-1   O        *
*  NUPTT    total N uptake by crop from soil          kg N ha-1   O        *
*  NLIVT    Amount of N in living crop organs         kg N ha-1   O        *
*  NLOSST   Amount of N in dead crop organs incl. harvests kg N ha-1 O     *
*  YCH      indicates year change                     -           I        * 
*                                                                          *
*  File usage: LINTER                                                      *
*------------------------------------------------------------------------
       SUBROUTINE CROPG(ICROP,ISOIL,INITI,IOPT, IDAY,IDEM,IDEMERG,IDPL,
     $    IDFLOW,IDHALT,PL,TERMIN,EMERG,TMIN,TMAX,AVRAD,CO,TRANRF,RDMSO,
     $	  DAYLP,WLVG, WLVDD, WRT, WRE,GRASS,YIELD,TADRW,TILLER,
     $       RD,RDMCR,RR,RDM,LAI,
     $        CFET,DEPNR,IAIRDU,TSUM,DVS,DVSEND,TSULP,FINT,TPARINT,TPAR,
     $        TSUML,NNI,NMINT,NMIN,NUPTT,NLIVT,NLOSST,YCH)


       IMPLICIT REAL (A-Z)
       INTEGER ICROP,IDAY,IDEM,IDEMERG,IDPL,IDFLOW,IDHALT,IDSL
       INTEGER ILCO,IOPT,IAIRDU,ILDTSM,ILSLA,ILSSA,ILKDIF,ILRUE
	 INTEGER ILTMPF,ILTMNF,ILFR,ILFL,ILFS,ILFO,ILRDRL, ILRDD,ISOIL
	 INTEGER ILRDRR,ILRDRS,ILNMXL,ILPHOT,ILFERT,ILNRFT,ILFRT,INH
	 INTEGER IMOPT

       LOGICAL TERMIN,INITI,PL,FLOW,EMERG,YCH
       REAL DTSMTB (30), SLATB (30),SSATB (30),KDIFTB (30)  
       REAL TMPFTB (30), TMNFTB (30), RUETB (30), FLTB (30), RDRLTB (30) 
       REAL COTB (30), RDRRTB (30), RDRSTB (30), NMXLV (30), PHOTTB (30)
	 REAL FSTB (30), FOTB (30), FERTAB(30), NRFTAB(30), VERNRT(30)
	 REAL FRTB (30), MNDAT (30), FRTRB (30), RDDTB (30) 

       COMMON /CROPOUT/ WLVG1,WLVDD1,WRE1,GRASS1,YIELD1,TADRW1
	 COMMON /CROPOUT/ TPARINT1, RUEC
	 COMMON /CROPOUT/ NUPTT1,NLIV1,NLOSS1,NROOT1,TRANRF1,NNI1


*----------------------------------------------------
*      INITIALIZATION
*----------------------------------------------------

       IF (INITI) THEN
         FLOW= .FALSE.
         EMERG= .FALSE.
       ELSE
         GOTO 200
       ENDIF

*---- Reading crop data from file CROPG.DAT
       
       IF (ICROP .EQ. 1) CALL RDINIT (14,0, 'CROPG1.DAT')
       IF (ICROP .EQ. 2) CALL RDINIT (14,0, 'CROPG2.DAT')
       IF (ICROP .EQ. 3) CALL RDINIT (14,0, 'CROPG3.DAT')
       IF (ICROP .EQ. 4) CALL RDINIT (14,0, 'CROPG4.DAT')
       IF (ICROP .EQ. 5) CALL RDINIT (14,0, 'CROPG5.DAT')
       IF (ICROP .EQ. 6) CALL RDINIT (14,0, 'CROPG6.DAT')

      
!    Read initial states and parameter values
	call RDSINT ('IDSL', IDSL)
      CALL RDSREA ('TSUM1',TSUM1)
      CALL RDSREA ('TSUM2',TSUM2)

      CALL RDSREA ('DVSI', DVSI)
      CALL RDSREA ('DVSEND',DVSEND)
      CALL RDSREA ('TDWI',TDWI)
      CALL RDSREA ('RGRLAI',RGRLAI)
      
      CALL RDSREA ('SPA',SPA)
      CALL RDAREA ('DTSMTB',DTSMTB,30,ILDTSM)      
      CALL RDAREA ('SLATB',SLATB,30,ILSLA)
      CALL RDAREA ('SSATB',SSATB,30,ILSSA)
      CALL RDSREA ('TBASE',TBASE)
	CALL RDAREA ('KDIFTB',KDIFTB,30,ILKDIF)
      CALL RDAREA ('RUETB',RUETB,30,ILRUE)
      CALL RDAREA ('TMPFTB',TMPFTB,30,ILTMPF)
      CALL RDAREA ('TMNFTB',TMNFTB,30,ILTMNF)
      CALL RDAREA ('COTB',COTB,30,ILCO)
      CALL RDAREA ('FRTB',FRTB,30,ILFR)
      CALL RDAREA ('FLTB',FLTB,30,ILFL)
      CALL RDAREA ('FSTB',FSTB,30,ILFS)
      CALL RDAREA ('FOTB',FOTB,30,ILFO)
      CALL RDSREA ('RDRL',RDRL)
      CALL RDAREA ('RDRLTB',RDRLTB,30,ILRDRL)
      call RDSREA ('RDRSHM', RDRSHM)
      call RDSREA ('RDRNS', RDRNS)
      CALL RDAREA ('RDRRTB',RDRRTB,30,ILRDRR)
      CALL RDAREA ('RDRSTB',RDRSTB,30,ILRDRS)
      CALL RDSREA ('CFET',CFET)
      CALL RDSREA ('DEPNR',DEPNR)
      CALL RDSINT ('IAIRDU',IAIRDU)
      CALL RDSREA ('RDI',RDI)
      CALL RDSREA ('RRI',RRI)
      CALL RDSREA ('RDMCR',RDMCR)
      CALL RDSREA ('DVSDLT', DVSDLT)
      CALL RDSREA ('DVSNLT', DVSNLT)
      CALL RDSREA ('DVSNT', DVSNT)
      CALL RDSREA ('TBASEM',TBASEM)
      CALL RDSREA ('TEFFMX',TEFFMX)
      CALL RDSREA ('TSUMEM',TSUMEM)
      call RDSREA ('FNTRT', FNTRT)
      call RDSREA ('FRNX', FRNX)
      call RDSREA ('LAICR', LAICR)
      call RDSREA ('LRNR', LRNR)
      call RDSREA ('LSNR', LSNR)
      call RDSREA ('NLAI', NLAI)
      call RDSREA ('NLUE', NLUE)
      call RDSREA ('NMAXSO', NMAXSO)
      call RDSREA ('NPART', NPART)
      call RDSREA ('NFIXF', NFIXF)
      call RDSREA ('NSLA', NSLA)
      call RDSREA ('RNFLV', RNFLV)
      call RDSREA ('RNFRT', RNFRT)
      call RDSREA ('RNFST', RNFST)
      call RDSREA ('TCNT', TCNT) 
      call RDAREA ('NMXLV', NMXLV, 30, ILNMXL)
      call RDAREA ('PHOTTB', PHOTTB, 30, ILPHOT)
      call RDAREA ('RDDTB', RDDTB, 30, ILRDD)
      call RDAREA ('FRTRB', FRTRB, 30, ILFRT)
      call RDSREA ('TILLI', TILLI)      
	call RDSREA ('WREI', WREI)
	call RDSREA ('TMBAS1', TMBAS1)	 

      CLOSE (14, STATUS= 'DELETE')

*     Read management data from file MANAGG.DAT
      Call RDINIT (15,0, 'MANAGG.DAT')
      call RDAREA ('FERTAB', FERTAB, 30, ILFERT) 
      call RDAREA ('NRFTAB', NRFTAB, 30, ILNRFT)
      call RDAREA ('MNDAT', MNDAT, 30, INH)
	      
	call RDSREA ('NMINS', NMINS)
	call RDSREA  ('RTMINS', RTMINS)
	call RDSREA ('CLAI', CLAI)
	call RDSINT ('IMOPT', IMOPT)
	call RDSREA ('CWGHT', CWGHT)

      CLOSE (15, STATUS= 'DELETE')

*---- Reading soil data from file SOILG.DAT
       
       IF (ISOIL .EQ. 1) CALL RDINIT (21,0, 'SOILG1.DAT')
       IF (ISOIL .EQ. 2) CALL RDINIT (21,0, 'SOILG2.DAT')
       IF (ISOIL .EQ. 3) CALL RDINIT (21,0, 'SOILG3.DAT')
       IF (ISOIL .EQ. 4) CALL RDINIT (21,0, 'SOILG4.DAT')
       IF (ISOIL .EQ. 5) CALL RDINIT (21,0, 'SOILG5.DAT')
       IF (ISOIL .EQ. 6) CALL RDINIT (21,0, 'SOILG6.DAT')
  
	call RDSREA ('SOITMI', SOITMI)
      CLOSE (21, STATUS= 'DELETE')


*----- initialization of state variables

       TDW= TDWI
	 DVS= DVSI
	 DVR= 0.
	 NMINI=NMINS
	 NMIN=NMINI
	 WRTI=  LINT (FRTB, ILFR, DVSI) * TDW
	 WRT= WRTI
	 TAGB= TDW- WRT 
	 WLVGI= LINT (FLTB, ILFL, DVSI) * TAGB
	 WLVG= WLVGI
	 SLAI= LINT(SLATB,ILSLA,DVSI)
	 LAII= WLVGI * SLAI 
	 LAI= LAII
	 RLAI= 0.
	 WLVD= 0.
	 WSTI=  LINT (FSTB, ILFS, DVSI) * TAGB
	 WST= WSTI 
       WSOI=  LINT (FOTB, ILFO, DVSI) * TAGB
	 WSO= WSOI
	 WSTD= 0.
	 WRTD= 0.
	 RWLVG= 0.
	 RWST=0.
	 RWRT=0.
	 RWSO= 0.
	 DLV=0.
	 DRST=0.
	 DRRT= 0.
       GTW= 0.
       TSULP= 0. 
       TSUM= 0.
	 TSUML= 0.
       DTSULP= 0. 
       DTSUM= 0.
	 DTSUML= 0.
	 VERN= 0.
	 RVERNR= 0.
	 DVRED= 1.
       TPAR= 0.
       TPARINT= 0.
       PAR= 0.
       PARINT= 0. 
	 ATN= 0. 
	 GTSUM= 0.  
	 RD= RDI 
	 DELT= 1. 
	 NMINT= 0.
	 NMAXLVI= LINT (NMXLV, ILNMXL, DVSI)
       NMAXSTI= LSNR * NMAXLVI
	 NMAXRTI= LRNR * NMAXLVI
       ANLVI= NMAXLVI * WLVGI
	 ANSTI= NMAXSTI * WSTI
	 ANRTI= NMAXRTI * WRTI
	 ANLV= ANLVI
	 ANST= ANSTI
	 ANRT= ANRTI
	 ANSOI= 0.
	 ANSO= ANSOI
	 NLOSSL= 0.
	 NLOSSR= 0.
	 NLOSSS= 0.
	 NUPTT= 0.
	 RNMINS= 0.
	 RNMINT= 0.
	 CTRAN= 0.
       CNNI= 0.
	 INCUT= 0.
	 DAHA= 0.
	 RDAHA= 0.
	 RSOITM= 0.
	 SOITMP= SOITMI
	 TILLER= TILLI
	 DTIL= 0.
	 GRASS= 0.
	 CWLVG= CLAI/SLAI
	 HARV= 0.
	 WRE= WREI
	 RRE= 0.
	 GRE= 0.
	 DRE= 0.
	 TADRW= 0.
       YIELD= 0.
       LENGTH= 0.

       DAYPL= REAL(IDPL)
	 TRANRF= 1.
	 NNI= 1.
       
       IF (.NOT. PL) THEN
          DAYEM= REAL(IDEM)
	    EMERG= .TRUE.
       ELSE
          DAYEM= 999.
       ENDIF


200     CONTINUE

        


        IF (TERMIN) GOTO 999


*---------------------------------------------------------
*      INTEGRATION
*---------------------------------------------------------

*----- Temperature sums (C.d) from sowing/planting (P) and from emergence with and without daylength effect 

       TSULP= INTGRL(TSULP, DTSULP, 1.)
       TSUM= INTGRL(TSUM, DTSUM, 1.)
	 TSUML= INTGRL(TSUML, DTSUML, 1.)

*----- Vernalisation effect accumulated
       VERN= INTGRL(VERN, RVERNR, 1.)

*----- Start of flowering
*       IF (.NOT. FLOW .AND. TSUML .GE. TSUM1)  IDFLOW= IDAY
*       IF (.NOT. FLOW .AND. TSUML .GE. TSUM1)  FLOW= .TRUE.

*---- Total photosynthetically active radiation (MJ/m2)
      TPAR= INTGRL(TPAR, PAR, 1.)

*---- Total intercepted radiation (MJ/m2)
      TPARINT= INTGRL(TPARINT, PARINT, 1.)

*---- Dry weights of total biomass, living crop organs, and total above-ground living biomass (kg DM/ha)
	GTSUM= INTGRL(GTSUM, GTW, 1.)
      WLVG= INTGRL(WLVG, RWLVG, 1.) 
	WST=  INTGRL(WST, RWST, 1.)
	WRT=  INTGRL(WRT,RWRT, 1.)
	WSO=  INTGRL(WSO, RWSO, 1.)
	TAGBG= WLVG+WST+WSO
      
*---- Dry weights of dead crop organs and total above-ground biomass incl. dead crop organs and harvests (kg DM/ha)
	WSTD= INTGRL(WSTD,DRST, 1.)
	WRTD= INTGRL(WRTD, DRRT, 1.)
      WLVD= INTGRL(WLVD, DLV, 1.)


*----- Rooting depth and Leaf area index
       RD= INTGRL(RD, RR, 1.)
	 LAI= INTGRL(LAI, RLAI, 1.)

*     Soil mineral N  and Total mineral N available from both fertiliser and soil (kg N ha-1)
       NMIN= INTGRL(NMIN, RNMINS, 1.)
	 NMINT= INTGRL(NMINT, RNMINT,1.) 

*----- Total N uptake by crop over time (kg N ha-1) from soil and by biological fixation
       NUPTT= INTGRL(NUPTT, NUPTR, 1.)

*-----Actual N amount in various living organs and total living N amount(kg N ha-1)
      ANLV =  INTGRL (ANLV,RNLV, 1.)
      ANST =  INTGRL (ANST,RNST, 1.)
      ANRT =  INTGRL (ANRT,RNRT, 1.)
      ANSO =  INTGRL (ANSO,RNSO, 1.)
	NLIVT= ANLV+ANST+ANRT+ANSO

*-----N losses from leaves, roots and stems due to senescence and total N loss (kg N ha-1)
      NLOSSL=  INTGRL(NLOSSL, RNLDLV, 1.)
	NLOSSR=  INTGRL(NLOSSR, RNLDRT, 1.)
	NLOSSS=  INTGRL(NLOSSS, RNLDST, 1.)
      NLOSST=  NLOSSL+NLOSSR+NLOSSS

*---- total N in living and dead roots
      NROOT= ANRT + NLOSSR

*---- New state variables for grassland
*----

*----  Hypothetical development stage with TSUM1 (600 C*d) taken from routine TILSUB
       DVS= TSUML/TSUM1

*----  Soil temperature (C)
       SOITMP= INTGRL(SOITMP, RSOITM, 1.)

*---   Days after harvest (d)
       DAHA= INTGRL(DAHA, RDAHA, 1.)

*---   Number of tilers (m-2)
       TILLER= INTGRL(TILLER, DTIL, 1.)

*---- Dry weight of cutted green leaves  (kg ha-1)
       GRASS= INTGRL( GRASS, HARV, 1.)

*----  dry weight of dead leaves and of total above ground living and dead biomass (kg ha-1)
       WLVDD= WLVD - GRASS
	 TAGB= TAGBG + WLVDD + WSTD

*----  Harvestable leaf weight (kg ha-1) in the field
       HRVBL= WLVG -CWLVG

*---- Dry weight of reserves (i.e. storage carbohydrates)
      WRE= INTGRL(WRE, RRE, 1.)

*---- Total above ground living and dead dry weight plus harvest (kg ha-1)
       TADRW= GRASS + TAGB

*---- Harvestable leaf weight in the field plus the harvested amounts of grass (kg ha-1) 
       YIELD= GRASS + MAX(0., HRVBL)

*----  Length of leaves (cm)
       LENGTH= INTGRL(LENGTH, LERA2, 1.)

*----  Calculating specific leaf area index (ha kg-1)
       SLAINT= LAI/MAX(0.1, WLVG)

*-----------------------------------------------------------------------------
*      RATE CALCULATIONS
*-----------------------------------------------------------------------------


*------ Weather calculations

*------ Daily photosynthetically active radiation (PAR, MJ/m2)
        PAR= AVRAD/1.0E6 * 0.50

*------ Average daily temperature (C)
        TMPA= 0.5 *(TMIN + TMAX)
        DAY= REAL(IDAY)

*------ date of planting/sowing
        IF (PL) THEN
          IF (.NOT. YCH) PUSHPL= INSW(DAY - DAYPL, 0., 1.)
	    IF (YCH) PUSHPL= INSW(365. + DAY - DAYPL, 0., 1.)
        ELSE
          PUSHPL= 0.
        ENDIF

        IF (PL .AND. (TSULP .GE. TSUMEM) .AND. (.NOT. EMERG)) THEN
           IDEMERG= IDAY
           DAYEM= REAL(IDEMERG)
           EMERG= .TRUE.
        ENDIF 

*-----  Reduction of development rate until flowering by  day length

        DVRED= LINT(PHOTTB, ILPHOT, DAYLP)

        IF (IDSL .EQ. 1 .AND. .NOT. FLOW) THEN
	      RDAYL= DVRED 
        ELSE
           RDAYL= 1.
        ENDIF

*       emergence date set
        IF (.NOT. PL) IDEMERG= IDEM

*------ Change in temperature sums from sowing/planting (P) and from emergence for crop development without and with day length and vernal. effects
        DTSULP= LIMIT(0., TEFFMX-TBASEM, TMPA-TBASEM) * PUSHPL
        IF (.NOT. YCH) PUSHEM= INSW(DAY - DAYEM, 0., 1.)
        IF (YCH) PUSHEM=INSW(365. + DAY-DAYEM, 0., 1.)
        DTSU=  MAX(0.,LINT(DTSMTB, ILDTSM, TMPA))
	  DTSUM= DTSU * PUSHEM
        DTSUML= DTSU * PUSHEM * RDAYL
		   

*------ Plant growth
        
*------ Radiation use efficiency as dependent on development stage (g DM MJ-1)
        RUE= LINT(RUETB,ILRUE,DVS)
        
*------ Correction of radiation use efficiency for change in atmospheric CO2 concentration (-)
        RCO= LINT(COTB,ILCO,CO)
        
*------ Reduction of radiation use efficiency for non-optimal day-time temperatures, 
*        for low minimum temperature and high radiation level (in MJ m-2 d-1)
        DTEMP= TMAX - 0.25*(TMAX-TMIN)
        RTMP= LINT(TMPFTB,ILTMPF,SOITMP) * LINT(TMNFTB,ILTMNF,TMIN)
	  RRDD= LINT(RDDTB,ILRDD,AVRAD/1.E6)
*------ Correction of RUE for both non-optimal temperatures, atmospheric CO2 and high radiation level
	  RTMCO= RTMP * RCO * RRDD
        
!      Calling the subroutine for translocatable N in leaves, stem, roots and
!      storage organs (kg N ha-1)
       CALL NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT,RNFLV,RNFST,RNFRT,FNTRT,
     $        ATNLV,ATNST,ATNRT,ATN)
     
!    * Total vegetative living above-ground biomass (kg DM ha-1)
        TBGMR =WLVG+WST
    
!      N concentration (kg N kg-1 DM) in the living leaves, stem (if applic.), roots and storage
!      organs (if aplic.)
       NFLV = ANLV/ NOTNUL(WLVG)
       NFST = ANST/ NOTNUL(WST)
       NFRT = ANRT/ NOTNUL(WRT)
       NFSO = ANSO/ NOTNUL(WSO)
    
!    Total N in vegetative living above-ground biomass (kg N ha-1)
       NUPGMR = ANLV + ANST
     
!    Fertilizer N application (kg N ha-1 d-1) and its recovery fraction (-)---------------------------------*
       FERTN  = LINT (FERTAB,ILFERT, DAY)
       NRF    = LINT (NRFTAB,ILNRFT, DAY)
       FERTNS = FERTN * NRF

!      Check on N balance
       NBALAN = ABS(NUPTT+(ANLVI+ANSTI+ANRTI+ANSOI)-(ANLV
     $          +ANST+ANRT+ANSO+NLOSSL+NLOSSR+NLOSSS))  

      IF (NBALAN .GE. 1.) STOP
     $    ' nitrogen balance NBALAN not 0, program aborted'
   
!    Total leaf weight, both green and dead (kg DM ha-1)
       WLV    = WLVG + WLVDD
  
!     Relative death rate of roots and stems (if applic.)(d-1)
       RDRRT = LINT(RDRRTB, ILRDRR, DVS)
	 RDRST = LINT(RDRSTB, ILRDRS, DVS)
   
!     Total N in living above-ground crop organs (kg N ha-1)   
       NTAG   = ANLV   + ANST   + ANSO
   
!     Relative death rate of leaves due to senescence/ageing as dependent on mean daily temperature (d-1)
       RDRTMP = LINT(RDRLTB,ILRDRL,TMPA)
  
!     Maximum N concentration in the leaves, from which the N conc. in the
!     stem and roots are derived, as a function of development stage (kg N kg-1 DM)  
       NMAXLV = LINT (NMXLV, ILNMXL, DVS)
     
!    * N concentration in above-ground living biomass incl. possibly seeds (kg N kg-1 DM)
	 NTAC  = NTAG/TAGBG
    
!    N supply to the storage organs (kg N ha-1 d-1)
       NSUPSO = INSW (DVS-DVSNT,0.,ATN/TCNT)

!      N concentration in total vegetative living above-ground biomass  (kg N kg-1 DM) 
       NFGMR  = NUPGMR/NOTNUL(TBGMR)
     
!    * Residual N concentration in total vegetative living above-ground biomass  (kg N kg-1 DM) 
       NRMR   = (WLVG*RNFLV+WST*RNFST)/NOTNUL(TBGMR)
  
!     Nitrogen uptake limiting factor (-) at low moisture conditions in the
!     rooted soil layer before DVSNLT. If DVS becomes larger than DVSNLT, there is no
!     N uptake from the soil
       NLIMIT = INSW(DVS-DVSNLT, INSW(TRANRF-0.01,0.,1.) , 0.0)
     

!     Biomass partitioning functions under non-stressed situations (-)
       FRTWET = LINT(FRTB,ILFR, DVS )
       FLVT   = LINT(FLTB,ILFL, DVS )
       FSTT   = LINT(FSTB,ILFS, DVS)
       FSOT   = LINT(FOTB,ILFO, DVS )
 
*----  Carbon balance check  
       CBALAN = ABS(GTSUM + (WRTI+WLVGI+WSTI+WSOI)-(WLVG+WST+WSO+
     $          WRT+WLVD+WRTD+WSTD))	

       IF (CBALAN .GE. 1.) STOP
     $    ' carbon balance CBALAN not 0, program aborted'
 
    
*----  Maximum N concentrations in stems (if applic.) and roots (kg N kg-1 DM)
       NMAXST = LSNR * NMAXLV
       NMAXRT = LRNR * NMAXLV

*----- Root growth (cm d-1)
       IF (EMERG) RR = MIN(RRI * INSW( TRANRF-0.01, 0., 1. ),  RDM-RD)
        
!    Calling the subroutine for calculating optimal nitrogen concentrations in leaves
!    and stems
       CALL NOPTM(FRNX,NMAXLV,NMAXST, NOPTLV,NOPTST)
    
!    Optimal amount of N in vegetative above-ground living biomass and its N concentration
       NOPTS = NOPTST* WST
       NOPTL = NOPTLV* WLVG 
       NOPTMR = (NOPTL+ NOPTS)/NOTNUL(TBGMR)
    
!    Calling the subroutine for calculating the Nitrogen Nutrition Index (NNI)
       CALL NNINDX(DAY,DAYEM,EMERG,NFGMR,NRMR,NOPTMR,NNI)
!    With potential and water limited conditions there is no N stress and NNI is set to 1   
       IF (IOPT .EQ. 1 .or. IOPT .EQ. 2) NNI= 1.


!     Biomass partitioning functions under water-stress situation (-)
       FRT = FRTWET * LINT(FRTRB,ILFRT, TRANRF)
       FLV   = FLVT 
       FST   = FSTT
       FSO   = FSOT


       FCHECK  = ABS(FRT + (FLV+ FST+ FSO) * (1.- FRT) -1.)

       IF (FCHECK .GE. 0.05) STOP
     $ ' assimilate allocation check over crop organs FCHECK '
     $'not 0, program aborted'

*      Call to subroutine for grassland management options
       CALL MOWING (IMOPT,INCUT,MNDAT,INH,DAY,WLVG, 
     $  CWGHT,CWLVG,DAHA,RDAHA,HARV)
*
*        Temperature dependent leaf appearance rate, according to
*        (Davies and Thomas, 1983), soil temperature (SOITMP)is used as
*        driving force which is estimated from a 10 day running
*        average
         REDTMP= LINT(TMPFTB,ILTMPF,SOITMP)
         LEAFNE=   FCNSW(REDTMP, 0., 0.,SOITMP * 0.01 )

*        Leaf elongation rate affected by temperature
*        cm day-1 tiller-1
         LERA= FCNSW(TMPA-TMBAS1, 0., 0.,
     $    0.83*LOG(MAX(TMPA, 3.))-0.8924 )
*
         LERA2 = INSW (HARV-0.1, LERA, -LENGTH)

*      Maximum site filling new buds (FSMAX) decreases due
*      to low nitrogen contents, Van Loo and Schapendonk (1992)
*      Theoretical maximum tillering size = 0.693
        FSMAX = NNI*0.693

        CALL TILSUB (TILLER,FSMAX,LAI,LAICR,DAHA,LEAFNE,TSUML,REDTMP, 
     $                DTIL)

*        Rate of sink limited leaf growth, unit of TILLER is tillers m-2,
*        1.0E-8 is conversion from cm-2 to ha-1, ha leaf ha ground-1 d-1
         GLAISI = (TILLER * 1.0E4 * (LERA * 0.3)) * 1.0E-8

    
!    Calling the subroutine for total growth rate (kg DM ha-1 d-1)
       KDIF= LINT(KDIFTB,ILKDIF,DVS)
       CALL GROWTH(DAY,EMERG,PAR,KDIF,NLUE,LAI,RUE,RTMCO,TRANRF,FINT,
     $	 NNI,HARV,PARINT,GTWSO1)
  
       GTWSO2 = GTWSO1+WRE/DELT
         DRE   = WRE/DELT

!        Specific Leaf area(ha/kg)
       SLA = LINT (SLATB,ILSLA,DVS)*EXP(-NSLA * (1.-NNI))

*        Conversion to total sink limited carbon demand,
*        kg DM ha-1 ground area d-1
         GTWSI= FCNSW(HARV,GLAISI * (1./SLA) * (1./FLV) * (1./(1.-FRT)),
     $      GLAISI * (1./SLA) * (1./FLV) * (1./(1.-FRT)), 0.)
*
*        Actual growth switches between sink- and source limitation
*        (more or less dry matter formed than can be stored)
         GRE= FCNSW(GTWSO2-GTWSI,0.,0., GTWSO2-GTWSI)
         GTW= FCNSW(GTWSO2-GTWSI,GTWSO2,GTWSO2, GTWSI)

*        Change in reserves
         RRE = GRE-DRE


!     Water-Nitrogen stress factor
       RNW = MIN(TRANRF,NNI)
     
!     cumulative values for TRANRF and NNI over growth period
       IF (EMERG) CTRAN= CTRAN + TRANRF
	 IF (EMERG) CNNI= CNNI + NNI
      
    
!    Calling the subroutine for relative death rate of leaves.
       CALL DEATHL(DAY,EMERG,DVS,DVSDLT,RDRTMP,RDRSHM,RDRL,TRANRF,LAI,
     $	 LAICR,WLVG,RDRNS,NNI,SLA,SLAINT,
     $     RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,DLAIS,DLAINS,DLAI,HARV)
     
!    ** Leaf growth
       GLV    = FLV * GTW * (1-FRT)
    
!    Calling the subroutine for calculating the daily increase of leaf area index (m2 m-2 d-1).
       DTEFF= MAX(0., TMPA-TBASE)
       CALL GLA(DAY,EMERG,DTEFF,LAII,RGRLAI,DELT,SLA,LAI,GLV,NLAI,
     $      DVS,TRANRF,NNI,GLAI)
    
!    Calling the subroutine for N loss due to death of leaves,stems and roots (kg N ha-1 d-1)
       CALL RNLD (DVS,WRT,WST,RDRRT,RDRST,RNFLV,DLV,RNFRT,
     $ RNFST,DRRT,DRST,RNLDLV,RNLDRT,RNLDST,NFGMR,HARV)
  
!    Net rate of change of Leaf area (m2 leaf area m-2 d-1)
       RLAI   = GLAI - DLAI
    
!    Calling the subroutine for calculating the relative growth rate of roots, leaves, stem (if applic.)
!    and storage organs (if applic.) (kg ha-1 d-1)
       CALL RELGR(DAY,DAYEM,EMERG,GTW,FLV,FRT,FST,
     $ FSO,DLV,DRRT,DRST,RWLVG,RWRT,RWST,RWSO)
    
!    Calling the subroutine for N demand of leaves, roots, stems and storage
!    organs (if applic.) (kg N ha-1 d-1)
       CALL NDEMND(NMAXLV,NMAXST,NMAXRT,NMAXSO,WLVG,WST,WRT,WSO,
     $   ANLV,ANST,ANRT,ANSO,TCNT,NDEML,NDEMS,
     $   NDEMR,NDEMSO)
     
!        Total Nitrogen demand (kg N ha-1)
       NDEMTO = MAX (0.0,(NDEML + NDEMS + NDEMR))
    
!    Rate of N uptake in grains (if applic.) (kg N ha-1 d-1)
       RNSO =  AMIN1 (NDEMSO,NSUPSO)
      
	IF (EMERG) THEN
!        Total N uptake (kg N ha-1 d-1) from soil 
         NUPTR = (MAX (0., MIN (NDEMTO, NMINT))* NLIMIT)/DELT
	   
!        No N limitation for optimal and water limited production
	   IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) NUPTR=
     $    (MAX (0., NDEMTO)* NLIMIT ) / DELT 
      ELSE
	  NUPTR= 0.
	END IF    

!    Calling the subroutine for calculating N translocated from leaves, stem, and roots (kg N ha-1 d-1)
       CALL NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN, RNTLV,RNTST,RNTRT)
    
!    Calling the subroutine to compute the partitioning of the total
!    N uptake rate (NUPTR) over the leaves, stem and roots (kg N ha-1 d-1)
       CALL RNUSUB(DAY,DAYEM,EMERG,NDEML,NDEMS,NDEMR,NUPTR,
     $ NDEMTO,RNULV,RNUST,RNURT)
 
!     Soil N supply (g N m-2 d-1) through mineralization during crop growth
       IF (EMERG) RNMINS  = -MAX(0.,MIN( RTMINS * NMINI * NLIMIT, NMIN))
!     Change in total inorganic N in soil as function of fertilizer
!     input, soil N mineralization and crop uptake.
       RNMINT = FERTNS/DELT -NUPTR - RNMINS
     
*----Rate of change of N in crop organs   
       RNST = RNUST-RNTST-RNLDST
       RNRT = RNURT-RNTRT-RNLDRT
       RNLV = RNULV-RNTLV-RNLDLV       

*---- Soil temperature change
      RSOITM= (TMPA -SOITMP) / 10.

999      CONTINUE

* ----- Data for summary output

        IF (TERMIN) THEN
        WRITE (*, '(//2A)') ' Crop development completed'

        WLVG1=WLVG
	  WLVDD1=WLVDD
	  WRE1=WRE
	  GRASS1=GRASS
	  YIELD1= YIELD
	  TADRW1=TADRW
        TPARINT1= TPARINT
        IDHALT= IDAY
        NUPTT1= NUPTT
	  NLIV1= NLIVT
	  NLOSS1= NLOSST
	  NROOT1= NROOT
	    IF (IDAY .GT. IDEMERG) THEN
	    TRANRF1= CTRAN/(IDAY -IDEMERG)
	    NNI1= CNNI/(IDAY-IDEMERG)
	    ELSE
	    TRANRF1= CTRAN/ (IDAY+365-IDEMERG)
          NNI1= CNNI/(IDAY+365-IDEMERG)
	    END IF

*----- calculated radiation use efficiency in g above-ground D.M./radiation intercepted in MJ PAR
        RUEC= (TADRW/10.)/TPARINT
        END IF


*---- Finish conditions (not used for grass; LAI kept at minimal value of 0.08 to allow regrowth)
*     IF (FINT .LT. 0.05 .AND.  TAGB .GT. 200.) TERMIN= .TRUE.


         RETURN
         END
