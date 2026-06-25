*---------------------------------------------------------------------*
*                                                                     *
*     PROGRAM LINGRAN.FOR                                             *
*     Author : joost wolf                                             *
*     Date of last revision:  March 2012                              *
*                                                                     *
*     Purpose: This model is the LINGRA fst model but in FORTRAN      *
*              It simulates the growth of a perennial rye grass       *
*              (L. Perenne) as function of intercepted radiation,     *
*               temperature and light use efficiency (i.e.            *
*              (potential production) and also the water and N limited*
*               productions. In these situations soil water (free     *
*               drainage) and simple nitrogen balances are simulated  *
*               and also the effects of water and nitrogen supply on  *
*               crop growth.                                          *
*---------------------------------------------------------------------*

      PROGRAM LINGRAN

      IMPLICIT REAL (A-Z)
      INTEGER ISYR,ISYEAR,IUWE,IIYR,INYEAR,IDPL,IDEM,IFINIT,IGAP,ISOIL
      INTEGER ICROP,IOPT,IOUT,IDAY,IYEAR,IWAR,IRRI,IDEMERG,IDFLOW
      INTEGER IDHALT,IAIRDU

      CHARACTER RUNNAME*5,REMARK*80,STATR*5,CONTIN*1,WTRDIR*80,CONTW*1
      COMMON /ITIME/ ISYR, ISYEAR

      LOGICAL TERMIN,INITI,TERMY,PL,EMERG,YCH,CWAVT      
      
*-----iunit number for weather
      DATA IUWE/65/

*---- initialize program
*     LOOP for different weather stations

      
1      CALL INITP (IIYR,INYEAR,IDPL,IDEM,IFINIT,IGAP,REMARK,RUNNAME,
     $  STATR,CONTIN,CONTW,SENSP,SENSV,SENSW,SENSR,SENST,CO,
     $  ISOIL,ICROP,IOPT,IRRI,IOUT,WTRDIR)


       TERMY= .FALSE.

*---- LOOP for a number of years

      IYEAR = IIYR

5     CONTINUE
            
      TERMIN = .FALSE.

	IF (IYEAR .NE. IIYR .AND. (CONTW .EQ.'Y'.OR. CONTW .EQ.'y')) THEN    
	 CWAVT= .TRUE.
	ELSE
	 CWAVT= .FALSE.
	ENDIF	 

      IF (IDPL .LE. 0) THEN
        IDAY= IDEM
        PL= .FALSE.
      ELSE 
        IDAY= IDPL
        PL= .TRUE.
      ENDIF
      
      ISYR = IYEAR 
      ISYEAR= IYEAR + 1900
      INITI = .TRUE.

*---- LOOP for one growth period (year)
        
	 YCH= .FALSE.
        
10       CONTINUE

*------- read weather data
         CALL WEATHR (IWAR,WTRDIR,STATR,IYEAR,IDAY,IGAP,LONGIE,LATIN,
     $                ALTI,TMIN,TMAX,DTR,RAIN,VAP,WIND)


* ----- to prevent temperature effect on relative humidity
          
          TMPA = (TMIN + TMAX)/2.
          SVAP1  = 6.10588 * EXP (17.32491*TMPA/(TMPA+238.102))
          VAP= AMIN1(VAP,SVAP1)
          RH= VAP/SVAP1

* ------  sensitivity analyses
        
           TMIN= TMIN + SENST
           TMAX= TMAX + SENST  
         
            TMPA = (TMIN + TMAX)/2.
           SVAP2  = 6.10588 * EXP (17.32491*TMPA/(TMPA+238.102))
           VAP= RH*SVAP2
         
               
           VAP= AMIN1(SVAP2, VAP * SENSV)
           WIND= WIND * SENSW
           DTR= DTR * SENSR
           RAIN= RAIN * SENSP
          
*------- calculate daylength

         CALL ASTRO (IDAY,LATIN,DAYL,DAYLP,SINLD,COSLD)                    


*------- calculate potential soil evaporation and crop transpiration

         CALL PENMAN (IDAY,DAYL,SINLD,COSLD,ALTI,TMIN,TMAX,
     $                DTR,WIND,VAP,CO,E0,ES0,ETC,AVRAD)


*------- calculate soil water balance  

         CALL WATBALS(ICROP,ISOIL,INITI,IOPT, IRRI,TERMIN,EMERG,IDAY,
     $              ES0,ETC,RAIN, FINT,DEPNR,RD,RDMCR,RR,RDM,
     $              CFET,IAIRDU,SMACT,TTRANS,TDRAIN,TRAIN,TESOIL,TRUNOF,
     $             TIRR,TRANRF,RUNFR,WTOT, WTOTL,WAVT,WAVTL,CWAVT,CONTW)


*------- when finish conditions are reached (TERMIN .TRUE.)

*-------- calculate crop growth
         
         IF (TERMIN) THEN
         CALL CROPG(ICROP,ISOIL,INITI,IOPT, IDAY,IDEM,
     $    IDEMERG,IDPL,IDFLOW,IDHALT,PL,TERMIN,EMERG,TMIN,TMAX,AVRAD,CO,
     $	 TRANRF,RDMSO,DAYLP,WLVG, WLVDD, WRT, WRE,GRASS,YIELD,TADRW,
     $      TILLER,RD,RDMCR,RR,RDM,LAI,
     $        CFET,DEPNR,IAIRDU,TSUM,DVS,DVSEND,TSULP,FINT,TPARINT,TPAR,
     $        TSUML,NNI,NMINT,NMIN,NUPTT,NLIVT,NLOSST,YCH)


*------- daily output
           CALL DAILOUT (IOPT,IDAY,IOUT,IDEM,IDPL,IDEMERG,
     $                   WLVG, WLVDD, WRT,WRE,GRASS,YIELD,TADRW,
     $                   TILLER,TPARINT,TPAR,REMARK,DVS,RUNNAME,LAI,
     $                   SMACT,TTRANS,TDRAIN,TRAIN,TESOIL,TRUNOF,TIRR,
     $                   TRANRF,TERMIN,STATR,ISOIL,ICROP, IRRI,RUNFR,CO,
     $                   WTOT,WAVT,TSUML,NNI,NMINT,NMIN,NUPTT,
     $                   NLIVT,NLOSST)

*------- yearly output
           CALL STOUT (RUNNAME,REMARK,STATR,ISOIL,ICROP,TERMY,
     $         IDEMERG,IDPL,IDFLOW,IDHALT,CO,IOPT)

         GOTO 20

         ENDIF


*-------- calculate crop growth

         CALL CROPG(ICROP,ISOIL,INITI,IOPT, IDAY,IDEM,
     $    IDEMERG,IDPL,IDFLOW,IDHALT,PL,TERMIN,EMERG,TMIN,TMAX,AVRAD,CO,
     $	 TRANRF,RDMSO,DAYLP,WLVG, WLVDD, WRT, WRE,GRASS,YIELD,TADRW,
     $      TILLER,RD,RDMCR,RR,RDM,LAI,
     $        CFET,DEPNR,IAIRDU,TSUM,DVS,DVSEND,TSULP,FINT,TPARINT,TPAR,
     $        TSUML,NNI,NMINT,NMIN,NUPTT,NLIVT,NLOSST,YCH)


   
*------- daily output 

           CALL DAILOUT (IOPT,IDAY,IOUT,IDEM,IDPL,IDEMERG,
     $                   WLVG, WLVDD, WRT,WRE,GRASS,YIELD,TADRW,
     $                   TILLER,TPARINT,TPAR,REMARK,DVS,RUNNAME,LAI,
     $                   SMACT,TTRANS,TDRAIN,TRAIN,TESOIL,TRUNOF,TIRR,
     $                   TRANRF,TERMIN,STATR,ISOIL,ICROP, IRRI,RUNFR,CO,
     $                   WTOT,WAVT,TSUML,NNI,NMINT,NMIN,NUPTT,
     $                   NLIVT,NLOSST)
 
         IF (TERMIN) GOTO 15   


*------- update daynumber and year

         CALL TIMER (IDAY,IYEAR,IFINIT,TERMIN,INITI,DVS,DVSEND,YCH)

         INITI = .FALSE.                
 
 
15      GOTO 10

20      CONTINUE
      
*-----  For permanent grassland the year number is changed in the TIMER routine 
*       at day 365, so no updating of year number is needed here
*-----  For temperate grassland the new grass crop may start in a new year	     
        IF (YCH) THEN
          IYEAR= IYEAR 
	  ELSE
	    IYEAR= IYEAR + 1
	  ENDIF         

        IF (IYEAR .NE. IIYR + INYEAR)  GO TO 5

*------- when finish conditions are reached for total number of years (TERMY= .TRUE.)

         TERMY= .TRUE.

*------- statistical analysis
         CALL STOUT (RUNNAME,REMARK,STATR,ISOIL,ICROP,TERMY,
     $         IDEMERG,IDPL,IDFLOW,IDHALT,CO,IOPT)


*------- runs for new site or years ?

         IF (CONTIN .EQ. 'Y' .OR. CONTIN .EQ. 'y') GOTO 1
      END
