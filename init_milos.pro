;+
; NAME:
;   INIT_MILOS
;
; AUTHOR:
;   D. Orozco Suarez
;      National Astronomical Observatory of Japan,
;      2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;      d.orozco@nao.ac.jp
;
;   J.C. del Toro Iniesta
;      Instituto de Astrofisica de Andalucia (CSIC)
;      Apdo de Correos 3004, 18080 Granada, SPAIN
;      jti@iaa.es
;
; PURPOSE:
;   SETUP OF MILOS CODE. MANDATORY TO EXECUTE TO USE MILOS.
;   IT CREATES THE COMMON BLOCK QUANTEN AND ALSO COMPILE ALL NECESSARY ROUTINES BY MILOS.
;
; CATEGORY:
;   Spectropolarimetric fitting
;
; CALLING SEQUENCE:
;   INIT_MILOS, LINES, WLI
;
; INPUTS:
;        LINES: Line of interest. String.
;
; OUTPUTS:
;   WLI: A scalar or array with the central wavelength(s) of the line(s)
;      in Angstroms
;      Dimensions: (number of lines, central wavelengths)
;      Examples: [1, 6302.5]; [2, 6301.5, 6302.5]
;      WLI IS THE FIRST INPUT PARAMETER NEEDED BY MILOS.PRO
;
; KEYWORD PARAMETERS:
;
; COMMON BLOCKS:
;
; CALLED ROUTINES
;      CREATE_NC, TYPE
;
; EXAMPLE:
;      IDL> INIT_MILOS,'6301',WLI
;
; USAGE:
;     Run IDL> INIT_MODEL with no options to see Usage and AVAILABLE SPECTRAL LINES.
;
;
; MODIFICATION HISTORY:
;   First beta version created, D Orozco Suarez (DOS), 2008
;   Added some error handling, 24 Feb, 2008
;   Changed eta ceto from 0.342466 to 0.3266, 21 Jan, 2010. DOS
;-

pro init_milos,LINES,WLI,quiet=quiet,not_normalize=not_normalize

;Check input parameters
IF N_PARAMS() NE 2 THEN BEGIN
	print,' '
	print,' ------------  ATTENTION: CALLING ERROR -------------------------'
	print,' It seems that you are trying to execute INIT_MILOS with less or '
	print,'   more that 2 parameters. Usage example:'
    print,'      IDL> INIT_MILOS,"6301",WL  '
	print,' AVAILABLE LINES: '
    print,'   Label = "6301"     , Line  ->  # Fe 6301.5 A               #    '
    print,'   Label = "6302"     , Line  ->  # Fe 6302.5 A               #'
    print,'   Label = "63016302" , Line  ->  # Fe 6301.5 A + Fe 6302.5 A #'
    print,'   Label = "6767"     , Line  ->  # Fe 6767.79 A              #'
    print,'   Label = "5250.2"   , Line  ->  # Fe 5250.208 A             #'
    print,'   Label = "5250.6"   , Line  ->  # Fe 5250.208 y 645 A       #'
    print,'   Label = "5250"     , Line  ->  # Fe 5250.645 A             #'
    print,'   Label = "6173"     , Line  ->  # Fe 6173.3356 A            #'
    print,'   Label = "10830"    , Line  ->  # 10830 HELIUM LINE TRIPLET  '
    print,'                                    (SL COUPLING)             #'
	print,' ----------------------------------------------------------------'
	print,' '
	return
ENDIF
a=type(lines)
IF a ne 7 then begin
	print,' '
	print,' ------------  ATTENTION: CALLING ERROR -------------------------'
	print,' INPUT parameter LINES is not of required type '
	print,'   LINES:  string.'
	print,' ----------------------------------------------------------------'
	print,' '
	return
ENDIF

;
; Compilation of several needed routines
;
RESOLVE_ROUTINE,'me_der',/either
RESOLVE_ROUTINE,'mil_sinrf',/either
RESOLVE_ROUTINE,'fvoigt',/either
RESOLVE_ROUTINE,'filtro',/either
RESOLVE_ROUTINE,'create_nc',/either
RESOLVE_ROUTINE,'type',/either
RESOLVE_ROUTINE,'quanten',/either
RESOLVE_ROUTINE,'sign',/either
RESOLVE_ROUTINE,'weights_init',/either
RESOLVE_ROUTINE,'check_param',/either
RESOLVE_ROUTINE,'lm_mils',/either
RESOLVE_ROUTINE,'covarm',/either

;
;
;


;*** SPECTRAL LINES ***

;### Fe 6173. ###;

CASE LINES OF
'6301': BEGIN
print,'### Fe 6301.5 ###'
           WL=6301.515d0
           slo=2 & llo=1 & jlo=2
           sup=2 & lup=2 & jup=2
           DATA=[1.,slo,llo,jlo,sup,lup,jup]
        END
'6302': BEGIN
print,'### Fe 6302.5 ###'
           WL=6302.494d0
           slo=2 & llo=1 & jlo=1
           sup=2 & lup=2 & jup=0
           DATA=[1.,slo,llo,jlo,sup,lup,jup]
        END
'63016302': BEGIN
print,'### Fe 6301.5 + Fe 6302.5 ###'
           wl=[6301.5012d0 , 6302.4936d0]
           slo1=2 & llo1=1 & jlo1=2
           sup1=2 & lup1=2 & jup1=2
           slo2=2 & llo2=1 & jlo2=1
           sup2=2 & lup2=2 & jup2=0
           DATA=[2.,slo1,llo1,jlo1,sup1,lup1,jup1,1.,$
                  slo2,llo2,jlo2,sup2,lup2,jup2,0.3266];0.342466]
            END
'5250': BEGIN
print,'### Fe 5250.2 + Fe 5250.6 ###'
           wl=[5250.208d0 , 5250.645d0]
           slo1=2 & llo1=2 & jlo1=0
           sup1=3 & lup1=2 & jup1=1
           slo2=2 & llo2=1 & jlo2=2
           sup2=2 & lup2=1 & jup2=3
           DATA=[2.,slo1,llo1,jlo1,sup1,lup1,jup1,1.,$
                  slo2,llo2,jlo2,sup2,lup2,jup2,0.3]
            END
'5250.2': BEGIN
print,'### Fe 5250.208 ###'
           WL=5250.208d0
           slo=2 & llo=2 & jlo=0
           sup=3 & lup=2 & jup=1
           DATA=[1.,slo,llo,jlo,sup,lup,jup]
            END
'5250.6': BEGIN
print,'### Fe 5250.6 ###'
           WL=5250.645d0
           slo=2 & llo=1 & jlo=2
           sup=2 & lup=1 & jup=3
           DATA=[1.,slo,llo,jlo,sup,lup,jup]
            END
'6173': BEGIN
print,'### Fe 6173.3 ###'
           WL = 6173.3500 ; C_MILOS 6173.3356d0
           slo=2 & llo=1 & jlo=1
           sup=2 & lup=2 & jup=0
           DATA=[1.,slo,llo,jlo,sup,lup,jup]
        END
'6767': BEGIN
print,'### Fe 6767.79 ###'
           WL=6767.79d0
           slo=0 & llo=0 & jlo=0
           sup=1 & lup=1 & jup=1
           DATA=[1.,slo,llo,jlo,sup,lup,jup]
        END
'10830': BEGIN
print,'### HELIUM LINE TRIPLET (SL COUPLING) ###'
           WL = [ 10830.2501d0, 10830.3397d0, 10829.0911d0]
           slo0=1 & llo0=0 & jlo0=1
           sup0=1 & lup0=1 & jup0=1
           slo1=1 & llo1=0 & jlo1=1
           sup1=1 & lup1=1 & jup1=2
           slo2=1 & llo2=0 & jlo2=1
           sup2=1 & lup2=1 & jup2=0
           DATA=[3.,slo0,llo0,jlo0,sup0,lup0,jup0,0.3333,$
                slo1,llo1,jlo1,sup1,lup1,jup1,0.556,slo2,llo2,jlo2,sup2,lup2,jup2,0.1111]
        END
				'10830+si': BEGIN
				print,'### HELIUM LINE TRIPLET (SL COUPLING) + Silicon ###'
				           WL = [10827.089d0,10830.2501d0, 10830.3397d0, 10829.0911d0]
									 slo=1 & llo=1 & jlo=2
									 sup=1 & lup=1 & jup=2
				           slo0=1 & llo0=0 & jlo0=1
				           sup0=1 & lup0=1 & jup0=1
				           slo1=1 & llo1=0 & jlo1=1
				           sup1=1 & lup1=1 & jup1=2
				           slo2=1 & llo2=0 & jlo2=1
				           sup2=1 & lup2=1 & jup2=0
				           DATA=[4.,slo,llo,jlo,sup,lup,jup,1,slo0,llo0,jlo0,sup0,lup0,jup0,0.3333,$
				                slo1,llo1,jlo1,sup1,lup1,jup1,0.556,slo2,llo2,jlo2,sup2,lup2,jup2,0.1111]
				        END
								'si': BEGIN
								print,'### Silicon ###'
								           WL = [10827.089d0]
													 slo=1 & llo=1 & jlo=2
													 sup=1 & lup=1 & jup=2
								           DATA=[1.,slo,llo,jlo,sup,lup,jup,1]
								        END
        ;3=FE 1      15648.515       1.0         5.426    -0.669  7D 1.0- 7D 1.0  0.229  2.7289e-14
		'15648': BEGIN
print,'### Fe 15648 ###'
           WL=15648.515d0
           slo=3 & llo=2 & jlo=1
           sup=3 & lup=2 & jup=1
           DATA=[1.,slo,llo,jlo,sup,lup,jup]
        END

ELSE: BEGIN
	print,' '
	print,' ------------  ATTENTION: LINE NOT FOUND  -----------------------'
	print,' You set LINES = ', LINES
	print,' AVAILABLE LINES: '
    print,'   Label = "6301"     , Line  ->  # Fe 6301.5 A               #    '
    print,'   Label = "6302"     , Line  ->  # Fe 6302.5 A               #'
    print,'   Label = "63016302" , Line  ->  # Fe 6301.5 A + Fe 6302.5 A #'
    print,'   Label = "6767"     , Line  ->  # Fe 6767.79 A              #'
    print,'   Label = "5250.2"   , Line  ->  # Fe 5250.208 A             #'
    print,'   Label = "5250.6"   , Line  ->  # Fe 5250.645 A             #'
    print,'   Label = "6173"     , Line  ->  # Fe 6173.3356 A            #'
    print,'   Label = "10830"    , Line  ->  # 10830 HELIUM LINE TRIPLET  '
    print,'                                    (SL COUPLING)             #'
    print,' ----------------------------------------------------------------'
	print,' '
	return
        END
ENDCASE

CREATE_NC,DATA,not_normalize=not_normalize

WLI = FLTARR ( N_ELEMENTS ( WL ) + 1)
WLI[0] = N_ELEMENTS ( WL )
IF N_ELEMENTS(WL) GT 1 THEN WLI[1:*] = WL ELSE WLI[1]=WL

return
end
