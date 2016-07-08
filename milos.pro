;+
; NAME:
;   MILOS
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
; VERSIONS:
;   For current and updated versions go to "hinode.nao.ac.jp/~orozco/MILOS" or
;
; PURPOSE:
;   Synthesis and inversion of Stokes profiles under the Milne-Eddington approximation
;
; CATEGORY:
;   Spectropolarimetric fitting
;
; CALLING SEQUENCE:
;   MILOS, WLI, AXIS, MODEL, STOKESPROF[, YFIT=yfit][, ERR=err][,$
;      CHISQR=chisqr][,ITER=iter][,SLIGHT=slight][,TOPLIM=toplim][,$
;      TRIPLET=triplet][,QUIET=quiet][,MITER=miter][,$
;      WEIGHT=weight][,FIX=fix][,SIGMA=sigma][,$
;      SYNTHESIS=synthesis][,INVERSION=inversion][,RFS=rfs][,ILAMBDA=ilambda][,$
;      FILTER=filter][,NOISE=noise][,POL=pol][,GETSHI=getshi][,DOPLOT=doplot][,AC_RATIO=ac_ratio]
;
; INPUTS:
;   WLI: A scalar or array with the central wavelength(s) of the line(s)
;      in Angstroms
;      Dimensions: (number of lines, central wavelengths)
;      Examples: [1, 6302.5]; [2, 6301.5, 6302.5]
;   AXIS: Wavelength axis in Angstroms. Just one for all the lines in the sample
;      It should be a double precision array
;   MODEL: Array with the 11 model-atmosphere parameters
;      eta0 = line-to-continuum absorption coefficient ratio
;      B = magnetic field strength       [Gauss]
;      vlos = line-of-sight velocity     [km/s]
;      dopp = Doppler width              [Angstroms]
;      aa = damping parameter
;      gm = magnetic field inclination   [deg]
;      az = magnetic field azimuth       [deg]
;      S0 = source function constant
;      S1 = source function gradient
;      mac = macroturbulent velocity     [km/s]
;      alpha = filling factor of the magnetic component [0->1]
;   STOKESPROF: Array with the input (if the /INVERSION keyword is set) Stokes profiles
;      Dimensions: (N_ELEMENTS(AXIS), 4)
;   SLIGHT: Array with the stray-light profile
;      Dimensions: (N_ELEMENTS(AXIS), 4)
;   TOPLIM: Optional minimum relative difference between two succesive merit-function values
;   MITER: Maximum number of iterations allowed
;   WEIGHT: Array with arbitrary weights for the Stokes profiles
;      Dimensions: (4) (for individual Stokes profiles)
;      Dimensions: (N_ELEMENTS(AXIS), 4) (as many as the profiles)
;         It can be used, for example, to neglect telluric lines in the fits
;         Also, this keyword is very useful when dealing with filter data that
;            has only 4-5 wavelegth samples. In order to convolve the spectra with
;            an instrumental filter profile, the synthesized Stokes vector has to be of, at least,
;            12 wavelength samples. Less number of samples will increase the errors in the
;            convolution process significantly.
;   FIX: Array with as many (integers) zeroes or ones as parameters to be inverted
;      Dimensions: (11)
;      Example: FIX = [1,1,1,1,1,1,1,1,1,0,0] means "invert everything but mac and alpha"
;   SIGMA: Array with noise values for the profiles
;      Dimensions: (4)
;   ILAMBDA: Initial value for the Levenberg-Marquardt's fudge parameter
;   FILTER: Width (in mA) of a Gaussian filter profile or array with the instrument
;      profile samples. Dimensions: N_ELEMENTS(LMB)
;   NOISE: Scalar with the noise values for the Stokes profiles (Ignored when SIGMA is set)
;   POL: Polarization threshold. If all Q, U, and V are below POL, then the
;      vector magnetic field is set to zero
;
; KEYWORD PARAMETERS:
;   SYNTHESIS: Set this keyword to carry out just a synthesis
;   INVERSION: Set this keyword to carry out an inversion
;   TRIPLET: Set this keyword to do the calculations for the effective triplet
;   QUIET: Set this keyword to cancel any printed message on the screen
;   DOPLOT: Set this keyword to plot the fit
;   MU: Scalar containing the cosine of the heliocentric angle
;   PARLIMITS: Limitations on the model parameter changes can be imposed throught this keywork
;            N-Array of 4 elements [PAR_INDEX,SET,LOWERLIMIT,UPPERLIMIT,COMPONENT].
;     	PAR_INDEX: scalar value indicating the model parameter (index) to change
;       SET: 1-> ON, 0 -> OFF
;       LOWERLIMIT,UPPERLIMIT: Min and Max limits for the selected model parameter
;            Example:
;                 PARLIMITS = [2, 1, -5d0, 5d0]
;					change limits for velocity (index 2, set ON) to vary between -5 and 5 km/s
;                 PARLIMITS = [0, 0, 0, 0]
;					change limits for eta0 (index 0), free variation set -> OFF
;                 PARLIMITS = [[9, 1, 1d0, 2d0],[1, 1, 0d0, 600d0]]
;					change limits for macrovelocity (index 9), to vary between 1 and 2 km/s
;					AND for magnetic field strength (index 1), to vary between 0 and 600 gauss
;                COMPONENT  = 0 -> first component, 1-> 2nd component, etc. NO!!
;            See CHECK_PARAM subroutine for DEFAULT LIMITS
;   VARLIMITS: The maximum change to be made in the model parameter can be imposed with VARLIMITS
;              N-Array of 4 elements [PAR_INDEX,SET,MAXLIMIT_LEFT,MAXLIMIT_RIGHT,COMPONENT].
;     	PAR_INDEX: scalar value indicating the model parameter (index) to change
;       SET: 1-> ON, 0 -> OFF
;       MAXLIMIT_LEFT,MAXLIMIT_RIGHT: Left and right maximum change to be made
;            Example:
;                 VARLIMITS = [2, 1, -5d0, 5d0]
;					Limits the change to a maximum value of 5 km (for velocity, index 2, set ON)
;                    in both directions
;            The usage is similar to PARLIMITS
;   AC_RATIO: This keywork is for coupling the eta0 parameters of the 630 nm lines
;		using a polinomial function. More information in "Simultaneous ME inversion of the Fe I 630 nm lines"
;		Orozco Su\'arez, D., Bellot Rubio, L.R., and Del Toro Iniesta, J.C., A&ARN, inpress, (2010)
;
; OUTPUTS:
;   STOKESPROF: Array with the output (if the /SYNTHESIS keyword is set) Stokes profiles.
;      Dimensions: (N_ELEMENTS(AXIS), 4)
;   YFIT: Array with the fit Stokes profiles.
;      Dimensions: (N_ELEMENTS(AXIS), 4)
;   ERR: Array with errors in the model parameters.
;      Dimensions: N_ELEMENTS(MODEL)
;   CHISQR: Scalar with the chi^2 merit function of the inversion
;   ITER: Counter with the last iteration number
;   RFS: Array with the response functions to model parameter perturbations
;      Dimensions: (N_ELEMENTS(AXIS), 11, 4)
;   GETSHI: Array with chi^2 or each Stokes parameter.
;      Dimensions: (4)
;
; COMMON BLOCKS:
;   QUANTIC
;
; EXAMPLE:
;   Initialize MILOS code
;      IDL> INIT_MILOS,'63016302',WL

;   INIT MODEL PARAMETER FOR INVERSION OF SYNTHESIS
;      IDL> S0=0.2d0
;      IDL> S1=0.8d0
;      IDL> Eta0=12.5d0
;      IDL> Strength=1200d0
;      IDL> Gamma=20d0
;      IDL> Azimuth=20d0
;      IDL> Vlos=0.25d0  ;km/s
;      IDL> Macro=2.5d0
;      IDL> Lambdadopp=0.09d0
;      IDL> Damp=0.09d0
;      IDL> Alpha=0d0
;      IDL>  INIT_MODEL=[Eta0,Strength,Vlos,Lambdadopp,Damp,Gamma,Azimuth,S0,S1,Macro,Alpha]
;
;   Define an axis
;      IDL> Init_landa = Wl(1) - 0.3           ;Initial lambda
;      IDL> Step = 5d0 / 1000d0                ;Step in mA
;      IDL> Points=150.                        ;Number of wavelength points
;      IDL> AXIS = Init_lambda + Dindgen (Points) * Step
;
;   Example of synthesis
;      IDL> Milos, WL, AXIS, INIT_MODEL, STOKESPROF, /Synthesis
;
;   Example of how to make synthesis and get the RFs
;      IDL> Milos, WL, AXIS, INIT_MODEL, STOKESPROF, RFS = rfs
;
;   Example of inversion: asummed INIT_MODEL is different
;      IDL> FIX = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.]
;      IDL> WEIGHT = [1.,1.,1.,1.]
;      IDL> Milos, WL, AXIS, INIT_MODEL, STOKESPROF, YFIT=yfit,$
;           WEIGHT=weight,FIX=fix,/INVERSION
;
; NOTES:
;      Milne-Eddington inverting more than one spectral line in a row is a little risky.
;      It increases considerably the number of free parameters. An excellent opportunity
;      is found when two (or more) lines are very close in wavelength (in order to use the
;      same source function) and belong to the same species multiplet (in order for both
;      lines to be described with the same thermodynamic parameters plus a relative strength
;      for the lines) as is the case for the well-known Fe I pair of lines @630 nm.
;
;      The stray-light profile is usually assumed to be non polarized. Hence, it should
;      have the last three dimensions equal to zero. The program is however prepared
;      with the more general case of a polarized stray-light profile.
;
;      The filter profile must have the same spectral resolution than the Stokes profiles.
;
;      *IMPORTANT*
;      MILOS DOES NOT WORK IF NOT INITIALIZED CALLING INIT_MILOS
;      Refer to MILOS documentation or to init_milos routine header for further info.
;
; CALLED ROUTINES
;   MIL_SINRF: Synthesizes the Stokes profiles without calculating response functions
;   LM_MILS: Carries out the inversion through the Levenberg-Marquardt algorithm
;   ME_DER: Synthesizes the Stokes profiles and the response functions
;   TYPE
;
; REFERENCES:
;   Orozco Suarez, D., and del Toro Iniesta, J.C., A&A 462, 1137 (2007)
;   Orozco Suarez, D., PhD Thesis, 2008 (ISBN: )
;   Orozco Suarez, D., Bellot Rubio, L.R., and del Toro Iniesta, J.C., A&ARN, inpress (2010)

; MODIFICATION HISTORY:
;   First beta version created, D Orozco Suarez (DOS) and J.C. del Toro Iniesta (JTI), 2007
;   First beta documentation version, JTI, June 2008
;   Major changes: improved documentation and code. Feb 2009. DOS
;      Changed TOPRINT keyword. It is now called QUIET. 24 Feb
;      Added MU keyword. Heliocentric angle cos(mu). 24 Feb
;      Moved part of the code from LM_MILS.pro to execute DOPLOT keyword here. 24 Feb.
;      Added error handling. 26 Feb
;      Added parameter constrain keyworks PARLIMITS and VARLIMITS. 27 Feb
;   Fixed bug in var. howmany. 21 Jan, 2010. DOS
;   Added keywork AC_RATIO. 21 Jan, 2010. DOS
;   Added local straylight weight in merit function (See A. Asensio Ramos 2010) Aug, 2011. DOS
;   Added two components in Nov. 2011. DOS
;   Added numerical derivatives in Dec. 2011. DOS
;
;-
; Original MILOS, Copyright (C) 2004-2011,  D. Orozco Suarez
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-


pro MILOS, WLI, AXIS, MODEL, STOKESPROF, YFIT=yfit, ERR=err,$
    CHISQR=chisqr,ITER=iter,SLIGHT=slight,TOPLIM=toplim,$
    TRIPLET=triplet,QUIET=quiet,MITER=miter,$
    WEIGHT=weight,FIX=fix,SIGMA=sigma,$
    SYNTHESIS=synthesis,INVERSION=inversion,RFS=rfs,ILAMBDA=ilambda,$
    FILTER=filter,NOISE=noise,POL=pol,GETSHI=getshi,DOPLOT=doplot,MU=mu,$
	PARLIMITS=parlimits,VARLIMITS=varlimits,AC_RATIO=ac_ratio,MLOCAL=MLOCAL,$
	N_COMPONENTS=n_components,numerical=numerical,iter_info = iter_info,$
  use_svd_cordic = use_svd_cordic,ipbs=ipbs

COMMON QUANTIC,C_N

;Error handling
;on_error,2

;check whether the COMMON block quantic exist!

a=type(c_n)
if a eq 0 then begin
	print,' '
	print,' ---------------------  ATTENTION -------------------------------'
	print,' It seems that the commom block QUANTIC has not been definet jet.'
	print,' Check whether INIT_MILOS has been executed before runing MILOS.'
	print,' Check README file for more info on how to run MILOS.'
	print,' ----------------------------------------------------------------'
	print,' '
	return
endif

;check minimum of input parameters
nparams = N_PARAMS()
if nparams lt 4 then begin
	print,' '
	print,' ------------------- CALLING  ERROR -----------------------------'
	print,' Not enought input parameters:'
	print,'    MILOS, WLI, AXIS, MODEL, STOKESPROF'
	print,' ----------------------------------------------------------------'
	print,' '
	return
endif

;if the code does not stop until here everything should be ok'
;before continuing we check whether variable sizes correspond each'

;WLI has to be 2,3... elements vector depending on QUANTEN
c2 = size (C_N)
if (wli(0) gt c2(1)) or (wli(0) gt n_elements(wli)-1) then begin
	print,' '
	print,' --------------------- INPUT  ERROR ----------------------------------------'
	print,' Input parameter WLI does not have the correct dimensions      '
	print,'    Required: Array ( Number of lines , central wavelength 1, c. w. 2, ... )                                       '
	print,'    Current input: ', wli
	print,' ---------------------------------------------------------------------------'
	print,' '
	return
endif

if not(keyword_set(n_components)) then n_comp=1 else n_comp = n_components

;MODEL has to be 11*c elements vector
c1 = size (MODEL)
if (c1(0) ne 1) or (c1(1) ne 11*n_comp) then begin
	print,' '
	print,' --------------------- INPUT  ERROR -----------------------------'
	print,' Input parameter MODEL does not have the correct dimensions      '
	print,'    Required: Array ( 11 , n_comp )                                       '
	print,'    Current input: Array ( ' , c1(0:1) , ' ) '
	print,' ----------------------------------------------------------------'
	print,' '
	return
endif

; -----------------------------------------------------------------------
; Synthesis of the Stokes profiles without calculating response functions
; -----------------------------------------------------------------------

  if keyword_set(synthesis) then begin

    mil_sinrf,MODEL,wli,AXIS,STOKESPROF,TRIPLET=triplet,SLIGHT=slight,$
    FILTER=filter,MU=mu,AC_RATIO=ac_ratio,N_COMP = n_comp,ipbs=ipbs

  if keyword_set(doplot) then begin
    !p.multi=[0,2,2]
    plot,AXIS,STOKESPROF(*,0),thick=2,title='I/Ic'
    plot,AXIS,STOKESPROF(*,3),thick=2,title='V/Ic'
    plot,AXIS,STOKESPROF(*,1),thick=2,title='Q/Ic'
    plot,AXIS,STOKESPROF(*,2),thick=2,title='U/Ic'
  endif

; -----------------------------------------------------------------------
; Inversion of Stokes profiles
; -----------------------------------------------------------------------

  endif else if keyword_set(inversion) then begin

     pn = [1,0,1,1,1,0,0,1,1,0,1]
 if not(keyword_set(fix)) then begin
 	fix = fltarr(11*n_comp)
 	fix[0:8] = 1d0
 	if n_comp gt 1 then for i=1, n_comp-1 do fix[indgen(11)+i*11] = pn
 endif else begin
    cc = size(fix)
    if (cc(0) ne 1) or (cc(1) ne 11*n_comp) then begin
 	print,'error in fix par'
 	return
 	endif
 endelse

;AXIS and STOKESPROF has to be same column dimension
c1 = size (AXIS)
c2 = size (STOKESPROF)
if (c1(1) ne c2(1)) or (c2(2) ne 4) then begin
	print,' '
	print,' --------------------- INPUT  ERROR ------------------------------------------'
	print,' Input parameters AXIS and STOKESPROF do not have the correct dimensions      '
	print,'    Required: AXIS -> Array ( wavelength samples )                            '
	print,'    Required: STOKESPROF -> Array ( wavelength samples , 4 )                  '
	print,'    Current input: AXIS -> Array ( ' , c1(1) , ' ) '
	print,'    Current input: STOKESPROF -> Array ( ' , c2(1) , ' , ' , c2(2) , ' ) '
	print,' -----------------------------------------------------------------------------'
	print,' '
	return
endif

if keyword_set(slight) then begin
c3 = size (SLIGHT)
if (c2(1) ne c3(1)) or (c2(2) ne c3(2) and (c3(0) ne 1)) then begin
	print,' '
	print,' --------------------- INPUT  ERROR ------------------------------------------'
	print,' Input parameters SLIGHT and STOKESPROF do not have the same dimensions      '
	print,'    Required: SLIGHT -> Array ( wavelength samples , 4)                            '
	print,'    Required: STOKESPROF -> Array ( wavelength samples , 4 )                  '
	print,'    Current input: SLIGHT -> Array ( ' , c3(1) , ' , ' , c3(2) , ' ) '
	print,'    Current input: STOKESPROF -> Array ( ' , c2(1) , ' , ' , c2(2) , ' ) '
	print,' -----------------------------------------------------------------------------'
	print,' '
	return
endif
endif

;SET PLIMITS
  index = indgen(n_comp)*11
	PLIMITS = replicate( {set:0, limits:[0d0,0d0]} , 11* n_comp) ;define variable plimits
	PLIMITS[*].SET = 1  ; Activate all limits
	PLIMITS[index].LIMITS = [1d0,2500d0]  ; Eta0
	PLIMITS[index+1].LIMITS = [0d0,4500d0]  ; Magnetic field
	PLIMITS[index+2].LIMITS = [-100d0,100d0]  ; Velocity (km/s)
	PLIMITS[index+3].LIMITS = [1d-4,6d-1]  ; Doppler width
	PLIMITS[index+4].LIMITS = [1d-4,1d1]  ; Damping
	PLIMITS[index+5].LIMITS = [0d0,180d0]  ; Inclination
	PLIMITS[index+6].LIMITS = [0d0,360d0]  ; Azimuth
	PLIMITS[index+7].LIMITS = [1d-4,1d1]  ; S0
	PLIMITS[index+8].LIMITS = [1d-4,1d1]  ; S1
	PLIMITS[index+9].LIMITS = [0d0,4d0]  ; Macroturbulence
	PLIMITS[index+10].LIMITS = [0d0,1d0]  ; Magnetic filling factor

if keyword_set(PARLIMITS) then begin
	Sz = size(PARLIMITS)
	if Sz(1) ne 4 then begin
		print,' '
		print,' --------------------- INPUT  ERROR -----------------------------'
		print,' Input parameter PARLIMITS does not have the correct dimensions      '
		print,'    Required: Array ( 4 )                                       '
		print,'    Current input: Array ( ' , c1(1) , ' ) '
		print,' ----------------------------------------------------------------'
		print,' '
		return
	endif
	check_size = Size(PARLIMITS)
  If check_size(0) eq 1 then Howmany = 1 else Howmany = check_size(2)
	FOR i=0, howmany -1 do begin
		PLIMITS[PARLIMITS(0,i)].LIMITS = [ PARLIMITS(2,i) ,  PARLIMITS(3,i) ]
		PLIMITS[PARLIMITS(0,i)].SET = PARLIMITS(1,i)
	ENDFOR

ENDIF

    VLIMITS = replicate( {set:0, limits:[0d0,0d0]} , 11* n_comp) ;define variable vlimits
    VLIMITS[index+1].SET = 1
    VLIMITS[index+1].LIMITS = [-1000d0,1000d0]
    VLIMITS[index+5].SET = 1
    VLIMITS[index+5].LIMITS = [-90d0,90d0]
    VLIMITS[index+6].SET = 1
    VLIMITS[index+6].LIMITS = [-90d0,90d0]
    VLIMITS[index+9].SET = 1
    VLIMITS[index+9].LIMITS = [0d0,4d0] ;MACRO

if keyword_set(VARLIMITS) then begin
	Sz = size(VARLIMITS)
	if Sz(1) ne 4 then begin
		print,' '
		print,' --------------------- INPUT  ERROR -----------------------------'
		print,' Input parameter VARLIMITS does not have the correct dimensions      '
		print,'    Required: Array ( 4 )                                       '
		print,'    Current input: Array ( ' , c1(1) , ' ) '
		print,' ----------------------------------------------------------------'
		print,' '
		return
	endif
	check_size = Size(VARLIMITS) ;FIXED Bug PARLIMITS -> VARLIMITS (24-Oct-2015)
	If check_size(0) eq 1 then Howmany = 1 else Howmany = check_size(2)
	FOR i=0, howmany -1 do begin
		VLIMITS[VARLIMITS(0,i)].LIMITS = [ VARLIMITS(2,i) ,  VARLIMITS(3,i) ]
		VLIMITS[VARLIMITS(0,i)].SET = VARLIMITS(1,i)
	ENDFOR

ENDIF

    LM_MILS,Wli, AXIS, STOKESPROF, MODEL, YFIT, ERR,$
      CHISQR,ITER,SLIGHT=SLIGHT,TOPLIM=toplim,$
      TRIPLET=triplet,QUIET=quiet,MITER=miter,$
      WEIGHT=weight,FIX=fix,SIGMA=sigma,$
      FILTER=filter,ILAMBDA=ilambda,NOISE=noise,POL=pol,$
	  GETSHI=getshi,MU=mu,PLIMITS=plimits,VLIMITS=vlimits,$
	  AC_RATIO=ac_ratio,MLOCAL=MLOCAL,N_COMP=n_comp,numerical=numerical,$
	  iter_info = iter_info,use_svd_cordic = use_svd_cordic,ipbs=ipbs

  if keyword_set(doplot) then begin
    !p.multi=[0,2,2]
    plot,AXIS,STOKESPROF(*,0),thick=2,title='I/Ic',/ynoz
    oplot,AXIS,yfit(*,0),thick=3,color=4,line=2
    plot,AXIS,STOKESPROF(*,3),thick=2,title='V/Ic'
    oplot,AXIS,yfit(*,3),thick=3,color=4,line=2
    plot,AXIS,STOKESPROF(*,1),thick=2,title='Q/Ic'
    oplot,AXIS,yfit(*,1),thick=3,color=4,line=2
    plot,AXIS,STOKESPROF(*,2),thick=2,title='U/Ic'
    oplot,AXIS,yfit(*,2),thick=3,color=4,line=2
  endif

; -----------------------------------------------------------------------
; Synthesis of the Stokes profiles and of the response functions (derivatives for inversion)
; -----------------------------------------------------------------------

  endif else if ARG_PRESENT(RFS) then begin

    ME_DER,MODEL,WLI,AXIS,STOKESPROF,RFS,TRIPLET=triplet,ipbs=ipbs,$
    SLIGHT=slight,FILTER=filter,MU=mu,AC_RATIO=ac_ratio,N_COMP=n_comp,numerical=numerical

  if keyword_set(doplot) then begin
    !p.multi=[0,2,2]
    plot,AXIS,STOKESPROF(*,0),thick=2,title='I/Ic',/ynoz
    plot,AXIS,STOKESPROF(*,3),thick=2,title='V/Ic'
    plot,AXIS,STOKESPROF(*,1),thick=2,title='Q/Ic'
    plot,AXIS,STOKESPROF(*,2),thick=2,title='U/Ic'
  endif

  endif else begin

	print,'no option has been set: /inversion, /synthesis, or rfs = rfs'

	return

  endelse

end
