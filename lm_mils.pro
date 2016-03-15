;+
; NAME:
;	LM_MILS
;
; AUTHOR:
;	D. Orozco Su�rez
;						National Astronomical Observatory of Japan,�
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
;	J.C. del Toro Iniesta
;						Instituto de Astrof�sica de Andaluc�a (CSIC)
;						Apdo de Correos 3004, 18080 Granada, SPAIN
;						jti@iaa.es
;
; PURPOSE: Core routine of the MILOS procedure. It performs a non-linear least squares fit
;          using the Levenberg-Marquardt algorithm
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: LM_MILS, WLI, AXIS, STOKESPROF, p_i, yfit[, err, chisqf,iter,SLIGHT=slight,TOPLIM=toplim,$
;    TRIPLET=triplet, QUIET=quiet,MITER=miter, WEIGHT=weight,FIX=fix,SIGMA=sigma,$
;    FILTER=filter, ILAMBDA=ilambda, NOISE=noise, POL=pol, GETSHI=getshi, MU=mu,$
;    PLIMITS=plimits,VLIMITS=vlimits,N_COMPONENTS=n_components,$
;    AC_RATIO=ac_ratio,MLOCAL=mlocal]
;
; INPUTS (mandatory):
;           WLI: A scalar or array with the central wavelength(s) of the line(s)
;                in Angstroms
;                Dimensions: (number of lines, central wavelengths)
;                Examples: [1, 6302.5]; [2, 6301.5, 6302.5]
;           AXIS: Wavelength axis in Angstroms. Just one for all the lines in the sample
;                 It should be a double precision array
;           STOKESPROF: Array with the input Stokes profiles
;              Dimensions: (N_ELEMENTS(AXIS), 4)
;          P_I: Array with the 11*c model-atmosphere parameters
;               eta0 = line-to-continuum absorption coefficient ratio
;               B = magnetic field strength       [Gauss]
;               vlos = line-of-sight velocity     [km/s]
;               dopp = Doppler width              [Angstroms]
;               aa = damping parameter
;               gm = magnetic field inclination   [deg]
;               az = magnetic field azimuth       [deg]
;               S0 = source function constant
;               S1 = ource function gradient
;               mac = macroturbulent velocity     [km/s]
;               alpha = filling factor of the magnetic component [0->1]
;               Whenever n_components gt 1 p_i dimension shoul be scaled accordingly
; INPUTS (optional):
;           SLIGHT: Array with the stray-light profile
;                   Dimensions: (N_ELEMENTS(AXIS), 4)
;           TOPLIM: Optional minimum relative difference between two succesive merit-function values
;           MITER: Maximum number of iterations allowed
;           WEIGHT: Array with arbitrary weights for the Stokes profiles
;                   Dimensions: (N_ELEMENTS(AXIS), 4) (as many as the profiles)
;                   It can be used, for example, to neglect telluric lines in the fits
;           FIX: Array with as many (integers) zeroes or ones as parameters to be inverted
;                Dimensions: (11)
;                Example: [1,1,1,1,1,1,1,1,1,0,0] means "invert everything but mac and alpha"
;                Default: all set to 1 but macro and filling factor
;           SIGMA: Array with noise values for the profiles (standard deviations)
;                  Dimensions: (4)
;           ILAMBDA: Initial value for the Levenberg-Marquardt's fudge parameter
;           FILTER: Width (in mA) of a Gaussian filter profile or array with the instrument
;                   profile samples. Dimensions: N_ELEMENTS(AXIS)
;           NOISE: Scalar with the noise values for the Stokes profiles
;           POL: Polarization threshold. If all Q, U, and V are below POL, then the
;                vector magnetic field is set to zero
;           MU: Scalar containing the cosine of the heliocentric angle
;           PLIMITS: Limitations on the model parameter changes can be imposed throught this keywork
;                    Array of structures. Dimension (11).
;                    Each model parameter is associated with one element of the array
;                    The parameter can be bounded on the lower/upper side.
;            	PLIMITS.SET: scalar value indicating whether the parameter constraint in on/off
;                            1-> on, 0-> off
;               PLIMITS.LIMITS: two element array with the lower/upper sides
;               DEFAULT LIMITS ARE:
;           	PLIMITS = replicate( {set:0, limits:[0d0,0d0]} , 11) ;define variable plimits
;           	PLIMITS[*].SET = 1  ; Activate all limits
;               PLIMITS[0].LIMITS = [1d0,2500d0]  ; Eta0
;               PLIMITS[1].LIMITS = [0d0,4500d0]  ; Magnetic field
;               PLIMITS[2].LIMITS = [-20d0,20d0]  ; Velocity (km/s)
;               PLIMITS[3].LIMITS = [1d-4,6d-1]  ; Doppler width
;               PLIMITS[4].LIMITS = [1d-4,1d1]  ; Damping
;               PLIMITS[5].LIMITS = [0d0,180d0]  ; Inclination
;               PLIMITS[6].LIMITS = [0d0,180d0]  ; Azimuth
;               PLIMITS[7].LIMITS = [1d-4,1d1]  ; S0
;               PLIMITS[8].LIMITS = [1d-4,1d1]  ; S1
;               PLIMITS[9].LIMITS = [0d0,4d0]  ; Macroturbulence
;               PLIMITS[10].LIMITS = [0d0,1d0]  ; Magnetic filling factor
;           VLIMITS: The maximum change to be made in the model parameter can be imposed with VLIMITS
;                    Array of structures. Dimension (11).
;                    Each model parameter is associated with one element of the array
;                    The maximum change can be set in both directions
;						VLIMITS[*].LIMITS(0) > DELTA < VLIMITS[*].LIMITS(1)
;            	VLIMITS.SET: scalar value indicating whether a parameter is limited or not
;                            1-> on, 0-> off
;               VLIMITS.LIMITS: two element array with the maximum variation when decreasing/increasing
;                               the model parameter.
;               DEFAULT LIMITS ARE:
;           	VLIMITS = replicate( {set:0, limits:[0d0,0d0]} , 11) ;define variable vlimits
;           	VLIMITS[1,2,5,6].SET = 1  ; Activate four limits
;               VLIMITS[1].LIMITS = [-300d0,300d0]  ; Magnetic field
;               VLIMITS[2].LIMITS = [-10d0,10d0]  ; Velocity (km/s)
;               VLIMITS[5].LIMITS = [-15d0,15d0]  ; Inclination
;               VLIMITS[6].LIMITS = [-15d0,15d0]  ; Azimuth
;			N_COMP : Number of components (besides stray-light)
;           AC_RATIO: See Orozco Suarez et al 2010. Lines coupling.
;           M_LOCAL: Number of pixels taken to average the stray light porfile (box)
; KEYWORD PARAMETERS:
;           TRIPLET: Set this keyword to do the calculations for the effective triplet
;           QUIET: Set this keyword to cancel any printed message on the screen
;
; OUTPUTS:
;           P_I: Array with the 11*c model-atmosphere parameters (overwritten)
;           YFIT: Array with the fit Stokes profiles.
;                       Dimensions: (N_ELEMENTS(AXIS), 4)
;           ERR: Array with errors in the model parameters.
;                       Dimensions: N_ELEMENTS(MODEL)
;                       Assuming the true reduced chi-squared value is unity.
;                       ERR = SQRT(COVAR(DIAG)*OCHISQR*4*NPOINTS/NFREE)
;           CHISQF: Scalar with the chi^2 merit function of the inversion
;           ITER: Counter with the last iteration number
;           GETSHI: Array with chi^2 or each Stokes parameter.
;                       Dimensions: (4)
;
; COMMON BLOCKS:
;
; NOTES:
;
; CALLED ROUTINES
;           WEIGHTS_INIT: Gives initial values to the Stokes profile weights
;           ME_DER: Synthesizes the Stokes profiles and their derivatives with respect to
;                   model parameters
;           COVARM: Calculates the gradient and approximate-Hessian (covariance) matrices
;           MIL_SVD: Inverts the covariance matrix through the singular-value-decomposition
;                    method
;           CHECK_PARAM: Checks that parameters are within limits
;
; MODIFICATION HISTORY:
; 	First beta version created, D Orozco Su�rez (DOS) and J.C. del Toro Iniesta (JTI), 2007
;   First beta documentation version, JTI, June 2008
;	Improved documentation. 24 Feb 2009. DOS
;	Changed TOPRINT keyword. It is now called QUIET. 24 Feb, 2009. DOS
;	Removed DOPLOT keyword. It is now in milos.pro. 24 Feb, 2009. DOS
;   Added keywork PLIMITS=plimits to pass to check_model procedure. 27 Feb, 2009. DOS
;   Added keywork VLIMITS=vlimits to constraint Delta. 27 Feb, 2009. DOS
;   Added local straylight weight in merit function (See A. Asensio Ramos 2010) Aug, 2011. DOS
;   Added two components in Nov. 2011. DOS
;
;-

pro LM_MILS, WLI, AXIS, STOKESPROF, p_i, yfit, err, chisqf,iter,slight=slight,toplim=toplim,$
    triplet=triplet, QUIET=quiet,miter=miter, weight=weight,fix=fix,sigma=sigma,$
    filter=filter, ilambda=ilambda, noise=noise, pol=pol, getshi=getshi, $
	PLIMITS=plimits,VLIMITS=vlimits,MU=mu,AC_RATIO=ac_ratio,MLOCAL=mlocal,$
	N_COMP=n_comp,numerical=numerical,iter_info = iter_info,use_svd_cordic = use_svd_cordic,$
  ipbs=ipbs

; Enviromental parameters
  prt=keyword_set(QUIET)
  if not(keyword_set(miter)) then miter=50  ; Maximum number of iterations
  if not(keyword_set(ilambda)) then ilambda=0.1d0;0.1d0 ; Assumes close to global minimum
  if not(keyword_set(toplim)) then toplim=1d-12 ; Very low, It may stop by MITER
;  if not(keyword_set(slight)) then print,'No stray light';slight=0 ; No stray light
;  if not(keyword_set(filter)) then print,'No inst. filter profile'
  if not(prt) then print,'triplet='+string(keyword_set(triplet))
  if not(prt) then print,'Iterations='+string(miter)
  if not(keyword_set(noise)) then noise=1d-3 ; Default noise level
  if not(keyword_set(n_comp)) then n_comp=1d0  ; Maximum number of iterations
;  n_comp = n_comp*1.

  nterms = N_ELEMENTS(P_I)*1d0
  NFIXSS = total(fix)
  flambda=ILAMBDA   ; Initial fudge parameter
  DIAG = INDGEN(NTERMS)*(NTERMS+1) ;SUBSCRIPTS OF DIAGONAL ELEMENTS
  FIXED=FIX*1d0
  NAXIS=N_ELEMENTS(AXIS)
  valid=WHERE(STOKESPROF ne 0.,HAY) ; If Y = 0 then the samples are not valid (e.g. blends)
  nfree=DOUBLE(N_ELEMENTS(VALID)-NFIXSS)  ; Number of degrees of freedom
  npoints=(SIZE(STOKESPROF))(1)
  iter_info = {lmb:dblarr(miter+1),iter:0,citer:0,conv_crit:dblarr(miter+1),Params_stored:fltarr(nterms,miter+1)}

  if NFREE le 0 then begin
	PRINT,'Not enough points'
	return
  endif
  P_M=DBLARR(NTERMS) ;New model aprameters

  goodc=0
  max_stored=Miter+1
  Params_stored=fltarr(nterms,max_stored)
  conv_crit = intarr(max_stored)
  ;TO BE IMPLEMENTED
  ;   ;INITIALIZATION
  ;   p_i(6)=atan(max(STOKESPROF(*,1))/max(STOKESPROF(*,2)))/2.*180/!dpi
  ;   signo=total(STOKESPROF(0:NPOINTS/2,3))-total(STOKESPROF(NPOINTS/2+1:*,3))
  ;   if signo lt 0 then p_i(5)=100. else p_i(5)=10.

  ;INITIALIZING WEIGHTS
  WEIGHTS_INIT,NAXIS,STOKESPROF,w,sig,weight=weight,SIGMA=sigma,NOISE=noise
  CLANDA=0 ;SOME DEFINITIONS TO CONTROL THE FUDGE PARAMETER
  REPITE=1
  pillado=0
  ITER=0  ;DEFINING FIRST ITER
  LANDA_STORE=DBLARR(MITER+1) ;Array to store lambda variation with iter

  ;POLARIZATION THRESHOLD (If Stokes q, u, and v are below it, set B, GAMMA, and AZI to cero
  index=indgen(n_comp)
  IF keyword_set(pol) then begin
    limp=pol
    IF (max(abs(STOKESPROF(*,1))) lt limp) and  (max(abs(STOKESPROF(*,2))) lt limp) $
      and (max(abs(STOKESPROF(*,3))) lt limp) then begin
      p_i(1)=0. & p_i(5)=0. & p_i(6)=0
      fixed(11*index+1)=0 & fixed(11*index+5)=0 & fixed(11*index+6)=0
    endif
  endif

  ;Compute the ME RFs and the initial synthetic Stokes profiles for INIT MODEL
  ME_DER,p_i,wli,AXIS,yfit,pder,TRIPLET=triplet,SLIGHT=slight,FILTER=filter,$
    AC_RATIO=ac_ratio,N_COMP=n_comp,numerical=numerical,ipbs=ipbs

  ;The derivatives with respect to fixed parameters are set to cero
  fxx=where(fixed eq 1)
  for I=0,nterms-1 DO PDER(*,I,*)=PDER(*,I,*)*FIXED(I)

; Take into account the local stray light in the merit function as in Asensio Ramos and Manso Sainz, 2010.
If keyword_set(mlocal) then begin
	;mlocal is the number of averaged pixels
	rho = dblarr(4)
	drho = dblarr(4)
	rho[0] = reform( 1./6. * total( (stokesprof(0:5,0) - mean(stokesprof(0:5,0)) ) * (slight(0:5,0) - mean(slight(0:5,0)) ) ) )
	rho[1:3] = 0d0
	drho[0] = 2d0*(1d0 - p_i(10))/(mlocal*1d0)*sig[0] - 2*rho[0]
	term1 = 1d0 + (1d0 - p_i(10))^2d0/(mlocal*1d0) ;1-(1-alpha)^2/M
	sigo = sig
	sig = sigo*term1 - 2d0*(1d0 - p_i(10))*rho
endif

  ;Computes the gradient, the covariance and the initial merit function value
  COVARM,W,SIG,YFIT,STOKESPROF,PDER,NTERMS,NFREE,BETA,ALPHA,CHISQR,drho=dhro
;  flambda = CHISQR
;  flambda = 1.
;  Jope = max(alpha(diag))*1.e-8
;  if jope gt 1 then cuantos = (n_elements(digits_get(float(Jope)))) else flambda=0.1
;  for i=0,cuantos do flambda = flambda*10 ;< 1000000.
;  print,flambda

  ;Store old chisqr value
  OCHISQR=CHISQR

  ;SOME DEBUGGING STAFF
  ;   flambda=max(ALPHA(DIAG))*1d-6
  ;   ya=1 & ll=1.
  ;   while (ya) do begin
  ;       if ceil(flambda/ll) eq 1. then begin
  ;           flambda=1.>double(ll)<10000.
  ;           ya=0
  ;           endif
  ;       ll=ll*10
  ;   endwhile
  ;print,flambda

   P_M=P_I
   yfitt = yfit

  ;BEGIND THE ITERATION LOOP FOR THE LM OPTIMIZATION ALGORITHM
REPEAT BEGIN

	iter_info.lmb[ITER] = FLAMBDA

    IF (FLAMBDA GT 1d25) OR (FLAMBDA LT 1d-25) THEN CLANDA=1 ;Cond to Flambda !!!!!!!!!!
    COVAR=ALPHA
    COVAR(DIAG)=COVAR(DIAG)*(1.+FLAMBDA)
    BETAD=BETA

;       INVERT COVARIANCE MATRIX TO FIND NEW MODEL PARAMETERS.

    MIL_SVD,COVAR,BETAD,DELTA,LOSW,use_svd_cordic = use_svd_cordic

	;CHECK MAX VARIATION OF DELTA AND SET NEW MODEL P_M

    for i=0,11*n_comp-1 do IF VLIMITS(i).SET EQ 1 THEN $
	DELTA(i) = VLIMITS(i).LIMITS(0) > DELTA(i) < VLIMITS(i).LIMITS(1)

    P_M(fxx) = P_I(fxx) - DELTA(fxx) ;NEW PARAMS

    ;CHECK PARAMETERS
    CHECK_PARAM,P_M,PLIMITS=plimits,n_comp=n_comp

	;EVALUATE FUNCTION
    mil_sinrf,p_m,wli,AXIS,yfit,triplet=triplet,slight=slight,filter=filter,$
      AC_RATIO=ac_ratio,n_comp=n_comp,ipbs=ipbs

	;new chisqr
    CHISQR = TOTAL(TOTAL((YFIT-STOKESPROF)^2d0*W,1,/double)/SIG^2d0,/double)/NFREE

;print,delta
;TEST!!!!
goto,notest
    COVARM,W,SIG,YFIT,STOKESPROF,PDER,NTERMS,NFREE,BETA,ALPHA,CHISQRB,drho=dhro
    COVAR=ALPHA
    COVAR(DIAG)=COVAR(DIAG)*(1.+FLAMBDA)
    BETAD=BETA
    MIL_SVD,COVAR,BETAD,DELTA,LOSW,use_svd_cordic = use_svd_cordic
    for i=0,11*n_comp-1 do IF VLIMITS(i).SET EQ 1 THEN $
  	 DELTA(i) = VLIMITS(i).LIMITS(0) > DELTA(i) < VLIMITS(i).LIMITS(1)
    P_M(fxx) = P_M(fxx) - DELTA(fxx)
    CHECK_PARAM,P_M,PLIMITS=plimits,n_comp=n_comp
    mil_sinrf,p_m,wli,AXIS,yfit,triplet=triplet,slight=slight,filter=filter,$
      AC_RATIO=ac_ratio,n_comp=n_comp,ipbs=ipbs
    CHISQR = TOTAL(TOTAL((YFIT-STOKESPROF)^2d0*W,1,/double)/SIG^2d0,/double)/NFREE
notest:
;TEST!!!!
;print,delta

      ;****CONVERGENCE CONDITION *****
    IF CHISQR-OCHISQR lt 0. then begin ;;FIT GOT BETTER
      ;****CONVERGENCE CONDITION *****
      conv_crit[ITER] = 1
      iter_info.conv_crit[ITER] = conv_crit[ITER]

      if (ABS((OCHISQR-CHISQR)*100./CHISQR) LT TOPLIM) OR (CHISQR lt 0.0001) $
        THEN CLANDA = 1 ;Stoping criteria

      ;DECREASE FLAMBDA BY FACTOR OF 10

;help,delta,betad,OCHISQR,CHISQR,flambda,transpose(flambda*Delta - BETAD),(0.5*DELTA##transpose(flambda*Delta - BETAD));
;;	rho = abs((OCHISQR-CHISQR)/(0.5*DELTA##transpose(flambda*Delta - BETAD)))
;;	rho = rho[0] * 100.
	;print, 'rho',rho

;;pbeta=2d0;6d0
;;pp=3d0
;rho = findgen(1000)/999.
;plot,rho,1-(pbeta-1)*(2*rho-1)^pp
;print, 'fidge', ((pbeta-1)*(2*rho-1)^pp+10)
;flambda = flambda / ((pbeta-1)*(2*rho-1)^pp+10) ;*1-(pbeta-1)*(2*rho-1)^pp

      ;FIT GOT BETTER SO DECREASE FLAMBDA BY FACTOR OF 10
      flambda=flambda/10d0

            ;RNOISE_PROBLEM
;      IF FLAMBDA LT 1 THEN BEGIN
;      II = 1d0
;      REPEAT II = II * 10. UNTIL FLAMBDA-ROUND(FLAMBDA*II) LT 0
;      II = II*10d0
;      FLAMBDA = ROUND(FLAMBDA*II)/II
;      ENDIF

      ;flambda=flambda/(10d0);+(OCHISQR-CHISQR)/CHISQR)
      ;      ;SOME DEBUGGING STAFF
      ;if repite then begin
      ;    if flambda gt 1.e-4 then begin
      ;        flambda=flambda*0.1d0
      ;    endif else begin
      ;        flambda=flambda*0.5d0
      ;    endelse
      ;endif

      ;store new model parameters
      P_I(fxx)=P_M(fxx)
      ;store parameters for studing the evolution of them (opt mode)
      Params_stored(*,goodc) = p_i
      iter_info.Params_stored[*,ITER] = P_I
      ;Store the number if iterations that converged (opt mode)
      goodc = goodc + 1

	 ; if (goodc mod 10 eq 0)  then begin
	 ;  	acc,goodc - 1,fxx,params_stored,pnew
	 ;  	CHECK_PARAM,Pnew,PLIMITS=plimits
	 ;   p_i = pnew
	 ; endif

      ;info
      IF not(prt) THEN PRINT,'ITERATION =',ITER,' , CHISQR =',CHISQR,$
        '  CONVERGE - LAMBDA= ',flambda,(OCHISQR-CHISQR)/CHISQR

      ;compute new RFs and Stokes profiles
      me_der,p_i,wli,AXIS,yfit,pder,triplet=triplet,slight=slight,filter=filter,$
        AC_RATIO=ac_ratio,n_comp=n_comp,numerical=numerical,ipbs=ipbs
      for I=0,nterms-1 DO PDER(*,I,*)=PDER(*,I,*)*FIXED(I)

      ;Local stray light option (modification of chi2)
      If keyword_set(mlocal) then begin
        drho[0] = 2d0*(1d0 - p_i(10))/(mlocal*1d0)*sig[0] - 2*rho[0]
        term1 = 1d0 + (1d0 - p_i(10))^2d0/(mlocal*1d0) ;1-(1-alpha)^2/M
        sig = sigo*term1 - 2d0*(1d0 - p_i(10))*rho
      endif

      COVARM,W,SIG,YFIT,STOKESPROF,PDER,NTERMS,NFREE,BETA,ALPHA,OCHISQR,drho=dhro

      ;****CONVERGENCE CONDITION *****
    ENDIF ELSE BEGIN ;ASSUMES FIT GOT WORSE
      ;****CONVERGENCE CONDITION *****
      conv_crit[ITER] = 0
      iter_info.conv_crit[ITER] = conv_crit[ITER]
      iter_info.Params_stored[*,ITER] = P_I

      ;increase flambda by a factor 10

      flambda=flambda*10d0

;    IF (FLAMBDA GT 1d25) OR (FLAMBDA LT 1d-25) THEN CLANDA=1 ;Cond to Flambda !!!!!!!!!!

      ;RNOISE_PROBLEM
;      IF FLAMBDA LT 1 THEN BEGIN
;      II = 1d0
;      REPEAT II = II * 10. UNTIL FLAMBDA-ROUND(FLAMBDA*II) LT 0
;      II = II*100d0
;      FLAMBDA = ROUND(FLAMBDA*II)/II
;      ENDIF

       ;SOME DEBUGGING STAFF
      ;if flambda le 1.e-3 then begin
      ;    flambda=100d0*flambda
      ;endif  else begin
      ;    if flambda lt 1.e3 then begin
      ;        flambda=10d0*flambda
      ;    endif else begin
      ;        flambda=2d0*flambda
      ;    endelse
      ;endelse


      ;info
      IF not(prt) THEN PRINT,'ITERATION =',ITER,' , CHISQR =',OCHISQR,$
        '  NOT CONVERGE - LAMBDA= ',flambda,(OCHISQR-CHISQR)/CHISQR

    ENDELSE

	;Lambda storage to avoid oscillations in lambda, in such a case the inversion is finished
    LANDA_STORE(ITER)=FLAMBDA

    if (ITER gt 5) then begin
      if (conv_crit(ITER) eq conv_crit(ITER-2)) $
        AND (conv_crit(ITER-1) eq conv_crit(ITER-3)) $
        AND (conv_crit(ITER) ne conv_crit(ITER-1))then $
        flambda = flambda + flambda*1d-1
    endif

    ITER=ITER+1

ENDREP UNTIL (ITER GT MITER) OR (CLANDA EQ 1)

	iter_info.iter = ITER
	iter_info.citer = goodc

  ;END OF MAIL LOOP AND INVERSION

;print,'ITER ' , ITER, 'goodc ',goodc, ' CHI ', OCHISQR
  ;ERRORS

  COVAR=1./ALPHA
  ERR = SQRT(COVAR(DIAG)*OCHISQR*4*NPOINTS/NFREE) ; Almeida

  ;keep last merit function value
  CHISQF=OCHISQR

  ;Stokes I, Q, U, and V merit function values.
  getshi=fltarr(4)
  for ii=0,3 do getshi(ii)=TOTAL(TOTAL((YFIT(*,ii)-STOKESPROF(*,ii))^2d0*W(*,ii),1,/double)/SIG(ii)^2d0,/double)/NFREE

;  plot,LANDA_STORE;
;wait,0.2
END
