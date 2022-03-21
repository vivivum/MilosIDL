;+
; NAME:
;	LM_MILS
;
; AUTHORS:
;   D. Orozco Suarez
;   orozco@iaa.es
;	and
;   J.C. del Toro Iniesta
; ADDRESS:
;      Instituto de Astrofisica de Andalucia (CSIC)
;      Apdo de Correos 3004, 18080 Granada, SPAIN
;
; PURPOSE: Core routine of the MILOS procedure. It performs a non-linear least squares fit
;          using the Levenberg-Marquardt algorithm
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: LM_MILS, WLI, AXIS, STOKESPROF, p_i, yfit[, err, chisqf,iter,SLIGHT=slight,TOPLIM=toplim,$
;    TRIPLET=triplet, QUIET=quiet,MITER=miter, WEIGHT=weight,FIX=fix,SIGMA=sigma,
;    FILTER=filter, ILAMBDA=ilambda, NOISE=noise, POL=pol, GETSHI=getshi, MU=mu,
;    PLIMITS=plimits,VLIMITS=vlimits,N_COMP=n_comp,
;    AC_RATIO=ac_ratio,MLOCAL=mlocal,numerical=numerical,iter_info = iter_info,
;    use_svd_cordic = use_svd_cordic,ipbs=ipbs,crosst = crosst,OTHERS = OTHERS,
;    nlte=nlte,CHISQR_LIMIT=CHISQR_LIMIT

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
;   Added keyword NLTE in Nov-2020 and modified program to deal with nlte DOS
;-

pro LM_MILS, WLI, AXIS, STOKESPROF, p_i, yfit, err, chisqf,iter,slight=slight,toplim=toplim,$
    triplet=triplet, QUIET=quiet,miter=miter, weight=weight,fix=fix,sigma=sigma,$
    filter=filter, ilambda=ilambda, noise=noise, pol=pol, getshi=getshi, $
	  PLIMITS=plimits,VLIMITS=vlimits,MU=mu,AC_RATIO=ac_ratio,MLOCAL=mlocal,$
	  N_COMP=n_comp,numerical=numerical,iter_info = iter_info,use_svd_cordic = use_svd_cordic,$
    ipbs=ipbs,crosst = crosst,OTHERS = OTHERS,nlte=nlte,C_LIMIT=C_LIMIT,saverfs=saverfs

    ; Enviromental parameters
    prt = keyword_set(QUIET)
    if not(keyword_set(miter)) then miter=50  ; Maximum number of iterations
    if not(keyword_set(ilambda)) then ilambda=10D0;0.1d0 ; Assumes close to global minimum
    if not(keyword_set(toplim)) then toplim=1D-12 ; Very low, It may stop by MITER
    if not(keyword_set(C_LIMIT)) then C_LIMIT=1D-12 ; Very low, It may stop by MITER
    if not(keyword_set(llimit)) then llimit=[1e-32,1e+32]; Very low, It may stop by MITER
    ; if not(keyword_set(slight)) then print,'No stray light';slight=0 ; No stray light
    ; if not(keyword_set(filter)) then print,'No inst. filter profile'
    if not(prt) then print,'triplet='+string(keyword_set(triplet))
    if not(prt) then print,'Iterations='+string(miter)
    if not(keyword_set(noise)) then noise=1D-3 ; Default noise level
    if not(keyword_set(n_comp)) then n_comp=1D0  ; Maximum number of iterations
    IF KEYWORD_SET(NLTE) THEN NPAR = 15 ELSE NPAR = 11

    NTERMS = N_ELEMENTS(P_I)*1D0
    NFIXSS = TOTAL(FIX)
    FLAMBDA = ILAMBDA   ; Initial fudge parameter
    DIAG = INDGEN(NTERMS)*(NTERMS+1) ;SUBSCRIPTS OF DIAGONAL ELEMENTS
    FIXED = FIX*1D0
    NAXIS = N_ELEMENTS(AXIS)
    VALID = WHERE(STOKESPROF ne 0.,HAY) ; If Y = 0 then the samples are not valid (e.g. blends)
    NFREE = (N_ELEMENTS(VALID)-NFIXSS)*1D0  ; Number of degrees of freedom
    NPOINTS = (SIZE(STOKESPROF))(1)
    ITER_INFO = {lmb:dblarr(miter+1),iter:0,citer:0,conv_crit:dblarr(miter+1),$
                Params_stored:dblarr(nterms,miter+1),CHISQR:dblarr(miter+1)}

    if NFREE le 0 then begin
        PRINT,'Not enough points'
        STOP
        RETURN
    endif
    P_M = DBLARR(NTERMS) ;New model aprameters

    goodc = 0

    ;TO BE IMPLEMENTED
    ;   ;INITIALIZATION
    ;   p_i(6)=atan(max(STOKESPROF(*,1))/max(STOKESPROF(*,2)))/2.*180/!dpi
    ;   signo=total(STOKESPROF(0:NPOINTS/2,3))-total(STOKESPROF(NPOINTS/2+1:*,3))
    ;   if signo lt 0 then p_i(5)=100. else p_i(5)=10.

    ;INITIALIZING WEIGHTS
    WEIGHTS_INIT,NAXIS,STOKESPROF,W,SIG,WEIGHT=WEIGHT,SIGMA=SIGMA,NOISE=NOISE

    CLANDA = 0 ;SOME DEFINITIONS TO CONTROL THE FUDGE PARAMETER
    ITER = 0  ;DEFINING FIRST ITER
    LANDA_STORE = DBLARR(MITER+1) ;Array to store lambda variation with iter

    ;POLARIZATION THRESHOLD (If Stokes q, u, and v are below it, set B, GAMMA, and AZI to cero
    IF KEYWORD_SET(pol) then begin
        INDEX = INDGEN(N_COMP)
        LIMP = POL
        IF (MAX(ABS(STOKESPROF[*,1])) lt LIMP) AND (MAX(ABS(STOKESPROF[*,2])) lt LIMP) $
          AND (MAX(ABS(STOKESPROF[*,3])) lt LIMP) THEN BEGIN
            P_I[1] = 0D0 & P_I[5] = 0D0 & P_I[6] = 0D0
            FIXED[NPAR*INDEX+1] = 0D0 & FIXED[NPAR*INDEX+5] = 0D0 & FIXED[NPAR*INDEX+6] = 0D0
        ENDIF
    ENDIF

    ;Compute the ME RFs and the initial synthetic Stokes profiles for INIT MODEL
    ME_DER,P_I,WLI,AXIS,YFIT,PDER,TRIPLET=triplet,SLIGHT=slight,FILTER=filter,$
      AC_RATIO=ac_ratio,N_COMP=n_comp,NUMERICAL=numerical,IPBS=ipbs,crosst=crosst,nlte=nlte,saverfs=saverfs

    ;The derivatives with respect to fixed parameters are set to cero
    fxx = where(fixed eq 1)
    for I=0,nterms-1 DO PDER(*,I,*)=PDER(*,I,*)*FIXED(I)

    ; Take into account the local stray light in the merit function as in Asensio Ramos and Manso Sainz, 2010.
    If keyword_set(mlocal) then begin ;only valid out of NLTE
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

    ;Store old chisqr value
    OCHISQR=CHISQR

    P_M = P_I

    ;BEGIN THE ITERATION LOOP FOR THE LM OPTIMIZATION ALGORITHM
    REPEAT BEGIN

        iter_info.lmb[ITER] = FLAMBDA
        iter_info.CHISQR[ITER] = CHISQR
        ;IF (FLAMBDA GT 1d25) OR (FLAMBDA LT 1d-25) THEN CLANDA=1 ;Cond to Flambda !!!!!!!!!!
        COVAR=ALPHA
        COVAR(DIAG)=COVAR(DIAG)*(1D0+FLAMBDA)
        BETAD=BETA

        ;INVERT COVARIANCE MATRIX TO FIND NEW MODEL PARAMETERS.

        MIL_SVD,COVAR,BETAD,DELTA,LOSW,use_svd_cordic = use_svd_cordic

        ;CHECK MAX VARIATION OF DELTA AND SET NEW MODEL P_M

        for i=0,NPAR*n_comp-1 do IF VLIMITS(i).SET EQ 1 THEN $
        DELTA(i) = VLIMITS(i).LIMITS(0) > DELTA(i) < VLIMITS(i).LIMITS(1)


        P_M(fxx) = P_I(fxx) - DELTA(fxx) ;NEW PARAMS
    
        ;CHECK PARAMETERS
        CHECK_PARAM,P_M,NPAR,PLIMITS=plimits,n_comp=n_comp

        ;EVALUATE FUNCTION
        mil_sinrf,p_m,wli,AXIS,yfit,triplet=triplet,slight=slight,filter=filter,$
          AC_RATIO=ac_ratio,n_comp=n_comp,ipbs=ipbs,crosst=crosst,nlte=nlte

        ;new chisqr
        CHISQR = TOTAL(TOTAL((YFIT-STOKESPROF)^2d0*W,1,/double)/SIG^2d0,/double)/NFREE

        ;****CONVERGENCE CONDITION *****
        IF CHISQR-OCHISQR lt 0. then begin ;;FIT GOT BETTER
        ;****CONVERGENCE CONDITION *****

            if (ABS((OCHISQR-CHISQR)*100./CHISQR) LT TOPLIM) THEN CLANDA = 1
            if (CHISQR lt C_LIMIT) THEN CLANDA = 1
               ;Stoping criteria

            ;FIT GOT BETTER SO DECREASE FLAMBDA BY FACTOR OF 10
            flambda=0.1*flambda

            ;store new model parameters
            P_I(FXX) = P_M(FXX)
            ;store parameters for studing the evolution of them (opt mode)
            iter_info.Params_stored[*,ITER] = P_I
            iter_info.conv_crit[ITER] = 1.
            goodc = goodc + 1

            ;info
            IF not(prt) THEN PRINT,'ITERATION =',ITER,' , CHISQR =',CHISQR,$
              '  CONVERGE     - LAMBDA= ',flambda,(OCHISQR-CHISQR)/CHISQR

            ;compute new RFs and Stokes profiles
            me_der,p_i,wli,AXIS,yfit,pder,triplet=triplet,slight=slight,filter=filter,$
              AC_RATIO=ac_ratio,n_comp=n_comp,numerical=numerical,ipbs=ipbs,crosst=crosst,nlte=nlte,saverfs=saverfs

            for I=0,nterms-1 DO PDER(*,I,*)=PDER(*,I,*)*FIXED(I)

            ;Local stray light option (modification of chi2)
            If keyword_set(mlocal) then begin
              drho[0] = 2d0*(1d0 - p_i(10))/(mlocal*1d0)*sig[0] - 2*rho[0]
              term1 = 1d0 + (1d0 - p_i(10))^2d0/(mlocal*1d0) ;1-(1-alpha)^2/M
              sig = sigo*term1 - 2d0*(1d0 - p_i(10))*rho
            endif

            COVARM,W,SIG,YFIT,STOKESPROF,PDER,NTERMS,NFREE,BETA,ALPHA,OCHISQR,drho=dhro

        ENDIF ELSE BEGIN ;ASSUMES FIT GOT WORSE

            if (flambda GT llimit[1]) THEN CLANDA = 1
            if (flambda LT llimit[0]) THEN CLANDA = 1

            ;store parameters for studing the evolution of them (opt mode)
            iter_info.Params_stored[*,ITER] = P_I
            iter_info.conv_crit[ITER] = 0.

            ;increase flambda by a factor 10
            flambda=10.*flambda

            ;info
            IF not(prt) THEN PRINT,'ITERATION =',ITER,' , CHISQR =',OCHISQR,$
              '  NOT CONVERGE - LAMBDA= ',flambda,(OCHISQR-CHISQR)/CHISQR
        ENDELSE

        ITER=ITER+1

    ENDREP UNTIL (ITER GT MITER) OR (CLANDA EQ 1)

    iter_info.citer = goodc
    iter_info.iter = iter

    ;END OF MAIN LOOP

    COVAR=1./ALPHA
    ERR = SQRT(COVAR(DIAG)*OCHISQR*4*NPOINTS/NFREE) ; Almeida

    ;keep last merit function value
    CHISQF=OCHISQR

    ;Stokes I, Q, U, and V merit function values.
    getshi=fltarr(4)
    for ii=0,3 do getshi(ii)=TOTAL(TOTAL((YFIT(*,ii)-STOKESPROF(*,ii))^2d0*W(*,ii),1,/double)/SIG(ii)^2d0,/double)/NFREE

END
