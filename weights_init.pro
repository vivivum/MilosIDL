;+
; NAME: 
;	WEIGHTS_INIT
;
; AUTHOR: 
;	D. Orozco Suarez  	
;						National Astronomical Observatory of Japan, 
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: It reformat the array of weights and evaluates sigma for the Stokes profiles (if not given)
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: WEIGHTS_INIT,NAXIS,PROFILES,SIGMA=sigma,WEIGHTS=weights,W,SIG,NOISE=noise
;
; INPUTS:
;   NAXIS: Number of wavelength samples
;   STOKESPROF: Array with the input (if the /INVERSION keyword is set) Stokes profiles
;      Dimensions: (N_ELEMENTS(AXIS), 4)
;   WEIGHTS: Array with arbitrary weights for the Stokes profiles
;      Dimensions: (4) (for individual Stokes profiles)
;      Dimensions: (N_ELEMENTS(AXIS), 4) (as many as the profiles)
;   SIGMA: Array with noise values for the profiles 
;      Dimensions: (4)
;   NOISE: Scalar with the noise values for the Stokes profiles (Ignored when SIGMA is set)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   W: Array with arbitrary weights for the Stokes profiles
;      Dimensions: (N_ELEMENTS(AXIS), 4) (as many as the profiles)
;   SIG: Array with noise values for the profiles 
;
; COMMON BLOCKS: 
;
; NOTES:
;     In principle, it would be posible to calculate the sigma values automatically. However, 
;        this function is still under development (weights.pro).
;     
; CALLED ROUTINES
;    WEIGHTS: (commented)
;
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Suarez, Feb, 2009
; 	Bug fixed (axis -> naxis): DOS, March, 2010
;
;-

PRO WEIGHTS_INIT,NAXIS,STOKESPROF,SIGMA=SIGMA,WEIGHTS=WEIGHTS,W,SIG,NOISE=noise

  ;INITIALIZES THE SIGMA AND WEIGHTS FOR THE PROFILES
  ;THE WEIGHTS ARE USEFULL TO SELECT ONLY A FEW WAVELENGTH POSITIONS

  W=DBLARR(NAXIS,4)
  IF KEYWORD_SET(WEIGHTS) THEN BEGIN
    WT=WEIGHTS
    IF N_ELEMENTS(WT) GT 4 THEN BEGIN
      ;IF N_ELEMENTS(WT) EQ LLN THEN BEGIN
      ;  FOR I=0,3 DO W(*,i)=WT
      ;ENDIF ELSE BEGIN
        FOR I=0,3 DO W(*,i)=WT(*,i)
      ;ENDELSE
    ENDIF ELSE BEGIN
      FOR I=0,3 DO W(*,i)=WT(i)
    ENDELSE
  ENDIF ELSE BEGIN
    W(*,*)=1d0
  ENDELSE
  
  IF NOT(KEYWORD_SET(SIGMA)) THEN BEGIN
	;Indented to calculate sigma automatically from the stokes profiles (under development)
    ;WEIGHTS, STOKESPROF, SIG, W,noise=noise
    sig=DBLARR(4) & SIG[*]= noise
  ENDIF ELSE BEGIN
    SIG=DBLARR(4) & SIG[*]= SIGMA
  ENDELSE
  
END
