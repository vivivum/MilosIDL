;+
; NAME:
;	COVARM
;
; AUTHOR:
;	D. Orozco Suarez
;						National Astronomical Observatory of Japan,�
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: It calculates the gradient and covariance matrixes, and the merit function
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: COVARM,W,SIG,YFIT,STOKESPROF,PDER,NTERMS,NFREE,BETA,ALPHA,CHISQR
;
; INPUTS:
;			W: Weights for Stokes profiles. Dimensions: (N_ELEMENTS(AXIS), 4 )
;           SIG: Standar deviation of Stokes profiles. Dimensions: (4)
;           STOKESPROF: Array with the input Stokes profiles. Dimensions: (N_ELEMENTS(AXIS), 4)
;           YFIT: Array with the synthetic Stokes profiles. Dimensions: (N_ELEMENTS(AXIS), 4)
;           PDER: Response functions of Stokes profiles. Dimensions:  (N_ELEMENTS(AXIS), 4, 11)
;           NTERMS: Number of parameters to fit. Dimensions:  Scalar
;           NFREE: Degrees of freedom. Dimensions:  Scalar
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;           BETA: Gradient of merit function. Dimensions: (NTERMS)
;           ALPHA: approximate-Hessian (covariance) matrix. Dimensions:  (NTERMS, NTERMS)
;           CHISQF: Scalar with the chi^2 merit function of the inversion. Dimensions: scalar
;
; COMMON BLOCKS:
;
; NOTES:
;
; CALLED ROUTINES:
;
; CALLING ROUTINES:  LM_MILS
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Suarez Feb, 2008
;
;-

PRO COVARM,W,SIG,YFIT,STOKESPROF,PDER,NTERMS,NFREE,BETA,ALPHA,CHISQR,drho=dhro

  BT = DBLARR(NTERMS,4)
  AP = DBLARR(NTERMS,NTERMS,4)
  UNOS = DBLARR(NTERMS)+1D

if not(keyword_set(dhro)) then begin

  for I=0,3 do begin
    BT(*,I) = (W(*,I)*(YFIT(*,I)-STOKESPROF(*,I))#PDER(*,*,I))/SIG(I)^2
    AP(*,*,I)=(TRANSPOSE(PDER(*,*,I))#(W(*,I)$
      #UNOS*PDER(*,*,I)))/SIG(I)^2
  endfor

endif else begin

  for I=0,3 do begin
    BT(*,I) = (W(*,I)*(YFIT(*,I)-STOKESPROF(*,I))#PDER(*,*,I))/SIG(I)^2d0
    BT(10,I) = BT(10,I) +  1D/2D* total( W(*,I) * (YFIT(*,I)-STOKESPROF(*,I))^2d0 * drho(I) / (SIG(I)^2d0)^2d0 )
    AP(*,*,I)=(TRANSPOSE(PDER(*,*,I))#(W(*,I)#UNOS*PDER(*,*,I)))/SIG(I)^2d0
  endfor

  endelse

  BETA = TOTAL(BT,2,/double)
  ALPHA = TOTAL(AP,3,/double)
  CHISQR = TOTAL(TOTAL((YFIT-STOKESPROF)^2*W,1,/double)/SIG^2,/double)/NFREE

END
