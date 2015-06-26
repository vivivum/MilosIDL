;+
; NAME: 
;	MIL_SVD
;
; AUTHOR: 
;	D. Orozco Suarez  	
;						National Astronomical Observatory of Japan, 
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: It calculates the inverse of the covariance matrix and 
;          change to be made in the model parameter
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: MIL_SVD,H,BETA,DELTA
;
; INPUTS:
;        H: Approximate-Hessian (covariance) matrix. Dimension (11,11)
;        BETA: Gradient. Dimension (11)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS: 
;        DELTA: change to be made in the model parameter. Dimension (11)
;
; COMMON BLOCKS: 
;
; NOTES:
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Suarez, Feb 2009
;
;-

PRO MIL_SVD,H,BETA,DELTA,w

  R=(SIZE(H))(1)
  EPSILON=1d-18  ;threshold to reject small eigenvalues (very small, this is usefull in other situations)
  TOP=1d0  
  cual=where(finite(h) eq 0,hay)
  if hay gt 0 then h(cual)=0.

  ;GOTO,ALTERNATIVA
  SVDC,H,W,U,V,/double  ;SVD decomposition
  zz=dblarr(r,r)
  For j=0,r-1 do begin
    zz(j,j) = (abs(w(j)) GT EPSILON*TOP) ? 1D0/w(j) : 0D0  
  endfor
  DELTA=(V##ZZ##TRANSPOSE(U))#BETA ;INVERSE MATRIX and DELTA
  
  RETURN
  ALTERNATIVA:

  SVDC,H,W,U,V,/double  ;SVD decomposition
  ;LA INVERSA SERIA V##(1/W)##Transpose(U)
  ;LA DIRECTA SERIA U##W##U
  ;LA NUEVA: INVERT(V##W^2##TRANSPOSE(V))##V##W##TRANSPOSE(V)
  WW=dblarr(r,r)
  For j=0,r-1 do WW(j,j) = W(j)  
  MM = V##WW^2.##TRANSPOSE(V)
  SVDC,MM,W,UU,VV,/double  ;SVD decomposition
  zz=dblarr(r,r)
  For j=0,r-1 do zz(j,j) = (abs(w(j)) GT EPSILON*TOP) ? 1D0/w(j) : 0D0  
  INVERSE_H = (VV##ZZ##TRANSPOSE(UU) )##V##WW##TRANSPOSE(U)
  
  DELTA=INVERSE_H#BETA ;INVERSE MATRIX and DELTA
  
END
