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
;   Added SVD cordic in 2015 by David Orozco (Experimental).
;-

PRO MIL_SVD,H,BETA,DELTA,W,use_svd_cordic = use_svd_cordic

  R=(SIZE(H))(1)
  EPSILON=1D-40  ;threshold to reject small eigenvalues (very small, this is usefull in other situations)
  TOP=1d0
  cual=where(finite(h) eq 0,hay)
  if hay gt 0 then h(cual)=0.

  if keyword_set(use_svd_cordic) then begin
    if use_svd_cordic eq 2 then begin
      ;compile once
      use_svd_cordic = 1
      INCLUDE=STREGEX(!MAKE_DLL.CC, '-I[^ ]+', /EXTRACT)
     ; Build the sharable library, using the CC keyword to specify gcc:
     MAKE_DLL, 'mil_svd_c', 'mil_svd_c','mil_svd_c', INPUT_DIRECTORY='./', /verbose,EXTRA_CFLAGS='-shared'
     spawn,'cp -v /Users/orozco/.idl/idl/compile_dir-118-idl_8_4-darwin-x86_64-m64-f64/mil_svd_c.so .',ff
   endif

  D1 =  h[0:R-2,0:R-2]*1.0D
  D2 =  beta(0:R-2)*1.0D
  D3 = DBLARR(10) & D3[*] = 1.0D
  I = 9
  DELTA = DBLARR(R)
;print,'i: ',i
;print,'d1: ',d1
;print,'d2: ',d2
;print,'d3: ',d3
ejecuta = CALL_EXTERNAL('mil_svd_c.so', 'mil_svd_c', d1,d2,I,D3, /CDECL)
;print,'i: ',i
;print,'d1: ',d1
;print,'d2: ',d2
;print,'d3: ',D3
;print,'d3 (reverse): ',D3[REVERSE(SORT(D3))]

DELTA(0:R-2) = D3

goto,noc

  SVDC,H,W,U,V,/double  ;SVD decomposition
  zz=dblarr(r,r)
  For j=0,r-1 do begin
    zz(j,j) = (abs(w(j)) GT EPSILON*TOP) ? 1D0/w(j) : 0D0
  endfor
  DELTA_SVD=(V##ZZ##TRANSPOSE(U))#BETA ;INVERSE MATRIX and DELTA

  print,'DELTA_SVD: ',DELTA_SVD[0:R-2]
  print,'DELTA_CORDIC: ',DELTA
noc:


endif else begin

  ;LA_SVD,H,W,U,V,/double  ;SVD decomposition
  SVDC,H,W,U,V,/double  ;SVD decomposition
  zz=dblarr(r,r)
  For j=0,r-1 do begin
    zz(j,j) = (abs(w(j)) GT EPSILON*TOP) ? 1D0/w(j) : 0D0
  endfor
  DELTA=(V##ZZ##TRANSPOSE(U))#BETA ;INVERSE MATRIX and DELTA

endelse

return

;OTROS tests

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
