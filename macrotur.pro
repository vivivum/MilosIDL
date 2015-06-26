;+
; NAME: 
;	MACROTUR
;
; AUTHOR: 
;	D. Orozco Suarez  	
;						National Astronomical Observatory of Japan,Ê
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
;	J.C. del Toro Iniesta  
;						Instituto de Astrofisica de Andalucia (CSIC)
;						Apdo de Correos 3004, 18080 Granada, SPAIN 
;						jti@iaa.es
;
; PURPOSE: It calculates a Gaussian of a given width (in km/s) at given wavelengths
;          to take macroturbulence into account through convolution
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: MACROTUR,WIDTH,AXIS,LAMBDA[,deriv=deriv]
;
; INPUTS:
;           WIDTH: A scalar with the macroturbulent velocity in km/s
;           AXIS: Wavelength axis in Angstroms. Just one for all the lines in the sample
;                It should be a double precision array
;           LAMBDA: Central wavelength of the line (double precision constant)
;
; KEYWORD PARAMETERS:
;           DERIV: Set this keyword to output the derivative with respect to
;                  the macroturbulent velocity
;
; OUTPUTS:
;
; COMMON BLOCKS: 
;
; NOTES:
;           The wavelength axis is assumed to have an odd number of samples
;
; CALLED ROUTINES
;
; MODIFICATION HISTORY:
; 	Added documentation: Jose Carlos del Toro Iniesta, June, 2008
;
;-


FUNCTION MACROTUR,WIDTH,AXIS,LAMBDA,deriv=deriv


vlight=double(2.99792458E+5)   ; speed of light (cm/s)

centro=AXIS(n_elements(AXIS)/2) ; central wavelength

ILd=LAMBDA*WIDTH/vlight ; Gaussian sigma in Angstroms

Term=((AXIS-centro)/ILd)^2d0/2d0 ; exponent

Ande=where(Term LT 1.e30,loai) ; Nan otherwise

IF loai GT 0 THEN mtb=exp(-Term[Ande]) ELSE mtb=exp(-Term)

cte=total(mtb)
mtb=mtb/cte ; normalization in area

; The derivative case

IF keyword_set(Deriv) THEN BEGIN

mtb2=mtb/WIDTH*(((AXIS-centro)/ILd)^2d0-1d0)

return,mtb2

ENDIF ELSE return, mtb

end

