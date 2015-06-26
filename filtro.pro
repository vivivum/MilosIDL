;+
; NAME: 
;	FILTRO
;
; AUTHOR: 
;	D. Orozco Suarez  	
;						National Astronomical Observatory of Japan,Ê
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: It calculates a Gaussian of a given full width at half maximum (FWHM, in mA) 
;          at given wavelengths to take the spectral smearing into account through convolution
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: filtro,FWHM,AXIS
;
; INPUTS:
;           FWHM: A scalar with the spectrograph PSF's
;           AXIS: Wavelength axis in Angstroms. Just one for all the lines in the sample
;                It should be a double precision array
;
; KEYWORD PARAMETERS:
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

FUNCTION filtro,FWHM,AXIS

lambda=AXIS(n_elements(AXIS)/2) ; central wavelength

; Conversion from FWHM to Gaussian sigma (1./(2*sqrt(2*alog2)))
sigma=FWHM*0.42466090d0/1000d0 ; in Angstroms

Term=((AXIS-lambda)/sigma)^2d0/2d0

Ande=where(Term LT 1.e30,loai)

IF loai GT 0 THEN filtro=exp(-Term[Ande]) ELSE filtro=exp(-Term)
filtro=filtro/total(filtro)  ; normalization in area

return,filtro

end

