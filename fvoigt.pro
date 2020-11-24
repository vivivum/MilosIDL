;+
; NAME:
;	FVOIGT
;
; AUTHOR:
;	D. Orozco Suarez
;						National Astronomical Observatory of Japan,�
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
;	J.C. del Toro Iniesta
;						Instituto de Astrofisica de Andalucia (CSIC)
;						Apdo de Correos 3004, 18080 Granada, SPAIN
;						jti@iaa.es
;
; PURPOSE: It calculates the Voigt and Faraday-Voigt functions for a given
;          damping parameter and at given wavelengths
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: fvoigt,DAMP,VV,H,F
;
; INPUTS:
;           DAMP: A scalar with the damping parameter
;           VV: Wavelength axis usually in Doppler units. Just one for
;               all the lines in the sample. It should be a double precision array
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;           H: Voigt function
;           F: Faraday-Voigt function
;
; COMMON BLOCKS:
;
; NOTES:
;           A rational approximation to the complex error function is used
;           after Hui, Armstrong, and Wray(1978, JQSRT 19, 509). H and F are
;           the real and imaginary parts of such function, respectively.
;           The procedure is inspired on that in SIR (Ruiz Cobo & del Toro
;           Iniesta 1992, ApJ 398, 385). On its turn, that routine was taken
;           from modifications by A. Wittmann (1986) to modifications by S.K.
;           Solanki (1985) to an original FORTRAN routine written by J.W. Harvey
;           and A. Nordlund.
;
;           The result is exactly the same as in SIR when the damping is not
;           zero (a non realistic case in which the line would be purely Gaussian)
;
; CALLED ROUTINES
;
; MODIFICATION HISTORY:
; 	Added documentation: Jose Carlos del Toro Iniesta, June, 2008
;
;-

pro fvoigt,DAMP,VV,H,F

	a=[122.607931777104326d0, 214.382388694706425d0, 181.928533092181549d0,$
	   93.155580458138441d0, 30.180142196210589d0, 5.912626209773153d0,$
	   0.564189583562615d0]

	b=[122.60793177387535d0, 352.730625110963558d0, 457.334478783897737d0, $
   	   348.703917719495792d0, 170.354001821091472d0, 53.992906912940207d0, $
	   10.479857114260399d0,1.d0]

	z=dcomplex(damp,-abs(vv))

  Z=((((((A(6)*Z+A(5))*Z+A(4))*Z+A(3))*Z+A(2))*Z+A(1))*Z+A(0))/$
          (((((((Z+B(6))*Z+B(5))*Z+B(4))*Z+B(3))*Z+B(2))*Z+B(1))*Z+B(0))
          
	h=double(z)
	f=sign(vv)*imaginary(z)*0.5

end
