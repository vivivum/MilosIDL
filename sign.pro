;+
; NAME: 
;	SIGN
;
; AUTHOR: 
;	D. Orozco Su‡rez  	
;						National Astronomical Observatory of Japan,Ê
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: Return the sign of scalar
;
; CALLING SEQUENCE: value = sign(scalar)
;
; OUTPUTS: 1 if positive, -1 if negative
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Su‡rez, 24 Feb, 2009
;
;-

function sign,f
	return, 2*(f ge 0) - 1
end
