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
; PURPOSE: 
; 		Returns the type of X as a long integer, in the (0,11) range.
; 		  0 Undefined
; 		  1 Byte
;		  2 Integer
;		  3 Long integer
;		  4 Float
; 		  5 Double precision
;		  6 Complex number
; 	 	  7 String
;		  8 Structure
;		  9   Double complex
;		  10 Pointer
;		  11 Object reference
;
; CALLING SEQUENCE: value = type(argument)
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Su‡rez, 24 Feb, 2009
;
;-

Function Type, x

; ++++
    dum = size(x)
    return, dum(dum(0) + 1)
; ++++

end

