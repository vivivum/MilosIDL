;+
; NAME: 
;		WRITEP
; AUTHOR: 
;	D. Orozco Suárez  	
;						National Astronomical Observatory of Japan, 
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
; 
; PURPOSE:
;    Write a profile archive with five columns
;           Wavelength,I,Q,U,V
;
; CALLING SEQUENCE:
;    Writep,'archive', wlaxis ,I [,Q,U,V]
;
; ARGUMENTS
;
;    archive    The profile archive (IN)
;    Wlaxis     The wavelength axis     (IN)
;    I,Q,U,V    If arguments Q,U,V are not present    (IN)
;               I is a matrix of dimension (nuber of points,stokes)
;               If they are present, each parametes is stored in 
;               each variable
;
; COMMENTS
;
;    Really intuitive procedure, (Attention)
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Suárez, June, 2008
;
;-

pro WRITEP,arc,l,i,q,u,v

lines=n_elements(l)

datos=dblarr(lines,5)

if n_params() eq 6 then begin

datos(*,0)=l
datos(*,1)=i
datos(*,2)=q
datos(*,3)=u
datos(*,4)=v

endif else begin

datos(*,0)=l
datos(*,1:4)=i

endelse

get_lun,unit

openw,unit,arc
for j=0,lines-1 do printf,unit,datos(j,*),FORMAT='(F16.4,4E12.4)'
free_lun,unit

end
