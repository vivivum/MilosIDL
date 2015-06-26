;+
; NAME: 
;		READP
; AUTHOR: 
;	D. Orozco Suárez  	
;						National Astronomical Observatory of Japan, 
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
; 
; PURPOSE:
;    Read a profile archive with five columns
;           Wavelength,I,Q,U,V
;
; CALLING SEQUENCE:
;    readp,'archive', wlaxis ,I [,Q,U,V]
;
; ARGUMENTS
;
;    archive    The profile archive (IN)
;    Wlaxis     The wavelength axis     (OUT)
;    I,Q,U,V    If arguments Q,U,V are not present    (OUT)
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

pro readp,arc,l,i,q,u,v

get_lun,unit

openr,unit,arc
lines=0
while not eof(unit) do begin
   readf,unit,l0,i0,q0,u0,v0
   lines=lines+1
endwhile
free_lun,unit

get_lun,unit

datos=dblarr(lines,5)
openr,unit,arc
for j=0,lines-1 do begin
    readf,unit,l0,i0,q0,u0,v0;,FORMAT='(5E12.4)'
    datos(j,*)=[l0,i0,q0,u0,v0]
endfor
free_lun,unit

if n_params() eq 6 then begin

l=datos(*,0)
i=datos(*,1)
q=datos(*,2)
u=datos(*,3)
v=datos(*,4)

endif else begin

l=datos(*,0)
i=datos(*,1:4)

endelse

end

