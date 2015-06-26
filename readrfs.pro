;+
; NAME: 
;		READRFS
; AUTHOR: 
;	D. Orozco Suárez  	
;						National Astronomical Observatory of Japan, 
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
; 
; PURPOSE:
;    Read response function from an archive
;
; CALLING SEQUENCE:
;    READRFS, 'archive', dim, Wlaxis, dat, rfsi, rfsq, rfsu, rfsv
;
; ARGUMENTS
;
;    archive    The profile archive (IN)
;    dim        The wavelength axis dimension (IN)
;    Wlaxis     The wavelength axis  (OUT)
;    dat        The Stokes vector (Dim,[I,Q,U,V])
;    rfsi       Response function of I to respect the parameters
;    rfsq       Response function of I to respect the parameters
;    rfsu        Response function of I to respect the parameters
;    rfsv       Response function of I to respect the parameters
;          Dimension - (dim,parameter)
;
; COMMENTS
;    Really intuitive procedure, (Attention)
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Suárez, June, 2008
;
;-

pro READRFS,arc,dim,l,datos,rfsi,rfsq,rfsu,rfsv

datos=dblarr(dim,5)
rfsi=dblarr(dim,11)
rfsq=dblarr(dim,11)
rfsu=dblarr(dim,11)
rfsv=dblarr(dim,11)

get_lun,unit

openr,unit,arc
for j=0,dim-1 do begin
    readf,unit,l0,i0,q0,u0,v0;,FORMAT='(5E12.4)'
    datos(j,*)=[l0,i0,q0,u0,v0]
endfor

for j=0,dim-1 do begin
    readf,unit,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11
    rfsi(j,*)=[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]
endfor
for j=0,dim-1 do begin
    readf,unit,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11
    rfsq(j,*)=[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]
endfor
for j=0,dim-1 do begin
    readf,unit,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11
    rfsu(j,*)=[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]
endfor
for j=0,dim-1 do begin
    readf,unit,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11
    rfsv(j,*)=[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]
endfor

free_lun,unit


l=datos(*,0)
datos=datos(*,1:4)

end
