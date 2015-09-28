;+
; NAME:
;	ipbs
;
; AUTHOR:
;	D. Orozco Suarez
;
; PURPOSE:
;;
; CATEGORY:
;	Splitting
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
;
; COMMON BLOCKS:
;			QUANTIC,C_N: store the results in this common block
;
; EXAMPLE:
;
; NOTES:
;			*IMPORTANT*
;			MILOS DOES NOT WORK IF NOT INITIALIZED WITH CREATE_NC. SEE INIT_MILOS.
;
; ROUTINES CALLING COMMON BLOCK:
;
; CALLED ROUTINES:
;
; MODIFICATION HISTORY:
; 	First version created, D Orozco Suarez (DOS), 2015
;-

PRO ipbs,MF,LD

COMMON QUANTIC,C_N

a=type(c_n)
if a then c_n=0

;this is for helium only
WL = [ 10830.2501d0, 10830.3397d0, 10829.0911d0]
slo0=1 & llo0=0 & jlo0=1
sup0=1 & lup0=1 & jup0=1
slo1=1 & llo1=0 & jlo1=1
sup1=1 & lup1=1 & jup1=2
slo2=1 & llo2=0 & jlo2=1
sup2=1 & lup2=1 & jup2=0
DAT=[3.,slo0,llo0,jlo0,sup0,lup0,jup0,0.3333,$
     slo1,llo1,jlo1,sup1,lup1,jup1,0.556,slo2,llo2,jlo2,sup2,lup2,jup2,0.1111]

LINES=DAT(0)     ;1,2,3 number of lines

SELECT=INDGEN(LINES)*7+1
Sl=DAT(SELECT)
Ll=DAT(SELECT+1)
Jl=DAT(SELECT+2)
Su=DAT(SELECT+3)
Lu=DAT(SELECT+4)
Ju=DAT(SELECT+5)
If Lines gt 1 then Fos=Dat(SELECT+6)

;ASSIGN THE VARIABLES TO THE STRUCTURE
FOR I=0,LINES-1 DO BEGIN
    N_PI=2.*MIN([jl(I),ju(I)])+1
    N_SIG=jl(I)+ju(I)
    ;QUANTEN COMPUTES THE CUANTEM NUMBERS
    QUANTEN,SL(I),SU(I),LL(I),LU(I),JL(I),JU(I),$
          NUB,NUP,NUR,WEB,WEP,WER,GLO,GUP,GEF,/not_normalize
    C_N(I).N_PI=N_PI
    C_N(I).N_SIG=N_SIG
    C_N(I).GL=GLO
    C_N(I).GU=GUP
    C_N(I).GEFF=GEF
    If Lines gt 1 then C_N(I).FO=Fos(I)
    FOR J=0,C_N(I).N_SIG-1 DO NUB[J]=NUB[J]+ipbs_shift('N_SIG_B',J,I,MF)/(MF*wl(I)^2D0/LD)/4.6686411D-13*1e-2
    FOR J=0,C_N(I).N_PI-1 DO NUP[J]=NUP[J]+ipbs_shift('N_PI',J,I,MF)/(MF*wl(I)^2D0/LD)/4.6686411D-13*1e-2
    FOR J=0,C_N(I).N_SIG-1 DO NUR[J]=NUR[J]+ipbs_shift('N_SIG_R',J,I,MF)/(MF*wl(I)^2D0/LD)/4.6686411D-13*1e-2
    FOR J=0,C_N(I).N_SIG-1 DO WEB[J]=WEB[J]+ipbs_strength('N_SIG_B',J,I,MF)
    FOR J=0,C_N(I).N_PI-1 DO WEP[J]=WEP[J]+ipbs_strength('N_PI',J,I,MF)
    FOR J=0,C_N(I).N_SIG-1 DO WER[J]=WER[J]+ipbs_strength('N_SIG_R',J,I,MF)
    WEB = WEB/total(WEB)
    WEP = WEP/total(WEP)
    WER = WER/total(WER)
    C_N(I).NUB(0:N_SIG-1)=NUB
    C_N(I).NUP(0:N_PI-1)=NUP
    C_N(I).NUR(0:N_SIG-1)=NUR
    C_N(I).WEB(0:N_SIG-1)=WEB
    C_N(I).WEP(0:N_PI-1)=WEP
    C_N(I).WER(0:N_SIG-1)=WER

ENDFOR

END

PRO TEST_IPBS

common quantic,c_n

init_milos,'10830',wl

nub_m = fltarr(3,3,100)
nup_m = fltarr(3,3,100)
nur_m = fltarr(3,3,100)
nub_o = fltarr(3,3,100)
nup_o = fltarr(3,3,100)
nur_o = fltarr(3,3,100)
web_m = fltarr(3,3,100)
wep_m = fltarr(3,3,100)
wer_m = fltarr(3,3,100)
web_o = fltarr(3,3,100)
wep_o = fltarr(3,3,100)
wer_o = fltarr(3,3,100)

B = findgen(100)*20.
LD = 0.1

for i=0,99 do begin
  init_milos,'10830',wl;,/not_normalize
  for j=0,2 do begin
    nub_o(j,*,i) = c_n(j).nub*(B[i]*wl(j+1)^2D0/LD)*4.6686411D-13
    nup_o(j,*,i) = c_n(j).nup*(B[i]*wl(j+1)^2D0/LD)*4.6686411D-13
    nur_o(j,*,i) = c_n(j).nur*(B[i]*wl(j+1)^2D0/LD)*4.6686411D-13
  web_o(j,*,i) = c_n(j).web
  wep_o(j,*,i) = c_n(j).wep
  wer_o(j,*,i) = c_n(j).wer
  endfor
    ipbs,B[i],LD
    for j=0,2 do begin
      nub_m(j,*,i) = c_n(j).nub*(B[i]*wl(j+1)^2D0/LD)*4.6686411D-13
      nup_m(j,*,i) = c_n(j).nup*(B[i]*wl(j+1)^2D0/LD)*4.6686411D-13
      nur_m(j,*,i) = c_n(j).nur*(B[i]*wl(j+1)^2D0/LD)*4.6686411D-13
    web_m(j,*,i) = c_n(j).web
    wep_m(j,*,i) = c_n(j).wep
    wer_m(j,*,i) = c_n(j).wer
    endfor
endfor

colores
!p.multi=[0,2,2]
plot,b,nub_m(2,0,*),line=2,thick=2,yrange=[-2.5,2.5],/nodata
oplot,b,nub_m(2,0,*),line=2,thick=2,color=2
oplot,b,nub_o(2,0,*),color=2
oplot,b,nur_m(2,0,*),line=2,thick=2,color=3
oplot,b,nur_o(2,0,*) ,color=3
oplot,b,nup_m(2,0,*),line=2,thick=2,color=4
oplot,b,nup_o(2,0,*),color=4

plot,b,nub_m(1,0,*),line=2,thick=2,yrange=[-2.5,2.5],/nodata
oplot,b,nub_m(1,0,*),line=2,thick=2,color=2
oplot,b,nub_m(1,1,*),line=2,thick=2,color=2
oplot,b,nub_m(1,2,*),line=2,thick=2,color=2
oplot,b,nub_o(1,0,*),color=2
oplot,b,nub_o(1,1,*),color=2
oplot,b,nub_o(1,2,*),color=2

oplot,b,nur_m(1,0,*),line=2,thick=2,color=3
oplot,b,nur_m(1,1,*),line=2,thick=2,color=3
oplot,b,nur_m(1,2,*),line=2,thick=2,color=3
oplot,b,nur_o(1,0,*) ,color=3
oplot,b,nur_o(1,1,*) ,color=3
oplot,b,nur_o(1,2,*) ,color=3

oplot,b,nup_m(1,0,*),line=2,thick=2,color=4
oplot,b,nup_m(1,1,*),line=2,thick=2,color=4
oplot,b,nup_m(1,2,*),line=2,thick=2,color=4
oplot,b,nup_o(1,0,*),color=4
oplot,b,nup_o(1,1,*),color=4
oplot,b,nup_o(1,2,*),color=4

plot,b,nub_m(0,0,*),line=2,thick=2,yrange=[-2.5,2.5],/nodata
oplot,b,nub_m(0,0,*),line=2,thick=2,color=2
oplot,b,nub_m(0,1,*),line=2,thick=2,color=2
oplot,b,nub_o(0,0,*),color=2
oplot,b,nub_o(0,1,*),color=2

oplot,b,nur_m(0,0,*),line=2,thick=2,color=3
oplot,b,nur_m(0,1,*),line=2,thick=2,color=3
oplot,b,nur_o(0,0,*) ,color=3
oplot,b,nur_o(0,1,*) ,color=3

oplot,b,nup_m(0,0,*),line=2,thick=2,color=4
oplot,b,nup_m(0,1,*),line=2,thick=2,color=4
oplot,b,nup_m(0,2,*),line=2,thick=2,color=4
oplot,b,nup_o(0,0,*),color=4
oplot,b,nup_o(0,1,*),color=4
oplot,b,nup_o(0,2,*),color=4

STOP
!p.multi=[0,2,2]
plot,b,web_m(2,0,*),line=2,thick=2,yrange=[0.99,1.01],/nodata,ystyle=1
oplot,b,web_m(2,0,*),line=2,thick=2,color=2
oplot,b,web_o(2,0,*),color=2
oplot,b,wer_m(2,0,*),line=2,thick=2,color=3
oplot,b,wer_o(2,0,*) ,color=3
oplot,b,wep_m(2,0,*),line=2,thick=2,color=4
oplot,b,wep_o(2,0,*),color=4

plot,b,web_m(1,0,*),line=2,thick=2,yrange=[0.05,0.7],/nodata,ystyle=1
oplot,b,web_m(1,0,*),line=2,thick=2,color=2
oplot,b,web_m(1,1,*),line=2,thick=2,color=2
oplot,b,web_m(1,2,*),line=2,thick=2,color=2
oplot,b,web_o(1,0,*),color=2
oplot,b,web_o(1,1,*),color=2
oplot,b,web_o(1,2,*),color=2

oplot,b,wer_m(1,0,*),line=2,thick=2,color=3
oplot,b,wer_m(1,1,*),line=2,thick=2,color=3
oplot,b,wer_m(1,2,*),line=2,thick=2,color=3
oplot,b,wer_o(1,0,*) ,color=3
oplot,b,wer_o(1,1,*) ,color=3
oplot,b,wer_o(1,2,*) ,color=3

oplot,b,wep_m(1,0,*),line=2,thick=2,color=4
oplot,b,wep_m(1,1,*),line=2,thick=2,color=4
oplot,b,wep_m(1,2,*),line=2,thick=2,color=4
oplot,b,wep_o(1,0,*),color=4
oplot,b,wep_o(1,1,*),color=4
oplot,b,wep_o(1,2,*),color=4

plot,b,web_m(0,0,*),line=2,yrange=[-0.1,0.6],/nodata,ystyle=1
oplot,b,web_m(0,0,*),line=2,thick=2,color=2
oplot,b,web_m(0,1,*),line=2,thick=2,color=2
oplot,b,web_o(0,0,*),color=2
oplot,b,web_o(0,1,*),color=2

oplot,b,wer_m(0,0,*),line=2,thick=2,color=3
oplot,b,wer_m(0,1,*),line=2,thick=2,color=3
oplot,b,wer_o(0,0,*) ,color=3
oplot,b,wer_o(0,1,*) ,color=3

oplot,b,wep_m(0,0,*),line=2,thick=2,color=4
oplot,b,wep_m(0,1,*),line=2,thick=2,color=4
oplot,b,wep_m(0,2,*),line=2,thick=2,color=4
oplot,b,wep_o(0,0,*),color=4
oplot,b,wep_o(0,1,*),color=4
oplot,b,wep_o(0,2,*),color=4



stop
END
