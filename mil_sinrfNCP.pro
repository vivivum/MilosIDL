pro mil_sinrfNCP,PARAM,L0,LMB,SPECTRA,triplete=triplete

common N_CUANTICOS,nub,nup,nur,web,wep,wer,glo,gup,mlo,mup,geff
common HELIO,AH

;Hay que poner n_cuanticos como una entrada

;parametros
;E0,MF,VL,LD,A,GM,AZI,B1,B2,MC
B0=param(7)
B1=param(8)
E0=param(0)
MF=param(1)
GM=param(5)
AZI=param(6)
VL=param(2)
LD=param(3)
A=param(4)
MC=param(9)

;A RADIANES
CC=!dpi/180d0
AZI=AZI*CC
GM=GM*CC

NUML=n_elements(LMB)
SPECTRA=DBLARR(NUML,4)
Vlight=double(2.99792458E+5)   ;light speed (cm/s)
RR=5.641895836E-1

;doppler velocity (for voigt function)
u=(LMB-L0)/LD
;frecuency shift for v line of sight
ulos=(VL*L0)/(vlight*LD)

IF NOT KEYWORD_SET(TRIPLETE) THEN BEGIN
	;CASO MULTIPLETE
	;*************************************************
	nubB=(MF*L0^2./ld)*(4.6686411e-13*nub) ;(Spliting)
	nupB=(MF*L0^2./ld)*(4.6686411e-13*nup) ;(Spliting)
	nurB=(MF*L0^2./ld)*(4.6686411e-13*nur) ;(Spliting)
	;*************************************************
	FI_P=DBLARR(NUML) & FI_B=FI_P & FI_R=FI_P
	     SHI_P=FI_P & SHI_B=FI_P & SHI_R=FI_P

    ;CENTRAL COMPONENT
    FOR i=0, n_elements(wep)-1 do BEGIN
    	UU=u-ulos-nupB(i)
    	fvoigt, A, UU, H, F ;P
    	FI_P=FI_P+wep(i)*H & SHI_P=SHI_P+wep(i)*F
    ENDFOR
    FI_P=FI_P*RR & SHI_P=SHI_P*RR

    ;BLUE COMPONENT
    FOR i=0,n_elements(web)-1 do begin
    	UU=u-ulos-nubB(i)
    	fvoigt,A,UU, H, F ;B
    	FI_B=FI_B+web(i)*H & SHI_B=SHI_B+web(i)*F
    ENDFOR
    FI_B=FI_B*RR & SHI_B=SHI_B*RR
    ;RED COMPONENT
    FOR i=0,n_elements(wer)-1 do begin
    	UU=u-ulos-nurB(i)
    	fvoigt,A,UU, H, F ;R
    	FI_R=FI_R+wer(i)*H & SHI_R=SHI_R+wer(i)*F
    ENDFOR
    FI_R=FI_R*RR & SHI_R=SHI_R*RR

ENDIF ELSE BEGIN

	;CASO TRIPLETE
	SHIF=(MF*L0^2./ld)*(4.6686411e-13*geff)
    ;CENTRAL COMPONENT
    fvoigt, A, U-ULOS, H, F
    FI_P=H*RR & SHI_P=F*RR
    dH_u=2.*A*F-2.*(U-ULOS)*H
    dF_u=2.*(RR-A*H-(U-ULOS)*F)
    ;D1_C=dH_u*RR*(-L0/LD)
    ;D2_C=dF_u*RR*(-L0/LD) ;U
    D1_C=0.;
    D2_C=0.;

    ;BLUE COMPONENT
    fvoigt, A, U-ULOS+SHIF, H, F
    FI_B=H*RR & SHI_B=F*RR
    dH_u=2.*A*F-2.*(U-ULOS)*H
    dF_u=2.*(RR-A*H-(U-ULOS)*F)
    ;D1_B=dH_u*RR*(-L0/LD)
    ;D2_B=dF_u*RR*(-L0/LD) ;U
    D1_B=dH_u*RR*(SHIF/MF)
    D2_B=dF_u*RR*(SHIF/MF)

    ;RED COMPONENT
    fvoigt, A, U-ULOS-SHIF, H, F
    FI_R=H*RR & SHI_R=F*RR
    dH_u=2.*A*F-2.*(U-ULOS)*H
    dF_u=2.*(RR-A*H-(U-ULOS)*F)
    ;D1_R=dH_u*RR*(-L0/LD)
    ;D2_R=dF_u*RR*(-L0/LD) ;U
    D1_R=dH_u*RR*(-SHIF/MF)
    D2_R=dF_u*RR*(-SHIF/MF) ;B


ENDELSE

;DISPERSION PROFILES  (NUML)
ETAI=1.+E0*[FI_P*sin(GM)^2.+(FI_B+FI_R)*(1.+cos(GM)^2.)/2.]/2.
ETAQ=E0*[FI_P-(FI_B+FI_R)/2.]*sin(GM)^2.*cos(2.*AZI)/2.
ETAU=E0*(FI_P-(FI_B+FI_R)/2.)*sin(GM)^2.*sin(2.*AZI)/2.
ETAV=E0*(FI_R-FI_B)*cos(GM)/2.
RHOQ=E0*(SHI_P-(SHI_B+SHI_R)/2.)*sin(GM)^2.*cos(2.*AZI)/2.
RHOU=E0*(SHI_P-(SHI_B+SHI_R)/2.)*sin(GM)^2.*sin(2.*AZI)/2.
RHOV=E0*(SHI_R-SHI_B)*cos(GM)/2.

DELTA=ETAI^2.*(ETAI^2.-ETAQ^2.-ETAU^2.-ETAV^2.+RHOQ^2.+RHOU^2.+RHOV^2.)-$
      (ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)^2.
DELTAI=1./DELTA

SPECTRA(*,0)=B0+B1*DELTAI*ETAI*(ETAI^2.+RHOQ^2.+RHOU^2.+RHOV^2.)*AH
SPECTRA(*,1)=(-B1*DELTAI*(ETAI^2.*ETAQ+ETAI*(ETAV*RHOU-ETAU*RHOV)+$
  RHOQ*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*AH
SPECTRA(*,2)=(-B1*DELTAI*(ETAI^2.*ETAU+ETAI*(ETAQ*RHOV-ETAV*RHOQ)+$
  RHOU*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*AH
SPECTRA(*,3)=(-B1*DELTAI*(ETAI^2.*ETAV+ETAI*(ETAU*RHOQ-ETAQ*RHOU)+$
  RHOV*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*AH

;macroturbulencia
;we make gauss profile
G = FGAUSS(MC,LMB,L0)
FFTG=FFT(G, 1, /double)
;CONVOLUTION

;FOR I=0,3 DO BEGIN
;    IF abs(mean(SPECTRA(*,I))) LE 1.e-08 THEN GOTO, jump
;    SPECTRA(*,I)=FFT(FFT(SPECTRA(*,I), 1, /double)*FFTG, -1, /double)
;    SPECTRA(*,I)=shift(SPECTRA(*,I),-NUML/2)
;    Jump:
;ENDFOR


;USING matrizes
UNITY = [[1., 0., 0., 0.], $
         [0., 1., 0., 0.], $
         [0., 0., 1., 0.], $
         [0., 0., 0., 1.]]
UNIT = transpose([1., 0, 0, 0])

    MFHI = [[[ETAI],[ETAQ],[ETAU],[ETAV]],$
          [[ETAQ],[ETAI],[RHOV],[-RHOU]],$
          [[ETAU],[-RHOV],[ETAI],[RHOQ]],$
          [[ETAV],[RHOU],[-RHOQ],[ETAI]]]

sPECT=DBLARR(100,4)
FOR I=0,99 DO SPECT(I,*)=B1*AH*(INVERT(REFORM(MFHI(I,*,*)))##UNIT)+B0*UNIT

pert=DBLARR(100,4)

DEI=E0*(D1_C*sin(GM)^2.+(D1_B+D1_R)*(1.+cos(GM)^2.)/2.)/2.
DEQ=E0*(D1_C-(D1_B+D1_R)/2.)*sin(GM)^2.*cos(2.*AZI)/2.
DEU=E0*(D1_C-(D1_B+D1_R)/2.)*sin(GM)^2.*sin(2.*AZI)/2.
DEV=E0*(D1_R-D1_B)*cos(GM)/2.
DRQ=E0*(D2_C-(D2_B+D2_R)/2.)*sin(GM)^2.*cos(2.*AZI)/2.
DRU=E0*(D2_C-(D2_B+D2_R)/2.)*sin(GM)^2.*sin(2.*AZI)/2.
DRV=E0*(D2_R-D2_B)*cos(GM)/2.

;NUMERICA

;    DMFHI = [[[DEI],[DEQ],[DEU],[DEV]],$
;          [[DEQ],[DEI],[DRV],[-DRU]],$
;          [[DEU],[-DRV],[DEI],[DRQ]],$
;          [[DEV],[DRU],[-DRQ],[DEI]]]


DMFHII=MFHI
FOR I=0,99 DO DMFHII(I,*,*)=INVERT(REFORM(MFHI(I,*,*)))
DMFHIII=DMFHII

II=[0,-1,1]
JJ=[0,-1,1]

FOR I=0,3 DO BEGIN
FOR J=0,3 DO BEGIN
DMFHIII(*,I,J)=DERIV(U,DMFHII(*,I,J))*SHIF/MF
;DMFHIII(*,I,J)=DERIV(U,DMFHII(*,I,J))*(1./LD)
ENDFOR
ENDFOR


WA=0.2
BB=500.

;FOR I=0,99 DO PERT(I,*)=-B1*L0/LD*WA/VLIGHT*iNVERT(REFORM(MFHI(I,*,*)))##$
;(REFORM(DMFHIII(I,*,*)))##UNIT
FOR I=0,99 DO PERT(I,*)=-B1*BB*iNVERT(REFORM(MFHI(I,*,*)))##$
(REFORM(DMFHIII(I,*,*)))##UNIT

SPECTOT=SPECTRA+PERT

SET_PLOT,'PS'
DEVICE,FILENAME='MILNE.PS'
!P.MULTI=[0,2,2]
PLOT,SPECTOT(*,0)
OPLOT,SPECTRA(*,0),LINE=2,THICK=2
PLOT,SPECTOT(*,1)
OPLOT,SPECTRA(*,1),LINE=2,THICK=2
PLOT,SPECTOT(*,2)
OPLOT,SPECTRA(*,2),LINE=2,THICK=2
PLOT,SPECTOT(*,3)
OPLOT,SPECTRA(*,3),LINE=2,THICK=2
DEVICE,/CLOSE
stop

end
