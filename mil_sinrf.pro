;+
; NAME:
;	MIL_SINRF
;
; AUTHOR:
;	D. Orozco Su�rez
;						National Astronomical Observatory of Japan,�
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
;	J.C. del Toro Iniesta
;						Instituto de Astrof�sica de Andaluc�a (CSIC)
;						Apdo de Correos 3004, 18080 Granada, SPAIN
;						jti@iaa.es
;
; PURPOSE: Synthesis of Stokes profiles without calculating response functions
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: MIL_SINRF, PARAM, WL, LMB, SPECTRA, triplet=triplet,$
;                  mu=mu,slight=slight,filter=filter
;
; INPUTS:
;           PARAM: Array with the 11 model-atmosphere parameters
;               eta0 = line-to-continuum absorption coefficient ratio
;               B = magnetic field strength       [Gauss]
;               vlos = line-of-sight velocity     [km/s]
;               dopp = Doppler width              [Angstroms]
;               aa = damping parameter
;               gm = magnetic field inclination   [deg]
;               az = magnetic field azimuth       [deg]
;               S0 = source function constant
;               S1 = source function gradient
;               mac = macroturbulent velocity     [km/s]
;               alpha = filling factor of the magnetic component [0->1]
;           WL: A scalar or array with the central wavelength(s) of the line(s) in Angstroms
;                Dimensions: (number of lines, central wavelengths)
;                Examples: [1, 6302.5]; [2, 6301.5, 6302.5]
;           LMB: Wavelength axis in Angstroms. Just one for all the lines in the sample
;                It should be a double precision array
;           SLIGHT: Array with the stray-light profile
;                       Dimensions: (N_ELEMENTS(LMB), 4)
;           MU: Scalar containing the cosine of the heliocentric angle
;           FILTER: Width (in mA) of a Gaussian filter profile or array with the instrument
;                   profile samples. Dimensions: N_ELEMENTS(LMB)
;
; KEYWORD PARAMETERS:
;           TRIPLET: Set this keyword to do the calculations for the effective triplet
;
; OUTPUTS:
;           SPECTRA: Array with the input (if the /INVERSION keyword is set) Stokes profiles
;                       Dimensions: (N_ELEMENTS(AXIS), 4)
;
; COMMON BLOCKS:
;           QUANTIC: Contains all the atomic data relevant to the spectral lines
;
; NOTES:
;           The stray-light profile is usually assumed to be non polarized. Hence, it should
;           have the last three dimensions equal to zero. The program is however prepared
;           with the more general case of a polarized stray-light profile.
;
;           The filter profile must have the same spectral resolution than the Stokes profiles.
;
; CALLED ROUTINES
;           FVOIGT: Calculates the Voigt and Faraday-Voigt functions
;           MACROTUR: Evaluates a Gauss function
;           FILTRO: Outputs a Gaussian filter instrumental profile in case that only the
;                   width is provided
;
; MODIFICATION HISTORY:
; 	First beta version created, D Orozco Su�rez (DOS) and J.C. del Toro Iniesta (JTI), 2007
;   First beta documentation version, JTI, June 2008
;   Added keywork AC_RATIO. 21 Jan, 2010. DOS
;-

pro MIL_SINRF,PARAM,WL,LMB,SPECTRA,TRIPLET=triplet,MU=mu,SLIGHT=slight,$
	FILTER=filter,AC_RATIO=ac_ratio,N_COMP=n_comp

COMMON QUANTIC,C_N

lines=wl(0)

IF NOT(KEYWORD_SET(mu)) THEN mu=1d0
IF NOT(KEYWORD_SET(AC_RATIO)) THEN AC_RATIO=0

; Some initial settings

NUML=N_ELEMENTS(LMB)
SPECTRA=DBLARR(NUML,4)
SPECTRA_comp=DBLARR(NUML,4,n_comp)
vlight=2.99792458D+5   ; speed of light (cm/s)
RR=5.641895836D-1
CC=!dpi/180d0  ; conversion factor to radians

alpha = param(10)
fill_fractions = dblarr(n_comp)
if n_comp gt 1 then begin
	fill_fractions(1:*)=param(11*(indgen(n_comp-1)+1)+10)
	fill_fractions(0) = 1d0 - total(fill_fractions(1:*))
endif else fill_fractions = 1d0

for k=0,n_comp-1 do begin ;loop in components

;Model atmosphere parameters E0,MF,VL,LD,A,GM,AZI,B1,B2,MC,ALPHA

E00=param(11*k)
MF=param(11*k+1)
VL=param(11*k+2)
LD=param(11*k+3)
A=param(11*k+4)
GM=param(11*k+5)
AZI=param(11*k+6)
B0=param(11*k+7)
B1=param(11*k+8)
MC=param(11*k+9)

ETAI=1D0 & ETAQ=0D0 & ETAU=0D0 & ETAV=0D0  ; propagation matrix elements
RHOQ=0D0 & RHOU=0D0 & RHOV=0D0
VNULO=DBLARR(NUML)

AZI=AZI*CC
GM=GM*CC

sinis=sin(GM)^2d0
cosis=cos(GM)^2d0
cosi=cos(GM)
sina=sin(2D0*AZI)
cosa=cos(2D0*AZI)

; Loop in spectral lines

FOR IL=0,lines-1 DO BEGIN

	;ratio
	ratio = C_N(IL).Fo
	; Ratio 6301. / 6302.  = 2.92000
	; 1 / ratio = 0.342466 =  1. / 2.92000

	if (ratio ne 1) and (ac_ratio) then begin
		coeffa = [1.01702d-4, -9.33108d-3, 5.41421d-1, 3.11270d-2 ]
		coeffb = [2.26080d-1, 4.11138d0 ]
		if (E00 le 38) then begin

			E0 = coeffa(0)*E00^3d0 + coeffa(1)*E00^2d0 + $
	    	    coeffa(2)*E00 + coeffa(3)

		endif else if (E00 lt 300) then begin

			E0 = coeffb(0)*E00 + coeffb(1)

		endif

	endif else begin
    ; Line strength
	    E0=E00*ratio  ; C_N(IL).Fo is the relative line strength
	endelse

	;print,E00,ratio,E00*ratio,E0

    ; Wavelengths in units of the Doppler width
    u=(LMB-wl(IL+1))/LD

    ;frequency shift for v line of sight
    ulos=(VL*wl(IL+1))/(vlight*LD)

    IF NOT(KEYWORD_SET(triplet)) THEN BEGIN

        ; GENERAL MULTIPLET CASE
        ;***************************************************
        ; Zeeman splittings
        nubB=(MF*wl(IL+1)^2D0/ld)*(4.6686411D-13*C_N(IL).nub) ; sigma_b components
        nupB=(MF*wl(IL+1)^2D0/ld)*(4.6686411D-13*C_N(IL).nup) ; pi components
        nurB=(MF*wl(IL+1)^2D0/ld)*(4.6686411D-13*C_N(IL).nur) ; sigma_r components
        ;***************************************************

        ; Absorption and dispersion profiles initialization

        FI_P=VNULO & FI_B=VNULO & FI_R=VNULO
        SHI_P=VNULO & SHI_B=VNULO & SHI_R=VNULO

        ; ABSORPTION AND DISPERSION PROFILES

        ; PI COMPONENTS

        FOR i=0,C_N(IL).N_PI-1 DO BEGIN

            ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
            UU=u-ulos-nupB(i)

            fvoigt, A, UU, H, F
            FI_P=FI_P+C_N(IL).wep(i)*H
            SHI_P=SHI_P+C_N(IL).wep(i)*F

        ENDFOR
        SHI_P=2d0*SHI_P

        ; SIGMA BLUE COMPONENTS

        FOR i=0,C_N(IL).N_SIG-1 DO BEGIN

            ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
            UU=u-ulos-nubB(i)

            fvoigt, A, UU, H, F
            FI_B=FI_B+C_N(IL).web(i)*H
            SHI_B=SHI_B+C_N(IL).web(i)*F

        ENDFOR
        SHI_B=2d0*SHI_B

        ; SIGMA RED COMPONENTS

        FOR i=0,C_N(IL).N_SIG-1 DO BEGIN

            ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
            UU=u-ulos-nurB(i)

            fvoigt,A,UU, H, F   ;R
            FI_R=FI_R+C_N(IL).wer(i)*H
            SHI_R=SHI_R+C_N(IL).wer(i)*F

        ENDFOR
        SHI_R=2d0*SHI_R

    ENDIF ELSE BEGIN

        ; TRIPLET CASE

        SHIF=(MF*wl(IL+1)^2D0/LD)*(4.6686411D-13*C_N(IL).GEFF)

        ; PI COMPONENT

        fvoigt, A, U-ULOS, H, F
        FI_P=H & SHI_P=F
        SHI_P=2d0*SHI_P

        ; SIGMA BLUE COMPONENT

        fvoigt, A, U-ULOS+SHIF, H, F
        FI_B=H & SHI_B=F
        SHI_B=2d0*SHI_B

        ; SIGMA RED COMPONENT

        fvoigt, A, U-ULOS-SHIF, H, F
        FI_R=H & SHI_R=F
        SHI_R=2d0*SHI_R

    ENDELSE

    ; Adding the necessary geometry and line strength

    ETAI=ETAI+E0/2D0*(FI_P*sinis+(FI_B+FI_R)*(1D0+cosis)/2D0)
    ETAQ=ETAQ+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*cosa
    ETAU=ETAU+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*sina
    ETAV=ETAV+E0/2D0*(FI_R-FI_B)*cosi
    RHOQ=RHOQ+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*cosa
    RHOU=RHOU+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*sina
    RHOV=RHOV+E0/2D0*(SHI_R-SHI_B)*cosi

ENDFOR

;Etaq = etaq*(-1d0)
;Etau = etau*(-1d0)
;rhoq = rhoq*(-1d0)
;Etau = etau*(-1d0)

; THE STOKES PROFILES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

DELTA=ETAI^2D0*(ETAI^2D0-ETAQ^2D0-ETAU^2D0-ETAV^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)-$
      (ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)^2D0
DELTAI=1D0/DELTA

SPECTRA(*,0)=(B0+B1*DELTAI*ETAI*(ETAI^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)*MU)
SPECTRA(*,1)=(-B1*DELTAI*(ETAI^2D0*ETAQ+ETAI*(ETAV*RHOU-ETAU*RHOV)+$
  RHOQ*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*MU
SPECTRA(*,2)=(-B1*DELTAI*(ETAI^2D0*ETAU+ETAI*(ETAQ*RHOV-ETAV*RHOQ)+$
  RHOU*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*MU
SPECTRA(*,3)=(-B1*DELTAI*(ETAI^2D0*ETAV+ETAI*(ETAU*RHOQ-ETAQ*RHOU)+$
  RHOV*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*MU

; CONVOLUTIONS WITH THE MACROTURBULENT VELOCITY AND THE INSTRUMENTAL PROFILE

; Taking care of the wavelength centering. Odd number of samples

edge=(NUML mod 2)+1
axis=LMB(0:NUML-edge)
numln=n_elements(axis)

IF MC gt 0.001 THEN BEGIN  ; The macroturbulent velocity

; We assume a macroturbulent Gauss profile

    G = MACROTUR(MC,LMB(0:NUML-edge),Wl(1))
    FFTG=FFT(G, -1, /double)

    FOR I=0,3 DO BEGIN
        IF abs(mean(SPECTRA(*,I))) GT 1.e-25 THEN BEGIN  ; Just in case...
           FFTS=FFT(SPECTRA(0:NUML-edge,I), -1, /double)
           ; The convolution
           SPECTRA(0:NUML-edge,I)=DOUBLE(FFT(FFTS*FFTG, 1, /double))*NUMLn
           ;;;;;;;;;;;;;;;;;
           SPECTRA(0:NUML-edge,I)=shift(SPECTRA(0:NUML-edge,I),-NUML/2)
        ENDIF
    ENDFOR

ENDIF

; Instrumental profile
IF keyword_set(filter) THEN BEGIN
    DIMF=n_elements(filter)
    IF dimf GT 1 THEN BEGIN  ; The full filter profile is provided in this case

    ;CHECK IF FILTER SIZE IS LARGER THAN SPECTUM
    IF NUML lt DIMF then begin
    	;EXTEND THE PROFILES

    	SPECTRA_EXT = DBLARR(DIMF,4)
    	SPECTRA_EXT (DIMF/2-NUML/2:DIMF/2+NUML/2, *) = SPECTRA(0:NUML-1,*)
    	SPECTRA_EXT (0:DIMF/2-NUML/2-1, 0) = SPECTRA_EXT (DIMF/2-NUML/2, 0)
    	SPECTRA_EXT (DIMF/2+NUML/2+1:*, 0) = SPECTRA_EXT (DIMF/2+NUML/2, 0)
    	ENDIF ELSE SPECTRA_EXT = SPECTRA

        FFTF=FFT(filter, -1, /double)/total(filter)*dimf; Normalizing in area...
        SH=(WHERE(filter eq max(filter)))(0)
        FOR I=0,3 DO BEGIN
            IF abs(mean(SPECTRA_EXT(*,I))) GT 1.e-25 THEN BEGIN  ; Just in case...
               FFTS=FFT(SPECTRA_EXT(*,I), -1, /double)
               ;The convolution
               SPECTRA_EXT(*,I)=DOUBLE(FFT(FFTS*FFTF, 1, /double))
               ;;;;;;;;;;;;;;;;;
               SPECTRA_EXT(*,I)=shift(SPECTRA_EXT(*,I),-SH)
            ENDIF
        ENDFOR
        SPECTRA(*, *) = SPECTRA_EXT(DIMF/2-NUML/2:DIMF/2+NUML/2,*)

    ENDIF ELSE BEGIN ; Just the width is provided; the filter is assumed to be Gaussian
        filt=filtro(filter,LMB(0:NUML-edge))
        FFTF=FFT(filt,-1,/double)
        FOR I=0,3 DO BEGIN
            IF abs(mean(SPECTRA(*,I))) GT 1.e-25 THEN BEGIN  ; Just in case...
               FFTS=FFT(SPECTRA(0:NUML-edge,I), -1, /double)
               ;The convolution
               SPECTRA(0:NUML-edge,I)=DOUBLE(FFT(FFTS*FFTF, 1, /double))*numln
                ;;;;;;;;;;;;;;;;;
               SPECTRA(0:NUML-edge,I)=shift(SPECTRA(0:NUML-edge,I),-NUML/2)
            ENDIF
        ENDFOR
    ENDELSE
ENDIF

Spectra_comp(*,*,k) =  spectra

endfor ;COMPONENTS

;adding components

spectra(*) = 0

if n_comp gt 1 then begin
for k=0,n_comp-1 do begin
spectra = spectra + Spectra_comp(*,*,k)*fill_fractions(k)
endfor
endif else spectra = Spectra_comp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; ADDING THE STRAY-LIGHT PROFILE

IF keyword_set(slight) THEN BEGIN
    SPECTRA=SPECTRA*alpha+slight*(1D0-alpha)
ENDIF

END


;A ver. Lo que hay que calcular es ��[ eta0_6301 / eta0_6302 ]

;El eta 0 es el coeficiente de absorci�n de linea / por el coef. de abs. de continuo.
;Cuando divides los dos, te queda una parrafada de constantes entre las que hay que destacar la fuerza de oscilador y la poblaci�n del nivel inferior de la transici�n. Los potenciales de excitaci�n son muy parecidos: 3.654 frente a 3.686. Los niveles inferiores de las transiciones son los mismos. Te queda entonces solo el cociente de las fuerzas de oscilador.

;Log(gf) 6301.5 = -0.75 �(ojo, dependiendo de la literatura, puede variar un poco...)
;Log(gf) 6302.5 = -1.236
;g = 5 en las dos l�neas (2J+1)
;F_oscilador_6301 = 0.0356
;F_oscilador_6302 = 0.0116

;EL coeficiente de absorcion de linea es proporcional a la F_oscilador, osea que 6301 es mas profunda.
