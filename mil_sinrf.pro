;+
; NAME:
;	MIL_SINRF
;
; AUTHORS:
;   D. Orozco Suarez
;   orozco@iaa.es
;	and
;   J.C. del Toro Iniesta
; ADDRESS:
;      Instituto de Astrofisica de Andalucia (CSIC)
;      Apdo de Correos 3004, 18080 Granada, SPAIN
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
;   Added keywork CROSST for crostalk calculations. Nov, 2016. DOS
;   THIS ONE work in NLTE. Nov, 2020. DOS -> A. Dorantes PhD thesis
;-

pro MIL_SINRF,PARAM,WL,LMB,SPECTRA,TRIPLET=triplet,MU=mu,SLIGHT=slight,$
	FILTER=filter,AC_RATIO=ac_ratio,N_COMP=n_comp,ipbs=ipbs,crosst = crosst,NLTE=NLTE

COMMON QUANTIC,C_N

lines=wl(0)

IF NOT(KEYWORD_SET(mu)) THEN mu=1d0
IF NOT(KEYWORD_SET(AC_RATIO)) THEN AC_RATIO=0
IF (KEYWORD_SET(ipbs)) THEN he=1 else he=0

; Some initial settings

NUML=N_ELEMENTS(LMB)
SPECTRA=DBLARR(NUML,4)
SPECTRA_comp=DBLARR(NUML,4,n_comp)
vlight=2.99792458D+5   ; speed of light (cm/s)
RR=5.641895836D-1
CC=!dpi/180d0  ; conversion factor to radians

;alpha = param(10)
;ojo, ahora hay 4 mas
IF KEYWORD_SET(NLTE) THEN BEGIN
    ALPHA = PARAM[14] 
    NPAR = 15
ENDIF ELSE BEGIN
    ALPHA = PARAM[10]
    NPAR = 11
ENDELSE
;npar = 11 ;indice de nparams + 1 LTE
;npar = 15 ;indice de nparams + 1 NLTE

;FILL fraction in NLTE not tested yet
    FILL_FRACTIONS = FLTARR(N_COMP)
;IF N_COMP EQ 2 then begin
IF N_COMP GE 2 then begin
    FILL_FRACTIONS = FLTARR(N_COMP)
    FILL_FRACTIONS[1:*] = PARAM[NPAR*(INDGEN(N_COMP-1)+1)+NPAR-1]
    FILL_FRACTIONS[0] = 1d0 - TOTAL(FILL_FRACTIONS[1:*])
;ENDIF ELSE IF n_comp gt 2 THEN BEGIN
;    FILL_FRACTIONS[1:*] = PARAM[NPAR*(INDGEN(N_COMP-1)+1)+NPAR-1]
;    FILL_FRACTIONS[0] = PARAM[NPAR-1]
ENDIF ELSE FILL_FRACTIONS = 1D

FOR k = 0,N_COMP-1 DO BEGIN ;loop in components

;Model atmosphere parameters E0,MF,VL,LD,A,GM,AZI,B1,B2,MC,ALPHA
E00=param(npar*k)
MF=param(npar*k+1)
VL=param(npar*k+2)
LD=param(npar*k+3)
A=param(npar*k+4)
GM=param(npar*k+5)
AZI=param(npar*k+6)
B0=param(npar*k+7)
B1=param(npar*k+8)
MC=param(npar*k+9)
IF KEYWORD_SET(NLTE) THEN BEGIN
    A1  = PARAM(npar*k+10)
    ap1 = PARAM(npar*k+11)
    A2  = PARAM(npar*k+12)
    ap2 = PARAM(npar*k+13)
ENDIF

;recalculo los nñumeros cuanticos para paschen back
if KEYWORD_SET(IPBS) then IPBS,MF,LD

    ;SPECTRA is still valid. I have to add the additional terms which depends on for more parameters
    ;these parameters are followed after the straylight component so no 2-comp available for this mode right now

    ;first I need to calculate (α1 * µ * II + C)^−1
    ;for this reason ETAI will be first ME and then MENLTE

;****************************************************
;****************************************************
; FIRST PART OF EQN 9.116 (ME)
;****************************************************
;****************************************************

ETAI=1D0 & ETAQ=0D0 & ETAU=0D0 & ETAV=0D0  ; propagation matrix elements
RHOQ=0D0 & RHOU=0D0 & RHOV=0D0
VNULO=DBLARR(NUML)
;->
IF keyword_set(NLTE) THEN BEGIN
    ;NOW the MATRIX C is different. Instead of ETAI = 1 + ETAI -> ap1*MU + 1 + ETAI
    ETAI1_NLTE = 1D0 + ap1*MU
    ETAI2_NLTE = 1D0 + ap2*MU
ENDIF
;-<

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
            UU=u-ulos-(nupB(i))
            fvoigt, A, UU, H, F
            FI_P=FI_P+(C_N(IL).wep(i))*H
            SHI_P=SHI_P+(C_N(IL).wep(i))*F

        ENDFOR
        SHI_P=2d0*SHI_P

        ; SIGMA BLUE COMPONENTS

        FOR i=0,C_N(IL).N_SIG-1 DO BEGIN

            ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
            UU=u-ulos-(nubB(i))
            fvoigt, A, UU, H, F
            FI_B=FI_B+C_N(IL).web(i)*H
            SHI_B=SHI_B+C_N(IL).web(i)*F

        ENDFOR
        SHI_B=2d0*SHI_B

        ; SIGMA RED COMPONENTS

        FOR i=0,C_N(IL).N_SIG-1 DO BEGIN

            ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
            UU=u-ulos-(nurB(i))
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
    IF keyword_set(NLTE) THEN BEGIN
        ;->
        ETAI1_NLTE=ETAI1_NLTE+E0/2D0*(FI_P*sinis+(FI_B+FI_R)*(1D0+cosis)/2D0)
        ETAI2_NLTE=ETAI2_NLTE+E0/2D0*(FI_P*sinis+(FI_B+FI_R)*(1D0+cosis)/2D0)
        ;-<
    ENDIF
ENDFOR

;Etaq = etaq*(-1d0)
;Etau = etau*(-1d0)
;rhoq = rhoq*(-1d0)
;Etau = etau*(-1d0)

; THE STOKES PROFILES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Inverse of MAtrix C  EQN 9.91

DELTA=ETAI^2D0*(ETAI^2D0-ETAQ^2D0-ETAU^2D0-ETAV^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)-$
      (ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)^2D0
DELTAI=1D0/DELTA

;->
IF keyword_set(NLTE) THEN BEGIN
    DELTA1_NLTE=ETAI1_NLTE^2D0*(ETAI1_NLTE^2D0-ETAQ^2D0-ETAU^2D0-ETAV^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)-$
        (ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)^2D0
    DELTAI1_NLTE=1D0/DELTA1_NLTE

    DELTA2_NLTE=ETAI2_NLTE^2D0*(ETAI2_NLTE^2D0-ETAQ^2D0-ETAU^2D0-ETAV^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)-$
        (ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)^2D0
    DELTAI2_NLTE=1D0/DELTA2_NLTE
ENDIF
;<-

SPECTRA(*,0)=(B0+B1*DELTAI*ETAI*(ETAI^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)*MU)
SPECTRA(*,1)=(-B1*DELTAI*(ETAI^2D0*ETAQ+ETAI*(ETAV*RHOU-ETAU*RHOV)+$
  RHOQ*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*MU
SPECTRA(*,2)=(-B1*DELTAI*(ETAI^2D0*ETAU+ETAI*(ETAQ*RHOV-ETAV*RHOQ)+$
  RHOU*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*MU
SPECTRA(*,3)=(-B1*DELTAI*(ETAI^2D0*ETAV+ETAI*(ETAU*RHOQ-ETAQ*RHOU)+$
  RHOV*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)))*MU

IF keyword_set(NLTE) THEN BEGIN
    ;->
    SPECTRA(*,0)= SPECTRA(*,0) + $
                A1*(1d0 - ap1*MU*DELTAI1_NLTE*ETAI1_NLTE*$
                (ETAI1_NLTE^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)) - $
                A2*(1d0 - (1+ap2*MU)*DELTAI2_NLTE*ETAI2_NLTE*$
                (ETAI2_NLTE^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0))

    SPECTRA(*,1)= SPECTRA(*,1) + $
                A1*ap1*MU*DELTAI1_NLTE*(ETAI1_NLTE^2D0*ETAQ+ETAI1_NLTE*(ETAV*RHOU-ETAU*RHOV) + $
                RHOQ*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)) - $
                A2*(1+ap2*MU)*DELTAI2_NLTE*(ETAI2_NLTE^2D0*ETAQ+$
                ETAI2_NLTE*(ETAV*RHOU-ETAU*RHOV)+RHOQ*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV))

    SPECTRA(*,2)= SPECTRA(*,2) + $
                A1*ap1*MU*DELTAI1_NLTE*(ETAI1_NLTE^2D0*ETAU+$
                ETAI1_NLTE*(ETAQ*RHOV-ETAV*RHOQ)+$
                RHOU*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV))-$
                A2*(1+ap2*MU)*DELTAI2_NLTE*(ETAI2_NLTE^2D0*ETAU+$
                ETAI2_NLTE*(ETAQ*RHOV-ETAV*RHOQ)+$
                RHOU*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV))
    SPECTRA(*,3)= SPECTRA(*,3) + $
                A1*ap1*MU*DELTAI1_NLTE*(ETAI1_NLTE^2D0*ETAV+$
                ETAI1_NLTE*(ETAU*RHOQ-ETAQ*RHOU)+$
                RHOV*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV))-$
                A2*(1+ap2*MU)*DELTAI2_NLTE*(ETAI2_NLTE^2D0*ETAV+$
                ETAI2_NLTE*(ETAU*RHOQ-ETAQ*RHOU)+$
                RHOV*(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV))
ENDIF
;-<

;****************************************************
; CONVOLUTIONS WITH THE MACROTURBULENT VELOCITY AND THE INSTRUMENTAL PROFILE

; Taking care of the wavelength centering. Odd number of samples

edge=(NUML mod 2)+1
axis=LMB(0:NUML-edge)
numln=n_elements(axis)

IF MC gt 0.01 THEN BEGIN  ; The macroturbulent velocity
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


    IF keyword_set(crosst) THEN BEGIN
        IF keyword_set(NLTE) THEN BEGIN
            MatrixCrosst = [[1,0,0,0],[PARAM[15],1,0,0],[PARAM[16],0,1,0],[PARAM[17],0,0,1]]
        ENDIF ELSE BEGIN
        ;assumes 1 component !!!! It uses param to add crosstalk values
        ;derivatives are numerically calculated!!!!!
            MatrixCrosst = [[1,0,0,0],[PARAM[11],1,0,0],[PARAM[12],0,1,0],[PARAM[13],0,0,1]]
        ENDELSE
        FOR i=0,NUML-1 DO SPECTRA[i,*] = Matrixcrosst##SPECTRA[i,*]
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
