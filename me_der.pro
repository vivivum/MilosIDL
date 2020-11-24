;+
; NAME:
;	ME_DER
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
; PURPOSE: Synthesis of Stokes profiles with calculating response functions
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: ME_DER, PARAM, WL, LMB, SPECTRA,D_SPECTRA triplet=triplet,$
;                  mu=mu,slight=slight,filter=filter
;
; INPUTS:
;           PARAM: Array with the 11 model-atmosphere PARAMeters
;               eta0 = line-to-continuum absorption coefficient ratio
;               B = magnetic field strength       [Gauss]
;               vlos = line-of-sight velocity     [km/s]
;               dopp = Doppler width              [Angstroms]
;               aa = damping PARAMeter
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
;           SPECTRA: Array with the output Stokes profiles
;                       Dimensions: (N_ELEMENTS(AXIS), 4)
;           D_SPECTRA: Array with the output Stokes profiles Response Functions
;                       Dimensions: (N_ELEMENTS(AXIS), 4, N_ELEMENTS(PARAM))
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
;			Numerical = 2 computes both RFs for comparison
;
; CALLED ROUTINES
;           FVOIGT: Calculates the Voigt and Faraday-Voigt functions
;           MACROTUR: Evaluates a Gauss function
;           FILTRO: Outputs a Gaussian filter instrumental profile in case that only the
;                   width is provided
;
; MODIFICATION HISTORY:
; 	First beta version created, D Orozco Suárez (DOS) and J.C. del Toro Iniesta (JTI), 2007
;   First beta documentation version, JTI, June 2008
;   Improved documentation, DOS, 24 Feb, 2009
;   2015 n_comp
;   Added keywork CROSST for crostalk calculations. Nov, 2016. DOS
;   THIS ONE work in NLTE. July, 2019. DOS -> A. Dorantes PhD (only numerically)
;   Missing are the RFs to filter profile in NLTE!!!! It has to be implemented
;   Missing the Analytical response to alpha1 and alpha2 which are calculated numerically
;       Note that RF to macro is different when the previous RFs are added analytically
;   NOT TESTED FOR MORE THAN 2 COMPONENTS YET
;-

pro ME_DER,PARAM,WL,LMB,SPECTRA,D_SPECTRA,TRIPLET=triplet,SLIGHT=slight,$
           FILTER=filter,MU=mu,AC_RATIO=ac_ratio,N_COMP=n_comp,NUMERICAL=numerical,$
           IPBS=ipbs,crosst=crosst,nlte=nlte,saverfs=saverfs


COMMON QUANTIC,C_N

IF keyword_set(crosst) THEN NUMERICAL = 1

IF KEYWORD_SET(NUMERICAL) THEN D_N = NUMERICAL ELSE D_N = 0

IF (D_N EQ 0) OR (D_N EQ 2) then begin  ;ANALITIC RESPONSE FUNCTIONS

    IF NOT(KEYWORD_SET(MU)) THEN MU=1D
    IF NOT(KEYWORD_SET(AC_RATIO)) THEN AC_RATIO=0D

    LINES = WL[0]

    IF KEYWORD_SET(NLTE) THEN BEGIN
        ALPHA = PARAM[14] 
        NPAR = 15
        NTERMS = 15
    ENDIF ELSE BEGIN
        ALPHA = PARAM[10]
        NPAR = 11
        NTERMS = 11
    ENDELSE

    RR = 5.641895836D-1
    CC = !DPI/180D  ; conversion factor to radians
    NUML = N_ELEMENTS(LMB)
    D_SPECTRA = DBLARR(NUML,NTERMS,4) ;E0,MF,VL,LD,A,GM,AZI,DB,MC
    D_SPECTRA_COMP = DBLARR(NUML,NTERMS,4,N_COMP)
    SPECTRA = DBLARR(NUML,4)
    SPECTRA_COMP = DBLARR(NUML,4,N_COMP)
    VLIGHT = 2.99792458D+5   ; speed of light (cm/s)

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

        ;Model atmosphere PARAMeters E0,MF,VL,LD,A,GM,AZI,B1,B2,MC,ALPHA

        E00 = PARAM[npar*k]
        MF  = PARAM[npar*k+1]
        VL  = PARAM[npar*k+2]
        LD  = PARAM[npar*k+3]
        A   = PARAM[npar*k+4]
        GM  = PARAM[npar*k+5]
        AZI = PARAM[npar*k+6]
        B0  = PARAM[npar*k+7]
        B1  = PARAM[npar*k+8]
        MC  = PARAM[npar*k+9]
        IF KEYWORD_SET(NLTE) THEN BEGIN
            A1  = PARAM(npar*k+10)
            ap1 = PARAM(npar*k+11)
            A2  = PARAM(npar*k+12)
            ap2 = PARAM(npar*k+13)
        ENDIF

        if KEYWORD_SET(IPBS) then IPBS,MF,LD

        ; propagation matrix elements derivatives
        D_EI = DBLARR(NUML,7) ;D_ETAI
        D_EQ = DBLARR(NUML,7) ;D_ETAQ
        D_EU = DBLARR(NUML,7) ;D_ETAU
        D_EV = DBLARR(NUML,7) ;D_ETAV
        D_RQ = DBLARR(NUML,7) ;D_RHOQ
        D_RU = DBLARR(NUML,7) ;D_RHOU
        D_RV = DBLARR(NUML,7) ;D_RHOV

        ETAI = 1D0 & ETAQ = 0D0 & ETAU = 0D0 & ETAV = 0D0   ; propagation matrix elements
        RHOQ = 0D0 & RHOU = 0D0 & RHOV = 0D0
        VNULO = DBLARR(NUML)
        DNULO = DBLARR(NUML,4,3)
        FI_P = VNULO & FI_B = VNULO & FI_R = VNULO
        SHI_P = VNULO & SHI_B = VNULO & SHI_R = VNULO
        DFI = DNULO  ;U,A,B,LD
        DSHI = DNULO  ;U,A,B,LD
        ;->
        IF keyword_set(NLTE) THEN BEGIN
            ;NOW the MATRIC C is different. Instead of ETAI = 1 + ETAI -> ETAI = ap1*MU + 1 + ETAI
            ETAI1_NLTE = 1D0 + ap1*MU
            ETAI2_NLTE = 1D0 + ap2*MU
        ENDIF
        ;-<

        AZI = AZI*CC
        GM = GM*CC

        sinis=sin(GM)^2d0
        cosis=cos(GM)^2d0
        cosi=cos(GM)
        sina=sin(2D0*AZI)
        cosa=cos(2D0*AZI)

        sinda=cos(2D0*azi)*2D0*cc
        cosda=-sin(2D0*azi)*2D0*cc
        sindi=cos(gm)*sin(gm)*2D0*cc
        cosdi=-sin(gm)*cc

        ; Loop in spectral lines

        For IL=0,Lines-1 do begin

            ;ratio
            ratio = C_N(IL).Fo
            ; Ratio 6301. / 6302.  = 2.92000
            ; 1 / ratio = 0.342466 =  1. / 2.92000   ;teoretically 0.5444

            if (ratio ne 1) and (ac_ratio) then begin
                coeffa = [1.01702d-4, -9.33108d-3, 5.41421d-1, 3.11270d-2 ]
                coeffb = [2.26080d-1, 4.11138d0 ]
                if (E00 le 38) then begin

                    E0 = coeffa(0)*E00^3d0 + coeffa(1)*E00^2d0 + $
                        coeffa(2)*E00 + coeffa(3)
                    DE0 = 3.*coeffa(0)*E00^2d0 + 2*coeffa(1)*E00 + coeffa(2)

                endif else if (E00 lt 300) then begin

                    E0 = coeffb(0)*E00 + coeffb(1)
                    DE0 = coeffb(0)

                endif

            endif else begin
            ; Line strength
                E0=E00*ratio  ; C_N(IL).Fo is the relative line strength
                DE0=ratio  ; C_N(IL).Fo is the relative line strength
            endelse

            ; Wavelengths in units of the Doppler width
            U = (LMB-WL[IL+1])/LD

            ;frecuency shift for v line of sight
            ULOS = (VL*WL[IL+1])/(VLIGHT*LD)

            IF NOT(KEYWORD_SET(TRIPLET)) THEN BEGIN

                ; GENERAL MULTIPLET CASE
                ;***************************************************
                ; Zeeman splittings
                nubB=(wl(IL+1)^2D/ld)*(4.6686411D-13*C_N(IL).nub) ; sigma_b components
                nupB=(wl(IL+1)^2D/ld)*(4.6686411D-13*C_N(IL).nup) ; pi components
                nurB=(wl(IL+1)^2D/ld)*(4.6686411D-13*C_N(IL).nur) ; sigma_r components
                ;*************************************************

                ; Absorption and dispersion profiles initialization

                FI_P=VNULO & FI_B=VNULO & FI_R=VNULO
                SHI_P=VNULO & SHI_B=VNULO & SHI_R=VNULO

                ; ABSORPTION AND DISPERSION PROFILES

                ; PI COMPONENTS

                FOR i=0,C_N(IL).N_PI-1 do BEGIN

                    ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
                    UU=u-ulos-nupB(i)*MF

                    fvoigt, A, UU, H, F
                    FI_P = FI_P + C_N(IL).wep(i) * H   &     SHI_P=SHI_P+C_N(IL).wep(i)*F

                    ;Voigt and Voigt-Faraday Functions derivatives
                    dH_u = 4D0 * A * F - 2D0 * UU * H
                    dF_u = RR-A*H-2D0*UU*F
                    dH_a = -2D0*dF_u
                    dF_a = dH_u/2D0

                    ;numericamente es bastante más preciso y es que no es que esten mal las derivadas, sino que la fvoigt tiene error
                    HH=1e-3
                    fvoigt, A, UU - HH, H_i, F_i
                    fvoigt, A, UU + HH, H_d, F_d
                    dH_u = (H_d - H_i)/2./HH
                    dF_u = (F_d - F_i)/2./HH
                    fvoigt, A - HH, UU, H_i, F_i
                    fvoigt, A + HH, UU, H_d, F_d
                    dH_a = (H_d - H_i)/2./HH
                    dF_a = (F_d - F_i)/2./HH

                    ;Absortion and dispersion profiles derivatives
                    DFI(*,0,0)=DFI(*,0,0)+C_N(IL).wep(i)*dH_u*(-NUPB(i))  ;B
                    DSHI(*,0,0)=DSHI(*,0,0)+C_N(IL).wep(i)*dF_u*(-NUPB(i))  ;B
                    DFI(*,1,0)=DFI(*,1,0)+C_N(IL).wep(i)*dH_u*(-wl(IL+1))/(vlight*LD)  ;VLOS
                    DSHI(*,1,0)=DSHI(*,1,0)+C_N(IL).wep(i)*dF_u*(-wl(IL+1))/(vlight*LD) ;VLOS
                    DFI(*,2,0)=DFI(*,2,0)+C_N(IL).wep(i)*(dH_u*(-UU/LD));(-UU/LD))
                    DSHI(*,2,0)=DSHI(*,2,0)+C_N(IL).wep(i)*(dF_u*(-UU/LD));(-UU/LD))  ;B,U,LD,A
                    DFI(*,3,0)=DFI(*,3,0)+C_N(IL).wep(i)*dH_a ;a
                    DSHI(*,3,0)=DSHI(*,3,0)+C_N(IL).wep(i)*dF_a ;a

                ENDFOR
                SHI_P = 2d0*SHI_P
                DSHI(*,*,0) = 2d0*DSHI(*,*,0)
                ; SIGMA BLUE COMPONENTS

                FOR i=0,C_N(IL).N_SIG-1 DO BEGIN

                    ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
                    UU=u-ulos-nubB(i)*MF

                    fvoigt, A, UU, H, F
                    FI_B=FI_B+C_N(IL).web(i)*H & SHI_B=SHI_B+C_N(IL).web(i)*F

                    ;Voigt and Voigt-Faraday Functions derivatives
                    dH_u=4D0*A*F-2D0*UU*H
                    dF_u=RR-A*H-2D0*UU*F
                    dH_a=-2D0*dF_u
                    dF_a=dH_u/2D0

                    HH=1e-3
                    fvoigt, A, UU - HH, H_i, F_i
                    fvoigt, A, UU + HH, H_d, F_d
                    dH_u = (H_d - H_i)/2./HH
                    dF_u = (F_d - F_i)/2./HH
                    fvoigt, A - HH, UU, H_i, F_i
                    fvoigt, A + HH, UU, H_d, F_d
                    dH_a = (H_d - H_i)/2./HH
                    dF_a = (F_d - F_i)/2./HH

                    ;Absortion and dispersion profiles derivatives
                    DFI(*,1,1)=DFI(*,1,1)+C_N(IL).web(i)*dH_u*(-wl(IL+1))/(vlight*LD)
                    DSHI(*,1,1)=DSHI(*,1,1)+C_N(IL).web(i)*dF_u*(-wl(IL+1))/(vlight*LD)
                    DFI(*,3,1)=DFI(*,3,1)+C_N(IL).web(i)*dH_a
                    DSHI(*,3,1)=DSHI(*,3,1)+C_N(IL).web(i)*dF_a
                    DFI(*,0,1)=DFI(*,0,1)+C_N(IL).web(i)*dH_u*(-NUBB(i))
                    DSHI(*,0,1)=DSHI(*,0,1)+C_N(IL).web(i)*dF_u*(-NUBB(i))
                    DFI(*,2,1)=DFI(*,2,1)+C_N(IL).web(i)*(dH_u*(-UU/LD))
                    DSHI(*,2,1)=DSHI(*,2,1)+C_N(IL).web(i)*(dF_u*(-UU/LD))

                ENDFOR
                SHI_B=2d0*SHI_B
                DSHI(*,*,1)=2d0*DSHI(*,*,1)

                ; SIGMA RED COMPONENTS

                FOR i=0,C_N(IL).N_SIG-1 DO BEGIN

                    ; Wavelength axis shifted by Doppler and Zeeman, in Doppler width units
                    UU=u-ulos-nurB(i)*MF

                    fvoigt,A,UU, H, F ;R
                    FI_R=FI_R+C_N(IL).wer(i)*H & SHI_R=SHI_R+C_N(IL).wer(i)*F

                    ;Voigt and Voigt-Faraday Functions derivatives
                    dH_u=4D0*A*F-2D0*UU*H
                    dF_u=RR-A*H-2D0*UU*F
                    dH_a=-2D0*dF_u
                    dF_a=dH_u/2D0

                    HH=1e-3
                    fvoigt, A, UU - HH, H_i, F_i
                    fvoigt, A, UU + HH, H_d, F_d
                    dH_u = (H_d - H_i)/2./HH
                    dF_u = (F_d - F_i)/2./HH
                    fvoigt, A - HH, UU, H_i, F_i
                    fvoigt, A + HH, UU, H_d, F_d
                    dH_a = (H_d - H_i)/2./HH
                    dF_a = (F_d - F_i)/2./HH

                    ;Absortion and dispersion profiles derivatives
                    DFI(*,1,2)=DFI(*,1,2)+C_N(IL).wer(i)*dH_u*(-wl(IL+1))/(vlight*LD)
                    DSHI(*,1,2)=DSHI(*,1,2)+C_N(IL).wer(i)*dF_u*(-wl(IL+1))/(vlight*LD)
                    DFI(*,3,2)=DFI(*,3,2)+C_N(IL).wer(i)*dH_a
                    DSHI(*,3,2)=DSHI(*,3,2)+C_N(IL).wer(i)*dF_a
                    DFI(*,0,2)=DFI(*,0,2)+C_N(IL).wer(i)*dH_u*(-NURB(i))
                    DSHI(*,0,2)=DSHI(*,0,2)+C_N(IL).wer(i)*dF_u*(-NURB(i))
                    DFI(*,2,2)=DFI(*,2,2)+C_N(IL).wer(i)*(dH_u*(-UU/LD))
                    DSHI(*,2,2)=DSHI(*,2,2)+C_N(IL).wer(i)*(dF_u*(-UU/LD))
                ENDFOR
                SHI_R = 2d0*SHI_R
                DSHI(*,*,2) = 2d0*DSHI(*,*,2)

            ENDIF ELSE BEGIN

                ; TRIPLET CASE

                SHIF=(wl(IL+1)^2D0/LD)*(4.6686411D-13*C_N(IL).GEFF)

                ; PI COMPONENT

                fvoigt, A, U-ULOS, H, F
                FI_P=H & SHI_P=F
                SHI_P=2d0*SHI_P

                dH_u=4D0*A*F-2D0*(U-ULOS)*H
                dF_u=RR-A*H-2D0*(U-ULOS)*F
                dH_a=-2D0*dF_u
                dF_a=dH_u/2D0

                DFI(*,1,0)=dH_u*(-wl(IL+1))/(vlight*LD)
                DFI(*,3,0)=dH_a ;A
                DSHI(*,1,0)=2d0*dF_u*(-wl(IL+1))/(vlight*LD)
                DSHI(*,3,0)=2d0*dF_a ;A
                DFI(*,0,0)=0. ;B
                DSHI(*,0,0)=0. ;B
                DFI(*,2,0)=(dh_u*(-(U-ULOS)/LD)) ;LD
                DSHI(*,2,0)=2d0*(df_u*(-(U-ULOS)/LD)) ;LD

                ; SIGMA BLUE COMPONENT

                fvoigt, A, U-ULOS+SHIF*MF, H, F
                FI_B=H & SHI_B=F
                SHI_B=2d0*SHI_B

                dH_u=4D0*A*F-2D0*(U-ULOS+SHIF*MF)*H
                dF_u=RR-A*H-2D0*(U-ULOS+SHIF*MF)*F
                dH_a=-2D0*dF_u
                dF_a=dH_u/2D0

                DFI(*,1,1)=dH_u*(-wl(IL+1))/(vlight*LD)
                DFI(*,3,1)=dH_a ;A
                DSHI(*,1,1)=2d0*dF_u*(-wl(IL+1))/(vlight*LD)
                DSHI(*,3,1)=2d0*dF_a ;A
                DFI(*,0,1)=dh_u*(SHIF) ;B
                DSHI(*,0,1)=2d0*df_u*(SHIF) ;B
                DFI(*,2,1)=(dh_u*(-(U-ULOS+SHIF*MF)/LD)) ;LD
                DSHI(*,2,1)=2d0*(df_u*(-(U-ULOS+SHIF*MF)/LD)) ;LD

                ; SIGMA RED COMPONENT

                fvoigt, A, U-ULOS-SHIF*MF, H, F
                FI_R=H & SHI_R=F
                SHI_R=2d0*SHI_R

                dH_u=4D0*A*F-2D0*(U-ULOS-SHIF*MF)*H
                dF_u=RR-A*H-2D0*(U-ULOS-SHIF*MF)*F
                dH_a=-2D0*dF_u
                dF_a=dH_u/2D0

                DFI(*,0,2)=dh_u*(-SHIF) ;B
                DSHI(*,0,2)=2d0*df_u*(-SHIF) ;B
                DFI(*,1,2)=dH_u*(-wl(IL+1))/(vlight*LD)
                DSHI(*,1,2)=2d0*dF_u*(-wl(IL+1))/(vlight*LD)
                DFI(*,2,2)=(dh_u*(-(U-ULOS-SHIF*MF)/LD)) ;LD
                DSHI(*,2,2)=2d0*(df_u*(-(U-ULOS-SHIF*MF)/LD)) ;LD
                DFI(*,3,2)=dH_a ;A
                DSHI(*,3,2)=2d0*dF_a ;A

            ENDELSE

            ; Adding the necessary geometry and line strength

            ETAIN=E0/2D0*(FI_P*sinis+(FI_B+FI_R)*(1D0+cosis)/2D0)
            ETAQN=E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*cosa
            ETAUN=E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*sina
            ETAVN=E0/2D0*(FI_R-FI_B)*cosi
            RHOQN=E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*cosa
            RHOUN=E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*sina
            RHOVN=E0/2D0*(SHI_R-SHI_B)*cosi

            ETAI=ETAI+ETAIN
            ETAQ=ETAQ+ETAQN
            ETAU=ETAU+ETAUN
            ETAV=ETAV+ETAVN
            RHOQ=RHOQ+RHOQN
            RHOU=RHOU+RHOUN
            RHOV=RHOV+RHOVN

            IF keyword_set(NLTE) THEN BEGIN
                ;ETAI1N_NLTE=E0/2D0*(FI_P*sinis+(FI_B+FI_R)*(1D0+cosis)/2D0)
                ;ETAI2N_NLTE=E0/2D0*(FI_P*sinis+(FI_B+FI_R)*(1D0+cosis)/2D0)
                ;ETAI1_NLTE=ETAI1_NLTE+ETAI1N_NLTE
                ;ETAI2_NLTE=ETAI2_NLTE+ETAI2N_NLTE

                ;ETAI_NLTE=E0/2D0*(FI_P*sinis+(FI_B+FI_R)*(1D0+cosis)/2D0)
                ;ETAI1_NLTE=ETAI1_NLTE+ETAI_NLTE
                ;ETAI2_NLTE=ETAI2_NLTE+ETAI_NLTE

                ETAI1_NLTE=ETAI1_NLTE+ETAIN
                ETAI2_NLTE=ETAI2_NLTE+ETAIN
            ;ENDIF

            ;IF keyword_set(NLTE) THEN BEGIN

                ; Derivatives with respect eta0

                ; D_EI(*,0)=D_EI(*,0)+ETAIN/E0*DE0;C_N(IL).Fo
                ; D_EQ(*,0)=D_EQ(*,0)+ETAQN/E0*DE0;C_N(IL).Fo
                ; D_EU(*,0)=D_EU(*,0)+ETAUN/E0*DE0;C_N(IL).Fo
                ; D_EV(*,0)=D_EV(*,0)+ETAVN/E0*DE0;C_N(IL).Fo
                ; D_RQ(*,0)=D_RQ(*,0)+RHOQN/E0*DE0;C_N(IL).Fo
                ; D_RU(*,0)=D_RU(*,0)+RHOUN/E0*DE0;C_N(IL).Fo
                ; D_RV(*,0)=D_RV(*,0)+RHOVN/E0*DE0;C_N(IL).Fo

                ; ;;SAME because d(ap1*MU + 1 + ETAI)/d(e0) = d(ETAI)/d(e0)  !!!
                ; ;D_EI_NLTE(*,0) = D_EI_NLTE(*,0)+ETAIN/E0*DE0;C_N(IL).Fo
                ; ;D_EI_NLTE(*,0) = D_EI(*,0)+ETAIN/E0*DE0 

                ; ; Derivatives with respect B, VL, LDOPP, and A

                ; D_EI(*,1:4)=D_EI(*,1:4)+E0*(DFI(*,*,0)*sinis+(DFI(*,*,1)+DFI(*,*,2))*(1D0+cosis)/2D0)/2D0
                ; D_EQ(*,1:4)=D_EQ(*,1:4)+E0*(DFI(*,*,0)-(DFI(*,*,1)+DFI(*,*,2))/2D0)*sinis*cosa/2D0
                ; D_EU(*,1:4)=D_EU(*,1:4)+E0*(DFI(*,*,0)-(DFI(*,*,1)+DFI(*,*,2))/2D0)*sinis*sina/2D0
                ; D_EV(*,1:4)=D_EV(*,1:4)+E0*(DFI(*,*,2)-DFI(*,*,1))*cosi/2D0
                ; D_RQ(*,1:4)=D_RQ(*,1:4)+E0*(DSHI(*,*,0)-(DSHI(*,*,1)+DSHI(*,*,2))/2D0)*sinis*cosa/2D0
                ; D_RU(*,1:4)=D_RU(*,1:4)+E0*(DSHI(*,*,0)-(DSHI(*,*,1)+DSHI(*,*,2))/2D0)*sinis*sina/2D0
                ; D_RV(*,1:4)=D_RV(*,1:4)+E0*(DSHI(*,*,2)-DSHI(*,*,1))*cosi/2D0

                ; ; Derivatives with respect inclination

                ; D_EI(*,5)=D_EI(*,5)+E0/2D0*(FI_P*sindi+(FI_B+FI_R)*cosi*cosdi)
                ; D_EQ(*,5)=D_EQ(*,5)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sindi*cosa
                ; D_EU(*,5)=D_EU(*,5)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sindi*sina
                ; D_EV(*,5)=D_EV(*,5)+E0/2D0*(FI_R-FI_B)*cosdi
                ; D_RQ(*,5)=D_RQ(*,5)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sindi*cosa
                ; D_RU(*,5)=D_RU(*,5)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sindi*sina
                ; D_RV(*,5)=D_RV(*,5)+E0/2D0*(SHI_R-SHI_B)*cosdi

                ; ; Derivatives with respect azimuth

                ; D_EQ(*,6)=D_EQ(*,6)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*cosda
                ; D_EU(*,6)=D_EU(*,6)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*sinda
                ; D_RQ(*,6)=D_RQ(*,6)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*cosda
                ; D_RU(*,6)=D_RU(*,6)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*sinda

            ENDIF ;ELSE BEGIN
                
                        ; Derivatives with respect eta0

                D_EI(*,0)=D_EI(*,0)+ETAIN/E0*DE0;C_N(IL).Fo
                D_EQ(*,0)=D_EQ(*,0)+ETAQN/E0*DE0;C_N(IL).Fo
                D_EU(*,0)=D_EU(*,0)+ETAUN/E0*DE0;C_N(IL).Fo
                D_EV(*,0)=D_EV(*,0)+ETAVN/E0*DE0;C_N(IL).Fo
                D_RQ(*,0)=D_RQ(*,0)+RHOQN/E0*DE0;C_N(IL).Fo
                D_RU(*,0)=D_RU(*,0)+RHOUN/E0*DE0;C_N(IL).Fo
                D_RV(*,0)=D_RV(*,0)+RHOVN/E0*DE0;C_N(IL).Fo

                ; Derivatives with respect B, VL, LDOPP, and A

                D_EI(*,1:4)=D_EI(*,1:4)+E0*(DFI(*,*,0)*sinis+(DFI(*,*,1)+DFI(*,*,2))*(1D0+cosis)/2D0)/2D0
                D_EQ(*,1:4)=D_EQ(*,1:4)+E0*(DFI(*,*,0)-(DFI(*,*,1)+DFI(*,*,2))/2D0)*sinis*cosa/2D0
                D_EU(*,1:4)=D_EU(*,1:4)+E0*(DFI(*,*,0)-(DFI(*,*,1)+DFI(*,*,2))/2D0)*sinis*sina/2D0
                D_EV(*,1:4)=D_EV(*,1:4)+E0*(DFI(*,*,2)-DFI(*,*,1))*cosi/2D0
                D_RQ(*,1:4)=D_RQ(*,1:4)+E0*(DSHI(*,*,0)-(DSHI(*,*,1)+DSHI(*,*,2))/2D0)*sinis*cosa/2D0
                D_RU(*,1:4)=D_RU(*,1:4)+E0*(DSHI(*,*,0)-(DSHI(*,*,1)+DSHI(*,*,2))/2D0)*sinis*sina/2D0
                D_RV(*,1:4)=D_RV(*,1:4)+E0*(DSHI(*,*,2)-DSHI(*,*,1))*cosi/2D0

                ; Derivatives with respect inclination

                D_EI(*,5)=D_EI(*,5)+E0/2D0*(FI_P*sindi+(FI_B+FI_R)*cosi*cosdi)
                D_EQ(*,5)=D_EQ(*,5)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sindi*cosa
                D_EU(*,5)=D_EU(*,5)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sindi*sina
                D_EV(*,5)=D_EV(*,5)+E0/2D0*(FI_R-FI_B)*cosdi
                D_RQ(*,5)=D_RQ(*,5)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sindi*cosa
                D_RU(*,5)=D_RU(*,5)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sindi*sina
                D_RV(*,5)=D_RV(*,5)+E0/2D0*(SHI_R-SHI_B)*cosdi

                ; Derivatives with respect azimuth

                D_EQ(*,6)=D_EQ(*,6)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*cosda
                D_EU(*,6)=D_EU(*,6)+E0/2D0*(FI_P-(FI_B+FI_R)/2D0)*sinis*sinda
                D_RQ(*,6)=D_RQ(*,6)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*cosda
                D_RU(*,6)=D_RU(*,6)+E0/2D0*(SHI_P-(SHI_B+SHI_R)/2D0)*sinis*sinda

            ;ENDELSE
        ENDFOR

        ; THE STOKES PROFILES AND THEIR RESPONSE FUNCTIONS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; Are Stokes profiles here normalized to continuum intentity? NO

        ; Some terms
        GP1=(ETAI^2D0-ETAQ^2D0-ETAU^2D0-ETAV^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)
        GP2=(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)
        DT=ETAI^2D0*GP1-GP2^2D0
        DTI=1D0/DT
        GP3=(ETAI^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)
        GP4=ETAI^2D0*ETAQ+ETAI*(ETAV*RHOU-ETAU*RHOV)
        GP5=ETAI^2D0*ETAU+ETAI*(ETAQ*RHOV-ETAV*RHOQ)
        GP6=ETAI^2D0*ETAV+ETAI*(ETAU*RHOQ-ETAQ*RHOU)

        ;Stokes Profiles
        SPECTRA(*,0)=B0+B1*DTI*ETAI*GP3*MU
        SPECTRA(*,1)=(-B1*DTI*(GP4+RHOQ*GP2))*MU
        SPECTRA(*,2)=(-B1*DTI*(GP5+RHOU*GP2))*MU
        SPECTRA(*,3)=(-B1*DTI*(GP6+RHOV*GP2))*MU

        IF keyword_set(NLTE) THEN BEGIN
            GP1_NLTE1=(ETAI1_NLTE^2D0-ETAQ^2D0-ETAU^2D0-ETAV^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)
            GP2_NLTE1=(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)
            DT_NLTE1=ETAI1_NLTE^2D0*GP1_NLTE1-GP2_NLTE1^2D0
            DTI_NLTE1=1D0/DT_NLTE1
            GP3_NLTE1=(ETAI1_NLTE^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)
            GP4_NLTE1=ETAI1_NLTE^2D0*ETAQ+ETAI1_NLTE*(ETAV*RHOU-ETAU*RHOV)
            GP5_NLTE1=ETAI1_NLTE^2D0*ETAU+ETAI1_NLTE*(ETAQ*RHOV-ETAV*RHOQ)
            GP6_NLTE1=ETAI1_NLTE^2D0*ETAV+ETAI1_NLTE*(ETAU*RHOQ-ETAQ*RHOU)

            GP1_NLTE2=(ETAI2_NLTE^2D0-ETAQ^2D0-ETAU^2D0-ETAV^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)
            GP2_NLTE2=(ETAQ*RHOQ+ETAU*RHOU+ETAV*RHOV)
            DT_NLTE2=ETAI2_NLTE^2D0*GP1_NLTE2-GP2_NLTE2^2D0
            DTI_NLTE2=1D0/DT_NLTE2
            GP3_NLTE2=(ETAI2_NLTE^2D0+RHOQ^2D0+RHOU^2D0+RHOV^2D0)
            GP4_NLTE2=ETAI2_NLTE^2D0*ETAQ+ETAI2_NLTE*(ETAV*RHOU-ETAU*RHOV)
            GP5_NLTE2=ETAI2_NLTE^2D0*ETAU+ETAI2_NLTE*(ETAQ*RHOV-ETAV*RHOQ)
            GP6_NLTE2=ETAI2_NLTE^2D0*ETAV+ETAI2_NLTE*(ETAU*RHOQ-ETAQ*RHOU)

            SPECTRA(*,0)= SPECTRA(*,0) + $
                        A1*(1d0 - ap1*MU*DTI_NLTE1*ETAI1_NLTE*GP3_NLTE1) - $
                        A2*(1d0 - (1+ap2*MU)*DTI_NLTE2*ETAI2_NLTE*GP3_NLTE2)
            SPECTRA(*,1)= SPECTRA(*,1) + $
                        A1*ap1*MU*DTI_NLTE1*(GP4_NLTE1+RHOQ*GP2_NLTE1) - $
                        A2*(1+ap2*MU)*DTI_NLTE2*(GP4_NLTE2+RHOQ*GP2_NLTE2)
            SPECTRA(*,2)= SPECTRA(*,2) + $
                        A1*ap1*MU*DTI_NLTE1*(GP5_NLTE1+RHOU*GP2_NLTE1)-$
                        A2*(1+ap2*MU)*DTI_NLTE2*(GP5_NLTE2+RHOU*GP2_NLTE2)
            SPECTRA(*,3)= SPECTRA(*,3) + $
                        A1*ap1*MU*DTI_NLTE1*(GP6_NLTE1+RHOV*GP2_NLTE1)-$
                        A2*(1+ap2*MU)*DTI_NLTE2*(GP6_NLTE2+RHOV*GP2_NLTE2)

        ENDIF


        ; Loop to compute the Response Functions corresponding to E0, MF, VL, Doppl, A, GAMMA, and AZI

        FOR ITER=0,6 DO BEGIN

            DGP1=2D0*(ETAI*D_EI(*,ITER)-ETAQ*D_EQ(*,ITER)-ETAU*D_EU(*,ITER)-ETAV*D_EV(*,ITER)+$
                RHOQ*D_RQ(*,ITER)+RHOU*D_RU(*,ITER)+RHOV*D_RV(*,ITER))
            DGP2=RHOQ*D_EQ(*,ITER)+ETAQ*D_RQ(*,ITER)+RHOU*D_EU(*,ITER)+ETAU*D_RU(*,ITER)+$
                RHOV*D_EV(*,ITER)+ETAV*D_RV(*,ITER)

            D_DT=2D0*ETAI*D_EI(*,ITER)*GP1+ETAI^2D0*DGP1-2.*GP2*DGP2
            DGP3=2D0*(ETAI*D_EI(*,ITER)+RHOQ*D_RQ(*,ITER)+RHOU*D_RU(*,ITER)+RHOV*D_RV(*,ITER))

            DGP4=D_EI(*,ITER)*(2D0*ETAI*ETAQ+ETAV*RHOU-ETAU*RHOV)+ETAI^2D0*D_EQ(*,ITER)+$
                ETAI*(RHOU*D_EV(*,ITER)+ETAV*D_RU(*,ITER)-RHOV*D_EU(*,ITER)-ETAU*D_RV(*,ITER))

            DGP5=D_EI(*,ITER)*(2.*ETAI*ETAU+ETAQ*RHOV-ETAV*RHOQ)+ETAI^2.*D_EU(*,ITER)+$
                ETAI*(RHOV*D_EQ(*,ITER)+ETAQ*D_RV(*,ITER)-RHOQ*D_EV(*,ITER)-ETAV*D_RQ(*,ITER))

            DGP6=D_EI(*,ITER)*(2D0*ETAI*ETAV+ETAU*RHOQ-ETAQ*RHOU)+ETAI^2D0*D_EV(*,ITER)+$
                ETAI*(RHOQ*D_EU(*,ITER)+ETAU*D_RQ(*,ITER)-RHOU*D_EQ(*,ITER)-ETAQ*D_RU(*,ITER))

            D_SPECTRA(*,ITER,0)=B1*MU*((D_EI(*,ITER)*GP3+ETAI*DGP3)*DT-D_DT*ETAI*GP3)/DT^2D0
            D_SPECTRA(*,ITER,1)=-B1*MU*((DGP4+D_RQ(*,ITER)*GP2+RHOQ*DGP2)*DT-D_DT*(GP4+RHOQ*GP2))/DT^2D0
            D_SPECTRA(*,ITER,2)=-B1*MU*((DGP5+D_RU(*,ITER)*GP2+RHOU*DGP2)*DT-D_DT*(GP5+RHOU*GP2))/DT^2D0
            D_SPECTRA(*,ITER,3)=-B1*MU*((DGP6+D_RV(*,ITER)*GP2+RHOV*DGP2)*DT-D_DT*(GP6+RHOV*GP2))/DT^2D0

            IF keyword_set(NLTE) THEN BEGIN

                DGP1_NLTE1=2D0*(ETAI1_NLTE*D_EI(*,ITER)-ETAQ*D_EQ(*,ITER)-ETAU*D_EU(*,ITER)-ETAV*D_EV(*,ITER)+$
                    RHOQ*D_RQ(*,ITER)+RHOU*D_RU(*,ITER)+RHOV*D_RV(*,ITER))
                DGP2_NLTE1=RHOQ*D_EQ(*,ITER)+ETAQ*D_RQ(*,ITER)+RHOU*D_EU(*,ITER)+ETAU*D_RU(*,ITER)+$
                    RHOV*D_EV(*,ITER)+ETAV*D_RV(*,ITER)
                D_DT_NLTE1=2D0*ETAI1_NLTE*D_EI(*,ITER)*GP1_NLTE1+ETAI1_NLTE^2D0*DGP1_NLTE1-2.*GP2_NLTE1*DGP2_NLTE1
                DGP3_NLTE1=2D0*(ETAI1_NLTE*D_EI(*,ITER)+RHOQ*D_RQ(*,ITER)+RHOU*D_RU(*,ITER)+RHOV*D_RV(*,ITER))
                DGP4_NLTE1=D_EI(*,ITER)*(2D0*ETAI1_NLTE*ETAQ+ETAV*RHOU-ETAU*RHOV)+ETAI1_NLTE^2D0*D_EQ(*,ITER)+$
                    ETAI1_NLTE*(RHOU*D_EV(*,ITER)+ETAV*D_RU(*,ITER)-RHOV*D_EU(*,ITER)-ETAU*D_RV(*,ITER))
                DGP5_NLTE1=D_EI(*,ITER)*(2.*ETAI1_NLTE*ETAU+ETAQ*RHOV-ETAV*RHOQ)+ETAI1_NLTE^2.*D_EU(*,ITER)+$
                    ETAI1_NLTE*(RHOV*D_EQ(*,ITER)+ETAQ*D_RV(*,ITER)-RHOQ*D_EV(*,ITER)-ETAV*D_RQ(*,ITER))
                DGP6_NLTE1=D_EI(*,ITER)*(2D0*ETAI1_NLTE*ETAV+ETAU*RHOQ-ETAQ*RHOU)+ETAI1_NLTE^2D0*D_EV(*,ITER)+$
                    ETAI1_NLTE*(RHOQ*D_EU(*,ITER)+ETAU*D_RQ(*,ITER)-RHOU*D_EQ(*,ITER)-ETAQ*D_RU(*,ITER))

                DGP1_NLTE2=2D0*(ETAI2_NLTE*D_EI(*,ITER)-ETAQ*D_EQ(*,ITER)-ETAU*D_EU(*,ITER)-ETAV*D_EV(*,ITER)+$
                    RHOQ*D_RQ(*,ITER)+RHOU*D_RU(*,ITER)+RHOV*D_RV(*,ITER))
                DGP2_NLTE2=RHOQ*D_EQ(*,ITER)+ETAQ*D_RQ(*,ITER)+RHOU*D_EU(*,ITER)+ETAU*D_RU(*,ITER)+$
                    RHOV*D_EV(*,ITER)+ETAV*D_RV(*,ITER)
                D_DT_NLTE2=2D0*ETAI2_NLTE*D_EI(*,ITER)*GP1_NLTE2+ETAI2_NLTE^2D0*DGP1_NLTE2-2.*GP2_NLTE2*DGP2_NLTE2
                DGP3_NLTE2=2D0*(ETAI2_NLTE*D_EI(*,ITER)+RHOQ*D_RQ(*,ITER)+RHOU*D_RU(*,ITER)+RHOV*D_RV(*,ITER))
                DGP4_NLTE2=D_EI(*,ITER)*(2D0*ETAI2_NLTE*ETAQ+ETAV*RHOU-ETAU*RHOV)+ETAI2_NLTE^2D0*D_EQ(*,ITER)+$
                    ETAI2_NLTE*(RHOU*D_EV(*,ITER)+ETAV*D_RU(*,ITER)-RHOV*D_EU(*,ITER)-ETAU*D_RV(*,ITER))
                DGP5_NLTE2=D_EI(*,ITER)*(2.*ETAI2_NLTE*ETAU+ETAQ*RHOV-ETAV*RHOQ)+ETAI2_NLTE^2.*D_EU(*,ITER)+$
                    ETAI2_NLTE*(RHOV*D_EQ(*,ITER)+ETAQ*D_RV(*,ITER)-RHOQ*D_EV(*,ITER)-ETAV*D_RQ(*,ITER))
                DGP6_NLTE2=D_EI(*,ITER)*(2D0*ETAI2_NLTE*ETAV+ETAU*RHOQ-ETAQ*RHOU)+ETAI2_NLTE^2D0*D_EV(*,ITER)+$
                    ETAI2_NLTE*(RHOQ*D_EU(*,ITER)+ETAU*D_RQ(*,ITER)-RHOU*D_EQ(*,ITER)-ETAQ*D_RU(*,ITER))

                D_SPECTRA(*,ITER,0) = D_SPECTRA(*,ITER,0) - $
                            A1 * ap1 * MU * ((-1.) * D_DT_NLTE1*ETAI1_NLTE*GP3_NLTE1 / DT_NLTE1^2. + $
                            DTI_NLTE1*D_EI[*,ITER]*GP3_NLTE1 + DTI_NLTE1*ETAI1_NLTE*DGP3_NLTE1) + $ 
                            A2 * (1+ap2*MU) * ((-1.) * D_DT_NLTE2*ETAI2_NLTE*GP3_NLTE2 / DT_NLTE2^2.+ $
                            DTI_NLTE2*D_EI[*,ITER]*GP3_NLTE2 + DTI_NLTE2*ETAI2_NLTE*DGP3_NLTE2)

                D_SPECTRA(*,ITER,1)= D_SPECTRA(*,ITER,1) + $
                            A1 * ap1 * MU * ( (-1.) * D_DT_NLTE1*(GP4_NLTE1+RHOQ*GP2_NLTE1) / DT_NLTE1^2. + $
                            DTI_NLTE1 * (DGP4_NLTE1 + D_RQ[*,ITER]*GP2_NLTE1 + RHOQ*DGP2_NLTE1) ) - $
                            A2 * (1+ap2*MU) * ( (-1.) * D_DT_NLTE2*(GP4_NLTE2+RHOQ*GP2_NLTE2)  / DT_NLTE2^2.+ $
                            DTI_NLTE2 * (DGP4_NLTE2 + D_RQ[*,ITER]*GP2_NLTE2 + RHOQ*DGP2_NLTE2) )

                D_SPECTRA(*,ITER,2)= D_SPECTRA(*,ITER,2) + $
                            A1 * ap1 * MU * ( (-1.) * D_DT_NLTE1*(GP5_NLTE1+RHOU*GP2_NLTE1) / DT_NLTE1^2. + $
                            DTI_NLTE1 * (DGP5_NLTE1 + D_RU[*,ITER]*GP2_NLTE1 + RHOU*DGP2_NLTE1) ) - $
                            A2 * (1+ap2*MU) * ( (-1.) * D_DT_NLTE2*(GP5_NLTE2+RHOU*GP2_NLTE2)  / DT_NLTE2^2.+ $
                            DTI_NLTE2 * (DGP5_NLTE2 + D_RU[*,ITER]*GP2_NLTE2 + RHOU*DGP2_NLTE2) )

                D_SPECTRA(*,ITER,3)= D_SPECTRA(*,ITER,3) + $
                            A1 * ap1 * MU * ( (-1.) * D_DT_NLTE1*(GP6_NLTE1+RHOV*GP2_NLTE1) / DT_NLTE1^2. + $
                            DTI_NLTE1 * (DGP6_NLTE1 + D_RV[*,ITER]*GP2_NLTE1 + RHOV*DGP2_NLTE1) ) - $
                            A2 * (1+ap2*MU) * ( (-1.) * D_DT_NLTE2*(GP6_NLTE2+RHOV*GP2_NLTE2)  / DT_NLTE2^2.+ $
                            DTI_NLTE2 * (DGP6_NLTE2 + D_RV[*,ITER]*GP2_NLTE2 + RHOV*DGP2_NLTE2) )

            ENDIF

        ENDFOR 

        ;Response Functions with respect B0 and B1
        D_SPECTRA(*,7,0)=1D0 & D_SPECTRA(*,8,0)=DTI*ETAI*GP3*MU
        D_SPECTRA(*,7,1)=0D0 & D_SPECTRA(*,8,1)=-DTI*(GP4+RHOQ*GP2)*MU
        D_SPECTRA(*,7,2)=0D0 & D_SPECTRA(*,8,2)=-DTI*(GP5+RHOU*GP2)*MU
        D_SPECTRA(*,7,3)=0D0 & D_SPECTRA(*,8,3)=-DTI*(GP6+RHOV*GP2)*MU

        D_SPECTRA(*,6,0)=0D0  ;AZIMUTH STOKES I
        D_SPECTRA(*,6,3)=0D0  ;AZIMUTH STOKES V

        IF keyword_set(NLTE) THEN BEGIN

            ;Response Functions with respect A1, A2, ap1, ap2 
            print,A1,ap1,A2,ap2

            ;A1
            D_SPECTRA(*,10,0) = (1d0 - ap1*MU*DTI_NLTE1*ETAI1_NLTE*GP3_NLTE1)
            D_SPECTRA(*,10,1) = ap1*MU*DTI_NLTE1*(GP4_NLTE1+RHOQ*GP2_NLTE1)
            D_SPECTRA(*,10,2) = ap1*MU*DTI_NLTE1*(GP5_NLTE1+RHOU*GP2_NLTE1)
            D_SPECTRA(*,10,3) = ap1*MU*DTI_NLTE1*(GP6_NLTE1+RHOV*GP2_NLTE1)
            ;A2
            D_SPECTRA(*,12,0) = -(1d0 - (1+ap2*MU)*DTI_NLTE2*ETAI2_NLTE*GP3_NLTE2)
            D_SPECTRA(*,12,1) = -(1+ap2*MU)*DTI_NLTE2*(GP4_NLTE2+RHOQ*GP2_NLTE2)
            D_SPECTRA(*,12,2) = -(1+ap2*MU)*DTI_NLTE2*(GP5_NLTE2+RHOU*GP2_NLTE2)
            D_SPECTRA(*,12,3) = -(1+ap2*MU)*DTI_NLTE2*(GP6_NLTE2+RHOV*GP2_NLTE2)
            ;NOT CORRECT.(BELOW) I RUN 11 and 13 NUMERICALLY
            ;ap1 

            ; D_SPECTRA(*,11,0) = -A1*MU*DTI_NLTE1*ETAI1_NLTE*GP3_NLTE1 + $
            ;                     A1*ap1*MU*MU*(DTI_NLTE1*ETAI1_NLTE*GP3_NLTE1)^2 
            ; D_SPECTRA(*,11,1) = -A1*MU*DTI_NLTE1*(GP4_NLTE1+RHOQ*GP2_NLTE1) - $
            ;                     A1*ap1*MU*MU*(DTI_NLTE1*(GP4_NLTE1+RHOQ*GP2_NLTE1))^2*ap1^2 
            ; D_SPECTRA(*,11,2) = -A1*MU*DTI_NLTE1*(GP5_NLTE1+RHOU*GP2_NLTE1) - $
            ;                     A1*ap1*MU*MU*(DTI_NLTE1*(GP5_NLTE1+RHOU*GP2_NLTE1))^2*ap1^2
            ; D_SPECTRA(*,11,3) = -A1*MU*DTI_NLTE1*(GP6_NLTE1+RHOV*GP2_NLTE1) - $
            ;                     A1*ap1*MU*MU*(DTI_NLTE1*(GP6_NLTE1+RHOV*GP2_NLTE1))^2*ap1^2  

            D_SPECTRA(*,11,0) = -A1*MU*DTI_NLTE1*ETAI1_NLTE*GP3_NLTE1 - $ ;primer termino
                                A1*ap1*MU*DTI_NLTE1*GP3_NLTE1 -  $;segundo termino
                                A1*ap1*MU*DTI_NLTE1*ETAI1_NLTE*2*ETAI1_NLTE + $ ;tercer termino
                                A1*ap1*MU*DTI_NLTE1^2*ETAI1_NLTE*GP3_NLTE1*(2d0*ETAI1_NLTE*GP1_NLTE1 + $
                                ETAI1_NLTE^2D0*2*ETAI1_NLTE)
            D_SPECTRA(*,11,1) = A1*MU*DTI_NLTE1*(GP4_NLTE1+RHOQ*GP2_NLTE1) + $
                                A1*MU*ap1*DTI_NLTE1*(2.*ETAI1_NLTE*ETAQ) + $
                                A1*MU*ap1*DTI_NLTE1*(ETAV*RHOU-ETAU*RHOV) - $ ;tercer termino
                                A1*MU*ap1*(GP4_NLTE1+RHOQ*GP2_NLTE1)*DTI_NLTE1^2*(2d0*ETAI1_NLTE*GP1_NLTE1 + $
                                ETAI1_NLTE^2D0*2*ETAI1_NLTE)
            D_SPECTRA(*,11,2) = A1*MU*DTI_NLTE1*(GP5_NLTE1+RHOU*GP2_NLTE1) + $ ;primer termino
                                A1*MU*ap1*DTI_NLTE1*(2.*ETAI1_NLTE*ETAU) + $ ;segundo termino
                                A1*MU*ap1*DTI_NLTE1*(ETAQ*RHOV-ETAV*RHOQ) - $ ;tercer termino
                                A1*MU*ap1*(GP5_NLTE1+RHOU*GP2_NLTE1)*DTI_NLTE1^2*(2d0*ETAI1_NLTE*GP1_NLTE1 + $
                                ETAI1_NLTE^2D0*2*ETAI1_NLTE)
            D_SPECTRA(*,11,3) = A1*MU*DTI_NLTE1*(GP6_NLTE1+RHOV*GP2_NLTE1) + $ ;primer termino
                                A1*MU*ap1*DTI_NLTE1*(2.*ETAI1_NLTE*ETAV) + $ ;segundo termino
                                A1*MU*ap1*DTI_NLTE1*(ETAU*RHOQ-ETAQ*RHOU) - $ ;tercer termino
                                A1*MU*ap1*(GP6_NLTE1+RHOV*GP2_NLTE1)*DTI_NLTE1^2*(2d0*ETAI1_NLTE*GP1_NLTE1 + $
                                ETAI1_NLTE^2D0*2*ETAI1_NLTE)

            ;ap2     
            ; D_SPECTRA(*,13,0) = A2*MU*DTI_NLTE2*ETAI2_NLTE*GP3_NLTE2 - $
            ;                     A2*(1+ap2*MU)*MU*(DTI_NLTE2*ETAI2_NLTE*GP3_NLTE2)^2 
            ; D_SPECTRA(*,13,1) = A2*MU*DTI_NLTE2*(GP4_NLTE2+RHOQ*GP2_NLTE2) - $
            ;                     A2*(1+ap2*MU)*MU*(DTI_NLTE2*(GP4_NLTE2+RHOQ*GP2_NLTE2))^2                      
            ; D_SPECTRA(*,13,2) = A2*MU*DTI_NLTE2*(GP5_NLTE2+RHOU*GP2_NLTE2) - $
            ;                     A2*(1+ap2*MU)*MU*(DTI_NLTE2*(GP5_NLTE2+RHOU*GP2_NLTE2))^2 
            ; D_SPECTRA(*,13,3) = A2*MU*DTI_NLTE2*(GP6_NLTE2+RHOV*GP2_NLTE2) - $
            ;                     A2*(1+ap2*MU)*MU*(DTI_NLTE2*(GP6_NLTE2+RHOV*GP2_NLTE2))^2 

            D_SPECTRA(*,13,0) = A2*MU*DTI_NLTE2*ETAI2_NLTE*GP3_NLTE2 + $ ;primer termino
                                A2*(1+ap2*MU)*MU*DTI_NLTE2*GP3_NLTE2 +  $;segundo termino
                                A2*(1+ap2*MU)*MU*DTI_NLTE2*ETAI2_NLTE*2*ETAI2_NLTE - $ ;tercer termino
                                A2*(1+ap2*MU)*MU*DTI_NLTE2^2*ETAI2_NLTE*GP3_NLTE2*(2d0*ETAI2_NLTE*GP1_NLTE2 + $
                                ETAI2_NLTE^2D0*2*ETAI2_NLTE)
            D_SPECTRA(*,13,1) = -A2*MU*DTI_NLTE2*(GP4_NLTE2+RHOQ*GP2_NLTE2) - $
                                A2*MU*(1+ap2*MU)*DTI_NLTE2*(2.*ETAI2_NLTE*ETAQ) - $
                                A2*MU*(1+ap2*MU)*DTI_NLTE2*(ETAV*RHOU-ETAU*RHOV) + $ ;tercer termino
                                A2*MU*(1+ap2*MU)*(GP4_NLTE2+RHOQ*GP2_NLTE2)*DTI_NLTE2^2*(2d0*ETAI2_NLTE*GP1_NLTE2 + $
                                ETAI2_NLTE^2D0*2*ETAI2_NLTE)
            D_SPECTRA(*,13,2) = -A2*MU*DTI_NLTE2*(GP5_NLTE2+RHOU*GP2_NLTE2) - $ ;primer termino
                                A2*MU*(1+ap2*MU)*DTI_NLTE2*(2.*ETAI2_NLTE*ETAU) - $ ;segundo termino
                                A2*MU*(1+ap2*MU)*DTI_NLTE2*(ETAQ*RHOV-ETAV*RHOQ) + $ ;tercer termino
                                A2*MU*(1+ap2*MU)*(GP5_NLTE2+RHOU*GP2_NLTE2)*DTI_NLTE2^2*(2d0*ETAI2_NLTE*GP1_NLTE2 + $
                                ETAI2_NLTE^2D0*2*ETAI2_NLTE)
            D_SPECTRA(*,13,3) = -A2*MU*DTI_NLTE2*(GP6_NLTE2+RHOV*GP2_NLTE2) - $ ;primer termino
                                A2*MU*(1+ap2*MU)*DTI_NLTE2*(2.*ETAI2_NLTE*ETAV) - $ ;segundo termino
                                A2*MU*(1+ap2*MU)*DTI_NLTE2*(ETAU*RHOQ-ETAQ*RHOU) + $ ;tercer termino
                                A2*MU*(1+ap2*MU)*(GP6_NLTE2+RHOV*GP2_NLTE2)*DTI_NLTE2^2*(2d0*ETAI2_NLTE*GP1_NLTE2 + $
                                ETAI2_NLTE^2D0*2*ETAI2_NLTE)

        ;LAS CALCULO NUMERICAMENTE LA 11 y la 13 PORQUE NO ME SALE ANALITICAMENTE. ALGO DEBE ESTAR PASANDO.
        dom = "false" & goto,jum 
        ;dom = "true"
            H=0.001
            PARAMN=PARAM
            IP = [11,13]
            FOR Ipc = 0,n_elements(IP)-1 do begin
                IF (ABS(PARAM(IP[Ipc])) gt 1e-2) THEN PARAMN(IP[Ipc])=PARAM(IP[Ipc])*(1.+H) $
                    ELSE PARAMN(IP[Ipc])=PARAM(IP[Ipc])+H
                mil_sinrf,PARAMn,wl,lmb,yfitD,triplet=triplet,slight=slight,filter=filter,$
                        AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
                IF (ABS(PARAM(IP[Ipc])) gt 1e-2) THEN PARAMN(IP[Ipc])=PARAM(IP[Ipc])*(1.-H) $
                    ELSE PARAMN(IP[Ipc])=PARAM(IP[Ipc])-H
                mil_sinrf,PARAMn,wl,lmb,yfitI,triplet=triplet,slight=slight,filter=filter,$
                        AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
                FOR J=0,3 DO D_SPECTRA(*,IP[Ipc],J)=(YFITD(*,J)-YFITI(*,J))/(2.*(PARAM(IP[Ipc])-PARAMN(IP[Ipc])))
            ENDFOR
            jum:
        ENDIF

        ;stop
        ; CONVOLUTIONS WITH THE MACROTURBULENT VELOCITY AND THE INSTRUMENTAL PROFILE

        ; Taking care of the wavelength centering. Odd number of samples

        edge=(NUML mod 2)+1
        axis=LMB(0:NUML-edge)
        numln=n_elements(axis)

        ;MACRO NO FUNCIONA CON 2 COMPONENTES!!!!!!!!!!!

        IF MC gt 0.01 THEN BEGIN  ; The macroturbulent velocity

            ; We assume a macroturbulent Gauss profile
            G1 = MACROTUR(MC,LMB(0:NUML-edge),wl(1))         ; Gauss function
            G2 = MACROTUR(MC,LMB(0:NUML-edge),wl(1),/deriv)  ; Derivative of the Gauss F. with respect macro

            ; some nonesense calculus
            g1=g1;*(lmb(1)-lmb(0))
            g2=g2;*(lmb(1)-lmb(0))

            ; Fast Fourier Transform
            FFTG1=FFT(G1, -1, /double)
            FFTG2=FFT(G2, -1, /double)

            ; Convolution of Stokes profiles and RF to MACROVELOCITY
            FOR ITER=0,3 DO BEGIN
                IF abs(mean(SPECTRA(*,ITER))) GT 1.e-25 THEN BEGIN  ; Just in case...
                    FFTS=FFT(SPECTRA(0:NUML-edge,ITER), -1, /double)
                    ; The convolution of macrovelocity Response Function
                    D_SPECTRA(0:NUML-edge,9,ITER)=FFT(FFTS*FFTG2, 1, /double)*NUMLn
                    D_SPECTRA(0:NUML-edge,9,ITER)=shift(D_SPECTRA(0:NUML-edge,9,ITER),-NUML/2)
                    ; The convolution of the Stokes profiles
                    SPECTRA(0:NUML-edge,ITER)=FFT(FFTS*FFTG1, 1, /double)*NUMLn
                    SPECTRA(0:NUML-edge,ITER)=shift(SPECTRA(0:NUML-edge,ITER),-NUML/2)
                ENDIF
            ENDFOR

            ; Convolution of Stokes profiles Response Functions (except that for the macrovelocity)
            ;dd=d_spectra
            FOR PAR=0,3 DO BEGIN  ;loop in Stokes profiles
                FOR ITER=0,6 DO BEGIN    ;loop in model PARAMeters
                    If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin  ; Just in case...
                        FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                        ; The convolution of Response Functions
                        D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTG1,1,/double)*NUMLn
                        D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                    ENDIF
                ENDFOR
                ;Parameter 7 is missing because it is a constant)
                ITER = 8
                If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin  ; Just in case...
                    FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                    ; The convolution of Response Functions
                    D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTG1,1,/double)*NUMLn
                    D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                ENDIF

                IF keyword_set(NLTE) THEN BEGIN
                ;OJO. COMENTO PORQUE la 11 y 13 estan numericas arriba lo que significa que no hay que meterles la macro!!!!!
                    if dom ne 'true' then begin
                        FOR ITER=10,13 DO BEGIN    ;loop in model PARAMeters
                            ;If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin  ; Just in case...
                                FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                                ; The convolution of Response Functions
                                D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTG1,1,/double)*NUMLn
                                D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                            ;ENDIF
                        ENDFOR
                    endif else begin
                        ITER = 10 ;***************
                        If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin  ; Just in case...
                            FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTG1,1,/double)*NUMLn
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                        ENDIF
                        ITER = 12
                        If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin  ; Just in case...
                            FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTG1,1,/double)*NUMLn
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                        ENDIF
                    endelse
                ENDIF
            ENDFOR
        ENDIF

        IF KEYWORD_SET(filter) THEN BEGIN ;COMPONENTS SHOULD BE ONE AT THE END!!!!
            DIMF=n_elements(filter)
            IF keyword_set(NLTE) THEN print,'NLTE WITH FILTER NOT IMPLEM;EMNTED YET'
            If dimf gt 1 then begin  ; The full filter profile is provided in this case

                ;CHECK IF FILTER SIZE IS LARGER THAN SPECTUM
                IF NUML lt DIMF then begin
                ;EXTEND THE PROFILES

                    SPECTRA_EXT = DBLARR(DIMF,4)
                    SPECTRA_EXT (DIMF/2-NUML/2:DIMF/2+NUML/2, *) = SPECTRA(0:NUML-1,*)
                    SPECTRA_EXT (0:DIMF/2-NUML/2-1, 0) = SPECTRA_EXT (DIMF/2-NUML/2, 0)
                    SPECTRA_EXT (DIMF/2+NUML/2+1:*, 0) = SPECTRA_EXT (DIMF/2+NUML/2, 0)

                    D_SPECTRA_EXT = DBLARR(DIMF,11,4)
                    D_SPECTRA_EXT (DIMF/2-NUML/2:DIMF/2+NUML/2, *,*) = D_SPECTRA(0:NUML-1,*,*)
                    D_SPECTRA_EXT (0:DIMF/2-NUML/2-1, 8,0) = D_SPECTRA_EXT (DIMF/2-NUML/2, 8,0)
                    D_SPECTRA_EXT (DIMF/2+NUML/2+1:*, 8,0) = D_SPECTRA_EXT (DIMF/2+NUML/2, 8,0)
                ENDIF ELSE BEGIN
                        D_SPECTRA_EXT = SPECTRA
                        SPECTRA_EXT = SPECTRA
                ENDELSE

                FFTF=FFT(filter, -1, /double)/total(filter)*dimf; Normalizing in area...
                SH=(WHERE(filter eq max(filter)))(0)

                ; Convolution of the Stokes profiles

                FOR I=0,3 DO BEGIN
                    IF abs(mean(SPECTRA_EXT(*,I))) GT 1.e-25 THEN BEGIN  ; Just in case...
                        FFTS=FFT(SPECTRA_EXT(*,I), -1, /double)
                        ;The convolution
                        SPECTRA_EXT(*,I)=DOUBLE(FFT(FFTS*FFTF, 1, /double))
                        SPECTRA_EXT(*,I)=shift(SPECTRA_EXT(*,I),-SH)
                    ENDIF
                ENDFOR
                SPECTRA(*, *) = SPECTRA_EXT(DIMF/2-NUML/2:DIMF/2+NUML/2,*)

                ; Convolution of the Response Functions (except for Stray light and S0)
                ; NLTE PARAMETERS MISSING
                FOR PAR=0,3 DO BEGIN
                    FOR ITER=0,6 DO BEGIN
                        If abs(Mean(D_SPECTRA_EXT(*,ITER,PAR))) GT 1.e-25 then begin ; Just in case
                            FFTD=FFT(D_SPECTRA_EXT(*,ITER,PAR),-1,/double)
                            ;The convolution
                            D_SPECTRA_EXT(*,ITER,PAR)=FFT(FFTD*FFTF,1,/double)
                            D_SPECTRA_EXT(*,ITER,PAR)=shift(D_SPECTRA_EXT(*,ITER,PAR),-SH)
                        ENDIF
                    ENDFOR
                    FOR ITER=8,9 DO BEGIN  ;Y AHORA LO TENGO AQUI SIN COMENTAR!!!!!!!!!!!!!
                        If abs(Mean(D_SPECTRA_EXT(*,ITER,PAR))) GT 1.e-25 then begin ; Just in case
                            FFTD=FFT(D_SPECTRA_EXT(*,ITER,PAR),-1,/double)
                            ;The convolution
                            D_SPECTRA_EXT(*,ITER,PAR)=FFT(FFTD*FFTF,1,/double)
                            D_SPECTRA_EXT(*,ITER,PAR)=shift(D_SPECTRA_EXT(*,ITER,PAR),-SH)
                        ENDIF
                    ENDFOR
                ENDFOR
                D_SPECTRA(*, *, *) = D_SPECTRA_EXT(DIMF/2-NUML/2:DIMF/2+NUML/2, *, *)

            ENDIF ELSE BEGIN  ; Just the width is provided; the filter is assumed to be Gaussian

                filt=filtro(filter,LMB(0:NUML-edge))
                FFTF=FFT(filt,-1,/double)

                ; Convolution of the Stokes profiles
                FOR ITER=0,3 DO BEGIN
                    IF abs(mean(SPECTRA(*,ITER))) GT 1.e-25 THEN BEGIN ; Just in case
                        FFTS=FFT(SPECTRA(0:NUML-edge,ITER), -1, /double)
                            ;The convolution
                        SPECTRA(0:NUML-edge,ITER)=FFT(FFTS*FFTF, 1, /double)*NUMLn
                        SPECTRA(0:NUML-edge,ITER)=shift(SPECTRA(0:NUML-edge,ITER),-NUML/2)
                    ENDIF
                ENDFOR

                ; Convolution of the Response Functions (except for Stray light and S0)

                FOR PAR=0,3 DO BEGIN
                    FOR ITER=0,6 DO BEGIN
                        If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin ; Just in case
                            FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                            ;The convolution
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTF,1,/double)*NUMLn
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                        ENDIF
                    ENDFOR
                    FOR ITER=8,9 DO BEGIN
                                    If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin ; Just in case
                            FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                            ;The convolution
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTF,1,/double)*NUMLn
                            D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                        ENDIF
                    ENDFOR
                    IF keyword_set(NLTE) THEN BEGIN
                        FOR ITER=10,13 DO BEGIN    ;loop in model PARAMeters
                            If abs(Mean(D_SPECTRA(*,ITER,PAR))) GT 1.e-25 then begin  ; Just in case...
                                FFTD=FFT(D_SPECTRA(0:NUML-edge,ITER,PAR),-1,/double)
                                ; The convolution of Response Functions
                                D_SPECTRA(0:NUML-edge,ITER,PAR)=FFT(FFTD*FFTF,1,/double)*NUML
                                D_SPECTRA(0:NUML-edge,ITER,PAR)=shift(D_SPECTRA(0:NUML-edge,ITER,PAR),-NUML/2)
                            ENDIF
                        ENDFOR
                    ENDIF
                ENDFOR
            ENDELSE
        ENDIF

        SPECTRA_COMP(*,*,k) =  SPECTRA
        D_SPECTRA_COMP(*,*,*,k) = D_SPECTRA

    ENDFOR ;COMPONENTS


    ;adding components

    REPLICATE_INPLACE, SPECTRA, 0D0
    REPLICATE_INPLACE, D_SPECTRA, 0D0

    if n_comp gt 1 then begin
        if keyword_set(NLTE) then print,'NLTE NOT READY FOR MORE THAN 1 COMPONENT'
        for k=0,n_comp-1 do begin
            spectra = spectra + Spectra_comp(*,*,k)*fill_fractions(k)
        endfor
        ;la derivada no se suma!
        for k=0,n_comp-1 do begin
            D_SPECTRA(*,11*k:11*(k+1)-1,*) =  D_SPECTRA_comp(*,*,*,k)*fill_fractions(k)
        endfor
        ;Derivative with respefct fill fractions
        if n_comp eq 2 then begin
            D_SPECTRA(*,11+10,*) = Spectra_comp(*,*,1) - Spectra_comp(*,*,0)
        endif else if n_comp gt 2 then begin
            for k=0,n_comp-1 do D_SPECTRA(*,11*k+10,*) = Spectra_comp(*,*,k)
        endif

    endif else begin
        spectra = temporary(Spectra_comp)
        d_spectra = temporary(D_SPECTRA_comp)
    endelse
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; ADDING THE STRAY-LIGHT PROFILE

    IF KEYWORD_SET(slight) THEN BEGIN
        ;Response Functions
        D_SPECTRA=D_SPECTRA*alpha

        ;Magnetic filling factor Response Function
        D_SPECTRA(*,NPAR-1,*)=SPECTRA-slight

        ;Stokes profiles
        SPECTRA=SPECTRA*alpha+slight*(1D0-alpha)

    ENDIF

ENDIF ;END IF ANALYTICAL derivatives

If (D_N EQ 1) OR (D_N EQ 2) then begin  ;NUMERICAL RESPONSE FUNCTIONS

    H=1d-6
    mil_sinrf,PARAM,wl,lmb,yfit,triplet=triplet,slight=slight,filter=filter,$
            AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
    pder=DBLARR(N_ELEMENTS(LMB),N_ELEMENTS(PARAM),4)
    FOR I=0,N_ELEMENTS(PARAM)-1 DO BEGIN
        PARAMN=PARAM
        IF (ABS(PARAM(I)) gt 1e-2) THEN PARAMN(I)=PARAM(I)*(1.+H) ELSE PARAMN(I)=PARAM(I)+H
        ;PARAMN(I)=PARAM(I)*(1.+H)
        mil_sinrf,PARAMn,wl,lmb,yfitD,triplet=triplet,slight=slight,filter=filter,$
                AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
        IF (ABS(PARAM(I)) gt 1e-2) THEN PARAMN(I)=PARAM(I)*(1.-H) ELSE PARAMN(I)=PARAM(I)-H
        ;PARAMN(I)=PARAM(I)*(1.-H)
        mil_sinrf,PARAMn,wl,lmb,yfitI,triplet=triplet,slight=slight,filter=filter,$
                AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
        FOR J=0,3 DO PDER(*,I,J)=(YFITD(*,J)-YFITI(*,J))/(2.*(PARAM(I)-PARAMN(I)))
    ENDFOR
    ; FOR I=0,N_ELEMENTS(PARAM)-1 DO BEGIN
    ; ;    f(a−2h)−8f(a−h)+8f(a+h)−f(a+2h)
    ;     PARAMN=PARAM
    ;     PARAMN(I)=PARAM(I)-2*H
    ;     mil_sinrf,PARAMN,wl,lmb,yfit1,triplet=triplet,slight=slight,filter=filter,$
    ;             AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
    ;     PARAMN(I)=PARAM(I)-H
    ;     mil_sinrf,PARAMN,wl,lmb,yfit2,triplet=triplet,slight=slight,filter=filter,$
    ;             AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
    ;     PARAMN(I)=PARAM(I)+H
    ;     mil_sinrf,PARAMN,wl,lmb,yfit3,triplet=triplet,slight=slight,filter=filter,$
    ;             AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
    ;     PARAMN(I)=PARAM(I)+2*H
    ;     mil_sinrf,PARAMN,wl,lmb,yfit4,triplet=triplet,slight=slight,filter=filter,$
    ;             AC_RATIO=ac_ratio,n_comp=n_comp,crosst=crosst,nlte=nlte
    ;    FOR J=0,3 DO PDER(*,I,J)=(YFIT1(*,J)-8*YFIT2(*,J)+8*YFIT3(*,J)-YFIT4(*,J))/(12.*H)
    ; ENDFOR
ENDIF

If (D_N EQ 2) THEN BEGIN  ;ANALITIC RESPONSE FUNCTIONS
    SPECTRA_N = YFIT
    D_SPECTRA_N = PDER

    label = ["E0","MF","VL","LD","A","GM","AZI","B1","B2","MC","A1","alpha1","A2","alpha2","ALPHA"]
    window,1,xsize=1400,ysize=400
    if keyword_set(saverfs) then begin
        setpscf,filename='rfs-a.ps' ,aspect_ratio=2
        cht = 0.3
    endif else cht = 1
    if keyword_set(nlte) then !p.multi=[0,15,4] else !p.multi=[0,11,4] 
    FOR J=0,nterms*n_comp-1 DO BEGIN
        PLOT,D_SPECTRA(*,j,0),TITLE='I '+label[j],charsize=cht
        OPLOT,D_SPECTRA_N(*,j,0),line=2,thick=2
    ENDFOR
    ;window,2
    ;!p.multi=[0,4,3]
    FOR J=0,nterms*n_comp-1 DO BEGIN
        PLOT,D_SPECTRA(*,j,1),TITLE='Q '+label[j],charsize=cht
        OPLOT,D_SPECTRA_N(*,j,1),line=2,thick=2
    ENDFOR
    ;window,3
    ;!p.multi=[0,4,3]
    FOR J=0,nterms*n_comp-1 DO BEGIN
        PLOT,D_SPECTRA(*,j,2),TITLE='U '+label[j],charsize=cht
        OPLOT,D_SPECTRA_N(*,j,2),line=2,thick=2
    ENDFOR
    ;window,4
    ;!p.multi=[0,4,3]
    FOR J=0,nterms*n_comp-1 DO BEGIN
        PLOT,D_SPECTRA(*,j,3),TITLE='V '+label[j],charsize=cht
        OPLOT,D_SPECTRA_N(*,j,3),line=2,thick=2
    ENDFOR
    if keyword_set(saverfs) then endpscf else pause
ENDIF ELSE If (D_N EQ 1) THEN BEGIN
  SPECTRA = YFIT
  D_SPECTRA = PDER
ENDIF

END
