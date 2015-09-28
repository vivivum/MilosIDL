;+
; NAME:
;	CHECK_PARAM
;
; AUTHOR:
;	D. Orozco Suarez
;						National Astronomical Observatory of Japan, 
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: Check the new model parameters' limits
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: CHECK_PARAM,P_M,PLIMITS=plimits
;
; INPUTS:
;           P_M: 11-element float or double array containing the new model parameters
;
; KEYWORD PARAMETERS:
;   PLIMITS: Limitations on the model parameter changes can be imposed throught this keywork
;            Array of structures. Dimension (11).
;            	PLIMITS.SET: scalar value indicating whether the parameter limitation in on/off
;                            1-> on, 0-> off
;               PLIMITS.LIMITS: two element array with the min/max limits for the parameter
;               DEFAULT LIMITS ARE:
;               PLIMITS[*].SET = 1  ; Activate all limits
;               PLIMITS[0].LIMITS = [1d0,2500d0]  ; Eta0
;               PLIMITS[1].LIMITS = [0d0,4500d0]  ; Magnetic field
;               PLIMITS[2].LIMITS = [-20d0,20d0]  ; Velocity (km/s)
;               PLIMITS[3].LIMITS = [1d-4,6d-1]  ; Doppler width
;               PLIMITS[4].LIMITS = [1d-4,1d1]  ; Damping
;               PLIMITS[5].LIMITS = [0d0,180d0]  ; Inclination
;               PLIMITS[6].LIMITS = [0d0,180d0]  ; Azimuth
;               PLIMITS[7].LIMITS = [1d-4,1d1]  ; S0
;               PLIMITS[8].LIMITS = [1d-4,1d1]  ; S1
;               PLIMITS[9].LIMITS = [0d0,4d0]  ; Macroturbulence
;               PLIMITS[10].LIMITS = [0d0,1d0]  ; Magnetic filling factor
;
; OUTPUTS:
;           P_M (overwritten)
;
; COMMON BLOCKS: none
;
; NOTES:
;
; CALLED ROUTINES:
;
; MODIFICATION HISTORY:
; 	Added documentation: David Oozco Suarez (DOS), Feb, 2009
;   Added plimits keyword to allow user-interaction with CHECK_PARAM: 27 Feb, 2009 DOS
;
;-


PRO  CHECK_PARAM,P_M,PLIMITS=plimits,n_comp=n_comp

;ch=[0,0,0,0,0,0,0,0,0,0]
;lb=['eta0','magnet','vlos','landadopp','aa','gamma','azi','B1','B2','macro','alpha']

  ;Magnetic field
  if (P_M(1) lt 0) and (PLIMITS(1).SET EQ 1) then begin
	P_M(1) = (-1d0)*P_M(1)
	P_M(5) = 180d0 - P_M(5)
;    ch(1)=1
  endif
  ;Inclination
  If PLIMITS(5).SET EQ 1 THEN BEGIN
  If P_M(5) LT PLIMITS[5].LIMITS(0)  THEN begin
    P_M(5) = -P_M(5)
;    ch(5)=1
  endif
  If P_M(5) GT PLIMITS[5].LIMITS(1)  THEN begin
    P_M(5) = 360d0-P_M(5)
;    ch(5)=1
  endif
  ENDIF

  ;azimuth
  If PLIMITS(6).SET EQ 1 THEN BEGIN
  If P_M(6) LT 0. THEN  begin
    P_M(6) = 180d0 + P_M(6)
;    ch(6)=1
  endif
  If P_M(6) GT PLIMITS[6].LIMITS(1) THEN  begin
    P_M(6) = P_M(6) - 180d0
;    ch(6)=1
  endif
  ENDIF

  for i=0,4 do IF PLIMITS(i).SET EQ 1 THEN $
	P_M(i) = PLIMITS(i).LIMITS(0) > P_M(i) < PLIMITS(i).LIMITS(1)
  for i=7,10 do IF PLIMITS(i).SET EQ 1 THEN $
	P_M(i) = PLIMITS(i).LIMITS(0) > P_M(i) < PLIMITS(i).LIMITS(1)

if n_comp ge 2 then begin

for j = 1,n_comp do begin

  ;Magnetic field
  if (P_M(1+11*(j-1)) lt 0) and (PLIMITS(1+11*(j-1)).SET EQ 1) then begin
	P_M(1+11*(j-1)) = (-1d0)*P_M(1+11*(j-1))
	P_M(5+11*(j-1)) = 180d0 - P_M(5+11*(j-1))
;    ch(1)=1
  endif
  ;Inclination
  If PLIMITS(5+11*(j-1)).SET EQ 1 THEN BEGIN
  If P_M(5+11*(j-1)) LT PLIMITS[5+11*(j-1)].LIMITS(0)  THEN begin
    P_M(5+11*(j-1)) = -P_M(5+11*(j-1))
;    ch(5)=1
  endif
  If P_M(5+11*(j-1)) GT PLIMITS[5+11*(j-1)].LIMITS(1)  THEN begin
    P_M(5+11*(j-1)) = 360d0-P_M(5+11*(j-1))
;    ch(5)=1
  endif
  ENDIF

  ;azimuth
  If PLIMITS(6+11*(j-1)).SET EQ 1 THEN BEGIN
  If P_M(6+11*(j-1)) LT 0. THEN  begin
    P_M(6+11*(j-1)) = 180d0 + P_M(6+11*(j-1))
;    ch(6)=1
  endif
  If P_M(6+11*(j-1)) GT PLIMITS[6+11*(j-1)].LIMITS(1) THEN  begin
    P_M(6+11*(j-1)) = P_M(6+11*(j-1)) - 180d0
;    ch(6)=1
  endif
ENDIF

  for i=0,4 do IF PLIMITS(i+11*(j-1)).SET EQ 1 THEN $
	P_M(i+11*(j-1)) = PLIMITS(i+11*(j-1)).LIMITS(0) > P_M(i+11*(j-1)) < PLIMITS(i+11*(j-1)).LIMITS(1)
  for i=7,10 do IF PLIMITS(i+11*(j-1)).SET EQ 1 THEN $
	P_M(i+11*(j-1)) = PLIMITS(i+11*(j-1)).LIMITS(0) > P_M(i+11*(j-1)) < PLIMITS(i+11*(j-1)).LIMITS(1)

endfor

endif

END
