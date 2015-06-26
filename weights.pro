;+
; NAME: 
;	WEIGHTS
;
; AUTHOR: 
;	D. Orozco Suarez  	
;						National Astronomical Observatory of Japan, 
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: It calculates the sigma automatically from the Stokes profiles
;
; CATEGORY: Spectropolarimetric fitting
;
; CALLING SEQUENCE: weights, Stokesprof, Sigma, w, noise=noise
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMON BLOCKS: 
;
; NOTES:
;          *NOT USED* IT IS UNDER DEVELOPMENT * 
;     
; CALLED ROUTINES
;
; MODIFICATION HISTORY:
; 	Added documentation: David Orozco Suarez, Feb, 2009
;
;-

PRO weights, Stokesprof, Sigma, w, noise=noise

  Sigma=dblarr(4)
  
  ;Efficiency=[1d0,0.57735027,0.57735027,0.57735027]
  EFFICIENCY=[1d0,0.5d0,0.5d0,0.5d0]
  ;OPCION
  
  valid=where(stokesprof(*,0)*w ne 0.,hay)
  
  Imin=Min(ABS(Stokesprof(valid,0)))
  ;Imax=Max(abs(1d0-PROF(valid,0)))
  Qmax=Max(ABS(Stokesprof(valid,1)))
  Umax=Max(ABS(Stokesprof(valid,2)))
  Vmax=Max(ABS(Stokesprof(valid,3)))
  
  NORM=MIN([Imin,Qmax,Umax,Vmax])
  Sigma(0)=Noise*Imin/Norm/Eficiencia(0)
  ;NORM=MIN([Imax,Qmax,Umax,Vmax])
  ;Sigma(0)=Noise*Imax/norm/Eficiencia(0)
  Sigma(1)=Noise*Qmax/Norm/Eficiencia(1)
  Sigma(2)=Noise*Umax/Norm/Eficiencia(2)
  Sigma(3)=Noise*Vmax/Norm/Eficiencia(3)
  
End
