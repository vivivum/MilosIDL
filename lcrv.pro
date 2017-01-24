function LCVR,theta, delta, muestra=muestra

deltar = delta*!dpi/180d0
thethar = theta*!dpi/180d0
LCVRM = dblarr(4,4)
LCVRM[0] = 1

c2 = cos(2*thethar)
s2 = sin(2*thethar)

if keyword_set(muestra) then begin
  LCVRMS = strarr(4,4);dblarr(4,4)
LCVRMS[0] = '1'
LCVRMS[*] = '0'
LCVRMs[1,1] = 'c2^2d0+s2^2d0*cos(deltar)'
LCVRMs[1,2] = 'c2*s2*(1d0-cos(deltar))'
LCVRMs[1,3] = 's2*sin(deltar)'
LCVRMs[2,1] = 'c2*s2*(1d0-cos(deltar))'
LCVRMs[2,2] = 's2^2d0+c2^2d0*cos(deltar)'
LCVRMs[2,3] = '-c2*sin(deltar)'
LCVRMs[3,1] = '-s2*sin(deltar)'
LCVRMs[3,2] = 'c2*sin(deltar)'
LCVRMs[3,3] = 'cos(deltar)'
print,LCVRMs

endif

LCVRM[1,1] = c2^2d0+s2^2d0*cos(deltar)
LCVRM[1,2] = c2*s2*(1d0-cos(deltar))
LCVRM[1,3] = s2*sin(deltar)
LCVRM[2,1] = c2*s2*(1d0-cos(deltar))
LCVRM[2,2] = s2^2d0+c2^2d0*cos(deltar)
LCVRM[2,3] = -c2*sin(deltar)
LCVRM[3,1] = -s2*sin(deltar)
LCVRM[3,2] = c2*sin(deltar)
LCVRM[3,3] = cos(deltar)


return, LCVRM

end

function LR_Defocus,x,y,ro,xo=xo,yo=yo

if not(keyword_set(xo)) and not(keyword_set(yo)) then begin
xo = 0
yo = 0
endif

theta=atan(x-xo,y-yo)*180./!pi
delta=ro*((x-xo)^2.+(y-yo)^2.)*180./!pi
print,theta,delta
value = lcvr(theta,delta)
return,value

end

function LDiatenuation,x,y,ro,xo=xo,yo=yo,neto=neto

if not(keyword_set(xo)) and not(keyword_set(yo)) then begin
xo = 0
yo = 0
endif

theta=atan(x-xo,y-yo)*180./!pi
delta=ro*((x-xo)^2.+(y-yo)^2.)*180./!pi
print,theta,delta
value = lcvr(theta,delta)

if keyword_set(neto) then begin
 print,1
endif else return,value

end

function MDepolarization,beta

M = fltarr(4,4)
M(0) = 1
M(1,1) = 1/2.+1/2.*beta
M(2,2) = 1/2.+1/2.*beta
M(3,3) = beta

return,M

end

function pmp,delt1,delt2,ventana=ventana

MM = dblarr(4,4)
;delt1 = [225d0,225d0,315d0,315d0]
;delt2 = [234d0,125.26d0,54.74d0,305.26d0]
the1 = 0d0
the2 = 45d0
PL=[[[[1,1,0,0],[1,1,0,0],[0,0,0,0],[0,0,0,0]]]]

if keyword_set(ventana) then begin
  MV = LCVR(the1,ventana)
  endif else begin
  MV = fltarr(4,4)
  for i=0,3 do MV(i,i) = 1.
endelse

For j=0,3 do begin
  M0 = PL##(LCVR(the2,delt2[j])##LCVR(the1,delt1[j])##MV)
  MM[*,j] = M0[*,0]
endfor

return,mm

end
