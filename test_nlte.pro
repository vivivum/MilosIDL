pro test_nlte

;on_error,2
 
;run different tests

Text = ' '
test_milos_1_nlte
PRINT,''
PRINT,'Single fit, check that the fit is ok'
PRINT,''
PRINT,'----------------------------------------'
PRINT,''
PRINT,'Next text: inversion of 1000 random profiles with noise'
PRINT,'Press any key to continue or type exit to stop'
read,Text
if (Text eq 'exit') or (Text eq 'EXIT') then return
PRINT,''

stop

end

pro colors

tvlct,rr,gg,bb,/get
r=[0,255,0,255,0,255,0,0  ,255,255,0  ,175,180]
g=[0,255,0,0,255,0,255,235,0  ,148,126,0  ,180]
b=[0,0,255,0,0,255,255,228,200,0  ,0  ,201,180]
r=[r,r,r,r,r,r,r,r,r,r,r]
g=[g,g,g,g,g,g,g,g,g,g,g]
b=[b,b,b,b,b,b,b,b,b,b,b]
tvlct,[[r],[g],[b]]
return
end

pro test_milos_1_nlte

loadct,0
colors

;************************PHYSICS CONSTANTS (INITIAL)******************
B1=0.2 & B2=0.8
eta0=20d0
Magnet=1500.
GAMMA=60.
AZI=0.
vlos=0.  ;km/s
MACRO=3.
LANDADOPP=0.09
aa=0.05
alfa=0.
A1 = 4.*B1
A2 = 4.*B1
ap1 = 5.
ap2 = 10.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,A1,ap1,A2,ap2,alfa]
;**********************************************************
init_milos,'5250.6',wl

;wavelength axis
Init_landa=Wl(1)-0.4
step=5d0
Points=150.
STEP=STEP/1000d0
axis=Init_landa+findgen(Points)*step

milos, Wl, axis, init_model, y3, /synthesis, /doplot, /nlte
Init_model[12:13] = 0
milos, Wl, axis, init_model, y2, /synthesis, /doplot, /nlte
Init_model[10:11] = 0
milos, Wl, axis, init_model, y1, /synthesis, /doplot, /nlte
milos, Wl, axis, init_model[0:10], y4, /synthesis, /doplot
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,A1,ap1,A2,ap2,alfa]

!p.multi=[0,2,2]
plot,y1(*,0)/y1(0)
oplot,y2(*,0)/y2(0),line=2
oplot,y3(*,0)/y3(0),line=3
oplot,y4(*,0)/y4(0),line=4,thick=2
plot,y1(*,1)/y1(0)
oplot,y2(*,1)/y2(0),line=2
oplot,y3(*,1)/y3(0),line=3
oplot,y4(*,1)/y4(0),line=4,thick=2
plot,y1(*,2)/y1(0)
oplot,y2(*,2)/y2(0),line=2
oplot,y3(*,2)/y3(0),line=3
oplot,y4(*,2)/y4(0),line=4,thick=2
plot,y1(*,3)/y1(0)
oplot,y2(*,3)/y2(0),line=2
oplot,y3(*,3)/y3(0),line=3
oplot,y4(*,3)/y4(0),line=4,thick=2

PRINT,'CHECK NLTE PROFILES'
pause
PRINT,''

S0 = init_model(7)
S1 = init_model(8)
A1 = init_model(10)
ap1 = init_model(11)
A2 = init_model(12)
ap2 = init_model(13)
tau = dindgen(1000)/999.
Sc = S0 + S1*tau + A1*exp(-ap1*tau)
Sl = S0 + S1*tau + A1*exp(-ap1*tau) - A2*exp(-ap2*tau)
plot,tau,Sc
oplot,tau,Sl,line=2

PRINT,'CHECK SOURCE'
pause
PRINT,''

Magnet=500.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,A1,ap1,A2,ap2,alfa]
print,A1,0,A2,ap2
milos, Wl, axis, init_model, y1, rfs=rfs1,numerical=0,/nlte
milos, Wl, axis, init_model, y2, rfs=rfs2,numerical=1,/nlte

; window,0,xsize=2000,ysize=800
; !p.multi=[0,15,4]
; FOR I = 0,3 DO FOR J=0,15-1 DO PLOT,rfs(*,j,i),TITLE='RFS to '+strtrim(string(J),1)

window,1
colores
!p.multi=[0,2,2]
FOR I = 0,3 DO begin
    PLOT,y1(*,i)-y2(*,i),TITLE='Stokes to '+strtrim(string(I),1)
    ;oPLOT,y2(*,i),line=2,thick=2,color=2
endfor
window,2
; !p.multi=[0,4,4]
; FOR I = 0,3 DO begin
;   FOR J=10,13 DO BEGIN
;     PLOT,rfs1(*,j,i),TITLE='RFS to '+strtrim(string(J),1),thick=2,/nodata,$
;        yrange=[max([rfs1(*,j,i),rfs2(*,j,i)]),min([rfs1(*,j,i),rfs2(*,j,i)])]
;     oPLOT,rfs1(*,j,i),thick=2,color=4 ;Analitica GREEN
;     oPLOT,rfs2(*,j,i),line=2,thick=2,color=3 ;numerica red-dashed
;   endfor
; endfor
window,3,xsize=1200,ysize=800
!p.multi=[0,6,4]
FOR I = 0,3 DO begin
  FOR J=1,6 DO BEGIN
    PLOT,rfs1(*,j,i)-rfs2(*,j,i),TITLE='RFS to '+strtrim(string(J),1),thick=2,chars=2.5
  endfor
endfor

STOP
;endpscf

pause

milos, Wl, axis, init_model, y, rfs=rfs, /doplot,/nlte,numerical=2,/saverfs

PRINT,'CHECK RFSs NLTE'
pause

INVS:

OLD=INIT_MODEL
;Init model inversion
B1=0.3 & B2=0.7
eta0=6d0
Magnet=100.
GAMMA=90.
AZI=60.
vlos=0.25  ;km/s
MACRO=3.
LANDADOPP=0.01
aa=0.03
alfa=0.
A1 = 0.8
A2 = 0.8
ap1 = 2.
ap2 = 8.
INIT=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,A1,ap1,A2,ap2,alfa]
;A1 = 4.*0.2
;A2 = 4.*0.2
;ap1 = 5.
;ap2 = 10.

fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,1.,1.,1., 0.]
weight=[1.,1.,1.,1.]

;TIC,/PROFILER
init_model = init
milos, wl, axis, init_model, y, yfit=yfit1,fix=fix,/inversion,miter=150,/doplot,ilambda=10.,/nlte,iter_info=into1,TOPLIM=1e-5,numerical=1
init_model = init
milos, wl, axis, init_model, y, yfit=yfit2,fix=fix,/inversion,miter=150,/doplot,ilambda=10.,/nlte,iter_info=into2,TOPLIM=1e-5
;init_model = init
;milos, wl, axis, init_model, y, yfit=yfit3,fix=fix,/inversion,miter=150,/doplot,ilambda=100.,/nlte,iter_info=into2,TOPLIM=1e-5
;TOC
;PROFILER, /REPORT, /CODE_COVERAGE

window,1
colores
!p.multi=[0,2,2]
FOR I = 0,3 DO begin
    PLOT,y(*,i),TITLE='Stokes to '+strtrim(string(I),1)
    oPLOT,yfit1(*,i),line=3,thick=1,color=4
    oPLOT,yfit2(*,i),line=2,thick=3,color=2
endfor
stop

!p.multi=0
colores
plot,into1.CHISQR,/ylog,thick=2
oplot,into2.chisqr,color=2,thick=2

Print,'Old model:' , OLD
Print,'New model:' , init_model

PRINT,'CHECK FIT'
pause

ftsmio,data,5170,6,xlam=xlam,sdir='data/',dir='data/'

init_milos,'5173',wl
y = fltarr(3000,4)
y [*,0] = data / 10000.

fix=[1.,0.,1.,1.,1.,0.,0.,1.,1.,0., 1.,1.,1.,1., 0.]
INIT_MODEL=[eta0,0.,vlos,landadopp,aa,0.,0.,B1,B2,0.,A1,ap1,A2,ap2,0.]
weight = fltarr(3000,4)
weight[*,0] = 1/y[*,0] 
weight[280:560,0] = 0
weight[670:980,0] = 0
weight[1700:2040,0] = 0
milos, wl, xlam, init_model, y, yfit=yfit,$
  fix=fix,/inversion,miter=100,noise=1.d-3,ilambda=10.,toplim=1e-25,/nlte,weight = weight
  print,init_model
!p.multi=[0,2,2]
plot,xlam,y[*,0]
oplot,xlam,yfit[*,0],line=2,color=3
plot,xlam,(y[*,0]-yfit[*,0])*100.,ytitle="%",yrange=[-3,3]
taU = findgen(1000)/999.
Sc = init_model[0] + init_model[1]*tau + init_model[10]*exp(-init_model[11]*tau)
Sl = Sc - init_model[12]*exp(-init_model[13]*tau)
plot,tau,Sc,/ynoz,xrange=[0,0.5],xstyle=1
oplot,tau,Sl,line=2,color=3

pause

A1 = 0.2
A2 = 0.4
ap1 = 0.001
ap2 = 10

ftsmio,data,8530,24,xlam=xlam,sdir='data/',dir='data/'
xlamn = dindgen(1000)*(xlam(11999)-xlam(0))/1000.+xlam(0)
data = interpol(data,xlam,xlamn)
xlam = xlamn

nlan = n_elements(xlam)
init_milos,'8542',wl
y = fltarr(nlan,4)
y [*,0] = data / 10000.

fix=[1.,0.,1.,1.,1.,0.,0.,1.,1.,0., 1.,1.,1.,1., 0.]
INIT_MODEL=[eta0,0.,vlos,landadopp,aa,0.,0.,B1,B2,0.,A1,ap1,A2,ap2,0.]
weight = fltarr(nlan,4)
weight[*,0] = 1.;/y[*,0] 
milos, wl, xlam, init_model, y, yfit=yfit,$
  fix=fix,/inversion,miter=100,noise=1.d-3,ilambda=10.,toplim=1e-25,/nlte,weight = weight
  print,init_model
!p.multi=[0,2,2]
plot,xlam,y[*,0]
oplot,xlam,yfit[*,0],line=2,color=3
plot,xlam,(y[*,0]-yfit[*,0])*100,ytitle='%',yrange=[-3,3]
taU = findgen(1000)/999.
Sc = init_model[0] + init_model[1]*tau + init_model[10]*exp(-init_model[11]*tau)
Sl = Sc - init_model[12]*exp(-init_model[13]*tau)
plot,tau,Sc,/ynoz,xrange=[0,0.5],xstyle=1
oplot,tau,Sl,line=2,color=3

pause

init_model[1] = 1500
init_model[5] = 45
init_model[6] = 45

milos, Wl, xlam[420:580], init_model, y, rfs=rfs, /doplot,/nlte,numerical=2,/saverfs
milos, Wl, xlam, init_model, y, rfs=rfs, /doplot,/nlte,numerical=2,/saverfs

stop
end
