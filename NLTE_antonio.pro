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

pro NLTE_antonio

loadct,0
colors

;************************PHYSICS CONSTANTS (INITIAL)******************
B1=0.2 & B2=0.8
eta0=20d0
Magnet=1500.
GAMMA=60.
AZI=0.
vlos=0.  ;km/s
MACRO=0.
LANDADOPP=0.09
aa=0.05
alfa=0.
A1 = 4.*B1
A2 = 4.*B1
ap1 = 5.
ap2 = 10.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,A1,ap1,A2,ap2,alfa]
;**********************************************************
init_milos,'5173',wl  ;INICIALIZA MILOS CON MgI

;Hacemos una sintesis de prueba

;wavelength axis
Init_landa=Wl(1)-0.4
step=5d0
Points=150.
STEP=STEP/1000d0
axis=Init_landa+findgen(Points)*step

milos, Wl, axis, init_model, y, /synthesis, /nlte

PRINT,'CHECK NLTE PROFILES'

;plot profile and source
!p.multi=[0,2,3]
plot,axis,y[*,0],title='Stokes I',xtitle='Wavelength [A]',ytitle='Intensity'
plot,axis,y[*,1],title='Stokes Q',xtitle='Wavelength [A]',ytitle='Stokes Q/Ic'
plot,axis,y[*,2],title='Stokes U',xtitle='Wavelength [A]',ytitle='Stokes U/Ic'
plot,axis,y[*,3],title='Stokes V',xtitle='Wavelength [A]',ytitle='Stokes V/Ic'

S0 = init_model(7)
S1 = init_model(8)
A1 = init_model(10)
ap1 = init_model(11)
A2 = init_model(12)
ap2 = init_model(13)
tau = dindgen(1000)/999.
Sc = S0 + S1*tau + A1*exp(-ap1*tau)
Sl = S0 + S1*tau + A1*exp(-ap1*tau) - A2*exp(-ap2*tau)
plot,tau,Sc,title='Source function',xtitle='opacity'
oplot,tau,Sl,line=2

PRINT,'CHECK SOURCE'
pause
PRINT,'NEXT DEVUELVE LAS FUNCIONES DE RESPUESTA EN rfs'

Magnet=500.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,A1,ap1,A2,ap2,alfa]
milos, Wl, axis, init_model, y, rfs=rfs1,/nlte

window,3,xsize=1200,ysize=800
!p.multi=[0,14,4]
FOR I = 0,3 DO begin
  FOR J=0,12 DO BEGIN
    PLOT,rfs1(*,j,i),TITLE='RFS to '+strtrim(string(J),1),thick=2,chars=0.5
  endfor
endfor

pause

;INVERSION
ORIGINAL_model = INIT_MODEL
;Init model inversion
B1=0.3 & B2=0.7
eta0=6d0
Magnet=100.
GAMMA=90.
AZI=60.
vlos=0.25  ;km/s
MACRO=0.
LANDADOPP=0.01
aa=0.03
alfa=0.
A1 = 0.8
A2 = 0.8
ap1 = 2.
ap2 = 8.
INIT=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,A1,ap1,A2,ap2,alfa]

fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,1.,1.,1.,0.] ;que ajusta y que no
weight=[1.,1.,1.,1.]

init_model = init
milos, wl, axis, init_model, y, yfit=yfit1,fix=fix,/inversion,miter=150,/doplot,ilambda=10.,/nlte,iter_info=into2,TOPLIM=1e-5

window,1
colores
!p.multi=[0,2,2]
FOR I = 0,3 DO begin
    PLOT,y(*,i),TITLE='Stokes to '+strtrim(string(I),1)
    oPLOT,yfit1(*,i),line=3,thick=1,color=4
endfor

pause

;PRUEBA CON 517 FTS

ftsmio,data,5170,6,xlam=xlam,sdir='data/',dir='data/'

init_milos,'5173',wl
y = fltarr(3000,4)
y [*,0] = data / 10000.

fix=[1.,0.,1.,1.,1.,0.,0.,1.,1.,0., 1.,1.,1.,1., 0.] ;no ajusta el campo y por eso B,gamma y azi == 0
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

;PRUEBA CON 8530 FTS

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

stop
end
