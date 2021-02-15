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

pro NLTE_noise

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

stop

testp = 10000
b=randomu(seed,testp)*1500.
g=randomu(seed,testp)*180.
a=randomu(seed,testp)*180.
v=randomu(seed,testp)*4.-2.
ch=fltarr(testp)
fm=fltarr(testp,15)
err=fltarr(testp,15)

;Landa inicial
Init_landa=Wl(1)-1.5
;Step in mA
step=30d0
step=step/1000d0
;Samples
Points=100.

eje=Init_landa+Dindgen(Points)*step

fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0., 1.,1.,1.,1., 0.]
init_syn = init_model
init_syn[1] = 10
init_syn[5] = 10
init_syn[6] = 10

weight = fltarr(Points,4)
weight[*,*] = 1.
start = SYSTIME(/SECONDS)
gg = 0
for i=0,testp-1 do begin

    init_syn(1)=b(i)
    init_syn(2)=v(i)
    init_syn(5)=g(i)
    init_syn(6)=a(i)

    milos, Wl, eje, init_syn, y,/synthesis,/nlte

	y = y + randomn(Points,4,100)*1.e-3

    INIT_MODEL=init_syn

    milos, wl, eje, init_model, y, CHISQR=chi,yfit=yfit,$
      sigma=sigma,fix=fix,/inversion,miter=100,/doplot,/quiet,Err=rr,$
      iter_info = iter_info,/nlte,ilambda=10.,toplim=1e-25,weight = weight
	fm(i,*)=init_model
	err(i,*)=rr
	;wait,0.01
    print,'CHISQR VALUE: ',chi, ' Profile: ', i
    ch(i)=chi
    if i eq 0 then iter_info_r = iter_info
    if i gt 0 then iter_info_r = [iter_info_r,iter_info]

endfor

elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

;setpsc,filename='test2.ps'
!p.multi=[0,2,2]
plot,ch,yrange=[0,10],title=palabra(elapsed_time)
plot,err(*,1)
plot,err(*,0)
plot,err(*,5)
;endps
;spawn,'open test2.ps',dd
stop
chl=where(ch gt 1*3.,nl)
print,''
print,'------------------------------------------------------------------------'
print,'Numer of points which the program did not converge: ',nl/1000.*100., ' %'
print,'------------------------------------------------------------------------'

;histogauss,iter_info_r.citer/float(iter_info_r.iter),s1

stop

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
