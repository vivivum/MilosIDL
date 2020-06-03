pro test

;on_error,2

;run different tests

;goto,jump

Text = ' '
test_milos_1
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

test_milos_2
PRINT,''
PRINT,'The chisqr values have to be around 1'
print,'The percentage of points shuld be less than 3% (0.1%?)'
PRINT,''
PRINT,'----------------------------------------'
PRINT,''
PRINT,'Next text: inversion of 1024 real profiles.'
PRINT,'Press any key to continue or type exit to stop'
read,Text
if (Text eq 'exit') or (Text eq 'exit') then return
PRINT,''

test_milos_3,chi
PRINT,''
!p.multi=0
PLOT,CHI,Yrange=[0,50]
PRINT,'----------------------------------------'
print,'Check that all chi values are around 5-10 with a large oscillation up to 20-30'
print,'few of them (6-7) are larger than 50 (noise profiles are also fit)'
print,'----------------------------------------'


test_milos_5

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


pro test_milos_1

loadct,0
colors

;************************PHYSICS CONSTANTS (INITIAL)******************
B1=0.2 & B2=0.8
eta0=12.5d0
Magnet=1200.
GAMMA=20.
AZI=20.
vlos=0.25  ;km/s
MACRO=2.5
LANDADOPP=0.09
aa=0.09
alfa=0.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
;**********************************************************
init_milos,'5250.6',wl

;wavelength axis
Init_landa=Wl(1)-0.4
step=5d0
Points=150.
STEP=STEP/1000d0
axis=Init_landa+Dindgen(Points)*step

milos, Wl, axis, init_model, y, /synthesis, /doplot

WAIT,1

milos, Wl, axis, init_model, y, rfs=rfs, /doplot

OLD=INIT_MODEL
;Init model inversion
B1=0.3 & B2=0.7
eta0=6d0
Magnet=100.
GAMMA=90.
AZI=60.
vlos=0.25  ;km/s
MACRO=2.5
LANDADOPP=0.01
aa=0.03
alfa=0.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]

fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
weight=[1.,10.,10.,1.]


INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
milos, wl, axis, init_model, y, yfit=yfit,$
  fix=fix,/inversion,miter=100,noise=1.d-3,/doplot,ilambda=10.,toplim=1e-25

Print,'Old model:' , OLD
Print,'New model:' , init_model

end

pro test_milos_2

init_milos,'5250.6',wl

B1=0.2 & B2=0.8		;cociente B1/Bo terminos de la funcion de Planck b0+b1t
eta0=6.5d0	        ;cociente entre coef. abs. de la linea y el continuo
Magnet=1200.               ;campo magn�tico que no es cero por problemas con svd
GAMMA=20.
AZI=20.
vlos=0.25  ;km/s
LANDADOPP=0.03
aa=0.03
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,0d0,1d0]
;**********************************************************
INIT_SYN=INIT_MODEL
INIT=INIT_MODEL

b=randomu(seed,1000)*3000.
g=randomu(seed,1000)*180.
a=randomu(seed,1000)*180.
v=randomu(seed,1000)*4.-2.

;Landa inicial
Init_landa=Wl(1)-0.4
;Step in mA
step=10d0
;Samples
Points=75.

STEP=STEP/1000d0
fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
sigma=[1.,1.,1.,1.]/1000.

eje=Init_landa+Dindgen(Points)*step

ch=fltarr(1000)

fm=fltarr(1000,11)
err=fltarr(1000,11)

start = SYSTIME(/SECONDS)
gg = 0
for i=0,999 do begin

    init_syn(1)=b(i)
    init_syn(2)=v(i)
    init_syn(5)=g(i)
    init_syn(6)=a(i)

    milos, Wl, eje, init_syn, y,/synthesis

	y = y + randomn(Points,4,100)*1.e-3

    INIT_MODEL=INIT

    milos, wl, eje, init_model, y, chisqr=chisqr,yfit=yfit,$
      sigma=sigma,fix=fix,/inversion,miter=100,/doplot,/quiet,Err=rr,$
      iter_info = iter_info
	fm(i,*)=init_model
	err(i,*)=rr
	;wait,0.01
    print,'CHISQR VALUE: ',chisqr, ' Profile: ', i
    ch(i)=chisqr
    if i eq 0 then iter_info_r = iter_info
    if i gt 0 then iter_info_r = [iter_info_r,iter_info]

endfor

elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

setpsc,filename='test2.ps'
!p.multi=[0,2,2]
plot,ch,yrange=[0,10],title=palabra(elapsed_time)
plot,err(*,1)
plot,err(*,0)
plot,err(*,5)
endps
spawn,'open test2.ps',dd

chl=where(ch gt 1*3.,nl)
print,''
print,'------------------------------------------------------------------------'
print,'Numer of points which the program did not converge: ',nl/1000.*100., ' %'
print,'------------------------------------------------------------------------'

histogauss,iter_info_r.citer/float(iter_info_r.iter),s1

stop
end

pro test_milos_3,chi

init_milos,'63016302',wl

restore,'data/wavelength_axis.sav'
laxis=landas(5:94)

S0=0.2
S1=0.8
eta0=20.d0
Magnet=300.
GM=15.
AZ=15.
vlos=0.1  ;km/s
MACRO=0d0
LANDADOPP=0.029
aa=0.8
ffactor=1d0

INIT_MODEL=double([eta0,magnet,vlos,landadopp,aa,gm,az,S0,S1,macro,ffactor])

fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
sigma=[0.003,0.00035,0.00035,0.0003]

weight=[1.,1.,1.,1.]

profile=dblarr(90,4)

file='data/SP4D20070227_002010.0C.fits'

result=dblarr(1024,11)
error=dblarr(1024,11)
chi=fltarr(1024)
chii=fltarr(1024,4)

!p.multi=[0,2,2]
colors

readl1_sbsp,file,tmp_data,hdr

average=total(tmp_data(*,*,0),2)/1024.
conti=mean(average(1:6))

start = SYSTIME(/SECONDS)
for i=0,1023 do begin
        profile(*,0)=reform(tmp_data(5:94,i,0))/conti
        profile(*,1)=reform(tmp_data(5:94,i,1))/conti
        profile(*,2)=reform(tmp_data(5:94,i,2))/conti
        profile(*,3)=reform(tmp_data(5:94,i,3))/conti
        init=init_model

            milos, wl, laxis, init, profile, yfit=yfit,fix=fix,/inversion,$
              chisqr=chisqr,miter=400,err=err,filter=30d0,$
              sigma=sigma,/doplot,weight=weight,/quiet,getshi=getshi,/ac_ratio

            result(i,*)=init
            error(i,*)=err
            chi(i)=chisqr
            chii(i,*)=getshi
            print,'CHISQR VALUE: ',chisqr, ' Profile: ', i
endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

setpsc,filename='test.ps'
!p.multi=0
plot,chi,yrange=[1,30],title=palabra(elapsed_time)
endps
spawn,'open test.ps',dd

end

pro test_milos_4

init_milos,'6301',wl

B1=0.2 & B2=0.8		;cociente B1/Bo terminos de la funcion de Planck b0+b1t
eta0=6.5d0	        ;cociente entre coef. abs. de la linea y el continuo
Magnet=1200.               ;campo magn�tico que no es cero por problemas con svd
GAMMA=20.
AZI=0.
vlos=0.25  ;km/s
LANDADOPP=0.03
aa=0.03
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,0d0,1d0]
;**********************************************************
INIT_SYN=INIT_MODEL
INIT=INIT_MODEL

b=findgen(100)*10.
g=findgen(180)
a=0.
v=0.

;Landa inicial
Init_landa=Wl(1)-0.4
;Step in mA
step=10d0
;Samples
Points=75.

STEP=STEP/1000d0

eje=Init_landa+Dindgen(Points)*step

pr = fltarr(75,4,100,180)
prn = fltarr(75,4,100,180)

for i=0,99 do begin
for j=0,179 do begin

    init_syn(1)=b(i)
    init_syn(5)=g(j)

    milos, Wl, eje, init_syn, y,/synthesis
    pr(*,*,i,j) = y

endfor
    print,' Profile: ', i
endfor

ampV = fltarr(100,180)
ampQ = fltarr(100,180)

for i=0,99 do for j=0,179 do ampV(i,j) = max(abs(pr(*,3,i,j)))
for i=0,99 do for j=0,179 do ampQ(i,j) = max( [abs(pr(*,1,i,j)),abs(pr(*,2,i,j))])

STOP
end

pro test_milos_5

print,'Two COMPONENTS'

loadct,0
colors

;************************PHYSICS CONSTANTS (INITIAL)******************
B1=0.2 & B2=0.8
eta0=12.5d0
Magnet=200.
GAMMA=20.
AZI=20.
vlos=0.25  ;km/s
MACRO=0.;2.5
LANDADOPP=0.06
aa=0.09
alfa=0.
MODEL1=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
MODEL2=[eta0,100,-6,landadopp,aa,130.,azi,B1,B2,macro,0.3]
INIT_MODEL=[MODEL1,MODEL2]
;**********************************************************
init_milos,'5250.2',wl

;wavelength axis
Init_landa=Wl(1)-0.4
step=5d0
Points=150.
STEP=STEP/1000d0
axis=Init_landa+Dindgen(Points)*step

milos, Wl, axis, init_model, y, /synthesis, /doplot,n_comp=2

WAIT,1

milos, Wl, axis, init_model, y, rfs=rfs, /doplot,n_comp=2

WAIT,1
OLD=INIT_MODEL
;Init model inversion
B1=0.3 & B2=0.7
eta0=6d0
Magnet=100.
GAMMA=40.
AZI=60.
vlos=0.25  ;km/s
MACRO=0.
LANDADOPP=0.06
aa=0.09
alfa=0.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
MODEL1=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
MODEL2=[eta0,200,-5,landadopp,aa,170.,azi,B1,B2,macro,0.6]
INIT_MODEL=[MODEL1,MODEL2]

fix1=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
fix2=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.]
fix=[fix1,fix2]
weight=[1.,1.,1.,1.]

milos, wl, axis, init_model, y, yfit=yfit,$
  fix=fix,/inversion,miter=30,noise=1.d-6,/doplot,n_comp=2,/numerical

Print,'Old model:' , OLD
Print,'New model:' , init_model

stop
end

pro test_milos_6,chi

init_milos,'63016302',wl

restore,'data/wavelength_axis.sav'
laxis=landas(5:94)

S0=0.2
S1=0.8
eta0=20.d0
Magnet=300.
GM=15.
AZ=15.
vlos=0.1  ;km/s
MACRO=0d0
LANDADOPP=0.029
aa=0.8
ffactor=0.3d0

MODEL1=[eta0,magnet,vlos,landadopp,aa,gm,az,S0,S1,macro,ffactor]*1d0
INIT_MODEL=model1
fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.]

MODEL1=[eta0,magnet,vlos,landadopp,aa,gm,az,S0,S1,macro,ffactor]*1d0
MODEL2=[eta0,0,0.,landadopp,aa,0.,0.,s0,s1,macro,0.3]*1d0
INIT_MODEL=[model1,model2]
fix1=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.]
fix2=[1.,0.,1.,1.,1.,0.,0.,1.,1.,0.,1.]
fix=[fix1,fix2]

goto,jump

MODEL1=[eta0,magnet,vlos,landadopp,aa,gm,az,S0,S1,macro,ffactor]*1d0;
MODEL2=[eta0,500,2.,landadopp,aa,170.,30.,s0,s1,macro,0.3]*1d0
INIT_MODEL=[model1,model2]
fix1=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
fix2=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.]
fix=[fix1,fix2]

MODEL1=[eta0,magnet,vlos,landadopp,aa,gm,az,S0,S1,macro,ffactor]*1d0
MODEL2=[eta0,500,2.,landadopp,aa,170.,30.,s0,s1,macro,0.3]*1d0
MODEL3=[eta0,0,0.,landadopp,aa,0.,0.,s0,s1,macro,0.3]*1d0
INIT_MODEL=[model1,model2,model3]
fix1=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
fix2=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.]
fix3=[1.,0.,1.,1.,1.,0.,0.,1.,1.,0.,1.]
fix=[fix1,fix2,fix3]
jump:

sigma=[0.003,0.00035,0.00035,0.0003]

weight=[1.,1.,1.,1.]

profile=dblarr(90,4)

file='data/SP4D20070227_002010.0C.fits'

result=dblarr(1024,11*2)
error=dblarr(1024,11*2)
chi=fltarr(1024)
chii=fltarr(1024,4)

!p.multi=[0,2,2]
colors

readl1_sbsp,file,tmp_data,hdr

average=total(tmp_data(*,*,0),2)/1024.
conti=mean(average(1:6))

for i=0,1023 do begin
        profile(*,0)=reform(tmp_data(5:94,i,0))/conti
        profile(*,1)=reform(tmp_data(5:94,i,1))/conti
        profile(*,2)=reform(tmp_data(5:94,i,2))/conti
        profile(*,3)=reform(tmp_data(5:94,i,3))/conti
        init=init_model
        diff=reform(total(tmp_data(5:94,0>(i-5):(i+5)<1023,*),2))/( (i+5)<1023 - 0>(i-5) +1.)/conti
		diff(*,1:3) = 0.

            milos, wl, laxis, init, profile, yfit=yfit,fix=fix,/inversion,$
              chisqr=chisqr,miter=100,err=err,filter=30d0,slight=diff,$
              sigma=sigma,/doplot,weight=weight,getshi=getshi,/ac_ratio,n_comp=2,/quiet;,/numerical;,/quiet

            result(i,*)=init
            error(i,*)=err
            chi(i)=chisqr
            chii(i,*)=getshi
            print,'CHISQR VALUE: ',chisqr, ' Profile: ', i
endfor

end

pro test_milos_7

  RESOLVE_ROUTINE, ['milos','mil_svd','mil_sinrf','me_der',$
  'lm_mils','check_param','covarm','create_nc',$
  'fvoigt','init_milos','quanten','weights',$
  'weights_init']
  RESOLVE_routine, ['filtro','macrotur'], /IS_FUNCTION

  Profiler, /SYSTEM & Profiler

loadct,0
colors

;************************PHYSICS CONSTANTS (INITIAL)******************
B1=0.2 & B2=0.8
eta0=12.5d0
Magnet=1200.
GAMMA=20.
AZI=20.
vlos=0.25  ;km/s
MACRO=0;2.5
LANDADOPP=0.09
aa=0.09
alfa=0.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
;**********************************************************
init_milos,'5250.6',wl

;wavelength axis
Init_landa=Wl(1)-0.4
step=5d0
Points=150.
STEP=STEP/1000d0
axis=Init_landa+Dindgen(Points)*step

milos, Wl, axis, init_model, y, /synthesis, /doplot

WAIT,1

milos, Wl, axis, init_model, y, rfs=rfs, /doplot

OLD=INIT_MODEL
;Init model inversion
B1=0.3 & B2=0.7
eta0=6d0
Magnet=100.
GAMMA=90.
AZI=60.
vlos=2.25  ;km/s
MACRO=0.;1.8
LANDADOPP=0.01
aa=0.03
alfa=0.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]

fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
weight=[1.,1.,1.,1.]

INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
milos, wl, axis, init_model, y, yfit=yfit,ilambda = 10.,$
  fix=fix,/inversion,miter=100,noise=1.d-6,/doplot,iter_info = iter_info,/use_svd_cordic

Print,'Old model:' , OLD
Print,'New model:' , init_model

Profiler, /REPORT
Profiler,/clear,/system
Profiler,/clear,/reset
stop
end


pro test_milos_cordic

init_milos,'5250.6',wl

B1=0.2 & B2=0.8		;cociente B1/Bo terminos de la funcion de Planck b0+b1t
eta0=6.5d0	        ;cociente entre coef. abs. de la linea y el continuo
Magnet=1200.               ;campo magn�tico que no es cero por problemas con svd
GAMMA=20.
AZI=20.
vlos=0.25  ;km/s
LANDADOPP=0.03
aa=0.03
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,0d0,1d0]
;**********************************************************
INIT_SYN=INIT_MODEL
INIT=INIT_MODEL

b=randomu(seed,1000)*3000.
g=randomu(seed,1000)*180.
a=randomu(seed,1000)*180.
v=randomu(seed,1000)*4.-2.

;Landa inicial
Init_landa=Wl(1)-0.4
;Step in mA
step=10d0
;Samples
Points=75.

STEP=STEP/1000d0
fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
sigma=[1.,1.,1.,1.]/1000.

eje=Init_landa+Dindgen(Points)*step

ch=fltarr(1000,2)

fm=fltarr(1000,11,2)
err=fltarr(1000,11)

start = SYSTIME(/SECONDS)
for i=0,99 do begin

    init_syn(1)=b(i)
    init_syn(2)=v(i)
    init_syn(5)=g(i)
    init_syn(6)=a(i)

    milos, Wl, eje, init_syn, y,/synthesis

	y = y + randomn(Points,4,100)*1.e-3

    INIT_MODEL=INIT

    milos, wl, eje, init_model, y, chisqr=chisqr,yfit=yfit,$
      sigma=sigma,fix=fix,/inversion,miter=100,/quiet,/doplot,Err=rr,$
      iter_info = iter_info
	fm(i,*,0)=init_model
	err(i,*)=rr
	;wait,0.01
    print,'CHISQR VALUE: ',chisqr, ' Profile: ', i
    ch(i,0)=chisqr

endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

start = SYSTIME(/SECONDS)
for i=0,99 do begin

    init_syn(1)=b(i)
    init_syn(2)=v(i)
    init_syn(5)=g(i)
    init_syn(6)=a(i)

    milos, Wl, eje, init_syn, y,/synthesis

	y = y + randomn(Points,4,100)*1.e-3

    INIT_MODEL=INIT

    milos, wl, eje, init_model, y, chisqr=chisqr,yfit=yfit,$
      sigma=sigma,fix=fix,/inversion,miter=100,/quiet,/doplot,Err=rr,$
      iter_info = iter_info,use_svd_cordic=2
	fm(i,*,1)=init_model
	err(i,*)=rr
	;wait,0.01
    print,'CHISQR VALUE: ',chisqr, ' Profile: ', i
    ch(i,1)=chisqr

endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

STOP

end

pro test_crosst

loadct,0
colors

;************************PHYSICS CONSTANTS (INITIAL)******************
B1=0.2 & B2=0.8
eta0=12.5d0
Magnet=1200.
GAMMA=20.
AZI=20.
vlos=0.25  ;km/s
MACRO=2.5
LANDADOPP=0.09
aa=0.09
alfa=0.
e1 = 0.001
e2 = 0.005
e3 = 0.005
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa,e1,e2,e3]

;**********************************************************
init_milos,'5250.6',wl

;wavelength axis
Init_landa=Wl(1)-0.4
step=5d0
Points=150.
STEP=STEP/1000d0
axis=Init_landa+Dindgen(Points)*step

milos, Wl, axis, init_model, y, /synthesis, /doplot,/crosst

WAIT,1

milos, Wl, axis, init_model, y, rfs=rfs, /doplot,/crosst

OLD=INIT_MODEL
;Init model inversion
B1=0.3 & B2=0.7
eta0=12d0
Magnet=100.
GAMMA=90.
AZI=60.
vlos=0.25  ;km/s
MACRO=2.5
LANDADOPP=0.01
aa=0.03
alfa=0.
INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]

fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,1.,1.]
weight=[1.,10.,10.,1.]


INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa,0.002,0.002,0.002]
milos, wl, axis, init_model, y, yfit=yfit,$
  fix=fix,/inversion,miter=100,noise=1.d-6,/doplot,ilambda=10.,toplim=1e-25,/crosst

  Print,'New model:' , init_model
Print,'Old model:' , OLD

end
