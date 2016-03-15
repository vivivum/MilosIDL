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

pro test_sfunct

colors
init_milos,'6173',wl

INIT_MODEL=[6.5d0,10d0,0d0,0.03d0,0.03d0,25d0,25d0,0.2d0,0.8d0,0d0,1d0]
INIT = INIT_MODEL

b=randomu(seed,1000)*3000.
g=randomu(seed,1000)*180.
a=randomu(seed,1000)*180.
v=randomu(seed,1000)*4.-2.
s0 = randomu(seed,1000)*0.4+0.1
s1 = randomu(seed,1000)*0.4+0.6
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
s0e=fltarr(1000)
s1e=fltarr(1000)
INIT_SYN = FLTARR(1000,11)

start = SYSTIME(/SECONDS)

for i=0,999 do begin

    INIT_SYN[I,*]=INIT
    init_syn(I,1)=b(i)
    init_syn(I,2)=v(i)
    init_syn(I,5)=g(i)
    init_syn(I,6)=a(i)
    init_syn(I,7)=s0(i)
    init_syn(I,8)=s1(i)

    milos, Wl, eje, reform(init_syn[i,*]), y,/synthesis

    for j=0,0 do begin
	  y = y + randomn(Points,4,100)*1.e-3

    INIT_MODEL = INIT
    INIT_MODEL[7] = y[0]*0.30d0
    INIT_MODEL[8] = y[0]*0.70d0
    s0e[i] = y[0]*0.30d0
    s1e[i] = y[0]*0.70d0
    if b[i] gt 2000. then INIT_MODEL[1] = 2000. ;Fuente error 1

    LIMITS = [[0,1,-1.,1.],[3,1,-0.010,0.010]]  ;Fuente error 2

    milos, wl, eje, init_model, y, chisqr=chisqr,yfit=yfit,$
      sigma=sigma,fix=fix,/inversion,miter=50,/quiet,/doplot,Err=rr,$
      iter_info = iter_info,ilambda=10.,VARLIMITS = LIMITS




      IF chisqr GT 10 THEN STOP
endfor

    if i eq 0 then iter = iter_info else iter = [iter,iter_info]

    fm(i,*)=init_model
    err(i,*)=rr
    print,'CHISQR VALUE: ',chisqr, ' Profile: ', i
    ch(i)=chisqr
endfor

elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

plot,iter[*].PARAMS_STORED[8,*],psym=1

;setpsc,filename='test2.ps'
!p.multi=[0,2,2]
plot,ch,yrange=[0,10];,title=palabra(elapsed_time)
plot,err(*,1)
plot,err(*,0)
plot,err(*,5)
;endps
;spawn,'open test2.ps',dd

stop
end
