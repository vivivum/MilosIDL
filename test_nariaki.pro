pro test_nariaki,chi

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
chirmean = 0.
start = SYSTIME(/SECONDS)
dorepeat = 1
init=init_model

for i=0,1023 do begin
        profile(*,0)=reform(tmp_data(5:94,i,0))/conti
        profile(*,1)=reform(tmp_data(5:94,i,1))/conti
        profile(*,2)=reform(tmp_data(5:94,i,2))/conti
        profile(*,3)=reform(tmp_data(5:94,i,3))/conti

            milos, wl, laxis, init, profile, yfit=yfit,fix=fix,/inversion,$
              chisqr=chisqr,miter=100,err=err,filter=30d0,$
              sigma=sigma,/doplot,weight=weight,/quiet,getshi=getshi,/ac_ratio

            result(i,*)=init
            error(i,*)=err
            chi(i)=chisqr
            chii(i,*)=getshi
            chirmean = (chirmean + chisqr )/2.
            print,'CHISQR VALUE: ',chisqr, ' Profile: ', i, chirmean
            if (chisqr gt 1.2*chirmean) and (dorepeat ne 0) then begin
                i = i - 1
                dorepeat = 0
                init [5] = 134.
                init [6] = 144.                 
                init [1] = 1300. 
                init [2] = -0.1                
                print,'repeat'
            endif else begin
                dorepeat = 1
                init = init_model
            endelse
endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

plot,chi,yrange=[1,30],title=palabra(elapsed_time)

stop
end
