pro invf

;Tomamos el eje de longitudes de onda
;restore,'../deep-mode/wavelength_axis.sav'
restore,'/Users/orozco/Desktop/CalculosPaper/wavelength_axis.sav'
laxis=landas(5:94)
restore,'average.sav'
contI=mean(average[0:5])

;************************PHISYC CONSTANTS (INITIAL)********
B0=0.2 & B1=0.8         ;cociente B1/Bo terminos de la funcion de Planck b0+b1t
eta0=7.5d0      ;cociente entre coef. abs. de la linea y el continuo
Magnet=300.
GM=15.
AZ=15.
vlos=0.1  ;km/s
MACRO=0.
LANDADOPP=0.029
aa=0.78
model1=[eta0,300.,0.1,landadopp,aa,30,30,B0,B1,macro,1]*1d0
model2=[7.5d0,0.,-0.1,0.06d0,0.78d0,0.,0.,B0,B1,macro,0.8]*1d0
init_model=[model1,model2]
;************************PHISYC CONSTANTS (INITIAL)********

init_milos,'63016302',wl

dir='/Users/orozco/Desktop/CalculosPaper/center'
;dir='../data/sot/Center-to-limb/center'
files=file_search(dir+'/SM*.fits')
n_files=n_elements(files)

;readl1_sbsp,files,tmp_data,hdr 
;stop
;average=fltarr(112)
;np=0.
;for i=0,n_files-1 do begin    &$
;    readl1_sbsp,files(i),dat,hdr    &$
;    for j=0,1023 do begin    &$
;        if max(abs(dat(*,j,3))) lt 60 then begin    &$
;            average=average+dat(*,j,0)    &$
;            np=np+1.    &$
;        endif    &$
;    endfor    &$
;print,i   &$
;endfor    
;average=average/np
;save,filename='average_dm.sav',average

;Los resultados seran 
result=dblarr(n_files,1024,11*2)
error=dblarr(n_files,1024,11*2)
chi=fltarr(n_files,1024,4)

fix1=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
fix2=[1.,0.,1.,1.,1.,0.,0.,1.,1.,0.,1.]
fix=[fix1,fix2]
sig=[0.0003,0.0003,0.0003,0.0003]
weight=[1.,20.,20.,10.]
;weight=[1.,1.,1.,1.]

;LIMITE INVERSION 6e3

profile=dblarr(90,4)
diff=dblarr(90,4,1024)
;A VER COMO CALCULO YO LA LUZ DIFUSA....


;for i=0,n_files-1 do begin      ;Bucle principal
;progressbar,/init
for i=0,n_files-1 do begin      ;Bucle principal
   
    readl1_sbsp,files(i),tmp_data,hdr 
    ;tmp_data=float(readfits(files(i),hdr))

    ;for j=0,1023 do begin
    ;    ii=(j-3)>0
    ;    ff=(j+3)<1023
    ;    diff(*,*,j)=reform(total(tmp_data(5:94,ii:ff,*),2))/float(ff-ii)
    ;endfor 
    ;diff=diff/conti
    ;diff(*,1:3,*)=0.
    for j=0,1023 do begin
        
        profile(*,0)=reform(tmp_data(5:94,j,0))/conti 
        profile(*,1)=reform(tmp_data(5:94,j,1))/conti 
        profile(*,2)=reform(tmp_data(5:94,j,2))/conti 
        profile(*,3)=reform(tmp_data(5:94,j,3))/conti 
              
        init=init_model
        ;init(10)=0.95
        
        ;milos, wl, laxis, init, profile, yfit=yfit,fix=fix,/inversion,$
        ;  ilambda=1.,chisqr=chisqr,miter=300,err=err,toplim=1.e-30,$
        ;  sigma=sig,weight=weight,filter=25d0,/quiet,slight=diff(*,*,j),getshi=chisqr_new;,/doplot
        milos, wl, laxis, init, profile, yfit=yfit,fix=fix,/inversion,$
          chisqr=chisqr,miter=100,err=err,filter=30d0,toplim=1.e-30,$
          sigma=sig,weight=weight,/quiet,getshi=chisqr_new,n_comp=2,/doplot,/ac_ratio,VARLIMITS = [11, 1, -5d0, 5d0];,/numerical

        result(i,j,*)=init
        error(i,j,*)=err
        chi(i,j,*)=chisqr_new

        print,chisqr,'px: ',i,j 

    endfor 
;        progressbar,progress=[0,n_files-1,i],title='High S/N '
endfor

save,filename='deep_mode_result_center_2c2.sav',result,chi,error,/compress ;pesos
 ;progressbar,/destroy

STOP

;QUIERO CALCULAR LOS FLUJOS.....

campo=result(0:25,*,1)
incli=result(0:25,*,5)
alfa=result(0:25,*,10)

PRINT,'FLUJO MEDIO TOTAL: ',mean(campo * alfa)
PRINT,'FLUJO MEDIO LONGITUDINAL: ',mean(campo*cos(incli*!pi/180.) * alfa)
PRINT,'FLUJO MEDIO TRANSVERSAL: ',mean(campo*sin(incli*!pi/180.) * alfa)
PRINT,'FLUJO MEDIO LONGITUDINAL (SIN SIGNO): ',mean(abs(campo*cos(incli*!pi/180.)) * alfa)

;LA COSA ES QUE HAY RUIDO, QUE VIENE A SER 0.0003
;COJO EN 5,4.5,4,3,2,1 SIGMAS O MUCHAS SIGMAS
;USO LA SENHAL MAXIMA....


restore,'maximos_dm.sav'
sma=(findgen(8)/2.+1.)*0.0003
flong=fltarr(8)
ftrans=fltarr(8)
ftot=fltarr(8)
mask=maximos(0:25,*,0)
for kk=0,7 do begin     &$
mask(*,*)=0.    &$
    ss=sma[kk]     &$
    where2,maximos(0:25,*,0) gt ss,x1,y1     &$
    where2,maximos(0:25,*,1) gt ss,x2,y2     &$
    where2,maximos(0:25,*,2) gt ss,x3,y3     &$
    mask( (x1) , (y1) )  = 1.   &$
    mask( (x2) , (y2) )  = 1.   &$
    mask( (x3) , (y3) )  = 1.   &$
    flong(kk)=mean( abs(campo*cos(incli*!pi/180.)) * alfa * mask)     &$
    ftrans(kk)=mean( campo*sin(incli*!pi/180.) * alfa * mask)     &$
    ftot(kk)=mean( campo * alfa * mask)     &$
endfor



;maximos=fltarr(727,1024,3)
;for i=0,n_files-1 do begin     &$
;    readl1_sbsp,files(i),tmp_data,hdr     &$
;    for j=0,1023 do begin    &$
;        maximos(i,j,0)=max(abs(tmp_data(5:94,j,1)/conti))    &$
;        maximos(i,j,1)=max(abs(tmp_data(5:94,j,2)/conti))    &$
;        maximos(i,j,2)=max(abs(tmp_data(5:94,j,3)/conti))    &$
;    endfor     &$
;    print,i     &$
;endfor
;save,filename='maximos_dm.sav',maximos



ruido=fltarr(727,1024,3)
for i=0,n_files-1 do begin     &$ 
    readl1_sbsp,files(i),tmp_data,hdr     &$
    for j=0,1023 do begin    &$
        ruido(i,j,0)=(tmp_data(3,j,1)/tmp_data(3,j,0))    &$
        ruido(i,j,1)=(tmp_data(3,j,2)/tmp_data(3,j,0))    &$
        ruido(i,j,2)=(tmp_data(3,j,3)/tmp_data(3,j,0))    &$
    endfor     &$
    print,i     &$
endfor
print,stddev(ruido(*,*,0))


ss=fltarr(727,1024,3)
for i=0,n_files-1 do begin    &$
    readl1_sbsp,files(i),dat,hdr    &$
    for j=0,1023 do begin    &$
        ss(i,j,0)=stddev((dat(2:12,j,1)/dat(3,j,0)))    &$
        ss(i,j,1)=stddev((dat(2:12,j,2)/dat(3,j,0)))    &$
        ss(i,j,2)=stddev((dat(2:12,j,3)/dat(3,j,0)))    &$
    endfor    &$
print,i   &$
endfor    
;save,filename='ruidos_dm.sav',maximos



pdf,0,1000,20,campo,cc,pp
plot,cc,pp,/ylog,yrange=[0.00001,0.01]

pdf,0,180,10,incli,cc,pp

pdf,0,1,0.01,alfa,cc,pp


setpsc,/landscape
tvframe,maximos(*,*,0)<0.01,/asp,/bar,/sample
tvframe,maximos(*,*,1)<0.01,/asp,/bar,/sample
tvframe,maximos(*,*,2)<0.01,/asp,/bar,/sample
endps

STOP
END
