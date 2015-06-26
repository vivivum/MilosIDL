pro v_guess,landa,landa_ref,wt,profile,vlos, beff, azi, aaa, ic, id,vf

iprofile = profile(*,0)
q=profile(*,1)
u=profile(*,2)
v=profile(*,3)

;***************
;Velocidad
;***************

Vlight=2.99792458d18 ;A/s

;Si lo busco
goto,jump
i1=(where(abs((landa-landa_ref)+0.090) eq min(abs((landa-landa_ref)+0.090))))(0)
i2=(where(abs((landa-landa_ref)+0.030) eq min(abs((landa-landa_ref)+0.030))))(0)
i3=(where(abs((landa-landa_ref)-0.030) eq min(abs((landa-landa_ref)-0.030))))(0)
i4=(where(abs((landa-landa_ref)-0.090) eq min(abs((landa-landa_ref)-0.090))))(0)
i5=(where(abs((landa-landa_ref)-0.200) eq min(abs((landa-landa_ref)-0.200))))(0)
laxis=(landa-landa_ref)([i1,i2,i3,i4,i5])/1000.

lprof=[iprofile(i1),iprofile(i2),iprofile(i3),iprofile(i4),iprofile(i5)]
;si se lo doy a pelo
jump:

pp=where(wt(*,0) eq 1)
i1=pp(0)
i2=pp(1)
i3=pp(2)
i4=pp(3)
i5=pp(4)
laxis=(landa-landa_ref)([i1,i2,i3,i4,i5])/1000D0
lprof=[iprofile(i1),iprofile(i2),iprofile(i3),iprofile(i4),iprofile(i5)]

cond=(lprof(0)+lprof(1)-lprof(2)-lprof(3))
if cond gt 0. then alpha=cond/(lprof(0)-lprof(2)) else $
  alpha=cond/(lprof(3)-lprof(1))
vf=alpha
id=sqrt(2.*((lprof(0)-lprof(2))^2.+(lprof(1)-lprof(3))^2.))
ic=2.*lprof(4)/2.+mean(lprof(0:3))        
aaa=atan(lprof(0)+lprof(1)-lprof(2)-lprof(3),$
         lprof(0)-lprof(1)-lprof(2)+lprof(3))
;vlos=2.*Vlight*0.060/!dpi/landa_ref*aaa*1.d-13
vlos=2.*Vlight*0.040/!dpi/landa_ref*aaa*1.d-13

;***************
;Campo
;***************

gamma=0.
beff=0.
azi=0.
noise=5.e-4
ff=1.

a = n_elements(iprofile) 

dera=deriv(laxis, iprofile) 
Va=v(0:3)
C1=-4.67d-13*g*landa_ref^2.*dera(0:3)
C2=total(c1^2.,/double)
beff=total(Va*C1,/double)/c2/ff

maa = where(dera eq max(dera))
mia = where(dera eq min(dera))
maa=maa[0]
mia=mia[0]

;***************
;B*sin(gamma)*ff
if beff lt 0 then gamma=135.
if beff gt 0 then gamma=45.
;***************

;***************
;AZIMUTH (Va bien)
;***************

uu=total(u)
qq=total(q)

if qq gt 3.*noise then begin
    if qq gt 0. then x0=0
    if (qq lt 0) and (uu gt 0) then x0=!pi/2.
    if (qq lt 0) and (uu lt 0) then x0=-!pi/2.
    azi=0.5*atan(uu/qq)+x0
endif else begin
    if uu gt 0 then azi=!pi/4. else azi=-!pi/4.
endelse
azi=azi*180./!pi


end



