;pro milos_int_v5, wl, axis, imodel,fmodel, yy, chisqr=chisqr,yfited=yfited,sigma=sigma,fix=fix,toplim=toplim,$
;	ilambda=ilambda,miter=miter,Err=Err,filter=filter,weight=weight,badpixel=badpixel,simple=simple,$
;	doplot=doplot,icontlimit=icontlimit

pro milos_int_v5, YY,VARS,fmodel,yfited,Error,chisqr,slight=slight,mlocal=mlocal,doplot=doplot

wl = vars.wl
axis = vars.axis
imodel = vars.init_model 
sigma = vars.sigma 
fix = vars.fix 
ilambda = vars.ilambda 
miter = vars.miter 
filter = vars.filter
badpixel = vars.badpixel
simple = vars.simple 
icontlimit = vars.icontlimit
weight = vars.weight 

sz=size(yy)
s_x = sz(1) 
s_y = sz(2) 
fmodel = dblarr(s_x,s_y,11 )
error = dblarr(s_x,s_y,11)
chisqr = dblarr(s_x,s_y)
sw = size(weight)
yfited = dblarr( s_x,s_y, sw(1), sz(4))

;milos_ex =    'milos,wl, axis, init_model, reform(yy(x,y,*,*)),chisqr=chisqrt,yfit=yfit,/quiet,$
;              sigma=sigma,fix=fix,/inversion,ilambda=ilambda,miter=miter,Err=Error,filter=filter,$
;              weight=weight,VARLIMITS = [[7, 1, -0.05, 0.05],[8, 1, -0.05, 0.05],[3, 1, -0.01, 0.01]],$
;              PARLIMITS = [0, 1, 1d0, 150d0]'

if simple eq 1 then begin
    for x = 0 , s_x-1   DO BEGIN
        for y = 0 , s_y-1   DO BEGIN
            init_model = imodel
            if max(yy(x,y,*,0)) gt icontlimit then begin
            if keyword_set(slight) then begin
            	if keyword_set(mlocal) then begin
            		stray = reform(slight(x,y,*,*))
            	endif else stray = slight
            endif
            	
            milos,wl, axis, init_model, reform(yy(x,y,*,*)),chisqr=chisqrt,yfit=yfit,$
              sigma=sigma,fix=fix,/inversion,ilambda=ilambda,miter=miter,Err=Err,filter=filter,$
              weight=weight,VARLIMITS = [[7, 1, -0.1, 0.1],[8, 1, -0.1, 0.1],[3, 1, -0.005, 0.005]],$
              PARLIMITS = [0, 1, 1d0, 150d0],/quiet,doplot=doplot,slight=stray,mlocal=mlocal
            ;print,x,y
            fmodel (x ,y,*)= init_model
            error (x ,y,*)= err
            chisqr (x ,y)=  chisqrt
            yfited (x ,y, *, *)=  yfit
            endif else begin
            fmodel (x ,y,*)= 0
            err (x ,y,*)= 0
            chisqr (x ,y)=  0
            yfited (x ,y, *, *)=  0
            
            endelse
        endfor
    endfor
endif else begin
    
    imodel1=[[imodel],[imodel],[imodel]]
    imodel1(1,0)=300.
    imodel1(1,1)=800.
    imodel1(1,2)=1300.
    for x = 0 , s_x-1   DO BEGIN
        for y = 0 , s_y-1   DO BEGIN
            for im = 0, 2 do begin
             if max(yy(x,y,*,0)) gt icontlimit then begin
                 if keyword_set(slight) then begin
            	if keyword_set(mlocal) then begin
            		stray = reform(slight(x,y,*,*))
            	endif else stray = slight
            endif

                init_model = imodel1(*,im)
                milos,wl, axis, init_model, reform(yy(x,y,*,*)),chisqr=chisqrt,yfit=yfit,/quiet,$
                  sigma=sigma,fix=fix,/inversion,ilambda=ilambda,miter=miter,Err=Err,filter=filter,$
                  weight=weight,VARLIMITS = [[7, 1, -0.1, 0.1],[8, 1, -0.1, 0.1],[3, 1, -0.005, 0.005],[0, 1, -10d0, 10d0]],$
                  PARLIMITS = [[0, 1, 1d0, 250d0] , [3, 1, 0.01d0, 0.06d0], [4, 1, 0.01d0, 5d0]],doplot=doplot,toplim=toplim,slight=stray,mlocal=mlocal
                ;print,chisqr(x,y),chisqrt
                if (im eq 0) OR (chisqrt lt chisqr(x,y)) then begin
                    fmodel (x ,y,*)= init_model
                    error (x ,y,*)= err
                    chisqr (x ,y)=  chisqrt
                    yfited (x ,y, *, *)=  yfit
                endif
                            endif else begin
            fmodel (x ,y,*)= 0
            error (x ,y,*)= 0
            chisqr (x ,y)=  0
            yfited (x ,y, *, *)=  0
            
            endelse

            endfor
        endfor
    endfor
    
endelse


if badpixel eq 1 then begin
    temp = badconversionfitfix(fmodel,0.2)
    init_block = temp(*,*,0:10)
    
    for x = 0 , s_x-1  DO BEGIN
        for y = 0 , s_y-1  DO BEGIN
                        if max(yy(x,y,*,0)) gt icontlimit then begin
            init_model = reform(init_block(x,y,*))
            milos,wl, axis, init_model, reform(yy(x,y,*,*)),chisqr=chisqrt,yfit=yfit,/quiet,$
              sigma=sigma,fix=fix,/inversion,ilambda=ilambda*10.,miter=miter+30.,Err=Err,filter=filter,$
              weight=weight,VARLIMITS = [[7, 1, -0.05, 0.05],[8, 1, -0.05, 0.05],[3, 1, -0.005, 0.005],[0, 1, -10d0, 10d0]],$
              PARLIMITS = [[0, 1, 1d0, 250d0], [3, 1, 0.01d0, 0.06d0], [4, 1, 0.01d0, 5d0]],doplot=doplot,slight=stray,mlocal=mlocal
            ;print,chisqr(x,y),chisqrt
            if chisqrt lt chisqr(x,y) then begin
                fmodel (x ,y,*)= init_model
                error(x ,y,*)= err
                chisqr (x ,y)=  chisqrt
                yfited (x ,y, *, *)=  yfit
                print,'better'
            endif else print,'worse'
                        endif else begin
            fmodel (x ,y,*)= 0
            error (x ,y,*)= 0
            chisqr (x ,y)=  0
            yfited (x ,y, *, *)=  0
            
            endelse

        endfor
    endfor
endif

end
