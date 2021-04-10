pro check_ajustes,data,landas,result,ch,landa,order=order,$
        n_comp=n_comp,nlte=nlte,zoom=zoom,color=color,yrange=yrange,si = si,filter=filter

if not(keyword_set(yrange)) then yrange='[-0.1,0.1]'
if not(keyword_set(si)) then yrangeI = '[0.2,1.2]'
if not(keyword_set(n_comp)) then n_comp=1
if not(keyword_set(filter)) then filter=0
if not(keyword_set(order)) then order='xylp'
if not(keyword_set(zoom)) then zoom=1
if not(keyword_set(color)) then color=0
ll = n_elements(landas)
case order of
  'xylp': begin
    p1 = 'plot,landas,data(x,y,*,0) ,/ynoz,yrange='+yrangeI
    p2 = 'plot,landas,data(x,y,*,1),yrange='+yrange
    p3 = 'plot,landas,data(x,y,*,2),yrange='+yrange
    p4 = 'plot,landas,data(x,y,*,3),yrange='+yrange
    s=size(result)
    end
  'xypl': begin
    p1 = 'plot,landas,data(x,y,0,*) ,/ynoz,yrange='+yrangeI
    p2 = 'plot,landas,data(x,y,1,*),yrange='+yrange
    p3 = 'plot,landas,data(x,y,2,*),yrange='+yrange
    p4 = 'plot,landas,data(x,y,3,*),yrange='+yrange
    s=size(result)
    end
    'lxpy': begin
      p1 = 'plot,landas,data(*,x,0,y) ,/ynoz,yrange='+yrangeI
      p2 = 'plot,landas,data(*,x,1,y),yrange='+yrange
      p3 = 'plot,landas,data(*,x,2,y),yrange='+yrange
      p4 = 'plot,landas,data(*,x,3,y),yrange='+yrange
      s=size(result)
      end
  else: print,'Fack'
endcase

loadct,color
s=size(result)
window,0,xsize=s(1)*zoom,ysize=s(2)*zoom
tvscl,congrid(result(*,*,7)+result(*,*,8),s(1)*zoom,s(2)*zoom),0
window,1,xsize=500,ysize=500
wset,0
colores

init_milos,landa,wl

goto,jjjj
xlast=-1
repeat begin
CROSSM,x,y,[s(1),s(2),0,1,1],BUT=ib,xy_window=1
if ib eq 1 then begin
    ;genera fit
    initi=reform(result(x,y,*))
    milos, wl,landas,initi, pr,/synthesis,n_comp=n_comp,filter=filter
    xlast=x
    wset,1
    !p.multi=[0,2,2]
    e=execute(p1)
    oplot,landas,pr(*,0),color=4
    e=execute(p2)
    oplot,landas,pr(*,1),color=4
    e=execute(p3)
    oplot,landas,pr(*,2),color=4
    e=execute(p4)
    oplot,landas,pr(*,3),color=4
    !p.multi=0
    wset,0
    print,reform(result(x,y,*))
    print,x,y
endif
endrep until ib eq 2

jjjj:

!mouse.button=0
xb=0
yb=0
device,set_graphics=6
plots,/dev,[xb,xb]*zoom,[0,s(2)]*zoom
plots,/dev,[0,s(1)]*zoom,[yb,yb]*zoom
while !mouse.button ne 4 do begin
    cursor,x,y,2,/device
    x = x / zoom
    y = y / zoom
    if ( ((x ne xb) or (y ne yb)) and ((x ge 0) and (x lt s(1)-1) and (y ge 0) and (y lt s(2)-1))) then begin
        plots,/dev,[xb,xb]*zoom,[0,s(2)]*zoom
        plots,/dev,[0,s(1)]*zoom,[yb,yb]*zoom
        plots,/dev,[x,x]*zoom,[0,s(2)]*zoom
        plots,/dev,[0,s(1)]*zoom,[y,y]*zoom
        device,set_graphics=3
 	initi=reform(result(x,y,*))
        milos, wl,landas,initi, pr,/synthesis,n_comp=n_comp,nlte=nlte,filter=filter
  	xlast=x
        wset,1
        colores
        !p.multi=[0,2,2]
        e=execute(p1)
    	oplot,landas,pr(*,0),color=4
        e=execute(p2)
    	oplot,landas,pr(*,1),color=4
        plots,[landas[0],landas[ll-1]],[-0.003,-0.003],line=3
        plots,[landas[0],landas[ll-1]],[0.003,0.003],line=3
        
        e=execute(p3)
    	oplot,landas,pr(*,2),color=4
        plots,[landas[0],landas[ll-1]],[-0.003,-0.003],line=3
        plots,[landas[0],landas[ll-1]],[0.003,0.003],line=3
        e=execute(p4)
    	oplot,landas,pr(*,3),color=4
        plots,[landas[0],landas[ll-1]],[-0.003,-0.003],line=3
        plots,[landas[0],landas[ll-1]],[0.003,0.003],line=3

        print,'x: ',text(x),' y: ',text(y)
        print,'B (B): ',text(result(x,y,1)),' Inc (G): ',text(result(x,y,5)),' Azi (A): ',text(result(x,y,6)),' Vdop (D): ',text(result(x,y,3))
        print, 'damp: (T)',text(result(x,y,4)),' eta0 (H): ',text(result(x,y,0)),' VLOS (V): ',text(result(x,y,2)),' CHISQR (I): ',text(ch(x,y))
        ;print, 'A1: (E)',text(result(x,y,10)),' appha1 (R): ',text(result(x,y,11)),' A2 (T): ',text(result(x,y,12)),' alpha2: (Y)',text(result(x,y,13))
        xb=x & yb=y
        loadct,4,/silent
        wset,0
        device,set_graphics=6
     endif
     case strupcase(get_kbrd(0)) of
' ': print,'gg'
'C':	begin
        device,set_graphics=3
        tvscl,congrid(reform(data(0,*,*,0)),s(1)*zoom,s(2)*zoom)>0.2
        device,set_graphics=6
	endcase
'I':	begin
        device,set_graphics=3
        tvscl,congrid(ch,s(1)*zoom,s(2)*zoom)<1300,0
        device,set_graphics=6
	endcase
'B':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,1),s(1)*zoom,s(2)*zoom)<300,0
        device,set_graphics=6
	endcase
'H':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,0),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'T':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,4),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'G':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,5),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'A':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,6),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'D':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,3),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'V':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,2),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'E':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,10),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'R':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,11),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'T':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,12),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
'Y':	begin
        device,set_graphics=3
        tvscl,congrid(result(*,*,13),s(1)*zoom,s(2)*zoom),0
        device,set_graphics=6
	endcase
  else: aa=0;print,'nothing'
endcase

    if !mouse.button eq 2 then begin
       print,'HHHH'
    endif
endwhile
device,set_graphics=3




;if ib eq 4 then begin
;    fitea,wl,landas,reform(tmp_data(*,y,*)),yfit,chinew,fitpar,pol=noise*4.,stray=straylight
;    initi=reform(result(x,y,*))
;    milos,wl,landas,initi,pr,/synthesis,slight=straylight_local(*,*,y)
;    xlast=x
;    wset,1
;    !p.multi=[0,4,1]
;    plot,tmp_data(*,y,0),/ynoz
;    oplot,yfit(*,0),color=3
;    oplot,pr(*,0),color=4
;    plot,tmp_data(*,y,1)
;    oplot,yfit(*,1),color=3
;    oplot,pr(*,1),color=4
;    plot,tmp_data(*,y,2)
;    oplot,yfit(*,2),color=3
;    oplot,pr(*,2),color=4
;    plot,tmp_data(*,y,3)
;    oplot,yfit(*,3),color=3
;    oplot,pr(*,3),color=4
;    !p.multi=0
;    wset,0
;    print,reform(result(x,y,*))
;    print,fitpar
;    print,chi(x,y),'New CHI:',chinew
;    print,x,y
;endif

end
