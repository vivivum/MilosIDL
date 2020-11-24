pro set_style,ppt=ppt,ps=ps,back=back

common set_style_var,var_p_thick,var_x_thick,var_y_thick,var_p_charthick,alreadyset

if keyword_set(back) then begin
    !p.thick=var_p_thick 
    !x.thick=var_x_thick
    !y.thick=var_y_thick
    !p.charthick=var_p_charthick
    alreadyset = 0
endif else begin
   if not(keyword_set(alreadyset)) then begin
    var_p_thick=!p.thick
    var_x_thick=!x.thick
    var_y_thick=!y.thick
    var_p_charthick=!p.charthick
    endif
    if keyword_set(ppt) then begin
        !p.thick=5
        !x.thick=5
        !y.thick=5
        !p.charthick=5
        alreadyset = 1
    endif else if keyword_set(ps) then begin
        !p.thick=2
        !x.thick=2
        !y.thick=2
        !p.charthick=2
        alreadyset = 1
    endif else print,'no fmt selected'
endelse
end

pro setpscf, filename=filename,_extra=_extra,aspect_ratio=aspect_ratio,mydevice = mydevice

  mydevice = !D.NAME

if not keyword_set(filename) then filename = 'idl.ps'

if not keyword_set(aspect_ratio) then  aspect_ratio=1.5
xsize=18
ysize=xsize/aspect_ratio

set_plot, 'PS' 
!p.font=0
device, filename=filename,/color,LANGUAGE_LEVEL=2,bits=8,_extra=_extra,/encapsulated, /helvetica
device, xsize=xsize, ysize=ysize

end

pro endpscf,mydevice = mydevice
device, /close
if keyword_set(mydevice) then set_plot,mydevice else set_plot, 'X'
!p.font=-1
end
