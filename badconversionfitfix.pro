function badconversionfitfix,data,level

on_error,2

nl = level lt 1 ? level : 0.2

;ss=size(gram.b_long)
;nx=ss(1) & ny=ss(2)
;tags=[1,2,3,4,5,6] ; will need to correct b_long -> b_fill, not others
;easytags=[1,4,5,6]

a=data
snx=(size(data))(1) & sny=(size(data))(2)
mask=bytarr(snx,sny)
maskfinal=bytarr(snx,sny)

; correct spikes
ele=[0,1,2,3,4,5,6,7,8,9,10]
for i=0,10 do begin
aa=a(*,*,ele(i))
saa=size(aa)
nloop=0
repeat begin
	if saa(0) gt 1 then begin
   xa = abs((aa-shift(aa,1,0))*(aa-shift(aa,-1,0)))
   ya = abs((aa-shift(aa,0,1))*(aa-shift(aa,0,-1)))
   endif else return,[[[a]],[[maskfinal]]]
   xm = max(xa) & ym = max(ya)
   xa = xa / xm & ya = ya / ym 
   xx = where(xa gt nl,countxx)
   yy = where(ya gt nl,countyy)
   if (countxx+countyy gt 5.) then begin
	mask([xx,yy])=1 
   ;create smooth image
   sm = smooth(aa,7,/edge_truncate)
   bb = aa * (1 - mask) + sm * mask
   aa([xx,yy]) = bb([xx,yy])
;   ;create smooth image 
;   sm = smooth(aa,7,/edge_truncate)
;   bb = aa * (1 - mask) + sm * mask
;   aa([xx,yy]) = bb([xx,yy])
	endif else begin
	print, 'No pixels!!!',ELE(I)
	endelse
   nloop = nloop + 1
endrep until nloop eq 2
a(*,*,ele(i)) = aa
maskfinal = maskfinal + mask
endfor

maskfinal = maskfinal < 1.

return,[[[a]],[[maskfinal]]]

end

pro ttt

;.r badconversionfitfix  
b=badconversionfitfix(a,0.2)

window,0,xsize=500,ysize=500
window,1,xsize=500,ysize=500


wset,0
erase
tvscl,b(*,*,11)                 
wset,1
erase
tvscl,b(*,*,5)                     
blink,[0,1]

end