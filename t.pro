pro t

a = 0.09
d = findgen(100)/9.-5.5

;& if i mod 100 eq 0 then print,i 

start = SYSTIME(/SECONDS)
for i=0,30000. do begin & fvoigt,a,d,h,f & endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time
RR=5.641895836D-1
dH_u=4D0*A*F-2D0*d*H
dF_u=RR-A*H-2D0*d*F

start = SYSTIME(/SECONDS)
for i=0,30000. do begin & dummy = fft_shift(h,randomu(seed,1)+1.) & endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

start = SYSTIME(/SECONDS)
for i=0,30000. do begin & dummy = interpol([h,f],[d,d],[d+0.234,d+0.234]) & endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time

start = SYSTIME(/SECONDS)
for i=0,30000. do begin & mid=where(d eq min(abs(d))) & fvoigt,a,d(0:mid),h,f & h=temporary([h,reverse(h)]) & & f=temporary([f,reverse(f)]) & endfor
elapsed_time = SYSTIME(/SECONDS) - start
print,elapsed_time


stop

end
