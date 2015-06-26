pro acc,n,fxx,par,pnew

pnew = par[*,0]

wt = 1d0 ;/ par[fxx,n]
Coef1 = par[fxx,n] - 2d0*par[fxx,n-1] + par[fxx,n-2]
A1 = Total( Coef1^2d0 * wt) 
Coef2 =  (par[fxx,n] - par[fxx,n-1] - par[fxx,n-2] + par[fxx,n-3])
B1 = Total( Coef2 * wt * Coef1 ) 
C1 = Total( Coef1 * (par[fxx,n] - par[fxx,n-1]) * wt )
A2 = B1
B2 = Total( Coef2^2d0*wt )
C2 = Total( Coef2 * (par[fxx,n] - par[fxx,n-1]) * wt )

Coef3 = A1*B2 - A2 * B1

if coef3 eq 0 then begin
x_new = par[fxx,n] 
endif else  begin

a = (C1*B2 - C2*B1 ) / Coef3
b = (C2*A1 - C1*A2) / Coef3
x_new = (1 - a - b)*par[fxx,n] + a*par[fxx,n-1] + b*par[fxx,n-2]

endelse

pnew(fxx) = x_new

end
