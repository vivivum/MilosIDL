PRO cog, profile, eje, g, lda, blong, vlos, vlosm

;this soubroutine calculates the magnetic field with
;the COG technique
; David Orozco - Goettingen (2003)

i=profile(*,0)
q=profile(*,1)
u=profile(*,2)
v=profile(*,3)

;lda = 6301.5012
;g = 1.66667

;eje = laxis1 / 1000. + 6301.5012  ;EN A

t1 = eje * ( profile(0) - ( profile(*,0) + profile(*,3) ) )
It1 = total(t1) / (eje(2) - eje(1))
t2 =  ( profile(0) - ( profile(*,0) + profile(*,3) ) )
It2 = total(t2) / (eje(2) - eje(1))

I_mas = It1 / It2

t1 = eje * ( profile(0) - ( profile(*,0) - profile(*,3) ) )
It1 = total(t1) / (eje(2) - eje(1))
t2 =  ( profile(0) - ( profile(*,0) - profile(*,3) ) )
It2 = total(t2) / (eje(2) - eje(1))

I_menos = It1 / It2


blong = (I_mas - I_menos )/2d0 / 4.67d-13/g/lda^2.

vlosm = (I_mas + I_menos )/2d0

t1 = eje * ( profile(0) - profile(*,0) )
It1 = total(t1) / (eje(2) - eje(1))
t2 =  ( profile(0) - profile(*,0)  )
It2 = total(t2) / (eje(2) - eje(1))

vlos = (lda - it1/it2 ) * 2.99792458D+5 / lda * (-1d0)


end
