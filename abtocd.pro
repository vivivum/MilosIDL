pro abtocd,a,rescale=rescale

;scaling the parameters

  amin=[0.001,0.,-50.,0.0001,0.0001,0.,0.,0.0001,0.0001,0.,0.]
  amax=[500.,3500.,50.,0.6000,10.,180.,360.,4.,4.,4.,1.]

  if keyword_set(rescale) then begin
    a[*]=a*(amax-amin)+amin
  endif else begin
    a[*]=(a-amin)/(amax-amin)
  endelse

end
