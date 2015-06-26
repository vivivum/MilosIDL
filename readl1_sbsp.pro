pro readl1_sbsp,in_file,dat,hdr

;  routine to read in and properly normalize the level1 SP data

;	INPUT:		in_file	= the fits input file

;	OUTPUTS:	dat	= the floating point Stokes data
;			hdr	= the level1 fits header


n_files=n_elements(in_file)

if n_files gt 1 then begin

if n_files gt 100 then print,'MEMORY CAN BE FILLEDDDDD'

dat = float(readfits(in_file(0),hdr))
s=size(dat)

dat=fltarr(n_files,s(1),s(2),s(3))

for i=0,n_files-1 do begin

dat_tmp=float(readfits(in_file(i),hdr))

temp = dat_tmp(*,*,0)
whr = where(temp lt 0., countwrap)
if countwrap ne 0 then temp(whr) = temp(whr) + 65536.
dat_tmp(0,0,0) = temp
bitshft = sxpar(hdr,'SPBSHFT')
print,bitshft
switch bitshft of
    3:      dat_tmp(*,*,1:2) = 2.*dat_tmp(*,*,1:2)
    2:      dat_tmp(*,*,3) = 2.*dat_tmp(*,*,3)
    1:      dat_tmp(*,*,0) = 2.*dat_tmp(*,*,0)
endswitch

dat(i,*,*,*)=dat_tmp

endfor

endif else begin

dat = float(readfits(in_file,hdr))

;  adjust for intensity wraparound
temp = dat(*,*,0)
whr = where(temp lt 0., countwrap)
if countwrap ne 0 then temp(whr) = temp(whr) + 65536.
dat(0,0,0) = temp
;  adjust for intensity wraparound TWICE RETURN 0 IF NOTHING FOUND
;temp = dat(*,*,0)
;whr = where(temp lt 0., countwrap)
;if countwrap ne 0 then temp(whr) = temp(whr) + 65536. else print,'NOTHING TO WRAPBACK'
;dat(0,0,0) = temp


;  bit shifting for SP data
        bitshft = sxpar(hdr,'SPBSHFT')
;  account for bit shifting of the data by multiplying the appropriate
;  Stokes images by 2
        switch bitshft of
                3:      dat(*,*,1:2) = 2.*dat(*,*,1:2)
                2:      dat(*,*,3) = 2.*dat(*,*,3)
                1:      dat(*,*,0) = 2.*dat(*,*,0)
        endswitch

    endelse

return
end
