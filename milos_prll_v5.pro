pro MILOS_prll_v5,THREADS=threads,Block_data=block_data,xrange=xrange,yrange=yrange,$
                  final_model=final_model,SBLOCKS=sblocks,nCONT=ncont,sav_file=sav_file,$
                  badpixel=badpixel,simple=simple,max_iter=max_iter,init_lambda=init_lambda,$
                  not_normalize=not_normalize,waxis = waxis,INITIAL_MODEL=INITIAL_MODEL,IMAX_SAVE=imax_save,$
                  filter_width=filter_width,weights=weights,fix=fix,sigma=sigma,icontlimit=icontlimit,line=line,muestreo=muestreo,$
                  slight=slight,mlocal=mlocal,ipbs=ipbs,n_components=n_components

;Block_data = variable containing the data [X, Y, Wavelengths, Stokes profiles], or a save file (ending in .sav)
;line = spectral line (see milos)

;THREADS: Max number of CPUs. if 1 it does not use IDL_BRIGDE

;XRANGE and YRANGE are useful if one wants to inver only apiece of
;data. Default, all pixels are inverted. Example: xrange=[0,100],yrange[0,200]
;IMPORTANT: when selected destroy Block_data and replaces with the region

;SBLOCKS is the dimension of the square blocks (in pixels) in to which
;the program will devide the image to be process in the different
;threads, default= 25 px

;Icont: average continuum of the data to which the programs normalize
;the data. By default, the program normalize the data to the average
;continuum taking all pixels into account.

;not_normalize To avoid normalization

;sav_file: sav filename to store the results

;final_model: results from the inversion

;badpixel: if set look for bad pixels and rerun inversion
;simple: if set, only one inversion is carried out. This programs
;invert three times each pixel

;max_iter: maximum number of iterations. Default 30. (see MILOS)
;init_lambda: initial lambda parameter. Default 1. (see MILOS)
;filter_width: filter width (defaul 85ma)

;fix... see milos
;sigma... see milos

;TESTS.....
;IDL> milos_prll_v5,xrange=[50,150],yrange=[50,150],/badpixel,imax_sav='../Data/reduc_rnr_163_208.save',max_iter=100.,init_lambda=10d0,sav_file='test1'
;$cp milos_int_v5.pro test1_pro

;milos_prll_v5,xrange=[50,900],yrange=[50,900],/badpixel,imax_sav='../early_paper_data/reduc_rnr_163_208.save',max_iter=100.,init_lambda=10d0,sav_file='test_confiltro'

;milos_prll_v5,xrange=[50,900],yrange=[50,900],/badpixel,imax_sav='../Data/nuevos_110110/reduc_rnr_163_208.save',max_iter=100.,init_lambda=10d0,sav_file='test_confiltro'

;***************************************************
print, '--------------------------------------------------     '
print, ' Welcome to Milos in Parallel Wrapper '
print, '   (Beta version)                           '
print, '  D. Orozco 22/Jul/2011 d.orozco@nao.ac.jp '
print, '  Modified D. Orozco 27/Jun/2014 dorozco@iac.es '
print, '---------------------------------------------------     '
;****************************************************

;READ INPUT IMaX DATA. The program assumes IID name for IMaX data and all data in block input

if keyword_set(IMaX_save) then begin
	;save_img=OBJ_NEW('IDL_Savefile', '../Data/reduc_163_346.save')
	check=file_test(IMaX_save)
	if check eq 0 then message,"IMaX_save file does not exist" else print,'reading... '+IMaX_save
	restore,IMaX_save,/verbose
	;save_img=OBJ_NEW('IDL_Savefile', IMaX_save)
	;var_name=save_img->Names()
	;dimensions=save_img->Size(var_name(0),/dimensions)
    Print,'Restoring data...'
	restore,'../Data/reduc_163_346_r.save',/verbose
	;IID             FLOAT     = Array[936, 936, 5, 4]
	block_data = fltarr(936, 936, 12, 4)
	IID(*,*,*,0) = IIDN(*,*,*,0)
	IID(*,*,*,3) = IIDN(*,*,*,1)*1.8
	IIDN=0
endif


; DIMENSIONS OF DATA
;the program assumes N wavelengths samples and 4 Stokes profiles

sz = size(block_data)
if sz(0) ne 4 then begin
print,'Image dimensions not valid'
return
endif

npr_x = sz(1) ;X size
npr_y = sz(2) ;Y size
n_landas = sz(3)

;INVERSION RANGE
;we take the data we want to analyze only

IF keyword_set(xrange) and keyword_set(yrange) then begin
	 if xrange(0) lt 0 or xrange(1) gt npr_x then message,"Xrange exceeds data dimensions... stopped"
	 if yrange(0) lt 0 or yrange(1) gt npr_y then message,"Yrange exceeds data dimensions... stopped"
	 Block_data = temporary(Block_data(xrange(0):xrange(1),yrange(0):yrange(1),*,*))
endif else if keyword_set(xrange) and not(keyword_set(yrange)) then begin
 	 if xrange(0) lt 0 or xrange(1) gt npr_x then message,"Xrange exceeds data dimensions... stopped"
	 Block_data = temporary(Block_data(xrange(0):xrange(1),0:npr_y-1,*,*))
endif else if not(keyword_set(xrange)) and keyword_set(yrange) then begin
 	 if yrange(0) lt 0 or yrange(1) gt npr_y then message,"Yrange exceeds data dimensions... stopped"
	 Block_data = temporary(Block_data(0:npr_x-1,yrange(0):yrange(1),*,*))
endif else begin
	 xrange = [0,npr_x-1] & yrange=[0,npr_y-1]
endelse

;Update dimensionss_x,s_y

sz = size(Block_data)
npr_x = sz(1)
npr_y = sz(2)
npr_total = npr_x * npr_y ;total number of pixels
if not(keyword_set(sblocks)) then blocks = 25 else blocks=sblocks  ;Size of the blocks

;Normalization of the data to average continuum
if not(keyword_set(not_normalize)) then begin
	Print,'Normalization to averaged continuum...'
	if not(keyword_set(ncont)) then begin
		Ic = mean (Block_data(*,*,1:5,0))
		Block_data = Block_data / Ic
		Print,'Averaged continuum...(counts): ',Ic
    endif else begin
		Block_data = temporary(Block_data) / ncont
		Print,'Averaged continuum...(counts) (user defined): ',ncont
    endelse
endif

;**********************************************************
;VARIABLE DEFINITION  INVERSION PART
;**********************************************************

if not(keyword_set(line)) then print,'Not line selected...'

init_milos,line,wl

if not(keyword_set(icontlimit)) then icontlimit=0
;fitting parameters
if not(keyword_set(fix)) then fix=[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.]
;stdev of the data
if not(keyword_set(sigma)) then sigma=[1.,1.,1.,1.]/1000.
;maximum iterations
if keyword_set(max_iter) then miter = max_iter else miter=30
if keyword_set(init_lambda) then ilambda=init_lambda else ilambda=1d0
if not(keyword_set(filter_width)) then filter=0d0 else filter = filter_width

if not(keyword_set(waxis)) then begin
print,'landas not provided...'
return
endif
axis = waxis
landas = n_elements(axis)
;landas y n_landas debe coincidir, sino, necesito muestreo
if (n_elements(landas) ne n_elements(n_landas)) and not(keyword_set(muestreo)) then STOP

if not(keyword_set(muestreo)) then muestreo = !NULL

wt = fltarr(landas,4)
if not(keyword_set(weights)) then begin
if muestreo NE !NULL then begin
	wt(muestreo,0) = 1.0
	wt(muestreo,1) = 1.0
	wt(muestreo,2) = 1.0
	wt(muestreo,3) = 1.0
endif else begin
	wt(*,0) = 1.0
	wt(*,1) = 4.0
	wt(*,2) = 4.0
	wt(*,3) = 2.0
endelse
endif else wt = weights

if not(keyword_set(INITIAL_MODEL)) then begin

;initial model
S0=0.2d0 & S1=0.95d0
eta0=4.5d0
Magnet=1200.
GAMMA=45d0
AZI=45d0
vlos=0.25d0  ;km/s
MACRO=0d0
LANDADOPP=0.03d0
aa=0.09d0
alfa=0d0

INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,S0,S1,macro,alfa]

endif else INIT_MODEL = INITIAL_MODEL

;definition of output parameters

FINAL_chisqr=fltarr(npr_x,npr_y)
FINAL_MODEL=fltarr(npr_x,npr_y,11)
ERRORS=fltarr(npr_x,npr_y,11)
;FINAL_FIT=fltarr(npr_x,npr_y,n_landas,4)

;definition of some other parameters
INIT_MODEL_BLOCK = fltarr(BLOCKS,BLOCKS,11)

if keyword_set(badpixel) then badpixel = 1 else badpixel = 0
if keyword_set(simple) then simple = 1 else simple = 0

;**********************************************************
;VARIABLE DEFINITION  INVERSION PART END
;**********************************************************

;Number of threads
nbparallel = !CPU.TPOOL_NTHREADS ; 1 thread per core
if keyword_set(threads) then begin
    nbparallel = nbparallel < threads
endif else begin
    nbparallel = !CPU.TPOOL_NTHREADS ; 1 thread per core
    threads = nbparallel
endelse


vars = {wl:fltarr(2), axis:fltarr(landas), init_model:fltarr(11), sigma:fltarr(4) ,fix:fltarr(11),$
	ilambda:0.0, miter:0.0,filter:0.0,badpixel:0,simple:0,icontlimit:0.0,weight:fltarr(landas,4)}
vars.wl = wl
vars.axis = axis
vars.init_model = init_model
vars.sigma = sigma
vars.fix = fix
vars.ilambda = ilambda
vars.miter = miter
vars.filter = filter
vars.badpixel = badpixel
vars.simple = simple
vars.icontlimit = icontlimit
vars.weight = wt

IF THREADS GT 1 then begin

;**********************************************************
;PARALLEL PART
;**********************************************************

;Startup file
startup_file = '@' + PREF_GET('IDL_STARTUP')

bridges = ptrarr(nbparallel)

print, 'Available number of threads: ', nbparallel

;Initialize child processes

for ipar=0L,nbparallel-1 do begin
    bridges[ipar] = ptr_new(obj_new('IDL_IDLBridge'))
    print,'Built thread number ',ipar
endfor

;oBridge[i] = obj_new('IDL_IDLBridge', Callback='bridgeFunctionCallback')

;  oBridge[i].setProperty, userData=0

;setup IDL and milos in each child

for ipar=0L,nbparallel-1 do begin
    print,'Initiallizing IDL/PATH in thread... ',ipar
    (*bridges[ipar])->Execute, startup_file
    if (*bridges[ipar])->Status() then message,'Thread busy. Come back again later.',/info
endfor

for ipar=0L,nbparallel-1 do begin
    print,'Initiallizing MILOS in thread ',ipar
	(*bridges[ipar])->SetVar, "line", line
    (*bridges[ipar])->Execute, "init_milos,''+line+'',wl"
	print,(*bridges[ipar])->GetVar( "wl")
endfor

for ipar=0L,nbparallel-1 do begin
 (*bridges[ipar])->SetVar, "landas", landas
 (*bridges[ipar])->Execute,"vars ={wl:fltarr(2), axis:fltarr(landas), init_model:fltarr(11), sigma:fltarr(4) ,fix:fltarr(11),ilambda:0.0, miter:0.0,filter:0.0,badpixel:0,simple:0,icontlimit:0.0,weight:fltarr(landas,4)}"
 (*bridges[ipar])->SetVar, "vars.wl", vars.wl
 (*bridges[ipar])->SetVar, "vars.axis", vars.axis
 (*bridges[ipar])->SetVar, "vars.init_model", vars.init_model
 (*bridges[ipar])->SetVar, "vars.sigma", vars.sigma
 (*bridges[ipar])->SetVar, "vars.fix", vars.fix
 (*bridges[ipar])->SetVar, "vars.ilambda", vars.ilambda
 (*bridges[ipar])->SetVar, "vars.miter", vars.miter
 (*bridges[ipar])->SetVar, "vars.filter", vars.filter
 (*bridges[ipar])->SetVar, "vars.badpixel", vars.badpixel
 (*bridges[ipar])->SetVar, "vars.simple", vars.simple
 (*bridges[ipar])->SetVar, "vars.icontlimit", vars.icontlimit
 (*bridges[ipar])->SetVar, "vars.weight", vars.weight
endfor

;INPUT DATA FOR milos_int_v5.pro

;Send vars to bridges
;for ipar=0L,nbparallel-1 do (*bridges[ipar])->SetVar, "vars", vars
;for ipar=0L,nbparallel-1 do (*bridges[ipar])->setproperty, userdata=vars

running = 1
stats =intarr(nbparallel)
stats[*] = 0
nblock_x = ceil(npr_x*1./blocks) ;+ 1 ;number of blocks in x direction
nblock_y = ceil(npr_y*1./blocks) ;+ 1 ;number of blocks in y direction
failed=0L
iprof = 0L
iprof_x = 0L
iprof_y = 0L
once = 1

print,'Total blocks: ', strtrim(string(nblock_x * nblock_y),2)
print,'Blocks in x and y dim: ',' [',strtrim(string(nblock_x),2),',',strtrim(string(nblock_y),2),']'
print,'Size of the blocks: ',blocks


while (running) do begin ;MAIN LOOP

;send a profile

if once then begin
for ipar=0L,nbparallel-1 do begin
	if stats[ipar] eq 0 then begin
            toop_x = (blocks*(iprof_x+1)-1)<(npr_x-1) - blocks*iprof_x
            toop_y = (blocks*(iprof_y+1)-1)<(npr_y-1) - blocks*iprof_y
            ely = Block_data(blocks*iprof_x:(blocks*(iprof_x+1)-1)<(npr_x-1), blocks*iprof_y:(blocks*(iprof_y+1)-1)<(npr_y-1),*,*)
            ;en caso de que el muestreo sea diferente, este entra en muestreo
            if muestreo NE !NULL then begin
            y = fltarr(toop_x+1 , toop_y+1,landas,4)
            for ii = 0, n_landas-1 do y(*,*,muestreo(ii),*) = ely(*,*,ii,*)
            endif else y = ely
            if keyword_set(slight) then begin
            	if not(keyword_set(mlocal)) then begin
                	(*bridges[ipar])->SetVar, "slight", slight
            	endif else begin
                	slight_block = slight(blocks*iprof_x:(blocks*(iprof_x+1)-1)<(npr_x-1), blocks*iprof_y:(blocks*(iprof_y+1)-1)<(npr_y-1),*,*)
                	(*bridges[ipar])->SetVar, "slight", slight_block
                	(*bridges[ipar])->SetVar, "mlocal", mlocal
                endelse
            endif
            (*bridges[ipar])->SetVar, "y", y
            (*bridges[ipar])->SetVar, "iprof_x", iprof_x ;x block id
            (*bridges[ipar])->SetVar, "iprof_y", iprof_y ;y block id
            (*bridges[ipar])->Execute, 'milos_int_v5,Y,VARS,final_model,yfited,Err,chisqr,slight=slight,mlocal=mlocal',/nowait

        print,'Processing block ',strtrim(string(iprof_x+iprof_y*nblock_x+1),2),$
          ', [x,y] = [',strtrim(string(iprof_x),2),',',strtrim(string(iprof_y),2),'], ',$
          ' of ', strtrim(string(nblock_x * nblock_y),2)

		iprof = iprof + 1 ;update block
                iprof_x = iprof mod nblock_x
                iprof_y = iprof/nblock_x
	endif
endfor
once = 0
endif

for ipar=0L,nbparallel-1 do begin
    stats[ipar] = (*bridges[ipar])->Status()
    if stats[ipar] eq 2 then begin

        which_iprof_x = (*bridges[ipar])->GetVar('iprof_x') ;; get the profile id
        which_iprof_y = (*bridges[ipar])->GetVar('iprof_y') ;; get the profile id
        var1 = (*bridges[ipar])->GetVar("final_model") ;; get the result
        var2 = (*bridges[ipar])->GetVar("Err") ;; get the result
        var3 = (*bridges[ipar])->GetVar("chisqr") ;; get the result
        ;var4 = (*bridges[ipar])->GetVar("yfited") ;; get the result

		print,'Thread ',strtrim(string(ipar),2),' finished.',' Block ',$
                  strtrim(string(which_iprof_x+which_iprof_y*nblock_x+1),2),$
                  ', [x,y] = [',strtrim(string(which_iprof_x),2),',',$
                  strtrim(string(which_iprof_y),2),'], ',' of ', $
                  strtrim(string(nblock_x * nblock_y),2),'. Avr. Chi^2 VALUE: ',mean(var3)

        de_x = blocks*which_iprof_x<(npr_x-1)
        hasta_x = (blocks*(which_iprof_x+1)-1)<(npr_x-1)
		de_y = blocks*which_iprof_y<(npr_y-1)
		hasta_y = (blocks*(which_iprof_y+1)-1)<(npr_y-1)

		final_model(de_x:hasta_x,de_y:hasta_y,*) = var1
		errors(de_x:hasta_x,de_y:hasta_y,*) = var2
		final_chisqr(de_x:hasta_x,de_y:hasta_y,*) = var3
		;final_fit(de_x:hasta_x,de_y:hasta_y,*,*) = var4

		if iprof lt nblock_x*nblock_y then begin ; set new profile


                    de_x = blocks*iprof_x<(npr_x-1)
                    hasta_x = (blocks*(iprof_x+1)-1)<(npr_x-1)
                    de_y = blocks*iprof_y<(npr_y-1)
                    hasta_y = (blocks*(iprof_y+1)-1)<(npr_y-1)

		toop_x = (blocks*(iprof_x+1)-1)<(npr_x-1) - blocks*iprof_x
		toop_y = (blocks*(iprof_y+1)-1)<(npr_y-1) - blocks*iprof_y
	               ely = Block_data(de_x:hasta_x,de_y:hasta_y,*,*)
            ;en caso de que el muestreo sea diferente, este entra en muestreo
			if muestreo NE !NULL then begin
        	     y = fltarr(toop_x+1 , toop_y+1,landas,4)
        	    for ii = 0, n_landas-1 do y(*,*,muestreo(ii),*) = ely(*,*,ii,*)
           	 endif else y = ely
            if keyword_set(slight) then begin
            	if not(keyword_set(mlocal)) then begin
                	(*bridges[ipar])->SetVar, "slight", slight
            	endif else begin
                	slight_block = slight(blocks*iprof_x:(blocks*(iprof_x+1)-1)<(npr_x-1), blocks*iprof_y:(blocks*(iprof_y+1)-1)<(npr_y-1),*,*)
                	(*bridges[ipar])->SetVar, "slight", slight_block
                	(*bridges[ipar])->SetVar, "mlocal", mlocal
                endelse
            endif
            (*bridges[ipar])->SetVar, "y", y
            (*bridges[ipar])->SetVar, "iprof_x", iprof_x ;x block id
            (*bridges[ipar])->SetVar, "iprof_y", iprof_y ;y block id
            (*bridges[ipar])->Execute, 'milos_int_v5,Y,VARS,final_model,yfited,Err,chisqr,slight=slight,mlocal=mlocal',/nowait
;            (*bridges[ipar])->Execute, 'milos_int_v5,Y,VARS,final_model,yfited,Err,chisqr',/nowait

        print,'Processing block ',strtrim(string(iprof_x+iprof_y*nblock_x+1),2),$
          ', [x,y] = [',strtrim(string(iprof_x),2),',',strtrim(string(iprof_y),2),'], ',$
          ' of ', strtrim(string(nblock_x * nblock_y),2)

		iprof = iprof + 1
                iprof_x = iprof mod nblock_x
                iprof_y = iprof/nblock_x

		endif else print,'Resuming thread ',strtrim(string(ipar),2)
       endif
       if stats[ipar] eq 3 then begin
             print,'ERROR in ', ipar
       endif
endfor

if total(stats) eq 0 then running=0

endwhile ;END MAIL LOOP

message,/info, "Parallel computation done!"
message,/info, " Destroying bridges........."
wait,10
for ipar=0L,nbparallel-1 do print,(*bridges[ipar])->status()

if not(keyword_set(sav_file)) then begin
message,/info, " saving results in results.save"
save,filename='results.save',Block_data,final_model,errors,final_chisqr ;sin macro a fijo, New CI, constraits
;save,filename='results_perfiles.save',final_fit,axis,wt
    endif else begin
message,/info, " saving results in "+sav_file
save,filename=sav_file+'.save',Block_data,final_model,errors,final_chisqr ;sin macro a fijo, New CI, constraits
save,filename=sav_file+'_perfiles.save',final_fit,axis,wt
endelse

if not(keyword_set(sav_file)) then begin
message,/info, " saving results in results.save"
save,filename='results.save',Block_data,final_model,errors,final_chisqr ;sin macro a fijo, New CI, constraits
save,filename='results_perfiles.save',final_fit,axis,wt
    endif else begin
message,/info, " saving results in "+sav_file
save,filename=sav_file+'.save',Block_data,final_model,errors,final_chisqr ;sin macro a fijo, New CI, constraits
save,filename=sav_file+'_perfiles.save',final_fit,axis,wt
endelse

for ipar=0L,nbparallel-1 do begin
    stat = (*bridges[ipar])->status()
    if stat eq 0 then obj_destroy, (*bridges[ipar]) else print,ipar,' thread not destroyed'
    print, " .... Bridge ",ipar,' destroyed...'
    wait,5
endfor

ENDIF ELSE BEGIN

    nblock_x = ceil(npr_x*1./blocks) ;+ 1 ;number of blocks in x direction
    nblock_y = ceil(npr_y*1./blocks) ;+ 1 ;number of blocks in y direction
    iprof = 0L
    iprof_x = 0L
    iprof_y = 0L
    while (iprof lt nblock_x*nblock_y) do begin ;MAIN LOOP

            print,'Processing block ',strtrim(string(iprof_x+iprof_y*nblock_x+1),2),$
              ', [x,y] = [',strtrim(string(iprof_x),2),',',strtrim(string(iprof_y),2),'], ',$
              ' of ', strtrim(string(nblock_x * nblock_y),2)

            de_x = blocks*iprof_x<(npr_x-1)
            hasta_x = (blocks*(iprof_x+1)-1)<(npr_x-1)
            de_y = blocks*iprof_y<(npr_y-1)
            hasta_y = (blocks*(iprof_y+1)-1)<(npr_y-1)

            toop_x = (blocks*(iprof_x+1)-1)<(npr_x-1) - blocks*iprof_x
            toop_y = (blocks*(iprof_y+1)-1)<(npr_y-1) - blocks*iprof_y
	        ely = Block_data(de_x:hasta_x,de_y:hasta_y,*,*)
            ;en caso de que el muestreo sea diferente, este entra en muestreo
			if muestreo NE !NULL then begin
        	     y = fltarr(toop_x+1 , toop_y+1,landas,4)
        	    for ii = 0, n_landas-1 do y(*,*,muestreo(ii),*) = ely(*,*,ii,*)
           	 endif else y = ely

            if keyword_set(slight) then begin
                	slight_block = slight(blocks*iprof_x:(blocks*(iprof_x+1)-1)<(npr_x-1), blocks*iprof_y:(blocks*(iprof_y+1)-1)<(npr_y-1),*,*)
            endif

            ;milos_int_v5, wl, axis, init_model, fmodel, y, chisqr=chisqr,$
            ;  yfited=yfited,sigma=sigma,fix=fix,ilambda=ilambda,miter=miter,Err=Err,$
            ;  filter=filter,weight=wt,badpixel=badpixel,simple=simple
            milos_int_v5,Y,VARS,fmodel,yfited,Err,chisqr,slight=slight_block,mlocal=mlocal,/doplot
            print,'Task finished.',' Block ',$
              strtrim(string(iprof_x+iprof_y*nblock_x+1),2),$
              ', [x,y] = [',strtrim(string(iprof_x),2),',',$
              strtrim(string(iprof_y),2),'], ',' of ', $
              strtrim(string(nblock_x * nblock_y),2),'. Avr. Chi^2 VALUE: ',mean(chisqr)


            final_model(de_x:hasta_x,de_y:hasta_y,*) = fmodel
            errors(de_x:hasta_x,de_y:hasta_y,*) = err
            final_chisqr(de_x:hasta_x,de_y:hasta_y,*) = chisqr
            ;final_fit(de_x:hasta_x,de_y:hasta_y,*,*) = yfited

            iprof = iprof + 1   ;update block
            iprof_x = iprof mod nblock_x
            iprof_y = iprof/nblock_x

if not(keyword_set(sav_file)) then begin
message,/info, " saving results in results.save"
save,filename='results.save',Block_data,final_model,errors,final_chisqr ;sin macro a fijo, New CI, constraits
save,filename='results_perfiles.save',final_fit,axis,wt
    endif else begin
message,/info, " saving results in "+sav_file
save,filename=sav_file+'.save',Block_data,final_model,errors,final_chisqr ;sin macro a fijo, New CI, constraits
save,filename=sav_file+'_perfiles.save',final_fit,axis,wt
endelse

endwhile


ENDELSE

end
