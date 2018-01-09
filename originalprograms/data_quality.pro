pro data_quality, path=path,fits_file=fits_file,first,shift_x, shift_y,size_x,size_y,nameplot=nameplot,namegal=namegal

; input can be arrays
;Input
;Fits_file: input acube_obj.fits file for the save file that comes out of metallicity_map.pro
;first : 1 or 2. whether the line is the first line (ran with fitspec_general) or the second line (ran with fitspec_intensity) Usually Ha and OIII are first lines while Hb and NII are second lines
;shift_x/y: x(y)_shift_ha/ or x(y)_shift hb in metallicity map.pro - The shift that would align Hb and Ha images
;nameplot : plot titles in the outputfiles

set_plot,'ps'
device, file='output/metallicity/'+namegal+'_dataqual.eps', encapsulated =1, /color
device, xsize=40, ysize=18
!p.multi = [0,4,2]

for ii=0,n_Elements(fits_file)-1 do begin
input = fits_file(ii)
if input eq 'no' then goto, skip
xshift = shift_x(ii)
yshift = shift_y(ii)
name_pos = strpos(input,'/',/reverse_search)
name = strmid(input, name_pos+1)

if first(ii) eq 1 then begin
   name = strmid(name, 0,strlen(name)-15) ;delete '_acube_obj.fits'
   savefile = path+name+'.sav'
endif else if first(ii) eq 2 then begin
   name = strmid(name, 0,strlen(name)-26) ;delete '_secondline_acube_obj.fits'
   savefile = path+name+'_2ndline.sav'
endif else print,'something is wrong with your input ''first'' parameter'

restore, savefile
sz = size(fitcube_obj)
fitnew = findgen(size_x,size_y,sz(3))
cubenew = findgen(size_x,size_y,sz(3))
acubenew = findgen(size_x,size_y,6)
outputnew = findgen(size_x,size_y,3)

;1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Shift images to the output of metallicty map
for xx=0,size_x-1 do begin for yy=0, size_y-1 do begin
fitnew(xx,yy,*) = fitcube_obj(xx+xshift,yy+yshift,*)
cubenew(xx,yy,*) = finalcube_obj(xx+xshift,yy+yshift,*)
if first(ii) eq 1. then acubenew(xx,yy,*) = acube_obj(xx+xshift,yy+yshift,*)
if first(ii) eq 2. then outputnew(xx,yy,*) = output_obj(xx+xshift,yy+yshift,*)
endfor
endfor

;2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Make graphs of total signals over all x,y for each wavelengths
sumfit = fltarr(sz(3))
sumcube = fltarr(sz(3))

for jj=0,sz(3)-1 do begin
   sumfit(jj) = total(fitnew(*,*,jj))
   sumcube(jj) = total(cubenew(*,*,jj))
endfor

plot,wl,sumcube,xtitle='wavelength',ymin=0.,title=nameplot(ii)
oplot,wl,sumfit,linestyle=2,thick=2
legend,['data','fit'],linestyle=[0,2]
;3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Make histogram of fit significance
if first(ii) eq 1. then begin
   acubenew(where(acubenew(*,*,4) eq 0.)) = 1./0.
   plothist, acubenew(*,*,4), xtitle='fit significance',/nan
endif

if first(ii) eq 2. then begin
   outputnew(where(outputnew(*,*,2) eq 0.)) = 1./0.
   if where(finite(outputnew(*,*,2)) ne [-1]) then plothist, outputnew(*,*,2),xtitle='fit significance',/nan else begin
      blankx = fltarr(10)
      blanky = fltarr(10)
      plot,blankx,blanky,title='All values are not significant (nsig = Nan)'
   endelse
; do not plot histogram if all values are NAN and plot blank plot instead
endif
skip: continue
endfor

device,/close

;stop
end
