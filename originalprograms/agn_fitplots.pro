pro agn_fitplots,path=path,fits_file=fits_file,first,shift_x, shift_y,size_x,size_y,nameplot=nameplot,namegal=namegal,agnpointsxy,N2divHa,O3divHb,N2divHa_err,O3divHb_err,agnpointsinn2divha
  common BPTshare, goodagn
;Inputs
;Fits_file: input acube_obj.fits file for the save file that comes out of metallicity_map.pro Order = ha,nii,oiii,hb
;first : 1 or 2. whether the line is the first line (ran with fitspec_general) or the second line (ran with fitspec_intensity) Usually Ha and OIII are first lines while Hb and NII are second lines
;shift_x/y: x(y)_shift_ha/ or x(y)_shift hb in metallicity map.pro - The shift that would align Hb and Ha images
;nameplot : plot titles in the outputfiles
;AGNpointsxy: pixels with values indicated as AGN. The coordinates of the pixel are according to the intersected map of Ha and Hb.agnpointsxy(0,*) are x values while agnpointsxy(1,*) are y values.

set_plot,'ps'

;Test if a directory for this galaxy existed. If not, create a new directory name namegal_agntest
foldflag = file_test('/scr2/nichal/workspace/output/metallicity/'+namegal+'_agntest',/directory)
if foldflag eq 0 then file_mkdir,'/scr2/nichal/workspace/output/metallicity/'+namegal+'_agntest'

;Make data nice, shifting, and store all lines data
for ii=0,n_Elements(fits_file)-1 do begin
   input = fits_file(ii)
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

;Shift images to the output of metallicty map
   for xx=0,size_x-1 do begin for yy=0, size_y-1 do begin
      fitnew(xx,yy,*) = fitcube_obj(xx+xshift,yy+yshift,*)
      cubenew(xx,yy,*) = finalcube_obj(xx+xshift,yy+yshift,*)
      if first(ii) eq 1. then acubenew(xx,yy,*) = acube_obj(xx+xshift,yy+yshift,*)
      if first(ii) eq 2. then outputnew(xx,yy,*) = output_obj(xx+xshift,yy+yshift,*)
   endfor
endfor

if ii eq 0 and first(ii) eq 1. then begin 
   fitha = fitnew
   cubeha = cubenew
   fitsigha = acubenew(*,*,4)
endif else if ii eq 1 and first(ii) eq 2. then begin
   fitnii = fitnew
   cubenii = cubenew
   fitsignii = outputnew(*,*,2)
endif else if ii eq 2 and first(ii) eq 1.then begin
   fitoiii = fitnew
   cubeoiii = cubenew
   fitsigoiii = acubenew(*,*,4)
endif else if ii eq 3 and first(ii) eq 2.then begin
   fithb = fitnew
   cubehb = cubenew
   fitsighb = outputnew(*,*,2)
endif else begin
   print, 'boo'
   stop
endelse
endfor


;Do plots for every indicated AGN points

goodAGN = fltarr(n_elements(agnpointsxy)/2)

for ii = 0,(n_elements(agnpointsxy)/2)-1 do begin
   xpos = agnpointsxy(0,ii)
   ypos = agnpointsxy(1,ii)

   device, file='/scr2/nichal/workspace/output/metallicity/'+namegal+'_agntest/'+'x'+strtrim(xpos,2)+'y'+strtrim(ypos,2)+'.eps', encapsulated =1, /color
   device, xsize=30, ysize=18
   !p.multi = [0,2,2]

 
;2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Make graphs of the fits
plot, wl, cubeha(xpos,ypos,*),xtitle = 'wavelength',title='(x='+strtrim(xpos,2)+',y='+strtrim(ypos,2)+')'
oplot, wl, fitha(xpos,ypos,*),linestyle=2,thick=2
legend,['Ha data','fit'],linestyle=[0,2]
xyouts,3000,15400,'fit significance:'+string(fitsigha(xpos,ypos),format='(F5.2)'),/device

plot, wl, cubenii(xpos,ypos,*),xtitle = 'wavelength',title='(x='+strtrim(xpos,2)+',y='+strtrim(ypos,2)+')'
oplot, wl, fitnii(xpos,ypos,*),linestyle=2,thick=2
legend,['NII data','fit'],linestyle=[0,2]
xyouts,18000,15400,'fit significance:'+string(fitsignii(xpos,ypos),format='(F5.2)'),/device

plot, wl, cubeoiii(xpos,ypos,*),xtitle = 'wavelength',title='(x='+strtrim(xpos,2)+',y='+strtrim(ypos,2)+')'
oplot, wl, fitoiii(xpos,ypos,*),linestyle=2,thick=2
legend,['OIII data','fit'],linestyle=[0,2]
xyouts,3000,6400,'fit significance:'+string(fitsigoiii(xpos,ypos),format='(F5.2)'),/device

plot, wl, cubehb(xpos,ypos,*),xtitle = 'wavelength',title='(x='+strtrim(xpos,2)+',y='+strtrim(ypos,2)+')'
oplot, wl, fithb(xpos,ypos,*),linestyle=2,thick=2
legend,['Hb data','fit'],linestyle=[0,2]
xyouts,18000,6400,'fit significance:'+string(fitsighb(xpos,ypos),format='(F5.2)'),/device

logn2divhastr = string(alog10(n2divha(ii)),format='(F5.2)')
logo3divhbstr = string(alog10(o3divhb(ii)),format='(F5.2)')
stringbpt = 'log[NII/Ha]= '+logn2divhastr+'  log[OIII/Hb]= '+logo3divhbstr
xyouts,2,5,stringbpt,/device
device,/close

if finite(fitsignii(xpos,ypos)) eq 1 and finite(fitsighb(xpos,ypos)) eq 1 then goodagn(ii)= 1.

endfor

goodagn = where(goodagn eq 1)

;stop
end
