pro fitstifconver
;Find the linear parameters for the conversion between fits file and tif files used in matlab relensing
tiffiles = file_search('*_test.tif')
fitsfiles = repstr(tiffiles,'_test.tif','.fits')
print, tiffiles,fitsfiles



;stop
setplot,14
!p.multi = [0,3,3]
!p.charsize = 1.5
window, 0, xsize = 1200, ysize=700
linvararr = fltarr(2,n_elements(tiffiles))
for ii=0,n_elements(tiffiles)-1 do begin
   tif = read_tiff(tiffiles(ii))
   fits = readfits(fitsfiles(ii))
   plot, tif,fits,psym=1,title=tiffiles(ii)

   x0 = mode(tif)
   x1 = max(tif)
   y0 = mode(fits)
   y1 = max(fits)
   print, x0,y0,x1,y1
   slope = (y1-y0)/(x1-x0)
   intercept = y0 - x0*slope
   linvar = [intercept,slope]
;   stop
;   linvar  = linfit(tif,fits)
   linvararr(*,ii) = linvar
   print,tiffiles(ii), linvar
   oplot, !x.crange,!x.crange*linvar(1)+linvar(0),color=20
   legend,string(linvar)  
   ;stop
endfor
save, tiffiles,linvararr,filename='fitstifconver.sav'
stop
end
