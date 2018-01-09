pro makeprettymaps

file_mkdir,'prettymaps'
;make the image plane pics with 0 pixels to nan
imagefiles = file_search('*interp.fits')
outname = repstr(imagefiles,'.fits','_pretty.fits')
outname = 'prettymaps/'+outname

for ii=0, n_elements(imagefiles)-1 do begin
   img = readfits(imagefiles(ii),hdr)
   img(where(img eq 0.)) = 1./0.
   writefits,outname(ii),img,hdr
endfor

; make the nan dots in source plane gone by averaging pixels around it
;remove the detection less than 5.
sourcehadetection = readfits('sourceha_detection.fits')
badpix = where(sourcehadetection lt 5.)

sourcefiles = file_search('source*.fits')
outname = repstr(sourcefiles,'.fits','_pretty.fits')
outname = 'prettymaps/'+outname
for ii=0, n_Elements(sourcefiles)-1 do begin
   img = readfits(sourcefiles(ii),hdr)
   img(where(abs(img) le 1.e-5)) = 1./0.
   img(badpix) = 1./0.
   if sourcefiles(ii) eq 'sourcekinematic.fits' then begin
      img(where(img eq 0)) = 1./0.
      img = img+45.
   endif
   nanpix = where(finite(img) eq 0.)
   xyind = array_indices(img,nanpix)
   xind = xyind(0,*)
   yind = xyind(1,*)
   xind = reform(xind)
   yind = reform(yind)

   goodnanpix = where(xind ge 95 and xind le 116 and yind ge 88 and yind le 109)
   xind = xind(goodnanpix)
   yind = yind(goodnanpix)
   countchange = 0.
   for jj=0, n_elements(xind)-1 do begin
      x = xind(jj)
      y = yind(jj)
      surroundpix = [img(x-1,y-1),img(x-1,y),img(x-1,y+1),img(x,y-1),img(x,y+1),img(x+1,y-1),img(x+1,y),img(x+1,y+1)]
      goodpix = where(finite(surroundpix) eq 1.)
      goodpixcount = n_elements(goodpix)
      
      if goodpixcount ge 6. then begin
         img(x,y) = mean(surroundpix(goodpix))
         countchange = countchange+1.
      endif
   endfor
   print, countchange, outname(ii)
   img = img[95:118,88:111]
   ;rebin into 2" per pixel
   img = rebin(img,12,12)
   writefits,outname(ii),img,hdr

   
endfor
end
