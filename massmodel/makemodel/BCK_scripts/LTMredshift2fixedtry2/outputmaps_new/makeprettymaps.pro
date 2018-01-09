pro makeprettymaps

file_mkdir,'prettymaps'

hadetect=readfits('sourceha_detection.fits')
baddetect = where(hadetect lt 10.)

; make the nan dots in source plane gone by averaging pixels around it
sourcefiles = file_search('source*.fits')
outname = repstr(sourcefiles,'.fits','_pretty.fits')
outname2 = repstr(sourcefiles,'.fits','_pretty_maingal.fits')
outname = 'prettymaps/'+outname
outname2 = 'prettymaps/'+outname2

for ii=0, n_Elements(sourcefiles)-1 do begin
   img = readfits(sourcefiles(ii),hdr)
   img(baddetect) = 1./0.
   img(where(abs(img) le 1.e-5)) = 1./0.
   if sourcefiles(ii) eq 'sourcekinematic.fits' then img(where(img eq 0)) = 1./0.
   nanpix = where(finite(img) eq 0.)
   xyind = array_indices(img,nanpix)
   xind = xyind(0,*)
   yind = xyind(1,*)
   xind = reform(xind)
   yind = reform(yind)

   goodnanpix = where(xind ge 684 and xind le 770 and yind ge 581 and yind le 666)
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
;crop 
   imgnew = img[680:795,572:682]
   writefits,outname(ii),imgnew,hdr

;extract main galaxy
   slope=(105.5-3.)/(9.45-96.)
   detection = readfits('prettymaps/sourceha_detection_pretty.fits')
   indices=array_indices(detection,indgen(n_elements(detection)))
   x=reform(indices(0,*))
   y=reform(indices(1,*))

   badpix = where(y gt slope*(x-9.75)+105.5 or detection lt 9.)
   imgnew(badpix)=1./0.
   imgnew = imgnew[0:85,0:85]
   writefits,outname2(ii),imgnew,hdr
endfor
end
