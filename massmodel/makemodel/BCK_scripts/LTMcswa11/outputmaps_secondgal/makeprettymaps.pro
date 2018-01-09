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
sourcefiles = file_search('source*.fits')
outname = repstr(sourcefiles,'.fits','_pretty.fits')
outname = 'prettymaps/'+outname
for ii=0, n_Elements(sourcefiles)-1 do begin
   
   img = readfits(sourcefiles(ii),hdr)
   if sourcefiles(ii) eq 'sourcekinematic.fits' then img(where(img eq 0.))=1./0.
   if sourcefiles(ii) eq 'sourcen2index.fits' then img(where(img eq 0.))=1./0.

   nanpix = where(finite(img) eq 0.)
   xyind = array_indices(img,nanpix)
   xind = xyind(0,*)
   yind = xyind(1,*)
   xind = reform(xind)
   yind = reform(yind)

   goodnanpix = where(xind ge 138 and xind le 174 and yind ge 126 and yind le 161)
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

    for jj=0, n_elements(xind)-1 do begin
      x = xind(jj)
      y = yind(jj)
      surroundpix = [img(x-1,y-1),img(x-1,y),img(x-1,y+1),img(x,y-1),img(x,y+1),img(x+1,y-1),img(x+1,y),img(x+1,y+1)]
      goodpix = where(finite(surroundpix) eq 1.)
      goodpixcount = n_elements(goodpix)
      
      if goodpixcount ge 5. then begin
         img(x,y) = mean(surroundpix(goodpix))
         countchange = countchange+1.
      endif
   endfor

   print, countchange, outname(ii)
   writefits,outname(ii),img,hdr
endfor
end
