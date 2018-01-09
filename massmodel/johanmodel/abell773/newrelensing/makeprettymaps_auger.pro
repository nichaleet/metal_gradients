pro makeprettymaps_auger,img,outimg,xmin,xmax,ymin,ymax,minsurround
  outimg=img
  nanpix = where(finite(img) eq 0.)
  xyind = array_indices(img,nanpix)
  xind = xyind(0,*)
  yind = xyind(1,*)
  xind = reform(xind)
  yind = reform(yind)

  goodnanpix = where(xind ge xmin and xind le xmax and yind ge ymin and yind le ymax)
  xind = xind(goodnanpix)
  yind = yind(goodnanpix)
  countchange = 0.
  for jj=0, n_elements(xind)-1 do begin
     x = xind(jj)
     y = yind(jj)
     surroundpix = [img(x-1,y-1),img(x-1,y),img(x-1,y+1),img(x,y-1),img(x,y+1),img(x+1,y-1),img(x+1,y),img(x+1,y+1)]
     goodpix = where(finite(surroundpix) eq 1.)
     goodpixcount = n_elements(goodpix)      
     if goodpixcount ge minsurround then begin
        outimg(x,y) = mean(surroundpix(goodpix))
        countchange = countchange+1.
     endif
  endfor
  print, countchange  
end

