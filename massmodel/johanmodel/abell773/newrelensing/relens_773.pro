pro relens_773,imagefile

outputfile = repstr(imagefile,'.fits','_newsourceplane.fits')
outputfile_interp = repstr(imagefile,'.fits','_newsourceplane_interp.fits')
image = readfits(imagefile,header_im)
imsize = size(image)
model = readfits('a773_model.fits',header_mo)
Ha_im = model(*,*,0)
x_source = model(*,*,1) 
y_source = model(*,*,2)

minx = min(x_source(where(finite(x_source))))
miny = min(y_source(where(finite(y_source))))

res = 0.04
print, 'resolution',1./res
;use the same resolution in both direction
xsize = round((max(x_source(where(finite(x_Source))))-minx)/res)
ysize = round((max(y_source(where(finite(y_Source))))-miny)/res)

print, 'new source plane size ',xsize,ysize

stop
newim = fltarr(xsize,ysize,imsize(3))
newim = (newim+1.)/0.  ;nan array

for ii=0,imsize(3)-1 do begin
   countim = fltarr(xsize,ysize)   
   current_im = image(*,*,ii)
   goodpixels = where(finite(current_im))
   ind = array_indices(current_im,goodpixels)
   xind_im = reform(ind(0,*))  ;xindices of image plane
   yind_im = reform(ind(1,*))
   newxarr = []
   newyarr = []
   for jj=0, n_Elements(xind_im)-1 do begin
      newx = model(xind_im(jj),yind_im(jj),1)    ;xindices of source plane in auger's map
      newy = model(xind_im(jj),yind_im(jj),2)    ;yindices of source plane in auger's map
      if finite(newx) eq 0 or finite(newy) eq 0 then goto, skip
      newx = fix((newx-minx)/res)
      newy = fix((newy-miny)/res)
                                ;if finite(newim(newx,newy,ii)) eq 0 then newim(newx,newy,ii) = current_im(xind_im(jj),yind_im(jj)) else 
      countim(newx,newy) = countim(newx,newy)+1.
      ncount = countim(newx,newy)
      if ncount eq 1. then newim(newx,newy,ii) = current_im(xind_im(jj),yind_im(jj))
      if ncount gt 1. and finite(newim(newx,newy,ii)) eq 1  then begin
         newim(newx,newy,ii) =  (current_im(xind_im(jj),yind_im(jj))+(ncount-1.)*newim(newx,newy,ii))/ncount
         print,'the image plane pixels were map to the same source plane pixel'
      endif

      newxarr = [newxarr,newx]
      newyarr = [newyarr,newy]
      skip: continue
   endfor
endfor

writefits, outputfile,newim,header_im


newim_interp = newim
;interpolation
for ii=0, imsize(3)-1 do begin
   current_im = reform(newim(*,*,ii))
   makeprettymaps_Auger,current_im,outimg,233,258,233,248,4
   current_im = outimg
   makeprettymaps_Auger,current_im,outimg,233,258,233,248,5
   newim_interp(*,*,ii) = outimg
endfor

writefits,outputfile_interp,newim_interp, header_im
;stop
end
