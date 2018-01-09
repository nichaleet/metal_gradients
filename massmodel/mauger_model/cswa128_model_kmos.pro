pro cswa128_model_kmos
model = readfits('cswa128_SDSSJ1956+5950_deflections.fits',hdr)
model = model[*,0:73,*]
newmodel = fltarr(33,37,5)
for i=0,4 do  newmodel[*,*,i] = rebin(reform(model[*,*,i]),33,37)
writefits,'cswa128_deflection_kmos.fits',newmodel

imagefile='/scr2/nichal/workspace/output/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_KMOS_acube.fits'
outputfile = repstr(imagefile,'.fits','_sourceplane.fits')
outputfile_interp = repstr(imagefile,'.fits','_sourceplane_interp.fits')
image = readfits(imagefile,header_im)
imsize = size(image)
model = newmodel
model(where(model eq 0.)) = 1./0. 
x_source = model(*,*,1)
y_source = model(*,*,2)

sort_ind_x = bsort(x_Source,x_source_sort)
sort_ind_y = bsort(y_source,y_Source_sort)
x_diff_arr = fltarr(n_elements(x_source_sort)-1)
y_diff_arr = fltarr(n_elements(y_source_sort)-1)

for ii=0,n_elements(x_diff_arr)-1 do x_Diff_arr(ii) = x_source_sort(ii+1)-(x_source_sort(ii))
for ii=0,n_elements(y_diff_arr)-1 do y_Diff_arr(ii) = y_source_sort(ii+1)-(y_source_sort(ii))

x_diff_arr  = x_diff_arr(where(x_diff_arr ne 0. and finite(x_diff_arr) eq 1))
y_diff_arr  = y_diff_arr(where(y_diff_arr ne 0. and finite(x_diff_arr) eq 1))

statx = summary(x_diff_arr) ; min, 1st quartile, median, 3rd quartile, max
staty = summary(y_diff_arr)

res_x = statx(2)
res_y = staty(2)

minx = min(x_source(where(x_source ne 0.)))
miny = min(y_source(where(y_source ne 0.)))
res = min(res_x,res_y)
res = 1.
print, 'resolution',1./res
;use the same resolution in both direction
xsize = round((max(x_source(where(finite(x_Source))))-minx)/res)+1
ysize = round((max(y_source(where(finite(y_Source))))-miny)/res)+1

print, 'new source plane size ',xsize,ysize

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
   makeprettymaps_Auger,current_im,outimg,92,104,86,98,5
   newim_interp(*,*,ii) = outimg
endfor

writefits,outputfile_interp,newim_interp, header_im
stop
end

stop
end
