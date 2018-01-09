pro magnification165_area
model = readfits('cswa165_SDSSJ0105+0144_deflections.fits',header_mo)
model(where(model eq 0.)) = 1./0. 
x_source = model(*,*,1)
y_source = model(*,*,2)

minx = min(x_source(where(x_source ne 0.)))
miny = min(y_source(where(y_source ne 0.)))
res = .125   
print, 'resolution',1./res   ;times of the original osiris scale
;use the same resolution in both direction
xsize = round((max(x_source(where(finite(x_Source))))-minx)/res)
ysize = round((max(y_source(where(finite(y_Source))))-miny)/res)

print, 'new source plane size ',xsize,ysize

newim = fltarr(xsize,ysize)+1./0.

current_im = model[*,*,0]
goodpixels = where(finite(current_im))
ind = array_indices(current_im,goodpixels)
xind_im = reform(ind(0,*))      ;xindices of image plane
yind_im = reform(ind(1,*))
for jj =0, n_elements(xind_im)-1 do begin
   newx = model(xind_im(jj),yind_im(jj),1) ;xindices of source plane in auger's map
   newy = model(xind_im(jj),yind_im(jj),2) ;yindices of source plane in auger's map
   if finite(newx) eq 0 or finite(newy) eq 0 then goto, skip
   newx = fix((newx-minx)/res)
   newy = fix((newy-miny)/res)
   newim(newx,newy)= model(xind_im(jj),yind_im(jj),0)
   skip: continue
endfor
writefits,'test_before.fits',newim
makeprettymaps_Auger,newim,outimg,0,29,1,21,4
newim = outimg
makeprettymaps_Auger,newim,outimg,0,29,1,21,4
newim = outimg
writefits,'test.fits',outimg

  ;Magnification in model region
  map = readfits('cswa165_SDSSJ0105+0144_deflections.fits')
  mapx = map[*,*,1]
  mapy = map[*,*,2]
  remove,where(mapx eq 0.),mapx,mapy
  plot,mapx,mapy,psym=1,xrange=[min(mapx)-2,max(mapx)+2],yrange=[min(mapy)-2,max(mapy)+2]
  
  im = map[*,*,1]
  im_detect = where(im ne 0.,areaim)
  areaim = float(areaim)

  source = outimg
  source_detect = where(finite(source),areasource)
  areasource = float(areasource)
  print,areaim, areasource
  mag = areaim/areasource*8.^2
  print, 'model magnification by area', mag

  ;Magnification in Ha detection region
  im=readfits('/scr2/nichal/workspace/output/cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr_acube.fits') ;
  imdetect = im[*,*,4]
  goodimpix = where(finite(imdetect),ngoodimpix)
  source = readfits('/scr2/nichal/workspace/output/cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr_acube_sourceplane_interp.fits');
  sourcedetect = source[*,*,4]
  goodimsource = where(finite(sourcedetect),ngoodimsource)
  ;print, ngoodimpix,ngoodimsource
  mag = float(ngoodimpix)/float(ngoodimsource)*8.^2
  print,'magnification by area where Ha is detected:',mag

    ;flux weighted magnification
  imHa = im[*,*,1]
  Flux_im = total(imHa(where(finite(imHa))))*0.1^2
  sourceHa = source[*,*,1]
  Flux_source = total(sourceHa(where(finite(sourceHa))))*0.0125^2
  mag = flux_im/flux_source
  print,'magnifcation by Ha flux:', mag

     ;flux weighted magnification2
  map = readfits('cswa165_SDSSJ0105+0144_deflections.fits')
  imHa = map[*,*,0]
  imHa = imHa[22:60,30:45]
  Flux_im = total(imHa(where(finite(imHa))))*0.1^2
  sourceHa = outimg
  Flux_source = total(sourceHa(where(finite(sourceHa))))*0.0125^2
  mag = flux_im/flux_source
  print,'magnifcation by Ha flux:', mag

stop
end
