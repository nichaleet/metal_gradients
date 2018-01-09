pro cropsourceplanevel
path = '/scr2/nichal/workspace/output/'
files = ['abell773_Ha_Kc3_mosaic_sky_230hr_acube_sourceplane.fits','cswa128_Ha_Kn2_mosaic_sky_130hr_acube_sourceplane_interp.fits','cswa165_Ha_Kn2_mosaic_sky_1hr_acube_sourceplane_interp.fits','cswa20_Ha_Hn2_mosaic_scaledsky_100hr_acube_sourceplane_interp.fits']
files = path+files

outfile = ['abell773.fits','cswa128.fits','cswa165.fits','cswa20.fits']
outfile_disp = ['abell773_disp.fits','cswa128_disp.fits','cswa165_disp.fits','cswa20_disp.fits']
outfile_velerr = ['abell773_velerr.fits','cswa128_velerr.fits','cswa165_velerr.fits','cswa20_velerr.fits']
outfile_Ha = ['abell773_Ha.fits','cswa128_Ha.fits','cswa165_Ha.fits','cswa20_Ha.fits']

;x1 = [136,233,4,5]
;x2 = [176,257,30,31]
;y1 = [18,223,0,8]
;y2 = [30,249,21,45]
xc = [156,245,17,18]
yc = [24,236,10,26]
width = [40,26,26,36] 
hwidth = width/2

for i=0,n_elements(files)-1 do begin
   img = readfits(files(i))
   size_file = size(img)

   x1 = xc(i)-hwidth(i)
   x2 = xc(i)+hwidth(i)
   y1 = yc(i)-hwidth(i)
   y2 = yc(i)+hwidth(i)

   if x1 lt 0 then x1=0
   if y1 lt 0 then y1=0
   if x2 ge size_file(1) then x2 = size_file(1)-1
   if y2 ge size_file(2) then y2 = size_file(2)-1

   imgvel = img(*,*,0)
   imgvel = imgvel(x1:x2,y1:y2)
   
   newsize = size(imgvel)
   xappend = (width(i)-newsize(1))/2
   if xappend ge 1 then begin
     xblank = fltarr(xappend,newsize(2))+1./0.
     imgvel =[xblank,imgvel,xblank]
  endif

   newsize = size(imgvel)
   yappend = (width(i)-newsize(2))/2
   if yappend ge 1 then begin
      yblank = fltarr(newsize(1),yappend)+1./0.
      imgvel = [[yblank],[imgvel],[yblank]]
   endif
   help,imgvel
   writefits,outfile(i),imgvel

; Now doing velocity dispersion
   imgdisp = img(*,*,2)
   imgdisp = imgdisp(x1:x2,y1:y2)
   
   newsize = size(imgdisp)
   xappend = (width(i)-newsize(1))/2
   if xappend ge 1 then begin
     xblank = fltarr(xappend,newsize(2))+1./0.
     imgdisp =[xblank,imgdisp,xblank]
  endif

   newsize = size(imgdisp)
   yappend = (width(i)-newsize(2))/2
   if yappend ge 1 then begin
      yblank = fltarr(newsize(1),yappend)+1./0.
      imgdisp = [[yblank],[imgdisp],[yblank]]
   endif
   writefits,outfile_disp(i),imgdisp

; Now doing velocity uncertainty
   imgvelerr = img(*,*,8)
   imgvelerr = imgvelerr(x1:x2,y1:y2)
   
   newsize = size(imgvelerr)
   xappend = (width(i)-newsize(1))/2
   if xappend ge 1 then begin
     xblank = fltarr(xappend,newsize(2))+1./0.
     imgvelerr =[xblank,imgvelerr,xblank]
  endif

   newsize = size(imgvelerr)
   yappend = (width(i)-newsize(2))/2
   if yappend ge 1 then begin
      yblank = fltarr(newsize(1),yappend)+1./0.
      imgvelerr = [[yblank],[imgvelerr],[yblank]]
   endif
   writefits,outfile_velerr(i),imgvelerr

; Now doing Halpha
   imgHa = img(*,*,1)
   imgHa = imgHa(x1:x2,y1:y2)
   
   newsize = size(imgHa)
   xappend = (width(i)-newsize(1))/2
   if xappend ge 1 then begin
     xblank = fltarr(xappend,newsize(2))+1./0.
     imgHa =[xblank,imgHa,xblank]
  endif

   newsize = size(imgHa)
   yappend = (width(i)-newsize(2))/2
   if yappend ge 1 then begin
      yblank = fltarr(newsize(1),yappend)+1./0.
      imgHa = [[yblank],[imgHa],[yblank]]
   endif
   writefits,outfile_Ha(i),imgHa

endfor

end
