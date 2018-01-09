pro cropsourceplanemetal
path = '/scr2/nichal/workspace/output/metallicity/'
files = ['Abell773_metallicity_sourceplane.fits','CSWA128_metallicitynew_sourceplane_interp.fits','CSWA165_metallicitynew_sourceplane_interp.fits','CSWA20_metallicitynew_sourceplane_interp.fits']


files = path+files

outfile = ['abell773.fits','cswa128.fits','cswa165.fits','cswa20.fits']
outfile_err = ['abell773err.fits','cswa128err.fits','cswa165err.fits','cswa20err.fits']
outfile_m08 = ['abell773m08.fits','cswa128m08.fits','cswa165m08.fits','cswa20m08.fits']
outfile_m08err = ['abell773m08err.fits','cswa128m08err.fits','cswa165m08err.fits','cswa20m08err.fits']
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

   imgvel = img(*,*,0)  ;NII/Ha
   imgvel = imgvel(x1:x2,y1:y2)
   ;Make it into square
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
;==============================
   imgvel = img(*,*,3) ;NII/Ha _err
   imgvel = imgvel(x1:x2,y1:y2)
   ;Make it into square
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
   writefits,outfile_err(i),imgvel
;=================================
   imgvel = img(*,*,2)   ;Bayesian metal
   imgvel = imgvel(x1:x2,y1:y2)
   ;Make it into square
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
   writefits,outfile_m08(i),imgvel
;===================================
   imgvel = img(*,*,5)
   imgvel = imgvel(x1:x2,y1:y2)
   ;Make it into square
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
   writefits,outfile_m08err(i),imgvel
endfor

end
