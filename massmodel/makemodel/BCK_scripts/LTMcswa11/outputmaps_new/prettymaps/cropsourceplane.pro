pro cropsourceplane
file = file_search('source*.fits')
print,file
file_mkdir,'crop'
hadetect=readfits('sourceha_detection_pretty.fits')
bad = where(hadetect lt 5.)
for i=0,n_elements(file)-1 do begin
   img = readfits(file(i),hdr)
   img(bad) = 1./0. 
   imgnew = img(140:172,130:162)
 
      ;fix the velocity dispersion (subtract the intrinsic dispersion)
   if file(i) eq 'sourceveldisp.fits' then begin
      vel = imgnew(where(finite(imgnew)))
      velsq = vel^2-50.^2
      velsq(where(velsq lt 0.)) = 0.
      velnew = sqrt(velsq)
      imgnew(where(finite(imgnew))) = velnew
   endif
   if file(i) eq 'sourceha_pretty.fits' then imgnew(where(imgnew ge 0.09))=1./0.

   writefits,'crop/'+file(i),imgnew
   print, 'yay'


endfor
stop
end
