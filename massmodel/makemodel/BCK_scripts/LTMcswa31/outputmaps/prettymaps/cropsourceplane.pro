
pro cropsourceplane
file = file_search('source*.fits')
print,file
stop
file_mkdir,'crop'
imgha = readfits('sourceha_pretty.fits')
imgha = imgha(90:125,80:115)
badpix = where(finite(imgha) eq 0.)

for i=0,n_elements(file)-1 do begin
   img = readfits(file(i),hdr)
   imgnew = img(90:125,80:115)
   imgnew(where(imgnew eq 0)) = 1./0.
   imgnew(badpix) = 1./0.
      ;fix the velocity dispersion (subtract the intrinsic dispersion)
   if file(i) eq 'sourceveldisp.fits' then begin
      vel = imgnew(where(finite(imgnew)))
      velsq = vel^2-50.^2
      velsq(where(velsq lt 0.)) = 0.
      velnew = sqrt(velsq)
      imgnew(where(finite(imgnew))) = velnew
   endif
   ;if file(i) eq 'sourcekinematic_pretty.fits' then imgnew(where(imgnew ge 0))=1./0.

   writefits,'crop/'+file(i),imgnew
   print, 'yay'


endfor
stop
end
