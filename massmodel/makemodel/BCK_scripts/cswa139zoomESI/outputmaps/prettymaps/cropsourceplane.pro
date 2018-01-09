pro cropsourceplane
file = file_search('source*.fits')
print,file
stop
file_mkdir,'crop'
for i=0,n_elements(file)-1 do begin
   img = readfits(file(i),hdr)
   imgnew = img(145:180,235:270)
      ;fix the velocity dispersion (subtract the intrinsic dispersion)
   if file(i) eq 'sourceveldisp_pretty.fits' then begin
      vel = imgnew(where(finite(imgnew)))
      velsq = vel^2-50.^2
      velsq(where(velsq lt 0.)) = 0.
      velnew = sqrt(velsq)
      imgnew(where(finite(imgnew))) = velnew
   endif

   writefits,'crop/'+file(i),imgnew
   print, 'yay'


endfor
stop
end
