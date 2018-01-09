pro cropsourceplane
for ii=0,10 do begin
   file = file_search('cswa159_'+strtrim(string(fix(ii)),2)+'source*.fits')
   print,file
   hadetect=readfits('cswa159_'+strtrim(string(fix(ii)),2)+'sourceha_detection_pretty.fits')
   bad = where(hadetect lt 5.)
   for i=0,n_elements(file)-1 do begin
      img = readfits(file(i),hdr)
      img(bad)=1./0.
      imgnew = img(27:55,102:130)
                                ;fix the velocity dispersion (subtract the intrinsic dispersion)
      if file(i) eq 'cswa159_'+strtrim(string(fix(ii)),2)+'sourceveldisp_pretty.fits' then begin
         vel = imgnew(where(finite(imgnew)))
         velsq = vel^2-50.^2
         velsq(where(velsq lt 0.)) = 0.
         velnew = sqrt(velsq)
         imgnew(where(finite(imgnew))) = velnew
      endif
      
      writefits,'crop/'+file(i),imgnew
   endfor
endfor
end
