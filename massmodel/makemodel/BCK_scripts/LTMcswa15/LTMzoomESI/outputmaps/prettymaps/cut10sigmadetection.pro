pro cut10sigmadetection
;find where sigma detection is less than 10

detection = readfits('crop/sourceha_detection_pretty.fits')
badpix = where(detection le 10. or finite(detection) eq 0.)


file = file_search('crop/source*.fits')
print,file
stop
file_mkdir,'crop/10sigmadetection'

for i=0,n_elements(file)-1 do begin
   img = readfits(file(i),hdr)
   img(badpix)=1./0.
   namefile = strsplit(file(i),'/',/extract)
   namefile=namefile(1)
   writefits,'crop/10sigmadetection/'+namefile,img
   print, 'finished cutting file:',file(i)
endfor
stop
end
