pro cut5sigmadetection

detection = readfits('crop/sourceha_detection_pretty.fits')

badpix = where(detection lt 5. or finite(detection) eq 0)

file = file_search('crop/source*.fits')
print,file
stop
file_mkdir,'crop/5sigmadetection'

for i=0,n_elements(file)-1 do begin
   img = readfits(file(i),hdr)
   img(badpix)=1./0.
   namefile = strsplit(file(i),'/',/extract)
   namefile=namefile(1)
   writefits,'crop/5sigmadetection/'+namefile,img
   print, 'finished cutting file:',file(i)
endfor
stop
end
