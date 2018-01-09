pro cutmaingalaxy
;Extract only the main (bigger) galaxy by a straight line

slope=(105.5-3.)/(9.45-96.)


detection = readfits('crop/sourceha_detection_pretty.fits')
indices=array_indices(detection,indgen(n_elements(detection)))
x=reform(indices(0,*))
y=reform(indices(1,*))

badpix = where(y gt slope*(x-9.75)+105.5 or detection lt 9.)


file = file_search('crop/source*.fits')
print,file
stop
file_mkdir,'crop/maingalaxy'

for i=0,n_elements(file)-1 do begin
   img = readfits(file(i),hdr)
   img(badpix)=1./0.
   img = img[0:85,0:85]
   namefile = strsplit(file(i),'/',/extract)
   namefile=namefile(1)
   writefits,'crop/maingalaxy/'+namefile,img
   print, 'finished cutting file:',file(i)
endfor
stop
end
