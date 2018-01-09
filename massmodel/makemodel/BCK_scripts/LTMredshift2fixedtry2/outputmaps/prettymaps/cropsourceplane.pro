pro cropsourceplane
file = file_search('source*pretty.fits')
print,file
stop
file_mkdir,'crop'
for i=0,n_elements(file)-1 do begin
   img = readfits(file(i),hdr)
   imgnew = img(675:785,570:680)
   writefits,'crop/'+file(i),imgnew
   print, 'yay'
endfor
stop
end
