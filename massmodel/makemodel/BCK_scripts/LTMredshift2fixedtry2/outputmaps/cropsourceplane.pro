pro cropsourceplane
file = file_search('source*pretty.fits')
print,file
stop
file_mkdir,'crop'
for i=0,6 do begin
   img = readfits(file(i),hdr)
   imgnew = img(680:777,573:680)
   writefits,'crop/'+file(i),imgnew
   print, 'yay'
endfor
stop
end
