pro positive_starspec

files=file_search('*.fits')
for nfile=0,n_elements(files)-1 do begin
   file = files(nfile)
   spec = readfits(file,hdr)
   writefits,file,abs(spec),hdr
   plot,abs(spec)
   wait,1.5
endfor

end
