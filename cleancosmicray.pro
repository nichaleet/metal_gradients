pro cleancosmicray,threshold
;clean all the fits file in the folder. Make all the pixels above the threshold value to 0. and put all the files in a new folder called 'clean'
file_mkdir,'clean'
files = file_search('*.fits')
if n_elements(files) ne -1 then begin
   for ii=0,n_elements(files)-1 do begin
      cube = readfits(files(ii),header)
      badpix = where(abs(cube) ge threshold)
      cube(badpix) = 0.
      writefits,'clean/'+files(ii),cube,header
   endfor
endif
end
