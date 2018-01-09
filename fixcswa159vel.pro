pro fixcswa159vel
;somehow the interpolation of the velocity map make an excessive rim
files= file_Search('/scr2/nichal/workspace/output/forlensmodel/CSWA159/*_interp.fits')
imref= readfits('/scr2/nichal/workspace/output/forlensmodel/CSWA159/CSWA159_kinematic_interp.fits')
ind = array_indices(imref,findgen(n_elements(imref)))
y = ind(1,*)
bad1 = where(imref gt -35. and y ge 168)
bad2 = where(abs(imref) le 12. and y lt 168)
bad=[bad1,bad2]
for n=0,n_Elements(files)-1 do begin
   im = readfits(files(n),hdr)
   im(bad) = 0.
   writefits,files(n),im,hdr
endfor

end
