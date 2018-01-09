pro cropsourceplane
file_mkdir,'crop'

for i=0,19 do begin
hadetect = readfits('sourcehadetection_random'+strtrim(string(i),2)+'_pretty.fits')
sizefits = size(hadetect,/dimensions)
xind     = rebin(indgen(sizefits(0)),sizefits(0),sizefits(1))
yind     = rebin(transpose(indgen(sizefits(1))),sizefits(0),sizefits(1))
good     = where(hadetect ge 5. and finite(hadetect))
xc       = median(xind(good))
yc       = median(yind(good))
bad      = where(hadetect lt 5. and finite(hadetect))

files    = ['sourcehadetection_random'+strtrim(string(i),2)+'_pretty.fits','sourcekinematic_random'+strtrim(string(i),2)+'_pretty.fits','sourcekinematic_err_random'+strtrim(string(i),2)+'_pretty.fits']
  for ii=0,n_elements(files)-1 do begin
     img = readfits(files(ii),hdr)
     img(bad)=1./0.
     imgnew = img(xc-23:xc+23,yc-23:yc+23)
     imgnew[*,0:16] = 1./0.
     imgnew[*,30:46] = 1./0.
     writefits,'crop/'+files(ii),imgnew
  endfor
endfor

stop
end
