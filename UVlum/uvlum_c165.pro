pro uvlum_c165
magnification = 42.
z = 2.13
;Note on changing flux to SDSS magnitude
;https://www.sdss3.org/dr10/algorithms/magnitudes.php#asinh
;sdss parameters
b = 0.9e-10


;load sdss uband catalog
pos1=569
cat = rsex('cswa165_g.cat')
if abs(cat(pos1).x_image-1328) gt 10. then stop  
flux_iso = cat(pos1).flux_iso
flux_iso_err = cat(pos1).fluxerr_iso

errmag = 2.5*abs(flux_iso_err/flux_iso/alog(10.))
mag_std  = 22.5-2.5*alog10(flux_iso)
mag_sdss = -2.5/alog(10)*(asinh(flux_iso*1.e-9/(2.*b))+alog(b))

if abs(mag_std-mag_sdss) gt 0.1 then stop,'You need to really calibrate the magnitude. Poor you' else mag = 0.5*(mag_std+mag_sdss)

print, 'magnitude of the arc is ',mag
if 0.5*abs(mag_std-mag_sdss) gt errmag then errmag = 0.5*abs(mag_std-mag_sdss)

mag_source = mag+2.5*alog10(magnification)
print, 'magnitude of the source is ', mag_source

;Calculate absolute magnitude
dist = lumdist(z)*1.e6 ;parsecs
UVabmag = mag_source-5.*(alog10(dist)-1.)+2.5*alog10(1.+z)
print,'UV Absolute mag',Uvabmag,errmag
stop


stop
end
