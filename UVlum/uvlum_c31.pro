pro uvlum_c31
  b = 1.4e-10
  z = 1.49
  magnification = 3.3
  im   = readfits('cswa31sdss_u.fits')
  mask = readfits('cswa31sdss_u_mask.fits')
  background = readfits('check_c31_background.fits')

  im_sub = im-background
  
  goodpix = where(mask eq 1. and im_sub gt 0.)
  
  flux_iso = total(im_sub(goodpix))
  flux_iso_err = sqrt(total(im_sub(goodpix)))

  mag_std = 22.5-2.5*alog10(flux_iso)
  mag_sdss = -2.5/alog(10)*(asinh(flux_iso*1.e-9/(2.*b))+alog(b))
  errmag = 2.5*abs(flux_iso_err/flux_iso/alog(10.))

  if abs(mag_std-mag_sdss) gt 0.1 then stop,'You need to really calibrate the magnitude. Poor you' else mag = 0.5*(mag_std+mag_sdss)

  print, 'magnitude of the arc is ',mag
  if 0.5*abs(mag_std-mag_sdss) gt errmag then errmag = 0.5*abs(mag_std-mag_sdss)

  mag_source = mag+2.5*alog10(magnification)
  print, 'magnitude of the source is ', mag_source

;Calculate absolute magnitude
  dist = lumdist(z)*1.e6        ;parsecs
  UVabmag = mag_source-5.*(alog10(dist)-1.)+2.5*alog10(1.+z)
  print,'UV Absolute mag',Uvabmag,errmag
  stop
end
