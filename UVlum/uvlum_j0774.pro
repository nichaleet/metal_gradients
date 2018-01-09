pro uvlum_j0774
  z             = 2.21
  magnification = 16.
  magnification_err=3.
  cat = rsex('J0774.cat')
  pos1 = where(cat.number eq 847)
  pos2 = where(cat.number eq 855)
  
  totflux    = total(cat(pos1).flux_iso+cat(pos2).flux_iso)
  totfluxerr = sqrt(total((cat(pos1).fluxerr_iso)^2+(cat(pos2).fluxerr_iso)^2))
  
  mag_ref  = cat(pos1).mag_iso
  flux_ref = cat(pos1).flux_iso

  mag      = mag_ref-2.5*alog10(totflux/flux_ref)
  magerr   = 2.5*0.4343*totfluxerr/totflux

  print, 'magnitude of the arc is ',mag,magerr

  mag_source = mag+2.5*alog10(magnification)
  magerr = sqrt(magerr^2+(2.5*magnification_err/magnification*0.4343)^2)
  print, 'magnitude of the source is ', mag_source,magerr

;Calculate absolute magnitude
  dist = lumdist(z)*1.e6        ;parsecs
  UVabmag = mag_source-5.*(alog10(dist)-1.)+2.5*alog10(1.+z)
  print,'UV Absolute mag',Uvabmag,magerr
  stop
end
