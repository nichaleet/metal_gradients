pro uvlum_j1148
  z             = 2.38
  magnification = 10.3
  magnification_err= 5
  cat = rsex('j1148.cat')
  pos1 = where(cat.number eq 1238)
  
  
  mag      = cat(pos1).mag_iso
  magerr   = cat(pos1).magerr_iso

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
