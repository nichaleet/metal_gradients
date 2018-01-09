pro uvlum_c159
;CSWA159
;SDSS Gband = 22.22 +/- 0.1 ;SDSS10
z = 2.30
UVmag = 22.22
UVmagerr = 0.1
;divide by magnification
magnification = 6.2
magnificationerr = 0.2
UVmag_delensed = UVmag+2.5*alog10(magnification)
UVmagerr_delensed = sqrt((2.5*0.4343*magnificationerr/magnification)^2+UVmagerr^2)
;Calculate absolute magnitude
dist = lumdist(z)*1.e6 ;parsecs
UVabmag = UVmag_delensed-5.*(alog10(dist)-1)+2.5*alog10(1.+z)
print,'UV Ab mag',Uvabmag,UVmagerr_delensed
stop
end
