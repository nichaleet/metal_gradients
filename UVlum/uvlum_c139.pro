pro uvlum_c139
;CSWA139
;SDSS Gband = 22.709 +/- 0.08 ;SDSS10
z = 2.54
UVmag = 23.29
UVmagerr = 0.24
;divide by magnification
magnification = 8.7
magnificationerr = 0.2
UVmag_delensed = UVmag+2.5*alog10(magnification)
UVmagerr_delensed = sqrt((2.5*0.4343*magnificationerr/magnification)^2+UVmagerr^2)
;Calculate absolute magnitude
dist = lumdist(z)*1.e6 ;parsecs
UVabmag = UVmag_delensed-5.*(alog10(dist)-1)+2.5*alog10(1.+z)
print,'UV Ab mag',Uvabmag,UVmagerr_delensed
stop
end
