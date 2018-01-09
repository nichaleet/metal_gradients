pro uvlum_c15

clight = 3e10

;CSWA15
z=2.16
lambda = 1500.e-8
catfile = 'cswa15_B.cat'
imfile = 'cswa15_B_small_match_osirisRES.fits'
catfits_file = 'check_c15.fits'
magfile = 'cswa15_magnification.fits'
backgroundfile = 'check_c15_background.fits'

;read files
cat = rsex(catfile)
catfits = readfits(catfits_file)
im = readfits(imfile)
mag = readfits(magfile)
background = readfits(backgroundfile)

id = 54
sex_roi = where(catfits eq id)  ;get the region

;get the values from catalog
pos = where(cat.number eq id,cpos)
if cpos ne 1 then stop
sex_flux = cat[pos].flux_iso
sex_fluxerr = cat[pos].fluxerr_iso
sex_mag = cat[pos].mag_auto 

;demagnification
sex_delensed_flux_tot = total((im(sex_roi)-background(sex_roi))/mag(sex_roi))
sex_delensed_fluxerr_tot = sex_fluxerr/total(mag(sex_roi))

print,'Delensed fluxerr (sex) = ',sex_delensed_fluxerr_tot
UVmag = -2.5*alog10(sex_delensed_flux_tot/sex_flux)+sex_mag
UVmagerr = 2.5*0.4343*sex_delensed_fluxerr_tot/sex_delensed_Flux_tot
print,'UV mag', UVmag,UVmagerr

;Calculate absolute magnitude
dist = lumdist(z)*1.e6 ;parsecs
UVabmag = UVmag-5.*(alog10(dist)-1.)+2.5*alog10(1+z)
print,'UV Ab mag',Uvabmag,UVmagerr
stop
end

