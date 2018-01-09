pro uvlum_c28

clight = 3e10

;CSWA28
;The measured magnitudes of compact galaxies from cswa19_B.cat are about 0.2-0.3 mag dimmer than g band SDSS magnitudes
z=2.09
lambda = 1500.e-8
catfile = 'cswa28_B.cat'
imfile = 'cswa28_450_small.fits'
catfits_file = 'check_c28.fits'
magfile = 'cswa28_magnification.fits'
backgroundfile = 'check_c28_background.fits'

;read files
cat = rsex(catfile)
catfits = readfits(catfits_file)
im = readfits(imfile)
mag = readfits(magfile)
background = readfits(backgroundfile)

;There are 3 images in this arc
;One is above y = 597
;Another is between y=524 and 597
;Another is below 524

;array of indices
sizearr = size(im,/dimensions)
sizex = sizearr(0)
sizey = sizearr(1)
xarr = rebin(findgen(sizex),sizex,sizey)
yarr = rebin(transpose(findgen(sizey)),sizex,sizey)

;calculate total flux from sextractor
ids = [177,119,116,117,120,121,115,118]
sex_roi = []
sex_flux_arr = fltarr(n_elements(ids))+99
sex_fluxerr_arr = fltarr(n_elements(ids))+99
sex_mag_arr = fltarr(n_elements(ids))+99

for i=0,n_elements(ids)-1 do begin
   current_sex_roi = where(catfits eq ids(i) and yarr gt 597,croi)
;   current_sex_roi = where(catfits eq ids(i) and yarr gt 524 and yarr lt 597,croi)
;   current_sex_roi = where(catfits eq ids(i) and yarr lt 524,croi)
   if croi eq 0 then goto, skip
   sex_roi = [sex_roi,current_sex_roi]
   pos = where(cat.number eq ids(i),cpos)
   if cpos ne 1 then stop
   sex_Flux_arr(i) = cat[pos].flux_iso
   sex_fluxerr_arr(i) = cat[pos].fluxerr_iso
   sex_mag_arr(i) = cat[pos].mag_auto
   skip:continue
endfor

badelement = where(sex_flux_arr eq 99,ctbad)
if ctbad ne 0 then remove, badelement,sex_flux_arr,sex_fluxerr_arr,sex_mag_arr,ids

background_tot = total(background(sex_roi))
sex_flux_tot = total(im(sex_roi))-background_tot ;This one should be the same as flux_ISO
sex_fluxerr_tot = sqrt(total(sex_fluxerr_Arr^2))

print, 'Total Flux from sextractor',sex_flux_tot

sex_mag_tot = -2.5*alog10(sex_flux_tot/sex_flux_arr)+sex_mag_arr
sex_magerr_tot = 2.5*0.4343*sex_fluxerr_tot/sex_Flux_tot
print,'Lensed UV magnitude(sex)', sex_mag_tot
print, 'mean lensed UV magnitude(sex)', mean(sex_mag_tot),'+/-',sex_magerr_Tot

;demagnification
sex_delensed_flux_tot = total((im(sex_roi)-background(sex_roi))/mag(sex_roi))
sex_delensed_fluxerr_tot = sex_fluxerr_tot/total(mag(sex_roi))

print,'Delensed flux = ',sex_delensed_flux_tot
print,'Delensed fluxerr  = ',sex_delensed_fluxerr_tot
UVmag = -2.5*alog10(sex_delensed_flux_tot/sex_flux_arr)+sex_mag_arr
UVmagerr = 2.5*0.4343*sex_delensed_fluxerr_tot/sex_delensed_Flux_tot
print,'UV mag', UVmag,UVmagerr
UVmag = median(UVmag)

;Calculate absolute magnitude
dist = lumdist(z)*1.e6 ;parsecs
UVabmag = UVmag-5.*(alog10(dist)-1)+2.5*alog10(1+z)
print,'UV Ab mag',Uvabmag,UVmagerr
stop
end

