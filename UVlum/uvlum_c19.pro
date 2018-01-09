pro uvlum_c19

clight = 3e10

;CSWA19
;The measured magnitudes of compact galaxies from cswa19_B.cat are about 0.2-0.3 mag dimmer than g band SDSS magnitudes
z=2.03
lambda = 4773.7e-8
catfile = 'cswa19_B.cat'
imfile = 'cswa19_hst_wfc3_f475w_drz.fits'
roifile = 'cswa19_roi.fits'
catfits_file = 'check_c19.fits'
magfile = 'cswa19_magnification.fits'
backgroundfile = 'check_c19_background.fits'

;read files
cat = rsex(catfile)
catfits = mrdfits(catfits_file,1)
im = mrdfits(imfile,1)
im_roi = readfits(roifile)
mag = readfits(magfile)
background = mrdfits(backgroundfile,1)

;calculate total flux from sextractor
roi = where(im_roi eq 1.)
ids = (catfits(roi))[uniq((catfits(roi)),sort(catfits(roi)))]
ids = ids(where(ids ne 0))

sex_roi = []
sex_flux_arr = fltarr(n_elements(ids))
sex_fluxerr_arr = fltarr(n_elements(ids))
sex_mag_arr = fltarr(n_elements(ids))

for i=0,n_elements(ids)-1 do begin
   current_sex_roi = where(catfits eq ids(i),croi)
   if croi eq 0 then stop
   sex_roi = [sex_roi,current_sex_roi]

   pos = where(cat.number eq ids(i),cpos)
   if cpos ne 1 then stop
   sex_Flux_arr(i) = cat[pos].flux_iso
   sex_fluxerr_arr(i) = cat[pos].fluxerr_iso
   sex_mag_arr(i) = cat[pos].mag_auto 
endfor

background_tot = total(background(sex_roi))
sex_flux_tot = total(im(sex_roi))-background_tot ;This one should be similar to flux_ISO
hand_flux_tot = total(im(roi))-total(background(roi))
flux_diff = abs(hand_flux_tot-sex_flux_tot)
sex_fluxerr_tot = sqrt(total(sex_fluxerr_Arr^2))
if flux_diff gt sex_fluxerr_tot then sex_fluxerr_tot = flux_diff

print, 'Total Flux from sextractor',sex_flux_tot,total(sex_flux_arr) ;These two should be the same!
print, 'Total Flux from hand',hand_flux_tot

sex_mag_tot = -2.5*alog10(sex_flux_tot/sex_flux_arr)+sex_mag_arr
sex_magerr_tot = 2.5*0.4343*sex_fluxerr_tot/sex_Flux_tot
print,'Lensed UV magnitude(sex)', sex_mag_tot
print, 'mean lensed UV magnitude(sex)', mean(sex_mag_tot),'+/-',sex_magerr_Tot

;demagnification
hand_delensed_flux_tot = total((im(roi)-background(roi))/mag(roi))
sex_delensed_flux_tot = total((im(sex_roi)-background(sex_roi))/mag(sex_roi))
sex_delensed_fluxerr_tot = sex_fluxerr_tot/total(mag(sex_roi))
if abs(hand_delensed_flux_tot-sex_delensed_flux_tot) gt sex_delensed_fluxerr_tot then sex_delensed_fluxerr_tot = abs(hand_delensed_flux_tot-sex_delensed_flux_tot)
print,'Delensed flux (hand,sex) = ',hand_Delensed_flux_tot,sex_delensed_flux_tot
print,'Delensed fluxerr (sex) = ',sex_delensed_fluxerr_tot
UVmag = -2.5*alog10(sex_delensed_flux_tot/sex_flux_arr)+sex_mag_arr
UVmagerr = 2.5*0.4343*sex_delensed_fluxerr_tot/sex_delensed_Flux_tot
print,'UV mag', UVmag,UVmagerr
UVmag = median(UVmag)

;Calculate absolute magnitude
dist = lumdist(z)*1.e6 ;parsecs
UVabmag = UVmag-5.*(alog10(dist)-1)
print,'UV Ab mag',Uvabmag,UVmagerr
stop
end
