pro cswa20ana

window, 0,xsize=1000,ysize=800
!p.multi = [0,2,2]
z=1.4330
setplot, 14
;1) do H alpha image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa20_Ha_Hn2_handmosaic_scaledsky_1hr.fits",header)
image(where(abs(image) gt .1)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[322:333,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.05,0.07,0.1,0.2,0.3]
writefits,'hasum.fits',ha_sum

;To plot wavelengths vs intensity map. Take only galaxy region. 


ha_spectrum = fltarr(n_elements(image[*,0,0]))
cube = readfits('/scr2/nichal/workspace/output/cswa20_Ha_tlc_Hn2_handmosaic_scaledsky_1hr_acube.fits')
ha=cube[*,*,0]
goodha=where(finite(ha))
for i =0, n_elements(image[*,0,0])-1 do begin
   im_i = image[i,*,*]
   ha_spectrum[i] = total(im_i(goodha))
endfor
;for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,19:29,23:25])+total(image[i,39:53,30:40])
wavelength = 1.532+0.0002*findgen(391)


;fit Halpha spectrum
good = where(wavelength gt 1.5945 and wavelength lt 1.606)
wl = wavelength(good)
ha = ha_spectrum(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[-0.2,2.],position=[.1,0.55,.45,.95],title='Halpha'
fit = gaussfit(wl,ha,param,nterms=4,sigma=sigma)
oplot, wl,fit,color=50,thick =2
redshift =  param(1)/.656461-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print,'parameter a,b,c,d:', param
print, redshift,veldisp,intrinsic_veldisp
areaHa = abs(param(0)*param(2))*sqrt(2.*!pi)
error_areaHa = areaHa*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)*sqrt(2.*!pi)

restore,'/scr2/nichal/workspace/Calibrate_Ha/flux_calibration.sav'
pos1 = where(flux_calibration.name eq 'cswa20')
lamb_res = sxpar(header,'CDELT1')/1000.
conversionfactor = flux_calibration(pos1).calibfactor
conversionfactor_err = flux_calibration(pos1).calibfactor_err
totHaflux = areaHa*conversionfactor/lamb_res ;erg/s/cm^2
tothaflux_err = totHaflux*sqrt(error_areaHa^2/areaHa^2+conversionfactor_err^2/conversionfactor^2) ;erg/s/cm^2
print,'total Ha flux ',totHaflux,'erg/s/cm^2',tothaflux_err
stop

;fit NII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,34:53,21:29])
good = where(wavelength gt 1.599 and wavelength lt 1.606)
wl = wavelength(Good)
NII = ha_spectrum(good)
weight = weight(good)
print, .658523*(redshift+1.)
fit2ndline,param(2),.658523*(redshift+1.),wl,NII,weight,fit,area,sigma_area
oplot,wl,fit,color=20,thick=2
plot, wl,NII,psym=10,xtitle='micron',ytitle = 'flux',/noerase,position=[.6,0.55,.95,.95],title='NII'
oplot,wl,fit,color=20,thick=2
areaNII = area
error_areaNII = sigma_area

print, 'log(NII/Ha) = ',alog10(areaNII/areaHa)
print, 'error = ',0.4343*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)

print, '(NII/Ha) = ',(areaNII/areaHa)
print, 'error = ',(areaNII/areaHa)*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)

;OIII image

window, 1,xsize=1000,ysize=800
!p.multi = [0,2,2]
z=1.433
setplot, 14
;1) do OIII image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa20_Hb_Jn1_handmosaic_sky_075hr.fits",header)
image(where(abs(image) gt .1)) = 0.

for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[292:304,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.03,0.04,0.05,0.1,0.2]
writefits,'oiiisum.fits',ha_sum
;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,40:55,28:38])+total(image[i,22:28,20:25])
wavelength = 1.174+0.00014987*findgen(388)


;fit OIII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,22:38,29:36])
good = where(wavelength gt 1.215 and wavelength lt 1.222)
wl = wavelength(good)
ha = ha_spectrum(good)
weightHa = weight(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[-0.2,2.5],position=[.1,0.5,.45,.9],title='[OIII]'
weightha(where(wl ge 1.565 and wl le 1.5655)) = 0.
error = 1./sqrt(weightHa)
fit = gaussfit(wl,ha,param,nterms=4,measure_errors=error,sigma=sigma)
oplot, wl,fit,color=50,thick =2
;oplot, wl,weightHa/max(weightHa),color=200,psym=10
redshift =  param(1)/.5006843-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print, param
print, redshift,veldisp,intrinsic_veldisp
areaOIII = abs(param(0)*param(2))*sqrt(2*!pi)
error_areaOIII = areaOIII*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)*sqrt(2.*!pi)

;fit Hb spectrum

good = where(wavelength gt 1.178 and wavelength lt 1.186)
wl = wavelength(Good)
NII = ha_spectrum(good)
weightNII = weight(good)
fit2ndline,param(2),0.4861363*(redshift+1.),wl,NII,weightNII,fit

plot,wl,NII,psym=10,xtitle='nm',ytitle = 'flux',/noerase,yrange=[-0.2,2.5],position=[.6,0.5,.95,.9],title='Hb'
oplot,wl,fit,color=20,thick=2
areaHb = area
error_areaHb = sigma_area
print, 'log(OIII/Hb) = ',alog10(areaOIII/areaHb)
print, 'error = ',0.4343*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)

print, '(OIII/Hb) = ',(areaOIII/areaHb)
print, 'error = ',(areaOIII/areaHb)*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)

end
