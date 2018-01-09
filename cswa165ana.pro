pro cswa165ana
window, 0,xsize=1000,ysize=500
!p.multi = [0,2,1]
z=2.1282
setplot, 14
;1) do H alpha image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa165_Ha_Kn2_handmosaic_sky_1hr.fits",header)
image(where(abs(image) gt .5)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[64:74,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,1.], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.015,0.02,0.03,0.04,0.05]
writefits,'hasum.fits',ha_sum

;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
cube = readfits('/scr2/nichal/workspace/output/cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr_acube.fits')
ha=cube[*,*,0]
goodha=where(finite(ha))
for i =0, n_elements(image[*,0,0])-1 do begin
   im_i = image[i,*,*]
   ha_spectrum[i] = total(im_i(goodha))
endfor
;for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,25:36,33:40])+total(image[i,46:56,33:40])
wavelength = 2.036+0.00025*findgen(421)


;fit Halpha spectrum
good = where(wavelength lt 2.056 and wavelength gt 2.0505)
wl = wavelength(good)
ha = ha_spectrum(good)
plot, wavelength,ha_spectrum,psym=10,xtitle='micron',ytitle = 'flux',/noerase,xrange=[2.045,2.070],yrange=[-0.2,3.],position=[.55,0.15,.9,.9],title='Halpha'
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
pos1 = where(flux_calibration.name eq 'cswa159')
lamb_res = sxpar(header,'CDELT1')/1000.
conversionfactor = flux_calibration(pos1).calibfactor
conversionfactor_err = flux_calibration(pos1).calibfactor_err
totHaflux = areaHa*conversionfactor/lamb_res ;erg/s/cm^2
tothaflux_err = totHaflux*sqrt(error_areaHa^2/areaHa^2+conversionfactor_err^2/conversionfactor^2) ;erg/s/cm^2
print,'total Ha flux ',totHaflux,'erg/s/cm^2',tothaflux_err
stop ;get 14-18 e-17 erg/s/cm2
;fit NII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,8:20,30:40])
good = where(wavelength gt 2.055 and wavelength lt 2.064)
wl = wavelength(Good)
NII = ha_spectrum(good)
weight = weight(good)
fit2ndline,param(2),.658523*(redshift+1.),wl,NII,weight,fit,area,sigma_area
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
z=2.1282
setplot, 14
;1) do OIII image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa165_Hb_Hbb_mosaic_scaledsky_1hr.fits",header)
image(where(abs(image) gt .08)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[467:473,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.015,0.02,0.03,0.04,0.05]
writefits,'oiiisum.fits',ha_sum
;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,43:55,12:22])+total(image[i,20:38,14:19])
wavelength = 1.473+0.0002*findgen(1651)


;fit Halpha spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,5:13,13:19])
good = where(wavelength gt 1.563 and wavelength lt 1.570)
wl = wavelength(good)
ha = ha_spectrum(good)
weightHa = weight(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[-0.2,2.5],position=[.1,0.5,.45,.9],title='[OIII]'
weightha(where(wl ge 1.565 and wl le 1.5655)) = 0.
error = 1./sqrt(weightHa)
fit = gaussfit(wl,ha,param,nterms=4,measure_errors=error,sigma=sigma)
oplot, wl,fit,color=50,thick =2
oplot, wl,weightHa/max(weightHa),color=200,psym=10
redshift =  param(1)/.5006843-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print, param
print, redshift,veldisp,intrinsic_veldisp
areaOIII = abs(param(0)*param(2))*sqrt(2*!pi)
error_areaOIII = areaOIII*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)*sqrt(2.*!pi)

;fit Hb spectrum

good = where(wavelength gt 1.517 and wavelength lt 1.524)
wl = wavelength(Good)
NII = ha_spectrum(good)
weightNII = weight(good)
fit2ndline,param(2),0.4861363*(redshift+1.),wl,NII,weightNII,fit,area,sigma_are

plot,wl,NII,psym=10,xtitle='nm',ytitle = 'flux',/noerase,yrange=[-0.2,2.5],position=[.6,0.5,.95,.9],title='Hb'
oplot,wl,fit,color=20,thick=2
areaHb = area
error_areaHb = sigma_area
print, 'log(OIII/Hb) = ',alog10(areaOIII/areaHb)
print, 'error = ',0.4343*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)

print, '(OIII/Hb) = ',(areaOIII/areaHb)
print, 'error = ',(areaOIII/areaHb)*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)
end
