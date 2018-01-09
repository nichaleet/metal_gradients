pro cswa31ana
window, 0,xsize=1000,ysize=500
!p.multi = [0,2,1]
z=1.487
setplot, 39
;1) do H alpha image
!p.thick=1
;image = readfits('/scr2/nichal/workspace/reduced_data/mosaic/cswa31_Ha_Hn3_handmosaic_scaledsky_3hr.fits',header)
image = readfits('/scr2/nichal/workspace/reduced_data/mosaic/cswa31_Ha_tlc_Hn3_handmosaic_scaledsky_3hr.fits',header)

image(where(abs(image) gt .5)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[190:198,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,1.], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.015,0.02,0.03,0.04,0.05]
writefits,'hasum.fits',ha_sum

;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
cube = readfits('/scr2/nichal/workspace/output/cswa31_Ha_tlc_Hn3_handmosaic_scaledsky_3hr_acube.fits')
ha=cube[*,*,0]
goodha=where(finite(ha))
for i =0, n_elements(image[*,0,0])-1 do begin
   im_i = image[i,*,*]
   ha_spectrum[i] = total(im_i(goodha))
endfor
;for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,10:28,35:50])+total(image[i,29:55,30:42])
wavelength=getwl_filter('Hn3')

;fit Halpha spectrum
good = where(wavelength gt 1.625 and wavelength lt 1.643)
wl = wavelength(good)
ha = ha_spectrum(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[-0.2,3.],position=[.55,0.15,.9,.9],title='Halpha'

;somehow the fit for wl and ha just gave the peak at the wrong peak. so we have to modify some input here;
forhafit = where(wavelength gt 1.630 and wavelength lt 1.635)
wlfine = wavelength(forhafit)
hafine = ha_spectrum(forhafit)

fit = gaussfit(wlfine,hafine,param,nterms=4,sigma=sigma)

oplot, wlfine,fit,color=50,thick =2
redshift =  param(1)/.656461-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print,'parameter a,b,c,d:', param
print, redshift,veldisp,intrinsic_veldisp
areaHa = abs(param(0)*param(2))*sqrt(2.*!pi)
error_areaHa = areaHa*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)*sqrt(2.*!pi)

restore,'/scr2/nichal/workspace/Calibrate_Ha/flux_calibration.sav'
pos1 = where(flux_calibration.name eq 'cswa31')
lamb_res = sxpar(header,'CDELT1')/1000.
conversionfactor = flux_calibration(pos1).calibfactor
conversionfactor_err = flux_calibration(pos1).calibfactor_err
totHaflux = areaHa*conversionfactor/lamb_res ;erg/s/cm^2
tothaflux_err = totHaflux*sqrt(error_areaHa^2/areaHa^2+conversionfactor_err^2/conversionfactor^2) ;erg/s/cm^2
print,'total Ha flux ',totHaflux,'erg/s/cm^2',tothaflux_err
stop


;fit NII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,8:20,30:40])
good = where(wavelength gt 1.634 and wavelength lt 1.643)
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
setplot, 14
;1) do OIII image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa31_Hb_Jbb_pipelinemosaic_scaledsky_030hr.fits",header)
image(where(abs(image) gt .08)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[467:473,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.015,0.02,0.03,0.04,0.05]
writefits,'oiiisum.fits',ha_sum
;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,5:40,5:15])
wavelength = getwl_filter('Jbb')
;plot,wavelength,ha_spectrum,psym=10,/noerase,yrange=[-0.5,1.0],position=
;vline,(z+1)*[656.4,658.5]

;fit Halpha spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,40:60,5:15])
good = where(wavelength gt 1.243 and wavelength lt 1.250)
wl = wavelength(good)
ha = ha_spectrum(good)
weightHa = weight(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[-0.2,2.5],position=[.1,0.5,.45,.9],title='[OIII]'
error = 1./sqrt(weightHa)
fit = gaussfit(wl,ha,param,nterms=4,measure_errors=error,sigma=sigma)
oplot, wl,fit,color=50,thick =2
oplot, wl,weightHa/max(weightHa),color=200,psym=10
redshift =  param(1)/.500824-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print, param
print, redshift,veldisp,intrinsic_veldisp
areaOIII = abs(param(0)*param(2))*sqrt(2*!pi)
error_areaOIII = areaOIII*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)*sqrt(2.*!pi)

;fit Hb spectrum

good = where(wavelength gt 1.205 and wavelength lt 1.213)
wl = wavelength(Good)
NII = ha_spectrum(good)
weightNII = weight(good)
fit2ndline,param(2),0.486269*(redshift+1.),wl,NII,weightNII,fit,area,sigma_are

plot,wl,NII,psym=10,xtitle='nm',ytitle = 'flux',/noerase,position=[.6,0.5,.95,.9],title='Hb'
oplot,wl,fit,color=20,thick=2
areaHb = area
error_areaHb = sigma_area
print, 'log(OIII/Hb) = ',alog10(areaOIII/areaHb)
print, 'error = ',0.4343*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)

print, '(OIII/Hb) = ',(areaOIII/areaHb)
print, 'error = ',(areaOIII/areaHb)*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)
end
