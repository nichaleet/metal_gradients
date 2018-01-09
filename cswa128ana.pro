pro cswa128ana
window, 0,xsize=1000,ysize=800
!p.multi = [0,2,2]
z=2.226
setplot, 14
;1) do H alpha image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa128_Ha_Kn2_handmosaic_sky_130hr.fits",header)

image(where(abs(image) gt 1.)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[320:340,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,1.], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.05,0.07,0.1,0.2,0.3]
writefits,'hasum.fits',ha_sum

;To plot wavelengths vs intensity map. Take only galaxy region. 
cubesize = size(image,/dimensions)
ha_spectrum = fltarr(cubesize(1),cubesize(2))
magnification2d = readfits('/scr2/nichal/workspace/massmodel/mauger_model/cswa128_magnification_idl.fits')
magnification_cube = fltarr(cubesize(0),cubesize(1),cubesize(2))
for i=0,cubesize(0)-1 do magnification_cube[i,*,*]=magnification2d


cube = readfits('/scr2/nichal/workspace/output/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube.fits')
ha=cube[*,*,0]
goodha=where(finite(ha))
for i =0, n_elements(image[*,0,0])-1 do begin
   im_i = image[i,*,*]
   ha_spectrum[i] = total(im_i(goodha))
endfor
;for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,13:34,29:39])+total(image[i,40:62,27:36])
wavelength = 2.036+2.5e-4*findgen(421)


;fit Halpha spectrum
good = where(wavelength lt 2.127 and wavelength gt 2.113)
wl = wavelength(good)
ha = ha_spectrum(good)
secondgal_model = [1.40841,2.45099,3.44813,4.25331,4.43877,4.39465,3.44813,2.45099,1.40841]-0.48
secondgal_wl = [wl[17:25]]
ha_firstgal = ha
ha_firstgal[17:25] = ha[17:25]-secondgal_model
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[0.,6.],position=[.55,0.5,.9,.9],title='Halpha'
oplot,wl,ha_firstgal,linestyle=2
fit = gaussfit(wl,ha_firstgal,param,nterms=4,sigma=sigma)
oplot, wl,fit,color=50,thick =2
redshift =  param(1)/.656461-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print,'parameter a,b,c,d:', param
print, redshift,veldisp,intrinsic_veldisp
areaHa = abs(param(0)*param(2))*sqrt(2.*!pi)
error_areaHa = areaHa*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)*sqrt(2.*!pi)

;find total Ha flux
ha_secondgal = ha-fit
fit2 = gaussfit(wl,ha_secondgal,param2,nterms=4,sigma=sigma2)
oplot,wl,fit+fit2,color=200,thick=2
area2 = abs(param2(0)*param2(2))*sqrt(2.*!pi)
area2_err = area2*sqrt((sigma2(0)/param2(0))^2+(sigma2(2)/param2(2))^2)*sqrt(2.*!pi)
restore,'/scr2/nichal/workspace/Calibrate_Ha/flux_calibration.sav'
pos = where(flux_calibration.name eq 'cswa128')
conversionfactor = flux_calibration(pos).calibfactor
conversionfactor_err = flux_calibration(pos).calibfactor_err
totHaflux = (areaHa+area2)*conversionfactor/2.5e-4 ;erg/s/cm^2
tothaflux_err = sqrt(error_areaHa^2+area2_err^2)*conversionfactor/2.5e-4 ;erg/s/cm^2
print,'total Ha flux ',totHaflux,'erg/s/cm^2',tothaflux_err
stop

;fit NII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,6:10,32:42])
good = where(wavelength gt 2.121 and wavelength lt 2.127)
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

stop
;OIII image

window, 1,xsize=1000,ysize=800
!p.multi = [0,2,2]
z=2.226
setplot, 14
;1) do OIII image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa128_Hb_Hbb_mosaic_scaledsky_130hr.fits",header)
image(where(abs(image) gt .08)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[707:719,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.03,0.04,0.05,0.1,0.2]
writefits,'oiiisum.fits',ha_sum
;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,12:30,12:19])+total(image[i,40:56,10:16])
wavelength = 1.473+0.0002*findgen(1651)


;fit OIII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,6:10,10:25])
good = where(wavelength gt 1.610 and wavelength lt 1.620)
wl = wavelength(good)
ha = ha_spectrum(good)
weightHa = weight(good)
plot, wl,ha,psym=10,xtitle='nm',ytitle = 'flux',/noerase,yrange=[-0.2,2.5],position=[.1,0.5,.45,.9],title='[OIII]'
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

good = where(wavelength gt 1.564 and wavelength lt 1.572)
wl = wavelength(Good)
NII = ha_spectrum(good)
weightNII = weight(good)
fit2ndline,param(2),0.4861363*(redshift+1.),wl,NII,weightNII,fit,area,sigma_area

plot,wl,NII,psym=10,xtitle='nm',ytitle = 'flux',/noerase,yrange=[-0.2,2.5],position=[.6,0.5,.95,.9],title='Hb'
oplot,wl,fit,color=20,thick=2
areaHb = area
error_areaHb = sigma_area
print, 'log(OIII/Hb) = ',alog10(areaOIII/areaHb)
print, 'error = ',0.4343*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)

print, '(OIII/Hb) = ',(areaOIII/areaHb)
print, 'error = ',(areaOIII/areaHb)*sqrt((error_areaHb/areaHb)^2+(error_areaOIII/areaOIII)^2)
end
