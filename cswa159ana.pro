pro cswa159ana

window, 0,xsize=1000,ysize=800
!p.multi = [0,2,2]
z=2.30
setplot, 12
;1) do H alpha image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa159_Ha_tlc_Kc3_pipelinemosaic_sky_2hr.fits",header)

image(where(abs(image) gt 1.)) = 0.
;scale = abs(stddev(image[where(image ne 0. and finite(image) eq 1)]))
;image = image/scale
image(204,44,37) = 0.

;;add noise
;noise = randomn(10,size(image,/dimension))*0.5421*0.
;plot,image[*,40,36],psym=10
;image = image+noise
;oplot,image[*,40,36],color=50,psym=10
;pause=get_kbrd()
;erase

for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[173:183,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.01,0.02,0.03,0.2,0.3]
writefits,'hasum.fits',ha_sum

;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
cube = readfits('/scr2/nichal/workspace/output/cswa159_Ha_tlc_Kc3_pipelinemosaic_sky_2hr_acube.fits')
ha=cube[*,*,0]
goodha=where(finite(ha))
for i =0, n_elements(image[*,0,0])-1 do begin
   im_i = image[i,*,*]
   ha_spectrum[i] = total(im_i(goodha))
endfor
;for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,25:60,33:40])
wavelength = (sxpar(header,'crval1')+sxpar(header,'cdelt1')*findgen(sxpar(header,'naxis1')))/1000. ;in micron

;;calculate noise
;cleanregion = where(wavelength gt 2.155 and wavelength lt 2.162)
;noise = stdev(ha_spectrum(cleanregion))
;print, 'noise = ',noise
;noise_per_pixel = noise/sqrt(float(n_elements(goodha)))
;print, 'noise per pixel = ',noise_per_pixel*scale, noise_per_pixel,'(no scale)'
;pause=get_kbrd() 

;fit Halpha spectrum
good = where(wavelength gt 2.155 and wavelength lt 2.1745)
wl = wavelength(good)
ha = ha_spectrum(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,position=[.1,0.55,.45,.95],title='Halpha'
fit = gaussfit(wl,ha,param,nterms=4,sigma=sigma) ; flux = aExp[(x-b)^2/2c^2]+d
oplot, wl,fit,color=50,thick =2
redshift =  param(1)/.656461-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print,'parameter a,b,c,d:', param
print, sigma
print,'z, disp, disp_intrinsic: ', redshift,veldisp,intrinsic_veldisp
areaHa = abs(param(0)*param(2))*sqrt(2.*!pi)
error_areaHa = areaHa*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)

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
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,6:20,32:45])
good = where(wavelength gt 2.1704 and wavelength lt 2.1735)
good = where(wavelength gt 2.168 and wavelength lt 2.175)

wl = wavelength(Good)
NII = ha_spectrum(good)
weight = weight(good)
print, .658523*(redshift+1.)
width=param(2)
fit2ndline,width,.658523*(redshift+1.),wl,NII,weight,fit,area,sigma_area
oplot,wl,fit,color=20,thick=2

plot, wl,NII,psym=10,xtitle='micron',ytitle = 'flux',/noerase,position=[.6,0.55,.95,.95],title='NII'
oplot,wl,fit,color=20,thick=2
plot, wl,weight,psym=10,xtitle='micron',ytitle = 'weight',/noerase,position=[.6,0.05,.95,.45]
areaNII = area
error_areaNII = sigma_area

print, 'log(NII/Ha) = ',alog10(areaNII/areaHa)
print, 'error = ',0.4343*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)

print, '(NII/Ha) = ',(areaNII/areaHa)
print, 'error = ',(areaNII/areaHa)*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)

totalnoise = sqrt(width/(wavelength(1)-wavelength(0))*5.)*noise
print, 'SN of NII = ',areaNII/totalnoise/(wavelength(1)-wavelength(0))
stop


end
