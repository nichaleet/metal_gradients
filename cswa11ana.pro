pro cswa11ana

window, 0,xsize=1000,ysize=800
!p.multi = [0,2,2]
z=1.41
setplot, 12
;1) do H alpha image
!p.thick=1
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa11_Ha_Hn2_handmosaic_scaledsky_430hr.fits",header)
image(where(abs(image) gt .1)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[248:254,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.01,0.02,0.03,0.2,0.3]
writefits,'hasum.fits',ha_sum


;;For the left galaxy (shorter lambda)
;
;;To plot wavelengths vs intensity map. Take only galaxy region. 
;ha_spectrum = fltarr(n_elements(image[*,0,0]))
;for i =0, n_elements(image[*,0,0])-1 do begin
;   ha_spectrum[i] = total(image[i,7:63,44:52])
;endfor
;wavelength = (sxpar(header,'crval1')+sxpar(header,'cdelt1')*findgen(sxpar(header,'naxis1')))/1000. ;in micron
;
;
;;fit Halpha spectrum
;good = where(wavelength gt 1.578 and wavelength lt 1.589)
;wl = wavelength(good)
;ha = ha_spectrum(good)
;plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[0.,5.],position=[.1,0.55,.45,.95],title='Halpha'
;fit = gaussfit(wl,ha,param,nterms=4,sigma=sigma) ; flux = aExp[(x-b)^2/2c^2]+d
;oplot, wl,fit,color=50,thick =2
;redshift =  param(1)/.65628-1.
;veldisp = param(2)/param(1)*3.e5
;intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
;print,'parameter a,b,c,d:', param
;print, sigma
;print,'z, disp, disp_intrinsic: ', redshift,veldisp,intrinsic_veldisp
;areaHa = abs(param(0)*param(2))*sqrt(2.*!pi)
;error_areaHa = areaHa*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)
;
;;fit NII spectrum
;weight = wavelength*0.
;for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,50:60,40:50])
;good = where(wavelength gt 1.5835 and wavelength lt 1.589)
;wl = wavelength(Good)
;NII = ha_spectrum(good)
;weight = weight(good)
;print, .658345*(redshift+1.)
;
;fit2ndline,param(2),.658345*(redshift+1.),wl,NII,weight,fit,area,sigma_area
;oplot,wl,fit,color=20,thick=2
;
;plot, wl,NII,psym=10,xtitle='micron',ytitle = 'flux',/noerase,position=[.6,0.55,.95,.95],title='NII'
;oplot,wl,fit,color=20,thick=2
;areaNII = area
;error_areaNII = sigma_area
;
;print, 'log(NII/Ha) = ',alog10(areaNII/areaHa)
;print, 'error = ',0.4343*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)
;
;print, '(NII/Ha) = ',(areaNII/areaHa)
;print, 'error = ',(areaNII/areaHa)*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)


;For the right galaxy

window, 1,xsize=1000,ysize=800
!p.multi = [0,2,2]
z=1.41
setplot, 12


;To plot wavelengths vs intensity map. Take only galaxy region. 
cube = readfits('/scr2/nichal/workspace/output/cswa11_Ha_tlc_Hn2_handmosaic_scaledsky_430hr_secondgal_acube.fits')
ha = cube[*,*,0]
goodha = where(finite(ha))
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do begin
   im_i = image[i,*,*]
   ha_spectrum[i] = total(im_i(goodha))
endfor
wavelength = (sxpar(header,'crval1')+sxpar(header,'cdelt1')*findgen(sxpar(header,'naxis1')))/1000. ;in micron


;fit Halpha spectrum
good = where(wavelength gt 1.576 and wavelength lt 1.595)
wl = wavelength(good)
ha = ha_spectrum(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[0.,1.5],position=[.1,0.55,.45,.95],title='Halpha'
error = sqrt(abs(ha))
bad = where(wl gt 1.58 and wl le 1.5812)
skymodel = [0.400675,0.856361,1.24982,1.45918,1.45918,1.24982,0.856361,0.400675]-0.26
skywl    = wl[20:27]
ha_sub_sky = ha
for i=0,n_elements(skywl)-1 do begin
   skyreg = where(wl eq skywl(i))
   if skyreg eq -1 then stop
   ha_sub_sky(skyreg) =  ha_sub_sky(skyreg)-skymodel(i)
endfor
oplot,wl,ha_sub_sky,linestyle=2   
fit = gaussfit(wl,ha_sub_sky,param,nterms=4,sigma=sigma) ; flux = aExp[(x-b)^2/2c^2]+d
oplot, wl,fit,color=50,thick =2
redshift =  param(1)/.656461-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print,'parameter a,b,c,d:', param
print, sigma
print,'z,vel disp, intrinsic vel disp: ', redshift,veldisp,intrinsic_veldisp
areaHa = abs(param(0)*param(2))*sqrt(2.*!pi)
error_areaHa = areaHa*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)

restore,'/scr2/nichal/workspace/Calibrate_Ha/flux_calibration.sav'
pos1 = where(flux_calibration.name eq 'cswa11feb')
pos2 = where(flux_calibration.name eq 'cswa11mar')
conversionfactor = mean([flux_calibration(pos1).calibfactor,flux_calibration(pos2).calibfactor])
conversionfactor_err = 0.5*sqrt((flux_calibration(pos1).calibfactor_err)^2+(flux_calibration(pos2).calibfactor_err)^2)
totHaflux = areaHa*conversionfactor/2.e-4 ;erg/s/cm^2
tothaflux_err = totHaflux*sqrt(error_areaHa^2/areaHa^2+conversionfactor_err^2/conversionfactor^2) ;erg/s/cm^2
print,'total Ha flux ',totHaflux,'erg/s/cm^2',tothaflux_err
stop

;fit NII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,50:60,40:50])
good = where(wavelength gt 1.5855 and wavelength lt 1.589)
wl = wavelength(good)
NII = ha_spectrum(good)
weight = weight(good)
print,'NII wavelength: ', .658523*(redshift+1.)
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

end
