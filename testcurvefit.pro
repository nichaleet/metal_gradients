
pro singlet, x, par, f, pder

  x0 = par[0]
  aa = par[1]
  ww = par[2]
  cont = par[3]

  f = (aa/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-x0)^2/ww^2) + cont

  pder = fltarr(n_elements(x),n_elements(par)) ; no value returned.

end

pro testcurvefit
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa11_Ha_Hn2_mosaic_scaledsky_430hr.fits",header)
image(where(abs(image) gt .1)) = 0.
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[248:254,*,*],1)
imdisp, ha_sum,position=[.05,.05,.5,.5], out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.01,0.02,0.03,0.2,0.3]
writefits,'hasum.fits',ha_sum


;For the left galaxy

;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,18:34,37:47])
wavelength = (sxpar(header,'crval1')+sxpar(header,'cdelt1')*findgen(sxpar(header,'naxis1')))/1000. ;in micron


;fit Halpha spectrum
good = where(wavelength gt 1.578 and wavelength lt 1.585)
wl = wavelength(good)
ha = ha_spectrum(good)

plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[0.,2.],position=[.1,0.55,.45,.95],title='Halpha'
a=[1.5817,0.0016,0.0008,0.]
fita=[1,1,1,1]
fit = curvefit(wl,ha,weight,a,sigmafit,chisq=1,fita=fit_a,function_name='singlet',/noderivative)
oplot, wl,fit,color=50,thick =2

stop
end
