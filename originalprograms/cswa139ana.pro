pro doublet,x,par,f

zz   = par[0]
aa1 = par[1]
aa2 = par[2]
ww  = par[3]
cont = par[4]
f = ((aa1/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-656.4*(1+zz))^2/ww^2))+((aa2/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-658.5*(1+zz))^2/ww^2))+cont

end

pro doubletoiii,x,par,f

zz   = par[0]
aa1 = par[1]
aa2 = par[2]
ww  = par[3]
cont = par[4]
f = ((aa1/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-500.8*(1+zz))^2/ww^2))+((aa2/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-486.2*(1+zz))^2/ww^2))+cont

end

pro cswa139ana
z=2.54
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa139_Ha_Kc5_mosaic_sky_130hr.fits",header)
sz= size(image)
for i =0, sz[1]-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[120:154,*,*],1)

imdisp, ha_sum, out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.04,0.05,0.08,0.1,0.2]

;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, sz[1]-1 do ha_spectrum[i] = total(image[i,16:40,40:48])
wavelength = findgen(sz[1])*0.25+2292
plot,wavelength,ha_spectrum,psym=1

xmins = 16
xmaxs = 40
ymins = 40
ymaxs = 48
lineha = 656.4*(1.+z)
lineNII = 658.5*(1.+z)
range = 10
sigmalim = 5. 
wres = wavelength[1]-wavelength[0] ; in nm per pixel

ok = where(wavelength ge lineha-range and wavelength le lineNII+range)
finalcube = image[ok,*,*]
;rescale final cube 
; NOTE - IDL fit routine does not work with very small numbers
;        (i.e. cgs flux units ~1e-19), so make sure cube values are of
;        order unity before fitting spectra.
; Need to scale flux values to order unity
scale = abs(stddev(finalcube[where(finalcube)]))
finalcube = finalcube/scale

wl = wavelength[ok]
finalsz = size(finalcube)

weight = fltarr(finalsz[1])
for i = 0,finalsz[1]-1 do weight[i] = 1/variance(finalcube[i,xmins:xmaxs,ymins:ymaxs])

signal = 3.*stddev(finalcube[where(finalcube)])
offset = median(finalcube[where(finalcube)])

fitcube = 1./fltarr(finalsz[1],finalsz[2],finalsz[3]) ;holds best fit
acube   = 1./fltarr(finalsz[2],finalsz[3],5) ;holds best parameters
aerrorcube = 1./fltarr(finalsz[2],finalsz[3],4) ;holds sigma values
count = 0

for ii=xmins,xmaxs do for jj=ymins,ymaxs do if finalcube[0,ii,jj] then begin
   ; only look at pixels where zero-th wl bin is non zero
   x = wl
   y = finalcube[*,ii,jj]
   av = (moment(y))[0]
   chi0 = total(weight*(y-av)^2)

   a = [z,signal,signal,1.,offset]    ;[redshift,height1,height2,width,y offset]   
   fit = curvefit(x,y,weight,a,sigmafit,chisq = chi2,function_name = "doublet",fita = [1,1,1,1,1],/noderivative)
   
   chi2 = total(weight*(y-fit)^2) ;chisquare of line fit
   chi2_var = chi0+1.-chi2-n_elements(a) ;variance of the line fit
   nsigma = sqrt(chi2_var)  ; significance of the line fit in standard deviations

   if nsigma gt sigmalim then begin
      plot, x, y, yrange=[min(y),max(y)],psym=1
      oplot, x, fit, color=255
      print, stddev(y), median(1/sqrt(weight))
      fitcube[*,ii,jj] = fit
      acube[ii,jj,*] = [a[0:3], nsigma]
      aerrorcube[ii,jj,*] = sigmafit[0:3]
      print, 'fit found with ', nsigma, ' sigma at ', ii, jj
      print, 'chi^2: no line, line, sigma(line) = ', chi0, chi2, nsigma
      count = count+1
      ;wait, 0.5
   endif 
endif

print, count,'pixels out of',(xmaxs-xmins+1)*(ymaxs-ymins+1),'pixels have sigma fit >', sigmalim

; Turn fit parameters into actual physical properties
acube[where(acube eq 0)] = 1/0.
aerrorcube[where(aerrorcube eq 0)] = 1/0
; Convert linewidth in nm to rest-frame km/s
acube[*,*,3] = abs(acube[*,*,3])/lineHA* 3d5   ; sigma in km/s
aerrorcube[*,*,3] = abs(aerrorcube[*,*,3]) / lineHA * 3d5
; Append the measured redshift, z, to the 'acube' data
z_meas = acube[*,*,0] 
acube = [[[acube]], [[z_meas]]]
; Convert redshift to rest-frame km/s
acube[*,*,0] = (((acube[*,*,0] +1 )/(z+1))-1) * 3d5
aerrorcube[*,*,0] = aerrorcube[*,*,0] * 3d5
; Fit amplitude is multiplied by scale factor
acube[*,*,1:2] = acube[*,*,1:2] * scale
aerrorcube[*,*,1:2] = aerrorcube[*,*,1:2] * scale
; Fit amplitude is integral of (flux/pixel)*micron, so multiply by
; pixels/micron to get total flux in erg/s/cm^2
acube[*,*,1:2] = acube[*,*,1:2] / wres
aerrorcube[*,*,1:2] = aerrorcube[*,*,1:2] / wres

;get header parameters(and add them to new fits files)
RA = sxpar(header,'RA')  ;RA at spatial [0,0] in mosaic ;format='(f14.10)'
DEC = sxpar(header,'DEC');DEC at spatial [0,0] in mosaic ;format='(f14.10)'
scale = sxpar(header,'CDELT2') ;scale
PA_Spec = sxpar(header,'PA_SPEC') ;position angle of spectrograph on sky
PA_IMAG = sxpar(header,'PA_IMAG') ; position angle of imager on sky
NAXIS_datacube = sxpar(header,'NAXIS*')
CRVAL2 = sxpar(header, 'CRVAL2') ;[deg] R.A. at reference pixel       
CRVAL3 = sxpar(header, 'CRVAL3') ;[deg] DEC at reference pixel       
CRPIX2 = sxpar(header, 'CRPIX2') ; Reference pixel location        
CRPIX3 = sxpar(header,'CRPIX3')  ;Reference pixel location        

mkhdr,header_acube,acube
;paraname = ['RA','DEC','scale','PA_SPEC','PA_IMAG','NAXIS_C1','NAXIS_C2','NAXIS_C3','CRVAL2','CRVAL3','CRPIX2','CRPIX3']

;writefits, 'output/cswa139_acube.fits', acube,header_acube
writefits, 'output/cswa139_Ha_acube.fits', acube

;2) do OII/Hb image

image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa139_Hb_Hbb_mosaic_scaledsky_130hr.fits",header)
for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the OIII wavelengths
oiii_sum = total(image[1502:1517,*,*],1)

imdisp, oiii_sum, out_pos=pos1

contour, oiii_sum, pos=pos1, /xs, /ys, /noerase,levels = [0.07,0.09,0.11]

;To plot wavelengths vs intensity map. Take only galaxy region. 
oiii_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do oiii_spectrum[i] = total(image[i,40:45,18:23])
wavelength = findgen(1651)*0.2+1473
plot,wavelength,oiii_spectrum,psym=1
vline,(z+1)*[486.2,500.8]

xmins = 30
xmaxs = 55
ymins = 15
ymaxs = 30
lineoiii = 500.8*(1.+z)
lineb2 = 486.2*(1.+z)
range = 10
sigmalim = 3. 
wres = wavelength[1]-wavelength[0] ; in nm per pixel

ok = where(wavelength ge lineb2-range and wavelength lt lineoiii+range)
finalcube = image[ok,*,*]
wl = wavelength[ok]

junk1 = where(wl lt lineb2-3.)
junk2 = where(wl gt lineb2+3. and wl le lineoiii-3.)
junk3 = where(wl gt lineoiii+3.)
;finalcube[[junk1,junk2,junk3],*,*]=0.

;rescale final cube 
; NOTE - IDL fit routine does not work with very small numbers
;        (i.e. cgs flux units ~1e-19), so make sure cube values are of
;        order unity before fitting spectra.
; Need to scale flux values to order unity
scale = abs(stddev(finalcube[where(finalcube)]))
finalcube = finalcube/scale


finalsz = size(finalcube)

weight = fltarr(finalsz[1])
for i = 0,finalsz[1]-1 do weight[i] = 1/variance(finalcube[i,xmins:xmaxs,ymins:ymaxs])

signal = 3.*stddev(finalcube[where(finalcube)])
offset = median(finalcube[where(finalcube)])

fitcube = 1./fltarr(finalsz[1],finalsz[2],finalsz[3]) ;holds best fit
acube   = 1./fltarr(finalsz[2],finalsz[3],5) ;holds best parameters
aerrorcube = 1./fltarr(finalsz[2],finalsz[3],4) ;holds sigma values
count = 0

for ii=xmins,xmaxs do for jj=ymins,ymaxs do if finalcube[0,ii,jj] then begin
   ; only look at pixels where zero-th wl bin is non zero
   x = wl
   y = finalcube[*,ii,jj]
   av = (moment(y))[0]
   chi0 = total(weight*(y-av)^2)

   a = [z,signal,signal,1.,offset]    ;[redshift,height1,height2,width,y offset]   
   fit = curvefit(x,y,weight,a,sigmafit,chisq = chi2,function_name = "doubletoiii",fita = [1,1,1,1,1],/noderivative)
   
   chi2 = total(weight*(y-fit)^2) ;chisquare of line fit
   chi2_var = chi0+1.-chi2-n_elements(a) ;variance of the line fit
   nsigma = sqrt(chi2_var)  ; significance of the line fit in standard deviations

   if nsigma gt sigmalim then begin
      plot, x, y, yrange=[min(y),max(y)],psym=1,title='('+String(ii)+String(jj)+')'
      oplot, x, fit, color=255
      print, stddev(y), median(1/sqrt(weight))
      fitcube[*,ii,jj] = fit
      acube[ii,jj,*] = [a[0:3], nsigma]
      aerrorcube[ii,jj,*] = sigmafit[0:3]
      print, 'fit found with ', nsigma, ' sigma at ', ii, jj
      print, 'chi^2: no line, line, sigma(line) = ', chi0, chi2, nsigma
      count = count+1
      wait, 0.5
   endif 
endif

print, count,'pixels out of',(xmaxs-xmins+1)*(ymaxs-ymins+1),'pixels have sigma fit >', sigmalim

; Turn fit parameters into actual physical properties
acube[where(acube eq 0)] = 1/0.
aerrorcube[where(aerrorcube eq 0)] = 1/0
; Convert linewidth in nm to rest-frame km/s
acube[*,*,3] = abs(acube[*,*,3])/lineHA* 3d5   ; sigma in km/s
aerrorcube[*,*,3] = abs(aerrorcube[*,*,3]) / lineHA * 3d5
; Append the measured redshift, z, to the 'acube' data
z_meas = acube[*,*,0] 
acube = [[[acube]], [[z_meas]]]
; Convert redshift to rest-frame km/s
acube[*,*,0] = (((acube[*,*,0] +1 )/(z+1))-1) * 3d5
aerrorcube[*,*,0] = aerrorcube[*,*,0] * 3d5
; Fit amplitude is multiplied by scale factor
acube[*,*,1:2] = acube[*,*,1:2] * scale
aerrorcube[*,*,1:2] = aerrorcube[*,*,1:2] * scale
; Fit amplitude is integral of (flux/pixel)*micron, so multiply by
; pixels/micron to get total flux in erg/s/cm^2
acube[*,*,1:2] = acube[*,*,1:2] / wres
aerrorcube[*,*,1:2] = aerrorcube[*,*,1:2] / wres

;get header parameters(and add them to new fits files)
RA = sxpar(header,'RA')  ;RA at spatial [0,0] in mosaic ;format='(f14.10)'
DEC = sxpar(header,'DEC');DEC at spatial [0,0] in mosaic ;format='(f14.10)'
scale = sxpar(header,'CDELT2') ;scale
PA_Spec = sxpar(header,'PA_SPEC') ;position angle of spectrograph on sky
PA_IMAG = sxpar(header,'PA_IMAG') ; position angle of imager on sky
NAXIS_datacube = sxpar(header,'NAXIS*')
CRVAL2 = sxpar(header, 'CRVAL2') ;[deg] R.A. at reference pixel       
CRVAL3 = sxpar(header, 'CRVAL3') ;[deg] DEC at reference pixel       
CRPIX2 = sxpar(header, 'CRPIX2') ; Reference pixel location        
CRPIX3 = sxpar(header,'CRPIX3')  ;Reference pixel location        

mkhdr,header_acube,acube
;paraname = ['RA','DEC','scale','PA_SPEC','PA_IMAG','NAXIS_C1','NAXIS_C2','NAXIS_C3','CRVAL2','CRVAL3','CRPIX2','CRPIX3']

;writefits, 'output/cswa139_acube.fits', acube,header_acube
writefits, 'output/cswa139_OII_acube.fits', acube


;Calculate metallicity
Halpha = acube[*,*,1]
NII = acube[*,*,2]
;1) N2 method




end
