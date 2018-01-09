; FITSPEC_INTENSITY
;Nicha editted april 2014
; IDL function to fit spectra of individual spaxels in a data cube with Gaussian line profiles. Designed to generate emission line maps in conjunction with "fitspec_general.pro".
;
; CALLING SEQUENCE:
; FITSPEC_INTENSITY, CUBE, WL, LINE, Z, DISP, NBIN, XMINS, XMAXS, YMINS, YMAXS
;
; INPUT:
; CUBE          - data cube. Dimensions must be [lambda,x,y].
; WL            - wavelength array corresponding to the lambda dimension of CUBE.
; LINE          - rest frame wavelength of the spectroscopic line to be fit, in same units as WL (e.g. 0.6563 for Halpha in microns).
; Z             - grid of redshifts for each spaxel in the data cube. Same dimensions as the data cube's [x,y].
; DISP          - grid of velocity dispersion for each spaniel in the data cube. Same dimensions as the data cube's [x,y], with same units as wl and line.
; NBIN          - grid of binning factor to use for each pixel position. Same dimensions as the data cube's [x,y]. The data cube is binned in square boxes of size 2*n+1 centered on the pixel being fit.
; XMINSÉYMAXS - x,y min/max pixel values of a blank sky region, used to measure the sky spectrum. This determines the Gaussian weights for fitting each spaxel. The sky region used is FINALCUBE[*,XMINS:XMAXS,YMINS:YMAXS].
;
; OUTPUT:
; Writes a FITS file to disk.
; output.fits    - contains two-dimensional arrays of [best-fit amplitude, 1-sigma uncertainty in amplitude, significance level of fit (# of standard deviations)].
;
; EXAMPLE:
; This program is designed to be used in conjunction with "fitspec_general.pro", where FITSPEC_GENERAL is used to fit the kinematics (precise redshift and line width) for a strong emission line, and FITSPEC_INTENSITY is then used to measure the intensity of a fainter emission line using the same kinematics.
; To fit Halpha (6562.79 Angstroms) emission at redshift z=2.20 in an OSIRIS Kn2-band data cube, requiring a 5-sigma line detection, with no binning:
; IDL> kcube = readfits('s130911_a029001_Kn2_100.fits')
; IDL> wl = 2.036 + 0.00025*findgen(421)
; IDL> fitspec_general, kcube, wl, .656279, 2.20, 5.0, 7,13,29,42, 0, 1, 0, 500, 500
; And to apply the best-fit Halpha kinematics to [N II] 6583.45 Angstrom line:
; IDL> kcube = readfits('s130911_a029001_Kn2_100.fits')
; IDL> wl = 2.036 + 0.00025*findgen(421)
; IDL> ha_fit = readfits('acube.fits')
; IDL> z_cube = ha_fit[*,*,5]
; IDL> disp_cube = (ha_fit[*,*,2] / 2.998d5) * z_cube * .658345   ; convert velocity dispersion from km/s to microns
; IDL> bin_cube = ha_fit[*,*,3]
; IDL> fitspec_intensity, kcube, wl, .658345, z_cube, disp_cube, bin_cube, 7,13,29,42
;
; ----------------------------------------------------------------

pro singlet, x, par, f, pder

x0 = par[0]
aa = par[1]
ww = par[2]
cont = par[3]

f = (aa/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-x0)^2/ww^2) + cont

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


; ----------------------------------------------------------------

pro fitspec_intensity_kmos, cube, wl, line, z, disp, nbin, xmins, xmaxs, ymins, ymaxs, name=name,header,xmin,xmax,ymin,ymax
COMMON cube, acube,aerrorcube
window,0,xsize=400,ysize=800
!p.multi=[0,1,2]

finalcube = transpose(cube, [1,2,0])
cube = finalcube

; Parameters (may be tweaked)
wres = wl[1]-wl[0] ; wavelength resolution (in micron/pixel) of data
range = 0.004    ; delta wavelength range to fit
range_low = range
range_high = range
sz = size(cube)


linewave = line * (1. + median(z[where(finite(z))]))   ; median wavelength

if name eq 'cswa165_Ha_Kn2_handmosaic_sky_1hr' then begin
   z(34,34) = 2.12824
   range_low = 0.0025 ;change the range because there is a bright sky line at 2.057 that increases the background level in line fitting and resulting in faulty negative NII lines. 
endif 

if name eq 'cswa19_Ha_Kn1_mosaic_scaledsky_4hr' then range_high = 0.002

linewaves = line * (1. + z)   ; array of wavelength values

; cut cube in wavelength for speed

ok = where(wl ge linewave-range_low and wl le linewave+range_high) 

finalcube = cube[*,*,ok]
wl = wl[ok]
sz = size(finalcube)
stop
; NOTE - IDL fit routine does not work with very small numbers
;        (i.e. cgs flux units ~1e-19), so make sure cube values are of
;        order unity before fitting spectra.
; Need to scale flux values to order unity
scale = abs(stddev(finalcube[where(finalcube)]))
finalcube = finalcube / scale

;restore the fits of the first line (either Ha or OIII) to pit it below the second line fit in comparison
restore,'fitstruct.sav'

; Determine weights for spectrum fitting. Use Gaussian weights
; w = 1/variance.
; Note that these are weights for one pixel, so when binning the
; variance V becomes V/(number of pixels binned)
; Therefore use weight = w*(number of pixels binned)
print, 'Calculating weights from sky spectrum'
w = fltarr(sz[3])
for i=0,sz[3]-1 do w[i] = 1/variance(finalcube[xmins:xmaxs,ymins:ymaxs,i])
plot, wl,w

oplot, wl,w,linestyle=2, color=255
plot,[0,1],[1,1],psym=1 ;nothing. Just to fill the plot space so that the next ones come in pairs
;stop

; Create variable in which to write the intensity and uncertainty.
; Store 3 values at each point: intensity, error on intensity, and
; significance of fit.
output = make_array(sz[1],sz[2],3)
fitcube = fltarr(sz[1],sz[2],sz[3]) ; holds best fit 

; Do the fit!
print, 'Fitting spectra'
fitnegativestr = []
fitpositivestr = []
for j=ymin,ymax-1 do for i=xmin,xmax-1 do begin 
   if finite(z[i,j]) then begin
                                ; Only look at pixels where a redshift has been measured
      x = wl
      y = fltarr(sz[3])
      nsum = (nbin[i,j]-1.)/2.  ; binning factor
      
                                ; Determine the spectrum to be fit
      for k=0,sz[3]-1 do y[k] = mean(finalcube[(i-nsum)>0:((i+nsum)<(sz[1]-1)),$
                                             (j-nsum)>0:((j+nsum)<(sz[2]-1)),k])
      
      weight = w * (2*nsum+1)^2. ; adjust variance for binning

      av = (moment(y))[0]
      chi0 = total(weight*(y - av)^2) ; chi-square of fit with no line
      df0 = n_elements(y) - 1.        ; degrees of freedom for chi0
      sigma = sqrt(2.*df0)            ; standard dev of chi-square distribution
      print, 'chi-square, mean, and sigma (no line fit):', chi0, df0, sigma
      if n_elements([chi0,df0,sigma]) gt 3 then stop
    
    ; initial guess of fit parameters:
      a = [linewaves[i,j], 3.*stddev(y)*disp[i,j], disp[i,j], median(y)]
      fit_a = [0,1,0,1]         ; fit only the intensity and continuum
    
      fit = curvefit(x,y,weight,a,sigmafit,chisq=chi2, fita=fit_a, $
                   function_name="singlet",/noderivative)

      chi2 = total(weight*(y - fit)^2)           ; chi-square of line fit
      chi0 = total(weight*(y-a(3))^2)
      chi2_var = chi0+1. - chi2   ; variance of the line fit
      nsigma = sqrt(chi2_var)   ; significance of line fit, in standard devs
      print, 'chi^2: no line, line, sigma(line) = ', chi0, chi2, nsigma
    
    ;if i ge 20 and i le 23 and j ge 33 and j le 35 then stop

    ; Write the best-fit results to the 'output' array
      fitcube[i,j,*] = fit
      output[i,j,*] = [a[1], sigmafit[1], nsigma] ; [intensity, error on intensity, significance]

      plot, x, y*scale, yrange=[min(y*scale),max(y*scale)],psym=10,title=string(i)+','+string(j)
      oplot, x, fit*scale, color=255
      oplot,[linewaves(i,j),linewaves(i,j)],!y.crange,linestyle=2

      ;Plot the first line in comparison
      position = where(fitstruct.x eq i and fitstruct.y eq j)
      npositions = n_elements(position)
      if npositions gt 1 then position = position(npositions-1)
      plot, fitstruct(position).wl,fitstruct(position).data,psym=10,title=string(fitstruct(position).sigma)
      oplot,fitstruct(position).wl,fitstruct(position).fit,color=255
      print, 'redshift=',z(i,j)
      ;read, doesntmatter,prompt='press any number to continue'
      ;wait,.5
      if a[1] lt 0. then begin
         fitnegative={x:i,y:j,wl:x,data:y,fit:fit,sigma:nsigma}
         fitnegativestr = [fitnegativestr,fitnegative]
         ;stop                   ; stop if the line is negative
      endif

      if a[1] gt 0. then begin
         fitpositive={x:i,y:j,wl:x,data:y,fit:fit,sigma:nsigma}
         fitpositivestr = [fitpositivestr,fitpositive]
         ;stop                   ; stop if the line is negative
      endif

   endif
endfor
save,fitnegativestr,filename='negative_2ndline.sav'
save,fitpositivestr,filename='positive_2ndline.sav'

; Convert best-fit intensity into physical units
   ;set bad fit values to Nan
output[where(output eq 0)] = 1/0.
   ; Fit amplitude is multiplied by scale factor
output[*,*,0] = output[*,*,0] * scale
output[*,*,1] = output[*,*,1] * scale
   ; Fit amplitude is integral of (flux/pixel)*micron, so multiply by
   ; pixels/micron to get total flux in erg/s/cm^2
output[*,*,0] = output[*,*,0] / wres
output[*,*,1] = output[*,*,1] / wres
acube = [[[acube]],[[output[*,*,0]]]]
acube = [[[acube]],[[output[*,*,2]]]]

aerrorcube = [[[aerrorcube]],[[output[*,*,1]]]]
stop

end

; ----------------------------------------------------------------
