
; FITSPEC_GENERAL
; IDL function to fit spectra of individual spaxels in a data cube with Gaussian line profiles. Designed to generate emission line maps.
;
; CALLING SEQUENCE:
; FITSPEC_GENERAL, FINALCUBE, WL, LINE, Z, SIGMALIM, XMINS, XMAXS, YMINS, YMAXS, NMIN, NI, NMAX, VMAX, WMAX
;
; INPUT:
; FINALCUBE     - data cube. Dimensions must be [lambda,x,y].
; WL            - wavelength array (microns) corresponding to the lambda dimension of FINALCUBE.
; LINE          - rest frame wavelength of the spectroscopic line to be fit, in microns (e.g. 0.6563 for Halpha).
; Z             - approximate redshift of LINE.
; SIGMALIM      - minimum significance required for acceptable fit (standard deviations).
; XMINSÉYMAXS - x,y min/max pixel values of a blank sky region, used to measure the sky spectrum. This determines the Gaussian weights for fitting each spaxel. The sky region used is FINALCUBE[*,XMINS:XMAXS,YMINS:YMAXS].
; NMIN        - initial binning factor. Use 0 for no initial binning.
; NI          - increase in binning factor in each iteration.
; NMAX        - highest allowed binning factor. Use 0 for no binning.
; VMAX        - maximum allowed absolute velocity offset from Z (km/s). Best-fit lines with centroids offset by a value >VMAX will be rejected.
; WMAX        - maximum allowed velocity dispersion (km/s). Best-fit lines with dispersions >WMAX will be rejected.
;
; OUTPUT:
; Writes two FITS files to disk.
; acube.fits      - contains two-dimensional arrays of best-fit [velocity centroid (km/s), amplitude (same units as FINALCUBE), velocity dispersion (km/s), binning factor, fit significance (# of standard deviations), redshift].
; aerrorcube.fits - contains two-dimensional arrays of 1-sigma uncertainty on the best-fit [velocity, amplitude, velocity dispersion].
;
; EXAMPLE:
; To fit Halpha emission at redshift z=2.20 in an OSIRIS Kn2-band data cube, requiring a 5-sigma line detection, with no binning:
; IDL> kcube = readfits('s130911_a029001_Kn2_100.fits')
; IDL> wl = 2.036 + 0.00025*findgen(421)
; IDL> fitspec_general, kcube, wl, .6563, 2.2275, 5.0, 7,13,29,42, 0, 1, 0, 500, 500,'s130911_a029001_Kn2_100'
;
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

pro ha_nii, x, par, f, pder

; Fit Halpha and NII lines in the same spectrum with same redshift, linewidth
; x = wavelength
; par = [x0(halpha), A(halpha), A(NII), linewidth, continuum]

xha = par[0]              ; Halpha line
xn2 = par[0]*.6583/.6563  ; NII line
aha = par[1]
an2 = par[2]
ww  = par[3]
cont= par[4]

f = (aha/ww/sqrt(2.0*!pi))*exp(-0.5*(x-xha)^2/ww^2) + (an2/ww/sqrt(2.0*!pi))*exp(-0.5*(x-xn2)^2/ww^2) + cont

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


; ----------------------------------------------------------------
pro fitspec_general, finalcube, wl, line, z, sigmalim, xmins, xmaxs, ymins, ymaxs, nmin, ni, nmax, vmax, wmax, name=name,imgymins,imgymaxs,header

;;; Note, the instrumental resolution is not accounted for in
;;; determining line widths.

; Fit lines in a data cube using iterative binning. Adapted from fitspec_0712.
; finalcube = data cube. Dimensions must be [lambda,x,y]
; wl = wavelength in microns
; line = line wavelength in microns
; z = redshift of line
; sigmalim = significance required for acceptable fit (standard deviations)
; xmins...ymaxs = x,y min/max pixel values of blank sky region, used
; for calculating Gaussian weights for fitting.
; nmin = initial binning factor
; ni = increase in binning factor in each iteration
; nmax = highest allowed binning factor
; vmax = maximum allowed velocity offset (km/s)
; wmax = maximum allowed velocity dispersion (km/s)
; name = name of output files
; imgymins/max = y of positive image (the number of image that you click to move) 




; Transpose the input data cube
; (for historical reasons, the script was originally written to take an input cube with dimensions [x,y,lambda])
finalcube = transpose(finalcube, [1,2,0])

; Parameters (may be tweaked)
wres = wl[1]-wl[0] ; wavelength resolution (in micron/pixel) of data
range = .005       ; delta wavelength in which to fit emission
sz = size(finalcube)
xmin = 0      ; minimum x pixel value for extracting spectra
xmax = sz[1]-1
ymin = 0
ymax = sz[2]-1

linewave = line * (1.+z)

; NOTE - IDL fit routine does not work with very small numbers
;        (i.e. cgs flux units ~1e-19), so make sure cube values are of
;        order unity before fitting spectra.
; Need to scale flux values to order unity
scale = abs(stddev(finalcube[where(finalcube)]))
finalcube = finalcube / scale

; Subtract any continuum in the spectrum.
;sz = size(finalcube)
;for j=0,sz[2]-1 do for i=0,sz[1]-1 do finalcube[i,j,*] = $
;  finalcube[i,j,*] - median(finalcube[i,j,*])

; cut cube in wavelength for speed
ok = where(wl ge linewave-range and wl le linewave+range)
finalcube = finalcube[*,*,ok]
wl = wl[ok]
sz = size(finalcube)

; Determine weights for spectrum fitting. Use Gaussian weights
; w = 1/variance.
; Note that these are weights for one pixel, so when binning the
; variance V becomes V/(number of pixels binned)
; Therefore use weight = w*(number of pixels binned)
print, 'Calculating weights from sky spectrum'
w = fltarr(sz[3])
for i=0,sz[3]-1 do w[i] = 1/variance(finalcube[xmins:xmaxs,ymins:ymaxs,i])

fitcube = fltarr(sz[1],sz[2],sz[3]) ; holds best fit 
acube = fltarr(sz[1],sz[2],5)   ; holds best parameters
aerrorcube = fltarr(sz[1],sz[2],3) ; holds +- 1 sigma values for velocity, intensity, line width

; first guess parameters
linewidth = 0.0003
signal = 3.*stddev(finalcube)
offset = median(finalcube[where(finalcube)])
alast = [linewave, signal*linewidth*sqrt(2*!pi), linewidth, offset]


print, 'Fitting spectra'
for j=ymin,ymax-1 do for i=xmin,xmax-1 do if finalcube[i,j,0] then begin
    ; Only look at pixels where zero-th wavelength bin is non-zero
    x = wl
    y = fltarr(sz[3])
    var = fltarr(sz[3])

    nsum = nmin    ; initial binning
    nsigma = 0
    while nsum le nmax AND nsigma lt sigmalim do begin
;        print, 'nsum = ', nsum, ' and nsigma = ', nsigma, ' at pixel ', [i,j]
        ; Iteratively bin until good enough fit is found
        ; bin 2*nsum+1 by 2*nsum+1 pixels
        for k=0,sz[3]-1 do y[k] = mean(finalcube[(i-nsum)>0:((i+nsum)<(sz[1]-1)),$
                                                 (j-nsum)>0:((j+nsum)<(sz[2]-1)),k])

        weight = w * (2*nsum+1)^2.  ; adjust variance for binning
;        weight[*] = 1/variance(y)    ; assumes variance is equal everywhere

        av = (moment(y))[0]
        chi0 = total(weight*(y - av)^2)  ; chi-square of fit with no line
        df0 = n_elements(y) - 1.    ; degrees of freedom for chi0
        sigma = sqrt(2.*df0) ; standard dev of chi-square distribution
        print, 'chi-square, mean, and sigma (no line fit):', chi0, df0, sigma

        a = alast
        fit = curvefit(x,y,weight,a,sigmafit,chisq=chi2, $
                       function_name="singlet",/noderivative, fita=[1,1,1,1])
                       
;        guess = [signal, linewave, linewidth, av]
;        fit = gaussfit(x,y,a,chisq=chi2,estimates=guess,nterms=4)
;        print, chi2
        chi2 = total(weight*(y - fit)^2) ; chi-square of line fit
        chi2_var = chi0+1. - chi2 - n_elements(a)    ; variance of the line fit
        nsigma = sqrt(chi2_var)    ; significance of line fit, in standard devs

;        print, 'chi^2: no line, line, sigma(line) = ', chi0, chi2, nsigma
;        print, 'Line fit significance is ', nsigma, ' sigma'

;        plot, x, y
;        oplot, x, fit, color=255
;        wait, 0.1

 
        if (nsigma gt sigmalim) then begin
            plot, x, y, yrange=[min(y),max(y)]
            oplot, x, fit, color=255
            print, stddev(y), median(1/sqrt(weight))
            fitcube[i,j,*] = fit
            acube[i,j,*] = [a[0:2], 2*nsum+1., nsigma]
;            acube[i,j,*] = [a[0:2], 2*nsum+1, chi2]
            aerrorcube[i,j,*] = sigmafit[0:2]
;            print, 'best parameters = ', a
            print, 'fit found with ', nsigma, ' sigma at ', i, j
            print, 'chi^2: no line, line, sigma(line) = ', chi0, chi2, nsigma
        endif
        
        nsum = nsum+ni
    endwhile
endif


; Turn fit parameters into actual physical properties
acube[where(acube eq 0)] = 1/0.   ; set bad fit values to NaN
aerrorcube[where(aerrorcube eq 0)] = 1/0.
; Convert sigma in micron to rest-frame km/s
acube[*,*,2] = abs(acube[*,*,2]) / linewave * 3d5   ; sigma in km/s
aerrorcube[*,*,2] = abs(aerrorcube[*,*,2]) / linewave * 3d5
; Append the measured redshift, z, to the 'acube' data
z_meas = acube[*,*,0] / line - 1.0
acube = [[[acube]], [[z_meas]]]
; Convert velocity in micron to rest-frame km/s
acube[*,*,0] = (acube[*,*,0] - linewave) / linewave * 3d5
aerrorcube[*,*,0] = aerrorcube[*,*,0] / linewave * 3d5
; Fit amplitude is multiplied by scale factor
acube[*,*,1] = acube[*,*,1] * scale
aerrorcube[*,*,1] = aerrorcube[*,*,1] * scale
; Fit amplitude is integral of (flux/pixel)*micron, so multiply by
; pixels/micron to get total flux in erg/s/cm^2
acube[*,*,1] = acube[*,*,1] / wres
aerrorcube[*,*,1] = aerrorcube[*,*,1] / wres


; Impose some constraints for legitimate fits
;vmax = 200
;wmax = 200
for i=xmin,xmax do for j=ymin,ymax do begin
    ; limit velocity range to \pm vmax km/s
    if abs(acube[i,j,0]) gt vmax then begin
        acube[i,j,*]=1/0.
        aerrorcube[i,j,*] = 1/0.
    endif
    ; require amplitude > 0 for emission line
    if acube[i,j,1] lt 0 then begin
        acube[i,j,*]=1/0.
        aerrorcube[i,j,*] = 1/0.
    endif
    ; limit line width range to < wmax km/s
    if acube[i,j,2] gt wmax then begin
        acube[i,j,*]=1/0.
        aerrorcube[i,j,*] = 1/0.
    endif
endfor



; Make acube_obj (only galaxy portion of the whole image)
acube_obj = fltarr(sz[1],imgymaxs-imgymins+1,6)
aerrorcube_obj = fltarr(sz[1],imgymaxs-imgymins+1,3)
fitcube_obj = fltarr(sz[1],imgymaxs-imgymins+1,sz[3])
finalcube_obj = fltarr(sz[1],imgymaxs-imgymins+1,sz[3])
for i=0,sz[1]-1 do for j=0,imgymaxs-imgymins do  begin
   acube_obj(i,j,*) = acube(i,imgymins+j,*)
   aerrorcube_obj(i,j,*) = aerrorcube(i,imgymins+j,*)
   fitcube_obj(i,j,*) = fitcube(i,imgymins+j,*)
   finalcube_obj(i,j,*) = finalcube(i,imgymins+j,*)
   
endfor
save,wl,finalcube_obj,acube_obj,aerrorcube_obj,fitcube_obj,file='output/'+name+'.sav'
;=============================================================
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
sxaddpar,header_acube,'RA',RA,'RA at spatial [0,0] in mosaic'
sxaddpar,header_acube,'DEC',dec,'DEC at spatial [0,0] in mosaic'
sxaddpar,header_acube,'scale',scale,'plate scale in degree per pixel'
sxaddpar,header_acube,'PA_SPEC',PA_SPEC,'position angle of spectrograph on sky'
sxaddpar,header_acube,'PA_IMAG',PA_IMAG,'position angle of imager on sky'
sxaddpar,header_acube,'NAXIS_C1',naxis_datacube(0),'sizes of original mosaic datacube'
sxaddpar,header_acube,'NAXIS_C2',naxis_datacube(1),'sizes of original mosaic datacube'
sxaddpar,header_acube,'NAXIS_C3',naxis_datacube(2),'sizes of original mosaic datacube'
sxaddpar,header_acube,'LINE',line,'rest frame wavelength of the spectroscopic line fit in this image'
sxaddpar,header_acube,'SIGMALIM',sigmalim,'minimum significance required for acceptable fit (standard deviations)'
sxaddpar,header_acube,'CRVAL2',CRVAL2,'[deg] R.A. at reference pixel '
sxaddpar,header_acube,'CRVAL3',CRVAL3,'[deg] DEC at reference pixel '
sxaddpar,header_acube,'CRPIX2',CRPIX2,'Reference pixel location'
sxaddpar,header_acube,'CRPIX3',CRPIX3,'Reference pixel location'

mkhdr,header_acube_obj,acube_obj
sxaddpar,header_acube_obj,'RA',RA,'RA at spatial [0,0] in mosaic'
sxaddpar,header_acube_obj,'DEC',dec,'DEC at spatial [0,0] in mosaic'
sxaddpar,header_acube_obj,'scale',scale,'plate scale in degree per pixel'
sxaddpar,header_acube_obj,'PA_SPEC',PA_SPEC,'position angle of spectrograph on sky'
sxaddpar,header_acube_obj,'PA_IMAG',PA_IMAG,'position angle of imager on sky'
sxaddpar,header_acube_obj,'NAXIS_C1',naxis_datacube(0),'sizes of originall mosaic datacube'
sxaddpar,header_acube_obj,'NAXIS_C2',naxis_datacube(1),'sizes of original mosaic datacube'
sxaddpar,header_acube_obj,'NAXIS_C3',naxis_datacube(2),'sizes of original mosaic datacube'

sxaddpar,header_acube_obj,'pix_y1',imgymins,'y pixel# in acube.fits of y=0 in this image'
sxaddpar,header_acube_obj,'pix_y2',imgymaxs,'y pixel# in acube.fits of y=NAXIS2 in this image'
sxaddpar,header_acube_obj,'LINE',line,'rest frame wavelength of the spectroscopic line fit in this image'
sxaddpar,header_acube_obj,'SIGMALIM',sigmalim,'minimum significance required for acceptable fit (standard deviations)'
sxaddpar,header_acube_obj,'CRVAL2',CRVAL2,'[deg] R.A. at reference pixel '
sxaddpar,header_acube_obj,'CRVAL3',CRVAL3,'[deg] DEC at reference pixel '
sxaddpar,header_acube_obj,'CRPIX2',CRPIX2,'Reference pixel location'
sxaddpar,header_acube_obj,'CRPIX3',CRPIX3,'Reference pixel location'
;===================================================


writefits, 'output/'+name+'_acube.fits', acube,header_acube
writefits, 'output/'+name+'_aerrorcube.fits', aerrorcube
writefits, 'output/'+name+'_acube_obj.fits', acube_obj,header_acube_obj
writefits, 'output/'+name+'_aerrorcube_obj.fits', aerrorcube_obj

end

; ----------------------------------------------------------------
