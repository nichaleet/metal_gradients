; CSWA11FITSPEC_GENERAL
;Tucker's original
;Nicha modified oct 2014
; IDL function to fit spectra of individual spaxels in a data cube with Gaussian line profiles. Designed to generate emission line maps for cswa11 with 2 merging components. I.e. fit 2 gaussians of Ha at the same time
;
; CALLING SEQUENCE:
; FITSPEC_GENERAL, finalcube, wl, line, z, sigmalim, xmins, xmaxs, ymins, ymaxs, nmin, ni, nmax, vmax, wmax, name=name,xmin,xmax,ymin,ymax,header
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
; acube.fits      - contains two-dimensional arrays of best-fit [velocity centroid (km/s), amplitude (same units as FINALCUBE), velocity dispersion (km/s),velocity centroid (km/s) of the 2nd line, amplitude of the 2nd line (same units as FINALCUBE), velocity dispersion (km/s)of the 2nd line, binning factor, fit significance (# of standard deviations), redshift of the 1st line, redshift of the 2nd line].
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

pro doublet, x, par, f, pder

x0 = par[0]
aa = par[1]
ww = par[2]
cont = par[3]

x0_1 = par[4]
aa_1 = par[5]
ww_1 = par[6]


f = (aa/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-x0)^2/ww^2)+(aa_1/ww_1/sqrt(2.0*!pi)) * exp(-0.5*(x-x0_1)^2/ww_1^2) + cont

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


pro cswa11fitspec_general, finalcube, wl, line, z, sigmalim, xmins, xmaxs, ymins, ymaxs, nmin, ni, nmax, vmax, wmax, name=name,xmin,xmax,ymin,ymax,header,header_acube,crpix2,crpix3
Common cube,acube,aerrorcube

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
; xmin/xmax = minimum/maximum x pixel value for extracting spectra
; ymin/ymax = minimum/maximum y pixel value for extracting spectra



; Transpose the input data cube
; (for historical reasons, the script was originally written to take an input cube with dimensions [x,y,lambda])
finalcube = transpose(finalcube, [1,2,0])

; Parameters (may be tweaked)
wres = wl[1]-wl[0] ; wavelength resolution (in micron/pixel) of data
range = .004       ; delta wavelength in which to fit emission original = 0.005
sz = size(finalcube)


linewave = 1.581
linewave1 = 1.580
linewave2 = 1.582 

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

acube = fltarr(sz[1],sz[2],8)   ; holds best parameters
aerrorcube = fltarr(sz[1],sz[2],6) ; holds +- 1 sigma values for velocity, intensity, line width

; first guess parameters
linewidth = 0.0004
signal = 3.*stddev(finalcube)
offset = median(finalcube[where(finalcube)])

;alast = [linewave, signal*linewidth*sqrt(2*!pi), linewidth, offset]
alast = [1.580,0.5,linewidth,offset,1.582,0.3,linewidth]
plot, wl,w,ytitle='weight',xtitle='micron'
;stop

oplot, wl,w,linestyle=2,color=255
stop
print, 'Fitting spectra'

;fitstruct is a variable to hold a fit data
fitstruct = []

for j=ymin,ymax-1 do for i=xmin,xmax-1 do if finalcube[i,j,0] then begin
    ; Only look at pixels where zero-th wavelength bin is non-zero
    x = wl
    y = fltarr(sz[3])
    var = fltarr(sz[3])
    ww= w

    nsum = nmin    ; initial binning
    nsigma = 0
    sigmalimnow = sigmalim
    
    while nsum le nmax AND nsigma lt sigmalimnow+2.*nsum do begin
;        print, 'nsum = ', nsum, ' and nsigma = ', nsigma, ' at pixel ', [i,j]
        ; Iteratively bin until good enough fit is found
        ; bin 2*nsum+1 by 2*nsum+1 pixels
        for k=0,n_Elements(y)-1 do y[k] = mean(finalcube[(i-nsum)>0:((i+nsum)<(sz[1]-1)),$
                                                 (j-nsum)>0:((j+nsum)<(sz[2]-1)),k])


        weight = ww * (2*nsum+1)^2.  ; adjust variance for binning
;        weight[*] = 1/variance(y)    ; assumes variance is equal everywhere

        av = (moment(y))[0]
        chi0 = total(weight*(y - av)^2)  ; chi-square of fit with no line
        df0 = n_elements(y) - 1.    ; degrees of freedom for chi0
        sigma = sqrt(2.*df0) ; standard dev of chi-square distribution
        print, 'chi-square, mean, and sigma (no line fit):', chi0, df0, sigma

        a = alast
        fita=[1,1,1,1,1,1,1]

        fit = curvefit(x,y,weight,a,sigmafit,chisq=chi2, $
                       function_name="doublet",/noderivative, fita=fita)
       
        chi2 = total(weight*(y - fit)^2) ; chi-square of line fit
        chi2_var = chi0-1. - chi2 + n_elements(a)    ; variance of the line fit
        nsigma = sqrt(chi2_var)    ; significance of line fit, in standard devs

        ;if (nsigma gt sigmalimnow+5.*nsum) then begin
        if (nsigma gt sigmalimnow) then begin
            plot, x, y*scale, yrange=[min(y*scale),max(y*scale)],psym=10,title='First line '+String(i)+String(j)+string(nsum)+string(fix(nsigma))
            oplot, x, fit*scale, color=255
            currentfit = {x:i,y:j,wl:x,data:y*scale,fit:fit*scale,sigma:nsigma}
            fitstruct = [fitstruct,currentfit]
               
            print, stddev(y), median(1/sqrt(weight))
            acube[i,j,*] = [a[0:2],a[4:6], 2*nsum+1., nsigma]
            aerrorcube[i,j,*] = [sigmafit[0:2],sigmafit[4:6]]
            print, 'best parameters = ', a
            print, 'fit found with ', nsigma, ' sigma at ', i, j
            print, 'chi^2: no line, line, sigma(line) = ', chi0, chi2, nsigma
            print, 'nbins = ' ,2*nsum+1
            ;if nsum gt 0 then wait, 0.5
            ;wait,1 
           endif
         ;To check why they don't fit properly for specific galaxies
        if a(1)*scale/wres gt .02 and abs(a(2))/linewave*3d5 lt 200. and abs((a(0)-linewave)/linewave*3d5) lt 300. and j lt 44 then begin
           print, a(1)*scale/wres
           print, abs((a(0)-linewave)/linewave*3d5)
           stop
        endif
        nsum = nsum+ni
     endwhile
 endif

save,fitstruct,filename='fitstruct.sav'


; Turn fit parameters into actual physical properties
acube[where(acube eq 0)] = 1./0.   ; set bad fit values to NaN
aerrorcube[where(aerrorcube eq 0)] = 1./0.
; Convert sigma in micron to rest-frame km/s
acube[*,*,2] = abs(acube[*,*,2]) / linewave * 3d5   ; sigma in km/s
aerrorcube[*,*,2] = abs(aerrorcube[*,*,2]) / linewave * 3d5
acube[*,*,5] = abs(acube[*,*,5]) / linewave * 3d5   ; sigma in km/s
aerrorcube[*,*,5] = abs(aerrorcube[*,*,5]) / linewave * 3d5

; Append the measured redshift, z, to the 'acube' data
z_meas = acube[*,*,0] / line - 1.0
acube = [[[acube]], [[z_meas]]]  ;frame number 8 (start from 0)
z_meas = acube[*,*,3] / line - 1.0
acube = [[[acube]], [[z_meas]]]  ;frame number 9

; Convert velocity in micron to rest-frame km/s
acube[*,*,0] = (acube[*,*,0] - linewave) / linewave * 3d5
aerrorcube[*,*,0] = aerrorcube[*,*,0] / linewave * 3d5
acube[*,*,3] = (acube[*,*,3] - linewave) / linewave * 3d5
aerrorcube[*,*,3] = aerrorcube[*,*,3] / linewave * 3d5
; Fit amplitude is multiplied by scale factor
acube[*,*,1] = acube[*,*,1] * scale
aerrorcube[*,*,1] = aerrorcube[*,*,1] * scale
acube[*,*,4] = acube[*,*,4] * scale
aerrorcube[*,*,4] = aerrorcube[*,*,4] * scale
; Fit amplitude is integral of (flux/pixel)*micron (area under gaussian curve), so multiply by
; pixels/micron to get total flux in erg/s/cm^2
acube[*,*,1] = acube[*,*,1] / wres
aerrorcube[*,*,1] = aerrorcube[*,*,1] / wres
acube[*,*,4] = acube[*,*,4] / wres
aerrorcube[*,*,4] = aerrorcube[*,*,4] / wres


; Impose some constraints for legitimate fits
vmax = 400
wmax = 200
for i=xmin,xmax do for j=ymin,ymax do begin
    ; limit velocity range to \pm vmax km/s
   if abs(acube[i,j,0]) gt vmax then begin
      acube[i,j,0:2]=1/0.
      acube[i,j,8]=1/0.
      aerrorcube[i,j,0:2] = 1/0.
   endif
   if abs(acube[i,j,3]) gt vmax then begin
      acube[i,j,3:5]=1/0.
      acube[i,j,9]=1/0.
      aerrorcube[i,j,3:5] = 1/0.
   endif
   ; if both amplitudes are negative then everything is nan
   if acube[i,j,1] lt 0 and acube[i,j,4] lt 0 then begin
      acube[i,j,*]=1/0.
      aerrorcube[i,j,*] = 1/0.
    endif
    ; require amplitude > 0 for the first emission line
   if acube[i,j,1] lt 0 then begin
      acube[i,j,0:2]=1/0.
      acube[i,j,8]=1/0.
      aerrorcube[i,j,0:2] = 1/0.
   endif
; require amplitude > 0 for the second emission line
    if acube[i,j,4] lt 0 then begin
        acube[i,j,3:5]=1/0.
        acube[i,j,9]=1/0.
        aerrorcube[i,j,3:5] = 1/0.
    endif
    ; limit line width range to < wmax km/s
    if acube[i,j,2] gt wmax then begin
       acube[i,j,0:2]=1/0.
       acube[i,j,8]=1/0.
       aerrorcube[i,j,0:2] = 1/0.
     endif
    if acube[i,j,5] gt wmax then begin
        acube[i,j,3:5]=1/0.
        acube[i,j,9]=1/0.
        aerrorcube[i,j,3:5] = 1/0.
    endif
endfor


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
;CRPIX2 = sxpar(header, 'CRPIX2') ; Reference pixel location        
;CRPIX3 = sxpar(header,'CRPIX3')  ;Reference pixel location   
CDELT2 = sxpar(header, 'CDELT2') ; Pixel scale in degree/pixel
CDELT3 = sxpar(header, 'CDELT3') ; Pixel scale in degree/pixel
PC2_2    =  sxpar(header, 'PC2_2')  ;/RA, Dec axes rotated by  degr.           
PC2_3    =  sxpar(header, 'PC2_3')  ;/RA, Dec axes rotated by  degr.           
PC3_2    =  sxpar(header, 'PC3_2')  ;/RA, Dec axes rotated by  degr.           
PC3_3    =  sxpar(header, 'PC3_3')  ;/RA, Dec axes rotated by  degr.           


mkhdr,header_acube,acube
paraname = ['RA','DEC','scale','PA_SPEC','PA_IMAG','NAXIS_C1','NAXIS_C2','NAXIS_C3','LINE','SIGMALIM','WCSAXES','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2','CD1_1','CD1_2','CD2_1','CD2_2','EQUINOX']
angle = (Pa_spec-90.)/180.*!dpi

angle = [cos(angle),sin(angle),-1.*sin(angle),cos(angle)]*scale
paraval = [ra,dec,scale,pa_spec,pa_imag,naxis_datacube,line,sigmalim,2,crval2,crval3,crpix2,crpix3,cdelt2,cdelt3,angle(0),angle(1),angle(2),angle(3),2000.]
comment = ['RA at spatial [0,0] in mosaic','DEC at spatial [0,0] in mosaic','plate scale in degree per pixel','position angle of spectrograph on sky','position angle of imager on sky','sizes of original mosaic datacube','sizes of original mosaic datacube','sizes of original mosaic datacube','rest frame wavelength of the spectroscopic line fit in this image','minimum significance required for acceptable fit (standard deviations)','Number of axes in WCS system','[deg] R.A. at reference pixel ','[deg] DEC at reference pixel ','Reference pixel location','Reference pixel location','Pixel scale in degree/pixel','Pixel scale in degree/pixel','RA, Dec axes rotated by  degr','RA, Dec axes rotated by  degr','RA, Dec axes rotated by  degr','RA, Dec axes rotated by  degr','J2000 coordinates']
hdr_paraname_str = ['CTYPE1','CTYPE2','CUNIT1','CUNIT2','RADESYS']
hdr_paraname_str_val = ['RA---TAN','DEC--TAN','deg','deg','FK5']
for ii=0,n_elements(hdr_paraname_Str)-1 do sxaddpar,header_acube,hdr_paraname_str(ii),hdr_paraname_str_val(ii)
for ii=0,n_elements(paraname)-1 do sxaddpar, header_acube,paraname[ii],paraval[ii],comment[ii]
;stop

;writefits, 'output/'+name+'_acube.fits', acube,header_acube
;writefits, 'output/'+name+'_aerrorcube.fits', aerrorcube
print, 'finished running fit spec general'
end

; ----------------------------------------------------------------
