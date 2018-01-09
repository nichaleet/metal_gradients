pro run_fitspec_general,input=input,filter=filter,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,run_fitspec_int_index,secondline,fwhm_gaussian,xmin,xmax,ymin,ymax,crpix2,crpix3

Common cube, acube,aerrorcube

; To run fitspec_general (and) fitspec_intensity

;INPUT:
; LINE          - rest frame wavelength of the spectroscopic line to be fit, in microns (e.g. 0.6563 for Halpha).
; Z             - approximate redshift of LINE.
; SIGMALIM      - minimum significance required for acceptable fit (standard deviations).
; XMINS,YMAXS - x,y min/max pixel values of a blank sky region, used to measure the sky spectrum. This determines the Gaussian weights for fitting each spaxel. The sky region used is FINALCUBE[*,XMINS:XMAXS,YMINS:YMAXS]. y ranges from 0 to 33
; run_fit_spec_int_index - equal to 1 if want to run
;                          fitspec_intensity.pro 
; secondline     - rest frame wavelength of the spectroscopic line o  metal element (to be ratioed with the first line)
; fwhm_gaussian  - fwhm to smooth data with filter_image (#pixels to smooth data with)
; xmin, xmax, ymin,ymax - frames for the arc
; crpix2, crpix3 - the pixels in KECK image of CRVALS that match the lensing image e.g. from HST. If don't know then put 32,19

;OUTPUT: 
;x-y image of 9 things: 
;0) velocity according to the firstling line in comparing to the input redshift
;1) amplitude of the first line
;2) velocity dispersion
;3) number of bins
;4) first line fit significance
;5) redshift
;6) amplitude of the second line
;7) second line fit significance                  
;8) 1sigma error in velocity
;9) 1 sigma error in the first line amplitude
;10)1 sigma error in velocity dispersion
;11) 1 sigma error in the second line amplitude

;Example:
;file = "/scr2/nichal/workspace/reduced_data/mosaic/cswa19_Ha_Kn1_mosaic_scaledsky_4hr.fits"
;run_fitspec_general,input=file, filter ='Kn1',.6562801,2.0332,5.,30,40,18,23,1,0.658345,3.,0,60,15,40,26,44
 
;Note: xymins xymaxs, the x and y in the program are actually
;opposite to  when viewed in datacube. y is frame# in datacube. x is y
;dimension in ds9

incube = readfits(input,header)

if filter eq 'Jbb' then wlarr =1.180+0.00015*findgen(1574)else $
if filter eq 'Jn1' then wlarr =1.174+0.00014987*findgen(388)else $
if filter eq 'Jn2' then wlarr =1.228+0.00015*findgen(408)else $
if filter eq 'Hbb' then wlarr =1.473+0.0002*findgen(1651)else $
if filter eq 'Hn1' then wlarr =1.466+0.0002*findgen(376)else $
if filter eq 'Hn2' then wlarr =1.532+0.0002*findgen(391)else $
if filter eq 'Hn3' then wlarr =1.594+0.0002*findgen(411)else $
if filter eq 'Hn4' then wlarr =1.652+0.0002*findgen(426)else $
if filter eq 'Kn1' then wlarr =1.955+0.00025*findgen(401)else $
if filter eq 'Kn2' then wlarr =2.036+0.00025*findgen(421)else $
if filter eq 'Kc5' then wlarr =2.292+0.00025*findgen(465)else $
if filter eq 'Kc3' then wlarr =2.121+0.00025*findgen(433)else begin 
   print, 'no filter match'   
   stop
endelse

name_pos = strpos(input,'/',/reverse_search)
name = strmid(input, name_pos+1)
name = strmid(name, 0,strlen(name)-5)  ;delete '.fits'

;Clean finalcube (set the outliers to 0.)
clean_treshold = 1.
bad_pixels = where(abs(incube) gt clean_treshold)
if bad_pixels(0) ne -1 then begin
   incube[bad_pixels] = 0.
   ind = array_indices(incube,bad_pixels)
endif

;Below are for fixing the glitches that go right across the galaxy

if name eq 'macs0717_Hb_Hbb_mosaic_scaledsky_1hr' then begin
   for wi=1518,1519 do incube(wi,42,20)=0.5*(incube(wi,41,20)+incube(wi,43,20))
   for wi=1524,1526 do incube(wi,45,20)=0.5*(incube(wi,44,20)+incube(wi,46,20))
   for wi=1520,1524 do incube(wi,43,19)=0.5*(incube(wi,42,19)+incube(wi,44,19))
   for wi=1522,1524 do incube(wi,44,17)=0.5*(incube(wi,43,17)+incube(wi,45,17))
   ;incube(1521,43,18) = 0.017
endif

;smooth the data
if fwhm_gaussian ne 0 then for i=0,n_elements(incube[*,0,0])-1 do incube[i,*,*] = filter_image(reform(incube[i,*,*]),fwhm_gaussian=fwhm_gaussian)
writefits,'checkcube.fits',incube,header

finalcube = incube
wl=wlarr
if name eq 'cswa31_Hb_tlc_Jbb_handmosaic_scaledsky_030hr' then nmax = 5 else nmax = 0

fitspec_general,finalcube,wl,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,0,1,nmax,350.,200.,name=name,xmin,xmax,ymin,ymax,header,header_acube,crpix2,crpix3
;fitspec_general, finalcube, wl, line, z, sigmalim, xmins, xmaxs, ymins, ymaxs, nmin, ni, nmax, vmax, wmax, name=name,xmin,xmax,ymin,ymax,header,header_acube,crpix2,crpix3

if run_fitspec_int_index eq 1 then begin
   z_cube    = acube[*,*,5]
   disp_cube = (acube[*,*,2] / 2.998d5) * z_cube * line ; convert velocity dispersion from km/s to microns
   bin_cube  = acube[*,*,3]
   wl=wlarr
   fitspec_intensity, incube, wl,secondline,z_cube,disp_cube,bin_cube,xmins,xmaxs,ymins,ymaxs,name=name,header,xmin,xmax,ymin,ymax
endif

output_cube = [[[acube]],[[aerrorcube]]]
acubefile='/scr2/nichal/workspace/output/'+name+'_acube.fits'
writefits, acubefile, output_cube,header_acube

;make point spread function center on the center of the field. This will be used in the model of the kinematics
sizecube=size(output_cube)
sizex=sizecube(1)
sizey=sizecube(2)
sizepsf=min(sizex,sizey)
fwhm=sizepsf/3.
array = PSF_GAUSSIAN( Npixel=sizepsf, FWHM=[fwhm,fwhm], /NORMAL )
arraypsf=extend_array(array,sizex,sizey)
writefits,'/scr2/nichal/workspace/output/psf_files/'+name+'_psf.fits',arraypsf,header_acube


;Analyze 2nd line detection
seclineflux = output_cube[*,*,6]
uncertainty  = output_cube[*,*,11]
goodpix = where(finite(seclineflux))
negative_pix = where(seclineflux lt 0.)
total_2ndline = total(seclineFlux(goodpix))
print, 'total of 2nd line flux: ',total_2ndline
negative_upperlim = where(seclineflux+uncertainty lt 0.)
plothist, seclineflux(goodpix),title='NII flux',/autobin
print, 'Out of', n_Elements(goodpix),'pixels,',n_elements(negative_pix),'are negative.', n_Elements(negative_upperlim),' have negative upper limit.'
end

