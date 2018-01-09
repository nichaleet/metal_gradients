pro cswa11run_fitspec_general,input=input,filter=filter,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,run_fitspec_int_index,secondline,fwhm_gaussian,xmin,xmax,ymin,ymax,crpix2,crpix3

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
;[0)velocity centroid (km/s) of the first galaxy 
;[1]amplitude (same units as FINALCUBE) of the first galaxy
;[2]velocity dispersion (km/s)of the first galaxy
;[3]velocity centroid (km/s) of the 2nd galaxy 
;[4]amplitude of the 2nd galaxy (same units as FINALCUBE)
;[5]velocity dispersion (km/s)of the 2nd line
;[6]binning factor
;[7]fit significance (# of standard deviations)
;[8]redshift of the 1st galaxy
;[9]redshift of the 2nd galaxy
;[10] amplitude of the NII(or Hb) of the first galaxy
;[11] amplitude of the NII(or Hb) of the second galaxy
;[12] NII line fit significance
;[13] 1 sigma error in velocity of the first galaxy
;[14] 1 sigma error in Ha or OIII of the first galaxy
;[15] 1 sigma error in velocity dispersion of the first galaxy           ;[16] 1 sigma error in velocity of the second galaxy
;[17] 1 sigma error in Ha or OIII of the second galaxy
;[18] 1 sigma error in velocity dispersion of the second galaxy          ;[19] 1 sigma error in the NII line amplitude of the first galaxy
;[20] 1 sigma error in the NII line amplitude of the first galaxy

;Example:
;cswa11run_fitspec_general,input=file,filter ='Hn2',.6564,1.409,5.,50,60,40,50,1,0.6585,3.,15,55,30,55,32,23

 
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
clean_treshold = .1
bad_pixels = where(abs(incube) gt clean_treshold)
if bad_pixels(0) ne -1 then begin
   incube[bad_pixels] = 0.
   ind = array_indices(incube,bad_pixels)
endif

;smooth the data
if fwhm_gaussian ne 0 then for i=0,n_elements(incube[*,0,0])-1 do incube[i,*,*] = filter_image(reform(incube[i,*,*]),fwhm_gaussian=fwhm_gaussian)
writefits,'checkcube.fits',incube,header

finalcube = incube
wl=wlarr
nmax = 0

cswa11fitspec_general,finalcube,wl,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,0,1,nmax,500.,500.,name=name,xmin,xmax,ymin,ymax,header,header_acube,crpix2,crpix3

;fitspec_general, finalcube, wl, line, z, sigmalim, xmins, xmaxs, ymins, ymaxs, nmin, ni, nmax, vmax, wmax, name=name,xmin,xmax,ymin,ymax,header,header_acube,crpix2,crpix3

;Now 
;acube has 10 frames: [1)velocity centroid (km/s), 2)amplitude (same units as FINALCUBE), 3)velocity dispersion (km/s),4)velocity centroid (km/s) of the 2nd line, 5)amplitude of the 2nd line (same units as FINALCUBE), 6)velocity dispersion (km/s)of the 2nd line, 7)binning factor, 8)fit significance (# of standard deviations), 9)redshift of the 1st line, 10)redshift of the 2nd line].
;aerrorcube has 6 frames: [1] 1sigma error in velocity [2] 1 sigma error in the first line amplitude [3] 1 sigma error in velocity dispersion [4,5,6] the same as [1,2,3] but of the second galaxy

if run_fitspec_int_index eq 1 then begin
   z_cube1    = acube[*,*,8]
   disp_cube1 = (acube[*,*,2] / 2.998d5) * z_cube1 * line ; convert velocity dispersion from km/s to microns
   z_cube2    = acube[*,*,9]
   disp_cube2 = (acube[*,*,5] / 2.998d5) * z_cube2 * line ; convert velocity dispersion from km/s to microns
   bin_cube  = acube[*,*,6]
   wl=wlarr
   cswa11fitspec_intensity, incube, wl,secondline,z_cube1,z_cube2,disp_cube1,disp_cube2,bin_cube,xmins,xmaxs,ymins,ymaxs,name=name,header,xmin,xmax,ymin,ymax
endif

output_cube = [[[acube]],[[aerrorcube]]]
writefits, '/scr2/nichal/workspace/output/'+name+'_acube.fits', output_cube,header_acube

;Separate the output into 2 galaxies like the output of run_fitspec_General.pro
output_cube1 = [[[output_cube[*,*,0]]],[[output_cube[*,*,1]]],[[output_cube[*,*,2]]],[[output_cube[*,*,6]]],[[output_cube[*,*,7]]],[[output_cube[*,*,8]]],[[output_cube[*,*,10]]],[[output_cube[*,*,12]]],[[output_cube[*,*,13]]],[[output_cube[*,*,14]]],[[output_cube[*,*,15]]],[[output_cube[*,*,19]]]]
output_cube1[*,0:39,*] = 1./0.
output_cube1[*,55:82,*]= 1./0. 

output_cube2 = [[[output_cube[*,*,3]]],[[output_cube[*,*,4]]],[[output_cube[*,*,5]]],[[output_cube[*,*,6]]],[[output_cube[*,*,7]]],[[output_cube[*,*,9]]],[[output_cube[*,*,11]]],[[output_cube[*,*,12]]],[[output_cube[*,*,16]]],[[output_cube[*,*,17]]],[[output_cube[*,*,18]]],[[output_cube[*,*,20]]]]

output_cube2[55:67,*,*] = 1./0.
output_cube2[0:15,*,*] = 1./0.
output_cube2[*,55:82,*]= 1./0.

writefits, '/scr2/nichal/workspace/output/'+name+'_firstgal_acube.fits', output_cube1,header_acube
writefits, '/scr2/nichal/workspace/output/'+name+'_secondgal_acube.fits', output_cube2,header_acube

stop
;Analyze 2nd line detection
seclineflux = output_cube[*,*,6]
uncertainty  = output_cube[*,*,11]
goodpix = where(finite(seclineflux))
negative_pix = where(seclineflux lt 0.)
total_2ndline = total(seclineFlux(goodpix))
print, 'total of 2nd line flux: ',total_2ndline
negative_upperlim = where(seclineflux+uncertainty lt 0.)
;plothist, seclineflux(goodpix),title='NII flux',/autobin
print, 'Out of', n_Elements(goodpix),'pixels,',n_elements(negative_pix),'are negative.', n_Elements(negative_upperlim),' have negative upper limit.'
end

