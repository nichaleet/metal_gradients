pro run_fitspec_general,input=input,filter=filter,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,run_fitspec_int_index,secondline,smooth,imgymins,imgymaxs

; To run fitspec_general (and) fitspec_intensity

;INPUT:
; LINE          - rest frame wavelength of the spectroscopic line to be fit, in microns (e.g. 0.6563 for Halpha).
; Z             - approximate redshift of LINE.
; SIGMALIM      - minimum significance required for acceptable fit (standard deviations).
; XMINS,YMAXS - x,y min/max pixel values of a blank sky region, used to measure the sky spectrum. This determines the Gaussian weights for fitting each spaxel. The sky region used is FINALCUBE[*,XMINS:XMAXS,YMINS:YMAXS]. y ranges from 0 to 33
; fwhm_gaussian  - fwhm to smooth data with filter_image (#pixels to smooth data with)
; run_fit_spec_int_index - equal to 1 if want to run
;                          fitspec_intensity.pro 
; secondline     - rest frame wavelength of the spectroscopic line o  metal element (to be ratioed with the first line)
                 

;Example: 
;file ='/scr2/nichal/osiris/drs/output/scaledskysubtracted/s130911_a023_mosaic_Hbb_100.fits'
;run_fitspec_general,input=file,filter='Hbb',.4861,2.3,5.,29,42,7,11,1,0.5006,2.,11,22

;Note: xymins xymaxs, the x and y in the program are actually
;opposite to  when viewed in datacube. y is frame# in datacube. x is y
;dimension in ds9

finalcube = readfits(input,header)

if filter eq 'Jbb' then wl =1.180+0.00015*findgen(1574)else $
if filter eq 'Jn2' then wl =1.228+0.00015*findgen(408)else $
if filter eq 'Hbb' then wl =1.473+0.0002*findgen(1651)else $
if filter eq 'Hn1' then wl =1.466+0.0002*findgen(376)else $
if filter eq 'Hn2' then wl =1.532+0.0002*findgen(391)else $
if filter eq 'Hn3' then wl =1.594+0.0002*findgen(411)else $
if filter eq 'Hn4' then wl =1.652+0.0002*findgen(426)else $
if filter eq 'Kn1' then wl =1.955+0.00025*findgen(401)else $
if filter eq 'Kn2' then wl =2.036+0.00025*findgen(421)else $
if filter eq 'Kc3' then wl =2.121+0.00025*findgen(433)else begin 
   print, 'no filter match'   
   stop
   stop
endelse


;Clean finalcube (set the outliers to 0.)
clean_treshold_op = 10.*stdev(finalcube[where(abs(finalcube) le 1.)]) ;10sigma
clean_treshold = 1.
print, 'Now limit for outliers is ', clean_treshold
print, 'want to change to ''clean_treshold = 10.*stdev(finalcube[where(abs(finalcube) le 1.)]) ;10sigma''(',clean_treshold_op,') Input .cont for no. yes then copy and paste.'
stop
bad_pixels = where(abs(finalcube) gt clean_treshold)
finalcube[bad_pixels] = 0.

;smooth the data
if smooth gt 1 then fwhm = smooth-1. else fwhm = 1.
for i=0,n_elements(finalcube[*,0,0])-1 do finalcube[i,*,*] = filter_image(reform(finalcube[i,*,*]),smooth=smooth, fwhm_gaussian=fwhm)

name_pos = strpos(input,'/',/reverse_search)
name = strmid(input, name_pos+1)
name = strmid(name, 0,strlen(name)-5)  ;delete '.fits'


fitspec_general,finalcube,wl,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,0,1,0,500.,500.,name=name,imgymins,imgymaxs,header

if run_fitspec_int_index eq 1 then run_fitspec_intensity,input=input,filter=filter,secondline,xmins,xmaxs,ymins,ymaxs,smooth,header,imgymins,imgymaxs

end

