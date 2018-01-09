pro meanclipcswa159
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.3 arc CSWA 159.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s130911_a027001_Kc3_100.fits')  ; "A" position
c2 = readfits('s130911_a027002_Kc3_100.fits')       ; "B" position
  ; Second night
  ; Sep's ref position is 4 pixels to the right of Dec's ref position
c3 = readfits('s141207_a007002_Kc3_100.fits',hdr) ; "A" position
c4 = readfits('s141207_a007001_Kc3_100.fits') ; "B" position 
c5 = readfits('s141207_a007004_Kc3_100.fits') ; "A" position
c6 = readfits('s141207_a007003_Kc3_100.fits') ; "B" position 
c7 = readfits('s141207_a008002_Kc3_100.fits') ; "A" position
c8 = readfits('s141207_a008001_Kc3_100.fits') ; "B" position 

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data cubes are offset by 24 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
; Sep's ref position is 4 pixels to the right of Dec's ref position
cube = make_array(433,70,75,8)
cube[*,4:69,24:74,0] = c1
cube[*,4:69,0:50,1] = c2
cube[*,0:65,24:74,2] = c3
cube[*,0:65,0:50,3] = c4
cube[*,0:65,24:74,4] = c5
cube[*,0:65,0:50,5] = c6
cube[*,0:65,24:74,6] = c7
cube[*,0:65,0:50,7] = c8

  ; Remove outliers -- this is an efficient cosmic ray rejection
cube[where(abs(cube) ge 0.2)] = 0.
  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data
  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'c159_kc3_unwrapped_crrej.fits'

writefits,'cswa159_Ha_Kc3_handmosaic_sky_2hr.fits', cube_mean, hdr
 
end
