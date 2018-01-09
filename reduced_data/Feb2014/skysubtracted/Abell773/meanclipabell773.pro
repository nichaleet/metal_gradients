pro meanclipabell773
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.3 arc Abell773.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
  c1 = readfits('s130303_a008001_Kc3_100.fits', hdr) ; "A" position
  c2 = readfits('s130303_a008002_Kc3_100.fits')      ; "B" position
  c3 = readfits('s130303_a008003_Kc3_100.fits')      ; "A" position
  c4 = readfits('s130303_a008004_Kc3_100.fits')      ; "B" position

  ;second night
  c5 = readfits('s140222_a006001_Kc3_100.fits')     ; "A" position
  c6 = readfits('s140222_a006002_Kc3_100.fits')     ; "B" position
  c7 = readfits('s140222_a006003_Kc3_100.fits')     ; "A" position
  c8 = readfits('s140222_a006004_Kc3_100.fits')     ; "B" position
  c9 = readfits('s140222_a007001_Kc3_100.fits')     ; "A" position
  c10 = readfits('s140222_a007002_Kc3_100.fits')    ; "B" position

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially. Original size is 66x51
cube = make_array(433,89,80,10)
cube[*,0:65,29:79,0] = c1
cube[*,0:65,5:55,1] = c2
cube[*,0:65,29:79,0] = c3
cube[*,0:65,5:55,1] = c4
cube[*,23:88,24:74,0] = c5
cube[*,23:88,0:50,1] = c6
cube[*,23:88,24:74,0] = c7
cube[*,23:88,0:50,1] = c8
cube[*,23:88,24:74,0] = c9
cube[*,23:88,0:50,1] = c10


  ; Remove outliers -- this is an efficient cosmic ray rejection
cube[where(abs(cube) ge 0.1)] = 0.
  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data
  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'a773_kc3_unwrapped_crrej.fits'

writefits,'abell773_Ha_Kc3_handmosaic_sky_230hr.fits', cube_mean, hdr
 
end
