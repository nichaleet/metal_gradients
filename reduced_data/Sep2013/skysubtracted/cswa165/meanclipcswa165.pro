pro meanclipcswa165
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of arc CSWA 165.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s130911_a029001_Kn2_100.fits', hdr)  ; "A" position
c2 = readfits('s130911_a029002_Kn2_100.fits')       ; "B" position
c3 = readfits('s130911_a029003_Kn2_100.fits')  ; A
c4 = readfits('s130911_a029004_Kn2_100.fits')  ; B


  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data cubes are offset by 24 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(421,66,75,4)
cube[*,*,24:74,0] = c1
cube[*,*,0:50,1] = c2
cube[*,*,24:74,2] = c3
cube[*,*,0:50,3] = c4

  ; Remove outliers -- this is an efficient cosmic ray rejection
cube[where(abs(cube) ge 0.5)] = 0.
  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data
  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'c165_kn2_unwrapped_crrej.fits'

writefits,'cswa165_Ha_Kn2_handmosaic_sky_1hr.fits', cube_mean, hdr
 
end
