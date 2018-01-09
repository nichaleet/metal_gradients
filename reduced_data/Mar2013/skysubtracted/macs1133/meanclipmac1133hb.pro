pro meanclipmac1133hb
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of mac1133hb

;;; OIII + Hb
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First exposure
c1 = readfits('s130304_a021001_Jn2_100.fits', hdr)  ; "A" position
c2 = readfits('s130304_a021002_Jn2_100.fits')       ; "B" position
c3 = readfits('s130304_a021003_Jn2_100.fits')  ; A
c4 = readfits('s130304_a021004_Jn2_100.fits')  ; B


  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data cubes are offset by 26 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(408,66,77,4)
cube[*,*,0:50,0] = c1
cube[*,*,26:76,1] = c2
cube[*,*,0:50,0] = c3
cube[*,*,26:76,1] = c4




  ; Remove outliers -- this is an efficient cosmic ray rejection
cube[where(abs(cube) ge 0.05)] = 0.
  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data
  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'm1133_Jn2_unwrapped_crrej.fits'

writefits,'macs1133_Hb_Jn2_handmosaic_sky_1hr.fits', cube_mean, hdr
 
end
