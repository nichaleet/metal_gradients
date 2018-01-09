pro meanclipmacs1133mar
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of Macs1133Feb2014 data.

;;; Ha + [NII]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s130303_a016001_Hn4_100.fits', hdr)  ; "A" position
c2 = readfits('s130303_a016002_Hn4_100.fits')       ; "B" position
c3 = readfits('s130303_a016003_Hn4_100.fits')  ; A
c4 = readfits('s130303_a016004_Hn4_100.fits')  ; B
c5 = readfits('s130303_a017001_Hn4_100.fits')  ; A
c6 = readfits('s130303_a017002_Hn4_100.fits')  ; B
  ; Second night
c7  = readfits('s130304_a022001_Hn4_100.fits')  ; A
c8  = readfits('s130304_a022002_Hn4_100.fits')  ; B
c9  = readfits('s130304_a022003_Hn4_100.fits')  ; A
c10 = readfits('s130304_a022004_Hn4_100.fits')  ; B

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data cubes are offset by 26 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
  ; original size is 66x51 pixels
cube = make_array(426,68,77,10)

cube[*,2:67,0:50,0] = c1
cube[*,2:67,26:76,1] = c2
cube[*,2:67,0:50,2] = c3
cube[*,2:67,26:76,3] = c4
cube[*,2:67,0:50,4] = c5
cube[*,2:67,26:76,5] = c6
cube[*,0:65,0:50,6] = c7
cube[*,0:65,26:76,7] = c8
cube[*,0:65,0:50,8] = c9
cube[*,0:65,26:76,9] = c10

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
unwrap, smooth(cube_mean,[1,3,3]), 'm1133Mar_Hn4_unwrapped_crrej.fits'

writefits,'macs1133Mar_Ha_Hn4_handmosaic_scaledsky_230hr.fits', cube_mean, hdr
 
end
