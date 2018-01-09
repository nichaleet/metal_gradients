pro meanclipabell773
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.3 arc Abell773

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First exposure
c1 = readfits('s130303_a006001_tlc_Hn3_100.fits', hdr)  ; "A" position
c2 = readfits('s130303_a006002_tlc_Hn3_100.fits')       ; "B" position
c3 = readfits('s130303_a006003_tlc_Hn3_100.fits')  ; A
c4 = readfits('s130303_a006004_tlc_Hn3_100.fits')  ; B
c5 = readfits('s130303_a007001_tlc_Hn3_100.fits')  ; A
c6 = readfits('s130303_a007002_tlc_Hn3_100.fits')  ; B
c7 = readfits('s130303_a007003_tlc_Hn3_100.fits')  ; A 
c8 = readfits('s130303_a007004_tlc_Hn3_100.fits')  ; B

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data from both nights are almost exactly aligned (within ~0.1")!
; Data cubes are offset by 24 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(411,66,75,8)
cube[*,*,0:50,0] = c1
cube[*,*,24:74,1] = c2
cube[*,*,0:50,2] = c3
cube[*,*,24:74,3] = c4
cube[*,*,0:50,4] = c5
cube[*,*,24:74,5] = c6
cube[*,*,0:50,6] = c7
cube[*,*,24:74,7] = c8



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
unwrap, smooth(cube_mean,[1,3,3]), 'a773_tlc_Hn3_unwrapped_crrej.fits'

writefits,'abell773_Hb_tlc_Hn3_handmosaic_sky_2hr.fits', cube_mean, hdr
 
end
