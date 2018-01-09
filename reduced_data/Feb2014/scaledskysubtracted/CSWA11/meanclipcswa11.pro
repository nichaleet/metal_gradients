pro meanclipcswa11
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the arc CSWA 20.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s130303_a010001_Hn2_100.fits',hdr)  ; "A" position  
c2 = readfits('s130303_a010002_Hn2_100.fits')       ; "B" position  
c3 = readfits('s130303_a010003_Hn2_100.fits')  ; A  
c4 = readfits('s130303_a010004_Hn2_100.fits')  ; B
  ; second night
c5 = readfits('s130304_a009001_Hn2_100.fits')  ; A  
c6 = readfits('s130304_a009002_Hn2_100.fits')  ; B  
c7 = readfits('s130304_a009003_Hn2_100.fits')  ; A  
c8 = readfits('s130304_a009004_Hn2_100.fits')  ; B  
c9 = readfits('s130304_a010001_Hn2_100.fits')  ; A  
c10 = readfits('s130304_a010002_Hn2_100.fits')  ; B
c11 = readfits('s130304_a010003_Hn2_100.fits')  ; A
c12 = readfits('s130304_a010004_Hn2_100.fits')  ; B
  ; third night
c13 = readfits('s140223_a006001_Hn2_100.fits')  ; B
c14 = readfits('s140223_a006002_Hn2_100.fits')  ; A
c15 = readfits('s140223_a006003_Hn2_100.fits')  ; B
c16 = readfits('s140223_a006004_Hn2_100.fits')  ; A
c17 = readfits('s140223_a007001_Hn2_100.fits')  ; B
c18 = readfits('s140223_a007002_Hn2_100.fits')  ; A

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))

; The two exposure offsets about ~0.14" from each other
; The original cube size is 66*51
; Data cubes are offset by 30 or 24 pixels along the short axis (axis 2)
                                ; Combine all frames into a mega-array
                                ; Shift each frame to align them spatially (see configuration in Nicha's logbook page 60)
cube = make_array(391,68,83,18)
cube[*,2:67,0:50,0] = c1
cube[*,2:67,30:80,1] = c2
cube[*,2:67,0:50,2] = c3
cube[*,2:67,30:80,3] = c4

cube[*,0:65,8:58,4] = c5
cube[*,0:65,32:82,5] = c6
cube[*,0:65,8:58,6] = c7
cube[*,0:65,32:82,7] = c8
cube[*,0:65,8:58,8] = c9
cube[*,0:65,32:82,9] = c10
cube[*,0:65,8:58,10] = c11
cube[*,0:65,32:82,11] = c12

cube[*,2:67,30:80,12] = c13
cube[*,2:67,0:50,13] = c14
cube[*,2:67,30:80,14] = c15
cube[*,2:67,0:50,15] = c16
cube[*,2:67,30:80,16] = c17
cube[*,2:67,0:50,17] = c18



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
unwrap, smooth(cube_mean,[1,3,3]), 'c11_Hn2_unwrapped_crrej.fits'

writefits,'cswa11_Ha_Hn2_handmosaic_sky_430hr.fits', cube_mean, hdr
 
end
