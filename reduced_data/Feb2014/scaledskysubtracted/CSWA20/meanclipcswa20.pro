pro meanclipcswa20
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the arc CSWA 20.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First exposure
c1 = readfits('s140224_a021001_Hn2_100.fits')  ; "A" position
c2 = readfits('s140224_a021002_Hn2_100.fits',hdr)       ; "B" position
 ; Second exposure
c3 = readfits('s140224_a024001_Hn2_100.fits')  ; A
c4 = readfits('s140224_a024002_Hn2_100.fits')  ; B

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))

; The two exposure offsets about ~0.14" from each other
; The original cube size is 66*51
; Data cubes are offset by 20 pixels along the short axis (axis 2)
                                ; Combine all frames into a mega-array
                                ; Shift each frame to align them spatially (see configuration in Nicha's logbook page 61)
cube = make_array(391,67,72,4)
cube[*,0:65,21:71,0] = c1
cube[*,0:65,1:51,1] = c2
cube[*,1:66,20:70,2] = c3
cube[*,1:66,0:50,3] = c4

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
unwrap, smooth(cube_mean,[1,3,3]), 'c20_kn2_unwrapped_crrej.fits'

writefits,'cswa20_Ha_Hn2_handmosaic_scaledsky_1hr.fits', cube_mean, hdr
 
end
