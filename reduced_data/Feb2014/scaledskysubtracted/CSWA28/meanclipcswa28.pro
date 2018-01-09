pro meanclipcswa28
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.164 arc CSWA 15.

;;; [NII] + Ha
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s130303_a022001_Kn1_100.fits', hdr)  ; "A" position
c2 = readfits('s130303_a022002_Kn1_100.fits')       ; "B" position

  ; Second night
c3 = readfits('s140223_a023001_Kn1_100.fits')  ; A
c4 = readfits('s140223_a023002_Kn1_100.fits')  ; B
c5 = readfits('s140223_a023003_Kn1_100.fits')  ; A
c6 = readfits('s140223_a023004_Kn1_100.fits')  ; B
  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data from both nights are not exactly aligned. The second night is 3 pixels to the right (in x direction) and 1 pixel above (in y direction) of the first night
; Data cubes are offset by 20 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(401,69,72,10)
cube[*,0:65,20:70,0] = c1
cube[*,0:65,0:50,1] = c2
cube[*,3:68,21:71,2] = c3
cube[*,3:68,1:51,3] = c4
cube[*,3:68,21:71,4] = c5
cube[*,3:68,1:51,5] = c6

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
unwrap, smooth(cube_mean,[1,3,3]), 'c28_kn1_unwrapped_crrej.fits'

writefits,'cswa28_Ha_Kn1_handmosaic_sky_130hr.fits', cube_mean, hdr
 
end
