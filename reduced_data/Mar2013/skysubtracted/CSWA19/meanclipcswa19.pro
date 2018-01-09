pro meanclipcswa19
 
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.03 arc CSWA19Hbeta

;;; [OIII] and Hbeta
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First exposure
c1 = readfits('s130304_a018001_Hn1_100.fits')  ; "A" position
c2 = readfits('s130304_a018002_Hn1_100.fits',hdr)       ; "B" position
c3 = readfits('s130304_a018003_Hn1_100.fits')  ; A
c4 = readfits('s130304_a018004_Hn1_100.fits')  ; B

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data from both nights are almost exactly aligned (within ~0.1")!
; Data cubes are offset by 24 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(376,66,71,4)
cube[*,*,20:70,0] = c1
cube[*,*,0:50,0] = c2
cube[*,*,20:70,1] = c3
cube[*,*,0:50,0] = c4



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
unwrap, smooth(cube_mean,[1,3,3]), 'cswa19_Hn1_unwrapped_crrej.fits'

writefits,'cswa19_Hb_Hn1_handmosaic_sky_1hr.fits', cube_mean, hdr
 
end
