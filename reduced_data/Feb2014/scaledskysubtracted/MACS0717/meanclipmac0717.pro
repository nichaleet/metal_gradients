pro meanclipmac0717
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.55 arc MACS0717.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s140224_a005001_Kc5_100.fits')  ; "A" position
c2 = readfits('s140224_a005002_Kc5_100.fits',hdr)       ; "B" position
c3 = readfits('s140224_a005003_Kc5_100.fits')  ; A
c4 = readfits('s140224_a005004_Kc5_100.fits')  ; B
c5 = readfits('s140224_a006001_Kc5_100.fits')  ; A
c6 = readfits('s140224_a006002_Kc5_100.fits')  ; B

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))

; Data cubes are offset by 20 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(465,66,71,6)
cube[*,*,20:70,0] = c1
cube[*,*,0:50,1] = c2
cube[*,*,20:70,2] = c3
cube[*,*,0:50,3] = c4
cube[*,*,20:70,4] = c5
cube[*,*,0:50,5] = c6

  ; Remove outliers -- this is an efficient cosmic ray rejection
cube[where(abs(cube) ge 0.08)] = 0.
  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data
  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'macs0717_kc5_unwrapped_crrej.fits'

writefits,'macs0717_Ha_Kc5_handmosaic_scaledsky_130hr.fits', cube_mean, hdr
 
end
