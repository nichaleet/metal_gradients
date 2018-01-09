pro meanclipcswa31
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.164 arc CSWA 15.

;;; [NII] + Ha
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s141206_a013001_tlc_Hn3_100.fits', hdr)  ; "A" position
c2 = readfits('s141206_a013002_tlc_Hn3_100.fits')       ; "B" position
c3 = readfits('s141206_a013003_tlc_Hn3_100.fits')  ; A

c4 = readfits('s141206_a015001_tlc_Hn3_100.fits')  ; B
c5 = readfits('s141206_a015002_tlc_Hn3_100.fits')  ; A

c6 = readfits('s141206_a016001_tlc_Hn3_100.fits')  ; B
c7 = readfits('s141206_a016002_tlc_Hn3_100.fits')  ; A
c8 = readfits('s141206_a016003_tlc_Hn3_100.fits')  ; B
c9 = readfits('s141206_a016004_tlc_Hn3_100.fits')  ; A

c10 = readfits('s141206_a017001_tlc_Hn3_100.fits')  ; B
c11 = readfits('s141206_a017002_tlc_Hn3_100.fits')  ; A
c12 = readfits('s141206_a017003_tlc_Hn3_100.fits')  ; B

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data from all the frames are not exactly aligned. The 13th frame is 1 pixels to the right (in x direction) and 1 pixel above (in y direction) of the other frames

; Data cubes are offset by 20 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(411,67,72,12)
cube[*,1:66,21:71,0] = c1
cube[*,1:66,1:51,1]  = c2
cube[*,1:66,21:71,2] = c3

cube[*,0:65,20:70,3] = c4
cube[*,0:65,0:50,4]  = c5
cube[*,0:65,20:70,5] = c6
cube[*,0:65,0:50,6]  = c7
cube[*,0:65,20:70,7] = c8
cube[*,0:65,0:50,8]  = c9
cube[*,0:65,20:70,9] = c10
cube[*,0:65,0:50,10]  = c11
cube[*,0:65,20:70,11] = c12

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
unwrap, smooth(cube_mean,[1,3,3]), 'c31_tlc_Hn3_unwrapped_crrej.fits'

writefits,'cswa31_Ha_tlc_Hn3_handmosaic_scaledsky_3hr.fits', cube_mean, hdr
 
end
