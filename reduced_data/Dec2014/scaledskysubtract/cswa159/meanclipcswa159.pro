pro meanclipcswa159
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.30 arc CSWA 159.

;;;[OIII], Hb
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s141207_a009001_tlc_Hn3_100.fits', hdr)  ; "A" position
c2 = readfits('s141207_a009002_tlc_Hn3_100.fits')       ; "B" position
c3 = readfits('s141207_a009003_tlc_Hn3_100.fits')  ; A
c4 = readfits('s141207_a009004_tlc_Hn3_100.fits')  ; B

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data from all the frames are not exactly aligned. The 13th frame is 1 pixels to the right (in x direction) and 1 pixel above (in y direction) of the other frames

; Data cubes are offset by 24 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(411,66,75,4)
cube[*,0:65,24:74,0] = c1
cube[*,0:65,0:50,1]  = c2
cube[*,0:65,24:74,2] = c3
cube[*,0:65,0:50,3] = c4


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
unwrap, smooth(cube_mean,[1,3,3]), 'c31_tlc_Hn3_unwrapped_crrej.fits'

writefits,'cswa159_Hb_tlc_Hn3_handmosaic_scaledsky_1hr.fits', cube_mean, hdr
writefits,'/scr2/nichal/workspace/reduced_data/mosaic/cswa159_Hb_tlc_Hn3_handmosaic_scaledsky_1hr.fits', cube_mean, hdr
 
end
