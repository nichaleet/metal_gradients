pro meanclipcswa31hb
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.164 arc CSWA 15.

;;; [NII] + Ha
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s141206_a020001_Jbb_100.fits', hdr)  ; "A" position
c2 = readfits('s141207_a024001_Jbb_100.fits')       ; "B" position


  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(1574,67,20,2) ;original size for each frame = 64x19
cube[*,3:66,1:19,0] = c1
cube[*,0:63,0:18,1]  = c2


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
unwrap, smooth(cube_mean,[1,3,3]), 'c31_Hn3_unwrapped_crrej.fits'

writefits,'cswa31_Hb_Jbb_handmosaic_scaledsky_030hr.fits', cube_mean, hdr
 
end
