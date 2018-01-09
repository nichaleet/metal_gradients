pro meanclipcswa139
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the z=2.54 arc CSWA 139.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First run
c1 = readfits('cswa139_Hb_tlc_Hbb_pipelinemosaic_scaledsky_feb_130hr.fits',hdr)  ; "A" position
 
  ; Second run
c2 = readfits('cswa139_Hb_tlc_Hbb_pipelinemosaic_scaledsky_dec_130hr.fits') ;B position

;Each data cube is of the size 465x66x71
  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data cubes are offset by 20 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(1651,100,45,2)
cube[*,2:85,0:34,0] = c1
cube[*,0:97,8:44,1] = c2


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
unwrap, smooth(cube_mean,[1,3,3]), 'c139_Hbb_tlc_unwrapped_crrej.fits'

writefits,'cswa139_Hb_tlc_Hbb_handmosaic_scaledsky_3hr.fits', cube_mean, hdr
 
end
