pro meanclipcswa20hb
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of the arc CSWA 20Hb.

;;; [OIII] + Hbeta
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First exposure
c1 = readfits('s140223_a026001_tlc_Jn1_100.fits')  ; "A" position
c2 = readfits('s140223_a026002_tlc_Jn1_100.fits',hdr)       ; "B" position
c3 = readfits('s140223_a026003_tlc_Jn1_100.fits')  ; A

  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))

; The original cube size is 66*51
; Data cubes are offset by 20 pixels along the short axis (axis 2)
; Combine all frames into a mega-array
cube = make_array(388,66,71,3)
cube[*,0:65,20:70,0] = c1
cube[*,0:65,0:50,1] = c2
cube[*,0:65,20:70,2] = c3

  ; Remove outliers -- this is an efficient cosmic ray rejection
cube[where(abs(cube) ge 0.1)] = 0.
                                ; Remove the traces that are clear to be cosmic rays that lies on galaxy image but have the values less than the threshold
cube[41:46,53,28,0]=(cube[41:46,52,28,0]+cube[41:46,54,28,0]+cube[41:46,53,27,0]+cube[41:46,53,29,0])/4.
cube[45:50,30,28,0]=(cube[45:50,29,28,0]+cube[45:50,31,28,0]+cube[45:50,30,29,0]+cube[45:50,30,27,0])/4.
;print,cube[41:46,53,28,0]
;print,cube[45:50,30,28,0]
;stop

  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data

;There is still a streak of cosmic ray line across the OIII emission line.
cube_mean[290:300,22,18] = (cube_mean[290:300,23,18]+cube_mean[290:300,21,18]+cube_mean[290:300,22,17]+cube_mean[290:300,22,19])/4.
;cube_mean[41:46,53,28]=0.
;cube_mean[45:50,30,28]=0.

  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'c20_tlc_Jn1_unwrapped_crrej.fits'

writefits,'cswa20_Hb_tlc_Jn1_handmosaic_sky_075hr.fits', cube_mean, hdr
 
end
