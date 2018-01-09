pro meanclipcswa165
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate 

; IDL script to examine OSIRIS data of arc CSWA 165.

;;; H alpha + [N II]
; Add individual reduced data cubes in such a way that cosmic rays can be rejected from each frame.
  ; First night
c1 = readfits('s130911_a029001_tlc_Kn2_100.fits', hdr)  ; "A" position
c2 = readfits('s130911_a029002_tlc_Kn2_100.fits')       ; "B" position
c3 = readfits('s130911_a029003_tlc_Kn2_100.fits')  ; A
c4 = readfits('s130911_a029004_tlc_Kn2_100.fits')  ; B
  ; Second night
;c5=readfits('s141207_a011001_tlc_Kn2_100.fits') ;A
;c6=readfits('s141207_a011002_tlc_Kn2_100.fits') ;B
;c7=readfits('s141207_a011003_tlc_Kn2_100.fits') ;A
;c8=readfits('s141207_a011004_tlc_Kn2_100.fits') ;B

c5=readfits('s141207_a011001_Kn2_100.fits') ;A
c6=readfits('s141207_a011002_Kn2_100.fits') ;B
c7=readfits('s141207_a011003_Kn2_100.fits') ;A
c8=readfits('s141207_a011004_Kn2_100.fits') ;B
  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data cubes are offset by 24 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(421,68,77,8)
cube[*,2:67,26:76,0] = c1
cube[*,2:67,2:52,1]  = c2
cube[*,2:67,26:76,2] = c3
cube[*,2:67,2:52,3]  = c4

cube[*,0:65,24:74,0] = c5
cube[*,0:65,0:50,1]  = c6
cube[*,0:65,24:74,2] = c7
cube[*,0:65,0:50,3]  = c8


  ; Remove outliers -- this is an efficient cosmic ray rejection
cube[where(abs(cube) ge 0.5)] = 0.
  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data
  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'c165_kn2_unwrapped_crrej.fits'

writefits,'cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr.fits', cube_mean, hdr
 
end
