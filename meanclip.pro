pro meanclipcswa15
  
;This program will stack Osiris data cube with (lambda,x,y) coordinate of all the files in the folder

; IDL script to examine OSIRIS data of the z=2.164 arc CSWA 15.

;;; H alpha + [N II]

; Mosaiced data cubes, reduced by Nicha:
cube = readfits('cswa15_Ha_Kn2_mosaic_sky_230hr.fits')
;unwrap, cube, 'c15_kn2_unwrapped.fits'
unwrap, smooth(cube,[1,3,3]), 'c15_kn2_unwrapped.fits'


; Add individual reduced data cubes in such a way that cosmic rays can be
rejected from each frame.
  ; First night
c1 = readfits('s140223_a019001_Kn2_100.fits', hdr)  ; "A" position
c2 = readfits('s140223_a019002_Kn2_100.fits')       ; "B" position
c3 = readfits('s140223_a019003_Kn2_100.fits')  ; A
c4 = readfits('s140223_a019004_Kn2_100.fits')  ; B
c5 = readfits('s140223_a020001_Kn2_100.fits')  ; A
c6 = readfits('s140223_a020002_Kn2_100.fits')  ; B
  ; Second night
c7 = readfits('s140224_a019001_Kn2_100.fits')  ; A (within ~1 pixel)
c8 = readfits('s140224_a019002_Kn2_100.fits')  ; B
c9 = readfits('s140224_a019003_Kn2_100.fits')  ; A
c10= readfits('s140224_a019004_Kn2_100.fits')  ; B
  ; Wavelength scale
wl = sxpar(hdr,'crval1') + sxpar(hdr,'cdelt1')*findgen(sxpar(hdr,'naxis1'))


; Data from both nights are almost exactly aligned (within ~0.1")!
; Data cubes are offset by 24 pixels along the short axis (axis 2)
  ; Combine all frames into a mega-array
  ; Shift each frame to align them spatially
cube = make_array(421,66,75,10)
cube[*,*,24:74,0] = c1
cube[*,*,0:50,1] = c2
cube[*,*,24:74,2] = c3
cube[*,*,0:50,3] = c4
cube[*,*,24:74,4] = c5
cube[*,*,0:50,5] = c6
cube[*,*,24:74,6] = c7
cube[*,*,0:50,7] = c8
cube[*,*,24:74,8] = c9
cube[*,*,0:50,9] = c10
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
unwrap, smooth(cube_mean,[1,3,3]), 'c15_kn2_unwrapped_crrej.fits'
 
