pro mosaic_mac1133handmosaic

;Below is for hand mosaic files
file1='/scr2/nichal/workspace/reduced_data/Mar2013/scaledskysubtracted/macs1133/macs1133Mar_Ha_Hn4_handmosaic_scaledsky_230hr.fits'
file2='/scr2/nichal/workspace/reduced_data/Feb2014/scaledskysubtracted/MACS1133/macs1133Feb_Ha_Hn4_handmosaic_scaledsky_130hr.fits'

cube1 = readfits(file1,head1)
cube2 = readfits(file2,head2)

;Below will give the coordinate of head's crpix in head1 coordinate
pixmatch,sxpar(head1,'PA_SPEC'),0.1,0.1,sxpar(head1, 'CRVAL2'),sxpar(head2, 'CRVAL2'),sxpar(head1, 'CRVAL3'),sxpar(head2, 'CRVAL3'),sxpar(head1, 'CRPIX2'),sxpar(head1, 'CRPIX3'),sxpar(head2, 'CRPIX2'),sxpar(head2, 'CRPIX3'),x_of_ref2_in1,y_of_ref2_in1
;Found that (32,28) of image2 is (38,50) in image1.

;The 2 images are up-side down to each other so have to invert both
;size of image1(mar) is 68x77
;size of image2(feb) is 66x77

;2.5 hours of the first file
cube =make_array(426,70,78,8)
cube[*,0:67,0:76,0]=cube1
cube[*,0:67,0:76,1]=cube1
cube[*,0:67,0:76,2]=cube1
cube[*,0:67,0:76,3]=cube1
cube[*,0:67,0:76,4]=cube1

;1.5 hours of the second file
for ii=0,65 do for jj=0,76 do cube[*,69-ii,77-jj,5]=cube2[*,ii,jj]
cube[*,*,*,6]=cube[*,*,*,5]
cube[*,*,*,7]=cube[*,*,*,5]

  ; Take the average of all non-rejected data
tot = total(cube,4)
nframes_4d = cube
nframes_4d[where(cube)] = 1.
nframes = total(nframes_4d,4)
cube_mean = tot/nframes  ; average flux per pixel, with bad values rejected
cube_mean[where(nframes eq 0)] = 0.  ; pixels with no data
  ; Unwrap this data cube and write the result to a FITS file
unwrap, smooth(cube_mean,[1,3,3]), 'm1133_unwrapped_crrej.fits'

writefits,'macs1133_Ha_Hn4_handmosaic_scaledsky_4hr.fits', cube_mean, head1
 

end
