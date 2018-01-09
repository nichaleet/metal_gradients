pro crop_c159_fits

file = ['frame-i-008116-6-0087.fits','cswa159_sdss_r_corr.fits','frame-g-008116-6-0087.fits']
outfile = ['sdss_i.fits','sdss_r.fits','sdss_g.fits']
;current angle = 75.627 EON
rot_angle = 360.-75.627
for i =0,2 do begin
   im = readfits(file(i),header)
   costheta = sxpar(header,'CD1_1')
   sintheta = sxpar(header,'CD1_2')
   oldscale = sqrt(costheta^2+sintheta^2) ;degree
   angle = acos(costheta/oldscale)
   angle = 180.+(angle*180./!pi)
   hrot,im,header,newim,newhd,angle,1877,517,1
   hextract,newim,newhd,940,1079,680,819 ;140*140 pix
   newscale = 0.1
   newsize  = round(140*oldscale*3600./newscale)
   ;stop
   hrebin,newim,newhd,imout,hdout,newsize,newsize
   writefits,outfile(i),imout,hdout
endfor

;Align images to r band image
Rim = readfits('sdss_r.fits',Rhd)
Iim = readfits('sdss_i.fits',Ihd)
Gim = readfits('sdss_g.fits',Ghd)

hastrom,Iim,Ihd,newIim,newIhd,Rhd,missing=0.
hastrom,Gim,Ghd,newGim,newGhd,Rhd,missing=0.

writefits,'sdss_i.fits',newIim,newIhd
writefits,'sdss_g.fits',newGim,newGhd
stop
end
