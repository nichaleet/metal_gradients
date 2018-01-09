

;;Below is for /massmodel/BCK_script/cswa139zoomESI/ (for the version that we increase the resolution of ESI to match the OSIRIS

pro interp_map_139,file0
;ESI image = /scr2/nichal/workspace/imaging_data/cswa139/cswa139_ESI_R_match_osirisRES.fits  
;The ESI image rotate  166.5 degree from OSIRIS image clockwise. Size 464x464. The scale is 0.1542" per pixel.

oldim = readfits(file0,oldhd) ; read this frame
;orient the Osiris image with the ESI
hrot,oldim,oldhd,-1,-1,193.,25,44,1,missing=0. ;now oldim and oldhd has a rotation that can be aligned with the ESI image

;make the number of pixels in osiris frame match with the ESI
;(35,35) of rotated Osiris is 252,188  of ESI
;(1,1) = (218,154)
newim = make_Array(464,464,/float,value=0.)
newim[217:285,153:225]=oldim

;fix the header
newcrpix1 = sxpar(oldhd,'crpix1')
newcrpix2 = sxpar(oldhd,'crpix2')
fxaddpar,oldhd,'crpix1',newcrpix1
fxaddpar,oldhd,'crpix2',newcrpix2
fxaddpar,oldhd,'naxis1',464
fxaddpar,oldhd,'naxis2',464

filename = strmid(file0,0,strpos(file0,'.fits'))
writefits,filename+'_interp.fits',newim,oldhd
;write_tiff,filename+'_interp.tif',newim
;stop
end
