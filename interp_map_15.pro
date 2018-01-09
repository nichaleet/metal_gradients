;;Below is for /massmodel/BCK_script/LTMcswa15/LTM (for the version that we lower the resolution of OSIRIS to match the ESI
;pro interp_map_15,file0
  
;;The ESI image has a Y axis flipped. Size 552x552. The scale is 0.15" per pixel.

;oldim = readfits(file0,oldhd) ; read this frame
;;orient the Osiris image with the ESI
;hrot,oldim,oldhd,-1,-1,185,-1,-1,1,missing=0. ;now oldim and oldhd has a rotation that can be aligned with the ESI image

;;expand the scale of the osiris to match the ESI (Osiris is
;;0.1"/pix. ESI is 0.15"/pix)
;newsize = size(oldim)*2./3.
;hcongrid,oldim,oldhd,out=[newsize(1),newsize(2)],cubic=-0.5,/half

;;make the number of pixels in osiris frame match with the ESI
;;(0.5,50.5) of Osiris is 252.125 and 257.25 of ESI
;;(0.5,0.5) of Osiris is 252.125 and 207.25 of ESI
;newim = make_Array(552,552,/float,value=0.)
;newim[252:295,207:256]=oldim

;;fix the header
;newcrpix1 = sxpar(oldhd,'crpix1')
;newcrpix2 = sxpar(oldhd,'crpix2')
;fxaddpar,oldhd,'crpix1',newcrpix1
;fxaddpar,oldhd,'crpix2',newcrpix2
;fxaddpar,oldhd,'naxis1',552
;fxaddpar,oldhd,'naxis2',552

;filename = strmid(file0,0,strpos(file0,'.fits'))
;writefits,filename+'_interp.fits',newim,oldhd
;write_tiff,filename+'_interp.tif',newim
;stop
;end

;;Below is for /massmodel/BCK_script/LTMcswa15/LTMzoomESI (for the version that we increase the resolution of ESI to match the OSIRIS

pro interp_map_15,file0
  
;The ESI image has a Y axis flipped. Size 552x552. The scale is 0.15" per pixel.

oldim = readfits(file0,oldhd) ; read this frame
;orient the Osiris image with the ESI
hrot,oldim,oldhd,-1,-1,185,-1,-1,1,missing=0. ;now oldim and oldhd has a rotation that can be aligned with the ESI image


;make the number of pixels in osiris frame match with the ESI
;(0.5,0.5) of Osiris is 378.5 and 308.5 of ESI
newim = make_Array(828,828,/float,value=0.)
newim[377:442,307:381]=oldim

;fix the header
newcrpix1 = sxpar(oldhd,'crpix1')
newcrpix2 = sxpar(oldhd,'crpix2')
fxaddpar,oldhd,'crpix1',newcrpix1
fxaddpar,oldhd,'crpix2',newcrpix2
fxaddpar,oldhd,'naxis1',828
fxaddpar,oldhd,'naxis2',828

filename = strmid(file0,0,strpos(file0,'.fits'))
writefits,filename+'_interp.fits',newim,oldhd
;write_tiff,filename+'_interp.tif',newim
;stop
end
