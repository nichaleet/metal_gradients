pro interp_map,file0,fileref


oldim = readfits(file0,oldhd)
refim = readfits(fileref,refhd)

;oldim(where(finite(oldim) eq 0)) = 0.
hastrom,oldim,oldhd,newim,newhd,refhd,missing=0.

filename = strmid(file0,0,strpos(file0,'.fits'))
writefits,filename+'_interp.fits',newim,newhd
write_tiff,filename+'_interp.tif',newim
;stop
end
