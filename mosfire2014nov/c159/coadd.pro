pro coadd

c1=readfits('m141127_0167.fits',hdr) 
c2=readfits('m141127_0168.fits',hdr)
c3=readfits('m141127_0169.fits',hdr)
c4=readfits('m141127_0170.fits',hdr)


imA = (c1+c3)/2.
imB = (c2+c4)/2.
c_add= imA-imB
gal = c_add(968:1027,982:1015)
writefits,'c159_all.fits',c_add
writefits,'c159.fits',gal

c1_reduced = readfits('rectified_c159long_c159_star_K_A-B_0167-0168.fits',hdr)
c2_reduced = readfits('rectified_c159long_c159_star_K_A-B_0169-0170.fits')
writefits,'c159_reduced.fits',c1_reduced+c2_reduced,hdr

end
