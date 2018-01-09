pro coadd

c1=readfits('m141127_0310.fits',hdr) 
c2=readfits('m141127_0311.fits',hdr)
c3=readfits('m141127_0312.fits',hdr)
c4=readfits('m141127_0313.fits',hdr)


imA = (c1+c3)/2.
imB = (c2+c4)/2.
c_add= imA-imB
meanclip,c_add,meanim
c_add(where(abs(c_add) ge 80)) = meanim
gal = c_add(1610:1690,935:1025)
writefits,'c164_all.fits',c_add
writefits,'c164.fits',gal

;coadd the new data(new slit position with ImA only)
c1=readfits('m141129_0346.fits',hdr) 
c2=readfits('m141129_0347.fits',hdr)
c3=readfits('m141129_0348.fits',hdr)
c4=readfits('m141129_0349.fits',hdr)

imA = (c1+c3)/2.
imB = (c2+c4)/2.
c_add= imA-imB
meanclip,c_add,meanim
c_add(where(abs(c_add) ge 100)) = meanim
galpos = c_add(1600:1690,935:1025)
galneg = c_add(1594:1684,1047:1137)
writefits,'c164im1_all.fits',c_add
writefits,'c164im1.fits',galpos-galneg





stop
end
