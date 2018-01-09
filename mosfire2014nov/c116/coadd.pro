pro coadd

c1=readfits('m141126_0182.fits',hdr) 
c2=readfits('m141126_0183.fits',hdr)
c3=readfits('m141126_0184.fits',hdr)
c4=readfits('m141126_0185.fits',hdr)
c5=readfits('m141126_0186.fits',hdr)
c6=readfits('m141126_0187.fits',hdr)

imA = (c1+c3+c5)/3.
imB = (c2+c4+c6)/3.
c_add= imA-imB
gal = c_add(1066:1137,955:1025)
writefits,'c116.fits',gal

end
