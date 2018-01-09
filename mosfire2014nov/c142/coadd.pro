pro coadd

c1=readfits('m141127_0303.fits',hdr) 
c2=readfits('m141127_0304.fits')
c3=readfits('m141127_0305.fits')
c4=readfits('m141127_0306.fits')
c5=readfits('m141127_0307.fits')
c6=readfits('m141127_0308.fits')


imA = (c1+c3+c5)/3.
imB = (c2+c4+c6)/3.
c_add= imA-imB
gal = c_add(884:938,975:1035)
gal(42,25)=0. ;there is a cosmic ray here.
gal_no_skyline = gal
gal_no_skyline(24:31,*)=3.7
writefits,'c142_all.fits',c_add
writefits,'c142.fits',gal
writefits,'c142_noskyline.fits',gal_no_skyline
;stop
end
