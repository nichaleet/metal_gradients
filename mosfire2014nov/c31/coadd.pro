pro coadd

c1=readfits('m141126_0320.fits',hdr) 
c2=readfits('m141126_0321.fits')
c3=readfits('m141126_0322.fits')
c4=readfits('m141126_0323.fits')
c5=readfits('m141126_0324.fits')
c6=readfits('m141126_0325.fits')

imA = (c1+c3+c5)/3.
imB = (c2+c4+c6)/3.
c_add= imA-imB
gal = c_add(1000:1100,940:1100)
meanclip,c_add,meanim
gal(59:60,110)=meanim
writefits,'c31.fits',gal

end
