pro coadd

c1=readfits('m141128_0194.fits',hdr) 
c2=readfits('m141128_0195.fits',hdr)
c3=readfits('m141128_0196.fits',hdr)
c4=readfits('m141128_0197.fits',hdr)


imA = (c1+c3)/2.
imB = (c2+c4)/2.
c_add= imA-imB
meanclip,c_add,meanim
gal = c_add(312:367,945:1005)
gal(5,48:49)=meanim
gal(44,24) = meanim
gal(19,56:57)=meanim
gal(18:19,47)=meanim
gal(where(gal ge 300))=meanim
writefits,'c163_all.fits',c_add
writefits,'c163.fits',gal

c1_reduced = readfits('rectified_c163slit_c163_star_K_A-B_0194-0195.fits',hdr)
c2_reduced = readfits('rectified_c163slit_c163_star_K_A-B_0196-0197.fits')
creduced=c1_reduced+c2_reduced
meanclip,creduced,meanim
creduced(where(abs(creduced) ge 5.))= meanim
creduced(497:498,974:975)=meanim
creduced(500:501,969)=meanim
writefits,'c163_reduced.fits',creduced,hdr
stop
end
