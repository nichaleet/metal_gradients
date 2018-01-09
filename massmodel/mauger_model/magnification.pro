pro magnification

imsource= readfits('/scr2/nichal/workspace/output/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube_sourceplane_interp.fits')
source_Detection = imsource(*,*,4)
imim = readfits('/scr2/nichal/workspace/output/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube.fits')
im_detection = imim(*,*,4)
source_Res = 0.1/2.5 ;arcsec per pix
im_res = 0.1 ;arcsec per pix
area_source = n_Elements(where(finite(source_Detection)))*source_res
area_im= n_Elements(where(finite(im_detection)))*im_res
mag=area_im/area_source
print, mag
stop
end
