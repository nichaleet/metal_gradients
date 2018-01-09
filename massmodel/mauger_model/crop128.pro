pro crop128
path = '/scr2/nichal/workspace/output/'
hacube = readfits(path+'cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube_sourceplane_interp.fits',hahdr)
hbcube = readfits(path+'cswa128_Hb_tlc_Hbb_pipelinemosaic_scaledsky_130hr_acube_aligned_sourceplane_interp.fits',hbhdr)
metalcube = readfits(path+'metallicity/CSWA128_metallicitynew_sourceplane_interp.fits',metalhdr)


ha_vel = hacube(*,*,0)
badha = where(ha_vel gt -62.) 
hb_vel = hbcube(*,*,0)
badhb = where(hb_vel gt -60.)

for ii=0,11 do begin
   framenow = hacube(*,*,ii)
   framenow(badha) = 1./0.
   framenow(242:245,243:247) = 1./0.
   hacube(*,*,ii) = framenow

   framenow = hbcube(*,*,ii)
   framenow(badhb) = 1./0.
   hbcube(*,*,ii) = framenow

   if ii le 6 then begin
      framenow = metalcube(*,*,ii)
      framenow(badha) = 1./0.
      framenow(badhb) = 1./0.
      framenow(242:245,243:247) = 1./0.
      metalcube(*,*,ii) = framenow
   endif
endfor

writefits,path+'cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube_sourceplane_interp_cropped.fits',hacube,hahdr
writefits,path+'cswa128_Hb_tlc_Hbb_pipelinemosaic_scaledsky_130hr_acube_aligned_sourceplane_interp_cropped.fits',hbcube,hbhdr
writefits,path+'metallicity/CSWA128_metallicitynew_sourceplane_interp_cropped.fits',metalcube,metalhdr
end
