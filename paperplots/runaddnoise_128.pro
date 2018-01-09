pro runaddnoise_128

for i=0,10 do begin
   file='/scr2/nichal/workspace/output/addnoise/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube_noise'+strtrim(string(i),2)+'_sourceplane_interp.fits'
   name = 'cswa128_noise'+strtrim(string(i),2)
   ha_acubetoeps,file,246,238,34,'kpc', .336,name
endfor
end
