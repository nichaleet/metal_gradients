pro relens_128_addnoise
for i=0,10 do begin
file = '/scr2/nichal/workspace/output/addnoise/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube_noise'+strtrim(string(i),2)+'.fits'
print,file
relens_128,file
endfor

end
