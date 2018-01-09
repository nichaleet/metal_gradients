pro runaddnoise_bpt
slope=findgen(4,11)
sigma=slope
namearr=[]
for i=0,10 do begin
   name = 'cswa128_noise'+strtrim(string(i),2)
   path = '/scr2/nichal/workspace/output/chopped_sourceplane/cswa128_noise'+strtrim(string(i),2)+'/'
   bptmetalanalysis_addnoise,name,path,name+'_sourceha_detection.fits',name+'_sourcenii.fits',name+'_sourcenii_err.fits',name+'_sourceha.fits',name+'_sourceha_err.fits',name+'_sourceoiii.fits',name+'_sourceoiii_err.fits',name+'_sourcehb.fits',name+'_sourcehb_err.fits',param,sigma_param,angle=150.,xc=13,yc=16,pixelscale=0.336,slitwidth=1.,inc=57.,pos=210.,detectlim=0.0058,restrict_rmax=5.
   slope[*,i]=param
   sigma[*,i]=sigma_param
   namearr=[namearr,name]
endfor
save,namearr,slope,sigma,filename='/scr2/nichal/workspace/output/bptmetalanalysis_result/addnoise/slope.sav'
end
