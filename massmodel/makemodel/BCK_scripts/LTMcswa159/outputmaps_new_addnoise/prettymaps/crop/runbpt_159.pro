pro runbpt_159

for i=0,10 do begin
   name ='cswa159_'+strtrim(string(i),2)
   path = '/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa159/outputmaps_new_addnoise/prettymaps/crop/'
   bptmetalanalysis_159,name,path, name+'sourceha_detection_pretty.fits',name+'sourcenii_pretty.fits',name+'sourcenii_err_pretty.fits',name+'sourceha_pretty.fits',name+'sourceha_err_pretty.fits','sourceoiii_pretty.fits','sourceoiii_err_pretty.fits','sourcehb_pretty.fits','sourcehb_err_pretty.fits',outparams,angle=5.,xc=13,yc=17,pixelscale=0.418205,slitwidth=1.5,inc=13,pos=175.
   outparams = [[outparams],[outparams]]
   stop
endfor

end
