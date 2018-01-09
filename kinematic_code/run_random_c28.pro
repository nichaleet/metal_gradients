pro run_random_c28
  openw,1,'cswa28result.dat'
  printf,1,'CSWA28 Kineticmodel Result:'
  close,1
  bestmodel_path = '/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMredshift2fixed_cswa28/outputmaps_new/prettymaps/crop/'
  bestmodel_file = 'sourcekinematic_pretty.fits'
  bestmodel_errfile = 'sourcekinematic_err_pretty.fits'

  random_path = '/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMredshift2fixed_cswa28/random_model/outputmaps/prettymaps/crop/'

  kernelfile='/scr2/nichal/workspace/psf_sourceplane/cswa28_psf_fit.fits'
  kernel = readfits(kernelfile)
  kernel = kernel[5:11,7:14]
  kinematicmodel,bestmodel_path+bestmodel_file,bestmodel_path+bestmodel_errfile,bestmodel_path+'sourceha_detection_pretty.fits',kernel,45.,90.,13,16,15,23,'cswa28',vcmax=130.,vmeanmin=-50.,vmeanmax=50.,rtmax=5.,pixscale=0.424

  for i=0,19 do begin
     openw,1,'cswa28result.dat',/append
     printf,1,'CSWA28 random',strtrim(string(i),2)
     close,1
     velfile = random_path+'sourcekinematic_random'+strtrim(string(i),2)+'_pretty.fits'
     velerrfile = random_path+'sourcekinematic_err_random'+strtrim(string(i),2)+'_pretty.fits'
     hadetect=random_path+'sourcehadetection_random'+strtrim(string(i),2)+'_pretty.fits'
     ha = readfits(hadetect)
     xyind = array_indices(ha,where(finite(ha)))
     xcmin = min(xyind[0,*])
     xcmax = xcmin+4.
     ycmax = max(xyind[1,*])
     ycmin = ycmax-10
     kinematicmodel,velfile,velerrfile,hadetect,kernel,45.,90.,xcmin,xcmax,ycmin,ycmax,'cswa28_random'+strtrim(string(i),2),vcmax=130.,vmeanmin=-50.,vmeanmax=50.,rtmax=5.,pixscale=0.424
    
  endfor
end
