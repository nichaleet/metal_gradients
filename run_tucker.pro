pro run_tucker
  names  = ['0744','1148']
  pixarcsec = [0.005,0.01] ;arcsec per pixel
  z = [2.21,2.38]
  xc = [54.,58]
  yc = [50.,60]
  angle = [63.3,93.6]
  angle_bpt = [116.7,86.4]
  incline = [71.,32.]
  restricted_rmax=[1.2,2.0]
  angdist = lumdist(z)/(1.+z)^2               ;mpc
  angsize = 1./3600.*!dpi/180.*1.e3*angdist   ;kpc per arcsec
  pixscale = pixarcsec*angsize ;kpc per pixel

  path =  '/scr2/nichal/workspace/reduced_data/Tuckers/'

;  for ii=0,n_elements(names)-1 do begin
;  name = names[ii]
;  rotationcurve_fix,path,name+'_vel.fits',name+'_velerr.fits',name+'_veldisp.fits',name+'_veldisperr.fits',name+'_ha_detect.fits',xc[ii],yc[ii],angle[ii],pixscale[ii],name+'_NII_Ha.fits',name+'_NII_Ha_err.fits','Bay1','Bay2',name
;  endfor


  cd,'/scr2/nichal/workspace/bptmetalanalysis/'
  for ii=0,n_elements(names)-1 do begin
     name=names[ii]
     bptmetalanalysis,name,path,name+'_hadetect.fits',name+'_NII.fits',name+'_NI_err.fits',name+'_Ha.fits',name+'_Haerr.fits','nope','nope','nope','nope',angle=angle[ii],xc=xc[ii],yc=yc[ii],pixelscale=pixscale[ii],slitwidth=3.,inc=incline[ii],pos=angle_bpt[ii],minhadetect=10.,restrict_rmax=restricted_rmax[ii],detectlim=0.004
  endfor
end
