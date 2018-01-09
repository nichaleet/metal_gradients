function broadgauss,wl,zB,vdispB,NII_Ha,SIIa_Ha,SIIb_SIIa,HaB
  COMMON LINEWAVE, Ha,Hb,NII,OIII,SIIa,SIIb,clight
 ;Broad Region
  SIIaWl   = SIIa*(zB+1.)
  SIIbWl   = SIIb*(zB+1.)
  HaWl     = Ha*(zB+1.)
  NIIWl    = NII*(zB+1.)
  sigmaHa  = vdispB/clight*hawl
  sigmaNII  = vdispB/clight*NIIwl
  sigmaSIIa = vdispB/clight*SIIawl
  sigmaSIIb  = vdispB/clight*SIIbwl
  broadspec = HaB*exp(-0.5*(wl-hawl)^2/sigmaha^2)+HaB*Nii_Ha*exp(-0.5*(wl-NIIwl)^2/sigmaNII^2)+HaB*SIIa_Ha*exp(-0.5*(wl-SIIawl)^2/sigmaSIIa^2)+HaB*SIIa_Ha*SIIb_SIIa*exp(-0.5*(wl-SIIbwl)^2/sigmaSIIb^2)
  return,broadspec
end

function narrowgauss,wl,zN,vdispN,NII_Ha,SIIa_Ha,SIIb_SIIa,HaN
  COMMON LINEWAVE, Ha,Hb,NII,OIII,SIIa,SIIb,clight
 ;Broad Region
  ;Narrow Region
  SIIaWl   = SIIa*(zN+1.)
  SIIbWl   = SIIb*(zN+1.)
  HaWl     = Ha*(zN+1.)
  NIIWl    = NII*(zN+1.)
  sigmaHa  = vdispN/clight*hawl
  sigmaNII  = vdispN/clight*NIIwl
  sigmaSIIa = vdispN/clight*SIIawl
  sigmaSIIb  = vdispN/clight*SIIbwl
  narrowspec = HaN*exp(-0.5*(wl-hawl)^2/sigmaha^2)+HaN*Nii_Ha*exp(-0.5*(wl-NIIwl)^2/sigmaNII^2)+HaN*SIIa_Ha*exp(-0.5*(wl-SIIawl)^2/sigmaSIIa^2)+HaN*SIIa_Ha*SIIb_SIIa*exp(-0.5*(wl-SIIbwl)^2/sigmaSIIb^2)
  return,narrowspec
end

function twogaussAGN,wl,a
COMMON LINEWAVE, Ha,Hb,NII,OIII,SIIa,SIIb,clight
;pi.parname = ['    const', '  zN','   ZB', 'vdispN','vdispB','NII_Ha','SIIa_Ha','SIIb_SIIa']
  const  = a[0]
  zN     = a[1]
  zB     = a[2]
  vdispN = a[3]
  vdispB = a[4]
  Nii_Ha = a[5]
  SIIa_Ha   = a[6]
  SIIb_SIIa = a[7]
  HaN = a[8]
  HaB = a[9]
  broadspec=broadgauss(wl,zB,vdispB,NII_Ha,SIIa_Ha,SIIb_SIIa,HaB)
  narrowspec=narrowgauss(wl,zN,vdispN,NII_Ha,SIIa_Ha,SIIb_SIIa,HaN)
  specstruct = const+broadspec+narrowspec
  return,specstruct
end

function twogaussHa,wl,a
COMMON LINEWAVE, Ha,Hb,NII,OIII,SIIa,SIIb,clight
;pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
  const = a[0]
  z     = a[1]
  vdisp = a[2]
  HHa   = a[3]
  HNII   = a[4]
  HaWl   = Ha*(z+1.)
  NIIWl   = NII*(z+1.)
  sigmaHa = vdisp/clight*hawl
  sigmaNII  = vdisp/clight*NIIwl
  specstruct = const+HHa*exp(-0.5*(wl-hawl)^2/sigmaha^2)+HNII*exp(-0.5*(wl-NIIwl)^2/sigmaNII^2)
  return,specstruct
end

function twogaussSII,wl,a
COMMON LINEWAVE, Ha,Hb,NII,OIII,SIIa,SIIb,clight
;pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
  const = a[0]
  z     = a[1]
  vdisp = a[2]
  HSIIa   = a[3]
  HSIIb   = a[4]
  SIIaWl   = SIIa*(z+1.)
  SIIbWl   = SIIb*(z+1.)
  sigmaSIIa = vdisp/clight*SIIawl
  sigmaSIIb  = vdisp/clight*SIIbwl
  specstruct = const+HSIIa*exp(-0.5*(wl-SIIawl)^2/sigmaSIIa^2)+HSIIb*exp(-0.5*(wl-SIIbwl)^2/sigmaSIIb^2)
  return,specstruct
end

function main,filehacube,filehadetect,filebptmap,filehadetect_source,z,clustername,filter,lowersig=lowersig
  COMMON LINEWAVE, Ha,Hb,NII,OIII,SIIa,SIIb,clight
  COMMON PATH, pathbptout,pathcube
;look for SiII lines and broad line region from the identified potential AGN region
  cube3D = readfits(filehacube,cube3Dhdr)
  cube3Derr = sqrt(cube3D)
  Hadetect_image = readfits(filehadetect)
  Hadetect_image = Hadetect_image[*,*,4]
  bptmap = readfits(filebptmap)
  hadetect_source = readfits(filehadetect_source)
  
  ;AGN region in source plane
  bptreg = where(bptmap eq 0.5)
  !p.multi=[0,1,1]
  plothist,hadetect_source(bptreg),bin=5,xtitle='significance'
  if ~keyword_set(lowersig) then read,lowersig,prompt='Input lower detection significance: '

  ;check size
  sizecube = size(cube3d,/dimensions)
  sizeim   = size(hadetect_image,/dimensions)
  if sizecube[1] ne sizeim[0] or sizecube[2] ne sizeim[1] then stop

  ;AGN region in image plane
  bptreg_im = where(Hadetect_image gt lowersig and finite(Hadetect_image))
  npix = n_elements(bptreg_im)
  toshow = fltarr(sizeim)+1
  toshow(where(finite(hadetect_image))) = 3
  toshow(bptreg_im) = 7
  window,0,xsize=500,ysize=800,title=clustername
  cgLoadCT, 33, CLIP=[30,255]
  !p.multi=[0,1,4]
  !p.font = 0
  position =   [0.1, 0.1, 0.9, 0.9]
  cgImage,toshow , Stretch=1,/Axes, XTitle='x',YTitle='y',XRange=[0,sizeim[0]-1], YRange=[0,sizeim[1]-1],font=0,position=position

  ;take the spectra from AGN region
  ;Make SPEC1D and Plot
  spec1d = fltarr(n_elements(cube3D[*,0,0]))
  spec1derr = spec1d
  for i =0, n_elements(spec1d)-1 do begin
     Imi    = cube3D[i,*,*]
     ImiErr = cube3DErr[i,*,*]
     spec1d[i] = total(imi(bptreg_im))
     ;spec1derr[i] = sqrt(total((Imierr(bptreg_im))^2))
     spec1derr[i] = stdev(Imi(bptreg_im))*sqrt(npix)
  endfor
  wl = getwl_filter(filter)
  Cgplot,wl,spec1d,color=fsc_color('black'),xrange=minmax(wl),/nodata,title='Fitting Gassian for one width component'
  oplot,wl,spec1d+spec1derr,color=fsc_color('rose')
  oplot,wl,spec1d-spec1derr,color=fsc_color('rose')
  oplot,wl,spec1d,color=fsc_color('black')

  ;FIT gaussian to Ha,NII line
  pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)
  pi.value = double([median(spec1d),z,200.0,max(spec1d),1])
  pi[0].limits = [-1,1]
  pi[1].limits = [z-0.1,z+0.1]
  pi[2].limits = [0.0, 1000.0]
  pi[3].limits = minmax(spec1d)
  pi[4].limits = minmax(spec1d)
  pi.step = double([0.1, 0.001, 25.0,0.001,0.001])
  pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
  pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)']
  mask = bytarr(n_elements(wl))+1
  RegOff = where(wl lt Ha*(1.+z)*(1.-1500./clight) or wl gt NII*(1+z)*(1+1500./clight))
  ;RegOff = [RegOff,where(abs(spec1d-median(spec1d)) gt 2.*stdev(spec1d))]
  mask(RegOff) = 0
  RegIn  = [where(wl gt Ha*(1.+z)*(1.-400./clight) and wl lt Ha*(1.+z)*(1.+400./clight)),where(wl gt NII*(1.+z)*(1.-400./clight) and wl lt OIII*(1.+z)*(1.+400./clight))]
  mask(RegIn) = 1
  mask[0,(where(~finite(spec1d)))] = 0
  won = where(mask eq 1)
  oplot,wl(won),spec1d(won),color=fsc_color('dark green')
  ;fix peak Ha
  pi[3].fixed = 1
  params = MPFITFUN('twogaussHa',wl(won),spec1d(won),spec1derr(won),parinfo=pi,/nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit)
  wlmodel = wl(where(wl gt min(wl(won)) and wl lt max(wl(won))))
  model = twogaussHa(wlmodel,params)
  oplot,wlmodel,model,color=fsc_color('purple')
  z = params[1]
  vdisp = params[2]
  xyouts,mean(!x.crange),0.7*params[3],'Halpha dispersion='+string(vdisp,format='(I4)')+'km/s',color=fsc_color('black')
  ;calculate NII/Ha of the central region
  NII_Ha = params[4]*NII/(params[3]*Ha)
  NII_Haerr = NII_Ha*sqrt((perror[4]/params[4])^2+(perror[3]/params[3])^2)
  print,'NII_Ha =',NII_Ha,NII_Haerr

  ;FIT gaussian to SII lines
  pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)
  pi.value = double([median(spec1d),z,200.0,1,1])
  pi[0].limits = [-1,1]
  pi[1].limits = [z-0.1,z+0.1]
  pi[2].limits = [0.0, 1000.0]
  pi[3].limits = minmax(spec1d)
  pi[4].limits = minmax(spec1d)
  pi.step = double([0.1, 0.001, 25.0,0.001,0.001])
  pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
  pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)']
  mask = bytarr(n_elements(wl))+1
  RegOff = where(wl lt SIIa*(1.+z)*(1.-1500./clight) or wl gt SIIb*(1+z)*(1+1500./clight))
  mask(RegOff) = 0
  mask[0,(where(~finite(spec1d)))] = 0
  won = where(mask eq 1)
  oplot,wl(won),spec1d(won),color=fsc_color('dark green')
  paramsSII = MPFITFUN('twogaussSII',wl(won),spec1d(won),spec1derr(won),parinfo=pi,/nocatch, bestnorm=bestnorm, dof=dof, perror=perrorSII, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit)
  wlmodel = wl(where(wl gt min(wl(won)) and wl lt max(wl(won))))
  model = twogaussSII(wlmodel,paramsSII)
  oplot,wlmodel,model,color=fsc_color('purple')
  zSII = paramsSII[1]
  vdispSII = paramsSII[2]
  print,'Ha vdisp = '+strtrim(string(vdisp,format='(F7.2)'),2)+'+-'+strtrim(string(perror[2],format='(F7.2)'),2)+'km/s'
  print,'SII vdisp = '+strtrim(string(vdispSII,format='(F7.2)'),2)+'+-'+strtrim(string(perrorSII[2],format='(F7.2)'),2)+'km/s'
  ;calculate SII_Ha
  SII_Ha = (paramsSII[3]*SIIa+paramsSII[4]*SIIb)*vdispSII/(params[3]*Ha*vdisp)
  A = (paramsSII[3]*SIIa+paramsSII[4]*SIIb)
  B = (params[3]*Ha)
  C = vdispSII/vdisp
  ErrA = sqrt((perrorSII[3]*SIIa)^2+(perrorSII[4]*SIIb)^2)
  ErrB = Ha*perror[3]
  ErrC = vdispSII/vdisp*sqrt((perrorSII[2]/vdispSII)^2+(perror[2]/vdisp)^2)
  SII_HaErr = A*C/B*sqrt((A^2*ErrC^2+C^2*ErrA^2)/(A*C)^2+(ErrB/B)^2)
  print, 'SII_Ha =',SII_Ha,SII_HaErr
  print, pi.parname
  print, params
  print, paramsSII
  vline,(z+1)*SIIa
  vline,(z+1)*SIIb

  
  ;Fit gaussian AGN with broad and narrow line component simultaneously
  spec1dN = spec1d/params[3] 
  spec1dNerr = spec1derr/params[3]
  Cgplot,wl,spec1dN,color=fsc_color('black'),xrange=minmax(wl),title='Fitting Gassian for two width components'
  pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 10)
  pi.value = double([median(spec1dN),z,z,100.0,500.,0.1,0.05,1.,1.,0.1])
  pi[0].limits = [-1,1]
  pi[1].limits = [z-0.1,z+0.1]
  pi[2].limits = [z-0.1,z+0.1]
  pi[3].limits = [0.0, 1000.0]
  pi[4].limits = [100., 2000.0]
  pi[5].limits = [0.,1]
  pi[6].limits = [0,1.]
  pi[7].limits = [0,5.]
  pi[8].limits = [0,1.]
  pi[9].limits = [0,1.]
  pi.step = double([0.1, 0.001,0.001, 25.,25.,0.001,0.001,.001,0.01,0.01])
  pi.parname = ['    const', '  zN','  zB', 'vdispN','vdispB','NII_Ha','SIIa_Ha','SIIb_SIIa','HaN','HaB']
  pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)','(D6.3)']
  mask = bytarr(n_elements(wl))
  RegIN1 = where(wl gt Ha*(1.+z)*(1.-1500./clight) and wl lt NII*(1+z)*(1+1500./clight)) 
  won1 = regin1
  RegIN2 = where(wl gt SIIa*(1.+z)*(1.-1500./clight) and wl lt SIIb*(1+z)*(1+1500./clight))
  won2 = regin2
  RegIN = [RegIn1,RegIn2]
  mask(RegIn) = 1
  won = where(mask eq 1)
  oplot,wl(won1),spec1dN(won1),color=fsc_color('dark green')
  oplot,wl(won2),spec1dN(won2),color=fsc_color('dark green')
  paramsAGN = MPFITFUN('twogaussAGN',wl(won),spec1dN(won),spec1dNerr(won),parinfo=pi,/nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit)
  wlmodel = wl(won1)
  model = twogaussAGN(wlmodel,paramsAGN)
  oplot,wlmodel,model,color=fsc_color('purple')
  wlmodel = wl(won2)
  model = twogaussAGN(wlmodel,paramsAGN)
  oplot,wlmodel,model,color=fsc_color('purple')

  Cgplot,wl,spec1dN,color=fsc_color('black'),xrange=minmax(wl),title='Fitting Gassian for two width components:Show separate components'
  oplot,wl(won1),spec1dN(won1),color=fsc_color('dark green')
  oplot,wl(won2),spec1dN(won2),color=fsc_color('dark green')
  wlmodel = wl(won1)
  modelB = broadgauss(wlmodel,paramsAGN[2],paramsAGN[4],paramsAGN[5],paramsAGN[6],paramsAGN[7],paramsAGN[9])
  modelN = narrowgauss(wlmodel,paramsAGN[1],paramsAGN[3],paramsAGN[5],paramsAGN[6],paramsAGN[7],paramsAGN[8])
  oplot,wlmodel,paramsAGN[0]+modelB,color=fsc_color('blue')
  oplot,wlmodel,paramsAGN[0]+modelN,color=fsc_color('red')
  wlmodel = wl(won2)
  modelB = broadgauss(wlmodel,paramsAGN[2],paramsAGN[4],paramsAGN[5],paramsAGN[6],paramsAGN[7],paramsAGN[9])
  modelN = narrowgauss(wlmodel,paramsAGN[1],paramsAGN[3],paramsAGN[5],paramsAGN[6],paramsAGN[7],paramsAGN[8])
  oplot,wlmodel,paramsAGN[0]+modelB,color=fsc_color('blue')
  oplot,wlmodel,paramsAGN[0]+modelN,color=fsc_color('red')

  vdispb = paramsAGN[4]
  vdispn = paramsAGN[3]
  xyouts,max(wl(won1)),0.85,'Broad dispersion='+string(vdispb,format='(I4)')+'km/s',color=fsc_color('black')
  xyouts,max(wl(won1)),0.7,'Narrow dispersion='+string(vdispn,format='(I4)')+'km/s',color=fsc_color('black')

  ;contribution
  Habroad_HaNarrow = paramsAGN[4]*paramsAGN[9]/paramsAGN[3]/paramsAGN[8]*100.
  xyouts,max(wl(won1)),0.55,'Ha Broad Contribution ='+string(Habroad_HaNarrow,format='(F5.2)')+'%',color=fsc_color('black')

  
  stop
return, 4
end

pro agnanalysis
  COMMON LINEWAVE, Ha,Hb,NII,OIII,SIIa,SIIb,clight
  COMMON PATH, pathbptout,pathcube
  Ha = .656461
  Hb = .486269
  NII= .658523
  OIII= .500824
  SIIa = .671829
  SIIb = .6732267
  clight = 299792.458
  pathout = '/scr2/nichal/workspace/output/'
  pathcube= '/scr2/nichal/workspace/reduced_data/mosaic/'
  pathbptout='/scr2/nichal/workspace/output/bptmetalanalysis_result/'

;cswa15
  filebptmap = pathbptout+'cswa15_bptmap.fits'
  filehadetect_source = '/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa15/LTMzoomESI/outputmaps_new/prettymaps/crop/sourceha_detection_pretty.fits'
  filehadetect = pathout+'cswa15_Ha_tlc_Kn2_handmosaic_sky_230hr_acube.fits'
  filehacube = pathcube+'cswa15_Ha_tlc_Kn2_handmosaic_sky_230hr.fits'
  filehacube = pathcube+'cswa15_Ha_Kn2_handmosaic_sky_230hr.fits'
  out = main(filehacube,filehadetect,filebptmap,filehadetect_source,2.16,'cswa15','Kn2',lowersig=20)

end
