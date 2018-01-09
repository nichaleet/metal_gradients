function twogaussHa,wl,a
COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
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

function twogaussHb,wl,a
COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
;pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
  const = a[0]
  z     = a[1]
  vdisp = a[2]
  Hhb   = a[3]
  Ho3   = a[4]
  hbwl   = Hb*(z+1.)
  o3wl   = OIII*(z+1.)
  sigmahb = vdisp/clight*hbwl
  sigmao3  = vdisp/clight*o3wl
  specstruct = const+Hhb*exp(-0.5*(wl-hbwl)^2/sigmahb^2)+Ho3*exp(-0.5*(wl-o3wl)^2/sigmao3^2)
  return,specstruct
end

function calculate_reddening,filehacube,filehadetect,filehbcube,filehbdetect,calibha,calibhb,z,filterha,filterhb,stoptag=stoptag
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  !p.multi=[0,1,2]
  HaCube  = readfits(fileHaCube,HdrHa)
  HaDetect= readfits(fileHaDetect,HdrHaDetect)
  HaCubeErr = sqrt(abs(HaCube))
  gooddetect = where(finite(HaDetect[*,*,0]))
  npixdetect = float(n_Elements(gooddetect))

  ;Make SPEC1D and Plot
  spec1d = fltarr(n_elements(HaCube[*,0,0]))
  spec1derr = spec1d
  for i =0, n_elements(spec1d)-1 do begin
     Imi    = HaCube[i,*,*]
     ImiErr = HaCubeErr[i,*,*]
     spec1d[i] = total(imi(gooddetect))
     ;spec1derr[i] = sqrt(total((Imierr(gooddetect))^2))
     spec1derr[i] = stdev(Imi(gooddetect))*sqrt(npixdetect)
  endfor
  wl = getwl_filter(filterHa)
  !p.multi=[0,1,2]
  Cgplot,wl,spec1d,color=fsc_color('black'),title='Ha'
  oplot,wl,spec1d+spec1derr,color=fsc_color('rose')
  oplot,wl,spec1d-spec1derr,color=fsc_color('rose')
  oplot,wl,spec1d,color=fsc_Color('black')
  ;FIT gaussian to Ha and NII line
  pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)
  pi.value = double([median(spec1d),z,200.0,0.5,0.5])
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
  mask(RegOff) = 0
  won = where(mask)
  oplot,wl(won),spec1d(won),color=fsc_color('dark green')
  params = MPFITFUN('twogaussHa',wl(won),spec1d(won),spec1derr(won),parinfo=pi,/nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit)
  model = twogaussHa(wl(won),params)
  oplot,wl(won),model,color=fsc_color('red')
  z = params[1]
  vdisp = params[2]
  sigmaHa = Ha*(z+1.)*vdisp/clight
  sigmaHaerr = sigmaHa*sqrt((perror[2]/vdisp)^2+(perror[1]/z)^2)
  pixres = wl(1)-wl(0) ; micron/pixel
  Haflux = sqrt(2.*!pi)*params[3]*sigmaHa*calibHa.conv_factor/pixres ;erg/s/cm^2
 ;The integration integrates ADU/s/keck_area/pixel over micron so have to divide by micron/pixel
  Hafluxerr = Haflux*sqrt((sigmaHaerr/sigmaHa)^2+(perror[3]/params[3])^2+(calibHa.conv_factorerr/calibha.conv_factor)^2)
  print, 'Ha Flux '+string(Haflux,format='(E10.3)')+'+-'+string(Hafluxerr,format='(E10.3)')
  vline,(z+1)*Ha
  wait,1
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  HbCube  = readfits(fileHbCube,HdrHb)
  HbDetect= readfits(fileHbDetect,HdrHbDetect)
  HbCubeErr = sqrt(abs(HbCube))
  gooddetect = where(finite(HbDetect[*,*,0]))
  npixdetect = float(n_Elements(gooddetect))
  ;Make SPEC1D and Plot
  spec1d = fltarr(n_elements(HbCube[*,0,0]))
  spec1derr = spec1d
  for i =0, n_elements(spec1d)-1 do begin
     Imi    = HbCube[i,*,*]
     ImiErr = HbCubeErr[i,*,*]
     spec1d[i] = total(imi(gooddetect))
     ;spec1derr[i] = sqrt(total((Imierr(gooddetect))^2))
     spec1derr[i] = stdev(Imi(gooddetect))*sqrt(npixdetect)
  endfor
  wl = getwl_filter(filterHb)
  Cgplot,wl,spec1d,color=fsc_color('black'),title='Hb'
  oplot,wl,spec1d+spec1derr,color=fsc_color('rose')
  oplot,wl,spec1d-spec1derr,color=fsc_color('rose')
  oplot,wl,spec1d,color=fsc_Color('black')
  vline,(z+1)*Hb
  vline,(z+1)*OIII
  ;FIT gaussian to Hb and OIII line
  pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)
  pi.value = double([median(spec1d),z,vdisp,0.3,1])
  pi[0].limits = [-1,1]
  pi[1].limits = [z-0.1,z+0.1]
  pi[2].limits = [vdisp-50.,vdisp+50.]
  pi[3].limits = minmax(spec1d)
  pi[4].limits = minmax(spec1d)
  pi.fixed = [0,0,1,0,0]
  pi.step = double([0.1, 0.001, 25.0,0.001,0.001])
  pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
  pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)']
  mask = bytarr(n_elements(wl))+1
  RegOff = where(wl lt Hb*(1.+z)*(1.-1500./clight) or wl gt OIII*(1+z)*(1+1500./clight))
  RegOff = [RegOff,where(abs(spec1d-median(spec1d)) gt 2.*stdev(spec1d))]
  mask(RegOff) = 0
  RegIn  = [where(wl gt Hb*(1.+z)*(1.-200./clight) and wl lt Hb*(1.+z)*(1.+200./clight)),where(wl gt OIII*(1.+z)*(1.-200./clight) and wl lt OIII*(1.+z)*(1.+200./clight))]
  mask(RegIn) = 1
  won = where(mask)
  won1 = won
  if filehbcube eq '/scr2/nichal/workspace/reduced_data/mosaic/cswa159_Hb_tlc_Hn3_handmosaic_scaledsky_1hr.fits' then begin
     skyline1 = where(wl gt 1.6024 and wl lt 1.604)
     middle = where(wl gt 1.61); and wl lt 1.64)
     mask([skyline1,middle]) = 0
     pi[3].limits = [0.1,0.4]
     won = where(mask)
  endif
  if filehbcube eq '/scr2/nichal/workspace/reduced_data/mosaic/cswa128_Hb_tlc_Hbb_pipelinemosaic_scaledsky_130hr.fits' then begin
     mask(where(wl gt 1.585)) = 0
     won = where(mask)
  endif
  oplot,wl(won),spec1d(won),color=fsc_color('dark green')
  paramsHb = MPFITFUN('twogaussHb',wl(won),spec1d(won),spec1derr(won),parinfo=pi,/nocatch, bestnorm=bestnorm, dof=dof, perror=perrorHb, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit)
  model = twogaussHb(wl(won1),paramsHb)
  oplot,wl(won1),model,color=fsc_color('red')
  z = paramsHb[1]
  vdisp = paramsHb[2]
  sigmaHb = Hb*(z+1.)*vdisp/clight;micron
  sigmaHberr = sigmaHb*sqrt((perrorHb[2]/vdisp)^2+(perrorHb[1]/z)^2)
  pixres = wl(1)-wl(0) ; micron/pixel
  Hbflux = sqrt(2.*!pi)*paramsHb[3]*sigmaHb*calibHb.conv_factor/pixres ;erg/s/cm^2
 ;The integration integrates ADU/s/keck_area/pixel over micron so have to divide by micron/pixel
  Hbfluxerr = Hbflux*sqrt((sigmaHberr/sigmaHb)^2+(perrorHb[3]/paramsHb[3])^2+(calibHb.conv_factorerr/calibha.conv_factor)^2)
  print, 'Hb Flux '+string(Hbflux,format='(E10.3)')+'+-'+string(Hbfluxerr,format='(E10.3)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Calculate E(B-V)
EBV    = 1.965*alog10(Haflux/Hbflux/2.86)
EBVerr = 1.965*sqrt((Hafluxerr/Haflux)^2+(Hbfluxerr/Hbflux)^2)/alog(10.)/2.86
AHa    = 3.33*EBV
AHaerr = 3.33*EBVerr
print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
print,'Extinction in Ha ='+string(AHa,format='(F5.2)')+'+-'+string(AHaerr,format='(F5.2)')+' mag'
print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
;print, params,perror,paramshb,perrorhb
if keyword_Set(stoptag) then stop
wait,1
return,[AHa,AHaerr]
end
