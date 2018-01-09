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

function calculate_SFR,filehacube,filehadetect,calibHa,z,mu,filter,AHa,AHaerr
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  Mpctocm = 3.085678e+24
  !p.multi=[0,1,1]
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
  wl = getwl_filter(filter)
  Cgplot,wl,spec1d,color=fsc_color('black'),title='Ha fitting'
  oplot,wl,spec1d+spec1derr,color=fsc_color('rose')
  oplot,wl,spec1d-spec1derr,color=fsc_color('rose')
  ;FIT gaussian to Ha and NII line
  pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)
  pi.value = double([median(spec1d),z,200.0,max(spec1d),0.5])
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
  RegOff = [RegOff,where(abs(spec1d-median(spec1d)) gt 2.*stdev(spec1d))]
  mask(RegOff) = 0
  RegIn  = [where(wl gt Ha*(1.+z)*(1.-400./clight) and wl lt Ha*(1.+z)*(1.+400./clight)),where(wl gt NII*(1.+z)*(1.-400./clight) and wl lt OIII*(1.+z)*(1.+400./clight))]
  mask(RegIn) = 1
  mask(0) = 0
  won = where(mask)
  oplot,wl(won),spec1d(won),color=fsc_color('dark green')
  params = MPFITFUN('twogaussHa',wl(won),spec1d(won),spec1derr(won),parinfo=pi,/nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit)
  model = twogaussHa(wl,params)
  oplot,wl,model,color=fsc_color('red')
  z = params[1]
  vdisp = params[2]
  sigmaHa = Ha*(z+1.)*vdisp/clight
  sigmaHaerr = sigmaHa*sqrt((perror[2]/vdisp)^2+(perror[1]/z)^2)
  pixres = wl(1)-wl(0) ; micron/pixel
  Haflux = sqrt(2.*!pi)*params[3]*sigmaHa*calibHa.conv_factor/pixres ;erg/s/cm^2
 ;The integration integrates ADU/s/keck_area/pixel over micron so have to divide by micron/pixel
  Hafluxerr = Haflux*sqrt((sigmaHaerr/sigmaHa)^2+(perror[3]/params[3])^2+(calibHa.conv_factorerr/calibha.conv_factor)^2)

  ;correct for extinction and magnification
  HafluxCorrected = Haflux*10.^(Aha/2.5)/mu
  HafluxCorrectedErr = HafluxCorrected*sqrt((alog(10)*AHaerr/2.5)^2+(HafluxErr/Haflux)^2)

  ;calculate SFR
  dist = lumdist(z)*Mpctocm
  SFR = 7.9e-42*HafluxCorrected*4.*!dpi*dist^2
  SFRerr = 7.9e-42*HafluxCorrectedErr*4.*!dpi*dist^2
  print, 'Ha Flux '+string(HafluxCorrected,format='(E10.3)')+'+-'+string(HafluxCorrectederr,format='(E10.3)')
  print, 'SFR '+string(SFR,format='(F6.2)')+'+-'+string(SFRerr,format='(F6.2)')
  ;stop
  return,{SFR:SFR,SFRerr:SFRerr,Haflux:Hafluxcorrected,HafluxErr:HafluxCorrectederr}
end
