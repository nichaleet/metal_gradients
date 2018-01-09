function cswa19
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=2.0325
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Mar2013/skysubtracted/FS26'
  filestarha = path+'/s130304_a026_mosaic_Kn1_100.fits'
  filestarhb = path+'/s130304_a027001__mosaic_Hn1_100.fits'
  ;calculate_fluxconv,[filestarha,filestarhb],name='FS26',cluster='c19',filters=['Kn1','Hn1'],photmags=[7.93,8.127],mags_err=[.009,.01],z=z,linewl=[Ha,Hb],outthres=[5000.,500.],calibfactors=calibfactors,stoptag=0
  calculate_fluxconv_2mass,[filestarha,filestarhb],name='FS26',cluster='c19',filters=['Kn1','Hn1'],photmags=[7.99,8.148],mags_err=[.024,.036],z=z,linewl=[Ha,Hb],outthres=[5000.,500.],calibfactors=calibfactors,stoptag=0
  ;calculate flux conversion of Ha with tt star
  filettstar='/scr2/nichal/workspace/reduced_data/Mar2013/skysubtracted/cswa19tt/s130303_a011001__mosaic_Kn1_100.fits'
  calculate_fluxconv_2mass,[filettstar],name='cswa19_tt',cluster='c19',filters=['Kn1'],photmags=[12.502],mags_err=[.05],z=z,linewl=[Ha],outthres=[5.],calibfactors=calibfactortt,stoptag=0
  calibfactors=[calibfactors,calibfactortt]
  ;calculate reddening
  filehacube   = pathcube+'cswa19_Ha_tlc_Kn1_pipelinemosaic_scaledsky_3hr.fits'
  filehadetect = pathout+'cswa19_Ha_tlc_Kn1_pipelinemosaic_scaledsky_3hr_acube.fits'
  filehbcube   = pathcube+'cswa19_Hb_tlc_Hn1_handmosaic_sky_1hr.fits'
  filehbdetect = pathout+'cswa19_Hb_tlc_Hn1_handmosaic_sky_1hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[0],calibfactors[1],2.03,'Kn1','Hn1')
  save,calibfactors,filename='outfluxconv/cswa19_calibfactors.sav'
  ;calculate SFR
  AHaArr = [0.5,0.5]
  second = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,4.3,'Kn1',AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,4.3,'Kn1',AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa15
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=2.16
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/FS26'
  filestarha = path+'/s140223_a032001_tlc_Kn2_100.fits'
  ;not using mosaic because the 2nd frame is too faint. close to sunrise
  filestarhb = path+'/s140223_a030001__mosaic_Hbb_100.fits'
  filestarhb = path+'/s140223_a030002_tlc_Hbb_100.fits'
  ;calculate_fluxconv,[filestarha,filestarhb],name='FS26',cluster='c15',filters=['Kn2','Hbb'],photmags=[7.93,8.127],mags_err=[.009,.01],z=z,linewl=[Ha,Hb],outthres=[20.,1000.],calibfactors=calibfactors,stoptag=0
  calculate_fluxconv_2mass,[filestarha,filestarhb],name='FS26',cluster='c15',filters=['Kn2','Hbb'],photmags=[7.99,8.148],mags_err=[.024,.036],z=z,linewl=[Ha,Hb],outthres=[20.,1000.],calibfactors=calibfactors,stoptag=0
  ;calculate flux conversion of Ha with tt star
  filettstar='/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/cswa15tt/s140223_a018001__mosaic_Kn2_100.fits'
  calculate_fluxconv_2mass,[filettstar],name='cswa15_tt',cluster='c15',filters=['Kn2'],photmags=[14.193],mags_err=[0.05],z=z,linewl=[Ha,Hb],outthres=[0.25],calibfactors=calibfactortt,fixcenter=1,centerpos=[32,19],stoptag=0
  calibfactors=[calibfactors,calibfactortt]
  ;calculate reddening
  filehacube   = pathcube+'cswa15_Ha_tlc_Kn2_handmosaic_sky_230hr.fits'
  filehadetect = pathout+'cswa15_Ha_tlc_Kn2_handmosaic_sky_230hr_acube.fits'
  filehbcube   = pathcube+'cswa15_Hb_tlc_Hbb_pipelinemosaic_scaledsky_1hr.fits'
  filehbdetect = pathout+'cswa15_Hb_tlc_Hbb_pipelinemosaic_scaledsky_1hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[0],calibfactors[1],z,'Kn2','Hbb')
  save,calibfactors,filename='outfluxconv/cswa15_calibfactors.sav'
  ;calculate SFR
  ;AHaArr=[0.5,0.5]
  second = calculate_SFR(filehacube,filehadetect,calibfactors[0],z,7.6,'Kn2',AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[0],z,7.6,'Kn2',AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa11
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  COMMON EXTINCT, meanAHaArr 

  z=1.41
  mu=1.9
  ;calculate flux conversion of Ha with tt star
  filettstar_mar='/scr2/nichal/workspace/reduced_data/Mar2013/skysubtracted/cswa11tt/s130304_a008001__mosaic_Hn2_100.fits'
  filettstar_feb='/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/CSWA11tt/s140223_a005001__mosaic_Hn2_100.fits'
  calculate_fluxconv_2mass,[filettstar_mar,filettstar_feb],name='cswa11_tt',cluster='c11',filters=['Hn2','Hn2'],photmags=[14.201,14.201],mags_err=[0.042,0.042],z=z,linewl=[Ha,Ha],outthres=[2,2],calibfactors=calibfactortt
  calculate_fluxconv_2mass,[filettstar_feb],name='cswa11_tt',cluster='c11',filters=['Hn2'],photmags=[14.201],mags_err=[0.042],z=z,linewl=[Ha],outthres=[2],calibfactors=calibfactortt_pipe,usepipeline=1
  calibfactors=[calibfactortt,calibfactortt_pipe]
  save,calibfactors,filename='outfluxconv/cswa11_calibfactors.sav'
  ;calculate SFR
  ;AHaArr = [0.5,0.5]
  AHaArr = meanAHaArr
  print, 'Using Extinction '+AHaArr
  filehacube   = pathcube+'cswa11_Ha_tlc_Hn2_handmosaic_scaledsky_430hr.fits'
  filehadetect = pathout+'cswa11_Ha_tlc_Hn2_handmosaic_scaledsky_430hr_secondgal_acube.fits'
  second = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,'Hn2',AHaarr[0],AHaarr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,'Hn2',AHaarr[0],AHaarr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  ;stop
  return,SFRarr
end

function cswa20
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=1.43
  mu=14
  filterha='Hn2'
  filterhb='Jn1'
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/'
  filestarha = path+'FS27/s140224_a026001__mosaic_Hn2_100.fits'
  filestarhb = path+'FS26/s140223_a027001__mosaic_Jn1_100.fits'
  calculate_fluxconv_2mass,[filestarha],name='FS27',cluster='c20',filters=[filterha],photmags=[13.186],mags_err=[.028],z=z,linewl=[Ha],outthres=[10.],calibfactors=calibfactora,stoptag=0
  calculate_fluxconv_2mass,[filestarhb],name='FS26',cluster='c20',filters=[filterhb],photmags=[8.828],mags_err=[.019],z=z,linewl=[Hb],outthres=[250.],calibfactors=calibfactorb,stoptag=0
  calibfactorFS=[calibfactora,calibfactorb]
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/cswa20tt'
  filestarha = path+'/s140224_a020001__mosaic_Hn2_100.fits'
  ;not using mosaic because the 2nd frame is too faint. close to sunrise
  filestarhb = path+'/s140223_a024002_tlc_Jn1_100.fits'
  calculate_fluxconv_2mass,[filestarha,filestarhb],name='c20tt',cluster='c20',filters=[filterha,filterhb],photmags=[15.495,15.763],mags_err=[.128,.077],z=z,linewl=[Ha,Hb],outthres=[3.,0.5],calibfactors=calibfactortt,fixpeak=1,floor=[1,1],stoptag=0
  calibfactors=[calibfactortt[0],calibfactorFS[1]]
  ;calculate reddening
  filehacube   = pathcube+'cswa20_Ha_tlc_Hn2_handmosaic_scaledsky_1hr.fits'
  filehadetect = pathout+'cswa20_Ha_tlc_Hn2_handmosaic_scaledsky_1hr_acube.fits'
  filehbcube   = pathcube+'cswa20_Hb_tlc_Jn1_handmosaic_sky_075hr.fits'
  filehbdetect = pathout+'cswa20_Hb_tlc_Jn1_handmosaic_sky_075hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[0],calibfactors[1],z,filterha,filterhb)
  save,calibfactors,filename='outfluxconv/cswa20_calibfactors.sav'
  ;calculate SFR
  ;AHaArr=[0,0.5]
  second = calculate_SFR(filehacube,filehadetect,calibfactors[0],z,mu,filterha,AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[0],z,mu,filterha,AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa28
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  COMMON EXTINCT, meanAHaArr 
  z=2.09
  mu=9.3
  ;calculate flux conversion of Ha with tt star
  filettstar='/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/cswa28tt/s140223_a022001__mosaic_Kn1_100.fits'
  fileFSstar='/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/FS26/s140223_a031001__mosaic_Kn1_100.fits'
  calculate_fluxconv_2mass,[filettstar],name='cswa28_tt',cluster='c28',filters=['Kn1'],photmags=[12.397],mags_err=[0.023],z=z,linewl=[Ha],outthres=[30],calibfactors=calibfactortt,stoptag=0
  calculate_fluxconv_2mass,[filettstar],name='cswa28_tt_pipeline',cluster='c28',filters=['Kn1'],photmags=[12.397],mags_err=[0.023],z=z,linewl=[Ha],outthres=[30],calibfactors=calibfactortt_pipe,usepipeline=1,stoptag=0
  calculate_fluxconv_2mass,[fileFSstar],name='FS26',cluster='c28',filters=['Kn1'],photmags=[8.148],mags_err=[0.036],z=z,linewl=[Ha],outthres=[200],calibfactors=calibfactorFS,stoptag=0

  calibfactors = [calibfactortt,calibfactortt_pipe,calibfactorFS]
  save,calibfactors,filename='outfluxconv/cswa28_calibfactors.sav'
  ;calculate SFR
  filehacube   = pathcube+'cswa28_Ha_tlc_Kn1_handmosaic_scaledsky_130hr.fits'
  filehadetect = pathout+'cswa28_Ha_tlc_Kn1_handmosaic_scaledsky_130hr_acube.fits'
  ;main result
  AHaArr = [0.5,0.5]
  AHaArr =meanAHaArr
  main = calculate_SFR(filehacube,filehadetect,calibfactors[1],z,mu,'Kn1',AHaArr[0],AHaArr[1])
  second = calculate_SFR(filehacube,filehadetect,calibfactors[1],z,mu,'Kn1',AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa31
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  COMMON EXTINCT, meanAHaArr 

  z=1.486
  mu=3.3
  ;calculate flux conversion of Ha with tt star
  filettstar='/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/cswa31tt/s141206_a012001__mosaic_Hn3_100.fits'
  fileFSstar='/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/FS19/s141207_a028001__mosaic_Hn3_100.fits'
  calculate_fluxconv_2mass,[filettstar,filettstar],name='cswa31_tt',cluster='c31',filters=['Hn3','Hn3'],photmags=[13.071,12.853],mags_err=[0.033,0.07],z=z,linewl=[Ha,Ha],outthres=[30,30],calibfactors=calibfactortt,stoptag=0
  calculate_fluxconv_2mass,[filettstar,filettstar],name='cswa31_tt_pipe',cluster='c31',filters=['Hn3','Hn3'],photmags=[13.071,12.853],mags_err=[0.033,0.07],z=z,linewl=[Ha,Ha],outthres=[30,30],calibfactors=calibfactortt_pipe,usepipeline=1,stoptag=0
  calculate_fluxconv_2mass,[fileFSstar],name='FS19',cluster='c31',filter=['Hn3'],photmags=[13.695],mags_err=[0.022],z=z,linewl=[Ha],outthres=[1],calibfactors=calibfactorFS,stoptag=0,/fwtmtag
;The tt star is double stars
;The first conv is from the main star. The 2nd conv is from combining two stars.
  calibfactors = [calibfactortt,calibfactortt_pipe,calibfactorFS]
  save,calibfactors,filename='outfluxconv/cswa31_calibfactors.sav'
  ;calculate SFR
  filehacube   = pathcube+'cswa31_Ha_tlc_Hn3_handmosaic_scaledsky_3hr.fits'
  filehadetect = pathout+'cswa31_Ha_tlc_Hn3_handmosaic_scaledsky_3hr_acube.fits'
  ;main result
  AHa_arr=[0.5,0.5]
  AHa_Arr =meanAHaArr
  main   = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,'Hn3',AHa_arr[0],AHa_Arr[1])
  second = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,'Hn3',AHa_arr[0],AHa_Arr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]

  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa128
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=2.229
  mu = 10
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Sep2013/skysubtracted/FS26'
  filestarha = path+'/s130911_a016001__mosaic_Kn2_100.fits'
  filestarhb = path+'/s130911_a019001__mosaic_Hbb_100.fits'
  calculate_fluxconv_2mass,[filestarha,filestarhb],name='FS26',cluster='c128',filters=['Kn2','Hbb'],photmags=[7.99,8.148],mags_err=[.024,.036],z=z,linewl=[Ha,Hb],outthres=[1500.,1000.],calibfactors=calibfactora,stoptag=0,/fwhmtag
  calculate_fluxconv_2mass,[filestarha,filestarhb],name='FS26',cluster='c128',filters=['Kn2','Hbb'],photmags=[7.99,8.148],mags_err=[.024,.036],z=z,linewl=[Ha,Hb],outthres=[1500.,1000.],calibfactors=calibfactorb,stoptag=0,/fwtmtag

  ;calculate flux conversion of Ha with tt star
  filettstar='/scr2/nichal/workspace/reduced_data/Sep2013/ttstars/s130911_a020001__mosaic_Kn2_100.fits'
  calculate_fluxconv_2mass,[filettstar],name='cswa128_tt',cluster='c128',filters=['Kn2'],photmags=[14.493],mags_err=[0.072],z=z,linewl=[Ha],outthres=[2],calibfactors=calibfactortt,stoptag=0,/fwtmtag
  calibfactors=[calibfactora[0],calibfactorb[1],calibfactortt]
  ;calculate reddening
  filehacube   = pathcube+'cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr.fits'
  filehadetect = pathout+'cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube.fits'
  filehbcube   = pathcube+'cswa128_Hb_tlc_Hbb_pipelinemosaic_scaledsky_130hr.fits'
  filehbdetect = pathout+'cswa128_Hb_tlc_Hbb_pipelinemosaic_scaledsky_130hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[0],calibfactors[1],z,'Kn2','Hbb',/stoptag)
  save,calibfactors,filename='outfluxconv/cswa128_calibfactors.sav'
  ;calculate SFR
  ;AHaArr = [0,0.5]
  second = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,'Kn2',AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,'Kn2',AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa139
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=2.54
  filterha='Kc5'
  filterhb='Hbb'
  mu = 9.7
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/FS27'
  filestarha = path+'/s140224_a029001__mosaic_Kc5_100.fits'
  path = '/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/FS19'
  filestarhb = path+'/s141207_a026001__mosaic_Hbb_100.fits'
  calculate_fluxconv_2mass,[filestarha,filestarhb],name='FS27,FS19',cluster='c139',filters=[filterha,filterhb],photmags=[13.140,13.695],mags_err=[.035,.022],z=z,linewl=[Ha,Hb],outthres=[5,5.],calibfactors=calibfactors,stoptag=0,/fwtmtag
  ;calculate flux conversion of Ha with tt star
  pathtt='/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/cswa139tt/'
  filettstarHa = pathtt+'s141207_a015001__mosaic_Kc5_100.fits'
  filettstarHb = pathtt+'s141207_a019001__mosaic_Hbb_100.fits'
  calculate_fluxconv_2mass,[filettstarHa,filettstarHb],name='cswa139_tt',cluster='c139',filters=[filterha,filterhb],photmags=[14.892,14.872],mags_err=[0.103,0.059],z=z,linewl=[Ha,Hb],outthres=[1.,1],calibfactors=calibfactortt,stoptag=0,/fwhmtag
  calibfactors=[calibfactors,calibfactortt]
  ;calculate reddening
  filehacube   = pathcube+'cswa139_Ha_tlc_Kc5_handmosaic_sky_3hr.fits'
  filehadetect = pathout+'cswa139_Ha_tlc_Kc5_handmosaic_sky_3hr_acube.fits'
  filehbcube   = pathcube+'cswa139_Hb_tlc_Hbb_pipelinemosaic_scaledsky_dec_130hr.fits'
  filehbdetect = pathout+'cswa139_Hb_tlc_Hbb_pipelinemosaic_scaledsky_dec_130hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[2],calibfactors[3],z,filterha,filterhb)
  save,calibfactors,filename='outfluxconv/cswa139_calibfactors.sav'
  ;calculate SFR
  ;AHaArr=[0.5,0.5]
  second = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,filterha,AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,filterha,AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa159
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=2.30
  filterha='Kc3'
  filterhb='Hn3'
  mu = 4.6
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Sep2013/skysubtracted/FS9'
  filestarha = path+'/s130911_a033001__mosaic_Kc3_100.fits'
  path = '/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/FS19'
  filestarhb = path+'/s141207_a028001__mosaic_Hn3_100.fits'
  calculate_fluxconv_2mass,[filestarhb],name='FS19',cluster='c159',filters=[filterhb],photmags=[13.695],mags_err=[.022],z=z,linewl=[Hb],outthres=[5.],calibfactors=calibfactorHb,stoptag=0
  calculate_fluxconv_2mass,[filestarha],name='FS9',cluster='c159',filters=[filterha],photmags=[8.223],mags_err=[.021],z=z,linewl=[Ha],outthres=[300],calibfactors=calibfactorHa,stoptag=0
  
  ;calculate flux conversion of Ha with tt star
  filettstar='/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/cswa159tt/s141207_a005001__mosaic_Kc3_100.fits'
  calculate_fluxconv_2mass,[filettstar],name='cswa159_tt',cluster='c159',filters=[filterha],photmags=[12.293],mags_err=[0.023],z=z,linewl=[Ha],outthres=[20.],calibfactors=calibfactortt,stoptag=0
  calibfactors=[calibfactorHa,CalibfactorHb,calibfactortt]
  ;calculate reddening
  filehacube   = pathcube+'cswa159_Ha_tlc_Kc3_pipelinemosaic_sky_2hr.fits'
  filehadetect = pathout+'cswa159_Ha_tlc_Kc3_pipelinemosaic_sky_2hr_acube.fits'
  filehbcube   = pathcube+'cswa159_Hb_tlc_Hn3_handmosaic_scaledsky_1hr.fits'
  filehbdetect = pathout+'cswa159_Hb_tlc_Hn3_handmosaic_scaledsky_1hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[2],calibfactors[1],z,filterha,filterhb,/stoptag)
  save,calibfactors,filename='outfluxconv/cswa139_calibfactors.sav'
  AHaArr[1]=AHaArr[1]*2. ;the Hb line is on the skyline
  ;calculate SFR
  AHaArr=[0.5,0.5]
  print, AHaArr
  second = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,filterha,AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,filterha,AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function cswa165
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=2.13
  filterha='Kn2'
  filterhb='Hbb'
  mu = 42
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Sep2013/skysubtracted/FS9'
  filestarha1 = path+'/s130911_a032001__mosaic_Kn2_100.fits'
  calculate_fluxconv_2mass,[filestarha1],name='FS9',cluster='c165',filters=[filterha],photmags=[8.223],mags_err=[.021],z=z,linewl=[Ha],outthres=[1000],calibfactors=calibfactorHa1,stoptag=0
  path = '/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/FS19'
  filestarha2 = path+'/s141207_a029001__mosaic_Kn2_100.fits'
  filestarhb = path+'/s141207_a026001__mosaic_Hbb_100.fits'
  calculate_fluxconv_2mass,[filestarha2,filestarhb],name='FS19',cluster='c139',filters=[filterha,filterhb],photmags=[13.833,13.695],mags_err=[.047,.022],z=z,linewl=[Ha,Hb],outthres=[1,5.],calibfactors=calibfactors,stoptag=0
  ;calculate flux conversion of Ha with tt star
  pathtt='/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/cswa165tt/'
  filettstarHa = pathtt+'s141207_a010_mosaic_Kn2_100.fits'
  filettstarHb = pathtt+'s141207_a013_mosaic_Hbb_100.fits'
  calculate_fluxconv_2mass,[filettstarHa,filettstarHb],name='cswa165_tt',cluster='c165',filters=[filterha,filterhb],photmags=[13.471,13.451],mags_err=[0.032,0.025],z=z,linewl=[Ha,Hb],outthres=[10.,5],calibfactors=calibfactortt,stoptag=0
  calibfactors=[calibfactorHa1,calibfactors,calibfactortt]
  ;calculate reddening
  filehacube   = pathcube+'cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr.fits'
  filehadetect = pathout+'cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr_acube.fits'
  filehbcube   = pathcube+'cswa165_Hb_tlc_Hbb_pipelinemosaic_scaledsky_130hr.fits'
  filehbdetect = pathout+'cswa165_Hb_tlc_Hbb_pipelinemosaic_scaledsky_130hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[1],calibfactors[2],z,filterha,filterhb,/stoptag)
  save,calibfactors,filename='outfluxconv/cswa139_calibfactors.sav'
  ;calculate SFR
  ;AHaArr=[0.5,0.5]
  second = calculate_SFR(filehacube,filehadetect,calibfactors[3],z,mu,filterha,AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[3],z,mu,filterha,AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

function abell773
  COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
  COMMON PATH, pathout,pathcube
  z=2.30
  filterha='Kc3'
  filterhb='Hn3'
  mu = 20.3
  ;flux conversion with standard stars
  path = '/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/FS26'
  filestarha = path+'/s140222_a013_mosaic_Kc3_100.fits'
  path = '/scr2/nichal/workspace/reduced_data/Mar2013/skysubtracted/FS11'
  filestarhb = path+'/s130304_a004_mosaic_Hn3_100.fits'
  calculate_fluxconv_2mass,[filestarha,filestarhb],name='FS26,FS11',cluster='a773',filters=[filterha,filterhb],photmags=[7.99,11.311],mags_err=[.024,.024],z=z,linewl=[Ha,Hb],outthres=[300,10.],calibfactors=calibfactors,stoptag=0
  ;calculate flux conversion of Ha with tt star
  pathtt='/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/a773tt/'
  filettstarHa = pathtt+'s140222_a005001__mosaic_Kc3_100.fits'
  pathtt='/scr2/nichal/workspace/reduced_data/Mar2013/skysubtracted/a773tt/'
  filettstarHb = pathtt+'s130303_a005_mosaic_Hn3_100.fits'
  calculate_fluxconv_2mass,[filettstarHa,filettstarHb],name='a773tt',cluster='c139',filters=[filterha,filterhb],photmags=[10.668,10.735],mags_err=[0.019,0.028],z=z,linewl=[Ha,Hb],outthres=[15.,15],calibfactors=calibfactortt,stoptag=0
  calibfactors=[calibfactors,calibfactortt]
  ;calculate reddening
  filehacube   = pathcube+'abell773_Ha_tlc_Kc3_handmosaic_sky_330hr.fits'
  filehadetect = pathout+'abell773_Ha_tlc_Kc3_handmosaic_sky_330hr_acube.fits'
  filehbcube   = pathcube+'abell773_Hb_tlc_Hn3_pipelinemosaic_sky_2hr.fits'
  filehbdetect = pathout+'abell773_Hb_tlc_Hn3_pipelinemosaic_sky_2hr_acube.fits'
  AHaArr = calculate_reddening(filehacube,filehadetect,filehbcube,filehbdetect,calibfactors[2],calibfactors[3],z,filterha,filterhb)
  save,calibfactors,filename='outfluxconv/cswa139_calibfactors.sav'
  ;calculate SFR
  AHaArr=[0.5,0.5]
  second = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,filterha,AHaArr[0],AHaArr[1])
  ;main result
  main = calculate_SFR(filehacube,filehadetect,calibfactors[2],z,mu,filterha,AHaArr[0],AHaArr[1])
  if main.sfr gt second.sfr then begin
     high=main
     low =second
  endif else begin
     high = second
     low  = main
  endelse 
  SFRarr = [low.sfr-low.sfrerr,main.sfr,high.sfr+high.sfrerr]
  Hafluxarr = [low.haflux-low.hafluxerr,main.haflux,high.haflux+high.hafluxerr]
  print, 'FINAL SFR:',SFRarr
  print, 'FINAL Haflux:',Hafluxarr
  return,SFRarr
end

pro runcalculate_reddening
COMMON LINEWAVE, Ha,Hb,NII,OIII,clight
COMMON PATH, pathout,pathcube
COMMON EXTINCT, meanAHaArr 
  Ha = .656461
  Hb = .486269
  NII= .658523
  OIII= .500824
  clight = 299792.458
  ;find average AHA
  AHaAll    = [0.52,0.46,0.27,1.96,0.38,0.44,1.03,2.5,0.53,0.67,2.94,0.3]
  AHaAllErr = [0.1,0.17,0.34,0.6,0.52,0.95,0.5,0.25,1.2,1.3,1.05,0.12]
  meanerr,AHaAll,AHaAllErr,meanAHa,sigmam,meanAHaErr,sigmas
  meanAHaArr = [meanAHa,meanAHaErr]
  pathout = '/scr2/nichal/workspace/output/'
  pathcube= '/scr2/nichal/workspace/reduced_data/mosaic/'
  ;SFR19=cswa19()
  ;SFR15=cswa15()
  ;SFR11=cswa11()
  ;SFR20=cswa20()
  ;SFR28 =cswa28()
  ;SFR31 = cswa31()
  SFR128 = cswa128()
  ;SFR139 = cswa139()
  ;SFR159 = cswa159()
  ;SFR165 = cswa165()
  ;SFR773 = abell773()
end
