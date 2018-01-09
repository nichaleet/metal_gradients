pro run_half_light
  namearr=[]
  radiusarr=[]

  radius=half_light('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa11/outputmaps_secondgal/prettymaps/crop/sourceha_pretty.fits',15,12,64,16,0.626)
  namearr=[namearr,'cswa11']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa15/LTMzoomESI/outputmaps/prettymaps/crop/10sigmadetection/sourceha_pretty.fits',26,13,26,194,0.422)
  namearr=[namearr,'cswa15']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMredshift2fixedtry2/outputmaps_handmosaic/prettymaps/crop/maingalaxy/sourceha_pretty.fits',44,43,49,156,0.167)
  namearr=[namearr,'cswa19']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/output/chopped_sourceplane/cswa20/cswa20_sourceha.fits',15,20,40,162,0.172)
  namearr=[namearr,'cswa20']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMredshift2fixed_cswa28/outputmaps_handmosaic/prettymaps/crop/5sigmadetection/sourceha_pretty.fits',6,12,23,129,0.424)
  namearr=[namearr,'cswa28']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa31/outputmaps/prettymaps/crop/sourceha_pretty.fits',13.5,21,25,150,0.628)
  namearr=[namearr,'cswa31']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/output/chopped_sourceplane/cswa128/cswa128_sourceha.fits',13,16,57.,150.,0.336)
  namearr=[namearr,'cswa128']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/cswa139zoomESI/outputmaps/prettymaps/crop/sourceha_pretty.fits',19,15,24,297,0.40975)
  namearr=[namearr,'cswa139']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa159/outputmaps/prettymaps/crop/sourceha_pretty.fits',13,17,13,175,0.418205)
  namearr=[namearr,'cswa159']
  radiusarr=[radiusarr,radius]
  
  radius=half_light('/scr2/nichal/workspace/output/chopped_sourceplane/cswa165/cswa165_sourceha.fits',14,12.5,57.,218.,0.106)
  namearr=[namearr,'cswa165']
  radiusarr=[radiusarr,radius]

  radius=half_light('/scr2/nichal/workspace/output/chopped_sourceplane/Abell773/Abell773_sourceha.fits',20,23,22,144,0.159)
  namearr=[namearr,'Abell773']
  radiusarr=[radiusarr,radius]
  stop
end
