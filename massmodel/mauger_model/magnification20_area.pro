pro magnification20_area
  imcube=readfits('/scr2/nichal/workspace/output/cswa20_Ha_tlc_Hn2_handmosaic_scaledsky_1hr_acube.fits')
  im = imcube[*,*,4]
  im_detect = where(im gt 7. and finite(im),areaim)
  areaim = areaim*0.1^2

  sourcecube = readfits('/scr2/nichal/workspace/output/cswa20_Ha_tlc_Hn2_handmosaic_scaledsky_1hr_acube_sourceplane_interp.fits')
  source = sourcecube[*,*,4]
  source_detect = where(source gt 7. and finite(source),areasource)
  areasource = areasource*0.02^2

  mag = areaim/areasource
  print,'By area:', mag

  ;flux weighted
  sourceha = sourcecube[*,*,1]
  flux_source = total(sourceha(where(finite(sourceha) and source gt 7.)))*0.02^2
  imha = imcube[*,*,1]
  flux_im = total(imha(where(finite(imha) and im gt 7.)))*0.1^2
  mag = flux_im/flux_source
  print,'By flux:', mag
  stop
end
