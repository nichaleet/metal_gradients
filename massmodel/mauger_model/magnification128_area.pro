pro magnification128_area
  imcube=readfits('/scr2/nichal/workspace/output/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube.fits') ;
  im = imcube[*,*,4]
  im_detect = where(im gt 5. and finite(im),areaim)
  areaim = areaim*0.1^2

  sourcecube = readfits('/scr2/nichal/workspace/output/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube_sourceplane_interp.fits');
  source = sourcecube[*,*,4]
  source_detect = where(source gt 5. and finite(source),areasource)
  areasource = areasource*0.04^2

  mag = areaim/areasource
  print,'By area:', mag

  ;flux weighted
  sourceha = sourcecube[*,*,1]
  flux_source = total(sourceha(where(finite(sourceha))))*0.04^2
  imha = imcube[*,*,1]
  flux_im = total(imha(where(finite(imha))))*0.1^2
  mag = flux_im/flux_source
  print,'By flux:', mag
  stop
end
