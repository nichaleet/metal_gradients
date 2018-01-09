pro magnification128
  im=readfits('cswa128_SDSSJ1956+5950_deflections.fits')
  sx=im[*,*,1]                                          
  sy=im[*,*,2]
  leng=size(sx,/dimensions)
  ix = rebin(findgen(leng(0)),leng(0),leng(1))
  iy = rebin(transpose(findgen(leng(1))),leng(0),leng(1))
  ax = ix-sx
  ay = iy-sy
  da_x_dy = pdiv(ax,2)
  da_x_dx = pdiv(ax,1)
  da_y_dy = pdiv(ay,2)
  da_y_dx = pdiv(ay,1)                                
  poisson_ALL=da_x_dx+da_y_dy                                              
  magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx*da_y_dy-da_x_dy*da_y_dx)) 
  writefits,'cswa128_magnification_idl.fits',magnification_ALL
  imref=readfits('/scr2/nichal/workspace/output/cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_acube.fits') ;
  detection=imref[*,*,4]        ;
  good1=where(finite(detection) and detection ge 5. and magnification_all gt 1.) ;
  mag=mean(magnification_ALL(good1))
  print,'By mag map:', mag
  
  ;flux weighted
  Haflux = imref[*,*,1]
  good2 = where(finite(Haflux) and finite(magnification_all) and magnification_all lt 1.e6 and magnification_all gt 1.)
  mag = total(Haflux(good2))/total(Haflux(good2)/magnification_all(good2))
  print,'By mag mag weighted by Ha flux:', mag
  
  ;total Ha luminosity
  tot_ha = total(Haflux(good2))

  stop
end
