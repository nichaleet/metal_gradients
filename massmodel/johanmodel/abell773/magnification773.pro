pro magnification773
  im=readfits('a773_model.fits')
  ax=im[*,*,1]*10. ;change arcsec to 0.1 arcsec    
  ay=im[*,*,2]*10.
  da_x_dy = pdiv(ax,2)
  da_x_dx = pdiv(ax,1)
  da_y_dy = pdiv(ay,2)
  da_y_dx = pdiv(ay,1)                                
  poisson_ALL=da_x_dx+da_y_dy                                              
  magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx*da_y_dy-da_x_dy*da_y_dx)) 
  writefits,'a773_magnification_idl.fits',magnification_ALL
  imref=readfits('/scr2/nichal/workspace/output/abell773_Ha_tlc_Kc3_handmosaic_sky_330hr_acube.fits')
  detection=imref[*,*,4]        ;
  good1=where(finite(detection) and detection ge 5. and magnification_all gt 1.) ;
  mag=mean(magnification_ALL(good1))
  print,'By mag map:', mag
  
  ;flux weighted
  Haflux = imref[*,*,1]
  good2 = where(finite(Haflux) and finite(magnification_all) and magnification_all lt 1.e6 and magnification_all gt 1.)
  mag = total(Haflux(good2))/total(Haflux(good2)/magnification_all(good2))
  print,'By mag mag weighted by Ha flux:', mag
  
  stop
end
