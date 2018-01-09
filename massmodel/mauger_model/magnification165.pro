pro magnification165
  im=readfits('cswa165_SDSSJ0105+0144_deflections.fits')
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
  writefits,'cswa165_magnification_idl.fits',magnification_ALL
  lengmag = size(magnification_All,/dimensions)

  imref=readfits('/scr2/nichal/workspace/output/cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr_acube.fits') 
  detection=imref[*,*,4] 
  detection = detection[0:lengmag(0)-1,0:lengmag(1)-1]
  good=where(finite(detection) and detection ge 5. and finite(magnification_all) and magnification_all ge 1.) 
  mag=mean(magnification_ALL(good))
  print, mag

  ;flux weighted
  Haflux = im[*,*,0]
  good = where(finite(Haflux) and finite(magnification_all) and magnification_all ge 1.)
  mag = total(Haflux(good))/total(Haflux(good)/magnification_all(good))
  print, mag
end
