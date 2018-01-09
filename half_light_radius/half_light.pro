function half_light,haim,xc,yc,inc,theta,pixscale
  inc = inc*!pi/180.
  theta = theta*!pi/180.
  haim = readfits(haim)
  sizeim = size(haim,/dimensions)
  x = rebin(findgen(sizeim(0)),sizeim(0),sizeim(1))
  y = rebin(transpose(findgen(sizeim(1))),sizeim(0),sizeim(1))
  Rsquare = ((x-xc)*cos(theta)-(y-yc)*sin(theta))^2+(((x-xc)*sin(theta)+(y-yc)*cos(theta))/sin(inc))^2
  R = sqrt(Rsquare)
  writefits,'r.fits',r
  haim(where(~finite(haim))) = 0.
  tot_light = total(haim)
  
  light_ratio = 0.
  npix =1
  while light_ratio lt 0.5 do begin
     light_ratio_old = light_ratio
     lightin = total(haim(where(R lt npix)))
     light_ratio = lightin/tot_light
     npix = npix+1.
  endwhile
  ;interpolate
  fraction = (0.5-light_ratio_old)/(light_ratio-light_ratio_old)
  final_npix = fraction+npix-2.
  half_light_radius = final_npix*pixscale
  print,half_light_radius

  return,half_light_radius
end
