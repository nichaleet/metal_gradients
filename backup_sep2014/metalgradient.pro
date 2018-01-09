pro metalgradient,metalmap,xmid,ymid,x,y,n_slits,slit_size,slitwidth,radius_metal,metal,metal_err,linparam
  set_plot,'x'
  metal = fltarr(n_slits)
  metal_err = fltarr(n_slits)
  radius = sqrt((x-xmid)^2+(y-ymid)^2)
  for i=0,n_slits-1 do begin
     radius_in = i*slitwidth
     radius_out = (i+1.)*slitwidth
     pix_in_slit = where(radius ge radius_in and radius le radius_out)
     good = metalmap(pix_in_slit)
     if pix_in_slit(0) ne -1 then begin 
        print,'There are', n_Elements(where(finite(good) eq 1)),'pixels with velocity in this slit'
        metal(i) = mean(good(where(finite(good) eq 1)))
        metal_err(i) = stddev(good(where(finite(good) eq 1)))
     endif else print,'There is no pixel in this slit.'
  endfor
  radius_metal = findgen(n_slits)*slit_size
  baddata = where(metal eq 0.)
  if baddata(0) ne -1 then remove,badata, radius_metal,metal,metal_err
  ploterror,radius_metal,metal,metal_err,xtitle='radius from center(kpc)',ytitle='log(O/H)+12',psym=1  
  linparam = linfit(radius_metal,metal,chisq=chisq,measure_errors=metal_err,sigma=sigma)
  oplot,!x.crange,!x.crange*linparam(1)+linparam(0),color=60
  print, 'metal gradient =',linparam(1),'+/-',sigma(1)
;stop
end
