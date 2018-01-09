pro metalgradient_oct26,metalmap,metalerrmap,angle,xmid,ymid,x,y,n_slits,slit_size,slitwidth,radius_metal,metal,metal_err,linparam
  set_plot,'x'
  ;major axis
  ;equation for major axis is y=mx+c st. m=slope_major
  m_major = tan(angle/180.*!pi)  
  c_major = ymid-m_major*xmid
  ;distance of a point(x,y) to this major axis is |-mx+y-c|/sqrt(m^2+1)
  dist = abs((-1.*m_major*x)+y-c_major)/sqrt(m_major^2+1.)
  
  ;Next are the angle and equations for the slit dividers 
  angle_perpend = angle-90.  ;angle of the little slit
  slope = tan(angle_perpend/180.*!pi)  

  ;ref points for each little slit dividers
  xref = fltarr(n_slits+1)+xmid  
  yref_up = ymid+findgen(n_slits+1)*slitwidth; for the slits above the ref points
  yref_down = ymid-findgen(n_slits+1)*slitwidth; for the slits below the ref points
; y = slope*x+intercept
; intercept = y-slope*x
  intercept_up = yref_up-slope(0)*xref
  intercept_down = yref_down-slope(0)*xref
 
  metal = fltarr(n_slits)
  metal_err = fltarr(n_slits)

  for i=0,n_slits-1 do begin
     pix_in_slit_up = where(y le slope*x+intercept_up(i+1) and y ge slope*x+intercept_up(i) and dist le 4.) ; The slit width is 8 pixels    
     pix_in_slit_down = where(y ge slope*x+intercept_down(i+1) and y le slope*x+intercept_down(i) and dist le 4.) ; The slit width is 8 pixels 
     pix_in_slit = SetUnion(pix_in_slit_up,pix_in_slit_up)
     good = metalmap(pix_in_slit)
     gooderr = metalerrmap(pix_in_slit)
     if pix_in_slit(0) ne -1 then begin 
        print,'There are', n_Elements(where(finite(good) eq 1)),'pixels with velocity in this slit'
        metalin = good(where(finite(good)))
        metalerrin = gooderr(where(finite(good)))
        meanerr,metalin,metalerrin,meanmetal,sigmam,sigmad
        metal(i) = mean(metalin)
        metal_err(i) = sqrt(sigmam^2*1./sqrt(n_Elements(metalin)-1.)*total((metalin-meanmetal)^2/metalerrin^2))
       ; metal_err(i)=sigmad/sqrt(float(n_elements(metalin)))
        if mean(metalin) eq 7.0 then metal(i) = 0.
        ;if n_elements(metalin) le 3. then metal(i) = 0.
        ;stop
     endif else print,'There is no pixel in this slit.'
  endfor
  radius_metal = findgen(n_slits)*slit_size
  baddata = where(metal eq 0. or finite(metal) eq 0)
  if baddata(0) ne -1 then remove,baddata, radius_metal,metal,metal_err
  ploterror,radius_metal,metal,metal_err,xtitle='radius from center(kpc)',ytitle='log(O/H)+12',psym=1  
  metal_err(where(metal_err lt 0.001)) = mean(metal_err)
  linparam = linfit(radius_metal,metal,chisq=chisq,measure_errors=metal_err,sigma=sigma)
  oplot,!x.crange,!x.crange*linparam(1)+linparam(0),color=60
  print, 'metal gradient =',linparam(1),'+/-',sigma(1)
;stop
end
