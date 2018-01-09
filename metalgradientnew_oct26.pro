pro funcPP04,x,A,f,pder

gradient = A[0]
Zc = A[1]

f = 10.^((gradient*x+Zc-8.90)/0.57)

pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro metalgradientnew_oct26,angle,metalmap,metalmaperr,xmid,ymid,x,y,n_slits,slit_size,slitwidth,radius_N2,N2,N2_err,A
;This one selects bins from slit in a similar manner to the rotation curve. 
  ;The inputmetalmap is map of N2 = NII/Ha
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
 
  ;output arrays
  N2 = fltarr(n_slits)
  N2_err = fltarr(n_slits)


  for i=0,n_slits-1 do begin
     pix_in_slit_up = where(y le slope*x+intercept_up(i+1) and y ge slope*x+intercept_up(i) and dist le 4.) ; The slit width is 8 pixels    
     pix_in_slit_down = where(y ge slope*x+intercept_down(i+1) and y le slope*x+intercept_down(i) and dist le 4.) ; The slit width is 8 pixels 
     pix_in_slit = SetUnion(pix_in_slit_up,pix_in_slit_up)
     
     good = metalmap(pix_in_slit)
     gooderr = metalmaperr(pix_in_slit)
     if pix_in_slit(0) ne -1 then begin 
        print,'There are', n_Elements(where(finite(good) eq 1)),'pixels with velocity in this slit'
        N2in = good(where(finite(good) eq 1))
        N2errin = gooderr(where(finite(good) eq 1))
        meanerr,N2in,N2errin,meanN2,sigmam,sigmad
        if meann2 gt 0.5 then meann2 = 0.
        if meann2 lt -0.2 then meann2 = 0.
        N2(i) = meanN2
        N2_err(i) = sigmam;sqrt(sigmam^2*1./sqrt(n_Elements(N2in)-1.)*total((N2in-meanN2)^2/N2errin^2))
;sigmad/sqrt(float(n_elements(N2in)))
        
     endif else print,'There is no pixel in this slit.'
  endfor
  radius_N2 = findgen(n_slits)*slit_size
  baddata = where(N2 eq 0.)
  if baddata(0) ne -1 then remove,baddata, radius_N2,N2,N2_err
  ploterror,radius_N2,N2,N2_err,xtitle='radius from center(kpc)',ytitle='[NII]/Ha',psym=1  

  ;Fit the function
  A = [0.,8.9] ;guess parameter(gradient and Z center)
  weight = 1./N2_err^2
  x=radius_N2
  y=N2
  N2fit = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04',/noderivative)  
  print, 'metal gradient =',A(0),'+/-',sigmaA(0)
  print, 'central metal =', A(1),'+/-',sigmaA(1)
;stop
end
