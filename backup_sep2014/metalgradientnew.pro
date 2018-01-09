pro funcPP04,x,A,f,pder

gradient = A[0]
Zc = A[1]

f = 10.^((gradient*x+Zc-8.90)/0.57)

pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro metalgradientnew,metalmap,metalmaperr,xmid,ymid,x,y,n_slits,slit_size,slitwidth,radius_N2,N2,N2_err,A

  ;The inputmetalmap is map of N2 = NII/Ha
  set_plot,'x'
  N2 = fltarr(n_slits)
  N2_err = fltarr(n_slits)
  radius = sqrt((x-xmid)^2+(y-ymid)^2)
  for i=0,n_slits-1 do begin
     radius_in = i*slitwidth
     radius_out = (i+1.)*slitwidth
     pix_in_slit = where(radius ge radius_in and radius le radius_out)
     good = metalmap(pix_in_slit)
     gooderr = metalmaperr(pix_in_slit)
     if pix_in_slit(0) ne -1 then begin 
        print,'There are', n_Elements(where(finite(good) eq 1)),'pixels with velocity in this slit'
        N2in = good(where(finite(good) eq 1))
        N2errin = gooderr(where(finite(good) eq 1))
        meanerr,N2in,N2errin,meanN2,sigmam,sigmad
        N2(i) = meanN2
        N2_err(i) = sigmam
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
;stop
end
