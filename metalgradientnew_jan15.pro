pro funcPP04,x,A,f,pder

gradient = A[0]
Zc = A[1]

f = 10.^((gradient*x+Zc-8.90)/0.57)

pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end


pro metalgradientnew_jan15,angle,xref,yref,xmid,ymid,x,y,n_slits,midpoint,N2map,N2maperr,distance,N2_distance,N2_distance_err,radius,N2_radius,N2_radius_err,slit_size,name,A,sigmaA

    ;PP04 method to get the slope of the metallicity from N2 index
    ;This version only take the pixels in a slit along a major axis 
    ;to calculate.

  slitmap = N2map ;this map will show where the slits are
  ;major axis
  ;equation for major axis is y=mx+c st. m=slope_major
  m_major = tan(angle/180.*!pi)  
  c_major = ymid-m_major*xmid
  ;distance of a point(x,y) to this major axis is |-mx+y-c|/sqrt(m^2+1)
  dist = abs((-1.*m_major*x)+y-c_major)/sqrt(m_major^2+1.)
  
  ;Next are the angle and equations for the slit dividers 
  angle_perpend = angle-90.  ;angle of the little slit
  slope = tan(angle_perpend/180.*!pi)  

; y = slope*x+intercept
; intercept = y-slope*x
  intercept = yref-slope(0)*xref

  N2_distance = fltarr(n_slits)
  N2_distance_err = fltarr(n_slits)
  for i=0,n_slits-1 do begin
     pix_in_slit = where(y le slope*x+intercept(i+1) and y ge slope*x+intercept(i) and dist le 1.5) ; The slit width is ~3 pixels     
     if pix_in_slit(0) ne -1 then begin
        slitmap(pix_in_slit)=i*100.
        goodN2 = N2map(pix_in_slit) 
        goodN2err = N2maperr(pix_in_slit)
        
        thegood = where(finite(goodN2) eq 1 and goodN2 gt 0.015)
        print,'There are', n_Elements(thegood),'pixels with velocity in this slit'
        simplemean = mean(goodN2(thegood))
        meanerr,goodN2(thegood),goodN2err(thegood),wmean,sigmam,sigmad
        print,simplemean,wmean,' N2index'

        if finite(wmean) eq 1 then begin 
           N2_distance(i) = wmean 
           N2_distance_err(i) = sigmad
        endif else begin
           N2_distance(i)=simplemean
            N2_distance_err(i) = stddev(goodN2(thegood))
         endelse         
        if N2_distance(i) le 0.05 then stop
     endif else print,'There is no pixel in this slit.'
  endfor
  oldslitmap=readfits('/scr2/nichal/workspace/slitmaps/'+name+'_slitmap.fits')
  writefits,'/scr2/nichal/workspace/slitmaps/'+name+'_slitmap.fits',[[[oldslitmap]],[[slitmap]],[[N2map]]]
 
  distance = findgen(n_slits)*slit_size
  badpix=where(finite(N2_distance_err) eq 0 or N2_distance_err eq 0.) ;This will remove all the subslit with 1 pixel in it.
  if badpix(0) ne -1 then remove,badpix, distance,N2_distance,N2_distance_err
  
  ploterror,distance,N2_distance,N2_distance_err,xtitle='distance(kpc)',ytitle='N2 ([NII]/Ha)',psym=1,title='Best angle='+string(angle)
  stop
  
; Now get the N2 index as a function of radius along major axis (N2_radius)
  distance_indices = round(distance/slit_size)
  print,distance_indices,'(distance indices)'
  
  n_radii = max([max(distance_indices)-midpoint,midpoint])+1
  radius = findgen(n_radii)*slit_size ;kpc
  N2_radius= fltarr(n_radii)
  N2_radius_err= fltarr(n_radii)

  for i=0,n_radii-1 do begin
     good_radii = where(abs(distance_indices-midpoint) eq i)
     print,'good slit to be average', good_radii
     ;stop
     if good_radii(0) ne -1 and n_elements(good_radii) gt 1. then begin
        N2_radius(i) = mean(N2_distance(good_radii))
        good_error = N2_distance_err(good_radii)
        good_error_sq = good_error^2
        N2_radius_err(i) = sqrt(total(good_error_sq))/n_elements(good_radii)
     endif 
     if good_radii(0) ne -1 and n_elements(good_radii) eq 1. then begin
        N2_radius(i) = N2_distance(good_radii)
        N2_radius_err(i) = N2_distance_err(good_radii)
     endif
     if good_radii(0) eq -1 then begin
        N2_radius(i) = 1./0.
        N2_radius_err(i) = 1./0.
     endif
  endfor
  badelements=where(finite(N2_radius) eq 0.)
  if badelements(0) ne -1 then begin
     remove, badelements, N2_radius,N2_radius_err,radius
     print,'some elements have been removed'
     stop
  endif

  ploterror,radius,N2_radius,N2_radius_err,xtitle='radius(kpc)',ytitle='N2 ([NII]/Ha)',psym=1,title='Best angle='+string(angle)


;Fit the function
  A = [0.,8.9] ;guess parameter(gradient and Z center)
  weight = 1./N2_radius_err^2
  x=radius
  y=N2_radius
  N2fit = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04',/noderivative)  
  print, 'metal gradient =',A(0),'+/-',sigmaA(0)
  print, 'central metal =', A(1),'+/-',sigmaA(1)

if name eq 'cswa128' then begin ;the last 2 values are from interaction/merging
   print,'if remove the last 2 elements then'
   radius_short=radius
   N2_radius_short = N2_radius
   N2_radius_err_short = N2_radius_err
   remove,[5,6],radius_short,N2_radius_short,N2_radius_err_short
   A = [0.,8.9] ;guess parameter(gradient and Z center)
   weight = 1./N2_radius_err_short^2
   x=radius_short
   y=N2_radius_short
   N2fit = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04',/noderivative)  
   print, 'metal gradient =',A(0),'+/-',sigmaA(0)
   print, 'central metal =', A(1),'+/-',sigmaA(1)

endif


stop
end
