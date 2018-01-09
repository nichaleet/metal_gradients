pro N2distance_slitmap,name,outdir,n_slits,slitwidth,pixelscale,xc,yc,angle,N2,N2_err,N2_distance,N2_distance_err,distance,N2mapslit,midpoint,slit_indices,slit_size

pix_size = pixelscale ; in kpc per pixel
xmid = xc
ymid = yc
if angle ge 45. and angle le 135. then begin
   xref = fltarr(n_slits+1)+xmid    ; xref = xmid
   yref = findgen(n_slits+1)*slitwidth ;yref = 0,slitwidt,2slitwidth,...
   slit_size = slitwidth*pix_size*sin(angle*!pi/180.)
endif
if angle gt 135. then begin
   yref = fltarr(n_slits+1)+ymid       ; yref = xymid
   xref = reverse(findgen(n_slits+1)*slitwidth) ; xref = n*slitwidth,(n-1)*slitwidth,...,slitwidth,0
   slit_size = -1.*slitwidth*pix_size*cos(angle*!pi/180.)
endif
if angle lt 45. then begin
   yref = fltarr(n_slits+1)+ymid                ; yref = xmid
   xref = findgen(n_slits+1)*slitwidth ; xref = n*slitwidth,(n-1)*slitwidth,...,slitwidth,0
   slit_size = slitwidth*pix_size*cos(angle*!pi/180.)
endif
;stop
ind = array_indices(N2,findgen(n_elements(N2)-1))
x = ind(0,*)
y = ind(1,*)

;making subslits
m_major = tan(angle/180.*!pi)  
c_major = ymid-m_major*xmid
;distance of a point(x,y) to this major axis is |-mx+y-c|/sqrt(m^2+1)
dist = abs((-1.*m_major*x)+y-c_major)/sqrt(m_major^2+1.)
;Next are the angle and equations for the slit dividers 
angle_perpend = angle-90.       ;angle of the little slit
slope = tan(angle_perpend/180.*!pi)  
; y = slope*x+intercept
; intercept = y-slope*x
intercept = yref-slope(0)*xref

N2_distance = fltarr(n_slits)
N2_distance_err = fltarr(n_slits)
N2mapslit = N2
slit_indices = findgen(n_slits)

for i=0,n_slits-1 do begin
   pix_in_slit = where(y le slope*x+intercept(i+1) and y ge slope*x+intercept(i) and dist le 1.5) ; The slit width is ~3 pixels     
   if pix_in_slit(0) ne -1 then begin
      goodN2 = N2(pix_in_slit) 
      goodN2err = N2_err(pix_in_slit)    
      N2mapslit(pix_in_slit)=i*100
      thegood = where(finite(goodN2) eq 1 and goodN2+goodn2err ge 0.)
      if thegood(0) ne -1 then begin
         negativept = where(goodn2 lt 0. and goodn2+goodn2err gt 0.)
         if negativept(0) ne -1 then begin
            goodn2(negativept) = 0.
            goodn2err(negativept) = goodn2(negativept)+goodn2err(negativept)
         endif
         print,'There are', n_Elements(thegood),'pixels with finite values in this slit'
         simplemean = mean(goodN2(thegood))
         meanerr,goodN2(thegood),goodN2err(thegood),wmean,sigmam,sigmad
         print,'subslit',i,' N2index=',simplemean,wmean
         if finite(wmean) eq 1 then begin 
            N2_distance(i) = wmean 
            N2_distance_err(i) = sigmam
         endif else stop,'Stop: Error with taking meanerr'
      endif else begin
         n2_distance(i)=1./0.
         n2_distance_err(i)=1./0.
      endelse
   endif else print,'There is no pixel in this slit.'
endfor
;write N2 map 
writefits,outdir+name+'_N2map.fits',[[[N2]],[[N2_err]],[[N2mapslit]]]
;find the subslit in which the center lies
  midindex = N2mapslit[xc,yc]/100.
  midpoint = where(slit_indices eq midindex)
  if midpoint(0) eq -1 or n_elements(midpoint) gt 1. then stop,'ERROR:CANNOT FIND THE LOCATION OF THE INPUT CENTER'
  midpoint = slit_indices(midpoint)
  midpoint = midpoint(0)
;Remove all the subslit with 1 pixel in it.
  badslit=where(finite(N2_distance_err) eq 0 or N2_distance_err eq 0.)
  if badslit(0) ne -1 then remove,badslit,N2_distance,N2_distance_err,slit_indices
;correct the distance
distance = (slit_indices-midpoint)*slit_size


end
