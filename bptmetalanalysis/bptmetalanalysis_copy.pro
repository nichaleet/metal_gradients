pro funcPP04,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((gradient*x+Zc-8.90)/0.57)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcSteidel,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((gradient*x+Zc-8.62)/0.36)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcPP04_O3N2,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((8.73-(gradient*x+Zc))/0.32)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcSteidel_O3N2,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((8.66-(gradient*x+Zc))/0.28)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro bptmetalanalysis,name,path,nii,nii_err,ha,ha_err,oiii,oiii_err,hb,hb_err,angle=angle,xc=xc,yc=yc,pixelscale=pixelscale,slitwidth=slitwidth

outdir = '/scr2/nichal/workspace/output/bptmetalanalysis_result/'
;set plot for eps maps
set_plot, 'ps'
page_height = 15.
page_width = 15.
plot_left = 3.
plot_bottom = 2.
xsize = 10.
ysize = 10.

;Define greek letters
alphaLetter = '!9' + String("141B) + '!X'  
betaletter  = '!9' + String("142B) + '!X'

;load data
nii     = readfits(path+nii)
nii_err = readfits(path+nii_err)
ha      = readfits(path+ha)
ha_err  = readfits(path+ha_err)
if file_search(path+oiii) ne '' then begin
   N2onlyflag = 0
   oiii    = readfits(path+oiii)
   oiii_err= readfits(path+oiii_err)
   hb      = readfits(path+hb)
   hb_err  = readfits(path+hb_err)
endif else begin
   N2onlyflag = 1
   oiii    = fltarr(size(ha,/dimension))+1./0.
   oiii_err= oiii
   hb      = oiii
   hb_err  = oiii
endelse
;fix the bad nii
bad = where(nii gt 1.)
nii(bad)     = 1./0.
nii_err(bad) = 1./0.
ha(bad)      = 1./0.
ha_err(bad)  = 1./0.

;calculate the n2, o3hb, and o3n2 indices
n2       = nii/ha
n2_err   = abs(n2)*sqrt((nii_err/nii)^2+(ha_err/ha)^2)
o3hb     = oiii/hb
o3hb_err = abs(o3hb)*sqrt((oiii_err/oiii)^2+(hb_err/hb)^2)

o3n2     = o3hb/n2
o3n2_Err = abs(o3n2)*sqrt((o3hb_err/o3hb)^2+(n2_err/n2)^2)

;remove the o3n2 points where both n2 and o3hb are negative
bad=where(n2 lt 0. and o3hb lt 0.)
if bad(0) ne -1 then begin
   o3n2(bad)=1./0.
   o3n2_err(bad)=1./0.
endif

;write N2 map and O3N2 map
writefits,outdir+name+'_N2map.fits',[[[N2]],[[N2_err]]]
if N2onlyflag eq 0 then begin
   writefits,outdir+name+'_O3N2map.fits',[[[O3N2]],[[O3N2_err]]]

;plot BPT diagram
   loadct,2
   set_plot,'ps'
   psname=outdir+name+'_bpt.eps'
   device, filename = psname, $
           xsize = 10, $
           ysize = 8, $
           xoffset = 0, $
           yoffset = 0, $
           scale_factor = 1.0,/encapsulated,/color

   ploterror,n2,o3hb,n2_err[*],o3hb_err[*],type=3,errcolor=240,xtitle='[NII]/H'+alphaletter,ytitle='[OIII]/H'+betaletter,font=0.,xrange=[1.e-3,0.8],yrange=[1,60.],psym=1

;Find out where is significantly in AGN region
   o3hb_lowerbound = o3hb-o3hb_err
   o3hb_lowerbound(where(o3hb lt 0.)) = -1./0.
   n2_lowerbound = n2-n2_err
   n2_lowerbound(where(n2 lt 0.)) = -1./0.
   finite_Elements = where(o3hb gt 0. and n2 gt 0.)
   bptmap = fltarr(size(ha,/dimension))+1./0.
   bptmap(where(finite(ha))) = 0.5 ; values for where the map is  
   agnflag = []                 ; to keep the points where it's likely an agn
   for i=0,n_Elements(finite_elements)-1 do begin
      element = finite_elements(i)
      if o3hb(element) gt 10.^(0.66/(alog10(n2(element))-0.31)+1.14) then bptmap(element)=0.5
      if o3hb_lowerbound(element) gt 10.^(0.66/(alog10(n2_lowerbound(element))-0.31)+1.14) then begin ;greater than ysteidel
         bptmap(element) = 1.
         agnflag = [agnflag,element]
      endif
   endfor
;oplot the errorbar in bpt diagram for the points that is likely AGN
   If n_elements(agnflag) gt 0. then begin
      oploterror,n2(agnflag),o3hb(agnflag),n2_err(agnflag),o3hb_err(agnflag),errcolor=100,psym=1
   endif
;oplot every points with x dots over the color bars
   loadct, 12
   oplot,n2,o3hb,psym=7,thick=2
;overplot kewley and steidel
   logx=findgen(78)/100.*5.-4.
   ykewley = 10.^( 0.61/(logx+0.08)+1.1)
   ysteidel = 10.^(0.66/(logx-0.31)+1.14)
   oplot,10.^logx,ykewley,color=15
   oplot,10.^logx,ysteidel,color=100
   device,/close

;plot bptmap
   writefits,outdir+name+'_bptmap.fits',bptmap

;plot bptmap to .eps
   psname = outdir+name+'_bptmap.eps'
   device, filename = psname, xsize = page_width,ysize = page_height, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated
   bptmap(where(finite(bptmap) eq 0.)) = 0.
   data = bytscl(bptmap, min = 0., max = 1.)
   tam = size(data, /dimensions)
   cgloadct, 3, ncolors = 256, bottom = 0, clip = [0,256], /reverse
   tvlct, redvector, greenvector, bluevector, /get
   cgimage,data, position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height] 
   device,/close
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Metalicity;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;1 The PP04,Steidel N2 index ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
sizeim = size(N2,/dimension)
width = sizeim(1)
n_slits = fix(1.4*width/slitwidth) ;number of subslits
if name eq 'cswa165' then n_slits = fix(2.*width/slitwidth)
pix_size = pixelscale ; in kpc per pixel
xmid = xc
ymid = yc
if angle ge 45. and angle le 135. then begin
   xref = fltarr(n_slits+1)+xmid    ; xref = xmid
   yref = findgen(n_slits+1)*slitwidth ;yref = 0,slitwidt,2slitwidth,...
endif
if angle gt 135. then begin
   yref = fltarr(n_slits+1)+ymid       ; yref = xymid
   xref = reverse(findgen(n_slits+1)*slitwidth) ; xref = n*slitwidth,(n-1)*slitwidth,...,slitwidth,0
endif
if angle lt 45. then begin
   yref = fltarr(n_slits+1)+ymid                ; yref = xmid
   xref = findgen(n_slits+1)*slitwidth ; xref = n*slitwidth,(n-1)*slitwidth,...,slitwidth,0
endif
;stop
ind = array_indices(N2,findgen(n_elements(N2)-1))
x = ind(0,*)
y = ind(1,*)
slit_size = slitwidth*pix_size*cos((angle-90.)*!pi/180.)

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
            N2_distance_err(i) = sigmad
         endif else begin
            N2_distance(i)=simplemean
            N2_distance_err(i) = stddev(goodN2(thegood))
         endelse 
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

;Get the N2 index as a function of radius along major axis (N2_radius)
n_radii = max([max(slit_indices)-midpoint,midpoint-min(slit_indices)])+1
radius = findgen(n_radii)*slit_size ;kpc
N2_radius= fltarr(n_radii)+1./0.
N2_radius_err= fltarr(n_radii)+1./0.
for i=0,n_radii-1 do begin
   good_radii = where(abs(slit_indices-midpoint) eq i)
   print,'good slit to be average', slit_indices(good_radii)
   if good_radii(0) ne -1 then begin
      goodslit_ind = slit_indices(good_radii)
      N2_inslit = []
      N2err_inslit = []
      n_detects=0.
      for j=0,n_elements(good_radii)-1 do begin
         good_pix_in_slit = where(N2mapslit eq goodslit_ind(j)*100. and finite(N2) eq 1 and N2+N2_err ge 0.)
         n_detects = n_detects+n_elements(where(N2mapslit eq goodslit_ind(j)*100. and finite(N2) eq 1))
         if good_pix_in_slit(0) ne -1 then begin
            N2_candidates = N2(good_pix_in_slit)
            N2err_candidates = N2_err(good_pix_in_slit)
            negative_N2 = where(N2_candidates lt 0.)
            if negative_N2(0) ne -1 then begin
               N2err_candidates(negative_N2)=N2_candidates(negative_N2)+N2err_candidates(negative_N2)
               N2_candidates(negative_N2) = 0.
            endif
            
            N2_inslit =[N2_inslit,N2_candidates]
            N2err_inslit = [N2err_inslit,N2err_candidates]
         endif
      endfor
                                
      if n_elements(N2_inslit) ge 1. and n_detects gt 1 then begin
;do not take the subslit where there is only 1 pixel of h_alpha detection (a lonely island)
         N2_radius(i) = mean(N2_inslit)
         good_error_sq = N2err_inslit^2
         N2_radius_err(i) = sqrt(total(good_error_sq))/n_elements(good_error_sq) ;normal error propagation
      endif else begin
         N2_radius(i) = 1./0.
         N2_radius_err(i) = 1./0.
      endelse
   endif
endfor

badelements=where(finite(N2_radius) eq 0. or n2_radius lt 0.001)
;if name eq 'cswa19_maingal' then badelements=where(finite(N2_radius) eq 0. or n2_radius lt 0.002)
if badelements(0) ne -1 then begin
   ;stop,'Stop. Some elements have been removed. You may want to check.'
   remove, badelements, N2_radius,N2_radius_err,radius     
endif

;plot N2 as a function of distance 
loadct,2
set_plot,'ps'
psname=outdir+name+'_N2distance.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
if min(N2_distance) lt 1.e-3 then yrange = [1.e-2,max(N2_distance)] else yrange = [min(N2_distance),max(N2_distance)] 

ploterror,distance,N2_distance,N2_distance_err,xtitle='distance(kpc)',ytitle='[NII]/H'+alphaletter,psym=7,type=1,font=0,errcolor=200,yrange=yrange,hatlength=!D.X_VSIZE/50.
  
device,/close

;plot N2 as a function of radius
psname=outdir+name+'_N2radius.eps'
device, filename = psname,xsize = 10,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

if min(N2_radius) lt 1.e-3 then yrange = [1.e-3,max(N2_radius)] else yrange = [min(N2_radius),max(N2_radius)] 

ploterror,radius,N2_radius,N2_radius_err,xtitle='radius(kpc)',ytitle='[NII]/H'+alphaletter,psym=7,type=1,font=0,errcolor=200,yrange=yrange,xmargin=[10,7],ystyle=8

;save the radius values
radius_of_n2 = radius

;Fit the PP04 function
A = [0.,8.90+0.57*alog10(N2_radius(0))]  ;guess parameter(gradient and Z center)
if A(1) lt 7. then A(1) = 7.2
weight = 1./N2_radius_err^2
x=radius
y=N2_radius
N2fit = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04',/noderivative)
print, 'N2 PP04:'
print, 'metal gradient =',A(0),'+/-',sigmaA(0)
print, 'central metal =', A(1),'+/-',sigmaA(1)

;make the fitted plot
xfit = !x.crange
funcPP04,xfit,A,yfit,pder
oplot,xfit,yfit,color=5

;add the second y axis in (12+logO/H)
metalrange=8.90+0.57*!y.crange
Axis,yaxis=1,yrange=metalrange,ylog=0,ytitle='12+log(O/H)!LPP04',font=0,ystyle=1

;Append the metal gradient values to the N2method gradient file
openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
printf,1, ''
printf,1, name
printf,1, systime()
printf,1, 'N2 PP04'
printf,1, 'metal gradient =',A(0),'+/-',sigmaA(0)
printf,1, 'central metal =', A(1),'+/-',sigmaA(1)
close,1

;Fit the Steidel function
N2fit = curvefit(x,y,weight,A,sigmaA,function_name='funcSteidel',/noderivative)
print, 'N2 Steidel:'
print, 'metal gradient =',A(0),'+/-',sigmaA(0)
print, 'central metal =', A(1),'+/-',sigmaA(1)

;make the fitted plot
funcSteidel,xfit,A,yfit,pder
oplot,xfit,yfit,color=5
device,/close

;Append the metal gradient values to the N2method gradient file
openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
printf,1, 'N2 Steidel'
printf,1, 'metal gradient =',A(0),'+/-',sigmaA(0)
printf,1, 'central metal =', A(1),'+/-',sigmaA(1)
close,1

if N2onlyflag eq 0 then begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;2 The PP04&Steidel O3N2 index ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   O3N2_distance = fltarr(n_slits)+1./0.
   O3N2_distance_err = fltarr(n_slits)+1./0.
   O3N2mapslit = O3N2
   slit_indices = findgen(n_slits)
   for i=0,n_slits-1 do begin
      pix_in_slit = where(N2mapslit eq i*100.)
      O3N2mapslit(pix_in_slit) = i*100.
      if pix_in_slit(0) ne -1 then begin
         goodO3N2 = O3N2(pix_in_slit) 
         goodO3N2err = O3N2_err(pix_in_slit)    
         thegood = where(finite(goodO3N2) eq 1 and goodO3N2+goodo3n2err ge 0.)
         if thegood(0) ne -1 then begin
            negativept = where(goodo3n2 lt 0. and goodo3n2+goodo3n2err gt 0.)
            if negativept(0) ne -1 then begin
               goodo3n2err(negativept) = goodO3n2(negativept)+goodO3n2err(negativept)
               goodo3n2(negativept) = 0.
            endif
            print,'There are', n_Elements(thegood),'pixels with finite values in this slit'
            simplemean = mean(goodO3N2(thegood))
            meanerr,goodO3N2(thegood),goodO3N2err(thegood),wmean,sigmam,sigmad
            print,'subslit',i,' N2index=',simplemean,wmean
            if finite(wmean) eq 1 then begin 
               O3N2_distance(i) = wmean 
               O3N2_distance_err(i) = sigmad
            endif else begin
               O3N2_distance(i)=simplemean
               O3N2_distance_err(i) = stddev(goodO3N2(thegood))
            endelse 
         endif
      endif else print,'There is no pixel in this slit.'
   endfor
;write O3N2 map 
   writefits,outdir+name+'_O3N2map.fits',[[[O3N2]],[[O3N2_err]],[[O3N2mapslit]]]
;Remove all the subslit with 1 pixel in it.
   badslit=where(finite(O3N2_distance_err) eq 0 or O3N2_distance_err eq 0.)
   if badslit(0) ne -1 then remove,badslit,O3N2_distance,O3N2_distance_err,slit_indices
;correct the distance
   distance = (slit_indices-midpoint)*slit_size
   
;Get the O3N2 index as a function of radius along major axis (O3N2_radius)
   n_radii = max([max(slit_indices)-midpoint,midpoint-min(slit_indices)])+1
   radius = findgen(n_radii)*slit_size ;kpc
   O3N2_radius= fltarr(n_radii)+1./0.
   O3N2_radius_err= fltarr(n_radii)+1./0.
   for i=0,n_radii-1 do begin
      good_radii = where(abs(slit_indices-midpoint) eq i)
      print,'good slit to be average', slit_indices(good_radii)
      if good_radii(0) ne -1 then begin
         goodslit_ind = slit_indices(good_radii)
         O3N2_inslit = []
         O3N2err_inslit = []
         for j=0,n_elements(good_radii)-1 do begin
            good_pix_in_slit = where(N2mapslit eq goodslit_ind(j)*100. and finite(O3N2) eq 1 and O3N2+O3N2_err ge 0.)
            if good_pix_in_slit(0) ne -1 then begin
               O3N2_candidates = O3N2(good_pix_in_slit)
               O3N2err_candidates = O3N2_err(good_pix_in_slit)
               negative_O3N2 = where(O3N2_candidates lt 0.)
               if negative_O3N2(0) ne -1 then begin
                  O3N2err_candidates(negative_O3N2)=O3N2_candidates(negative_O3N2)+O3N2err_candidates(negative_O3N2)
                  O3N2_candidates(negative_O3N2) = 0.
               endif
               O3N2_inslit =[O3N2_inslit,O3N2_candidates]
               O3N2err_inslit = [O3N2err_inslit,O3N2err_candidates]
            endif
         endfor
         if n_elements(O3N2_inslit) gt 1. then begin
            o3n2modified = where(o3n2_inslit eq 0)
            if o3n2modified(0) ne -1 then begin
               O3N2_radius(i) = mean(O3N2_inslit)
               good_error_sq = O3N2err_inslit^2
               O3N2_radius_err(i) = sqrt(total(good_error_sq))/n_elements(good_error_sq) ;normal error propagation. Because the new uncertainty of the modified one will be smaller than usual. Then the weighted mean will give more faulty weight to these modified ones.
            endif else begin
               meanerr,o3n2_inslit,o3n2err_inslit,wmean,sigmam,sigmad
               o3n2_radius(i) = wmean
               o3n2_radius_err(i) = sigmam
            endelse
         endif else begin
            O3N2_radius(i) = 1./0.
            O3N2_radius_err(i) = 1./0.
         endelse
      endif
   endfor
   
   badelements=where(finite(O3N2_radius) eq 0.)
   if badelements(0) ne -1 then begin
                                ;stop,'Stop: Some elements have been removed. You may want to check.'
      remove, badelements, O3N2_radius,O3N2_radius_err,radius     
   endif
   
;plot O3N2 as a function of distance 
   loadct,2
   set_plot,'ps'
   psname=outdir+name+'_O3N2distance.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   
   ploterror,distance,O3N2_distance,O3N2_distance_err[*],xtitle='distance(kpc)',ytitle='([OIII]/H'+betaletter+')/([NII]/H'+alphaletter+')',psym=7,type=1,font=0,errcolor=200,yrange=[min(O3N2_distance),max(O3N2_distance)],HATLENGTH=!D.X_VSIZE/50.
   
   device,/close

;plot O3N2 as a function of radius
   psname=outdir+name+'_O3N2radius.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   ploterror,radius,O3N2_radius,O3N2_radius_err[*],xtitle='radius(kpc)',ytitle='([OIII]/H'+betaletter+')/([NII]/H'+alphaletter+')',psym=7,type=1,font=0,errcolor=200,yrange=[min(O3N2_radius),300],ystyle=8,xmargin=[10,7],HATLENGTH=!D.X_VSIZE/50.
   
;Fit the function
   A = [0.,8.66-0.28*alog10(O3N2_radius(0))]                 ;guess parameter(gradient and Z center)
   weight = 1./O3N2_radius_err^2
   x = radius
   y = O3N2_radius
   O3N2fit = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04_O3N2',/noderivative)
   print, 'O3N2 PP04:'
   print, 'metal gradient =',A(0),'+/-',sigmaA(0)
   print, 'central metal =', A(1),'+/-',sigmaA(1)

;make the fitted plot
   xfit = !x.crange
   funcPP04_O3N2,xfit,A,yfit,pder
   oplot,xfit,yfit,color=5
   
;Append the metal gradient values to the gradient file
   openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
   printf,1, 'O3N2 PP04'
   printf,1, 'metal gradient =',A(0),'+/-',sigmaA(0)
   printf,1, 'central metal =', A(1),'+/-',sigmaA(1)
   close,1

;Fit the Steidel function
   A = [0.,8.66-0.28*alog10(O3N2_radius(0))]                 ;guess parameter(gradient and Z center)
   x = radius
   y = O3N2_radius
   O3N2fit = curvefit(x,y,weight,A,sigmaA,function_name='funcSteidel_O3N2',/noderivative)
   print, 'O3N2 Steidel:'
   print, 'metal gradient =',A(0),'+/-',sigmaA(0)
   print, 'central metal =', A(1),'+/-',sigmaA(1)

;make the fitted plot
   xfit = !x.crange
   funcSteidel_O3N2,xfit,A,yfit,pder
   oplot,xfit,yfit,color=5

;add the second y axis in (12+logO/H)
   metalrange=8.66-0.28*!y.crange
   Axis,yaxis=1,yrange=metalrange,ylog=0,ytitle='12+log(O/H)!L S14',font=0,ystyle=1
   
   device,/close
   
;Append the metal gradient values to the gradient file
   openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
   printf,1, 'O3N2 Steidel'
   printf,1, 'metal gradient =',A(0),'+/-',sigmaA(0)
   printf,1, 'central metal =', A(1),'+/-',sigmaA(1)
   close,1
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;3 O3N2 Maiolino08: measure metalicity at each pixel;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
m08exist = file_search('/scr2/nichal/workspace/output/bptmetalanalysis_result/',+name+'_m08_metal.fits')
if m08exist eq '' then begin
   set_plot,'x'
   Bayesian_metal,name, OIII,hb,NII,Ha,OIII_err,hb_err,NII_err,ha_err,m08metal,median_m08metal,lowerbound_m08metal,upperbound_m08metal,typebayesian
   m08metal_err = (upperbound_m08metal-lowerbound_m08metal)/2.
   writefits,'/scr2/nichal/workspace/output/bptmetalanalysis_result/'+name+'_m08_metal.fits',[[[m08metal]],[[median_m08metal]],[[lowerbound_m08metal]],[[upperbound_m08metal]],[[m08metal_err]],[[typebayesian]]]
endif else begin
   m08 = readfits('/scr2/nichal/workspace/output/bptmetalanalysis_result/'+name+'_m08_metal.fits')
   m08metal            = m08[*,*,0]
   median_m08metal     = m08[*,*,1]
   lowerbound_m08metal = m08[*,*,2]
   upperbound_m08metal = m08[*,*,3]
   m08metal_err        = m08[*,*,4]
   typebayesian        = m08[*,*,5]
endelse
m08 = m08metal
m08_err = m08metal_err

m08_distance = fltarr(n_slits)+1./0.
m08_distance_err = fltarr(n_slits)+1./0.
m08mapslit = m08metal
slit_indices = findgen(n_slits)
for i=0,n_slits-1 do begin
   pix_in_slit = where(N2mapslit eq i*100.)
   m08mapslit(pix_in_slit) = i*100.
   if pix_in_slit(0) ne -1 then begin
      goodm08 = m08(pix_in_slit) 
      goodm08err = m08_err(pix_in_slit)
      goodtype   = typebayesian(pix_in_slit)
      thegood = where(goodtype mod 27 eq 0 or goodtype mod 8 eq 0 and goodtype  ne 0 and finite(goodm08) eq 1) ;take the pixels where either Nii or Hb is positively detected
      if thegood(0) ne -1 then begin
         print,'There are', n_Elements(thegood),'pixels with finite values in this slit'
         simplemean = mean(goodm08(thegood))
         meanerr,goodm08(thegood),goodm08err(thegood),wmean,sigmam,sigmad
         print,'subslit',i,' N2index=',simplemean,wmean
         if finite(wmean) eq 1 then begin 
            m08_distance(i) = wmean 
            m08_distance_err(i) = sigmad
         endif else begin
            m08_distance(i)=simplemean
            m08_distance_err(i) = stddev(goodm08(thegood))
         endelse 
      endif
   endif else print,'There is no pixel in this slit.'
endfor
;write m08 map 
writefits,outdir+name+'_m08map.fits',[[[m08]],[[m08_err]],[[m08mapslit]],[[typebayesian]]]
;Remove all the subslit with 1 pixel in it.
  badslit=where(finite(m08_distance_err) eq 0)
  if badslit(0) ne -1 then remove,badslit,m08_distance,m08_distance_err,slit_indices
;correct the distance
distance = (slit_indices-midpoint)*slit_size

;Get the m08 index as a function of radius along major axis (m08_radius)
n_radii = max([max(slit_indices)-midpoint,midpoint-min(slit_indices)])+1
radius = findgen(n_radii)*slit_size ;kpc
m08_radius= fltarr(n_radii)+1./0.
m08_radius_err= fltarr(n_radii)+1./0.
for i=0,n_radii-1 do begin
   good_radii = where(abs(slit_indices-midpoint) eq i)
   print,'good slit to be average', slit_indices(good_radii)
   if good_radii(0) ne -1 then begin
      goodslit_ind = slit_indices(good_radii)
      m08_inslit = []
      m08err_inslit = []
      n_detects=0.
      for j=0,n_elements(good_radii)-1 do begin
         good_pix_in_slit = where(typebayesian mod 27 eq 0 or typebayesian mod 8 eq 0 or typebayesian eq 8 and typebayesian ne 0. and N2mapslit eq goodslit_ind(j)*100. and finite(m08) eq 1)

         n_detects = n_detects+n_elements(where(N2mapslit eq goodslit_ind(j)*100. and finite(N2) eq 1)) ;Using finite(N2) is correct. This takes where Ha is detected. 
         if good_pix_in_slit(0) ne -1 then begin
            m08_inslit =[m08_inslit,m08(good_pix_in_slit)]
            m08err_inslit = [m08err_inslit,m08_err(good_pix_in_slit)]
         endif
      endfor
      
      if n_elements(m08_inslit) eq 1. and n_detects gt 1. then begin
         m08_radius(i) = mean(m08_inslit)
         good_error_sq = m08err_inslit^2
         m08_radius_err(i) = sqrt(total(good_error_sq))/n_elements(good_error_sq) ;normal error propagation
      endif else if n_elements(m08_inslit) gt 1. then begin 
         meanerr,m08_inslit,m08err_inslit,wmean,sigmam,sigmad
         m08_radius(i) = wmean
         m08_radius_err(i) = sigmam
      endif else begin
         m08_radius(i) = 1./0.
         m08_radius_err(i) = 1./0.
      endelse
   endif
endfor

badelements=where(finite(m08_radius) eq 0.)
if badelements(0) ne -1 then begin
   remove, badelements, m08_radius,m08_radius_err,radius     
endif

;plot m08 as a function of distance 
loadct,2
set_plot,'ps'
psname=outdir+name+'_m08distance.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

ploterror,distance,m08_distance,m08_distance_err,xtitle='distance(kpc)',ytitle='12+log(O/H)!LM08',psym=7,font=0,errcolor=200,hatlength=!D.X_VSIZE/50.
  
device,/close

;plot m08 as a function of radius
psname=outdir+name+'_m08radius.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,radius,m08_radius,m08_radius_err,xtitle='radius(kpc)',ytitle='12+log(O/H)!LM08',psym=7,font=0,errcolor=200,hatlength=!D.X_VSIZE/50.

;Fit the function
x = radius
y = m08_radius
A=linfit(radius,m08_radius,measure_errors=m08_radius_err,sigma=sigmaA)

print, 'O3N2 m08:'
print, 'metal gradient =',A(1),'+/-',sigmaA(1)
print, 'central metal =', A(0),'+/-',sigmaA(0)

;make the fitted plot
xfit = !x.crange
yfit = A(0)+A(1)*xfit
oplot,xfit,yfit,color=5
device,/close
;Append the metal gradient values to the gradient file
openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
printf,1, 'O3N2 M08'
printf,1, 'metal gradient =',A(1),'+/-',sigmaA(1)
printf,1, 'central metal =', A(0),'+/-',sigmaA(0)
close,1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;4 N2, O3N2 Maiolino08: measure metalicity at each subslit;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Already have radius_of_n2,n2_radius, n2_radius_err
;calculate the metal
;N2
Bayesian_metal_slit_n2,name,n2_radius,n2_radius_err,radius_of_n2,m08_n2_metal,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,m08_n2_metal_err

;plot m08 as a function of radius
psname=outdir+name+'_m08_n2_slit.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,radius,m08_n2_metal,m08_n2_metal_err,xtitle='radius(kpc)',ytitle='12+log(O/H)!LM08,N2',psym=7,font=0,errcolor=200,hatlength=!D.X_VSIZE/50.

;Fit the function
A=linfit(radius_of_n2,m08_n2_metal,measure_errors=m08_n2_metal_err,sigma=sigmaA)

print, 'N2 slit  m08:'
print, 'metal gradient =',A(1),'+/-',sigmaA(1)
print, 'central metal =', A(0),'+/-',sigmaA(0)

;make the fitted plot
xfit = !x.crange
yfit = A(0)+A(1)*xfit
oplot,xfit,yfit,color=5
device,/close
;Append the metal gradient values to the gradient file

openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
printf,1, 'N2 M08 Slit'
printf,1, 'metal gradient =',A(1),'+/-',sigmaA(1)
printf,1, 'central metal =', A(0),'+/-',sigmaA(0)
close,1

;O3N2
if N2onlyflag eq 0 then begin
;Need to find o3hb_radius
   o3hb_distance = fltarr(n_slits)+1./0.
   o3hb_distance_err = fltarr(n_slits)+1./0.
   o3hbmapslit = o3hb
   slit_indices = findgen(n_slits)
   for i=0,n_slits-1 do begin
      pix_in_slit = where(N2mapslit eq i*100.)
      o3hbmapslit(pix_in_slit) = i*100.
      if pix_in_slit(0) ne -1 then begin
         goodo3hb = o3hb(pix_in_slit) 
         goodo3hberr = o3hb_err(pix_in_slit)    
         thegood = where(finite(goodo3hb) eq 1 and goodo3hb ge 0.)
         if thegood(0) ne -1 then begin
            print,'There are', n_Elements(thegood),'o3hb pixels with positive values in this slit'
            simplemean = mean(goodo3hb(thegood))
            meanerr,goodo3hb(thegood),goodo3hberr(thegood),wmean,sigmam,sigmad
            print,'subslit',i,' N2index=',simplemean,wmean
            if finite(wmean) eq 1 then begin 
               o3hb_distance(i) = wmean 
               o3hb_distance_err(i) = sigmad
            endif else begin
               o3hb_distance(i)=simplemean
               o3hb_distance_err(i) = stddev(goodo3hb(thegood))
            endelse 
         endif
      endif else print,'There is no pixel in this slit.'
   endfor
;write o3hb map 
   writefits,outdir+name+'_o3hbmap.fits',[[[o3hb]],[[o3hb_err]],[[o3hbmapslit]]]
;Remove all the subslit with 1 pixel in it.
   badslit=where(finite(o3hb_distance_err) eq 0 or o3hb_distance_err eq 0.)
   if badslit(0) ne -1 then remove,badslit,o3hb_distance,o3hb_distance_err,slit_indices
;correct the distance
   distance = (slit_indices-midpoint)*slit_size
   
;Get the o3hb index as a function of radius along major axis (o3hb_radius)
   n_radii = max([max(slit_indices)-midpoint,midpoint-min(slit_indices)])+1
   radius = findgen(n_radii)*slit_size ;kpc
   o3hb_radius= fltarr(n_radii)+1./0.
   o3hb_radius_err= fltarr(n_radii)+1./0.
   for i=0,n_radii-1 do begin
      good_radii = where(abs(slit_indices-midpoint) eq i)
      print,'good slit to be average', slit_indices(good_radii)
      if good_radii(0) ne -1 then begin
         goodslit_ind = slit_indices(good_radii)
         o3hb_inslit = []
         o3hberr_inslit = []
         for j=0,n_elements(good_radii)-1 do begin
            good_pix_in_slit = where(N2mapslit eq goodslit_ind(j)*100. and finite(o3hb) eq 1 and o3hb ge 0.)
            if good_pix_in_slit(0) ne -1 then begin
               o3hb_candidates = o3hb(good_pix_in_slit)
               o3hberr_candidates = o3hb_err(good_pix_in_slit)
               o3hb_inslit =[o3hb_inslit,o3hb_candidates]
               o3hberr_inslit = [o3hberr_inslit,o3hberr_candidates]
            endif
         endfor
         if n_elements(o3hb_inslit) gt 1. then begin
            meanerr,o3hb_inslit,o3hberr_inslit,wmean,sigmam,sigmad
            o3hb_radius(i) = wmean
            o3hb_radius_err(i) = sigmam
         endif else begin
            o3hb_radius(i) = 1./0.
            o3hb_radius_err(i) = 1./0.
         endelse
      endif
   endfor
   
   badelements=where(finite(o3hb_radius) eq 0.)
   if badelements(0) ne -1 then begin
      remove, badelements, o3hb_radius,o3hb_radius_err,radius     
   endif
   
;plot o3hb as a function of distance 
   loadct,2
   set_plot,'ps'
   psname=outdir+name+'_o3hbdistance.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   
   ploterror,distance,o3hb_distance,o3hb_distance_err[*],xtitle='distance(kpc)',ytitle='[OIII]/H'+betaletter,psym=7,type=1,font=0,errcolor=200,yrange=[min(o3hb_distance),max(o3hb_distance)],HATLENGTH=!D.X_VSIZE/50.
   
   device,/close

;plot o3hb as a function of radius
   psname=outdir+name+'_o3hbradius.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   ploterror,radius,o3hb_radius,o3hb_radius_err[*],xtitle='radius(kpc)',ytitle='[OIII]/H'+betaletter,psym=7,type=1,font=0,errcolor=200,yrange=[min(o3hb_radius),max(o3hb_radius)],xmargin=[10,7],HATLENGTH=!D.X_VSIZE/50.
device,/close   
   
radius_of_o3hb=radius
print,'check radius of n2 and o3hb'
print, radius_of_n2
print, radius_of_o3hb
;stop

;find the match of radius
match,radius_of_n2,radius_of_o3hb,sub_n2,sub_o3hb
radius = radius_of_n2(sub_n2)
o3hb_radius_match = o3hb_radius(sub_o3hb)
o3hb_radius_err_match = o3hb_radius_err(sub_o3hb)
n2_radius_match   = n2_radius(sub_n2)
n2_radius_err_match   = n2_radius_err(sub_n2)

;O3N2
Bayesian_metal_slit_o3n2,name,n2_radius_match,n2_radius_err_match,o3hb_radius_match,o3hb_radius_err_match,radius,m08_o3n2_metal,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,m08_o3n2_metal_err
 
;plot m08 as a function of radius
psname=outdir+name+'_m08_o3n2_slit.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,radius,m08_o3n2_metal,m08_o3n2_metal_err,xtitle='radius(kpc)',ytitle='12+log(O/H)!LM08,O3N2',psym=7,font=0,errcolor=200,hatlength=!D.X_VSIZE/50.

;Fit the function
A=linfit(radius,m08_o3n2_metal,measure_errors=m08_o3n2_metal_err,sigma=sigmaA)

print, 'O3N2 slit  m08:'
print, 'metal gradient =',A(1),'+/-',sigmaA(1)
print, 'central metal =', A(0),'+/-',sigmaA(0)

;make the fitted plot
xfit = !x.crange
yfit = A(0)+A(1)*xfit
oplot,xfit,yfit,color=5
device,/close
;Append the metal gradient values to the gradient file

openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
printf,1, 'O3N2 M08 Slit'
printf,1, 'metal gradient =',A(1),'+/-',sigmaA(1)
printf,1, 'central metal =', A(0),'+/-',sigmaA(0)
close,1
endif


stop
end
