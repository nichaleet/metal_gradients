

;This chunck is taken from bptmetalanalysis.pro 
;It calculate M08 for individual pixel then calculate metallicity as a function of radius/distance from measured metalicity.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;3 O3N2 Maiolino08: measure metalicity at each pixel;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;----------------- Calculate metallicity at each pixel -----------------------
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
            beep
            beep
            beep
            stop
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
      
      if n_elements(m08_inslit) ge 1. and n_detects gt 1. then begin
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
