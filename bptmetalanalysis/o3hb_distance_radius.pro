pro o3hb_distance_radius,o3hb,o3hb_err,n2mapslit,n_slits,midpoint,slit_size,o3hbmapslit,distance,o3hb_distance,o3hb_distance_err,radius,o3hb_radius,o3hb_radius_err

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
           ;simplemean = mean(goodo3hb(thegood))
           meanerr,goodo3hb(thegood),goodo3hberr(thegood),wmean,sigmam,sigmad
           ;print,'subslit',i,' o3hbindex=',simplemean,wmean
           if finite(wmean) eq 1 then begin 
              o3hb_distance(i) = wmean 
              o3hb_distance_err(i) = sigmam
           endif else begin
              stop,'ERROR: WMEAN is not working'
              ;o3hb_distance(i)=simplemean
              ;o3hb_distance_err(i) = stddev(goodo3hb(thegood))
           endelse 
        endif
     endif else print,'There is no pixel in this slit.'
  endfor

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
end
