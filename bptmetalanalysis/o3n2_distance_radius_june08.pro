pro o3n2_distance_radius,o3n2,o3n2_err,n2mapslit,n_slits,midpoint,slit_size,o3n2mapslit,distance,o3n2_distance,o3n2_distance_err,radius,o3n2_radius,o3n2_radius_err

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
           print,'There are', n_Elements(thegood),'pixels with finite values in this slit'
                                ;simplemean = mean(goodO3N2(thegood))
           meanerr,goodO3N2(thegood),goodO3N2err(thegood),wmean,sigmam,sigmad
                                ;print,'subslit',i,' N2index=',simplemean,wmean
           if finite(wmean) eq 1 then begin 
              O3N2_distance(i) = wmean 
              O3N2_distance_err(i) = sigmad
           endif else stop,'STOP:ERROR with meanerr'
        endif
     endif else print,'There is no pixel in this slit.'
  endfor

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
           good_pix_in_slit = where(N2mapslit eq goodslit_ind(j)*100. and finite(O3N2) eq 1 and O3N2 gt 0.)
           if good_pix_in_slit(0) ne -1 then begin
              O3N2_candidates = O3N2(good_pix_in_slit)
              O3N2err_candidates = O3N2_err(good_pix_in_slit)
              O3N2_inslit =[O3N2_inslit,O3N2_candidates]
              O3N2err_inslit = [O3N2err_inslit,O3N2err_candidates]
            endif
        endfor
        if n_elements(O3N2_inslit) gt 1. then begin
           ;meanerr,o3n2_inslit,o3n2err_inslit,wmean,sigmam,sigmad
           meanerr,o3n2_inslit,o3n2err_inslit*0.+1.,wmean,sigmam,sigmad
           o3n2_radius(i) = wmean
           o3n2_radius_err(i) = sigmam
        endif else begin
           O3N2_radius(i) = 1./0.
           O3N2_radius_err(i) = 1./0.
        endelse
     endif
  endfor
   
  badelements=where(finite(O3N2_radius) eq 0.)
  if badelements(0) ne -1 then begin
     remove, badelements, O3N2_radius,O3N2_radius_err,radius     
  endif
  
end

