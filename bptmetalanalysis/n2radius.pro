pro n2radius,name,outdir,N2mapslit,N2,N2_err,midpoint,slit_indices,slit_size,N2_radius,N2_radius_err,radius

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
         good_pix_in_slit = where(N2mapslit eq goodslit_ind(j)*100. and finite(N2) eq 1 and N2 gt -3.5)
         n_detects = n_detects+n_elements(where(N2mapslit eq goodslit_ind(j)*100. and finite(N2) eq 1) and N2 gt -3.5)
         if good_pix_in_slit(0) ne -1 then begin
            N2_candidates = N2(good_pix_in_slit)
            N2err_candidates = N2_err(good_pix_in_slit)
            N2_inslit =[N2_inslit,N2_candidates]
            N2err_inslit = [N2err_inslit,N2err_candidates]
         endif
      endfor
                                
      if n_elements(N2_inslit) ge 1. and n_detects gt 1 then begin
;do not take the subslit where there is only 1 pixel of h_alpha detection (a lonely island)
         meanerr,n2_inslit,n2err_inslit,wmean,sigmam,sigmad  ; weighted
         ;print,radius(i),'kpc :',n2_inslit
         ;print, 'mean',wmean,mean(n2_inslit)
         ;waitt = get_kbrd()
         N2_radius(i) = wmean
         N2_radius_err(i) = sigmam
      endif else begin
         N2_radius(i) = 1./0.
         N2_radius_err(i) = 1./0.
      endelse
   endif
endfor
;stop
badelements=where(finite(N2_radius) eq 0.)
;if name eq 'cswa19_maingal' then badelements=where(finite(N2_radius) eq 0. or n2_radius lt 0.002)
if badelements(0) ne -1 then begin
   ;stop,'Stop. Some elements have been removed. You may want to check.'
   remove, badelements, N2_radius,N2_radius_err,radius     
endif

end
