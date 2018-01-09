
function calB,RR_II,UU_BB,g,r
;transform SDSS u g r i to B magnitude
if RR_II lt 1.15 then begin
   if UU_BB lt 0 then B = g+0.33*(g-r)+0.2 else  B = g+0.39*(g-r)+0.21 
endif else B=99
return,B
end

pro photometry_c139

id_stars=[12,9,5,65,59,38,17,9,11,49]
sdss_u = [25.24,24.23,24.02,22.597,26.395,23.732,26.311,24.076,24.138,23.954]
sdss_g = [23.15,24.95,23.57,24.612,23.905,22.612,23.366,24.896,24.116,23.132]
sdss_r = [21.12,23.08,22.32,22.742,22.547,22.365,21.682,23.535,22.135,21.798]
sdss_i = [20.54,21.93,21.44,22.593,22.031,2,.687,21.320,22.254,21.511,21.276]
R_I = 0.930*(sdss_r-sdss_i)+0.259
U_B = 0.79*(sdss_u-sdss_g)-0.93

cat = rsex('cswa15_B.cat')
Breal = fltarr(n_elements(id_stars))
Bmeasure = Breal
Bmeasure_err = Breal

for i=0,n_Elements(id_stars)-1 do begin
   id = id_stars(i)
   pos = where(cat.number eq id, ct)
   if ct eq 0 then stop
   Breal(i)  = calB(R_I(i),U_B(i),sdss_g(i),sdss_r(i))
   Bmeasure(i) = cat(pos).mag_auto
   Bmeasure_err(i) = cat(pos).magerr_auto
endfor
bad = where(breal eq 99., ct)
if ct ne 0 then remove, bad, breal,bmeasure,bmeasure_err,id_stars

diff = breal-bmeasure
print, diff
print, 'median',median(diff)
print, 'mean',mean(diff)

;The lensed galaxy is id 54
pos = where(cat.number eq 54)
print, cat(pos).mag_auto
stop
end
