pro photometry_c28
id = [250,182,194,88,160,256,231]
real_iso=[22.8136,22.615,22.6129,22.6089,22.4641,22.4187,22.1322]
cat=rsex('cswa28_B.cat')
measure = real_iso*0.

for i=0, n_Elements(id)-1 do begin
   pos = where(cat.number eq id(i),ct)
   if ct ne 1 then stop
   measure(i)=cat[pos].mag_iso
endfor
diff = measure-real_iso
print,diff
stop
end
