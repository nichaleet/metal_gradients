pro get_bcg,cat,posx,posy
readcol,cat,v1,V2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28
for i=0,n_elements(posx)-1 do begin
x=posx(i)
y=posy(i)
dis = (v2-x)^2+(v3-y)^2
number =  where(dis eq min(dis))
el  = (v21(number)-v22(number))/(v21(number)+v22(number))
print, 'number',v1(number),' x',v2(number),' y',v3(number),' mag',v8(number),'el',el, ' angle',v23(number)
endfor

;3 Brightest gals
print, '3 Brightest galaxies'
galregion = where(abs(v2-posx(0)) lt 200. and abs(v3-posy(0)) lt 200. and v8 gt 20)
v8temp = v8(galregion)

bright = []
for ii=0,2 do begin
   bright = [bright,where(v8temp eq min(v8temp))]
   v8temp(bright)= 99.
endfor
print, 'number',v1(galregion(bright))
print, 'x:',v2(galregion(bright))
print, 'y:',v3(galregion(bright))
print, 'mag', v8(galregion(bright))
print, 'elipticity', (v21(galregion(bright))-v22(galregion(bright)))/(v21(galregion(bright))+v22(galregion(bright)))
print, 'angle', v23(galregion(bright))
;stop
end
