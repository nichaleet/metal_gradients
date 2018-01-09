pro matchobj 

cat = 'cswa15_B.cat'

readcol,cat,v1,V2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28

cat = 'cswa15_R.cat'
readcol,cat,vR1,VR2,vR3,vR4,vR5,vR6,vR7,vR8,vR9,vR10,vR11,vR12,vR13,vR14,vR15, $
            vR16,vR17,vR18,vR19,vR20,vR21,vR22,vR23,vR24,vR25,vR26,vR27,vR28
distarr=[]
xdiffarr=[]
ydiffarr=[]
for i=0,n_elements(v2)-1 do begin
   dist=sqrt((v2(i)-vR2)^2+(v3(i)-vR3)^2)
   mindist = min(dist)
   pos = where(dist eq mindist)
   print, mindist
   if mindist(0) le 3. then begin
      distarr=[distarr,mindist]
      xdiffarr=[xdiffarr,v2(i)-vR2(pos)]
      ydiffarr=[ydiffarr,v3(i)-vR3(pos)]
   endif

endfor
!p.multi=[0,2,2]
plothist,distarr,/autobin
plothist,xdiffarr,/autobin
plothist,ydiffarr,/autobin
help,distarr
print,'mean xB-xR',mean(xdiffarr)
print,'mean yB-yR',mean(ydiffarr)
stop

xshift = round(mean(xdiffarr))
yshift = round(mean(ydiffarr))

end
