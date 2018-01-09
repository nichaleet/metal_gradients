pro makeregion


;Write objects of which we want to get spectra.
;Numbers are in the F105 catalog
;READ IN CATALOT
cat = RSEX('cswa31_R.cat')
pos = where(cat.mag_auto le 18.)
mag = strtrim(cat.mag_auto,2)
;write_ds9_regionfile,cat(pos).alpha_j2000,cat(pos).delta_j2000,comment=mag(pos),filename='c31_R.reg',color='red'

obj1 =[882,897,637] ;Images of the arcs
obj2 =[1138,627,634]            ; Image of the arc but less important
obj3 = [838,620,895,840,851]   ;cluster members (I guessed)
obj4 = [925] ;lensed background objs
obj5 = [961,975,1081,1016,301,439,540,1099,116] ;random objects that seems interesting

listobj=list(obj1,obj2,obj3,obj4,obj5)
nameobj=['arcm','arc','cl','bg','rd']
importantvalue = [500,450,200,550,30]

openw,lun,'c31mask.txt',/get_lun
posarr=[]
namearr=[]
for ni = 0,n_Elements(listobj)-1 do begin
   objnow = listobj(ni)
   for ii=0,n_elements(objnow)-1 do begin
      id = objnow(ii)
      pos = where(cat.number eq id)
      pos = pos(0)
      posarr=[posarr,pos]
      namearr=[namearr,nameobj(ni)+strtrim(ii+1,2)]
      radec, cat(pos).alpha_J2000,cat(pos).delta_j2000,ihr,imin,xsec,ideg,imn,xsc
      mag = cat(pos).mag_auto
      printf,lun,nameobj(ni)+strtrim(ii+1,2)+' '+strtrim(importantvalue(ni),2)+' '+strtrim(mag,2)+' '+strtrim(ihr,2)+' '+strtrim(imin,2)+' '+strtrim(xsec,2)+' '+strtrim(ideg,2)+' '+strtrim(imn,2)+' '+strtrim(xsc,2)+' 2000.0 2000.0 0 0'
   endfor
endfor
free_lun,lun

write_ds9_regionfile,cat(posarr).alpha_j2000,cat(posarr).delta_j2000,comment=namearr,filename='objects.reg',color='green'

;2mass point source from the catalog in ds9
cat = read_ascii('stdstars.tsv',data_start=1,header=header)
cat = cat.field01
ra = cat(0,*)
dec= cat(1,*)
Jmag = cat(5,*)
write_ds9_regionfile,ra,dec,comment=strtrim(Jmag,2),filename='stdstars.reg',color='red'
radec,ra,dec,rahr,ramin,rasec,decdeg,decmin,decsec

;Write standard star in a format for mask file
openu,lun,'c31mask.txt',/get_lun,/append
for ii = 0L,n_elements(RAhr)-1 do printf,lun,'s'+strtrim(ii+1,2)+' -1 '+strtrim(Jmag(ii),2)+' '+strtrim(RAhr(ii),2)+' '+strtrim(RAmin(ii),2)+' '+strtrim(RAsec(ii),2)+' '+strtrim(Decdeg(ii),2)+' '+strtrim(Decmin(ii),2)+' '+strtrim(Decsec(ii),2)+' 2000.0 2000.0 0 0'
free_lun,lun

end
