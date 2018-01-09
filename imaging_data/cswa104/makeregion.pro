pro makeregion


;Write objects of which we want to get spectra.
;Numbers are in the F105 catalog
;READ IN CATALOT
cat = RSEX('cswa104_F105W.cat')
obj1 =[113,114] ;Images of the arcs
obj2 =[123] ; Image of the arc but less important
obj3 = [168,153,154,156,155,179,160,219,259,299] ;cluster members (I guessed)
obj4 = [260,238,247,246]                     ;lensed background objs
obj5 = [136,181,285,66,538,325,362,448,183,387,394,395,357,38,39,198,139,101,107,103] ;random objects that seems interesting

listobj=list(obj1,obj2,obj3,obj4,obj5)
nameobj=['arcm','arc','cl','bg','rd']
importantvalue = [500,400,200,100,30]

openw,lun,'c104mask.txt',/get_lun
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


;stdstsars.txt is generated from http://www2.keck.hawaii.edu/software/fcTool/findingChart.php

;readcol,'stdstars.txt',F='(A,I,I,F,I,I,F,F,A,A,A,A,A)',v1,Rahr,RAmin,RAsec,Decdeg,Decmin,Decsec,V2,V3,V4,mag,v4,v5
;mag = strsplit(mag,'=',/extract)
;mag = mag.ToArray()
;mag = mag(*,1)
;RA = RAhr*15.+RAmin*15./60.+RAsec*15./3600.
;Dec = Decdeg+Decmin/60.+Decsec/3600.
;write_ds9_regionfile,ra,dec,comment=mag,filename='stdstars.reg',color='red'

;2mass point source from the catalog in ds9
cat = read_ascii('SDSSptsource.tsv',data_start=1,header=header)
cat = cat.field01
ra = cat(0,*)
dec= cat(1,*)
Rmag = cat(9,*)
write_ds9_regionfile,ra,dec,comment=strtrim(Rmag,2),filename='stdstars.reg',color='red'
radec,ra,dec,rahr,ramin,rasec,decdeg,decmin,decsec

;Write standard star in a format for mask file
openu,lun,'c104mask.txt',/get_lun,/append
for ii = 0L,n_elements(RAhr)-1 do printf,lun,'s'+strtrim(ii+1,2)+' -1 '+strtrim(Rmag(ii),2)+' '+strtrim(RAhr(ii),2)+' '+strtrim(RAmin(ii),2)+' '+strtrim(RAsec(ii),2)+' '+strtrim(Decdeg(ii),2)+' '+strtrim(Decmin(ii),2)+' '+strtrim(Decsec(ii),2)+' 2000.0 2000.0 0 0'
free_lun,lun

end
