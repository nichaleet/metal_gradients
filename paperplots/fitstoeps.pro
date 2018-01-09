pro fitstoeps,file,scale,xtitle,ytitle,colorbartitle,directory_out,format=format
  file_mkdir,directory_out
  nameout = strmid(file,strpos(file,'/',/reverse_search)+1,strpos(file,'.fits'))+'.eps'  
  nameout = directory_out+'/'+nameout
  print, 'output is',nameout
  set_plot,'ps'
  device,filename=nameout,/encapsulated
  device,decomposed=1,color=1,bits_per_pixel=8
  device,xsize=6.5,ysize=7,/inches
  device,set_font='Helvetica',/TT_FONT

  loadct,39
  img = readfits(file)  
  
  size = size(img)
  xc = findgen(size[1])*scale
  yc = findgen(size[2])*scale
     
  good = img(where(finite(img) eq 1.))
  ;reject outliers
  sd = stdev(good)
  meanval = mean(good)
  outliers = where(img le meanval-5.*sd or img ge meanval+5.*sd and finite(img) eq 1.)
  print,'there values are considered outliers:',img(outliers)
  img(outliers) = 1./0.
  good = img(where(finite(img) eq 1.))
  stat = summary(good)      
                      
  !p.charsize=1.5

  ;for NII/Ha
  if colorbartitle eq '[NII]/Ha' then begin
     if stat(0) le -0.5 then stat(0) = -0.5
     if stat(4) ge 0.6 then stat(4) = 0.6
  endif
  ;for metal
  if colorbartitle eq '12+logO/H(M08)' then stat(0) = 7.

  range=[stat(0),stat(4)]

  rdisplay,img,xc,yc,range=range,xtitle=xtitle,ytitle=ytitle,ymargin=[5,7] 

  !p.charsize=1.5
  colourbar, range=range,position=[0.2, 0.85, 0.95, 0.9],title=colorbartitle,format=format,font=1
  device,/close
  set_plot,'x'
  !p.background = 255
  !p.color=0
  rdisplay,img,xc,yc,range=range,xtitle=xtitle,ytitle=ytitle,ymargin=[5,6] ;,maskevalue=9999.
  !p.charsize=1.5
  colourbar, range=range,position=[0.2, 0.85, 0.95, 0.9],title=colorbartitle,format=format,font=1
wait,2
end
