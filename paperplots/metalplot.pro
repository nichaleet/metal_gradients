pro metalplot

files = file_search('/scr2/nichal/workspace/output/cropmetal/*.fits')
print, files
 
scale = [0.16,0.336,0.106,0.167,0.17,0.425]   ;kpc per pixel  ;[0.019",0.04",]
;;(abell773,128,165,19,20,28)
name =['Abell773','CSWA128','CSWA165','CSWA19','CSWA20','CSWA28']

;setplot, 5
;for i=0, n_elements(files)-1 do begin
;img = readfits(files(i))
;vel = img(where(finite(img) eq 1.))
;size = size(img)
;xc = findgen(size[1])*scale(i)
;yc = findgen(size[2])*scale(i)
;img(where(finite(img) eq 0.)) = max(vel)*3.
;stat = summary(vel)

;!p.charsize=2.
;rdisplay,img,xc,yc,range=[stat(0),stat(4)],xtitle='kpc',ytitle=name(i),ymargin=[5,6]

;!p.charsize=1.5
;colourbar, range=[stat(0),stat(4)],position=[0.15, 0.9, 0.95, 0.95],title='velocity (km/s)'
;stop

;endfor

nameout ='metal_'+ name+'.eps'

for i=0, n_elements(files)-1 do begin
   set_plot,'ps'
   device,filename=nameout(i),/encapsulated
   device,decomposed=1,color=1,bits_per_pixel=8
   device,xsize=6.5,ysize=7,/inches
   loadct,39
   img = readfits(files(i))
   vel = img(where(finite(img) eq 1.))
   size = size(img)
   xc = findgen(size[1])*scale(i)
   yc = findgen(size[2])*scale(i)
   ;img(where(finite(img) eq 0.)) = 9999.
   stat = summary(vel)
   !p.charsize=1.5
   range=[stat(0),stat(4)]
   
   rdisplay,img,xc,yc,range=range,xtitle='kpc',ytitle=name(i),ymargin=[5,7];,maskevalue=9999.
   !p.charsize=1.5
   colourbar, range=range,position=[0.2, 0.85, 0.95, 0.9],title='12+log(O/H)(N2)'
   device,/close
set_plot,'x'
!p.background = 255
!p.color=0
rdisplay,img,xc,yc,range=range,xtitle='kpc',ytitle=name(i),ymargin=[5,6];,maskevalue=9999.
   !p.charsize=1.5
   colourbar, range=range,position=[0.2, 0.85, 0.95, 0.9],title='12+log(O/H)'
;wait,2
   
endfor
end
