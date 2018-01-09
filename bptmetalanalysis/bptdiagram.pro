pro bptdiagram,name,outdir,n2,n2_err,o3hb,o3hb_err,page_width,page_height,plot_left,plot_bottom,xsize,ysize
;Define greek letters
alphaLetter = '!9' + String("141B) + '!X'  
betaletter  = '!9' + String("142B) + '!X'
;Plot
loadct,0
set_plot,'ps'
psname=outdir+name+'_bpt.eps'
device, filename = psname, $
        xsize = 10, $
        ysize = 8, $
        xoffset = 0, $
        yoffset = 0, $
        scale_factor = 1.0,/encapsulated,/color

;ploterror,n2,o3hb,n2_err[*],o3hb_err[*],type=0,errcolor=240,xtitle='[NII]/H'+alphaletter,ytitle='[OIII]/H'+betaletter,font=0.,psym=1,xrange=[-3,0],yrange=[0,2]
;make a filled symbol
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL
plot,n2,o3hb,xtitle='[NII]/H'+alphaletter,ytitle='[OIII]/H'+betaletter,font=0.,xrange=[-3,0],yrange=[0,2.],psym=3,title=strupcase(name)
oplot,n2,o3hb,psym=8,color=100

loadct,2
;Find out where is significantly in AGN region
o3hb_lowerbound = o3hb-o3hb_err
n2_lowerbound = n2-n2_err
finite_Elements = where(finite(o3hb) and finite(n2))
bptmap = fltarr(size(n2,/dimension))+1./0.
bptmap(where(finite(n2))) = 1. ; values for where the map is  
agnflag = []                    ; to keep the points where it's likely an agn
sigarr = []
for i=0,n_Elements(finite_elements)-1 do begin
   element = finite_elements(i) 
                                ;this is going to skip the "bad points" that makes graph looks bad but these points are not likely an AGN since they are not in the middle
   if name eq 'cswa15' and o3hb(element) gt 1.3 then goto, skip
;   if o3hb_lowerbound(element) gt 0.61/(n2_lowerbound(element)-0.47)+1.19 then begin ;greater than ykewley
   if o3hb(element) gt 0.61/(n2(element)-0.47)+1.19 then begin ;greater than ykewley
      bptmap(element) = .5
      agnflag = [agnflag,element]
      nsig = (o3hb(element)-(0.61/(n2(element)-0.47)+1.19))/o3hb_err(element)
      print, 'significance:',nsig
      sigarr = [sigarr,nsig]
   endif
   skip: continue
endfor
stop
;oplot the errorbar in bpt diagram for the points that is likely AGN
If n_elements(agnflag) gt 0. then begin
;oploterror,n2(agnflag),o3hb(agnflag),n2_err(agnflag),o3hb_err(agnflag),errcolor=100,psym=1
X = [-1, 0, 1, 0, -1]
Y = [0, 1, 0, -1, 0]
USERSYM, X, Y,/fill
oplot,n2(agnflag),o3hb(agnflag),color=80,psym=8
endif

;oplot typical errorbar
oploterror,[-2.2],[0.3],[median(n2_err(where(finite(n2))))],[median(o3hb_err(where(finite(o3hb))))],errcolor=200,thick=2
;stop

;oplot every points with x dots over the color bars
loadct, 12
;oplot,n2,o3hb,psym=7,thick=2
;overplot kewley and steidel
logx=findgen(78)/100.*5.-4.
ykewley = (0.61/(logx-0.47)+1.19) ;kewley01 maximum starburst
ysteidel = (0.67/(logx-0.33)+1.13);locus of z=2
oplot,logx,ykewley,color=15,thick=2
oplot,logx,ysteidel,color=100,thick=2
device,/close

;plot bptmap
writefits,outdir+name+'_bptmap_all.fits',bptmap

;plot bptmap to .eps
psname = outdir+name+'_bptmap.eps'
device, filename = psname, xsize = page_width,ysize = page_height, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated
bptmap(where(finite(bptmap) eq 0.)) = 0.
data = bytscl(bptmap, min = 0., max = 1.)
tam = size(data, /dimensions)
cgloadct, 3, ncolors = 256, bottom = 0, clip = [0,256], /reverse
tvlct, redvector, greenvector, bluevector, /get
cgimage,data, position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height] 
device,/close

end
