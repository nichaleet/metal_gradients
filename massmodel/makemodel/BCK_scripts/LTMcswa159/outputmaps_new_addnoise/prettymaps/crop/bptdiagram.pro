pro bptdiagram,name,outdir,ha,n2,n2_err,o3hb,h2_err,o3hb_err,page_width,page_height,plot_left,plot_bottom,xsize,ysize
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

;ploterror,n2,o3hb,n2_err[*],o3hb_err[*],type=3,errcolor=240,xtitle='[NII]/H'+alphaletter,ytitle='[OIII]/H'+betaletter,font=0.,xrange=[1.e-3,0.8],yrange=[1,60.],psym=1
;make a filled symbol
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL
plot,n2,o3hb,xtitle='[NII]/H'+alphaletter,ytitle='[OIII]/H'+betaletter,font=0.,xrange=[1.e-3,0.8],yrange=[1,60.],psym=3,/xlog,/ylog
oplot,n2,o3hb,psym=8,color=100

loadct,2
;Find out where is significantly in AGN region
o3hb_lowerbound = o3hb-o3hb_err
o3hb_lowerbound(where(o3hb lt 0.)) = -1./0.
n2_lowerbound = n2-n2_err
n2_lowerbound(where(n2 lt 0.)) = -1./0.
finite_Elements = where(o3hb gt 0. and n2 gt 0.)
bptmap = fltarr(size(ha,/dimension))+1./0.
bptmap(where(finite(ha))) = 1. ; values for where the map is  
agnflag = []                    ; to keep the points where it's likely an agn
for i=0,n_Elements(finite_elements)-1 do begin
   element = finite_elements(i)
   if o3hb(element) gt 10.^(0.61/(alog10(n2(element))-0.47)+1.19) then bptmap(element)= 1.
   if o3hb_lowerbound(element) gt 10.^(0.61/(alog10(n2_lowerbound(element))-0.47)+1.19) then begin ;greater than ysteidel
      bptmap(element) = .5
      agnflag = [agnflag,element]
   endif
endfor
;oplot the errorbar in bpt diagram for the points that is likely AGN
If n_elements(agnflag) gt 0. then begin
                                ;oploterror,n2(agnflag),o3hb(agnflag),n2_err(agnflag),o3hb_err(agnflag),errcolor=100,psym=1
X = [-1, 0, 1, 0, -1]
Y = [0, 1, 0, -1, 0]
USERSYM, X, Y,/fill
oplot,n2(agnflag),o3hb(agnflag),color=80,psym=8
endif

;oplot typical errorbar
oploterror,[0.05],[9],[median(n2_err)*.75],[median(o3hb_err)],errcolor=200,thick=2

;oplot every points with x dots over the color bars
loadct, 12
;oplot,n2,o3hb,psym=7,thick=2
;overplot kewley and steidel
logx=findgen(78)/100.*5.-4.
ykewley = 10.^( 0.61/(logx-0.47)+1.19) ;kewley01 maximum starburst
ysteidel = 10.^(0.67/(logx-0.33)+1.13);locus of z=2
oplot,10.^logx,ykewley,color=15,thick=2
oplot,10.^logx,ysteidel,color=100,thick=2
device,/close

;plot bptmap
writefits,outdir+name+'_bptmap.fits',bptmap

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
