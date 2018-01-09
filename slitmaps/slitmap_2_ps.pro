pro slitmap_2_ps
; define page size (in cm)
page_height = 10
page_width = 10

; define position of plot
; bottom left corner (in cm)
plot_left = 1
plot_bottom = 1

; define plot size (in cm)
xsize = 7
ysize = 7

; use postscript output 
set_plot, 'ps'
; name the output file 
psname = 'plot.ps'

; open the postscript file
device, filename = psname, $
xsize = page_width, $
ysize = page_height, $
xoffset = 0, $
yoffset = 0, $
scale_factor = 1.0, $
/portrait 

; read FITS file
image = mrdfits('cswa31_slitmap.fits',0, header)
image = image(*,*,1)
; select a sub array
;image = image[100:800,100:800]
  
; scale the image
data = 255b - bytscl(image, min = -80., max = 40.)

tam = size(data, /dimensions)

stop
; load colour table
cgloadct, 33, ncolors = 256, bottom = 0, clip = [0,256], /reverse
tvlct, redvector, greenvector, bluevector, /get
   
; display background image   
cgimage, data, $
position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height] 

; load colour table
cgloadct, 0
tvlct, redvector, greenvector, bluevector, /get
   
; display axes   
scale=0.628 ;kpc per pixel
cgplot, [0], [0], $
xcharsize = 1, ycharsize = 1, $
thick = 2, $
xrange= [0, tam[0]]*scale, $
yrange= [0, tam[1]]*scale, $
xtitle = 'kpc',$
ytitle = 'kpc',$
xstyle = 1, ystyle = 1, $
/nodata, /normal, /noerase, $
position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height] 

; define coordinates of sources
xc = [14]
yc = [23]

; display black X at sources 
oplot, xc, yc, color = cgcolor("white"), psym = symcat(7), symsize = 2

; label two of the sources
;xyouts, 125.00, 397.0, '1', charsize = 1.0, color = cgcolor("white") 
;xyouts, 205.00, 440.0, '2', charsize = 1.0, color = cgcolor("white") 

; draw a line to indicate scale (values provided are fictitious)
;oplot, [680,680], [150,300], linestyle = 0, thick = 3, color = cgcolor("white") 
;xyouts, 675, 200, ' 1Kpc', charsize = 1., color = cgcolor("white"), orientation = 90 
 
 
; load colour table for colour bar
cgloadct, 33, ncolors = 256, bottom = 0, clip = [0,256], /reverse
tvlct, redvector, greenvector, bluevector, /get

; define position of colour bar   
plot_left = 1
plot_bottom = 8
xsize = 7
ysize = 0.5

; define colour bar 
cgcolorbar, divisions = 5, $
range = [max(image), min(image)], $
xminor = 1, $
xticklen = 0.1, format = '(I3)', $
palette= [[redvector],[greenvector],[bluevector]], $
position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height] ,$
annotatecolor = 'white', $
charsize = 0.8
  
; close postscript file
device,/close

; create a PDF file
cgps2pdf, psname, /delete_ps 


end
