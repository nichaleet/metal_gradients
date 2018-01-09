pro make_eps_map_colorbar,outname,name,data,page_width,page_height,colornumber
set_plot, 'ps'
device, filename = outname, $
xsize = page_width,ysize = page_height, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated

minval = min(image(where(finite(image))))
maxval = max(image(where(finite(image))))
data(where(finite(data) eq 0.)) = minval
data = bytscl(data, min = minval, max = maxval)
tam = size(data, /dimensions)
cgloadct, colornumber, ncolors = 256, bottom = 0, clip = [0,256], /reverse
tvlct, redvector, greenvector, bluevector, /get
cgimage,data, position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height] 
device,/close
