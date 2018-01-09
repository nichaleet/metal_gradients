pro annuli_calc,val,val_err,xc,yc,pos,inc,pixelscale,badpix,radius,radius_pixel,val_radius_pixel,val_radius_err_pixel,radius_annuli,val_radius_annuli,val_radius_err_annuli
;calculate the values(val) as a function of radius along annuli
;input  = value (2D map of values),xc yc (center of the galaxy),
;         pos=position angle,inc=inclination of the galaxy
;output = radius_pixel,val_radius_pixel,val_radius_err_pixel,radius_annuli,val_radius_annuli,val_radius_err_annuli

;Program
;1) find values and radius for all pixels
sizeim = size(val,/dimension)
x_dist = rebin(findgen(sizeim(0)),sizeim(0),sizeim(1))-xc
y_dist = rebin(transpose(findgen(sizeim(1))),sizeim(0),sizeim(1))-yc
radius = pixelscale*sqrt((x_dist*cos(pos)-y_dist*sin(pos))^2+(x_dist*sin(pos)+y_dist*cos(pos))^2/(sin(inc))^2)
badpix = [badpix,where(finite(val) eq 0.)]
radius(badpix)=1./0.
good = where(finite(radius))
radius_pixel = radius(good)
val_radius_pixel = val(good)
val_radius_err_pixel = val_err(good)

;2) Find values from annuli 
val_radius_annuli = []
val_radius_err_annuli = []
radius_annuli = []
n_Element_annuli =[]
step_annuli = max(radius_pixel)/pixelscale
For ii=0,step_annuli do begin
   goodpix = where(radius_pixel gt ii*pixelscale and radius_pixel le (ii+1.)*pixelscale)
   if n_elements(goodpix) gt 1. then begin
      radius_annuli = [radius_annuli,(ii+0.5)*pixelscale]
      ;meanerr,val_radius_pixel(goodpix),val_radius_err_pixel(goodpix),wmean,sigmam,sigmad  ;weighted mean
      meanerr,val_radius_pixel(goodpix),1.+fltarr(n_elements(goodpix)),wmean,sigmam,sigmad  ;un-weighted mean

      val_radius_annuli = [val_radius_annuli,wmean]
      val_radius_err_annuli = [val_radius_err_annuli,sigmam]
      n_Element_annuli = [n_element_annuli,n_Elements(goodpix)]
   endif
endfor

end

