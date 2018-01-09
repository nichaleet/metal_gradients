pro run_fitspec_intensity,input=input,filter=filter,line,xmins,xmaxs,ymins,ymaxs,smooth,header,imgymins,imgymaxs

finalcube = readfits(input)

if filter eq 'Jbb' then wl =1.180+0.00015*findgen(1574)else $
if filter eq 'Jn2' then wl =1.228+0.00015*findgen(408)else $
if filter eq 'Hbb' then wl =1.473+0.0002*findgen(1651)else $
if filter eq 'Hn1' then wl =1.466+0.0002*findgen(376)else $
if filter eq 'Hn2' then wl =1.532+0.0002*findgen(391)else $
if filter eq 'Hn3' then wl =1.594+0.0002*findgen(411)else $
if filter eq 'Hn4' then wl =1.652+0.0002*findgen(426)else $
if filter eq 'Kn1' then wl =1.955+0.00025*findgen(401)else $
if filter eq 'Kn2' then wl =2.036+0.00025*findgen(421)else $
if filter eq 'Kc3' then wl =2.121+0.00025*findgen(433)else begin 
   print, 'no filter match'   
   stop
endelse


;Clean finalcube (set the outliers to 0.)
;clean_treshold = 10.*stdev(finalcube[where(abs(finalcube) le 1.)]) ;5sigma
clean_treshold = 1.
bad_pixels = where(abs(finalcube) gt clean_treshold)
finalcube[bad_pixels] = 0.

;smooth the data
if smooth gt 1 then fwhm = smooth-1. else fwhm = 1.
for i=0,n_elements(finalcube[*,0,0])-1 do finalcube[i,*,*] = filter_image(reform(finalcube[i,*,*]),smooth=smooth,fwhm_gaussian=fwhm)
name_pos = strpos(input,'/',/reverse_search)
name = strmid(input, name_pos+1)
name = strmid(name, 0,strlen(name)-5)  ;delete '.fits'

acube_fit = readfits('output/'+name+'_acube.fits')
z_cube    = acube_fit[*,*,5]
disp_cube = (acube_fit[*,*,2] / 2.998d5) * z_cube * line   ; convert velocity dispersion from km/s to microns
bin_cube  = acube_fit[*,*,3]

fitspec_intensity, finalcube, wl,line,z_cube,disp_cube,bin_cube,xmins,xmaxs,ymins,ymaxs,name=name,header,imgymins,imgymaxs

end
