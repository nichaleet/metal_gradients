pro run_fitspec_kmos,input=input,filter=filter,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,run_fitspec_int_index,secondline,fwhm_gaussian,xmin,xmax,ymin,ymax,crpix2,crpix3

Common cube, acube,aerrorcube


incube = readfits(input,header)

if filter eq 'Jbb' then wlarr =1.180+0.00015*findgen(1574)else $
if filter eq 'Jn1' then wlarr =1.174+0.00014987*findgen(388)else $
if filter eq 'Jn2' then wlarr =1.228+0.00015*findgen(408)else $
if filter eq 'Hbb' then wlarr =1.473+0.0002*findgen(1651)else $
if filter eq 'Hn1' then wlarr =1.466+0.0002*findgen(376)else $
if filter eq 'Hn2' then wlarr =1.532+0.0002*findgen(391)else $
if filter eq 'Hn3' then wlarr =1.594+0.0002*findgen(411)else $
if filter eq 'Hn4' then wlarr =1.652+0.0002*findgen(426)else $
if filter eq 'Kn1' then wlarr =1.955+0.00025*findgen(401)else $
if filter eq 'Kn2' then wlarr =2.036+0.00025*findgen(421)else $
if filter eq 'Kc5' then wlarr =2.292+0.00025*findgen(465)else $
if filter eq 'Kc3' then wlarr =2.121+0.00025*findgen(433)else begin 
   print, 'no filter match'   
   stop
endelse

name_pos = strpos(input,'/',/reverse_search)
name = strmid(input, name_pos+1)
name = strmid(name, 0,strlen(name)-5)  ;delete '.fits'
name = name+'_KMOS'

;Clean finalcube (set the outliers to 0.)
clean_treshold = 1.
bad_pixels = where(abs(incube) gt clean_treshold)
if bad_pixels(0) ne -1 then begin
   incube[bad_pixels] = 0.
   ind = array_indices(incube,bad_pixels)
endif

;smooth the data
if fwhm_gaussian ne 0 then for i=0,n_elements(incube[*,0,0])-1 do incube[i,*,*] = filter_image(reform(incube[i,*,*]),fwhm_gaussian=fwhm_gaussian)
writefits,'checkcube.fits',incube,header

;delete the last row/column if size is odd number
sizecube=size(incube,/dimensions)
sizex= sizecube(1)
sizey= sizecube(2)
if sizex mod 2. eq 1. then incube=incube[*,0:sizex-2,*]
if sizey mod 2. eq 1. then incube=incube[*,*,0:sizey-2]
sizecube=size(incube,/dimensions)
;rebin into 0.2"X0.2" pixels
newsizex=sizecube(1)*0.5
newsizey=sizecube(2)*0.5
newcube=findgen(sizecube(0),newsizex,newsizey)
for i=0,sizecube(0)-1 do  newcube[i,*,*] = rebin(reform(incube[i,*,*]),newsizex,newsizey)
writefits,'/scr2/nichal/workspace/reduced_data/mosaic/'+name+'.fits',newcube,header

stop
finalcube = newcube
wl=wlarr
nmax=0.

fitspec_general,finalcube,wl,line,z,sigmalim,xmins,xmaxs,ymins,ymaxs,0,1,nmax,300.,200.,name=name,xmin,xmax,ymin,ymax,header,header_acube,crpix2,crpix3
;fitspec_general, finalcube, wl, line, z, sigmalim, xmins, xmaxs, ymins, ymaxs, nmin, ni, nmax, vmax, wmax, name=name,xmin,xmax,ymin,ymax,header,header_acube,crpix2,crpix3
;stop
if run_fitspec_int_index eq 1 then begin
   z_cube    = acube[*,*,5]
   disp_cube = (acube[*,*,2] / 2.998d5) * z_cube * line ; convert velocity dispersion from km/s to microns
   bin_cube  = acube[*,*,3]
   wl=wlarr
   fitspec_intensity_kmos, newcube, wl,secondline,z_cube,disp_cube,bin_cube,xmins,xmaxs,ymins,ymaxs,name=name,header,xmin,xmax,ymin,ymax
endif

output_cube = [[[acube]],[[aerrorcube]]]
;for cswa128
if name eq 'cswa128_Ha_tlc_Kn2_handmosaic_sky_130hr_KMOS' then begin
output_cube[0:18,12:13,*] = 1./0.
output_Cube[13,14,*]=1./0.
output_Cube[22,12,*]=1./0.
endif
cdelt1=sxpar(header_acube,'CDELT1')
sxaddpar,header_acube,'CDELT1',cdelt1*2.
sxaddpar,header_acube,'CDELT2',cdelt1*2.
sxaddpar,header_acube,'SCALE',cdelt1*2.
sxaddpar,header_acube,'CD1_1',sxpar(header_acube,'CD1_1')*2.
sxaddpar,header_acube,'CD1_2',sxpar(header_acube,'CD1_2')*2.
sxaddpar,header_acube,'CD2_1',sxpar(header_acube,'CD2_1')*2.
sxaddpar,header_acube,'CD2_2',sxpar(header_acube,'CD2_2')*2.
acubefile='/scr2/nichal/workspace/output/'+name+'_acube.fits'
writefits, acubefile, output_cube,header_acube

;make point spread function center on the center of the field. This will be used in the model of the kinematics
sizecube=size(output_cube)
sizex=sizecube(1)
sizey=sizecube(2)
sizepsf=min(sizex,sizey)
fwhm=3 ;0.6" = 3 pixels
array = PSF_GAUSSIAN( Npixel=sizepsf, FWHM=[fwhm,fwhm], /NORMAL )
arraypsf=extend_array(array,sizex,sizey)

writefits,'/scr2/nichal/workspace/output/psf_files/'+name+'_psf.fits',arraypsf,header_acube


;Analyze 2nd line detection
seclineflux = output_cube[*,*,6]
uncertainty  = output_cube[*,*,11]
goodpix = where(finite(seclineflux))
negative_pix = where(seclineflux lt 0.)
total_2ndline = total(seclineFlux(goodpix))
print, 'total of 2nd line flux: ',total_2ndline
negative_upperlim = where(seclineflux+uncertainty lt 0.)
plothist, seclineflux(goodpix),title='NII flux',/autobin
print, 'Out of', n_Elements(goodpix),'pixels,',n_elements(negative_pix),'are negative.', n_Elements(negative_upperlim),' have negative upper limit.'
end

