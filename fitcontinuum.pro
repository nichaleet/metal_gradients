pro fitcontinuum,filter=filter,filename=filename,filecube=filecube,miny=miny,maxy=maxy
setplot,12
!p.multi = [0,1,2]

image = readfits(filename,header)
;for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
wavelength  = getwl_filter(filter)
nspec       = n_elements(wavelength)

;sum to get spectrum of all pixels in the field
total_spectrum = fltarr(nspec)
sizeim =size(image)
if n_elements(miny) eq 0 then miny=0
if n_elements(maxy) eq 0 then maxy=sizeim(3)-1

for i=0,nspec-1 do total_spectrum[i] = total(image[i,*,miny:maxy])
total_spectrum=total_Spectrum-median(total_spectrum)
plot,wavelength,total_spectrum
;findout where the sky lines and OH are
ok='no'
percents=[7.5,92.5]
while ok ne 'y' do begin
   percents=percents/100.
   plot,wavelength,total_spectrum,yrange=[-2,2]
   badlines = percentiles(total_spectrum,value=[percents(0),percents(1)])
   for i=0,1 do oplot,!x.crange,[badlines(i),badlines(i)],color=100,linestyle=2
   read,ok,prompt='ok?(y/n)'
   if ok ne 'y' then begin
      read,percentlow,prompt='bottom % of data cut off:'
      read,percenttop,prompt='top % of data to cut off:'
      percents=[percentlow,percenttop]
   endif

endwhile
goodwl = where(total_spectrum gt badlines(0) and total_spectrum lt badlines(1))
badwl  = where(total_spectrum le badlines(0) or total_spectrum ge badlines(1))

if filename eq '/scr2/nichal/workspace/reduced_data/mosaic/cswa15_Ha_tlc_Kn2_handmosaic_sky_230hr.fits' then begin
   goodwl = where(total_spectrum gt badlines(0) and total_spectrum lt badlines(1) and wavelength lt 2.105)
   badwl  = where(total_spectrum le badlines(0) or total_spectrum ge badlines(1)or wavelength ge 2.105)

endif

ngoodwl=n_Elements(goodwl)

;calculate variance across the whole field for each wavelength
weight = fltarr(nspec)
for i=0,nspec-1 do weight[i]=1./variance(image[i,*,*])
weight = weight/max(weight)*0.5*max(total_spectrum) ;normallized
weight(badwl) = 1./0.
oplot,wavelength,weight,color=200


imsize = size(image)
continuummap = fltarr(imsize(2),imsize(3))
for ii=0,imsize(2)-1 do for jj=0,imsize(3)-1 do begin
   specreduced = image[goodwl,ii,jj]
   sigmareduced= 1./sqrt(weight(goodwl))
   meanerr,specreduced,sigmareduced,xmean
   continuummap(ii,jj) = xmean*ngoodwl
endfor
;continuummap2 = total(image[goodwl,*,*],1)
;imdisp,continuummap
;writefits,'test.fits',[[[continuummap]],[[continuummap2]]]

writefits,'continuummap.fits',continuummap

nytag='y'
read,nytag,prompt='Do you want to try to see the continuum?'
if nytag eq 'y' then begin
   read,xmin,prompt='xmin'
   read,xmax,prompt='xmax'
   read,ymin,prompt='ymin'
   read,ymax,prompt='ymax'
   total_spectrum = fltarr(nspec)
   sizeim =size(image)
   for i=0,nspec-1 do total_spectrum[i] = total(image[i,xmin:xmax,ymin:ymax])
   plot,wavelength,total_spectrum
   read,conval,prompt='Want to change the value where the peak is (0 if no, otthe values if yes)??'
   if conval ne 0 then continuummap[xmin:xmax,ymin:ymax]=conval
endif

;add to the current acube
blank = 0.*continuummap
cube=readfits(filecube,hdr)
sizecube = size(cube)
if sizecube(3) eq 12 then cube=[[[cube]],[[blank]],[[continuummap]]]
if sizecube(3) eq 13 then cube=[[[cube]],[[continuummap]]]
if sizecube(3) eq 14 then cube[*,*,13] = continuummap
writefits,filecube,cube,hdr

stop
end
