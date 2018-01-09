pro measure_psf_sourceplane

setplot,14
!p.multi=[0,1,2]

cswa11={name:'cswa11',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa11/outputmaps_secondgal/prettymaps/crop/sourceha_psf_pretty.fits',scale:0.626}

cswa15={name:'cswa15',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa15/LTMzoomESI/outputmaps/prettymaps/crop/sourceha_psf_pretty.fits',scale:0.422}
cswa19={name:'cswa19',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMredshift2fixedtry2/outputmaps_handmosaic/prettymaps/crop/sourceha_psf_pretty.fits',scale:0.167}

cswa20={name:'cswa20',file:'/scr2/nichal/workspace/output/chopped_sourceplane/cswa20/cswa20_sourceha_psf.fits',scale:0.172}

cswa20Hb={name:'cswa20Hb',file:'/scr2/nichal/workspace/output/chopped_sourceplane/cswa20/cswa20_sourceoiii_psf.fits',scale:0.172}

cswa28={name:'cswa28',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMredshift2fixed_cswa28/outputmaps_handmosaic/prettymaps/crop/sourceha_psf_pretty.fits',scale:0.424}

cswa31 = {name:'cswa31',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa31/outputmaps/prettymaps/crop/sourceha_psf_pretty.fits',scale:0.628}

cswa128 = {name:'cswa128',file:'/scr2/nichal/workspace/output/chopped_sourceplane/cswa128/cswa128_sourceha_psf.fits',scale:0.336}

cswa139 = {name:'cswa139',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/cswa139zoomESI/outputmaps/prettymaps/crop/sourceha_psf_pretty.fits',scale:0.40975}

cswa139Hb = {name:'cswa139Hb',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/cswa139zoomESI/outputmaps/prettymaps/crop/sourceoiii_psf_pretty.fits',scale:0.40975}
;IDL is stupid doesn't work here. I used matlab and it gave 0.09754x0.2785 pixels for Sqrt(2)*sigma--> FWHM=w/sqrt(2)*2.3548*0.40975 = 0.07*0.19 kpc

cswa159 = {name:'cswa159',file:'/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa159/outputmaps/prettymaps/crop/sourceha_psf_pretty.fits',scale:0.83641}
;Fit is bad. Guess from projection of raw image is about 3.3x1.5 pix=
cswa165 = {name:'cswa165',file:'/scr2/nichal/workspace/output/chopped_sourceplane/cswa165/cswa165_sourceha_psf.fits',scale:0.106}

cswa165Hb = {name:'cswa165Hb',file:'/scr2/nichal/workspace/output/chopped_sourceplane/cswa165/cswa165_sourceoiii_psf.fits',scale:0.106}

a773 = {name:'abell773',file:'/scr2/nichal/workspace/output/chopped_sourceplane/Abell773/Abell773_sourceha_psf.fits',scale:0.159}

filearr=[cswa11,cswa15,cswa19,cswa20,cswa20Hb,cswa28,cswa31,cswa128,cswa165,cswa165Hb,cswa139,cswa139Hb,cswa159,a773]

ngal = n_elements(filearr)
namearr    = strarr(ngal)
psf_fwhm_x  = fltarr(ngal)
psf_fwhm_y  = fltarr(ngal)
psf_theta = fltarr(ngal)


for i=0,ngal-1 do begin
window,0,xsize=500,ysize=800
   str  = filearr(i)
   file = str.file
   name = str.name
   scale = str.scale ;kpc per pix
   image= readfits(file,hdr)
   image(where(finite(image) eq 0)) = 0.
   psf=gauss2dfit(image,param,/tilt)

   print, param
   rdisplay,image,title=name
   psf(where(finite(psf) eq 0)) = param(0)
   rdisplay,psf
   print, name
   print, 'x_width (kpc) , y_width(kpc), rotation: ',param(2)*scale,scale*param(3),param(6)
   namearr(i)    = name
   psf_fwhm_x(i)  = param(2)*scale*2.3548
   psf_fwhm_y(i)  = param(3)*scale*2.3548
   psf_theta(i) = param(6)   
   writefits,'psf_sourceplane/'+name+'_psf_fit.fits',psf
   wait,1
endfor
save, namearr, psf_FWHM_x,psf_fwhm_y,psf_theta,filename='/scr2/nichal/workspace/output/psf/psf_variables_sourceplane.sav'


print, 'name, fwhm_x(kpc), fwhm_y(kpc),fwhm_x(pix),fwhm_y(pix), theta'
for i=0,n_elements(namearr)-1 do print, namearr(i),psf_FWHM_x(i),psf_fwhm_y(i),psf_FWHM_x(i)/scale,psf_fwhm_y(i)/scale,psf_theta(i)

stop
end
