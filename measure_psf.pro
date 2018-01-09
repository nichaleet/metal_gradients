pro measure_psf ; get the psf parameters of tt star for each arc at Halpha wavelength in image plane
setplot,14
!p.multi=[0,2,2]
smoothwidth = 2.
cswa11={z:1.41,name:'cswa11',file:'/scr2/nichal/WMKO/MAR2013/SPEC/ORP/s130304_a008001_Hn2_100.fits.gz',filter:'Hn2'}
cswa15={z:2.16,name:'cswa15',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140224_a018001_Kn2_100.fits',filter:'Kn2'}
cswa19={z:2.03,name:'cswa19',file:'/scr2/nichal/WMKO/MAR2013/SPEC/ORP/s130304_a015001_Kn1_100.fits.gz',filter:'Kn1'}
cswa20={z:1.43,name:'cswa20',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140224_a020001_Hn2_100.fits',filter:'Hn2'}
cswa20Hb={z:1.43,name:'cswa20Hb',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140223_a025001_Jn1_100.fits',filter:'Jn1'}
cswa28={z:2.09,name:'cswa28',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140223_a022001_Kn1_100.fits',filter:'Kn1'}
cswa31={z:1.49,name:'cswa31',file:'/scr2/nichal/WMKO/DEC2014/SPEC/ORP/s141206_a014001_Hn3_100.fits',filter:'Hn3'}
cswa31_1Hb={z:1.49,name:'cswa31_1Hb',file:'/scr2/nichal/WMKO/DEC2014/SPEC/ORP/s141206_a019001_Jbb_100.fits',filter:'Jbb'}
cswa31_2Hb={z:1.49,name:'cswa31_2Hb',file:'/scr2/nichal/WMKO/DEC2014/SPEC/ORP/s141207_a023002_Jbb_100.fits',filter:'Jbb'}
cswa128={z:2.22,name:'cswa128',file:'/scr2/nichal/WMKO/SEP2013/SPEC/ORP/s130911_a020001_Kn2_100.fits.gz',filter:'Kn2'}
cswa139={z:2.54,name:'cswa139',file:'/scr2/nichal/WMKO/DEC2014/SPEC/ORP/s141207_a015002_Kc5_100.fits',filter:'Kc5'}
cswa139febHb={z:2.54,name:'cswa139Hb',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140224_a013001_Hbb_100.fits',filter:'Hbb'}
cswa159={z:2.30,name:'cswa159',file:'/scr2/nichal/WMKO/DEC2014/SPEC/ORP/s141207_a005001_Kc3_100.fits',filter:'Kc3'}
cswa165={z:2.13,name:'cswa165',file:'/scr2/nichal/WMKO/DEC2014/SPEC/ORP/s141207_a010002_Kn2_100.fits',filter:'Kn2'}
cswa165Sep={z:2.13,name:'cswa165Sep',file:'/scr2/nichal/WMKO/SEP2013/SPEC/ORP/s130911_a028001_Kn2_100.fits.gz',filter:'Kn2'}
cswa165Hb={z:2.13,name:'cswa165Hb',file:'/scr2/nichal/WMKO/DEC2014/SPEC/ORP/s141207_a013002_Hbb_100.fits',filter:'Hbb'}
A773={z:2.30,name:'Abell773',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140222_a005001_Kc3_100.fits',filter:'Kc3'}
A773Hb={z:2.30,name:'Abell773Hb',file:'/scr2/nichal/WMKO/MAR2013/SPEC/ORP/s130303_a005001_Hn3_100.fits.gz',filter:'Hn3'}
M0717={z:2.55,name:'MACS0717',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140224_a004001_Kc5_100.fits',filter:'Kc5'}
M0717Hb={z:2.55,name:'MACS0717Hb',file:'/scr2/nichal/WMKO/FEB2014/SPEC/ORP/s140224_a007001_Hbb_100.fits',filter:'Hbb'}
M1133={z:1.55,name:'MACS1133',file:'/scr2/nichal/WMKO/MAR2013/SPEC/ORP/s130303_a015001_Hn4_100.fits.gz',filter:'Hn4'}
M1133Hb={z:1.55,name:'MACS1133Hb',file:'/scr2/nichal/WMKO/MAR2013/SPEC/ORP/s130304_a020001_Jn2_100.fits.gz',filter:'Jn2'}

filearr=[cswa11,cswa15,cswa19,cswa20,cswa20Hb,cswa28,cswa31,cswa31_1Hb,cswa31_2Hb,cswa128,cswa139,cswa139febHb,cswa159,cswa165,cswa165Sep,cswa165Hb,A773,A773Hb,M0717,M0717Hb,M1133,M1133Hb]
;filearr=[cswa139,cswa139febHb,cswa159,cswa165,cswa165Sep,cswa165Hb,A773,A773Hb,M0717,M0717Hb,M1133,M1133Hb]

ngal = n_elements(filearr)
namearr    = strarr(ngal)
psf_width  = fltarr(ngal)
psf_height = fltarr(ngal)
psf_floor  = fltarr(ngal)

for i=0,ngal-1 do begin
window,0,xsize=900,ysize=900
   str  = filearr(i)
   file = str.file
   name = str.name
   z    = str.z
   filter= str.filter

   wlarr = getwl_filter(filter)
   cube= readfits(file,hdr)
   dif = abs(wlarr-(.6564*(z+1.))) ;need to find the wl that corresponds to Ha
   if strmid(name,1,3,/reverse_offset) eq 'Hb' then  dif = abs(wlarr-(.5008*(z+1.))) ; for OIII line

   pos = where(dif eq min(dif))
   ;cutout glitches/outliers
   resistant_mean,cube,10.,meancube
   cube(where(abs(cube) ge 15.)) = meancube

   image = mean(cube(pos-10:pos+10,*,*),dimension=1)
   image = filter_image(image,fwhm_Gaussian=smoothwidth)
   if name eq 'cswa20Hb' then image = -image ;somehow the negative line is stronger

   ;cutout the negative image all pixels below 3xfirst quartile
   stat = summary(image)
   maskbar = stat(1)*3.
   if maskbar gt 0. then maskbar=-0.01
   maskedpix = where(image lt stat(1))
   mask = abs(image/image)
   mask(maskedpix)= 0.

   if name eq 'cswa139' then begin
      mask[*,36:50]=0      
   endif

   if name eq 'cswa15' then begin
      param = [0.,0.04,2.4,2.4,30,22,0]
      image(where(abs(image) gt 0.5)) = meancube
      ;image = filter_image(image,fwhm_Gaussian=3)
   endif 

   if name eq 'cswa20Hb' then mask[27:60,*]=0
  
   writefits,'/scr2/nichal/workspace/output/psf/'+name+'_tt_psf.fits',image

   sizecube = size(image)
   x = findgen(sizecube(1))
   y = findgen(sizecube(2))   

   if name eq 'cswa15' then begin
      ;psf=gauss2dfit(image,param,x,y,fita=fita,mask=mask,/negative)
      psf=gauss2dfit(image,param,x,y,fita = [1,1,1,1,0,0,1])
   endif else psf=gauss2dfit(image,param,mask=mask)

   if name eq 'cswa139' then begin ; we will use the gaussian fit for the 1D fit
      xflux = image[*,param(5)]
      xflux(where(xflux lt 0.))=0.
      pix   = findgen(n_elements(xflux))
      fit   = gaussfit(pix,xflux,param1d,nterms=4)
      dim   = size(image)
      npix  = min([dim(1),dim(2)])
      sigma = param1d(2)
      psfsq   = psf_gaussian(npixel=npix,st_dev=sigma,centroid=[param(4),param(5)],/normalize)
      sigma=sigma(0)
      psfsq  = 0.2*2.*!pi*sigma^2*psfsq
      psf  = fltarr(dim(1),dim(2))
      psf[0:npix-1,0:npix-1] = psfsq
      param(2)=param1d(2)
      param(3)=param1d(2)
   endif 

   print, param
   psf = psf-param(0)
   rdisplay,image,title=name
   psf(where(finite(psf) eq 0)) = param(0)
   rdisplay,psf
   print, name
   print, 'x_width, y_width, rotation: ',param(2),param(3),param(6)
   namearr(i)    = name
   psf_width(i)  = 0.5*(abs(param(2))+abs(param(3)))
   psf_height(i) = abs(param(1))
   psf_floor(i)  = param(0)
   

   ;Plot 1D gaussian
   xflux = image[*,param(5)]
   yflux = image[param(4),*]
   xfit = psf[*,param(5)]
   yfit = psf[param(4),*]
   plot,xflux,psym=10
   oplot,xfit,color=50
   oplot,!x.crange,[maskbar,maskbar],linestyle=2
   plot, yflux,psym=10
   oplot,yfit,color=50
   oplot,!x.crange,[maskbar,maskbar],linestyle=2
   if name eq 'cswa139' then stop
   wait,1
endfor


;according to WIKI, FWHM ~ 2.35482*width param
psf_FWHM = 2.35482*psf_width
;subtract the smooth param
psf_FWHM = sqrt(psf_FWHM^2-smoothwidth^2)
save, namearr, psf_FWHM,psf_height,psf_floor,filename='/scr2/nichal/workspace/output/psf/psf_variables.sav'

print,'name, fwhm, floor value'
for i=0,n_elements(namearr)-1 do print, namearr(i),psf_FWHM(i),psf_floor(i)
stop
end

