pro calculate_fluxconv,filestarcubes,name=name,cluster=cluster,filters=filters,photmags=mags,mags_err=mags_err,z=z,linewl=linewl,outthres=outthres,calibfactors=output,fixpeak=fixpeak,fixcenter=fixcenter,centerpos=centerpos,stoptag=stoptag
;Explanation:
  ;calculate the flux calibration for Ha and Hb using standard stars.
  ;since most of the tip-tilt stars weren't observed in Hb band, we use standard star to calibrate both Ha and Hb for reddening.
clight = 299792.458 ;km
clightcgs = clight*1.e5
hplanckcgs = 6.6260755e-27
nlines = n_elements(filestarcubes)
output = replicate({filestarcube:'',filter:'',starname:'',clustername:'',mag:0.,z:0.,linewl:0.,conv_factor:0.,conv_factorerr:0.,unitto:'erg/s/cm^2',unitfrom:'ADU/s/keck_area',gain:0.},nlines)

for i=0,nlines-1 do begin
   starcube = readfits(filestarcubes[i],hdr)
   starcubeerr = sqrt(abs(starcube))
   wl       = getwl_filter(filters[i])
   sizecube = size(starcube,/dimensions)

;remove outliers
   outlier = where(abs(starcube) gt outthres[i] or ~finite(Starcube))
   if outlier[0] ne -1 then begin
      starcube(outlier) = median(starcube(where(starcube ne 0 and finite(starcube))))
      starcubeerr(outlier) = median(starcubeerr(where(starcube ne 0 and finite(starcube))))
   endif

;collapse the cube in wavelength to see the psf
   star_2d = total(starcube,1)*0.
   star_2d_sigma = star_2d
   imsize = size(star_2d,/dimensions)
   for ii=0,imsize(0)-1 do begin
      for jj=0,imsize(1)-1 do begin
         meanclip,starcube[*,ii,jj],meanflux,sigma
         star_2d(ii,jj) = meanflux
         star_2d_sigma(ii,jj) = sigma
      endfor
   endfor

;Find stars
   ;find peak location in a smoothed image
   peakval = max(star_2d,peakloc)
   xarr = rebin(findgen(imsize(0)),imsize(0),imsize(1))
   yarr = rebin(transpose(findgen(imsize(1))),imsize(0),imsize(1))
   xc0 = xarr(peakloc)
   yc0 = yarr(peakloc)
   if keyword_Set(centerpos) then xc0=centerpos[0]
   if keyword_set(centerpos) then yc0=centerpos[1]
;Fit 2d Gaussian
   int_param = [0.,peakval,2.,2.,xc0,yc0,0]
   pi = replicate({fixed:0,limits:[0.,0.]},7)
   if keyword_set(fixpeak) then pi[1].fixed = 1   
   if keyword_set(fixcenter) then pi[4:5].fixed=1
   pi[4].limits=[xc0-5,xc0+5]
   pi[5].limits=[yc0-5,yc0+5]
   psf = mpfit2dpeak(star_2d,param,error=star_2d_sigma,estimates=int_param,parinfo=pi,status=status)
   xc = param[4]
   yc = param[5]
   print,'position of star:',xc,yc
   print,'FWHM(pixels):',param([2,3])*2.35482
   FWHM = round(param([2,3])*2.35482)
   FWTM = round(param([2,3])*4.29193)
   threesigma = round(param([2,3])*3.)
   print,'FWTM(pixels):',param([2,3])*4.29193,FWTM

;make star region mask
   xarr = xarr-xc
   yarr = yarr-yc
   ;Elipse equation: total(distance from the point to each focus) = 2a
   a = float(max(threesigma))
   b = float(min(threesigma))
   c = sqrt(a^2-b^2)
   if FWTM(0) gt FWTM(1) then rarr = sqrt((xarr+c)^2+yarr^2)+sqrt((xarr-c)^2+yarr^2) $
   else rarr = round(sqrt((yarr+c)^2+xarr^2)+sqrt((yarr-c)^2+xarr^2))

   mask = bytarr(imsize(0),imsize(1))
   mask(where(rarr le 2.*a)) = 1
   mask2d = mask
   mask3d = fltarr(sizecube[0],sizecube[1],sizecube[2])
   for jj=0,sizecube[0]-1 do mask3d[jj,*,*]=mask
   if outlier[0] ne -1 then mask3d(outlier) = 0

;make model star with gaussian peak
   sigmax = param[2]
   sigmay = param[3]
   modelflux = param[0]+param[1]*exp(-1.*(xarr^2/(2*sigmax^2)+yarr^2/(2*sigmay^2)))

;Plot
   window,0,xsize=800,ysize=700,title=name
   cgLoadCT, 33, CLIP=[30,255]
   !p.multi=[0,2,2]
   !p.font = 0
   position =   [0.1, 0.1, 0.9, 0.9]
   cgImage, star_2d, Stretch=1, MinValue=median(star_2d), MaxValue=max(modelflux)/10,/Axes, XTitle='x',YTitle='y',XRange=minmax(xarr), YRange=minmax(yarr),font=0,position=position
   tvellipse,FWHM[0]/2,FWHM[1]/2,0,0,/data
   tvellipse,FWTM[0]/2,FWTM[1]/2,0,0,/data
   tvellipse,threesigma[0],threesigma[1],0,0,/data

   cgImage, modelflux,Stretch=1, MinValue=median(star_2d), MaxValue=max(modelflux)/10.,/Axes, XTitle='x',YTitle='y',XRange=minmax(xarr), YRange=minmax(yarr),font=0,position=[0.1,0.1,0.9,0.9]
   tvellipse,FWHM[0]/2,FWHM[1]/2,0,0,/data
   tvellipse,FWTM[0]/2,FWTM[1]/2,0,0,/data
   tvellipse,threesigma[0],threesigma[1],0,0,/data
   writefits,'outfluxconv/'+cluster+'_'+name+'_'+filters[i]+'.fits',[[[star_2d]],[[modelflux]]]

;Calculate 1D spectra
   starcube_cleaned = starcube
   starcubeerr_cleaned = starcubeerr
   starcube_cleaned = starcube_cleaned-param[0] ;adjust floor level to zero
   starcube_cleaned(where(~mask3d)) = 0 
   starcubeerr_cleaned(where(~mask3d)) = 0
   spec1d    = total(total(starcube_cleaned,3),2)
   spec1derr = sqrt(total(total(starcubeerr_cleaned^2,3),2))
   cgPlot,wl,spec1d,xtitle='wavelength',ytitle='SPEC1D',color=fsc_color('black'),title='Nicha Reduction'
   oplot,wl,spec1d-spec1derr,color=fsc_color('rosy brown')
   oplot,wl,spec1d+spec1derr,color=fsc_color('rosy brown')
   vline,(z+1.)*linewl[i]
;search for odrfgui reduced spec1d
   file_odrf = strmid(filestarcubes[i],0,strpos(filestarcubes[i],'_mosaic'))
   if file_odrf eq '' then file_odrf = strmid(filestarcubes[i],0,strpos(filestarcubes[i],'_tlc'))
   file_odrf = file_search(file_odrf+'*1d*.fits')
   if file_odrf ne '' then begin
      odrf_1d = readfits(file_odrf[0])
      if total(finite(odrf_1d)) ne 0 then plot,wl,odrf_1d,title='pipeline reduction',color=fsc_color('black')
   endif

;Get flux at wavelength in num_res resolution elements
      ;The unit in reduced data cube is DN/sec
      ;The spectral resolution is about 3000 -> ~100 km/s
   num_Res  = 3.
   linewlnow= (z+1.)*linewl[i]
   wl_pix   = sxpar(hdr,'CDELT1')/1000. ;microns
   if wl_pix-(wl(1)-wl(0)) gt 1.e-10 then stop,'Stop: Problems with wavelength resolution.'
   pix_res = num_Res*round(linewlnow*100./clight/wl_pix) ;number of pixels per one resolution element
   if pix_res mod 2. ne 1. then pix_res = pix_res+1.
   minval = min(abs(wl-linewlnow),pos_line)
   pos_min = pos_line-0.5*(pix_res-1.)
   pos_max = pos_line+0.5*(pix_res-1.)
   print,'Collapsing flux in pixel [min,line,max]:',pos_min,pos_line,pos_max
   ;collapse spec1D in wavelength
   star_lineflux    = total(spec1d[pos_min:pos_max])
   star_linefluxerr = sqrt(total(spec1derr[pos_min:pos_max]))
   star_lineflux    = star_lineflux/(pix_res*wl_pix) ;ADU/second/micron
   star_linefluxerr = star_linefluxerr/(pix_res*wl_pix) ;ADU/second/micron

;Calculate Real Flux from stars
   ;THIS IS FOR FAINT STANDARD ONLY https://www2.keck.hawaii.edu/inst/nirc/UKIRTstds.html ;constant = -48.58 
   mag       = mags[i]
   magerr    = mags_err[i]
   fnu           = 10.^((mag+48.58)/(-2.5)) ;erg/s/cm^2/Hz
   fnuerr        = abs(fnu*alog(10.)*magerr/2.5)
   clight_micron = clight*1.e9
   flambda       = clight_micron/linewlnow^2*fnu ;erg/s/cm^2/micron 
   calibfactor = flambda/star_lineflux ;erg/s/cm^2 per ADU/second/telescope area
   calibfactorerr = calibfactor*sqrt((star_linefluxerr/star_lineflux)^2+(fnuerr/fnu)^2)
   gain = hplanckcgs*clightcgs/calibfactor/(linewlnow*1.e-4)

   outnow = output[i]
   outnow.filestarcube = filestarcubes[i]
   outnow.filter   = filters[i]
   outnow.starname = name
   outnow.clustername = cluster
   outnow.mag = mag
   outnow.z = z
   outnow.linewl = linewl[i]
   outnow.conv_factor = calibfactor
   outnow.conv_factorerr = calibfactorerr
   outnow.gain = gain
   output[i] = outnow
   help,outnow,/str
   wait,1
   if keyword_set(stoptag) and stoptag eq 1 then stop
endfor

end
