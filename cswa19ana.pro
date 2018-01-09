function twogauss,wl,a
;pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
  const = a[0]
  z     = a[1]
  vdisp = a[2]
  Hhb   = a[3]
  Ho3   = a[4]
  clight = 299792.458
  hbwl   = 486.269*(z+1.)
  o3wl   = 500.824*(z+1.)
  sigmahb = vdisp/clight*hbwl
  sigmao3  = vdisp/clight*o3wl
  specstruct = const+Hhb*exp(-0.5*(wl-hbwl)^2/sigmahb^2)+Ho3*exp(-0.5*(wl-o3wl)^2/sigmao3^2)
  return,specstruct
end

pro cswa19ana
clight = 299792.458
setplot, 14
!p.multi = [0,2,2]
z=2.30
;1) do H alpha image
;loadct, 22,ncolors=100

image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa19_Ha_Kn1_mosaic_scaledsky_4hr.fits",head)
image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa19_Ha_tlc_Kn1_pipelinemosaic_scaledsky_3hr.fits",head)


;get header parameters

RA     = sxpar(head,'RA')  ;RA at spatial [0,0] in mosaic ;format='(f14.10)'
DEC    = sxpar(head,'DEC');DEC at spatial [0,0] in mosaic ;format='(f14.10)'
PA_Spec = sxpar(head,'PA_SPEC') ;position angle of spectrograph on sky
PA_IMAG = sxpar(head,'PA_IMAG') ; position angle of imager on sky
NAXIS   = sxpar(head,'NAXIS*')

CRVAL1 = sxpar(head, 'CRVAL1') ;[nm] Wavelength at reference pixel 
CRVAL2 = sxpar(head, 'CRVAL2') ;[deg] R.A. at reference pixel       
CRVAL3 = sxpar(head, 'CRVAL3') ;[deg] DEC at reference pixel       

CRPIX1 = sxpar(head, 'CRPIX1') ; Reference pixel location        
CRPIX2 = sxpar(head, 'CRPIX2') ; Reference pixel location        
CRPIX3 = sxpar(head, 'CRPIX3')  ;Reference pixel location       

CDELT1 = sxpar(head, 'CDELT1') ; Wavelength scale in nm/channel
CDELT2 = sxpar(head, 'CDELT2') ; Pixel scale in degree/pixel
CDELT3 = sxpar(head, 'CDELT3') ; Pixel scale in degree/pixel
PC2_2   =  sxpar(head, 'PC2_2')  ;/RA, Dec axes rotated by 205.000 degr.           
PC2_3    =  sxpar(head, 'PC2_3')  ;/RA, Dec axes rotated by 205.000 degr.           
PC3_2    =  sxpar(head, 'PC3_2')  ;/RA, Dec axes rotated by 205.000 degr.           
PC3_3   =  sxpar(head, 'PC3_3')  ;/RA, Dec axes rotated by 205.000 degr.           

for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[140:151,*,*],1)
;imdisp, ha_sum,position=[.05,.05,.5,1.], out_pos=pos1
imdisp, ha_sum, out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.03,0.04,0.05,0.08]

;make fit file of the ha_map
ha_sum = transpose(ha_sum)
mkhdr,ha_sum_hdr,ha_sum
hdr_paranames=['CDELT1','CDELT2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','PC1_1','PC1_2','PC2_1','PC2_2']
hdr_val1 = [CDELT3,CDELT2,CRPIX3,CRPIX2,CRVAL2,CRVAL3,PC2_2,PC2_3,PC3_2,PC3_3]
hdr_paraname_str = ['CTYPE1','CTYPE2','CUNIT1','CUNIT2']
hdr_paraname_str_val = ['RA---TAN','DEC--TAN','deg','deg']

for ii=0,n_elements(hdr_paranames)-1 do begin 
   sxaddpar,ha_sum_hdr,hdr_paranames(ii),hdr_val1(ii)
endfor
for ii=0,n_elements(hdr_paraname_Str)-1 do begin 
   sxaddpar,ha_sum_hdr,hdr_paraname_str(ii),hdr_paraname_str_val(ii)
endfor
writefits,'cswa19.fits',ha_sum,ha_sum_hdr


;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
cube = readfits('/scr2/nichal/workspace/output/cswa19_Ha_tlc_Kn1_pipelinemosaic_scaledsky_3hr_acube.fits')
ha=cube[*,*,0]
goodha=where(finite(ha))
for i =0, n_elements(image[*,0,0])-1 do begin
   im_i = image[i,*,*]
   ha_spectrum[i] = total(im_i(goodha))
endfor
;for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,25:45,25:35])
wavelength = findgen(401)*0.00025+1.955
plot,wavelength[120:180],ha_spectrum[120:180],psym=1,/noerase

;find the peak of the ha_spectrum
wavepeak = wavelength(where(ha_spectrum eq max(ha_spectrum)))
z= wavepeak/656.4-1.
vline,(z[0]+1.)*[656.4,658.5]
Print, 'redshift is ',z
print, where(ha_spectrum eq max(ha_spectrum))


;fit Halpha spectrum
good = where(wavelength gt 1.980 and wavelength lt 2.0)
wl = wavelength(good)
ha = ha_spectrum(good)
plot, wl,ha,psym=10,xtitle='micron',ytitle = 'flux',/noerase,yrange=[-0.5,3.1],position=[.1,0.05,.45,.45],title='Halpha'
fit = gaussfit(wl,ha,param,nterms=4,sigma=sigma)
oplot, wl,fit,color=50,thick =2
redshift =  param(1)/.656461-1.
veldisp = param(2)/param(1)*3.e5
intrinsic_Veldisp = sqrt(veldisp^2-50.^2)
print,'parameter a,b,c,d:', param
print, redshift,veldisp,intrinsic_veldisp
areaHa = abs(param(0)*param(2))*sqrt(2.*!pi)
error_areaHa = areaHa*sqrt((sigma(0)/param(0))^2+(sigma(2)/param(2))^2)*sqrt(2.*!pi)

restore,'/scr2/nichal/workspace/Calibrate_Ha/flux_calibration.sav'
pos1 = where(flux_calibration.name eq 'cswa19')
lamb_res = sxpar(head,'CDELT1')/1000.
conversionfactor = flux_calibration(pos1).calibfactor
conversionfactor_err = flux_calibration(pos1).calibfactor_err
totHaflux = areaHa*conversionfactor/lamb_res ;erg/s/cm^2
tothaflux_err = totHaflux*sqrt(error_areaHa^2/areaHa^2+conversionfactor_err^2/conversionfactor^2) ;erg/s/cm^2
print,'total Ha flux ',totHaflux,'erg/s/cm^2',tothaflux_err
stop

;fit NII spectrum
weight = wavelength*0.
for i=0, n_Elements(weight)-1 do weight[i] = 1./variance(image[i,55:60,25:45])
good = where(wavelength gt 1.993 and wavelength lt 2.0)
wl = wavelength(Good)
NII = ha_spectrum(good)
weight = weight(good)
print, .658523*(redshift+1.)
fit2ndline,param(2),.658523*(redshift+1.),wl,NII,weight,fit,area,sigma_area
oplot,wl,fit,color=20,thick=2

plot, wl,NII,psym=10,xtitle='micron',ytitle = 'flux',/noerase,position=[.6,0.05,.95,.45],title='NII'
oplot,wl,fit,color=20,thick=2
areaNII = area
error_areaNII = sigma_area

print, 'log(NII/Ha) = ',alog10(areaNII/areaHa)
print, 'error = ',0.4343*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)

print, '(NII/Ha) = ',(areaNII/areaHa)
print, 'error = ',(areaNII/areaHa)*sqrt((error_areaHa/areaHa)^2+(error_areaNII/areaNII)^2)

;OIII Hb
image2 = readfits('/scr2/nichal/workspace/reduced_data/mosaic/cswa19_Hb_tlc_Hn1_handmosaic_sky_1hr.fits',head2)
;ALIGNED IMAGE WITH HA
imsize1 = size(image,/dimensions)
imsize2 = size(image2,/dimensions)

;To plot 1D spec Take only galaxy region. 
path='/scr2/nichal/workspace/output/'
detection= readfits(path+'cswa19_Hb_tlc_Hn1_handmosaic_sky_1hr_acube.fits',head2d2)
gooddetect = where(finite(detection[*,*,0]))
spec1d = fltarr(n_elements(image2[*,0,0]))
spec1derr = spec1d
for i =0, n_elements(spec1d)-1 do begin
   im_i = image2[i,*,*]
   spec1d[i] = total(im_i(gooddetect))
   spec1derr[i] = n_elements(gooddetect)*stdev(im_i(gooddetect))
endfor
wl = getwl_filter('Hn1')*1000.
!p.multi=[0,1,1]
plot,wl,spec1d
vline,3.034*[500.824,486.269]

;fit gaussian to two lines
znow = 2.034
pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)
pi.value = double([0.0,znow,200.0,1,1])
pi[0].limits = [-0.5,0.5]
pi[1].limits = [2.03,2.04]
pi[2].limits = [0.0, 1000.0]
pi[3].limits = minmax(spec1d)
pi[4].limits = minmax(spec1d)
pi.step = double([0.1, 0.001, 25.0,0.001,0.001])
pi.parname = ['    const', '  z', 'vdisp','peakval1','peakval2']
pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)']
mask = bytarr(n_elements(wl))+1
mask(where(wl gt 1503 and wl lt 1510)) = 0
won = where(mask)
params = MPFITFUN('twogauss',wl(won),spec1d(won),spec1derr(won),parinfo=pi,/nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit)
model = twogauss(wl,params)
oplot,wl,model,color=150
;calculate equivalent width
z      = params[1]
model_normal = model/params[0]
spec_normal  = spec1d/params[0]
hbwl   = 486.269*(z+1.)
o3wl   = 500.824*(z+1.)
vdisp  = params[2]
range  = [-1000,1000]/clight*o3wl+o3wl
won    = where(wl gt range[0] and wl lt range[1])
o3_eqw1 = int_tabulated(wl(won),1.-model_normal(won)) ;nm
range  = [-2.*vdisp,2*vdisp]/clight*o3wl+o3wl
won    = where(wl gt range[0] and wl lt range[1])
o3_eqw2 = int_tabulated(wl(won),1.-spec_normal(won)) ;nm

range  = [-1000,1000]/clight*hbwl+hbwl
won    = where(wl gt range[0] and wl lt range[1])
hb_eqw1 = int_tabulated(wl(won),1.-model_normal(won)) ;nm
range  = [-2.*vdisp,2*vdisp]/clight*hbwl+hbwl
won    = where(wl gt range[0] and wl lt range[1])
hb_eqw2 = int_tabulated(wl(won),1.-spec_normal(won)) ;nm

print,'EW OIII:',o3_eqw1,'+-',abs(o3_eqw1-o3_eqw2)
print,'EW Hb:',hb_eqw1,'+-',abs(hb_eqw1-hb_eqw2)

stop
end
