pro metalnew,type=type,name=name,path=path, Ha_file=Ha_file, OIII_file=OIII_file

;===============================================================
Ha_file   = path+Ha_file
OIII_file = path+OIII_file

;1) N2 metallicity
;From Pettini and Pagel, 2004
hacube = readfits(Ha_file,Ha_hdr)
halpha     = hacube[ *,*,1]
halpha_err = hacube[*,*,9]
NII        = hacube[*,*,6]
NII_err    = hacube[*,*,11]

N2 = NII/Halpha
N2_err = sqrt((NII_err/NII)^2+(halpha_err/halpha)^2)

PA_spec    = sxpar(ha_hdr,'PA_SPEC')
scale      = sxpar(ha_hdr,'scale') ; in degrees

;Now just to get the N2
N2metal = N2
O3N2metal = N2*0.
best_metal_arr=N2*0.
N2metal_err = N2_err
O3N2metal_err = N2*0.
Bayesian_metal_err = N2*0.
typebayesian = N2*0.

;*Combine all results into a fits image
   O3N2_fits = [[[N2metal]],[[O3N2metal]],[[best_metal_arr]],[[N2metal_err]],[[O3N2metal_err]],[[Bayesian_metal_err]],[[typebayesian]]]

;*writing fits file
o3n2_header = ha_hdr
hdr_para_names = ['frame1','frame2','frame3','frame4','frame5','frame6','frame7']
hdr_para_vals = [1/0,1/0,1/0,1/0,1/0,1/0,1/0]
comment1 = ['N2metal','O3N2metal','Bayesian metal','N2metal Error','O3N2metal-N2metal error','Bayesian metal error','Values in Bayesian method.[2,4,8] for NII''s [Non detection, upper limit detection, detection] [3,9,27] for OIII''s']
 
for ii=0,n_elements(hdr_para_names)-1 do begin 
   sxaddpar,o3n2_header,hdr_para_names(ii),hdr_para_vals(ii),comment1(ii)
endfor


writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_metallicitynew.fits',O3N2_fits,o3n2_header



end
