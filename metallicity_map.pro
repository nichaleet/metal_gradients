pro metallicity_map,type=type,name=name,path=path, Ha_file=Ha_file, OIII_file=OIII_file
common BPTshare, goodagn
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

;calculate N2metallicity
N2metal,halpha,halpha_err,NII,NII_err,N2metal,N2metal_err,upper_tags

PA_spec    = sxpar(ha_hdr,'PA_SPEC')
scale      = sxpar(ha_hdr,'scale') ; in degrees

;0000000000000000000000000000000000000000000000000
if type eq "N2" then begin  ;Making a fits file for N2 metallicity and BAyesian method only
;---------------------------------------------------------
;calculate Bayesian method metallicity
   ;Prepare blank arrays for Bayesian method(below)
   sizeha = size(halpha)
   new_OIII = fltarr(sizeha(1),sizeha(2))*1./0.
   new_hbeta= fltarr(sizeha(1),sizeha(2))*1./0.
   new_OIII_err= fltarr(sizeha(1),sizeha(2))*1./0.
   new_hbeta_err= fltarr(sizeha(1),sizeha(2))*1./0.
   Bayesian_metal,name, new_OIII,new_hbeta,NII,Halpha,new_OIII_err,new_hbeta_err,NII_err,halpha_err,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,typebayesian
   Bayesian_metal_err = (upperbound_metal_arr-lowerbound_metal_arr)/2.
;-----------------------------------------------------
;writing fits file
   N2fits  = [[[N2metal]],[[best_metal_arr]],[[N2metal_err]],[[Bayesian_metal_err]],[[upper_tags]]] 
   n2_header = ha_hdr
   hdr_para_names = ['frame1','frame2','frame3','frame4','frame5']
   hdr_para_vals = [1/0,1/0,1/0,1/0,1/0]
   comment = ['N2metal','Bayesian[NII]/Ha metal','N2metal Error','Bayesian metal error','1=detection,2=upper limit detection, 3=no detection']
   for ii=0,n_elements(hdr_para_names)-1 do sxaddpar,n2_header,hdr_para_names(ii),hdr_para_vals(ii),comment(ii)
   writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_N2_metallicity.fits',N2fits,n2_header
   
endif

;==============================================================
;2 O3N2 metallicity
;0000000000000000000000000000000000000000000000000
if type eq 'O3N2' then begin
   OIIIcube = readfits(OIII_file,OIII_hdr)

;--------------------------------------------
;* Align pixels from 2 sets. Align OIII to Halpha

   writefits,'test.fits',oiiicube[*,*,1],oiii_hdr
   oiii=readfits('test.fits',oiii1d_hdr)

   hbeta     = oiiicube[*,*,6]
   oiii_err  = oiiicube[*,*,9]
   hbeta_err = oiiicube[*,*,11] 
   velocity = oiiicube[*,*,0]
   vel_disp = oiiicube[*,*,2]
   bins     = oiiicube[*,*,3]
   oiiisig  = oiiicube[*,*,4]
   redshift = oiiicube[*,*,5]
   hbsig    = oiiicube[*,*,7]
   velsig   = oiiicube[*,*,8]
   dispsig  = oiiicube[*,*,10]

   hastrom,oiii,oiii1d_hdr,new_oiii,new_oiii_hdr,ha_hdr;,missing=1./0.
   hastrom,hbeta,oiii1d_hdr,new_hbeta,newhdr,ha_hdr,missing=1./0.
   hastrom,oiii_err,oiii1d_hdr,new_oiii_err,newhdr,ha_hdr,missing=1./0.
   hastrom,hbeta_err,oiii1d_hdr,new_hbeta_err,newhdr,ha_hdr,missing=1./0.
   hastrom,velocity,oiii1d_hdr,new_velocity,newhdr,ha_hdr,missing=1./0.
   hastrom,vel_disp,oiii1d_hdr,new_vel_disp,newhdr,ha_hdr,missing=1./0.
   hastrom,bins,oiii1d_hdr,new_bins,newhdr,ha_hdr,missing=1./0.
   hastrom,oiiisig,oiii1d_hdr,new_oiiisig,newhdr,ha_hdr,missing=1./0.
   hastrom,redshift,oiii1d_hdr,new_redshift,newhdr,ha_hdr,missing=1./0.
   hastrom,hbsig,oiii1d_hdr,new_hbsig,newhdr,ha_hdr,missing=1./0.
   hastrom,velsig,oiii1d_hdr,new_velsig,newhdr,ha_hdr,missing=1./0.
   hastrom,dispsig,oiii1d_hdr,new_dispsig,newhdr,ha_hdr,missing=1./0.
;Now all the new_param is aligned with the Halpha image

   new_OIIIcube = [[[new_velocity]],[[new_oiii]],[[new_vel_disp]],[[new_bins]],[[new_oiiisig]],[[new_redshift]],[[new_hbeta]],[[new_hbsig]],[[new_velsig]],[[new_oiii_err]],[[new_dispsig]],[[new_hbeta_err]]]

   namenewoiii = strmid(oiii_File,0,strpos(oiii_file,'.fits'))+'_aligned.fits'
   writefits,namenewoiii,new_OIIIcube,new_oiii_hdr
;This new_OIIIcube has the same size/position as Halpha image.The new fits file is Hbeta's original name+_aligned.fits

;--------------------------------------------------------------
;*calculate O3N2 metallicity
   O3N2metal,new_OIII,new_hbeta,NII,Halpha,new_OIII_err,new_hbeta_err,NII_err,halpha_err,O3N2metal,O3N2metal_err

;--------------------------------------------------------------
;*calculate metallicity with Bayesian method
   Bayesian_metal,name, new_OIII,new_hbeta,NII,Halpha,new_OIII_err,new_hbeta_err,NII_err,halpha_err,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,typebayesian
   Bayesian_metal_err = (upperbound_metal_arr-lowerbound_metal_arr)/2.
;--------------------------------------------------------------
;*making BPT plot
  bptanalysis,name, new_OIII,new_hbeta,NII,Halpha,new_OIII_err,new_hbeta_err,NII_err,halpha_err

;----------------------------------------------------------------
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


writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_metallicity.fits',O3N2_fits,o3n2_header


endif

end




