pro metalnew,type=type,name=name,path=path, Ha_file=Ha_file, OIII_file=OIII_file,outpath=outpath,silent=silent

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
N2_err = abs(N2)*sqrt((NII_err/NII)^2+(halpha_err/halpha)^2)
;stop

PA_spec    = sxpar(ha_hdr,'PA_SPEC')
scale      = sxpar(ha_hdr,'scale') ; in degrees


if type eq 'N2' then begin
;calculate metallicity by Bayesian method
   new_OIII = N2*0.+1./0.       ;In case want only N2 index in the analysis
   new_Hbeta = new_OIII
   new_OIII_err = new_OIII
   new_hbeta_err = new_OIII

   Bayesian_metal,name, new_OIII,new_hbeta,NII,Halpha,new_OIII_err,new_hbeta_err,NII_err,halpha_err,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,typebayesian,silent=silent
   Bayesian_metal_err = (upperbound_metal_arr-lowerbound_metal_arr)/2.

   sizefits = size(N2)
   O3N2 = fltarr(sizefits(1),sizefits(2))*0.+1./0.
   O3N2_err = O3N2
   best_metal_arr = best_metal_arr
   Bayesian_metal_err = (upperbound_metal_arr-lowerbound_metal_arr)/2.

;*Combine all results into a fits image
   N2_fits = [[[N2]],[[O3N2]],[[best_metal_arr]],[[N2_err]],[[O3N2_err]],[[Bayesian_metal_err]],[[typebayesian]]]

;*writing fits file
header = ha_hdr
hdr_para_names = ['frame1','frame2','frame3','frame4','frame5','frame6','frame7']
hdr_para_vals = [1/0,1/0,1/0,1/0,1/0,1/0,1/0]
comment1 = ['[NII]/Ha','[OIII]/Hb/[NII]/Ha','Bayesian metal','[NII]/Ha Error','OIII/Hb/NII/Ha error','Bayesian metal error','Values in Bayesian method.[2,4,8] for NII''s [Non detection, upper limit detection, detection] [3,9,27] for OIII''s']
 
for ii=0,n_elements(hdr_para_names)-1 do begin 
   sxaddpar,header,hdr_para_names(ii),hdr_para_vals(ii),comment1(ii)
endfor

if keyword_set(outpath) then outfile = outpath+name+'_metallicitynew.fits' else outfile = '/scr2/nichal/workspace/output/metallicity/'+name+'_metallicitynew.fits'
writefits, outfile,N2_fits,header

endif


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

   hastrom,oiii,oiii1d_hdr,new_oiii,new_oiii_hdr,ha_hdr,missing=1./0.
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
;*calculate metallicity with Bayesian method
new_OIII = new_OIII
   Bayesian_metal,name, new_OIII,new_hbeta,NII,Halpha,new_OIII_err,new_hbeta_err,NII_err,halpha_err,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,typebayesian
   Bayesian_metal_err = (upperbound_metal_arr-lowerbound_metal_arr)/2.

N2metal = N2
O3N2metal = (new_OIII/new_hbeta)/N2
N2metal_err = N2_err
O3N2metal_err = abs(O3N2metal)*sqrt((new_OIII_err/new_OIII)^2+(new_hbeta_err/new_hbeta)^2+(NII_err/NII)^2+(halpha_err/halpha)^2)


;*Combine all results into a fits image
   O3N2_fits = [[[N2metal]],[[O3N2metal]],[[best_metal_arr]],[[N2metal_err]],[[O3N2metal_err]],[[Bayesian_metal_err]],[[typebayesian]]]

;*writing fits file
o3n2_header = ha_hdr
hdr_para_names = ['frame1','frame2','frame3','frame4','frame5','frame6','frame7']
hdr_para_vals = [1/0,1/0,1/0,1/0,1/0,1/0,1/0]
comment1 = ['[NII]/Ha','[OIII]/Hb/[NII]/Ha','Bayesian metal','[NII]/Ha Error','OIII/Hb/NII/Ha error','Bayesian metal error','Values in Bayesian method.[2,4,8] for NII''s [Non detection, upper limit detection, detection] [3,9,27] for OIII''s']
 
for ii=0,n_elements(hdr_para_names)-1 do begin 
   sxaddpar,o3n2_header,hdr_para_names(ii),hdr_para_vals(ii),comment1(ii)
endfor

if keyword_set(outpath) then outfile = outpath+name+'_metallicitynew.fits' else outfile = '/scr2/nichal/workspace/output/metallicity/'+name+'_metallicitynew.fits'
writefits, outfile,O3N2_fits,o3n2_header
stop
endif

end
