pro upresolution
;Make the original ESI resolution of 0.1542"/pix to match the Osiris 0.1"/pix
;Input: cswa139_ESI_V.fits cswa139_ESI_R.fits cswa139_ESI_I.fits
;Output: cswa139_ESI_V_match_osirisRES.fits cswa139_ESI_R_match_osirisRES.fits cswa139_ESI_I_match_osirisRES.fits 
;Oldsize is 301*301 --> new size is 464*464

im_V = readfits('cswa139_ESI_V.fits',hdV)
im_R = readfits('cswa139_ESI_R.fits',hdR)
im_I = readfits('cswa139_ESI_I.fits',hdI)

oldsize = size(im_V)
hcongrid,im_V,hdV,out=[464,464],cubic=-0.5,/half
hcongrid,im_R,hdR,out=[464,464],cubic=-0.5,/half
hcongrid,im_I,hdI,out=[464,464],cubic=-0.5,/half

writefits,'cswa139_ESI_V_match_osirisRES.fits',im_V,hdV
writefits,'cswa139_ESI_R_match_osirisRES.fits',im_R,hdV
writefits,'cswa139_ESI_I_match_osirisRES.fits',im_I,hdV

end
