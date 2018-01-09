pro upresolution
;Make the original ESI resolution of 0.15"/pix to match the Osiris 0.1"/pix
;Input: cswa15_B_small_match.fits,cswa15_R_small_match.fits,cswa15_I_small_match.fits
;Output: cswa15_B_small_match_osirisRES.fits,cswa15_R_small_match_osirisRES.fits,cswa15_I_small_match_osirisRES.fits

im_B = readfits('cswa15_B_small_match.fits',hdB)
im_R = readfits('cswa15_R_small_match.fits',hdR)
im_I = readfits('cswa15_I_small_match.fits',hdI)

oldsize = size(im_B)
hcongrid,im_B,hdB,out=[oldsize(1)*1.5,oldsize(2)*1.5],cubic=-0.5,/half
hcongrid,im_R,hdR,out=[oldsize(1)*1.5,oldsize(2)*1.5],cubic=-0.5,/half
hcongrid,im_I,hdI,out=[oldsize(1)*1.5,oldsize(2)*1.5],cubic=-0.5,/half

writefits,'cswa15_B_small_match_osirisRES.fits',im_B,hdB
writefits,'cswa15_R_small_match_osirisRES.fits',im_R,hdB
writefits,'cswa15_I_small_match_osirisRES.fits',im_I,hdB

end
