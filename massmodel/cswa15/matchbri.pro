pro matchBRI
fits = readfits('cswa15_B_small.fits',header)
size = size(fits)
backgroundval = median(fits)
addrow = fltarr(size[1])+backgroundval
fits = [[addrow],[fits]]
fits = transpose(fits)
addrow = fltarr(size[2]+1)+backgroundval
fits = [[addrow],[fits]]
fits = transpose(fits)
sxaddpar,header,'BZERO',0
writefits,'cswa15_B_small_match.fits',fits,header
writefits,'checkmatch.fits',fits,header

fits = readfits('cswa15_I_small.fits',header)
size = size(fits)
backgroundval = median(fits)
addrow = fltarr(size[1])+backgroundval
fits = [[fits],[addrow]]
fits = transpose(fits)
addrow = fltarr(size[2]+1)+backgroundval
fits = [[fits],[addrow]]
fits = transpose(fits)
sxaddpar,header,'BZERO',0
writefits,'cswa15_I_small_match.fits',fits,header

fits = readfits('cswa15_R_small.fits',header)
size = size(fits)
backgroundval = median(fits)
addrow = fltarr(size[1])+backgroundval
fits = [[fits],[addrow]]
fits = transpose(fits)
addrow = fltarr(size[2]+1)+backgroundval
fits = [[fits],[addrow]]
fits = transpose(fits)
sxaddpar,header,'BZERO',0
writefits,'cswa15_R_small_match.fits',fits,header

end
