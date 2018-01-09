pro convert_hubble, fileref


hubbleim = 'cswa28_814_small_osiris_source.tif'
outname  = repstr(hubbleim,'tif','fits')
ref = readfits(fileref,header)


img = read_tiff(hubbleim)

fits = readfits('cswa28_814_small_osiris.fits')
slope = max(fits)/max(img)
print, slope

img = img*slope
img(where(img eq 0.)) = 1./0.
writefits,outname,img,header
end
