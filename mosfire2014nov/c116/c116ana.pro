pro c116ana
setplot,14
img=readfits('c116.fits')
gal = img(*,26:52)
sum = total(gal,2)
meanf = sum/27.
wl = findgen(n_elements(sum))*1.63
plot, wl,meanf,xtitle='angstrom'
stop
end
