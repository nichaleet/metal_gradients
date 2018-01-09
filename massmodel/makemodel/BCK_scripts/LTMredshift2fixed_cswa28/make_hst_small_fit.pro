pro make_hst_small_fit
newim = readfits('cswa28_814_small.fits',header)


newim(0:347,*) = 0.
newim(*,0:578) = 0.
newim(390:1000,*) = 0.
newim(*,641:1000) = 0.

interest = where(newim)
ind = array_indices(newim,interest)
xind = ind(0,*)
yind = ind(1,*)

line1 = where(yind gt (xind-348.)*40./31.+595.)  
line2 = where(yind lt (xind-363.)*4./3.+578.)
bad = [line1,line2]
newim(interest(bad)) = 0.

writefits,'cswa28_814_small_osiris.fits',newim,header
stop

end
