pro sfr_cal

Haflux = [49,31,18,3,5,14,14,1.6,9,0.4,0.6,1.3,16,7,18] ;10^-17erg/s/cm2
z = [1.41,2.16,2.03,1.43,2.09,1.49,2.22,2.54,2.30,2.13,2.30,2.21,2.2,2.38,2.0]
lumdist = lumdist(z) ;mpc
lumdist = lumdist*1.d6*3.08567758d18; cm
SFR = 7.9d-42*4.*!dpi*lumdist^2*Haflux*1.e-17

print,SFR
stop
end
