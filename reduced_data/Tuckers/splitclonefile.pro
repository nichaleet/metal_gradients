pro splitclonefile

cube=readfits('clone_o3n2ha_acube_source1.fits')
cube2 = readfits('clone_o3n2ha_aerrorcube_source1.fits')
cube = cube[40:120,75:155,*]
cube2 = cube2[40:120,75:155,*]
cube2(where(cube2 eq 0)) = 1./0.
cube(where(cube eq 0)) = 1./0.
name='clone'
writefits,name+'_vel.fits',cube[*,*,0] ;
writefits,name+'_Ha.fits',cube[*,*,1]  ;
writefits,name+'_NII_Ha.fits',cube[*,*,3]/cube[*,*,1] ;

writefits,name+'_velerr.fits',cube2[*,*,0]
writefits,name+'_Haerr.fits',cube2[*,*,1]
writefits,name+'_veldisperr.fits',cube2[*,*,2]

writefits,name+'_NII.fits',cube[*,*,3] ;
writefits,name+'_NII_err.fits',cube2[*,*,3]

writefits,name+'_hadetect.fits',cube[*,*,6] ;
writefits,name+'_veldisp.fits',cube[*,*,2] ;

end
