pro run_tucker
  names  = ['0744','1148','1038','clone']
  fwhm   = [[.8,.3],[.9,.6],[1.6,0.4],[3.,0.5]] ;kpc
  theta0 = [118.,70.,300.,350.]
  cenxy  =[[[22.,49,35,45]+15.],[40,82,50,70.],[[20,50,74,100]+20.],[32,47,35,48]]
  pixarcsec = [0.005,0.01,0.01,0.02] ;arcsec per pixel
  z = [2.21,2.38,2.20,2.00]
  angle_psf=[350.,315.,160.,210.]
  
  angdist = lumdist(z)/(1.+z)^2  ;mpc
  angsize = 1./3600.*!dpi/180.*1.e3*angdist ;kpc per arcsec
  pixscale = pixarcsec*angsize ;kpc per pixel
  fwhm = fwhm/rebin(transpose(pixscale),2,n_elements(names)) ;pixel
  
  path =  '/scr2/nichal/workspace/reduced_data/Tuckers/'

  for i =3,3 do begin ;n_elements(names)-1 do begin
     name=names[i]
     velfile= path+name+'_vel.fits'
     velerrfile = path+name+'_velerr.fits'
     if file_test(velerrfile) ne 1 then begin
        print,'uncertainty file does not exist'
        vel = readfits(velfile)
        temp = vel
        temp(where(finite(temp))) = 5.
        if name eq '0744' then begin
           temp(where(vel lt -200.)) = 100.
           temp[29,53]=100.
        endif 
        writefits,path+name+'_velerr_temp.fits',temp
        velerrfile = path+name+'_velerr_temp.fits'
     endif
     hadetectfile = path+name+'_hadetect.fits'
     if file_test(hadetectfile) ne 1 then begin
        temp = readfits(velfile)
        temp(where(finite(temp))) = 10.
        writefits,path+name+'_hadetec_temp.fits',temp
        hadetectfile = path+name+'_hadetec_temp.fits'
     endif 
     
     ;create PSF gaussian
     npixkernel = 2.*max(fwhm[*,i])
     kernel = psf_gaussian_tilt(nx=npixkernel,ny=npixkernel,fwhm=fwhm[*,i],tilt=angle_psf(i))
     writefits,path+name+'_fwhm_kernel.fits',kernel
     kinematicmodel,velfile,velerrfile,hadetectfile,kernel,25.,theta0[i],cenxy[0,i],cenxy[1,i],cenxy[2,i],cenxy[3,i],name,pixscale=pixscale[i],rtmax=100.,vcmax=1000.

  endfor
end
