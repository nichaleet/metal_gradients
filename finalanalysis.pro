pro finalanalysis

;make graphs between metal gradients and other parameters

name = ['cswa15','cswa20','cswa28','cswa128','cswa165','Abell773','macs744','j1148','j1206']
newdata=[0,1,2,3,4,5]
olddata=[6,7,8]
name_merge  =['cswa19','j1038']

gradient     = [-0.057,0.278,0.021,0.028,-0.045,-0.615,-0.06,-0.28,-0.25]
gradient_err = [0.053,0.322,0.019,0.006,-0.045,0.265,0.04,0.05,0.06]
gradient_merge     = [-0.078,0.08]
gradient_merge_err = [0.059,0.03]

;NII_HA_ratio=[]

z = [2.16,1.43,2.09,2.23,2.13,2.30,2.21,2.38,2.00]
z_merge = [2.03,2.20]

vrot =[50.67,72.20,38.97,251.5,45.90,52.58,252.,148.,159.]
vrot_err =[3.45,26.95,8.92,15.72,11.52,5.33,33.,3.,38.]
vdisp =[29.1,59.91,59.94,50.35,82.64,47.21,89.,90.,104]
vdisp_err=[14.26,14.55,14.43,29.77,42.24,9.79,26.,33.,37.]

vrot_merge=[55.84,160.]
vrot_merge_err =[9.45,10.]

vdisp_merge=[73.46,82.]
vdisp_merge_err = [13.62,22.]

rotation=0.5*vrot/vdisp
rotation_err = rotation*sqrt(vrot_err^2/vrot^2+vdisp_err^2/vdisp^2)
rotation_merge = 0.5*vrot_merge/vdisp_merge
rotation_err_merge = rotation_merge*sqrt(vrot_merge_err^2/vrot_merge^2+vdisp_merge_err^2/vdisp_merge^2)


;redshift and gradient
gradient_z_old = errorplot(z(olddata),gradient(olddata),gradient_err(olddata),xtitle='redshift',ytitle='Metallicity gradient(PP04) dex/kpc',symbol='diamond',sym_color='purple',sym_filled=1,errorbar_color='purple',linestyle=6,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9],xrange=[2.5,1.3],yrange=[0.5,-1.],name='isolate galaxy(Jones13)')

gradient_z_new = errorplot(z(newdata),gradient(newdata),gradient_err(newdata),/current,overplot=1,symbol='circle',sym_color='purple',sym_filled=1,errorbar_color='purple',linestyle=6,thick=2,name='isolate galaxy(new)')

gradient_z_merge_old = errorplot([z_merge(1)],[gradient_merge(1)],[gradient_merge_err(1)],/current,overplot=1,symbol='diamond',sym_color='green',sym_filled=1,errorbar_color='green',linestyle=6,thick=2,name='merger galaxy(Jones13)')

gradient_z_merge_new = errorplot([z_merge(0)],[gradient_merge(0)],[gradient_merge_err(0)],/current,overplot=1,symbol='circle',sym_color='green',sym_filled=1,errorbar_color='green',linestyle=6,thick=2,name='merger galaxy(new)')

leg = LEGEND(TARGET=[gradient_z_old,gradient_z_new,gradient_z_merge_old,gradient_z_merge_new],position=[0.5,0.85],/normal)


gradient_z_merge_new.save,'/scr2/nichal/workspace/output/finalanalysis/gradient_redshift.png'
stop
gradient_z_merge_new.close



;rotation and gradient

rotation_gradient_old = errorplot(rotation(olddata),gradient(olddata),rotation_err(olddata),gradient_err(olddata),xtitle='$\Delta V/2\sigma$',ytitle='metallicity gradient(dex/kpc)',symbol='diamond',sym_color='purple',sym_filled=1,errorbar_color='purple',linestyle=6,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9],name='isolate galaxy(Jones13)',xrange=[0,3])

rotation_gradient_new = errorplot(rotation(newdata),gradient(newdata),rotation_err(newdata),gradient_err(newdata),/current,overplot=1,symbol='circle',sym_color='purple',sym_filled=1,errorbar_color='purple',linestyle=6,thick=2,name='isolate galaxy(new)')

rotation_gradient_merge_old = errorplot([rotation_merge(1)],[gradient_merge(1)],[rotation_err_merge(1)],[gradient_merge_err(1)],/current,overplot=1,symbol='diamond',sym_color='green',sym_filled=1,errorbar_color='green',linestyle=6,thick=2,name='merger galaxy(Jones13)')

rotation_gradient_merge_new = errorplot([rotation_merge(0)],[gradient_merge(0)],[rotation_err_merge(0)],[gradient_merge_err(0)],/current,overplot=1,symbol='circle',sym_color='green',sym_filled=1,errorbar_color='green',linestyle=6,thick=2,name='merger galaxy(new)')

leg = LEGEND(TARGET=[rotation_gradient_old,rotation_gradient_new,rotation_gradient_merge_old,rotation_gradient_merge_new],position=[.5,0.85],/normal)

stop

rotation_gradient_merge_new.save,'/scr2/nichal/workspace/output/finalanalysis/rotation_gradient.png'

rotation_gradient_merge_new.close


;rotation velocity and gradient

vrot_gradient_old = errorplot(vrot(olddata),gradient(olddata),vrot_err(olddata),gradient_err(olddata),xtitle='$\Delta V$',ytitle='metallicity gradient(dex/kpc)',symbol='diamond',sym_color='purple',sym_filled=1,errorbar_color='purple',linestyle=6,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9],name='isolate galaxy(Jones13)')

vrot_gradient_new = errorplot(vrot(newdata),gradient(newdata),vrot_err(newdata),gradient_err(newdata),/current,overplot=1,symbol='circle',sym_color='purple',sym_filled=1,errorbar_color='purple',linestyle=6,thick=2,name='isolate galaxy(new)')

vrot_gradient_merge_old = errorplot([vrot_merge(1)],[gradient_merge(1)],[vrot_merge_err(1)],[gradient_merge_err(1)],/current,overplot=1,symbol='diamond',sym_color='green',sym_filled=1,errorbar_color='green',linestyle=6,thick=2,name='merger galaxy(Jones13)')

vrot_gradient_merge_new = errorplot([vrot_merge(0)],[gradient_merge(0)],[vrot_merge_err(0)],[gradient_merge_err(0)],/current,overplot=1,symbol='circle',sym_color='green',sym_filled=1,errorbar_color='green',linestyle=6,thick=2,name='merger galaxy(new)')

leg = LEGEND(TARGET=[vrot_gradient_old,vrot_gradient_new,vrot_gradient_merge_old,vrot_gradient_merge_new],position=[.5,0.85],/normal)

stop

vrot_gradient_merge_new.save,'/scr2/nichal/workspace/output/finalanalysis/vrot_gradient.png'

vrot_gradient_merge_new.close

end
