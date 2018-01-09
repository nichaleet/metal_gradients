pro plotnoise
restore,'slope.sav'
noise = findgen(11)*0.2
slope_res = fltarr(4,11)
slope_res_err = fltarr(4,11)
slope_real = slope[*,0]
slope_real_err = sigma[*,0]
for i=0,10 do begin
   slope_res[*,i]=slope[*,i]-slope_real
   slope_res_err[*,i] = sqrt(sigma[*,i]^2+slope_real_err^2)
endfor


loadct,2
set_plot,'ps'
psname='cswa128slope_addnoise.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

ploterror,noise,slope_res[0,*],slope_res_err[0,*],ytitle='metallicity gradient residual(dex/kpc)',xtitle='Noise Added (x [NII] detection limit)',psym=3,type=0,font=0,errcolor=200,/nodata,xrange=[-0.1,2.1],xstyle=1
cgplot,noise,slope_res[0,*],psym=14,color=200,/overplot

oploterror,noise,slope_res[2,*],slope_res_err[2,*],ytitle='metallicity gradient residual(dex/kpc)',xtitle='Noise Added (x[NII] detection limit)',psym=3,type=0,font=0,errcolor=20,/nodata
cgplot,noise,slope_res[2,*],psym=14,color=20,/overplot
loadct,0
oplot,!x.crange,[0,0],color=100
device,/close

stop
end
