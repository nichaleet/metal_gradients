pro fig3

;group1 = ['cswa11','cswa15','cswa28','cswa159','Abell773'] ;rot supported
;group2 =['cswa19','cswa20','cswa31','cswa128','cswa139','cswa165']

name = ['cswa11','cswa15','cswa19both','cswa20','cswa28','cswa31','cswa128','cswa139','cswa159','cswa165','Abell773']

group1 = fltarr(11)
vc_sig = [1.84,0.85,3.55,2.15,1.62,2.17,8.63,4.27,2.31,5.89,3.91]
vc_sig_low =[1.12,0.48,2.81,1.12,0.80,1.80,4.83,1.90,1.07,3.08,2.17]
vc_sig_high =[2.52,1.20,5.76,3.03,2.19,2.56,12.43,6.51,3.48,8.49,5.30]

delv_2sig = [0.89,0.71,1.54,0.37,0.26,0.91,1.98,1.34,0.65,0.94,0.60]
delv_2sig_err = [0.19,0.06,0.30,0.09,0.12,0.08,0.44,0.36,0.17,0.26,0.28]

chisq = float([17,18,246,35,6.6,93,990,38,7.6,141,17])

group1(where(chisq lt 40.)) = 1. ;rotation is flagged as 1. 

DeltaLetter = '!9'+String("104B)+'!X'
sigmaLetter = '!9'+String("163B)+'!X'
ChiLetter = '!9'+String("166B)+'!X'
ChiLetter = '!9'+String("143B)+'!X'

set_plot,'ps'
loadct,11
psname='vc_delv_plot.eps'
device, filename = psname,xsize = 10,ysize = 12, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
multiplot,/init
multiplot,[1,2]

ploterror,chisq,delv_2sig,delv_2sig_err,psym=3,/xlog,ytitle=DeltaLetter+'V/2'+SigmaLetter,type=0,font=0,/nodata,xrange=[4,1200],yrange=[0.1,3],/ylog,xstyle=1,ystyle=1
cgplot,chisq,delv_2sig,color=200,psym=14,symsize=1.5,/overplot
loadct,0
cgplot,[1,10000],[0.4,0.4],linestyle=2,color=100,/overplot
cgplot,[20,20],[0.1,10],linestyle=2,color=100,/overplot

multiplot
loadct,11
plot,chisq,vc_sig,/xlog,psym=3,ytitle='V'+'!IC!N'+'/'+SigmaLetter,font=0,/nodata,xtitle=Chiletter+textoidl('_{red}^2'),xrange=[4,1200],yrange=[0.1,14],/ylog,xstyle=1,ystyle=1
errplot,chisq,vc_sig_low,vc_sig_high
cgplot,chisq,vc_sig,color=100,psym=14,symsize=1.5,/overplot
loadct,0
cgplot,[1,10000],[1.,1.],linestyle=2,color=100,/overplot
cgplot,[20,20],[0.1,10],linestyle=2,color=100,/overplot

device,/close
stop
end

;ploterror,delv_2sig(rot),vc_sig(rot),delv_2sig_err(rot),vc_sig(rot)*0.,psym=3,xmargin=[10,7],xtitle=DeltaLetter+'V/2'+SigmaLetter,ytitle='V'+'!IC!N'+'/'+SigmaLetter,type=0,font=0,/nodata,xrange=[0,2.5],yrange=[0,12]
;errplot,delv_2sig(rot),vc_sig_low(rot),vc_sig_high(rot)
;cgplot,delv_2sig(rot),vc_sig(rot),color=200,psym=14,symsize=1.5,/overplot

;oploterror,delv_2sig(chao),vc_sig(chao),delv_2sig_err(chao),vc_sig(chao)*0.,psym=3
;errplot,delv_2sig(chao),vc_sig_low(chao),vc_sig_high(chao)
;cgplot,delv_2sig(chao),vc_sig(chao),color=200,psym=4,symsize=1.5,/overplot,thick=2
;loadct,0
;oplot,!x.crange,[1,1],linestyle=5,color=100
;oplot,[0.4,0.4],!y.crange,linestyle=5,color=100
;;oplot,!x.crange,2.5*!x.crange,color=100
