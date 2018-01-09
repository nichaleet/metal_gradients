
pro metal_redshift
;make redshift vs metal gradient plot
;make total metallicity vs metal gradient plot
;make SFR vs metal gradient plot

;GLASS object
nameG14=['glass_arc4'];
zG14 = [1.855];
gradientG14=[-0.05];
gradient_errG14 = [0.05];
interactingG14=[1]

;Old objects
name_old=['M0744','S1038','S1148','S1206'];
zold=[2.21,2.20,2.38,2.00];
gradientOld=[-0.06,0.08,-0.28,-0.25];
gradient_errOld=[0.04,0.03,0.05,0.06];
interactingOld =[0,1,0,0]
N2HaOld = [0.45,0.04,0.11,0.22]
N2HaOlderr=[0.11,0.07,0.04,0.03]
delv_2sigmaOld = [1.42,0.98,0.82,0.76]
delv_2sigma_errOld = [0.45,0.27,0.30,0.32];
HafluxOld = [1.3,16.,7.,18.]*1.d-17 
Haflux_errOld = [0.3,2.,3.,1.]*1.d-17
HalumOld=fluxtolum(HafluxOld,zold)/1.d42
Halum_errOld = fluxtolum(Haflux_errOld,zold)/1.d42
SFRold      = [5.4,38,210,68] 
SFRoldupper = [15.3,75,377,132]
SFRoldlower = [3.6,23,43,44]
MassOld = [9.4,9.1,9.9,10.1]
MassOldUpper = [9.8,9.3,10.1,10.3]
MassOldLower = [9.3,9.0,9.6,9.9]

;New objects
name=['cswa11','cswa15','cswa19','cswa20','cswa28','cswa31','cswa128','cswa139','cswa159','cswa165','a733'];
z=[1.41,2.16,2.03,1.43,2.09,1.49,2.22,2.54,2.30,2.13,2.30] ;
;gradient    =[-0.07,-0.04,-0.01,-0.15,0.11,0.02,-0.10,-0.01,-0.01,0.03,-0.13]; pseudo slit
;gradient_err=[0.02 ,0.01 ,0.01 ,0.05 ,0.03,0.01,0.02 ,0.01 ,0.01,0.01,0.18]; pseudo slit
gradient     =[-0.11,-0.03,0.02,0.05,0.04,0.00,-0.04,0.04,-0.01,0.01,-0.05]
gradient_err =[0.02,0.01,0.01,0.02,0.03,0.01,0.01,0.01,0.01,0.02,0.02]
N2Ha =[0.14,0.03,0.09,0.05,0.09,0.33,0.12,0.17,0.08,0.26,0.17];
N2Ha_err = [0.01,0.01,0.01,0.01,0.01,0.06,0.03,0.08,0.01,0.06,0.07]
chisq = [17,18,9.6,35,6.6,93,989,38,7.6,141,17];
interacting = [0,0,1,1,0,1,1,1,0,1,0]
interacting(where(chisq gt 20.)) = 1
interacting(where(chisq lt 20.)) = 0
interacting(where(name eq 'cswa19'))=1
delv_2sigma = [0.89,0.71,0.57,0.37,0.26,0.91,1.98,1.34,0.65,0.94,0.60]
delv_2sigma_err = [0.19,0.06,0.05,0.09,0.12,0.08,0.44,0.36,0.17,0.26,0.07];

Haflux = [49.,31.,18.,3.,5.,14.,14.,1.6,9.,0.4,0.6]*1.d-17
Haflux_err = [5.,3.,2.,1.,1.,3.,3.,0.3,1.,0.1,0.1]*1.d-17
Halum      = fluxtolum(Haflux,z)/1.d42
Halum_err  = fluxtolum(Haflux_err,z)/1.d42
halflight_radius = [4.54165,3.48052,2.52721,1.67163,7.28319,7.33211,2.44646,4.24517,4.16522,0.618414,2.31872] ;kpc
SFRupper = [189.,47,70,9,23,69,321,49,100,9,37]
SFR      = [99.,42,59,6,12,36,250,33,53,7,30]
SFRlower = [11,37,48,4,1,4,108,9,6,2.8,23]
mass = [-99,-99,9.6,-99,9.4,-99,-99,-99,-99,-99,9.16]
massupper = [-99,-99,10.0,-99,9.8,-99,-99,-99,-99,-99,9.37]
masslower = [-99,-99,9.5,-99,9.1,-99,-99,-99,-99,-99,8.95]
;%Maciel et al 2003
zM03 = [1.4,0.65,0.116];
gradientM03 = [-0.1164,-0.0720,-0.0648];
gradient_errM03 =[0.0127,0.0057,0.0168];
interactingM03=[0,0,0]

;Yuan2011
zY11 = [1.5]
gradientY11 = [-0.16]
gradient_errY11 = [0.02]
interactingY11=[0]

;S12
zS12=[0.8425,1.455,1.4608,1.4625,1.4471,1.4858,2.2418]
gradientS12=[-0.037,-0.019,0.06,-0.027,-0.031,-0.087,-0.024]
gradientS12up = [0.03,0.019,0.017,0.010,0.016,0.032,0.012]
gradientS12low =[-0.058,-0.04,-0.004,-0.018,-0.014,-0.006,-0.012]
interactingS12=[0,0,0,0,1,0,1]
N2HaS12 = [999,0.43,0.1,0.27,0.13,0.6,0.6] ;the first gal redshift is too low to be compared
N2HaS12err=[0,0.05,0.1,0.03,.04,0.1,0.05]
HafluxS12 = [999,12.,11,10,16,13,16]*1.d-17 
Haflux_errS12 = [0.,1,1,1,2,1,1]*1.d-17
HalumS12=fluxtolum(HafluxS12,zS12)/1.d42
Halum_errS12 = fluxtolum(Haflux_errS12,zold)/1.d42
SFRs12 = [9999.,8,7,6,10,8,27]
SFRS12err = haflux_errs12/Hafluxs12*SFRs12
;simulation G13
ziso = [0,0.1,0.2,0.4,0.6,0.8,1.1,1.5,1.75,2.,2.3]
gradientiso = [-0.035,-0.05,-0.055,-0.062,-0.078,-0.09,-0.12,-0.24,-0.32,-0.22,-0.285]
zint = [0,0.1,0.2,0.4,0.6,0.8,1.1,1.5,1.75,2.]
gradientint = [-0.035,-0.03,-0.035,-0.032,-0.03,-0.02,-0.018,-0.025,-0.018,-0.02]


;PLOT
loadct,2
set_plot,'ps'
psname='z_gradient.eps'
device, filename = psname,xsize = 12,ysize = 9, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.thick=1.5
;M03
ploterror,zm03,gradientm03,gradient_errm03,psym=3,errcolor=14,yrange=[-0.35,0.12],xrange=[2.55,0.3],xmargin=[10,7],xtitle='Redshift',ytitle='Metallicity gradient (dex/kpc)',type=0,font=0,/nodata
int=where(interactingm03 eq 1,cint)
iso=where(interactingm03 eq 0,ciso)
if ciso ne 0 then cgplot,zm03(iso),gradientm03(iso),psym=14,color=fsc_color('navy'),/overplot,symsize=1.2
if cint ne 0 then cgplot,zm03(int),gradientm03(int),psym=4,color=fsc_color('navy'),/overplot,symsize=1.2
oplot,zm03,gradientm03,color=fsc_color('navy')

;S12

errplot,zs12,gradients12+gradients12up,gradients12+gradients12low,color=fsc_color('purple')
int=where(interactings12 eq 1)
iso=where(interactings12 eq 0)
cgplot,zs12(iso),gradients12(iso),psym=14,color=fsc_color('purple'),/overplot,symsize=1.2
cgplot,zs12(int),gradients12(int),psym=4,color=fsc_color('purple'),/overplot,symsize=1.2
;Yuan11
oploterror,zy11,gradienty11,gradient_erry11,psym=3,errcolor=fsc_color('dodger blue'),/nodata
cgplot,zy11,gradienty11,psym=14,color=fsc_color('dodger blue'),/overplot,symsize=1.2
;G14
oploterror,zg14,gradientg14,gradient_errg14,psym=3,errcolor=fsc_color('deep pink'),/nodata
cgplot,zg14,gradientg14,psym=6,color=fsc_color('deep pink'),/overplot,symsize=1.

;Old data
oploterror,zold,gradientold,gradient_errold,psym=3,/nodata,errcolor=fsc_color('lime green')
int=where(interactingold eq 1)
iso=where(interactingold eq 0)
cgplot,zold(iso),gradientold(iso),psym=15,color=fsc_color('lime green'),/overplot,symsize=1.
cgplot,zold(int),gradientold(int),psym=6,color=fsc_color('lime green'),/overplot,symsize=1.

;mydata
oploterror,z,gradient,gradient_err,psym=3,errcolor=fsc_color('dark red'),/nodata
int=where(interacting eq 1)
iso=where(interacting eq 0)
cgplot,z(iso),gradient(iso),psym=15,color=fsc_color('dark red'),/overplot,symsize=1.
cgplot,z(int),gradient(int),psym=6,color=fsc_color('dark red'),/overplot,symsize=1.

;simulation
oplot,ziso,gradientiso,color=fsc_color('turquoise'),thick=2.5
oplot,zint,gradientint,color=fsc_color('turquoise'),linestyle=2,thick=2.5

;legend
al_Legend,['L15','J13','J15','Y11','S12','M03'],psym=[15,15,6,14,14,14],color=fsc_color(['dark red','lime green','deep pink','dodger blue','purple','navy']),position=[0.15,-0.08],/right,box=0,font=0
al_Legend,['','','','',''],psym=[6,6,3,3,4],position=[0.,-0.08],color=[fsc_color('dark red'),fsc_color('lime green'),fsc_color('deep pink'),fsc_color('dodger blue'),200],/right,box=0,font=0,linestyle=[6,6,-1,-1,6]
al_Legend,[''],psym=[14],position=[-0.05,-0.265],color=[fsc_color('navy blue')],/right,box=0,font=0,linestyle=[0]
al_Legend,['','G13 normal feedback','G13 enhanced feedback'],psym=[0,0,0],color=fsc_color(['navy blue','turquoise','turquoise']),linestyle=[0,0,2],linsize=0.4,position=[-0.06,-0.264],thick=3,/right,box=0,font=0

device,/close
;----------------------------------------------------------------------
;Plot SFR and stellar mass
psname='SFR_mass.eps'
device, filename = psname,xsize = 12,ysize = 9, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
goodmass = where(mass ne -99)
massplot = [massold,mass(goodmass)]
massplotupper = [massoldupper,massupper(goodmass)]
massplotlower = [massoldlower,masslower(goodmass)]
SFRplot = [SFRold,SFR(goodmass)]
SFRplotupper = [SFRoldupper,SFRupper(goodmass)]
SFRplotlower = [SFRoldlower,SFRlower(goodmass)]
intplot = [0,1,0,0,1,0,0]
cgplot,massplot,SFRplot,psym=3,xmargin=[10,3],ymargin=[5,3],ytitle='SFR(M'+textoidl('_{sun}')+'/year)',xtitle='Log(M) [M'+textoidl('_{sun}')+']',font=0,/nodata,xrange=[8.5,10.5],xstyle=1,charsize=1,/ylog
cgErrPlot,SFRplot[0:3],massplotupper[0:3],massplotlower[0:3],/horizontal,color='lime green',thick=1.5
cgErrPlot,massplot[0:3],SFRplotupper[0:3],SFRplotlower[0:3],color='lime green',thick=1.5
cgErrPlot,SFRplot[4:6],massplotupper[4:6],massplotlower[4:6],/horizontal,color='dark red',thick=1.5
cgErrPlot,massplot[4:6],SFRplotupper[4:6],SFRplotlower[4:6],color='dark red',thick=1.5
cgPlot,massplot([5,6]),SFRplot([5,6]),color='dark red',psym=15,/overplot
cgPlot,massplot[4],SFRplot[4],color='dark red',psym=6,/overplot
cgplot,massplot([0,2,3]),SFRplot([0,2,3]),color='lime green',psym=15,/overplot
cgPlot,massplot[1],SFRplot[1],color='lime green',psym=6,/overplot
;oplot with main sequence from Behroozi2013b
msequence = [9,9.5,10,10.5]
ceq = [3.135,2.169,1.873,1.129]*1.e-9
aeq = [-1.028,-0.993,-1.219,-1.426]
beq = [-0.06,-0.08,-0.023,-0.083]
z0eq = [1.,1.,1.,1.]
SSFRsequence = ceq/[10.^(aeq*(2.-z0eq))+10.^(beq*(2.-z0eq))]
SFRsequence = SSFRsequence*10.^msequence
SFRsequenceupper = SFRsequence+10.^[0.22,0.22,0.24,0.35]
SFRsequencelower = SFRsequence-10.^[0.22,0.22,0.24,0.35]
cgErrPlot,msequence,SFRsequenceupper,SFRsequencelower,color='blue',thick=1.5
cgPlot,msequence,SFRsequence,psym=14,color='blue',/overplot
al_Legend,['L15','J13','Behroozi 2013b'],psym=[15,15,14],color=fsc_color(['dark red','lime green', 'blue']),position=[10.4,5.],/right,box=0,font=0,thick=2
al_Legend,['',''],psym=[6,6],position=[10.5,5],color=fsc_color(['dark red','lime green']),/right,box=0,font=0,thick=2

device,/close
SFRall = [SFRold,SFR]
SFRallErr = [(SFRoldupper-SFRoldlower)/2.,(SFRupper-SFRlower)/2.]
meanerr,SFRall,SFRallerr,SFRwmean,SFRmeanErr,SFRdev,SFRwdev
SFRmedian = median(SFRall)
SFRmean = mean(SFRall)
print,'SFR'
print,'mean wmean meanerr median stdev weithed stdev'
print,SFRmean,SFRwmean,SFRmeanErr,SFRmedian,SFRdev,SFRwdev
stop
;-----------------------------------------------------------------------
;Plot total metallicity vs gradient with S12
alphaLetter = '!9' + String("141B) + '!X'  
;betaletter  = '!9' + String("142B) + '!X'
loadct,2
psname='totalmetal_gradient_withS12.eps'
device, filename = psname,xsize = 12,ysize = 9, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,N2Ha,gradient,N2Ha_err,gradient_err,psym=3,errcolor=fsc_color('dark red'),yrange=[-0.35,0.12],xrange=[0.02,0.65],xmargin=[10,7],xtitle='total [NII]/H'+alphaLetter,ytitle='Metallicity gradient (dex/kpc)',type=0,font=0,/nodata,errthick=1.5,/xlog,xstyle=5
int=where(interacting eq 1)
iso=where(interacting eq 0)
cgplot,N2Ha(iso),gradient(iso),psym=15,/overplot,color=fsc_color('dark red')
cgplot,N2Ha(int),gradient(int),psym=6,/overplot,color=fsc_color('dark red'),thick=2

oploterror,N2Haold,gradientold,N2Haolderr,gradient_errold,psym=3,errcolor=fsc_color('lime green'),errthick=1.5
int=where(interactingold eq 1)
iso=where(interactingold eq 0)
cgplot,N2Haold(iso),gradientold(iso),psym=15,color=fsc_color('lime green'),/overplot
cgplot,N2Haold(int),gradientold(int),psym=6,color=fsc_color('lime green'),/overplot,thick=2

oploterror,N2HaS12,gradientS12,N2HaS12err,gradientS12*0.,psym=3,errcolor=fsc_color('purple'),errthick=1.5
errplot,N2HaS12,gradients12+gradients12up,gradients12+gradients12low,color=fsc_color('purple'),thick=1.5
int=where(interactingS12 eq 1)
iso=where(interactingS12 eq 0)
cgplot,N2HaS12(iso),gradientS12(iso),psym=15,color=fsc_color('purple'),/overplot
cgplot,N2HaS12(int),gradientS12(int),psym=6,color=fsc_color('purple'),/overplot,thick=2
plotsym,6,2,/fill
oplot,[N2HaS12(2)],[gradientS12(2)],psym=8,thick=2,color=fsc_color('purple')

al_Legend,['L15','J13','S12'],psym=[15,15,15],color=fsc_color(['dark red','lime green', 'purple']),position=[0.52,-0.27],/right,box=0,font=0,thick=2
al_Legend,['','',''],psym=[6,6,6],position=[0.6,-0.27],color=fsc_color(['dark red','lime green', 'purple']),/right,box=0,font=0,thick=2

xtitle='total [NII]/H'+alphaLetter
;add x axis
for i=0,1 do begin 
   axis,XAXIS=i,/XLOG,/XSTYLE,XRANGE=10.^!X.CRANGE, $
        XTICKFORMAT=i?'(A1)':'',XTICKS=3, $
        XTICKV=[.02,0.1,.5],XMINOR=0,XTITLE=i?'':xtitle,font=0
endfor 
device,/close
;Find correlation coefficient
x_all =[N2Ha,N2Haold,N2HaS12]
y_all = [gradient,gradientold,gradientS12]
y_allerr = [gradient_err,gradient_errold,0.5*(gradientS12up-gradientS12low)]
x_all = alog10(x_all)
int_all = [interacting,interactingold,interactingS12]
result = regress(x_all,y_all,sigma=sigma,const=const,measure_errors=y_allerr,correlation=correlation)
print,'All galaxies', correlation,correlate(x_all,y_all)

iso = where(int_all eq 0)
int = where(int_all eq 1)
x = x_all(iso)
y = y_all(iso)
y_err = y_allerr(iso)
result = regress(x,y,sigma=sigma,const=const,measure_errors=y_err,correlation=correlation)
print,'isolated galaxies', correlation,correlate(x,y)

x = x_all(int)
y = y_all(int)
y_err = y_allerr(int)
result = regress(x,y,sigma=sigma,const=const,measure_errors=y_err,correlation=correlation)
print,'interacting galaxies', correlation,correlate(x,y)


;Plot SFR vs gradient with S12
;alphaLetter = '!9' + String("141B) + '!X'  
loadct,33
psname='SFR_gradient_withS12.eps'
device, filename = psname,xsize = 12,ysize = 9, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
cgplot,SFR,gradient,psym=3,yrange=[-0.35,0.12],xmargin=[10,3],ymargin=[5,3],xtitle='Star Formation Rate (M'+textoidl('_{sun}')+'/year)',ytitle='Metallicity gradient (dex/kpc)',font=0,/nodata,xrange=[1,400],/xlog,xstyle=1,charsize=1
;plot Tucker
cgErrPlot,gradientold,SFRoldlower,SFRoldupper,/horizontal,color='lime green',thick=1.5
cgErrPlot,SFRold,gradientold-gradient_errold,gradientold+gradient_errold,color='lime green',thick=1.5
cgplot,sfrold(iso),gradientold(iso),psym=15,color='lime green',/overplot
cgplot,sfrold(int),gradientold(int),psym=6,color='lime green',/overplot,thick=2
;Plot S12
cgErrPlot,SFRs12,gradients12low,gradients12up,color='purple',thick=1.5
cgErrPlot,gradients12,SFRs12-sfrs12err,sFRs12+sfrs12err,color='purple',thick=1.5,/horizontal
int=where(interactingS12 eq 1)
iso=where(interactingS12 eq 0)
cgplot,SFRS12(iso),gradientS12(iso),psym=15,color='purple',/overplot
cgplot,SFRS12(int),gradientS12(int),psym=6,color='purple',/overplot,thick=2

;Plot my data
cgErrPlot,gradient,SFRlower,SFRupper,/horizontal,color=fsc_color('dark red'),thick=1.5
cgErrPlot,SFR,gradient-gradient_err,gradient+gradient_err,color=fsc_color('dark red'),thick=1.5
int=where(interacting eq 1)
iso=where(interacting eq 0)
cgplot,SFR(iso),gradient(iso),psym=15,/overplot,color=fsc_color('dark red')
cgplot,SFR(int),gradient(int),psym=6,/overplot,color=fsc_color('dark red'),thick=2
int=where(interactingold eq 1)
iso=where(interactingold eq 0)

al_Legend,['L15','J13','S12'],psym=[15,15,15],color=fsc_color(['dark red','lime green','purple']),position=[4,-0.25],/right,box=0,font=0,thick=2
al_Legend,['','',''],psym=[6,6,6],position=[5.2,-0.25],color=fsc_color(['dark red','lime green','purple']),/right,box=0,font=0,thick=2
device,/close


;Plot Delv_sigma vs Gradeint
DeltaLetter = '!9'+String("104B)+'!X'
sigmaLetter = '!9'+String("163B)+'!X'
psname='vel_gradient_annuli.eps'
device, filename = psname,xsize = 14,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,delv_2sigma,gradient,delv_2sigma_err,gradient_err,psym=0,xrange=[0,3],yrange=[-.5,0.4],xmargin=[10,7],xtitle=DeltaLetter+'V/2'+SigmaLetter,ytitle='Metallicity gradient (dex/kpc)',type=0,font=0,/nodata,errthick=1.5,color=fsc_color('black'),errcolor=fsc_color('dark red')
int=where(interacting eq 1)
iso=where(interacting eq 0)
cgplot,delv_2sigma(iso),gradient(iso),psym=15,/overplot,color=fsc_color('dark red'),size=2
cgplot,delv_2sigma(int),gradient(int),psym=6,/overplot,color=fsc_color('dark red'),thick=2,size=2

oploterror,delv_2sigmaold,gradientold,delv_2sigma_errold,gradient_errold,errcolor=fsc_color('lime green'),psym=3
int=where(interactingold eq 1)
iso=where(interactingold eq 0)
cgplot,delv_2sigmaold(iso),gradientold(iso),psym=15,/overplot,color=fsc_color('lime green'),size=2
cgplot,delv_2sigmaold(int),gradientold(int),psym=6,/overplot,color=fsc_color('lime green'),thick=2,size=2

al_Legend,['L15','J13'],psym=[15,15],color=fsc_color(['dark red','lime green']),position=[2.9,-0.4],/right,box=0,font=0,thick=2
al_Legend,['',''],psym=[6,6],position=[3.,-0.4],color=fsc_color(['dark red','lime green']),/right,box=0,font=0,thick=2
device,/close

;gradient    =[-0.07,-0.04,-0.01,-0.15,0.11,0.02,-0.10,-0.01,-0.01,0.03,-0.13]; pseudo slit
;gradient_err=[0.02 ,0.01 ,0.01 ,0.05 ,0.03,0.01,0.02 ,0.01 ,0.01,0.01,0.18]; pseudo slit
;psname='vel_gradient_pseudoslit.eps'
;device, filename = psname,xsize = 14,ysize = 8, $
;        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
;ploterror,delv_2sigma,gradient,delv_2sigma_err,gradient_err,psym=0,xrange=[0,3],yrange=[-.5,0.4],xmargin=[10,7],xtitle=DeltaLetter+'V/2'+SigmaLetter,ytitle='Metallicity gradient (dex/kpc)',type=0,font=0,/nodata,errthick=1.5,errcolor=fsc_color('dark red')
;int=where(interacting eq 1)
;iso=where(interacting eq 0)
;cgplot,delv_2sigma(iso),gradient(iso),psym=14,/overplot,color=fsc_color('dark red'),size=2
;cgplot,delv_2sigma(int),gradient(int),psym=4,/overplot,color=fsc_color('dark red'),thick=2,size=2

;oploterror,delv_2sigmaold,gradientold,delv_2sigma_errold,gradient_errold,errcolor=fsc_color('lime green'),psym=3
;int=where(interactingold eq 1)
;iso=where(interactingold eq 0)
;cgplot,delv_2sigmaold(iso),gradientold(iso),psym=14,/overplot,color=fsc_color('lime green'),size=2
;cgplot,delv_2sigmaold(int),gradientold(int),psym=4,/overplot,color=fsc_color('lime green'),thick=2,size=2
;device,/close

;Find correlation coefficient
x_all =[delv_2sigma,delv_2sigmaold]
y_all = [gradient,gradientold]
y_allerr = [gradient_err,gradient_errold]
int_all = [interacting,interactingold]
result = regress(x_all,y_all,sigma=sigma,const=const,measure_errors=y_allerr,correlation=correlation)
print,'All galaxies', correlation,correlate(x_all,y_all)

iso = where(int_all eq 0)
int = where(int_all eq 1)
x = x_all(iso)
y = y_all(iso)
y_err = y_allerr(iso)
result = regress(x,y,sigma=sigma,const=const,measure_errors=y_err,correlation=correlation)
print,'isolated galaxies', correlation,correlate(x,y)
x = x_all(int)
y = y_all(int)
y_err = y_allerr(int)
result = regress(x,y,sigma=sigma,const=const,measure_errors=y_err,correlation=correlation)
print,'interacting galaxies', correlation,correlate(x,y)
stop

;Plot half light radius vs gradient
alphaLetter = '!9' + String("141B) + '!X'  
betaletter  = '!9' + String("142B) + '!X'
loadct,2
psname='halflight_gradient.eps'
device, filename = psname,xsize = 12,ysize = 9, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

ploterror,halflight_Radius,gradient,halflight_radius*0.1,gradient_err,psym=3,errcolor=fsc_color('dark red'),yrange=[-0.35,0.12],xmargin=[10,7],xtitle='Half light radius (kpc)',ytitle='Metallicity gradient (dex/kpc)',type=0,font=0,/nodata,errthick=1.5,xrange=[0.1,10.],xstyle=1
int=where(interacting eq 1)
iso=where(interacting eq 0)
cgplot,halflight_radius(iso),gradient(iso),psym=15,/overplot,color=fsc_color('dark red')
cgplot,halflight_radius(int),gradient(int),psym=6,/overplot,color=fsc_color('dark red'),thick=2

al_Legend,['L15','J13'],psym=[15,15],color=[75,10],position=[9.5,-0.3],/right,box=0,font=0,thick=2
al_Legend,['',''],psym=[6,6],position=[10,-0.3],color=[75,10],/right,box=0,font=0,thick=2
device,/close
stop
end


