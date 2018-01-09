pro ana_radius

inc     = [33.6,32.,35.,36,43,33,33,43,31,33,33,38,32,43,47,40,32,40,35,36,34,43]
inc_err = [10.,10.,10,9,9,13,13,9,10,13,10,10,10,9,7,8,10,9,10,12,10,10]
pos     = [152,152.,148,157,158,205,158,149,164,151,211,153,167,158,210,153,210,162,152,148,158]
pos_err = [7.,10.,10,10,9,30,9,10,15,9,40,10,30,10,40,10,35,11,8,9,9]

!p.multi = [0,1,2]
modelnum = indgen(19)
ploterror,modelnum,inc,inc_err,xtitle='model number', ytitle='inclination',xrange=[-1,20]
ploterror,modelnum,pos,pos_err,xtitle='model number', ytitle='position angle',xrange=[-1,20]
stop
end
