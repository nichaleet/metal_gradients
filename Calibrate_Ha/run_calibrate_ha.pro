pro run_calibrate_ha
c128 = {starfile:'/scr2/nichal/workspace/reduced_data/Sep2013/ttstars/s130911_a020001__mosaic_Kn2_100.fits',name:'cswa128',filter:'Kn2',mag:14.493,z:2.225}
c11mar = {starfile:'/scr2/nichal/workspace/reduced_data/Mar2013/skysubtracted/cswa11tt/s130304_a008001__mosaic_Hn2_100.fits',name:'cswa11mar',filter:'Hn2',mag:14.201,z:1.41}
c11feb = {starfile:'/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/CSWA11tt/s140223_a005001__mosaic_Hn2_100.fits',name:'cswa11feb',filter:'Hn2',mag:14.201,z:1.41}
c15={starfile:'/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/cswa15tt/s140223_a018001__mosaic_Kn2_100.fits',name:'cswa15',filter:'Kn2',mag:14.193,z:2.16}
c19={starfile:'/scr2/nichal/workspace/reduced_data/Mar2013/skysubtracted/cswa19tt/s130303_a011001__mosaic_Kn1_100.fits',name:'cswa19',filter:'Kn1',mag:12.502,z:2.0327}
c20={starfile:'/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/cswa20tt/s140224_a020001__mosaic_Hn2_100.fits',name:'cswa20',filter:'Hn2',mag:15.495,z:1.43274}
c28={starfile:'/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/cswa28tt/s140223_a022001__mosaic_Kn1_100.fits',name:'cswa28',filter:'Kn1',mag:12.397,z:2.09234}
c31={starfile:'/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/cswa31tt/s141206_a012001__mosaic_Hn3_100.fits',name:'cswa31',filter:'Hn3',mag:13.071,z:1.48737}
c139 = {starfile:'/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/cswa139tt/s141207_a015001__mosaic_Kc5_100.fits',name:'cswa139',filter:'Kc5',mag:14.892,z:2.54213}
c159 = {starfile:'/scr2/nichal/workspace/reduced_data/Dec2014/skysubtract/cswa159tt/s141207_a005001__mosaic_Kc3_100.fits',name:'cswa159',filter:'Kc3',mag:12.293,z:2.29863}
c165 = {starfile:'/scr2/nichal/workspace/reduced_data/Sep2013/skysubtracted/cswa165tt/s130911_a028001__mosaic_Kn2_100.fits',name:'cswa165',filter:'Kn2',mag:13.471,z:2.1278}
a773 = {starfile:'/scr2/nichal/workspace/reduced_data/Feb2014/skysubtracted/a773tt/s140222_a005001__mosaic_Kc3_100.fits',name:'abell773',filter:'Kc3',mag:10.668,z:2.3027}
list = [c31,c20,a773]
list = [c128]
flux_calibration = replicate({name:'name',calibfactor:0D,calibfactor_err:0D},n_elements(list))
for i=0,n_elements(list)-1 do begin
   calibrate_ha,list[i].starfile,name=list[i].name,filter=list[i].filter,photmag=list[i].mag,z=list[i].z,calibfactor=calibfactor
   flux_calibration[i].name = list[i].name
   flux_calibration[i].calibfactor = calibfactor[0]
   flux_calibration[i].calibfactor_err = calibfactor[1]
endfor
save,flux_calibration,filename='flux_calibration_c128.sav'
stop
end
