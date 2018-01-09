function calculate_numberdensity_profile;
i_x=1500;
i_y=1500;
ff=load('model_colors_temp.txt');
   count=0;
  tot(1:3000/20)=0; 
  dist(1:3000/20)=0; 
  area_fiducial(1:3000/20)=0;
  dens(1:3000/20)=0;
for j=20:20:3000-20
    count=count+1;
    dist(count)=j * 0.05; %in arcsec
    if j==20
   area_fiducial(count)=pi*dist(count)^2;
    elseif j>20
        area_fiducial(count)=pi*dist(count)^2-area_fiducial(count-1);    
        end
  for i=1:length(ff(:,1))
      if sqrt((ff(i,2)-3.1792240e+03)^2+(ff(i,3)-2.8254380e+03)^2)<j+20 && sqrt((ff(i,2)-3.1792240e+03)^2+(ff(i,3)-2.8254380e+03)^2)<j ;
      tot(count)=tot(count)+1;
     
      end
  end
end
dens=tot./area_fiducial;
 KappaVsR(:,1)=(dist);
 KappaVsR(:,2)=(dens);
 KappaVsR(:,3)=(tot);
 KappaVsR(:,4)=(area_fiducial);
 KappaVsR(:,5)=log10(dens);

 plot(KappaVsR(:,1),KappaVsR(:,2),'black');
 
  save 'densplot.txt' KappaVsR -ascii