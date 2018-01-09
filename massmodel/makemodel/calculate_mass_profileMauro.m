function calculate_mass_profileMauro;
load transformedPIEMD1
% [da_x_dy,da_x_dx]=gradient(alpha_x_ALL);
% [da_y_dy,da_y_dx]=gradient(alpha_y_ALL);
% poisson_ALL=da_x_dx+da_y_dy;
% magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
%  
 d_rad=1
total_mass(1:3000/1)=0;
counta(1:3000/1)=0;
sum_mass=0;
poisson_ALL=poisson_ALL*0.5;
 for m=d_rad:2997-d_rad
    for n=d_rad:2997-d_rad
        rad1=sqrt((m+i_x-2470)^2+(n+i_y-2482)^2);
       
       bin=floor(rad1/d_rad+1);
        total_mass(bin)=total_mass(bin)+poisson_ALL(m,n);
        
        counta(bin)=counta(bin)+1;
         if (bin <=90)
            sum_mass=sum_mass+poisson_ALL(m,n)-0.9155;
        end
    end
 end
 sum_mass
 for i=1:3000/1
     total_mass(i);
     surface_density(i)=total_mass(i)/counta(i);
     radius(i)=i*0.065;
 end
 KappaVsR(:,1)=(radius);
 KappaVsR(:,2)=(surface_density);
 KappaVsR(:,3)=log10(radius);
 KappaVsR(:,4)=log10(surface_density);
 KappaVsR(:,5)=counta;
 
%plot((1:2500/20),(surface_density/2));
 plot(KappaVsR(:,1),KappaVsR(:,4),'black');
 sum_mass
  save 'KappaVsR_tryPIEMD1.txt' KappaVsR -ascii