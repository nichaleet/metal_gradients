function calculate_mass_profile;
load transformedTry1BCG1RedApr16Gaussian
[da_x_dy,da_x_dx]=gradient(alpha_x_ALL);
[da_y_dy,da_y_dx]=gradient(alpha_y_ALL);
poisson_ALL=da_x_dx+da_y_dy;

 d_rad=1
total_mass(1:leng/1)=0;
counta(1:leng/1)=0;
sum_mass=0;
poisson_ALL=poisson_ALL*0.5;
 for m=d_rad:leng-d_rad-2
    for n=d_rad:leng-d_rad-2
        rad1=sqrt((m+i_x-2501.343)^2+(n+i_y-2500.412)^2);
        
       bin=floor(rad1/d_rad+1);
        total_mass(bin)=total_mass(bin)+poisson_ALL(m,n);
        counta(bin)=counta(bin)+1;
%          if (bin <=90)
%             sum_mass=sum_mass+poisson_ALL(m,n)-0.9155;
%         end
    end
 end
 sum_mass
 for i=1:leng/1
     total_mass(i);
     surface_density(i)=total_mass(i)/counta(i);
     radius(i)=i*pix_scale;
 end
 KappaVsR(:,1)=(radius);
 KappaVsR(:,2)=(surface_density);
 KappaVsR(:,3)=log10(radius);
 KappaVsR(:,4)=log10(surface_density);
 KappaVsR(:,5)=counta;
 
%plot((1:2500/20),(surface_density/2));
 plot(KappaVsR(:,1),KappaVsR(:,4),'blue');
 sum_mass
  save 'KappaVsR_try.txt' KappaVsR -ascii