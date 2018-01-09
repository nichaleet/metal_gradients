function add_galaxies_by_click(cat_name)

display('loading everything this may take up to 5 minutes...')
% temp_im=imread('stiff.tif');
% 
% for m=1:3000
%  for n=1:3000
%    new_im(m,n,:)=temp_im((3000+2)-n,m,:);
%      
%  end
% end
i_x=1250;
i_y=1250;
plot_circles_like_regions('stiff2500.tif',cat_name);
g=load(cat_name);
length_g=length(g(:,1));
beg=length_g;
f=load('temp_cat_allIR.txt');
[x,y]=ginput;
x=x+i_x
y=(2500-y)+i_y
lines=length(f(:,1));
length(x)
for bla=1:length(x)
for i=1:lines
    if (x(bla)<f(i,2)+5 && x(bla)>f(i,2)-5 && (y(bla))<f(i,3)+5 && (y(bla))>f(i,3)-5)
        length_g=length_g+1;
        g(length_g,:)=f(i,:); %*10^-2
    end
end
end
length_g
display('850 fluxes are:')
g(beg+1:length_g,7)
save 'model_colors_new.txt' g -ascii
plot_circles_like_regions('stiff2500.tif','model_colors_new.txt');