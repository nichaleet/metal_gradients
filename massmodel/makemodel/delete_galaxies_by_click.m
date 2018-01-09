function delete_galaxies_by_click(cat_name)

display('loading everything this may take up to 5 minutes...')

i_x=800;
i_y=800;
plot_circles_like_regions('stiff.tif',cat_name);
g=load(cat_name);

[x,y]=ginput;
x=x+i_x
y=(3000-y)+i_y
lines=length(g(:,1));
length(x)
for bla=1:length(x)
for i=1:lines
    if (x(bla)<g(i,2)+10 && x(bla)>g(i,2)-10 && (y(bla))<g(i,3)+10 && (y(bla))>g(i,3)-10)
   display('OK')
        g(i,:)=0;
    end
end
end

save 'model_colors_new.txt' g -ascii
%clean_temp_cat_new_from_zeros
length(g(:,1))
plot_circles_like_regions('stiff.tif','model_colors_new.txt');