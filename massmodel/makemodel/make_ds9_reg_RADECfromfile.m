function make_ds9_reg_RADECfromfile(filename)
f=load(filename);
%circlestr='circle('
fid = fopen('model_colors_cswa19v1toedit.reg', 'w');

% print a title, followed by a blank line
%fprintf(fid, '# ReadMe File\n\n # Constructed and supplied by Adi Zitrin');

for i=1:length(f(:,1))
fprintf(fid, 'circle(%f, %f, 20) # text={}\n', f(i,2), f(i,3));

end
fclose(fid);
%save 'MembersRXJ2129.reg' new_file -ascii