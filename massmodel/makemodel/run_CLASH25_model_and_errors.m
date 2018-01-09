function run_CLASH25_model_and_errors
%This function is the "main" function
%user has to define (many) inputs:
user_name='Adi Zitrin'
GPU_OR_CPU='CPU';
cluster_name='cswa15'; %make sure teher are no spaces in name: eg 'Abell 383' prohibited, 'Abell383' or 'A383' are good. 
%(read the readme file for file input definitions 
%for the multiple images and member galaxies formats) %%comment: I did not
%yet write this file

galaxies_file=load('model_colors_templim24updatedAdi.txt'); % member galaxy photometry file: col 2,3 should be x and y in the original image pixels, and col 7 should be the flux
images_file=load('imagev1.txt'); %multiple image file: col1,2: x and y in original image pixels, col3: dls/ds rationm relative to a dls/ds==1 for the "main_source_z" 
%value which user will input below, col 4: no. of images in each system, col 5: no. of systems in file, col 6: parity if applicable (either 1, -1, or 0 if not relevant), cols 7,8 are not used but should be for future 
%purposes: col7: arc name digits only (e.g. 1.1), col8: redshift estimate of that system (the one that corresponds to the dls/ds relative ratio specified in col 3). 
pix_scale=0.15; %pixel scale of original image
shapes_file=load('blank.txt');%WL/ ellipticity shapes file, col1,2 should be the RA and DEC offsets from the BCG (or the first galaxy in the galaxies file), col 5,6: e1 and e2 (NEEDS CHECKING WITH ADI)
%Also, if the user does not want to use WL and only wants to use SL, they
%should load here an empty file with their name of choice (e.g.
%load('WLnullfile.txt')


Wresolution=1; %factor of lower resolution wrt to pixel scale (used only for 
%first stage of minimization, result will be in original pixel scale),
%usually for HST a value of 4,5,10, or 20 works well here (the lower value
%the slower the first part chains are),

%if spline interpolation is chosen below, resolution will be always /4
GalaxyRepresentation='PowerLaw'; %representations of each galaxy; 
%current legal entries are: 'PowerLaw' or 'PIEMD' ; option 'PIEMD' is
%currently disabled for the LTM part so only use PowerLaw
SplineOrGaussian='Gaussian'; %representation of DM component; 
%current legit options are: 
% 'Spline' 2D spline interpolation (with an external shear addition)
% 'Gaussian' (with an external shear addition)
% 'eGaussian' elliptical Gaussian kernel (no extr. shear) % comment: not
% available yet for the LTM
% 'eNFW' elliptical NFW % comment: not
% available yet for the LTM, use ../PIEMDeNFW library instead if you wish
% to use that one
% 'PE_NFW' pseudo elliptical NFW % comment: not
% available yet 
No_of_Bright_Free=0; %no. of brightest galaxies to optimize weight, user can specify galaxies' index below, or else code will take the XX brightest ones
No_of_redshift_Free=0;%8; %no. of systems with only photometric redshifts to be optimized by the model, indexes have to be specified below
No_of_Parity_Force=1; %no. of systems to force critical curves pass exactly between its specified images

cluster_redshift=0.306; % cluster redshift
main_source_z=2.162; %main_system redshift, this will correspond to dls/ds==1 in images file.
M_star=-21.08; %relevant only for galaxy PIEMD representation, not for PowerLaw LTM
save_file_name='MCMC_TryTRYTRY'; % ignore this, do not touch

%BCG ellipticity parameters:
BCGel=0; % (b-a)/(a+b) form e.g. SExtractor  %these are ignored if choosing to leave free below...
BCGPA=0; % as in SExtractor (ie in degrees!)
%BCG ellipse parameter priors:
BCGell=0; BCGelh=0.5; %BCGeli=0.2;
BCGPAl=-60; BCGPAh=60; %BCGPAi=0;
%Should I leave these values free to be optimized by the MC? (1==yes, 0==no) :
BCGelFreeflag=0;
BCGPAFreeflag=0; %important note: if the BCG is *not* left free, make sure the priors above are wide enough to contain the BCG fixed el and PA, 
%even though the priors are not used in practice - since they will anyhow
%be randomized 

%As above, second brightest free ellipticity?:
BCGel2=0; % (b-a)/(a+b)  %these are ignore if choosing to leave free...
BCGPA2=0; % as in SExtractor (ie in degrees!)
%2nd BCG ellipticity parameters, if applicable:
BCGell2=0; BCGelh2=0.5; %BCGeli=0.2;
BCGPAl2=-60; BCGPAh2=60; %BCGPAi=0;
BCGelFreeflag2=0;
BCGPAFreeflag2=0; 
%this function saves the input parameters later used for naming saved
%files:
try_save_ignore(cluster_name,cluster_redshift,main_source_z,SplineOrGaussian)

%%%% now define priors
ql=1.25; qh=1.35; qi=1;%range of power law values and initial value (the latter is ignored), values between 1 and 1.5 are OK/enough.
if strcmp(SplineOrGaussian,'Spline')==1
sl=4; sh=24; si=24; %range of smoothing values and initial value (the latter is ignored), if Spline, , values between 2 and 24 are OK/enough. (typical best-fit value is 10 or 12)
elseif strcmp(SplineOrGaussian,'Gaussian')==1
sl=8; sh=30; si=500; %range of smoothing values and initial value (the latter is ignored), if Gaussian, , values between few dozens and few hundrads are OK/enough. (typical best-fit value is ~250)
end
k_newl=0.1; k_newh=0.2; k_newi=0.42; %range of overall scaling values and initial% value (the latter is ignored), values between 0.5 and 4 are OK/enough
%
k_gall=0.1; k_galh=0.2; k_gali=0.03;  %range of relative galaxies weight-wrt-DM values and initial% value (the latter is ignored), values between few percents (eg 0.03) and few dozen percents (eg 0.6) are OK/enough
gammal=0; gammah=0.1; gammai=0.14; %range of external shear strength values and initial% value (the latter is ignored) (values between 0 and 0.3 are usually enough/OK)
phil=0; phih=pi; phii=pi/2;% %range of external shear PA values, should be between 0 and pi unless user prefers to guess the direction a priori
rcl=0; rch=80; rci=0; % core radius of the BCG, in pixels, between 0 and 200 is more than enough
rcl2=0; rch2=80; rci=0; % core radius 2nd brightest
CoreFree=0; % mark 1 if you want it to be a free parameter, otherwise it will be set to zero
CoreFree2=0; % mark 1 if you want it to be a free parameter, otherwise it will be set to zero

BCGl=0.3; BCGh=5; BCGi=4; % range of possibel weight for the BCG if left freely weighted (eg between 0.5 to 5 times the original photometry value)
galsl=0.3; galsh=5; galsi=1; % range of possibel weight for the other galaxies left freely weighted, (eg between 0.5 to typically 2-3 times the original photometry value)
%galsl2=0.8; galsh2=1.2; galsi2=1; %ignore this line
redsl=-0.3; redsh=0.1; redsi=0; % %range of values to add to the specified in column 3 of images files, i.e. range of dls/ds values that can be added to the relative dls/ds of each system that was left to be freely scaled. This will be later automatically translated to redshift
%comment: above reds should be in dls/ds if working regularly, or in
%redshift if for cosmo!
CosmoFree=0; %leave cosmology free? (always choose no, ==0, unless you are working on constraining the cosmology with lensing...)
% priors ranges for cosmology
if CosmoFree==1
wai=0; %constant wa
w_0l=-2; w_0h=0;  w_0i=-1;
omegaMl=0; omegaMh=1; omegaMi=0.3;
omegaLl=0; omegaLh=1; omegaLi=0.3;
omegaLi=1-omegaMi; % flat universe
end
H_zeroi=70; % constant H_zero
KTmc=180; %boltzman; for first chains
%other definitions:
%Nmcmc=max_steps;
mean_sample_redshift=1.14;  %%ATTENZIONE parametro!
relat_weight_SL=1; % weight SL in chi2 minimization   %%ATTENZIO%N%E parametro!
relat_weight_WL=0; % weight WL in chi2 minimization   %%ATTENZIO%N%E parametro!
fiducial_chi2_mult=1; % just because we work in reduced chi2
jumpi=60; %salto  %%%ATTENZIONE parametro!!!
DOFnormalize=0;% flag whether to normalize according to the DOF of the SL and WL regimes separately, usually choose no, in which case the SL and WL chi^2's will be simply added

%decide if this is for publication and you would like to make images of the
%results, at least in this stage (can be run independently later):
prepare_images_flag=1; %in general choose 0, 1 if you have supplied 
% %define input for prepare images scripts
% if prepare_images_flag==1
%    image_to_plot_on='a383.png';
%    multiple_images_list_toplot='';
%    multiple_images_list_totable=''; % same as above but in ra/dec format
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_many_err_models=0;%this parameters is ==1 saves all models generated for the error estimation; not recommended to use more than 100 error models if ==1
save('save_many_err_models.mat','save_many_err_models');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now caulculating effective number of parameters - User does not have to
%touch this at all:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CosmoFree==1
Nparams=17+No_of_Bright_Free+No_of_redshift_Free;
num_param=Nparams-9;%since Ho wa and lamba known, and asusming bcgs fixed ;
if BCGelFreeflag==1
    num_param=num_param+1;
end
if BCGPAFreeflag==1
    num_param=num_param+1;
end
if BCGelFreeflag2==1
    num_param=num_param+1;
end
if BCGPAFreeflag2==1
    num_param=num_param+1;
end
if CoreFree==1
    num_param=num_param+1;
end
if CoreFree2==1
    num_param=num_param+1;
end
else
Nparams=12+No_of_Bright_Free+No_of_redshift_Free;
num_param=Nparams-6;
if BCGelFreeflag==1
    num_param=num_param+1;
end
if BCGPAFreeflag==1
    num_param=num_param+1;
end
if BCGelFreeflag2==1
    num_param=num_param+1;
end
if BCGPAFreeflag2==1
    num_param=num_param+1;
end
if CoreFree==1
    num_param=num_param+1;
end
if CoreFree2==1
    num_param=num_param+1;
end
end

Nparams 
% global number of params
num_param 
%effective number of params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_data1=galaxies_file;
if No_of_Bright_Free>=1
for i=1:No_of_Bright_Free
maxtemp=max(file_data1(:,7));%7th column is the flux
maxtemp
Ind_Free_Gal(i)=find(file_data1(:,7)==maxtemp);
Ind_Free_Gal(i)
file_data1(Ind_Free_Gal(i),:)=0;
end
%or define manually, with the BCG first, eg:
elseif No_of_Bright_Free==0
  Ind_Free_Gal(1)=0;
  
end

%Ind_Free_Gal(1:3)=[1 2 3];
%define index for redshift free systems, index should be system number!
if No_of_redshift_Free>=1
   Ind_Free_Reds(1:4)=[1 3 4 6];
%eg:
%Ind_Free_Reds(1:4)=[4 5 6 7];
% Ind_Free_Reds(1:9)=[3 4 5 6 7 8 12 13 14];
elseif No_of_redshift_Free==0
   Ind_Free_Reds(1)=0;%
end
if No_of_Parity_Force>=1  % if used
    Ind_Parity_Force(1:2)=[1 3]; %recall to add parity sign to images file column 6!
else
    Ind_Parity_Force(1)=[0];
end

i_x=300 %FOV start point on x axis wrt original image 
i_y=650 %FOV start point on y axis wrt original image 
leng=550 %length of FOV on each axis
save MCinput
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%we suggest that before running the long chain, check ByEye to see that it
%fits...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define first short chains control parameters
m_chains=round(Nparams/5);  %number of chains
n_steps_per_chain=round(330*Nparams/10);
run_multiple_MCs(m_chains,n_steps_per_chain)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % still to do - check that (a) converged thus far and (b) that acceptance
% % % rate is reasonable, otherwise adjust parameters
% % %-generalize resolution
% % %-parity choice
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %run one long final chain, until convergence is achieved,
clear all
max_steps=300;
KTmc2=125;
WresolutionOne=1;
run_final_MC_chain(max_steps,KTmc2,WresolutionOne)
save step2
%then, refine around the best fit and re-check for onvergence
%(or maybe simply take errors from before refinement but after convergence)
max_steps2=300;
KTmc3=5;
desired_acceptedsteps2=2000;
%WresolutionOne=2;
run_final_MC_chainFiner(max_steps2,KTmc3,desired_acceptedsteps2,WresolutionOne)

%then, run examination by eye, measure r_e and mass etc., and save the best
%model, including alpha, magnification, and time delays also
WresolutionTwo=1;
[rms_best,chi2_WL_best,ndf_WL_best,chi2red_WL_best,chi2red_SL_best,ndf_SL_best,num_images,num_systems,num_param,num_pixels_shear]=MCstep_PowerLaw_and_GaussianOrSplineByEyeSave(WresolutionTwo)
 save('rms_best.mat','rms_best','chi2_WL_best','ndf_WL_best','chi2red_WL_best','chi2red_SL_best','ndf_SL_best','num_images','num_systems','num_param','num_pixels_shear');
 ref_fits_file='cswa15_R_small.fits';
save_bestModel_to_fits(ref_fits_file)
%calculate error matrices (and input file names)
how_many_toerr=100; % the error on error estimation goes like 1/sqrt(how_manys_toerr); ~500-1000 should be enough
WresolutionThree=1;
calculate_Map_errorsRand(how_many_toerr,ref_fits_file,WresolutionThree) %(but needs revision! also make draw random...?)
calculate_Map_errorsRandDefs(how_many_toerr,ref_fits_file,WresolutionThree) %(but needs revision! also make draw random...?)

%still to do: change to 68/2 instead of std, include profile separately?
prepare_images_flag=1;
if prepare_images_flag==1
   preparefiles_and_errors_general(WresolutionThree,how_many_toerr,ref_fits_file) 
end
%calculate errors on parameters (and input file names)
extract_1sigma_errs
% 
% %check redshift prediciton for systems and for candidates and errs 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % still to do - make 2d plots of resulting params
% %also would be nice to include/ expand in future
% %- WL!
% %- include gaussian and cosmological params
% %- choose red sequence
% %- full gpu make (if faster), and cpu version
% %-reminimize/refine with larger FOV
% %-BCg ellip
% %-PIEMDeNFW through multiple chains
% %-solve ra/dec shift when interpolating...maybe by finalizing with an orig
% %resolution chain (also more accurtae solution like this...)?
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % calculate evidence
% %calculate_evidence
% %input all to a readme file (to do)
% %include in package also ds9 region files...
% 
% %get_critarea_andcausticarea % but needs revision!
%choose redshift: taking for now redshift 2
[refirst,final_re]=get_critarea_andcausticarea_z(2)
% raedme:
write_readme_file
Post_processing2_addParamsHeader
mkdir 2PostOnline
copyfile('*.fits', '2PostOnline/');
copyfile('Params_*.txt', '2PostOnline/');
copyfile('ReadMe.txt', '2PostOnline/');
zip('2PostOnline','2PostOnline');
