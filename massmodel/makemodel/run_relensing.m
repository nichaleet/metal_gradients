function run_relensing
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
mkdir('outputmaps')
% maketifmaps('CSWA15_Ha_interp.fits')
% maketifmaps('CSWA15_Hb_interp.fits')
% maketifmaps('CSWA15_kinematic_interp.fits')
% maketifmaps('CSWA15_NII_interp.fits')
% maketifmaps('CSWA15_OIII_interp.fits')
% maketifmaps('CSWA15_N2index_interp.fits')
% maketifmaps('CSWA15_veldisp_interp.fits')
% maketifmaps('CSWA15_metaltype_interp.fits')
% maketifmaps('CSWA15_Bayesianmetal_interp.fits')
% maketifmaps('cswa15_O3N2index_interp.fits')

transform_imagefromMCspeedup('cswa15_stiff.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourceHa.tif')


transform_imagefromMCspeedup('CSWA15_Ha_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourceha.tif')

transform_imagefromMCspeedup('CSWA15_Hb_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourcehb.tif')

transform_imagefromMCspeedup('cswa15_kinematic_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourcekinematic.tif')

transform_imagefromMCspeedup('cswa15_N2index_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourcen2index.tif')

transform_imagefromMCspeedup('cswa15_NII_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourcenii.tif')

transform_imagefromMCspeedup('cswa15_veldisp_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourceveldisp.tif')

transform_imagefromMCspeedup('cswa15_metaltype_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourcemetaltype.tif')

transform_imagefromMCspeedup('cswa15_Bayesianmetal_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourceBayesianmetal.tif')

transform_imagefromMCspeedup('cswa15_OIII_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourceoiii.tif')

transform_imagefromMCspeedup('cswa15_O3N2index_interp.tif')
lensing_tool_real_color_hr_zoom(1,567,880,60)
movefile('source.tif','outputmaps/sourceo3n2index.tif')

end

