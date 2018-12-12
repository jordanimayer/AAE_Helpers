%%%%%
% Jordan Mayer
% AAE 532
%
% save_all_views:
%   Save four different views of a 3-D plot as files: 3D, XY, XZ, and
%   YZ.
%
% Inputs:
%   fig: figure to save
%   fileName: suffix for file names; note all files will be .jpg with
%             prefices _3D, _XY, _XZ, _YZ respectively
%%%%%

function [] = save_all_views(fig, fileName)

    fileName_3d = [fileName '_3D.jpg'];
    fileName_xy = [fileName '_XY.jpg'];
    fileName_xz = [fileName '_XZ.jpg'];
    fileName_yz = [fileName '_YZ.jpg'];
    
    saveas(fig, fileName_3d);
    view(0, 90);
    saveas(fig, fileName_xy);
    view(0, 0);
    saveas(fig, fileName_xz);
    view(90, 0);
    saveas(fig, fileName_yz);
end
