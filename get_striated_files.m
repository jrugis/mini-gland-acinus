%
% get the striated file names
%
% J.Rugis
% 27.09.21
%
%

function [fnames] = get_striated_files()
  fnames = dir("*sCell*.ply");
end
