function [ras2xyz, c] = read_vg(fid)
%
% read and parse volume geometry information at
% the end of a MRI surface file.
% the format of vol geom is inferred from freesurfer/utils/transform.c
%

% ignore the beginning part of volume geometry 
fgets(fid);
fgets(fid);
fgets(fid);
fgets(fid);

% local to world transformation
ras2xyz = zeros(3, 3);
for i=1:3
  ras2xyz(i, :) = sscanf(fgets(fid), '%*s   = %e %e %e')';
end

% translation vector
c = sscanf(fgets(fid), '%*s   = %e %e %e')';