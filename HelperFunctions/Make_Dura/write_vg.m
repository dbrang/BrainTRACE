function write_vg(fid, volsize, voxelsize, ras2xyz, c)
%
% write volume geometry information to fid.
% the format of volume geometry information is inferred 
% from freesurfer/utils/transform.c
%

fwrite(fid, 20, 'int32');
fprintf(fid, 'valid = 1  # volume info valid\n');
fprintf(fid, 'filename = \n');
fprintf(fid, 'volume = %d %d %d\n', volsize);
fprintf(fid, 'voxelsize = %.15e %.15e %.15e\n', voxelsize);
fprintf(fid, 'xras   = %.15e %.15e %.15e\n', ras2xyz(1, :));
fprintf(fid, 'yras   = %.15e %.15e %.15e\n', ras2xyz(2, :));
fprintf(fid, 'zras   = %.15e %.15e %.15e\n', ras2xyz(3, :));
fprintf(fid, 'cras   = %.15e %.15e %.15e\n', c);
