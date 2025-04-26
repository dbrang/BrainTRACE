function write_pointset(fname, points)
%
% write a point set file readable by Freeview
%
fp = fopen(fname, 'wb');

for i=1:size(points, 1)
  fprintf(fp, '%f %f %f\n', points(i, :));
end

fprintf(fp, 'info\n');
fprintf(fp, 'numpoints %d\n', size(points, 1));
fprintf(fp, 'useRealRAS 1');

fclose(fp);
