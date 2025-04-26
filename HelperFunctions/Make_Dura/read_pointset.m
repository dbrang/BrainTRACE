function points = read_pointset(fname)
%
% read a point set file produced by Freeview
%
fp = fopen(fname, 'r');

point = sscanf(fgets(fp), '%f', 3)';
points = [];

while size(point, 2) == 3
  points = [points; point];
  point = sscanf(fgets(fp), '%f', 3)';
end

fclose(fp);