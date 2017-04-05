function [M] = readbin(filename)
% [M] = readbin2(filename)
%
% read various binary format matrix files
% M - single matrix or array of matrices
%
% use this instead of readbin

fid = fopen(filename,'r');

magic = 'MATRIX ARRAY ';
[a,c] = fread(fid,length(magic),'char');
if c == length(magic) & char(a) == magic',
  % cell array of matrices
  a = fgetl(fid);
  [nmat,count,errmsg,nextindex] = sscanf(a,'%d');
  if count ~= 1,
    error(['error parsing ' filename]);
  end
  M = cell(nmat,1);
  for i = 1:nmat
    [sz,c] = fread(fid,2,'int32');
    [M{i},c] = fread(fid,sz','double');
  end;
  fclose(fid);
else
  % single matrix
  fclose(fid);
  fid = fopen(filename,'r');
  [sz,c] = fread(fid,2,'int32');
  [M,c] = fread(fid,sz','double');
  fclose(fid);
end

