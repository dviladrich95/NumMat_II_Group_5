% read a mesh from triangle
% sample code by Dirk Peschka

% Note: Code probably needs to be changed if extra marker are introduced

function [x,y,npoint,nelement,e2p,idp,ide] = readtria(fname)

fnode = strcat(fname,'.node');
fele  = strcat(fname,'.ele');

fprintf('read mesh file: %s \n',fname);

% read .node file
fid = fopen(fnode,'r');
mesh_header = fscanf(fid,'%i',4);

npoint = mesh_header(1); % number of points
ndim   = mesh_header(2); % dimension
npattr = mesh_header(3); % number of point attributes
nbound = mesh_header(4); % number of boundary markers

x   = zeros(npoint,1);
y   = zeros(npoint,1);
idp = zeros(npoint,1);

for i=1:npoint
    dump  =fscanf(fid,'%i',1); % index running from 1..npoint
    x(i)  =fscanf(fid,'%f',1); % x-position
    y(i)  =fscanf(fid,'%f',1); % y-position
    idp(i)=fscanf(fid,'%i',1); % boundary marker
end
fclose(fid);

% read .ele file
fid = fopen(fele,'r');
mesh_header = fscanf(fid,'%i',3);

nelement = mesh_header(1); % number of points
nphi     = mesh_header(2); % dimension
neattr   = mesh_header(3); % number of element attributes

e2p = zeros(nelement,nphi);
ide = zeros(nelement,neattr);

for k=1:nelement
    dump    =fscanf(fid,'%i',1); % index running from 1..nelement
    for j=1:nphi
        e2p(k,j)=fscanf(fid,'%i',1); % element-to-point 1..nphi
    end
    for j=1:neattr
        ide(k,j)=fscanf(fid,'%i',1); % element-marker 1..neattr
    end
end
fclose(fid);

npoint = max(max(e2p(:,1:3))); % use npoint for P1 points