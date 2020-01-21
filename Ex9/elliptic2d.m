%% FEM MATLAB sample code
[x,y,np,ne,e2,idp,ide] = readtria('disc');% read mesh from file
localtoglobal2DP1 = e2;                   % could it be this simple? why?
int = ~(idp==1);                          % select points without Dirichlet bc
nphi = 3;

%% build matrices
ii = zeros(ne,nphi^2); % sparse i-index
jj = zeros(ne,nphi^2); % sparse j-index
aa = zeros(ne,nphi^2); % entry of Galerkin matrix
bb = zeros(ne,nphi^2); % entry in mass-matrix (to build rhs)

%% build global from local
for k=1:ne             % loop over elements
    [Fdet,Finv] = generatetransformation2D(k,e2,x,y); % compute trafo

    % build local matrices (mass, stiffness, ...)
    sloc = localstiff2D(Fdet,Finv); % element stiffness matrix
    mloc = localmass2D(Fdet);       % element mass matrix

    % compute i,j indices of the global matrix
    dofs = localtoglobal2DP1(k,:);
    ii( k,: ) = [dofs([1 2 3 1 2 3 1 2 3])]; % local-to-global
    jj( k,: ) = [dofs([1 1 1 2 2 2 3 3 3])]; % local-to-global

    % compute a(i,j) values of the global matrix
    aa( k,: ) = sloc(:);
    bb( k,: ) = mloc(:);
end
% create sparse matrices
A=sparse(ii(:),jj(:),aa(:));
M=sparse(ii(:),jj(:),bb(:));

% build rhs and take into account Dirichlet bcs, solve, plot
rhs = M*ones(np,1);
u   = zeros(np,1);
u(int) = A(int,int) \ rhs(int);
% plotting
trisurf(e2,x,y,u)