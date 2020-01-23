%% FEM MATLAB sample code
N    =  32;                         % 1D: number of points
nphi =   2;                         % P1: 2 local basis functions
x = linspace(0,1,N);                % array of points
[ne,np,e1] = generateelements1D(x); % generate mesh
localtoglobal1DP1 = e1;             % could it be this simple? why?

%% build matrices
ii = zeros(ne,nphi^2); % sparse i-index
jj = zeros(ne,nphi^2); % sparse j-index
aa = zeros(ne,nphi^2); % entry of Galerkin matrix
bb = zeros(ne,nphi^2); % entry in mass-matrix (to build rhs)

%% build global from local
for k=1:ne             % loop over elements
    [Fdet,Finv] = generatetransformation1D(k,e1,x); % compute trafo
    
    % build local matrices (mass, stiffness, ...)
    sloc = localstiff1D(Fdet,Finv); % element stiffness matrix
    mloc = localmass1D(Fdet);       % element mass matrix
    
    % compute i,j indices of the global matrix
    dofs = localtoglobal1DP1(k,:);
    ii( k,: ) = [dofs(1) dofs(2) dofs(1) dofs(2)]; % local-to-global
    jj( k,: ) = [dofs(1) dofs(1) dofs(2) dofs(2)]; % local-to-global
    
    % compute a(i,j) values of the global matrix
    aa( k,: ) = sloc(:);
    bb( k,: ) = mloc(:);
end
% create sparse matrices
A=sparse(ii(:),jj(:),aa(:));
M=sparse(ii(:),jj(:),bb(:));

% build rhs and take into account Dirichlet bcs, solve, plot
rhs = M*ones(N,1);
u   = zeros(N,1);
u(2:end-1) = A(2:end-1,2:end-1) \ rhs(2:end-1);
plot(x,u,x,x.*(1-x)/2,'r.')