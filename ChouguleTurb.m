function [PHI,k11,k2_log,k3_log] = ChouguleTurb(alphaEps,GAMMA,L,zeta,varargin)
% PHI,k11,k2_log,k3_log] = ChouguleTurb(alphaEps,GAMMA,L,zeta,varargin)
% compute the two sided uniform shear model with stability correction by [1]  spectral tensor for wind turbulence.
% 
% INPUT
% 
% alphaEps = 1st constant of Mann spectral tensor [1 x 1]
% GAMMA = shear constant (2nd constant) [1 x 1]
% L = Integral length scale [1 x 1]
% zeta = Non-dimensional Obukhov length [1 x 1]
% varargin: It can be;
%  - N1: Number of points in the along-wind direction [1 x 1]
%  - N2: Number of points in the across-wind direction [1 x 1]
%  - N3: Number of points in the vertical-wind direction [1 x 1]
%  - k1min: min value of wavenumber for k1 [1 x 1]
%  - k2min: min value of wavenumber for k2 [1 x 1]
%  - k3min: min value of wavenumber for k3 [1 x 1]
%  - k1max = max value of wavenumber for k1 [1 x 1]
%  - k2max: max value of wavenumber for k2 [1 x 1]
%  - k3max: max value of wavenumber for k3 [1 x 1]
%  - Ninterp: Number of interpolation points for 2F1 approximation [ 1 x 1]
%  - Niter: Number of iteration for solving the Rapid-distortion theory 
%    equation with 4th order Runge-Kutta method [ 1 x 1]
% 
% OUTPUT
% 
% PHI : 5D Spectral tensor for the three wind components [2N1 x 2N2 x 2N3 x 4 x 4]
% k11: : Single-sided wavenumber [1,1 x N1]
% k2_log : two sided wavenumber  [1, 2 x N2]
% k3_log : two sided wavenumber  [1, 2 x N3]
% 
%  SYNTHAX:
% 
% [PHI,k2,k3,k11,k2_log,k3_log] = ChouguleTurb(alphaEps,GAMMA,L,zeta)
% [PHI,k2,k3,k11,k2_log,k3_log] = ChouguleTurb(alphaEps,GAMMA,L,zeta'N1',50)
% [PHI,k2,k3,k11,k2_log,k3_log] = ChouguleTurb(alphaEps,GAMMA,L,zeta,'Ninterp',100)
% 
% References:
% Chougule, A., Mann, J., Kelly, M., & Larsen, G. C. (2018).
% Simplification and validation of a spectral-tensor model for turbulence 
% including atmospheric stability. 
% Boundary-Layer Meteorology, 167(3), 371-397.
% 
% author: E Cheynet - UiB - last modified : 10-03-2020
% 
% see also MannCoherence.m fitMannTensor.m
%  
%% INPUT parser
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Niter',10); % 
p.addOptional('N1',30); % 
p.addOptional('N2',40); % 
p.addOptional('N3',40); % 
p.addOptional('k1min',log10(1/300)); % 
p.addOptional('k2min',-4); % 
p.addOptional('k3min',-4); % 
p.addOptional('k1max',log10(10)); % 
p.addOptional('k2max',log10(10)); % max value of wavenumber for k2
p.addOptional('k3max',log10(10)); % 
p.addOptional('Ninterp',200); % 
p.parse(varargin{:});
% check number of input: Number of outputs must be >=3 and <=6.
nargoutchk(0,4)
% shorthen the variables name
N1 = p.Results.N1 ;
N2 = p.Results.N2 ;
N3 = p.Results.N3 ;
k1min = p.Results.k1min ;
k2min = p.Results.k2min ;
k3min = p.Results.k3min ;
k1max = p.Results.k1max ;
k2max = p.Results.k2max ;
k3max = p.Results.k3max ;
Ninterp = p.Results.Ninterp ;
Niter = p.Results.Niter ;
%% construction of the wavenumber vectors
k1_log=[-fliplr(logspace(k1min,k1max,N1)),logspace(k1min,k1max,N1)];
k2_log=[-fliplr(logspace(k2min,k2max,N2)),logspace(k2min,k2max,N2)];
k3_log=[-fliplr(logspace(k3min,k3max,N3)),logspace(k3min,k3max,N3)];
%% Get beta and k0
[k2,k1,k3]=meshgrid(k2_log,k1_log,k3_log);
k = sqrt(k1.^2+k2.^2+k3.^2); % k from previous step
ValToInterp = -(k(:).*L).^(-2);
x = sort(ValToInterp);
x = x(1:round(numel(ValToInterp)/Ninterp):end);
F = griddedInterpolant(x,hypergeom([1/3,17/6],4/3,x));
HYPERGEOM = F(ValToInterp);
beta = GAMMA.*(k(:).*L).^(-2/3).*(HYPERGEOM).^(-1/2);
beta = reshape(beta,2*N1,2*N2,2*N3);
k30 = k3 + beta.*k1;
k0 = sqrt(k1.^2+k2.^2+k30.^2); % k from previous step
%% Initial conditions and preallocation
NDOF = 4; % the tree velocity component  + temperature
dZ0 =  zeros(NDOF,NDOF,2*N1,2*N2,2*N3);
PHI = nan(2*N1,2*N2,2*N3,NDOF,NDOF);
for i1=1:2*N1
    for i2=1:2*N2
        for i3=1:2*N3
            dZ0(:,:,i1,i2,i3) = getDZ0(k0(i1,i2,i3),k1(i1,i2,i3),k2(i1,i2,i3),k30(i1,i2,i3),zeta,alphaEps,L); % isotropic covariance
        end
    end
end
dZ = dZ0;
%% Main loop
Fun = @(Y,A) A*Y; % RDT equation
h = beta/(Niter); % time step (non-dimensional)
for idt=1:Niter
    newK3 = k30 - beta.*k1.*idt/Niter;
    newK = sqrt(k1.^2+k2.^2+newK3.^2); % k from previous step
    [M] = getM(k1,k2,newK3,newK,zeta);
    for i1=1:2*N1
        for i2=1:2*N2
            for i3=1:2*N3
                [dZ(:,:,i1,i2,i3)] = RK4(Fun,dZ(:,:,i1,i2,i3),h(i1,i2,i3),M(:,:,i1,i2,i3));
            end
        end
    end
end

for i1=1:2*N1
    for i2=1:2*N2
        for i3=1:2*N3
            PHI(i1,i2,i3,:,:) = dZ(:,:,i1,i2,i3)*dZ(:,:,i1,i2,i3)';
        end
    end
end
k11 = k1_log(end-N1+1:end);
%% Nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Y] = RK4(Fun, Y, dt,A)
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A);
        k_2 = Fun(Y+0.5*dt*k_1,A);
        k_3 = Fun(Y+0.5*dt*k_2,A);
        k_4 = Fun(Y+k_3*dt,A);
        
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end
    function [DZ] = getDZ0(k0,k1,k2,k30,zeta,alphaEps,L)
        
        Beta1 = 0.8;
        alpha = 1.7; %  the (spectral) Kolmogorov constant
        
        if zeta<=0
            etaT = zeta./(1/zeta.*(1+16.*abs(zeta)).^(-1/4)-1);
        elseif zeta>0 && zeta <=1
            A = (zeta./(1+5.*zeta)).^2;
            B = 1-zeta./(1+5*zeta);
            etaT = A./B;
        end
        
        
        E = alphaEps.*L.^(5/3).*(L.*k0).^4./((1+(L.*k0).^2).^(17/6));
        
        Beta = Beta1./alpha;
        S_prime = Beta.*etaT.*(1+(k0.*L).^2)./((k0.*L).^2).*E;
        Z4 = sqrt(k0.^2.*S_prime./E);
        DZ1 = [0, k30, -k2, 0];
        DZ2 = [-k30, 0, k1, 0];
        DZ3 = [k2, -k1, 0, 0];
        DZ4 = [0, 0, 0, Z4];
        DZ = [DZ1;DZ2;DZ3;DZ4].*sqrt(E./(4*pi.*k0.^4));
    end
    function [M] = getM(k1,k2,k3,k,zeta)
        
        M = zeros(NDOF,NDOF,2*N1,2*N2,2*N3);
        
        if zeta >=0
            Ri = zeta./(1+5.*zeta);
        else
            Ri = zeta;
        end
        
        M(1,1,:,:,:) =  0    ;
        M(1,2,:,:,:) =  0    ;
        M(1,3,:,:,:) =  2*k1.^2./k.^2-1    ;
        M(1,4,:,:,:) =  -k1.*k3./k.^2    ;
        
        M(2,1,:,:,:) =  0    ;
        M(2,2,:,:,:) =  0    ;
        M(2,3,:,:,:) =  2*k1.*k2./k.^2    ;
        M(2,4,:,:,:) =  -k2.*k3./k.^2    ;
        
        M(3,1,:,:,:) =  0    ;
        M(3,2,:,:,:) =  0    ;
        M(3,3,:,:,:) =  2*k1.*k3./k.^2    ;
        M(3,4,:,:,:) =  -(k3.^2./k.^2-1)    ;
        
        M(4,1,:,:,:) =  0    ;
        M(4,2,:,:,:) =  0    ;
        M(4,3,:,:,:) =  -Ri    ;
        M(4,4,:,:,:) =  0    ;
    end
end
