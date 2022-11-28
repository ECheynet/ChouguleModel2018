function [GAMMA,L,alphaEps,zeta] = fitChougule(k11,Su,Sv,Sw,Suw,zeta0,varargin)
% [Gamma,L,alphaEps] =
% fitChougule(k11,Su,Sv,Sw,Suw,alphaEps,guess) fits the
% stability-corrected Uniform shear model [1] to measured 2-sided single 
% and cross-spectra. The fitting procedure can be done using 3 or 4-unknown parameters. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k11: single sided wavenumber in the along wind direction.
% Su, Sv, Sw, Suw : [1 x N1] single-point, single-sided non-normalized 
% wind spectra.
% varargin: it can be;
%  - N2: Number of points in the across-wind direction [1 x 1]
%  - N3: Number of points in the vertical-wind direction [1 x 1]
%  - k2min: min value of wavenumber for k2 [1 x 1]
%  - k3min: min value of wavenumber for k3 [1 x 1]
%  - k2max: max value of wavenumber for k2 [1 x 1]
%  - k3max: max value of wavenumber for k3 [1 x 1]
%  - tolX: tolerance for fitting procedure
%  - tolFun: tolerance for fitting procedure
%  - Ninterp: Number of interpolation points for 2F1 approximation [ 1 x 1]
%  - guess: [GAMMA,L] or [GAMMA,L,alphaEps] is the first guess of the 
% coefficients to be fitted. By default it is a [ 1 x 3 ] vector
%  - alphaEps: if unknown, it is empty by default. Otherwise, it is user's
%  defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alphaEps = 1st constant of Mann spectral tensor [1 x 1]
% GAMMA = shear constant (2nd constant) [1 x 1]
% L = Integral length scale [1 x 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author and version:
%  Etienne cheynet  - last modified:  16/04/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see also MannTurb.m MannCoherence.m

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT parser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
p.addOptional('Ninterp',100);
p.addOptional('tolX',1e-3);
p.addOptional('tolFun',1e-3);
p.addOptional('guess',[3,20,0.1]);
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
Niter = p.Results.Niter ;
guess = p.Results.guess ;
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Ninterp = p.Results.Ninterp ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVE NUMBER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double sided k2 and k3 with logarithmically spaced interval points
k2_log=[-fliplr(logspace(k2min,k2max,N2)),logspace(k2min,k2max,N2)];
k3_log=[-fliplr(logspace(k3min,k3max,N3)),logspace(k3min,k3max,N3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a vector kTot used for nlinfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(N1),    N1 = numel(k11);end


%% Check for incomplete data
if isempty(Su)
    flagU = nan;
else
    flagU = 1;
end
if isempty(Sv)
    flagV = nan;
else
    flagV = 1;
end
if isempty(Sw)
    flagW = nan;
else
    flagW = 1;
end
if isempty(Suw)
    flagUW = nan;
else
    flagUW = 1;
end


if N1~=numel(k11)
    dummyK = k11;
    k11 = logspace(log10(dummyK(1)),log10(dummyK(end)),N1);
    
    if isempty(Su)
        Su = nan(1,N1);
    else
        Su = interp1(dummyK,Su,k11,'cubic');
    end

    if isempty(Sv)
        Sv = nan(1,N1);
    else
        Sv = interp1(dummyK,Sv,k11,'cubic');
    end

    if isempty(Sw)
        Sw = nan(1,N1);
    else
        Sw = interp1(dummyK,Sw,k11,'cubic');
    end

    if isempty(Suw)
        Suw = nan(1,N1);
    else
        Suw = interp1(dummyK,Suw,k11,'cubic');
    end

else
    if isempty(Su), Su = nan(1,N1); end
    if isempty(Sv), Sv = nan(1,N1); end
    if isempty(Sw), Sw = nan(1,N1); end
    if isempty(Suw), Suw = nan(1,N1); end
end

%% Concatenate spectra
kTot = zeros(N1,3,3);
for ii=1:3
    for jj=1:3
        kTot(:,ii,jj) = k11;
    end
end

kTot = reshape(kTot,[],1); % is [9 x Ndk1,1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D to 1D transformation (numerical trick)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformation of 3D data into 1D
clear S
S(:,1,1)= k11(:).*Su(:);
S(:,2,2)= k11(:).*Sv(:);
S(:,3,3)= k11(:).*Sw(:);
S(:,1,3)= k11(:).*real(Suw(:));
S(:,3,1)=S(:,1,3);
S(:,1,2)= nan;
S(:,2,3)= nan;
S(:,2,1)=nan;
S(:,3,2)=nan;
S = reshape(S,[],1); % is [9 x Ndk1,1]

indNan = find(isnan(S));
S(indNan) = [];
kTot(indNan) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=optimset('TolX',tolX,'TolFun',tolFun,'Display','iter');
modelFun1 = @ChouguleTurb1; % transform a nested function into anonymous function

if isempty(zeta0)
    guess1 = [guess,0.05];
    minLim = [0,1,0,-0.1];
    maxLim = [6,200,2,1.5];
else
    guess1 = guess;
    minLim = [0,1,0];
    maxLim = [6,200,2];
end


Coeff = lsqcurvefit(@(para,kTot) modelFun1(para,kTot),guess1,kTot,S,minLim,maxLim,options);

GAMMA = abs(Coeff(1));
L = abs(Coeff(2));
alphaEps = abs(Coeff(3));
if isempty(zeta0)
    zeta = abs(Coeff(4));
else
    zeta = zeta0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [FM,k1] = ChouguleTurb1(para,kTot)
        kmax = max(kTot);
        kmin = min(kTot);
        GAMMA0 = abs(para(1));
        L0 = abs(para(2));
        alphaEps0 = abs(para(3));
        if isempty(zeta0)
            zeta1 = (para(4));
        else
            zeta1 = zeta0;
        end
        [PHI,k1,k2_log,k3_log] = ChouguleTurb(alphaEps0,GAMMA0,L0,zeta1,'Niter',Niter,'N1',N1,'Ninterp',Ninterp,'k1max',log10(kmax),'k1min',log10(kmin));
        FM0= squeeze(trapz(k3_log,trapz(k2_log,PHI,2),3));
        FM0 = FM0(end-N1+1:end,:,:);
        FM = nan(size(FM0));
        FM(:,1,1) = k1(:).*FM0(:,1,1)*flagU;
        FM(:,2,2) = k1(:).*FM0(:,2,2)*flagV;
        FM(:,3,3) = k1(:).*FM0(:,3,3)*flagW;
        FM(:,1,3) = k1(:).*FM0(:,1,3)*flagUW;
        FM(:,3,1) = k1(:).*FM0(:,3,1)*flagUW;
        FM = reshape(FM,[],1);
        indNan = find(isnan(FM));
        FM(indNan) = [];
    end
end