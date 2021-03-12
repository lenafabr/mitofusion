% Steady-state results for analytical CoG model
% by A Agrawal and EF Koslover, 2020

% input: structure containing optional parameters
% outputs:
% opt = actual parameters used, including defaults
% Pmx = density of anterograde-moving protein
% Mmx = density of retrograde-moving protein
% Sm = proteins in each discrete site
% regpos = discrete site positions
% xplot= x positions at which densities are calculated
% dx = binsize for x

function [opt,Pmx,Mmx,Sm,regpos,xplot,dx] = calcStopRestart_localtrans(options)
%%
% declare default parameters
opt.cm = 1; % protein per mitochondria
opt.ps = 0.4; % stopping probability at sinks
opt.v = 1; % velocity of moving mito
opt.kd = 0.06; % protein degradation rate
opt.kw = 0.1; % restarting prob at sinks
opt.L = 1; % total domain length
opt.M = 1500; % total number of mitochondria
opt.Nregions = 10; % number of regions the domain is divided into

opt.npt = 1000; %discretization while plotting

opt.maxS = 0; %use when setting a limit to how many can stop 

% local translation parameters
% if have certain %age 
opt.trans_frac = 0.2;
%opt.r = 0.05; %prefactor for local translation term

% copy over input options
if (exist('options')==1)
    opt = copyStruct(options, opt);
end

%opt.Nsinks = opt.Nregions; 
% redeclare for simple reading of code
del = opt.L/(opt.Nregions+1); % distance between regions
v = opt.v; % velocity of protein transport
kd = opt.kd;
kw = opt.kw;
ps = opt.ps;
L = opt.L;
%kp = unitprot; % production rate over velocity
%Nsinks = opt.Nsinks;
Nregions = opt.Nregions;
cm = opt.cm;
M = opt.M;
P0 = cm * M / (2*(L+opt.Nregions*v*ps/kw));

% derived known parameters
B = kw*ps/(2*(kd+kw)); %Beta, factor in equations
a = kd*del/v; % alpha, exponent used in equation

% calculate required rate of production to contribute to trans_frac
%r = opt.trans_frac * (L*kw+opt.Nregions*v*ps)*(kw+kd) /(v*ps*Nregions);
% below: set %age of total influx into the domain due to transport from soma 
r = opt.trans_frac * kw /(2*v*ps*Nregions); 
opt.r = r;

%% boundary conditions
% order of elements: P1 M1 P2 M2 ... Pn Mn
% Pn : positive moving mito in n region, Mn = neg moving
nMat = 2*(Nregions+1); %total number of elements
constmatrix = zeros(nMat,1); % column vector of rhs of matrix eqn
pmmatrix = zeros(nMat); %declare a zero matrix

% declare edge boundary conditions
% constant influx of proteins
pmmatrix(1,1) = 1;
constmatrix(1,1) = P0; 
localterm = opt.r*v*ps*M/(L*kw+opt.Nregions*v*ps)/(kd+kw); %local translation term
constmatrix(2:end-1) = 1/2*localterm*kw;



% iteratively declare rest of the elements, going region by region
% better way to do this?
rowstart = 2;
columnstart = 1;



for m = 1:Nregions
    % P eqn
    pmmatrix(rowstart,columnstart:columnstart+3) = [-(1- ps + B) * exp(-a),0,1,-B];
    % M eqn
    pmmatrix(rowstart+1,columnstart:columnstart+3) = [-B*exp(-a),exp(a),0,-(1- ps + B)];
    % redeclare rowstart and columnstart
    rowstart = rowstart+2;
    columnstart = columnstart+2;
end

% reflecting boundary condition
pmmatrix(end,end-1) = 1;
pmmatrix(end,end) = -exp(2*a);

% solution for parameter set
solnvector = pmmatrix\constmatrix;

%% function values

% function values by sink
Pm0 = solnvector(1:2:end-1);
Mm0 = solnvector(2:2:end);


Sm = v*ps/(kd+kw)* (Pm0(1:end-1) .* exp(-a) + Mm0(2:end)) + localterm ;

%%
ptsperregion = ceil((opt.npt-1)/(opt.Nregions+1));
dx = del/ptsperregion;
xreg = linspace(dx,del,ptsperregion);


xplot = zeros(1,ptsperregion*(opt.Nregions+1)+1);
Pmx = xplot;
Mmx = xplot;

xplot(1) = 0;
Pmx(1) = Pm0(1);
Mmx(1) = Mm0(1);

regind = 1 + ptsperregion*(1:Nregions);
%%
ind = 1;
for j = 1:Nregions + 1 % cycle over regions
    startind = 2 + (j-1)*ptsperregion;
    Pmx(startind:startind+ptsperregion-1) = Pm0(ind) * exp(-kd*xreg/v);
    Mmx(startind:startind+ptsperregion-1) = Mm0(ind) * exp(kd*xreg/v);
    xplot(startind:startind+ptsperregion-1) = xreg+del*(j-1);
    ind = ind + 1;
    
end

% site positions
regpos = xplot(regind);


end