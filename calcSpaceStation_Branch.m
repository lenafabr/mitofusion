% Steady-state results for analytical SS model
% by A Agrawal and EF Koslover, 2020

% input: structure containing optional parameters
% outputs:
% Smreg = health in each discrete sites
% opt = actual parameters used, including defaults
% Nm = number of mitochondria stationed in each discrete site
% Pmx = density of anterograde-moving protein
% Mmx = density of retrograde-moving protein
% regpos = discrete site positions
% xplot= x positions at which densities are calculated
% dx = binsize for x

function [Smreg,opt,Nm,Pmx,Mmx,regpos,xplot,dx] = calcSpaceStation_Branch(options)
%%
% declare default parameters
opt.pf = 0.3; % protein exchange probability
opt.v = 1; % velocity of moving mito
opt.kd = 0.1; % protein degradation rate
opt.L = 1; % total domain length

opt.Nsinks = 100; %number of sinks along the region
opt.Nregions = 4; % number of regions the domain is divided into
% for a branched network, this is just the number of regions along one
% trajectory of mitochondria. 
%Actual number of regions = (ndeg^(Nregions)-1)/(ndeg-1)

Nsinkperreg = opt.Nsinks/opt.Nregions; % make sure this is a whole number
opt.npt = 1000; %discretization while plotting
opt.M = 1500; % total number of mitochondria 
opt.sinkinreginput = 0;
% declare branching degree
opt.ndeg = 1;

% copy over input options
if (exist('options')==1)
    opt = copyStruct(options, opt);
end

%opt.rho = (opt.M - opt.Nsinks)/opt.L;
% declare rho as option
if (opt.ndeg == 1)
     opt.rho = (opt.M - opt.Nsinks)/opt.L;
else
     opt.rho = (opt.M - (opt.Nsinks/opt.Nregions)*(opt.ndeg^(opt.Nregions)-1)/(opt.ndeg-1))/opt.L;
end

%opt.rho 

% redeclare for simple reading of code
del = opt.L/(opt.Nregions+1); % distance between regions
v = opt.v; % velocity of protein transport
kd = opt.kd;
%kp = unitprot; % production rate over velocity
Nsinks = opt.Nsinks;
Nregions = opt.Nregions;
rho = opt.rho;
pf = opt.pf;
ndeg = opt.ndeg;


% derived known parameters
%B = pf * v * rho / (4 * (kd + v*rho*pf/2)); %Beta, factor in equations
a = kd*del/v; % alpha, exponent used in equation

%% position sinks in nregions
if (opt.sinkinreginput == 0) %if no predefined distribution exists
    clear sinkinreg
    sinkinreg = randi(Nregions,1,Nsinks);
    Nm = zeros(1,Nregions); % Nm gives the number of sink mito in a given region
%     lastsink = 1;
%     for i = 1:Nregions
%         Nm(i) = Nsinkperreg * ndeg^(i-1); 
%         sinkinreg(lastsink:lastsink+Nm(i)-1) = i;
%         lastsink = lastsink+Nm(i);
%     end
        
else
    sinkinreg = options.sinkinreginput;
end
sinkinreg = sort(sinkinreg);
for k = 1:Nregions
    Nm(k) = sum(sinkinreg == k);
end
% reassign total number of sinks as lastsink
%Nsinks = lastsink-1;
%opt.Nsinks = Nsinks;
%%
% get spacing between individual sinks
dellist = del * (sinkinreg(2:end) - sinkinreg(1:end-1)); 
alist = [sinkinreg(1) * del * kd /v  kd * dellist / v] ; % list of multipliers for individual sinks
% Nm = ones(1,Nregions); %default with 1 sink in every region

%% boundary conditions
% order of elements: P1 M1 P2 M2 ... Pn Mn
% Pn : positive moving mito in n region, Mn = neg moving
nMat = 2*(Nsinks+1); %total number of elements
constmatrix = zeros(nMat,1); % column vector of rhs of matrix eqn
pmmatrix = zeros(nMat); %declare a zero matrix

% declare edge boundary conditions
% constant influx of proteins
pmmatrix(1,1) = 1;
constmatrix(1,1) = rho/2; 



% iteratively declare rest of the elements, going region by region
% better way to do this?
rowstart = 2;
columnstart = 1;



for m = 1:Nsinks
    currReg = sinkinreg(m); % region which the sink belongs to
    rhoN = rho/ndeg^(currReg-1); % scaled rho for input
    B = pf * v *rhoN/ (4 * (kd + v*rhoN*pf/2)); %Beta, factor in equations
    if (alist(m) > 0)
        % P eqn
        pmmatrix(rowstart,columnstart:columnstart+3) = [-(1- (pf/2) + B*pf/2) * exp(-alist(m))/ndeg,0,1,-pf*B/2/ndeg];
        % M eqn
        pmmatrix(rowstart+1,columnstart:columnstart+3) = [-pf*B*exp(-alist(m))/2/ndeg,exp(alist(m))/ndeg,0,-(1- (pf/2) + pf*B/2/ndeg)];
    else
        % P eqn
        pmmatrix(rowstart,columnstart:columnstart+3) = [-(1- (pf/2) + B*pf/2) * exp(-alist(m)),0,1,-pf*B/2];
        % M eqn
        pmmatrix(rowstart+1,columnstart:columnstart+3) = [-pf*B*exp(-alist(m))/2,exp(alist(m)),0,-(1- (pf/2) + pf*B/2)];
    end
    % redeclare rowstart and columnstart
    rowstart = rowstart+2;
    columnstart = columnstart+2;
end

% reflecting boundary condition
pmmatrix(end,end-1) = 1;
pmmatrix(end,end) = -exp(2*(opt.L - sinkinreg(end)*del)*kd/v);

% solution for parameter set
solnvector = pmmatrix\constmatrix;

%% function values

% function values by sink
Pm0 = solnvector(1:2:end-1);
Mm0 = solnvector(2:2:end);
Sm = zeros(Nsinks,1);
%Sm = B*2/rho * (Pm0(1:end-1) .* exp(-alist)' + Mm0(2:end));

for m = 1:Nsinks
    currReg = sinkinreg(m); % region which the sink belongs to
    rhoN = rho/ndeg^(currReg-1);
    B = pf * v *rhoN / (4 * (kd + v*rhoN*pf/2));
    Sm(m) = B*2/rhoN * (Pm0(m) .* exp(-alist(m)) + Mm0(m+1));
end
    
Smreg = zeros(1,Nregions);

%%
for reg = 1:Nregions
    Smreg(reg) = sum(Sm' .* (sinkinreg == reg));
end 
nzNm = sum(Nm>0);
Pmreg = zeros(1,sum(Nm>0) + 1); %exists over regions where there are sinks
Mmreg = Pmreg;
Pmreg(1) = Pm0(1);
Mmreg(1) = Mm0(1);
ind = 2;
for reg = 1:Nregions
    if (Nm(reg)>0) % if there are sinks positioned
        Pmreg(ind) = sum(Pm0(2:end)' .* (sinkinreg == reg)) ./ sum(sinkinreg == reg);
        Mmreg(ind) = sum(Mm0(2:end)' .* (sinkinreg == reg))./ sum(sinkinreg == reg);
        ind = ind+1;
    end
end
%%
ptsperregion = ceil((opt.npt-1)/(opt.Nregions+1));
dx = del/ptsperregion;
xreg = linspace(dx,del,ptsperregion);


xplot = zeros(1,ptsperregion*(opt.Nregions+1)+1);
Pmx = xplot;
Mmx = xplot;

xplot(1) = 0;
Pmx(1) = Pmreg(1);
Mmx(1) = Mmreg(1);

regind = 1 + ptsperregion*(1:Nregions);

%%
ind = 1;
for j = 1:Nregions + 1 % cycle over regions 
    startind = 2 + (j-1)*ptsperregion;
    if ((j >1 && Nm(j-1) > 0) || (j == 1))
        Pmx(startind:startind+ptsperregion-1) = Pmreg(ind) * exp(-kd*xreg/v);
        Mmx(startind:startind+ptsperregion-1) = Mmreg(ind) * exp(kd*xreg/v);
        xplot(startind:startind+ptsperregion-1) = xreg+del*(j-1);
        ind = ind + 1;
    else % extrapolate from previous region
        Pmx(startind:startind+ptsperregion-1) =  Pmx(startind-1) * exp(-kd*xreg/v);
        Mmx(startind:startind+ptsperregion-1) =  Mmx(startind-1) * exp(kd*xreg/v);
        xplot(startind:startind+ptsperregion-1) = xreg+del*(j-1);
        
    end        
            
end

regpos = xplot(regind);



end