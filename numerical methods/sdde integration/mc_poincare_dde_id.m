function ip = mc_poincare_dde_id(h,g,tu,U,Ud,nPmax,dpmin)
%MC_POINCARE_DDE Extract Poincare section crossings from the output of
%Monte Carlo simulations (only approximate indices)
% Input:
%   h: Poincare section definition h(t,x,xd) = 0
%   g: inequality condition g(t,x,xd)>0 at zero crossings
%   tu: simulation time (equidistant mesh)
%   U: ensamble of generated stochastic solution signals
%   Ud: ensamble of generated stochastic delayed term signals
%   nPmax: maximum number of crossings detected per realization 
%       (default 100)
%   dpmin: minimum number of timesteps between subsequent crossings
%       (default 10)
% Output:
%   ip: indices of Poicare section crossings

if nargin < 6
    nPmax = 100; % defualt crossing count limit
end
if nargin < 7
    dpmin = 10; % defualt crossing distance limit
end

% Initialize output array
d = size(U,1);
n = size(U,3);
nt = size(Ud,4);
ip = NaN(nPmax,n);

% Go trajectory by trajectory
niPmax = 0;
for i = 1:n
    % Find zero crossings of h
    Ui = squeeze(U(:,:,i)); 
    Udi = reshape(squeeze(Ud(:,:,i,:)),d,[],nt);
    hi = h(tu,Ui,Udi);
    shi = sign(hi(1:end-1)).*sign(hi(2:end)); 
    pid = find(shi == -1); % detect sign changes in h
    if ~isempty(pid)
        gpi = g(tu(pid),Ui(:,pid),Udi(:,pid,:)); % inequlaity conditions
        pid(gpi<=0) = []; % delete zero crossings that do not satisfy g(x,xd)>0
        pid(diff(pid)<dpmin) = []; % delete zero crossings that are too close
        pid(nPmax+1:length(pid)) = []; % delete surplus crossings
    end
    if length(pid) > niPmax
        niPmax = length(pid);
    end
    ip(1:length(pid),i) = pid;

end

% Delete surplus columns
ip = ip(1:niPmax,:);

end