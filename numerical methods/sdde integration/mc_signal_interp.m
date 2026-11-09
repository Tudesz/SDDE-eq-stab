function Ue = mc_signal_interp(t,U,Ud,tau,th,dmap)
%MC_SIGNAL_INTERP Interpolate on the output of MC simulations
% Input:
%   t: time signal of the solution
%   U: ensamble of generated stochastic signals
%   Ud: ensamble of generated stochastic delayed term signals
%   tau: fixed, point time delays
%   th: querry points in history (th>-tau_max!)
%   dmap: boolean vector which marks with ones the n0 relevant delayed
%       state variables
% Output:
%   Ue: interpolated ensamble of simulation outputs

% Initialization
[tau_max,id_max] = max(tau);
Ud_max = Ud(:,:,:,id_max);
[~,N,n] = size(Ud_max);
de = length(find(dmap==1));
M = length(th);

% Interpolation on Ud (nearest neighbour, zeroth order)
dt = t(end)/(N-1);
Ue = zeros(de*M,N,n);
for i = 1:N
    for j = 1:M
        tdel = t(i) - th(j);
        if tdel > t(1)
            id = i-ceil(th(j)/dt);
            Ue((j-1)*de+(1:de),i,:) = U(dmap==1,id,:);
        else
            id = i-ceil((th(j)-tau_max)/dt);
            id = min([id,size(Ud_max,2)]);
            id = max([id,1]);
            Ue((j-1)*de+(1:de),i,:) = Ud_max(dmap==1,id,:);
        end
    end
    
end

end