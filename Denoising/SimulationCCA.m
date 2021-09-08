%% Simulation to test CCA 
% See appendix

%% set parameters

% basic parameters
Npix_out = 500; % dimensionality of 'out' (in one slice)
Npix_in = 500; % dimensionality of 'in' (in one slice)
Nsnd = 200; % number of sounds
Nrep = 4; % number of repeats
T = 20; % length of a trial
Tbase = 7; % baseline period before sound onset

% a bunch of parameters for data creation
Ncomp_a = 20; % number of artefact components
p_List = .0:.05:1.0; % loop over how much 'a_k' is going to be correlated to 's' (if high, high correlation)
noisesigma = 1; % noise scaling factor for voxels

% cca params
cca_params.pc2Keep = 250;
cc2Remove_List = 20;
cca_params.baseline_tps = 1:Tbase;

%% get data parts (anything that is randomized)

% shared components in response to sounds
Ncomp_u = 10; % dimensionality of the responses to sound
u = randn(T*Nsnd,Ncomp_u);

% movement components ('a_k')
n_s2a = randn(Ncomp_u,Ncomp_a); % resample weights everytime
% a is the reliable, sound-evoked part, b is non-reliable part
a_s = u*n_s2a;
a_s = repmat(a_s,[Nrep,1]);

b_k = randn(T*Nsnd*Nrep,Ncomp_a);

% get 'out' data parts
w_a2dout = randn(Ncomp_a,Npix_out);
Dout_n = randn(T*Nsnd*Nrep,Npix_out); % voxel noise

% get 'in' data parts
m_s2din = randn(Ncomp_u,Npix_in);
Din_s = u*m_s2din;
Din_s = repmat(Din_s,[Nrep,1]);

w_a2din = randn(Ncomp_a,Npix_in);
Din_n = randn(T*Nsnd*Nrep,Npix_in); % voxel noise

%% loop

NSE1 = nan(numel(p_List), numel(cc2Remove_List));
NSE2 = nan(numel(p_List), numel(cc2Remove_List));
NSE3 = nan(numel(p_List), numel(cc2Remove_List));
NSE_noiseceiling = nan(numel(p_List), numel(cc2Remove_List));

for idx_p = 1:numel(p_List)
    % loop over how much a_k is correlated with u
    p = p_List(idx_p);
    
    % assemble parts
    a_k = p*a_s + (1-p)*b_k;
    
    Dout_a = a_k*w_a2dout;

    Dout = Dout_a + noisesigma*Dout_n;
    Dout = reshape(Dout, [T,Nsnd,Nrep,Npix_out]);
    
    Din_a = a_k*w_a2din;
    
    Din = Din_s + Din_a + noisesigma*Din_n;
    Din = reshape(Din, [T,Nsnd,Nrep,Npix_in]);
    
    for idx_cc = 1:numel(cc2Remove_List)
        % loop over number of components to remove
        
        cca_params.cc2Remove = cc2Remove_List(idx_cc);
        
        disp([' p: ' num2str(p) '/ cc: ' num2str(cca_params.cc2Remove)])
        
        % perform CCA with recentering
        cca_params.recenter = 1;
        D_corr_r1 = CCAcorrection(Din,Dout,cca_params);
        
        % try vanilla CCA on non-centered data        
        cca_params.recenter = 0;
        D_corr_r0 = CCAcorrection(Din,Dout,cca_params);
        
        % get nse
        realSig = Din_s; % here put nsf*Din_n to have the lower floor=0
        NSE1(idx_p,idx_cc) = median(NSE(realSig,reshape(D_corr_r1,[T*Nsnd*Nrep,Npix_in]),1));
        NSE2(idx_p,idx_cc) = median(NSE(realSig,reshape(D_corr_r0,[T*Nsnd*Nrep,Npix_in]),1));
        NSE3(idx_p,idx_cc) = median(NSE(realSig,reshape(Din,[T*Nsnd*Nrep,Npix_in]),1));
        NSE_noiseceiling(idx_p,idx_cc) = median(NSE(realSig,Din_s + noisesigma*Din_n,1));
        
    end
end

%% plot

figure; hold all
plot(p_List, NSE1(:,idx_cc), 'LineWidth', 2)
plot(p_List, NSE2(:,idx_cc), 'LineWidth', 2)
plot(p_List, NSE3(:,idx_cc), 'LineWidth', 2)
plot(p_List, NSE_noiseceiling(:,idx_cc),'k--', 'LineWidth', 2)
xlabel('p')
ylabel('NSE')
legend({'rCCA','CCA','raw','noise ceiling'}, 'Location', 'NorthWest'); legend boxoff

%%
idx_p = 1;

figure; hold all
plot(cc2Remove_List,NSE1(idx_p,:),'Linewidth',2); 
plot(cc2Remove_List,NSE2(idx_p,:),'Linewidth',2);
plot(cc2Remove_List,NSE3(idx_p,:),'Linewidth',2); 
plot(cc2Remove_List,NSE_noiseceiling(idx_p,:),'k--', 'LineWidth', 2); 
vline(Ncomp_a,'k--')
xlabel('Number of components removed')
ylabel('Error')
title(['p: ' num2str(p_List(idx_p))])