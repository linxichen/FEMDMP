%% Housekeeping
clear
close all
clc
format long
addpath(genpath('./tools'))
addpath(genpath('./param'))

%% Set the stage
mypara;
nA = 10;
nK = 30;
nN = 100;
[P,lnAgrid] = rouwen(rrho,0,ssigma/sqrt(1-rrho^2),nA);
Anodes = exp(lnAgrid);
P = P';
min_lnA = lnAgrid(1); max_lnA = lnAgrid(end);
min_K = 500; max_K = 3000;
min_N = 0.8; max_N = 0.98;
damp_factor = 0.8;
maxiter = 10000;
tol = 1e-10;
options = optimoptions(@fsolve,'Display','none','Jacobian','off');

%% Grid creaton
Knodes = linspace(min_K,max_K,nK);
Nnodes = linspace(min_N,max_N,nN);
N = nA*nK*nN;
[Kmesh,Amesh,Nmesh] = meshgrid(Knodes,Anodes,Nnodes);

%% Encapsulate all parameters
param = [... 
 bbeta; % 1
 ggamma; % 2
 kkappa; % 3
 eeta; % 4
 rrho; %5
 ssigma; %6
 x; % 7
 aalpha; % 8
 ddelta; % 9
 xxi; % 10
 ttau; % 11
 z % 12
 ];

%% Precomputation and initial guess
if (exist('PEA_Em_FEM.mat','file'))
    load('PEA_Em_FEM.mat','EMHval','EMFval');
end
if isequal(size(EMHval),[nA nK nN]) ~= 1
    coeff_lnmh = zeros(4,1); coeff_lnmf = zeros(4,1);
    coeff_lnmh(1) = 2.247337592951108;
    coeff_lnmh(2) = -0.041544383160081;
    coeff_lnmh(3) = -0.607644008294761;
    coeff_lnmh(4) = -0.004314696290213;
    
    coeff_lnmf(1) = 2.351435745790115;
    coeff_lnmf(2) = 2.203515288267346;
    coeff_lnmf(3) = -0.364368568546649;
    coeff_lnmf(4) = -0.011952817385299;
    tot_stuff = zeros(N,1); ustuff = zeros(N,1);
    EMHval = zeros(nA,nK,nN); EMFval = EMHval;
    EMHval_temp = EMHval; EMFval_temp = EMFval;
    parfor i = 1:N
        [i_a,i_k,i_n] = ind2sub([nA,nK,nN],i);
        a = Anodes(i_a); k  = Knodes(i_k); n = Nnodes(i_n); %#ok<PFBNS>
        tot_stuff(i) = a*k^aalpha*n^(1-aalpha) + (1-ddelta)*k + z*(1-n);
        ustuff(i) = xxi*(1-n)^(1-eeta);
        EMHval(i) = exp([1 log(a) log(k) log(n)]*coeff_lnmh);
        EMFval(i) = exp([1 log(a) log(k) log(n)]*coeff_lnmf);
    end
end


%% Solve for SS
kss = k_ss;
nss = n_ss;

%% Iteration
diff = 10; iter = 0;
while (diff>tol && iter <= maxiter)
    %% Time iter step, uses endo grid technique
    parfor i = 1:N
        [i_a,i_k,i_n] = ind2sub([nA,nK,nN],i);
        n = Nnodes(i_n); k = Knodes(i_k);
        EMH = globaleval(k,n,Knodes,Nnodes,squeeze(EMHval(i_a,:,:)));
        EMF = globaleval(k,n,Knodes,Nnodes,squeeze(EMFval(i_a,:,:)));
        c = 1/(bbeta*EMH);
        q = kkappa/c/(bbeta*EMF);
        if q <= 0
            warning('q <= 0!!')
            q = 0;
            ttheta = 0;
            v = 0;
            kplus = tot_stuff(i) - c - kkappa*v;
            nplus = (1-x)*n;
        else
            ttheta = (q/xxi)^(1/(eeta-1));
            v = ttheta*(1-n);
            kplus = tot_stuff(i) - c - kkappa*v;
            nplus = (1-x)*n + xxi*v^eeta*(1-n)^(1-eeta);
        end

        EMH_hat = 0; EMF_hat = 0;
        for i_node = 1:nA
            aplus = Anodes(i_node);
            EMH_plus = globaleval(kplus,nplus,Knodes,Nnodes,squeeze(EMHval(i_node,:,:)));
            EMF_plus = globaleval(kplus,nplus,Knodes,Nnodes,squeeze(EMFval(i_node,:,:)));
            cplus = 1/(bbeta*EMH_plus);
            qplus = kkappa/cplus/(bbeta*EMF_plus);
            tthetaplus = (qplus/xxi)^(1/(eeta-1));
            EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
            EMF_hat = EMF_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) + (1-x)*kkappa/qplus - ttau*kkappa*tthetaplus )/cplus );
        end
        
        EMHval_temp(i) = EMH_hat;
        EMFval_temp(i) = EMF_hat;
    end
    
    %% Damped update
    EMHval_new = (1-damp_factor)*EMHval_temp+(damp_factor)*EMHval;
    EMFval_new = (1-damp_factor)*EMFval_temp+(damp_factor)*EMFval;
    
    %% Compute norm
    diff = norm([EMHval(:);EMFval(:)]-[EMHval_new(:);EMFval_new(:)],Inf);
    
    %% Update
    EMHval = EMHval_new;
    EMFval = EMFval_new;
    iter = iter+1;
    %% Display something
    iter
    diff
    coeff_lnmh;
    coeff_lnmf;

end;

%% Euler equation error
nk_ee = 10; nnn_ee = 10;
Kgrid = linspace(1100,1500,nk_ee);
Agrid = exp(lnAgrid);
Ngrid = linspace(0.9,0.97,nnn_ee);
EEerror_c = 999999*ones(nA,nk_ee,nnn_ee);
EEerror_v = 999999*ones(nA,nk_ee,nnn_ee);
cc = zeros(nA,nk_ee,nnn_ee);
vv = zeros(nA,nk_ee,nnn_ee);

for i_a = 1:nA
    a = Agrid(i_a);
    for i_k = 1:nk_ee
        k = Kgrid(i_k);
        for i_n = 1:nnn_ee
            n = Ngrid(i_n);
			tot_stuff = a*k^aalpha*n^(1-aalpha)+(1-ddelta)*k+z*(1-n);
			ustuff = xxi*(1-n)^(1-eeta);
            
            EMH = globaleval(k,n,Knodes,Nnodes,squeeze(EMHval(i_a,:,:)));
            EMF = globaleval(k,n,Knodes,Nnodes,squeeze(EMFval(i_a,:,:)));
            c = 1/(bbeta*EMH);
            q = kkappa/c/(bbeta*EMF);            
            
            if q <= 0
                warning('q <= 0!!')
                q = 0;
                ttheta = 0;
                v = 0;
                kplus = tot_stuff - c - kkappa*v;
                nplus = (1-x)*n;
            else
                ttheta = (q/xxi)^(1/(eeta-1));
                v = ttheta*(1-n);
                kplus = tot_stuff - c - kkappa*v;
                nplus = (1-x)*n + xxi*v^eeta*(1-n)^(1-eeta);
            end
            
            cc(i_a,i_k,i_n) = c;
            vv(i_a,i_k,i_n) = v;

			% Find expected mh, mf tomorrow if current coeff applies tomorrow
            EMH_hat = 0; EMF_hat = 0;
            for i_node = 1:nA
                aplus = Anodes(i_node);
                EMH_plus = globaleval(kplus,nplus,Knodes,Nnodes,squeeze(EMHval(i_node,:,:)));
                EMF_plus = globaleval(kplus,nplus,Knodes,Nnodes,squeeze(EMFval(i_node,:,:)));
                cplus = 1/(bbeta*EMH_plus);
                qplus = kkappa/cplus/(bbeta*EMF_plus);
                tthetaplus = (qplus/xxi)^(1/(eeta-1));
                EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
                EMF_hat = EMF_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) + (1-x)*kkappa/qplus - ttau*kkappa*tthetaplus )/cplus );
            end

			c_imp = 1/(bbeta*EMH_hat);
			q_imp = kkappa/c_imp/(bbeta*EMF_hat);
			ttheta_imp = (q_imp/xxi)^(1/(eeta-1));
			v_imp = ttheta_imp*(1-n);

            EEerror_c(i_a,i_k,i_n) = abs((c-c_imp)/c_imp);   
            EEerror_v(i_a,i_k,i_n) = abs((v-v_imp)/v_imp);  
        end
    end
end
EEerror_c_inf = norm(EEerror_c(:),inf)
EEerror_v_inf = norm(EEerror_v(:),inf)

EEerror_c_mean = mean(EEerror_c(:));
EEerror_v_mean = mean(EEerror_v(:));

%% Export results
mkdir('results')
h_c = figure;
plot(Kgrid,squeeze(EEerror_c(ceil(nA/2),:,ceil(nnn_ee/2))))
title('Euler Error of Consumption')
print(h_c,'-dpsc','./results/EEerror_c.eps')

h_v = figure;
plot(Kgrid,squeeze(EEerror_v(ceil(nA/2),:,ceil(nnn_ee/2))))
title('Euler Error of Vacancy')
print(h_v,'-dpsc','./results/EEerror_v.eps')

result_mf = @(k,n) globaleval(k,n,Knodes,Nnodes,squeeze(EMFval(1,:,:)));
h_EMF = figure;
ezsurf(result_mf,[Kgrid(1),Kgrid(end),Ngrid(1),Ngrid(end)])
print(h_EMF,'-dpsc','./results/EMF.eps')

save('PEA_Em_FEM.mat');
