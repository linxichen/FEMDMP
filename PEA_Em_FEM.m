%% Housekeeping
clear
close all
clc
format long
addpath(genpath('./tools'))
addpath(genpath('./param'))

%% Set the stage
mypara;
nA = 9;
nK = 50;
nN = 50;
[P,lnAgrid] = rouwen(rrho,0,ssigma/sqrt(1-rrho^2),nA);
Anodes = exp(lnAgrid);
P = P';
min_lnA = lnAgrid(1); max_lnA = lnAgrid(end);
min_K = 300; max_K = 4000;
min_N = 0.7; max_N = 0.99;
damp_factor = 0.0;
maxiter = 10000;
tol = 1e-6;
options = optimoptions(@fsolve,'Display','none','Jacobian','off');

%% Grid creaton
k_rate = (log(max_K)-log(min_K))/nK;
n_rate = (log(max_N)-log(min_N))/nN;
Knodes = zeros(1,nK); Nnodes = zeros(1,nN);
Knodes(1) = min_K; Nnodes(1) = min_N;
for i = 2:nK
    Knodes(i) = (1+k_rate)*Knodes(i-1);
end
for i = 2:nN
    Nnodes(i) = (1+n_rate)*Nnodes(i-1);
end
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
tot_stuff = zeros(N,1); ustuff = zeros(N,1);
EMHval = zeros(nA,nK,nN); EMFval = EMHval;
EMHval_temp = EMHval; EMFval_temp = EMFval;
parfor i = 1:N
	[i_a,i_k,i_n] = ind2sub([nA,nK,nN],i);
	a = Anodes(i_a); k  = Knodes(i_k); n = Nnodes(i_n); %#ok<PFBNS>
	tot_stuff(i) = a*k^aalpha*n^(1-aalpha) + (1-ddelta)*k + z*(1-n);
	ustuff(i) = xxi*(1-n)^(1-eeta);
end
if (exist('PEA_Em_FEM.mat','file'))
    load('PEA_Em_FEM.mat','EMHval','EMFval');
end
if isequal(size(EMHval),[nA nK nN]) ~= 1
	EMHval = zeros(nA,nK,nN); EMFval = EMHval;
	EMHval_temp = EMHval; EMFval_temp = EMFval;
    coeff_lnmh = zeros(4,1); coeff_lnmf = zeros(4,1);
    coeff_lnmh(1) = 2.247337592951108;
    coeff_lnmh(2) = -0.041544383160081;
    coeff_lnmh(3) = -0.607644008294761;
    coeff_lnmh(4) = -0.004314696290213;
    
    coeff_lnmf(1) = 2.351435745790115;
    coeff_lnmf(2) = 2.203515288267346;
    coeff_lnmf(3) = -0.364368568546649;
    coeff_lnmf(4) = -0.011952817385299;
    parfor i = 1:N
        [i_a,i_k,i_n] = ind2sub([nA,nK,nN],i);
        a = Anodes(i_a); k  = Knodes(i_k); n = Nnodes(i_n); %#ok<PFBNS>
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
            % warning('q <= 0!!')
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
			if qplus <= 0
				% warning('qplus <= 0!!')
				qplus = 0;
				tthetaplus = 0;
				vplus = 0;
				EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
				EMF_hat = EMF_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) )/cplus );
			else
				tthetaplus = (qplus/xxi)^(1/(eeta-1));
				vplus = tthetaplus*(1-nplus);
				EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
				EMF_hat = EMF_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) + (1-x)*kkappa/qplus - ttau*kkappa*tthetaplus )/cplus );
			end
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

end;

%% Euler equation error
nk_ee = 60; nnn_ee = 60;
Kgrid = linspace(0.5*k_ss,1.5*k_ss,nk_ee);
Agrid = exp(lnAgrid);
Ngrid = linspace(0.96*n_ss,1.04*n_ss,nnn_ee);
EEerror_c = 999999*ones(nA,nk_ee,nnn_ee);
EEerror_v = 999999*ones(nA,nk_ee,nnn_ee);
cc = zeros(nA,nk_ee,nnn_ee);
vv = zeros(nA,nk_ee,nnn_ee);
tthetattheta = zeros(nA,nk_ee,nnn_ee);
cc_dynare = cc;
vv_dynare = vv;
tthetattheta_dynare = tthetattheta; 

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
            cc_dynare(i_a,i_k,i_n) = exp( 2.130385+0.039519*(log(a)/rrho-0)+0.606879*(log(k)-log(k_ss))+0.005573*(log(n)-log(n_ss)) );
            vv(i_a,i_k,i_n) = v;
            vv_dynare(i_a,i_k,i_n) = exp( -2.899249+3.417972*(log(a)/rrho-0)+0.451375*(log(k)-log(k_ss))+(-17.928147)*(log(n)-log(n_ss)) );
            tthetattheta(i_a,i_k,i_n) = ttheta;
            tthetattheta_dynare(i_a,i_k,i_n) = exp( 0+3.417972*(log(a)/rrho-0)+0.451375*(log(k)-log(k_ss))+(-0.767653)*(log(n)-log(n_ss)) );

			% Find expected mh, mf tomorrow if current coeff applies tomorrow
            EMH_hat = 0; EMF_hat = 0;
            for i_node = 1:nA
                aplus = Anodes(i_node);
                EMH_plus = globaleval(kplus,nplus,Knodes,Nnodes,squeeze(EMHval(i_node,:,:)));
                EMF_plus = globaleval(kplus,nplus,Knodes,Nnodes,squeeze(EMFval(i_node,:,:)));
                cplus = 1/(bbeta*EMH_plus);
                qplus = kkappa/cplus/(bbeta*EMF_plus);
				if qplus <= 0
					% warning('qplus <= 0!!')
					qplus = 0;
					tthetaplus = 0;
					vplus = 0;
					EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
					EMF_hat = EMF_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) )/cplus );
				else
					tthetaplus = (qplus/xxi)^(1/(eeta-1));
					vplus = tthetaplus*(1-nplus);
					EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
					EMF_hat = EMF_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) + (1-x)*kkappa/qplus - ttau*kkappa*tthetaplus )/cplus );
				end
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

v_policy = figure;
plot(Ngrid,squeeze(vv(1,1,:)))
title('Vacancy policy at lowerest productivity and capital.')
print(v_policy,'-dpsc','./results/v_policy.eps')

c_policy = figure;
plot(Kgrid,squeeze(cc(ceil(nA/2),:,ceil(nnn_ee/2))),Kgrid,squeeze(cc_dynare(ceil(nA/2),:,ceil(nnn_ee/2))))
title('Consumption policies at SS.')
print(c_policy,'-dpsc','./results/c_policy.eps')
xlabel('Capital')

ttheta_policy = figure;
plot(Kgrid,squeeze(tthetattheta(ceil(nA/2),:,ceil(nnn_ee/2))),Kgrid,squeeze(tthetattheta_dynare(ceil(nA/2),:,ceil(nnn_ee/2))))
title('\theta around SS')
print(ttheta_policy,'-dpsc','./results/ttheta_policy.eps')
xlabel('Capital')

ttheta_policyN = figure;
plot(Ngrid,squeeze(tthetattheta(ceil(nA/2),ceil(nk_ee/2),:)),Ngrid,squeeze(tthetattheta_dynare(ceil(nA/2),ceil(nk_ee/2),:)))
title('\theta around SS')
print(ttheta_policyN,'-dpsc','./results/ttheta_policy2.eps')
xlabel('Employment')

ttheta_policyA = figure;
plot(Anodes,squeeze(tthetattheta(:,ceil(nk_ee/2),ceil(nnn_ee/2))),Anodes,squeeze(tthetattheta_dynare(:,ceil(nk_ee/2),ceil(nnn_ee/2))))
title('\theta around SS')
xlabel('Productivity')
print(ttheta_policyA,'-dpsc','./results/ttheta_policy3.eps')

save('PEA_Em_FEM.mat');
