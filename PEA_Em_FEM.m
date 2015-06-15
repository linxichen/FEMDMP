%% Housekeeping
clear
close all
clc
format long
addpath(genpath('./tools'))
addpath(genpath('./param'))

%% Set the stage
mypara;
nA = 15;
nK = 50;
nN = 50;
[P,lnAgrid] = rouwen(rrho,0,ssigma/sqrt(1-rrho^2),nA);
Anodes = exp(lnAgrid);
P = P';
min_lnA = lnAgrid(1); max_lnA = lnAgrid(end);
min_K = 30; max_K = 5000;
min_N = 0.1; max_N = 0.99;
damp_factor = 0.5;
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
        ttheta = (q/xxi)^(1/(eeta-1));
        v = ttheta*(1-n);
        kplus = tot_stuff(i) - c - kkappa*v;
        nplus = (1-x)*n + xxi*v^eeta*(1-n)^(1-eeta);
        
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
nk = 10; nnn = 10;
Kgrid = linspace(1100,1500,nk);
Agrid = exp(lnAgrid);
Ngrid = linspace(0.9,0.97,nnn);
EEerror_c = 999999*ones(nA,nk,nnn);
EEerror_v = 999999*ones(nA,nk,nnn);

for i_a = 1:nA
    a = Agrid(i_a);
	lna_cheby = lnAchebygrid(i_a);
    for i_k = 1:nk
        k = Kgrid(i_k);
        for i_n = 1:nnn
            n = Ngrid(i_n);
			tot_stuff = a*k^aalpha*n^(1-aalpha)+(1-ddelta)*k+z*(1-n);
			ustuff = xxi*(1-n)^(1-eeta);
			state = [log(a),log(k),log(n),tot_stuff,ustuff];
			lnk_cheby = -1 + 2*(log(k)-min_lnK)/(max_lnK-min_lnK);
			lnn_cheby = -1 + 2*(log(n)-min_lnN)/(max_lnN-min_lnN);
			lnEMH = ChebyshevND(degree,[lna_cheby,lnk_cheby,lnn_cheby])*coeff_lnmh;
			lnEMF = ChebyshevND(degree,[lna_cheby,lnk_cheby,lnn_cheby])*coeff_lnmf;
			c = 1/(bbeta*exp(lnEMH));
			q = kkappa/c/(bbeta*exp(lnEMF));
			v = (q/ustuff)^(1/(eeta-1));
			kplus = tot_stuff - c - kkappa*v;
			nplus = (1-x)*n + q*v;
			lnkplus = log(kplus); lnnplus = log(nplus);
			lnkplus_cheby = -1 + 2*(lnkplus-min_lnK)/(max_lnK-min_lnK);
			lnnplus_cheby = -1 + 2*(lnnplus-min_lnN)/(max_lnN-min_lnN);
			if (lnkplus_cheby < -1 || lnkplus_cheby > 1)
				lnkplus
				error('kplus out of bound')
			end
			if (lnnplus_cheby < -1 || lnnplus_cheby > 1)
				lnnplus_cheby
				lnnplus
				error('nplus out of bound')
			end

			% Find expected mh, mf tomorrow if current coeff applies tomorrow
			EMH_hat = 0;
			EMF_hat = 0;
			for i_node = 1:nA
				lnaplus = lnAgrid(i_node);
				lnaplus_cheby = -1 + 2*(lnaplus-min_lnA)/(max_lnA-min_lnA);
				if (lnaplus_cheby < -1 || lnaplus_cheby > 1)
					error('Aplus out of bound')
				end
				lnEMH_plus = ChebyshevND(degree,[lnaplus_cheby,lnkplus_cheby,lnnplus_cheby])*coeff_lnmh;
				lnEMF_plus = ChebyshevND(degree,[lnaplus_cheby,lnkplus_cheby,lnnplus_cheby])*coeff_lnmf;
				cplus = 1/(bbeta*exp(lnEMH_plus));
				qplus = kkappa/cplus/(bbeta*exp(lnEMF_plus));
				tthetaplus = (qplus/xxi)^(1/(eeta-1));
				EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*exp(lnaplus)*(kplus/nplus)^(aalpha-1))/cplus);
				EMF_hat = EMF_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*exp(lnaplus)*(kplus/nplus)^aalpha-z-ggamma*cplus) + (1-x)*kkappa/qplus - ttau*kkappa*tthetaplus )/cplus );
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

figure
plot(Kgrid,squeeze(EEerror_c(ceil(nA/2),:,ceil(nnn/2))))

figure
plot(Kgrid,squeeze(EEerror_v(ceil(nA/2),:,ceil(nnn/2))))

%% Implied policy functions and find wages
lnAgrid = csvread('../CUDA_VFI/results/Agrid.csv');
lnKgrid = csvread('../CUDA_VFI/results/Kgrid.csv');
lnNgrid = csvread('../CUDA_VFI/results/Ngrid.csv');
nA = length(lnAgrid);
nk = length(lnKgrid);
nnn = length(lnNgrid);

kk = zeros(nA,nk,nnn);
cc = kk;
vv = kk;
nn = kk;
ttheta_export = kk;
wage_export = kk;
cc_dynare = kk;
kk_dynare = kk;
nn_dynare = kk;
vv_dynare = kk;

mmummu = kk;
for i_k = 1:nk
    for i_n = 1:nnn
        for i_A = 1:nA
            state(1) = lnAgrid(i_A); A = state(1);
            state(2) = lnKgrid(i_k); k = state(2);
            state(3) = lnNgrid(i_n); n = state(3);
            EM = exp([1 log(state)]*[coeff_mh coeff_mf]);
            
            y = A*(k)^(aalpha)*(n)^(1-aalpha);
            c = (bbeta*EM(1))^(-1);
            ttheta = (kkappa/(c*xxi*bbeta*EM(2)))^(1/(eeta-1));
            v = ttheta*(1-n);
            mh = (1-ddelta+aalpha*y/k)/c;
            mf = ( (1-ttau)*((1-aalpha)*y/n-z-ggamma*c) + (1-x)*kkappa/xxi*ttheta^(1-eeta) - ttau*kkappa*ttheta )/c;
            w = ttau*A*k^(aalpha)*(1-aalpha)*n^(-aalpha) + (1-ttau)*(z+ggamma*c) + ttau*kkappa*ttheta;
    
            kk(i_A,i_k,i_n) = y - c +(1-ddelta)*k - kkappa*v + z*(1-nn(i_A,i_k,i_n));
            nn(i_A,i_k,i_n) = (1-x)*n + xxi*ttheta^(eeta)*(1-n);
            cc(i_A,i_k,i_n) = c;
            vv(i_A,i_k,i_n) = v;
            
            cc_dynare(i_A,i_k,i_n) = exp(2.111091 + 0.042424/rrho*log(lnAgrid(i_A))/ssigma + 0.615500*(log(lnKgrid(i_k))-log(k_ss)) + 0.014023*(log(lnNgrid(i_n))-log(n_ss)) );
            kk_dynare(i_A,i_k,i_n) = exp(7.206845 + 0.006928/rrho*log(lnAgrid(i_A))/ssigma + 0.997216*(log(lnKgrid(i_k))-log(k_ss)) + 0.005742*(log(lnNgrid(i_n))-log(n_ss)) );
            nn_dynare(i_A,i_k,i_n) = exp(-0.056639 + 0.011057/rrho*log(lnAgrid(i_A))/ssigma + 0.001409*(log(lnKgrid(i_k))-log(k_ss)) + 0.850397*(log(lnNgrid(i_n))-log(n_ss)) );
            
            % Export prices
            wage_export(i_A,i_k,i_n) = w;
            ttheta_export(i_A,i_k,i_n) = ttheta;
        end
    end
end
save('PEA_Em.mat');


i_mid_n = ceil(nnn/2);
i_mid_A = ceil(nA/2);
linewitdh=1.5;
figure
plot(lnKgrid,squeeze(kk(i_mid_A,:,i_mid_n)),lnKgrid,squeeze(kk_dynare(i_mid_A,:,i_mid_n)),'LineWidth',linewitdh)
axis('tight')
xlabel('k(t)')
ylabel('k(t+1)')
legend('Nonlinear','Linear')

figure
plot(lnKgrid,squeeze(nn(i_mid_A,:,i_mid_n)),lnKgrid,squeeze(nn_dynare(i_mid_A,:,i_mid_n)),'LineWidth',linewitdh)
axis('tight')
xlabel('k(t)')
ylabel('n(t+1)')
legend('Nonlinear','Linear')

figure
plot(lnKgrid,squeeze(cc(i_mid_A,:,i_mid_n)),lnKgrid,squeeze(cc_dynare(i_mid_A,:,i_mid_n)),'LineWidth',linewitdh)
axis('tight')
xlabel('k(t)')
ylabel('c(t)')
legend('Nonlinear','Linear')

figure
plot(lnKgrid,squeeze(wage_export(i_mid_A,:,i_mid_n)),'LineWidth',linewitdh)
axis('tight')
xlabel('k(t)')
ylabel('wage')
legend('Nonlinear')

figure
plot(lnKgrid,squeeze(ttheta_export(i_mid_A,:,i_mid_n)),'LineWidth',linewitdh)
axis('tight')
xlabel('k(t)')
ylabel('Tightness')
legend('Nonlinear')

%% Ergodic set where art thou?
    figure
    scatter3(Asim,ksim,nsim)
    xlabel('Productivity')
    ylabel('Capital')
    zlabel('Employment')

%% Dynamics
Aindex = ceil(nA/2);
figure
[Kmesh,Nmesh] = meshgrid(lnKgrid,lnNgrid);
DK = squeeze(kk(Aindex,:,:))-Kmesh';
DN = squeeze(nn(Aindex,:,:))-Nmesh';
quiver(Kmesh',Nmesh',DK,DN,2);
axis tight

%% Paths 1
T = 5000; scale = 0;
A = 0.6;
k1 = zeros(1,T); n1 = zeros(1,T);
k1(1) = 1100; n1(1) = 0.90;
for t = 1:T
    state = [A k1(t) n1(t)];
    EM = exp([1 log(state)]*[coeff_mh coeff_mf]);
    y = A*(k1(t))^(aalpha)*(n1(t))^(1-aalpha);
    c = (bbeta*EM(1))^(-1);
    ttheta = (kkappa/(c*xxi*bbeta*EM(2)))^(1/(eeta-1));
    v = ttheta*(1-n1(t));
    
    if t < T
    k1(t+1) = y - c +(1-ddelta)*k1(t) - kkappa*v;
    n1(t+1) = (1-x)*n1(t) + xxi*ttheta^(eeta)*(1-n1(t));
    end
end
xx = k1; y = n1;
u = [k1(2:end)-k1(1:end-1) 0];
v = [n1(2:end)-n1(1:end-1) 0];

figure
quiver(xx,y,u,v,scale,'Linewidth',0.3);



wage_export = wage_export(:);
ttheta_export = ttheta_export(:);
cc = cc(:);
kk = kk(:);
nn = nn(:);
dlmwrite('../CUDA_VFI/wage_export.csv',wage_export,'precision',16);
dlmwrite('../CUDA_VFI/ttheta_export.csv',ttheta_export,'precision',16);
dlmwrite('../CUDA_VFI/cPEA_export.csv',cc,'precision',16);
dlmwrite('../CUDA_VFI/kPEA_export.csv',kk,'precision',16);
dlmwrite('../CUDA_VFI/nPEA_export.csv',nn,'precision',16);


save('PEA_Em.mat');
