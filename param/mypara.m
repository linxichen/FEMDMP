% Follows Nir's Code
bbeta = 0.999; % discount rate
% ddelta = 0.025; % depreciation
% kkappa = 0.213; % cost of vacancy
aalpha = 0.33; % share of capital
xxi = 1.390; % matching efficiency
% ggamma = 0.399; % work disutility
eeta = 0.4; % share of vacancy in matching function
ttau = 1-eeta; % worker's bargaining power
x = 0.0081; % exo seperation rate
agg_u_target=x/(x+0.139); % we target this unemployment rate
Qss = bbeta; % steady state SDF

% Markov Process
ssigma = 0.03;
rrho = 0.95;
Abar = 1;

pars(1) = bbeta;
pars(2) = Abar;
pars(3) = aalpha;
pars(4) = eeta;
pars(5) = ttau;
pars(6) = Qss;
pars(7) = x;

% The set of momemts to be mathced
rr         = 0.4;     %Replacement ratio
z_HM       = 0.4;     % match the leisure value as a fraction of avg. wage as in Hall-Milgrom
agg_ttheta  = 1;     % Average tightness ratio
agg_jfr    = 0.139;   %avergae job finding rate
i_y_target = 0.2;    % This is now the target of inv/output ratio

% z = z_HM; % unemployment benefit

% Parameters that are targets
pars(8)   = i_y_target   ; 
pars(9)   = rr;
pars(10)   = z_HM;
pars(11)   = agg_ttheta;
pars(12)   = agg_jfr;
pars(13) = agg_u_target;

guess1 = [0.3 1.2   0.5 0.002]

options=optimset('MaxIter',2000,'MaxFunEvals',20000,'TolFun',1e-22,'Diagnostics','on', 'Display','on');
%[x,fval,flag]=fsolve('rbc_nir_ss_function',guess,options,rbc_pars);
[solution,fval,flag]=fsolve('DMP_RBC_ss',guess1,options,pars); % The objective function can't be found, switching to solve_rbc_ss_function_alt
% [solution,fval,flag]=fsolve('solve_rbc_ss_function_alt',guess1,options,pars); 
solution;

kkappa  = solution(1);
z      = solution(2);
ggamma  = solution(3);
ddelta  = solution(4);

r_ss = 1/bbeta - 1 + ddelta;
k_n_ss = (r_ss/(Abar*aalpha))^(1/(aalpha-1));
v_ss=agg_ttheta*agg_u_target;
n_ss=1-agg_u_target;
ttheta_ss=agg_ttheta;
mmu_ss=agg_jfr;
xxi = mmu_ss/agg_ttheta^eeta;
q_ss  = xxi*agg_ttheta^(eeta-1);
k_ss = k_n_ss*n_ss;
inv_ss = ddelta*k_ss;
y_ss = Abar*(k_ss^aalpha)*(n_ss^(1-aalpha));
c_ss = y_ss - ddelta*k_ss - kkappa*v_ss + z*(1-n_ss);
oomega_ss = ttau*Abar*(1-aalpha)*k_n_ss^aalpha + (1-ttau)*(z + ggamma*c_ss) + ttau*kkappa*ttheta_ss;
const1 = 1-(1-x)*bbeta;
Jn_ss = (Abar*(1-aalpha)*k_n_ss^aalpha - oomega_ss)/const1;
rep_ratio = z/oomega_ss;
leisure_value = (ggamma*c_ss)/oomega_ss;
i_y_ss = inv_ss/(y_ss);
lp_agg_ss=y_ss/n_ss;
moments = [
        rep_ratio - rr;
        leisure_value - z_HM;
        i_y_ss - i_y_target;
        kkappa -  (q_ss)*bbeta*(Jn_ss);
        ];
moments;

 
A_ss=Abar;u_ss=agg_u_target;SDF_ss=bbeta;invest_ss=i_y_ss;

ssigma_A = ssigma;
rrho_A = rrho;
save DMP_RBC_PARAS.mat ...
    lp_agg_ss rrho ssigma ssigma_A rrho_A z ttau ggamma kkappa bbeta aalpha xxi eeta x ddelta y_ss A_ss Abar k_ss n_ss mmu_ss q_ss ttheta_ss  u_ss SDF_ss oomega_ss r_ss Jn_ss c_ss invest_ss v_ss;


