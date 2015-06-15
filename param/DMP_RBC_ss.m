function moments=DMP_RBC_ss(solution,pars)


beta      = pars(1);
A        = pars(2);
alpha     = pars(3);
eta       = pars(4);
tau       = pars(5);
Qss       = pars(6);
x        = pars(7);
i_y_target= pars(8);
rr        = pars(9);
z_HM      = pars(10);
agg_theta = pars(11);
agg_jfr   = pars(12);
agg_u_target=pars(13);
 



kappa  = solution(1);
z      = solution(2);
gamma  = solution(3);
delta  = solution(4);
solution

 
theta_ss=agg_theta;


% Solving the steady state of the moel given a set of parameters

rss = 1/beta - 1 + delta;
k_n_ss = (rss/(A*alpha))^(1/(alpha-1));
 
v_ss=agg_theta*agg_u_target;
n_ss=1-agg_u_target;
 

% Defining q's and mu's
 
mu_ss=agg_jfr;
xi=mu_ss/agg_theta^eta;
q_ss  = xi*agg_theta^(eta-1);


% Solving for the rest of the model (see pdf notes)

k_ss = k_n_ss*n_ss;
inv_ss = delta*k_ss;
 

y_ss = A*(k_ss^alpha)*(n_ss^(1-alpha));
c_ss = y_ss - delta*k_ss - kappa*v_ss + z*(1-n_ss); % [Linxi:] add home production here


omega_ss = tau*A*(1-alpha)*k_n_ss^alpha + (1-tau)*(z + gamma*c_ss) + tau*kappa*theta_ss;
const1 = 1-(1-x)*beta;
Jn_ss = (A*(1-alpha)*k_n_ss^alpha - omega_ss)/const1;

 


% Replacement ratio as a function of the L wage
rep_ratio = z/omega_ss;
% The number that should match the Hall-Milgrom calibration
leisure_value = (gamma*c_ss)/omega_ss;
% investment output ratio
i_y_ss = inv_ss/(y_ss);



moments = [
        rep_ratio - rr;
        leisure_value - z_HM;
        i_y_ss - i_y_target;
        kappa -  (q_ss)*beta*(Jn_ss);
        % added this line to force z = 0
        ];
    moments
end
    

 

