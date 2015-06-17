//----------------------------------------------------------------
// 1. Declare variables
//----------------------------------------------------------------
var 

// endogenous variables
c
k           // variable 1
n           // variable 2
v           // variable 3
ttheta      // variable 4
y           // variable 5
mh          // 
mf
q
mmu

// exogenous variables
A          // variable 14
;




//----------------------------------------------------------------
// 2. Exogenous shocks
//----------------------------------------------------------------

varexo eps;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

bbeta ttau aalpha ddelta Abar kkappa ssigma rrho eeta xxi z ggamma x

;


//----------------------------------------------------------------
// 4. Calibration 
//----------------------------------------------------------------
load DMP_RBC_PARAS.mat;
for i=1:length(M_.params)
    deep_parameter_name = M_.param_names(i,:);
    eval(['M_.params(i)  = ' deep_parameter_name ' ;'])
end

//----------------------------------------------------------------
// 6. Model
//----------------------------------------------------------------

model;

  // 1. Definiton of log output
  y = A + aalpha*k(-1) + (1-aalpha)*n(-1);

  // 2. Definition of q
  exp(q) = xxi*(exp(ttheta))^(eeta-1);

  // 3. Definition of mmu
  exp(mmu) = exp(q)*exp(ttheta);

  // 4. LOM for TFP
  A = rrho*A(-1) + ssigma*eps;

  // 5. Def mh
  exp(mh) = 1/exp(c)*(1-ddelta+aalpha*exp(y)/exp(k(-1)));

  // 6. Def mf
  exp(mf) = 1/exp(c)*( (1-ttau)*((1-aalpha)*exp(y)/exp(n(-1))-z-ggamma*exp(c)) + (1-x)*kkappa/exp(q) - ttau*kkappa*exp(ttheta) );  

  // 7. HH Euler 
  1/exp(c) = bbeta/exp(c(+1))*(1-ddelta+aalpha*exp(y(+1))/exp(k));
  
  // 8. Firm Euler
  kkappa/(exp(c)*exp(q)) = bbeta/exp(c(+1))*( (1-ttau)*((1-aalpha)*exp(y(+1))/exp(n)-z-ggamma*exp(c(+1))) + (1-x)*kkappa/exp(q(+1)) - ttau*kkappa*exp(ttheta(+1)) );  

  // 9. Resource
  exp(y) + z*(1-exp(n(-1))) = exp(k) - (1-ddelta)*exp(k(-1)) + exp(c) + kkappa*exp(v);

  // 10. Defintion of tightness
  exp(ttheta) = exp(v)/(1-exp(n(-1)));

  // 11. Employment LOM
  exp(n) = (1-x)*exp(n(-1)) + exp(mmu)*(1-exp(n(-1)));

end;

//----------------------------------------------------------------
// 7. Computation
//----------------------------------------------------------------
initval;
k = log(k_ss); 
n = log(n_ss);     
c = log(c_ss);
ttheta = log(ttheta_ss);
v = log(v_ss);
y = log(y_ss); 
mh = log(1/2.3/0.98);
mf = log(0.21/2.3/0.98/(1.355*1.3^(0.28-1)));
q = log(q_ss);
mmu =  log(mmu_ss);
eps = 0;
A = log(Abar);
end;


steady;
check;

shocks;
var eps;
stderr 1;
end;

stoch_simul(order = 1,periods=250000,irf=200); % compute polices up to 1st order
save