var pih ${\pi_h}$ (long_name='Domestic inflation')
    x $x$ (long_name='Output gap')
    y $y$ (long_name='Output')
    ynat ${\bar y}$ (long_name='Natural output')
    rnat ${\bar r}$ (long_name='Natural interest rate')
    r $r$ (long_name='Nominal interest rate')
    s $s$ (long_name='Terms of trade')
    pi ${\pi}$ (long_name='CPI inflation')
    p $p$ (long_name='CPI level')
    ph ${p_h}$ (long_name='Domestic price level')
    e $e$ (long_name='Exchange rate')
    ystar ${y^*}$ (long_name='World output')
    pistar ${\pi^{*}}$ (long_name='World inflation')
    n ${n}$ (long_name='Employment')
    nx ${nx}$ (long_name='Net Exports')
    real_wage ${w-p}$ (long_name='Real Wage')
    a $a$ (long_name='Risk aversion')
    c $c$ (long_name='Domestic consumption')
    deprec_rate $\Delta e_t$ (long_name='Nominal depr. rate')
	de
    ;

varexo eps_star ${\varepsilon^{*}}$ (long_name='World output shock')
    eps_a ${\varepsilon^{a}}$ (long_name='World output shock');

parameters sigma $\sigma$ (long_name='risk aversion')
    eta $\eta$ (long_name='Substitution home foreign')
    gamma $\gamma$ (long_name='Substitution between foreign')
    phi $\varphi$ (long_name='Inverse Frisch elasticity')
    epsilon $\varepsilon$ (long_name='Elasticit of substitution')
    theta $\theta$ (long_name='Calvo parameter')
    beta $\beta$ (long_name='discount factor')
    alpha $\alpha$ (long_name='openness')
    phi_pi $\phi_\pi$ (long_name='Feedback Taylor rule inflation')
    rhoa $\rho_a$ (long_name='autocorrelation TFP')
    rhoy $\rho_y$ (long_name='autocorrelation foreign output')
    phi_x
    phi_e
    rho_i
;

% set deep parameters
sigma = 1;
eta = 1 ;
gamma = 1;
phi = 3;
epsilon = 6;
theta = 0.75;
beta  = 0.96;
alpha = 0.62;
phi_pi = 1.5;
phi_x= 0.5;
phi_e= 0;
rhoa = 0.87;                                                                
rhoy = 0.87;
rho_i = 0.99; 

model(linear);
//define parameter dependencies
//steady state real interest rate, defined below equation (11)
#rho  = beta^(-1)-1;
//defined below equation (27)
#omega = sigma*gamma+(1-alpha)*(sigma*eta-1);
//defined below equation (29)
#sigma_a =sigma/((1-alpha)+alpha*omega);

#Theta=(sigma*gamma-1)+(1-alpha)*(sigma*eta-1);
//defined below equation (32)
#lambda = (1-(beta*theta))*(1-theta)/theta;
//defined below equation (36)
#kappa_a =lambda*(sigma_a+phi);
//defined below equation (35)
#Gamma = (1+phi)/(sigma_a+phi);
#Psi = -Theta*sigma_a/(sigma_a+phi);

//1. Equation (37), IS Curve
x    = x(+1) - sigma_a^(-1)*(r - pih(+1) - rnat) ;                              
//2. Equation (36), Philips Curve
pih  = beta * pih(+1)+ kappa_a*x;                                                
//3. Equation below (37)
rnat = -sigma_a*Gamma*(1-rhoa)*a + alpha*sigma_a*(Theta+Psi)*(ystar(+1)-ystar);
//4. Equation (35), definition natural level of output
ynat = Gamma*a + alpha*Psi*ystar;                                                 
//5. Equation above (35), definition output gap
x    = y - ynat;                                                               
//6. Equation (29)
y = ystar + sigma_a^(-1)*s;
//7. Equation (14)
pi   = pih + alpha*(s-s(-1));
//8. Equation 15 (first difference)
s    = s(-1) + e - e(-1) + pistar - pih;
//9. Constant world inflation, see p.724 (Given constant world prices) 
pistar = 0;
//10. Equation (22), employment
y = a + n;
//11. Equation (31), net exports
nx = alpha*(omega/sigma-1)*s;
//12. Equation (27), defines consumption
y = c+alpha*omega/sigma*s;
//13. Above equation (11), defines real wage
real_wage = sigma*c+phi*n;

//10-12. Equations on p. 723, stochastic processes
a    = rhoa*a(-1) + eps_a;
ystar= rhoy*ystar(-1) + eps_star;

// interest rate rule
r = rho_i*r(-1)+(1-rho_i)*(phi_pi*pih + phi_x*x + phi_e*de);

  
//definition consumer price level
pi   = p - p(-1);
//definition domestic price level
pih  = ph - ph(-1);
//definition nominal depreciation rate of exchange rate
deprec_rate=e-e(-1);
de=e-e(-6);
end;

shocks;
var eps_a = 0.018;
var eps_star = 0.02;
end;

optim_weights;
pih -1*((1-alpha)/2)*(epsilon/((1-(beta*theta))*(1-theta)/theta));
x -1*((1-alpha)/2)*(1+phi);
end;

osr_params phi_pi phi_x phi_e;

osr_params_bounds;
phi_pi ,0 ,2;
phi_x ,-0.5 ,1;
phi_e ,-0.5 ,0.5;
end;
osr(opt_algo=1) pih x de;
