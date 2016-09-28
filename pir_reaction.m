function [R, J] = pir_reaction(Q, ...
                                 rho_c,rho_r, K, ...
                                 km, pa, EC50, alpha, ...
                                 eta_c, q_c, q_r)
                             
% Given Q=[c; r; p], this function calculates the reaction
% R=R(Q), and the Jacobian J=R'(Q). 

    c = Q(1);
    r = Q(2);
    p = Q(3);

%--------FUNCTIONS USED IN R AND J-----------

% I like the neatness of useing

% Compute Michaelis-Menten terms:
    beta = @(s) s /(km + s);
%   cell surface receptor saturation ratio, paracrine: beta(p) = p / (km + p); % formerly delta(p)
%   cell surface receptor saturation, paracrine and autocrine: beta(p+pa); % beta(p+pa) = (p + pa) / (km + p + pa);
    beta_p=beta(p);
    beta_ppa=beta(p+pa);
    beta_EC50=beta(EC50);
% Compute Downstream Reponse terms:
%     response = @(s) beta(s) / (beta(EC50) + beta(s)); % instead of gamma, using beta(EC50)
    response_p = beta_p/(beta_EC50 + beta_p);
    response_ppa = beta_ppa/(beta_EC50 + beta_ppa);
% Compute Derivatives of Michaelis_Menten and Response terms:
%     dbeta_dp = @(s) km/(km+s)^2; % dp_dppa  = dbeta_dp(p+pa); % = km / (km + p + pa)^2; % dp_dp  = dbeta_dp(p); % = km / (km + p)^2;
%     dresponse_dp = @(s) (beta(EC50)*km)/((km+s)^2*(beta(EC50)+beta(s))^2); % dp_gp  = (gamma * km) / ( (km + p)^2 * (gamma + beta_p)^2 ); % dp_gppa = (gamma * km) / ( (km+(p+pa))^2*(gamma+beta_ppa)^2);
    dbeta_p_dp = km/(km+p)^2;
    dresponse_p_dp = (beta_EC50*km)/((km+p)^2*(beta_EC50+beta_p)^2);
    dresponse_ppa_dp = (beta_EC50*km)/((km+p+pa)^2*(beta_EC50+beta_ppa)^2);
% Logistic Growth:
    cap = (1 - (c+r)/K);
    
%--------REACTION-----------      
    
% The reaction vector, R(Q)
    R = zeros(3,1);

%     Written in functional form - but makes code take longer (since
%     computed on each iteration..)
%     R(1) = rho_c * response(p+pa) * c * cap;    
%     R(2) = rho_r * response(p) * r * cap;
%     R(3) = eta_c * c - q_r * beta(p) * r - q_c * beta(p) * alpha * c; % was this way but think it's wrong: - q_c * beta(p+pa) * alpha * c (SCM, 27-May-2016)

    R(1) = rho_c * response_ppa * c * cap;    
    R(2) = rho_r * response_p * r * cap;
    R(3) = eta_c * c - q_r * beta_p * r - q_c * beta_p * alpha * c; % was this way but think it's wrong: - q_c * beta(p+pa) * alpha * c (SCM, 27-May-2016)

%------JACOBIAN------------

% The Jacobian matrix, J = R'(Q)
    J = zeros(3,3);
    
%     Written in functional form - but makes code take longer (since
%     computed on each iteration..)
% % First row
%     J(1,1) =  rho_c * response(p+pa) * (1 - (2*c +r)/K); % note: because c+r=K *at equilibruim,* this becomes = 1-r/K - (c+r =K)/K = -r/K
%     J(1,2) = -rho_c * response(p+pa) * c/K;
%     J(1,3) =  rho_c * c * cap * dresponse_dp(p+pa);
% % Second row
%     J(2,1) = -rho_r * response(p) * r/K;
%     J(2,2) =  rho_r * response(p) * (1- (2*r + c)/K);
%     J(2,3) =  rho_r * r * cap * dresponse_dp(p);
% % Third row
%     J(3,1) =  eta_c - q_c * beta(p) * alpha; % - q_c * beta_ppa * alpha; was this way, but think wrong - see comment on R(3)
%     J(3,2) = -q_r * beta(p);
%     J(3,3) = -q_r * r * dbeta_dp(p) - q_c * c * alpha * dbeta_dp(p); % - q_c * c * alpha * dbeta_dp(p+pa); consistent with note above, was this way, but think it's wrong

% First row
    J(1,1) =  rho_c * response_ppa * (1 - (2*c +r)/K); % note: because c+r=K *at equilibruim,* this becomes = 1-r/K - (c+r =K)/K = -r/K
    J(1,2) = -rho_c * response_ppa * c/K;
    J(1,3) =  rho_c * c * cap * dresponse_ppa_dp;
% Second row
    J(2,1) = -rho_r * response_p * r/K;
    J(2,2) =  rho_r * response_p * (1- (2*r + c)/K);
    J(2,3) =  rho_r * r * cap * dresponse_p_dp;
% Third row
    J(3,1) =  eta_c - q_c * beta_p * alpha; % - q_c * beta_ppa * alpha; was this way, but think wrong - see comment on R(3)
    J(3,2) = -q_r * beta_p;
    J(3,3) = -q_r * r * dbeta_p_dp - q_c * c * alpha * dbeta_p_dp; % - q_c * c * alpha * dbeta_dp(p+pa); consistent with note above, was this way, but think it's wrong
