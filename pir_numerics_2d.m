function Simdata = pir_numerics_2d(Dweight1, Dweight2, ...
    Origin1, Origin2, ...
    Terminus1, Terminus2, ...
    n1, n2, ind_nz, boundary, ...
    IC, species,...
    dayspersave, t_finish, radiussave, radius_final, P)


%Code to do primary numerics for solve

%Input parameters
%Dweight{1,2} = weights to affect diffusion coefficients differently for
%       gray and white matter
%Origin{1,2} = Holds appropriate indices for the origin of cells in x and y
%       direction respectively
%Terminus{1,2} = Holds appropriate indices for the terminus of cells in x
%       and y direction respecitvely
%[n1,n2] = dimensions of slice (i.e. number of pixels in x and y direction
%       of the currently used brain slice)
%IC = vector holding initial conditions for all variables (c,r,p)
%dayspersave = an indicator for the interval between saves according to
%       timestep. If this is set to 0, will instead base termination off of 
%       radius values
%t_finish = if dayspersave ~=0, will indicate the final time step for
%       simulation
%radiussave = if dayspersave==0, will be used to indicate how frequently
%       simulation will be saved based on tumor radius size (note
%       dayspersave is default)
%radius_final = analagous to t_finish, will determine final size of tumor
%       for simulation termination 
%P = struct holding all parameters relevant for simulation 


%-------------- INITIALIZATION, OR SOMETHING ---------------
% defining radii savePoints
if(dayspersave==0)
  savePoints = radiussave:radiussave:radius_final;
else
  savePoints = .1:.1:2;
end

% Total number of grid cells
N = n1*n2;

% Initialize arrays for saving data
if(dayspersave~=0)
 %values if saving by predefined timestep
 num2save = round(t_finish / dayspersave);
else
 %values if saving by tumor growth
 num2save = size(savePoints,2);
end

%Simdata.t = zeros(num2save, 1); 

%Simdata.c = zeros(num2save, N);
%Simdata.r = zeros(num2save, N);
%Simdata.p = zeros(num2save, N);

% Set step sizes:
dt = P.dt; %.5; % This is an input in run_pir_sim though. Half day time step in 1d PIR; is that a good choice for 2d case?
dx = P.h; %.01; % 1/375 in 1d.

%------PARAMETERS----------------------------
% below are the values passed in from run_pir_sim_newR0.m;
% the values given are the same in both 1D and 2D PIR models
% Dp = 0.0005; % cm^2/day; 
% k = (2.3)*10^8; % cells/cm^3; carrying capacity. Computed assuming cell radius = 
%     % 13.365 microns (ie, volume of 10^4 microns^3), per Peter Canoll.
% lambda = 10^(-8); % cm^3/cell; inverse of carrying capacity.
% epsilon_p = 0.175; %days^-1; 4.2 hr half-life; Lynch, et al.
% p_auto = 1; % ng/mL; autocrine PDGF bound to infected cells.
% km = 30; % ng/mL; half-max receptor activation;
% EC50= 10^0.5; %ng/mL; half max PDGF for downstream effects;
% gamma = EC50/(km+EC50); %set gamma = beta(EC50)
% alpha=1-(p_auto/(km+p_auto)); = 1 - beta(p_auto)
% DrBar = (2/3)*8.7e-5; %cm2/day;
% DcBar = DrBar;
% RhoCBar = log(2)/(18/24); % 1/day;
% RhoRBar = RhoCBar;
% nchi = 10^(-8);
% Qc=10^(-5.15);
% Qr=Qc;
% eta_c = 10^(-5); % ng/mL/day (range: 10^-(4 to 6)
D_p = P.D_p;
D_c = P.D_c;
D_r = P.D_r;
rho_c = P.RhoCBar; % log(2)/(18/24); % 1/day; 
rho_r = P.RhoRBar;
K = P.K; % carrying capacity, cells/mL
eta_c = P.eta_c; %ng/mL/day; PDGF secretion rate. % INPUT?
q_c = P.Qc; 
q_r = P.Qr;
km = P.km; % ng/mL; 
% gamma = P.gamma;
EC50 = P.EC50;
pa = P.p_auto; % ng/mL; 
alpha = P.alpha;

% Initial conditions:
c0 = IC.c0;
r0 = IC.r0;
p0 = IC.p0;

% these correspond to:
% O2a= 250/((4/3)*pi*ninitx^3); %cells/mL = cells/cc = cells/cm^3 
% c0 = 0.99*O2a; % cells/mL - ie, we're assuming that 50% of the cells are infected
% r0 = O2a-c0; % cells/cm^3
% p0 = 0; 
% n0 = (k*.35)-O2a; % cells/cm^3 - have left out of 2d PIR


%----- SET THE PDGF DIFFUSION OPERATOR -------------------
dfsn_pdgf = @(p) diffusion_op(p, D_p, dx, ...
    Dweight1, Dweight2, ...
    Origin1, Origin2, ...
    Terminus1, Terminus2);

%----- INITIALIZE THE UNKNOWNS ------------------------
c = c0;
r = r0;
p = p0;

%----- BEGIN TIME LOOP --------------------------------
t_dim = 0;
t_step = 0;
isave = 0;
rT1 = 0; %initial T1 radius

%for tstep = 1:numsteps
while(isave < num2save)
    t_dim = t_dim + P.dt;
    t_step = t_step + 1;

    %------ BEGIN DIFFUSION STEP -------

    % Compute the available space, and PDGF gradient, and interpolate the values to the interfaces.
%     S = 1 - (c+r)/K;
%     S1 = (S(Origin1) + S(Terminus1)) / 2;
%     S2 = (S(Origin2) + S(Terminus2)) / 2;
                                
    SR = (p./(km+p))./(EC50/(km+EC50)+p./(km+p)).*(1 - (c+r)/K);
    SR1 = (SR(Origin1)+SR(Terminus1))/2;
    SR2 = (SR(Origin2)+SR(Terminus2))/2;
                                
    SC = ((p+pa)./(km+p+pa))./(EC50/(km+EC50)+(p+pa)./(km+p+pa)).*(1 - (c+r)/K);
    SC1 = (SC(Origin1)+SC(Terminus1))/2;
    SC2 = (SC(Origin2)+SC(Terminus2))/2;

    % Set up the transduced cell diffusion operator.
    dfsn_c = @(u,D) diffusion_op(u, D, dx, ...
        SC1.*Dweight1, SC2.*Dweight2, ...
        Origin1, Origin2, ...
        Terminus1, Terminus2);
    
    % Set up the recruited cell diffusion operator.
    dfsn_r = @(u,D) diffusion_op(u, D, dx, ...
        SR1.*Dweight1, SR2.*Dweight2, ...
        Origin1, Origin2, ...
        Terminus1, Terminus2);

    % The system I'm solving with conjugate gradients is
    %     (I - dt*D)*du = dt*D(u),
    % where dfsn represents the diffusion operator.

    % A. Transduced cell diffusion
    % Set up the left-hand side
    b = dt * dfsn_c(c, D_c);
    % Solve the system with conjugate gradients
    [dc,flag] = pcg(@(c) c - dt*dfsn_c(c,D_c), b);
    % Update c
    c = c + dc;
    c = max(c,0);

    % B. Recruited cell diffusion
    % Set up the left-hand side
    b = dt * dfsn_r(r, D_r);
    % Solve the system with conjugate gradients
    [dr,flag] = pcg(@(r) r - dt*dfsn_r(r,D_r), b);
    % Update r
    r = r + dr;
    r = max(r,0);

    % C. TAF diffusion
    % Set up the left-hand side
    b = dt * dfsn_pdgf(p);
    % Solve the system with conjugate gradients
    [dp,flag] = pcg(@(p) p - dt*dfsn_pdgf(p), b);
    % Update p
    p = p + dp;

    
    %------ BEGIN REACTION STEP -------- 
    
    %------ SET THE HANDLE FOR THE REACTION FUNCTION ------
    
%     %Check for treatment and modify parameters as needed
%     if(P.Rx.treatment & r >= P.Rx.T1RxStart)
%          rxn = @(Q) pir_reaction(Q, ...
%             rho, K, ...
%             beta_treat, gamma_treat, alpha_h, ...
%             delta_c, delta_h, K_M_treat, q, lambda_a, omega, ...
%             mu_v, ...
%             alpha_n,...
%             Ktrans_low,Ktrans_hi, 0.1*Ktrans_hi);
%             %Kt,0.001*Ktrans_hi);
%     else
         rxn = @(Q) pir_reaction(Q, ...
                                 rho_c,rho_r, K, ...
                                 km, pa, EC50, alpha, ...
                                 eta_c, q_c, q_r);
%     end
    %------------------------------------------------------

    % The reaction time step is done using the Runge-Kutta TR-BDF2 method. 
    for i = 1:N
        Q = [c(i); r(i); p(i)];
        if norm([c(i); r(i); p(i)]) ~= 0
            
            % First step of TR-BDF2 is to solve
            %     Qs = Q(n) + dt/4 * (R(Q(n)) + R(Qs)),
            % where Qn is Q at the nth time step.  This is
            % done using Newton's method.
            
            Qn = Q;
            [Rn,Jn] = rxn(Qn);
            M = (Qn + dt/4*Rn);
            Qs = Qn;
%             if(i==18289)
%                 i
%             end
            for k=1:1000
                [Rs, Js] = rxn(Qs);
                dQs = (eye(3,3) - dt/4 * Js) \ -(Qs - dt/4*Rs - M );
                Qs = Qs + dQs;
                if norm(dQs) < 1e-4
                    break
                end
            end  % end Newton loop for step 1
            
            % Now we have Qs, so we can start the second stage of
            % TR-BDF2.  This is
            %     Q(n+1) = 1/3 * ( 4*Qs - Q(n) + dt*R(Q(n+1)) ).
            % This is another Newton's method, and we start with
            % initial guess Q(n+1) = Q(n).
            
            M = (4*Qs - Qn)/3;
            Qnp1 = Qs;
            for k=1:1000
                [Rnp1, Jnp1] = rxn(Qnp1);
                dQnp1 = (eye(3,3) - dt/3 * Jnp1) \ -(Qnp1 - dt/3*Rnp1 - M);
                if rcond(eye(3,3) - dt/3 * Jnp1)<10e-16
                    fprintf(['\n badly scaled at time step = ' num2str(t_step), ', i = ', num2str(i), ', dQnp1 = ', num2str(dQnp1(1)),', ', num2str(dQnp1(2)),', ', num2str(dQnp1(3)), '\n']);
                    break
                end;
                Qnp1 = Qnp1 + dQnp1;
                if norm(dQnp1) < 1e-4
                    break
                end
            end  % end Newton loop for step 2

            % Update c, r, and p.
            c(i) = Qnp1(1);
            r(i) = Qnp1(2);
            p(i) = Qnp1(3);

        end  % end norm(Q) if statement

    end  % end i loop
    %          outtime = toc;
    %          fprintf([num2str(outtime) ' seconds for reaction step.\n']);

    %------ END REACTION STEP ----------

    
    isave = isave+1;

    Simdata.t(isave) = t_dim; 
    Simdata.c(isave,:) = c;
    Simdata.r(isave,:) = r;
    Simdata.p(isave,:) = p;

    
    [rT1, rT2] = T2radius(Simdata,isave,K,dx); % replaced "r" or "rad" with rT1
    
    % save the data according to desired method (i.e. days or radius)
    
    %First check if want to save by days:
    if(dayspersave~=0 && mod(t_step,dayspersave)==0)
        Simdata.T1radius(isave,:) = rT1;
        Simdata.T2radius(isave,:) = rT2; % T2radius(Simdata,isave,K,dx);
        fprintf(['\nFinished day ' num2str(t_dim) ' with T1 radius = ' num2str(rT1) ' cm and T2 radius = ', num2str(rT2) ' cm.\n']); % replaced num2str(T2radius(Simdata,isave,K)) with num2str(rT1)
    elseif(dayspersave==0 && (t_dim == P.dt || rT1>savePoints(1)))
        savePoints = savePoints(2:end);
        [Simdata.T1radius(isave,:), Simdata.T2radius(isave,:)] = T2radius(Simdata,isave,K,dx);
        fprintf(['\nFinished day ' num2str(t_dim) ' with T1 radius = ' num2str(rT1) ' cm and T2 radius = ', num2str(rT2) ' cm.\n']); % replaced num2str(T2radius(Simdata,isave,K)) with num2str(rT1)
    else
        isave = isave-1;
    end

%     % Save to a file every 200 days
%     if (mod(t_dim, 200) < mod(t_dim - P.dt, 200)) && (isave < num2save)
%         save TempSave Simdata
%         clock;
%     end
    
    %Definitely kill if tumor gets above a certain radius (1 cm = entire rat brain)
    if strcmp(species,'mouse') && (rT2 > 1)
        break;
    elseif strcmp(species,'human') && (rT2 > 4.5)
        break;
    end
    
    %Also kill if we get to a ridiculous time point (100 days? 1 year? - rats die after 20 days with tumors, ~1000 if well) 
    if strcmp(species,'mouse') && (t_dim > 1000)
        break;
    elseif strcmp(species,'human') && (t_dim > 1000)
        break;
    end

end
%----- END TIME LOOP ----------------------------------
