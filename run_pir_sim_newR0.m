function savefile = run_pir_sim_newR0(eta_c,Dwg_factor,species,...
    dt,dayspersave,t_finish,...
    radiussave,radius_final,savefile,...
    center_ind1,center_ind2)

%This function is meant to run the two dimensional PIR model with spatially
%variable initial distribution of OPCs.
%For reproducibility, please create script to define parameters which calls
%this function. (i.e. use the script as an input file.)

% A good script to use: run_2dPIR_sim.m

%Input Parameters:
%eta_c = PDGF production rate (ng/cell/day) by transduced cells
%Dwg_factor = factor by which diffusion in white is faster than in gray
%species = mouse or human; to load correct atlas.
%dt = time step size (day)
%dayspersave = allows the option to write out the solution at various time
%      steps ("days" was in reference to when time steps were always 1 day)
%      if value of 0 is provided, will instead save according to tumor
%      radius increments as defined by additional input parameter, note
%      days is the default and will be used unless set to 0
%t_finish = if saving according to dayspersave, defines the final time step
%radiussave = analagous to dayspersave, will save simulation data according
%       to radius size instead of the time points. 
%radius_final = analagous to t_finish defines the radius at which to stop
%       the simulation
%savefile = string giving the name of the file in which simulation data
%       will be saved to
%center_ind{1,2} = indices corresponding to pixel location defining the
%       center of the tumor

% dayspersave=dayspersave/dt; % Since dt is in fraction of days for me - scm

tic

%Begin by defining all parameters and storing them in struct P:

  % ---------- PROGENITOR PARAMETERS ----------------------------------
  
      % Proliferation rate (1 / day)
        P.RhoCBar = log(2)/(18/24); %rho;  
        P.RhoRBar = log(2)/(18/24); %rho;
%         P.RhoCBar = log(2)/(48/24); %rho;  
%         P.RhoRBar = log(2)/(48/24); %rho;
        
      % Diffusivity (cm^2 / day)       
        P.D_c = (2/3)*8.7e-5; % D_c;
%         P.D_c = 1e-3; % D_c;
        P.D_r = P.D_c;
        P.Dwg = Dwg_factor;

      % Carrying capacity (cells / cc tissue)
        P.K = (2.3)*10^8;   
      % Computed assuming cell radius = 13.365 microns (ie, volume of 10^4 microns^3), per Peter Canoll.
      

  % ---------- PDGF PARAMETERS ----------------------------------

      % Diffusion coeff 
        P.D_p = 0.0005; % cm^2/day;

      % Secretion rate
        P.eta_c = eta_c; %10^(-4); %ng/mL/day; PDGF secretion rate. 

      % Consumption rate
        P.Qc = 10^(-5.15); %10^(-3.15) - with no tis2ecs factor... should be bigger, but then have stiffness?
        P.Qr = P.Qc;

      % Half-max Receptor Activation (Michaelis Menten)
        P.km = 30; % ng/mL; 

      % EC_50 - Half max for downstream effects
        P.EC50 = 10^0.5; %ng/mL; Pringle 1989.

      % Gamma - Scaling term for downstream phenomena affected by PDGF
        P.gamma = P.EC50/(P.km+P.EC50);

      % Autocrine PDGF bound to Infected (Transduced) Cells
        P.p_auto = 1; % ng/mL; 

      % Formerly "not dumb"
        P.alpha=1-(P.p_auto/(P.km+P.p_auto));

% ----------------------------------------------------------------  
  
% ---------- GRID PARAMETERS ----------------------------------
      % Time step (oughta make this dynamic someday, check (ahd 6-30-2011))
        P.dt = dt;
        
      % Grid spacing, set to be fixed to brain map spacing
      if strcmp(species,'mouse')
        P.h = .01; % 0.1 mm = 0.01 cm, 
      elseif strcmp(species,'human')
%         P.h=.1; % 1 mm = 0.1 cm
        P.h=0.025; %cm; for refined grid 
      end;
      
% ----------------------------------------------------------------  


% ---------- SIM PARAMETERS ----------------------------------
      % Tumor Origin (injection site) for this sim
        P.XYorigin = [center_ind1 center_ind2]; % start as input, but likely fixed
        
      % Slice Type to use for sims
        slice_type = 'coronal'; %do I want to pass this in instead? good to save in P though, for sure!
        P.slice_type = slice_type;
%         P.slice_num = slice_num;
        
      % Remember what species this is for:
        P.species = species;
      
      % Name to save simulation
        P.Simname = savefile;      
% ----------------------------------------------------------------  

fprintf('\nYou are running a simulation with the following parameters:\n')
P % PRINTS THE PARAMETERS
fprintf(['\nSaving here:\n' regexprep([pwd '/' savefile '.mat'],'\','\\\') '\n'])
fprintf('Simulating....')
% ---------------------------------------------------------------- 

% Procede to loading brain map, setting cell interfaces and initial conditions:
if strcmp(species,'human') % if want to use the human atlas
    load ./Brain.mat 
    grey=GrayMatter; 
    white=WhiteMatter;
    P.slice_num = 138; % picked this coronal slice because of ventricle spacing/orientation - matches mouse better; also interested
    % in slices 103 or 104 (center-ish).
    
    % Pick out the brain slice based on the slice number
    % and slice type.        
    if strfind(slice_type, 'sagittal')
      Slice = squeeze(GrayMatter(center_ind1,:,:) + Dwg_factor*WhiteMatter(center_ind1,:,:));
    elseif strfind(slice_type, 'coronal')
%       Slice = squeeze(GrayMatter(:,center_ind2,:) + Dwg_factor*WhiteMatter(:,center_ind2,:));
      Slice = squeeze(GrayMatter(:,P.slice_num,:) + Dwg_factor*WhiteMatter(:,P.slice_num,:));
    elseif strfind(slice_type, 'axial')
%       Slice = squeeze(GrayMatter(:,:,center_ind3) + Dwg_factor*WhiteMatter(:,:,center_ind3));
      Slice = squeeze(GrayMatter(:,:,center_ind2) + Dwg_factor*WhiteMatter(:,:,center_ind2));
    end
    
elseif strcmp(species,'mouse') % if want to use the mouse atlas.
    load 2D_mouse_brain.mat % mouse brain atlas on which to run PIR sims; we used http://brainatlas.mbi.ufl.edu
%     MouseBrain=grey(:,:)+white(:,:);
%     GwMouseBrain=grey(:,:)+Dwg_factor*white(:,:);
    Slice=grey(:,:)+Dwg_factor*white(:,:);
end;

%------ SET CELL INTERFACES ------
% 
%     % Pick out the brain slice based on the slice number 
%     % and slice type.        
%     if strfind(slice_type, 'axial')
%       Slice = squeeze(Dg_over_Deff*GrayMatter(:,:,center_ind3) + Dw_over_Deff*WhiteMatter(:,:,center_ind3));
%     elseif strfind(slice_type, 'sagittal')
%       Slice = squeeze(Dg_over_Deff*GrayMatter(center_ind1,:,:) + Dw_over_Deff*WhiteMatter(center_ind1,:,:));
%     elseif strfind(slice_type, 'coronal')
%       Slice = squeeze(Dg_over_Deff*GrayMatter(:,center_ind2,:) + Dw_over_Deff*WhiteMatter(:,center_ind2,:));
%     end
    
    % Measure the dimensions of the brain slice.
    n1 = size(Slice, 1);
    n2 = size(Slice, 2);
    
    if strcmp(species,'human')
        NewSlice=zeros(4*n1,4*n2);
        for i=1:n1
            for j=1:n2
                NewSlice(4*i-3:4*i,4*j-3:4*j)=Slice(i,j);
            end;
        end;
        clear Slice
        Slice=NewSlice;
        n1=size(NewSlice, 1);
        n2=size(NewSlice, 2);
    end;
    
    % Set up the cell interfaces. 
    [Dweight1, Dweight2, Origin1, Origin2, Terminus1, Terminus2, ind_nz,boundary] ...
          = set_interfaces(Slice);
    
%------------------------------------------------------------


%----------- INITIAL CONDITIONS -----------
  % Set up coordinate vectors for the spatial dimensions. 
  % only doing this to set the initial condition.
    x1 = linspace(0, (n1-1)*.1, n1);
    x2 = linspace(0, (n2-1)*.1, n2);

  % Now meshgrid it.  But I need to transpose the result, since I'm always
  % taking the first space dimension to correspond to the first index.
    [X1,X2] = meshgrid(x1,x2);
    X1 = X1';
    X2 = X2';

    if strcmp(species,'mouse')
        x1_center = x1(center_ind1);
        x2_center = x2(center_ind2);
    elseif strcmp(species,'human')
        x1_center = x1(4*center_ind1);
        x2_center = x2(4*center_ind2);
    end;
    
    ninitx = 0.03; %cm or 0.3 mm; radius of initial sphere affected by injection - NEED TO RELATE THIS TO X1,X2_CTR
    O2a = 250/((4/3)*pi*ninitx^3); %cells/mL = cells/cc = cells/cm^3 
    % code for human sims is breaking with O2a... trying a lower amount to
    % see if it goes okay.. or if it's the relative c0 vs r0.
%     O2a=1e3;
    
    a=X1-x1_center;
    b=X2-x2_center;
    
    c0=zeros(size(Slice));
    r0=zeros(size(Slice));
    
    if strcmp(species,'mouse')
        for i=1:length(Slice(:,1))
            for j=1:length(Slice(1,:))
                % Check if inside brain!:
                if grey(i,j)==1 %|| white(i,j)==1
                    % Then check to see if inside sphere of injection:
                    if abs(a(i,j))<ninitx && abs(b(i,j))<ninitx
                        c0(i,j) = 0.99*O2a; % cells/mL
                        r0(i,j) = O2a-c0(i,j);   % cells/mL
                    else
                        c0(i,j) = 0;   % cells/mL
                        r0(i,j) = O2a; % cells/mL
                    end;
                elseif white(i,j)==1
                    if abs(a(i,j))<ninitx && abs(b(i,j))<ninitx
                        c0(i,j) = 0.99*(8/3)*O2a; % cells/mL
                        r0(i,j) = (8/3)*O2a-c0(i,j);   % cells/mL
                    else
                        c0(i,j) = 0;   % cells/mL
                        r0(i,j) = (8/3)*O2a; % cells/mL
                    end;
                else % not in the brain; these need to be zero! 
                    c0(i,j) = 0;   % cells/mL
                    r0(i,j) = 0; % cells/mL
                end;
            end;
        end;
    elseif strcmp(species,'human')
        % Initially "infected" (transduced) cells set to small blip with
        % max = O2a (the baseline transduceable/recruitable OPCs)
%         [indx,indy]=find(sqrt((X1-x1_center).^2 + (X2-x2_center).^2)<ninitx*2);
%         c0(indx,indy) = .99*O2a;
%         c0 = .99*O2a * exp(-50*((X1-x1_center).^2 + (X2-x2_center).^2));
        c0 = .99* O2a * exp(-100*((X1-x1_center).^2 + (X2-x2_center).^2));
%         c0 = 1e5 * exp(-100*((X1-x1_center).^2 + (X2-x2_center).^2));
        c0 = c0 .* (Slice > 0);
        r0 = O2a - c0;
        r0 = r0 .* (Slice > 0);
    end;
    
    IC.c0 = c0(:);
    IC.r0=r0(:);
    IC.p0 = zeros(size(IC.c0)); %could also try making this nonzero..

    clear x1 x2 X1 X2 

%------------------------------------------


%-------------- CALL THE NUMERICS FUNCTION --------------
    Simdata = pir_numerics_2d(Dweight1, Dweight2, ...
                                Origin1, Origin2, ...
                                Terminus1, Terminus2, ...
                                n1, n2, ind_nz,boundary,...
                                IC, species,...
                                dayspersave, t_finish, radiussave, radius_final, P);
%--------------------------------------------------------


%---- SAVE THE PARAMETERS AND THE BRAIN SLICE -----------
    Simdata.Parameters = P;
    Simdata.Slice = Slice;
    Simdata.Simname = savefile;
    save(savefile, 'Simdata');
%--------------------------------------------------------
toc
