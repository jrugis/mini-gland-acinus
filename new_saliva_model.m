% This is my program to solve a 3d reaction-diffusion problem in an
% irregular domain using a tetrahedral mesh.

% The program here uses linear basis functions and follows closely the
% book by Larson and Bengzon, with bits from Gockenbach also.

% Modified in 2021 to use the new saliva secretion model and incorporate
% fluid flow. The secretion bits are contained in the file secretion.m,
% which is called by ode15s.

% First, run the program make_cell_mesh.m which takes the mesh file from John, changes a few names, and precalculates the 
% tet volumes and triangle areas.

% Results are viewed by running plot_results.m.

% A few important things to note. 1. We do the IPR fluxes as boundary
% fluxes on the apical membrane. 2. The parameters only work for cell 4. To do
% another cell you first have to adjust the parameters to make sure you get
% the correct cell volume at rest. There are many ways to do this, but here
% I fiddle with KCa densities to get the proper resting volume. There may be better
% ways.


clear all
close all
clc

load('my_acinar_mesh.mat');             % Loads the acinar mesh.
load('par.mat');                        % Load the parameters


par.VPLC = 0.0025; 
par.spatialchoose = 1;                  % 1 for spatially distributed IPR.
par.apicalKCadensity = 3;               % 1 for equal density of KCa on apical and basal. Can run from 0 to bigger than 1.

delt = 0.05; tend = 50; numt = floor(tend/delt); tim = (0:delt:tend);

[mass,stiff] = makematrices(p,np,tets,ntets);                % Get the mass and stiffness matrices
Amat = sparse(mass + delt*stiff);   % Amat contains only entries for c and ip, not for h, which doesn't diffuse 


% This is where we do the actual solve. We use a backward Euler, and every
% time step we need to call rhs to get the new load vector (which includes
% the c and ip3 dependence of the rhs).

%   First np components of u    -   c
%   Second np components of u   -   ip3
%   Third np components of u    -   h

% options and mass matrix for ode15s (the fluid secretion solve)
M = eye(11);
M(10,10)=0; M(11,11)=0;   % 2 DAEs for Va and Vb
options = odeset('Mass',M,'RelTol', 1e-11, 'AbsTol', 1e-11);                    

[IC,xnew] = initialise(np);         % xnew is the initialisation for the PDE, IC is the initial condition for the secretion ODE bit
sol(:,1) = xnew;                    % Initialise and store in solution matrix. calcium variables in sol. secretion variables in SSsol.
SSsol(:,1) = IC;

for i=1:numt

    % First step forward the secretion model by delt, using ode15s. c, ip
    % and h are not stepped forward in this step.

    % You need to calculate the Cl and K permeabilities each time step, as
    % they change when the calcium changes each time step.
    
    [PrCl,PrKa,PrKb] = getconductances(xnew,np,par,surftrilist,apicaltrilist,basaltrilist,...
                                         apicalarea,basalarea,apicaltriarea,basaltriarea);

    f_secretion = @(t,x) secretion(t,x,par,4,PrCl,PrKa,PrKb);           % always using cell 4 for now
    [~,SS] = ode15s(f_secretion, [0,delt], IC, options);                % step forward the secretion model (for given c, ip and h)

    % Now pull out the cell volume and its derivative, for use in the PDE
    % solve. Note that in the PDE we use the volume at the end of the time step, not
    % the beginning. Although this makes no difference that I can tell.
    
    % Also note that h is stepped forward using a forward Euler on every
    % grid point, even though it's only used on the apical membrane to
    % calculate the boundary flux. This is inefficient, but easier to
    % program.
    
    cellvol = SS(end,4);
    dum = secretion(0,SS(end,:),par,4,PrCl,PrKa,PrKb);
    cellvolprime = dum(4);                                              % the volume derivative is the 4th component of the ODE

    xold = xnew;
    react = get_load(xold,np,par,tets,ntets,tetvol,cellvol,cellvolprime,dist_to_apical,dist_to_basal,apicaltrilist,apicaltriarea);                        % make the load vector
    xnew(1:2*np) = Amat\(mass*xold(1:2*np)+ delt*react(1:2*np));        % c and ip are solved by matrix inversion, as they diffuse.
    xnew(2*np+1:3*np) = xold(2*np+1:3*np) + delt*react(2*np+1:3*np);    % Forward Euler for h. 


    % reset initial conditions, ready for the next loop. 
    IC = SS(end,:); 

    % Now store the solution at the end of each time step, for later plotting  
    sol(:,i+1) = xnew;
    SSsol(:,i+1) = SS(end,:);

    % Just to keep track of where you are
%    if (tim(i)==floor(tim(i))) 
        display(tim(i)) 
%    end

end

save(strcat('output_KCa',num2str(par.apicalKCadensity),'_VPLC',num2str(par.VPLC),'.mat'))
%save new_saliva_output.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load vector
function out=get_load(u,np,par,tets,ntets,tetvol,cellvol,cellvolprime,dist_to_apical,dist_to_basal,apicaltrilist,apicaltriarea)

%   First np components of u    -   c
%   Second np components of u   -   ip3
%   Third np components of u    -   h

c=u(1:np);
ip=u(np+1:2*np);
h=u(2*np+1:3*np);

load_c=zeros(np,1);
load_ip=zeros(np,1);

for K = 1:ntets
    loc2glb = tets(K,1:4);
    
    % get the reactions at each of the tetrahedron nodes. Each r will have
    % three components; the reaction terms for c, ip and h.
    r1 = getrhsreactions(c(loc2glb(1)),ip(loc2glb(1)),h(loc2glb(1)),par,loc2glb(1),cellvol,cellvolprime,dist_to_apical,dist_to_basal);
    r2 = getrhsreactions(c(loc2glb(2)),ip(loc2glb(2)),h(loc2glb(2)),par,loc2glb(2),cellvol,cellvolprime,dist_to_apical,dist_to_basal);
    r3 = getrhsreactions(c(loc2glb(3)),ip(loc2glb(3)),h(loc2glb(3)),par,loc2glb(3),cellvol,cellvolprime,dist_to_apical,dist_to_basal);
    r4 = getrhsreactions(c(loc2glb(4)),ip(loc2glb(4)),h(loc2glb(4)),par,loc2glb(4),cellvol,cellvolprime,dist_to_apical,dist_to_basal);

    b_cal = [r1(1); r2(1); r3(1); r4(1)]/4*tetvol(K);      % local element load vector.
    b_ip = [r1(2); r2(2); r3(2); r4(2)]/4*tetvol(K);
    load_c(loc2glb) = load_c(loc2glb) + b_cal;          % add element loads to global load vector
    load_ip(loc2glb) = load_ip(loc2glb) + b_ip;
end

% Now add in the boundary conditions. These are obtained by integrating the
% IPR flux over the apical triangles, and adding these integrals to the
% appropriate entries in the load vector.

for K = 1:size(apicaltrilist,1)
    loc2glb = apicaltrilist(K,1:3);
    r1 = getbndyreact(c(loc2glb(1)),ip(loc2glb(1)),h(loc2glb(1)),par,cellvol);
    r2 = getbndyreact(c(loc2glb(2)),ip(loc2glb(2)),h(loc2glb(2)),par,cellvol);
    r3 = getbndyreact(c(loc2glb(3)),ip(loc2glb(3)),h(loc2glb(3)),par,cellvol);
    
    bndy_cal = [r1; r2; r3]/3*apicaltriarea(K);
    load_c(loc2glb) = load_c(loc2glb) + bndy_cal;
end

out(1:np,1)             =    load_c;       
out(np+1:2*np,1)        =    load_ip;
 
% The variables with no diffusion can be done differently. No loop over
% the triangles is needed, just a simpler loop over the nodes, but then
% time must be stepped forward using a forward Euler, not with a solve
% using the mass matrix.

% This is inefficient, as it steps forward h on all the grid points, even
% though h is only used on the apical membrane. But easier to program.
% 
  for i=1:np
      reactions=getrhsreactions(c(i),ip(i),h(i),par,i,cellvol,cellvolprime,dist_to_apical,dist_to_basal);
      out(2*np+i,1)   =   reactions(3);                     
  end

end
%-----------------------------------------

%% mass and stiffness matrices
function [mass,stiff] = makematrices(p,np,tets,ntets)

Dc=5; Dp=280;

stiffc=sparse(np,np);
stiffp=sparse(np,np);
small_mass=sparse(np,np);


small_mass = sparse(np,np); % allocate mass matrix
for K = 1:ntets % loop over elements
    loc2glb = tets(K,1:4); % local-to-global map
    P1 = p(loc2glb(1),:); P2 = p(loc2glb(2),:); P3 = p(loc2glb(3),:);   P4 = p(loc2glb(4),:);   % mostly just for convenience
    TT = [1 P1; 1 P2; 1 P3; 1 P4];
    volume = abs(det(TT))/6;
    
    % first the mass matrix
    MK = [2 1 1 1;
    1 2 1 1;
    1 1 2 1;
    1 1 1 2]/20*volume; % element mass matrix
    small_mass(loc2glb,loc2glb) = small_mass(loc2glb,loc2glb) + MK; % add element masses to M
    
    % then the stiffness matrix
    dum = inv(TT'); b = dum(:,2); c = dum(:,3); d = dum(:,4);   % This method of getting the coeffs follows Gockenbach
    Ac = Dc*(b*b' + c*c' + d*d')*volume;                        % element stiffness matrix
    Ap = Dp*(b*b' + c*c' + d*d')*volume;
    stiffc(loc2glb,loc2glb) = stiffc(loc2glb,loc2glb) + Ac;     % add element stiffnesses to global stiffness matrix
    stiffp(loc2glb,loc2glb) = stiffp(loc2glb,loc2glb) + Ap;
end


% Since h doesn't diffuse, it is solved directly by the forward Euler, and
% doesn't need entries in the mass and stiffness matrices.
Z=zeros(np);
mass  = sparse([ small_mass Z; Z small_mass]);
stiff = sparse([ stiffc Z; Z stiffp]);

end

%% Initialisation and reaction terms

function [IC,x] = initialise(np)

% Initialise the solution vector of the PDE
% You have to be careful with the structure of the solution vector
%
%   First np components of u    -   c
%   Second np components of u   -   ip
%   Third np components of u    -   h

c0  = 0.07;
ip0 = 0;
h0  = 0.619; 
x(1:np,1)=c0;    
x(np+1:2*np,1)=ip0;  
x(2*np+1:3*np,1)=h0; 

% Now initialise the ODE variables of the secretion model.

Nal_0     = 116.7; %116.719380654555;                       % Na in the lumen
Kl_0      = 6.8; %6.60670400306032;                         % K in the lumen
Cll_0     = Nal_0+Kl_0; %123.6; % 123.32608465962;          % Cl in the lumen
w_0       = 1123.8; % 1131.87051332195;                     % cell volume
Na_0      = 27.213; % 26.9628344160249;                     % Na in the cell
K_0       = 117.28; % 117.360283287559;                     % K in the cell
Cl_0      = 50.0; % 50.3594085578432;                       % Cl in the cell
HCO_0     = 10.63; % 10.6614665560422;                      % bicarbonate in the cell
H_0       = 0.0001425; %0.00014270477477718;                % H ions in the cell (determines the pH)
Va_0      = -58.46; %-59.6686463447929;                     % apical membrane potential
Vb_0      = -54.81; % -55.9511234948445;                    % basal membrane potential

IC = [Nal_0,Kl_0,Cll_0,w_0,Na_0,K_0,Cl_0,HCO_0,H_0,Va_0,Vb_0];   % The initial condition
        
        
end

%------------------------------------------------
function reactions = getrhsreactions(c,ip,h,par,point,vol,volprime,dist_to_apical,dist_to_basal)  

% nodepoint is the point (as an (x,y,z) vector) at which we are calculating the reaction
% terms. We can use this to define spatial heterogeneity.

% vol is the volume calculate from the secretion model, and volprime is its
% derivative. These are used in the volume correction terms in the ODE.

ce = (1/par.gamma) * ( ( par.w0{4} / vol )* par.ct - c );        % closed-cell calcium model

% IPR
H_inf = par.K_h^4 / ( par.K_h^4 + c^4 );
TAU = par.tau_max*par.K_tau^4/(par.K_tau^4+c^4);

% serca
Jserca = par.V_p*(c^2-par.K_bar*ce^2)/(c^2 + par.K_p^2);

% Turn on PLC only close to the basal region. Can comment this out mostly,
% as IP3 diffuses so fast there is no spatial distribution at all.

% if (dist_to_basal(point)>1 && par.spatialchoose==1) 
%     vplc = 0;
% else
    vplc=par.VPLC;
% end

% Now define the ODE RHS. Note how they include the volume factor
% corrections. But the only calcium body reaction is the serca pump, as the
% IPR is a boundary flux only.
reactions(1) = (- Jserca)*par.w0{4}/vol - c*volprime/vol;
reactions(2) = (vplc-par.V_5K*ip)*par.w0{4}/vol - ip*volprime/vol;
reactions(3) = (H_inf - h) / TAU;

end

%------------------------------------------------
function IPRflux = getbndyreact(c,ip,h,par,vol)  

% function for calculating the IPR boundary flux

ce = (1/par.gamma) * ( ( par.w0{4} / vol )* par.ct - c );        % closed-cell calcium model

phi_c = c^4 / ( c^4 + par.K_c^4 );
phi_p = ip^2 / ( par.K_p^2 + ip^2);
phi_p_down = par.K_p^2 / ( par.K_p^2 + ip^2);
H_inf = par.K_h^4 / ( par.K_h^4 + c^4 );
TAU = par.tau_max*par.K_tau^4/(par.K_tau^4+c^4);
beta = phi_p * phi_c * h;
alpha = phi_p_down * ( 1 - phi_c * H_inf );
po = beta/(beta+par.k_beta*(beta+alpha));

IPRflux = (par.k_f*po )*(ce-c)*par.w0{4}/vol;

end

%% conductances
function [PrCl,PrKa,PrKb]=getconductances(u,np,par,surftrilist,apicaltrilist,basaltrilist,...
                                         apicalarea,basalarea,apicaltriarea,basaltriarea)

% Ca2+ activated channels. These are spatially dependent, so have to be
% determined by integrating over the relevant surface triangles. All
% triangle areas and total areas are precalculated when the mesh is made,
% for efficiency.

% The scaling of the KCa apical and basal densities is a bit strange, and
% very ad hoc at the minute. Basically, if you change one you need to
% change the other also so that the resting cell volume comes out
% correctly. I do this by keeping the total resting KCa conductance constant.

c = u(1:np);

totalarea = apicalarea + basalarea;
total_KCa = totalarea/4;                            % determined by fiddling to get nice resting volume

apical_TMEM_density = 1;
apical_KCa_density = par.apicalKCadensity;
basal_KCa_density = (total_KCa - apicalarea*apical_KCa_density)/basalarea;

PrCl = 0;    % initialise
PrKa = 0;
for i=1:size(apicaltrilist,1)
    pts = surftrilist(apicaltrilist(i),1:3);
    cav = mean(c(pts));
    PrCl = PrCl + apicaltriarea(i) / ( 1 + ( par.KCaCC / cav )^par.eta1 ) / apicalarea;
    PrKa = PrKa + apicaltriarea(i) / ( 1 + ( par.KCaKC/cav )^par.eta2 ) / apicalarea;
end
PrCl = apical_TMEM_density * PrCl;
PrKa = apical_KCa_density * PrKa;


% Then integrate the KCa current over the basal region
PrKb = 0;
for i=1:size(basaltrilist,1)
    pts = surftrilist(basaltrilist(i),1:3);
    cav = mean(c(pts));  
    PrKb = PrKb + basaltriarea(i) / ( 1 + ( par.KCaKC/cav )^par.eta2 ) / basalarea;     % Total KCa conductance
end
PrKb = basal_KCa_density * PrKb;

end


