% lepe_vmodel.m

% Using the LEPE method (arxiv: 2301.06135) to assess the behaviour of
% coherences and populations in transient dynamics for the V model in the
% basis including the gound state, |1>, and the '+' and '-' superpositions
% of the two excited states: |+/->= (|2> +/- |3>)/sqrt(2), where |2> and
% |3> are nearly degenerate.

% System part of the interaction Hamiltonian in this basis: 
% S = [[0,1,0],[1,0,0],[0,0,0]] - jumps upwards and downwards between |1>
% and |+>. In this basis, all coherences involving leve
% l 1 and the real 
% part of the coherence rho_(+-) decouple from the rest of the dynamics,
% leaving a homogeneous system of three equations.

% Matthew Gerry, October 2023

% Variables used in the system of equations:
% P1 = (rho_(++) - exp(-b*nu)*rho_(11))/2 - Polarization between the
%                                           excited and ground states
%                                           weighted by temperature.
% P2 = (rho_(++) - rho_(--))/2            - Polarization between the two
%                                           exited states.
% cI                                      - Imaginary part of the coherence
%                                           rho_(+-).
% At steady state (Gibbs state), P1=P2=cI=0.


%%%% CONSTANTS %%%%

% Values we will substitute in for numerics, plotting later on
a = 0.02;
T = 1;
b_val = 1/T;
nu_val = 1;
k_val = a*nu_val/(exp(b_val*nu_val) - 1); % Spectral density times BE dist
D_val = 0.001;
time = logspace(-1,6,1000);

f_val = exp(-b_val*nu_val); Z_val = 1 + 2*f_val;

%%% LEPE - SYMBOLIC MANIPULATION %%%

syms k D b nu t

assume(b > 0)
assume(nu > 0)
assume(D > 0)
assume(k > 0)

f = exp(-b*nu); % Boltzmann factor for splitting
Z = 1 + 2*f; % Partition function

% Eigenvalues of the Liouvillian to leading order in D (done by hand)
l1 = -(1 + 2*f)*D^2/((1 + f)*k);
l2 = -k/2; l3 = -(1 + f)*k;
lambda = [l1, l2, l3];

Lambda = [[1,1,1];
          [l1, l2, l3];
          [l1^2, l2^2, l3^2]];

% Matrix for expressing derivatives of the P1 in terms of elements of the initial state
B1 = [[1, 0, 0];
      [-k*(1 + f), 0, -D/2];
      [k^2*(1 + f)^2, -D^2/2, k*D/4]];

% Analogous matrix for P2
B2 = [[0, 1, 0];
      [-k, 0, -D];
      [k^2*(1 + f), -D^2, k*D]];

% And for the imaginary part of the coherence
B3 = [[0, 0, 1];
      [0, D, -k/2];
      [-k*D, -k*D/2, k^2/4-D^2]];


% Initial state in terms of the matrix elements in the 1,+,- basis
% rho_init_11 = 0;
% rho_init_pp = 0.64;
% cI_init = 0.48;
rho_init_11 = 1;
rho_init_pp = 0;
cI_init = 0;

rho_init_mm = 1 - rho_init_11 - rho_init_pp;


% Initial state in the format of x = [P1; P2; cI]
P1_init = (rho_init_pp - exp(-b*nu)*rho_init_11)/2;
P2_init = (rho_init_pp - rho_init_mm)/2;
x_init = [P1_init; P2_init; cI_init];

% Components of x and their derivatives at steady state (all zero)
x1_inf = [0; 0; 0]; x2_inf = [0; 0; 0]; x3_inf = [0; 0; 0];


% Solve the LEPE matrix-vector equation to get the coefficients for each
% variable when written as a linear combination of exponentially decaying
% terms (with rates given by the Liouvillian eigenvalues).
a1_general = Lambda\(B1*x_init - x1_inf);
a1 = simplify(subs(a1_general,D,0)); % Limit of small Delta

a2_general = Lambda\(B2*x_init - x2_inf);
a2 = simplify(subs(a2_general,D,0)); % Limit of small Delta

a3_general = Lambda\(B3*x_init - x3_inf);
a3 = simplify(subs(a3_general,D,0)); % Limit of small Delta

% Plug into expressions for these quantities in terms of the Liouvillian
% eigenvalues
P1 = sum(exp(lambda*t).*a1');
P2 = sum(exp(lambda*t).*a2');
cI = sum(exp(lambda*t).*a3');

% Density operator matrix elements in terms of our variables:
rho11 = (2*(P2 - 2*P1) + 1)/Z;
rhopp = (2*(P1 + f*P2) + f)/Z;
rhomm = (f + 2*P1 - 2*(1+f)*P2)/Z;


%%%% UNIFIED QME - FULLY NUMERICAL %%%%
% Compare results of just solving the UQME numerically with LEPE
L = [[-k_val*(1 + f_val), 0, -D_val/2];
     [-k_val, 0, -D_val];
     [0, D_val, -k_val/2]]; % Liouvillian

% Numeric version of initial state
P1_init_u = (rho_init_pp - exp(-b_val*nu_val)*rho_init_11)/2;
P2_init_u = (rho_init_pp - rho_init_mm)/2;
x_init_u = [P1_init_u; P2_init_u; cI_init];

% Pre-allocate array to store state as a function of time
x_u = zeros(length(x_init_u),length(time));
% Numerically solve the UQME
for ii=1:length(time)
    x_u(:,ii) = expm(L*time(ii))*x_init_u;
end

rho11_u = (2*(x_u(2,:) - 2*x_u(1,:)) + 1)/Z_val;
rhopp_u = (2*(x_u(1,:) + f_val*x_u(2,:)) + f_val)/Z_val;
rhomm_u = (f_val + 2*x_u(1,:) - 2*(1+f_val)*x_u(2,:))/Z_val;
cI_u = x_u(3,:);

%%%% PLOTTING %%%%

% Timescales to include in plots
tau1 = (1 + f_val)*k_val/(Z_val*D_val^2);
tau2 = 1/(k_val*(1 + f_val));

% Plug in numbers and evaluate as a function of time
rho11_plot = subs(rho11, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
rho11_plot = double(subs(rho11_plot,t,time));

rhopp_plot = subs(rhopp, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
rhopp_plot = double(subs(rhopp_plot,t,time));

rhomm_plot = subs(rhomm, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
rhomm_plot = double(subs(rhomm_plot,t,time));

cI_plot = subs(cI, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
cI_plot = double(subs(cI_plot,t,time));

% Plot the matrix elements
figure(1)
subplot(1,1,1,xscale="log")
hold on; box on
plot(nu_val*time, rho11_u, DisplayName="$\rho_{11}$", linewidth=1.5)
ana11 = plot(nu_val*time, rho11_plot, '--k', DisplayName="$\rho_{11}$");
ana11.Annotation.LegendInformation.IconDisplayStyle = "off";
plot(nu_val*time, rhopp_u, DisplayName="$\rho_{++}$", linewidth=1.5);
anapp = plot(nu_val*time, rhopp_plot, '--k', DisplayName="$\rho_{++}$");
anapp.Annotation.LegendInformation.IconDisplayStyle = "off";
plot(nu_val*time, rhomm_u, DisplayName="$\rho_{--}$", linewidth=1.5)
anamm = plot(nu_val*time, rhomm_plot, '--k', DisplayName="$\rho_{--}$");
anamm.Annotation.LegendInformation.IconDisplayStyle = "off";
taulines = xline([tau1,tau2],'--',["$\tau_1$","$\tau_2$"],Interpreter='latex',fontsize=14,LabelOrientation='horizontal',HandleVisibility="off");

xlim([nu_val*min(time), nu_val*max(time)])
xlabel("$\nu t$", Interpreter="latex")
legend(Location="northeast", Interpreter="latex")
set(gca, fontsize=14)
hold off

% subplot(1,2,2,xscale="log"); hold on; box on
% plot(time, cI_u, DisplayName="Im[$\rho_{+-}$]", linewidth=1.5)
% anapm = plot(time, cI_plot, '--k', DisplayName="Im[$\rho_{+-}$]");
% anapm.Annotation.LegendInformation.IconDisplayStyle = "off";
% xlim([nu_val*min(time), nu_val*max(time)])
% xlabel("$\nu t$", Interpreter="latex")
% legend(Location="northeast", Interpreter="latex")
% set(gca, fontsize=14)
% hold off

%% Same thing, but in the energy eigenbasis

% System part of the interaction Hamiltonian in this basis: 
% S = [[0,1,1],[1,0,0],[1,0,0]]/sqrt(2)

% In this case, if we start with rho_(22) = rho_(33), they remain equal at
% all times. We consider that case here, and eliminate the variable
% representing the population difference between them. However, the real
% part of the coherence couples to the rest of the dynamics so we still
% have three variables.

% Matthew Gerry, October 2023

% Variables used in the system of equations:
% P1 = (rho_(22) - exp(-b*nu)*rho_(11))/2 - Polarization between the
%                                           excited and ground states
%                                           weighted by temperature.
% cR                                      - Real part of the coherence
%                                           rho_(32).
% cI                                      - Imaginary part of the coherence
%                                           rho_(32).
% At steady state (Gibbs state), P1=P2=cI=0.


%%%% CONSTANTS %%%%

% Values we will substitute in for numerics, plotting later on
a = 0.02;
T = 1;
b_val = 1/T;
nu_val = 1;
k_val = a*nu_val/(exp(b_val*nu_val) - 1); % Spectral density times BE dist
D_val = 0.001;
time = logspace(-1,6,1000);

f_val = exp(-b_val*nu_val); Z_val = 1 + 2*f_val;

%%% LEPE - SYMBOLIC MANIPULATION %%%

syms k D b nu t

assume(b > 0)
assume(nu > 0)
assume(D > 0)
assume(k > 0)

f = exp(-b*nu); % Boltzmann factor for splitting
Z = 1 + 2*f; % Partition function

% Eigenvalues of the Liouvillian to leading order in D (done by hand)
l1 = -(1 + 2*f)*D^2/((1 + f)*k);
l2 = -k/2; l3 = -(1 + f)*k;
lambda = [l1, l2, l3];

Lambda = [[1,1,1];
          [l1, l2, l3];
          [l1^2, l2^2, l3^2]];

% Matrix for expressing derivatives of the P1 in terms of elements of the initial state
B1 = [[1, 0, 0];
      [-k*Z/2, -k*Z/2, 0];
      [(k*Z)^2/4 + k^2*Z/4, (k*Z)^2/4 + k^2*Z/4, -k*Z*D/2]];

% Analogous matrix for the real part of the coherence
B2 = [[0, 1, 0];
      [-k/2, -k/2, D];
      [Z*k^2/4 + k^2/4, k^2*Z/2 + k^2/4 - D^2, -k*D]];

% And for the imaginary part of the coherence
B3 = [[0, 0, 1];
      [0, -D, -k/2];
      [k*D/2, k*D, k^2/4 - D^2]];


% Initial state in terms of the matrix elements in the 1,+,- basis
% rho_init_11 = 0;
% rho_init_pp = 0.64;
% cI_init = 0.48;
rho_init_11 = 1;
cR_init = 0;
cI_init = 0;

% Initial state in the format of x = [P1; cR; cI]
P1_init = Z*(1 - rho_init_11)/2 - f;
x_init = [P1_init; cR_init; cI_init];

% Components of x and their derivatives at steady state (all zero)
x1_inf = [0; 0; 0]; x2_inf = [0; 0; 0]; x3_inf = [0; 0; 0];


% Solve the LEPE matrix-vector equation to get the coefficients for each
% variable when written as a linear combination of exponentially decaying
% terms (with rates given by the Liouvillian eigenvalues).
a1_general = Lambda\(B1*x_init - x1_inf);
a1 = simplify(subs(a1_general,D,0)); % Limit of small Delta

a2_general = Lambda\(B2*x_init - x2_inf);
a2 = simplify(subs(a2_general,D,0)); % Limit of small Delta

a3_general = Lambda\(B3*x_init - x3_inf);
a3 = simplify(subs(a3_general,D,0)); % Limit of small Delta

% Plug into expressions for these quantities in terms of the Liouvillian
% eigenvalues
P1 = sum(exp(lambda*t).*a1');
cR = sum(exp(lambda*t).*a2');
cI = sum(exp(lambda*t).*a3');

% Density operator matrix elements in terms of our variables:
rho11 = 1 - 2*(P1 + f)/Z;
rho22 = (1 - rho11)/2;
rho33 = rho22; % This fact is initial state dependent


%%%% UNIFIED QME - FULLY NUMERICAL %%%%
% Compare results of just solving the UQME numerically with LEPE
L = [[-k_val*Z_val/2, -k_val*Z_val/2, 0];
     [-k_val/2, -k_val/2, D_val];
     [0, -D_val, -k_val/2]]; % Liouvillian

% Numeric version of initial state
P1_init_u = Z_val*(1 - rho_init_11)/2 - f_val;
x_init_u = [P1_init_u; cR_init; cI_init];

% Pre-allocate array to store state as a function of time
x_u = zeros(length(x_init_u),length(time));
% Numerically solve the UQME
for ii=1:length(time)
    x_u(:,ii) = expm(L*time(ii))*x_init_u;
end

% Re-derive the density matrix elements
rho11_u = 1 - 2*(x_u(1,:) + f_val)/Z_val;
rho22_u = (1 - rho11_u)/2;
cR_u = x_u(2,:);
cI_u = x_u(3,:);


%%%% PLOTTING %%%%

% Timescales to include in plots
tau1 = (1 + f_val)*k_val/(Z_val*D_val^2);
tau2 = 1/(k_val*(1 + f_val));

% Plug in numbers and evaluate as a function of time
rho11_plot = subs(rho11, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
rho11_plot = double(subs(rho11_plot,t,time));

rho22_plot = subs(rho22, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
rho22_plot = double(subs(rho22_plot,t,time));

cR_plot = subs(cR, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
cR_plot = double(subs(cR_plot,t,time));

cI_plot = subs(cI, [k,D,b,nu],[k_val,D_val,b_val,nu_val]);
cI_plot = double(subs(cI_plot,t,time));

% Plot the matrix elements

figure(1)
subplot(1,2,1,xscale="log"); hold on; box on
plot(nu_val*time, rho11_u, DisplayName="$\rho_{11}$", LineWidth=1.5)
curve1 = plot(nu_val*time, rho11_plot, '--k', DisplayName="Analytic (perturbative)");
curve1.Annotation.LegendInformation.IconDisplayStyle = "off";
plot(nu_val*time, rho22_u, DisplayName="$\rho_{22}$", LineWidth=1.5)
curve2 = plot(nu_val*time, rho22_plot, '--k');
curve2.Annotation.LegendInformation.IconDisplayStyle = "off";
taulines = xline([tau1,tau2],'--',["$\tau_1$","$\tau_2$"],Interpreter='latex',fontsize=14,LabelOrientation='horizontal',HandleVisibility="off");
% plot(nu_val*time, rhomm_plot, DisplayName="$\rho_{--}$")
% plot(nu_val*time, rhomm_u, '--r', DisplayName="$\rho_{--}$ (numerical)")
xticks([1e-1,1e1,1e3,1e5])
xlim([nu_val*min(time), nu_val*max(time)])
xlabel("$\nu t$", Interpreter="latex")
legend(Location="west", Interpreter="latex", fontsize=14)
text(0.3,0.9,"(a)",Interpreter="latex",FontSize=14)
set(gca, fontsize=14)
hold off

subplot(1,2,2,xscale="log"); hold on; box on
plot(time, cR_u, DisplayName="$\rho_{32}^R$", color=[0.56,0.17,0.92], LineWidth=1.5)
curve3 = plot(time, cR_plot, '--k', DisplayName="Analytic (perturbative)");
curve3.Annotation.LegendInformation.IconDisplayStyle = "off";
plot(time, cI_u, DisplayName="$\rho_{32}^I$", color=[0.10,0.70,0.12], LineWidth=1.5)
curve4 = plot(time, cI_plot, '--k');
curve4.Annotation.LegendInformation.IconDisplayStyle = "off";
taulines = xline([tau1,tau2],'--',["$\tau_1$","$\tau_2$"],Interpreter='latex',fontsize=14,LabelOrientation='horizontal',HandleVisibility="off");
xticks([1e-1,1e1,1e3,1e5])
xlim([nu_val*min(time), nu_val*max(time)])
xlabel("$\nu t$", Interpreter="latex")
legend(Location="west", Interpreter="latex", fontsize=14)
text(0.3,0.13,"(b)",Interpreter="latex",FontSize=14)
set(gca, fontsize=14)
hold off















