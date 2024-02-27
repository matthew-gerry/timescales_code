% lepe_2ls.m

% Using the LEPE method (Phys Rev E 108, 014130, 2023) to assess the
% behaviour of coherences and populations in transient dynamics for a 
% spin-boson model in the high temperature limit (nearly degenerate states).

% System operator in the interaction Hamiltonian Z + X or:
%  [[1,1]
%   [1,-1]]/sqrt(2)
% so as to give nontrivial behaviour with respect to coherences (not pure
% decoherence or equal probabilities to transition to the two different
% basis states - which turns out to be roughly the same thing in the high-
% temp limit).

% Matthew Gerry, July 2023

%%% LEPE - SYMBOLIC MANIPULATION %%%

syms k D t 

% Eigenvalues of the Liouvillian to leading order in D (done by hand)
l1 = -D^2/(8*k); l2 = -4*k-D/sqrt(2); l3 = -4*k+D/sqrt(2);
lambda = [l1, l2, l3];

Lambda = [[1,1,1];
          [l1, l2, l3];
          [l1^2, l2^2, l3^2]];

% Matrix for expressing derivatives of the 'polarization' P=0.5*(rho_11-rho_00)
% in terms of elements of the initial state
B1 = [[1,0,0];
      [-2*k, -2*k, 0];
      [8*k^2, 8*k^2, -2*k*D]];

% Matrix for expressing derivatives of the real part of the coherence in
% terms of elements of the initial state
B2 = [[0, 1, 0];
      [-2*k, -2*k, D];
      [8*k^2, 8*k^2-D^2, -6*k*D]];

% And for the imag part of the coherence
B3 = [[0,0,1];
      [0, -D, -4*k];
      [2*k*D, 6*k*D, 16*k^2-D^2]];

% Initial state in the format of [0.5*(rho_11-rho_0); Re(rho_10); Im(rho_10)]
x0 = [-1/2; 0; 0];
x1_inf = [0; 0; 0]; % Polarization and its derivatives at long times (everything goes to zero)
x2_inf = [0; 0; 0]; % Real part of coherence and its derivatives at long times
x3_inf = [0; 0; 0]; % Imag part of coherence and its derivatives at long times

% Solve for polarization and real and imag parts of coherence using LEPE
c1_general = Lambda\(B1*x0 - x1_inf);
c1 = simplify(subs(c1_general,D,0)); % Limit of small Delta

c2_general = Lambda\(B2*x0 - x2_inf);
c2 = simplify(subs(c2_general,D,0)); % Limit of small Delta

c3_general = Lambda\(B3*x0 - x3_inf);
c3 = simplify(subs(c3_general,D,0)); % Limit of small Delta

% Plug into expressions for these quantities in terms of the Liouvillian
% eigenvalues
P = sum(exp(lambda*t).*c1');
rho11 = 1/2 + P; % Population of the excited state
rho10_real = sum(exp(lambda*t).*c2'); % Real part of the coherence
rho10_imag = sum(exp(lambda*t).*c3'); % Imag part of the coherence

% Figure out the values of population and coherence in the transient regime
rho11_transient = subs(rho11,D,0);
rho11_transient = subs(rho11_transient,t,Inf);
rho10_real_transient = subs(rho10_real,D,0);
rho10_real_transient = subs(rho10_real_transient,t,Inf);

% Substitute in values for plotting
a = 0.02;
T = 1;
k_val = a*T;
D_val = 0.001;
time = logspace(-2,7,1000);

% Substitute values into symbolic expressions and convert to double
rho11_plot = subs(rho11, [k,D],[k_val,D_val]);
rho11_plot = double(subs(rho11_plot,t,time));

rho10_real_plot = subs(rho10_real,[k,D],[k_val,D_val]);
rho10_real_plot = double(subs(rho10_real_plot,t,time));

rho10_imag_plot = subs(rho10_imag,[k,D],[k_val,D_val]);
rho10_imag_plot = double(subs(rho10_imag_plot,t,time));

% Also plot the handwritten expression I have in my notes as a sanity check
rho10_real_expression = 0.25*exp(-D_val^2*time/(8*k_val)) - 0.25*exp(-4*k_val*time).*cosh(D_val*time/sqrt(2));
rho10_imag_expression = 0.25*sqrt(2)*exp(-4*k_val*time).*sinh(D_val*time/sqrt(2));

%%% FULL UQME - FULLY NUMERICAL %%%
% For comparison to analytic results from LEPE
% This is easy to do since the equations of motion for x (simplified way to
% express the state) are still homogeneous for this model

% Liovillian
L = -2*k_val*[[1,1,0];
              [1,1,0];
              [0,0,2]] + D_val*[[0,0,0];
                                [0,0,1]
                                [0,-1,0]];

% Pre-allocate array to store state as a function of time
x = zeros(length(x0),length(time));
for ii=1:length(time)
    x(:,ii) = expm(L*time(ii))*x0;
end

% Get components of the state from the result
rho10_real_uqme = x(2,:);
rho10_imag_uqme = x(3,:);
rho11_uqme = 1/2 + x(1,:);

%%% PLOTTING %%%

% Timescales to include on the plots
tau1 = 4*k_val/(D_val^2);
tau2 = 1/(2*k_val);

% Plot each element of the density operator as desired
subplot(2,1,1,xscale="log"); hold on; box on
semilogx(time, rho11_uqme, linewidth=1.5, DisplayName="Unified QME")
semilogx(time, rho11_plot, '--k', DisplayName="Analytic")
xlim([min(time), max(time)])
ylim([-0.01, 0.56])
xticks([1,1e3,1e6])
% xlabel("$t$",Interpreter="latex")
ylabel("$\rho_{11}(t)$",Interpreter="latex")
taulines = xline([tau1,tau2],'--',["$\tau_1$","$\tau_2$"],Interpreter='latex',fontsize=14,LabelOrientation='horizontal');
text(0.05,0.5,"(a)",Interpreter="latex",FontSize=14)
% legend(location='northwest',Interpreter='latex')
set(gca,fontsize=14)
hold off

subplot(2,1,2,xscale="log"); hold on; box on
semilogx(time,rho10_real_uqme, color=[0.64,0.08,0.18], Linewidth=1.5, DisplayName="Unified QME")
semilogx(time,rho10_real_plot, '--k', linewidth=1,DisplayName="Analytic")
% semilogx(time,rho10_real_expression, '--k')
xlim([min(time), max(time)])
ylim([-0.01, 0.3])
xticks([1,1e3,1e6])
xlabel("$t$",Interpreter="latex")
ylabel("$\rho_{10}^R(t)$",Interpreter="latex")
taulines = xline([tau1,tau2],'--',["$\tau_1$","$\tau_2$"],Interpreter='latex',fontsize=14,LabelOrientation='horizontal');
% legend(location='northwest',Interpreter='latex')
text(0.05,0.26,"(b)",Interpreter="latex",FontSize=14)
set(gca,fontsize=14)
hold off

% subplot(3,1,3,xscale="log"); hold on; box on
% semilogx(time,rho10_imag_plot, linewidth=1.5)
% semilogx(time,rho10_imag_expression, '--k')
% semilogx(time,rho10_imag_uqme,'--r')
% xlim([min(time), max(time)])
% ylim([-0.01, 0.026])
% xticks([1,1e3,1e6])
% xlabel("$t$",Interpreter="latex")
% ylabel("Im$[\rho_{10}(t)]$",Interpreter="latex")
% set(gca,fontsize=14)
% hold off


%% DOING IT ALL AGAIN IN A NEW BASIS
% We have changed to the basis where the interaction looks like Z and the 
% system Hamiltonian in turn looks like X + Z

% This changes the Liouvillian eigenvalues and very much changes the
% equations for the different aspects of the state

syms k D t O22 O23
assume(k, 'real')
assume(D, 'real')

% Eigenvalues of the Liouvillian to leading order in D (done by hand)
l1 = -D^2/(4*k); l2 = -2*k + O22; l3 = -2*k + O23;
% To avoid a singularity, we include the next order contributions to l2 and
% l3 but remain agnostic as to their value. We know they should be of order
% ~D which is small, so they'll be taken to zero later.

lambda = [l1, l2, l3];

Lambda = [[1,1,1];
          [l1, l2, l3];
          [l1^2, l2^2, l3^2]];

% Matrix for expressing derivatives of the 'polarization' P=0.5*(rho_11-rho_00)
% in terms of elements of the initial state
B1 = [[1,0,0];
      [0, 0, D/sqrt(2)];
      [-D^2/2, D^2/2, -sqrt(2)*k*D]];

% Matrix for expressing derivatives of the real part of the coherence in
% terms of elements of the initial state
B2 = [[0, 1, 0];
      [0, -2*k, -D/sqrt(2)];
      [D^2/2, 4*k^2-D^2/2, 2*sqrt(2)*k*D]];

% And the imag part of the coherence
B3 = [[0, 0, 1];
      [-D/sqrt(2), D/sqrt(2),-2*k];
      [sqrt(2)*k*D, -2*sqrt(2)*k*D, -D^2 + 4*k^2]];

% Initial state in the format of [0.5*(rho_00-rho_11); Re(rho_10); Im(rho_10)]
x0 = [-1/(2*sqrt(2)); -1/(2*sqrt(2)); 0]; % Corresponding to the system's ground state

% x0 = [1/2; 0; -1i/(2*sqrt(2))]; % Unphysical intial state used in sigma_z correlation function calculation

% Steady state elements and derivatives
x1_inf = [0; 0; 0]; % Polarization and its derivatives at long times (everything goes to zero)
x2_inf = [0; 0; 0]; % Real part of coherence and its derivatives at long times
x3_inf = [0; 0; 0]; % Imag part of coherence and its derivatives at long times

% Solve for polarization and real part of coherence using LEPE
c1_general = Lambda\(B1*x0 - x1_inf);
c1 = simplify(subs(c1_general,D,0)); % Limit of small Delta

c2_general = Lambda\(B2*x0 - x2_inf);
c2 = simplify(subs(c2_general,D,0)); % Limit of small Delta

c3_general = Lambda\(B3*x0 - x3_inf);
c3 = simplify(subs(c3_general,D,0)); % Limit of small Delta

% Plug into expressions for these quantities in terms of the Liouvillian
% eigenvalues
P = sum(exp(lambda*t).*c1');

rho10_real = sum(exp(lambda*t).*c2'); % Real part of the coherence
rho10_imag = sum(exp(lambda*t).*c3'); % Real part of the coherence

% Now drop the unknown next-order contributions to the larger eigenvalues -
% it is sufficient to set just one to zero, we choose O22
P = simplify(subs(P, O22, 0));
rho10_real = simplify(subs(rho10_real, O22, 0));
rho10_imag = simplify(subs(rho10_imag, O22, 0));

rho11 = 1/2 - P; % Population of the |+> state
rho00 = 1/2 + P; % Population of the |-> state

% Figure out the values of population and coherence in the transient regime
rho11_transient = subs(rho11,D,0);
rho11_transient = subs(rho11_transient,t,Inf);
rho10_real_transient = subs(rho10_real,D,0);
rho10_real_transient = subs(rho10_real_transient,t,Inf);

% Substitute in values for plotting
a = 0.02;
T = 1;
k_val = a*T;
D_val = 0.001;
time = logspace(-2,7,1000);

% Substitute values into symbolic expressions and convert to double
rho00_plot = subs(rho00, [k,D],[k_val,D_val]);
rho00_plot = double(subs(rho00_plot,t,time));

rho10_real_plot = subs(rho10_real,[k,D],[k_val,D_val]);
rho10_real_plot = double(subs(rho10_real_plot,t,time));

rho10_imag_plot = subs(rho10_imag,[k,D],[k_val,D_val]);
rho10_imag_plot = double(subs(rho10_imag_plot,t,time));


%%% FULL UQME - FULLY NUMERICAL %%%
% For comparison to analytic results from LEPE
% This is easy to do since the equations of motion for x (simplified way to
% express the state) are still homogeneous for this model

% Liovillian
L = -2*k_val*[[0,0,0];
              [0,1,0];
              [0,0,1]] + (D_val/sqrt(2))*[[0,0,1];
                                          [0,0,-1]
                                          [-1,1,0]];

% Pre-allocate array to store state as a function of time
x = zeros(length(x0),length(time));

for ii=1:length(time) % Solve the UQME numerically
    x(:,ii) = expm(L*time(ii))*x0;
end

% Get components of the state from the result
rho10_real_uqme_a = x(2,:);
rho10_imag_uqme_a = x(3,:);
rho11_uqme_a = 1/2 - x(1,:);
rho00_uqme_a = 1/2 + x(1,:);

%%% PLOTTING %%%

% Timescales to include on the plots
tau1 = 4*k_val/(D_val^2);
tau2 = 1/(2*k_val);

% Plot each element of the density operator in the transformed basis
subplot(2,1,1,xscale="log"); hold on; box on
semilogx(time, rho00_uqme_a, linewidth=1.5, DisplayName="Unified QME")
semilogx(time, rho00_plot, '--k', DisplayName="Analytic")
taulines = xline([tau1,tau2],'--',["$\tau_1$","$\tau_2$"],Interpreter='latex',fontsize=14,LabelOrientation='horizontal');
xlim([min(time), max(time)])
xticks([1,1e3,1e6])
ylim([0.0, 1.1])
% xlabel("$t$",Interpreter="latex")
ylabel("$\rho_{--}(t)$",Interpreter="latex")
% legend(location='southwest',Interpreter='latex')
text(0.05,0.96,"(a)",Interpreter="latex",FontSize=14)
set(gca,fontsize=14)
hold off

subplot(2,1,2,xscale="log"); hold on; box on
semilogx(time,rho10_real_uqme_a, color=[0.64,0.08,0.18], linewidth=1.5,DisplayName="Unified QME")
semilogx(time,rho10_real_plot, '--k', DisplayName="Analytic")
taulines = xline([tau1,tau2],'--',["$\tau_1$","$\tau_2$"],Interpreter='latex',fontsize=14,LabelOrientation='horizontal');
xlim([min(time), max(time)])
xticks([1,1e3,1e6])
ylim([-0.4,0.1])
xlabel("$t$",Interpreter="latex")
ylabel("$\rho_{+-}^R(t)$",Interpreter="latex")
% legend(location='northwest',Interpreter='latex')
text(0.05,0.042,"(b)",Interpreter="latex",FontSize=14)
set(gca,fontsize=14)
hold off

% subplot(3,1,3,xscale="log"); hold on; box on
% semilogx(time,rho10_imag_plot, linewidth=1.5)
% semilogx(time,rho10_imag_uqme_a,'--r')
% xlim([min(time), max(time)])
% xlabel("$t$",Interpreter="latex")
% ylabel("$\tilde{c}^I(t)$",Interpreter="latex")
% set(gca,fontsize=14)
% hold off

%% Investigate the equivalence of the results in the two bases (numerical only)

% Construct density matrices in the two bases based on numerical results of
% the previous two sections (must run both previous sections of code in
% order to run this section)

% Eigenbasis
rho_mat = zeros(2,2,length(time)); % Pre-allocate time-series for matrix
% Populate with numerically calculated values
rho_mat(1,1,:) = 1-rho11_uqme;
rho_mat(1,2,:) = rho10_real_uqme-1i*rho10_imag_uqme;
rho_mat(2,1,:) = rho10_real_uqme+1i*rho10_imag_uqme;
rho_mat(2,2,:) = rho11_uqme;

% Repeat steps for alternate basis
rho_mat_a = zeros(2,2,length(time)); 
% Populate with numerically calculated values
rho_mat_a(1,1,:) = 1-rho11_uqme_a;
rho_mat_a(1,2,:) = rho10_real_uqme_a-1i*rho10_imag_uqme_a;
rho_mat_a(2,1,:) = rho10_real_uqme_a+1i*rho10_imag_uqme_a;
rho_mat_a(2,2,:) = rho11_uqme_a;


% Prepare change of basis matrix for transform
V = (1i/sqrt(2))*[[-(sqrt(1-1/sqrt(2))), sqrt(1+1/sqrt(2))];
                [1/(sqrt(2)*sqrt(1-1/sqrt(2))), 1/(sqrt(2)*sqrt(1+1/sqrt(2)))]];

% Pre-allocate time series of matrices to be obtained by transforming from
% the alternate basis back to the eigenbasis
rho_mat_a_eig = zeros(2,2,length(time));
% Transform back
for ii=1:length(time)
    rho_mat_a_eig(:,:,ii) = V*rho_mat_a(:,:,ii)/V;
end

% Reshape data for plotting
rho11_a_eig = reshape(real(rho_mat_a_eig(2,2,:)),[1,length(time)]);
rho10_real_a_eig = reshape(real(rho_mat_a_eig(2,1,:)),[1,length(time)]);
rho10_imag_a_eig = -reshape(imag(rho_mat_a_eig(1,2,:)),[1,length(time)]);

% Plot each matrix element in the eigenbasis as calculated using the UQME,
% alongside the same element as calculated using the UQME in the
% transformed basis and then transformed back, to analyze the effect of
% making the leading-order approximation before or after the transformation
figure(1)
subplot(1,3,1,xscale="log"); hold on; box on
semilogx(time, rho11_a_eig, '--r')
semilogx(time, rho11_uqme, '-b')
hold off

subplot(1,3,2,xscale="log"); hold on; box on
semilogx(time, rho10_real_a_eig, '--r')
semilogx(time, rho10_real_uqme, '-b')
hold off

subplot(1,3,3,xscale="log"); hold on; box on
semilogx(time, rho10_imag_a_eig, '--r')
semilogx(time, rho10_imag_uqme, '-b')
hold off