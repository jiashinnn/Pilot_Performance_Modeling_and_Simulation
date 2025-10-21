clc
clear;


% Time vector
% duration
T = 50;
% change rate (0 to 1)
dt = 0.1;
% time steps
numSteps = T / dt;
% vectorization of time space
time = linspace(0, T, numSteps);

% Initialize arrays to store the variables over time
HR = zeros(1, numSteps); % Heart Rate
SQ = zeros(1, numSteps); % Sleep Quality
MC = zeros(1, numSteps); % Mission Complexity
ES = zeros(1, numSteps); % Environmental Stressor
EL = zeros(1, numSteps); % Experienced Level

% Initial values
for t = 1:numSteps
  HR(t) = 0.80;
  SQ(t) = 0.05;
  MC(t) = 0.85;
  EL_base(t) = 0.10;
  ES(t) = 0.88;
end

% Initialize arrays for the other variables
CL = zeros(1, numSteps);
SLs = zeros(1, numSteps);
PF = zeros(1, numSteps);
RT = zeros(1, numSteps);
SA = zeros(1, numSteps);
Ps = zeros(1, numSteps);
SLl = zeros(1, numSteps);
Pl = zeros(1, numSteps);

% Initial values for long-term variables
SLl(1) = 0.0;
Pl(1) = 0.0;

% Define parameters (tune these as needed)
alpha_CL = 0.3;
w_CL1 = 0.3;
w_CL2 = 0.3;
w_CL3 = 0.4;
beta_CL = 0.3;

w_SLs1 = 0.3;
w_SLs2 = 0.3;
w_SLs3 = 0.4;
lamda_SLs = 0.3;
alpha_SLs = 0.4;

alpha_PF = 0.5;
w_PF1 = 0.5;
w_PF2 = 0.5;

gamma_RT = 0.5;
w_RT1 = 0.5;
w_RT2 = 0.5;
w_RT3 = 0.5;
w_RT4 = 0.5;

beta_SA = 0.5;
w_SA1 = 0.3;
w_SA2 = 0.3;
w_SA3 = 0.3;

delta_EL = 0.5;

eta_Ps = 0.5;
w_Ps1 = 0.5;
w_Ps2 = 0.5;
w_Ps3 = 0.5;
w_Ps4 = 0.5;

CL(1) = (1 - alpha_CL * EL(1)) * ((w_CL1 * MC(1) + w_CL2 * ES(1) + w_CL3 * HR(1)) - beta_CL * SQ(1));

SLs(1) = ((w_SLs1 * HR(1) + w_SLs2 * ES(1) + w_SLs3 * MC(1)) - lamda_SLs * SQ(1)) * (1 - alpha_SLs * EL(1)) ;

PF(1) = alpha_PF * (w_PF1 * HR(1) + w_PF2 * (1 - SQ(1))) + (1 - alpha_PF) * SLl(1);

RT(1) = gamma_RT *( w_RT1*SQ(1) + w_RT2*SA(1)) + (1 - gamma_RT) * (1-(w_RT3*CL(1) + w_RT4*PF(1)));

SA(1) = (beta_SA *(w_SA1 * SQ(1) + w_SA2 * EL(1) + w_SA3 * Pl(1)) + (1 - beta_SA)*(1 - CL(1)));

Ps(1) = eta_Ps *(w_Ps1 * SA(1) + w_Ps2 * RT(1)) + (1- eta_Ps)*(1-(w_Ps3* SLs(1) + w_Ps4 * PF(1)));

EL(1) = delta_EL * EL_base(1) + (1 - delta_EL) * Pl(1);

for t = 2 : numSteps
    % Compute Cognitive Load (CL)
    CL(t) = (1 - alpha_CL * EL(t)) * ((w_CL1 * MC(t) + w_CL2 * ES(t) + w_CL3 * HR(t)) - beta_CL * SQ(t));

    % Compute Short-Term Stress Level (SLs)
    SLs(t) = ((w_SLs1 * HR(t) + w_SLs2 * ES(t) + w_SLs3 * MC(t)) - lamda_SLs*SQ(t)) * (1 - alpha_SLs*EL(t)) ;

    % Compute Physical Fatigue (PF)
    PF(t) = alpha_PF * (w_PF1 * HR(t) + w_PF2 * (1 - SQ(t))) + (1 - alpha_PF) * SLl(t-1);

    % Compute Reaction Time (RT)
    RT(t) = gamma_RT *( w_RT1*SQ(t) + w_RT2*SA(t)) + (1 - gamma_RT) * (1-(w_RT3*CL(t) + w_RT4*PF(t)));

    % Compute Situational Awareness (SA)
    SA(t) = beta_SA *(w_SA1 * SQ(t) + w_SA2 * EL(t) + w_SA3 * Pl(t-1))+((1- beta_SA)*(1- CL(t)));

    % Compute Short-Term Performance (Ps)
    Ps(t) = eta_Ps *(w_Ps1 * SA(t) + w_Ps2 * RT(t)) + (1- eta_Ps)*(1-(w_Ps3* SLs(t) + w_Ps4 * PF(t)));

    % Update Long-Term Stress Level (SLl)
    SLl(t) = SLl(t-1) + (SLs(t) - SLl(t-1)) * dt;
    SLl(t) = max(0, min(1, SLl(t))); % Ensure SLl stays within [0, 1]

    % Update Long-Term Performance (Pl)
    Pl(t) = Pl(t-1) + (Ps(t) - Pl(t-1)) * dt;
    Pl(t) = max(0, min(1, Pl(t))); % Ensure Pl stays within [0, 1]


    % Update Experienced Level (EL) for the next step
    EL(t) = delta_EL * EL_base(t) + (1 - delta_EL) * Pl(t-1);
end

% set maximum axis y - for graphs
minLimX = 0;
% set maximum axis y - for graphs
maxLimY = 1.2;


% Plotting
t = 1:numSteps;
figure;
plot(t, Pl, 'LineWidth', 2); % plotting results from Long-Term Performance
hold on;
plot(t, SLl ,'LineWidth', 2); % plotting results from Long-Term Stress Level
plot(time, SLs, 'LineWidth', 2); % plotting results from Short-Term Stress Level
plot(time, Ps, 'LineWidth', 2); % plotting results from Short-Term Performance
title('Level');
xlabel('Time');
legend('Long Performance', 'Long stress', 'Short Stress', 'Short Performance');
xlim([0 t(end)]);
ylim([minLimX maxLimY]);
grid on;

