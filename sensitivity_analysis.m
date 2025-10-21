% Sensitivity Analysis

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
  HR(t) = 1;
  SQ(t) = 0.222;
  MC(t) = 0.25;
  EL_base(t) = 1;
  ES(t) = 0.667;
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
SLl(1) = 0.5;
Pl(1) = 0.5;

% Define parameters (tune these as needed)

alpha_CL = 0.3; %el
w_CL1 = 0.3;    %mc
w_CL2 = 0.3;    %es
w_CL3 = 0.4;    %hr
beta_CL = 0.3;  %sq

w_SLs1 = 0.3;     %hr
w_SLs2 = 0.3;     %es
w_SLs3 = 0.4;     %mc
lamda_SLs = 0.3;  %sq
alpha_SLs = 0.4;  %el

alpha_PF = 0.5; %hr,pf,sq
w_PF1 = 0.5;    %hr
w_PF2 = 0.5;    %sq


gamma_RT = 0.5; %sq, sa
w_RT1 = 0.5;    %sq
w_RT2 = 0.5;    %sa
w_RT3 = 0.5;    %cl
w_RT4 = 0.5;    %pf

beta_SA = 0.5;  %sq,el,pl
w_SA1 = 0.3;    %sq
w_SA2 = 0.3;    %el
w_SA3 = 0.3;    %pl

delta_EL = 0.5; %pl

eta_Ps = 0.5; %sa,rt
w_Ps1 = 0.5;  %sa
w_Ps2 = 0.5;  %rt
w_Ps3 = 0.5;  %sl
w_Ps4 = 0.5;  %pf

% Sensitivity Analysis: Vary Sleep Quality (SQ)
SQ_values = linspace(0, 1, 5); % Vary SQ from 0 to 1 in 5 steps
numSQ = length(SQ_values);

% Initialize arrays to store the sensitivity analysis results
results_PF = zeros(numSQ, numSteps);
results_RT = zeros(numSQ, numSteps);
results_SA = zeros(numSQ, numSteps);
results_CL = zeros(numSQ, numSteps);
results_SLs = zeros(numSQ, numSteps);

% Perform sensitivity analysis
for i = 1:numSQ
    % Set the sleep quality for this iteration
    SQ = SQ_values(i) * ones(1, numSteps);

    % Reset initial conditions
    SLl(1) = 0.5;
    Pl(1) = 0.5;
    SLs(1) = 0;
    CL(1) = 0;
    EL(:) = 1;

% Compute initial condiitons
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

% Store the results for this SQ value
    results_PF(i, :) = PF;
    results_RT(i, :) = RT;
    results_SA(i, :) = SA;
    results_Ps(i, :) = Ps;
    results_CL(i, :) = CL;
    results_SLs(i, :) = SLs;
end

% Plotting the results
figure;
for i = 1:numSQ
    subplot(2, 3, 1);
    plot(time, results_PF(i, :), 'DisplayName', ['SQ = ' num2str(SQ_values(i))]);
    hold on;
    title('Physical Fatigue (PF)');
    xlabel('Time');
    ylabel('PF');
    legend;

    subplot(2, 3, 2);
    plot(time, results_RT(i, :), 'DisplayName', ['SQ = ' num2str(SQ_values(i))]);
    hold on;
    title('Reaction Time (RT)');
    xlabel('Time');
    ylabel('RT');
    legend;

    subplot(2, 3, 3);
    plot(time, results_SA(i, :), 'DisplayName', ['SQ = ' num2str(SQ_values(i))]);
    hold on;
    title('Situational Awareness (SA)');
    xlabel('Time');
    ylabel('SA');
    legend;

    subplot(2, 3, 4);
    plot(time, results_Ps(i, :), 'DisplayName', ['SQ = ' num2str(SQ_values(i))]);
    hold on;
    title('Cognitive Load (CL)');
    xlabel('Time');
    ylabel('CL');
    legend;

    subplot(2, 3, 5);
    plot(time, results_Ps(i, :), 'DisplayName', ['SQ = ' num2str(SQ_values(i))]);
    hold on;
    title('Short-Term Stress Level (SLs)');
    xlabel('Time');
    ylabel('SLs');
    legend;
end


