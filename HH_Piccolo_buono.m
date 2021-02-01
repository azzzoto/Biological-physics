function HH_Piccolo()

close all
clear all

% Define constants here
Vrest = -65.0 ; % resting potential, mV

V_Na = 50.0 ; % sodium Nernst potential, mV 
gNa_max = 120  ; % sodium max conductance, mS/cm2

V_K = -77.0; %potassium Nernst potential, mV
gK_max = 36; %potassium conductance, mS/cm2

V_L = -54.4; %leak Nernst potential
gL_max = 0.3; % leakage conductance, mS/cm2
Cm = 1; %membrane capacitance 1uF/cm2

% time stepping and stimulus related constants
deltaT = 0.001 ; % time step, millisec
tStart = -1.000 ; % start time, millisec
tEnd = 15.000 ; % end time, millisec
nStep = ceil((tEnd-tStart)/deltaT); % number of time steps
outputInterval = 20 ; % number of time steps between screen output
StimDur=2; % duration of the stimulus, millisec
IStim=100; % [uA/cm^2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial value of state variables
Vm = Vrest ; % membrane potential, mV
m = 0.05293 ; % initial value of m gate
h = 0.59612; %initial value for h gate
n = 0.31768; %initial value of n gate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage for variables to plot
plot_Vm = zeros(nStep,1) ;
plot_m = zeros(nStep, 1);
plot_h = zeros(nStep, 1);
plot_n = zeros(nStep, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print a heading for the screen output
disp('Hodgkin-Huxley model');
fprintf('\t Time \t Vm \t n \t\t m \t\t h\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the simulation
tNow = tStart ;
for iStep = 1:nStep
% Compute ion currents & stimulus current
JNa = gNa_max*m*m*m*h*(Vm-V_Na) ;
JK = gK_max*n*n*n*n*(Vm - V_K);
JL = gL_max*(Vm-V_L);

if( 0<=tNow && tNow<StimDur ) % apply stimulus current at t=0
Jm = IStim ;
else
Jm = 0 ;
end

% Compute gates' opening and closing rates
alpha_m = 0.1.*(Vm+40)/(1-exp(-(Vm+40)/10));
beta_m = 0.108*exp(-Vm/18);
m = m + deltaT*(alpha_m*(1-m)-beta_m*m);

alpha_n = (0.01*(55+Vm))/(-exp(-(55+Vm)/10)+1);
beta_n = 0.055*exp(-Vm/80);
n = n + deltaT*(alpha_n*(1-n) - beta_n*n);

alpha_h = 0.0027*exp(-Vm/20);
beta_h = 1/(exp(-(35+Vm)/10)+1);
h = h + deltaT*(alpha_h*(1-h) - beta_h*n);

% Compute change in state variables
deltaVm = -deltaT * (JNa + JK + JL - Jm)/Cm ;

% Record/display state variables & other values
plot_Vm(iStep) = Vm ;
plot_m(iStep) = m;
plot_n(iStep) = n;
plot_h(iStep) = h;
plot_JNa(iStep)= JNa;
plot_JK(iStep)=JK;
plot_time(iStep) = tNow;
if mod(iStep,outputInterval) == 0
fprintf('%8.2f %7.3f %7.5f %7.5f %7.5f\n', tNow, Vm, n, m, h) ;
end % if mod(tNow)

% Update state variables
Vm = Vm + deltaVm ; % new Vm = current Vm + change in Vm
tNow = tStart + iStep*deltaT ; % increment the time
end % for iStep

%DataNa(:,1) = plot_time;
%DataNa(:,2) = plot_Vm;
% Plot the gates, probabilities, currents, Vm, etc

plot(plot_time, plot_Vm); grid on ;
xlabel('time (ms)')
ylabel('Transmembrane Voltage(Vm), mV')
title('Time course of Action Potential')
opengl software;

figure();
plot(plot_time,gK_max*(plot_n).^4)
hold on;
plot(plot_time,gNa_max.*plot_h.*(plot_m).^3)
grid on;
title('Na and K conductances vrs time')
xlabel('time (ms)')
ylabel('Conductance (mS/cm^2)')
legend({'g_K','g_{Na}'},'Location', 'northeast')

figure();
plot(plot_time, plot_m); grid on ;
hold on
plot(plot_time, plot_h);
hold on
plot(plot_time, plot_n);
legend('m', 'h', 'n')
xlabel('time (ms)')
ylabel('gating activation/ inactivation probability')
title('gating probabilities (m, h, n) vrs time')

end


