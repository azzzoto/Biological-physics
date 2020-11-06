function GHK_Piccolo_Giovanni()

    clear;

    %Data, system and simulation parameters
    V=5000e-18; %m^3
    R=(3*V/(4*pi))^(1/3); %m
    A=4*pi*R^2; %m^2

    C_in=[12 155]; %Intramembrane concentration of [Na+ K+] in mM
    C_out=[145 4]; %Extramembrane concentration of [Na+ K+] in mM
    P=[0.2e-3 10e-3]; %Permeability of [Na+, K+] in mm/s
    capacitance=1e-2; %F/m^2

    %the valence numbers are [1 1] for [Na+ K+], not included in the
    %computation since they are superfluos

    Delta_V=-65e-3; %V
    q=1.6e-19; %fundamental charge in C
    N_A = 6.022e23; %avogadro number

    k_B=1.380649e-23; %Boltzmann constant[J/K]
    T=25; %°C, Standard Ambient Pressure
    kT=k_B*(273.15+T);

    dt = 1e-8; %time steps
    total_time = 0.025;
    N_steps = floor( total_time / dt );
    potential = zeros(N_steps,1);
    potential(1) = Delta_V;
    %flux = zeros(1,2);
    C = zeros(N_steps,2); %concentration
    C(1,:) = C_in;
    %current = zeros(N_steps,2);

    %Calculations of initial fluxes using factors to reduce the computational time
    exponential=exp(q*Delta_V/kT);
    prefactor= P.*q*Delta_V/(kT);
    initial_fluxes= prefactor.*((C_in.*exponential)-C_out)./(1-exponential);
    disp('Initial fluxes of:')
    disp('      Na+        K+');
    disp(initial_fluxes);

    %Factors to speed up the simulation cycle
    alpha_simulation = P.*q/kT;
    beta_simulation = q/kT;
    gamma_simulation = A * dt / V;
    %delta_simulation =  q*N_A*A;

    n=input('Enter 0 for the basic case; enter -1 for the case with the negative charge:  ');
        switch n
        case 0

        %GHK simulation
        %tic
        for i = 2:N_steps

            flux = alpha_simulation.*potential(i-1).*(C(i-1,:).*exp(beta_simulation*potential(i-1))-C_out)./(1-exp(beta_simulation*potential(i-1)));
            C(i,:) = C(i-1,:) + flux * gamma_simulation ;
            %current(i,:) = delta_simulation*flux ;
            potential(i) = potential(i-1) + q*N_A*dt*sum(flux)/capacitance;

        end      
        %toc

        disp('Final fluxes of: ');
        disp('      Na+        K+');
        disp(flux);

        disp('Final concentrations in [mM] of: ');
        disp('      Na+        K+');
        disp(C(N_steps,:));

        disp('Equilibrium potential [V]: ');
        disp(potential(N_steps));

        %total_current=sum(current,2); %I = I_K + I_Na

        t=0:dt:dt*(N_steps-1); %X axe

        %double plot (t,V) and (t,C)
        y1 = potential*1e3;
        y2 = C;
        tiledlayout(2,1);

        % Top plot
        ax1 = nexttile;
        plot(ax1,t,y1,'r');
        title(ax1,'Intracellular potential');
        xlabel('Time [s]');
        ylabel(ax1,'Potential [mV]');
        legend({'\Delta V'},'Location','southeast');

        % Bottom plot
        ax2 = nexttile;
        plot(ax2,t,y2);
        title(ax2,'Intracellullar concentrations');
        xlabel('Time [s]');
        ylabel(ax2,'Concentrations (mM)');
        legend({'Na+','K+'},'Location','southeast');

    case -1

        Q_pm = A*Delta_V*capacitance - sum( C_in ) * N_A * q * V;
        disp("Initial charge [pC]: ");
        disp(A*Delta_V*capacitance*1e12);

        additional_charge_factor=Q_pm/(A*capacitance);

        %GHK simulation
        %tic
        for i = 2:N_steps

            flux = alpha_simulation.*potential(i-1).*(C(i-1,:).*exp(beta_simulation*potential(i-1))-C_out)./(1-exp(beta_simulation*potential(i-1)));
            C(i,:) = C(i-1,:) + flux * gamma_simulation ;
            %current(i,:) = delta_simulation*flux ;
            potential(i) = sum(C(i,:))*q*N_A*V/(A*capacitance)+additional_charge_factor;

        end      
        %toc

        disp('Final fluxes of: ');
        disp('      Na+        K+');
        disp(flux);

        disp('Final concentrations in [mM] of: ');
        disp('      Na+        K+');
        disp(C(N_steps,:));

        disp("Final Charge [pC]:");
        disp((sum(C(N_steps,:))*q*N_A*V+Q_pm)*1e12 )
        
        disp('Equilibrium potential [V]: ');
        disp(potential(N_steps));

        %total_current=sum(current,2); %I = I_K + I_Na

        t=0:dt:dt*(N_steps-1); %X axe

        %double plot (t,V) and (t,C)
        y1 = potential*1e3;
        y2 = C;
        tiledlayout(2,1);

        % Top plot
        ax1 = nexttile;
        plot(ax1,t,y1,'r');
        title(ax1,'Intracellular potential');
        xlabel('Time [s]');
        ylabel(ax1,'Potential [mV]');
        legend({'\Delta V'},'Location','southeast');

        % Bottom plot
        ax2 = nexttile;
        plot(ax2,t,y2);
        title(ax2,'Intracellullar concentrations');
        xlabel('Time [s]');
        ylabel(ax2,'Concentrations (mM)');
        legend({'Na+','K+'},'Location','southeast');

        figure();
        plot(t,(sum(C,2)*q*N_A*V+Q_pm)*1e12,'m');
        title("Charge inside the cell");
        xlabel("Time [s]");
        ylabel("Charge [pC]");
        legend({'Q'},'Location','southeast');

    otherwise
        disp('Input number not recognized. Terminating process...');    
end


