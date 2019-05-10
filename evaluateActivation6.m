function [axons, Ve_out, Vm_out] = evaluateActivation6(amp1,amp2,ipi,pw,nElectrodes,UIAxes,gauge) 


L = 2.5E-4;              % Nodal gap (cm)
dX = 0.1;                % Delta X (cm)
rhoI = 0.055;            % Internal resistance (kOhmcm)
rhoE = 0.3;              % External resistance (kOhmcm)

Temp = 273+20.15;        % Temperature - K
R = 8.31451;             % Ideal gas constant - J/(mol K)
F = 96485.3;             % Faraday's constant - C/mol

% Channel parameters
pK = 1;                  % K permeability
pNa = 0.03;              % Na permeability
pCl = 0.1;               % Cl permeability
KIn = 400;               % K concentration inside
NaIn = 50;               % Na concentration inside
ClIn = 40;               % Cl K concentration inside
KOut = 20;               % K concentration outside
NaOut = 460;             % Na concentration outside
ClOut = 540;             % Cl concentration outside

cm = 1;                  % Capacitance - uF/cm2
EL = -54.4;              % Leak reversal potential - mV
gL = 0.3 * 10;           % Leak conductance - mS/cm2
gK = 36 * 10;            % K conductance - mS/cm2
gNa = 120 * 10;          % Na conductance - mS/cm2
T = 2;                   % Total time - ms
dT = 0.001;              % delta T
t = 0:dT:T;              % Time course
nodes = 10;              % Number of nodes

nerveR = 1;              % Radius of the nerve (mm)
elDist = 0.05;           % Electrode distance (cm)

T_us = length(t);

Vr = ones(nodes,1) .* (R * Temp / F * log((pK * KOut + pNa * NaOut + pCl * ClIn) / (pK * KIn + pNa * NaIn + pCl * ClOut)) * 1000);

ENa = -R * Temp / F * log(NaIn / NaOut) * 1000;      % Na reversal potential (mV)
EK = -R * Temp / F * log(KIn / KOut) * 1000;         % K reversal potential (mV)

margin = 0.05;                      % Epineurium thickness
totalAxons = 0;
resolution = 0.05;



Ist = zeros(nElectrodes, size(t,2));    % Stimulation current
for i=1:nElectrodes
    Ist(i,500 : 500 + pw(i)) = amp1(i);
    Ist(i,(500 + pw(i) + ipi(i)) : (500 + 2 * pw(i) + ipi(i))) = amp2(i);
end



% axonDist = zeros(nElectrodes,1);
theta = 0 : 2 * pi / nElectrodes : 2 * pi - 0.01;

% Calculate total no of axons to be evaluated
for y =  -nerveR + margin : resolution : nerveR - margin
     for x = -nerveR + margin : resolution : nerveR - margin
         if (x^2 + y^2) + margin > nerveR^2
            continue;
         else
            totalAxons = totalAxons + 1;
         end
     end
end
Ve_out = zeros(totalAxons, T_us);   % To return the external voltages of each axon across all time points evaluated
Vm_out = zeros(totalAxons, T_us);   % To return the membrane voltages of each axon across all time points evaluated
axons = zeros(totalAxons,4);  

count = 0;                          % For measuring the percentage progress

%%
for y =  -nerveR + margin : resolution : nerveR - margin
     for x = -nerveR + margin : resolution : nerveR - margin
        if (x^2 + y^2) + margin > nerveR^2
            continue;
        end
        
        % d = 10E-4;               % Fibre diameter (cm)
        d = normrnd(1.47E-4,0.07E-4);               % Fibre diameter (cm)
        count = count + 1;   
        gauge.Value = count / totalAxons * 100;
        % TextArea.Value = [num2str(count / totalPoints * 100) ' %'];
        drawnow;
        
        % Distance of current axon from each of the stimulation electrodes
        axonDist = (((nerveR / 10 + elDist) * sin(theta) - x / 10) .^ 2 + ((nerveR / 10 + elDist) * cos(theta) - y / 10) .^ 2) .^ 0.5;

        d2Vedx2 = zeros(nodes, size(t,2)); % Activating function
        d2Vndx2 = zeros(nodes, size(t,2)); % Reduced voltage (2nd spatial derivative)
        dVndt = zeros(nodes, size(t,2));  % dVn / dt
        Vn = zeros(nodes, size(t,2));     % Vn (reduced voltage) = Vi - Ve - Vr
        Ve = zeros(nodes, size(t,2));     % Vn (reduced voltage) = Vi - Ve - Vr
        
        v = zeros(nodes, size(t,2));      % Membrane volrage

        m = zeros(nodes, size(t,2));      % Parameter m
        n = zeros(nodes, size(t,2));      % Parameter n
        h = zeros(nodes, size(t,2));      % Parameter h
        gNaI = zeros(nodes, size(t,2));   % Instanstaneous Na conductance
        gKI = zeros(nodes, size(t,2));    % Instanstaneous K conductance

        IK = zeros(nodes, size(t,2));     % K current
        INa = zeros(nodes, size(t,2));    % Na current
        IL = zeros(nodes, size(t,2));     % Leak current

        n(:, 1) = alphan(Vr) ./ (alphan(Vr) + betan(Vr));     % Initial value of n (0.3)
        m(:, 1) = alpham(Vr) ./ (alpham(Vr) + betam(Vr));     % Initial value of m (0.005)
        h(:, 1) = alphah(Vr) ./ (alphah(Vr) + betah(Vr));     % Initial value of h
        v(:, 1) = Vr;                                         % Membrane is at rest initially
        Vn(:,1) = 0;

        %%
        for i=1:length(t)

            % Activating function
            for nodeInd = 1:nodes
                z = abs(nodeInd - nodes / 2) * dX;
                Ve(nodeInd, i) = sum((rhoE / (4 * pi)) * (axonDist' .^ 2 + z ^ 2) .^ (-1 / 2) .* Ist(:, i));
                d2Vedx2(nodeInd, i) = sum(rhoE .* ((z ^ 2 + axonDist' .^ 2) .^ (-5 / 2)) .* (2 * z ^ 2 - axonDist' .^ 2) .* Ist(:, i) / (4 * pi));
            end
            v(:, i) = Vn(:, i) + Vr(1,1);

            % Conductances
            gNaI(:, i) = gNa .* (m(:, i) .^ 3) .* h(:, i);    % Instantaneous Na conductance
            gKI(:, i) = gK .* (n(:, i) .^ 4);                 % Instantaneous K conductance

            % Currents
            INa(:, i) = gNaI(:, i) .* (v(:, i) - ENa);        % Na current 
            IL(:, i) = gL .* (v(:, i) - EL);                  % Leak current
            IK(:, i) = gKI(:, i) .* (v(:, i) - EK);           % K current
            Iionic = INa(:,i) + IK(:,i) + IL(:,i); 

            for nodeInd = 2:nodes-1
                d2Vndx2(nodeInd,i) = (Vn(nodeInd - 1,i) - 2 * Vn(nodeInd,i) + Vn(nodeInd + 1,i)) / (dX ^ 2);   
            end

            % Calculate increments
            dVndt(:, i) = ( ( d * dX / (4 * rhoI * L) ) .* (d2Vndx2(:, i) + d2Vedx2(:, i)) - Iionic(:)) / cm;       
            dndt = alphan(v(:, i)) .* (1 - n(:, i)) - betan(v(:, i)) .* n(:, i);
            dmdt = alpham(v(:, i)) .* (1 - m(:, i)) - betam(v(:, i)) .* m(:, i);
            dhdt = alphah(v(:, i)) .* (1 - h(:, i)) - betah(v(:, i)) .* h(:, i);

            % Update all variables
            if i < length(t)
                Vn(:, i+1) = Vn(:, i) + dVndt(:, i) .* dT;    
                Vn(end, i+1) = Vn(end-1, i+1);
                Vn(1, i+1) = Vn(2, i+1);
                m(:, i+1) = m(:, i) + dmdt .* dT;
                n(:, i+1) = n(:, i) + dndt .* dT;
                h(:, i+1) = h(:, i) + dhdt .* dT;
            end
        end
        
        Vm_out(count, :) = v(nodes/2, :);
        Ve_out(count, :) = Ve(nodes/2, :);
        
        % Store the x and y coordinates of the current axon being evaluated
        axons(count, 1) = x;
        axons(count, 2) = y; 
        
        % store the diameter of the axon
        axons(count, 4) = d;  
        
        % Check if the 2nd node from stim node has triggered an AP
        if sum(v(nodes / 2 + 2,:) > 20) > 0
            axons(count, 3) = 1;           
        end        
     end
end

for i=1:length(theta)
   hold(UIAxes,'on');
   scatter(UIAxes,(nerveR + elDist) * sin(theta(i)),(nerveR + elDist) * cos(theta(i)),[],[1 0 0],'filled'); 
end

