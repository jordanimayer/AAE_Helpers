% AAE 334 - HW 09
% Problem 1
%
% Jordan Mayer

% Declar variables
L = 1;      % length of nozzle (ends up cancelling out)
At = 1;     % throat area (ends up cancelling out)
Ae = 4.235;   % ratio of exit area to throat area
x = linspace(-L/2, L/2, 1001);
A = zeros(1,1001);  % area along length
gamma = 1.4;    % specific heat ratio for air
for k=1:1001
    A(k) = At + (Ae-At)*(2*x(k)/L)^2;
end
pbp01 = [0.2812, 0.30, 0.50, 0.70, 0.90, 0.9867408];
    % ratios of back pressure to stagnation pressure
d = zeros(1,6); % array to hold d/L values

% choose pressure ratio
for k=1
    
    % CALCULATE d
    
    p01 = 1;    % stagnation pressure before throat
                % (will cancel out)
    pe = pbp01(k);  % ratio of exit pressure to
                    % stagnation pressure
    
    % calculate exit Mach number (break up equation)
    Me = 2*(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1));
    Me = Me * ((1/pe)*(At/Ae))^2;
    Me = Me + 1;
    Me = Me^(1/2);
    Me = Me + (-1);
    Me = Me * (1/(gamma-1));
    Me = sqrt(Me);  % THIS is exit Mach number
    
    % find stagnation pressure downstream of shock
    p02 = pe * (1 + (gamma-1)/2 * Me^2)^(gamma/(gamma-1));
    
    % find pressure ratio, p02/p01
    p02p01 = p02/p01;
    
    % use shock jump conditions and iteration to find Mach number
    % upstream of the shock: M1 = f(p02/p01)
    M1 = 0;
    for M=linspace(0, 10, 10000);
        pp = (gamma+1)*M^2/((gamma-1)*M^2+2);
        pp = pp^(gamma/(gamma-1));
        pp = pp*(1+2*gamma/(gamma+1)*(M^2-1))^(-1/(gamma-1));
        % we now have a value for p02/p01 if M1 == M
        if (abs(pp-p02p01) < 0.001 && pp < p02p01)   % if our pressure ratio is correct
            M1 = M; % we found M1
            break;
        end
    end
    if (M1 == 0)    % check if M1 was found
        fprintf('Could not find M1\n');
    end
    
    % use isentropic relations to find area ratio at
    % the shock: As/A1star = f(M1)
    AsA1star = (2/(gamma+1))*(1 + (gamma-1)/2 * M1^2);
    AsA1star = AsA1star^((gamma+1)/(gamma-1));
    AsA1star = AsA1star/M1^2;
    AsA1star = sqrt(AsA1star);  % THIS is AsA1star
    As = AsA1star;  % if At = Astar = 1
    
    % solve for position of the shock
    d(k) = 0.5*sqrt((As/At - 1)/(Ae/At - 1));  % technically d/L
    
    % CALCULATE PLOT DATA
    
    pp = zeros(1,1001); % p/p01
    TT = zeros(1,1001); % T/T0
    AA = zeros(1,1001); % A(x)/At
    Mx = zeros(1,1001);
    
    for i=1:1001 % -1/2 <= x/L <= 0
        if (x(i) <= 0)
            % choose x and find A(x)
            AAstarsq = (A(i)/At)^2;

            % iterate to find M=f(A/A1star), subsonic, A1star=At
            Mtrue = 0;
            for M=linspace(0.000001, 1, 10^5);
                if (x(i) == 0)  % M = 0 at throat, don't waste time with iteration
                    Mtrue = 1;
                    break;
                end
                AAsq = (2/(gamma+1)*(1+(gamma-1)/2*M^2))^((gamma+1)/(gamma-1));
                AAsq = AAsq/M^2;
                % we now have a value for (A/Astar)^2 if Mtrue == M
                if (abs(AAsq - AAstarsq) < 0.01 && AAsq < AAstarsq)   
                    % if our area ratio is correct AND we are subsonic
                    Mtrue = M; % we found M
                    break;
                end
            end
            
            if (Mtrue == 0)    % check if M was found
                fprintf('Could not find Mtrue for first case\n');
            end
            M = Mtrue;

            % find p/p01, T/T01 from isentropic relations
            TT(i) = 1+(gamma-1)/2*M^2;
            TT(i) = 1/TT(i);
            pp(i) = TT(i)^(gamma/(gamma-1));
            
        
        elseif (x(i) <= d(k))
            AAstarsq = (A(i)/At)^2;

            % iterate to find M=f(A/A1star), subsonic, A1star=At
            Mtrue = 0;
            
            for M=linspace(1, 4, 10^5);
                AAsq = (2/(gamma+1)*(1+(gamma-1)/2*M^2))^((gamma+1)/(gamma-1));
                AAsq = AAsq/M^2;
                % we now have a value for (A/Astar)^2 if Mtrue == M
                if (abs(AAsq - AAstarsq) < 0.01 && AAsq > AAstarsq)   
                    % if our area ratio is correct AND we are supersonic
                    Mtrue = M; % we found M
                    break;
                end
            end
            
            
            if (Mtrue == 0)    % check if M was found
                fprintf('Could not find Mtrue for second case \n');
            end
            M = Mtrue;
            
            % find p/p01, T/T01 from isentropic relations
            TT(i) = 1+(gamma-1)/2*M^2;
            TT(i) = 1/TT(i);
            pp(i) = TT(i)^(gamma/(gamma-1));
           
        
        else
            % choose x and find A(x)
            % find A/A2star = A/A1star * A1star/A2star = A/At * p02/p01
            AAstarsq = (A(i)/At*p02/p01)^2;
            
            % iterate to find M=f(A/A2star), subsonic branch
            Mtrue = 0;
            
            for M=linspace(0.000001, 0.999999, 10^5);
                AAsq = (2/(gamma+1)*(1+(gamma-1)/2*M^2))^((gamma+1)/(gamma-1));
                AAsq = AAsq/M^2;
                % we now have a value for (A/Astar)^2 if Mtrue == M
                if (abs(AAsq - AAstarsq) < 0.01 && AAsq < AAstarsq)   
                    % if our area ratio is correct AND we are subsonic
                    Mtrue = M; % we found M
                    break;
                end
            end
            
             
            if (Mtrue == 0)    % check if M was found
                fprintf('Could not find Mtrue for third case\n');
            end
            M = Mtrue;
            
            % find p/p01=p/p02*p02/p01, T/T01 from isentropic relations
            TT(i) = 1+(gamma-1)/2*M^2;
            TT(i) = 1/TT(i);
            pp02 = TT(i)^(gamma/(gamma-1));
            
            pp(i) = pp02 * p02/p01;
        end
        AA(i) = A(i)/At;
        Mx(i) = M;
    end
    
    % CREATE PLOTS
    
    ttl = sprintf('Wind Tunnel Analysis for varying p_{b}/p_{01}');
    % plot Mach number
    figure(2); plot(x, Mx); title(ttl);
    xlabel('x'); ylabel('M'); hold on;
    
    % plot temperature ratio
    figure(3); plot(x, TT); title(ttl);
    xlabel('x'); ylabel('T/T_{0}'); hold on;
    
    % plot pressure ratio
    figure(4); plot(x, pp); title(ttl);
    xlabel('x'); ylabel('p/p_{01}'); hold on;
end

% plot area ratio profile
figure(1); plot(x, AA); title(ttl); grid on;
xlabel('x'); ylabel('A(x)/A_{t}'); ylim([0 4.5]);

% format plots
figure(2);grid on; legend('0.2812', '0.30', '0.50', '0.70', '0.90', '0.9867408');
figure(3);grid on; legend('0.2812', '0.30', '0.50', '0.70', '0.90', '0.9867408');
figure(4);grid on; legend('0.2812', '0.30', '0.50', '0.70', '0.90', '0.9867408');
