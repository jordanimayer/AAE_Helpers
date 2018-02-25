%this code was written on a caffine and sleep deprivation induced rampage

%constants and starting information (units in metric, m,N,kg etc)
L=13.22;    % length of aircraft, m     
span=12.8;  % wingspan, m
H=4.35;     % height of aircraft, m
W_max= 4950;    % maximum takeoff weight of aircraft, kgf
W_max = W_max * 9.81;   % convert W_max to Newtons
LF_pos = 3.8;   % maximum positive load factor
LF_neg = -1.7;  % most negative load factor
L_panel= 0.16+0.12+0.75*2+sqrt(0.06^2 + 0.35^2)*2+sqrt(0.08^2+0.40^2)*2;
    % sum of all panel lengths
F_pos=W_max*LF_pos*0.5; % each wing takes 1/2 of force
F_neg=W_max*LF_neg*0.5; % each wing takes 1/2 of force
%a = 0.75;
b = 2;
K=8.98;
A_closed=0;

%material Properties: [2024-Aluminum, Ti-6Al-2Sn-4Zr-2Mo-Titanium]
UTS = [469, 1010];  % ultimate stress, MPa
YS = [324, 990];    % yield stress, MPa
E = [73.1, 120];    % modulus of elasticity, GPa
G = [28, 45.5];     % shear modulus, GPa
den = [2780, 4540]; % density, kg/m^3

% convert to consistent units
UTS = UTS*10^6; % MPa --> Pa
YS = YS*10^6;   % MPa --> Pa
E = E*10^9;     % GPa --> Pa
G = G*10^9;     % GPa --> Pa

% SM = [28,45.5]; %GPa shear modulus
% den = [2780, 4540]; %kg/m^3

%geometry calculations
y_bar=[-0.08 -0.11 -0.11 -0.06 0.06 0.11 0.11 0.08];

%Factors to change
A_st = sqrt(abs(F_pos*y_bar*1.5/(YS(2)*sum(y_bar.^2))));     %test values note see image in folder as to array position for each element
%t = 0.001;
%num_string = 8;
%initial calculations
Izz = sum((y_bar.^2).*A_st);


%positve force calcultion
%shear flows
q = zeros(1,9);
dp = zeros(1, 9);
q(1)=0;
for i=2:9
    dp(i)=F_pos*y_bar(i-1)*(0.25*span)*A_st(i-1)/Izz;     %delta p in terms of Ai and Izz
    q(i)=dp(i)/(0.25*span)+q(i-1);              %q in terms of Ai
end
%Calculating critical tau
t = zeros(1,8);
for i=1:8
    t(i)= 1.5*sqrt(12*YS(1)*b^2/(pi^2*K*E(1)));   %check K, this is to show instability if Tau_top exceeds Tcrit the structure is unstable
end

%moment calculation
M0=2*q*A_closed;                             %I'm gonna break my rusty cage and run
%compression(buckling) on top tension(yielding) on bottom
%Top Buckling
Pcrit=pi^2*E(2)*Izz/(0.25*span)^2;            %in terms of Izz, Izz cancels with Izz from dpi
panic = zeros(1,5);
for i=2:5
    panic(i-1)= Pcrit > dp(i);                      %if this is not one long array of ones you're in trouble
    if panic(i-1)==0
        print('stringer failure')
    end
end
%bottom Yielding
for i=6:9
    panic(i-1)= YS(1)/1.5 > q(i)/t(i-1);                      %if this is not one long array of ones you're in trouble
    if panic(i-1)==0
        print('skin failure')
    end
end

%wing weight (this is to be minimized)
fprintf('y_bar = ');
disp(y_bar);
fprintf('dp = ');
disp(dp);
fprintf('A_st = ');
disp(A_st);
fprintf('q = ');
disp(q);
fprintf('t = ');
disp(t);
%wing weight, kg (this is to be minimized)
W_wing = den(2)*(sum(A_st.*span/2))+den(1)*sum(L_panel*t*span/2)

%negative force
A_st = sqrt(abs(F_neg*y_bar*1.5/(YS(2)*sum(y_bar.^2))));     %test values note see image in folder as to array position for each element
%t = 0.001;
%num_string = 8;
%initial calculations
Izz = sum((y_bar.^2).*A_st);
q(1)=0;
for i=2:9
    dp(i)=F_neg*y_bar(i-1)*(0.25*span)*A_st(i-1)/Izz;     %delta p in terms of Ai and Izz
    q(i)=dp(i)/(0.25*span)+q(i-1);              %q in terms of Ai
end
%Calculating critical tau 
for i=1:8
    t(i)= 1.5*sqrt(12*YS(1)*b^2/(pi^2*K*E(1)));   %check K, this is to show instability if Tau_top exceeds Tcrit the structure is unstable
end
%bottom Buckling
Pcrit=pi^2*E(2)*Izz/(0.25*span)^2;            %in terms of Izz, Izz cancels with Izz from dpi
for i=6:9
    panic(i-1)= Pcrit > dp(i);                      %if this is not one long array of ones you're in trouble
    if panic(i-1)==0
        fprintf('stringer failure');
    end
end
%top Yielding
for i=2:5
    panic(i-1)= YS(1)/1.5 > q(i)/t(i-1);                      %if this is not one long array of ones you're in trouble
    if panic(i-1)==0
        fprintf('skin failure');
    end
end

% print variables (order: counterclockwise from bottom-left)
fprintf('y_bar = ');
disp(y_bar);
fprintf('dp = ');
disp(dp);
fprintf('A_st = ');
disp(A_st);
fprintf('q = ');
disp(q);
fprintf('t = ');
disp(t);
%wing weight, kg (this is to be minimized)
W_wing = den(2)*(sum(A_st.*span/2))+den(1)*sum(L_panel*t*span/2)
