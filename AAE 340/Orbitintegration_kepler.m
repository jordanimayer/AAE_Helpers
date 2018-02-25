% Keplerian orbit integration

% integration initialization
dt=1.0;
t0=0.0;
tfinal=2*86400.0;
tnext=t0+dt;
options = odeset('RelTol',1E-12,'AbsTol',1E-12,'InitialStep',dt,'MaxStep',dt);
c=1;

% allocation
limit=tfinal/dt;
state=zeros(limit,6);
epoch=zeros(limit,1);



% orbit initialization
x0=[-28685970.39591214;28907860.51670547;6929070.112774;-2410.68202332;-1998.99410321;191.88163377]; % initial position and velocity
state(1,:)=x0; % save initial state
epoch(1)=t0; % save initial epoch


%  orbit integration
while tnext <=tfinal
    c=c+1; % counter
    [tn,xn] = ode45('orbitpde',[t0 tnext],x0,options);
    %size(stateout)
    t0=tn(end); 
    tnext=tnext+dt; % move to next time step
    x0=xn(end,1:6);
    epoch(c)=t0; % save epoch
    state(c,:)=x0; % save current state (position and velocity)
end

% make a plot
d=1:c;
figure(1)
plot3(state(d,1),state(d,2),state(d,3))
xlabel('meter')
ylabel('meter')
zlabel('meter')

