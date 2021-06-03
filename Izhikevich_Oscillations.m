clc; clear all; close all;

% The spike neuron code is from Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons    Inhibitory neurons

Ne=800;                 Ni=200; % neuron number (Ne = excitstory, Ni = inibitory
re=rand(Ne,1);          ri=rand(Ni,1); % row vector, excitateoy (re) & row vector inhibitory (ri)
a=[0.02*ones(Ne,1);     0.02+0.08*ri]; % time scale of recovery variable
b=[0.2*ones(Ne,1);      0.25-0.05*ri]; % sensitivety of the recovery variables
c=[-65+15*re.^2;        -65*ones(Ni,1)]; % after spike reset value
d=[8-6*re.^2;           2*ones(Ni,1)]; % after spike rest of recovery variable
S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)]; % excitatory and inhibitory weights

for i = 801:1000;
    for j = 801:1000;
    S(i,j) = 0;
end
end

v=-65*ones(Ne+Ni,1);    % Initial values of v
u=b.*v;                 % Initial values of u
firings=[];             % spike timings

for t=1:1000            % simulation of 1000 ms
  I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
  fired=find(v>=30);    % indices of spikes
  
  firings=[firings; t+0*fired,fired];
  v(fired)=c(fired);
  u(fired)=u(fired)+d(fired);

  I=I+sum(S(:,fired),2);
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
  u=u+a.*(b.*v-u);                 % stability
end;

x1 = firings(:,2);

inhib_i = [];

% runs through neuron number and choose > Ni - 1 (801 - 1000) for
% inhibitory neurons
for i = 1:length(firings);
    
    if x1(i) > Ni-1;
        inhib_i(i,2) = x1(i);
        inhib_i(i,1) = firings(i,1);
    end
end

% Plot Oscillation Results: %

inhib_no_1 = inhib_i(:,1);
inhib_no_2 = inhib_i(:,2);

%removes unwanted zeros
inhib_no_1(inhib_no_1==0) = [];
inhib_no_2(inhib_no_2==0) = [];

inhibs = [];

inhibs(:,1) = inhib_no_1;
inhibs(:,2) = inhib_no_2; % n x 2 inhibitory matrix

plot(inhibs(:,1), inhibs(:,2),'.r'); % plot inhibitory matrix

hold on

x2 = firings(:,2);

excit_e = [];

% runs through neuron number and choose < Ne (1 - 799) for
% excitatory neurons neurons
for i = 1:length(firings);
    
    if x2(i) < Ne;
        
        excit_e(i,2) = x2(i);
        excit_e(i,1) = firings(i,1);
        
    end
end

excit_no_1 = excit_e(:,1);
excit_no_2 = excit_e(:,2);

%removes unwanted zeros
excit_no_1(excit_no_1==0) = [];
excit_no_2(excit_no_2==0) = [];

excits = [];

excits(:,1) = excit_no_1;
excits(:,2) = excit_no_2; % n x 2 excitatory matrix

plot(excits(:,1), excits(:,2),'.b'); % plot excitatory oscillations

title('Neuron(s) Oscillations | EE, EI, IE | FS');
ylabel('Neuron Number'); xlabel('Time (mS)');
[~,obj_def] = legend({'Inhibitory Neurons', 'Excitatory Neurons'}, 'location', 'NorthEast', 'Fontsize', 10);
%choose required font size
objhl = findobj(obj_def, 'type', 'line'); %legend objects defined
set(objhl, 'Markersize', 30); %define si\e of marker