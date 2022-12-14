clear;

Ninstance = 1000;    % # problem instances
Nod = 1000;   % # O-D pairs
Ntransit = 100;   % # transit station candidates

Nangles = 91;   % # angle slices for polar angle
THlist = linspace(0, 2*pi, Nangles);

Xori = zeros(Ninstance, Nod);  % X position of origin
Yori = zeros(Ninstance, Nod);  % Y position of origin
Xdes = zeros(Ninstance, Nod);  % X position of destination
Ydes = zeros(Ninstance, Nod);  % Y position of destination
Xtra = zeros(Ninstance, Ntransit);  % X position of P&R
Ytra = zeros(Ninstance, Ntransit);  % Y position of P&R
Nneighbor = zeros(Ninstance, 1);    % Number of resident neighborhood
Ncbd = zeros(Ninstance, 1); % Number of CBD

Vcars = cell(Ninstance, 1);  % Utility of car
Vtransits = cell(Ninstance, 1);  % Utility of P&R

for inst = 1:Ninstance
    Nneighbor(inst) = randi(6)+4; % Number of resident neighborhood: 5~10
    Rneighbor = rand(1, Nneighbor(inst))*4 + 6;   % Radius coordinate of neighborhood center in [6, 10]
    % THneighbor = rand(1, Nneighbor(inst))*2*pi;   % Polar angle of neighborhood center in [0, 2pi]
    THneighbor = THlist(randperm(Nangles-1, Nneighbor(inst)));
    Xneighbor = Rneighbor.*cos(THneighbor); % X coordinate of neighborhood center
    Yneighbor = Rneighbor.*sin(THneighbor); % Y coordinate of neighborhood center

    indNeighbor = randi(Nneighbor(inst), 1, Nod);  % Randomly assign neighborhood
    Xori(inst, :) = rand(1, Nod)*2 - 1 + Xneighbor(indNeighbor);    % X position of origin: [-1, 1] from neighborhood center
    Yori(inst, :) = rand(1, Nod)*2 - 1 + Yneighbor(indNeighbor);    % Y position of origin: [-1, 1] from neighborhood center

    Ncbd(inst) = randi(3)+2;  % Number of CBD: 3~5
    Rcbd = rand(1, Ncbd(inst)) + 1;   % Radius coordinate of CBD in [1, 2]
    % THcbd = rand(1, Ncbd(inst))*2*pi; % Polar angle of CBD in [0, 2pi]
    THcbd = THlist(randperm(Nangles-1, Ncbd(inst)));
    Xcbd = Rcbd.*cos(THcbd);    % X coordinate of CBD
    Ycbd = Rcbd.*sin(THcbd);    % Y coordinate of CBD

    indCBD = randi(Ncbd(inst), 1, Nod);
    Xdes(inst, :) = Xcbd(indCBD);   % X position of destination: a CBD
    Ydes(inst, :) = Ycbd(indCBD);   % Y position of destination: a CBD

    Rtra = rand(1, Ntransit)*2 + 5; % Radius coordinate of transit in [5, 7]
    % THtra = rand(1, Ntransit)*2*pi; % Polar angle of transit in [0, 2pi]
    THtra = THlist(randi(Nangles-1, 1, Ntransit));
    Xtra(inst, :) = Rtra.*cos(THtra);    % X position of transit
    Ytra(inst, :) = Rtra.*sin(THtra);    % X position of transit
    
    Vcars{inst} = -sqrt((Xdes(inst, :)'-Xori(inst, :)').^2 + (Ydes(inst, :)'-Yori(inst, :)').^2);   % Utility of car = direct distance of O-D : Nod-by-1
    Vtransits{inst} = -sqrt((Xtra(inst, :)-Xori(inst, :)').^2 + (Ytra(inst, :)-Yori(inst, :)').^2) - sqrt((Xdes(inst, :)'-Xtra(inst, :)).^2 + (Ydes(inst, :)'-Ytra(inst, :)).^2);   % Utility of transit = distance of O-Transit-D : Nod-by-Ntransit
end

save(['loc', num2str(Ninstance) , 'Instances.mat']);

fig = figure;
fig.Position(3:4) = [1050, 420];
subplot(1, 2, 2);
plot(Xori(inst, :), Yori(inst, :), 'g*');
hold on;
grid on;
plot(Xdes(inst, :), Ydes(inst, :), 'bd');
plot(Xtra(inst, :), Ytra(inst, :), 'ro');
for i=1:Nneighbor(inst)
    % text(Xneighbor(i)-0.5, Yneighbor(i), ['NH ', num2str(i)], 'color', 'green');
    text(Xneighbor(i)-0.5, Yneighbor(i), ['NH ', num2str(i)], 'color', 'black');
end
legend('Origin', 'Destination', 'Candidate P&R', 'location', 'southeast');
subplot(1, 2, 1);
plot(Xori(inst, 1:100), Yori(inst, 1:100), 'g*');
hold on;
grid on;
plot(Xdes(inst, 1:100), Ydes(inst, 1:100), 'bd');
plot(Xtra(inst, 1:30), Ytra(inst, 1:30), 'ro');
for i=1:Nneighbor(inst)
    text(Xneighbor(i)-0.5, Yneighbor(i), ['NH ', num2str(i)], 'color', 'black');
end
legend('Origin', 'Destination', 'Candidate P&R', 'location', 'southeast');
% saveas(gcf, 'flpEx2');
