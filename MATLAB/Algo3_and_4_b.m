clc
clear
close all

rng_seed_val = 9;
tStart = tic;
% theta_arr = 0.1:0.1:180; % The theta_th values we are checking for
% theta_arr = [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]; 
% theta_arr = linspace(5, 180, 130);
theta_arr = [40, 60, 80, 90, 100];
% theta_arr = 60;
len_the = length(theta_arr); 

% For each theta we will simulate the network for 10k times
itr = 10000;


N_arr = [30, 40, 50, 60, 70, 80];
% N_arr = 80;
% N_arr = [30, 40, 50];
len_N_arr = length(N_arr);

acc_N30 = zeros(1, len_the);
acc_N40 = zeros(1, len_the);
acc_N50 = zeros(1, len_the);
acc_N60 = zeros(1, len_the);
acc_N70 = zeros(1, len_the);
acc_N80 = zeros(1, len_the);

time_N30_A3 = zeros(1, len_the);
time_N40_A3 = zeros(1, len_the);
time_N50_A3 = zeros(1, len_the);
time_N60_A3 = zeros(1, len_the);
time_N70_A3 = zeros(1, len_the);
time_N80_A3 = zeros(1, len_the);

time_N30_A4 = zeros(1, len_the);
time_N40_A4 = zeros(1, len_the);
time_N50_A4 = zeros(1, len_the);
time_N60_A4 = zeros(1, len_the);
time_N70_A4 = zeros(1, len_the);
time_N80_A4 = zeros(1, len_the);



for N_ind = 1:len_N_arr


for theta_ind = 1:len_the
% disp("The theta index is:")
N = N_arr(N_ind);
disp('N is:')
disp(N)
disp('Theta index is:')
disp(theta_ind)
    
for sim_ind = 1:itr

rng((sim_ind*theta_ind*N_ind+itr)*rng_seed_val);




cell_rad=75;% radius of the cell in m
% N=30;        % number of Nodes 
% N = N_arr(N_ind);
% plotting circle of radius=cell_rad
x=0:0.01:(2*pi);
x_c=cell_rad.*cos(x);
y_c=cell_rad.*sin(x);
% % % plot(x_c,y_c); hold on;

% scattering Nodes uniformly within circle of radius=cell_rad and EML Server at the centre of the cell

% % % scatter(0,0,30,'black','filled');hold on; % EML Server at the center of the cell
x_cb=10.*cos(x); % red circle of 10 for boundary around EML Server
y_cb=10.*sin(x); % red circle of 10 for boundary around EML Server
% % % plot(x_cb,y_cb,'r'); hold on;

r_tx=10+((cell_rad-10).*rand([1,N])); %+10 because no transmitter inside 10m of BS...radius of node 
theta_tx=(2*pi).*rand([1,N]); %theta of node
NumberLabel = cellstr(num2str([1:N]'));
x_tx=r_tx.*cos(theta_tx);
y_tx=r_tx.*sin(theta_tx);
% % % scatter(x_tx,y_tx,10,'blue','filled'); hold on;
% % % text(x_tx+1, y_tx+1, NumberLabel);
% % % legend('Cell Border','EML Station','Boundary of EML Server','Node','Location','NorthEastOutside');


% % % rbi_gen = rand(1,N);  %generating random RBI levels from uniform distribution, for N nodes
% rbi_gen = abs(normrnd(0.7,0.4,[1,N]));
% rbi_gen = abs(normrnd(0.7,0.5,[1,N]));
% rbi_gen = abs(normrnd(0.7,0.45,[1,N]));
rbi_gen = abs(normrnd(0.7,0.3,[1,N]));
for i = 1:length(rbi_gen)
    if rbi_gen(i) > 1
%         rbi_gen(i) = (randi([3,10]))/10;
          rbi_gen(i) = (randi([4,6]))/10;
    end
end

%Initialising Node struct
Node = repmat(struct('node_id',0, 'rbi', 0, 'x_cor', 0, 'y_cor', 0, 'status', 0, 'service_status', 0, 'service', 0,'angle_status', 0, 'Adjacency_vector' , zeros(1,N)),N,1);
%Using repmat instaed of preallocating for faster execution

% % % NSA = 0.6;

for i = 1:N
    Node(i).node_id = i;    %Assigning node id
    Node(i).rbi = rbi_gen(i);   %Assigning RBI level
    Node(i).x_cor = x_tx(i);   %Assigning X Coordinate 
    Node(i).y_cor = y_tx(i);   %Assigning Y Coordinate
    if Node(i).rbi < 0.3
        Node(i).status = 0;  % 0 = Energy Harvesting state
    else
        Node(i).status = 1;  % 1 = Active state
    end
    
    Node(i).min_sinr = -110; %The minimum detectable value of the received sinr such that the Tx is received, here it is taken arbitrarly as -110dB
    
% % %     if (Node(i).status == 1 && Node(i).rbi > 0.5 && NSA > 0.5)
% % %         Node(i).P_tx = 500*10^(-3);  % If Node is active and has RBI > 0.5, the transmit power is 500mW
% % %     elseif (Node(i).status == 1 && Node(i).rbi < 0.5 && Node(i).rbi > 0.3 && NSA > 0.5)
% % %          Node(i).P_tx = 250*10^(-3);  % If Node is active and has RBI < 0.5 and > 0.3, the transmit power is 250m
% % %     elseif (Node(i).status == 1 && Node(i).rbi > 0.5 && NSA < 0.5)
% % %         Node(i).P_tx = 500*10^(-3);
% % %     elseif (Node(i).status == 1 && Node(i).rbi < 0.5 && Node(i).rbi > 0.3 && NSA < 0.5)
% % %          Node(i).P_tx = 250*10^(-3);
% % %     end
    
    if Node(i).status == 1           % creating an Adjacency vector so as to create Adjacency matrix later
        Node(i).Adjacency_vector = ones(1,N);
        Node(i).Adjacency_vector(i) = 0;
    else
        Node(i).Adjacency_vector = zeros(1,N);
    end
    
end 

Adjacency_matrix = zeros(N,N);%Creating an Adjacency matrix for the network
for i = 1:N
    Adjacency_matrix(i,:) = Node(i).Adjacency_vector;
end

for i = 1:N
     if Node(i).status == 0   %If a node is in harvesting then making sure it does not recieve
        Adjacency_matrix(:,i) = 0;
     end
end


grap = digraph(Adjacency_matrix);
%figure;
% % % plot(grap);title('Graph of all active nodes');


% computing distance( in m) between Node k and Node l
% for all k and l from 1 to N(number of links)
d_m=zeros(N,N); % distance matrix
for k=1:N % LOOP for selecting transmitter 
    for l=1:N % LOOP for selecting receiver
        d_m(k,l)=sqrt((x_tx(k)-x_tx(l))^2 + (y_tx(k)-y_tx(l))^2);
    end
end

%to count the number of active nodes
No_active_nodes = 0;
for i = 1:N
    if Node(i).status ~= 0
        No_active_nodes = No_active_nodes + 1;
    end
end


% Assigning services to Nodes to simulate busy nodes and dynamic spectrum
% We take three services SP1 (highest priority, highest BW), SP2(medium
% priority and BW) and SP3 (Lowest priority and BW)
% Service status will tell if the active node is either Tx(1) or Rx(2) or
% idle(0)
% Servive will hold what kind of service it being transmitted, p1, p2, p3
% or none(0)

% norand = abs(normrnd(0, 1.1, [1,No_active_nodes])); 
% norand = abs(normrnd(0, 1.2, [1,No_active_nodes]));
% norand = abs(normrnd(0, 1.4, [1,No_active_nodes]));

norand = abs(normrnd(0, 0.4, [1,No_active_nodes]));
% norand = abs(normrnd(0,0.7,[1,No_active_nodes])); % Generates random values from normal dist. so that a small number of values are greater than 1 and these can be used as Tx nodes
%mean 0 and var 1 gives good amount of busy nodes for rng seed 3 to get NSA
%less than 0.5

% Now we write a piece of code to check if half or half - 1 or more than
% half of the values of norand are greater than 1, if yes then we modify
% the norand value such that only less than half - 1 number of values are
% greater than 1

countnor = 0;
for i = 1:No_active_nodes
    if norand(i) >= 1
       countnor = countnor + 1; 
    end
end

countcng = countnor;
if countnor >= (No_active_nodes/2) - 1 
   disp("Error, applying correction")
   for i = 1:No_active_nodes
      if norand(i) >= 1
         if countcng >= (No_active_nodes/2) - 1
            norand(i) = rand;
            countcng = countcng - 1;
         end
      end
   end
end

% Correction is done

j = 1;
for i = 1:N
    if Node(i).status == 1
        if norand(j) >= 1 
            Node(i).service_status = 1;
            % Now the followng code is to assign a service to the Tx node with
            % service p1 occuring with 50% probability, p2 with prob. 30% and p3 with prob. 20%
            rn = rand;
            if rn <= 0.5
                Node(i).service = 1;
            elseif rn > 0.5 && rn <= 0.8
                Node(i).service = 2;
            else 
                Node(i).service = 3;
            end
        end
        j = j+1;
    end
end

% Now we need to choose the destination for the nodes which are
% transmitting, for simplicity we assume dirct source to destination
% transmission
% We also create a service matrix will holds all the dirct transmissions
% and spectrum matrix will contain the BW used for the links from source (row number)
% and destination (column number)

% First we collect all the nodes which are capable of Rx into an array
collec_arr = zeros(1,N);
Already_tx_count = 0;
for i = 1:N
    if Node(i).status == 1 && Node(i).service_status == 0 % If node is active and idle
        collec_arr(1,i) = Node(i).node_id; 
    end
    if Node(i).status == 1 && Node(i).service_status == 1 % If node is active and idle
         Already_tx_count = Already_tx_count +1; % counts the number of already Tx nodes 
    end
end


collec_arr1 = nonzeros(collec_arr);
collec_arr_len = length(collec_arr1);

skip_itr2 = 0; % skip iteration count if number of already Tx nodes is greater then the avaliable nodes for reception
if Already_tx_count > collec_arr_len
    disp(ind);
    disp("This iteration has already Tx nodes greater than nodes avaliable for reception");
    skip_itr2 = skip_itr2 + 1;
    continue;
end
    
    
% Now we rearrange collec_arr1 according to random permutation to simulate
% random selection of target
collec_arr1 = collec_arr1(randperm(length(collec_arr1)));
% We also create a service matrix and a spectrum matrix
ser_mat = zeros(N,N);
spec_mat = zeros(N,N); % Holds the services in busy nodes
spec_mat_occp = zeros(N,N);   % Holds the spectrum usage in MHz in busy nodes

% Now we assign destinations
for i = 1:N
    if Node(i).status == 0
        continue;
    end
    if Node(i).service_status == 1
        Node(collec_arr1(1)).service_status = 2;
        ser_mat(i,collec_arr1(1)) = Node(i).service;
        spec_mat(i,collec_arr1(1)) = Node(i).service;
        collec_arr1(1,1) = 0;
        collec_arr1 = nonzeros(collec_arr1);
    end
end

% Plot of all occupied paths 
grap_occupied_paths = digraph(spec_mat);
% % % figure;
% % % plot(grap_occupied_paths);
% % % title('The graph of occupied paths among the active nodes');

% Now we create a service adjacency matrix and graph avaliable for new
% transmission
ser_adj = zeros(N,N);
for i = 1:N
    if Node(i).status == 1 && Node(i).service_status == 0 % If the node is active and idle
        ser_adj(i,:) = ones(1,N);
        ser_adj(i,i) = 0;
    end
end

for i = 1:N
    if Node(i).status == 0 || Node(i).service_status == 1  || Node(i).service_status == 2 % If node is asleep or Tx or Rx, making ser_adj such that it does not participate
        ser_adj(:,i) = 0;
    end
end

No_active_avaliable_nodes = 0;
for i = 1:N
    if Node(i).status ~= 0 && Node(i).service_status == 0
        No_active_avaliable_nodes = No_active_avaliable_nodes + 1;
    end
end

skip_itr = 0; % skip iteration if there are 0 or 1 nodes avaliable for new transmission
if No_active_avaliable_nodes == 0 || No_active_avaliable_nodes == 1
    disp(ind);
    disp("This iteration has 0 free nodes");
    skip_itr = skip_itr + 1;
    continue;
end


% Plot of all avaliable paths
grap_paths = digraph(ser_adj);
% % % figure;
% % % plot(grap_paths);
% % % title('All avaliable paths for new transmission');

% 4G max channel BW is 20Mhz
% 5G uses channel BW of upto 100Mhz
% We assume the total bandwidth avaliable for the cell is 1Ghz
% We allocate BW of 80Mhz for SP1, 50Mhz for SP2 and 30Mhz for SP3
% Now we calculate the NSA based on the spectrum matrix

cell_channel_bandwidth = 1000; % In Ghz
used_spectrum = 0; %in Mhz
for i = 1:N
    for j = 1:N
        if spec_mat(i,j) == 1
            used_spectrum = used_spectrum + 80;
            spec_mat_occp(i,j) = 80;
        elseif spec_mat(i,j) == 2
            used_spectrum = used_spectrum + 50;
            spec_mat_occp(i,j) = 50;
        elseif spec_mat(i,j) == 3
            used_spectrum = used_spectrum + 30;
            spec_mat_occp(i,j) = 30;
        elseif spec_mat(i,j) == 0
            used_spectrum = used_spectrum + 0;
            spec_mat_occp(i,j) = 0;
        end
    end
end

if (used_spectrum > cell_channel_bandwidth)
    disp('Error: Spectrum allocated for services exceeds maximum cell chanel bandwidth');
end

NSA = (cell_channel_bandwidth - used_spectrum)/cell_channel_bandwidth;

% Now we decide on the Ptx of nodes based on the RBI levels and NSA
% according to BCG matrix
Higher_Tx_power = 500*10^(-3);
Lower_Tx_power = 250*10^(-3);
for i = 1:N    
    if (Node(i).status == 1 && Node(i).rbi >= 0.5 && NSA >= 0.5)
% % %         Node(i).P_tx = 500*10^(-3);  % If Node is active and has RBI > 0.5, the transmit power is 500mW
        Node(i).P_tx = Higher_Tx_power;
    elseif (Node(i).status == 1 && Node(i).rbi < 0.5 && Node(i).rbi >= 0.3 && NSA >= 0.5)
% % %          Node(i).P_tx = 250*10^(-3);  % If Node is active and has RBI < 0.5 and > 0.3, the transmit power is 250m
         Node(i).P_tx = Lower_Tx_power;
    elseif (Node(i).status == 1 && Node(i).rbi >= 0.5 && NSA < 0.5)
% % %         Node(i).P_tx = 500*10^(-3);
        Node(i).P_tx = Higher_Tx_power;
    elseif (Node(i).status == 1 && Node(i).rbi < 0.5 && Node(i).rbi >= 0.3 && NSA < 0.5)
% % %          Node(i).P_tx = 250*10^(-3);
         Node(i).P_tx = Lower_Tx_power;
    end
end


f = 1000*10^9; %Carrier Frequency
n = 2; %Path Loss Component
grv = 3.1449;  %AWGN noise with 0 mean and 3.1449db variance
h_dB = zeros(N,N);
h_mag = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i == j
            h_dB(i,i) = NaN; 
            h_mag(i,i) = NaN;
            continue;
        end
        
        h_dB(i,j) = (20*(log10((4*pi*f)/(3*(10^8))))) + (10*n*log10(d_m(i,j))) + grv;   % Path Loss matrix in dB
        h_mag(i,j) = 10^(h_dB(i,j)/10);
    end
    
end

Ptx_inf = 100*10^(-3); %Interference from other nodes... Tx power is taken as 100mW
No_dBm = -174;   %AWGN Noise is taken as -174 dBm/Hz
B = 50*10^(9);   % Bandwidth is 50GHz 
No_W = (10^(No_dBm/10))*1000;
N_B = No_W*B;    %AWGN noise power

P_inf = zeros(1,N);
P_inf_dB = zeros(1,N);
for i = 1:N
    if Node(i).status == 0
        P_inf(i) = NaN;
        P_inf_dB(i) = NaN;
        continue;
    end
    for j = 1:N
        if j == i || Node(j).status == 0
            continue;
        end
       P_inf(i) = P_inf(i) + (Ptx_inf/h_mag(i,j)); 
       P_inf_dB(i) = 10*log10(P_inf(i));
    end
end



% % % Src = input('Enter Source Node ID \n');
% % % Des = input('Enter destination Node ID \n');
% % % set = 0;
% % % for i = 1:N
% % %     for j = 1:N
% % %         if ser_adj(i,j) == 1
% % %             Src = i;
% % %             Des = j;
% % %             set = 1;
% % %             break;
% % %         end
% % %     end
% % %     if set == 1
% % %         break;
% % %     end
% % % end

[r1,c1] = find(ser_adj); % Getting the row and coloumn numbers of all non zero values, i.e, all nodes between which new tx can occur 
ran_si = size(r1);
ran_tem = ran_si(1);
ran_ind = randi(ran_tem); % picking a random row and col number
Src = r1(ran_ind);
Des = c1(ran_ind);

if Node(Src).status == 0
    disp('Source is in Harvesting mode in iteration \n');
    disp(ind);
    return
elseif Node(Des).status == 0
    disp('Destination is in Harvesting mode in iteration \n');
    disp(ind);
    return
elseif Node(Src).service_status ~= 0
    disp('Source is already Transmitting or Receiving in iteration \n');
    disp(ind);
    disp(Src);
    disp(Des);
elseif Node(Des).service_status ~=0
    disp('Destination is already Transmitting or Receiving in iteration \n');
    disp(ind);
    disp(Src);
    disp(Des);
end


% Angle_status of a node will tell if the node lies within the cone or not
% if it is 1, the node lies inside the one else if it is 0, it will not

Node(Src).angle_status = 1; % Making the angle_status of the Src and Des nodes as 1
Node(Des).angle_status = 1;

% Now we create a N*1 matrix, ang_vec, which will hold all the absolute
% angles (0 to 180 degrees) of all the nodes with respect to the Src Des
% vector

% We also create a N*2 matrix, vec_mat which will hold the displacement vector or position vector from the
% Src to every other Node. Column 1 holds x^ and column 2 holds y^


tA3 = tic;


% Algorithm 3

%%% Finding best path between Src and Des

temp_SINR = zeros(N,N); % Contains all the Rx SINR values
w = zeros(N,N);
for i = 1:N
    if Node(i).status == 0 || Node(i).service_status ~= 0
        continue;
    end
    for j = 1:N
        if Node(j).status == 0 || Node(j).service_status ~= 0
            continue;
        end
        if i == j
            continue;
        end
        temp_SINR(i,j) = ((Node(i).P_tx)/(h_mag(i,j)))/(P_inf(j) + N_B); % Calculating the SINR from Node i to Node j
        w(i,j) = log2(1 + temp_SINR(i,j));
    end
end



W = zeros(N,N);
L_W = zeros(N,N);

% Initialization 

W(:,Src) = inf;
pred = zeros(N,N);
L_W(Src,:) = inf;
L_W(:,Src) = 0;

% Algorithm

for h = 1:(N-1)  
    for i = 1:N
        W((h+1),i) = W(h,i);
    end
  
    for i = 1:N
        for j = 1:N
            if i == j
                continue;
            end
            if min(W(h,i), w(i,j)) > W((h+1),j)
                W((h+1),j) = min(W(h,i), w(i,j));
                L_W((h+1),j) = L_W(h,i) + 1;
                pred((h+1),j) = i;
            end
            
        end
    end
end

temp_h = 0;
paths_cell = cell(N,1);
for i = N:-1:1
    temp_path = 0;
    path_counter = 0;
    if pred(i,Des) == 0
        continue;
    end
    temp_h = i;
    temp_path(1) = Des;
    path_counter = path_counter + 1;
    while temp_h ~= 1
        temp_path(end+1) = pred(temp_h, temp_path(path_counter));
        path_counter = path_counter + 1;
        temp_h = temp_h - 1;
    end
    temp_path = flip(temp_path);
    paths_cell(i) = mat2cell(temp_path, 1, length(temp_path));
    
end

temp_SE = zeros(N,1);
for i = 1:N
    if isinf(L_W(i,Des)) || ~any(L_W(i,Des))
        continue;
    end
    temp_SE(i,1) = W(i,Des)/(L_W(i,Des));
end

[max_SE, SE_ind] = max(temp_SE);
best_path_algo3 = cell2mat(paths_cell(SE_ind));

if N == 30 
    time_N30_A3(theta_ind) = time_N30_A3(theta_ind) + toc(tA3);
elseif N == 40
    time_N40_A3(theta_ind) = time_N40_A3(theta_ind) + toc(tA3);
elseif N == 50
    time_N50_A3(theta_ind) = time_N50_A3(theta_ind) + toc(tA3);
elseif N == 60
    time_N60_A3(theta_ind) = time_N60_A3(theta_ind) + toc(tA3);
elseif N == 70
    time_N70_A3(theta_ind) = time_N70_A3(theta_ind) + toc(tA3);
elseif N == 80
    time_N80_A3(theta_ind) = time_N80_A3(theta_ind) + toc(tA3);
end

% Algorithm 4


dis_th = 10; % The nodes within 12 m radius of Src and Des are also considered as within the cone regardless 
% of whether they are or not
theta_th = theta_arr(theta_ind);



main_vec = [0 0]; % Holds the Src to Des vector
main_vec(1,1) = x_tx(Des) - x_tx(Src);
main_vec(1,2) = y_tx(Des) - y_tx(Src);
main_vec_mag = sqrt((main_vec(1,1))^2 + (main_vec(1,2))^2);

vec_mat = zeros(N,2); % Holds the poistiuon vector of every node wrt Src node except Des node
ang_vec = zeros(N,1); % Holds the angles b/w the Src Des vector and position vection of all other nodes
for i = 1:N
    temp_ang = 0;
    if i == Src || i == Des % Skipping Src and Des nodes by maikng them nan
        vec_mat(i,1) = nan;
        vec_mat(i,2) = nan;
        ang_vec(i,1) = nan;
        continue;
    end
    vec_mat(i,1) = x_tx(i) - x_tx(Src);
    vec_mat(i,2) = y_tx(i) - y_tx(Src);
    temp_ang = ((vec_mat(i,1)*main_vec(1,1)) + (vec_mat(i,2)*main_vec(1,2)))/((main_vec_mag)*(sqrt((vec_mat(i,1))^2 + (vec_mat(i,2))^2)));
    ang_vec(i,1) = acosd(temp_ang);
    if ang_vec(i,1) <= theta_th
        Node(i).angle_status = 1;
    end
    dist_Src_node =  sqrt((x_tx(Src) - x_tx(i))^2 + (y_tx(Src) - y_tx(i))^2); % Distance between the node 'i' and the Src node
    dist_Des_node = sqrt((x_tx(Des) - x_tx(i))^2 + (y_tx(Des) - y_tx(i))^2); % Distance between the node 'i' and the Des node
    if dist_Src_node <= dis_th || dist_Des_node <= dis_th
         Node(i).angle_status = 1;
    end
end

tA4 = tic;

temp_SINR = zeros(N,N); % Contains all the Rx SINR values
w = zeros(N,N); 

% % % for i = 1:N
% % %     dist_Src_node =  sqrt((x_tx(Src) - x_tx(i))^2 + (y_tx(Src) - y_tx(i))^2); % Distance between the node 'i' and the Src node
% % %     dist_Des_node = sqrt((x_tx(Des) - x_tx(i))^2 + (y_tx(Des) - y_tx(i))^2); % Distance between the node 'i' and the Des node
% % %     if dist_Src_node <= dis_th || dist_Des_node <= dis_th
% % %          Node(i).angle_status = 1;
% % %     end
% % % end

% % % for i = 1:N
% % %     if Node(i).status == 0 || Node(i).service_status ~= 0 || Node(i).angle_status ~= 1
% % %         continue;
% % %     end
% % %     for j = 1:N
% % %         if Node(j).status == 0 || Node(j).service_status ~= 0 || Node(i).angle_status ~= 1
% % %             continue;
% % %         end
% % %         if i == j
% % %             continue;
% % %         end
% % %         temp_SINR(i,j) = ((Node(i).P_tx)/(h_mag(i,j)))/(P_inf(j) + N_B); % Calculating the SINR from Node i to Node j
% % %         w(i,j) = log2(1 + temp_SINR(i,j));
% % %     end
% % % end


N_new = zeros(1, N);
n_count = 1;
for i = 1:N
    if Node(i).status == 0 || Node(i).service_status ~= 0 || Node(i).angle_status ~= 1
        continue;
    end
   
    N_new(n_count) = i;
    n_count = n_count + 1;
    
end
N_new = nonzeros(N_new);
N_new_len = length(N_new);
temp_SINR = zeros(N_new_len,N_new_len); % Contains all the Rx SINR values
w = zeros(N_new_len,N_new_len);

for i = 1:N_new_len
    for j = 1:N_new_len
        if i == j
            continue;
        end
        temp_SINR(i,j) = ((Node(N_new(i)).P_tx)/(h_mag(N_new(i),N_new(j))))/(P_inf(N_new(j)) + N_B); % Calculating the SINR from Node i to Node j
        w(i,j) = log2(1 + temp_SINR(i,j));
     end
end


W = zeros(N_new_len,N_new_len);
L_W = zeros(N_new_len,N_new_len);

% % % Src_new = find(N_new == Src);
% % % Des_new = find(N_new == Des);

Src = find(N_new == Src);
Des = find(N_new == Des);

% Initialization 

W(:,Src) = inf;
pred = zeros(N_new_len,N_new_len);
L_W(Src,:) = inf;
L_W(:,Src) = 0;

% Algorithm

for h = 1:(N_new_len-1)  
    for i = 1:N_new_len
        W((h+1),i) = W(h,i);
    end
  
    for i = 1:N_new_len
        for j = 1:N_new_len
            if i == j
                continue;
            end
            if min(W(h,i), w(i,j)) > W((h+1),j)
                W((h+1),j) = min(W(h,i), w(i,j));
                L_W((h+1),j) = L_W(h,i) + 1;
                pred((h+1),j) = i;
            end
            
        end
    end
end

temp_h = 0;
paths_cell = cell(N_new_len,1);
for i = N_new_len:-1:1
    temp_path = 0;
    path_counter = 0;
    if pred(i,Des) == 0
        continue;
    end
    temp_h = i;
    temp_path(1) = Des;
    path_counter = path_counter + 1;
    while temp_h ~= 1
        temp_path(end+1) = pred(temp_h, temp_path(path_counter));
        path_counter = path_counter + 1;
        temp_h = temp_h - 1;
    end
    temp_path = flip(temp_path);
    paths_cell(i) = mat2cell(temp_path, 1, length(temp_path));
    
end

temp_SE = zeros(N_new_len,1);
for i = 1:N_new_len
    if isinf(L_W(i,Des)) || ~any(L_W(i,Des))
        continue;
    end
    temp_SE(i,1) = W(i,Des)/(L_W(i,Des));
end

[max_SEA4, SE_indA4] = max(temp_SE);
temp_best_path_algo4 = cell2mat(paths_cell(SE_indA4));
temp_len_best = length(temp_best_path_algo4);

best_path_algo4 = zeros(1,N_new_len);
for i = 1:temp_len_best
    best_path_algo4(i) = N_new(temp_best_path_algo4(i));
end
best_path_algo4 = nonzeros(best_path_algo4);
best_path_algo4 = best_path_algo4';

if N == 30 
    time_N30_A4(theta_ind) = time_N30_A4(theta_ind) + toc(tA4);
elseif N == 40
    time_N40_A4(theta_ind) = time_N40_A4(theta_ind) + toc(tA4);
elseif N == 50
    time_N50_A4(theta_ind) = time_N50_A4(theta_ind) + toc(tA4);
elseif N == 60
    time_N60_A4(theta_ind) = time_N60_A4(theta_ind) + toc(tA4);
elseif N == 70
    time_N70_A4(theta_ind) = time_N70_A4(theta_ind) + toc(tA4);
elseif N == 80
    time_N80_A4(theta_ind) = time_N80_A4(theta_ind) + toc(tA4);
end

if N == 30
    if length(best_path_algo3) == length(best_path_algo4)
        if best_path_algo3 == best_path_algo4
            acc_N30(theta_ind) = acc_N30(theta_ind) +  1;
        end
    end
elseif N == 40
    if length(best_path_algo3) == length(best_path_algo4)
        if best_path_algo3 == best_path_algo4
            acc_N40(theta_ind) = acc_N40(theta_ind) +  1;
        end
    end
elseif N == 50
    if length(best_path_algo3) == length(best_path_algo4)
        if best_path_algo3 == best_path_algo4
            acc_N50(theta_ind) = acc_N50(theta_ind) +  1;
        end
    end
elseif N == 60
    if length(best_path_algo3) == length(best_path_algo4)
        if best_path_algo3 == best_path_algo4
            acc_N60(theta_ind) = acc_N60(theta_ind) +  1;
        end
    end
elseif N == 70
    if length(best_path_algo3) == length(best_path_algo4)
        if best_path_algo3 == best_path_algo4
            acc_N70(theta_ind) = acc_N70(theta_ind) +  1;
        end
    end
elseif N == 80
    if length(best_path_algo3) == length(best_path_algo4)
        if best_path_algo3 == best_path_algo4
            acc_N80(theta_ind) = acc_N80(theta_ind) +  1;
        end
    end
end


end

end

end

avg_acc_N30 = acc_N30/itr;
avg_acc_N40 = acc_N40/itr;
avg_acc_N50 = acc_N50/itr;
avg_acc_N60 = acc_N60/itr;
avg_acc_N70 = acc_N70/itr;
avg_acc_N80 = acc_N80/itr;
y_acc = [avg_acc_N30; avg_acc_N40; avg_acc_N50; avg_acc_N60; avg_acc_N70; avg_acc_N80];

%%

avg_time_N30_A3 = time_N30_A3/itr;
avg_time_N40_A3 = time_N40_A3/itr;
avg_time_N50_A3 = time_N50_A3/itr;
avg_time_N60_A3 = time_N60_A3/itr;
avg_time_N70_A3 = time_N70_A3/itr;
avg_time_N80_A3 = time_N80_A3/itr;
y_time_A3 = [avg_time_N30_A3, avg_time_N40_A3, avg_time_N50_A3, avg_time_N60_A3, avg_time_N70_A3, avg_time_N80_A3];

avg_time_N30_A4 = time_N30_A4/itr;
avg_time_N40_A4 = time_N40_A4/itr;
avg_time_N50_A4 = time_N50_A4/itr;
avg_time_N60_A4 = time_N60_A4/itr;
avg_time_N70_A4 = time_N70_A4/itr;
avg_time_N80_A4 = time_N80_A4/itr;
y_time_A4 = [avg_time_N30_A4, avg_time_N40_A4, avg_time_N50_A4, avg_time_N60_A4, avg_time_N70_A4, avg_time_N80_A4];
x_time = 1:1:30;


figure;
hold on;
grid on;
title('Comparison of Algo 3 and Algo 4 in Acc');
xlabel('Number of Nodes');
ylabel('Accuracy');
bar_acc = bar(N_arr, y_acc);
set(bar_acc(1), 'FaceColor','b');
set(bar_acc(2), 'FaceColor','r');
set(bar_acc(3), 'FaceColor','g');
set(bar_acc(4), 'FaceColor','y');
set(bar_acc(5), 'FaceColor','c');
legend('\theta_{th} = 40', '\theta_{th} = 60', '\theta_{th} = 80', '\theta_{th} = 90', '\theta_{th} = 100');

figure;
hold on;
grid on;
title('Comparison of Algo 3 and Algo 4 in time');
xlabel('Number of Nodes');
ylabel('Time');
plot(x_time, y_time_A3,'b');
plot(x_time, y_time_A4,'r');
legend('Algo3', 'Algo4');

% % % figure;
% % % hold on;
% % % grid on;
% % % yyaxis left
% % % b = bar(N_arr, y_acc);
% % % xoff1=b(1).XOffset;                       % the offset of bars in each group (hidden property)
% % % dx1 = [xoff1, -xoff1, 5].';
% % % xoff2=b(2).XOffset;               
% % % dx2 = [xoff2, -xoff2, 5].';
% % % xoff3=b(3).XOffset;               
% % % dx3 = [xoff3, -xoff3, 5].';
% % % xoff4=b(4).XOffset;               
% % % dx4 = [xoff4, -xoff4, 5].';
% % % xoff5=b(5).XOffset;               
% % % dx5 = [xoff5, -xoff5, 5].';
% % % xoff6=b(6).XOffset;               
% % % dx6 = [xoff6, -xoff6, 5].';
% % % yyaxis right
% % % plot(1+dx1,avg_time_N30_A3,'b');
% % % plot(1+dx1,avg_time_N30_A4,'r');
% % % plot(1+dx2,avg_time_N40_A3,'b');
% % % plot(1+dx2,avg_time_N40_A4,'r');
% % % plot(1+dx3,avg_time_N50_A3,'b');
% % % plot(1+dx3,avg_time_N50_A4,'r');
% % % plot(1+dx4,avg_time_N60_A3,'b');
% % % plot(1+dx4,avg_time_N60_A4,'r');
% % % plot(1+dx5,avg_time_N70_A3,'b');
% % % plot(1+dx5,avg_time_N70_A4,'r');
% % % plot(1+dx6,avg_time_N80_A3,'b');
% % % plot(1+dx6,avg_time_N80_A4,'r');
% % % ylim([0 0.009])
% p = plot(N_arr,y_time_A3);



tEnd = toc(tStart);
disp("Total time taken for the program is:")
disp(tEnd)


