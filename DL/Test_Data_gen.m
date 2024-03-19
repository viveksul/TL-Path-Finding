%%....Main Script to Generate nodes and calculate the SE,EE and SINR margins in a simulated network...%%
clc
clear
close all

%rng(1);
tic
itrations = 10000;
rng_seed_val = 1003;
writematrix(rng_seed_val, "Rng_seed_val.txt")
for ind = 1:itrations
rng(ind*rng_seed_val);    
cell_rad=75;% radius of the cell in m
N=30;        % number of Nodes 
% plotting circle of radius=cell_rad
x=0:0.01:(2*pi);
x_c=cell_rad.*cos(x);
y_c=cell_rad.*sin(x);
% % % plot(x_c,y_c); hold on;

% scattering Nodes uniformly within circle of radius=cell_rad and EML Server at the centre of the cell

%scatter(0,0,30,'black','filled');hold on; % EML Server at the center of the cell
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
rbi_gen = abs(normrnd(0.7,0.4,[1,N]));
for i = 1:length(rbi_gen)
    if rbi_gen(i) > 1
        rbi_gen(i) = (randi([3,10]))/10;
    end
end

%Initialising Node struct
Node = repmat(struct('node_id',0, 'rbi', 0, 'x_cor', 0, 'y_cor', 0, 'status', 0, 'service_status', 0, 'service', 0, 'Adjacency_vector' , zeros(1,N)),N,1);
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
% % % figure;plot(grap);title('Graph of all active nodes');


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

% norand = abs(normrnd(0,0.8,[1,No_active_nodes])); % Generates random values from normal dist. so that a small number of values are greater than 1 and these can be used as Tx nodes
norand = abs(normrnd(0, 1.0, [1,No_active_nodes]));
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
    disp("This iteration has already Tx nodes greater than nodes avaliable for recepytion");
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


skip_itr = 0;
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
% % % for i = 1:N
% % %     for j = 1:N
% % %         if ser_adj(i,j) == 1
% % %             Src = i;
% % %             Des = j;
% % %         end
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
elseif Node(Des).service_status ~=0
    disp('Destination is already Transmitting or Receiving in iteration \n');
    disp(ind);
end

Src_Des_info = zeros(N,N);
Src_Des_info(Src,Des) = 1;

%for N=10...5,10
%for N=7....1,6
%for N=8...6,7
%for N=9...5,4
%for N=11...5,10
%for N=12...5,10 or 2,6 or 5,6
%for N=12 and var = 1.5... 6,9
% for N= 15 and var = 1 ... 4,3


if No_active_avaliable_nodes > 12
    paths = allpaths(grap_paths,Src,Des, 'MaxPathLength', 6);
else
    paths = allpaths(grap_paths,Src,Des);
    %paths = allpaths(grap_paths,4,3);
end
i = size(paths);
psize = i(1);



pat_mat = zeros(psize, No_active_nodes);

for i = 1:psize
    p = paths(i,:);
    p1 = cell2mat(p);
    pat_mat(i,1:length(p1)) = p1;
end
%pat_mat contains the paths in each row with zeros appended at end

%To perform calculations on each path
SINR = nan(psize,No_active_nodes); %SINR matrix which contains the rx_sinr vectors of each path in corresponding row
SINR_dB = nan(psize,No_active_nodes); %SINR matrix in dB
SINR_m_dB = nan(psize,No_active_nodes); %SINR margin matrix in dB
SE = zeros(psize,1);
SE_dB = zeros(psize,1);
EE = zeros(psize,1);
EE_dB = zeros(psize,1);
for i = 1:psize
    relay = pat_mat(i,:); %To select one path from the matrix
    relay = nonzeros(relay)'; %To remove zero padding at the end of each path
    L = length(relay); %To store length of the selected path
    %To calculate Rx_SINR for each node in each path
    sinr_vec = zeros(1,L);
    sinr_vec(1) = nan;
    for j = 2:L
        k = j-1;
        sinr_vec(j) =  ((Node(relay(k)).P_tx)/(h_mag(relay(k),relay(j))))/(P_inf(relay(j)) + N_B);%calculating the rx_sinr for each node in the path
    end
    SINR(i,1:length(sinr_vec)) = sinr_vec; % Putting the sinr vector of the path into the full sinr matrix
    SINR_dB(i,1:length(sinr_vec)) = 10*log10(sinr_vec); %Putting the sinr vector(in dB) of the path into the full sinr(dB) matrix
     
    
    %To calculate SE of the selected path
    SE(i,1) = (1/(L-1))*(log2(1 + min(sinr_vec)));
    SE_dB(i,1) = 10*log10(SE(i,1));
    
    %To calculate EE of the selected path
    Power_tx = 0;
    for j = 1:(L-1)
    Power_tx = Power_tx + Node(relay(j)).P_tx;
    end
    EE(i,1) = SE(i,1)/(Power_tx + 0.010);
    EE_dB(i,1) = 10*log10(EE(i,1));
  
end

%Calculating the SINR margin matrix in dB
for i = 1:psize
    for j = 2:No_active_nodes
        SINR_m_dB(i,j) = SINR_dB(i,j) - (-110);
    end
end

  
[max_SE, SE_ind] = max(SE_dB);
[max_EE, EE_ind] = max(EE_dB);

output_power = zeros(N,N);
best_path = cell2mat(paths(SE_ind));
var_s = size(best_path);
var_size = var_s(2);
for i = 1:(var_size-1)
    output_power(best_path(i),best_path(i+1)) = Node(best_path(i)).P_tx;
end

% % % figure; pl = plot(grap_paths);
% % % highlight(pl,paths{SE_ind},'EdgeColor','k','LineWidth',1.5,'NodeColor','k','MarkerSize',6)
% % % title('Path with highest SE')
% % % figure; pl = plot(grap_paths);
% % % highlight(pl,paths{EE_ind},'EdgeColor','magenta','LineWidth',1.5,'NodeColor','magenta','MarkerSize',6)
% % % title('Path with highest EE')

%%% Writing output in different way
out = zeros(1,N);
temp = nonzeros(pat_mat(SE_ind,:))';
out_size = size(temp);
for i = 1:(out_size(2)-1)
    out(1,temp(i)) = temp(i+1);
end

%%% Writing Source Destination in a different way
SrcDes = [Src, Des];

%%% Writing the channel matrices, spectrum matrices, rbi vector and source
%%% and destination information

%%% We also insert a header of 1:N array at the top of each file so that
%%% we can read in Jupyter notebooks with pandas properly
header = 1:N;
if ind == 1
    writematrix(header, "Test_Data_channel_matrices.csv");
    writematrix("","Test_Data_channel_matrices.csv", "writemode", "append");
    
    writematrix(h_dB, "Test_Data_channel_matrices.csv", "writemode", "append");
    writematrix("","Test_Data_channel_matrices.csv", "writemode", "append");
    
    
    writematrix(header, "Test_Data_spectrum_occpancy_matrix.csv");
    writematrix("", "Test_Data_spectrum_occpancy_matrix.csv", "writemode", "append");
    
    writematrix(spec_mat_occp, "Test_Data_spectrum_occpancy_matrix.csv", "writemode", "append");
    writematrix("", "Test_Data_spectrum_occpancy_matrix.csv", "writemode", "append");
    
    
    writematrix(header(1:1), "Test_Data_spectrum_matrix.csv");
    writematrix("", "Test_Data_spectrum_matrix.csv", "writemode", "append");
    
    writematrix(NSA, "Test_Data_spectrum_matrix.csv", "writemode", "append");
    writematrix("", "Test_Data_spectrum_matrix.csv", "writemode", "append");
    
    
    writematrix(header, "Test_Data_rbi_vectors.csv");
    writematrix("", "Test_Data_rbi_vectors.csv", "writemode", "append");
    
    writematrix(rbi_gen, "Test_Data_rbi_vectors.csv", "writemode", "append");
    writematrix("", "Test_Data_rbi_vectors.csv", "writemode", "append");
    
    
    writematrix(header, "Test_Data_Source_Destination_information_matrix.csv");
    writematrix("", "Test_Data_Source_Destination_information_matrix.csv", "writemode", "append");
    
    writematrix(Src_Des_info, "Test_Data_Source_Destination_information_matrix.csv", "writemode", "append");
    writematrix("", "Test_Data_Source_Destination_information_matrix.csv", "writemode", "append");
    
    
    writematrix(header(1:2), "Test_Data_Source_Destination_info_matrix(only Src Des).csv");
    writematrix("", "Test_Data_Source_Destination_info_matrix(only Src Des).csv", "writemode", "append");
    
    writematrix(SrcDes, "Test_Data_Source_Destination_info_matrix(only Src Des).csv", "writemode", "append");
    writematrix("", "Test_Data_Source_Destination_info_matrix(only Src Des).csv", "writemode", "append");
    
    
    writematrix(header, "Test_Data_output_power.csv");
    writematrix("", "Test_Data_output_power.csv", "writemode", "append");
    
    writematrix(output_power, "Test_Data_output_power.csv", "writemode", "append");
    writematrix("", "Test_Data_output_power.csv", "writemode", "append");
    
    
    writematrix(header, "Test_Data_output.csv");
    writematrix("", "Test_Data_output.csv", "writemode", "append");
    
    writematrix(out, "Test_Data_output.csv", "writemode", "append");
    writematrix("", "Test_Data_output.csv", "writemode", "append");
    
    
    writematrix(header, "Test_Data_output1.csv");
    writematrix("", "Test_Data_output1.csv", "writemode", "append");
    
    writematrix(pat_mat(SE_ind,:), "Test_Data_output1.csv", "writemode", "append");
    writematrix("", "Test_Data_output1.csv", "writemode", "append");

else    
writematrix(h_dB, "Test_Data_channel_matrices.csv", "writemode", "append");
writematrix("","Test_Data_channel_matrices.csv", "writemode", "append");

writematrix(spec_mat_occp, "Test_Data_spectrum_occpancy_matrix.csv", "writemode", "append");
writematrix("", "Test_Data_spectrum_occpancy_matrix.csv", "writemode", "append");

writematrix(NSA, "Test_Data_spectrum_matrix.csv", "writemode", "append");
writematrix("", "Test_Data_spectrum_matrix.csv", "writemode", "append");

writematrix(rbi_gen, "Test_Data_rbi_vectors.csv", "writemode", "append");
writematrix("", "Test_Data_rbi_vectors.csv", "writemode", "append");

writematrix(Src_Des_info, "Test_Data_Source_Destination_information_matrix.csv", "writemode", "append");
writematrix("", "Test_Data_Source_Destination_information_matrix.csv", "writemode", "append");

writematrix(SrcDes, "Test_Data_Source_Destination_info_matrix(only Src Des).csv", "writemode", "append");
writematrix("", "Test_Data_Source_Destination_info_matrix(only Src Des).csv", "writemode", "append");

writematrix(output_power, "Test_Data_output_power.csv", "writemode", "append");
writematrix("", "Test_Data_output_power.csv", "writemode", "append");

writematrix(out, "Test_Data_output.csv", "writemode", "append");
writematrix("", "Test_Data_output.csv", "writemode", "append");

writematrix(pat_mat(SE_ind,:), "Test_Data_output1.csv", "writemode", "append");
writematrix("", "Test_Data_output1.csv", "writemode", "append");
end

disp(ind);

end
toc
