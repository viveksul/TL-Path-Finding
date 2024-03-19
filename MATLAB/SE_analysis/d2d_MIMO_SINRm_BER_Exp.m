clc
clear
close all
rng(4);
tic;
cell_rad=75;% radius of the cell in m
N=11;        % number of Nodes 
Ntx4 = 4;
Nrx4 = 4;
Ntx16 = 16;
Nrx16 = 16;
Ntx32 = 32;
Nrx32 = 32;

c = 3e8;        % propagation speed
fc = 1000e9;      % carrier frequency
lambda = c/fc;  % wavelength


% plotting circle of radius=cell_rad
x=0:0.01:(2*pi);
x_c=cell_rad.*cos(x);
y_c=cell_rad.*sin(x);
plot(x_c,y_c); hold on;
% scattering Nodes uniformly within circle of radius=cell_rad and EML Server at the centre of the cell
scatter(0,0,30,'black','filled');hold on; % EML Server at the center of the cell
x_cb=10.*cos(x); % red circle of 10 for boundary around EML Server
y_cb=10.*sin(x); % red circle of 10 for boundary around EML Server
plot(x_cb,y_cb,'r'); hold on;
r_tx=10+((cell_rad-10).*rand([1,N])); %+10 because no transmitter inside 10m of BS...radius of node 
theta_tx=(2*pi).*rand([1,N]); %theta of node
NumberLabel = cellstr(num2str([1:N]'));
x_tx=r_tx.*cos(theta_tx);
y_tx=r_tx.*sin(theta_tx);
scatter(x_tx,y_tx,10,'blue','filled'); hold on;
text(x_tx+1, y_tx+1, NumberLabel);
legend('Cell Border','EML Station','Boundary of EML Server','Node','Location','NorthEastOutside');


% % % rbi_gen = rand(1,N);  %generating random RBI levels from uniform distribution, for N nodes
% rbi_gen = abs(normrnd(0.7,0.4,[1,N]));
rbi_gen = abs(normrnd(0.7,0.5,[1,N]));
for i = 1:length(rbi_gen)
    if rbi_gen(i) > 1
%         rbi_gen(i) = (randi([3,10]))/10;
          rbi_gen(i) = (randi([3,5]))/10;
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
    
    Node(i).min_sinr = -50; %The minimum detectable value of the received sinr such that the Tx is received, here it is taken arbitrarly as -110dB
    
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
figure;
plot(grap);title('Graph of all active nodes');


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
norand = abs(normrnd(0, 0.01, [1,No_active_nodes]));
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
%     disp(ind);
    disp("This iteration has already Tx nodes greater than nodes avaliable for reception");
    skip_itr2 = skip_itr2 + 1;
    return;
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
figure;
plot(grap_occupied_paths);
title('The graph of occupied paths among the active nodes');

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
%     disp(ind);
    disp("This iteration has 0 free nodes");
    skip_itr = skip_itr + 1;
    return;
end


% Plot of all avaliable paths
grap_paths = digraph(ser_adj);
figure;
plot(grap_paths);
title('All avaliable paths for new transmission');

% 4G max channel BW is 20Mhz
% 5G uses channel BW of upto 100Mhz
% We assume the total bandwidth avaliable for the cell is 1Ghz
% We allocate BW of 80Mhz for SP1, 50Mhz for SP2 and 30Mhz for SP3
% Now we calculate the NSA based on the spectrum matrix

cell_channel_bandwidth = 1000; % In Mhz
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
Higher_Tx_power = 1500*10^(-3);
Lower_Tx_power = 1000*10^(-3);
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

% % % [r1,c1] = find(ser_adj); % Getting the row and coloumn numbers of all non zero values, i.e, all nodes between which new tx can occur 
% % % ran_si = size(r1);
% % % ran_tem = ran_si(1);
% % % ran_ind = randi(ran_tem); % picking a random row and col number
% % % Src = r1(ran_ind);
% % % Des = c1(ran_ind);

Src = 7;
Des = 11;

if Node(Src).status == 0
    disp('Source is in Harvesting mode in iteration \n');
%     disp(ind);
    return
elseif Node(Des).status == 0
    disp('Destination is in Harvesting mode in iteration \n');
%     disp(ind);
    return
elseif Node(Src).service_status ~= 0
    disp('Source is already Transmitting or Receiving in iteration \n');
%     disp(ind);
    disp(Src);
    disp(Des);
elseif Node(Des).service_status ~=0
    disp('Destination is already Transmitting or Receiving in iteration \n');
%     disp(ind);
    disp(Src);
    disp(Des);
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


SE_SISO = zeros(psize, N);
SINR_SISO = nan(psize, N);
SINR_dB_SISO = nan(psize, N);
SINR_m_dB_SISO = nan(psize, N);
SE4 = zeros(psize, N);
SINR4 = nan(psize, N);
SINR4_dB = nan(psize, N);
SINR4_m_dB = nan(psize, N);
SE16 = zeros(psize, N);
SINR16_dB = nan(psize, N);
SINR16 = nan(psize, N);
SINR16_m_dB = nan(psize, N);
SE32 = zeros(psize, N);
SINR32 = nan(psize, N);
SINR32_dB = nan(psize, N);
SINR32_m_dB = nan(psize, N);
SE_SISO_path = zeros(psize, 1);
SE_SISO_path1 = zeros(psize, 1);
SE_SISO_dB = zeros(psize, 1);
SE4_path = zeros(psize, 1);
SE4_path1 = zeros(psize, 1);
SE4_dB = zeros(psize, 1);
SE16_path = zeros(psize, 1);
SE16_path1 = zeros(psize, 1);
SE16_dB = zeros(psize, 1);
SE32_path = zeros(psize, 1);
SE32_path1 = zeros(psize, 1);
SE32_dB = zeros(psize, 1);

EE_SISO = zeros(psize, 1);
EE_SISO_dB = zeros(psize, 1);
EE4 = zeros(psize, 1);
EE4_dB = zeros(psize, 1);
EE16 = zeros(psize, 1);
EE16_dB = zeros(psize, 1);
EE32 = zeros(psize, 1);
EE32_dB = zeros(psize, 1);

ber_sisomp = zeros(psize, N);
ber_mimomp4 = zeros(psize, N);
ber_mimomp16 = zeros(psize, N);
ber_mimomp32 = zeros(psize, N);
ber_sisomp_btlnk = zeros(psize, 1);
ber_mimomp4_btlnk = zeros(psize, 1);
ber_mimomp16_btlnk = zeros(psize, 1);
ber_mimomp32_btlnk = zeros(psize, 1);

relsiz = zeros(psize,1); %This matrix stores the length of each corresponding path, to be used in results

a_time = 0;
b_time = 0;



for k = 1:psize
    relay = pat_mat(k,:); %To select one path from the matrix
    relay = nonzeros(relay)'; %To remove zero padding at the end of each path
    L = length(relay); %To store length of the selected path
    relsiz(k,1) = L;
    
    
    Nframe = 1e3;
    Nbitperframe = 1e4;
    Nsamp = Nframe*Nbitperframe;
    ebn0_param = 240;
    Nsnr = numel(ebn0_param);
    Nscat = 8;
    rel_size = numel(relay);
    
    %C_mimo_cu4 = zeros(1,N);
    %C_mimo_ck4 = zeros(1,N);



    %C_mimo_cu16 = zeros(1, N);
    %C_mimo_ck16 = zeros(1,rel_size);



    %C_mimo_cu32 = zeros(1,N);
    %C_mimo_ck32 = zeros(1,N);



    %C_mimo_cu_siso = zeros(1,N);
    %C_mimo_ck_siso = zeros(1,N);
    
    for i = 1:rel_size
        if i == rel_size
            break;
        end
        
        
        txcenter = [x_tx(relay(i));y_tx(relay(i));0];
        rxcenter = [x_tx(relay(i+1));y_tx(relay(i+1));0];

        % SISO Channel

        [~,txang] = rangeangle(rxcenter,txcenter);
        [~,rxang] = rangeangle(txcenter,rxcenter);

        txsipos = [0;0;0];
        rxsopos = [0;0;0];

        % Path Loss 
        n = 2; %Path Loss Component
        grv = 3.1449;  %AWGN noise with 0 mean and 3.1449db variance
        h_dB =  (20*(log10((4*pi*fc)/(3*(10^8))))) + (10*n*log10(norm(rxcenter-txcenter))) + grv;   % Path Loss matrix in dB
        h_mag = 10^(h_dB/10);
        % h = randn + 1i*randn;
        % g = 1;  % gain for the path
        g = 1/h_mag;   
        
        
        sisochan = scatteringchanmtx(txsipos,rxsopos,txang,rxang,g);

        % Using BPSK modulation, the bit error rate (BER) for such a SISO channel can be plotted as
        Nsamp = 1e6;
        x = randi([0 1],Nsamp,1);
        
        txarray4 = phased.ULA('NumElements',Ntx4,'ElementSpacing',lambda/2);
        txmipos4 = getElementPosition(txarray4)/lambda;

        rxarray4 = phased.ULA('NumElements',Nrx4,'ElementSpacing',lambda/2);
        rxmopos4 = getElementPosition(rxarray4)/lambda;

        txarray16 = phased.ULA('NumElements',Ntx16,'ElementSpacing',lambda/2);
        txmipos16 = getElementPosition(txarray16)/lambda;

        rxarray16 = phased.ULA('NumElements',Nrx16,'ElementSpacing',lambda/2);
        rxmopos16 = getElementPosition(rxarray16)/lambda;

        txarray32 = phased.ULA('NumElements',Ntx32,'ElementSpacing',lambda/2);
        txmipos32 = getElementPosition(txarray32)/lambda;

        rxarray32 = phased.ULA('NumElements',Nrx32,'ElementSpacing',lambda/2);
        rxmopos32 = getElementPosition(rxarray32)/lambda;
        
        [txang,rxang,scatg,scatpos] = ...
        My_helperComputeRandomScatterer(txcenter,rxcenter,Nscat,fc);
%         My_helperComputeRandomScatterer(txcenter,rxcenter,Nscat,fc);
        
        % helperPlotSpatialMIMOScene(txsipos,rxsopos,...
        %     txcenter,rxcenter,scatpos);
        
        
        %%%
        %%% BER Calculation
        
        %%% SISO
        
        % SISO Multipath

        Nframe = 1e3;
        Nbitperframe = 1e3;
        Nsamp = Nframe*Nbitperframe;

        x = randi([0 1],Nbitperframe,1);
        x4 = randi([0 1],Nbitperframe,Ntx4);
        x16 = randi([0 1],Nbitperframe,Ntx16);
        x32 = randi([0 1],Nbitperframe,Ntx32);

        nerr_siso = zeros(1,Nsnr);
        nerr4 = zeros(Nrx4,Nsnr);
        nerr16 = zeros(Nrx16,Nsnr);
        nerr32 = zeros(Nrx32,Nsnr);

% % %         for m = 1:Nframe
% % %             
% % %             % SISO Multipath
% % %             
% % %             sisompchan = scatteringchanmtx(txsipos,rxsopos,txang,rxang,scatg);
% % %             wr_siso = sisompchan'/norm(sisompchan);
% % %             nerr_siso = nerr_siso+helperMIMOBER(sisompchan,x,ebn0_param,1,wr_siso);
% % %             
% % %             % 4x4 Multipath
% % %             
% % %             mimompchan4 = scatteringchanmtx(txmipos4,rxmopos4,txang,rxang,scatg);
% % %             [u4,s4,v4] = svd(mimompchan4);
% % %             wt_4 = u4(:,1)';
% % %             wr_4 = v4(:,1);
% % %             nerr4 = nerr4+helperMIMOBER(mimompchan4,x,ebn0_param,wt_4,wr_4);
% % %             
% % %             % 16x16 Multipath
% % %             
% % %             mimompchan16 = scatteringchanmtx(txmipos16,rxmopos16,txang,rxang,scatg);
% % %             [u16,s16,v16] = svd(mimompchan16);
% % %             wt16 = u16(:,1)';
% % %             wr16 = v16(:,1);
% % %             nerr16 = nerr16+helperMIMOBER(mimompchan16,x,ebn0_param,wt16,wr16);
% % %             
% % %             % 32x32 Multipath
% % %             
% % %             mimompchan32 = scatteringchanmtx(txmipos32,rxmopos32,txang,rxang,scatg);
% % %             [u32,s32,v32] = svd(mimompchan32);
% % %             wt32 = u32(:,1)';
% % %             wr32 = v32(:,1);
% % %             nerr32 = nerr32+helperMIMOBER(mimompchan32,x,ebn0_param,wt32,wr32);
% % %             
% % %         end
% % %         ber_sisomp(k,(i+1)) = nerr_siso/Nsamp;
% % %         ber_mimomp4(k,(i+1)) = nerr4/Nsamp;
% % %         ber_mimomp16(k,(i+1)) = nerr16/Nsamp;
% % %         ber_mimomp32(k,(i+1)) = nerr32/Nsamp;
        
       
        
        %%%

        C_mimo_ck4_temp = zeros(1,Nsnr);
        C_mimo_ck16_temp = zeros(1,Nsnr);
        C_mimo_ck32_temp = zeros(1,Nsnr);
        C_mimo_ck_siso_temp = zeros(1,Nsnr);
        
        Ntrial = 1000;
        
        for m = 1:Nsnr
            for n = 1:Ntrial
                
                % SISO
                
                % BER
                sisompchan = scatteringchanmtx(txsipos,rxsopos,txang,rxang,scatg);
                wr_siso = sisompchan'/norm(sisompchan);
                nerr_siso = nerr_siso+helperMIMOBER(sisompchan,x,ebn0_param,1,wr_siso);
                
                
                % Capacity
                N0 = db2pow(-ebn0_param(m));
                % [~,~,~,~,cu_siso] = diagbfweights(sisochan,1,N0,'uniform');
                [~,~,~,~,ck_siso] = diagbfweights(sisochan,1,N0,'waterfill');
                % C_mimo_cu_siso(m) = C_mimo_cu_siso(m)+cu_siso;
                C_mimo_ck_siso_temp(m) = C_mimo_ck_siso_temp(m)+ck_siso;


                % 4x4 BER
                
                mimompchan4 = scatteringchanmtx(txmipos4,rxmopos4,txang,rxang,scatg);
% % %                 [u4,s4,v4] = svd(mimompchan4);
% % %                 wt_4 = u4(:,1)';
% % %                 wr_4 = v4(:,1);
% % %                 nerr4 = nerr4+helperMIMOBER(mimompchan4,x,ebn0_param,wt_4,wr_4);
% % %                 [wp_4,wc_4] = diagbfweights(mimompchan4);
% % %                 nerr4 = nerr4+helperMIMOMultistreamBER(mimompchan4,x4,ebn0_param,wp_4,wc_4);
                
                
                
                
                % 4x4 Capacity
% % %                 mimompchan4 = scatteringchanmtx(txmipos4,rxmopos4,txang,rxang,scatg);
                N0 = db2pow(-ebn0_param(m));
                % [~,~,~,~,cu4] = diagbfweights(mimompchan4,1,N0,'uniform');
                [wp_4,wc_4,~,~,ck4] = diagbfweights(mimompchan4,1,N0,'waterfill');
                nerr4 = nerr4+helperMIMOMultistreamBER(mimompchan4,x4,ebn0_param,wp_4,wc_4);
                % C_mimo_cu4(m) = C_mimo_cu4(m)+cu4;
                C_mimo_ck4_temp(m) = C_mimo_ck4_temp(m)+ck4;


                % 16x16 BER
                
                mimompchan16 = scatteringchanmtx(txmipos16,rxmopos16,txang,rxang,scatg);
% % %                 [u16,s16,v16] = svd(mimompchan16);
% % %                 wt16 = u16(:,1)';
% % %                 wr16 = v16(:,1);
% % %                 nerr16 = nerr16+helperMIMOBER(mimompchan16,x,ebn0_param,wt16,wr16);
% % %                 [wp_16,wc_16] = diagbfweights(mimompchan16);
% % %                 nerr16 = nerr16+helperMIMOMultistreamBER(mimompchan16,x16,ebn0_param,wp_16,wc_16);
                
                
                
                
                % 16x16 Capacity
% % %                 mimompchan16 = scatteringchanmtx(txmipos16,rxmopos16,txang,rxang,scatg);
                N0 = db2pow(-ebn0_param(m));
                % [~,~,~,~,cu16] = diagbfweights(mimompchan16,1,N0,'uniform');
                [wp_16,wc_16,~,~,ck16] = diagbfweights(mimompchan16,1,N0,'waterfill');
                nerr16 = nerr16+helperMIMOMultistreamBER(mimompchan16,x16,ebn0_param,wp_16,wc_16);
                % C_mimo_cu16(m) = C_mimo_cu16(m)+cu16;
                C_mimo_ck16_temp(m) = C_mimo_ck16_temp(m)+ck16;


                % 32x32 BER
                
                mimompchan32 = scatteringchanmtx(txmipos32,rxmopos32,txang,rxang,scatg);
% % %                 [u32,s32,v32] = svd(mimompchan32);
% % %                 wt32 = u32(:,1)';
% % %                 wr32 = v32(:,1);
% % %                 nerr32 = nerr32+helperMIMOBER(mimompchan32,x,ebn0_param,wt32,wr32);
% % %                 [wp_32,wc_32] = diagbfweights(mimompchan32);
% % %                 nerr32 = nerr32+helperMIMOMultistreamBER(mimompchan32,x32,ebn0_param,wp_32,wc_32);
                
                
                
                
                
                % 32x32 Capacity
% % %                 mimompchan32 = scatteringchanmtx(txmipos32,rxmopos32,txang,rxang,scatg);
                N0 = db2pow(-ebn0_param(m));
                % [~,~,~,~,cu32] = diagbfweights(mimompchan32,1,N0,'uniform');
                [wp_32,wc_32,~,~,ck32] = diagbfweights(mimompchan32,1,N0,'waterfill');
                nerr32 = nerr32+helperMIMOMultistreamBER(mimompchan32,x32,ebn0_param,wp_32,wc_32);
                %C_mimo_cu32(m) = C_mimo_cu32(m)+cu32;
                C_mimo_ck32_temp(m) = C_mimo_ck32_temp(m)+ck32;
            end
            
            
            ber_sisomp(k,(i+1)) = nerr_siso/Nsamp;
% % %             ber_mimomp4(k,(i+1)) = nerr4/Nsamp;
            ber_mimompdiag4 = nerr4/Nsamp;
            ber_mimomp4(k,(i+1)) = ber_mimompdiag4(1,:).';
% % %             ber_mimomp16(k,(i+1)) = nerr16/Nsamp;
            ber_mimompdiag16 = nerr16/Nsamp;
            ber_mimomp16(k,(i+1)) = ber_mimompdiag16(1,:).';
% % %             ber_mimomp32(k,(i+1)) = nerr32/Nsamp;
            ber_mimompdiag32 = nerr32/Nsamp;
            ber_mimomp32(k,(i+1)) = ber_mimompdiag32(1,:).';
        
        end
        
        SE_SISO(k,i) = C_mimo_ck_siso_temp/Ntrial; % SE: 1st cell in each row gives the SE of link between 1st and 2nd node
        SINR_SISO(k,(i+1)) = (2^(SE_SISO(k,i))) - 1; % SINR: 2st cell in each row gives the Rx SINR of link between 1st and 2nd node at the 2nd node
        SINR_dB_SISO(k,(i+1)) = 10*log10(SINR_SISO(k,(i+1)));
        
        SINR_m_dB_SISO(k,(i+1)) = SINR_dB_SISO(k,(i+1)) - (-40);  % SINR threshold at Rx is -40 dB
        SE4(k,i) = C_mimo_ck4_temp/Ntrial;
        SINR4(k,(i+1)) = (2^(SE4(k,i))) - 1;
        SINR4_dB(k,(i+1)) = 10*log10(SINR4(k,(i+1)));
        SINR4_m_dB(k,(i+1)) = SINR4_dB(k,(i+1)) - (-40);
        SE16(k,i) = C_mimo_ck16_temp/Ntrial;
        SINR16(k,(i+1)) = (2^(SE16(k,i))) - 1;
        SINR16_dB(k,(i+1)) = 10*log10(SINR16(k,(i+1)));
        SINR16_m_dB(k,(i+1)) = SINR16_dB(k,(i+1)) - (-40);
        SE32(k,i) = C_mimo_ck32_temp/Ntrial;
        SINR32(k,(i+1)) = (2^(SE32(k,i))) - 1;
        SINR32_dB(k,(i+1)) = 10*log10(SINR32(k,(i+1)));
        SINR32_m_dB(k,(i+1)) = SINR32_dB(k,(i+1)) - (-40);
        
    end
    

    
    
    
    
    disp("Path No. done");
    disp(k);
    disp("Remaining paths are:")
    disp(psize - (k));
    a_time = toc;
    Act_time = a_time - b_time;
    b_time = a_time;
    disp("Time taken for this single path:")
    disp(Act_time);

    
end

time_to_cal_path_SE = toc;
disp("Time take to find SE of each link of all paths:")
disp(time_to_cal_path_SE);





for i = 1:psize
    relay = pat_mat(i,:); %To select one path from the matrix
    relay = nonzeros(relay)'; %To remove zero padding at the end of each path
    L = length(relay); %To store length of the selected path
    temp_SISO = nonzeros(SE_SISO(i,:));
    SE_SISO_path(i,:) = min(temp_SISO);
    SE_SISO_path1(i,:) = SE_SISO_path(i,:)/(L-1);
    SE_SISO_dB(i,:) = 10*log10(SE_SISO_path1(i,1));
    temp_4 = nonzeros(SE4(i,:));
    SE4_path(i,:) = min(temp_4);
    SE4_path1(i,:) = SE4_path(i,:)/(L-1);
    SE4_dB(i,:) = 10*log10(SE4_path1(i,1));
    temp_16 = nonzeros(SE16(i,:));
    SE16_path(i,:) = min(temp_16);
    SE16_path1(i,:) = SE16_path(i,:)/(L-1);
    SE16_dB(i,:) = 10*log10(SE16_path1(i,1));
    temp_32 = nonzeros(SE32(i,:));
    SE32_path(i,:) = min(temp_32);
    SE32_path1(i,:) = SE32_path(i,:)/(L-1); 
    SE32_dB(i,:) = 10*log10(SE32_path1(i,1));
    
    
    %To calculate EE of the selected path
    P_ckt = 0.010; % Circuit Power in W
    Power_tx = 0;
    for j = 1:(L-1)
    Power_tx = Power_tx + Node(relay(j)).P_tx + P_ckt;
    end
    %disp(Power_tx);
    EE_SISO(i,1) = SE_SISO_path1(i,1)/(Power_tx);
    EE_SISO_dB(i,1) = 10*log10(EE_SISO(i,1));
    EE4(i,1) = SE4_path1(i,1)/(Power_tx);
    EE4_dB(i,1) = 10*log10(EE4(i,1));
    EE16(i,1) = SE16_path1(i,1)/(Power_tx);
    EE16_dB(i,1) = 10*log10(EE16(i,1));
    EE32(i,1) = SE32_path1(i,1)/(Power_tx);
    EE32_dB(i,1) = 10*log10(EE32(i,1));
    
    
    %%% BER 
    
    % We put the max BER of a path as bottleneck of the path and put it
    % into ber_bottleneck
    ber_sisomp_btlnk(i,1) = max(ber_sisomp(i,:));  % For SISO multipath
    
    ber_mimomp4_btlnk(i,1) = max(ber_mimomp4(i,:)); % For 4x4 Multipath
    
    ber_mimomp16_btlnk(i,1) = max(ber_mimomp16(i,:)); % For 16x16 Multipath
    
    ber_mimomp32_btlnk(i,1) = max(ber_mimomp32(i,:)); % For 32x32 Multipath
   
end


%Result of EE,SE,SINR-m Vs No of intermediate relays
count = hist(relsiz,1:max(relsiz)); % Creates a vector where in each cell contains the total paths with length equal to column number
sm_SISO_mat = zeros((No_active_nodes),(psize)); %This matrix is used to store the sinr margins like ,row 3 contains the  avg sirn margins of all the paths that are of length 3 and so on... for a path with length 3 the avg sinr on that path is calculated and inserted into this matrix at row 3
sm4_mat = zeros((No_active_nodes),(psize)); %This matrix is used to store the sinr margins like ,row 3 contains the  avg sirn margins of all the paths that are of length 3 and so on... for a path with length 3 the avg sinr on that path is calculated and inserted into this matrix at row 3
sm16_mat = zeros((No_active_nodes),(psize)); %This matrix is used to store the sinr margins like ,row 3 contains the  avg sirn margins of all the paths that are of length 3 and so on... for a path with length 3 the avg sinr on that path is calculated and inserted into this matrix at row 3
sm32_mat = zeros((No_active_nodes),(psize)); %This matrix is used to store the sinr margins like ,row 3 contains the  avg sirn margins of all the paths that are of length 3 and so on... for a path with length 3 the avg sinr on that path is calculated and inserted into this matrix at row 3
se_SISO_mat = zeros((No_active_nodes),(psize)); %This matrix is used to store the SE value like ,row 3 contains the SE values of all the paths that are of length 3 and so on...
ee_SISO_mat = zeros((No_active_nodes),(psize)); %This matrix is used to store the EE value like ,row 3 contains the EE values of all the paths that are of length 3 and so on...
se4_mat = zeros((No_active_nodes),(psize));
ee4_mat = zeros((No_active_nodes),(psize));
se16_mat = zeros((No_active_nodes),(psize));
ee16_mat = zeros((No_active_nodes),(psize));
se32_mat = zeros((No_active_nodes),(psize));
ee32_mat = zeros((No_active_nodes),(psize));
ber_SISO_mat = nan((No_active_nodes),(psize)); % %This matrix is used to store the BER of SISO value like ,row 3 contains the BER  of SISO values of all the paths that are of length 3 and so on...
ber4_mat = nan((No_active_nodes),(psize));
ber16_mat = nan((No_active_nodes),(psize));
ber32_mat = nan((No_active_nodes),(psize));


for i = 1:psize
    relay = pat_mat(i,:); %To select one path from the matrix
    relay = nonzeros(relay)'; %To remove zero padding at the end of each path
    L = length(relay); %To store length of the selected path
    
    
    %sm_mat will contain many zeros in a row and the sirn margins will be
    %inserted at any(actuall ith column) column since I can't figure out
    %how to insert it into continuous columns
    
    %SINR_m_dB contains many NaNs in its rows, so we first select ith row
    %remove its NaNs the calulate its mean then insert it into the row of
    %sm_mat whre row number equals length of path
    
    % SINRm For SISO
    temp_row = SINR_m_dB_SISO(i,:); % Selecting ith row of the SINR_m_dB matirx
    
    temp_row = temp_row';% converting to column vector
    temp_row = temp_row(~isnan(temp_row))'; %Removing NaNs from vector
 
    sm_SISO_mat(L,i) = mean(temp_row); %Calculating the mean sinr margin of a path then inserting into the row number which is equal to the length of the path
    
    % SINRm for 4x4
    
    temp_row1 = SINR4_m_dB(i,:); % Selecting ith row of the SINR_m_dB matirx
    
    temp_row1 = temp_row1';% converting to column vector
    temp_row1 = temp_row1(~isnan(temp_row1))'; %Removing NaNs from vector
 
    sm4_mat(L,i) = mean(temp_row1); %Calculating the mean sinr margin of a path then inserting into the row number which is equal to the length of the path
    
    % SINRm for 16x16
    
    temp_row2 = SINR16_m_dB(i,:); % Selecting ith row of the SINR_m_dB matirx
    
    temp_row2 = temp_row2';% converting to column vector
    temp_row2 = temp_row2(~isnan(temp_row2))'; %Removing NaNs from vector
 
    sm16_mat(L,i) = mean(temp_row2); %Calculating the mean sinr margin of a path then inserting into the row number which is equal to the length of the path

    % SINRm for 32x32
    
    temp_row3 = SINR32_m_dB(i,:); % Selecting ith row of the SINR_m_dB matirx
    
    temp_row3 = temp_row3';% converting to column vector
    temp_row3 = temp_row3(~isnan(temp_row3))'; %Removing NaNs from vector
 
    sm32_mat(L,i) = mean(temp_row3); %Calculating the mean sinr margin of a path then inserting into the row number which is equal to the length of the path

    
    
    se_SISO_mat(L,i) = SE_SISO_dB(i,1); %Getting the SE of SISO of the path then inserting into the row number which is equal to the length of the path
    ee_SISO_mat(L,i) = EE_SISO_dB(i,1); %Getting the EE of SISO of the path then inserting into the row number which is equal to the length of the path
    se4_mat(L,i) = SE4_dB(i,1);
    ee4_mat(L,i) = EE4_dB(i,1);
    se16_mat(L,i) = SE16_dB(i,1);
    ee16_mat(L,i) = EE16_dB(i,1);
    se32_mat(L,i) = SE32_dB(i,1);
    ee32_mat(L,i) = EE32_dB(i,1);
    
    
    % For BER
    ber_SISO_mat(L,i) = ber_sisomp_btlnk(i,1); % Gets the bottleneck BER of SISO and then places it in the ber_SISO_mat, at the row equal to the length of the path
    ber4_mat(L,i) = ber_mimomp4_btlnk(i,1);
    ber16_mat(L,i) =  ber_mimomp16_btlnk(i,1);
    ber32_mat(L,i) =  ber_mimomp32_btlnk(i,1);
    
end

%%% Removing NaNs from BER matrices and getting mean of non NaNs
avg_ber_SISO = zeros(No_active_nodes,1);
avg_ber4 = zeros(No_active_nodes,1);
avg_ber16 = zeros(No_active_nodes,1);
avg_ber32 = zeros(No_active_nodes,1);

for i = 1:No_active_nodes
    
    % For SISO
    
    temp_row_BER_SISO = ber_SISO_mat(i,:); % Selecting ith row of the ber_mat matirx
    temp_row_BER_SISO = temp_row_BER_SISO'; % converting to column vector
    temp_row_BER_SISO = temp_row_BER_SISO(~isnan(temp_row_BER_SISO))'; %Removing NaNs from vector
    avg_ber_SISO(i,1) = mean(temp_row_BER_SISO); % Getting the mean of all non nans and storing in avg_ber
    
    % For 4x4
    
    temp_row_BER4 = ber4_mat(i,:); % Selecting ith row of the ber_mat matirx
    temp_row_BER4 = temp_row_BER4'; % converting to column vector
    temp_row_BER4 = temp_row_BER4(~isnan(temp_row_BER4))'; %Removing NaNs from vector
    avg_ber4(i,1) = mean(temp_row_BER4); % Getting the mean of all non nans and storing in avg_ber
    
    % For 16x16
    
    temp_row_BER16 = ber16_mat(i,:); % Selecting ith row of the ber_mat matirx
    temp_row_BER16 = temp_row_BER16'; % converting to column vector
    temp_row_BER16 = temp_row_BER16(~isnan(temp_row_BER16))'; %Removing NaNs from vector
    avg_ber16(i,1) = mean(temp_row_BER16); % Getting the mean of all non nans and storing in avg_ber
    
    % For 32x32
    
    temp_row_BER32 = ber32_mat(i,:); % Selecting ith row of the ber_mat matirx
    temp_row_BER32 = temp_row_BER32'; % converting to column vector
    temp_row_BER32 = temp_row_BER32(~isnan(temp_row_BER32))'; %Removing NaNs from vector
    avg_ber32(i,1) = mean(temp_row_BER32); % Getting the mean of all non nans and storing in avg_ber
    
    
    
end



%Creating avg_sm, avg_se, avg_ee of Lx1 dim, where each row of the vector holds the avg of all the nonzero values of the corresponding row, of sm_mat, se_mat and ee_mat respectively.
avg_sm_SISO = zeros(No_active_nodes,1);
avg_se_SISO = zeros(No_active_nodes,1);
avg_ee_SISO = zeros(No_active_nodes,1);
avg_sm4 = zeros(No_active_nodes,1);
avg_se4 = zeros(No_active_nodes,1);
avg_ee4 = zeros(No_active_nodes,1);
avg_sm16 = zeros(No_active_nodes,1);
avg_se16 = zeros(No_active_nodes,1);
avg_ee16 = zeros(No_active_nodes,1);
avg_sm32 = zeros(No_active_nodes,1);
avg_se32 = zeros(No_active_nodes,1);
avg_ee32 = zeros(No_active_nodes,1);


for i = 2:No_active_nodes
    avg_sm_SISO(i,1) = mean(nonzeros(sm_SISO_mat(i,:)));
    avg_se_SISO(i,1) = mean(nonzeros(se_SISO_mat(i,:)));
    avg_ee_SISO(i,1) = mean(nonzeros(ee_SISO_mat(i,:)));
    avg_sm4(i,1) = mean(nonzeros(sm4_mat(i,:)));
    avg_se4(i,1) = mean(nonzeros(se4_mat(i,:)));
    avg_ee4(i,1) = mean(nonzeros(ee4_mat(i,:)));
    avg_sm16(i,1) = mean(nonzeros(sm16_mat(i,:)));
    avg_se16(i,1) = mean(nonzeros(se16_mat(i,:)));
    avg_ee16(i,1) = mean(nonzeros(ee16_mat(i,:)));
    avg_sm32(i,1) = mean(nonzeros(sm32_mat(i,:)));
    avg_se32(i,1) = mean(nonzeros(se32_mat(i,:)));
    avg_ee32(i,1) = mean(nonzeros(ee32_mat(i,:)));
    
end
%%

% Plot of Avg SE,EE of path with given No. of intermediate Nodes Vs No. of intermediate relay nodes

figure;
hold on;
grid on;
title('Avg SE,EE of path with given No. of intermediate Nodes Vs No. of intermediate relay nodes');
xlabel('No. of Intermediate nodes in a relay path');
ylabel('Avg. SE,EE in dB');
xaxis = (0:1:(No_active_nodes-2)); 
plot(xaxis,avg_se_SISO(2:end),'bo-','LineWidth',2);
plot(xaxis,avg_ee_SISO(2:end),'rx-','LineWidth',2);
plot(xaxis,avg_se4(2:end),'g+:','LineWidth',2);
plot(xaxis,avg_ee4(2:end),'m*:','LineWidth',2);
plot(xaxis,avg_se16(2:end),'ks--','LineWidth',2);
plot(xaxis,avg_ee16(2:end),'cd--','LineWidth',2);
plot(xaxis,avg_se32(2:end),'r^-.','LineWidth',2);
plot(xaxis,avg_ee32(2:end),'bv-.','LineWidth',2);

%plot(linspace(0,No_active_nodes,No_active_nodes-1),avg_sm(2:end),'k*-');
legend('SE SISO','EE SISO', 'SE 4x4', 'EE 4x4', 'SE 16x16', 'EE 16x16', 'SE 32x32', 'EE 32x32');
hold off;

% Plot of Avg SINR-m of path with given No. of intermediate Nodes Vs No. of intermediate relay nodes

figure;
hold on;
grid on;
title('Avg SINR-m of path with given No. of intermediate Nodes Vs No. of intermediate relay nodes');
xlabel('No. of Intermediate nodes in a relay path');
ylabel('SINR-m in dB');
plot(xaxis,avg_sm_SISO(2:end),'bo-');
plot(xaxis,avg_sm4(2:end),'rx-');
plot(xaxis,avg_sm16(2:end),'m+-');
plot(xaxis,avg_sm32(2:end),'k*-');
legend('SISO', '4x4', '16x16', '32x32');
hold off;

% Plot of Avg BER Vs No. of intermediate relay nodes
figure;
hold on;
grid on;
title('Avg BER Vs No. of intermediate relay nodes')
xlabel('No. of Intermediate nodes in a relay path');
ylabel('BER');
plot(xaxis,avg_ber_SISO(2:end),'bo-');
plot(xaxis,avg_ber4(2:end),'rx-');
plot(xaxis,avg_ber16(2:end),'m+-');
plot(xaxis,avg_ber32(2:end),'k*-');
legend('SISO', '4x4', '16x16', '32x32');
hold off;