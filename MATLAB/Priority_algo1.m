clc
clear
close all

%rng(3);
tic;
R1 = 0; % Avg no. of nodes for condition 1
R10 = 0; % No. of times there are zero nodes for condition 1
R2 = 0; % Avg no. of nodes for condition 2
R20 = 0; % No. of times there are zero nodes for condition 2
R3 = 0; % Avg no. of nodes for condition 3
R30 = 0; % No. of times there are zero nodes for condition 3
count_err = 0; % No. of times the half way correction is applied
count_err1 = 0; % No. of times the iteration is skipped due to spectrum crowding
t1_waiting = 0; % Wait time for service 1
t2_waiting = 0; % Wait time for service 2
t3_waiting = 0; % Wait time for service 3


Avg_SE_P1 = 0; % Holds the Avg. SE of P1
Avg_SE_P2 = 0; % Holds the Avg. SE of P2
Avg_SE_P3 = 0; % Holds the Avg. SE of P3
Avg_SE_P1_dB = 0; % Holds the Avg. SE in dB of P1
Avg_SE_P2_dB = 0; % Holds the Avg. SE in dB of P2
Avg_SE_P3_dB = 0; % Holds the Avg. SE in dB of P3

Avg_EE_P1 = 0; % Holds the Avg. EE of P1
Avg_EE_P2 = 0; % Holds the Avg. EE of P2
Avg_EE_P3 = 0; % Holds the Avg. EE of P3
Avg_EE_P1_dB = 0; % Holds the Avg. EE in dB of P1
Avg_EE_P2_dB = 0; % Holds the Avg. EE in dB of P2
Avg_EE_P3_dB = 0; % Holds the Avg. EE in dB of P3


try_catch1 = 0; % Error indicator if the number of already Tx nodes is greater then the avaliable nodes for reception 
try_catch2 = 0;
try_catch3 = 0;

NP1 = 0; % No. of time P1 occurred
NP2 = 0; % No. of time P2 occurred
NP3 = 0; % No. of time P3 occurred

wait_err_in_P1 = 0; % Holds the number of times the wait exceeds 20 for P1 service
wait_err_in_P2 = 0; % Holds the number of times the wait exceeds 20 for P2 service
wait_err_in_P3 = 0; % Holds the number of times the wait exceeds 20 for P3 service

tic;
itrations = 15000;
% rng_seed_val = 77717;
% writematrix(rng_seed_val, "Rng_seed_val.txt")
for ind = 1 :itrations
rng(37*ind);
disp("Iteration no: \n");
disp(ind);
%tic2 = tic;
% rng(ind*rng_seed_val-5+ind);
t1_wait = 0; % Temp Wait time for service 1
t2_wait = 0; % Temp Wait time for service 2
t3_wait = 0; % Temp Wait time for service 3
done = 0; % Indicator of whether the service is completed or not, i.e, 0 = not done and 1 = done
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
% rbi_gen = abs(normrnd(0.7,0.4,[1,N]));
% rbi_gen = abs(normrnd(0.7,0.5,[1,N]));
% rbi_gen = abs(normrnd(0.7,0.45,[1,N]));
rbi_gen = abs(normrnd(0.8,0.3,[1,N]));
for i = 1:length(rbi_gen)
    if rbi_gen(i) > 1
%         rbi_gen(i) = (randi([3,10]))/10;
          rbi_gen(i) = (randi([7,9]))/10;
    end
end

%Initialising Node struct
Node = repmat(struct('node_id',0, 'rbi', 0, 'x_cor', 0, 'y_cor', 0, 'status', 0, 'service_status', 0,'time_stamp', 0, 'service', 0, 'Adjacency_vector' , zeros(1,N)),N,1);
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

Wh_counter = 0; % Tells if the while loop is executing more than one, implies waiting occured. If it is set to greater than 0 then waiting has occured

while done == 0
    if Wh_counter == 0
    
    


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

        
%         norand = abs(normrnd(0, 0.9, [1,No_active_nodes]));
        norand = abs(normrnd(0, 1, [1,No_active_nodes]));
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
                 if countcng >= (No_active_nodes/2) - 4
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
                        Node(i).time_stamp = Wh_counter;
                    elseif rn > 0.5 && rn <= 0.8
                        Node(i).service = 2;
                        Node(i).time_stamp = Wh_counter;
                    else 
                        Node(i).service = 3;
                        Node(i).time_stamp = Wh_counter;
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
        if Already_tx_count >= collec_arr_len
            disp(ind);
            disp("This iteration has already Tx nodes greater than nodes avaliable for reception");
            skip_itr2 = skip_itr2 + 1;
            t1_wait = 0;
            t2_wait = 0;
            t3_wait = 0;
            max_SE_P1_dB = 0;
            max_SE_P2_dB = 0;
            max_SE_P3_dB = 0;
            max_EE_P1_dB = 0;
            max_EE_P2_dB = 0;
            max_EE_P3_dB = 0;
            break;
%             continue;
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
            if (Node(i).status == 1 && Node(i).service_status == 0)  % If the node is active and idle
                ser_adj(i,:) = ones(1,N);
                ser_adj(i,i) = 0;
            end
        end

        for i = 1:N
            if Node(i).status == 0 || (Node(i).service_status == 1)  || (Node(i).service_status == 2) % If node is asleep or Tx or Rx, making ser_adj such that it does not participate
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
            t1_wait = 0;
            t2_wait = 0;
            t3_wait = 0;
            max_SE_P1_dB = 0;
            max_SE_P2_dB = 0;
            max_SE_P3_dB = 0;
            max_EE_P1_dB = 0;
            max_EE_P2_dB = 0;
            max_EE_P3_dB = 0;
            break;
%             continue;
        end


        % Plot of all avaliable paths
        grap_paths = digraph(ser_adj);
% % %         figure;
% % %         plot(grap_paths);
% % %         title('All avaliable paths for new transmission');

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
            count_err1 = count_err1 + 1;
            Wh_counter = Wh_counter + 1;
            continue;
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

        Src_Des_info = zeros(N,N);
        Src_Des_info(Src,Des) = 1;

        Node(Src).service_status = 1; % Setting the servive status of source to transmitting service
        Node(Des).service_status = 2; % Setting the servive status of destination to receiving service
        new_service = rand; % Randomly generating new service with 50% P1, 30% P2 and 20 % P3 chance.
        if new_service <= (1/3)
            Node(Src).service = 1;
            NP1 = NP1 + 1;
            disp("The new service is P1");
        elseif new_service > (1/3) && new_service <= (2/3)
            Node(Src).service = 2;
            disp("The new service is P2");
            NP2 = NP2 + 1;
        else 
            Node(Src).service = 3;
            disp("The new service is P3");
            NP3 = NP3 + 1;
        end

       
    else       % If Wh_counter is greater than 0
        
        % We assume it takes 1.5s, 1s and 0.5s to complete P1, P2 and P3
        % services respectively
        
        for i = 1:N
            if Node(i).service_status == 1 && i ~= Src
                if Node(i).service == 1
                    if (Wh_counter - Node(i).time_stamp) >= 1.5
                        Node(i).service_status = 0;
                        Node(i).service = 0;
                        for j = 1:N
                            if spec_mat(i,j) ~= 0
                                Node(j).service_status = 0;
                                Node(j).service = 0;
                            end
                        end
                    end
                elseif Node(i).service == 2
                    if (Wh_counter - Node(i).time_stamp) >= 1
                        Node(i).service_status = 0;
                        Node(i).service = 0;
                        for j = 1:N
                            if spec_mat(i,j) ~= 0
                                Node(j).service_status = 0;
                                Node(j).service = 0;
                            end
                        end
                    end
                elseif Node(i).service == 3
                    if (Wh_counter - Node(i).time_stamp) >= 0.5
                        Node(i).service_status = 0;
                        Node(i).service = 0;
                        for j = 1:N
                            if spec_mat(i,j) ~= 0
                                Node(j).service_status = 0;
                                Node(j).service = 0;
                            end
                        end
                    end    
                end
            end
        end
        
        temp_idle_nodes = 0; 
        for i = 1:N
            if Node(i).status ~= 0 && Node(i).service_status == 0
                temp_idle_nodes = temp_idle_nodes + 1;
            end
        end
        
        if temp_idle_nodes > 14
            norand_var = 1;
        else
            norand_var = exp(-Wh_counter/1.1);
        end
        norand = abs(normrnd(0, norand_var, [1,temp_idle_nodes]));
        countnor = 0;
        for i = 1:temp_idle_nodes
            if norand(i) >= 1
               countnor = countnor + 1; 
            end
        end
        
        count_Tx = 0;  % This will hold the no. of Tx nodes from previous while loop iteration
        for i = 1:N
            if i ~=Src || i ~= Des
                if Node(i).service == 1
                    count_Tx = count_Tx + 1;
                end
            end
        end
        
        
        countcng = countnor; 
        if (countnor + count_Tx) >= (temp_idle_nodes/2) - 1 
           disp("Error: The no. of services to be simulated is greater than half, applying correction")
           count_err = count_err + 1;
           for i = 1:temp_idle_nodes
              if norand(i) >= 1
                 if countcng >= (temp_idle_nodes/2) - 4
                    norand(i) = rand;
                    countcng = countcng - 1;
                 end
              end
           end
        end
        
        % Correction is done

        j = 1;
        Already_tx_count = 0;
        for i = 1:N
            if Node(i).status == 1 && Node(i).service_status == 0
                if norand(j) >= 1 
                    Already_tx_count = Already_tx_count +1; % counts the number of already Tx nodes
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
        %Already_tx_count = 0;
        
        for i = 1:N
            if Node(i).status == 1 && Node(i).service_status == 0 % If node is active and idle
                collec_arr(1,i) = Node(i).node_id; 
            end
    %         if Node(i).status == 1 && Node(i).service_status == 1 % If node is active and Tx
    %              Already_tx_count = Already_tx_count +1; % counts the number of already Tx nodes 
    %         end
        end
        
        collec_arr1 = nonzeros(collec_arr);
        collec_arr_len = length(collec_arr1);
        
        count_Tx_combined = 0;  % This will hold the no. of Tx nodes from previous and current while loop combined, whereas Already
        for i = 1:N
            if i ~=Src || i ~= Des
                if Node(i).service == 1
                    count_Tx_combined = count_Tx_combined + 1;
                end
            end
        end

        skip_itr2 = 0; % skip iteration count if number of already Tx nodes is greater then the avaliable nodes for reception
        if count_Tx_combined >= collec_arr_len

            disp("This iteration has already Tx nodes greater than nodes avaliable for reception");
            disp(ind);
            skip_itr2 = skip_itr2 + 1;
            t1_wait = 0;
            t2_wait = 0;
            t3_wait = 0;
            max_SE_P1_dB = 0;
            max_SE_P2_dB = 0;
            max_SE_P3_dB = 0;
            max_EE_P1_dB = 0;
            max_EE_P2_dB = 0;
            max_EE_P3_dB = 0;
            break;
%             continue;
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
            if Node(i).status == 1 && Node(i).service_status == 0 || (i == Src || i == Des) % If the node is active and idle or Src or Des
                ser_adj(i,:) = ones(1,N);
                ser_adj(i,i) = 0;
            end
        end

        for i = 1:N
            if (Node(i).status == 0 || (Node(i).service_status == 1)  || (Node(i).service_status == 2)) && i ~= Src && i ~=Des % If node is asleep or Tx or Rx, making ser_adj such that it does not participate
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
            t1_wait = 0;
            t2_wait = 0;
            t3_wait = 0;
            max_SE_P1_dB = 0;
            max_SE_P2_dB = 0;
            max_SE_P3_dB = 0;
            max_EE_P1_dB = 0;
            max_EE_P2_dB = 0;
            max_EE_P3_dB = 0;
            break;
%             continue;
        end
        
        % Plot of all avaliable paths
        grap_paths = digraph(ser_adj);
% % %         figure;
% % %         plot(grap_paths);
% % %         title('All avaliable paths for new transmission');

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
            count_err1 = count_err1 + 1;
            Wh_counter = Wh_counter + 1;
            continue;
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
        
    end
    
    
    
    if No_active_avaliable_nodes > 12
        paths = allpaths(grap_paths,Src,Des, 'MaxPathLength', 6);
        disp("No_active_avaliable_nodes > 12");
        disp(No_active_avaliable_nodes);
    else
        paths = allpaths(grap_paths,Src,Des, 'MaxPathLength', 6);
        disp("No_active_avaliable_nodes < 12");
        disp(No_active_avaliable_nodes);
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
    Avg_rbi = zeros(psize,1); % Store the avg rbi of the paths
    SE = zeros(psize,1);
    SE_dB = zeros(psize,1);
    EE = zeros(psize,1);
    EE_dB = zeros(psize,1);
    for i = 1:psize
        temp_rbi = 0;
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

        % To calculate the Avg. RBI of the path
        for j = 1:(L-1)
            temp_rbi = temp_rbi + Node(relay(j)).rbi; 
        end

        Avg_rbi(i,1) = temp_rbi/(L-1);

    end

    %Calculating the SINR margin matrix in dB
    for i = 1:psize
        for j = 2:No_active_nodes
            SINR_m_dB(i,j) = SINR_dB(i,j) - (-110);
        end
    end

  
    % [max_SE, SE_ind] = max(SE_dB);
    % [max_EE, EE_ind] = max(EE_dB);
    
    max_SE_P1_dB = 0;
    SE_ind_P1 = 0;
    
    max_SE_P2_dB = 0;
    SE_ind_P2 = 0;
    
    max_SE_P3_dB = 0;
    SE_ind_P3 = 0;
    
    
    max_EE_P1_dB = 0;
    EE_ind_P1 = 0;
    
    max_EE_P2_dB = 0;
    EE_ind_P2 = 0;
    
    max_EE_P3_dB = 0;
    EE_ind_P3 = 0;
    

    if Node(Src).service == 1 % If the new service is priority 1 service
        for i = 1:psize
           if Avg_rbi(i,1) < 0.7
                SE(i,1) = nan;
                SE_dB(i,1) = nan;
                EE(i,1) = nan;
                EE_dB(i,1) = nan;
            end
        end

        if any(SE_dB)
            if NSA >= 0.7
                [max_SE_P1_dB, SE_ind_P1] = max(SE_dB);
                [max_EE_P1_dB] = EE_dB(SE_ind_P1);
                deploy_path = cell2mat(paths(SE_ind_P1));
                done = 1;
                disp("Service P1 deployed on cond. 1");
            elseif NSA < 0.7 && NSA >= 0.5
                [max_SE_P1_dB, SE_ind_P1] = max(SE_dB);
                [max_EE_P1_dB] = EE_dB(SE_ind_P1);
                deploy_path = cell2mat(paths(SE_ind_P1));
                done = 1;
                disp("Service P1 deployed on cond. 2");
            else
                t1_wait = t1_wait + 1;
                Wh_counter = Wh_counter + 1;
                disp("Waiting occurred for P1");
                disp(Wh_counter);
            end
        else
            t1_wait = t1_wait + 1;
            Wh_counter = Wh_counter + 1;
            disp("Waiting occurred for P1");
            disp(Wh_counter)
        end

    elseif Node(Src).service == 2 % If the new service is priority 2 service
        for i = 1:psize
            if Avg_rbi(i,1) >= 0.7 || Avg_rbi(i,1) < 0.5
                SE(i,1) = nan;
                SE_dB(i,1) = nan;
                EE(i,1) = nan;
                EE_dB(i,1) = nan;
            end
        end

        if any(SE_dB)
            if NSA >= 0.5 && NSA < 0.7
                [max_SE_P2_dB, SE_ind_P2] = max(SE_dB);
                [max_EE_P2_dB] = EE_dB(SE_ind_P2);
                deploy_path = cell2mat(paths(SE_ind_P2));
                done = 1;
                disp("Service P2 deployed on cond. 1");
            elseif NSA >= 0.5
                [max_SE_P2_dB, SE_ind_P2] = max(SE_dB);
                [max_EE_P2_dB] = EE_dB(SE_ind_P2);
                deploy_path = cell2mat(paths(SE_ind_P2));
                done = 1;
                disp("Service P2 deployed on cond. 2");
            else
                t2_wait = t2_wait + 1;
                Wh_counter = Wh_counter + 1;
                disp("Waiting occurred for P2");
                disp(Wh_counter)
             end
        else
            t2_wait = t2_wait + 1;
            Wh_counter = Wh_counter + 1;
            disp("Waiting occurred for P2");
            disp(Wh_counter)
        end

    elseif Node(Src).service == 3 % If the new service is priority 3 service
        for i = 1:psize
            if t3_wait < 11
                if Avg_rbi(i,1) >= 0.5 || Avg_rbi(i,1) < 0.3
                    SE(i,1) = nan;
                    SE_dB(i,1) = nan;
                    EE(i,1) = nan;
                    EE_dB(i,1) = nan;
                end
            else
                if Avg_rbi(i,1) >= 0.7 || Avg_rbi(i,1) < 0.3
                    SE(i,1) = nan;
                    SE_dB(i,1) = nan;
                    EE(i,1) = nan;
                    EE_dB(i,1) = nan;
                end
            end
        end

        if any(SE_dB)
            if NSA >= 0.5
                [max_SE_P3_dB, SE_ind_P3] = max(SE_dB);
                [max_EE_P3_dB] = EE_dB(SE_ind_P3);
                deploy_path = cell2mat(paths(SE_ind_P3));
                done = 1;
                disp("Service P3 deployed on cond. 1");
            elseif NSA < 0.5
                [max_SE_P3_dB, SE_ind_P3] = max(SE_dB);
                [max_EE_P3_dB] = EE_dB(SE_ind_P3);
                deploy_path = cell2mat(paths(SE_ind_P3));
                done = 1;
                disp("Service P3 deployed on cond. 2");
            else
                t3_wait = t3_wait + 1;
                Wh_counter = Wh_counter + 1;
                disp("Waiting occurred for P3");
                disp(Wh_counter)
            end
        else
            t3_wait = t3_wait + 1;
            Wh_counter = Wh_counter + 1;
            disp("Waiting occurred for P3");
            disp(Wh_counter)
        end
    end
    
    if Wh_counter > 20
        if t1_wait ~= 0 && NP1 ~= 0
            NP1 = NP1 - 1;
            wait_err_in_P1 = wait_err_in_P1 + 1;
            disp("Wait error occurred in P1");
        elseif t2_wait ~= 0 && NP2 ~= 0
            NP2 = NP2 - 1;
            wait_err_in_P2 = wait_err_in_P2 + 1;
            disp("Wait error occurred in P2");
        elseif t3_wait ~= 0 && NP3 ~= 0
            NP3 = NP3 - 1;
            wait_err_in_P3 = wait_err_in_P3 + 1;
            disp("Wait error occurred in P3");
        end
        t1_wait = 0;
        t2_wait = 0;
        t3_wait = 0;
        max_SE_P1_dB = 0;
        max_SE_P2_dB = 0;
        max_SE_P3_dB = 0;
        max_EE_P1_dB = 0;
        max_EE_P2_dB = 0;
        max_EE_P3_dB = 0;
    break; 
    end
end

t1_waiting = t1_waiting + t1_wait;
t2_waiting = t2_waiting + t2_wait;
t3_waiting = t3_waiting + t3_wait;

Avg_SE_P1_dB = Avg_SE_P1_dB + max_SE_P1_dB;
Avg_SE_P2_dB = Avg_SE_P2_dB + max_SE_P2_dB;
Avg_SE_P3_dB = Avg_SE_P3_dB + max_SE_P3_dB;

Avg_EE_P1_dB = Avg_EE_P1_dB + max_EE_P1_dB;
Avg_EE_P2_dB = Avg_EE_P2_dB + max_EE_P2_dB;
Avg_EE_P3_dB = Avg_EE_P3_dB + max_EE_P3_dB;

end

%%

Avg_SE_P1_dB = Avg_SE_P1_dB/NP1;
Avg_EE_P1_dB = Avg_EE_P1_dB/NP1;

Avg_SE_P2_dB = Avg_SE_P2_dB/NP2;
Avg_EE_P2_dB = Avg_EE_P2_dB/NP2;

Avg_SE_P3_dB = Avg_SE_P3_dB/NP3;
Avg_EE_P3_dB = Avg_EE_P3_dB/NP3;

disp("No. of times P1 occurred");
disp(NP1);
disp("No. of times P2 occurred");
disp(NP2);
disp("No. of times P3 occurred");
disp(NP3);

xaxis = [1, 2, 3];
waiting_time = [t1_waiting/NP1, t2_waiting/NP2, t3_waiting/NP3];
Av_SE = [Avg_SE_P1_dB, Avg_SE_P2_dB, Avg_SE_P3_dB];
Av_EE = [Avg_EE_P1_dB, Avg_EE_P2_dB, Avg_EE_P3_dB];

figure;

yyaxis left
b = bar(xaxis,waiting_time);
yyaxis right
p_se = plot(xaxis,Av_SE);
p_se.Marker = 'o';
p_se.MarkerEdgeColor = 'r';
p_se.Color = 'r';
p_se.LineWidth = 3;
hold on;

p_ee = plot(xaxis,Av_EE);
p_ee.Marker = '*';
p_ee.MarkerEdgeColor = 'k';
p_ee.Color = 'k';
p_ee.LineWidth = 3;
hold on;


title('Waiting Time, Avg. SE and Avg. EE for different services priorities')
xlabel('Service Priorities')
yyaxis left
ylabel('Waiting Time in seconds')
yyaxis right
ylabel('dB')
legend('Wait Time','SE','EE');
toc;