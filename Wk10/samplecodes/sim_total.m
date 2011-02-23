%Activity is randomised in the Azimuthal (horizontal, j, x) direction.
function [Elev_phase, Elev_amp,Azim_phase, Azim_amp, connection_matrix,N] = sim_total(azim_ord,elev_ord,azim_spread,elev_spread)

% x = j and y = i for images,matrices

bar_width = 5; %number of pixels (width) that will be active at once
N = 250; %for NxN matrix
rad = 3; %number of pixels above and below in field that connect to a certain pixel in colliculus
num_connections = (2*rad+1)^2;
Elev_amp = zeros(N); %initialising
Elev_phase = zeros(N); %initialising
Azim_amp = zeros(N); %initialising
Azim_phase = zeros(N); %initialising
Coll_mixed = zeros(N); %initialising

sparse_i = zeros((2*rad+1)^2*N^2,1);
sparse_j = zeros((2*rad+1)^2*N^2,1);
sparse_s = ones((2*rad+1)^2*N^2,1);
disordered_azim = 0;
disordered_elev = 0;

%set whether x/y are ordered or disordered
if azim_ord == 'o'|| azim_ord == 'O'|| azim_ord == 0
    order_azim = 1:N;
else
    disordered_azim = 1;
end

if elev_ord == 'o'|| elev_ord == 'O'|| elev_ord == 0
    order_elev = 1:N;
else
    disordered_elev = 1;
end

%loop for creating connection matrix
for i = 1:N
    if disordered_elev == 0
        loc_i = order_elev(i);
        imin = max(loc_i-elev_spread,1);
        imax = min(loc_i+elev_spread,N);
        num_rows = length(imin:imax);
    end
    for j = 1:N %for each pixel
        if disordered_elev == 1
            loc_i = ceil(rand*N);
            imin = max(loc_i-elev_spread,1);
            imax = min(loc_i+elev_spread,N);
            num_rows = length(imin:imax);
        end
        
        if disordered_azim == 1
            loc_j = ceil(rand*N);
        else
            loc_j = order_azim(j);
        end
        jmin = max(loc_j-azim_spread,1);
        jmax = min(loc_j+azim_spread,N);
        num_connex = length(jmin:jmax); %find how far in j direction connections spread (will only differ for edges)
        %poss_conne_x = randperm(num_connex)+(jmin);
        %poss_conne_x = intersect(poss_connex,1:N);
        %poss_conne_y = randperm(num_rows)+(imin);
        %poss_conne_y = intersect(poss_connex,1:N);
        chosen_connex = zeros(num_connections,2);
        neuron_number = sub2ind([N,N],i,j);
        for k = 1:num_connections
            xx = ceil(rand*num_connex);
            yy = ceil(rand*num_rows);
            chosen_connex(k,1) = neuron_number;
            chosen_connex(k,2) = sub2ind([N,N],xx+jmin-1,yy+imin-1); %switched x,y to reflect mapping
        end
        
        
            
        num_added = length(chosen_connex);
        last_conn = find(sparse_i==0,1);
        sparse_i(last_conn:last_conn+num_added-1) =  chosen_connex(:,1);%add connection from neuron
        sparse_j(last_conn:last_conn+num_added-1) =  chosen_connex(:,2);%add connection to neuron
    end        
    i
end
%last_conn = find(sparse_i==0,1);
%sparse_i = sparse_i(1:last_conn-1);% get rid of extra zeros at the end from the initialisation
%sparse_j = sparse_j(1:last_conn-1);% get rid of extra zeros
%sparse_s = sparse_s(1:last_conn-1);
connection_matrix = sparse(sparse_i,sparse_j,sparse_s);%create sparse matrix                    

%Elevation direction main loop
%%For each row of pixels assume that is the bottom of the 'barwidth' pixel bar and
%%'stimulate' all neurons covered by the bar. 

for i = 1:N
    E_amp = zeros(N); %temp activity initialisation
    imin = max(i-bar_width+1,1);
    imax = i;
    num_rows = length(imin:imax);%how many rows are covered by the bar
    y_coord = repmat(imin:imax,1,N);
    y_coord = sort(y_coord); %create x coords for stimulated pixels
    x_coord = repmat(1:N,1,num_rows); %create y coords of stimulated pixels
    field_neurons = sub2ind([N,N],y_coord,x_coord); %assign x,y to neuron number in connection matrix
    activity = sum(connection_matrix(field_neurons,:)); %sum over COLL neurons of matrix only including active FIELD neurons
    active_neurons = find(activity); %find which COLL neurons have activity
    E_amp(active_neurons) = activity(active_neurons); %set temp activity to its current value for active neurons
    imagesc(E_amp);
    drawnow
    phase_max = find(E_amp>Elev_amp); %find which pixels have higher activity than ever before
    Elev_phase(phase_max) = i; %set their maximum phase to current row
    Elev_amp = max(Elev_amp,E_amp); % set max amp to the max between old max and current level
end
%loop for azimuthal direction
for j = 1:N
    A_amp = zeros(N); %temp activity
    jmin = max(j-bar_width+1,1);
    jmax = j;
    num_col = length(jmin:jmax);
    x_coord = repmat(jmin:jmax,1,N);
    x_coord = sort(x_coord);
    y_coord = repmat(1:N,1,num_col);
    field_neurons = sub2ind([N,N],y_coord,x_coord);
    activity = sum(connection_matrix(field_neurons,:));
    active_neurons = find(activity);
    A_amp(active_neurons) = activity(active_neurons);
    imagesc(A_amp);
    drawnow
    phase_max = find(A_amp>Azim_amp);
    Azim_phase(phase_max) = j;
    Azim_amp = max(Azim_amp,A_amp);
end
