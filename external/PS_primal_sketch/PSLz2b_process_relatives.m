function relative = PSLz4_process_relatives(relatives,current);
% relative = PSLz4_process_relatives(relatives,current)
%
% Process potential successors or predecessor of each point and choose the
% single one that is most likely to be the actual neighbor.
% Code is based on the Nevatia & Babu algorithm, but does not deal with 
% triple points.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

structure = relatives; expand_structure;

%%-------------------------------------------------------------------------
%% part 1:
%% find which orientations are fine with current linelet's  orientation:
%% generates in batch mode the neighboring orientations
%% into which the currnet linelets may branch. It is  assumed
%% orientations can change by 45 degrees max at each pixel
%%-------------------------------------------------------------------------

allowed = mod(repmat(current.orient',[3,1]) + repmat([-1:1]',[1,length(current.orient)])-1,4)+1;

%%-------------------------------------------------------------------------
%% and  find for each of the possible neigbors
%% check whether his orientation is within the set of allowable 
%% orientations 
%%-------------------------------------------------------------------------
ok = zeros(size(neigh_orient));

for k_potential = [1:3],        
    for k_neigh  =[1:3],        
        ok(k_neigh,:)= max(ok(k_neigh,:),allowed(k_potential,:) == neigh_orient(k_neigh,:));
        ok(k_neigh,:)= 1;
    end
end

%%-------------------------------------------------------------------------
%% part 2:
%% for each linelet check how many of its potential neighbors 
%% actually lie on the feature skeleton. Treat separately the
%% cases with double and single neighbors (we do not treat 
%% triple neighbors).
%%-------------------------------------------------------------------------
relative= -1*ones(1,length(neigh_ind));
nneighs = sum(neigh_active,1);
single_neigh = find(nneighs==1);
double_neigh = find(nneighs==2);

%%-------------------------------------------------------------------------
%% 2.1 single neighbor 
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
%% get `OK' index of the single active neighbor (batch mode code)
[closest_neighbor,dum] = find(neigh_active(:,single_neigh));
idx  = closest_neighbor' + (single_neigh-1)*3;
condition_single = ok(idx); 
continued = single_neigh(find(condition_single));

%%-------------------------------------------------------------------------
%% for those that are indeed ok, store the index of their 
%% single closest neighbor
relative(continued) = neigh_ind(idx(find(condition_single)));

%%-------------------------------------------------------------------------
%% 2.2 double neighbors 
%% case 1: 4-connected, similar orientation -> no fork
%% case 2: 4-connected, vertical orientations-> fork
%% case 3: 8-connected  -> fork 
%%-------------------------------------------------------------------------

% for ease, gather all data related to the 2-neighbor pairs
[neigh_active_d,neigh_ener_d,neigh_orient_d,dists_d]=  keep_points_double...
    (double_neigh,neigh_active,neigh_ener,neigh_orient,neigh_dist);

for k=1:length(double_neigh),
    %% examine separately every linelet's neighbors:    
    index = double_neigh(k);
    [neighs] = find(neigh_active_d(:,k));
    
    %% get the orientations, energies and distances
    %% of the two active neighbors of the current linelet
    orient(1) = neigh_orient_d(neighs(1),k);
    orient(2) = neigh_orient_d(neighs(2),k);

    ener_n(1) = neigh_ener_d(neighs(1),k);
    ener_n(2) = neigh_ener_d(neighs(2),k);

    dist_n(1) = dists_d(neighs(1),k);
    dist_n(2) = dists_d(neighs(2),k);

    %% are they on contiguous locations? 
    f_connected = (abs(neighs(2)-neighs(1))==1);
    
    %% do they have similar orientations
    siml_orient = ismember2(orient(1),mod(orient(2)+[-1:1]-1,4)+1);

    %% case 1: 4-connected, similar orientation -> no fork. e.g.:
    %% 0 l
    %% l l
    %% 0 0
    %% keep the neighbor that is closest (line will go through
    %% other neighbor afterwards, without changing route).
    if (f_connected)&(siml_orient)
        [mn,idx_min] = min(dist_n);
        rel = neighs(idx_min); 
    end

    %% case 2: 4-connected, vertical orientations-> fork
    %% 0 1
    %% 1 1
    %% 0 1
    %% not implemented yet

    
    %% case 3: 8-connected  -> fork
    %% 0 1
    %% 1 0
    %% 0 1
    %% keep the most active neighbor (naive heuristic)
    
    if (~siml_orient)|(~f_connected)
        [mx,idx_max] = max(ener_n);
        rel = neighs(idx_max);
    end    
    relative(index) = neigh_ind(rel,index);
end

%%-------------------------------------------------------------
%% internal functions
%%-------------------------------------------------------------
function res = ismember2(input,options),
res = (input == options(1))||(input==options(2));

function varargout = keep_points_double(pts,varargin);
for k=1:length(varargin)
    varargout{k} = varargin{k}(:,pts);
end