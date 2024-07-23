%% Figure 1 Code 

% Fig. 1c: spike-trajectory figure
% Fig. 1d: firing rate map
% Vertex score analysis/Firing field analysis (number, distance, size)

%% Session setting

idx_session =  4  ; % 1: Circle, 2: Triangle, 3: Square, 4: Hexagon

%% Data load

if idx_session == 1
    % Spike-trajectory map
    Fig1_beh_x_axis2 = Fig1_beh_x_axis_circ2 ; % Behavior x coordinates array
    Fig1_beh_y_axis2 = Fig1_beh_y_axis_circ2 ; % Behavior y coordinates array
    Fig1_spike_x_coordinates = Fig1_spike_x_coordinates_circ ; % Spike location - x coordinates
    Fig1_spike_y_coordinates = Fig1_spike_y_coordinates_circ ; % Spike location - y coordinates
    % Firing rate map
    Fig1_beh_x_axis = Fig1_beh_x_axis_circ ; % Behavior x coordinates array
    Fig1_beh_y_axis = Fig1_beh_y_axis_circ ; % Behavior y coordinates array
    Fig1_beh_t_axis = Fig1_beh_t_axis_circ ; % Behavior time array
    Fig1_beh_vDis_axis = Fig1_beh_vDis_axis_circ; % Behavior distance from vertex array
    Fig1_beh_spike_axis = Fig1_beh_spike_axis_circ  ; % Spike number array

elseif idx_session == 2
    % Spike-trajectory map
    Fig1_beh_x_axis2 = Fig1_beh_x_axis_tri2 ; % Behavior x coordinates array
    Fig1_beh_y_axis2 = Fig1_beh_y_axis_tri2 ; % Behavior y coordinates array
    Fig1_spike_x_coordinates = Fig1_spike_x_coordinates_tri ; % Spike location - x coordinates
    Fig1_spike_y_coordinates = Fig1_spike_y_coordinates_tri ; % Spike location - y coordinates
    % Firing rate map
    Fig1_beh_x_axis = Fig1_beh_x_axis_tri ; % Behavior x coordinates array
    Fig1_beh_y_axis = Fig1_beh_y_axis_tri ; % Behavior y coordinates array
    Fig1_beh_t_axis = Fig1_beh_t_axis_tri ; % Behavior time array
    Fig1_beh_vDis_axis = Fig1_beh_vDis_axis_tri; % Behavior distance from vertex array
    Fig1_beh_spike_axis = Fig1_beh_spike_axis_tri  ; % Spike number array
    % Vertex location
    x_vertex1 = -12 ;
    y_vertex1 = 21 ;
    x_vertex2 = 21;
    y_vertex2 = 3 ;
    x_vertex3 = -12 ;
    y_vertex3 = -16;

elseif idx_session == 3
    % Spike-trajectory map
    Fig1_beh_x_axis2 = Fig1_beh_x_axis_sqr2 ; % Behavior x coordinates array
    Fig1_beh_y_axis2 = Fig1_beh_y_axis_sqr2 ; % Behavior y coordinates array
    Fig1_spike_x_coordinates = Fig1_spike_x_coordinates_sqr ; % Spike location - x coordinates
    Fig1_spike_y_coordinates = Fig1_spike_y_coordinates_sqr ; % Spike location - y coordinates
    % Firing rate map
    Fig1_beh_x_axis = Fig1_beh_x_axis_sqr ; % Behavior x coordinates array
    Fig1_beh_y_axis = Fig1_beh_y_axis_sqr ; % Behavior y coordinates array
    Fig1_beh_t_axis = Fig1_beh_t_axis_sqr ; % Behavior time array
    Fig1_beh_vDis_axis = Fig1_beh_vDis_axis_sqr; % Behavior distance from vertex array
    Fig1_beh_spike_axis = Fig1_beh_spike_axis_sqr  ; % Spike number array
    % Vertex location
    x_vertex1 = -10 ;
    y_vertex1 = 19 ;
    x_vertex2 = 22;
    y_vertex2 = 19;
    x_vertex3 = 22;
    y_vertex3 = -13;
    x_vertex4 = -10;
    y_vertex4 = -13;
elseif idx_session == 4
    Fig1_beh_x_axis2 = Fig1_beh_x_axis_hex2 ; % Behavior x coordinates array
    Fig1_beh_y_axis2 = Fig1_beh_y_axis_hex2 ; % Behavior y coordinates array
    Fig1_spike_x_coordinates = Fig1_spike_x_coordinates_hex ; % Spike location - x coordinates
    Fig1_spike_y_coordinates = Fig1_spike_y_coordinates_hex ; % Spike location - y coordinates
    % Firing rate map
    Fig1_beh_x_axis = Fig1_beh_x_axis_hex ; % Behavior x coordinates array
    Fig1_beh_y_axis = Fig1_beh_y_axis_hex ; % Behavior y coordinates array
    Fig1_beh_t_axis = Fig1_beh_t_axis_hex ; % Behavior time array
    Fig1_beh_vDis_axis = Fig1_beh_vDis_axis_hex; % Behavior distance from vertex array
    Fig1_beh_spike_axis = Fig1_beh_spike_axis_hex  ; % Spike number array
    x_vertex1 = -1;
    y_vertex1 = 26;
    x_vertex2 = 20;
    y_vertex2 = 15;
    x_vertex3 = 20;
    y_vertex3 = -8;
    x_vertex4 = -1;
    y_vertex4 = -19;
    x_vertex5 = -19;
    y_vertex5 = -8;
    x_vertex6 = -19;
    y_vertex6 = 15;
end

%% Parameter setting

spacebin = 1.5  ; % spatial bin for firing rate map: 1.5 cm
Gaussian_smooting = 1 ; % Gaussian smoothing kernel (σ = 1) for firing rate analysis 
min_spk_num = 0 ; % minimum spike number for each spatial bin
min_t_span = 1 ; % minimum data point for each spatial bin
Firing_rate_threshold = 0.3 ; % threshold for firing rate map: 0.3 * max firing rate
min_field_size = 20/(spacebin*spacebin) ; % minimum firing field size: 20 cm^2

% x, y axis range for figure plotting
xlim_min = -30 ; 
xlim_max = xlim_min + 60 ; 
ylim_min = -30 ; 
ylim_max = ylim_min + 60 ;                                

%% Plot spike-trajectory map
figure() ;
% Plot behavior trajectory
L1 = plot(Fig1_beh_x_axis2 ,Fig1_beh_y_axis2 , 'linewidth', 1 ) ;  hold on;
L1.Color = [0.88 0.88 0.88 ] ; hold on;      
% Plot spikes
plot(Fig1_spike_x_coordinates, Fig1_spike_y_coordinates, 'o', 'markersize' , 3, 'markeredgecolor', 'r', 'markerfacecolor', 'r'); hold on ;
xlim([xlim_min , xlim_max ]);hold on;
ylim([ylim_min , ylim_max ]); hold on;
box off; axis off;


%% Firing rate map analysis

xsize = [-60 , 60 ];
ysize = [-60 , 60 ];
delta_t = Fig1_beh_t_axis(2) - Fig1_beh_t_axis(1) ;
xspace = [xsize(1):spacebin:xsize(end)]; % x coordiate for each spatial bin
yspace = [ysize(1):spacebin:ysize(end)]; % y coordiate for each spatial bin
zspace = zeros([(length(xspace) ), (length(yspace) )]); % mean firing rate for each spatial bin
zspace_spike = zeros([(length(xspace) ), (length(yspace) )]); % summed spike numbers for each spatial bin
zspace_time_span = zeros([(length(xspace) ), (length(yspace) )]); % time spent in each spatial bin
zspace_distance_vertex = zeros([(length(xspace) ), (length(yspace) )]);  % distance of each spatial bin from the closest vertex
                       
for x = 1:length(xspace)-1
    for y = 1:length(yspace)-1
        currentIdx = find(Fig1_beh_x_axis >= xspace(x) & Fig1_beh_x_axis< xspace(x+1) &...
            Fig1_beh_y_axis >= yspace(y) & Fig1_beh_y_axis < yspace(y+1));
        if ((length(currentIdx) > min_t_span)   )
            zspace_spike(x,y) = sum(Fig1_beh_spike_axis(currentIdx)) ;
            zspace_time_span(x,y) = length(currentIdx)*delta_t ;
            zspace_distance_vertex(x,y) =  mean(Fig1_beh_vDis_axis(currentIdx)) ;
            if zspace_spike(x,y) > min_spk_num
                zspace(x,y) = zspace_spike(x,y)/zspace_time_span(x,y)  ; % firing rate map
            end
        end
    end
end

%% Plot heatmap

figure() ;
xsize = [-60 , 60 ];
ysize = [-60 , 60 ];
xspace = [xsize(1):spacebin:xsize(end)];
yspace = [ysize(1):spacebin:ysize(end)];
[X, Y ] = meshgrid(xspace, yspace);

Filtered_zspace_figure = imgaussfilt(zspace' , 1.8  ,'FilterDomain','spatial', 'FilterSize', 5); % smoothed firing rate map
pcolor(X, Y, Filtered_zspace_figure) ;  hold on;
colormap(Fig1_jet) ;  hold on;
shading interp; hold on;
Max_c = max(max(Filtered_zspace_figure)) ;

xlim([xlim_min , xlim_max ]);hold on;
ylim([ylim_min , ylim_max ]); hold on;
box off; axis off;

%% Vertex score/Firing field analysis

if idx_session == 2  % Triangle 

Filtered_zspace = imgaussfilt(zspace' , Gaussian_smooting  ,'FilterDomain','spatial', 'FilterSize', 5); % Transposed
Max_firing_rate = max(max(Filtered_zspace)) ;

Temp_coordiation_matrix = [ ] ;
for y_idx = 1:size(Filtered_zspace, 2) % x coordinate
    for x_idx = 1:size(Filtered_zspace, 1) % y coordinate
        Temp_coordiation_matrix(end + 1, 1) = x_idx ; % y coordinate
        Temp_coordiation_matrix(end , 2) = y_idx ; %  x coordinate
    end
end

zspace_time_span = zspace_time_span' ; % Transposed
idx_time_span_ok = find( (zspace_time_span > min_t_span)   ) ;               

idx_firing_field_ok = find( (zspace_time_span > min_t_span) & (Filtered_zspace >(Firing_rate_threshold* Max_firing_rate) )  ) ;
Firing_field_logic = zeros(size(Filtered_zspace));
Firing_field_logic(idx_firing_field_ok) = 1 ;

Filtered_zspace_distance_vertex = zspace_distance_vertex' ;
Filtered_zspace_distance_vertex =  Filtered_zspace_distance_vertex.*  Filtered_zspace ;
Filtered_zspace_distance_vertex2 = Filtered_zspace_distance_vertex/mean( Filtered_zspace(idx_firing_field_ok) ) ;
idx_firing_field_ok2 = find( (zspace_time_span > min_t_span) & (Filtered_zspace >(Firing_rate_threshold* Max_firing_rate) )  ) ;
Firing_field_logic2 = zeros(size(Filtered_zspace));
Firing_field_logic2(idx_firing_field_ok2) = 1 ;


% Vertex1
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex1  )&(x_vertex1 < xspace(idxidx + 1)) )  )
        y_idx_vertex1 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex1  )&(y_vertex1 < yspace(idxidx + 1)) )  )
        x_idx_vertex1 =idxidx ; % Transposed matrix
    end
end
% vertex2
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex2  )&(x_vertex2 < xspace(idxidx + 1)) )  )
        y_idx_vertex2 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex2  )&(y_vertex2 < yspace(idxidx + 1)) )  )
        x_idx_vertex2 =idxidx ; % Transposed matrix
    end
end
% vertex3
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex3  )&(x_vertex3 < xspace(idxidx + 1)) )  )
        y_idx_vertex3 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex3  )&(y_vertex3 < yspace(idxidx + 1)) )  )
        x_idx_vertex3 =idxidx ; % Transposed matrix
    end
end

temp_num_field = 0 ;
temp_field_dis = [ ] ;
temp_field_size = [ ] ;

% Find firing field
% 1) Find the maximum firing spatial bin
% 2) Find all continuously adjacent spatial bins with firing rate above threshold (30% of maximum firing rate)
% 3) Remove already found spatial bins
% 4) Repeat this process until no spatial bin is found as firing fields

while ~isempty(find(Firing_field_logic2 > 0 ))
        
    idx_find = find(Firing_field_logic > 0 );
        To_search_idx = []  ;
    [ m , i ] = max(Filtered_zspace(idx_find))  ;
    To_search_idx(end + 1) = idx_find(i);
    Search_done_idx = [ ] ;
    
    while ~isempty(To_search_idx)
        
        Temp_kernel = ones(3,3) ;
        idx_search = 0 ;
        idx_search = To_search_idx(end);
        Temp_conv = zeros(size(Firing_field_logic)) ;
        Temp_conv(idx_search) = 1;
        Temp_neighbor_field = conv2(Temp_conv,Temp_kernel  , 'same' ) ;
        idx_neighbor = [];
        idx_neighbor = find(Temp_neighbor_field > 0 );
        idx_ttemp = [];
        idx_ttemp = Firing_field_logic(idx_neighbor) > 0  ;
        idx_found = [] ;
        idx_found = idx_neighbor(idx_ttemp) ;
        
        if ~isempty(idx_found)
            Search_done_idx(end + 1) = idx_search ;
        end
        
        To_search_idx(end) =  [] ;

        for idx_tempp = 1 : length(idx_found)
            if isempty(find(Search_done_idx == idx_found(idx_tempp) )  )
                if isempty(find(To_search_idx == idx_found(idx_tempp) )  )
                    if (pdist2([Temp_coordiation_matrix( idx_find(i), 1)  , Temp_coordiation_matrix( idx_find(i), 2)  ], [Temp_coordiation_matrix( idx_found(idx_tempp), 1)  , Temp_coordiation_matrix( idx_found(idx_tempp), 2)  ] )*spacebin <= 50 )
                        To_search_idx(end + 1) =  idx_found(idx_tempp) ;
                    end
                end
            end
        end

    end

    % Calculate distance of firing fields from the closest vertex
    if min_field_size <= length(Search_done_idx)
        temp_num_field = temp_num_field + 1 ;
        temp_field_size(end + 1) = length(Search_done_idx) * spacebin * spacebin ; % cm^2
        [ m , i ] =       max(Filtered_zspace(Search_done_idx))  ;                         
        d1 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex1, y_idx_vertex1] ) *spacebin) ;
        d2 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex2, y_idx_vertex2] ) *spacebin) ;
        d3 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex3, y_idx_vertex3] ) *spacebin) ;                 
        temp_field_dis(end + 1 , 1) = min( [ d1, d2 , d3  ] )     ;
    end

    Firing_field_logic(Search_done_idx) = 0;
    
    if isempty( find(Firing_field_logic > 0 ))
        break
    end

end

% Calculate the largest possible distance to vertex 
Largest_possible_dis = max(max(zspace_distance_vertex));

% Calculate vertex score (v) = 1 - d
% d is the distance from the (x, y) coordinate within the spatial receptive field with the peak firing rate to the nearest vertex, 
% which was normalized the largest possible distance to vertex in a given chamber. 

if ~isempty(temp_field_dis)  
    Temp_Vertex_Score = [ ] ;    
    for idx_temp_fieldN = [1: length(temp_field_dis(:,1))]
        % Temp_Vertex_Score(1,idx_temp_fieldN) = (1 - ( temp_field_dis(idx_temp_fieldN,1)/max(Temp_max_dis_array) ) ) ;
        Temp_Vertex_Score(1,idx_temp_fieldN) = (1 - ( temp_field_dis(idx_temp_fieldN,1)/Largest_possible_dis ) ) ;
    end

    disp('VScore: ')
    disp(mean( Temp_Vertex_Score( :)  , 'omitnan' ) )
    
    disp('The number of firing fields: ')
    disp( temp_num_field )

    disp('Distance from vertex (cm): ')
    mean( temp_field_dis( :)  , 'omitnan' ) 

    disp('Vertex field size (cm^2): ')
    mean( temp_field_size( :)  , 'omitnan' ) 

end

end

if idx_session == 3 % Square 

Filtered_zspace = imgaussfilt(zspace' , Gaussian_smooting  ,'FilterDomain','spatial', 'FilterSize', 5); % Transposed
Max_firing_rate = max(max(Filtered_zspace)) ;

Temp_coordiation_matrix = [ ] ;
for y_idx = 1:size(Filtered_zspace, 2) % x coordinate
    for x_idx = 1:size(Filtered_zspace, 1) % y coordinate
        Temp_coordiation_matrix(end + 1, 1) = x_idx ; % y coordinate
        Temp_coordiation_matrix(end , 2) = y_idx ; %  x coordinate
    end
end

zspace_time_span = zspace_time_span' ; % Transposed
idx_time_span_ok = find( (zspace_time_span > min_t_span)   ) ;               

idx_firing_field_ok = find( (zspace_time_span > min_t_span) & (Filtered_zspace >(Firing_rate_threshold* Max_firing_rate) )  ) ;
Firing_field_logic = zeros(size(Filtered_zspace));
Firing_field_logic(idx_firing_field_ok) = 1 ;

Filtered_zspace_distance_vertex = zspace_distance_vertex' ;
Filtered_zspace_distance_vertex =  Filtered_zspace_distance_vertex.*  Filtered_zspace ;
Filtered_zspace_distance_vertex2 = Filtered_zspace_distance_vertex/mean( Filtered_zspace(idx_firing_field_ok) ) ;
idx_firing_field_ok2 = find( (zspace_time_span > min_t_span) & (Filtered_zspace >(Firing_rate_threshold* Max_firing_rate) )  ) ;
Firing_field_logic2 = zeros(size(Filtered_zspace));
Firing_field_logic2(idx_firing_field_ok2) = 1 ;


% Vertex1
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex1  )&(x_vertex1 < xspace(idxidx + 1)) )  )
        y_idx_vertex1 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex1  )&(y_vertex1 < yspace(idxidx + 1)) )  )
        x_idx_vertex1 =idxidx ; % Transposed matrix
    end
end
% vertex2
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex2  )&(x_vertex2 < xspace(idxidx + 1)) )  )
        y_idx_vertex2 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex2  )&(y_vertex2 < yspace(idxidx + 1)) )  )
        x_idx_vertex2 =idxidx ; % Transposed matrix
    end
end
% vertex3
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex3  )&(x_vertex3 < xspace(idxidx + 1)) )  )
        y_idx_vertex3 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex3  )&(y_vertex3 < yspace(idxidx + 1)) )  )
        x_idx_vertex3 =idxidx ; % Transposed matrix
    end
end
% Vertex4
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex4  )&(x_vertex4 < xspace(idxidx + 1)) )  )
        y_idx_vertex4 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex4  )&(y_vertex4 < yspace(idxidx + 1)) )  )
        x_idx_vertex4 =idxidx ; % Transposed matrix
    end
end

temp_num_field = 0 ;
temp_field_dis = [ ] ;
temp_field_size = [ ] ;

% Find firing field
% 1) Find the maximum firing spatial bin
% 2) Find all continuously adjacent spatial bins with firing rate above threshold (30% of maximum firing rate)
% 3) Remove already found spatial bins
% 4) Repeat this process until no spatial bin is found as firing fields

while ~isempty(find(Firing_field_logic2 > 0 ))
        
    idx_find = find(Firing_field_logic > 0 );
    
    To_search_idx = []  ;
    [ m , i ] =       max(Filtered_zspace(idx_find))  ;
    To_search_idx(end + 1) = idx_find(i);
    
    Search_done_idx = [ ] ;
    
    while ~isempty(To_search_idx)
        
        Temp_kernel = ones(3,3) ;
        idx_search = 0 ;
        idx_search = To_search_idx(end);
        Temp_conv = zeros(size(Firing_field_logic)) ;
        Temp_conv(idx_search) = 1;
        Temp_neighbor_field = conv2(Temp_conv,Temp_kernel  , 'same' ) ;
        idx_neighbor = [];
        idx_neighbor = find(Temp_neighbor_field > 0 );
        idx_ttemp = [];
        idx_ttemp = Firing_field_logic(idx_neighbor) > 0  ;
        idx_found = [] ;
        idx_found = idx_neighbor(idx_ttemp) ;
        
        if ~isempty(idx_found)
            Search_done_idx(end + 1) = idx_search ;
        end
        
        To_search_idx(end) =  [] ;
        
        for idx_tempp = 1 : length(idx_found)
            if isempty(find(Search_done_idx == idx_found(idx_tempp) )  )
                if isempty(find(To_search_idx == idx_found(idx_tempp) )  )
                    if (pdist2([Temp_coordiation_matrix( idx_find(i), 1)  , Temp_coordiation_matrix( idx_find(i), 2)  ], [Temp_coordiation_matrix( idx_found(idx_tempp), 1)  , Temp_coordiation_matrix( idx_found(idx_tempp), 2)  ] )*spacebin <= 50 )
                        To_search_idx(end + 1) =  idx_found(idx_tempp) ;
                    end
                end
            end
        end

    end

    % Calculate distance of firing fields from the closest vertex
    if min_field_size <= length(Search_done_idx)
        temp_num_field = temp_num_field + 1 ;
        temp_field_size(end + 1) = length(Search_done_idx) * spacebin * spacebin ; % cm^2
        [ m , i ] =       max(Filtered_zspace(Search_done_idx))  ;                         
        d1 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex1, y_idx_vertex1] ) *spacebin) ;
        d2 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex2, y_idx_vertex2] ) *spacebin) ;
        d3 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex3, y_idx_vertex3] ) *spacebin) ;          
        d4 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex4, y_idx_vertex4] ) *spacebin) ;
        temp_field_dis(end + 1 , 1) =  min([ d1, d2 , d3  , d4 ] )   ;
    end

    Firing_field_logic(Search_done_idx) = 0;
    
    if isempty( find(Firing_field_logic > 0 ))
        break
    end

end

% Calculate the largest possible distance to vertex 
Largest_possible_dis = max(max(zspace_distance_vertex));

% Calculate vertex score (v) = 1 - d
% d is the distance from the (x, y) coordinate within the spatial receptive field with the peak firing rate to the nearest vertex, 
% which was normalized the largest possible distance to vertex in a given chamber. 
if ~isempty(temp_field_dis)  
    Temp_Vertex_Score = [ ] ;    
    for idx_temp_fieldN = [1: length(temp_field_dis(:,1))]
        Temp_Vertex_Score(1,idx_temp_fieldN) = (1 - ( temp_field_dis(idx_temp_fieldN,1)/Largest_possible_dis ) ) ;
    end

    disp('VScore: ')
    disp(mean( Temp_Vertex_Score( :)  , 'omitnan' ) )
    
    disp('The number of firing fields: ')
    disp( temp_num_field )

    disp('Distance from vertex (cm): ')
    mean( temp_field_dis( :)  , 'omitnan' ) 

    disp('Vertex field size (cm^2): ')
    mean( temp_field_size( :)  , 'omitnan' ) 

end

end

if idx_session == 4 % Hexagon 

Filtered_zspace = imgaussfilt(zspace' , Gaussian_smooting  ,'FilterDomain','spatial', 'FilterSize', 5); % Transposed
Max_firing_rate = max(max(Filtered_zspace)) ;

Temp_coordiation_matrix = [ ] ;
for y_idx = 1:size(Filtered_zspace, 2) % x coordinate
    for x_idx = 1:size(Filtered_zspace, 1) % y coordinate
        Temp_coordiation_matrix(end + 1, 1) = x_idx ; % y coordinate
        Temp_coordiation_matrix(end , 2) = y_idx ; %  x coordinate
    end
end

zspace_time_span = zspace_time_span' ; % Transposed
idx_time_span_ok = find( (zspace_time_span > min_t_span)   ) ;               

idx_firing_field_ok = find( (zspace_time_span > min_t_span) & (Filtered_zspace >(Firing_rate_threshold* Max_firing_rate) )  ) ;
Firing_field_logic = zeros(size(Filtered_zspace));
Firing_field_logic(idx_firing_field_ok) = 1 ;

Filtered_zspace_distance_vertex = zspace_distance_vertex' ;
Filtered_zspace_distance_vertex =  Filtered_zspace_distance_vertex.*  Filtered_zspace ;
Filtered_zspace_distance_vertex2 = Filtered_zspace_distance_vertex/mean( Filtered_zspace(idx_firing_field_ok) ) ;
idx_firing_field_ok2 = find( (zspace_time_span > min_t_span) & (Filtered_zspace >(Firing_rate_threshold* Max_firing_rate) )  ) ;
Firing_field_logic2 = zeros(size(Filtered_zspace));
Firing_field_logic2(idx_firing_field_ok2) = 1 ;

% Vertex1
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex1  )&(x_vertex1 < xspace(idxidx + 1)) )  )
        y_idx_vertex1 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex1  )&(y_vertex1 < yspace(idxidx + 1)) )  )
        x_idx_vertex1 =idxidx ; % Transposed matrix
    end
end
% vertex2
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex2  )&(x_vertex2 < xspace(idxidx + 1)) )  )
        y_idx_vertex2 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex2  )&(y_vertex2 < yspace(idxidx + 1)) )  )
        x_idx_vertex2 =idxidx ; % Transposed matrix
    end
end
% vertex3
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex3  )&(x_vertex3 < xspace(idxidx + 1)) )  )
        y_idx_vertex3 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex3  )&(y_vertex3 < yspace(idxidx + 1)) )  )
        x_idx_vertex3 =idxidx ; % Transposed matrix
    end
end
% Vertex4
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex4  )&(x_vertex4 < xspace(idxidx + 1)) )  )
        y_idx_vertex4 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex4  )&(y_vertex4 < yspace(idxidx + 1)) )  )
        x_idx_vertex4 =idxidx ; % Transposed matrix
    end
end
% Vertex5
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex5  )&(x_vertex5 < xspace(idxidx + 1)) )  )
        y_idx_vertex5 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex5  )&(y_vertex5 < yspace(idxidx + 1)) )  )
        x_idx_vertex5 =idxidx ; % Transposed matrix
    end
end
% Vertex6
for idxidx = 1 : ( length(xspace) - 1)
    if~isempty(find( (xspace(idxidx) <= x_vertex6  )&(x_vertex6 < xspace(idxidx + 1)) )  )
        y_idx_vertex6 =idxidx ; % Transposed matrix
    end
end
for idxidx = 1 : ( length(yspace) - 1)
    if~isempty(find( (yspace(idxidx) <= y_vertex6  )&(y_vertex6 < yspace(idxidx + 1)) )  )
        x_idx_vertex6 =idxidx ; % Transposed matrix
    end
end

temp_num_field = 0 ;
temp_field_dis = [ ] ;
temp_field_size = [ ] ;

% Find firing field
% 1) Find the maximum firing spatial bin
% 2) Find all continuously adjacent spatial bins with firing rate above threshold (30% of maximum firing rate)
% 3) Remove already found spatial bins
% 4) Repeat this process until no spatial bin is found as firing fields

while ~isempty(find(Firing_field_logic2 > 0 ))

    idx_find = find(Firing_field_logic > 0 );
    
    To_search_idx = []  ;
    [ m , i ] =       max(Filtered_zspace(idx_find))  ;
    To_search_idx(end + 1) = idx_find(i);
    
    Search_done_idx = [ ] ;
    
    while ~isempty(To_search_idx)
        
        Temp_kernel = ones(3,3) ;
        idx_search = 0 ;
        idx_search = To_search_idx(end);
        Temp_conv = zeros(size(Firing_field_logic)) ;
        Temp_conv(idx_search) = 1;
        Temp_neighbor_field = conv2(Temp_conv,Temp_kernel  , 'same' ) ;
        idx_neighbor = [];
        idx_neighbor = find(Temp_neighbor_field > 0 );
        idx_ttemp = [];
        idx_ttemp = Firing_field_logic(idx_neighbor) > 0  ;
        idx_found = [] ;
        idx_found = idx_neighbor(idx_ttemp) ;
        
        if ~isempty(idx_found)
            Search_done_idx(end + 1) = idx_search ;
        end
        
        To_search_idx(end) =  [] ;

        for idx_tempp = 1 : length(idx_found)
            if isempty(find(Search_done_idx == idx_found(idx_tempp) )  )
                if isempty(find(To_search_idx == idx_found(idx_tempp) )  )
                    if (pdist2([Temp_coordiation_matrix( idx_find(i), 1)  , Temp_coordiation_matrix( idx_find(i), 2)  ], [Temp_coordiation_matrix( idx_found(idx_tempp), 1)  , Temp_coordiation_matrix( idx_found(idx_tempp), 2)  ] )*spacebin <= 50 )
                        To_search_idx(end + 1) =  idx_found(idx_tempp) ;
                    end
                end
            end
        end

    end
    % Calculate distance of firing fields from the closest vertex
    if min_field_size <= length(Search_done_idx)
        temp_num_field = temp_num_field + 1 ;
        temp_field_size(end + 1) = length(Search_done_idx) * spacebin * spacebin ; % cm^2
        [ m , i ] =       max(Filtered_zspace(Search_done_idx))  ;                         
        d1 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex1, y_idx_vertex1] ) *spacebin) ;
        d2 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex2, y_idx_vertex2] ) *spacebin) ;
        d3 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex3, y_idx_vertex3] ) *spacebin) ;          
        d4 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex4, y_idx_vertex4] ) *spacebin) ;
        d5 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex5, y_idx_vertex5] ) *spacebin) ;
        d6 = (pdist2( [ Temp_coordiation_matrix( Search_done_idx(i)  , 1 )  ,  Temp_coordiation_matrix( Search_done_idx(i) , 2 ) ] ,  [x_idx_vertex6, y_idx_vertex6] ) *spacebin) ;
        temp_field_dis(end + 1 , 1) =  min( [ d1, d2 , d3  , d4  , d5  , d6 ] )    ;
    end

    Firing_field_logic(Search_done_idx) = 0;
    
    if isempty( find(Firing_field_logic > 0 ))
        break
    end

end

% Calculate the largest possible distance to vertex 
Largest_possible_dis = max(max(zspace_distance_vertex));

if ~isempty(temp_field_dis)  
    Temp_Vertex_Score = [ ] ;    
    for idx_temp_fieldN = [1: length(temp_field_dis(:,1))]
        Temp_Vertex_Score(1,idx_temp_fieldN) = (1 - ( temp_field_dis(idx_temp_fieldN,1)/Largest_possible_dis) ) ;
    end

    disp('VScore: ')
    disp(mean( Temp_Vertex_Score( :)  , 'omitnan' ) )
    
    disp('The number of firing fields: ')
    disp( temp_num_field )

    disp('Distance from vertex (cm): ')
    mean( temp_field_dis( :)  , 'omitnan' ) 

    disp('Vertex field size (cm^2): ')
    mean( temp_field_size( :)  , 'omitnan' ) 

end

end



