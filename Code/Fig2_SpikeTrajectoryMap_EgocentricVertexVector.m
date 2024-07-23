%% Figure 2 Code 

% Fig. 2: egocentric vector analysis
% Fig. 2d: spike-trajectory figure with head direction
% Fig. 2f: egocentric vertex-vector map

%% Session setting

idx_session =  3  ; % 1: Triangle, 2: Square, 3: Hexagon

%% Data load

if idx_session == 1
    % Spike-trajectory map
    Fig2_beh_x_axis2 = Fig2_beh_x_axis_tri2 ; % Behavior x coordinates array
    Fig2_beh_y_axis2 = Fig2_beh_y_axis_tri2 ; % Behavior y coordinates array
    Fig2_spike_x_coordinates = Fig2_spike_x_coordinates_tri ; % Spike location - x coordinates
    Fig2_spike_y_coordinates = Fig2_spike_y_coordinates_tri ; % Spike location - y coordinates
    Fig2_spike_vDeg_coordinates = Fig2_spike_vDeg_coordinates_tri ; % Spike location - corresponding egocentric degree from the closest vertex
    Fig2_spike_HD_coordinates = Fig2_spike_HD_coordinates_tri ; % Spike location - corresponding head direction
    % Egocentric vertex-vector analysis
    Fig2_beh_x_axis = Fig2_beh_x_axis_tri ; % Behavior x coordinates array
    Fig2_beh_y_axis = Fig2_beh_y_axis_tri ; % Behavior y coordinates array
    Fig2_beh_t_axis = Fig2_beh_t_axis_tri ; % Behavior time array
    Fig2_beh_vDis_axis = Fig2_beh_vDis_axis_tri; % Behavior distance from vertex array
    Fig2_beh_vDeg_axis = Fig2_beh_vDeg_axis_tri  ; % Egocentric degree from vertex array
    Fig2_beh_spike_axis = Fig2_beh_spike_axis_tri  ; % Spike number array
elseif idx_session == 2
    % Spike-trajectory map
    Fig2_beh_x_axis2 = Fig2_beh_x_axis_sqr2 ; % Behavior x coordinates array
    Fig2_beh_y_axis2 = Fig2_beh_y_axis_sqr2 ; % Behavior y coordinates array
    Fig2_spike_x_coordinates = Fig2_spike_x_coordinates_sqr ; % Spike location - x coordinates
    Fig2_spike_y_coordinates = Fig2_spike_y_coordinates_sqr ; % Spike location - y coordinates
    Fig2_spike_vDeg_coordinates = Fig2_spike_vDeg_coordinates_sqr ; % Spike location - corresponding egocentric degree from the closest vertex
    Fig2_spike_HD_coordinates = Fig2_spike_HD_coordinates_sqr ; % Spike location - corresponding head direction
    % Egocensqrc vertex-vector analysis
    Fig2_beh_x_axis = Fig2_beh_x_axis_sqr ; % Behavior x coordinates array
    Fig2_beh_y_axis = Fig2_beh_y_axis_sqr ; % Behavior y coordinates array
    Fig2_beh_t_axis = Fig2_beh_t_axis_sqr ; % Behavior time array
    Fig2_beh_vDis_axis = Fig2_beh_vDis_axis_sqr; % Behavior distance from vertex array
    Fig2_beh_vDeg_axis = Fig2_beh_vDeg_axis_sqr  ; % Egocensqrc degree from vertex array
    Fig2_beh_spike_axis = Fig2_beh_spike_axis_sqr  ; % Spike number array
elseif idx_session == 3
    % Spike-trajectory map
    Fig2_beh_x_axis2 = Fig2_beh_x_axis_hex2 ; % Behavior x coordinates array
    Fig2_beh_y_axis2 = Fig2_beh_y_axis_hex2 ; % Behavior y coordinates array
    Fig2_spike_x_coordinates = Fig2_spike_x_coordinates_hex ; % Spike location - x coordinates
    Fig2_spike_y_coordinates = Fig2_spike_y_coordinates_hex ; % Spike location - y coordinates
    Fig2_spike_vDeg_coordinates = Fig2_spike_vDeg_coordinates_hex ; % Spike location - corresponding egocentric degree from the closest vertex
    Fig2_spike_HD_coordinates = Fig2_spike_HD_coordinates_hex ; % Spike location - corresponding head direction
    % Egocenhexc vertex-vector analysis
    Fig2_beh_x_axis = Fig2_beh_x_axis_hex ; % Behavior x coordinates array
    Fig2_beh_y_axis = Fig2_beh_y_axis_hex ; % Behavior y coordinates array
    Fig2_beh_t_axis = Fig2_beh_t_axis_hex ; % Behavior time array
    Fig2_beh_vDis_axis = Fig2_beh_vDis_axis_hex; % Behavior distance from vertex array
    Fig2_beh_vDeg_axis = Fig2_beh_vDeg_axis_hex  ; % Egocenhexc degree from vertex array
    Fig2_beh_spike_axis = Fig2_beh_spike_axis_hex  ; % Spike number array
end

%% Parameter setting

degreebin = 15 ; %  Degree bin for egocentric vertex-vector map: 15^o
distancebin = 1.5  ; %  Distance bin for egocentric vertex-vector map: 1.5 cm
Gaussian_smooting = 0.7 ; % Gaussian smoothing kernel (σ = 1) for egocentric vertex-vector map
min_spk_num = 0 ; % minimum spike number for each spatial bin
min_t_span = 1 ; % minimum data point for each spatial bin
max_dis = 25 ; % maximum distance for egocentric vertex-vector map: 1.5 cm
temp_color_shift  = -10 ; % colormap setting
idx_pref_deg_range = 2;


% x, y axis range for figure plotting
xlim_min = -25 ; 
xlim_max = xlim_min + 60 ; 
ylim_min = -25 ; 
ylim_max = ylim_min + 60 ;     

%% Egocentric vertex-vector analysis

% Egocentric vertex-vector map

Dis_size = [0, (max_dis +distancebin) ];
Deg_size = [-180, (180 + degreebin) ];
delta_time = Fig2_beh_t_axis(2) - Fig2_beh_t_axis(1) ;
Dis_space   = [Dis_size(1):distancebin:Dis_size(end)];
Deg_space = [Deg_size(1):degreebin:Deg_size(end)];
zspace = zeros([(length(Dis_space) - 1 ) , (length(Deg_space)  - 1  ) ]); % mean firing rate for each distance x degree bin
zspace_spike = zeros([(length(Dis_space)  - 1) , (length(Deg_space)   - 1 ) ]); % summed spike numbers for each distance x degree bin
zspace_time_span = zeros([(length(Dis_space)  - 1 ) , (length(Deg_space)  - 1  ) ]); % time spent in each distance x degree bin

Indi_cell_vector = 0 ;

for x = 1:length(Dis_space)-1
    for y = 1:length(Deg_space)-1
        currentIdx = find(Fig2_beh_vDis_axis >= Dis_space(x) & Fig2_beh_vDis_axis< Dis_space(x+1) &...
            Fig2_beh_vDeg_axis >= Deg_space(y) & Fig2_beh_vDeg_axis < Deg_space(y+1));
        if (length(currentIdx) > min_t_span)
            zspace_spike(x,y) = sum(Fig2_beh_spike_axis(currentIdx)) ;
            zspace_time_span(x,y) = length(currentIdx)*delta_time ;
            if zspace_spike(x,y) > min_spk_num
                zspace(x,y) = zspace_spike(x,y)/zspace_time_span(x,y)  ;
            end
            Indi_cell_vector   = Indi_cell_vector + ( zspace(x,y) * exp(1j*deg2rad(Deg_space(y)))   )  ; % calculate egocentric vertex-vector
        end
    end
end

Vlength = 0 ; % Mean resultant vector length
Vdeg = 0 ; % Preferred degree (^o)
Vlength = abs(  Indi_cell_vector) / ( ( ( length(Dis_space) -1 ) * ( length(Deg_space) -1 )) * mean(mean( zspace(:,:)))  ) ;  % calculate mean resultant vector length
Vdeg = atan2d( imag( Indi_cell_vector),  real( Indi_cell_vector)) ;% calculate preferred degree
Vdis = 0 ; % Preferred distance (cm) 

Filtered_zspace = imgaussfilt(zspace' , Gaussian_smooting ) ;
Filtered_zspace = Filtered_zspace' ;

% Weibull distribution fitting to calculate preferred distance of egocentric vertex-vector
for y = 1:length(Deg_space)-1
    if ( Vdeg >= Deg_space(y)   ) && ( Vdeg< Deg_space(y + 1)   )
        a = Dis_space(1: end -1)' + 0.5  ;
        temp_range = circshift([1:(length(Deg_space)-1)] , -(y - 1) +idx_pref_deg_range ) ;
        b = mean(Filtered_zspace(:, temp_range(1:(2*idx_pref_deg_range + 1)) ) , 2 ) + 1 ;
        fit_weibull = fit(a, b, 'c*a*b*x^(b-1)*exp(-a*x^b)' , 'StartPoint', [0.01, 3 , 3 ] , 'Lower', [0 -Inf 0] , 'Exclude', b>= 50 )  ;
        for idx= 1 :  5
            fit_weibull = fit(a, b, 'c*a*b*x^(b-1)*exp(-a*x^b)', 'StartPoint', [fit_weibull.a, fit_weibull.b , fit_weibull.c ], 'Lower', [0 -Inf 0] , 'Exclude', b >= 50  )  ;
        end
        f = fittype( 'c*a*b*x^(b-1)*exp(-a*x^b)') ;
        c = cfit(f, fit_weibull.a , fit_weibull.b, fit_weibull.c ) ;
        a2 = [0:0.1:max_dis]'  ;
        y1 = feval(f, fit_weibull.a , fit_weibull.b, fit_weibull.c , a2)  ;
        [m_temp, i_temp] =  max(  y1 ) ;
        Vdis =  a2(i_temp) 
    end
end

disp('Mean resultant vector length: ')
disp( Vlength )

disp('Preferred degree (^o): ')
disp( Vdeg )

disp('Preferred distance (cm): ')
disp( Vdis )

%% Plot spike-trajectory map with color-coding of head dirction
    
figure()
% Plot behavior trajectory
L1 = plot(Fig2_beh_x_axis2 ,Fig2_beh_y_axis2 , 'linewidth', 1 ) ;  hold on;
L1.Color = [0.88 0.88 0.88 ] ; hold on;      

% Assigning spike plotting order depending on the corresponding egocentric degree from the closest vertex
Preferred_deg = Vdeg ;
temp_spike_vector_deg_order = zeros( length(Fig2_spike_x_coordinates) , 5 ) ;
temp_spike_vector_deg_order(:, 1) = [1 : length(Fig2_spike_x_coordinates)  ] ;
temp_spike_vector_deg_order(:, 2) = Fig2_spike_vDeg_coordinates' ;
temp_spike_vector_deg_order(:, 3) = Preferred_deg ;
temp_spike_vector_deg_order(:, 4) = rad2deg(angdiff(deg2rad(Fig2_spike_vDeg_coordinates') , deg2rad(Preferred_deg*ones(length(Fig2_spike_vDeg_coordinates) , 1)))) ;
temp_spike_vector_deg_order(:, 5) = abs(angdiff(deg2rad(Fig2_spike_vDeg_coordinates') , deg2rad(Preferred_deg*ones(length(Fig2_spike_vDeg_coordinates) , 1))))   ;
temp_spike_vector_deg_order = sortrows(temp_spike_vector_deg_order , 5, 'descend') ;


color_matrix = [hsv(722)] ;
color_matrix2 = circshift(color_matrix, 361 + temp_color_shift ) ;
idx_order_deg = [-180 : 0.5 : (180  + 0.5 ) ] ;
Temp_deg_idx_corrdinates = zeros(1, length(Fig2_spike_x_coordinates)) ;
for idx_deg = 1 : ( length(idx_order_deg) - 1)
    idx_temptemp = find(      ( idx_order_deg(idx_deg) <=  Fig2_spike_HD_coordinates  )  &   (   Fig2_spike_HD_coordinates  < idx_order_deg(idx_deg + 1) )         )        ;
    Temp_deg_idx_corrdinates(idx_temptemp) = idx_deg ;
end
% plot spike with corresponding head direction by color-coding
for idx_event = temp_spike_vector_deg_order(:,1)'
    color_temp = [ color_matrix2( Temp_deg_idx_corrdinates(idx_event) , 1 )    color_matrix2( Temp_deg_idx_corrdinates(idx_event), 2 )      color_matrix2( Temp_deg_idx_corrdinates(idx_event), 3 )     ] ;
    plot(Fig2_spike_x_coordinates(idx_event), Fig2_spike_y_coordinates(idx_event), 'o', 'markersize' , 3, 'markeredgecolor', color_temp, 'markerfacecolor', color_temp); hold on ;
end

xlim([xlim_min , xlim_max ]);hold on;
ylim([ylim_min , ylim_max ]); hold on;
box off; axis off;

%% Plot egocentric vertex-vector map

Dis_space2 = Dis_space ; 
Deg_space2 = Deg_space ; 
zspace2 = zspace ; 
Dis_space(end) = [] ;
Deg_space(end ) = [];
zspace(:, end ) = zspace(:, 1) ;
Filtered_zspace = imgaussfilt(zspace' , Gaussian_smooting );

figure();
[h, c] = polarPcolor( Dis_space,  Deg_space, Filtered_zspace, 'Ncircles', 1 , 'Nspokes', 1) ; hold on;
colormap(Fig2_jet) ;  hold on;
shading interp; hold on;
Max_c2 = max(max(Filtered_zspace)) 
caxis([0      Max_c2]);
ax = gca;
colorbar('off') ;
box off; axis off;
    










