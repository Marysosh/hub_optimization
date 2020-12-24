%получение данных из файлов csv 
china_airports = getting_data(chinaAirports);
eu_airports = getting_data(euAirpots);
airports_data = cat(1,china_airports, eu_airports);

%нормализация пассажиропотоков
normalized_passengers = airports_data(:,3);

longitude = airports_data(:,1);
latitude = airports_data(:,2);

%стартовые параметры
C=1;

%maxIter = 100;
hub_lat1 = round(rand + 43,7);
hub_start_lat = hub_lat1;

hub_long1 = round(rand + 60,7);
hub_start_long = hub_long1;

step = 0.05;

syms hub_long
syms hub_lat

exp1_cost_final = exp_cost_final(hub_long, hub_lat, normalized_passengers, C, airports_data(:,1:2));

%вычисление частных производных
dexp_lat = diff(exp1_cost_final, hub_long);
dexp_long = diff(exp1_cost_final, hub_lat);

% for i = 1:maxIter
hub_lat11 = hub_lat1;
hub_long11 = hub_long1;
i=0;
%применение метода градиентного спуска
while (sqrt(double(hub_lat11)^2 + double(hub_long11)^2)) > 0.005
    
    disp('Итерация:')
    disp(i)
    
    hub_lat_prev = hub_lat1;
    hub_long_prev = hub_long1;
    hub_lat11 = subs(dexp_lat,'hub_long',hub_long_prev);
    hub_lat11 = subs(hub_lat11,'hub_lat', hub_lat_prev);
    hub_lat1 = round(hub_lat_prev - step * hub_lat11,7);
    
    hub_long11 = subs(dexp_long,'hub_long', hub_long_prev);
    hub_long11 = subs(hub_long11, 'hub_lat', hub_lat_prev);
    hub_long1 = round(hub_long_prev - step * hub_long11,7);
    result = round(exp_cost(hub_lat1, hub_long1, normalized_passengers, C, airports_data(:,1:2)),7)
    
    i = i+1;
    
    %grad = double(sqrt(double(hub_lat11)^2 + double(hub_long11)^2))
end

%подсчет значения целевой функции
result = round(exp_cost(hub_lat1, hub_long1, normalized_passengers, C, airports_data(:,1:2)),7)

%построение координат точек на карте мира
load topo
topoR = georefcells(topolatlim,topolonlim,size(topo))
worldmap('world')
geoshow(topo,topoR,'DisplayType','texturemap')
demcmap(topo)

latitude_ch = latitude(1:20);
latitude_eu = latitude(21:40);
longitude_ch = longitude(1:20);
longitude_eu = longitude(21:40);

scatterm(latitude_ch,longitude_ch,'filled','MarkerFaceColor',[1 0 1])
scatterm(latitude_eu,longitude_eu,'filled','MarkerFaceColor',[1 1 0])
scatterm(double(hub_lat1),double(hub_long1),'filled','MarkerFaceColor',[1 0 0])

function [airports_data] = getting_data(airports)
    coordinates = table2array([airports(:,6) airports(:,7)]);
    passengers = table2array(airports(:,5));
    normalized_passengers = passengers(:,1) ./ sum(passengers(:,1));
    airports_data = cat(2,coordinates,normalized_passengers);
end

function [distances] = distances_to_point(point_long, point_lat, all_coordinates)
    EARTH_RADIUS = 6371.009;
    point_coordinates = [point_long point_lat];
    point_rad = deg2rad(point_coordinates(1,:));
    all_coordinates_rad = deg2rad(all_coordinates);
    
    sin_lat1  = sin(point_rad(:,1)); 
    cos_lat1 = cos(point_rad(:,1));
    sin_lat2  = sin(all_coordinates_rad(:,1)) ;  
    cos_lat2 = cos(all_coordinates_rad(:,1));
    
    delta_lng = all_coordinates_rad(:,2) - point_rad(:,2);
    cos_delta_lng = cos(delta_lng);
    sin_delta_lng = sin(delta_lng);
    
    d = atan2(sqrt((cos_lat2.*sin_delta_lng).^2 + (cos_lat1.*sin_lat2 - sin_lat1.*cos_lat2.*cos_delta_lng).^2),sin_lat1.*sin_lat2 + cos_lat1.*cos_lat2.*cos_delta_lng);
    distances = EARTH_RADIUS * d;
    
end

function [expected_cost] = exp_cost(hub_long, hub_lat, normalized_passengers, C, all_coordinates)
    distances = distances_to_point(hub_long, hub_lat, all_coordinates);
    expected_cost = sum(exp(C.*distances ./ 1000) .* normalized_passengers);
end

function [expected_cost] = exp_cost_final(hub_long, hub_lat, normalized_passengers, C, all_coordinates)
    EARTH_RADIUS = 6371.009;
    all_coordinates_rad = deg2rad(all_coordinates);
    expected_cost = sum(exp(C.*(EARTH_RADIUS .* atan2(sqrt((cos(all_coordinates_rad(:,1)).*sin(all_coordinates_rad(:,2) - deg2rad(hub_lat))).^2 + (cos(deg2rad(hub_long)).*sin(all_coordinates_rad(:,1)) - sin(deg2rad(hub_long)).*cos(all_coordinates_rad(:,1)).*cos(all_coordinates_rad(:,2) - deg2rad(hub_lat))).^2),sin(deg2rad(hub_long)).*sin(all_coordinates_rad(:,1)) + cos(deg2rad(hub_long)).*cos(all_coordinates_rad(:,1)).*cos(all_coordinates_rad(:,2) - deg2rad(hub_lat)))) ./ 1000) .* normalized_passengers);
end