function [ dist ] = latlon2dist(lat1,lon1,lat2,lon2)
%latlon2dist calculates distance between two points given by (lat,lon)

dlon = lon2 - lon1;
dlat = lat2 - lat1;
a = (sind(dlat/2))^2 + cosd(lat1) * cosd(lat2) * (sind(dlon/2))^2;
c = 2 * atan2( sqrt(a), sqrt(1-a) );
R=3961; % Radius of Earth in miles
dist = R * c; 

end

