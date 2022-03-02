function llh = xyz2llh(xyz_1)
%XYZ2LLH	Convert from ECEF cartesian coordinates to 
%               latitude, longitude and height.  WGS-84
%
%	llh = XYZ2LLH(xyz)	
%
%    INPUTS
%	xyz(1) = ECEF x-coordinate in meters
%	xyz(2) = ECEF y-coordinate in meters
%	xyz(3) = ECEF z-coordinate in meters
%
%    OUTPUTS
%	llh(1) = latitude in radians
%	llh(2) = longitude in radians
%	llh(3) = height above ellipsoid in meters
xyz_1=reshape(xyz_1,[],3);
x = xyz_1(:,1);
y = xyz_1(:,2);
z = xyz_1(:,3);

x2 = xyz_1(:,1).^2;
y2 = xyz_1(:,2).^2;
z2 = xyz_1(:,3).^2;

a = 6378137.0000;	% earth radius in meters
b = 6356752.3142;	% earth semiminor in meters
e =  0.0818;%sqrt (1-(b/a).^2);b2
b2=b.^2;e2 = e.^2;ep = e*(a/b);
r = sqrt(x2+y2);r2 = r.*r;
E2 = a^2 - b^2;
F = 54.*b2.*z2;
G = r2 + (1-e2).*z2 - e2.*E2;
c = (e2.*e2.*F.*r2)./(G.^3);
s = ( 1 + c + sqrt(c.*c + 2.*c) ).^(1/3);
P = F ./ (3 * (s+1./s+1).^2 .* G.^2);
Q = sqrt(1+2.*e2.*e2.*P);
ro = -(P.*e2.*r)./(1+Q) + sqrt((a.*a./2).*(1+1./Q) ...
    - (P.*(1-e2).*z2)./(Q.*(1+Q)) - P.*r2./2);
tmp = (r - e2*ro).^2;
U = sqrt( tmp + z2 );
V = sqrt( tmp + (1-e2)*z2 );
zo = (b2.*z)./(a*V);

height = U.*( 1 - b2./(a*V) );

lat = atan( (z + ep*ep*zo)./r );

temp = atan(y./x);

% if x >=0
%     long = temp;
% elseif (x < 0) & (y >= 0)
%     long = pi + temp;
% else
%     long = temp - pi;
%
long = temp;
long((x < 0) & (y >= 0))=pi + long((x < 0) & (y >= 0));
long((x < 0) & (y < 0))= long((x < 0) & (y< 0))- pi;
llh(:,1) = lat;
llh(:,2) = long;
llh(:,3) = height;
% 
% if size(xyz_1,2)<3
%     
%     xyz_1=xyz_1';
% end
% %
%  for i=1:size(xyz_1,1)
%      xyz=xyz_1(i,1:3);
% 	x = xyz(1);
% 	y = xyz(2);
% 	z = xyz(3);
% 	x2 = x^2;
% 	y2 = y^2;
% 	z2 = z^2;
% 
% 	a = 6378137.0000;	% earth radius in meters
% 	b = 6356752.3142;	% earth semiminor in meters	
% 	e = sqrt (1-(b/a).^2);
% 	b2 = b*b;
% 	e2 = e^2;
% 	ep = e*(a/b);
% 	r = sqrt(x2+y2);
% 	r2 = r*r;
% 	E2 = a^2 - b^2;
% 	F = 54*b2*z2;
% 	G = r2 + (1-e2)*z2 - e2*E2;
% 	c = (e2*e2*F*r2)/(G*G*G);
% 	s = ( 1 + c + sqrt(c*c + 2*c) )^(1/3);
% 	P = F / (3 * (s+1/s+1)^2 * G*G);
% 	Q = sqrt(1+2*e2*e2*P);
% 	ro = -(P*e2*r)/(1+Q) + sqrt((a*a/2)*(1+1/Q) ...
%                                 - (P*(1-e2)*z2)/(Q*(1+Q)) - P*r2/2);
% 	tmp = (r - e2*ro)^2;
% 	U = sqrt( tmp + z2 );
% 	V = sqrt( tmp + (1-e2)*z2 );
% 	zo = (b2*z)/(a*V);
% 
% 	height = U*( 1 - b2/(a*V) );
% 	
% 	lat = atan( (z + ep*ep*zo)/r );
% 
% 	temp = atan(y/x);
% 	if x >=0	
% 		long = temp;
% 	elseif (x < 0) & (y >= 0)
% 		long = pi + temp;
% 	else
% 		long = temp - pi;
% 	end
% 
% 	llh(1) = lat;
% 	llh(2) = long;
% 	llh(3) = height;
%     llh1(i,1:3)=llh;
% %     i=i+1;
%  end
