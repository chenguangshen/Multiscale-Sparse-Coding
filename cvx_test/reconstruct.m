% Total image pixels
R = 1024;

% Image height
H = sqrt(R);

% Max pixels per basis function
N = 4;
% Basis Function height
HB = sqrt(N);

% Number of basis functions
M = (sqrt(R)-1)^2; 

a = ones(N,M);
s = ones(M,1);
I = zeros(R,1);

for i=1:M
	%current column
	j = ((i-1)*H)+1;

	if (j+H+1 > R)
		q = j;
		j = mod(j,R-(H+1));
	end

	% last square in bottom right corner
	if j == 0
		j = ((H-2)*H)+(H-1);
	end

	I(j:j+1) = I(j:j+1) + a(1:2,i);
	I(j+H:j+H+1) = I(j+H:j+H+1) + a(3:4,i);
	%I(j+H+2:j+H+1) = I(j+H:j+H+1) + a(9:16,i);
end