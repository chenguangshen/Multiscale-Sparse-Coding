% Total image pixels
R = 4096;

% Image height
H = sqrt(R);

% Max pixels per basis function
N = 64;
% Basis Function height
HB = sqrt(N);

% Number of basis functions
M = (H-(HB-1))^2; 

a = ones(N,M);
s = ones(M,1);
I = zeros(R,1);

I = reshape(I,H,H);
a = reshape(a,HB,HB,M);

for i=1:M
	%current row
	j = ceil(i/(H-(HB-1)));
	
	%current column
	k = mod(i,H-(HB-1));
	if k == 0
		k = H-(HB-1);
	end

	I(j:j+HB-1,k:k+HB-1) = I(j:j+HB-1,k:k+HB-1) + a(:,:,i);
end