N = 4;
M = 9; 

R = 16; 

a = ones(N,M);
s = ones(M,1);
I = zeros(R,1);

b = reshape(a, 2, 2, M);

for i=1:M
	%current column
	j = ((i-1)*N)+1;

	if (j+5 > R)
		q = j;
		j = mod(j,11);
	end

	if j == 0
		j = ((3-1)*N)+3;;
	end

	I(j:j+1) = I(j:j+1) + b(:,1,i);
	I(j+4:j+5) = I(j+4:j+5) + b(:,2,i);
end