load IMAGES.mat

num_trials=10000;
batch_size=100;

[imsize imsize num_images]=size(IMAGES);
BUFF=4;

% N = Pixels in Input Image
N = 1024; % 32x32 Images
sz = sqrt(N);

% L = Pixels in each Basis Function
L = 256;  % 16x16 Images
szb = sqrt(L);

% Pseudo-overcompleteness Level
OC = 1;

% M = Total number of Basis Functions
M = OC*(sz-(szb-1))^2;

% Init image 
I = zeros(N,batch_size);

% Random Initiatialization of Basis Functions
Phi = randn(L,M);
Phi = Phi * diag(1./sqrt(sum(Phi.*Phi)));

% Learning rate
eta = 3.0/batch_size;

% Lambda
lambda = 0.1;

display_every = 100;
display_network(Phi);

for t=1:num_trials
	imi=ceil(num_images*rand);

	for i=1:batch_size
        r=BUFF+ceil((imsize-sz-2*BUFF)*rand);
        c=BUFF+ceil((imsize-sz-2*BUFF)*rand);
        I(:,i)=reshape(IMAGES(r:r+sz-1,c:c+sz-1,imi),N,1);
    end

    % Calculate coefficients with LCA
    ahat = sparsify(I,Phi,lambda);

    % Calculate 
    ihat = reshape(I,sz,sz,batch_size);
	phihat = reshape(Phi,szb,szb,M);

	for b=1:batch_size
		for i=1:M/OC
			%current row
			j = ceil(i/(sz-(szb-1)));

			%current column
			k = mod(i,sz-(szb-1));
			if k == 0
				k = sz-(szb-1);
			end

			for z=1:OC
				ihat(j:j+szb-1,k:k+szb-1,b) = ihat(j:j+szb-1,k:k+szb-1,b) + ahat(i,b)*phihat(:,:,i);
			end
		end
	end

	ihat = reshape(ihat,sz*sz,batch_size);
	Phi = reshape(phihat,szb*szb,M);

	% Calculate Residuals
	R = I-ihat;
	R = reshape(R,sz,sz,batch_size);

	% Update Basis Functions
	dPhi = zeros(L,M);

	for i=1:M/OC
		for b=1:batch_size
			%current row
			j = ceil(i/(sz-(szb-1)));

			%current column
			k = mod(i,sz-(szb-1));
			if k == 0
				k = sz-(szb-1);
			end

			dPhi(:,i) = dPhi(:,i) + reshape(R(j:j+szb-1,k:k+szb-1,b),szb*szb,1) * ahat(i,b)';
		end
	end

	%{
	for b=1:batch_size
        dPhihat = dPhihat + R(:,b) * ahat(:,b)'; %  learning rule here
    end
    
    dPhihat = dPhihat/batch_size;
    dPhihat = reshape(dPhihat,16,16);

   	for i=1:M
   		%current row
		j = ceil(i/(sz-(szb-1)));

		%current column
		k = mod(i,sz-(szb-1));
		if k == 0
			k = sz-(szb-1);
		end
    	
    	dPhi(:,i) = reshape(dPhihat(j:j+szb-1,k:k+szb-1),szb*szb,1); 
    end 
	%}
    Phi = Phi + eta*dPhi;
    Phi=Phi*diag(1./sqrt(sum(Phi.*Phi))); % normalize bases

    if (mod(t,display_every)==0)
    	fprintf('Trial %d \n', t);
        display_network(Phi);
    end

end




