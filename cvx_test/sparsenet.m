load IMAGES.mat

num_trials=10000;
batch_size=100;

[imsize imsize num_images]=size(IMAGES);
BUFF=4;

% N = Pixels in Input Image
N = 1024; % 32x32 Images
sz = sqrt(N);

% L = Starting # of pixels in each Basis Function
L = 64;  % 8x8 Images
szb = sqrt(L);
szb2 = szb+2;
L2 = szb2*szb2; 

% Border from which around to expand
B = 2;

% Pseudo-overcompleteness Level
OC = 1;

% M = Total number of Basis Functions
M = OC*(sz-(szb-1))^2;

% Init image 
I = zeros(N,batch_size);

% Random Initialization of Basis Functions
%Phi = randn(L,M);
% Zero Initialization of L2 Basis Functions
Phi2 = zeros(L2,M);
Phi = Phi * diag(1./sqrt(sum(Phi.*Phi)));

% Vector to record which level each basis function at
level = ones(M,1); % Start all basis functions at 1 

% Learning rate
eta = 3.0/batch_size;

% Lambda
lambda = 0.1;

display_every = 100;
display_network(Phi,Phi2);

% allow basis functions to converge before expanding again
expand_every = 400;

for t=1:num_trials
	imi=ceil(num_images*rand);

	for i=1:batch_size
        r=BUFF+ceil((imsize-sz-2*BUFF)*rand);
        c=BUFF+ceil((imsize-sz-2*BUFF)*rand);
        I(:,i)=reshape(IMAGES(r:r+sz-1,c:c+sz-1,imi),N,1);
    end

    % Calculate coefficients with LCA
    ahat = sparsify(I,Phi,Phi2,level,lambda);

    % Calculate 
    ihat = zeros(sz,sz,batch_size);
	phihat = reshape(Phi,szb,szb,M);
	phihat2 = reshape(Phi2,szb2,szb2,M);

	for b=1:batch_size
		for i=1:M/OC
			%current row
			j = ceil(i/(sz-(szb-1)));

			%current column
			k = mod(i,sz-(szb-1));
			if k == 0
				k = sz-(szb-1);
			end

			% Each OC "set" is contiguous
			for z=1:OC
				if level(i+(z-1)*M/OC) == 1
                	%Phi(j:j+szb-1,k:k+szb-1,i+(z-1)*M/OC) = Phihat(:,:,i+(z-1)*M/OC);
					ihat(j:j+szb-1,k:k+szb-1,b) = ihat(j:j+szb-1,k:k+szb-1,b) + ahat(i+(z-1)*M/OC,b)*phihat(:,:,i+(z-1)*M/OC);
            	else
                	p = j - max(0,(j+szb2)-sz); 
                	q = k - max(0,(k+szb2)-sz);
                	%Phi(p:p+szb2-1,q:q+szb2-1,i+(z-1)*M/OC) = Phihat2(:,:,i+(z-1)*M/OC);
					ihat(p:p+szb2-1,q:q+szb2-1,b) = ihat(p:p+szb2-1,q:q+szb2-1,b) + ahat(i+(z-1)*M/OC,b)*phihat2(:,:,i+(z-1)*M/OC);
            	end
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
	dPhi2 = zeros(L2,M);

	for i=1:M/OC
		for b=1:batch_size
			%current row
			j = ceil(i/(sz-(szb-1)));

			%current column
			k = mod(i,sz-(szb-1));
			if k == 0
				k = sz-(szb-1);
			end

			for z=1:OC
				if level(i+(z-1)*M/OC) == 1
                	%Phi(j:j+szb-1,k:k+szb-1,i+(z-1)*M/OC) = Phihat(:,:,i+(z-1)*M/OC);
					dPhi(:,i+(z-1)*M/OC) = dPhi(:,i+(z-1)*M/OC) + reshape(R(j:j+szb-1,k:k+szb-1,b),szb*szb,1) * ahat(i+(z-1)*M/OC,b)';
            	else
                	p = j - max(0,(j+szb2)-sz); 
                	q = k - max(0,(k+szb2)-sz);
                	%Phi(p:p+szb2-1,q:q+szb2-1,i+(z-1)*M/OC) = Phihat2(:,:,i+(z-1)*M/OC);
					dPhi2(:,i+(z-1)*M/OC) = dPhi2(:,i+(z-1)*M/OC) + reshape(R(p:p+szb2-1,q:q+szb2-1,b),szb2*szb2,1) * ahat(i+(z-1)*M/OC,b)';
            	end
			end
		end
	end

    Phi = Phi + eta*dPhi;
    Phi = Phi*diag(1./sqrt(sum(Phi.*Phi))); % normalize bases

    % need to be careful about 0s
    s = sum(Phi2.*Phi2);
    s(~s(:)) = 1;

    Phi2 = Phi2 + eta*dPhi2;
    Phi2 = Phi2*diag(1./sqrt(s)); % normalize bases


    % Expand basis functions if needed
   	if(mod(t,expand_every)==0)
   		phihat = reshape(Phi,szb,szb,M);
   		for i=1:M
   			% Fit amplitude envelope to basis function and find center
   			hf = hilbert(phihat(:,:,i));
   			[hfmax, hfmaxi] = max(abs(hf(:)));
			[yc, xc] = ind2sub(size(hf), hfmaxi);
   		
			if xc <= B || xc >= (szb - B) || yc <= B || yc >= (szb - B)
				level(i) = 2; % move to level 2
				Phi2(:,i) = reshape(padarray(phihat(:,:,i),[1 1]),szb2*szb2,1);
			end
   		end
   	end 

   	if (mod(t,display_every)==0)
    	fprintf('Trial %d \n', t);
        display_network(Phi,Phi2);
    end

    % Check if basis needs expanding
end




