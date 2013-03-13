% Polyphase decomposition in the frequency domain
% (matlab demo version only)
%
% (C) Alle Meije Wink (a.wink@vumc.nl)
%

M=32;               % signal length

x=0:M-1;            % time vector underlying the input signal
x=x(:);
x=x*(2*pi)/(M-1);

% types of functions
%x=sin(x);           % function f(t) describing the input signal
x=x.^2;              % function f(t) describing the input signal
%x=randn(M,1);       % function f(t) describing the input signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% go to the frequency domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=sqrt(-1);         % guarantee that i has the right value
X=fft(x);           % frequency representation of x
X=X(:);             % guarantee that X is a column vector

Q=4;                % downsampling factor (actually, the number of phases)
MQ=M/Q;             % length of each phase signal

k=(0:MQ-1)';        % time vector of each phase
phases=zeros(MQ,Q); % initialise output vector

%%%%%%%%%%%%%%%%%%%%
%%% downsampling %%%
%%%%%%%%%%%%%%%%%%%%

for j=0:Q-1         % cycle over phases
  Xjl=zeros(size(k));
  for l=0:Q-1       % cycle over frequency blocks
    expol=exp((2*pi*i*l*j)/Q); % complex exponential to account
			       % for the offset of each frequency block  
    Xjl=Xjl+expol*X(k+l*MQ+1); % encode the frequency block
  end
  expoj=exp((2*pi*i*k*j)/M);   % complex exponential to shift the
                               % right coefficient to index 0
  phases(:,j+1)=phases(:,j+1)+expoj.*Xjl; % shift the frequency block 
end
phases=phases/Q;      % account for adding together multiple freq. blocks

figure(1)
plot(real(ifft(phases)));

%%%%%%%%%%%%%%%%%%
%%% upsampling %%%
%%%%%%%%%%%%%%%%%%

k=(0:M-1)';           % time vector of reconstruction
recon=zeros(size(X)); % initialise reconstruction vector

for j=0:Q-1;
  expoj=exp(-(2*pi*i*k*j)/M);           % exponential to shift each
                                        % block back in the right position
  tmp=repmat(phases(:,j+1),Q,1).*expoj; % upsample and shift each block
  recon=recon+tmp;                      % add the usampled and shifted block
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% finally, go back to the time domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recon=real(ifft(recon));

figure(2)
plot([recon x]);

