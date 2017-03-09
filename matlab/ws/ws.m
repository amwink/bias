function smallworldcoeffs=ws(nv,preserve_degree,runningshow)
%
% simulates the generation of small-world networks using 
% the recipe by Watts and Strogatz (1998)
%
% syntax: swc = wattsstrogatz(nv,preserve_degree,runningshow)
%         
% inputs:  nv  -- 
%           => number of vertices of the graph (default:128)
%           => the number of edges is nv*2
%          preserve_degree --
%           => preserve node degrees (only let nodes swap connections)
%           => default: 1
%          runningshow --
%           => update the 4 plots after every re-labelling
%           => default: 1
%           => plot1: small worldness (r) path length (g) and clustering (b)
%           => plot2: connections, starting with ring-lattice
%           => plot3: connections matrix
%           => plot4: distance matrix
%
% output:  swc -- the small world coefficients after each re-labelling
%
% The final plot1 computes the normalised clustering, pathlength and small
% worldness based on the values of the final re-labelled graph. The optimal
% small worldness corresponds with the peak of the red line. Note that this 
% is usually after a relatively small number of re-labellings!
%
% written by Alle Meije Wink -- 3/3/2012 -- a.wink@vumc.nl
% if you use this code please give proper credit -- thanks
%

% use default arguments if none provided
if (nargin<3)
  runningshow=~nargout;
  if (nargin<2)
    preserve_degree=1;
    if (nargin<1)
      nv=128;
    end
  end
end

adj=eye(nv);                     % every vertex connected to itself 

adj = circshift(adj,-1)  + ...   % vertex connected to neighbour right
      circshift(adj,-2);         % vertex connected to neighbour 2 right

fa=find(adj);                    % find all original edges

% put vertices on circle in a good watts/strogatz way
for i=1:nv;                 
    xy(i,1)=sin(i/nv*2*pi); 
    xy(i,2)=cos(i/nv*2*pi);  
end  

% initialise some result vectors
relabellings=0;
clusterings=[];
pathlengths=[];
dist=[];

for f=0:length(fa)                          % cycle over all original (if not already re-labeleld) edges
        
  if (f)
    
  if (adj(fa(f)))                           % does edge still exist?
    
    [i,j]=ind2sub([nv nv],fa(f));           % find row and column: edge connects vertex i to j
    nli=[find(adj(i,:)) find(adj(:,i))'];   % neighbours of i
    
    empties=1:nv;
    empties([i nli])=[];
    
    if (~preserve_degree)
    
      % no degree preservation: we can choose an empty spot as i's new neighbour
      k=empties(ceil(rand*length(empties)));
      
      % now re-label
      adj(i,j)=0;
      adj(min(i,k),max(i,k))=1;

      fprintf('connecting vertex %04d to %04d instead of %04d\n',i,k,j);      
      relabellings=relabellings+1;
      
    else
      
      % degree preservation: we need to swap edges 
      % (i,j) and (i2,j2) to edges (i,j2) and (i2,j)
      nlj=[find(adj(j,:)) find(adj(:,j))']; % neighbours of j
      emptiesj=1:nv;
      emptiesj([j nlj])=[];
      
      empties=intersect(empties,emptiesj);  % points that are neighbour of neither i nor j      
      [i2s,j2s]=find(adj(empties,empties)); % any connections between them
      
      if (~isempty(i2s))
      
	i2=ceil(rand*length(i2s));
	j2=empties(j2s(i2));
	i2=empties(i2s(i2));
	
	adj(i,j)=0;
	adj(i2,j)=1;
	adj(i,j2)=1;
	adj(i2,j2)=0;
	
	fprintf('changing (%d,%d) and (%d,%d) to (%d,%d) and (%d,%d)\n',...
		i,j,i2,j2,i,j2,i2,j);
	relabellings=relabellings+1;
      
      end % if isempty
	
    end % if preserve_degree
      
	% /* from brain-connectivity-toolbox.net/bct/Home/functions/clustering_coef_bu.m            
	clus=zeros(nv,1);                    % clustering per vertex
	for v=1:nv
	  
	  nl=find(adj(v,:));             % find neighbours of vertex u
	  d=length(nl);                  % degree = |neighbours|
	  if (d>=2)                      % degree must be at least 2
	    la=adj(nl,nl);               % see how many neighbours of u are each other's neighbours
					 % (S is a |neighbours| x |neighbours| local adjacency matrix)
            clus(v)=sum(la(:))/(d^2-d);  % formula for C_i on en.wikipedia.org/wiki/Clustering_coefficient
	  end % if 
	  
	end % for v
	mC=mean(clus(:)); % average clustering
	
	%%
	%% compute average distance
	%%
	
	% /* from www.ee.columbia.edu/~marios/matlab/tips.pdf
	dist=1./adj;
	
	for c=1:nv
	  dist=min(dist,repmat(dist(:,c),[1 nv])+repmat(dist(c,:),[nv 1]));
    end % for c
	% from www.ee.columbia.edu/~marios/matlab/tips.pdf 
	mPL=mean(dist(dist~=inf)); % average shortest path
	
	% store this relabelling
	clusterings(relabellings)=mC;
	pathlengths(relabellings)=mPL;
    
  end % if (adj)
  
  end % if f
  
    if(runningshow)
      subplot(1,4,1)     
      if (prod(size(clusterings)))
      plot([clusterings(:)/max(clusterings) ...
	    pathlengths(:)/max(pathlengths) ...
	    (max(pathlengths)*clusterings(:))./(pathlengths(:)*max(clusterings))]) % see how the average clustering index changes
        set(gca,'xlim',[0 length(clusterings)])
      else
          cla;
      end
      subplot(1,4,2)
      gplot(adj,xy) % see the connections as a graph
      subplot(1,4,3)
      imagesc(adj); % see The Matrix
      subplot(1,4,4)
      if (prod(size(dist)))
        imagesc(dist,[0 nv/4]); % see the distance matrix
      else
          cla;
      end
      pause(.1)
    end
    
    %if (~f)
    %    
    %  for i=5:-1:1;
    %      fprintf('\t%d\r',i);
    %      pause(1);
    %  end       
    %  
    %end % if (f) -- used to record as an mpeg movie

end % for f
    
% normalise clusterings and pathlengths
clusterings=clusterings./clusterings(end)/nv*5; % normalise by value of completely re-labelled 'random' graph
pathlengths=pathlengths./pathlengths(end)/sqrt(nv)*5; % normalise by value of completely re-labelled 'random' graph
smallworldcoeffs=clusterings./pathlengths;

if (runningshow)
  subplot(1,4,1)
  semilogx([clusterings(:) pathlengths(:) smallworldcoeffs(:)]) % see how the average clustering index changes
  set(gca,'xlim',[1 length(clusterings)])
  subplot(1,4,2)
  gplot(adj,xy) % see the connections as a graph
  subplot(1,4,3)
  imagesc(adj); % see The Matrix
  subplot(1,4,4)
  imagesc(dist,[0 nv/4]); % see the distance matrix
end
 
return