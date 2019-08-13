clear; close all;

fi = @(ep,r)    exp(-(ep*r).^2);    % Gaussian; Radial function
u_Sol = @(x,t)  exp(-(20*(x-t)).^2);      % Analytic solution

n  = 30;                % Stencil size
ep = 20;                % Epsilon value
eu = 30;              % Shape of Gaussian initial condition
d = 4;                 %Degree of poynomials to attach
m = 5;                %Degree of PHS

xt = create_nodes_01;   % Create all the nodes. For changes to the node set see radius_variable_02.m

%Reorder node set
xtU = xt(xt(:,2)<=1 & xt(:,1)<=1,:);
[nt,~] = size(xtU);
xtU(1:nt,2) = xtU(1:nt,2);
xtU = sortrows(xtU,2);

u       = zeros(nt,1);  % Vector to hold computed solution
weights = zeros(nt,n);  % Array to hold computed weights for all stencils
p       = min(n/2-1,14);        % Number of trial weight sets for each node
W       = zeros(n,p,nt);% Array to hold all the temporary weight sets
        % Dimensions run over: (1) Weights, (2) Eval. points, (3) Stencils
ind2    = zeros(nt,n);  % Array to hold lower neighbors to node k
iBC     = zeros(nt,1);  % Array to mark nodes with given boundary values


for k = 1:nt            % Sequential loop over all nodes
                        % Create the weight sets to be optimized in a later
                        % parfor loop
    if xtU(k,1) < 1e-10 || xtU(k,2) < max(xtU(:,2))/10 || xtU(k,2)==1 || xtU(k,1)>(1-1e-10)  % Find k-values for bdy nodes
        u(k)   = u_Sol(xtU(k,1),xtU(k,2));     % Assign given boundary value
        iBC(k) = true;                      % Mark that this has been done
    else
        % Locate the n nearest neighbors below node k
        alpha = 1; %Wave speed and stencil streching
        if alpha>1
            sigma = 1/alpha;
            scaling = 1/sigma;
        else
            sigma = 1;
            scaling = 1;
        end
        %Finding neighbors
        [ind,dist]=knnsearch(xtU(xtU(:,2)<=xtU(k,2),:),xtU(k,:),'K',n,'Distance','mahalanobis','Cov',[1 0; 0 sigma^2]);
                          
        x = xtU(ind,1);  t = xtU(ind,2);  % Get stencil node locations x,t
        [~,indtemp] = sort(dist);   % Sort along row k to find nodes near
        %Center points
        xe      = x(1)+0.5*(x(2:p+1)-x(1));
        te      = t(1)+0.5*(t(2:p+1)-t(1));     
        
        % Calculate all the weights
        for P = 1:p
             w = RBF_TD_PHS_Transport_GitHub(x,t,xe(P),te(P),alpha,scaling,m,d);
             W(:,P,k) = w;
        end
        
        ind2(k,:) = ind(indtemp);
    end
end


for k = 1:nt            % Loop over all the nodes to find optimal comb.
% parfor k = 1:nt       % Use Parallel loop if cores available; 

    if ~iBC(k)          % If node not a boundary node, then optimize
        WI = W(:,:,k);
        
% l1decode version
        fact = 1./WI(1,:);
        W2 = bsxfun(@times,WI,fact);
        W2(:,2:p) = bsxfun(@minus,W2(:,2:p),W2(:,1));
        dv = l1decode_pd(zeros(p-1,1)/(p-1),W2(2:n,2:p),[],W2(2:n,1));
        wf = W2(:,1)-W2(:,2:p)*dv;
        weights(k,:) = wf';     % Store weights for actual time stepping       
        
    end
end

% Apply all the stencils for the actual time stepping
for k = 1:nt
    if ~iBC(k)
        u(k) = -weights(k,2:n)*u(ind2(k,2:n))/weights(k,1); % Appy stencil
    end
    
end

%Error in solution
er = u - u_Sol(xtU(:,1),xtU(:,2));
Error2 = norm(er,2)/sqrt(length(xtU(:,1)));
ErrorInf = norm(er,inf);

% Display result on original node set
figure
subplot(2,1,1)
gr = gray;
tri = delaunay(xtU(:,1),xtU(:,2)); trisurf(tri,xtU(:,1),xtU(:,2),u,'EdgeColor','none');
colormap(gr(end-20:end,:)); view([-30,40]); % Grey-to-white
hold on; plot3(xtU(:,1),xtU(:,2),u,'k.');     % Plot all points
xlabel('\itx'); ylabel('\itt');

% Edge on display
% figure
subplot(2,1,2)
plot3(xtU(:,1),xtU(:,2),u,'k.'); view([-45,0])    % View edge-on
xlabel('\itx'); ylabel('\itt');

%Display solution on original node set
figure
trisurf(tri,xtU(:,1),xtU(:,2),u);             % Show surface
colormap([1 1 1]); view([-16,60]);
xlabel('\itx'); ylabel('\itt');

figure
trisurf(tri,xtU(:,1),xtU(:,2),u_Sol(xtU(:,1),xtU(:,2)));             % Show surface
colormap([1 1 1]); view([-16,60]);
xlabel('\itx'); ylabel('\itt');

% Display error on original node set
figure
trisurf(tri,xtU(:,1),xtU(:,2),er);
colormap([1 1 1]); view([-16,60]); % zlim([-0.15,0.15])
xlabel('\itx'); ylabel('\itt');

% Display node set in the x,t-plane
figure
plot (xtU(:,1),xtU(:,2),'k.'); hold on
axis([0 1 0 1]); axis equal;
xlabel('\itx'); ylabel('\itt');
xd = xtU(logical(iBC),:);
plot (xd(:,1),xd(:,2),'ks');
box off
set(gca,'XAxisLocation','bottom')
set(gca,'YAxisLocation','left')
set(gca,'FontSize',18)
xlim([0,1])
ax = gca;
ax.TickLength = [0.025 0.1];

