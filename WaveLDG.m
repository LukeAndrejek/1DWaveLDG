% Implement a method from section 6.1 of the following paper:
% "Energy Conserving LDG Methods for Wave Propogation Problems"
clear all
% Polynomial degree
k = 2;
% Max time
tMax = 1;
% Max space
xMax = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use this section to input a desired Nx and dt/dx
Nx = 10;
dx = xMax/Nx;
Nt = 100 * Nx^2;
dt = tMax/Nt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Keep in mind that I_i = [(i-1)*delX, i*delX]
U = zeros((k+1)*Nx, Nt);
Q = zeros((k+1)*Nx, Nt);
% Note that polySine is the Taylor approximation of sin(pi*x)
trigDegree = 30;
polySine = [];
for p = 0:ceil((trigDegree-2)/2)
    polySine = [(-1)^p * pi^(2*p+1) / factorial(2*p+1), 0, polySine];
end
polyCosine = [1];
for p = 1:floor(trigDegree/2)
    polyCosine = [(-1)^p * pi^(2*p) / factorial(2*p), 0, polyCosine];
end

% Build matrices

% the ith row of polyBasis is ( (x-x_i)/delta(x_i) )^i
% we begin with x_i = 0
polyBasis = zeros(k+1,k+1);
for i = 1:k+1
    polyRowI = [1];
    if i > 1
        for j = 1:i-1
            polyRowI = conv(polyRowI, [1/dx -1/2]);
        end
    end
    polyBasis(i,:) = [zeros(1,k+1-length(polyRowI)) polyRowI];
end
oldPolyBasis = polyBasis;

M = zeros(k+1,k+1);
A = M;
B = M;
C = M;
D = M;
E = M;
for m = 0:k
    for l = 0:k
        M(m+1,l+1) = diff(polyval(polyint(conv(polyBasis(m+1,:),...
            polyBasis(l+1,:))),[0 dx]));
        A(m+1,l+1) = diff(polyval(polyint(conv(polyBasis(l+1,:),...
            polyder(polyBasis(m+1,:)))),[0 dx]));
        B(m+1,l+1) = 2^(-l-m);
        C(m+1,l+1) = (-1)^(m);
        D(m+1,l+1) = (-1)^(l+m);
        E(m+1,l+1) = (-1)^(l);
    end
end
G = M;
G(k+1,:) = ((-ones(1,k+1)) .^ (0:k)) .* ((2*ones(1,k+1)) .^ (0:-1:-k));
I = eye(Nx);
P = eye(Nx-1);
P = [zeros(1,Nx-1); P];
P = [P zeros(Nx,1)];
P(1,Nx) = 1;
blockM = kron(I, M);
block1 = kron(I, A-B) + kron(P, C.*B);
block2 = kron(I,-A-D.*B) + kron(P', E.*B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the values of U
for j = 1:Nt
    % Take care of initial conditions manually
    if j == 1
        for n = 1:Nx
            % Prepare matrix involving the initial condition
            
            % Set correct polynomial basis
            for l = 1:k+1
                polyRowI = [1];
                if l > 1
                    for m = 1:l-1
                        polyRowI = conv(polyRowI, [1/dx -n+1/2]);
                    end
                end
                polyBasis(l,:) = [zeros(1,k+1-length(polyRowI)) ...
                    polyRowI];
            end
            
            % Find U
            H = zeros(k+1,1);
            for l = 1:k
                H(l,1) = diff(polyval(polyint(conv(polySine,...
                    polyBasis(l,:))), [dx*(n-1) dx*n]));
            end
            H(k+1,1) = polyval(polySine, dx*(n-1));
            U((k+1)*n-k : (k+1)*n, 1) = G\H;
        end
        
        % Find Q
        Q(:,1) = blockM \ (block2*U(:,1));
        
    elseif j == 2
        
        % Find U
        for n = 1:Nx
            for l = 1:k+1
                polyRowI = [1];
                if l > 1
                    for m = 1:l-1
                        polyRowI = conv(polyRowI, [1/dx -n+1/2]);
                    end
                end
                polyBasis(l,:) = [zeros(1,k+1-length(polyRowI)) ...
                    polyRowI];
            end
            
            % Create polynomial in standard basis
            currentElement = (k+1)*n-k;
            Ustd = zeros(1,k+1);
            for l = 1:k+1
                Ustd = Ustd + U(currentElement+l-1) * polyBasis(l,:);
            end
            currentElement = (k+1)*n-k;
            if k > 0
                for p = 2:2:2*ceil(k/2)
                    for l = 1:k+1
                        % Take necessary derivatives of appropriate basis
                        R = polyBasis(l,:);
                        for deriv = 1:p
                            R = polyder(R);
                        end
                        R = [zeros(1,length(Ustd)-length(R)), flip(R)];
                        % Update standard basis polynomial
                        Ustd = Ustd + (U(currentElement+(k+1)-l, 1) ...
                            * dt^p / factorial(p)) * R;
                    end
                end
            end
            % Change from standard basis to finite element basis
            for s = k+1:-1:1
                a = 0;
                if s < k+1
                    for l = s+1:k+1
                        a = a + U(currentElement+l-1, 1) ...
                            * polyBasis(l,k+2-s);
                    end
                end
                U(currentElement+s-1, 2) = (Ustd(k+2-s) - a) ...
                    / polyBasis(s, k+2-s);
            end
        end
        % Find Q     
        Q(:,2) = blockM \ (block2*U(:,2));
        
    else
        % Run induction step
        U(:,j) = blockM \ (2*blockM*U(:,j-1) ...
            - blockM*U(:,j-2) - dt^2*block1*Q(:,j-1));
        Q(:,j) = blockM \ (block2*U(:,j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the result and compute error
plotFrames = 50;
plotGap = Nt/plotFrames;
X = [];
gaussPoints = ceil((k+1)/2) + 2;
plotPoints = k + 4;
exactPoints = 10;
for j = 1:Nt
    if j/plotGap == ceil(j/plotGap)
        [Xg, W] = gauss_points_and_weights(gaussPoints);
        X = [];
        Y = zeros(plotPoints*Nx, 1);
        L2error = 0;
        for n = 1:Nx
            % Set correct polynomial basis
            for l = 1:k+1
                polyRowI = [1];
                if l > 1
                    for m = 1:l-1
                        polyRowI = conv(polyRowI, [1/dx -n+1/2]);
                    end
                end
                polyBasis(l,:) = [zeros(1,k+1-length(polyRowI)) ...
                    polyRowI];
            end
            
            % Prepare for error computation
            aX = (n-1)*dx;
            bX = n*dx;
            Xe = ((bX-aX)/2)*Xg + ((bX+aX)/2)*ones(1,length(Xg));
            GQx = 0;
            X = [X (n-1)*dx+(1e-10)];
            X = [X (n-1)*dx+(dx/(plotPoints-1)) ...
                : dx/(plotPoints-1) : n*dx-(dx/(plotPoints-1))];
            X = [X n*dx-(1e-10)];
            
            % Obtain error data
            for m = 1:gaussPoints
                V = 0;
                for l = 1:k+1
                    V = V + U((k+1)*n+l-(k+1), j)...
                        * polyval(polyBasis(l,:),Xe(m));
                end
                GQx = GQx + W(m) * (V ...
                    - polyval( polyval(polyCosine,dt*(j-1)) ...
                    * polySine, Xe(m)) )^2;       
            end
            
            % Obtain plotting data
            for m = 1:plotPoints
                Y(plotPoints*n-m+1,1) = 0;
                currentElement = (k+1)*n-k;
                for l = 1:k+1
                    Y(plotPoints*n-m+1,1) = ...
                        Y(plotPoints*n-m+1,1) ...
                        + U(currentElement+l-1, j) ...
                        * polyval(polyBasis(l,:), ...
                        X(plotPoints*(n-1)+m));
                end
            end
            L2error = L2error + (dx/2)*GQx;
            % Flip Y
            Y(plotPoints*(n-1)+1:plotPoints*n,1) = flip( ...
                Y(plotPoints*(n-1)+1:plotPoints*n,1) );
        end
        L2error = sqrt(L2error);
        
        % Plot the function
        clf
        plot(X,Y)
        title(['t = ' num2str((j-1)*dt) ', err = ' num2str(L2error)])
        ylim([-1.1 1.1])
        hold on
        Z = 0:dx/exactPoints:xMax;
        plot(Z,polyval( polyval(polyCosine,dt*(j-1)) * polySine, Z))
        pause(0.001)
    end
end