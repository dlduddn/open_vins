%% Vectorization

function [xhat, N_eff] = sir(x0, P0, dSE2, z, sqrtQ, R, len, N)

    % Allocation
    stateDim = 3;
    xhat = zeros(stateDim, len);
    N_eff = zeros(1,len);
    
    % Initial Particle
    sqrtP0 = sqrt_psd(P0, stateDim);
    xi = repmat(x0, 1, N) + sqrtP0 * randn(stateDim, N);
    xi(3,:) = wrap_angle(xi(3,:));
    wi = (1/N) * ones(1,N);

    xhat(1:2,1) = sum(xi(1:2,:) .* wi, 2);
    xhat(3,1) = atan2(sum(sin(xi(3,:)) .* wi), sum(cos(xi(3,:)) .* wi));
    N_eff(1,1) = 1 / sum(wi.^2);
    
    % Particle filtering
    for t=2:len
        %% Propagation
        dX = dSE2(1, t-1);
        dY = dSE2(2, t-1);
        dYaw = dSE2(3, t-1);

        % dX, dY are body_{k-1}-frame controls. xi is in mapLocal frame.
        c = cos(xi(3,:));
        s = sin(xi(3,:));
        xi(1,:) = xi(1,:) + c .* dX - s .* dY;
        xi(2,:) = xi(2,:) + s .* dX + c .* dY;
        xi(3,:) = wrap_angle(xi(3,:) + dYaw);

        if any(sqrtQ(:))
            xi = xi + sqrtQ * randn(stateDim, N);
            xi(3,:) = wrap_angle(xi(3,:));
        end

        xhat(1:2,t) = sum(xi(1:2,:) .* wi, 2);
        xhat(3,t) = atan2(sum(sin(xi(3,:)) .* wi), sum(cos(xi(3,:)) .* wi));
        N_eff(1,t) = 1 / sum(wi.^2);

        %% Importance Weight Update & Normalization
        % % Likelihood evalutation
        % zhat = [atan2(xi(2,:),xi(1,:)); sqrt(xi(1,:).^2 + xi(2,:).^2)];
        % innovation = z(:,t) - zhat;
        % 
        % % likelihood = exp(-0.5 * sum((R \ innovation) .* innovation, 1)); % Fast
        % % likelihood = exp(-0.5 * sum((inv(R)*innovation) .* innovation, 1)); % Faster
        % 
        % % Update & Normalization
        % wi = likelihood;
        % wsum = sum(wi);
        % wi = wi./wsum;
        % 
        % % Estimate 
        % xhat(:,t) = sum(xi.*wi,2);
        % 
        % % Resampling
        % N_eff(1,t) = 1/sum(wi.^2);
        % if (1 / sum(wi.^2)) < N/2
        %     idx = sysresample(wi);
        %     xi = xi(:, idx);
        % end 
    end

end

function angle = wrap_angle(angle)
    angle = atan2(sin(angle), cos(angle));
end

function sqrtP = sqrt_psd(P, stateDim)
    if ~isequal(size(P), [stateDim, stateDim])
        error('P0 must be a %d x %d covariance matrix.', stateDim, stateDim);
    end

    P = (P + P') / 2;
    [V, D] = eig(P);
    d = max(diag(D), 0);
    sqrtP = V * diag(sqrt(d));
end
