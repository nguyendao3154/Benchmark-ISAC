%%
V_psd(:,:,1) = (V(:,:,1) + V(:,:,1)') / 2; % Ensure the matrix is symmetric
V_psd(:,:,1) = V(:,:,1) + 64 * eye(64); % Ensure it is positive definite if necessary
V_psd(:,:,2) = (V(:,:,2) + V(:,:,2)') / 2; % Ensure the matrix is symmetric
V_psd(:,:,2) = V(:,:,2) + 64 * eye(64); % Ensure it is positive definite if necessary
[eigvec(:,:,1), eigval(:,:,1)] = eig(V_psd(:,:,1));
[eigvec(:,:,2), eigval(:,:,2)] = eig(V_psd(:,:,2));
num_samples = 1000;
random_vectors = zeros(64, 2, num_samples);

% Generate random vectors from the Gaussian distribution
for i = 1:num_samples
    z = randn(64, 2); % Generate a standard normal random vector
    random_vectors(:, 1, i) = (eigvec(:,:,1) * sqrt(eigval(:,:,1)) * z(:,1))'; % Transform the vector
    random_vectors(:, 2, i) = (eigvec(:,:,2) * sqrt(eigval(:,:,2)) * z(:,2))'; % Transform the vector
end
for i = 1:num_samples
    random_vectors(:, 1, i) = random_vectors(:, 1, i) / norm(random_vectors(:, 1, i));
    random_vectors(:, 2, i) = random_vectors(:, 2, i) / norm(random_vectors(:, 2, i));
end

for i = 1:num_samples
    objective_values(i) = compute_capacity(DL.channel, abs(random_vectors(:,:,i)), angle(random_vectors(:,:,i))); % Example objective function
end
[best_OBJ, best_index] = max(objective_values);
best_solution = random_vectors(:, best_index);

result_SINR(1) = pow2db(DL_SINR(DL.channel, abs(random_vectors(:,:,best_index)), angle(random_vectors(:,:,best_index)), 1));
result_SINR(2) = pow2db(DL_SINR(DL.channel, abs(random_vectors(:,:,best_index)), angle(random_vectors(:,:,best_index)), 2));