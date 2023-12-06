function [res, F_normalized] = demo_main(KH, numClust, Y, r1Temp, r2Temp)
lambda1 = r1Temp;
lambda2 = r2Temp;
numSample = size(KH, 1);
numker = size(KH,3);

stream = RandStream.getGlobalStream;
reset(stream);

for p = 1 : numker
    H(:,:,p) = orth(rand(numSample,numClust));
    G(:,:,p) = H(:,:,p);
end

for p = 1 : numker
    R(:,:,p) = eye(numClust);
    W(:,:,p) = eye(numClust);
end
gamma = ones(numker,1)/(numker);
maxIter = 20;

HRGW = zeros(numSample,numClust);

flag = 1;
iter = 0;
while flag
    %% Update F
    for p = 1 : numker
        HRGW = HRGW + gamma(p)*(H(:,:,p)*R(:,:,p) + G(:,:,p)*W(:,:,p));
    end
    iter = iter +1;
    [U,~,V] = svd(lambda2 * HRGW,'econ');
    F = U * V';
    
    %% Update H
    for p = 1 : numker
        temp = KH(:,:,p) * G(:,:,p) + lambda1 * G(:,:,p) + lambda2 * gamma(p) * F * R(:,:,p)';
        [U,~,V] = svd(temp,'econ');
        H(:,:,p) = U * V';
    end
    
    %% Update G
    for p = 1 : numker
        temp = KH(:,:,p)' * H(:,:,p) + lambda1 * H(:,:,p) + lambda2 * gamma(p) * F * W(:,:,p)';
        [U,~,V] = svd(temp,'econ');
        G(:,:,p) = U * V';
    end
    
    %% Update R
    for p = 1 : numker
        temp = lambda2 * H(:,:,p)' * gamma(p) * F;
        [U,~,V] = svd(temp,'econ');
        R(:,:,p) = U * V';
    end
    
    %% Update W
    for p = 1 : numker
        temp = lambda2 * G(:,:,p)' * gamma(p) * F;
        [U,~,V] = svd(temp,'econ');
        W(:,:,p) = U * V';
    end
    
    %% Update gamma
    coef = zeros(1, numker);
    for p = 1 : numker
        coef(1,p) = trace(F' * (H(:,:,p) * R(:,:,p) + G(:,:,p) * W(:,:,p)));
    end
    gamma = coef/norm(coef,2);
    
    %% Cal obj
    obj_val = 0;
    for p = 1 : numker
        obj_val = obj_val + trace(H(:,:,p)' * KH(:,:,p) * G(:,:,p)) + lambda1*trace(H(:,:,p)'*G(:,:,p)) + ...
            lambda2 * trace(F'*gamma(p)*(H(:,:,p) * R(:,:,p) + G(:,:,p) * W(:,:,p)));
    end
    obj(iter) = obj_val;
    
    if (iter>2) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-6 || iter>maxIter)
        flag = 0;
    end
end
F_normalized = F ./ repmat(sqrt(sum(F.^2, 2)), 1,numClust);
for it = 1 : 10
    res_label = litekmeans(F_normalized, numClust, 'MaxIter',100, 'Replicates',50);
    result(it,:) = Clustering8Measure(Y, res_label);
end
res = [mean(result);std(result)];
end

