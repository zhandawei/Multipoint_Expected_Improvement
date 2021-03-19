function y = Infill_qEI(new_x,model,fmin)
warning off;
% reshape the input
if size(new_x,2) > size(model.S,2)
    new_x = reshape(new_x,size(model.S,2),[])';
end
% delete duplicated rows
new_x = unique(new_x,'rows');
% number of infill samples
q = size(new_x,1);
% predictions and covarince matrix
[u,s,~,Cov] = predictor(new_x,model);
mu = u';
sigma = Cov;

symetric_term = zeros(q,q);
nonsymetric_term = zeros(q,q);
pk = zeros(1,q);
first_term = zeros(1,q);
second_term = zeros(1,q);
if q == 1
    y = (fmin-u).*Gaussian_CDF((fmin-u)./s)+s.*Gaussian_PDF((fmin-u)./s);
else
    for k = 1:q
        % Sigma_k: covariance matrix of vector Z^{k}
        Sigma_k = sigma;
        for i = 1:q
            if i ~= k
                Sigma_k(i,k) = sigma(k,k)-sigma(k,i);
                Sigma_k(k,i) = Sigma_k(i,k);
            end
        end
        for i = 1:q
            for j = 1:q
                if i~=k && j~=k
                    if j>=i
                        Sigma_k(i,j) = sigma(i,j) + sigma(k,k) - sigma(k,i) - sigma(k,j);
                    else
                        Sigma_k(i,j) =   Sigma_k(j,i);
                    end
                end
            end
        end
        % mu_k: mean of vector Z^{k}
        mu_k = mu(k) - mu;
        mu_k(k) = mu(k);
        % vector b_k
        b_k = zeros(1,q);
        b_k(k) = fmin;
        % calcualte pk
        Sigma_k = Sigma_k + diag(ones(1,q)*1E-20);
        [~,p] = chol(Sigma_k);
        if p > 0
            [L, DMC, P] = modchol_ldlt(Sigma_k);
            Sigma_k = P'*L*DMC*L'*P;
        end
        %pk(k) = mvncdf(b_k - mu_k,zeros(1,q),Sigma_k);
        pk(k) = qsimvnv(200*q, Sigma_k, -inf*ones(q,1), (b_k - mu_k)');
        % replace NaN by 0
        if isnan(pk(k))
            pk(k) = 0;
        end
        first_term(k) = (fmin - mu(k))*pk(k);
        
        nonsymetric_term(k,:) = Sigma_k(:,k)';
        for i = 1:q
            if i >= k
                mik = mu_k(i);
                sigma_ii_k = Sigma_k(i,i);
                bik = b_k(i);
                phi_ik =  normpdf(bik,mik,sqrt(sigma_ii_k));
                % calculate c.i^(k)
                sigmai = Sigma_k(i,:)/Sigma_k(i,i);
                cik = (b_k - mu_k) - (b_k(i) - mu_k(i))*sigmai;
                cik(i) = [];
                % calculate sigma.i^(k)
                sigmaik = Sigma_k;
                for uu = 1:q
                    for vv = 1:q
                        if uu~=i && vv~=i
                            if Sigma_k(i,i) == 0
                                sigmaik(uu,vv) = Sigma_k(uu,vv);
                            else
                                sigmaik(uu,vv) = Sigma_k(uu,vv) - Sigma_k(uu,i)*Sigma_k(vv,i)/Sigma_k(i,i);
                            end
                        else
                            sigmaik(uu,vv) = 0;
                        end
                    end
                end
                sigmaik(i,:)=[];
                sigmaik(:,i)=[];
                % calculate phi_ik
                [~,p] = chol(sigmaik);
                if p > 0
                    [L, DMC, P] = modchol_ldlt(sigmaik);
                    sigmaik = P'*L*DMC*L'*P;
                end
                %Phi_ik = mvncdf(cik,zeros(1,q-1),sigmaik);
                Phi_ik = qsimvnv(200*(q-1), sigmaik, -inf*ones(q-1,1),cik');
                % replace NaN by 0
                if isnan(Phi_ik)
                    Phi_ik = 0;
                end
                symetric_term(k,i) = phi_ik*Phi_ik;
            else
                symetric_term(k,i) = symetric_term(i,k);
            end
        end
        second_term(k) = sum(nonsymetric_term(k,:).*symetric_term(k,:));
        
    end
    y = sum(first_term + second_term);
end
