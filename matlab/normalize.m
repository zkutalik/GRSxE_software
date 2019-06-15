function norm_matrix = normalize(matrix,dim)
warning off all

if nargin<2
    dim=1;
end

M = isnan(matrix);
mm = matrix;
mm(M==1) = 0;

sum1M = sum(1-M,dim);
sel = find(sum1M==0);
mea = sum(mm,dim)./sum1M;
mea(sel) = 0;

if dim==1
    for j=1:size(matrix,2)
        matrix(:,j) = matrix(:,j) - mea(j);
    end
    
%     matrix = matrix - repmat(mea,1,size(matrix,2));
else
%     matrix = matrix - repmat(mea,size(matrix,1),1);
    for j=1:size(matrix,1)
        matrix(j,:) = matrix(j,:) - mea(j);
    end
end    

mm = matrix.^2;
mm(M==1) = 0;
sel = find(sum1M<=1);
vari = sqrt(sum(mm,dim)./(sum1M-1));
vari(sel) = NaN;

if dim==1
%     matrix = matrix/repmat(vari,1,size(matrix,2));
    
    for j=1:size(matrix,2)
        matrix(:,j) = matrix(:,j)/vari(j);
    end
else
%     matrix = matrix/repmat(vari,size(matrix,1),1);

    for j=1:size(matrix,1)
        matrix(j,:) = matrix(j,:)/vari(j);
    end
end

norm_matrix = matrix;
