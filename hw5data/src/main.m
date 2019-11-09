function [] = main()
%MAIN Summary of this function goes here
%   Detailed explanation goes here


X = load('/Users/vamshimuthineni/Vamshi/ML/assign5/hw5data/digit/digit.txt');
Y = load('/Users/vamshimuthineni/Vamshi/ML/assign5/hw5data/digit/labels.txt');

X_Y = horzcat(X,Y);

k=2;
% this function performs k-means classification and returns within groupof
% sum of squares and pair-counting measures.
[ss,p1,p2,p3,iters] = k_means_clustering(X,Y,k);
print_outputs(ss,p1,p2,p3,k,iters);

k=4;
[ss,p1,p2,p3,iters] = k_means_clustering(X,Y,k);
print_outputs(ss,p1,p2,p3,k,iters);

k=6;
[ss,p1,p2,p3,iters] = k_means_clustering(X,Y,k);
print_outputs(ss,p1,p2,p3,k,iters);


end

function [] = print_outputs(ss,p1,p2,p3,k,iters)
    fprintf("\n-------------");
    fprintf("\n breaking for k=%d at: %d",k, iters);
    fprintf("\n Group sum of squares at k=%d is %f", k, ss);
    fprintf("\n Pair counting measures at k=%d is %f",k,p3);
    fprintf("\n-------------");
end

function sum = calculate_dist(X_i, C_i,d)
    sum = 0;
    for dim=1:d
        sum = sum + (X_i(1,dim) - C_i(1,dim))^2;
    end
    sum = sqrt(sum);
end

function [ss,p1,p2,p3,iters] = k_means_clustering(X,Y,k)
    [n,d] = size(X);

    centers = [];
    for i=1:k
        centers = [centers;X(i,:)];
    end
%     disp(size(centers));
    dists = [];
    clusters = zeros(n,1);
    old_clusters = zeros(n,1);
    for iters = 1:30
        old_clusters = clusters;
        for i=1:n
            for j=1:k
               dists(j,:) = calculate_dist(X(i,:),centers(j,:),d);
            end
        %     disp(dists);
            [~, min_index] = min(dists);
            clusters(i,:) = min_index;
        end
%         disp(size(clusters));
%         disp(clusters);
        
        if clusters == old_clusters | iters==20
            break;
        end

        new_center = zeros(size(centers));
        counts = zeros(k,1);
%         disp(size(centers));
        for i=1:n
            j = clusters(i,:);
            new_center(j,:) = new_center(j,:) + X(i,:);
            counts(j,:) = counts(j,:) + 1;
        end

        for i=1:k
            centers(i,:) = (new_center(i,:))/counts(i,:);
        end
    %     disp(centers);
    end

    ss = zeros(k,1);
    for i=1:n
        j = clusters(i,:);
        ss(j,:) = ss(j,:) + calculate_dist(X(i,:), centers(j,:),d);
    end
%     fprintf("displaying ss");
%     disp(size(ss));
%     disp(ss);
    ss = sum(sum(ss));
%     disp(ss);
%     fprintf("finished showing ss");
    p1 = 0;
    p2 = 0;
    same = 0;
    diff = 0;

    for i=1:n
        for j=i+1:n
            if Y(i) == Y(j)
                same = same + 1;
                if clusters(i,:) == clusters(j,:)
                    p1 = p1+1;
                end
            else
                diff = diff + 1;
                if clusters(i,:) ~= clusters(j,:)
                    p2 = p2+1;
                end
            end
        end
    end

    p3 = (p1/same + p2/diff)/2;
%     disp(p1/same);
%     disp(p2/diff);
%     disp(p3);
end
