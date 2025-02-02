
classdef k_means_cluster
    methods (Static)
        function [] = main()
            %MAIN Summary of this function goes here
            %   Detailed explanation goes here
            prob = 2;
            if prob == 3
                k_means_cluster.svm_rbf();
            else
                X = load('/Users/vamshimuthineni/Vamshi/ML/assign5/hw5data/digit/digit.txt');
                Y = load('/Users/vamshimuthineni/Vamshi/ML/assign5/hw5data/digit/labels.txt');

                X_Y = horzcat(X,Y);

                k=2;
                % this function performs k-means classification and returns within groupof
                % sum of squares and pair-counting measures.
                ks = [2,4,6];
                [n,d] = size(X);
                for k_i = 1:size(ks,2)
                    k = ks(k_i);
                    [centers,ss,clusters,iters] = k_means_cluster.k_means_clustering(X,k,'False');
                    [p1,p2,p3] = k_means_cluster.find_p3(n, Y, clusters);
                    k_means_cluster.print_outputs(ss,p1,p2,p3,k,iters);
                end
            end
        end
        
        function[] = prob_2_5_3_and_4()
                
            X = load('/Users/vamshimuthineni/Vamshi/ML/assign5/hw5data/digit/digit.txt');
            Y = load('/Users/vamshimuthineni/Vamshi/ML/assign5/hw5data/digit/labels.txt');


            ks = [1,2,3,4,5,6,7,8,9,10];
            [n,d] = size(X);
            p1s = [];
            p2s = [];
            p3s = [];
            ss_s = [];
            for k_i = 1:size(ks,2)
                k = ks(k_i);
                ss_sum = 0;
                p1_sum = 0;
                p2_sum = 0;
                p3_sum = 0;
                for iters = 1:10
                    [centers,ss,clusters,iters] = k_means_cluster.k_means_clustering(X, k, 'True');
                    [p1,p2,p3] = k_means_cluster.find_p3(n, Y, clusters);
                    ss_sum = ss_sum + ss;
                    p1_sum = p1_sum + p1;
                    p2_sum = p2_sum + p2;
                    p3_sum = p3_sum + p3;
                end
                p1s = [p1s,p1_sum/10];
                p2s = [p2s,p2_sum/10];
                p3s = [p3s,p3_sum/10];
                ss_s = [ss_s,ss_sum/10];
            end
            k_means_cluster.plot_graph_2_5_3(ks,ss_s)
            k_means_cluster.plot_graph_2_5_4(ks, p1s, p2s, p3s)
        end
        
        function [] = plot_graph_2_5_3(ks,ss_s)
            figure
            plot(ks,ss_s,'Color','blue')

            title('Plot of total within group sum of squares versus k')
            xlabel('k-values')
            ylabel('total within sum of squares(*e+08)')  
        end
        
        function plot_graph_2_5_4(ks, p1s, p2s, p3s)
            plot(ks, p1s, 'DisplayName','p1')
            hold on
            plot(ks, p2s, 'DisplayName','p2')
            plot(ks, p3s, 'DisplayName','p3')
            hold off
            title('Plot of p1,p2,p3 versus k')
            xlabel('k-values')
            ylabel('pair-counting measures')
            lgd = legend;
            lgd.NumColumns = 1;
        end
        
        
        function[] = cvv_svm_rbf()
            data = load('q3_variable_data');
            cvv_accuracy = svmtrain(data.trLbs, data.trD','-v 5 -t 2');
            fprintf("five-fold cross validation accuracy using rbf kernel %f", cvv_accuracy);
        end
        
        function [] = svm_rbf()
            data = load('q3_variable_data');
            Cs = [3000,5000,5000,7000,8000,10000,10000];
            gs = [0.7,0.7,0.9,0.9,0.9,1,2];
            for i = 1:size(Cs,2)
                S = sprintf('-q -v 5 -t 2 -c %d -g %d', Cs(i), gs(i));
                cv_accuracy = svmtrain(data.trLbs, data.trD',S);
                fprintf("cv_accuracy at C:%d, g:%d is %f", Cs(i),gs(i), cv_accuracy);
            
%             model = svmtrain(data.trLbs, data.trD','-c 10000 -g 2');
%             [predicted_label, accuracy, decision_values] = svmpredict(data.trLbs, data.trD', model);
%             fprintf("\n training accuracy is %f", accuracy);
%             trLbs = zeros(size(data.tstD,2),1);
%             [predicted_label] = svmpredict(data.tstIds', data.tstD', model);
            end
        end
        
        function[] = train_expX2()
            data = load('q3_variable_data');
            [trDK, tstDK] = k_means_cluster.cmpExpX2Kernel(data.trD', data.tstD', 2);
%             disp(size(trDK));
%             Cs = [9,10,11,12,13,14,15,16,17,18,19,20];
% %             Cs = [power(2,0),power(2,1),power(2,2),power(2,3),power(2,4),power(2,5),power(2,6), power(2,7),power(2,8),power(2,9),power(2,10),power(2,11)];
%             gs = [power(2,-6),power(2,-5),power(2,-4),power(2,-3),power(2,-2),power(2,-1),power(2,0),power(2,1),power(2,2),power(2,3),power(2,4),power(2,5)];
%             for i = 1:size(Cs,2)
%                 for j = 1:size(gs,2)
%                     S = sprintf('-q -v 5 -t 4 -c %d -g %d', Cs(i), gs(j));
%                     cv_accuracy = svmtrain(data.trLbs, [(1:size(trDK,1))',trDK], S);
%                     fprintf("\n cv_accuracy at C:%d, g:%d is %f", Cs(i),gs(j), cv_accuracy);
%                 end
%             end
    
            cv_accuracy = svmtrain(data.trLbs, [(1:size(trDK,1))',trDK],'-q -v 5 -t 4 -c 17 -g 1');
            disp(cv_accuracy);
            
            model = svmtrain(data.trLbs,[(1:size(trDK,1))',trDK],'-q -t 4 -c 17 -g 1');
            tstLbs = zeros(1600,1);
            [predicted_label, accuracy, decision_values] = svmpredict((1:size(tstDK,1))', [(1:size(tstDK,1))',tstDK], model);
        end
        
        function [trainDK, testDK] = cmpExpX2Kernel(trainD, testD, gamma)
            trainDK = zeros(size(trainD,1),size(trainD,1));
            for i=1:size(trainD,1)
                for j=1:size(trainD,1)
                    sum = 0;
                    for k=1:size(trainD,2)
                        sum = sum + power(trainD(i,k) - trainD(j,k), 2)/((trainD(i,k)+trainD(j,k) + eps));
                    end
                    trainDK(i,j) = exp((-1/gamma)*sum);
                end
            end
            
            testDK = zeros(size(testD));
            for i=1:size(testD,1)
                for j=1:size(trainD,1)
                    sum = 0;
                    for k=1:size(testD,2)
                        sum = sum + power(testD(i,k) - trainD(j,k), 2)/((testD(i,k)+trainD(j,k) + eps));
                    end
                    testDK(i,j) = exp((-1/gamma)*sum);
                end
            end
            
        end
        
        
        function [] = print_outputs(ss,p1,p2,p3,k,iters)
            fprintf("\n-------------");
            fprintf("\n breaking for k=%d at: %d",k, iters);
            fprintf("\n Group sum of squares at k=%d is %f", k, ss);
            fprintf("\n Pair counting measures at k=%d are (p1,p2) %f, %f",k,p1,p2);
            fprintf("\n Pair counting measure at k=%d is (p3) %f",k,p3);
            fprintf("\n-------------");
        end

        function sum = calculate_dist(X_i, C_i,d)
            sum = 0;
            for dim=1:d
                sum = sum + (X_i(1,dim) - C_i(1,dim))^2;
            end
            sum = sqrt(sum);
        end
        
        function [centers,ss,clusters,iters] = k_means_clustering(X, k, rand_center)
            [n,d] = size(X);
            if strcmp(rand_center, 'True')
                centers = X(randperm(size(X, 1), k), :);
            else
                centers = X(1:k,:);
            end
            
        %     disp(size(centers));
            dists = [];
            clusters = zeros(n,1);
            old_clusters = zeros(n,1);
            for iters = 1:30
                old_clusters = clusters;
                for i=1:n
                    for j=1:k
                       dists(j,:) = k_means_cluster.calculate_dist(X(i,:),centers(j,:),d);
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
                ss(j,:) = ss(j,:) + k_means_cluster.calculate_dist(X(i,:), centers(j,:),d);
            end
        %     fprintf("displaying ss");
        %     disp(size(ss));
        %     disp(ss);
            ss = sum(sum(ss));
        %     disp(ss);
        %     fprintf("finished showing ss");
        %     disp(p1/same);
        %     disp(p2/diff);
        %     disp(p3);
        end

        function [p1,p2,p3] = find_p3(n, Y, clusters)
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
            p1 = p1/same;
            p2 = p2/diff;
        end
    end
end