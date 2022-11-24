close all

% q3 = q2;
q3 = stableq2;

q4(:,:,:,:,:,1) = reshape(q3,pts,pts,pts,pts,pts);
q5 = zeros(size(q4));

for i = 1:length(hset)
    for j = 1:length(eccSet)
        for k = 1:length(dlaserSet)
            for l = 1:length(bEllSet)
                for m = 1:length(Lset)
%                     q5(m,l,k,1,i) = q4(m,l,k,1,i); %ideal circle
                    q5(m,l,k,j,i) = q4(m,l,k,j,i);
                end
            end
        end
    end
end


[top6,top6I] = maxk(q5(:),6); %look up top 6 values
% top7 = maxk(q2,7);
% out = find(q2 > top7(7)); %find the values less than top 6 values
out = top6I;
% out2 = q5(out);
% for i = 1:6
% out3(i) = find(q3 == out2(i),1) ;
% end
% [o1,o2,o3,o4,o5] = ind2sub([length(Lset),length(Lset),length(Lset),length(Lset),length(Lset)],out3);
[o1,o2,o3,o4,o5] = ind2sub([length(Lset),length(Lset),length(Lset),length(Lset),length(Lset)],out);

for i = 1:6
outliers6(i,:) = [o1(i),o2(i),o3(i),o4(i),o5(i)];
end
costsOut = q1(out);
outliers6

figure(22)
% sgtitle('Outlier (non turning point) Maxima')
for i = 1:6
    subplot(2,3,i)
    contourf(q1(:,:,outliers6(i,3),outliers6(i,4),outliers6(i,5)),[0:0.5:2.5 2.5:0.1:2.7 2.7:0.02:3],'ShowText','on')
    colorbar
%     colormap('jet')
    str = sprintf('L vs MajorAxis ||| Fixed dlaser = %0.2f | ecc = %0.2f | h = %0.2f ', dlaserSet(outliers6(i,3)),eccSet(outliers6(i,4)),hset(outliers6(i,5)));
    title(str)
    xlabel('Major Axis')
    xticks(1:length(bEllSet))
    xticklabels(num2cell(round(bEllSet),2))
    ylabel('L')
    yticks(1:length(Lset))
    yticklabels(num2cell(round(Lset),2))
    caxis([1.5 3])
    hold on
    scatter(outliers6(i,2),outliers6(i,1),75,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)
    hold off
end

% figure(33)
% [M,c] = contourf(q1(:,:,outliers6(3,3),outliers6(3,4),outliers6(3,5)),[0:0.5:1.5 2:0.1:2.5 2.5:0.05:2.7 2.7:0.03:3],'ShowText','on');
% c.LineWidth = 1;
% clabel(M,c,'FontSize',16);
% colorbar
% colormap('jet')
% str = sprintf('L vs MajorAxis | dlaser = %0.2f | ecc = %0.2f | h = %0.2f ', dlaserSet(outliers6(3,3)),eccSet(outliers6(3,4)),hset(outliers6(3,5)));
% title(str,'FontSize',40)
% xlabel('Major Axis','FontSize',32)
% xticks(1:length(bEllSet))
% xticklabels(num2cell(round(bEllSet),2))
% ylabel('L','FontSize',32)
% yticks(1:length(Lset))
% yticklabels(num2cell(round(Lset),2))
% caxis([2 3])
% hold on
% scatter(outliers6(3,2),outliers6(3,1),225,'MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',2)
% hold off

tabletest = cell(5,5);
combs = nchoosek(1:5,2);
for i = 1:length(combs)
    for j = 1:5
        tabletest{i,j} = {sprintf('outliers6(i,%d)',j)};
    end
    tabletest{i,combs(i,1)} = {sprintf(':')};
    tabletest{i,combs(i,2)} = {sprintf(':')};
end

titles = ["L","Major Axis","dlaser","eccentricity","h"];


%messy code but here is the jist
%foop is a string of permutations of the top6 like
%":,:,outliers6(i,3),outliers6(i,4),outliers6(i,5)"
%this is the first of the permutations to show all contours along all
%dimensions.  ie L vs Major axis, L vs dlaser, L vs ecc, etc

%eval forces the function to call as written ie
%q1(:,:,outliers6(i,3),outliers6(i,4),outliers6(i,5))
%the reshape is there because matlab is an ass about multi dimensions like
%q1(1,:,:,1,1) or q1(1,1,:,:,1) it starts trying to do 3 and 4 dimensional
%plots or something
k = 1;
 for i = 1:6
     tStart = tic;
     figure(i)
     for j = 1:10
     subplot(2,5,j)
     foop = str2mat(sprintf('%s,%s,%s,%s,%s', cell2mat(tabletest{j,1}),cell2mat(tabletest{j,2}),cell2mat(tabletest{j,3}),cell2mat(tabletest{j,4}),cell2mat(tabletest{j,5})));
%      B = reshape(eval(['q1(' foop ')']),pts,pts);
%      A = interp2(B,k,'bicubic');
%      contourf(A,[0:0.5:2.5 2.5:0.1:2.7 2.7:0.02:3],'ShowText','on')
     contourf(reshape(eval(['q1(' foop ')']),pts,pts),[0:0.5:2.5 2.5:0.1:2.7 2.7:0.02:3],'ShowText','on')
     colorbar
%      colormap('jet')
     ylabel(titles(combs(j,1)),'FontSize',16)
     xlabel(titles(combs(j,2)),'FontSize',16)
     str = sprintf('%s vs %s', titles(combs(j,1)),titles(combs(j,2)));
     title(str,'FontSize',16)
     caxis([1.5 3])
     set(gca,'XTick',[], 'YTick', []);
     hold on
     %      scatter(outliers6(i,combs(j,2))*2*k-1,outliers6(i,combs(j,1))*2*k-1,75,'MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',2)
     scatter(outliers6(i,combs(j,2)),outliers6(i,combs(j,1)),75,'MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',2)
     hold off
     end
     toc(tStart)
 end
 
%  i=1;
%  Bb = interpn(q1,1,'bilinear');
%  Ba = reshape(q1(outliers6(i,1),:,:,outliers6(i,4),outliers6(i,5)),pts,pts);
%  Aa = imresize(Ba,2,'bilinear');
%  Ca = interp2(Ba,1,'bicubic');
%  Cb = reshape(Bb(outliers6(i,1)*2*k-1,:,:,outliers6(i,4)*2*k-1,outliers6(i,5)*2*k-1),pts*2*k-1,pts*2*k-1);
%  figure(999)
%  subplot(2,2,1)
%  contourf(Ba,[0:0.5:2.5 2.5:0.1:2.7 2.7:0.02:3],'ShowText','on')
%  colorbar
%  subplot(2,2,2)
%  contourf(Aa,[0:0.5:2.5 2.5:0.1:2.7 2.7:0.02:3],'ShowText','on')
%  colorbar
%   subplot(2,2,3)
%  contourf(Ca,[0:0.5:2.5 2.5:0.1:2.7 2.7:0.02:3],'ShowText','on')
%  colorbar
%  subplot(2,2,4)
%  contourf(Cb,[0:0.5:2.5 2.5:0.1:2.7 2.7:0.02:3],'ShowText','on')
%  colorbar
 
 
% figure(55)
% for i = 1:length(outliers6)
% 
%     xlabel('Detector #','FontSize',20)
%     ylabel('Simulated Intensity','FontSize',20)
%     plot(profiles(out(i),:),'LineWidth',2)
%     hold on
%     title('Intensity Profile (top6)','FontSize',24)
%     xlim([0 320])
%     ylim([0 20000])
% end
%  legend('1st','2nd','3rd','4th','5th','6th')
 
