function customArea(X,Y,labels,alpha,colors)
% PLOTS A custom area stack
if nargin<5
     colors=[[38, 70, 83]/255,	          	
                [42, 157, 143]/255,	          	
                [233, 196, 106]/255,	          	
                [244, 162, 97]/255,	          	
                [231, 111, 81]/255,	          
                [144, 44, 20]/255    ];
end

if nargin<4
    alpha=0.9;
end
   

    s=size(Y);
    previous = zeros(s(1),1)';
    for k=1:s(2)
        jbfill(X,previous,previous+squeeze(Y(:,k))',colors(mod(k,length(colors)),:),[0.20 0.20 0.20],1,alpha,labels(k));
        %previous = max(squeeze(Y(:,k))'+previous,zeros(s(1),1)');
        previous = squeeze(Y(:,k))' + previous;
    end
    legend()
end