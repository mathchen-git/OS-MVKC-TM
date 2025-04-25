function NlzData = NormalizeData(normData,num_view,num_samp,data) %%% normalize data

%   NlzData = cell(num_view,1);
  if normData == 1
     for i = 1 :num_view
        for  j = 1:num_samp
            normItem = std(data{i}(j,:));
            if (0 == normItem)
                normItem = eps;
            end
            NlzData{i}(j,:) = (data{i}(j,:) - mean(data{i}(j,:)))/normItem;
        end
     end
  end
end

