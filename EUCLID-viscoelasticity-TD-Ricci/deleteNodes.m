function reference = deleteNodes(reference,missedNodes)
%   var is the missed node index
var = 0;
for i = 1:length(missedNodes)
    var = find(reference(:,1) == missedNodes(i));
    reference(var,:) = [];
end
end

