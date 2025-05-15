function [missedNodes] = clearNodes(nodes,missedNodes)
missingNode = 0;
cont = length(missedNodes) + 1;
for i = 1:length(nodes(:))
    if nodes(i) ~= missingNode
        while nodes(i) ~= missingNode
            missedNodes(cont) = missingNode;
            cont = cont + 1;
            missingNode = missingNode + 1;
        end
        missingNode = missingNode + 1;
    else
        missingNode = missingNode + 1;
    end
end
end

