function idx = logicalUnion(Vector1,Vector2)
%returns logical indexing for any time an element of vector 2 appears in
%vector 1
    idx = any(bsxfun(@eq,Vector1,Vector2),2);

end