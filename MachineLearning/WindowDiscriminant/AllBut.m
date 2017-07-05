function oc = AllBut(Array,But)
    oc = Array;
    oc(But,:) = [];
end