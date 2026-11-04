function v = od_pad3(x)
%OD_PAD3  A 2-vector padded with 1 to length 3 (R's c(size, 1)[1:3]).

    v = [double(x(:).'), 1];
    v = v(1:3);
end
