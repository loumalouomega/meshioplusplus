function v = od_rep3(x)
%OD_REP3  A scalar repeated three times (R's rep(size, 3)).

    v = repmat(double(x(:).'), 1, 3);
end
