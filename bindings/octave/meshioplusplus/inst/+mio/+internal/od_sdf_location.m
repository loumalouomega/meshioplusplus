function code = od_sdf_location(location)
%OD_SDF_LOCATION  The SDF sample location code: 'corner'/'point' -> 0,
%   'center'/'centre'/'cell' -> 1. These are the SDF locations, NOT the
%   point/cell/field data locations.

    names = {'corner', 'point', 'center', 'centre', 'cell'};
    codes = [0, 0, 1, 1, 1];
    hit = [];
    if ischar(location)
        hit = find(strcmp(location, names), 1);
    end
    if isempty(hit)
        error('meshioplusplus:argument', 'unknown location ''%s''', ...
              mio.internal.str(location, 'location'));
    end
    code = codes(hit);
end
