classdef (Abstract) Handle < handle
%HANDLE  Base class of the meshio++ native handles.
%   A handle owns one native object of the meshio++ C library. The object is
%   freed when the last MATLAB/Octave reference goes away (DELETE), or earlier
%   by an explicit RELEASE. Handles are created by mio.* functions, never
%   directly: the ID is an opaque key into a registry inside the MEX, and a
%   released, stale (loaded from a .mat file) or forged ID is rejected rather
%   than dereferenced.

    properties (SetAccess = private, Hidden = true)
        id = uint64(0);
    end

    methods
        function obj = Handle(id)
            if nargin > 0
                obj.id = uint64(id);
            end
        end

        function released = release(obj)
            %RELEASE  Free the native object now. Idempotent.
            released = false;
            for k = 1:numel(obj)
                released = meshioplusplus_mex('__release__', obj(k)) || released;
            end
        end

        function tf = is_open(obj)
            %IS_OPEN  True while the native object is live.
            tf = meshioplusplus_mex('__is_open__', obj);
        end

        function delete(obj)
            if obj.id ~= 0
                try
                    meshioplusplus_mex('__release__', obj);
                catch
                    % The MEX may already be gone during interpreter shutdown;
                    % its exit hook has freed everything by then.
                end
            end
        end
    end
end
