classdef PhaseVisualize < handle
%PHASEVISUALIZE  MATLAB wrapper for the phase_visualize C shared library.
%
%   pv = PhaseVisualize();                    % auto-locate library
%   pv = PhaseVisualize(libDir);              % explicit folder
%   pv = PhaseVisualize(libDir, headerPath);  % explicit folder + header
%
%   The library file must be named:
%       phase_visualize.dll   (Windows)
%       libphase_visualize.so (Linux)
%       libphase_visualize.dylib (macOS)
%
%   All array inputs/outputs are double, column-major (MATLAB native).
%
%   Methods:
%       unwrapped   = pv.unwrapGoldstein(wrapped_phase)
%       [c, s]      = pv.sincos(phase)
%       [gr, gc, m] = pv.gradient(phase)
%       coh         = pv.coherence(phase, window)
%       rgb         = pv.toRGB(phase, coherence)   % rows x cols x 3
%
%   Call pv.unload() (or 'clear pv') to release the library.

    properties (Constant, Access = private)
        LIB_NAME    = 'phase_visualize';
        HEADER_NAME = 'phase_visualize_matlab.h';
    end

    properties (Access = private)
        loadedHere logical = false;
    end

    methods
        function obj = PhaseVisualize(libDir, headerPath)
            if nargin < 1 || isempty(libDir)
                libDir = fileparts(mfilename('fullpath'));
            end
            if nargin < 2 || isempty(headerPath)
                headerPath = fullfile(libDir, obj.HEADER_NAME);
            end

            if ~libisloaded(obj.LIB_NAME)
                libFile = obj.libraryFile(libDir);
                if ~isfile(libFile)
                    error('PhaseVisualize:LibraryNotFound', ...
                          'Could not find shared library: %s', libFile);
                end
                if ~isfile(headerPath)
                    error('PhaseVisualize:HeaderNotFound', ...
                          'Could not find header file: %s', headerPath);
                end

                % Make sure the library's directory is on the OS search
                % path so dependent DLLs (e.g. OpenMP runtime) load.
                addpath(libDir);

                warning('off', 'MATLAB:loadlibrary:TypeNotFound');
                cleanup = onCleanup(@() warning('on', 'MATLAB:loadlibrary:TypeNotFound'));
                loadlibrary(libFile, headerPath, 'alias', obj.LIB_NAME);
                obj.loadedHere = true;
            end
        end

        function delete(obj)
            if obj.loadedHere && libisloaded(obj.LIB_NAME)
                unloadlibrary(obj.LIB_NAME);
            end
        end

        function unload(obj)
            if libisloaded(obj.LIB_NAME)
                unloadlibrary(obj.LIB_NAME);
            end
            obj.loadedHere = false;
        end

        % ------------------------------------------------------------ unwrap
        function unwrapped = unwrapGoldstein(obj, wrapped, maxBox)
            % unwrapGoldstein(wrapped)        - default search box (9)
            % unwrapGoldstein(wrapped, maxBox)- override Goldstein box cap
            %
            % maxBox is the dominant performance knob. Smaller values
            % (e.g. 5-9) are dramatically faster; larger values may resolve
            % more residue pairs at significant runtime cost. Pass 0 or omit
            % to use the built-in default.
            obj.checkMatrix(wrapped, 'wrapped');
            wrapped = double(wrapped);
            [rows, cols] = size(wrapped);
            if nargin < 3 || isempty(maxBox)
                maxBox = 0;
            end

            unwrapped = zeros(rows, cols);
            [status, ~, unwrapped] = calllib(obj.LIB_NAME, ...
                'phase_unwrap_goldstein_ex', wrapped, unwrapped, ...
                int32(rows), int32(cols), int32(maxBox));
            obj.checkStatus(status, 'phase_unwrap_goldstein_ex');
        end

        % ----------------------------------------------------------- sincos
        function [cosOut, sinOut] = sincos(obj, phase)
            phase = double(phase);
            sz = size(phase);
            n  = numel(phase);

            cosOut = zeros(sz);
            sinOut = zeros(sz);
            [status, ~, cosOut, sinOut] = calllib(obj.LIB_NAME, ...
                'phase_sincos_decompose', phase(:), cosOut(:), sinOut(:), ...
                int32(n));
            obj.checkStatus(status, 'phase_sincos_decompose');
            cosOut = reshape(cosOut, sz);
            sinOut = reshape(sinOut, sz);
        end

        % --------------------------------------------------------- gradient
        function [gradRow, gradCol, fringeMag] = gradient(obj, phase)
            obj.checkMatrix(phase, 'phase');
            phase = double(phase);
            [rows, cols] = size(phase);

            gradRow   = zeros(rows, cols);
            gradCol   = zeros(rows, cols);
            fringeMag = zeros(rows, cols);
            [status, ~, gradRow, gradCol, fringeMag] = calllib(obj.LIB_NAME, ...
                'phase_gradient', phase, gradRow, gradCol, fringeMag, ...
                int32(rows), int32(cols));
            obj.checkStatus(status, 'phase_gradient');
        end

        % --------------------------------------------------------- coherence
        function coh = coherence(obj, phase, window)
            obj.checkMatrix(phase, 'phase');
            if nargin < 3 || isempty(window), window = 3; end
            phase = double(phase);
            [rows, cols] = size(phase);

            coh = zeros(rows, cols);
            [status, ~, coh] = calllib(obj.LIB_NAME, ...
                'phase_coherence', phase, coh, ...
                int32(rows), int32(cols), int32(window));
            obj.checkStatus(status, 'phase_coherence');
        end

        % ------------------------------------------------------------ toRGB
        function rgb = toRGB(obj, phase, coherenceMap)
            obj.checkMatrix(phase, 'phase');
            phase = double(phase);
            [rows, cols] = size(phase);

            if nargin < 3 || isempty(coherenceMap)
                coherenceArg = libpointer('doublePtr');  % NULL
            else
                if ~isequal(size(coherenceMap), [rows, cols])
                    error('PhaseVisualize:SizeMismatch', ...
                          'coherence must match phase size [%d %d].', rows, cols);
                end
                coherenceArg = double(coherenceMap);
            end

            rOut = zeros(rows, cols);
            gOut = zeros(rows, cols);
            bOut = zeros(rows, cols);
            [status, ~, ~, rOut, gOut, bOut] = calllib(obj.LIB_NAME, ...
                'phase_to_rgb_hsv', phase, coherenceArg, rOut, gOut, bOut, ...
                int32(rows), int32(cols));
            obj.checkStatus(status, 'phase_to_rgb_hsv');
            rgb = cat(3, rOut, gOut, bOut);
        end
    end

    methods (Access = private, Static)
        function file = libraryFile(libDir)
            base = PhaseVisualize.LIB_NAME;
            if ispc
                file = fullfile(libDir, [base '.dll']);
            elseif ismac
                file = fullfile(libDir, ['lib' base '.dylib']);
            else
                file = fullfile(libDir, ['lib' base '.so']);
            end
        end

        function checkMatrix(A, name)
            if ~ismatrix(A) || ~isnumeric(A) || isempty(A)
                error('PhaseVisualize:BadInput', ...
                      '''%s'' must be a non-empty 2-D numeric matrix.', name);
            end
        end

        function checkStatus(status, fn)
            if status ~= 0
                error('PhaseVisualize:CallFailed', ...
                      '%s returned status %d.', fn, status);
            end
        end
    end
end
