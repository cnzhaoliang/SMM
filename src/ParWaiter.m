classdef ParWaiter < handle

    properties
        count = 0;
        totality = 0;
        h;
    end

    properties (Access = private)
        old = 0;
        step = 0.01;
    end

    methods

        function obj = ParWaiter(totality, step)

            obj.totality = totality;
            obj.step = step / 100;
            obj.h = waitbar(0, sprintf('计算进度: %.2f%%\n', 0));
            tic;

        end

        function obj = update(obj)

            obj.count = obj.count + 1;

            if (obj.count - obj.old) / obj.totality >= obj.step
                obj.old = obj.count;
                waitbar(obj.count / obj.totality, obj.h, sprintf('计算进度: %.2f%%\n', obj.count / obj.totality * 100));
            end

        end

        function obj = close(obj)

            fprintf('计算耗时: %.2fs\n', toc);
            obj.count = 0;

            if isvalid(obj.h)
                close(obj.h);
            end

        end

    end

end
