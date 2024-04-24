classdef ExpParam
%     传递一些实验中需要设定的参数
%     这里统一按照组进行设定，防止后面突然出现一堆参数
    
    properties
        selectArea % 设定 range focus Area 坐标
        cfarAzimuth
        fileName
        ldrFileName
        rawFileName
    end
    
    methods
        function obj = ExpParam(selectArea, cfarAzimuth, fileName)
            %EXPPARAM 构造此类的实例
            %   此处显示详细说明
            obj.selectArea = selectArea;
            obj.cfarAzimuth = cfarAzimuth;
            obj.fileName = fileName;
            obj.ldrFileName = sprintf('./%s/data/%s.000.ldr', fileName, fileName);
            obj.rawFileName = sprintf('./%s/data/%s.000.raw', fileName, fileName);
        end
    end
end

