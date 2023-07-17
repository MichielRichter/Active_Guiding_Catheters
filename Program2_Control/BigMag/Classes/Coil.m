classdef Coil
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        transform; %4x4 transformation matrix local w.r.t global
        transform_init;
        id;
    end
    
    methods 
        %Constructor with arguments:
        %Starting orientation R
        %Starting position p
        %Coil identification number id
        function obj = Coil(R,p,id)
            obj.transform_init = [R p; zeros(1,3) 1];
            obj.transform = [R p; zeros(1,3) 1];
            obj.id = id;
        end
        
        %Global to local position of the field
        function localvector = globalToLocalPosition(obj, globalvector)
            a = [globalvector; 1];
            b = obj.transform\a;
            localvector = b(1:3,1);
        end
        %Local to global position of the field
        function globalvector = localToGlobalField(obj, localvector)
            R = obj.transform(1:3,1:3);
            globalvector = R * localvector;
        end
        
        %Global to local position of the field
        function globaltensor = localToGlobalTensor(obj, localtensor)
            R = obj.transform(1:3,1:3);
            globaltensor = R * localtensor * R';
        end

    end
end

