classdef AircoreCoil < Coil
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods 
        function fieldB = getUnitFieldatLocation(obj,globalposition)
            %Get local position
            localposition = obj.globalToLocalPosition(globalposition);
            %Get local field
            localB = getMultipoleField_Aircore(localposition);
            %Return global field
            fieldB = obj.localToGlobalField(localB);
        end
        
        function gradB = getUnitGradientatLocation(obj,globalposition)
            %Get local position
            localposition = obj.globalToLocalPosition(globalposition);
            %Get local field
            localgradB = getMultipoleGradient_Aircore(localposition);
            %Return global field
            gradB = obj.localToGlobalTensor(localgradB);
        end
    end
    
end

