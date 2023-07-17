classdef IroncoreCoil < Coil
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
     properties
    end
    
    methods 
        function fieldB = getUnitFieldatLocation(obj,globalposition)
            %Get local position
            localposition = obj.globalToLocalPosition(globalposition);
            %Get local field
            localB = getMultipoleField_Ironcore(localposition);
            %Return global field
            fieldB = obj.localToGlobalField(localB);
        end
        
        function gradB = getUnitGradientatLocation(obj,globalposition)
            %Get local position
            localposition = obj.globalToLocalPosition(globalposition);
            %Get local field
            localgradB = getMultipoleGradient_Ironcore(localposition);
            %Return global field
            gradB = obj.localToGlobalTensor(localgradB);
        end
    end
end

