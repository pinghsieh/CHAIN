classdef wifi_node < matlab.mixin.SetGet     
    properties         
        is_CHAIN = 0
        is_active = 1
        %pred_id
        %succ_id
    end

    methods
        % Constructor
        function obj = wifi_node(c, a)
        	obj.is_CHAIN = c;
            obj.is_active = a;
            %obj.pred_id = p;
            %obj.succ_id = s;
        end
    end

end