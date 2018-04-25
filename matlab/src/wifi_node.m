classdef wifi_node < matlab.mixin.SetGet     
    properties         
        is_CHAIN = 0
        is_active = 0
        pred_id = 0
        succ_id = 0
    end

    methods
        % Constructor
        function obj = wifi_node(c, a)
        	obj.is_CHAIN = c;
            obj.is_active = a;
            %obj.pred_id = p; % id of the predecessor 
            %obj.succ_id = s; % id of the successor
        end
    end

end