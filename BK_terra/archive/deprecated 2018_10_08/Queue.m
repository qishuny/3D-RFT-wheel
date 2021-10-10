%%
classdef Queue < handle
    
    properties
        list = {};
    end

    %%
    methods
        % constructor synchronize with hmapCellContext at the beginning
        function obj = Queue()

        end
        
        function add(obj,element)
            
            if ~isempty(obj.list) % if this is not the first element
                % add it to the list if they are not the same
                for i=1:length(obj.list) % for it is not the same with my list
                    if element.isSame(obj.list{i})
                        return
                    end
                end
                % if it is not the same add it to the list
                obj.list{end+1} = element;
                
            elseif isempty(obj.list) %if this is first element just add it
                obj.list{end+1} = element;
            end
        end
        
        function clear(obj)
            % clear all list
            obj.list = {};
        end
        

    end
end