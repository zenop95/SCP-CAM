function CDM = read_CDM_XML(xmlcdm_file_path)
% Given an input XML-CDM file, the output is a struct containing its data
% -------------------------------------------------------------------------
% INPUT
% - xmlcdm_file_path: XML CDM File path
% -------------------------------------------------------------------------
% OUTPUT
% -	CDM: struct with fields:
%          - header.contents       : header contents from the CDM file
%          - body.metadata.comments: metadata comments from the CDM file
%          - body.metadata.contents: metadata contents from the CDM file
%          - body.data             : data from the CDM file
% -------------------------------------------------------------------------
% Author:   Marco Felice Montaruli, Politecnico di Milano, 14 December 2020
%           e-mail: marcofelice.montaruli@polimi.it

% Create output structure
CDM = struct();

% Read XML file
xml = xmlread(xmlcdm_file_path);

% Get the main node (<tdm>)
cdm_node = xml.getElementsByTagName('cdm').item(0);

if isempty(cdm_node)
    error('tdm node not found in xml');
end

% Get the <header> node
header_node = cdm_node.getElementsByTagName('header').item(0);
if isempty(header_node)
    error('tdm node not found in xml');
end

% Process <header> content
header_children = header_node.getChildNodes();
comment_index = 1;
for k = 0:header_children.getLength()-1
    item = header_children.item(k);
    item_name = string(item.getNodeName());
    item_text = string(item.getTextContent());
    if item_name == "#text"
        continue;
    end
    if item_name == "COMMENT"
        CDM.header.comments(comment_index) = item_text;
        comment_index = comment_index + 1;
    else
        CDM.header.contents.(item_name) = item_text;
    end
end

% Get <body> node
body_node = cdm_node.getElementsByTagName('body').item(0);
if isempty(body_node)
    error('body node not found in xml');
end

% Process body content
body_children = body_node.getChildNodes();

% Get <relative metadata> node
rel_metadata_node = body_children.getElementsByTagName('relativeMetadataData').item(0);
if isempty(rel_metadata_node)
    error('relativeMetadataData node not found in xml');
end

for k = 0:rel_metadata_node.getLength()-1
    item = rel_metadata_node.item(k);
    item_name = string(item.getNodeName());
    item_text = string(item.getTextContent());
    if item_name == "#text"
        continue;
    end
    if item_name == "relativeStateVector"
        %         rel_vect = str2double(split(item_text));
        %         item_text = rel_vect(~isnan(rel_vect));
        relativeStateVector_node = body_children.getElementsByTagName('relativeStateVector').item(0);
        
        for kk = 0:relativeStateVector_node.getLength()-1
            item_rel = relativeStateVector_node.item(kk);
            item_name_rel = string(item_rel.getNodeName());
            item_text_rel = string(item_rel.getTextContent());
            if item_name_rel == "#text"
                continue;
            end
            
            CDM.body.relativeMetadataData.contents.relativeStateVector.contents.(item_name_rel) = item_text_rel;
        end
        
    else
        CDM.body.relativeMetadataData.contents.(item_name) = item_text;
    end
end


% Get <relative metadata> node
obj_node1 = body_children.getElementsByTagName('segment').item(0);
obj_struct1 = object_info(obj_node1);


obj_node2 = body_children.getElementsByTagName('segment').item(1);
obj_struct2 = object_info(obj_node2);


CDM.body.object1 = obj_struct1;
CDM.body.object2 = obj_struct2;


end

function obj_struct = object_info(obj_node)

for k = 0:obj_node.getLength()-1
    item = obj_node.item(k);
    item_name = string(item.getNodeName());
    if item_name == "#text"
        continue;
    end
    if item_name == "metadata"
        
        segm_metedata_node = obj_node.getElementsByTagName('metadata').item(0);
        
        for kk = 0:segm_metedata_node.getLength()-1
            item_met = segm_metedata_node.item(kk);
            item_name_met = string(item_met.getNodeName());
            item_text_met = string(item_met.getTextContent());
            if item_name_met == "#text"
                continue;
            end
            
            obj_struct.metadata.contents.(item_name_met) = item_text_met;
        end
        
    elseif item_name == "data"
        data_node1 = obj_node.getElementsByTagName('data').item(0);
        for kk = 0:data_node1.getLength()-1
            item_data = data_node1.item(kk);
            item_name_data = string(item_data.getNodeName());
            
            if item_name_data == "#text"
                continue;
            end
            
            if item_name_data == "additionalParameters"
                
                iter_node = obj_node.getElementsByTagName('additionalParameters').item(0);
                
                for kh = 0:iter_node.getLength()-1
                    
                    
                    item_iter = iter_node.item(kh);
                    item_name_iter = string(item_iter.getNodeName());
                    item_text_iter = string(item_iter.getTextContent());
                    if item_name_iter == "#text"
                        continue;
                    end
                    
                    obj_struct.data.additionalParameters.contents.(item_name_iter) = item_text_iter;
                end
            elseif item_name_data == "stateVector"
                stateVector_node = obj_node.getElementsByTagName('stateVector').item(0);
                
                for kh = 0:stateVector_node.getLength()-1
                    item_iter = stateVector_node.item(kh);
                    item_name_iter = string(item_iter.getNodeName());
                    item_text_iter = string(item_iter.getTextContent());
                    if item_name_iter == "#text"
                        continue;
                    end
                    
                    obj_struct.data.stateVector.contents.(item_name_iter) = item_text_iter;
                end
            elseif item_name_data == "covarianceMatrix"
                covarianceMatrix_node = obj_node.getElementsByTagName('covarianceMatrix').item(0);
                
                for kh = 0:covarianceMatrix_node.getLength()-1
                    item_iter = covarianceMatrix_node.item(kh);
                    item_name_iter = string(item_iter.getNodeName());
                    item_text_iter = string(item_iter.getTextContent());
                    if item_name_iter == "#text"
                        continue;
                    end
                    
                    obj_struct.data.covarianceMatrix.contents.(item_name_iter) = item_text_iter;
                end
            end
        end
        
    end
end


end