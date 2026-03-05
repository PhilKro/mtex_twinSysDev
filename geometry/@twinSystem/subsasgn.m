function tS = subsasgn(tS,s,value)
% overloads subsasgn

if ~isa(tS,'twinSystem') && ~isempty(value)
  tS = value;

  tS.eta1.x = [];
  tS.eta1.y = [];
  tS.eta1.z = [];

  tS.eta2.x = [];
  tS.eta2.y = [];
  tS.eta2.z = [];

  tS.k1.x = [];
  tS.k1.y = [];
  tS.k1.z = [];

  tS.k2.x = [];
  tS.k2.y = [];
  tS.k2.z = [];

  tS.rotAxis.x = [];
  tS.rotAxis.y = [];
  tS.rotAxis.z = [];

  tS.CRSS = [];
  tS.twinType = [];
  tS.parent = [];
end

switch s(1).type
  
  case '()'
      
    if numel(s)>1, value =  builtin('subsasgn',subsref(tS,s(1)),s(2:end),value); end
      
    if isempty(value)
      tS.eta1 = subsasgn(tS.eta1,s(1),[]);
      tS.eta2 = subsasgn(tS.eta2,s(1),[]);
      
      tS.k1 = subsasgn(tS.k1,s(1),[]);
      tS.k2 = subsasgn(tS.k2,s(1),[]);

      tS.rotAxis = subsasgn(tS.rotAxis,s(1),[]);
     
      tS.CRSS = subsasgn(tS.CRSS,s(1),[]);
      tS.twinType = subsasgn(tS.twinType,s(1),[]);
      tS.parent = subsasgn(tS.parent,s(1),[]);
    else      
      tS.eta1 = subsasgn(tS.eta1,s(1),value.eta1);
      tS.eta2 = subsasgn(tS.eta2,s(1),value.eta2);
      
      tS.k1 = subsasgn(tS.k1,s(1),value.k1);
      tS.k2 = subsasgn(tS.k2,s(1),value.k2);

      tS.rotAxis = subsasgn(tS.rotAxis,s(1),value.rotAxis);
     
      tS.CRSS = subsasgn(tS.CRSS,s(1),value.CRSS);
      tS.twinType = subsasgn(tS.twinType,s(1),value.twinType);
      tS.parent = subsasgn(tS.parent,s(1),value.parent);
    end
  otherwise
    
    tS =  builtin('subsasgn',tS,s,value);    

end

end