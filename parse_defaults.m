function o = parse_defaults(default, varargin)

% % Define function properties and default values:
%   default = ...
%   struct( ...
%           'property_1' , 100, ...
%           'property_2' , 200 ...
%         );

% Handle variable {'property', value} list:
  if length(varargin) == 1 && isstruct(varargin{1})
    o = varargin{1};
  else
    o = struct(varargin{:});
  end;
  given_elements = fieldnames(o);
  valid_elements = fieldnames(default);
  for n=1:length(valid_elements)
    if ~ isfield(o, valid_elements{n})
      o = setfield(o, valid_elements{n}, ...
                   getfield(default, valid_elements{n}));
    end;
  end;
  for n=1:length(given_elements)
    if ~ any(strcmp(given_elements{n}, valid_elements))
      error(['invalid property: "' given_elements{n} '"']);
    end;
  end;

end

%Felipe Orduña <felipe.orduna@ccadet.unam.mx>
%[http://www.academicos.ccadet.unam.mx/felipe.orduna]
