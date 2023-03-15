classdef fortran_ascii < handle
%FORTRAN_ASCII FORTRAN ASCII format string
    properties (Dependent)
        %MANTISSA_SCALE Mantissa scale spec
        mantissa_scale char
        %DEC_PLACES Decimal places
        dec_places(1,1) int32
        %EXPONENT_CHAR Exponent character
        exponent_char(1,1) char
        %REPITITIONS Number of values per line
        repititions(1,1) int32
        %FIELD_WIDTH Width of the field
        field_width(1,1) int32
    end
    properties (SetAccess = immutable)
        %STR_FMT String format descriptor
        str_fmt
    end
    methods
        % Getters and setters
        function mantissa_scale = get.mantissa_scale(this)
            mantissa_scale = this.str_fmt(1:strfind(this.str_fmt,',')-1);
        end
        function dec_places = get.dec_places(this)
            dec_places = str2num(this.str_fmt(strfind(this.str_fmt,'.')+1:end));
        end
        function exponent_char = get.exponent_char(this)
            FMT_substr = this.str_fmt(strfind(this.str_fmt,',')+1 ...
                                     :strfind(this.str_fmt,'.')-1);
            exponent_char = FMT_substr(isletter(FMT_substr));
        end
        function repititions = get.repititions(this)
            FMT_substr = this.str_fmt(strfind(this.str_fmt,',')+1 ...
                                     :strfind(this.str_fmt,'.')-1);
            repititions = str2num(FMT_substr(1:find(isletter(FMT_substr))-1));
        end
        function field_width = get.field_width(this)
            FMT_substr = this.str_fmt(strfind(this.str_fmt,',')+1 ...
                                     :strfind(this.str_fmt,'.')-1);
            field_width = str2num(FMT_substr(find(isletter(FMT_substr))+1:end));
        end
        % Constructor
        function this = fortran_ascii(str_fmt)
        %FORTRAN_ASCII Class constructor
            
            FMT_substr = str_fmt(strfind(str_fmt,',')+1 ...
                                :strfind(str_fmt,'.')-1);
            ex_posn = isletter(FMT_substr);
            if ~any(ex_posn)
                error('Unsupported FORTRAN format');
            end
            this.str_fmt = str_fmt;
            
        end
    end
end

%#ok<*ST2NM>

