function [ colorsout badflag ] = bs2cs( bases, LigationSequence, initialbase, finalbase )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% LigationSequence is the order in which the ligations are obtained. For
% example, if we do primer N ligation 1, then primer N-1 ligation 1, Primer
% N-1 ligation 2, primer N-2 ligation 2, primer N-3 ligation 2, and primer
% N-4 ligation 2, in that order, LigationSequence is [2, 1, 6, 5, 4, 3].
badflag=0;

seq=[initialbase,bases,finalbase]; %initialbase is the last base of the primer

colors=[];
for i=1:(length(seq)-1)
    switch seq(i)
        case 'A'
            switch seq(i+1)
                case 'A'
                    colors(i)=1;
                case 'C'
                    colors(i)=2;
                case 'G'
                    colors(i)=3;
                case 'T'
                    colors(i)=4;
                case 'N'
                    colors(i)=0;
                    badflag=1;
            end
        case 'C'
            switch seq(i+1)
                case 'A'
                    colors(i)=2;
                case 'C'
                    colors(i)=1;
                case 'G'
                    colors(i)=4;
                case 'T'
                    colors(i)=3;
                case 'N'
                    colors(i)=0;
                    badflag=1;
            end
        case 'G'
            switch seq(i+1)
                case 'A'
                    colors(i)=3;
                case 'C'
                    colors(i)=4;
                case 'G'
                    colors(i)=1;
                case 'T'
                    colors(i)=2;
                case 'N'
                    colors(i)=0;
                    badflag=1;
            end
        case 'T'
            switch seq(i+1)
                case 'A'
                    colors(i)=4;
                case 'C'
                    colors(i)=3;
                case 'G'
                    colors(i)=2;
                case 'T'
                    colors(i)=1;
                case 'N'
                    colors(i)=0;
                    badflag=1;
            end
        case 'N'
            colors(i)=0;
            badflag=0;
    end
           
           

end
colorsout=colors(LigationSequence);     
end
