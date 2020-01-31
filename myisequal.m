function [res,reldevi,absdevi] = myisequal(A,B,tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        myisequal.m
% function [res,reldevi,absdevi] = myisequal(A,B,tol)
% compares two arrays with a tolerance
% res = logical output 1=equal = 2=non-equal
% reldevi = relative difference
% absdevi = absolute difference
% A = first array
% B = second array
% array dimension less than 4
% tol = tolerance with which to compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = 'max_relative';

if size(A)~=size(B)
    res = 0;
    reldevi = 1;
    absdevi = 1;
else
    
    switch method
        case 'original_dmitri'
            mA = max(max(max(max(abs(A)))));
            mB = max(max(max(max(abs(B)))));
            absdevi = max(max(max(max(abs(A-B)))));
            reldevi = 2*absdevi/(mA+mB);
        case 'max_relative'
            [reldevi ii] = max( 2.*abs(A-B) ./ (abs(A) + abs(B)) );
            absdevi = abs(A(ii)-B(ii));
        case 'max_absolute'
            [absdevi ii] = max(abs(A-B));
            reldevi = 2*absdevi/( abs(A(ii)) + abs(B(ii)) );
    end
    
    if reldevi<tol
        res = 1;
    else
        res = 0;
    end
    
end